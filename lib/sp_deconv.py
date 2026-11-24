"""
dsh-spatial — 反卷积（DESIGN.md §3 `deconv` 步骤 + §4 group=反卷积）

思路（文献）：
  * 非负最小二乘 / 阻尼加权最小二乘的 spot 反卷积（Tsoucas 2019 DWLS；SPOTlight 2021；
    Sturm 2019 对 Visium 反卷积方法的评测）。
  * 参考谱来自单细胞 h5ad（ctx["ref"] 懒加载）：每类细胞的平均表达谱（signature 矩阵）。

本模块只依赖 numpy/scipy/matplotlib + sp_common，不依赖 scanpy/squidpy。
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

import sp_common as sc

_GENERIC_TYPES = ("cell_type", "celltype", "cell_type_l1", "annotation", "clusters",
                  "leiden", "seurat_clusters", "louvain")


# ---------------------------------------------------------------------------
# ctx 小工具
# ---------------------------------------------------------------------------


def _p(params, key, default=None):
    value = (params or {}).get(key, default)
    return default if value in (None, "", "null") else value


def _p_int(params, key, default):
    try:
        return int(float(_p(params, key, default)))
    except (TypeError, ValueError):
        return int(default)


def _p_float(params, key, default):
    try:
        return float(_p(params, key, default))
    except (TypeError, ValueError):
        return float(default)


def _X(ctx):
    X = ctx.get("X")
    if X is None:
        raise sc.prerequisite("数据集没有表达矩阵")
    return np.asarray(X, dtype=np.float32)


def _genes(ctx):
    return [str(item) for item in (ctx.get("genes") or [])]


def _spots(ctx):
    return [str(item) for item in (ctx.get("spots") or [])]


def _info(ctx) -> pd.DataFrame:
    frame = ctx.get("sample_info")
    return frame if isinstance(frame, pd.DataFrame) else pd.DataFrame(index=pd.Index(_spots(ctx), name="spot"))


def _outdir(ctx, outdir) -> str:
    return sc.ensure_dir(outdir or ctx.get("outdir") or ".")


def _label_array(ctx, key=""):
    for candidate in ([key] if key else []) + ["domains", "clusters", ctx.get("cluster_key") or ""]:
        if not candidate:
            continue
        frame = _info(ctx)
        if candidate in frame.columns:
            return np.asarray([str(v) for v in frame[candidate].tolist()], dtype=object)
    labels = ctx.get("cluster")
    return np.asarray([str(v) for v in labels], dtype=object) if labels is not None else None


# ---------------------------------------------------------------------------
# 参考谱
# ---------------------------------------------------------------------------


def _ref_from_path(path: str, cell_type_col: str = "") -> dict:
    """直接读一个 h5ad 作为参考（步骤参数 ref 传路径时用，避免必须先 ref_attach）。"""
    import os

    if not path:
        raise sc.prerequisite("请提供参考 h5ad 路径")
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(path):
        raise sc.prerequisite(f"参考文件不存在：{path}")
    try:
        import anndata

        adata = anndata.read_h5ad(path)
    except Exception as error:  # noqa: BLE001
        raise sc.prerequisite(f"读取参考 h5ad 失败：{error}")
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return {"path": path, "X": np.asarray(X, dtype=np.float32),
            "genes": [str(g) for g in adata.var_names], "obs": adata.obs,
            "cell_type_col": cell_type_col or ""}


def _ref_columns(ref) -> tuple[np.ndarray, list[str], list[str]]:
    """从 ctx["ref"] 取 (X, genes, cell_types)。ref 可能是 LazyRef（dict 子类）。"""
    if ref is None:
        raise sc.prerequisite("还没有指定单细胞参考：先在「反卷积」页选一个 h5ad，或传 ref 参数")
    try:
        ref = ref.materialize() if hasattr(ref, "materialize") else ref
    except Exception as error:  # noqa: BLE001
        raise sc.prerequisite(f"读取参考失败：{error}")
    X = ref.get("X") if hasattr(ref, "get") else getattr(ref, "X", None)
    genes = ref.get("genes") if hasattr(ref, "get") else getattr(ref, "genes", None)
    if X is None or genes is None:
        err = (ref.get("error") if hasattr(ref, "get") else None) or "参考里没有表达矩阵"
        raise sc.prerequisite(f"参考 h5ad 不可用：{err}")
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=np.float32)
    genes = [str(g) for g in genes]
    obs = ref.get("obs") if hasattr(ref, "get") else None
    types = None
    if obs is not None and len(obs) == X.shape[0]:
        for col in [ref.get("cell_type_col")] + list(_GENERIC_TYPES):
            if col and str(col) in obs.columns:
                types = [str(v) for v in obs[str(col)].astype(str).tolist()]
                break
    if types is None:
        types = ["all"] * X.shape[0]
    return X, genes, types


def _signature(ref, min_genes: int = 50, cell_type_col: str = "") -> tuple[np.ndarray, list[str], list[str], int]:
    """返回 (signature[genes × types], genes, cell_types, 判定方式)。"""
    X, genes, types = _ref_columns(ref)
    uniq = [t for t in dict.fromkeys(types)]
    if len(uniq) <= 1 and cell_type_col:
        raise sc.prerequisite(f"参考里找不到细胞类型列「{cell_type_col}」：可用列见 /refs 的返回")
    index = {gene: i for i, gene in enumerate(genes)}
    mat = np.zeros((len(genes), len(uniq)), dtype=np.float32)
    for j, name in enumerate(uniq):
        sel = np.array([t == name for t in types])
        if not sel.any():
            continue
        mat[:, j] = X[sel].mean(axis=0)
    return mat, genes, uniq, len(uniq)


def _common_genes(ctx, sig_genes):
    genes = _genes(ctx)
    sig_index = {gene: i for i, gene in enumerate(sig_genes)}
    keep = [(i, sig_index[gene]) for i, gene in enumerate(genes) if gene in sig_index]
    if len(keep) < 10:
        raise sc.prerequisite(
            f"参考与数据集的基因交集只有 {len(keep)} 个：请确认参考与数据集的基因命名一致（symbol/Ensembl）")
    return keep


def _fit_nnls(spot_vec: np.ndarray, sig: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    """非负最小二乘（带可选权重）。"""
    from scipy.optimize import nnls

    a = sig if weights is None else sig * weights[:, None]
    b = spot_vec if weights is None else spot_vec * weights
    try:
        coef, _res = nnls(np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64))
    except Exception:  # noqa: BLE001 - 退化时退回最小二乘并截断负值
        coef = np.linalg.lstsq(np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64), rcond=None)[0]
        coef = np.clip(coef, 0, None)
    total = float(coef.sum())
    if total <= 0:
        return np.full(coef.shape, 1.0 / max(len(coef), 1))
    return coef / total


def _deconvolve(ctx, params) -> dict:
    X = _X(ctx)
    method = str(_p(params, "method", "nnls")).lower()
    min_genes = _p_int(params, "min_genes", 50)
    ref = ctx.get("ref")
    ref_path = str(_p(params, "ref", "") or "")
    if ref_path and (ref is None or not (hasattr(ref, "get") and ref.get("X") is not None)):
        ref = _ref_from_path(ref_path, str(_p(params, "cell_type_col", "") or ""))
    sig, sig_genes, types, n_types = _signature(ref, min_genes=min_genes,
                                               cell_type_col=str(_p(params, "cell_type_col", "") or ""))
    pairs = _common_genes(ctx, sig_genes)
    data_idx = np.array([p[0] for p in pairs], dtype=int)
    sig_idx = np.array([p[1] for p in pairs], dtype=int)
    if len(pairs) < min_genes:
        raise sc.prerequisite(f"参考与数据集的共同基因 {len(pairs)} 个 < min_genes={min_genes}：请降低阈值或换参考")
    sub = np.asarray(X[:, data_idx], dtype=np.float64)
    sig_mat = np.asarray(sig[sig_idx, :], dtype=np.float64)
    # 归一化 signature 列（每类细胞谱缩放到同一量级），避免库大小差异主导比例
    col_sum = sig_mat.sum(axis=0)
    col_sum[col_sum <= 0] = 1.0
    sig_mat = sig_mat / col_sum[None, :]
    weights = None
    if method in ("dwls", "wls", "weighted"):
        # 阻尼权重：每个基因在 signature 里的离散程度（Tsoucas 2019 做法）
        spread = sig_mat.std(axis=1) + 1e-3
        weights = np.sqrt(1.0 / spread)
        weights = weights / np.mean(weights)
    props = np.zeros((sub.shape[0], sig_mat.shape[1]), dtype=np.float32)
    for i in range(sub.shape[0]):
        props[i] = _fit_nnls(sub[i], sig_mat, weights).astype(np.float32)
    return {"proportions": props, "cell_types": list(types), "method": method,
            "n_genes": int(len(pairs)), "genes": [genes for genes in _genes(ctx) if genes],
            "n_types": int(n_types)}


# ---------------------------------------------------------------------------
# 图
# ---------------------------------------------------------------------------


def _fig_composition(ctx, params, outdir, props, types, name="deconv_bar", title="Estimated cell type composition"):
    labels = _label_array(ctx, str(_p(params, "group_col", "") or ""))
    order = None
    rows = []
    if labels is not None:
        order = list(dict.fromkeys([str(v) for v in labels.tolist()]))
        for value in order:
            sel = np.array([str(v) == value for v in labels])
            rows.append({"group": value, "n_spots": int(sel.sum()),
                         **{types[j]: float(props[sel, j].mean()) for j in range(len(types))}})
        xlabels = order
        means = np.array([[row[t] for t in types] for row in rows], dtype=float)
    else:
        means = props.mean(axis=0, keepdims=True)
        xlabels = ["all spots"]
    fig, ax = sc.new_fig(7.4, 4.4)
    bottom = np.zeros(means.shape[0], dtype=float)
    for j, t in enumerate(types):
        ax.bar(np.arange(means.shape[0]), means[:, j], bottom=bottom,
               color=sc.PALETTE[j % len(sc.PALETTE)], label=t, width=0.72)
        bottom += means[:, j]
    ax.set_xticks(np.arange(len(xlabels)))
    ax.set_xticklabels(xlabels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("mean proportion", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, frameon=False, ncol=2, loc="upper right")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return sc.save_fig(fig, outdir, name, title=None)


# ---------------------------------------------------------------------------
# 步骤
# ---------------------------------------------------------------------------


def step_deconv(ctx, params, outdir):
    out = _deconvolve(ctx, params)
    props, types = out["proportions"], out["cell_types"]
    outdir = _outdir(ctx, outdir)
    spots = _spots(ctx)
    frame = pd.DataFrame(props, columns=types)
    frame.insert(0, "spot", spots[:len(frame)])
    rows = frame.head(500).to_dict(orient="records")
    figures = [_fig_composition(ctx, params, outdir, props, types)]
    top_type = types[int(np.argmax(props.mean(axis=0)))] if types else (types[0] if types else "")
    if top_type:
        fig, ax = sc.new_fig(6.2, 5.0)
        sc.scatter_spatial(ax, ctx, values=props[:, types.index(top_type)], point_size=6.0,
                           legend_title=top_type)
        ax.set_title(f"Estimated proportion: {top_type}", fontsize=10)
        figures.append(sc.save_fig(fig, outdir, "deconv-map", title=None))
    bullets = [
        f"{out['method'].upper()} 反卷积：{len(spots)} 个 spot × {len(types)} 种细胞（共同基因 {out['n_genes']} 个）",
        "各类型平均比例：" + "，".join(f"{t}:{float(np.mean(props[:, j])):.2f}" for j, t in enumerate(types[:6])),
    ]
    return sc.result(
        summary={"title": "反卷积", "bullets": bullets,
                 "metrics": {"n_types": len(types), "n_genes": out["n_genes"], "method": out["method"]}},
        tables=[sc.table("deconv", "每个 spot 的细胞类型比例",
                         columns=[{"key": "spot", "label": "spot", "type": "string"}]
                                 + [{"key": t, "label": t, "type": "number"} for t in types],
                         rows=rows)],
        figures=figures,
        data_out=[{"kind": "deconv_csv", "path": sc.export_csv(frame, outdir, "deconv-proportions.csv", index=False),
                   "label": "反卷积比例表"}],
        store={"deconv": {"cell_types": types, "proportions": props, "method": out["method"]}},
    )


# ---------------------------------------------------------------------------
# 分析
# ---------------------------------------------------------------------------


def analysis_deconv_nnls(ctx, params, outdir):
    merged = dict(params or {})
    merged["method"] = "nnls"
    out = _deconvolve(ctx, merged)
    outdir = _outdir(ctx, outdir)
    props, types = out["proportions"], out["cell_types"]
    figures = [_fig_composition(ctx, params, outdir, props, types, name="deconv-nnls-bar",
                               title="NNLS deconvolution composition")]
    bullets = [f"NNLS 反卷积（Tsoucas 2019 / SPOTlight 2021 思路）：{len(types)} 种细胞，共同基因 {out['n_genes']} 个",
               "平均比例：" + "，".join(f"{t}:{float(np.mean(props[:, j])):.2f}" for j, t in enumerate(types[:6]))]
    return sc.result(summary={"title": "反卷积（NNLS）", "bullets": bullets,
                              "metrics": {"n_types": len(types), "n_genes": out["n_genes"]}},
                     tables=[sc.table("deconv_nnls", "平均比例",
                                      columns=[{"key": "cell_type", "label": "细胞类型", "type": "string"},
                                               {"key": "mean", "label": "平均比例", "type": "number"},
                                               {"key": "median", "label": "中位比例", "type": "number"}],
                                      rows=[{"cell_type": t, "mean": float(np.mean(props[:, j])),
                                             "median": float(np.median(props[:, j]))}
                                            for j, t in enumerate(types)])],
                     figures=figures,
                     store={"deconv_nnls": {"cell_types": types, "proportions": props}})


def analysis_deconv_dwls(ctx, params, outdir):
    merged = dict(params or {})
    merged["method"] = "dwls"
    out = _deconvolve(ctx, merged)
    outdir = _outdir(ctx, outdir)
    props, types = out["proportions"], out["cell_types"]
    figures = [_fig_composition(ctx, params, outdir, props, types, name="deconv-dwls-bar",
                               title="Dampened weighted least squares composition")]
    bullets = [f"阻尼加权最小二乘反卷积（Tsoucas 2019）：{len(types)} 种细胞，共同基因 {out['n_genes']} 个"]
    return sc.result(summary={"title": "反卷积（DWLS）", "bullets": bullets,
                              "metrics": {"n_types": len(types), "n_genes": out["n_genes"]}},
                     tables=[sc.table("deconv_dwls", "平均比例",
                                      columns=[{"key": "cell_type", "label": "细胞类型", "type": "string"},
                                               {"key": "mean", "label": "平均比例", "type": "number"}],
                                      rows=[{"cell_type": t, "mean": float(np.mean(props[:, j]))}
                                            for j, t in enumerate(types)])],
                     figures=figures,
                     store={"deconv_dwls": {"cell_types": types, "proportions": props}})


def _cached_props(ctx, params):
    store = ctx.get("results") or {}
    cell_type = str(_p(params, "cell_type", "") or "")
    for key in ("deconv", "deconv_nnls", "deconv_dwls"):
        entry = store.get(key)
        if isinstance(entry, dict) and "proportions" in entry:
            names = list(entry.get("cell_types") or [])
            props = np.asarray(entry["proportions"], dtype=np.float32)
            if cell_type and cell_type in names:
                return props[:, names.index(cell_type)], cell_type, names
            if not cell_type and names:
                return props[:, 0], names[0], names
    out = _deconvolve(ctx, params)
    names = out["cell_types"]
    if cell_type and cell_type in names:
        return out["proportions"][:, names.index(cell_type)], cell_type, names
    return out["proportions"][:, 0], names[0], names


def analysis_deconv_spatial_map(ctx, params, outdir):
    values, cell_type, names = _cached_props(ctx, params)
    outdir = _outdir(ctx, outdir)
    fig, ax = sc.new_fig(6.2, 5.0)
    sc.scatter_spatial(ax, ctx, values=values, point_size=_p_float(params, "point_size", 8.0),
                       legend_title=cell_type)
    ax.set_title(f"Estimated proportion: {cell_type}", fontsize=10)
    figure = sc.save_fig(fig, outdir, "deconv-spatial", title=None)
    bullets = [f"细胞类型「{cell_type}」在组织切片上的比例分布（均值 {float(np.mean(values)):.3f}）",
               f"可用类型：{'、'.join(names[:8])}"]
    return sc.result(summary={"title": "细胞类型空间分布", "bullets": bullets,
                              "metrics": {"cell_type": cell_type, "mean": float(np.mean(values))}},
                     figures=[figure])


def analysis_deconv_compare(ctx, params, outdir):
    a = _deconvolve(ctx, {**(params or {}), "method": "nnls"})
    b = _deconvolve(ctx, {**(params or {}), "method": "dwls"})
    outdir = _outdir(ctx, outdir)
    types = a["cell_types"]
    rows, corr = [], {}
    for j, t in enumerate(types):
        va = np.asarray(a["proportions"][:, j], dtype=float)
        vb = np.asarray(b["proportions"][:, j], dtype=float)
        if va.size < 3 or np.std(va) == 0 or np.std(vb) == 0:
            r = float("nan")
        else:
            r = float(np.corrcoef(va, vb)[0, 1])
        corr[t] = r
        rows.append({"cell_type": t, "pearson": r, "nnls_mean": float(va.mean()), "dwls_mean": float(vb.mean())})
    fig, ax = sc.new_fig(6.6, 4.2)
    values = [corr[t] for t in types]
    ax.bar(np.arange(len(types)), values, color=[sc.PALETTE[i % len(sc.PALETTE)] for i in range(len(types))],
           width=0.7)
    ax.set_xticks(np.arange(len(types)))
    ax.set_xticklabels(types, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Pearson r (NNLS vs DWLS)", fontsize=9)
    ax.set_ylim(min(0, min(v for v in values if np.isfinite(v)) - 0.05) if any(np.isfinite(v) for v in values) else 0, 1.02)
    ax.set_title("Consistency between deconvolution methods", fontsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    figure = sc.save_fig(fig, outdir, "deconv-compare", title=None)
    finite = [v for v in values if np.isfinite(v)]
    bullets = [f"两种方法一致性的中位相关系数 {np.median(finite):.3f}" if finite else "两种方法在多数细胞类型上无变异，无法比较"]
    return sc.result(summary={"title": "反卷积方法一致性", "bullets": bullets,
                              "metrics": {"median_r": float(np.median(finite)) if finite else None}},
                     tables=[sc.table("deconv_compare", "NNLS vs DWLS",
                                      columns=[{"key": "cell_type", "label": "细胞类型", "type": "string"},
                                               {"key": "pearson", "label": "Pearson r", "type": "number"},
                                               {"key": "nnls_mean", "label": "NNLS 均值", "type": "number"},
                                               {"key": "dwls_mean", "label": "DWLS 均值", "type": "number"}],
                                      rows=rows)],
                     figures=[figure])


def analysis_deconv_roi(ctx, params, outdir):
    store = ctx.get("results") or {}
    entry = store.get("deconv") or store.get("deconv_nnls") or store.get("deconv_dwls")
    if not isinstance(entry, dict) or "proportions" not in entry:
        raise sc.prerequisite("请先跑「反卷积」步骤（或 deconv_nnls 分析），再做区域组成差异")
    props = np.asarray(entry["proportions"], dtype=np.float32)
    types = list(entry.get("cell_types") or [])
    labels = _label_array(ctx, str(_p(params, "group_col", "") or ""))
    if labels is None:
        raise sc.prerequisite("找不到分组/域标签：请先设置分组列或跑空间域")
    labels = np.asarray([str(v) for v in labels.tolist()], dtype=object)
    uniq = list(dict.fromkeys(labels.tolist()))
    if len(uniq) < 2:
        raise sc.prerequisite("分组只有 1 个取值，无法比较")
    rows = []
    from scipy import stats

    for j, t in enumerate(types):
        groups = [props[labels == u, j] for u in uniq]
        groups = [g for g in groups if g.size > 0]
        if len(groups) < 2:
            continue
        try:
            stat, p = stats.kruskal(*groups)
        except Exception:  # noqa: BLE001
            stat, p = float("nan"), float("nan")
        row = {"cell_type": t, "stat": float(stat), "p": float(p)}
        for u, g in zip(uniq, groups):
            row[f"mean_{u}"] = float(np.mean(g))
        rows.append(row)
    try:
        pvals = [row["p"] for row in rows]
        adjusted = sc.bh_fdr(pvals)
        for row, q in zip(rows, adjusted):
            row["p_adj"] = float(q)
    except Exception:  # noqa: BLE001
        pass
    outdir = _outdir(ctx, outdir)
    figure = _fig_composition(ctx, params, outdir, props, types, name="deconv-roi-bar",
                              title="Cell type composition by group")
    bullets = [f"按「{'/'.join(uniq[:6])}」分组比较 {len(types)} 种细胞比例",
               "Kruskal-Wallis + BH 校正；p_adj<0.05 的类型：" +
               ("、".join(row["cell_type"] for row in rows if row.get("p_adj", 1) < 0.05) or "无")]
    return sc.result(summary={"title": "区域细胞组成差异", "bullets": bullets,
                              "metrics": {"n_types": len(types), "n_groups": len(uniq)}},
                     tables=[sc.table("deconv_roi", "各分组细胞比例（Kruskal-Wallis）",
                                      columns=[{"key": "cell_type", "label": "细胞类型", "type": "string"},
                                               {"key": "stat", "label": "统计量", "type": "number"},
                                               {"key": "p", "label": "p", "type": "number"},
                                               {"key": "p_adj", "label": "p_adj", "type": "number"}]
                                              + [{"key": f"mean_{u}", "label": f"{u} 均值", "type": "number"}
                                                 for u in uniq],
                                      rows=rows)],
                     figures=[figure])


def analysis_deconv_marker_corr(ctx, params, outdir):
    store = ctx.get("results") or {}
    entry = store.get("deconv") or store.get("deconv_nnls")
    if not isinstance(entry, dict) or "proportions" not in entry:
        raise sc.prerequisite("请先跑「反卷积」步骤，再比较比例与标记基因表达")
    types = list(entry.get("cell_types") or [])
    cell_type = str(_p(params, "cell_type", "") or "")
    if cell_type and cell_type not in types:
        raise sc.prerequisite(f"细胞类型「{cell_type}」不在反卷积结果里（可用：{'、'.join(types[:8])}）")
    if not cell_type:
        cell_type = types[0] if types else ""
    props = np.asarray(entry["proportions"], dtype=np.float32)
    j = types.index(cell_type)
    ratio = props[:, j]
    X = _X(ctx)
    genes = _genes(ctx)
    requested = [str(g) for g in (_p(params, "genes", []) or []) if str(g).strip()]
    index = {g: i for i, g in enumerate(genes)}
    requested = [g for g in requested if g in index]
    if not requested:
        # 取与比例相关性最高的 10 个基因
        corr = np.zeros(X.shape[1], dtype=np.float32)
        rc = ratio - ratio.mean()
        denom = np.linalg.norm(rc) + 1e-9
        for gi in range(X.shape[1]):
            col = np.asarray(X[:, gi], dtype=np.float32)
            cc = col - col.mean()
            corr[gi] = float(cc @ rc / (np.linalg.norm(cc) * denom + 1e-12))
        requested = [genes[i] for i in np.argsort(-corr)[:10]]
    rows = []
    rc = ratio - ratio.mean()
    for gene in requested:
        col = np.asarray(X[:, index[gene]], dtype=np.float32)
        cc = col - col.mean()
        r = float(cc @ rc / (np.linalg.norm(cc) * (np.linalg.norm(rc) + 1e-12) + 1e-12))
        rows.append({"gene": gene, "corr": r})
    outdir = _outdir(ctx, outdir)
    fig, ax = sc.new_fig(6.4, 4.2)
    ax.bar(np.arange(len(rows)), [row["corr"] for row in rows],
           color=[sc.PALETTE[i % len(sc.PALETTE)] for i in range(len(rows))], width=0.7)
    ax.set_xticks(np.arange(len(rows)))
    ax.set_xticklabels([row["gene"] for row in rows], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Pearson r with proportion", fontsize=9)
    ax.set_title(f"Marker consistency ({cell_type})", fontsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    figure = sc.save_fig(fig, outdir, "deconv-marker-corr", title=None)
    bullets = [f"「{cell_type}」比例与 {len(rows)} 个基因表达的相关性",
               "相关最高的基因：" + "、".join(row["gene"] for row in rows[:5])]
    return sc.result(summary={"title": "比例-标记基因一致性", "bullets": bullets,
                              "metrics": {"cell_type": cell_type, "n_genes": len(rows)}},
                     tables=[sc.table("deconv_marker", f"与「{cell_type}」比例相关的基因",
                                      columns=[{"key": "gene", "label": "基因", "type": "string"},
                                               {"key": "corr", "label": "相关系数", "type": "number"}],
                                      rows=rows)],
                     figures=[figure])


# ---------------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------------


def plot_deconv_bar(ctx, params=None, outdir=None):
    store = ctx.get("results") or {}
    entry = store.get("deconv") or store.get("deconv_nnls") or store.get("deconv_dwls")
    if not isinstance(entry, dict) or "proportions" not in entry:
        merged = dict(params or {})
        out = _deconvolve(ctx, merged)
        props, types = out["proportions"], out["cell_types"]
    else:
        props = np.asarray(entry["proportions"], dtype=np.float32)
        types = list(entry.get("cell_types") or [])
    return _fig_composition(ctx, params or {}, _outdir(ctx, outdir), props, types)


def plot_deconv_spatial(ctx, params=None, outdir=None):
    values, cell_type, names = _cached_props(ctx, params or {})
    outdir = _outdir(ctx, outdir)
    fig, ax = sc.new_fig(6.2, 5.0)
    sc.scatter_spatial(ax, ctx, values=values, point_size=_p_float(params, "point_size", 8.0),
                       legend_title=cell_type)
    ax.set_title(f"Estimated proportion: {cell_type}", fontsize=10)
    return sc.save_fig(fig, outdir, "deconv-spatial", title=None)


# ---------------------------------------------------------------------------
# 注册
# ---------------------------------------------------------------------------

STEPS = [
    {
        "id": "deconv", "title": "反卷积", "icon": "deconv", "requires": ["normalize"],
        "params": [
            {"key": "ref", "label": "参考 h5ad 路径（可留空用已挂载）", "type": "text", "default": ""},
            {"key": "method", "label": "方法", "type": "select", "default": "nnls",
             "options": [{"value": "nnls", "label": "NNLS 非负最小二乘"},
                         {"value": "dwls", "label": "DWLS 阻尼加权"}]},
            {"key": "cell_type_col", "label": "参考的细胞类型列", "type": "text", "default": ""},
            {"key": "min_genes", "label": "共同基因下限", "type": "number", "default": 50},
        ],
        "run": step_deconv,
    },
]

ANALYSES = [
    {
        "id": "deconv_nnls", "title": "反卷积（NNLS）", "group": "反卷积",
        "method": "非负最小二乘（Tsoucas 2019 DWLS 基础；SPOTlight 2021 思路；内置实现）",
        "needs": ["normalize"],
        "params": [{"key": "min_genes", "label": "共同基因下限", "type": "number", "default": 50},
                   {"key": "cell_type_col", "label": "细胞类型列", "type": "text", "default": ""}],
        "run": analysis_deconv_nnls,
    },
    {
        "id": "deconv_dwls", "title": "反卷积（加权 NNLS）", "group": "反卷积",
        "method": "阻尼加权最小二乘（Tsoucas 2019；内置实现）",
        "needs": ["normalize"],
        "params": [{"key": "min_genes", "label": "共同基因下限", "type": "number", "default": 50},
                   {"key": "cell_type_col", "label": "细胞类型列", "type": "text", "default": ""}],
        "run": analysis_deconv_dwls,
    },
    {
        "id": "deconv_spatial_map", "title": "细胞类型空间分布", "group": "反卷积",
        "method": "把比例映射到空间坐标（内置实现）",
        "needs": ["normalize"],
        "params": [{"key": "cell_type", "label": "细胞类型（留空取第一个）", "type": "text", "default": ""},
                   {"key": "point_size", "label": "点大小", "type": "number", "default": 8}],
        "run": analysis_deconv_spatial_map,
    },
    {
        "id": "deconv_compare", "title": "反卷积方法一致性", "group": "反卷积",
        "method": "NNLS vs DWLS 的逐类型 Pearson 相关（内置实现）",
        "needs": ["normalize"],
        "params": [],
        "run": analysis_deconv_compare,
    },
    {
        "id": "deconv_roi", "title": "区域细胞组成差异", "group": "反卷积",
        "method": "域/分组内比例 + Kruskal-Wallis + BH（内置实现）",
        "needs": ["normalize", "cluster"],
        "params": [{"key": "group_col", "label": "分组列（留空用域标签）", "type": "column", "default": None}],
        "run": analysis_deconv_roi,
    },
    {
        "id": "deconv_marker_corr", "title": "比例-标记基因一致性", "group": "反卷积",
        "method": "比例与基因表达的相关（内置实现）",
        "needs": ["normalize"],
        "params": [{"key": "cell_type", "label": "细胞类型", "type": "text", "default": ""},
                   {"key": "genes", "label": "基因列表（可留空）", "type": "genes", "default": []}],
        "run": analysis_deconv_marker_corr,
    },
]

PLOTS = {
    "deconv_bar": plot_deconv_bar,
    "deconv_spatial": plot_deconv_spatial,
}

__all__ = ["STEPS", "ANALYSES", "PLOTS"]
