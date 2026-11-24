"""
dsh-spatial — 基础流程步骤与分析（DESIGN.md §3 / §4 group=质控与整合）

本模块提供：
  STEPS   : qc / filter / normalize / hvg / pca / neighbors / cluster / embed / annotate
  ANALYSES: qc_summary / outlier_spots / slice_align / spatial_trajectory

全部用 numpy/scipy/sklearn 自实现（本机无 scanpy/umap-learn/leidenalg）：
  * 标准化        log1p(CPM-like)（扫描式文库大小校正，Luecken & Theis 2019 常规做法）
  * 高变基因      Seurat v3 的分箱 dispersion（Stuart 2019）或 Cell Ranger 的归一化 dispersion
  * PCA           z-score 后 SVD
  * 邻接图        表达 kNN（cosine）与空间图按 spatial_weight 融合（空间平滑正则，BayesSpace/GraphST 思想）
  * 聚类          networkx Louvain（Blondel 2008），resolution 控制粒度
  * 二维嵌入      sklearn TSNE（van der Maaten & Hinton 2008）或 PCA
  * 空间轨迹      扩散伪时序（Haghverdi 2016 DPT 思想：转移矩阵特征分解 + 扩散距离）
"""

from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd

from sp_common import (
    EngineError,
    PALETTE,
    add_image_background,
    bh_fdr,
    categorical_colors,
    coord_knn,
    ensure_dir,
    estimate_radius,
    export_csv,
    new_fig,
    palette,
    prerequisite,
    radius_graph,
    result,
    save_fig,
    scatter_spatial,
    spatial_weights,
    table,
)

# ---------------------------------------------------------------------------
# ctx 小工具
# ---------------------------------------------------------------------------


def _X(ctx) -> np.ndarray:
    X = ctx.get("X")
    if X is None:
        raise prerequisite("数据集没有表达矩阵")
    arr = np.asarray(X, dtype=np.float32)
    if arr.ndim != 2:
        raise EngineError("internal", f"表达矩阵维度异常：{arr.shape}")
    return arr


def _genes(ctx) -> list[str]:
    return [str(item) for item in (ctx.get("genes") or [])]


def _spots(ctx) -> list[str]:
    return [str(item) for item in (ctx.get("spots") or [])]


def _coords(ctx) -> np.ndarray:
    coords = ctx.get("coords")
    if coords is None:
        raise prerequisite("该数据集没有空间坐标：导 Visium 目录或提供 coords 坐标表后再跑空间分析")
    arr = np.asarray(coords, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] < 2 or arr.shape[0] != _X(ctx).shape[0]:
        raise prerequisite("空间坐标形状与表达矩阵不一致（spot 数对不上）")
    return arr[:, :2]


def _info(ctx) -> pd.DataFrame:
    frame = ctx.get("sample_info")
    if isinstance(frame, pd.DataFrame):
        return frame
    return pd.DataFrame(index=pd.Index(_spots(ctx), name="spot"))


def _column(ctx, key: str):
    if not key:
        return None
    frame = _info(ctx)
    if key in frame.columns:
        return np.asarray([str(item) for item in frame[key].tolist()], dtype=object)
    return None


def _label_columns(ctx) -> list[str]:
    frame = _info(ctx)
    out = []
    for col in frame.columns:
        series = frame[col]
        try:
            if series.nunique(dropna=True) <= max(60, int(len(series) * 0.2)) and series.nunique() > 1:
                out.append(str(col))
        except Exception:  # noqa: BLE001
            continue
    return out


def _store(ctx) -> dict:
    store = ctx.get("results")
    return store if isinstance(store, dict) else {}


def _sorted_labels(values) -> list[str]:
    labels = [str(item) for item in values]
    uniq = list(dict.fromkeys(labels))
    try:
        return sorted(uniq, key=lambda item: (0, float(item)) if _is_number(item) else (1, item))
    except Exception:  # noqa: BLE001
        return uniq


def _is_number(value) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def _numeric_params(params: dict, key: str, default=None):
    value = (params or {}).get(key, default)
    if value in (None, "", "null"):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _int_params(params: dict, key: str, default: int) -> int:
    value = _numeric_params(params, key, default)
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(default)


def _bool_params(params: dict, key: str, default: bool = False) -> bool:
    value = (params or {}).get(key, default)
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on", "是")
    return bool(value)


def _subsample(n: int, limit: int, seed: int = 0) -> np.ndarray:
    if n <= limit:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, size=int(limit), replace=False))


def _dense(value) -> np.ndarray:
    arr = value.toarray() if hasattr(value, "toarray") else value
    return np.asarray(arr, dtype=np.float32)


def _mt_mask(genes: list[str]) -> np.ndarray:
    return np.array([bool(g) and (str(g).upper().startswith("MT-") or str(g).startswith("mt-"))
                     for g in genes], dtype=bool)


def _spot_metrics(X: np.ndarray, genes: list[str]) -> pd.DataFrame:
    counts = np.asarray(X.sum(axis=1), dtype=np.float64)
    detected = np.asarray((X > 0).sum(axis=1), dtype=np.int64)
    mt = _mt_mask(genes)
    if mt.any() and counts.max() > 0:
        mt_counts = np.asarray(X[:, mt].sum(axis=1), dtype=np.float64)
        pct_mt = np.where(counts > 0, mt_counts / np.maximum(counts, 1e-9) * 100.0, 0.0)
    else:
        pct_mt = np.zeros(X.shape[0], dtype=np.float64)
    return pd.DataFrame({
        "spot": list(range(X.shape[0])),
        "counts": counts,
        "genes": detected,
        "pct_mt": pct_mt,
    })


def _gene_stats(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.asarray(X.mean(axis=0), dtype=np.float64)
    var = np.asarray(X.var(axis=0), dtype=np.float64)
    nonzero = np.asarray((X > 0).mean(axis=0), dtype=np.float64)
    return mean, var, nonzero


# ---------------------------------------------------------------------------
# 图：本模块的步骤自带图（kind 由 sp_plots 负责，这里只画步骤结果）
# ---------------------------------------------------------------------------


def _fig_spatial_metric(ctx, outdir: str, values, name: str, title: str,
                        cbar_label: str, point_size: float = 6.0) -> dict:
    fig, ax = new_fig(6.2, 5.0)
    scatter_spatial(ax, ctx, values=np.asarray(values, dtype=float), point_size=point_size,
                    legend_title=cbar_label)
    ax.set_title(title, fontsize=10)
    return save_fig(fig, outdir, name, title=None)


def _fig_hist(ctx, outdir: str, series: list[tuple[str, np.ndarray]], name: str,
              title: str, xlabel: str) -> dict:
    fig, ax = new_fig(6.6, 4.2)
    for index, (label, values) in enumerate(series):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            continue
        ax.hist(values, bins=48, alpha=0.55, color=PALETTE[index % len(PALETTE)], label=label,
                density=False)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("spots", fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return save_fig(fig, outdir, name, title=None)


def _fig_scatter(ctx, outdir: str, x, y, name: str, title: str, xlabel: str, ylabel: str,
                 color=None, label=None) -> dict:
    fig, ax = new_fig(6.2, 4.6)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    if color is not None and label is not None:
        labels = [str(item) for item in label]
        mapping = categorical_colors(labels)
        for value in mapping:
            sel = np.array([item == value for item in labels]) & finite
            ax.scatter(x[sel], y[sel], s=10, color=mapping[value], label=value, linewidths=0)
        ax.legend(fontsize=7, frameon=False, markerscale=2)
    else:
        ax.scatter(x[finite], y[finite], s=8,
                   color=PALETTE[0] if color is None else color, linewidths=0, alpha=0.8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return save_fig(fig, outdir, name, title=None)


def _fig_bar(ctx, outdir: str, labels: list[str], values, name: str, title: str,
             ylabel: str = "") -> dict:
    fig, ax = new_fig(6.8, 4.0)
    positions = np.arange(len(labels))
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(labels))]
    ax.bar(positions, np.asarray(values, dtype=float), color=colors, width=0.72)
    ax.set_xticks(positions)
    ax.set_xticklabels([str(item) for item in labels], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return save_fig(fig, outdir, name, title=None)


# ---------------------------------------------------------------------------
# 步骤：qc
# ---------------------------------------------------------------------------


def step_qc(ctx, params, outdir):
    X = _X(ctx)
    genes = _genes(ctx)
    spots = _spots(ctx)
    metrics = _spot_metrics(X, genes)
    min_counts = _numeric_params(params, "min_counts", None)
    min_genes = _numeric_params(params, "min_genes", None)
    max_pct_mt = _numeric_params(params, "max_pct_mt", 20.0)
    passed = np.ones(X.shape[0], dtype=bool)
    notes = []
    if min_counts is not None:
        passed &= metrics["counts"].to_numpy() >= float(min_counts)
        notes.append(f"counts ≥ {result_fmt(min_counts)}")
    if min_genes is not None:
        passed &= metrics["genes"].to_numpy() >= float(min_genes)
        notes.append(f"检出基因 ≥ {result_fmt(min_genes)}")
    if max_pct_mt is not None and _mt_mask(genes).any():
        passed &= metrics["pct_mt"].to_numpy() <= float(max_pct_mt)
        notes.append(f"线粒体比例 ≤ {result_fmt(max_pct_mt)}%")
    metrics["passed"] = passed

    frame = metrics.copy()
    frame.insert(1, "spot_name", [spots[i] if i < len(spots) else str(i) for i in range(len(frame))])
    rows = frame.to_dict(orient="records")

    figures = [
        _fig_spatial_metric(ctx, outdir, metrics["counts"].to_numpy(), "qc-counts-spatial",
                            "UMI counts per spot", "counts"),
        _fig_spatial_metric(ctx, outdir, metrics["pct_mt"].to_numpy(), "qc-mt-spatial",
                            "Mitochondrial fraction (%)", "%MT"),
        _fig_hist(ctx, outdir, [("counts", metrics["counts"].to_numpy())], "qc-counts-hist",
                  "Counts distribution", "UMI counts"),
    ]
    n_fail = int((~passed).sum())
    bullets = [
        f"{X.shape[0]} 个 spot，中位 counts {np.median(metrics['counts']):.0f}，"
        f"中位检出基因 {np.median(metrics['genes']):.0f} 个",
        f"线粒体比例中位 {np.median(metrics['pct_mt']):.2f}%"
        + ("（没有 MT- 基因，按 0 处理）" if not _mt_mask(genes).any() else ""),
    ]
    if n_fail:
        bullets.append(f"按当前阈值 {n_fail} 个 spot 会被过滤（{n_fail / max(X.shape[0], 1) * 100:.1f}%）")
    elif notes:
        bullets.append("按当前阈值没有 spot 被过滤")
    data_out = []
    try:
        path = export_csv(frame, outdir, "qc-spots.csv", index=False)
        data_out.append({"kind": "qc_csv", "path": path, "label": "spot 级质控表"})
    except Exception:  # noqa: BLE001
        pass
    payload = result(
        summary={"title": "空间质控", "bullets": bullets, "metrics": {
            "n_spots": int(X.shape[0]), "median_counts": float(np.median(metrics["counts"])),
            "median_genes": float(np.median(metrics["genes"])),
            "median_pct_mt": float(np.median(metrics["pct_mt"])),
            "fail": n_fail, "thresholds": notes,
        }},
        tables=[table("qc", "spot 级质控指标",
                      columns=[{"key": "spot_name", "label": "spot", "type": "string"},
                               {"key": "counts", "label": "counts", "type": "number"},
                               {"key": "genes", "label": "检出基因", "type": "number"},
                               {"key": "pct_mt", "label": "%MT", "type": "number"},
                               {"key": "passed", "label": "通过", "type": "boolean"}],
                      rows=rows)],
        figures=figures,
        data_out=data_out,
    )
    return payload


def result_fmt(value) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{number:.0f}" if abs(number - round(number)) < 1e-9 else f"{number:g}"


# ---------------------------------------------------------------------------
# 步骤：filter
# ---------------------------------------------------------------------------


def step_filter(ctx, params, outdir):
    X = _X(ctx)
    genes = _genes(ctx)
    spots = _spots(ctx)
    coords = _coords(ctx)
    count_arr = ctx.get("counts")
    info = _info(ctx).copy()

    min_counts = _numeric_params(params, "min_counts", 500)
    min_genes = _numeric_params(params, "min_genes", 200)
    max_pct_mt = _numeric_params(params, "max_pct_mt", 20)
    min_cells = _int_params(params, "min_cells", 3)

    metrics = _spot_metrics(X, genes)
    keep = np.ones(X.shape[0], dtype=bool)
    if min_counts is not None:
        keep &= metrics["counts"].to_numpy() >= float(min_counts)
    if min_genes is not None:
        keep &= metrics["genes"].to_numpy() >= float(min_genes)
    if max_pct_mt is not None and _mt_mask(genes).any():
        keep &= metrics["pct_mt"].to_numpy() <= float(max_pct_mt)
    if keep.sum() < max(10, int(0.1 * X.shape[0])):
        raise prerequisite(
            f"按当前阈值只剩 {int(keep.sum())} 个 spot（< 10%）：请放宽 min_counts/min_genes/max_pct_mt")

    X1 = X[keep]
    spots1 = [spot for spot, ok in zip(spots, keep) if ok]
    coords1 = coords[keep]
    info1 = info.loc[spots1] if len(info) == len(spots) else pd.DataFrame(index=pd.Index(spots1, name="spot"))
    info1.index = pd.Index(spots1, name="spot")

    detect = np.asarray((X1 > 0).sum(axis=0))
    gene_keep = detect >= max(1, int(min_cells))
    if gene_keep.sum() < 50:
        raise prerequisite(f"按 min_cells={min_cells} 过滤后只剩 {int(gene_keep.sum())} 个基因，请降低阈值")
    X2 = X1[:, gene_keep]
    genes2 = [gene for gene, ok in zip(genes, gene_keep) if ok]
    counts1 = ctx.get("counts")
    counts2 = None
    if counts1 is not None:
        counts2 = _dense(counts1[keep])[:, gene_keep]

    replace = {"X": X2, "genes": genes2, "spots": spots1, "coords": coords1,
               "sample_info": info1, "counts": counts2}
    bullets = [
        f"spot：{X.shape[0]} → {X2.shape[0]}（剔除 {int((~keep).sum())}）",
        f"基因：{X.shape[1]} → {X2.shape[1]}（剔除 {int((~gene_keep).sum())}）",
        f"阈值：counts≥{result_fmt(min_counts)}、检出基因≥{result_fmt(min_genes)}、"
        f"%MT≤{result_fmt(max_pct_mt)}、基因至少在 {min_cells} 个 spot 检出",
    ]
    figures = [
        _fig_spatial_metric(ctx, outdir, keep.astype(float), "filter-keep-spatial",
                            "Kept spots (1=kept)", "kept", point_size=6.0),
    ]
    return result(
        summary={"title": "过滤", "bullets": bullets, "metrics": {
            "n_spots_before": int(X.shape[0]), "n_spots_after": int(X2.shape[0]),
            "n_genes_before": int(X.shape[1]), "n_genes_after": int(X2.shape[1]),
        }},
        figures=figures,
        replace=replace,
    )


# ---------------------------------------------------------------------------
# 步骤：normalize
# ---------------------------------------------------------------------------


def step_normalize(ctx, params, outdir):
    X = _X(ctx)
    target_sum = _numeric_params(params, "target_sum", 10000.0) or 10000.0
    use_log = _bool_params(params, "log", True)
    method = str((params or {}).get("method") or "lognorm").strip().lower()
    kind = str(ctx.get("X_kind") or "unknown")
    if method == "none":
        return result(summary={"title": "标准化（跳过）",
                               "bullets": ["按参数跳过标准化：注意后续分析的可比性"],
                               "metrics": {"X_kind": kind}})
    if kind == "lognorm":
        return result(summary={"title": "标准化（已标准化）",
                               "bullets": ["数据已经是 lognorm 尺度，本次不重复标准化"],
                               "metrics": {"X_kind": kind, "skipped": True}})
    arr = np.asarray(X, dtype=np.float64)
    lib = np.asarray(arr.sum(axis=1), dtype=np.float64)
    lib[lib <= 0] = 1.0
    scaled = arr / lib[:, None] * float(target_sum)
    if use_log:
        scaled = np.log1p(scaled)
    scaled = scaled.astype(np.float32)
    before_median = float(np.median(lib))
    figures = [
        _fig_hist(ctx, outdir, [("library size (before)", lib)], "normalize-libsize",
                  "Library size before normalization", "library size"),
        _fig_scatter(ctx, outdir, lib, np.asarray(scaled.sum(axis=1)),
                     "normalize-total", "Total before vs after normalization",
                     "library size (before)", "total (after, log1p)"),
    ]
    bullets = [
        f"counts → {'log1p' if use_log else ''}(CPM-like)，目标总量 {result_fmt(target_sum)}",
        f"标准化前中位文库大小 {before_median:.0f}；输出 X_kind=lognorm",
    ]
    return result(
        summary={"title": "标准化", "bullets": bullets,
                 "metrics": {"target_sum": float(target_sum), "log": bool(use_log),
                             "median_lib": before_median, "X_kind": "lognorm"}},
        figures=figures,
        replace={"X": scaled, "X_kind": "lognorm"},
    )


# ---------------------------------------------------------------------------
# 步骤：hvg
# ---------------------------------------------------------------------------


def step_hvg(ctx, params, outdir):
    X = _X(ctx)
    genes = _genes(ctx)
    n_top = _int_params(params, "n_top_genes", 2000)
    flavor = str((params or {}).get("flavor") or "seurat").strip().lower()
    if str(ctx.get("X_kind")) == "counts":
        raise prerequisite("请先跑标准化（hvg 需要 lognorm 表达值）")
    mean, var, nonzero = _gene_stats(X)
    valid = (mean > 0) & (nonzero > 0.01)
    if valid.sum() < 10:
        raise prerequisite("可用基因太少（很多基因在绝大多数 spot 不表达）")
    disp = np.where(mean > 0, var / np.maximum(mean, 1e-9), 0.0)
    if flavor == "cell_ranger":
        bins = np.linspace(np.log1p(mean[valid].min()), np.log1p(mean[valid].max()) + 1e-9, 21)
        score = np.zeros_like(disp)
        idx_valid = np.where(valid)[0]
        assign = np.digitize(np.log1p(mean[idx_valid]), bins) - 1
        for group in np.unique(assign):
            sel = idx_valid[assign == group]
            if sel.size < 3:
                score[sel] = 0.0
                continue
            values = disp[sel]
            mu, sd = float(np.mean(values)), float(np.std(values) + 1e-12)
            score[sel] = (values - mu) / sd
    else:
        # Seurat v3：按平均表达分箱后对 dispersion 做组内 z-score
        idx_valid = np.where(valid)[0]
        values = np.log1p(mean[idx_valid])
        bins = np.quantile(values, np.linspace(0, 1, 21))
        bins[0] -= 1e-9
        bins[-1] += 1e-9
        assign = np.digitize(np.log1p(mean[idx_valid]), bins) - 1
        score = np.zeros_like(disp)
        for group in np.unique(assign):
            sel = idx_valid[assign == group]
            if sel.size < 3:
                score[sel] = 0.0
                continue
            sub = disp[sel]
            mu, sd = float(np.mean(sub)), float(np.std(sub) + 1e-12)
            score[sel] = (sub - mu) / sd
    order = np.argsort(-np.where(valid, score, -np.inf))
    n_top = int(max(50, min(n_top, int(valid.sum()))))
    chosen = [int(i) for i in order[:n_top] if valid[i]]
    table_rows = []
    for idx in order[:min(200, len(order))]:
        if not valid[idx]:
            continue
        table_rows.append({"gene": genes[idx], "mean": float(mean[idx]), "dispersion": float(disp[idx]),
                           "score": float(score[idx]), "pct": float(nonzero[idx] * 100)})
    fig, ax = new_fig(6.4, 4.6)
    sel = np.zeros(X.shape[1], dtype=bool)
    sel[chosen] = True
    ax.scatter(mean[~sel], np.log1p(disp[~sel]), s=5, color="#b8bcc4", linewidths=0, label="other")
    ax.scatter(mean[sel], np.log1p(disp[sel]), s=7, color=PALETTE[2], linewidths=0, label="HVG")
    ax.set_xscale("log")
    ax.set_xlabel("mean expression", fontsize=9)
    ax.set_ylabel("log1p(dispersion)", fontsize=9)
    ax.set_title(f"Highly variable genes (n={len(chosen)})", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    figure = save_fig(fig, outdir, "hvg-scatter", title=None)
    return result(
        summary={"title": "高变基因", "bullets": [
            f"选出 {len(chosen)} 个高变基因（flavor={flavor}）",
            "下一步 PCA 只在这些基因上计算（避免噪声基因主导）"],
            "metrics": {"n_hvg": len(chosen), "flavor": flavor}},
        tables=[table("hvg", "高变基因（前 200）",
                      columns=[{"key": "gene", "label": "基因", "type": "string"},
                               {"key": "mean", "label": "平均表达", "type": "number"},
                               {"key": "dispersion", "label": "dispersion", "type": "number"},
                               {"key": "score", "label": "标准化分数", "type": "number"},
                               {"key": "pct", "label": "检出率%", "type": "number"}],
                      rows=table_rows)],
        figures=[figure],
        store={"hvg": chosen, "hvg_genes": [genes[i] for i in chosen]},
    )


# ---------------------------------------------------------------------------
# 步骤：pca
# ---------------------------------------------------------------------------


def _hvg_indices(ctx) -> np.ndarray:
    store = _store(ctx)
    hvg = store.get("hvg")
    if isinstance(hvg, (list, tuple, np.ndarray)) and len(hvg) > 0:
        idx = np.asarray([int(item) for item in hvg if 0 <= int(item) < _X(ctx).shape[1]], dtype=int)
        if idx.size:
            return idx
    hvg_genes = store.get("hvg_genes") or []
    if hvg_genes:
        index = {gene: i for i, gene in enumerate(_genes(ctx))}
        idx = np.asarray([index[g] for g in hvg_genes if g in index], dtype=int)
        if idx.size:
            return idx
    return np.arange(_X(ctx).shape[1])


def _scaled_matrix(X: np.ndarray, idx: np.ndarray, scale: bool = True) -> np.ndarray:
    sub = np.asarray(X[:, idx], dtype=np.float64)
    if scale:
        mu = sub.mean(axis=0)
        sd = sub.std(axis=0)
        sd[sd <= 1e-9] = 1.0
        sub = (sub - mu) / sd
    return sub


def _pca_embedding(X: np.ndarray, n_pcs: int) -> tuple[np.ndarray, np.ndarray]:
    """返回 (scores, explained_variance_ratio)。"""
    arr = np.asarray(X, dtype=np.float64)
    n_pcs = int(max(2, min(n_pcs, min(arr.shape) - 1)))
    try:
        from sklearn.decomposition import PCA

        model = PCA(n_components=n_pcs, random_state=0)
        scores = model.fit_transform(arr)
        return np.asarray(scores, dtype=np.float32), np.asarray(model.explained_variance_ratio_, dtype=float)
    except Exception:  # noqa: BLE001 - sklearn 不可用时退回 SVD
        arr = arr - arr.mean(axis=0, keepdims=True)
        u, s, _vt = np.linalg.svd(arr, full_matrices=False)
        scores = u[:, :n_pcs] * s[:n_pcs]
        ratio = (s ** 2) / np.maximum((s ** 2).sum(), 1e-12)
        return np.asarray(scores, dtype=np.float32), np.asarray(ratio[:n_pcs], dtype=float)


def _pca_cached(ctx, n_pcs: int = 30, scale: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """取（或计算）PCA 坐标：优先 results["pca"]，其次现场算。"""
    store = _store(ctx).get("pca")
    if isinstance(store, dict):
        coords = store.get("coords")
        ratio = store.get("ratio")
        if isinstance(coords, np.ndarray) and coords.ndim == 2 and coords.shape[0] == _X(ctx).shape[0]:
            return np.asarray(coords, dtype=np.float32), np.asarray(ratio if ratio is not None else [], dtype=float)
    idx = _hvg_indices(ctx)
    mat = _scaled_matrix(_X(ctx), idx, scale=scale)
    return _pca_embedding(mat, n_pcs)


def step_pca(ctx, params, outdir):
    X = _X(ctx)
    n_pcs = _int_params(params, "n_pcs", 30)
    scale = _bool_params(params, "scale", True)
    if str(ctx.get("X_kind")) == "counts":
        raise prerequisite("请先跑标准化（PCA 需要 lognorm 表达值）")
    idx = _hvg_indices(ctx)
    mat = _scaled_matrix(X, idx, scale=scale)
    scores, ratio = _pca_embedding(mat, n_pcs)
    fig, ax = new_fig(6.2, 4.2)
    xs = np.arange(1, len(ratio) + 1)
    ax.bar(xs, ratio * 100, color=PALETTE[0], width=0.72)
    ax.plot(xs, np.cumsum(ratio) * 100, color=PALETTE[2], marker="o", markersize=3, linewidth=1.2)
    ax.set_xlabel("PC", fontsize=9)
    ax.set_ylabel("% variance", fontsize=9)
    ax.set_title("PCA scree / cumulative variance", fontsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    figure = save_fig(fig, outdir, "pca-scree", title=None)
    bullets = [
        f"在 {len(idx)} 个基因上做 PCA（scale={scale}），保留 {scores.shape[1]} 个主成分",
        f"PC1 解释 {ratio[0] * 100:.1f}%，前 5 个累计 {np.cumsum(ratio)[:5].sum() * 100:.1f}%",
    ]
    return result(
        summary={"title": "PCA", "bullets": bullets,
                 "metrics": {"n_pcs": int(scores.shape[1]), "pc1": float(ratio[0]),
                             "cum5": float(np.cumsum(ratio)[:5].sum())}},
        figures=[figure],
        store={"pca": {"coords": scores, "ratio": ratio}},
        coords_store={"pca": scores},
    )


# ---------------------------------------------------------------------------
# 步骤：neighbors（表达 kNN + 空间图融合）
# ---------------------------------------------------------------------------


def _knn_graph(X: np.ndarray, n_neighbors: int, metric: str = "cosine"):
    """kNN 邻接（对称化、权重=相似度），返回 scipy 稀疏矩阵（行归一化前的二值/权重图）。"""
    from scipy import sparse

    arr = np.asarray(X, dtype=np.float64)
    n = arr.shape[0]
    k = int(max(1, min(n_neighbors, n - 1)))
    if metric == "cosine":
        norm = np.linalg.norm(arr, axis=1, keepdims=True)
        norm[norm <= 1e-12] = 1.0
        unit = arr / norm
        sim = unit @ unit.T
    elif metric == "euclidean":
        sq = (arr ** 2).sum(axis=1)
        sim = -(sq[:, None] + sq[None, :] - 2 * arr @ arr.T)
    else:
        sim = arr @ arr.T
    np.fill_diagonal(sim, -np.inf)
    idx = np.argpartition(-sim, kth=k - 1, axis=1)[:, :k]
    rows = np.repeat(np.arange(n), k)
    cols = idx.ravel()
    weights = sim[rows, cols]
    weights = np.clip(weights, 0, None)
    graph = sparse.csr_matrix((weights, (rows, cols)), shape=(n, n))
    graph = graph.maximum(graph.T)
    return graph


def _normalize_graph(graph):
    from scipy import sparse

    g = graph.tocsr().astype(np.float64).copy()
    g.setdiag(0.0)
    g.eliminate_zeros()
    row = np.asarray(g.sum(axis=1)).ravel()
    row[row <= 0] = 1.0
    return sparse.diags(1.0 / row) @ g


def _fused_graph(ctx, params) -> object:
    n_neighbors = _int_params(params, "n_neighbors", 15)
    metric = str((params or {}).get("metric") or "cosine").strip().lower()
    weight = _numeric_params(params, "spatial_weight", 0.3)
    weight = float(np.clip(weight if weight is not None else 0.3, 0.0, 0.95))
    scores, _ratio = _pca_cached(ctx)
    expr_graph = _normalize_graph(_knn_graph(scores, n_neighbors, metric if metric in ("cosine", "euclidean") else "cosine"))
    if weight <= 0:
        return expr_graph, weight, "仅表达图"
    coords = _coords(ctx)
    radius = ctx.get("neighbor_radius") or estimate_radius(coords)
    space = _normalize_graph(radius_graph(coords, radius, k_fallback=max(3, n_neighbors // 2)))
    fused = (1.0 - weight) * expr_graph + weight * space
    return _normalize_graph(fused), weight, f"表达图 {1 - weight:.2f} + 空间图 {weight:.2f}"


def step_neighbors(ctx, params, outdir):
    n_neighbors = _int_params(params, "n_neighbors", 15)
    graph, weight, desc = _fused_graph(ctx, params)
    n = graph.shape[0]
    # 行归一化后行和为 1，度数要用非零元素个数
    degree = np.diff(graph.tocsr().indptr).astype(float)
    try:
        from scipy.sparse.csgraph import connected_components

        n_comp, labels = connected_components(graph, directed=False)
    except Exception:  # noqa: BLE001
        n_comp, labels = 1, np.zeros(n, dtype=int)
    bullets = [
        f"邻接图：{n} 个 spot，平均邻居 {degree.mean():.1f}，连通分量 {n_comp}",
        f"融合方式：{desc}（n_neighbors={n_neighbors}）",
    ]
    if n_comp > 1:
        bullets.append(f"注意：图不连通（{n_comp} 个分量），聚类可能被切成互不相连的块")
    figure = _fig_spatial_metric(ctx, outdir, np.asarray(labels, dtype=float), "neighbors-components",
                                f"Connected components (n={n_comp})", "component")
    return result(
        summary={"title": "邻接图", "bullets": bullets,
                 "metrics": {"n_neighbors": n_neighbors, "mean_degree": float(degree.mean()),
                             "components": int(n_comp), "spatial_weight": float(weight)}},
        figures=[figure],
        store={"graph": {"n_neighbors": n_neighbors, "spatial_weight": float(weight),
                         "components": int(n_comp)}},
    )


# ---------------------------------------------------------------------------
# 步骤：cluster（Louvain）
# ---------------------------------------------------------------------------


def _louvain(graph, resolution: float = 1.0) -> np.ndarray:
    """图上的 Louvain 社区检测（networkx）；不可用时退回谱聚类。"""
    n = graph.shape[0]
    try:
        import networkx as nx
        from scipy import sparse

        g = graph.tocoo()
        net = nx.Graph()
        net.add_nodes_from(range(n))
        for i, j, w in zip(g.row, g.col, g.data):
            if i < j and w > 0:
                net.add_edge(int(i), int(j), weight=float(w) + 1e-9)
        communities = nx.community.louvain_communities(net, weight="weight",
                                                      resolution=float(max(0.05, resolution)), seed=0)
        labels = np.full(n, -1, dtype=int)
        for index, group in enumerate(sorted(communities, key=lambda item: (-len(item), min(item)))):
            for node in group:
                labels[int(node)] = index
        # 孤立点（没有任何边）单独归一簇
        isolated = np.where(labels < 0)[0]
        if isolated.size:
            labels[isolated] = labels.max() + 1
        return labels
    except Exception:  # noqa: BLE001 - networkx 缺失或失败时退回谱聚类
        from sklearn.cluster import KMeans

        scores, _ratio = _pca_embedding(np.asarray(graph.toarray(), dtype=np.float64), 12)
        k = int(max(2, min(12, round(2 + 3 * resolution))))
        return KMeans(n_clusters=k, random_state=0, n_init=10).fit_predict(scores)


def _graph_summary(graph) -> dict:
    return {"shape": list(graph.shape)}


def step_cluster(ctx, params, outdir):
    entry = _store(ctx).get("graph")
    if not entry:
        # 没跑 neighbors 也允许：现场按默认参数构造融合图，避免卡住用户
        graph, _weight, _desc = _fused_graph(ctx, {"n_neighbors": 15, "spatial_weight": 0.3})
    else:
        graph, _weight, _desc = _fused_graph(ctx, {"n_neighbors": entry.get("n_neighbors", 15),
                                                  "spatial_weight": entry.get("spatial_weight", 0.3)})
    resolution = _numeric_params(params, "resolution", 1.0) or 1.0
    labels = _louvain(graph, float(resolution))
    labels_str = [str(int(item)) for item in labels]
    counts = pd.Series(labels_str).value_counts()
    order = _sorted_labels(counts.index.tolist())
    bullets = [
        f"Louvain 聚类（resolution={result_fmt(resolution)}）得到 {len(order)} 个簇",
        "簇大小：" + "，".join(f"{key}:{int(counts[key])}" for key in order[:10]),
    ]
    fig, ax = new_fig(6.0, 5.0)
    scatter_spatial(ax, ctx, labels=labels_str, point_size=7.0, legend_title="cluster")
    ax.set_title("Louvain clusters", fontsize=10)
    figure_spatial = save_fig(fig, outdir, "cluster-spatial", title=None)
    fig2, ax2 = new_fig(6.4, 3.8)
    positions = np.arange(len(order))
    ax2.bar(positions, [int(counts[key]) for key in order],
            color=[categorical_colors(order).get(key, PALETTE[0]) for key in order], width=0.7)
    ax2.set_xticks(positions)
    ax2.set_xticklabels(order, fontsize=8)
    ax2.set_ylabel("spots", fontsize=9)
    ax2.set_title("Cluster sizes", fontsize=10)
    for side in ("top", "right"):
        ax2.spines[side].set_visible(False)
    figure_bar = save_fig(fig2, outdir, "cluster-sizes", title=None)
    rows = [{"cluster": key, "n_spots": int(counts[key]),
             "pct": float(counts[key] / max(len(labels), 1) * 100)} for key in order]
    return result(
        summary={"title": "聚类", "bullets": bullets,
                 "metrics": {"n_clusters": len(order), "resolution": float(resolution)}},
        tables=[table("clusters", "簇大小",
                      columns=[{"key": "cluster", "label": "簇", "type": "string"},
                               {"key": "n_spots", "label": "spot 数", "type": "number"},
                               {"key": "pct", "label": "占比%", "type": "number"}],
                      rows=rows)],
        figures=[figure_spatial, figure_bar],
        cluster_store={"key": "clusters", "labels": labels_str},
    )


# ---------------------------------------------------------------------------
# 步骤：embed
# ---------------------------------------------------------------------------


def step_embed(ctx, params, outdir):
    layout = str((params or {}).get("layout") or "tsne").strip().lower()
    perplexity = _numeric_params(params, "perplexity", 30.0) or 30.0
    scores, _ratio = _pca_cached(ctx)
    n = scores.shape[0]
    limit = 20000
    index = _subsample(n, limit, seed=0)
    sub = scores[index]
    if layout == "pca":
        emb = sub[:, :2]
        layout = "pca"
    elif layout in ("tsne", "t-sne"):
        try:
            from sklearn.manifold import TSNE

            perp = float(min(max(5.0, perplexity), max(5.0, (sub.shape[0] - 1) / 3.0)))
            emb = TSNE(n_components=2, perplexity=perp, random_state=0, init="pca",
                       learning_rate="auto").fit_transform(sub)
            layout = "tsne"
        except Exception as error:  # noqa: BLE001
            raise prerequisite(f"t-SNE 失败（{error}）；可改用 layout=pca")
    else:
        emb = sub[:, :2]
        layout = "pca"
    coords_full = np.zeros((n, 2), dtype=np.float32)
    coords_full[index] = np.asarray(emb, dtype=np.float32)
    labels = _column(ctx, ctx.get("cluster_key") or "clusters")
    fig, ax = new_fig(6.0, 5.0)
    if labels is not None:
        labels = np.asarray(labels)[index]
        mapping = categorical_colors(_sorted_labels(labels))
        for value in mapping:
            sel = np.array([str(item) == value for item in labels])
            ax.scatter(coords_full[index][sel, 0], coords_full[index][sel, 1], s=6,
                       color=mapping[value], linewidths=0, label=value)
        ax.legend(fontsize=7, frameon=False, markerscale=2)
    else:
        ax.scatter(coords_full[index, 0], coords_full[index, 1], s=6, color=PALETTE[0], linewidths=0)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{layout.upper()} embedding", fontsize=10)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    figure = save_fig(fig, outdir, f"embedding-{layout}", title=None)
    bullets = [f"{layout.upper()} 二维嵌入完成（{len(index)} 个 spot"
               + ("，已降采样" if len(index) < n else "") + "）"]
    return result(
        summary={"title": "二维嵌入", "bullets": bullets,
                 "metrics": {"layout": layout, "n_points": int(len(index)),
                             "downsampled": bool(len(index) < n)}},
        figures=[figure],
        coords_store={layout: coords_full},
    )


# ---------------------------------------------------------------------------
# 步骤：annotate
# ---------------------------------------------------------------------------


def step_annotate(ctx, params, outdir):
    mapping = (params or {}).get("mapping") or {}
    key = str((params or {}).get("key") or "domain_anno").strip() or "domain_anno"
    base_key = str((params or {}).get("cluster_key") or ctx.get("cluster_key") or "clusters")
    if not isinstance(mapping, dict) or not mapping:
        raise prerequisite(
            f'请提供 mapping，例如 {{"0":"Tumor","1":"Stroma"}}（当前聚类列：{base_key}）')
    labels = _column(ctx, base_key)
    if labels is None:
        labels = ctx.get("cluster")
    if labels is None:
        raise prerequisite(f"找不到聚类列「{base_key}」：请先跑聚类")
    labels = [str(item) for item in labels]
    mapped = [str(mapping.get(item, mapping.get(str(item), ""))) or item for item in labels]
    values = mapped
    missing = sorted({item for item in labels if str(item) not in mapping})
    bullets = [f"已把 {len(set(labels))} 个簇映射为 {len(set(values))} 种注释（列名 {key}）"]
    if missing:
        bullets.append("未映射的簇保留原编号：" + "、".join(missing[:10]))
    fig, ax = new_fig(6.0, 5.0)
    scatter_spatial(ax, ctx, labels=values, point_size=7.0, legend_title=key)
    ax.set_title(f"Annotation ({key})", fontsize=10)
    figure = save_fig(fig, outdir, "annotate-map", title=None)
    return result(
        summary={"title": "区域注释", "bullets": bullets,
                 "metrics": {"n_labels": len(set(values)), "key": key}},
        figures=[figure],
        sample_columns=[{"key": key, "kind": "categorical"}],
        sample_values={key: values},
    )


# ---------------------------------------------------------------------------
# 分析：质控与整合
# ---------------------------------------------------------------------------


def analysis_qc_summary(ctx, params, outdir):
    X = _X(ctx)
    genes = _genes(ctx)
    metrics = _spot_metrics(X, genes)
    info = _info(ctx)
    rows = []
    for span in ("all",):
        rows.append({
            "span": span,
            "n_spots": int(X.shape[0]),
            "median_counts": float(np.median(metrics["counts"])),
            "min_counts": float(np.min(metrics["counts"])),
            "median_genes": float(np.median(metrics["genes"])),
            "median_pct_mt": float(np.median(metrics["pct_mt"])),
            "zero_spots": int((metrics["counts"] <= 0).sum()),
        })
    figures = [
        _fig_spatial_metric(ctx, outdir, metrics["counts"].to_numpy(), "qcsummary-counts",
                            "Counts per spot", "counts"),
        _fig_spatial_metric(ctx, outdir, metrics["genes"].to_numpy(), "qcsummary-genes",
                            "Detected genes per spot", "genes"),
    ]
    bullets = [f"{X.shape[0]} spot × {X.shape[1]} 基因；零表达 spot {int((metrics['counts'] <= 0).sum())} 个",
               f"中位 counts {np.median(metrics['counts']):.0f}，中位检出基因 {np.median(metrics['genes']):.0f}"]
    return result(summary={"title": "质控汇总", "bullets": bullets,
                           "metrics": {"n_spots": int(X.shape[0]), "n_genes": int(X.shape[1])}},
                  tables=[table("qcs", "质控汇总",
                                columns=[{"key": "span", "label": "范围", "type": "string"},
                                         {"key": "n_spots", "label": "spot", "type": "number"},
                                         {"key": "median_counts", "label": "counts 中位", "type": "number"},
                                         {"key": "min_counts", "label": "counts 最小", "type": "number"},
                                         {"key": "median_genes", "label": "基因中位", "type": "number"},
                                         {"key": "median_pct_mt", "label": "%MT 中位", "type": "number"},
                                         {"key": "zero_spots", "label": "零表达", "type": "number"}],
                                rows=rows)],
                  figures=figures)


def analysis_outlier_spots(ctx, params, outdir):
    k = _int_params(params, "k", 8)
    z_cut = _numeric_params(params, "z", 3.0) or 3.0
    scores, _ratio = _pca_cached(ctx)
    coords = _coords(ctx)
    idx, dist = coord_knn(coords, k)
    if idx.size == 0:
        raise prerequisite("spot 太少，无法做邻域平滑")
    sub = scores[:, : min(15, scores.shape[1])]
    neigh_mean = sub[idx].mean(axis=1)
    residual = np.linalg.norm(sub - neigh_mean, axis=1)
    z = (residual - np.mean(residual)) / (np.std(residual) + 1e-12)
    flagged = np.where(np.abs(z) >= float(z_cut))[0]
    bullets = [f"邻域表达残差 z≥{result_fmt(z_cut)} 的离群 spot：{len(flagged)} 个（k={k}）"]
    if len(flagged):
        bullets.append("这些 spot 可能是组织边缘、褶皱或低质量区域；建议看空间图确认")
    fig, ax = new_fig(6.2, 5.0)
    scatter_spatial(ax, ctx, values=z, point_size=6.0, legend_title="z-score")
    if len(flagged):
        ax.scatter(coords[flagged, 0], coords[flagged, 1], s=16, facecolors="none",
                   edgecolors="#d1495b", linewidths=0.7, label="outlier")
        ax.legend(fontsize=7, frameon=False)
    ax.set_title("Neighborhood expression residual", fontsize=10)
    figure = save_fig(fig, outdir, "outlier-spots", title=None)
    rows = [{"spot": _spots(ctx)[i] if i < len(_spots(ctx)) else str(i),
             "residual": float(residual[i]), "z": float(z[i])}
            for i in np.argsort(-np.abs(z))[:200]]
    return result(summary={"title": "空间离群点", "bullets": bullets,
                           "metrics": {"n_outliers": int(len(flagged)), "k": int(k)}},
                  tables=[table("outliers", "残差最大的 200 个 spot",
                                columns=[{"key": "spot", "label": "spot", "type": "string"},
                                         {"key": "residual", "label": "残差", "type": "number"},
                                         {"key": "z", "label": "z", "type": "number"}],
                                rows=rows)],
                  figures=[figure])


def _slice_values(ctx, column: str) -> np.ndarray:
    values = _column(ctx, column)
    if values is None:
        raise prerequisite(f"找不到切片列「{column}」：请在样本表里加一列（如 slice/sample）并设为切片列")
    return np.asarray([str(item) for item in values])


def analysis_slice_align(ctx, params, outdir):
    column = str((params or {}).get("slice_col") or "slice")
    slices = _slice_values(ctx, column)
    scores, _ratio = _pca_cached(ctx)
    uniq = _sorted_labels(slices)
    if len(uniq) < 2:
        raise prerequisite(f"切片列「{column}」只有 {len(uniq)} 个取值：至少需要 2 个切片才能做对齐诊断")
    rows = []
    for name in uniq:
        sel = np.array([item == name for item in slices])
        if sel.sum() == 0:
            continue
        sub = scores[sel][:, : min(5, scores.shape[1])]
        row = {"slice": name, "n_spots": int(sel.sum())}
        for pc in range(sub.shape[1]):
            row[f"PC{pc + 1}_mean"] = float(np.mean(sub[:, pc]))
        rows.append(row)
    fig, ax = new_fig(6.4, 4.6)
    mapping = categorical_colors(uniq)
    for name in uniq:
        sel = np.array([item == name for item in slices])
        ax.scatter(scores[sel, 0], scores[sel, 1], s=8, color=mapping[name], linewidths=0,
                   label=name, alpha=0.85)
    ax.set_xlabel("PC1", fontsize=9)
    ax.set_ylabel("PC2", fontsize=9)
    ax.set_title("Per-slice PCA overlap", fontsize=10)
    ax.legend(fontsize=8, frameon=False, markerscale=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    figure = save_fig(fig, outdir, "slice-align", title=None)
    centroids = np.array([[row.get(f"PC{pc + 1}_mean", 0.0) for pc in range(min(5, scores.shape[1]))]
                          for row in rows], dtype=float)
    spread = float(np.mean(np.linalg.norm(centroids - centroids.mean(axis=0), axis=1))) if len(rows) > 1 else 0.0
    bullets = [f"{len(uniq)} 个切片，各切片 PCA 重心平均偏移 {spread:.2f}（越小越一致）"]
    if spread > 1.0:
        bullets.append("切片间存在明显整体偏移：建议后续按切片分别标准化，或改用整合方法")
    return result(summary={"title": "多切片对齐诊断", "bullets": bullets,
                           "metrics": {"n_slices": len(uniq), "centroid_spread": spread}},
                  tables=[table("slices", "各切片 PCA 重心",
                                columns=[{"key": "slice", "label": "切片", "type": "string"},
                                         {"key": "n_spots", "label": "spot", "type": "number"}]
                                        + [{"key": f"PC{i}_mean", "label": f"PC{i} 均值", "type": "number"}
                                           for i in range(1, 6)],
                                rows=rows)],
                  figures=[figure])


def analysis_spatial_trajectory(ctx, params, outdir):
    start_gene = str((params or {}).get("start_gene") or "").strip()
    start_spot = _numeric_params(params, "start_spot", None)
    scores, _ratio = _pca_cached(ctx)
    sub = scores[:, : min(15, scores.shape[1])]
    try:
        graph, _weight, _desc = _fused_graph(ctx, {"n_neighbors": 15, "spatial_weight": 0.5})
    except EngineError:
        raise prerequisite("空间轨迹需要空间坐标")
    from scipy import sparse

    n = graph.shape[0]
    p = sparse.diags(np.ones(n)) * 0.5 + graph * 0.5
    p = _normalize_graph(p)
    try:
        from scipy.sparse.linalg import eigsh

        k = int(min(10, n - 2))
        values, vectors = eigsh(p.T, k=k, which="LM", v0=np.ones(n))
    except Exception:  # noqa: BLE001
        dense = np.asarray(p.todense())
        values, vectors = np.linalg.eig(dense.T)
        order = np.argsort(-np.abs(values))[: min(10, n)]
        values, vectors = np.real(values[order]), np.real(vectors[:, order])
    order = np.argsort(-np.abs(values))
    vectors = np.real(vectors[:, order])
    values = np.real(values[order])
    if vectors.shape[1] < 3:
        raise prerequisite("扩散图特征向量不足，无法计算伪时序")
    diffusion = vectors[:, 1:4]
    genes = _genes(ctx)
    if start_gene and start_gene in genes:
        gi = genes.index(start_gene)
        expr = np.asarray(_X(ctx)[:, gi], dtype=float)
        start = int(np.argmax(expr))
        start_desc = f"{start_gene} 表达最高的 spot (spot #{start})"
    elif start_spot is not None:
        start = int(np.clip(int(start_spot), 0, n - 1))
        start_desc = f"指定 spot #{start}"
    else:
        start = int(np.argmax(np.linalg.norm(diffusion - diffusion.mean(axis=0), axis=1)))
        start_desc = "扩散空间中最边缘的 spot（可传 start_gene 指定起点）"
    pt = np.linalg.norm(diffusion - diffusion[start], axis=1)
    pt = (pt - pt.min()) / (pt.max() - pt.min() + 1e-12)
    # 与伪时序相关性最高的基因
    X = np.asarray(_X(ctx), dtype=float)
    corr = np.zeros(X.shape[1])
    pt_c = pt - pt.mean()
    denom = np.linalg.norm(pt_c) + 1e-12
    for gi in range(X.shape[1]):
        col = X[:, gi]
        col_c = col - col.mean()
        corr[gi] = float(col_c @ pt_c / (np.linalg.norm(col_c) * denom + 1e-12))
    top = np.argsort(-np.abs(corr))[:20]
    rows = [{"gene": genes[i], "corr_with_pseudotime": float(corr[i])} for i in top]
    figure = _fig_spatial_metric(ctx, outdir, pt, "trajectory-spatial",
                                "Diffusion pseudotime", "pseudotime")
    fig2, ax2 = new_fig(6.4, 4.2)
    ax = fig2.add_subplot(111)
    bins = np.linspace(0, 1, 21)
    assign = np.clip(np.digitize(pt, bins) - 1, 0, len(bins) - 2)
    for rank, gi in enumerate(top[:4]):
        means = [float(np.mean(X[assign == b, gi])) if np.any(assign == b) else np.nan
                 for b in range(len(bins) - 1)]
        ax.plot(bins[:-1], means, marker="o", markersize=3, linewidth=1.2,
                color=PALETTE[rank % len(PALETTE)], label=genes[gi])
    ax.set_xlabel("pseudotime", fontsize=9)
    ax.set_ylabel("mean expression", fontsize=9)
    ax.set_title("Top genes along pseudotime", fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    figure2 = save_fig(fig2, outdir, "trajectory-genes", title=None)
    bullets = [f"扩散伪时序起点：{start_desc}",
               f"与伪时序最相关的基因：{'、'.join(genes[i] for i in top[:5])}"]
    return result(summary={"title": "空间轨迹", "bullets": bullets,
                           "metrics": {"start": start_desc, "n_spots": int(n)}},
                  tables=[table("trajectory", "与伪时序相关性最高的基因",
                                columns=[{"key": "gene", "label": "基因", "type": "string"},
                                         {"key": "corr_with_pseudotime", "label": "相关系数", "type": "number"}],
                                rows=rows)],
                  figures=[figure, figure2],
                  store={"trajectory": {"pseudotime": pt.tolist(), "start": start_desc}},
                  coords_store={"pseudotime": np.column_stack([pt, np.zeros_like(pt)])})


# ---------------------------------------------------------------------------
# 注册
# ---------------------------------------------------------------------------

STEPS = [
    {
        "id": "qc", "title": "空间质控", "icon": "qc", "requires": [],
        "params": [
            {"key": "min_counts", "label": "最少 counts", "type": "number", "default": None},
            {"key": "min_genes", "label": "最少检出基因", "type": "number", "default": None},
            {"key": "max_pct_mt", "label": "线粒体比例上限 %", "type": "number", "default": 20},
            {"key": "group_col", "label": "分组列", "type": "column", "default": None},
        ],
        "run": step_qc,
    },
    {
        "id": "filter", "title": "过滤", "icon": "filter", "requires": ["qc"],
        "params": [
            {"key": "min_counts", "label": "最少 counts", "type": "number", "default": 500},
            {"key": "min_genes", "label": "最少检出基因", "type": "number", "default": 200},
            {"key": "max_pct_mt", "label": "线粒体比例上限 %", "type": "number", "default": 20},
            {"key": "min_cells", "label": "基因最少检出 spot 数", "type": "number", "default": 3},
        ],
        "run": step_filter,
    },
    {
        "id": "normalize", "title": "标准化", "icon": "normalize", "requires": ["filter"],
        "params": [
            {"key": "target_sum", "label": "目标总量", "type": "number", "default": 10000},
            {"key": "log", "label": "log1p 变换", "type": "boolean", "default": True},
            {"key": "method", "label": "方法", "type": "select", "default": "lognorm",
             "options": [{"value": "lognorm", "label": "文库大小标准化 + log1p"},
                         {"value": "none", "label": "跳过"}]},
        ],
        "run": step_normalize,
    },
    {
        "id": "hvg", "title": "高变基因", "icon": "hvg", "requires": ["normalize"],
        "params": [
            {"key": "n_top_genes", "label": "高变基因数", "type": "number", "default": 2000},
            {"key": "flavor", "label": "算法", "type": "select", "default": "seurat",
             "options": [{"value": "seurat", "label": "Seurat v3 分箱 dispersion"},
                         {"value": "cell_ranger", "label": "Cell Ranger 归一化 dispersion"}]},
        ],
        "run": step_hvg,
    },
    {
        "id": "pca", "title": "PCA", "icon": "pca", "requires": ["hvg"],
        "params": [
            {"key": "n_pcs", "label": "主成分数", "type": "number", "default": 30},
            {"key": "scale", "label": "基因 z-score", "type": "boolean", "default": True},
        ],
        "run": step_pca,
    },
    {
        "id": "neighbors", "title": "邻接图", "icon": "neighbors", "requires": ["pca"],
        "params": [
            {"key": "n_neighbors", "label": "邻居数", "type": "number", "default": 15},
            {"key": "metric", "label": "表达距离", "type": "select", "default": "cosine",
             "options": [{"value": "cosine", "label": "cosine"}, {"value": "euclidean", "label": "euclidean"}]},
            {"key": "spatial_weight", "label": "空间图权重", "type": "number", "default": 0.3},
            {"key": "spatial_method", "label": "空间邻接", "type": "select", "default": "grid",
             "options": [{"value": "grid", "label": "半径邻接（自适应）"},
                         {"value": "knn", "label": "k 近邻"}]},
        ],
        "run": step_neighbors,
    },
    {
        "id": "cluster", "title": "聚类", "icon": "cluster", "requires": ["neighbors"],
        "params": [
            {"key": "resolution", "label": "resolution", "type": "number", "default": 1.0},
            {"key": "method", "label": "算法", "type": "select", "default": "louvain",
             "options": [{"value": "louvain", "label": "Louvain"}]},
        ],
        "run": step_cluster,
    },
    {
        "id": "embed", "title": "二维嵌入", "icon": "embed", "requires": ["pca"],
        "params": [
            {"key": "layout", "label": "布局", "type": "select", "default": "tsne",
             "options": [{"value": "tsne", "label": "t-SNE"}, {"value": "pca", "label": "PCA"}]},
            {"key": "perplexity", "label": "perplexity", "type": "number", "default": 30},
        ],
        "run": step_embed,
    },
    {
        "id": "annotate", "title": "区域注释", "icon": "annotate", "requires": ["cluster"],
        "params": [
            {"key": "mapping", "label": "簇→注释映射", "type": "text", "default": ""},
            {"key": "key", "label": "写入列名", "type": "text", "default": "domain_anno"},
        ],
        "run": step_annotate,
    },
]

ANALYSES = [
    {
        "id": "qc_summary", "title": "质控汇总", "group": "质控与整合",
        "method": "spot 级文库/检出/线粒体统计（内置实现）",
        "needs": ["coords"],
        "params": [],
        "run": analysis_qc_summary,
    },
    {
        "id": "outlier_spots", "title": "空间离群点", "group": "质控与整合",
        "method": "邻域表达残差 z-score（kNN 平滑，内置实现）",
        "needs": ["cluster"],
        "params": [
            {"key": "k", "label": "邻域大小", "type": "number", "default": 8},
            {"key": "z", "label": "z 阈值", "type": "number", "default": 3},
        ],
        "run": analysis_outlier_spots,
    },
    {
        "id": "slice_align", "title": "多切片对齐诊断", "group": "质控与整合",
        "method": "各切片 PCA 重心/重叠诊断（内置实现）",
        "needs": ["normalize"],
        "params": [
            {"key": "slice_col", "label": "切片列", "type": "column", "default": "slice"},
        ],
        "run": analysis_slice_align,
    },
    {
        "id": "spatial_trajectory", "title": "空间轨迹（扩散伪时序）", "group": "质控与整合",
        "method": "扩散图 + 伪时序（Haghverdi 2016 DPT 思想；内置实现）",
        "needs": ["cluster"],
        "params": [
            {"key": "start_gene", "label": "起点基因（可选）", "type": "text", "default": ""},
            {"key": "start_spot", "label": "起点 spot 序号（可选）", "type": "number", "default": None},
        ],
        "run": analysis_spatial_trajectory,
    },
]

__all__ = ["STEPS", "ANALYSES"]
