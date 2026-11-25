"""
dsh-spatial — 反卷积（DESIGN.md §3 `deconv` 步骤 + §4 group=反卷积）

思路（文献）：
  * 非负最小二乘 / 阻尼加权最小二乘的 spot 反卷积（Tsoucas 2019 DWLS；SPOTlight 2021；
    Sturm 2019 对 Visium 反卷积方法的评测）。
  * RCTD（Cable et al. 2022 Nat Biotechnol）走 R 包 spacexr 的官方实现：
    spacexr::Reference → spacexr::SpatialRNA → create.RCTD → run.RCTD，
    由 lib/r_rctd.R 执行、lib/sp_rbridge.py 负责 bundle/结果往返。这是**可选**路径
    （需要本机有 Rscript + spacexr），缺失时给 unsupported 并说明安装方式；
    nnls/dwls 是始终可用的内置实现（不经 R）。
  * 参考谱来自单细胞 h5ad（ctx["ref"] 懒加载）：每类细胞的平均表达谱（signature 矩阵）；
    RCTD 需要该参考的**原始 counts**（h5ad 的 layers["counts"] 或 raw/整数矩阵）。

本模块只依赖 numpy/scipy/matplotlib + sp_common（RCTD 路径额外用 sp_rbridge 调 R），
不依赖 scanpy/squidpy。
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

import sp_common as sc
import sp_rbridge as rbridge

_GENERIC_TYPES = ("cell_type", "celltype", "cell_type_l1", "annotation", "clusters",
                  "leiden", "seurat_clusters", "louvain")

#: RCTD 的 doublet_mode 取值（spacexr 官方）
_RCTD_MODES = ("full", "doublet", "multi")

#: log 值上界：非整数且 max < 该值 → 判为 log/标准化矩阵（RCTD 不接受）
_RCTD_LOG_MAX = 50.0

#: 从 params 透传给 spacexr::create.RCTD 的可选键（不给就用 spacexr 自己的默认值）
_RCTD_OPTIONAL_KEYS = ("umi_min", "umi_max", "counts_min", "cell_min", "max_multi_types",
                       "gene_cutoff", "fc_cutoff", "gene_cutoff_reg", "fc_cutoff_reg",
                       "confidence_threshold", "doublet_threshold")



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


def _p_flag(params, key, default=None):
    """布尔参数：兼容 true/false、"true"/"false"、1/0（前端可能传字符串）。"""
    value = _p(params, key, None)
    if value is None:
        return default
    if isinstance(value, str):
        low = value.strip().lower()
        if low in ("1", "true", "yes", "y", "on"):
            return True
        if low in ("0", "false", "no", "n", "off"):
            return False
    return bool(value)


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


# ---------------------------------------------------------------------------
# RCTD（spacexr 官方实现）：counts 形态检查 + 参考/空间端的原始 counts 读取
# ---------------------------------------------------------------------------


def _counts_sniff(mat, label: str) -> dict:
    """判断一个矩阵能不能当 RCTD 的原始 counts（非负、整数）。

    返回 {"ok": bool, "reason": str, "nonint_frac": float, "max": float, "min": float}。
    抽样判断（最多 20 万个元素），避免大矩阵白算一遍；**稀疏矩阵不 densify**
    （1M 细胞 × 30k 基因的 h5ad 稠密化会直接 OOM），只看非零值的整数性/非负性。
    """
    data = getattr(mat, "data", None)
    is_sparse = data is not None and hasattr(mat, "toarray")
    if is_sparse:
        values = np.asarray(data, dtype=np.float64)
        if values.size == 0:
            return {"ok": False, "reason": f"{label} 是空稀疏矩阵（没有非零值）",
                    "nonint_frac": 1.0, "max": 0.0, "min": 0.0}
        flat = values[:: max(1, values.size // 200_000)] if values.size > 200_000 else values
        vmax = float(np.max(flat))
        vmin = min(0.0, float(np.min(flat)))
        if not np.isfinite(flat).all():
            return {"ok": False, "reason": f"{label} 含 NA/Inf：RCTD 需要原始 counts",
                    "nonint_frac": 1.0, "max": 0.0, "min": 0.0}
    else:
        arr = np.asarray(mat, dtype=np.float64)
        if arr.size == 0:
            return {"ok": False, "reason": f"{label} 是空矩阵", "nonint_frac": 1.0,
                    "max": 0.0, "min": 0.0}
        flat = arr.ravel()
        if flat.size > 200_000:
            flat = flat[:: max(1, flat.size // 200_000)]
        if not np.isfinite(flat).all():
            return {"ok": False, "reason": f"{label} 含 NA/Inf：RCTD 需要原始 counts",
                    "nonint_frac": 1.0, "max": 0.0, "min": 0.0}
        vmax = float(np.max(flat))
        vmin = float(np.min(flat))
    nonint = float(np.mean(np.abs(flat - np.rint(flat)) > 1e-6))
    info = {"nonint_frac": nonint, "max": vmax, "min": vmin}
    if vmin < 0:
        return {"ok": False, "reason": f"{label} 含负值：RCTD 需要非负 counts", **info}
    if nonint <= 1e-6:
        return {"ok": True, "reason": "整数 counts", **info}
    hint = "（数值像 log/标准化后的结果）" if vmax < _RCTD_LOG_MAX else ""
    return {"ok": False,
            "reason": f"{label} 有 {nonint * 100:.1f}% 的值不是整数（max={vmax:.3g}）{hint}："
                      "RCTD 需要原始 counts",
            **info}


def _h5ad_ref_path(ref, params) -> str:
    """参考 h5ad 的路径：参数 ref 优先，其次 ctx["ref"]（LazyRef）的 path。"""
    path = str(_p(params, "ref", "") or "")
    if path:
        return os.path.abspath(os.path.expanduser(path))
    if hasattr(ref, "get"):
        return str(ref.get("path") or "")
    return ""


def _h5ad_count_sources(adata) -> list[tuple[str, object, list[str]]]:
    """列出 h5ad 里可能装原始 counts 的位置（按优先级）：layers['counts'] → raw.X → X。"""
    out: list[tuple[str, object, list[str]]] = []
    var_names = [str(g) for g in adata.var_names]
    layers = getattr(adata, "layers", None) or {}
    if "counts" in layers:
        out.append(("layers[\"counts\"]", layers["counts"], var_names))
    raw = getattr(adata, "raw", None)
    raw_x = getattr(raw, "X", None) if raw is not None else None
    if raw_x is not None:
        raw_genes = [str(g) for g in getattr(raw, "var_names", [])]
        if len(raw_genes) != int(np.shape(raw_x)[1] if np.ndim(raw_x) == 2 else len(var_names)):
            raw_genes = var_names[: int(np.shape(raw_x)[1])] if np.ndim(raw_x) == 2 else var_names
        out.append(("raw.X", raw_x, raw_genes))
    out.append(("X", adata.X, var_names))
    return out


def _slice_columns(mat, columns: list[int]):
    """按列取子矩阵（稀疏/稠密都支持；anndata 的 backed/层视图也能切）。"""
    try:
        return mat[:, columns]
    except Exception:  # noqa: BLE001 - 少见格式退回先 densify
        arr = np.asarray(mat.toarray() if hasattr(mat, "toarray") else mat, dtype=np.float32)
        return arr[:, columns]


def _rctd_dense_guard(n_rows: int, n_cols: int, label: str, warnings: list[str]) -> None:
    """RCTD 的 R 侧要 dense 矩阵：细胞×基因太大时先提示/拦截，避免直接 OOM。"""
    elements = int(n_rows) * int(n_cols)
    gb = elements * 8 / 1e9
    if elements > 500_000_000:
        raise sc.prerequisite(
            f"{label} {n_rows} 细胞 × {n_cols} 基因 稠密化约 {gb:.1f} GB，RCTD（spacexr）需要 dense counts 矩阵："
            "请先对参考做 HVG 筛选或下采样细胞（≤2 万）后重试")
    if elements > 200_000_000:
        warnings.append(f"{label} 稠密化约 {gb:.1f} GB（RCTD 要求 dense 矩阵）：内存紧张时先下采样参考细胞")


def _rctd_reference(ctx, params) -> dict:
    """RCTD 的参考端：原始整数 counts（细胞 × 基因）+ 每细胞的类型 + 细胞名。

    参考来源与 nnls/dwls 相同（ctx["ref"] 或 params["ref"]），但 RCTD 需要 counts：
    layers["counts"] → raw.X → X（X 必须本身像 counts，否则报 prerequisite）。
    取 counts 后会先按空间数据的基因集切列（少 densify 无关基因），再做稠密内存检查。
    """
    ref = ctx.get("ref")
    cell_type_col = str(_p(params, "cell_type_col", "") or "")
    if not cell_type_col and hasattr(ref, "get"):
        cell_type_col = str(ref.get("cell_type_col") or "")
    spot_genes = set(_genes(ctx))
    path = _h5ad_ref_path(ref, params)
    if path:
        if not os.path.exists(path):
            raise sc.prerequisite(f"参考文件不存在：{path}")
        try:
            import anndata

            adata = anndata.read_h5ad(path)
        except Exception as error:  # noqa: BLE001
            raise sc.prerequisite(f"读取参考 h5ad 失败：{error}")
        tried: list[str] = []
        for label, mat, genes in _h5ad_count_sources(adata):
            info = _counts_sniff(mat, f"参考 {label}")
            if not info["ok"]:
                tried.append(info["reason"])
                continue
            keep = [i for i, gene in enumerate(genes) if gene in spot_genes]
            if not keep:
                tried.append(f"参考 {label} 的 {len(genes)} 个基因与空间数据的基因名没有交集")
                continue
            warnings: list[str] = []
            if len(keep) < len(genes):
                warnings.append(
                    f"参考 {label} 有 {len(genes) - len(keep)} 个基因不在空间数据里，已按共同基因切掉")
                mat = _slice_columns(mat, keep)
                genes = [genes[i] for i in keep]
            _rctd_dense_guard(int(np.shape(mat)[0]), len(keep), f"参考 {label}", warnings)
            stub = {"path": path, "X": mat, "genes": genes, "obs": adata.obs,
                    "cell_type_col": cell_type_col}
            Xc, kept_genes, types = _ref_columns(stub)
            cells = [str(v) for v in adata.obs_names]
            if len(cells) != Xc.shape[0]:
                cells = [f"ref{i + 1}" for i in range(Xc.shape[0])]
            return {"counts": Xc, "genes": kept_genes, "types": types, "cells": cells,
                    "path": path, "source": label, "warnings": warnings}
        raise sc.prerequisite(
            "RCTD 需要原始 counts 参考（h5ad 的 counts 层或原始矩阵）："
            + "；".join(tried))
    # 没有路径：只能用已经 materialize 进 ctx 的参考矩阵
    if hasattr(ref, "get") and ref.get("X") is not None:
        X = ref.get("X")
        info = _counts_sniff(X, "参考 X")
        if not info["ok"]:
            raise sc.prerequisite(
                "RCTD 需要原始 counts 参考（h5ad 的 counts 层或原始矩阵）：" + info["reason"])
        Xc, genes, types = _ref_columns(ref)
        frames = ref.get("obs")
        cells = ([str(v) for v in frames.index] if frames is not None and len(frames) == Xc.shape[0]
                 else [f"ref{i + 1}" for i in range(Xc.shape[0])])
        return {"counts": Xc, "genes": genes, "types": types, "cells": cells,
                "path": str(ref.get("path") or ""), "source": "ctx.ref.X", "warnings": []}
    raise sc.prerequisite("还没有指定单细胞参考：先在「反卷积」页选一个 h5ad，或传 ref 参数")


def _rctd_spatial_counts(ctx, params) -> tuple[np.ndarray, str, list[str]]:
    """RCTD 的空间端 counts：优先 ctx["counts"]，其次 X 本身是整数 counts。

    返回 (counts[spot × gene], 来源说明, warnings)。counts 层不是原始 counts 时会先记录原因、
    再试整数 X（不静默换方法，只是换同一份数据的 counts 来源）。
    """
    spots = _spots(ctx)
    use_counts = _p_flag(params, "use_counts", None)
    warnings: list[str] = []
    reasons: list[str] = []
    counts = ctx.get("counts")
    if counts is not None and use_counts is not False:
        info = _counts_sniff(counts, "空间 counts")
        if info["ok"]:
            arr = np.asarray(counts, dtype=np.float32)
            if arr.shape[0] != len(spots):
                raise sc.prerequisite(
                    f"counts 行数（{arr.shape[0]}）与 spot 数（{len(spots)}）不一致")
            return arr, "ctx.counts", warnings
        reasons.append(info["reason"])
        if use_counts is True:
            raise sc.prerequisite(info["reason"] + "：RCTD 需要空间端原始 counts")
        warnings.append(f"{info['reason']}，改用整数表达矩阵 X 作为 counts 来源")
    elif use_counts is True:
        reasons.append("数据集没有 counts 层（params.use_counts=true）")
        raise sc.prerequisite(
            "params.use_counts=true 但没有可用的 counts：该数据集只有 log/标准化矩阵，"
            "请用 counts 方式重新导入（Visium/10x 的 counts 层）")
    X = _X(ctx)
    info = _counts_sniff(X, "空间表达矩阵 X")
    if info["ok"]:
        return np.asarray(X, dtype=np.float32), "X（整数，视为 counts）", warnings
    reasons.append(info["reason"])
    raise sc.prerequisite(
        "RCTD 需要空间端原始 counts（当前数据集只有 log/标准化矩阵）：" + "；".join(reasons)
        + "；请用 counts 方式重新导入，或改用 method=nnls/dwls")


def _rctd_method_label(mode: str, max_cores: int) -> str:
    return f"RCTD（spacexr，doublet_mode={mode}，max_cores={max_cores}）"


def _rctd_notes(where: str) -> list[str]:
    """RCTD 可用时给用户的操作提示（step 与分析共用）。"""
    try:
        if not rbridge.rctd_available():
            return []
        version = rbridge.spacexr_version()
    except Exception:  # noqa: BLE001 - 探测失败不影响主流程
        return []
    ver = f" {version}" if version else ""
    if where == "step":
        return [f"本机已装 R 包 spacexr（RCTD{ver}）：method=rctd 可得到主流 RCTD 结果"
                "（spacexr 官方实现，full/doublet/multi）"]
    return [f"本机已装 R 包 spacexr（RCTD{ver}）：可用「反卷积（RCTD）」分析（deconv_rctd）"
            "获得主流结果（spacexr 官方实现）"]


#: fit_status 取值（RCTD 后处理）：已拟合 / 缺失 / 被拒绝 / NaN 权重 / 全 0 权重
_FIT_OK = "ok"
_FIT_MISSING = "missing"
_FIT_REJECTED = "rejected"
_FIT_NAN = "nan_weights"
_FIT_ZERO = "zero_weights"


def _rctd_align_weights(weights, spots, results_df=None) -> dict:
    """把 RCTD 的 weights 表对齐到全部 spot（**保留 NA，不做均值/均匀分布填补**）。

    为什么不能补：create.RCTD 会按 UMI_min/counts_MIN 丢掉低质量 spot，doublet/multi 模式还会把
    「reject」的 spot 标出来；这些 spot 本身没有可靠的比例估计。用其余 spot 的平均值 / 均匀分布
    填补会把「未观测」伪装成「已观测」，下游的组成差异、空间图和均值全部被污染，覆盖率也被掩盖。

    约定：
      * ``proportions``：n_spots × n_types，未拟合/被拒绝/全 0 的行整行 NaN；
      * ``fit_status``：ok / missing / rejected / nan_weights / zero_weights；
      * ``spot_class``：RCTD 自带的 results_df.spot_class（没有则 "NA"）；
      * ``coverage``：n_spots / n_fitted / fraction_fitted 等覆盖率信息。
    """
    n = len(spots)
    if not isinstance(weights, pd.DataFrame) or weights.empty:
        raise sc.EngineError(
            "runtime", "RCTD 没有返回权重表（tables/weights.csv 缺失或为空）")
    cell_types = [str(col) for col in weights.columns]
    aligned = weights.reindex([str(spot) for spot in spots])
    arr = np.asarray(aligned.to_numpy(dtype=np.float64), dtype=np.float64)
    if arr.ndim != 2:
        arr = np.atleast_2d(arr)
    # 逐 spot 状态
    finite_row = np.isfinite(arr).all(axis=1)
    missing_mask = ~np.isfinite(arr).any(axis=1)
    row_sums = np.nansum(arr, axis=1)
    negative = finite_row & (np.nanmin(np.where(np.isfinite(arr), arr, np.inf), axis=1) < -1e-9)
    ok_mask = finite_row & ~negative & (row_sums > 0)
    zero_mask = finite_row & ~negative & (row_sums <= 0)
    props = np.full((n, len(cell_types)), np.nan, dtype=np.float32)
    if ok_mask.any():
        props[ok_mask] = (arr[ok_mask] / row_sums[ok_mask, None]).astype(np.float32)

    # spot_class（RCTD 的 results_df 只在 doublet/multi 模式生成；full 模式为空表）
    spot_class = ["NA"] * n
    if isinstance(results_df, pd.DataFrame) and "spot_class" in results_df.columns:
        try:
            series = results_df["spot_class"].reindex([str(spot) for spot in spots])
        except Exception:  # noqa: BLE001
            series = None
        if series is not None:
            for i, value in enumerate(series.tolist()):
                if value is None:
                    continue
                try:
                    if pd.isna(value):
                        continue
                except (TypeError, ValueError):
                    pass
                spot_class[i] = str(value)
    fit_status = np.where(ok_mask, _FIT_OK,
                          np.where(missing_mask, _FIT_MISSING,
                                   np.where(zero_mask, _FIT_ZERO, _FIT_NAN)))
    rejected_mask = np.array([cls.strip().lower() == "reject" for cls in spot_class], dtype=bool)
    fit_status = np.where(rejected_mask & (fit_status != _FIT_OK), _FIT_REJECTED, fit_status)
    fit_status = np.asarray(fit_status, dtype=object)

    n_fitted = int(ok_mask.sum())
    n_missing = int(np.sum(fit_status == _FIT_MISSING))
    n_rejected = int(np.sum(fit_status == _FIT_REJECTED))
    n_zero = int(np.sum(fit_status == _FIT_ZERO))
    n_nan = int(np.sum(fit_status == _FIT_NAN))
    by_class: dict[str, int] = {}
    for cls in spot_class:
        by_class[str(cls)] = int(by_class.get(str(cls), 0) + 1)
    coverage = {
        "n_spots": int(n),
        "n_fitted": n_fitted,
        "fraction_fitted": float(n_fitted / n) if n else 0.0,
        "n_missing": n_missing,
        "n_rejected": n_rejected,
        "n_zero_weights": n_zero,
        "n_nan_weights": n_nan,
        "by_spot_class": by_class,
    }
    warnings: list[str] = []
    if n_missing:
        warnings.append(
            f"RCTD 未返回 {n_missing} 个 spot 的权重（create.RCTD 按 UMI_min/counts_MIN 过滤）："
            "这些 spot 的比例保留 NA，不再用其余 spot 的均值填补")
    if n_rejected:
        warnings.append(f"{n_rejected} 个 spot 被 RCTD 标记为 reject（spot_class=reject）：比例保留 NA")
    if n_zero:
        warnings.append(f"{n_zero} 个 spot 的 RCTD 权重全为 0：比例保留 NA（不用均匀分布填补）")
    if n_nan:
        warnings.append(f"{n_nan} 个 spot 的 RCTD 权重含 NaN/负值：比例保留 NA")
    warnings.append(
        f"反卷积覆盖率：{n_fitted}/{n} 个 spot 有可靠比例（{coverage['fraction_fitted'] * 100:.1f}%）；"
        "下游汇总按可用 spot 计算（nanmean）")
    return {
        "proportions": props,
        "cell_types": cell_types,
        "fit_status": [str(v) for v in fit_status.tolist()],
        "spot_class": [str(v) for v in spot_class],
        "coverage": coverage,
        "warnings": warnings,
        "n_missing_aligned": n_missing,
    }


def _deconvolve_rctd(ctx, params, outdir=None) -> dict:
    """RCTD 反卷积：调 lib/r_rctd.R（spacexr 官方实现），返回与 _deconvolve 同构的结果。"""
    mode = str(_p(params, "mode", "full") or "full").lower()
    if mode not in _RCTD_MODES:
        raise sc.prerequisite(
            f"params.mode={mode} 不合法：只能是 {'/'.join(_RCTD_MODES)}（spacexr 的 doublet_mode）")
    max_cores = _p_int(params, "max_cores", 1)
    max_cores = max(1, min(max_cores, 4))
    min_genes = _p_int(params, "min_genes", 50)
    if not rbridge.rctd_available():
        raise sc.EngineError(
            "unsupported",
            "需要 R 包 spacexr（RCTD 官方实现）：" + rbridge.unavailable_reason()
            + "；按 INSTALL.md 安装（Bioconductor/GitHub），"
              "装好前请用 method=nnls/dwls 内置实现")
    ref = _rctd_reference(ctx, params)
    spot_counts, spot_source, spot_warnings = _rctd_spatial_counts(ctx, params)
    genes_spot = _genes(ctx)
    spots = _spots(ctx)
    ref_index = {gene: i for i, gene in enumerate(ref["genes"])}
    common = [gene for gene in genes_spot if gene in ref_index]
    if len(common) < min_genes:
        raise sc.prerequisite(
            f"参考与数据集的共同基因 {len(common)} 个 < min_genes={min_genes}：请确认两者基因命名一致"
            "（symbol/Ensembl），或降低 min_genes")
    spot_index = {gene: i for i, gene in enumerate(genes_spot)}
    ref_counts = np.asarray(ref["counts"], dtype=np.float64)[:, [ref_index[g] for g in common]]
    spot_sub = np.asarray(spot_counts, dtype=np.float64)[:, [spot_index[g] for g in common]]
    coords = ctx.get("coords")
    warnings: list[str] = list(ref.get("warnings") or []) + list(spot_warnings)
    if coords is None:
        warnings.append("数据集没有坐标：RCTD 用 (0,0) 占位坐标（full 模式的权重不依赖坐标）")
    r_params = {"mode": mode, "max_cores": max_cores, "min_genes": int(min_genes),
                "cell_type_col": str(_p(params, "cell_type_col", "") or "")}
    for key in _RCTD_OPTIONAL_KEYS:
        value = _p(params, key, None)
        if value is not None:
            r_params[key] = value
    payload = {
        "ref_counts": ref_counts,
        "ref_cells": [str(c) for c in ref["cells"]],
        "ref_types": [str(t) for t in ref["types"]],
        "spot_counts": spot_sub,
        "spot_names": [str(s) for s in spots],
        "genes": [str(g) for g in common],
        "coords": coords,
        "params": r_params,
    }
    rdir = os.path.join(_outdir(ctx, outdir), "rctd")
    run = rbridge.run_r("r_rctd.R", payload, rdir, timeout=rbridge.MAX_TIMEOUT)
    result = run.get("result") or {}
    weights = (run.get("tables") or {}).get("weights")
    if weights is None or not isinstance(weights, pd.DataFrame) or weights.empty:
        raise sc.EngineError(
            "runtime",
            "RCTD 跑完但没有返回权重表（tables/weights.csv 缺失或为空）："
            f"spacexr 版本 {run.get('spacexr_version') or '未知'}")
    results_df = (run.get("tables") or {}).get("results_df")
    if results_df is not None and not isinstance(results_df, pd.DataFrame):
        results_df = None
    if results_df is not None and results_df.empty:
        results_df = None
    aligned_out = _rctd_align_weights(weights, spots, results_df)
    arr = aligned_out["proportions"]
    cell_types = aligned_out["cell_types"]
    coverage = aligned_out["coverage"]
    n_missing = aligned_out["n_missing_aligned"]
    warnings.extend(aligned_out["warnings"])
    dropped = result.get("n_spots_dropped")
    if isinstance(dropped, int) and dropped > 0:
        warnings.append(f"create.RCTD 过滤掉 {dropped} 个低质量 spot（已计入覆盖率）")
    out = {
        "proportions": arr.astype(np.float32),
        "cell_types": cell_types,
        "fit_status": aligned_out["fit_status"],
        "spot_class": aligned_out["spot_class"],
        "coverage": coverage,
        "method": _rctd_method_label(mode, max_cores),
        "method_label": "RCTD",
        "n_genes": int(len(common)),
        "n_types": int(len(cell_types)),
        "genes": _genes(ctx),
        "warnings": warnings,
        "rctd": {
            "mode": mode,
            "max_cores": max_cores,
            "spacexr_version": run.get("spacexr_version") or result.get("spacexr_version"),
            "r_version": result.get("r_version"),
            "elapsed_s": result.get("elapsed_s"),
            "ms": run.get("ms"),
            "n_spots": result.get("n_spots"),
            "n_cells": result.get("n_cells"),
            "n_spots_dropped": dropped,
            "raw_rowsum_median": result.get("raw_rowsum_median"),
            "n_missing_aligned": n_missing,
            "coverage": coverage,
            "ref_source": ref.get("source"),
            "spot_counts_source": spot_source,
            "results_df": results_df,
        },
    }
    for item in (run.get("warnings") or []):
        out["warnings"].append(str(item))
    return out


def _deconvolve(ctx, params, outdir=None) -> dict:
    method = str(_p(params, "method", "nnls")).lower()
    if method in ("rctd", "spacexr"):
        return _deconvolve_rctd(ctx, params, outdir)
    X = _X(ctx)
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
    """细胞比例组成图：NaN（RCTD 未拟合/拒绝的 spot）自动排除，按可用 spot 求 nanmean。"""
    props = np.asarray(props, dtype=float)
    if props.ndim != 2 or props.shape[1] != len(types):
        raise sc.prerequisite("反卷积比例矩阵与细胞类型数不一致，无法画组成图")
    usable = np.isfinite(props).all(axis=1)
    if not usable.any():
        raise sc.prerequisite("所有 spot 的反卷积比例都是 NA（RCTD 全部未拟合），无法画组成图")
    labels = _label_array(ctx, str(_p(params, "group_col", "") or ""))
    order = None
    rows = []
    if labels is not None:
        order = list(dict.fromkeys([str(v) for v in labels.tolist()]))
        for value in order:
            sel = np.array([str(v) == value for v in labels]) & usable
            if not sel.any():
                continue
            rows.append({"group": value, "n_spots": int(sel.sum()),
                         **{types[j]: float(np.nanmean(props[sel, j])) for j in range(len(types))}})
        if not rows:
            raise sc.prerequisite("每个分组里都没有可用的反卷积比例（全为 NA），无法画组成图")
        xlabels = [row["group"] for row in rows]
        means = np.array([[row[t] for t in types] for row in rows], dtype=float)
    else:
        means = np.nanmean(props[usable], axis=0, keepdims=True)
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
    out = _deconvolve(ctx, params, outdir)
    props, types = out["proportions"], out["cell_types"]
    outdir = _outdir(ctx, outdir)
    spots = _spots(ctx)
    frame = pd.DataFrame(props, columns=types)
    frame.insert(0, "spot", spots[:len(frame)])
    coverage = out.get("coverage") or {}
    if out.get("fit_status"):
        frame.insert(1, "fit_status", out["fit_status"][:len(frame)])
    if out.get("spot_class"):
        frame.insert(2, "spot_class", out["spot_class"][:len(frame)])
    rows = frame.head(500).to_dict(orient="records")
    figures = [_fig_composition(ctx, params, outdir, props, types)]
    usable = np.isfinite(np.asarray(props, dtype=float)).all(axis=1)
    type_means = {types[j]: float(np.nanmean(np.asarray(props, dtype=float)[usable, j]))
                  for j in range(len(types))} if usable.any() else {t: float("nan") for t in types}
    top_type = max(type_means, key=lambda t: type_means[t]) if type_means else ""
    if top_type:
        fig, ax = sc.new_fig(6.2, 5.0)
        sc.scatter_spatial(ax, ctx, values=props[:, types.index(top_type)], point_size=6.0,
                           legend_title=top_type)
        ax.set_title(f"Estimated proportion: {top_type}", fontsize=10)
        figures.append(sc.save_fig(fig, outdir, "deconv-map", title=None))
    display = out.get("method_label") or str(out["method"]).upper()
    bullets = [
        f"{display} 反卷积：{len(spots)} 个 spot × {len(types)} 种细胞（共同基因 {out['n_genes']} 个）",
        "各类型平均比例：" + "，".join(f"{t}:{type_means[t]:.2f}" for t in list(types)[:6]),
    ]
    if coverage:
        bullets.append(
            f"反卷积覆盖率：{coverage.get('n_fitted')}/{coverage.get('n_spots')} 个 spot 有可靠比例"
            f"（{float(coverage.get('fraction_fitted') or 0) * 100:.1f}%）；"
            "未拟合/被拒绝的 spot 保留 NA（不填补），均值用 nanmean 计算")
    rctd = out.get("rctd") or {}
    if rctd:
        bullets.append(
            f"spacexr::run.RCTD（doublet_mode={rctd.get('mode')}，max_cores={rctd.get('max_cores')}）："
            f"{rctd.get('n_spots')} 个 spot 参与估计，耗时 {float(rctd.get('elapsed_s') or 0):.1f}s，"
            f"原始权重行和中位数 {float(rctd.get('raw_rowsum_median') or 0):.3f}"
            "（已按 spacexr::normalize_weights 归一化）")
        spot_class_values = out.get("spot_class") or []
        nontrivial = [v for v in spot_class_values if v not in ("NA", "", "None")]
        if nontrivial:
            series = pd.Series([str(v) for v in spot_class_values]).value_counts()
            bullets.append("spot 分类（results_df）：" + "，".join(
                f"{k}:{int(v)}" for k, v in series.items()))
        fit_values = out.get("fit_status") or []
        if fit_values:
            fit_counts = pd.Series([str(v) for v in fit_values]).value_counts()
            bullets.append("拟合状态：" + "，".join(f"{k}:{int(v)}" for k, v in fit_counts.items()))
    notes = ["方法说明：NNLS/DWLS 是内置实现；RCTD 走 R 包 spacexr 的官方实现"
             "（spacexr::create.RCTD/run.RCTD，Cable et al. 2022 Nat Biotechnol）"]
    if str(_p(params, "method", "nnls")).lower() not in ("rctd", "spacexr"):
        notes += _rctd_notes("step")
    table_columns = [{"key": "spot", "label": "spot", "type": "string"}]
    if out.get("fit_status"):
        table_columns.append({"key": "fit_status", "label": "拟合状态", "type": "string"})
    if out.get("spot_class"):
        table_columns.append({"key": "spot_class", "label": "RCTD spot_class", "type": "string"})
    table_columns += [{"key": t, "label": t, "type": "number"} for t in types]
    payload = sc.result(
        summary={"title": "反卷积", "bullets": bullets,
                 "metrics": {"n_types": len(types), "n_genes": out["n_genes"], "method": out["method"],
                             "coverage": coverage}},
        tables=[sc.table("deconv", "每个 spot 的细胞类型比例", table_columns, rows)],
        figures=figures,
        warnings=out.get("warnings") or [],
        notes=notes,
        data_out=[{"kind": "deconv_csv", "path": sc.export_csv(frame, outdir, "deconv-proportions.csv", index=False),
                   "label": "反卷积比例表"}],
        store={"deconv": {"cell_types": types, "proportions": props, "method": out["method"],
                          "fit_status": out.get("fit_status"), "spot_class": out.get("spot_class"),
                          "coverage": coverage}},
    )
    # summary.method 写真实实现（RCTD 时含 "RCTD（spacexr，doublet_mode=…）"）
    payload["summary"]["method"] = out["method"]
    return payload


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
               "平均比例：" + "，".join(f"{t}:{float(np.nanmean(props[:, j])):.2f}" for j, t in enumerate(types[:6]))]
    return sc.result(summary={"title": "反卷积（NNLS）", "bullets": bullets,
                              "metrics": {"n_types": len(types), "n_genes": out["n_genes"]}},
                     tables=[sc.table("deconv_nnls", "平均比例",
                                      columns=[{"key": "cell_type", "label": "细胞类型", "type": "string"},
                                               {"key": "mean", "label": "平均比例", "type": "number"},
                                               {"key": "median", "label": "中位比例", "type": "number"}],
                                      rows=[{"cell_type": t, "mean": float(np.nanmean(props[:, j])),
                                             "median": float(np.nanmedian(props[:, j]))}
                                            for j, t in enumerate(types)])],
                     figures=figures,
                     notes=["方法说明：本分析是内置非负最小二乘实现（不经 R）"]
                           + _rctd_notes("analysis"),
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
                                      rows=[{"cell_type": t, "mean": float(np.nanmean(props[:, j]))}
                                            for j, t in enumerate(types)])],
                     figures=figures,
                     store={"deconv_dwls": {"cell_types": types, "proportions": props}})


def analysis_deconv_rctd(ctx, params, outdir):
    """RCTD（spacexr 官方实现）反卷积。

    说明：本分析放在「反卷积」组、needs=["normalize"]。没有单细胞参考时**不抛错**，而是返回
    一个明确的「未执行」说明信封 —— 分析目录会被批量遍历（回归测试用 demo 跑全部分析），
    外部输入缺失不该让整批分析失败；但这里绝不产出数字、也绝不静默回退成 NNLS/DWLS。
    缺 R 包 spacexr 时仍然抛 unsupported（与契约一致，提示按 INSTALL.md 安装）。
    """
    outdir = _outdir(ctx, outdir)
    ref_available = bool(ctx.get("ref")) or bool(str(_p(params, "ref", "") or ""))
    if not ref_available:
        return sc.result(
            summary={"title": "反卷积（RCTD）· 缺少单细胞参考（未执行）",
                     "bullets": [
                         "未执行：RCTD 需要单细胞参考 h5ad（先 ref_attach 挂载，或在参数 ref 里填路径）",
                         "参考要求：原始 counts（整数矩阵，h5ad 的 layers[\"counts\"] 或 raw/整数 X）"
                         "+ 细胞类型列（cell_type/celltype/annotation/clusters，或用 cell_type_col 指定）；"
                         "细胞类型 ≥2 种、每类 ≥25 个细胞（spacexr 的 CELL_MIN_INSTANCE）",
                         "本分析不会回退成 NNLS/DWLS：需要内置实现请跑「反卷积（NNLS）」或「反卷积（加权 NNLS）」",
                     ],
                     "metrics": {"executed": 0}},
            warnings=["RCTD 未执行（缺单细胞参考）：没有产出任何比例，也没有换用其它方法"],
            notes=["分析目录会被批量遍历：缺外部输入时返回说明信封而不是抛错，"
                   "避免一次缺参考导致整批分析失败"])
    out = _deconvolve(ctx, {**(params or {}), "method": "rctd"}, outdir)
    props, types = out["proportions"], out["cell_types"]
    rctd = out.get("rctd") or {}
    coverage = out.get("coverage") or rctd.get("coverage") or {}
    figures = [_fig_composition(ctx, params, outdir, props, types, name="deconv-rctd-bar",
                                title="RCTD deconvolution composition")]
    rows = [{"cell_type": t, "mean": float(np.nanmean(props[:, j])),
             "median": float(np.nanmedian(props[:, j]))} for j, t in enumerate(types)]
    bullets = [
        f"{out['method']} 反卷积：{len(types)} 种细胞，共同基因 {out['n_genes']} 个",
        f"spacexr::create.RCTD + run.RCTD（doublet_mode={rctd.get('mode')}，"
        f"max_cores={rctd.get('max_cores')}，spacexr {rctd.get('spacexr_version') or '?'}）："
        f"{rctd.get('n_spots')} 个 spot、耗时 {float(rctd.get('elapsed_s') or 0):.1f}s",
        "平均比例（仅可用 spot，nanmean）：" + "，".join(
            f"{t}:{float(np.nanmean(props[:, j])):.2f}" for j, t in enumerate(types[:6])),
    ]
    if coverage:
        bullets.append(
            f"反卷积覆盖率：{coverage.get('n_fitted')}/{coverage.get('n_spots')} 个 spot 有可靠比例"
            f"（{float(coverage.get('fraction_fitted') or 0) * 100:.1f}%）；"
            f"未拟合 {coverage.get('n_missing', 0)}、reject {coverage.get('n_rejected', 0)}、"
            f"全 0 权重 {coverage.get('n_zero_weights', 0)} 个 spot 的比例保留 NA，不做任何填补")
    spot_class_values = out.get("spot_class") or []
    nontrivial = [v for v in spot_class_values if v not in ("NA", "", "None")]
    if nontrivial:
        series = pd.Series([str(v) for v in spot_class_values]).value_counts()
        bullets.append("spot 分类（RCTD results_df）：" + "，".join(
            f"{k}:{int(v)}" for k, v in series.items()))
    fit_values = out.get("fit_status") or []
    if fit_values:
        fit_counts = pd.Series([str(v) for v in fit_values]).value_counts()
        bullets.append("拟合状态：" + "，".join(f"{k}:{int(v)}" for k, v in fit_counts.items()))
    weight_frame = pd.DataFrame(props, index=pd.Index(_spots(ctx)[:len(props)], name="spot"),
                                columns=types)
    weight_frame.insert(0, "fit_status", (out.get("fit_status") or ["NA"] * len(weight_frame))[:len(weight_frame)])
    weight_frame.insert(1, "spot_class", (out.get("spot_class") or ["NA"] * len(weight_frame))[:len(weight_frame)])
    data_out = [{"kind": "deconv_rctd_weights", "path": sc.export_csv(
        weight_frame, outdir, "deconv-rctd-weights.csv", index=True),
        "label": "RCTD 每个 spot 的细胞类型比例 + fit_status + spot_class（未拟合为 NA）"}]
    if results_df is not None and not results_df.empty:
        data_out.append({"kind": "deconv_rctd_results_df", "path": sc.export_csv(
            results_df, outdir, "deconv-rctd-results_df.csv", index=True),
            "label": "RCTD 的 spot 分类表（results_df）"})
    payload = sc.result(
        summary={"title": "反卷积（RCTD）", "bullets": bullets,
                 "metrics": {"n_types": len(types), "n_genes": out["n_genes"],
                             "method": out["method"], "mode": rctd.get("mode"),
                             "elapsed_s": rctd.get("elapsed_s"),
                             "n_spots_dropped": rctd.get("n_spots_dropped"),
                             "coverage": coverage}},
        tables=[sc.table("deconv_rctd", "平均比例（RCTD 归一化权重）",
                         columns=[{"key": "cell_type", "label": "细胞类型", "type": "string"},
                                  {"key": "mean", "label": "平均比例", "type": "number"},
                                  {"key": "median", "label": "中位比例", "type": "number"}],
                         rows=rows)],
        figures=figures,
        warnings=out.get("warnings") or [],
        notes=["方法说明：RCTD 由 R 包 spacexr 官方实现（Cable et al. 2022 Nat Biotechnol）执行；"
               "只做 doublet_mode 切换与权重归一化（spacexr::normalize_weights），不重算、不近似",
               "回退条件：缺 Rscript/spacexr → unsupported（按 INSTALL.md 安装）；参考不是原始 counts、"
               "共同基因不足、运行报错 → 给真实原因，不自动改换方法；需要内置实现请跑 deconv_nnls/deconv_dwls"],
        data_out=data_out,
        store={"deconv_rctd": {"cell_types": types, "proportions": props,
                               "method": out["method"], "mode": rctd.get("mode"),
                               "fit_status": out.get("fit_status"),
                               "spot_class": out.get("spot_class"),
                               "coverage": coverage}})
    # summary.method 写真实实现（含 "RCTD（spacexr，doublet_mode=…）"）
    payload["summary"]["method"] = out["method"]
    return payload


def _cached_props(ctx, params):
    """取比例列（RCTD 未拟合 spot 为 NaN，原样保留；画图时散点自动跳过）。"""
    store = ctx.get("results") or {}
    cell_type = str(_p(params, "cell_type", "") or "")
    for key in ("deconv", "deconv_nnls", "deconv_dwls", "deconv_rctd"):
        entry = store.get(key)
        if isinstance(entry, dict) and "proportions" in entry:
            names = list(entry.get("cell_types") or [])
            props = np.asarray(entry["proportions"], dtype=np.float32)
            if props.ndim != 2 or props.shape[1] != len(names) or not names:
                continue
            if cell_type and cell_type in names:
                return props[:, names.index(cell_type)], cell_type, names
            if not cell_type:
                return props[:, 0], names[0], names
    out = _deconvolve(ctx, params)
    names = out["cell_types"]
    props = np.asarray(out["proportions"], dtype=np.float32)
    if cell_type and cell_type in names:
        return props[:, names.index(cell_type)], cell_type, names
    return props[:, 0], names[0], names


def analysis_deconv_spatial_map(ctx, params, outdir):
    values, cell_type, names = _cached_props(ctx, params)
    outdir = _outdir(ctx, outdir)
    fig, ax = sc.new_fig(6.2, 5.0)
    sc.scatter_spatial(ax, ctx, values=values, point_size=_p_float(params, "point_size", 8.0),
                       legend_title=cell_type)
    ax.set_title(f"Estimated proportion: {cell_type}", fontsize=10)
    figure = sc.save_fig(fig, outdir, "deconv-spatial", title=None)
    finite = np.isfinite(values)
    bullets = [f"细胞类型「{cell_type}」在组织切片上的比例分布（均值 {float(np.nanmean(values)):.3f}；"
               f"可用 spot {int(finite.sum())}/{values.size}）",
               f"可用类型：{'、'.join(names[:8])}"]
    if not finite.all():
        bullets.append(f"{int((~finite).sum())} 个 spot 的比例为 NA（RCTD 未拟合/被拒绝），图中不显示")
    return sc.result(summary={"title": "细胞类型空间分布", "bullets": bullets,
                              "metrics": {"cell_type": cell_type,
                                          "mean": float(np.nanmean(values)),
                                          "n_finite": int(np.isfinite(values).sum()),
                                          "n_spots": int(values.size)}},
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
        mask = np.isfinite(va) & np.isfinite(vb)
        if mask.sum() < 3 or np.std(va[mask]) == 0 or np.std(vb[mask]) == 0:
            r = float("nan")
        else:
            r = float(np.corrcoef(va[mask], vb[mask])[0, 1])
        corr[t] = r
        rows.append({"cell_type": t, "pearson": r, "nnls_mean": float(np.nanmean(va)),
                     "dwls_mean": float(np.nanmean(vb))})
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


_REGION_COLUMNS = ("domains", "domain", "spatial_domain", "domain_id", "region",
                   "clusters", "cluster")
_GROUP_COLUMNS = ("group", "condition", "treatment", "group_id", "分组")


def _deconv_group_labels(ctx, n: int, group_col: str = ""):
    """spot 级分组标签（group_col 显式指定 > ctx.group_key > 常见列名）。"""
    frame = _info(ctx)
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        return None, ""
    candidates = [group_col] if group_col else []
    ctx_group = str(ctx.get("group_key") or "")
    if ctx_group:
        candidates.append(ctx_group)
    candidates.extend(_GROUP_COLUMNS)
    for cand in candidates:
        if not cand or cand not in frame.columns:
            continue
        try:
            series = frame[cand]
            values = []
            for value in series.tolist():
                if value is None:
                    values.append("NA")
                    continue
                try:
                    if pd.isna(value):
                        values.append("NA")
                        continue
                except (TypeError, ValueError):
                    pass
                values.append(str(value).strip() or "NA")
        except Exception:  # noqa: BLE001
            continue
        if len(values) == n and any(v not in ("NA", "") for v in values):
            return np.asarray(values, dtype=object), str(cand)
    return None, ""


def _deconv_region_labels(ctx, n: int):
    """区域标签（域优先，其次聚类）：返回 (labels, key)。"""
    frame = _info(ctx)
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        return np.asarray(["all"] * n, dtype=object), "all"
    for cand in _REGION_COLUMNS:
        if cand not in frame.columns:
            continue
        try:
            values = [str(v) for v in frame[cand].tolist()]
        except Exception:  # noqa: BLE001
            continue
        if len(values) == n and len(set(values)) >= 1:
            return np.asarray(values, dtype=object), str(cand)
    return np.asarray(["all"] * n, dtype=object), "all"


def _deconv_units(ctx, props: np.ndarray, group_col: str = ""):
    """把 spot 级比例聚合成 供体/切片 × 区域 单位（NA spot 不参与均值）。

    返回 dict：units=[{block, region, group, n_spots, n_used, values[]}]、分组顺序、
    block_key、region_key、coverage、warnings。若每组独立供体/切片 < 2 → inference=descriptive。
    """
    n = int(props.shape[0])
    props = np.asarray(props, dtype=float)
    usable_rows = np.isfinite(props).all(axis=1)
    blocks, block_names, block_key, block_warns = sc.spatial_blocks(ctx)
    block_id = np.zeros(n, dtype=int)
    for i, idx in enumerate(blocks):
        block_id[np.asarray(idx, dtype=int)] = i
    groups, group_key = _deconv_group_labels(ctx, n, group_col)
    if groups is None:
        raise sc.prerequisite("找不到分组/条件列（group/condition）：请先设置分组列或指定 group_col")
    regions, region_key = _deconv_region_labels(ctx, n)
    warnings: list[str] = list(block_warns or [])
    grouped: dict[tuple[int, str], list[int]] = {}
    for i in range(n):
        grouped.setdefault((int(block_id[i]), str(regions[i])), []).append(i)
    units: list[dict] = []
    dropped = 0
    for (bi, region), idx in sorted(grouped.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
        idx_arr = np.asarray(idx, dtype=int)
        used = idx_arr[usable_rows[idx_arr]]
        if used.size < 3:
            dropped += 1
            continue
        unit_groups = list(dict.fromkeys([str(groups[i]) for i in idx_arr]))
        if len(unit_groups) > 1:
            warnings.append(
                f"单位 {block_names[bi]}/{region} 内混有多个分组（{unit_groups[:3]}）：按首个 spot 归类")
        units.append({"block": str(block_names[bi]), "region": str(region),
                      "group": unit_groups[0] if unit_groups else "NA",
                      "n_spots": int(idx_arr.size), "n_used": int(used.size),
                      "values": np.nanmean(props[used], axis=0)})
    if dropped:
        warnings.append(f"{dropped} 个单位可用 spot < 3（或全为 NA），已从单位级比较中剔除")
    if not units:
        raise sc.prerequisite("没有可用 spot 数 ≥ 3 的 供体/切片×区域 单位，无法做区域组成比较")
    group_order = list(dict.fromkeys([u["group"] for u in units]))
    group_order = [g for g in group_order if g not in ("NA", "")] or group_order
    n_used_total = int(usable_rows.sum())
    coverage = {
        "n_spots": n,
        "n_usable": n_used_total,
        "n_spots_used": n_used_total,
        "n_na": int(n - n_used_total),
        "n_spots_na": int(n - n_used_total),
        "fraction_usable": float(n_used_total / n) if n else 0.0,
    }
    # 推断前提：每组独立供体/切片 ≥ 2；同一供体的多个区域不是独立重复
    block_counts = {g: len({u["block"] for u in units if u["group"] == g}) for g in group_order}
    inference = "unit_level" if (len(group_order) >= 2 and all(v >= 2 for v in block_counts.values())) else "descriptive"
    return {"units": units, "group_order": group_order, "block_key": block_key,
            "region_key": region_key, "group_key": group_key, "block_counts": block_counts,
            "coverage": coverage, "inference": inference, "warnings": warnings}


def analysis_deconv_roi(ctx, params, outdir):
    """区域细胞组成差异：按 供体/切片 × 区域 单位检验（不再把 spot 当独立样本）。

    为什么不再按 spot 做 Kruskal-Wallis：同一张切片的 spot 不是独立生物学重复，
    spot 级检验的 p 值会被样本量放大到 1e-18 量级（空间伪重复）。正确做法：
    每个 供体/切片×区域 单位取该单位内 spot 比例的均值作为一个观测，单位之间再检验。
    NA（RCTD 未拟合/被拒绝的 spot）不参与均值，并报告覆盖率。每组独立供体/切片 < 2 时
    只输出描述性（descriptive）汇总，不做推断（inferential）。
    """
    params = dict(params or {})
    store = ctx.get("results") or {}
    entry = None
    for key in ("deconv", "deconv_nnls", "deconv_dwls", "deconv_rctd"):
        candidate = store.get(key)
        if isinstance(candidate, dict) and "proportions" in candidate:
            entry = candidate
            break
    if not isinstance(entry, dict) or "proportions" not in entry:
        raise sc.prerequisite("请先跑「反卷积」步骤（或 deconv_nnls 分析），再做区域组成差异")
    props = np.asarray(entry["proportions"], dtype=np.float32)
    types = list(entry.get("cell_types") or [])
    if props.ndim != 2 or props.shape[1] != len(types):
        raise sc.prerequisite("反卷积比例矩阵与细胞类型数不一致，无法比较区域组成")
    group_col = str(_p(params, "group_col", "") or "")
    units_info = _deconv_units(ctx, props, group_col)
    units = units_info["units"]
    group_order = units_info["group_order"]
    if len(group_order) < 2:
        raise sc.prerequisite(f"分组只有 {len(group_order)} 个有效取值，无法比较")
    inference = units_info["inference"]
    warnings: list[str] = list(units_info["warnings"])
    coverage = units_info["coverage"]
    region_key = units_info["region_key"]
    block_key = units_info["block_key"]
    from scipy import stats

    rows: list[dict] = []
    for j, t in enumerate(types):
        groups = [[u["values"][j] for u in units if u["group"] == g] for g in group_order]
        row: dict = {"cell_type": t}
        if inference == "unit_level" and all(len(g) >= 2 for g in groups):
            try:
                stat, p = stats.kruskal(*groups)
            except Exception:  # noqa: BLE001
                stat, p = float("nan"), float("nan")
        else:
            stat, p = float("nan"), float("nan")
        row.update({"stat": float(stat), "p": float(p)})
        for g, values in zip(group_order, groups):
            row[f"n_units_{g}"] = int(len(values))
            row[f"mean_{g}"] = float(np.mean(values)) if values else float("nan")
        rows.append(row)
    if inference == "unit_level":
        try:
            adjusted = sc.bh_fdr([row["p"] for row in rows])
            for row, q in zip(rows, adjusted):
                row["p_adj"] = float(q)
        except Exception:  # noqa: BLE001
            for row in rows:
                row["p_adj"] = float("nan")
    else:
        for row in rows:
            row["p_adj"] = float("nan")
    outdir = _outdir(ctx, outdir)
    figure = _fig_composition(ctx, params, outdir, props, types, name="deconv-roi-bar",
                              title="Cell type composition by group")
    n_units = len(units)
    unit_desc = f"{block_key}×{region_key}"
    group_text = "、".join(
        f"{g}={int(np.sum([u['group'] == g for u in units]))}" for g in group_order)
    if inference == "unit_level":
        bullets = [
            f"单位级（{unit_desc}）：{n_units} 个单位，分组 {group_text}",
            "Kruskal-Wallis + BH 校正；p_adj<0.05 的类型：" +
            ("、".join(row["cell_type"] for row in rows if row.get("p_adj", 1) < 0.05) or "无"),
        ]
        notes = [
            f"方法：按 供体/切片×区域 单位取比例均值，再对单位做 Kruskal-Wallis（scipy）+ BH；"
            f"单位={unit_desc}",
            "不按 spot 做检验：同一张切片的 spot 不是独立重复（空间伪重复）",
        ]
    else:
        block_text = "、".join(
            f"{g} 独立供体/切片={units_info['block_counts'].get(g, 0)} 个" for g in group_order)
        bullets = [
            f"仅描述性汇总：{n_units} 个单位（{unit_desc}），{block_text}",
            "每组独立供体/切片 < 2，不做推断（inferential），不报告 p 值",
        ]
        notes = [
            f"方法：描述性（descriptive）汇总——单位={unit_desc} 的比例均值；每组独立供体/切片不足 2 个，"
            "空间伪重复风险，禁止做推断",
        ]
    if coverage["n_na"]:
        bullets.append(
            f"覆盖率：{coverage['n_usable']}/{coverage['n_spots']} 个 spot 有可用比例"
            f"（{coverage['fraction_usable'] * 100:.1f}%）；{coverage['n_na']} 个 NA spot 不参与单位均值")
        warnings.append("反卷积比例含 NA（RCTD 未拟合/被拒绝的 spot）：单位均值只用可用 spot，并报告覆盖率")
    columns = [{"key": "cell_type", "label": "细胞类型", "type": "string"},
               {"key": "stat", "label": "统计量", "type": "number"},
               {"key": "p", "label": "p", "type": "number"},
               {"key": "p_adj", "label": "p_adj", "type": "number"}]
    for g in group_order:
        columns.append({"key": f"mean_{g}", "label": f"{g} 均值", "type": "number"})
        columns.append({"key": f"n_units_{g}", "label": f"{g} 单位数", "type": "number"})
    metrics = {
        "n_types": len(types),
        "n_groups": len(group_order),
        "inference": inference,
        "n_units": int(n_units),
        "unit": unit_desc,
        "coverage": coverage,
    }
    if len(group_order) == 2:
        metrics["n_case_units"] = int(np.sum([u["group"] == group_order[0] for u in units]))
        metrics["n_ctrl_units"] = int(np.sum([u["group"] == group_order[1] for u in units]))
        metrics["n_case_blocks"] = int(units_info["block_counts"].get(group_order[0], 0))
        metrics["n_ctrl_blocks"] = int(units_info["block_counts"].get(group_order[1], 0))
    return sc.result(
        title="区域细胞组成差异（单位级）" if inference == "unit_level" else "区域细胞组成差异（描述性）",
        bullets=bullets,
        metrics=metrics,
        tables=[sc.table("deconv_roi", "各分组细胞比例（单位级）", columns, rows)],
        figures=[figure],
        warnings=warnings,
        notes=notes,
        store={"deconv_roi": {"inference": inference, "n_units": int(n_units),
                              "coverage": coverage}},
    )


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
    ratio_all = props[:, j]
    usable = np.isfinite(ratio_all)
    if usable.sum() < 3:
        raise sc.prerequisite(
            f"「{cell_type}」只有 {int(usable.sum())} 个 spot 有可用比例（其余为 NA），无法做相关分析")
    ratio = ratio_all[usable]
    X = _X(ctx)
    if X.shape[0] == ratio_all.shape[0]:
        X = np.asarray(X, dtype=np.float32)[usable]
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
    if int((~usable).sum()):
        bullets.append(f"{int((~usable).sum())} 个 spot 的比例为 NA（RCTD 未拟合/被拒绝），相关分析只用可用 spot")
    return sc.result(summary={"title": "比例-标记基因一致性", "bullets": bullets,
                              "metrics": {"cell_type": cell_type, "n_genes": len(rows),
                                          "n_spot_used": int(usable.sum()),
                                          "n_spot_na": int((~usable).sum())}},
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
             "options": [{"value": "nnls", "label": "NNLS 非负最小二乘（内置）"},
                         {"value": "dwls", "label": "DWLS 阻尼加权（内置）"},
                         {"value": "rctd", "label": "RCTD（spacexr 官方实现，需 R 包）"}]},
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
        "id": "deconv_rctd", "title": "反卷积（RCTD）", "group": "反卷积",
        "method": "RCTD（spacexr 官方实现：create.RCTD/run.RCTD；Cable et al. 2022 Nat Biotechnol；需 R）",
        "needs": ["normalize"],
        "params": [{"key": "ref", "label": "参考 h5ad 路径（可留空用已挂载）", "type": "text", "default": ""},
                   {"key": "cell_type_col", "label": "细胞类型列", "type": "text", "default": ""},
                   {"key": "mode", "label": "doublet_mode", "type": "select", "default": "full",
                    "options": [{"value": "full", "label": "full（全混合，输出比例）"},
                                {"value": "doublet", "label": "doublet（单细胞/双细胞分类）"},
                                {"value": "multi", "label": "multi（多细胞类型组合）"}]},
                   {"key": "min_genes", "label": "共同基因下限", "type": "number", "default": 50}],
        "run": analysis_deconv_rctd,
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
