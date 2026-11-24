#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-spatial — 空间统计模块（``group="空间统计"``）。

本文件提供 9 个空间统计分析与 5 个可复用绘图入口，**全部为内置实现**
（只依赖 numpy / scipy / pandas / matplotlib），不依赖 squidpy、SPARK、SPARK-X
（本机未安装，见 AGENT_BRIEF §1）。

方法与文献
----------
=====================  ==========================================================
分析 id                方法与文献
=====================  ==========================================================
``morans_svg``         单变量 Moran's I（Moran 1950）+ 正态近似与标签置换
``gearys_svg``         Geary's C（Geary 1954）+ 标签置换经验 p
``svg_permutation``    纯标签置换经验 p（SPARK-X 思想，Sun et al. 2021）
``cross_moran``        双变量 Moran's I（Wartenberg 1985）
``hotspot``            Getis-Ord Gi*（Getis & Ord 1992），kNN 邻域
``spatial_gradient``   沿轴投影等频分箱 + 线性趋势检验
``spatial_corr_pairs`` Pearson 相关 + 邻域乘积一致性（空间共表达）
``ripley_k``           Ripley's K / L 函数（Ripley 1977）+ CSR 随机点包络
``gene_autocorr_scan`` 全基因向量化 Moran's I 扫描（稀疏矩阵乘法）
=====================  ==========================================================

空间权重的约定
--------------
* 统一使用 :func:`sp_common.spatial_weights`（半径邻接 + 行标准化）；Moran's I /
  Geary's C 的行标准化版本中 ``S0`` = 权重总和 = 非孤点个数，公式取
  ``I = (n / S0) * (z' W z) / (z' z)``（Cliff & Ord 1981 的行标准化形式）。
* 半径来自 ``ctx["neighbor_radius"]``（import 时估计），可用 ``neighbor_radius``
  参数覆盖；缺失时用 ``sp_common.estimate_radius`` 估计。
* 行和为 0 的孤点：只参与方差、不参与协方差（会在 ``warnings`` 里说明）。

输出约定（DESIGN.md §2.5）
--------------------------
* 结果一律走 ``sp_common.result()``；每张表都完整写入 ``data_out`` CSV，
  JSON 里 ``rows`` 由 ``sp_common.table()`` 截断到 500 行；
* 图存 ``outdir``、白底、dpi=140、**图上文字全部英文**；
* 缺前置 → :func:`sp_common.prerequisite`（中文）；参数非法 → ``bad_request``；
* 随机过程一律 ``np.random.default_rng(0)``（同数据同参数同结果）。
"""

from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd

from sp_common import (
    EngineError,
    bad_request,
    bh_fdr,
    data_file,
    ensure_dir,
    estimate_radius,
    export_csv,
    jsonable,
    new_fig,
    palette,
    prerequisite,
    result,
    save_fig,
    scatter_spatial,
    spatial_weights,
    table,
)

try:  # 仅用于大数据的等价加速（kNN / 半径图），不可用时退回 sp_common 的 O(n^2) 实现
    from scipy.spatial import cKDTree  # noqa: F401
    from scipy import sparse as _sparse
except Exception:  # pragma: no cover - scipy 一定存在，纯防御
    cKDTree = None
    _sparse = None

ANALYSIS_GROUP = "空间统计"

# 默认随机种子：所有置换 / 模拟都用它，保证可复现
DEFAULT_SEED = 0

# 单次分析的置换预算（矩阵乘法规模），避免大样本上长时间卡住
PERM_BUDGET = 2.5e9

# 统计量公式里使用的最小方差，避免除零
EPS = 1e-12

# 维度上限：用 O(n^2) 工具（coord_knn / radius_graph / estimate_radius）前的阈值
DENSE_LIMIT = 12000

# 高表达/低表达配色（热点图、相关网络、Ripley 曲线等固定语义）
HOT_COLOR = "#c0392b"
COLD_COLOR = "#2e6da4"
NS_COLOR = "#c9ced6"


# ===========================================================================
# 参数归一化（UI 可能把数字/布尔/列表传成字符串）
# ===========================================================================

def _p_str(value, default: str = "") -> str:
    if value is None:
        return default
    text = str(value).strip()
    return text if text else default


def _p_float(value, default: float) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    if not math.isfinite(number):
        return float(default)
    return number


def _p_int(value, default: int) -> int:
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return int(default)


def _p_bool(value, default: bool = False) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(value)
    text = str(value).strip().lower()
    if text in ("1", "true", "yes", "y", "on", "是"):
        return True
    if text in ("0", "false", "no", "n", "off", "否", ""):
        return False
    return bool(default)


def _p_list(value) -> list[str]:
    """把 ``"A,B C"`` / ``["A","B"]`` 统一成字符串列表。"""
    if value is None:
        return []
    if isinstance(value, (list, tuple, set, np.ndarray, pd.Index, pd.Series)):
        return [str(v).strip() for v in value if str(v).strip()]
    text = str(value).strip()
    if not text:
        return []
    out = []
    for chunk in text.replace(";", ",").replace("|", ",").replace("\t", ",").split(","):
        for item in chunk.split():
            item = item.strip()
            if item:
                out.append(item)
    return out


def _clip_p(pvals) -> np.ndarray:
    """p 值 clip 到 [0, 1]（NaN 保留为 NaN）。"""
    arr = np.asarray(pvals, dtype=float)
    return np.clip(arr, 0.0, 1.0)


# ===========================================================================
# ctx 校验 / 视图构造
# ===========================================================================

def _require_ctx(ctx) -> dict:
    if not isinstance(ctx, dict) or not ctx:
        raise prerequisite("缺少分析上下文（ctx），无法进行空间统计")
    return ctx


def _view(ctx: dict, *, need_x: bool = True) -> dict:
    """校验 ctx 并构造只读视图：``coords / X / genes / spot_idx / n / warnings``。

    * 坐标缺失 / 维度错误 / spot 过少 → ``prerequisite``（中文原因）；
    * 坐标有非有限值 → 剔除对应 spot 并写 warning（X 同步剔除）；
    * 表达矩阵含非有限值 → 置 0 并写 warning（不静默吞掉）。
    """
    ctx = _require_ctx(ctx)
    warnings: list[str] = []

    coords_raw = ctx.get("coords")
    if coords_raw is None:
        raise prerequisite("该数据集没有空间坐标（coords），无法进行空间统计")
    coords_all = np.asarray(coords_raw, dtype=float)
    if coords_all.ndim != 2 or coords_all.shape[1] < 2:
        raise prerequisite("空间坐标维度不正确（应为 spot × 2 的矩阵）")
    coords_all = coords_all[:, :2]
    n_raw = coords_all.shape[0]
    if n_raw < 3:
        raise prerequisite(f"空间坐标只有 {n_raw} 个 spot，无法进行空间统计")

    finite_coord = np.isfinite(coords_all[:, 0]) & np.isfinite(coords_all[:, 1])
    spot_idx = np.where(finite_coord)[0]
    if not finite_coord.all():
        warnings.append(f"{int((~finite_coord).sum())} 个 spot 的坐标含缺失值，已从空间统计中排除")
    coords = coords_all[finite_coord]
    n = int(coords.shape[0])
    if n < 5:
        raise prerequisite(f"有效空间坐标只有 {n} 个 spot（至少需要 5 个），无法进行空间统计")

    X = None
    genes: list[str] = []
    if need_x:
        X_raw = ctx.get("X")
        if X_raw is None:
            raise prerequisite("缺少表达矩阵 X，无法计算空间统计")
        X_all = np.asarray(X_raw)
        if X_all.ndim != 2:
            raise prerequisite("表达矩阵维度不正确（应为 spot × 基因）")
        if X_all.shape[0] != n_raw:
            raise prerequisite(f"表达矩阵行数（{X_all.shape[0]}）与坐标数（{n_raw}）不一致")
        if X_all.shape[1] == 0:
            raise prerequisite("表达矩阵没有任何基因列")
        X = X_all if finite_coord.all() else X_all[finite_coord]
        if not np.isfinite(X).all():
            bad = int((~np.isfinite(np.asarray(X, dtype=float))).sum())
            warnings.append(f"表达矩阵含 {bad} 个非有限值（NaN/Inf），已按 0 处理")
            X = np.nan_to_num(np.asarray(X, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
        genes = [str(g) for g in (ctx.get("genes") or [])]
        if len(genes) != X.shape[1]:
            genes = [f"GENE{i + 1}" for i in range(X.shape[1])]
            warnings.append("ctx.genes 与表达矩阵列数不一致，已使用占位基因名")

    x_kind = str(ctx.get("X_kind") or "unknown")
    if need_x and x_kind == "counts":
        warnings.append("当前表达矩阵是原始 counts，空间自相关建议先做 normalize（结果可能被测序深度驱动）")

    return {
        "ctx": ctx,
        "coords": coords,
        "coords_all": coords_all,
        "X": X,
        "genes": genes,
        "spot_idx": spot_idx,
        "n": n,
        "n_raw": n_raw,
        "warnings": warnings,
        "X_kind": x_kind,
        "neighbor_radius": _p_float(ctx.get("neighbor_radius"), 0.0),
        "image": ctx.get("image") or {},
        "spots": [str(s) for s in (ctx.get("spots") or [])],
    }


def _spot_names(view: dict) -> np.ndarray:
    """与视图 spot 顺序对齐的 spot 名称（缺失时用序号）。"""
    spots = view["spots"]
    idx = view["spot_idx"]
    if spots and max(idx) < len(spots):
        return np.asarray([spots[i] for i in idx], dtype=object)
    return np.asarray([f"spot{i}" for i in idx], dtype=object)


def _gene_lookup(genes: list[str]):
    exact: dict[str, int] = {}
    lowered: dict[str, int] = {}
    for i, gene in enumerate(genes):
        exact.setdefault(gene, i)
        lowered.setdefault(gene.lower(), i)
    return exact, lowered


def _find_gene(view: dict, name: str) -> int:
    target = _p_str(name)
    exact, lowered = _gene_lookup(view["genes"])
    if target in exact:
        return int(exact[target])
    hit = lowered.get(target.lower())
    if hit is not None:
        return int(hit)
    return -1


def _gene_hint(view: dict, limit: int = 6) -> list[str]:
    """基因缺失时给几个候选（按方差取前几个），帮助用户纠正拼写。"""
    try:
        X = np.asarray(view["X"])
        if X.shape[1] <= limit:
            return list(view["genes"])
        var = np.var(np.asarray(X, dtype=np.float64), axis=0)
        order = np.argsort(var)[::-1][:limit]
        return [view["genes"][i] for i in order]
    except Exception:  # noqa: BLE001
        return list(view["genes"])[:limit]


def _require_gene(view: dict, name: str) -> int:
    """取单个基因列号；缺失 → prerequisite（带候选）。"""
    target = _p_str(name)
    if not target:
        raise prerequisite("该分析需要指定基因参数（例如 gene=\"EPCAM\"）")
    idx = _find_gene(view, target)
    if idx < 0:
        raise prerequisite(
            f"基因「{target}」不在数据集中（共 {len(view['genes'])} 个基因）",
            candidates=_gene_hint(view),
        )
    return idx


def _auto_gene(view: dict, *, why: str, min_expr: float = 0.0) -> tuple[int, str | None]:
    """未指定基因时自动挑表达方差最大的基因（返回 idx 与提示语）。"""
    if not view["genes"]:
        raise prerequisite("数据集中没有可用的基因")
    idx, warn = _pick_variance_genes(view, n_top=1, min_expr=min_expr)
    if not idx:
        raise prerequisite("数据集中没有方差 > 0 的基因（表达矩阵可能为常数），无法挑选基因")
    return idx[0], f"未指定基因，已自动选择方差最大的「{view['genes'][idx[0]]}」（{why}）"


def _pick_variance_genes(view: dict, n_top: int, min_expr: float) -> tuple[list[int], str | None]:
    """按表达方差降序挑候选基因；优先保留平均表达 ≥ min_expr 的基因。"""
    X = np.asarray(view["X"], dtype=np.float64)
    var = np.var(X, axis=0)
    mean = np.mean(X, axis=0)
    usable = np.isfinite(var) & (var > EPS)
    cand = np.where(usable & (mean >= float(min_expr)))[0]
    warn = None
    if cand.size == 0:
        cand = np.where(usable)[0]
        if cand.size:
            warn = f"没有基因的平均表达 ≥ {min_expr:g}，已放宽到全部非恒定基因"
    if cand.size == 0:
        return [], None
    order = cand[np.argsort(var[cand])[::-1]]
    return [int(i) for i in order[: max(1, int(n_top))]], warn


def _constant_genes(view: dict, idx: list[int]) -> tuple[list[int], list[str]]:
    """剔除常数（全 0 / 无方差）基因，返回 (保留, 剔除名)。"""
    X = np.asarray(view["X"], dtype=np.float64)
    keep: list[int] = []
    dropped: list[str] = []
    for j in idx:
        col = X[:, int(j)]
        if not np.isfinite(col).all() or float(np.var(col)) <= EPS:
            dropped.append(view["genes"][int(j)])
        else:
            keep.append(int(j))
    return keep, dropped


def _resolve_gene_set(
    view: dict,
    params: dict,
    *,
    n_top: int,
    min_expr: float,
    default_n: int,
    key: str = "genes",
) -> tuple[list[int], list[str], str | None]:
    """解析基因参数：显式列表优先，否则按方差取前 n_top 个。返回 (idx, names, warn)。"""
    requested = _p_list(params.get(key)) or _p_list(params.get("gene"))
    warn = None
    if requested:
        exact, lowered = _gene_lookup(view["genes"])
        idx: list[int] = []
        missing: list[str] = []
        for name in requested:
            pos = exact.get(name)
            if pos is None:
                pos = lowered.get(name.lower())
            if pos is None:
                missing.append(name)
            elif int(pos) not in idx:
                idx.append(int(pos))
        if not idx:
            raise prerequisite(
                f"指定的 {len(requested)} 个基因都不在数据集中（共 {len(view['genes'])} 个基因）",
                candidates=_gene_hint(view),
            )
        if missing:
            warn = "以下基因不在数据集中，已跳过：" + "、".join(missing[:8])
        idx, dropped = _constant_genes(view, idx)
        if dropped:
            warn = (warn + "；" if warn else "") + "以下基因表达为常数（全 0 或无方差），已跳过：" + "、".join(dropped[:8])
        if not idx:
            raise prerequisite("指定基因的表达全是常数（无方差），无法计算空间自相关")
        return idx, [view["genes"][i] for i in idx], warn

    idx, pick_warn = _pick_variance_genes(view, n_top or default_n, min_expr)
    if pick_warn:
        warn = pick_warn
    if not idx:
        raise prerequisite("没有可用基因（全部为常数或平均表达过低），无法计算空间自相关")
    return idx, [view["genes"][i] for i in idx], warn


# ===========================================================================
# 空间权重
# ===========================================================================

def _neighbor_radius(view: dict, params: dict) -> float:
    """空间邻接半径：参数覆盖 > ctx["neighbor_radius"] > 自动估计。"""
    radius = _p_float(params.get("neighbor_radius"), 0.0)
    if radius > 0:
        return radius
    radius = _p_float(view.get("neighbor_radius"), 0.0)
    if radius > 0:
        return radius
    coords = view["coords"]
    if coords.shape[0] <= DENSE_LIMIT or cKDTree is None:
        return float(estimate_radius(coords))
    tree = cKDTree(coords)
    dist, _ = tree.query(coords, k=2)
    nn = dist[:, 1]
    nn = nn[np.isfinite(nn) & (nn > 0)]
    if nn.size == 0:
        return 100.0
    return float(max(np.median(nn) * 1.2, 1e-6))


def _radius_graph_large(coords: np.ndarray, radius: float) -> "_sparse.csr_matrix":
    """大数据集的半径邻接（cKDTree 等价实现，避免 O(n^2) 距离矩阵）。"""
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=float(radius), output_type="ndarray")
    n = coords.shape[0]
    if pairs.size:
        rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
        cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    else:
        rows = np.empty(0, dtype=int)
        cols = np.empty(0, dtype=int)
    data = np.ones(rows.size, dtype=float)
    graph = _sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    # 与 sp_common.radius_graph 一致：邻居不足 6 个的点补足 kNN
    degree = np.asarray(graph.sum(axis=1)).ravel()
    need = np.where(degree < min(6, n - 1))[0]
    if need.size:
        k = int(min(6, n - 1))
        dist, idx = tree.query(coords[need], k=k + 1)
        dist = np.atleast_2d(dist)
        idx = np.atleast_2d(idx)
        extra: dict[int, list[int]] = {}
        graph = graph.tolil()
        for row_i, point in enumerate(need):
            for col_j in np.atleast_1d(idx[row_i])[1:]:
                graph[int(point), int(col_j)] = 1.0
        graph = graph.tocsr()
    return graph


def _spatial_weights(view: dict, params: dict) -> tuple["_sparse.csr_matrix", float, dict]:
    """构造行标准化空间权重矩阵，返回 (W, radius, 信息 dict)。"""
    coords = view["coords"]
    radius = _neighbor_radius(view, params)
    info: dict = {"radius": float(radius), "fallback": False}
    n = coords.shape[0]
    if n <= DENSE_LIMIT or cKDTree is None or _sparse is None:
        W = spatial_weights(coords, radius)
    else:
        graph = _radius_graph_large(coords, radius)
        row = np.asarray(graph.sum(axis=1)).ravel()
        row[row == 0] = 1.0
        W = _sparse.diags(1.0 / row) @ graph
        info["fallback"] = True
    W = _sparse.csr_matrix(W, dtype=float)
    W.setdiag(0.0)
    W.eliminate_zeros()
    return W, radius, info


def _weight_stats(W, n: int) -> dict:
    """行标准化权重的 S0/S1/S2、Moran's I 的期望与正态近似方差（Cliff & Ord 1981）。"""
    W = _sparse.csr_matrix(W, dtype=float)
    rowsum = np.asarray(W.sum(axis=1)).ravel()
    colsum = np.asarray(W.sum(axis=0)).ravel()
    S0 = float(rowsum.sum())
    sym = _sparse.csr_matrix(W + W.T)
    S1 = 0.5 * float(np.asarray(sym.multiply(sym).sum()))
    S2 = float(np.sum((rowsum + colsum) ** 2))
    n_isolated = int(np.sum(rowsum <= 0))
    expected = -1.0 / (n - 1) if n > 1 else 0.0
    if S0 <= 0:
        variance = float("nan")
    else:
        denom = (n * n - 1.0) * S0 * S0
        variance = (n * n * S1 - n * S2 + 3.0 * S0 * S0) / denom - expected * expected
        if not math.isfinite(variance) or variance <= 0:
            variance = float("nan")
    return {
        "n": int(n),
        "S0": float(S0),
        "S1": float(S1),
        "S2": float(S2),
        "expected": float(expected),
        "variance": float(variance),
        "n_isolated": n_isolated,
        # rowsum 是标准化后的行和（恒为 1），邻居数要用非零元素个数
        "mean_neighbors": float(np.asarray(W.getnnz(axis=1)).mean()) if n else 0.0,
    }


def _submatrix(view: dict, idx: list[int]) -> np.ndarray:
    """按列取子矩阵并转 float64（只复制需要的列）。"""
    X = view["X"]
    return np.ascontiguousarray(np.asarray(X)[:, np.asarray(idx, dtype=int)], dtype=np.float64)


# ===========================================================================
# 核心统计量
# ===========================================================================

def _moran_vector(X_sub: np.ndarray, W, n: int, S0: float) -> np.ndarray:
    """向量化 Moran's I：``I = (n / S0) * (z' W z) / (z' z)``（z 为按列中心化）。"""
    if S0 <= 0:
        return np.full(X_sub.shape[1], np.nan)
    Z = X_sub - X_sub.mean(axis=0, keepdims=True)
    ss = np.einsum("ij,ij->j", Z, Z)
    WZ = np.asarray(W @ Z)
    num = np.einsum("ij,ij->j", Z, WZ)
    with np.errstate(divide="ignore", invalid="ignore"):
        I = np.where(ss > EPS, (float(n) / float(S0)) * num / ss, np.nan)
    return I


def _geary_vector(X_sub: np.ndarray, W, n: int, S0: float) -> np.ndarray:
    """向量化 Geary's C：``C = (n-1)/(2 S0) * Σw_ij(x_i-x_j)^2 / Σ(x_i-x̄)^2``。"""
    if S0 <= 0:
        return np.full(X_sub.shape[1], np.nan)
    Z = X_sub - X_sub.mean(axis=0, keepdims=True)
    ss = np.einsum("ij,ij->j", Z, Z)
    rowsum = np.asarray(W.sum(axis=1)).ravel().astype(float)
    colsum = np.asarray(W.sum(axis=0)).ravel().astype(float)
    quad = ((rowsum + colsum)[:, None] * (Z * Z)).sum(axis=0)
    cross = np.einsum("ij,ij->j", Z, np.asarray(W @ Z))
    edges = np.maximum(quad - 2.0 * cross, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        C = np.where(ss > EPS, (float(n) - 1.0) / (2.0 * float(S0)) * edges / ss, np.nan)
    return C


def _perm_budget(perms: int, n_genes: int, W, warn: list[str], label: str) -> tuple[int, int]:
    """限制置换规模：返回 (有效置换次数, 参与置换的基因数)，必要时写 warning。"""
    perms = max(0, int(perms))
    n_genes = max(1, int(n_genes))
    nnz = int(getattr(W, "nnz", 0)) or 1
    max_genes = n_genes
    if perms * n_genes * nnz > PERM_BUDGET:
        max_genes = max(1, int(PERM_BUDGET // max(perms * nnz, 1)))
        max_genes = min(max_genes, n_genes, 60)
        warn.append(
            f"{label}：置换开销过大（{perms} 次 × {n_genes} 基因），"
            f"仅对前 {max_genes} 个基因做置换检验"
        )
    if perms > 0:
        est = perms * max_genes * nnz
        if est > PERM_BUDGET:
            new_perms = max(10, int(perms * PERM_BUDGET / max(est, 1)))
            if new_perms < perms:
                warn.append(f"{label}：置换次数由 {perms} 降为 {new_perms}（控制计算量）")
                perms = new_perms
    return perms, max_genes


def _moran_permutation(X_sub: np.ndarray, W, n: int, S0: float, observed: np.ndarray,
                       perms: int, rng) -> tuple[np.ndarray, np.ndarray]:
    """标签置换零分布（Moran's I）。返回 (单侧经验 p, 零分布矩阵 perms×g)。"""
    g = X_sub.shape[1]
    null = np.empty((perms, g), dtype=float)
    for b in range(perms):
        perm = rng.permutation(n)
        Z = X_sub[perm]
        Z = Z - Z.mean(axis=0, keepdims=True)
        ss = np.einsum("ij,ij->j", Z, Z)
        num = np.einsum("ij,ij->j", Z, np.asarray(W @ Z))
        with np.errstate(divide="ignore", invalid="ignore"):
            null[b] = np.where(ss > EPS, (float(n) / float(S0)) * num / ss, np.nan)
    # 单侧（正空间自相关）经验 p：置换分布中 ≥ 观测值的比例
    ge = np.nansum(null >= observed[None, :], axis=0)
    pvals = (ge + 1.0) / (float(perms) + 1.0)
    return _clip_p(pvals), null


def _geary_permutation(X_sub: np.ndarray, W, n: int, S0: float, observed: np.ndarray,
                       perms: int, rng) -> tuple[np.ndarray, np.ndarray]:
    """标签置换零分布（Geary's C）。返回 (单侧经验 p, 零分布矩阵)。"""
    g = X_sub.shape[1]
    null = np.empty((perms, g), dtype=float)
    for b in range(perms):
        perm = rng.permutation(n)
        null[b] = _geary_vector(X_sub[perm], W, n, S0)
    le = np.nansum(null <= observed[None, :], axis=0)
    pvals = (le + 1.0) / (float(perms) + 1.0)
    return _clip_p(pvals), null


# ===========================================================================
# 绘图小工具（图上一律英文）
# ===========================================================================

def _point_size(view: dict, params: dict, default: float = 6.0) -> float:
    """按 spot 数自适应点大小（避免大盘点糊成一片）。"""
    requested = _p_float(params.get("point_size"), 0.0)
    if requested > 0:
        return requested
    n = max(view["n"], 1)
    if n > 40000:
        return 1.6
    if n > 15000:
        return 3.0
    if n > 6000:
        return 4.5
    return float(default)


def _spatial_ctx(view: dict, use_image: bool) -> dict:
    """给 sp_common.scatter_spatial 用的最小 ctx。"""
    return {"coords": view["coords"], "image": view["image"] if use_image else {}}


def _style_map_axis(ax, title: str | None = None) -> None:
    if title:
        ax.set_title(title, fontsize=9.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def _render_svg_map(view: dict, idx: list[int], names: list[str], outdir: str,
                    params: dict) -> dict:
    """Top SVG 空间表达图（2 列网格）。"""
    import matplotlib.pyplot as plt

    k = len(idx)
    ncol = 2 if k > 1 else 1
    nrow = int(math.ceil(k / ncol))
    width = 4.9 * ncol
    height = 4.4 * nrow
    fig, axes = plt.subplots(nrow, ncol, figsize=(width, height), squeeze=False)
    use_image = _p_bool(params.get("image"), False)
    ps = _point_size(view, params)
    X = view["X"]
    for pos in range(nrow * ncol):
        ax = axes[pos // ncol][pos % ncol]
        if pos >= k:
            ax.axis("off")
            continue
        j = int(idx[pos])
        values = np.asarray(X[:, j], dtype=float)
        try:
            scatter_spatial(ax, _spatial_ctx(view, use_image), values=values, point_size=ps,
                            cmap="viridis", legend_title="expression")
        except EngineError:
            raise
        except Exception as error:  # noqa: BLE001
            raise EngineError("internal", f"绘制「{names[pos]}」空间图失败：{error}")
        _style_map_axis(ax, names[pos])
    fig.suptitle("Top spatially variable genes", fontsize=11, y=0.998)
    return save_fig(fig, outdir, "svg_map", "Top SVG expression maps")


def _render_moran_plot(frame: pd.DataFrame, outdir: str, title: str = "Moran's I ranking") -> dict:
    """Moran's I 排序条形图（含期望值参考线、显著性着色）。"""
    import matplotlib.pyplot as plt

    plot = frame.dropna(subset=["moran_i"]).copy()
    if plot.empty:
        raise prerequisite("没有可绘制的 Moran's I（全部基因计算失败）")
    plot = plot.sort_values("moran_i", ascending=True)
    max_bars = 80
    if len(plot) > max_bars:
        plot = plot.tail(max_bars)
    k = len(plot)
    height = max(3.2, 0.19 * k + 1.5)
    fig, ax = plt.subplots(figsize=(7.2, min(height, 22.0)))
    ypos = np.arange(k)
    sig = plot["p_adj"].fillna(1.0).to_numpy() < 0.05 if "p_adj" in plot else np.zeros(k, dtype=bool)
    colors = ["#c0392b" if s else "#8c9bab" for s in sig]
    ax.barh(ypos, plot["moran_i"].to_numpy(), color=colors, height=0.72)
    ax.set_yticks(ypos)
    ax.set_yticklabels([str(g) for g in plot["gene"]], fontsize=6.4)
    expected = float(plot["expected"].iloc[0]) if "expected" in plot else float("-inf")
    if math.isfinite(expected):
        ax.axvline(expected, color="#2b2f36", linewidth=0.9, linestyle="--",
                   label=f"E[I] = {expected:.4f}")
        ax.legend(loc="lower right")
    ax.axvline(0.0, color="#c9ced6", linewidth=0.7)
    ax.set_xlabel("Moran's I")
    ax.set_title(title)
    ax.grid(axis="y", visible=False)
    return save_fig(fig, outdir, "moran_plot", "Moran's I ranking")


def _render_geary_plot(frame: pd.DataFrame, outdir: str) -> dict:
    """Geary's C 排序图（1-C 越大越聚集）。"""
    import matplotlib.pyplot as plt

    plot = frame.dropna(subset=["geary_c"]).copy()
    if plot.empty:
        raise prerequisite("没有可绘制的 Geary's C（全部基因计算失败）")
    plot["signal"] = 1.0 - plot["geary_c"]
    plot = plot.sort_values("signal", ascending=True).tail(60)
    k = len(plot)
    fig, ax = plt.subplots(figsize=(7.2, min(max(3.2, 0.19 * k + 1.5), 22.0)))
    ypos = np.arange(k)
    sig = plot["p_adj"].fillna(1.0).to_numpy() < 0.05 if "p_adj" in plot else np.zeros(k, dtype=bool)
    ax.barh(ypos, plot["signal"].to_numpy(), color=["#2e6da4" if s else "#8c9bab" for s in sig], height=0.72)
    ax.set_yticks(ypos)
    ax.set_yticklabels([str(g) for g in plot["gene"]], fontsize=6.4)
    ax.axvline(0.0, color="#2b2f36", linewidth=0.8, linestyle="--", label="E[C] = 1")
    ax.set_xlabel("1 - Geary's C  (higher = more spatial clustering)")
    ax.set_title("Geary's C ranking")
    ax.legend(loc="lower right")
    ax.grid(axis="y", visible=False)
    return save_fig(fig, outdir, "geary_plot", "Geary's C ranking")


def _render_permutation(null: np.ndarray, observed: np.ndarray, gene: str, outdir: str,
                        stat_name: str = "Moran's I") -> dict:
    """置换零分布直方图 + 观测值参考线。"""
    import matplotlib.pyplot as plt

    values = np.asarray(null, dtype=float).ravel()
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise prerequisite("置换零分布为空，无法绘图")
    fig, ax = new_fig(7.0, 4.2)
    ax.hist(values, bins=min(60, max(12, int(math.sqrt(values.size)))), color="#8fb4d9",
            edgecolor="white", linewidth=0.5, density=True)
    obs = float(np.asarray(observed, dtype=float).ravel()[0])
    if math.isfinite(obs):
        ax.axvline(obs, color="#c0392b", linewidth=1.8, label=f"observed = {obs:.4f}")
    ax.set_xlabel(f"{stat_name} under label permutation")
    ax.set_ylabel("Density")
    ax.set_title(f"Permutation null: {gene}")
    ax.legend()
    return save_fig(fig, outdir, "permutation_null", "Permutation null distribution")


def _hotspot_categories(z: np.ndarray, z_thresh: float) -> np.ndarray:
    cats = np.full(z.shape[0], "not significant", dtype=object)
    cats[z > z_thresh] = "hot spot"
    cats[z < -z_thresh] = "cold spot"
    return cats


def _render_hotspot_map(view: dict, cats: np.ndarray, values: np.ndarray, gene: str,
                        outdir: str, params: dict) -> dict:
    """Getis-Ord Gi* 热点图（hot 红 / cold 蓝 / 不显著灰）。"""
    import matplotlib.pyplot as plt

    use_image = _p_bool(params.get("image"), False)
    ps = _point_size(view, params)
    fig, ax = new_fig(6.6, 5.6)
    if use_image:
        from sp_common import add_image_background

        add_image_background(ax, _spatial_ctx(view, True))
    coords = view["coords"]
    order = [("not significant", NS_COLOR, 0.55, ps * 0.75),
             ("cold spot", COLD_COLOR, 0.95, ps),
             ("hot spot", HOT_COLOR, 0.95, ps)]
    for label, color, alpha, size in order:
        sel = np.array([c == label for c in cats])
        if not sel.any():
            continue
        ax.scatter(coords[sel, 0], coords[sel, 1], s=size, c=color, alpha=alpha,
                   linewidths=0, label=f"{label} (n={int(sel.sum())})")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    _style_map_axis(ax, f"Getis-Ord Gi*: {gene}")
    ax.legend(loc="upper right", fontsize=8, markerscale=1.8, title="Gi* class", title_fontsize=8)
    return save_fig(fig, outdir, "hotspot_map", "Getis-Ord Gi* hotspot map")


def _render_gradient_plot(bins: pd.DataFrame, gene: str, axis_label: str, slope: float,
                          pval: float, outdir: str) -> dict:
    """空间梯度：分箱均值 ± SE + 线性趋势。"""
    import matplotlib.pyplot as plt

    fig, ax = new_fig(7.2, 4.4)
    t = bins["t_center"].to_numpy(dtype=float)
    mean = bins["mean"].to_numpy(dtype=float)
    se = bins["se"].fillna(0.0).to_numpy(dtype=float)
    ax.errorbar(t, mean, yerr=se, fmt="o-", color="#2e6da4", markersize=4, linewidth=1.4,
                capsize=2.5, ecolor="#8fb4d9", elinewidth=1.0, label="binned mean ± SE")
    if len(bins) >= 2 and math.isfinite(slope):
        xs = np.linspace(np.nanmin(t), np.nanmax(t), 50)
        base = float(np.nanmean(mean - slope * (t - np.nanmean(t))))
        ax.plot(xs, base + slope * xs, color="#c0392b", linewidth=1.3, linestyle="--",
                label=f"linear trend (slope={slope:.3g}, p={pval:.3g})")
    ax.set_xlabel(axis_label)
    ax.set_ylabel("expression")
    ax.set_title(f"Spatial gradient: {gene}")
    ax.legend(fontsize=8)
    return save_fig(fig, outdir, "gradient_plot", "Spatial gradient")


def _render_corr_pairs(pairs: pd.DataFrame, genes: list[str], outdir: str) -> dict:
    """空间共表达：基因网络（top 对）+ Pearson–邻域一致性散点。"""
    import matplotlib.pyplot as plt

    top = pairs.head(min(28, len(pairs)))
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.8))

    ax = axes[0]
    nodes = list(dict.fromkeys(list(top["gene_a"]) + list(top["gene_b"])))
    if not nodes:
        raise prerequisite("没有可绘制的基因对")
    n_nodes = len(nodes)
    angles = np.linspace(0, 2 * math.pi, max(n_nodes, 3), endpoint=False)
    pos = {name: (math.cos(a), math.sin(a)) for name, a in zip(nodes, angles)}
    strength = {name: 0.0 for name in nodes}
    coh = top["neighbor_coherence"].fillna(0.0).to_numpy(dtype=float)
    lo, hi = float(np.nanmin(coh)), float(np.nanmax(coh))
    span = max(hi - lo, 1e-6)
    for (_, row), value in zip(top.iterrows(), coh):
        a, b = str(row["gene_a"]), str(row["gene_b"])
        if a not in pos or b not in pos:
            continue
        weight = 0.6 + 3.2 * (value - lo) / span
        color = "#c0392b" if value >= 0 else "#2e6da4"
        ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color=color,
                linewidth=weight, alpha=0.55, zorder=1)
        strength[a] += abs(value)
        strength[b] += abs(value)
    max_strength = max(strength.values()) if strength else 1.0
    for name in nodes:
        x, y = pos[name]
        size = 90 + 420 * (strength[name] / max(max_strength, 1e-9))
        ax.scatter([x], [y], s=size, color="#4e79a7", edgecolor="white", linewidth=0.8, zorder=2)
        ha = "left" if x >= 0 else "right"
        ax.text(x * 1.13, y * 1.13, name, fontsize=6.2, ha=ha, va="center")
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Top spatial co-expression pairs (edge ∝ neighbourhood coherence)", fontsize=9.5)

    ax = axes[1]
    sc = ax.scatter(pairs["pearson"], pairs["neighbor_coherence"], s=12,
                    c=pairs["neighbor_coherence"], cmap="coolwarm", linewidths=0, alpha=0.85)
    try:
        fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02).ax.tick_params(labelsize=7)
    except Exception:  # noqa: BLE001
        pass
    for _, row in top.head(6).iterrows():
        ax.annotate(f"{row['gene_a']}–{row['gene_b']}", (row["pearson"], row["neighbor_coherence"]),
                    fontsize=5.8, xytext=(3, 3), textcoords="offset points")
    ax.axhline(0.0, color="#c9ced6", linewidth=0.7)
    ax.axvline(0.0, color="#c9ced6", linewidth=0.7)
    ax.set_xlabel("Pearson r (global co-expression)")
    ax.set_ylabel("Neighbourhood product coherence")
    ax.set_title("Global vs spatial co-expression", fontsize=9.5)
    return save_fig(fig, outdir, "spatial_corr_network", "Spatial co-expression network")


def _render_ripley(curves: dict, outdir: str) -> dict:
    """Ripley K/L 曲线：每个簇 K(r) 与 CSR 均值，L(r)-r 偏差与包络。"""
    import matplotlib.pyplot as plt

    radii = np.asarray(curves["radii"], dtype=float)
    items = curves["clusters"]
    if not items:
        raise prerequisite("没有可计算 Ripley's K 的簇（每簇至少需要 5 个 spot）")
    colors = palette(len(items))
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.2))
    ax = axes[0]
    for color, item in zip(colors, items):
        ax.plot(radii, item["k"], color=color, linewidth=1.7, label=f"{item['label']} (n={item['n']})")
        ax.plot(radii, item["k_csr_mean"], color=color, linewidth=0.9, linestyle="--", alpha=0.75)
    ax.set_xlabel("radius r")
    ax.set_ylabel("Ripley's K(r)")
    ax.set_title("Ripley's K: observed (solid) vs CSR (dashed)")
    ax.legend(fontsize=7)

    ax = axes[1]
    for color, item in zip(colors, items):
        ax.plot(radii, item["l_dev"], color=color, linewidth=1.7, label=str(item["label"]))
    if items:
        first = items[0]
        ax.fill_between(radii, first["l_dev_lo"], first["l_dev_hi"], color="#c9ced6", alpha=0.32,
                        linewidth=0, label=f"CSR 95% envelope ({first['label']})")
    ax.axhline(0.0, color="#2b2f36", linewidth=0.9)
    ax.set_xlabel("radius r")
    ax.set_ylabel("L(r) - r  (observed - CSR)")
    ax.set_title("Spatial clustering deviation")
    ax.legend(fontsize=7)
    return save_fig(fig, outdir, "ripley", "Ripley's K / L curves")


def _render_autocorr_scan(frame: pd.DataFrame, expected: float, outdir: str) -> dict:
    """全基因扫描：I 分布直方图 + I vs -log10(FDR)。"""
    import matplotlib.pyplot as plt

    values = frame["moran_i"].dropna().to_numpy(dtype=float)
    if values.size == 0:
        raise prerequisite("全基因扫描没有有效结果，无法绘图")
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6))
    ax = axes[0]
    bins = int(min(80, max(20, math.sqrt(values.size) * 1.6)))
    ax.hist(values, bins=bins, color="#8fb4d9", edgecolor="white", linewidth=0.4)
    ax.axvline(expected, color="#2b2f36", linewidth=1.0, linestyle="--",
               label=f"E[I] = {expected:.4f}")
    ax.set_xlabel("Moran's I")
    ax.set_ylabel("Number of genes")
    ax.set_title("Distribution of Moran's I (all genes)")
    ax.legend(fontsize=8)

    ax = axes[1]
    y = -np.log10(np.clip(frame["p_adj"].to_numpy(dtype=float), 1e-300, 1.0))
    ax.scatter(frame["moran_i"], y, s=8, c="#4e79a7", alpha=0.6, linewidths=0)
    sig = frame["p_adj"] < 0.05
    ax.scatter(frame.loc[sig, "moran_i"], y[sig.to_numpy()], s=10, c="#c0392b", alpha=0.8, linewidths=0,
               label=f"FDR < 0.05 (n={int(sig.sum())})")
    ax.axhline(-math.log10(0.05), color="#2b2f36", linewidth=0.8, linestyle="--")
    ax.axvline(0.0, color="#c9ced6", linewidth=0.7)
    ax.set_xlabel("Moran's I")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Spatially variable genes")
    ax.legend(fontsize=8)
    return save_fig(fig, outdir, "autocorr_scan", "Genome-wide spatial autocorrelation scan")


def _render_cross_moran_map(view: dict, ia: int, ib: int, values_b_lag: np.ndarray,
                            names: tuple[str, str], bivar: dict, outdir: str,
                            params: dict) -> dict:
    """双变量 Moran：两基因空间图 + Moran 散点（za vs 邻居平均 zb）。"""
    import matplotlib.pyplot as plt

    X = view["X"]
    use_image = _p_bool(params.get("image"), False)
    ps = _point_size(view, params)
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.0))
    for ax, j, name in ((axes[0], ia, names[0]), (axes[1], ib, names[1])):
        scatter_spatial(ax, _spatial_ctx(view, use_image),
                        values=np.asarray(X[:, j], dtype=float), point_size=ps,
                        cmap="viridis", legend_title="expression")
        _style_map_axis(ax, name)

    ax = axes[2]
    za = np.asarray(bivar["za"], dtype=float)
    lag = np.asarray(values_b_lag, dtype=float)
    ax.scatter(za, lag, s=max(4.0, ps * 0.9), c="#4e79a7", alpha=0.55, linewidths=0)
    if za.size >= 3:
        slope, intercept = np.polyfit(za, lag, 1)
        xs = np.linspace(float(np.min(za)), float(np.max(za)), 50)
        ax.plot(xs, intercept + slope * xs, color="#c0392b", linewidth=1.4,
                label=f"slope = {float(bivar['moran_ab']):.3f}")
        ax.legend(fontsize=8)
    ax.axhline(0.0, color="#c9ced6", linewidth=0.7)
    ax.axvline(0.0, color="#c9ced6", linewidth=0.7)
    ax.set_xlabel(f"z({names[0]})")
    ax.set_ylabel(f"neighbour mean z({names[1]})")
    ax.set_title("Bivariate Moran scatter")
    return save_fig(fig, outdir, "cross_moran_map", "Bivariate Moran spatial maps")


# ===========================================================================
# 1. morans_svg — Moran's I
# ===========================================================================

_MORAN_COLUMNS = [
    {"key": "gene", "label": "基因", "type": "string"},
    {"key": "moran_i", "label": "Moran's I", "type": "number"},
    {"key": "expected", "label": "E[I]", "type": "number"},
    {"key": "z", "label": "z 值", "type": "number"},
    {"key": "p", "label": "解析 p", "type": "number"},
    {"key": "p_perm", "label": "置换 p", "type": "number"},
    {"key": "p_adj", "label": "FDR", "type": "number"},
    {"key": "mean", "label": "平均表达", "type": "number"},
    {"key": "n_spots", "label": "spot 数", "type": "number"},
]


def _moran_frame(view: dict, idx: list[int], names: list[str], W, radius: float,
                 stats: dict, perms: int, rng, warn: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    """计算 Moran's I 表（解析 p + 可选置换 p），返回 (表, 置换零分布)。"""
    X_sub = _submatrix(view, idx)
    n = view["n"]
    S0 = float(stats["S0"])
    if S0 <= 0:
        raise prerequisite("空间权重矩阵全为 0（spot 坐标过密或邻接半径过小），无法计算空间自相关")
    observed = _moran_vector(X_sub, W, n, S0)
    expected = float(stats["expected"])
    variance = float(stats["variance"])
    if math.isfinite(variance) and variance > 0:
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (observed - expected) / math.sqrt(variance)
        pvals = 2.0 * _norm_sf(np.abs(z))
    else:
        z = np.full(observed.shape, np.nan)
        pvals = np.ones(observed.shape)
        warn.append("正态近似方差不可用（权重结构退化），解析 z/p 记为缺失；请参考置换 p")
    pvals = np.where(np.isfinite(observed), _clip_p(pvals), np.nan)

    perms_eff, perm_genes = _perm_budget(perms, len(idx), W, warn, "Moran's I")
    p_perm = np.full(observed.shape, np.nan)
    null = np.zeros((0, len(idx)))
    if perms_eff > 0 and perm_genes > 0:
        p_sub, null = _moran_permutation(X_sub[:, :perm_genes], W, n, S0, observed[:perm_genes],
                                         perms_eff, rng)
        p_perm[:perm_genes] = p_sub
        if perm_genes < len(idx):
            warn.append(f"置换 p 只覆盖前 {perm_genes} 个基因，其余基因 p_perm 为空")

    means = X_sub.mean(axis=0)
    frame = pd.DataFrame({
        "gene": names,
        "moran_i": observed,
        "expected": expected,
        "z": z,
        "p": pvals,
        "p_perm": p_perm,
        "mean": means,
        "n_spots": int(n),
        "neighbor_radius": float(radius),
    })
    adj_source = np.where(np.isfinite(p_perm), p_perm, pvals)
    frame["p_adj"] = bh_fdr(adj_source)
    frame = frame.sort_values("moran_i", ascending=False, na_position="last").reset_index(drop=True)
    return frame, null


def _norm_sf(x) -> np.ndarray:
    """标准正态生存函数（避免依赖 scipy.stats 的导入开销/异常）。"""
    from scipy import stats

    return np.asarray(stats.norm.sf(np.asarray(x, dtype=float)), dtype=float)


def run_morans_svg(ctx: dict, params: dict, outdir: str) -> dict:
    """空间自相关（Moran's I）：Moran 1950 + 正态近似 / 标签置换（内置实现）。

    行标准化邻接权重（半径邻接）下 ``I = (n/S0)·(z'Wz)/(z'z)``；期望 ``-1/(n-1)``，
    方差用 Cliff & Ord (1981) 正态近似；``permutations>0`` 时另给单侧标签置换经验 p
    ``(k+1)/(n+1)``（检验正空间自相关）。
    """
    ensure_dir(outdir)
    params = dict(params or {})
    view = _view(ctx)
    warn = list(view["warnings"])
    n_top = max(1, min(2000, _p_int(params.get("n_top"), 20)))
    min_expr = max(0.0, _p_float(params.get("min_expr"), 0.1))
    perms = max(0, min(5000, _p_int(params.get("permutations"), 100)))
    idx, names, pick_warn = _resolve_gene_set(view, params, n_top=n_top, min_expr=min_expr, default_n=20)
    if pick_warn:
        warn.append(pick_warn)
    W, radius, winfo = _spatial_weights(view, params)
    if winfo.get("fallback"):
        warn.append("spot 数较多，空间权重改用 cKDTree 等价半径图（与 sp_common.spatial_weights 同规则）")
    stats = _weight_stats(W, view["n"])
    if stats["n_isolated"] > 0:
        warn.append(
            f"{stats['n_isolated']} 个 spot 没有任何邻居（权重行和为 0），只参与方差、不参与协方差，可能低估 Moran's I"
        )
    frame, null = _moran_frame(view, idx, names, W, radius, stats, perms, np.random.default_rng(DEFAULT_SEED), warn)
    frame_path = export_csv(frame, outdir, "morans_svg", index=False)

    figures: list[dict] = []
    top_rows = frame.dropna(subset=["moran_i"]).head(4)
    if len(top_rows):
        try:
            figures.append(_render_svg_map(view, [int(np.where(np.asarray(names) == g)[0][0]) for g in top_rows["gene"]],
                                           [str(g) for g in top_rows["gene"]], outdir, params))
        except Exception as error:  # noqa: BLE001
            warn.append(f"svg_map 绘制失败：{error}")
    try:
        figures.append(_render_moran_plot(frame, outdir))
    except Exception as error:  # noqa: BLE001
        warn.append(f"moran_plot 绘制失败：{error}")

    sig = frame["p_adj"].fillna(1.0) < 0.05
    n_sig_pos = int(((frame["p_adj"].fillna(1.0) < 0.05) & (frame["moran_i"] > 0)).sum())
    bullets = [
        f"对 {len(frame)} 个基因计算 Moran's I（行标准化半径邻接，半径 {radius:.3g}，"
        f"{view['n']} 个 spot，平均邻居 {stats['mean_neighbors']:.2f}）",
        f"FDR<0.05 且 I>0 的显著空间高变基因：{n_sig_pos} 个"
        + (f"，最强为「{frame.iloc[0]['gene']}」（I={frame.iloc[0]['moran_i']:.3f}）" if len(frame) else ""),
    ]
    if perms > 0:
        bullets.append(f"置换检验 {perms} 次（单侧正自相关经验 p=(k+1)/(n+1)，seed=0）")
    metrics = {
        "n_genes": int(len(frame)),
        "n_significant": n_sig_pos,
        "top_gene": str(frame.iloc[0]["gene"]) if len(frame) else None,
        "top_moran_i": float(frame.iloc[0]["moran_i"]) if len(frame) else None,
        "n_spots": int(view["n"]),
        "neighbor_radius": float(radius),
        "permutations": int(perms),
        "n_isolated": int(stats["n_isolated"]),
    }
    table_obj = table("svg", "空间高变基因（Moran's I）", _MORAN_COLUMNS, frame,
                      note="I>0 表示空间聚集；p 为解析正态近似双尾，p_perm 为单侧标签置换经验 p；p_adj 用二者中可得者做 BH")
    payload = result(
        title="空间自相关（Moran's I）",
        bullets=bullets,
        metrics=metrics,
        tables=[table_obj],
        figures=figures,
        warnings=warn,
        notes=[f"方法：Moran's I（Moran 1950），行标准化邻接权重，内置实现；半径 {radius:.3g}。"],
        data_out=[data_file("svg_csv", frame_path, "Moran's I 表")],
    )
    payload["store"] = {
        "morans_svg": {
            "genes": [str(g) for g in frame["gene"].head(200)],
            "table": jsonable(frame.head(500)),
            "params": {"n_top": n_top, "permutations": int(perms), "min_expr": min_expr},
            "neighbor_radius": float(radius),
        }
    }
    return payload


# ===========================================================================
# 2. gearys_svg — Geary's C
# ===========================================================================

_GEARY_COLUMNS = [
    {"key": "gene", "label": "基因", "type": "string"},
    {"key": "geary_c", "label": "Geary's C", "type": "number"},
    {"key": "expected", "label": "E[C]", "type": "number"},
    {"key": "z", "label": "置换 z", "type": "number"},
    {"key": "p_perm", "label": "置换 p", "type": "number"},
    {"key": "p_adj", "label": "FDR", "type": "number"},
    {"key": "mean", "label": "平均表达", "type": "number"},
    {"key": "n_spots", "label": "spot 数", "type": "number"},
]


def run_gearys_svg(ctx: dict, params: dict, outdir: str) -> dict:
    """空间自相关（Geary's C）：Geary 1954 + 标签置换经验 p（内置实现）。

    ``C = (n-1)/(2S0) · Σw_ij(x_i-x_j)² / Σ(x_i-x̄)²``；E[C]=1，C<1 表示相邻
    spot 表达接近（空间聚集）。p 用标签置换零分布（单侧 C 偏小）。
    """
    ensure_dir(outdir)
    params = dict(params or {})
    view = _view(ctx)
    warn = list(view["warnings"])
    n_top = max(1, min(2000, _p_int(params.get("n_top"), 20)))
    min_expr = max(0.0, _p_float(params.get("min_expr"), 0.1))
    perms = max(0, min(5000, _p_int(params.get("permutations"), 99)))
    idx, names, pick_warn = _resolve_gene_set(view, params, n_top=n_top, min_expr=min_expr, default_n=20)
    if pick_warn:
        warn.append(pick_warn)
    W, radius, winfo = _spatial_weights(view, params)
    if winfo.get("fallback"):
        warn.append("spot 数较多，空间权重改用 cKDTree 等价半径图（与 sp_common.spatial_weights 同规则）")
    stats = _weight_stats(W, view["n"])
    if stats["S0"] <= 0:
        raise prerequisite("空间权重矩阵全为 0（spot 坐标过密或邻接半径过小），无法计算空间自相关")
    if stats["n_isolated"] > 0:
        warn.append(f"{stats['n_isolated']} 个 spot 没有邻居（权重行和为 0），只参与方差、不参与协方差")

    X_sub = _submatrix(view, idx)
    n = view["n"]
    S0 = float(stats["S0"])
    observed = _geary_vector(X_sub, W, n, S0)
    perms_eff, perm_genes = _perm_budget(perms, len(idx), W, warn, "Geary's C")
    z = np.full(observed.shape, np.nan)
    p_perm = np.full(observed.shape, np.nan)
    if perms_eff > 0 and perm_genes > 0:
        p_sub, null = _geary_permutation(X_sub[:, :perm_genes], W, n, S0, observed[:perm_genes],
                                         perms_eff, np.random.default_rng(DEFAULT_SEED))
        p_perm[:perm_genes] = p_sub
        with np.errstate(invalid="ignore"):
            mean_null = np.nanmean(null, axis=0)
            sd_null = np.nanstd(null, axis=0, ddof=1) if null.shape[0] > 1 else np.full(perm_genes, np.nan)
            z[:perm_genes] = np.where(sd_null > EPS, (observed[:perm_genes] - mean_null) / sd_null, np.nan)
        if perm_genes < len(idx):
            warn.append(f"置换 p 只覆盖前 {perm_genes} 个基因，其余基因 p_perm 为空")
    else:
        warn.append("permutations=0：Geary's C 只给统计量，不给经验 p")

    frame = pd.DataFrame({
        "gene": names,
        "geary_c": observed,
        "expected": 1.0,
        "z": z,
        "p_perm": p_perm,
        "mean": X_sub.mean(axis=0),
        "n_spots": int(n),
        "neighbor_radius": float(radius),
    })
    frame["p_adj"] = bh_fdr(frame["p_perm"])
    frame = frame.sort_values("geary_c", ascending=True, na_position="last").reset_index(drop=True)
    frame_path = export_csv(frame, outdir, "gearys_svg", index=False)

    figures: list[dict] = []
    try:
        figures.append(_render_geary_plot(frame, outdir))
    except Exception as error:  # noqa: BLE001
        warn.append(f"geary_plot 绘制失败：{error}")
    if len(frame) and frame["p_perm"].notna().any():
        try:
            top = frame.dropna(subset=["p_perm", "geary_c"]).iloc[0]
            # 为直方图重算一次零分布（仅当置换基因数不多时，成本可接受）
            if perms_eff > 0 and perm_genes > 0:
                sub = _submatrix(view, idx[:perm_genes])
                _, null = _geary_permutation(sub, W, n, S0, _geary_vector(sub, W, n, S0), perms_eff,
                                             np.random.default_rng(DEFAULT_SEED))
                pos = names.index(str(top["gene"]))
                if pos < perm_genes:
                    figures.append(_render_permutation(null[:, pos], np.array([float(top["geary_c"])]),
                                                       str(top["gene"]), outdir, stat_name="Geary's C"))
        except Exception as error:  # noqa: BLE001
            warn.append(f"置换零分布图绘制失败：{error}")

    n_sig = int(((frame["p_adj"].fillna(1.0) < 0.05) & (frame["geary_c"] < 1.0)).sum())
    bullets = [
        f"对 {len(frame)} 个基因计算 Geary's C（行标准化半径邻接，半径 {radius:.3g}，{n} 个 spot）",
        f"C<1 且 FDR<0.05 的显著空间聚集基因：{n_sig} 个"
        + (f"，最小 C 为「{frame.iloc[0]['gene']}」（C={frame.iloc[0]['geary_c']:.3f}）" if len(frame) else ""),
        f"经验 p 来自 {perms} 次标签置换（单侧 C 偏小，p=(k+1)/(n+1)，seed=0）",
    ]
    metrics = {
        "n_genes": int(len(frame)),
        "n_significant": n_sig,
        "min_geary_c": float(frame["geary_c"].min()) if len(frame) else None,
        "min_gene": str(frame.iloc[0]["gene"]) if len(frame) else None,
        "n_spots": int(n),
        "neighbor_radius": float(radius),
        "permutations": int(perms),
    }
    payload = result(
        title="空间自相关（Geary's C）",
        bullets=bullets,
        metrics=metrics,
        tables=[table("geary", "空间高变基因（Geary's C）", _GEARY_COLUMNS, frame,
                      note="C<1 表示相邻 spot 表达接近（空间聚集），E[C]=1；p 为单侧标签置换经验 p")],
        figures=figures,
        warnings=warn,
        notes=[f"方法：Geary's C（Geary 1954），行标准化邻接权重，内置实现；半径 {radius:.3g}。"],
        data_out=[data_file("geary_csv", frame_path, "Geary's C 表")],
    )
    payload["store"] = {"gearys_svg": {"genes": [str(g) for g in frame["gene"].head(200)]}}
    return payload


# ===========================================================================
# 3. svg_permutation — 纯标签置换
# ===========================================================================

def run_svg_permutation(ctx: dict, params: dict, outdir: str) -> dict:
    """空间置换检验：纯标签置换经验 p（SPARK-X 思想，Sun et al. 2021；内置实现）。

    不依赖解析分布：把表达值在 spot 间随机重排 ``permutations`` 次，统计
    ``I_perm ≥ I_obs`` 的比例得到单侧经验 p ``(k+1)/(n+1)``，对候选基因做 BH。
    """
    ensure_dir(outdir)
    params = dict(params or {})
    view = _view(ctx)
    warn = list(view["warnings"])
    perms = max(1, min(5000, _p_int(params.get("permutations"), 200)))
    n_top = max(1, min(2000, _p_int(params.get("n_top"), 20)))
    min_expr = max(0.0, _p_float(params.get("min_expr"), 0.1))
    idx, names, pick_warn = _resolve_gene_set(view, params, n_top=n_top, min_expr=min_expr, default_n=20)
    if pick_warn:
        warn.append(pick_warn)
    W, radius, winfo = _spatial_weights(view, params)
    if winfo.get("fallback"):
        warn.append("spot 数较多，空间权重改用 cKDTree 等价半径图（与 sp_common.spatial_weights 同规则）")
    stats = _weight_stats(W, view["n"])
    if stats["S0"] <= 0:
        raise prerequisite("空间权重矩阵全为 0（spot 坐标过密或邻接半径过小），无法做置换检验")
    perms_eff, perm_genes = _perm_budget(perms, len(idx), W, warn, "标签置换")
    if perm_genes < len(idx):
        idx = idx[:perm_genes]
        names = names[:perm_genes]
    X_sub = _submatrix(view, idx)
    n = view["n"]
    S0 = float(stats["S0"])
    observed = _moran_vector(X_sub, W, n, S0)
    p_perm, null = _moran_permutation(X_sub, W, n, S0, observed, perms_eff, np.random.default_rng(DEFAULT_SEED))
    frame = pd.DataFrame({
        "gene": names,
        "moran_i": observed,
        "expected": float(stats["expected"]),
        "p_perm": p_perm,
        "mean": X_sub.mean(axis=0),
        "n_spots": int(n),
    })
    frame["p_adj"] = bh_fdr(frame["p_perm"])
    frame = frame.sort_values(["p_perm", "moran_i"], ascending=[True, False]).reset_index(drop=True)
    frame_path = export_csv(frame, outdir, "svg_permutation", index=False)

    figures: list[dict] = []
    if null.shape[0] > 0 and len(frame):
        try:
            figures.append(_render_permutation(null[:, 0], observed[:1], str(frame.iloc[0]["gene"]),
                                               outdir, stat_name="Moran's I"))
        except Exception as error:  # noqa: BLE001
            warn.append(f"置换零分布图绘制失败：{error}")
    try:
        plot_frame = frame.copy()
        plot_frame["p_adj"] = plot_frame["p_perm"]
        figures.append(_render_moran_plot(plot_frame, outdir, title="Permutation test: Moran's I"))
    except Exception as error:  # noqa: BLE001
        warn.append(f"置换检验条形图绘制失败：{error}")

    n_sig = int((frame["p_adj"] < 0.05).sum())
    bullets = [
        f"对 {len(frame)} 个基因做 {perms_eff} 次标签置换（{n} 个 spot，半径 {radius:.3g}）",
        f"经验 p<0.05 的基因：{n_sig} 个"
        + (f"，最小 p={frame.iloc[0]['p_perm']:.4f}（「{frame.iloc[0]['gene']}」）" if len(frame) else ""),
        "经验 p=(k+1)/(n+1)，k 为置换分布中 ≥ 观测 Moran's I 的次数（单侧正自相关）",
    ]
    metrics = {
        "n_genes": int(len(frame)),
        "n_significant": n_sig,
        "permutations": int(perms_eff),
        "min_p_perm": float(frame["p_perm"].min()) if len(frame) else None,
        "n_spots": int(n),
        "neighbor_radius": float(radius),
    }
    payload = result(
        title="空间置换检验",
        bullets=bullets,
        metrics=metrics,
        tables=[table("svg_perm", "空间置换检验（经验 p）",
                      [{"key": "gene", "label": "基因", "type": "string"},
                       {"key": "moran_i", "label": "Moran's I", "type": "number"},
                       {"key": "expected", "label": "E[I]", "type": "number"},
                       {"key": "p_perm", "label": "置换 p", "type": "number"},
                       {"key": "p_adj", "label": "FDR", "type": "number"},
                       {"key": "mean", "label": "平均表达", "type": "number"},
                       {"key": "n_spots", "label": "spot 数", "type": "number"}],
                      frame, note="p_perm 为单侧经验 p=(k+1)/(n+1)；FDR 为 BH")],
        figures=figures,
        warnings=warn,
        notes=["方法：标签置换经验 p（SPARK-X 思想，Sun 2021；内置实现），不依赖解析分布。"],
        data_out=[data_file("permutation_csv", frame_path, "置换检验表")],
    )
    payload["store"] = {"svg_permutation": {"genes": [str(g) for g in frame["gene"].head(200)],
                                            "permutations": int(perms_eff)}}
    return payload


# ===========================================================================
# 4. cross_moran — 双变量 Moran's I
# ===========================================================================

def run_cross_moran(ctx: dict, params: dict, outdir: str) -> dict:
    """双变量空间自相关：双变量 Moran's I（Wartenberg 1985；内置实现）。

    ``I_ab = (1/S0)·Σ_i z_a,i·(Wz_b)_i``（z 为标准化表达，W 行标准化邻接），
    即"a 的取值与 b 的邻域均值"的相关；同时给出 I_ba 与标签置换 p。
    """
    ensure_dir(outdir)
    params = dict(params or {})
    view = _view(ctx)
    warn = list(view["warnings"])
    gene_a = _p_str(params.get("gene_a"))
    gene_b = _p_str(params.get("gene_b"))
    if not gene_a or not gene_b:
        raise prerequisite("双变量 Moran's I 需要同时指定 gene_a 与 gene_b")
    ia = _require_gene(view, gene_a)
    ib = _require_gene(view, gene_b)
    if ia == ib:
        raise bad_request("gene_a 与 gene_b 不能是同一个基因（单变量 Moran's I 请用 morans_svg）")
    perms = max(0, min(5000, _p_int(params.get("permutations"), 199)))
    W, radius, winfo = _spatial_weights(view, params)
    if winfo.get("fallback"):
        warn.append("spot 数较多，空间权重改用 cKDTree 等价半径图（与 sp_common.spatial_weights 同规则）")
    stats = _weight_stats(W, view["n"])
    if stats["S0"] <= 0:
        raise prerequisite("空间权重矩阵全为 0（spot 坐标过密或邻接半径过小），无法计算双变量 Moran's I")
    n = view["n"]
    S0 = float(stats["S0"])
    X = np.asarray(view["X"], dtype=np.float64)
    a = X[:, ia]
    b = X[:, ib]
    sd_a = float(np.std(a))
    sd_b = float(np.std(b))
    if sd_a <= EPS or sd_b <= EPS:
        raise prerequisite("gene_a 或 gene_b 的表达为常数（无方差），无法计算双变量 Moran's I")
    za = (a - float(np.mean(a))) / sd_a
    zb = (b - float(np.mean(b))) / sd_b
    Wzb = np.asarray(W @ zb).ravel()
    Wza = np.asarray(W @ za).ravel()
    moran_ab = float(np.sum(za * Wzb) / S0)
    moran_ba = float(np.sum(zb * Wza) / S0)
    pearson_r = float(np.mean(za * zb))
    p_perm = float("nan")
    if perms > 0:
        rng = np.random.default_rng(DEFAULT_SEED)
        null = np.empty(perms, dtype=float)
        for t in range(perms):
            perm = rng.permutation(n)
            Zp = zb[perm]
            Zp = (Zp - Zp.mean()) / (Zp.std() if Zp.std() > 0 else 1.0)
            null[t] = float(np.sum(za * np.asarray(W @ Zp).ravel()) / S0)
        ge = int(np.sum(np.abs(null) >= abs(moran_ab)))
        p_perm = float(_clip_p(np.array([(ge + 1.0) / (perms + 1.0)]))[0])
        perm_lo, perm_hi = float(np.min(null)), float(np.max(null))
    else:
        null = np.zeros(0)
        perm_lo, perm_hi = float("nan"), float("nan")

    frame = pd.DataFrame([{
        "gene_a": view["genes"][ia],
        "gene_b": view["genes"][ib],
        "moran_ab": moran_ab,
        "moran_ba": moran_ba,
        "pearson_r": pearson_r,
        "p_perm": p_perm,
        "permutations": int(perms),
        "n_spots": int(n),
        "neighbor_radius": float(radius),
    }])
    frame_path = export_csv(frame, outdir, "cross_moran", index=False)

    figures: list[dict] = []
    try:
        figures.append(_render_cross_moran_map(
            view, ia, ib, Wzb, (view["genes"][ia], view["genes"][ib]),
            {"za": za, "moran_ab": moran_ab}, outdir, params))
    except Exception as error:  # noqa: BLE001
        warn.append(f"双基因空间图绘制失败：{error}")
    if perms > 0 and null.size:
        try:
            figures.append(_render_permutation(null, np.array([moran_ab]), f"{view['genes'][ia]} x {view['genes'][ib]}",
                                               outdir, stat_name="Bivariate Moran's I"))
        except Exception as error:  # noqa: BLE001
            warn.append(f"置换零分布图绘制失败：{error}")

    direction = "正" if moran_ab > 0 else "负"
    bullets = [
        f"双变量 Moran's I（{view['genes'][ia]} × {view['genes'][ib]}）：I_ab={moran_ab:.4f}，"
        f"I_ba={moran_ba:.4f}，表达相关 r={pearson_r:.4f}",
        f"方向：{direction}相关" + (f"；标签置换 {perms} 次双侧 p={p_perm:.4f}" if perms > 0 else ""),
        f"含义：{view['genes'][ia]} 高表达的 spot 周围 {view['genes'][ib]} 也偏高（或相反）",
    ]
    metrics = {
        "moran_ab": moran_ab,
        "moran_ba": moran_ba,
        "pearson_r": pearson_r,
        "p_perm": p_perm,
        "n_spots": int(n),
        "neighbor_radius": float(radius),
        "permutations": int(perms),
        "perm_range": [perm_lo, perm_hi],
    }
    payload = result(
        title="双变量空间自相关",
        bullets=bullets,
        metrics=metrics,
        tables=[table("cross_moran", "双变量 Moran's I",
                      [{"key": "gene_a", "label": "基因 A", "type": "string"},
                       {"key": "gene_b", "label": "基因 B", "type": "string"},
                       {"key": "moran_ab", "label": "I(A,B)", "type": "number"},
                       {"key": "moran_ba", "label": "I(B,A)", "type": "number"},
                       {"key": "pearson_r", "label": "Pearson r", "type": "number"},
                       {"key": "p_perm", "label": "置换 p", "type": "number"},
                       {"key": "n_spots", "label": "spot 数", "type": "number"}],
                      frame, note="I_ab = 标准化 a 与 b 的邻域均值的平均乘积（Wartenberg 1985）")],
        figures=figures,
        warnings=warn,
        notes=[f"方法：双变量 Moran's I（Wartenberg 1985），行标准化邻接权重，内置实现；半径 {radius:.3g}。"],
        data_out=[data_file("cross_moran_csv", frame_path, "双变量 Moran 表")],
    )
    payload["store"] = {"cross_moran": {"gene_a": view["genes"][ia], "gene_b": view["genes"][ib],
                                        "moran_ab": moran_ab, "p_perm": p_perm}}
    return payload


# ===========================================================================
# 5. hotspot — Getis-Ord Gi*
# ===========================================================================

def _knn_indices(coords: np.ndarray, k: int, warn: list[str]) -> np.ndarray:
    """kNN 索引（不含自身）。小样本用 sp_common.coord_knn，大样本用 cKDTree 等价实现。"""
    import sp_common

    n = coords.shape[0]
    k = int(max(1, min(k, n - 1)))
    if n <= DENSE_LIMIT or cKDTree is None:
        idx, _ = sp_common.coord_knn(coords, k)
        return np.asarray(idx, dtype=int)
    warn.append("spot 数较多，Gi* 邻域改用 cKDTree 的等价 kNN 实现（结果与 coord_knn 一致）")
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k + 1)
    idx = np.atleast_2d(np.asarray(idx, dtype=int))
    return idx[:, 1:]


def _gi_star(values: np.ndarray, knn_idx: np.ndarray) -> np.ndarray:
    """Getis-Ord Gi*（含自身权重）：返回每个 spot 的 z 值。"""
    n = values.shape[0]
    k = knn_idx.shape[1]
    x = np.asarray(values, dtype=float)
    total_w = float(k + 1)
    total_w2 = float(k + 1)
    xbar = float(np.mean(x))
    var = float(np.mean(x * x) - xbar * xbar)
    if var <= EPS:
        return np.full(n, np.nan)
    local_sum = x + np.sum(x[knn_idx], axis=1)
    denom = math.sqrt(var) * math.sqrt(max((n * total_w2 - total_w * total_w) / max(n - 1, 1), 0.0))
    if denom <= EPS:
        return np.full(n, np.nan)
    return (local_sum - xbar * total_w) / denom


def _hotspot_core(ctx: dict, params: dict, *, allow_auto: bool = True) -> dict:
    """热点计算内核（分析 / PLOTS 共用）。"""
    view = _view(ctx)
    warn = list(view["warnings"])
    gene = _p_str(params.get("gene"))
    if not gene:
        alt = _p_list(params.get("genes"))
        gene = alt[0] if alt else ""
    k = max(1, min(50, _p_int(params.get("k"), 8)))
    z_thresh = abs(_p_float(params.get("z_thresh"), 1.96))
    if not gene:
        if not allow_auto:
            raise prerequisite("热点检测需要指定基因参数（gene）")
        gi, hint = _auto_gene(view, why="热点检测")
        gene = view["genes"][gi]
        if hint:
            warn.append(hint)
    gi = _require_gene(view, gene)
    values = np.asarray(view["X"][:, gi], dtype=float)
    if float(np.var(values)) <= EPS:
        raise prerequisite(f"基因「{gene}」的表达为常数（无方差），无法计算 Getis-Ord Gi*")
    knn = _knn_indices(view["coords"], k, warn)
    z = _gi_star(values, knn)
    if not np.isfinite(z).any():
        raise prerequisite(f"基因「{gene}」的 Gi* 统计量全为缺失（表达方差为 0），无法检测热点")
    cats = _hotspot_categories(z, z_thresh)
    valid = np.isfinite(z)
    return {"view": view, "warn": warn, "gene": gene, "k": k, "z_thresh": z_thresh,
            "values": values, "z": z, "cats": cats, "valid": valid, "knn": knn}


def run_hotspot(ctx: dict, params: dict, outdir: str) -> dict:
    """热点检测：Getis-Ord Gi*（Getis & Ord 1992；内置实现）。

    用 ``coord_knn`` 的 k 邻域 + 自身构造二值权重，按
    ``z_i = (Σ_j w_ij x_j − x̄ Σ_j w_ij) / (S·sqrt((nΣw² − (Σw)²)/(n−1)))``
    计算每个 spot 的 z；``z > z_thresh`` 为热点、``z < −z_thresh`` 为冷点。
    """
    ensure_dir(outdir)
    params = dict(params or {})
    core = _hotspot_core(ctx, params)
    view, warn = core["view"], list(core["warn"])
    gene, z, cats = core["gene"], core["z"], core["cats"]
    valid = core["valid"]
    n_hot = int(np.sum(cats == "hot spot"))
    n_cold = int(np.sum(cats == "cold spot"))
    n_ns = int(np.sum(cats == "not significant"))

    names = _spot_names(view)
    spot_frame = pd.DataFrame({
        "spot": names,
        "x": view["coords"][:, 0],
        "y": view["coords"][:, 1],
        "expression": core["values"],
        "gi_star_z": z,
        "class": cats,
    })
    spot_frame = spot_frame.sort_values("gi_star_z", ascending=False).reset_index(drop=True)
    summary_frame = pd.DataFrame([
        {"class": "hot spot", "n_spots": n_hot, "percent": 100.0 * n_hot / max(view["n"], 1)},
        {"class": "cold spot", "n_spots": n_cold, "percent": 100.0 * n_cold / max(view["n"], 1)},
        {"class": "not significant", "n_spots": n_ns, "percent": 100.0 * n_ns / max(view["n"], 1)},
    ])
    summary_frame["mean_gi_z"] = [
        float(np.nanmean(z[cats == lab])) if np.any(cats == lab) and np.isfinite(z[cats == lab]).any() else np.nan
        for lab in summary_frame["class"]
    ]
    summary_frame["mean_expression"] = [
        float(np.nanmean(core["values"][cats == lab])) if np.any(cats == lab) else np.nan
        for lab in summary_frame["class"]
    ]
    summary_path = export_csv(summary_frame, outdir, "hotspot_summary", index=False)
    spots_path = export_csv(spot_frame, outdir, "hotspot_spots", index=False)

    figures: list[dict] = []
    try:
        figures.append(_render_hotspot_map(view, cats, core["values"], gene, outdir, params))
    except Exception as error:  # noqa: BLE001
        warn.append(f"热点图绘制失败：{error}")
    try:
        fig, ax = new_fig(7.4, 3.4)
        ax.hist(z[valid], bins=40, color="#8fb4d9", edgecolor="white", linewidth=0.4)
        ax.axvline(core["z_thresh"], color=HOT_COLOR, linewidth=1.2, linestyle="--", label="hot threshold")
        ax.axvline(-core["z_thresh"], color=COLD_COLOR, linewidth=1.2, linestyle="--", label="cold threshold")
        ax.set_xlabel("Getis-Ord Gi* z-score")
        ax.set_ylabel("Number of spots")
        ax.set_title(f"Gi* distribution: {gene}")
        ax.legend(fontsize=8)
        figures.append(save_fig(fig, outdir, "hotspot_z_hist", "Gi* z-score distribution"))
    except Exception as error:  # noqa: BLE001
        warn.append(f"Gi* 分布图绘制失败：{error}")

    bullets = [
        f"基因「{gene}」（k={core['k']} 邻域 + 自身）：热点 {n_hot} 个 spot、冷点 {n_cold} 个"
        f"（共 {view['n']} 个，阈值 |z|>{core['z_thresh']:.2f}）",
        f"Gi* z 范围 {float(np.nanmin(z)):.2f} ~ {float(np.nanmax(z)):.2f}；"
        f"热点平均表达 {float(np.nanmean(core['values'][cats == 'hot spot'])):.3g}"
        if n_hot else f"Gi* z 范围 {float(np.nanmin(z)):.2f} ~ {float(np.nanmax(z)):.2f}；未检出热点",
    ]
    metrics = {
        "gene": gene,
        "n_hot": n_hot,
        "n_cold": n_cold,
        "n_not_sig": n_ns,
        "k": int(core["k"]),
        "z_thresh": float(core["z_thresh"]),
        "max_z": float(np.nanmax(z)),
        "min_z": float(np.nanmin(z)),
        "n_spots": int(view["n"]),
    }
    payload = result(
        title="热点检测（Getis-Ord Gi*）",
        bullets=bullets,
        metrics=metrics,
        tables=[
            table("hotspot_summary", "热点/冷点汇总",
                  [{"key": "class", "label": "类别", "type": "string"},
                   {"key": "n_spots", "label": "spot 数", "type": "number"},
                   {"key": "percent", "label": "占比 %", "type": "number"},
                   {"key": "mean_gi_z", "label": "平均 Gi* z", "type": "number"},
                   {"key": "mean_expression", "label": "平均表达", "type": "number"}],
                  summary_frame),
            table("hotspot_spots", "Top Gi* spot",
                  [{"key": "spot", "label": "spot", "type": "string"},
                   {"key": "x", "label": "x", "type": "number"},
                   {"key": "y", "label": "y", "type": "number"},
                   {"key": "expression", "label": "表达", "type": "number"},
                   {"key": "gi_star_z", "label": "Gi* z", "type": "number"},
                   {"key": "class", "label": "类别", "type": "string"}],
                  spot_frame.head(200), note="完整 spot 级结果见 hotspot_spots.csv"),
        ],
        figures=figures,
        warnings=warn,
        notes=["方法：Getis-Ord Gi*（Getis & Ord 1992），kNN 邻域二值权重（含自身），内置实现。"],
        data_out=[data_file("hotspot_csv", summary_path, "热点汇总"),
                  data_file("hotspot_csv", spots_path, "spot 级 Gi* 结果")],
    )
    payload["store"] = {"hotspot": {"gene": gene, "n_hot": n_hot, "n_cold": n_cold}}
    return payload


# ===========================================================================
# 6. spatial_gradient — 沿轴梯度
# ===========================================================================

def _gradient_core(ctx: dict, params: dict, *, allow_auto: bool = True) -> dict:
    """空间梯度内核（分析 / PLOTS 共用）。"""
    view = _view(ctx)
    warn = list(view["warnings"])
    gene = _p_str(params.get("gene"))
    if not gene:
        alt = _p_list(params.get("genes"))
        gene = alt[0] if alt else ""
    if not gene:
        if not allow_auto:
            raise prerequisite("空间梯度需要指定基因参数（gene）")
        gi, hint = _auto_gene(view, why="空间梯度")
        gene = view["genes"][gi]
        if hint:
            warn.append(hint)
    gi = _require_gene(view, gene)
    axis = _p_str(params.get("axis"), "xy").lower()
    angle = _p_float(params.get("angle"), 0.0)
    n_bins = max(2, min(200, _p_int(params.get("bins"), 20)))
    coords = view["coords"]
    x = coords[:, 0].astype(float)
    y = coords[:, 1].astype(float)
    if axis in ("x", "xxx", "horizontal"):
        t = x
        axis_label = "x coordinate"
    elif axis in ("y", "vertical"):
        t = y
        axis_label = "y coordinate"
    elif axis in ("xy", "projection", "angle", "diagonal"):
        theta = math.radians(angle)
        t = x * math.cos(theta) + y * math.sin(theta)
        axis_label = f"projection (angle={angle:g} deg)"
    else:
        raise bad_request(f"axis 参数「{axis}」不支持（可选 x / y / xy）")
    trange = float(np.nanmax(t) - np.nanmin(t))
    if not math.isfinite(trange) or trange <= EPS:
        raise prerequisite("沿该轴的投影坐标没有变化（所有 spot 在同一位置），无法计算空间梯度")
    values = np.asarray(view["X"][:, gi], dtype=float)
    if float(np.var(values)) <= EPS:
        raise prerequisite(f"基因「{gene}」的表达为常数（无方差），无法计算空间梯度")

    # 等频分箱（按投影坐标排序后均分，避免空箱）
    order = np.argsort(t, kind="stable")
    ranks = np.empty(len(t), dtype=int)
    ranks[order] = np.arange(len(t))
    bins = np.minimum((ranks * n_bins) // len(t), n_bins - 1)
    rows = []
    for b in range(n_bins):
        sel = bins == b
        if not sel.any():
            continue
        vals = values[sel]
        sd = float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0
        se = sd / math.sqrt(vals.size) if vals.size > 0 else np.nan
        rows.append({
            "bin": int(b + 1),
            "n_spots": int(sel.sum()),
            "t_center": float(np.mean(t[sel])),
            "mean": float(np.mean(vals)),
            "sd": sd,
            "se": float(se),
        })
    bin_frame = pd.DataFrame(rows)
    if bin_frame.empty:
        raise prerequisite("分箱后没有有效 spot，无法计算空间梯度")

    from scipy import stats as _stats

    lin = _stats.linregress(t, values)
    slope = float(lin.slope)
    intercept = float(lin.intercept)
    pval = float(_clip_p(np.array([lin.pvalue]))[0])
    r_value = float(lin.rvalue)
    bin_frame["fitted"] = intercept + slope * bin_frame["t_center"]
    return {
        "view": view,
        "warn": warn,
        "gene": gene,
        "axis": axis,
        "axis_label": axis_label,
        "angle": float(angle),
        "bins": bin_frame,
        "slope": slope,
        "intercept": intercept,
        "p": pval,
        "r": r_value,
        "r2": float(r_value * r_value),
        "t": t,
        "values": values,
        "n_bins": int(len(bin_frame)),
    }


def run_spatial_gradient(ctx: dict, params: dict, outdir: str) -> dict:
    """空间梯度：沿指定轴投影 + 等频分箱均值±SE + 线性趋势（内置实现）。

    组织梯度分析的常用做法：把 spot 投影到指定方向（``axis=x/y/xy``，``angle``
    为 xy 投影角度），按等频分箱统计表达均值与标准误，并用
    ``scipy.stats.linregress`` 给斜率与双侧 p。
    """
    ensure_dir(outdir)
    params = dict(params or {})
    core = _gradient_core(ctx, params)
    view, warn = core["view"], list(core["warn"])
    frame = core["bins"]
    frame_path = export_csv(frame, outdir, "spatial_gradient", index=False)

    figures: list[dict] = []
    try:
        figures.append(_render_gradient_plot(frame, core["gene"], core["axis_label"],
                                            core["slope"], core["p"], outdir))
    except Exception as error:  # noqa: BLE001
        warn.append(f"梯度图绘制失败：{error}")
    try:
        fig, ax = new_fig(7.4, 3.4)
        order = np.argsort(core["t"])
        ax.scatter(core["t"][order], core["values"][order], s=max(3.0, _point_size(view, params) * 0.5),
                   c="#c9ced6", alpha=0.6, linewidths=0, label="spots")
        ax.plot(frame["t_center"], frame["mean"], "o-", color="#2e6da4", markersize=4,
                linewidth=1.4, label="binned mean")
        xs = np.linspace(float(np.min(core["t"])), float(np.max(core["t"])), 50)
        ax.plot(xs, core["intercept"] + core["slope"] * xs, color="#c0392b", linestyle="--",
                linewidth=1.2, label="linear fit")
        ax.set_xlabel(core["axis_label"])
        ax.set_ylabel("expression")
        ax.set_title(f"Gradient of {core['gene']}")
        ax.legend(fontsize=8)
        figures.append(save_fig(fig, outdir, "gradient_spots", "Spot-level gradient"))
    except Exception as error:  # noqa: BLE001
        warn.append(f"spot 级梯度散点图绘制失败：{error}")

    direction = "升高" if core["slope"] > 0 else "降低"
    bullets = [
        f"基因「{core['gene']}」沿 {core['axis_label']} 的线性趋势：slope={core['slope']:.4g}/单位，"
        f"r={core['r']:.3f}，p={core['p']:.4g}",
        f"方向：沿该轴表达呈{direction}趋势" if core["p"] < 0.05 else "方向：线性趋势不显著",
        f"等频分箱 {core['n_bins']} 个（每箱均值 ± SE 已写入 CSV）",
    ]
    metrics = {
        "gene": core["gene"],
        "axis": core["axis"],
        "angle": core["angle"],
        "slope": core["slope"],
        "r": core["r"],
        "r2": core["r2"],
        "p": core["p"],
        "n_bins": core["n_bins"],
        "n_spots": int(view["n"]),
    }
    payload = result(
        title="空间梯度",
        bullets=bullets,
        metrics=metrics,
        tables=[table("gradient", "分箱表达统计",
                      [{"key": "bin", "label": "分箱", "type": "number"},
                       {"key": "n_spots", "label": "spot 数", "type": "number"},
                       {"key": "t_center", "label": "投影坐标中心", "type": "number"},
                       {"key": "mean", "label": "均值", "type": "number"},
                       {"key": "sd", "label": "标准差", "type": "number"},
                       {"key": "se", "label": "标准误", "type": "number"},
                       {"key": "fitted", "label": "线性拟合值", "type": "number"}],
                      frame, note="按投影坐标等频分箱；斜率/ p 来自全部 spot 的线性回归")],
        figures=figures,
        warnings=warn,
        notes=["方法：沿轴投影 + 等频分箱均值±SE + 线性趋势（内置实现；组织梯度常做法）。"],
        data_out=[data_file("gradient_csv", frame_path, "空间梯度分箱表")],
    )
    payload["store"] = {"spatial_gradient": {"gene": core["gene"], "slope": core["slope"], "p": core["p"]}}
    return payload


# ===========================================================================
# 7. spatial_corr_pairs — 空间共表达基因对
# ===========================================================================

def _corr_pairs_core(ctx: dict, params: dict) -> dict:
    """基因对相关 + 邻域乘积一致性内核。"""
    view = _view(ctx)
    warn = list(view["warnings"])
    min_expr = max(0.0, _p_float(params.get("min_expr"), 0.1))
    top_pairs = max(1, min(500, _p_int(params.get("top_pairs"), 30)))
    max_genes = max(2, min(400, _p_int(params.get("max_genes"), 60)))
    requested = _p_list(params.get("genes"))
    if requested:
        idx, names, pick_warn = _resolve_gene_set(view, params, n_top=max_genes, min_expr=min_expr,
                                                  default_n=max_genes)
        if pick_warn:
            warn.append(pick_warn)
        if len(idx) > max_genes:
            warn.append(f"指定了 {len(idx)} 个基因，超过 max_genes={max_genes}，只取前 {max_genes} 个")
            idx, names = idx[:max_genes], names[:max_genes]
    else:
        idx, warn_pick = _pick_variance_genes(view, max_genes, min_expr)
        if warn_pick:
            warn.append(warn_pick)
        if not idx:
            raise prerequisite("没有可用基因（全部为常数或平均表达过低），无法计算空间共表达")
        names = [view["genes"][i] for i in idx]
    if len(idx) < 2:
        raise prerequisite("至少需要 2 个非恒定基因才能计算基因对相关")
    W, radius, winfo = _spatial_weights(view, params)
    if winfo.get("fallback"):
        warn.append("spot 数较多，空间权重改用 cKDTree 等价半径图（与 sp_common.spatial_weights 同规则）")
    if float(_weight_stats(W, view["n"])["S0"]) <= 0:
        raise prerequisite("空间权重矩阵全为 0（spot 坐标过密或邻接半径过小），无法计算邻域一致性")

    X = _submatrix(view, idx)
    mean = X.mean(axis=0)
    sd = X.std(axis=0)
    sd = np.where(sd > EPS, sd, 1.0)
    Z = (X - mean) / sd
    n = view["n"]
    g = Z.shape[1]
    # 全局 Pearson：z 标准化后即为均值内积
    corr = (Z.T @ Z) / float(n)
    np.fill_diagonal(corr, 1.0)

    pairs_i, pairs_j = [], []
    chunk = 256
    pair_list = [(i, j) for i in range(g) for j in range(i + 1, g)]
    if not pair_list:
        raise prerequisite("基因数不足，无法构造基因对")
    ii = np.array([p[0] for p in pair_list], dtype=int)
    jj = np.array([p[1] for p in pair_list], dtype=int)
    coh = np.empty(len(pair_list), dtype=float)
    mprod = np.empty(len(pair_list), dtype=float)
    for start in range(0, len(pair_list), chunk):
        stop = min(start + chunk, len(pair_list))
        prod = Z[:, ii[start:stop]] * Z[:, jj[start:stop]]  # (n, chunk)
        lag = np.asarray(W @ prod)                          # (n, chunk)
        pm = prod.mean(axis=0)
        lm = lag.mean(axis=0)
        psd = prod.std(axis=0)
        lsd = lag.std(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            cc = np.where((psd > EPS) & (lsd > EPS),
                          ((prod * lag).mean(axis=0) - pm * lm) / (psd * lsd), np.nan)
        coh[start:stop] = cc
        mprod[start:stop] = pm
    r_values = corr[ii, jj]
    # Pearson 相关检验 p（t 分布）
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(np.abs(r_values) < 1.0,
                         r_values * np.sqrt((n - 2) / np.maximum(1e-12, 1.0 - r_values * r_values)), np.inf)
    from scipy import stats as _stats

    pvals = _clip_p(2.0 * _stats.t.sf(np.abs(tstat), max(n - 2, 1)))
    frame = pd.DataFrame({
        "gene_a": [names[i] for i in ii],
        "gene_b": [names[j] for j in jj],
        "pearson": r_values,
        "p_pearson": pvals,
        "neighbor_coherence": coh,
        "mean_product": mprod,
        "n_spots": int(n),
    })
    frame["p_adj"] = bh_fdr(frame["p_pearson"])
    frame = frame.sort_values(["neighbor_coherence", "pearson"], ascending=[False, False],
                              na_position="last").reset_index(drop=True)
    return {"view": view, "warn": warn, "frame": frame, "top_pairs": top_pairs,
            "n_genes": len(idx), "radius": radius, "genes": names}


def run_spatial_corr_pairs(ctx: dict, params: dict, outdir: str) -> dict:
    """空间共表达基因对：Pearson 相关 + 邻域乘积一致性（内置实现）。

    * ``pearson``：两基因在 spot 间的全局 Pearson 相关（BH FDR）；
    * ``neighbor_coherence``：spot 上两基因标准化表达的乘积与其邻域平均乘积的相关，
      衡量"相邻 spot 上两基因同时高表达"是否一致（空间共表达）。
    """
    ensure_dir(outdir)
    params = dict(params or {})
    core = _corr_pairs_core(ctx, params)
    view, warn = core["view"], list(core["warn"])
    frame = core["frame"]
    full_path = export_csv(frame, outdir, "spatial_corr_pairs", index=False)
    top = frame.head(core["top_pairs"])

    figures: list[dict] = []
    try:
        figures.append(_render_corr_pairs(frame, core["genes"], outdir))
    except Exception as error:  # noqa: BLE001
        warn.append(f"共表达网络图绘制失败：{error}")

    strong = int((frame["neighbor_coherence"].fillna(0.0) > 0.5).sum())
    bullets = [
        f"对 {core['n_genes']} 个基因的 {len(frame)} 个基因对计算全局相关与邻域乘积一致性"
        f"（半径 {core['radius']:.3g}，{view['n']} 个 spot）",
        f"邻域一致性 > 0.5 的基因对：{strong} 个；最强一对为"
        f"「{top.iloc[0]['gene_a']}–{top.iloc[0]['gene_b']}」（coherence={top.iloc[0]['neighbor_coherence']:.3f}）"
        if len(top) else "没有可报告的基因对",
        "全局相关高但邻域一致性低 = 表达相关但空间不同步；两者都高提示空间共表达模块",
    ]
    metrics = {
        "n_genes": int(core["n_genes"]),
        "n_pairs": int(len(frame)),
        "n_coherence_gt_0.5": strong,
        "top_pair": f"{top.iloc[0]['gene_a']}–{top.iloc[0]['gene_b']}" if len(top) else None,
        "top_coherence": float(top.iloc[0]["neighbor_coherence"]) if len(top) else None,
        "neighbor_radius": float(core["radius"]),
        "n_spots": int(view["n"]),
    }
    payload = result(
        title="空间共表达基因对",
        bullets=bullets,
        metrics=metrics,
        tables=[table("corr_pairs", "空间共表达基因对（top）",
                      [{"key": "gene_a", "label": "基因 A", "type": "string"},
                       {"key": "gene_b", "label": "基因 B", "type": "string"},
                       {"key": "pearson", "label": "Pearson r", "type": "number"},
                       {"key": "p_pearson", "label": "P 值", "type": "number"},
                       {"key": "p_adj", "label": "FDR", "type": "number"},
                       {"key": "neighbor_coherence", "label": "邻域一致性", "type": "number"},
                       {"key": "mean_product", "label": "平均乘积", "type": "number"}],
                      top, note="按邻域一致性降序；完整基因对见 spatial_corr_pairs.csv")],
        figures=figures,
        warnings=warn,
        notes=["方法：Pearson 相关 + 邻域乘积一致性（内置实现）；邻域为行标准化半径邻接。"],
        data_out=[data_file("corr_pairs_csv", full_path, "基因对相关表")],
    )
    payload["store"] = {"spatial_corr_pairs": {"top_pairs": jsonable(top.head(50)),
                                               "n_genes": int(core["n_genes"])}}
    return payload


# ===========================================================================
# 8. ripley_k — Ripley's K / L
# ===========================================================================

def _cluster_labels(ctx: dict, key: str, spot_idx: np.ndarray) -> tuple[np.ndarray, str]:
    """取 spot 级聚类标签（与视图 spot 顺序对齐），返回 (标签数组, 列名)。"""
    key = _p_str(key) or _p_str(ctx.get("cluster_key")) or "clusters"
    info = ctx.get("sample_info")
    spots = [str(s) for s in (ctx.get("spots") or [])]
    labels = None
    if isinstance(info, pd.DataFrame) and key in info.columns:
        series = info[key]
        try:
            if spots and len(info) == len(spots):
                # 行顺序与 spots 一致（本插件的 sample_info 以 spot 为行）
                labels = series.astype(str).tolist()
            elif spots and len(set(map(str, info.index)) & set(spots)) > 0:
                labels = series.reindex(spots).astype(str).tolist()
            else:
                labels = series.astype(str).tolist()
        except Exception:  # noqa: BLE001
            labels = None
    if labels is None and ctx.get("cluster") is not None and (
        key in ("cluster", "clusters") or key == _p_str(ctx.get("cluster_key"))
    ):
        labels = [str(v) for v in ctx["cluster"]]
    if labels is None:
        alt = _p_str(ctx.get("cluster_key"))
        if alt and alt != key and isinstance(info, pd.DataFrame) and alt in info.columns:
            return _cluster_labels(ctx, alt, spot_idx)
        raise prerequisite(f"找不到聚类列「{key}」，请先运行 cluster 步骤或用 cluster_key 指定列名")

    arr = np.asarray(labels, dtype=object)
    n_total = len(np.asarray(ctx.get("coords")))
    if arr.size == n_total:
        arr = arr[np.asarray(spot_idx, dtype=int)]
    elif arr.size == len(spot_idx):
        arr = arr[np.arange(len(spot_idx))]
    else:
        raise prerequisite(f"聚类列「{key}」长度（{arr.size}）与坐标数（{n_total}）不匹配")
    arr = np.asarray([str(v) if v is not None else "NA" for v in arr], dtype=object)
    return arr, key


def _ripley_core(ctx: dict, params: dict) -> dict:
    """Ripley's K 内核：每簇 K(r) 与 CSR 随机点包络。"""
    view = _view(ctx, need_x=False)
    warn = list(view["warnings"])
    labels, key = _cluster_labels(ctx, _p_str(params.get("cluster_key"), "clusters"), view["spot_idx"])
    radii_n = max(4, min(200, _p_int(params.get("radii"), 20)))
    n_sim = max(0, min(999, _p_int(params.get("n_sim"), 19)))
    coords = view["coords"]
    n = view["n"]
    xspan = float(np.max(coords[:, 0]) - np.min(coords[:, 0]))
    yspan = float(np.max(coords[:, 1]) - np.min(coords[:, 1]))
    area = xspan * yspan
    if not math.isfinite(area) or area <= EPS:
        raise prerequisite("空间坐标在一个维度上没有跨度（退化成一条线），无法定义 Ripley's K 的窗口面积")
    lo = np.array([np.min(coords[:, 0]), np.min(coords[:, 1])])
    hi = np.array([np.max(coords[:, 0]), np.max(coords[:, 1])])
    uniq, counts = np.unique(labels.astype(str), return_counts=True)
    order = np.argsort(counts)[::-1]
    uniq, counts = uniq[order], counts[order]
    if uniq.size == 0:
        raise prerequisite("聚类标签为空，无法计算 Ripley's K")
    max_clusters = 12
    if uniq.size > max_clusters:
        warn.append(f"共 {uniq.size} 个簇，只对 spot 数最多的 {max_clusters} 个计算 Ripley's K")
        uniq, counts = uniq[:max_clusters], counts[:max_clusters]
    r_max = 0.5 * math.sqrt(area)
    r_min = r_max / max(8.0, radii_n / 2.0)
    radii = np.linspace(r_min, r_max, radii_n)
    rng = np.random.default_rng(DEFAULT_SEED)
    items: list[dict] = []
    skipped: list[str] = []
    max_points = 4000
    for label, count in zip(uniq, counts):
        pts = coords[labels.astype(str) == label]
        if pts.shape[0] > max_points:
            sel = rng.choice(pts.shape[0], size=max_points, replace=False)
            pts = pts[np.sort(sel)]
            warn.append(f"簇「{label}」spot 数 {int(count)} 过大，Ripley's K 内部下采样到 {max_points} 个")
        m = pts.shape[0]
        if m < 5:
            skipped.append(str(label))
            continue
        from scipy.spatial.distance import pdist

        d = np.sort(pdist(pts))
        counts_r = np.searchsorted(d, radii, side="right").astype(float)
        k_obs = area * counts_r / float(m * (m - 1))
        k_sims = np.empty((max(n_sim, 0), radii_n), dtype=float)
        for s in range(max(n_sim, 0)):
            sim = rng.uniform(lo, hi, size=(m, 2))
            ds = np.sort(pdist(sim))
            cs = np.searchsorted(ds, radii, side="right").astype(float)
            k_sims[s] = area * cs / float(m * (m - 1))
        if max(n_sim, 0) > 0:
            k_mean = k_sims.mean(axis=0)
            l_sims = np.sqrt(np.maximum(k_sims, 0.0) / math.pi) - radii[None, :]
            l_lo = np.percentile(l_sims, 2.5, axis=0)
            l_hi = np.percentile(l_sims, 97.5, axis=0)
        else:
            k_mean = np.zeros(radii_n)
            l_lo = l_hi = np.full(radii_n, np.nan)
        l_obs = np.sqrt(np.maximum(k_obs, 0.0) / math.pi) - radii
        items.append({
            "label": str(label),
            "n": int(m),
            "k": k_obs,
            "k_csr_mean": k_mean,
            "l_obs": l_obs,
            "l_csr_mean": np.sqrt(np.maximum(k_mean, 0.0) / math.pi) - radii,
            "l_dev": l_obs - (np.sqrt(np.maximum(k_mean, 0.0) / math.pi) - radii),
            "l_dev_lo": l_lo,
            "l_dev_hi": l_hi,
        })
    if not items:
        raise prerequisite("每个簇的 spot 数都少于 5，无法计算 Ripley's K（请先增加聚类或改用更大的簇）")
    if skipped:
        warn.append("以下簇 spot 数 < 5，已跳过：" + "、".join(skipped[:8]))
    return {"view": view, "warn": warn, "labels": labels, "key": key, "radii": radii,
            "area": area, "items": items, "r_min": float(radii[0]), "r_max": float(radii[-1]),
            "n_sim": int(n_sim)}


def run_ripley_k(ctx: dict, params: dict, outdir: str) -> dict:
    """Ripley's K/L 函数：Ripley 1977 + 完全空间随机（CSR）包络（内置实现）。

    每簇 ``K(r) = A/(m(m-1))·Σ_{i≠j} 1(d_ij ≤ r)``，``L(r)=sqrt(K(r)/π)-r``；
    包络来自 ``n_sim`` 次在观测窗口（外接矩形）内均匀随机撒点（CSR）的同规模模拟，
    观测与模拟用同一窗口，因此边界偏差可比。
    """
    ensure_dir(outdir)
    params = dict(params or {})
    core = _ripley_core(ctx, params)
    view, warn, items = core["view"], list(core["warn"]), core["items"]
    radii = core["radii"]
    rows = []
    for item in items:
        for t, radius in enumerate(radii):
            rows.append({
                "cluster": item["label"],
                "n_spots": int(item["n"]),
                "r": float(radius),
                "k_obs": float(item["k"][t]),
                "k_csr_mean": float(item["k_csr_mean"][t]),
                "l_obs": float(item["l_obs"][t]),
                "l_dev": float(item["l_dev"][t]),
                "l_dev_lo": float(item["l_dev_lo"][t]) if np.isfinite(item["l_dev_lo"][t]) else np.nan,
                "l_dev_hi": float(item["l_dev_hi"][t]) if np.isfinite(item["l_dev_hi"][t]) else np.nan,
            })
    frame = pd.DataFrame(rows)
    frame_path = export_csv(frame, outdir, "ripley_k", index=False)

    figures: list[dict] = []
    try:
        figures.append(_render_ripley({"radii": radii, "clusters": items}, outdir))
    except Exception as error:  # noqa: BLE001
        warn.append(f"Ripley 曲线绘制失败：{error}")

    dev_table = frame.groupby("cluster", as_index=False).agg(
        n_spots=("n_spots", "max"), r_max=("r", "max"),
        l_dev_max=("l_dev", "max"), l_dev_min=("l_dev", "min"))
    dev_table = dev_table.sort_values("l_dev_max", ascending=False).reset_index(drop=True)
    clustered = dev_table[dev_table["l_dev_max"] > 0]["cluster"].astype(str).tolist()
    bullets = [
        f"对 {len(items)} 个簇计算 Ripley's K/L（窗口面积 {core['area']:.3g}，"
        f"r ∈ [{core['r_min']:.3g}, {core['r_max']:.3g}]，{core['n_sim']} 次 CSR 模拟）",
        "L(r)-r 超出 CSR 包络且 > 0 的簇（空间聚集）：" + ("、".join(clustered[:8]) if clustered else "无"),
        "说明：窗口取 spot 外接矩形，观测与 CSR 用同一窗口，边界偏差可比；K 值偏小不代表无聚集",
    ]
    metrics = {
        "n_clusters": int(len(items)),
        "cluster_key": core["key"],
        "r_min": core["r_min"],
        "r_max": core["r_max"],
        "n_sim": core["n_sim"],
        "window_area": float(core["area"]),
        "clustered_clusters": clustered,
    }
    payload = result(
        title="Ripley's K/L 函数",
        bullets=bullets,
        metrics=metrics,
        tables=[
            table("ripley_curves", "Ripley's K/L 曲线（逐簇逐半径）",
                  [{"key": "cluster", "label": "簇", "type": "string"},
                   {"key": "n_spots", "label": "spot 数", "type": "number"},
                   {"key": "r", "label": "r", "type": "number"},
                   {"key": "k_obs", "label": "K(r)", "type": "number"},
                   {"key": "k_csr_mean", "label": "CSR 平均 K", "type": "number"},
                   {"key": "l_obs", "label": "L(r)", "type": "number"},
                   {"key": "l_dev", "label": "L(r)-r 偏差", "type": "number"},
                   {"key": "l_dev_lo", "label": "包络下界", "type": "number"},
                   {"key": "l_dev_hi", "label": "包络上界", "type": "number"}],
                  frame, note="l_dev = 观测 L(r)-r − CSR 平均 L(r)-r；包络为 CSR 2.5%~97.5%"),
            table("ripley_summary", "逐簇偏差汇总",
                  [{"key": "cluster", "label": "簇", "type": "string"},
                   {"key": "n_spots", "label": "spot 数", "type": "number"},
                   {"key": "r_max", "label": "最大 r", "type": "number"},
                   {"key": "l_dev_max", "label": "L 偏差最大值", "type": "number"},
                   {"key": "l_dev_min", "label": "L 偏差最小值", "type": "number"}],
                  dev_table),
        ],
        figures=figures,
        warnings=warn,
        notes=["方法：Ripley's K/L（Ripley 1977）+ CSR 随机点包络，内置实现。"],
        data_out=[data_file("ripley_csv", frame_path, "Ripley K/L 表")],
    )
    payload["store"] = {"ripley_k": {"clusters": [str(i["label"]) for i in items],
                                     "r_max": core["r_max"], "n_sim": core["n_sim"]}}
    return payload


# ===========================================================================
# 9. gene_autocorr_scan — 全基因向量化扫描
# ===========================================================================

def _col_mean_var(X, warn: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """分块计算列均值/方差（float64 累加，避免整表复制大矩阵）。"""
    n, g = X.shape
    rows = max(1, int(4_000_000 // max(g, 1)))
    s = np.zeros(g, dtype=np.float64)
    s2 = np.zeros(g, dtype=np.float64)
    for start in range(0, n, rows):
        block = np.asarray(X[start:start + rows], dtype=np.float64)
        s += block.sum(axis=0)
        s2 += np.einsum("ij,ij->j", block, block)
    mean = s / max(n, 1)
    var = np.maximum(s2 / max(n, 1) - mean * mean, 0.0)
    return mean, var


def run_gene_autocorr_scan(ctx: dict, params: dict, outdir: str) -> dict:
    """全基因空间自相关扫描：向量化 Moran's I（内置实现）。

    对每个合格基因用稀疏矩阵乘法一次性计算 ``Wz`` 得到 Moran's I，配正态近似
    z / p 与 BH FDR；``max_genes`` 限制扫描基因数（超过则按方差取前 N 个）。
    """
    ensure_dir(outdir)
    params = dict(params or {})
    view = _view(ctx)
    warn = list(view["warnings"])
    max_genes = max(1, min(60000, _p_int(params.get("max_genes"), 2000)))
    min_expr = max(0.0, _p_float(params.get("min_expr"), 0.1))
    W, radius, winfo = _spatial_weights(view, params)
    if winfo.get("fallback"):
        warn.append("spot 数较多，空间权重改用 cKDTree 等价半径图（与 sp_common.spatial_weights 同规则）")
    stats = _weight_stats(W, view["n"])
    if stats["S0"] <= 0:
        raise prerequisite("空间权重矩阵全为 0（spot 坐标过密或邻接半径过小），无法扫描空间自相关")
    if stats["n_isolated"] > 0:
        warn.append(f"{stats['n_isolated']} 个 spot 没有邻居（权重行和为 0），只参与方差、不参与协方差")

    X_all = np.asarray(view["X"])
    n, total_genes = X_all.shape
    mean, var = _col_mean_var(X_all, warn)
    usable = np.isfinite(var) & (var > EPS)
    expr_ok = usable & (mean >= float(min_expr))
    cand = np.where(expr_ok)[0]
    relaxed = False
    if cand.size == 0:
        cand = np.where(usable)[0]
        relaxed = True
        if cand.size:
            warn.append(f"没有基因的平均表达 ≥ {min_expr:g}，已放宽到全部非恒定基因")
    if cand.size == 0:
        raise prerequisite("没有非恒定基因（表达矩阵全部为常数），无法扫描空间自相关")
    if cand.size > max_genes:
        order = cand[np.argsort(var[cand])[::-1]][:max_genes]
        warn.append(f"合格基因 {int(cand.size)} 个，超过 max_genes={max_genes}，按方差取前 {max_genes} 个扫描")
        cand = np.sort(order)
    if relaxed:
        warn.append("本次扫描按放宽后的基因集合进行，min_expr 实际未生效")

    n_scan = int(cand.size)
    S0 = float(stats["S0"])
    expected = float(stats["expected"])
    variance = float(stats["variance"])
    I = np.full(n_scan, np.nan)
    chunk = max(8, int(8_000_000 // max(n, 1)))
    for start in range(0, n_scan, chunk):
        stop = min(start + chunk, n_scan)
        block = np.ascontiguousarray(X_all[:, cand[start:stop]], dtype=np.float64)
        I[start:stop] = _moran_vector(block, W, n, S0)
    if math.isfinite(variance) and variance > 0:
        with np.errstate(invalid="ignore"):
            z = (I - expected) / math.sqrt(variance)
        pvals = 2.0 * _norm_sf(np.abs(z))
    else:
        warn.append("正态近似方差不可用（权重结构退化），z/p 记为缺失")
        z = np.full(n_scan, np.nan)
        pvals = np.ones(n_scan)
    pvals = np.where(np.isfinite(I), _clip_p(pvals), np.nan)
    frame = pd.DataFrame({
        "gene": [view["genes"][int(j)] for j in cand],
        "moran_i": I,
        "expected": expected,
        "z": z,
        "p": pvals,
        "mean": mean[cand],
        "var": var[cand],
        "n_spots": int(n),
    })
    frame["p_adj"] = bh_fdr(frame["p"])
    frame = frame.sort_values("moran_i", ascending=False, na_position="last").reset_index(drop=True)
    frame_path = export_csv(frame, outdir, "gene_autocorr_scan", index=False)

    figures: list[dict] = []
    try:
        figures.append(_render_autocorr_scan(frame, expected, outdir))
    except Exception as error:  # noqa: BLE001
        warn.append(f"扫描分布图绘制失败：{error}")

    n_sig = int(((frame["p_adj"].fillna(1.0) < 0.05) & (frame["moran_i"] > 0)).sum())
    top = frame.dropna(subset=["moran_i"])
    bullets = [
        f"扫描 {n_scan} 个基因（数据集合计 {total_genes} 个，min_expr={min_expr:g}），"
        f"{n} 个 spot，行标准化半径邻接（半径 {radius:.3g}）",
        f"FDR<0.05 且 I>0 的显著空间高变基因：{n_sig} 个"
        + (f"，最强「{top.iloc[0]['gene']}」（I={top.iloc[0]['moran_i']:.3f}）" if len(top) else ""),
        f"I 中位数 {float(frame['moran_i'].median()):.4f}，E[I]={expected:.4f}",
    ]
    metrics = {
        "n_genes_scanned": n_scan,
        "n_genes_total": int(total_genes),
        "n_significant": n_sig,
        "median_moran_i": float(frame["moran_i"].median()),
        "top_gene": str(top.iloc[0]["gene"]) if len(top) else None,
        "top_moran_i": float(top.iloc[0]["moran_i"]) if len(top) else None,
        "neighbor_radius": float(radius),
        "n_spots": int(n),
    }
    payload = result(
        title="全基因空间自相关扫描",
        bullets=bullets,
        metrics=metrics,
        tables=[table("autocorr_scan", "全基因 Moran's I 扫描",
                      [{"key": "gene", "label": "基因", "type": "string"},
                       {"key": "moran_i", "label": "Moran's I", "type": "number"},
                       {"key": "expected", "label": "E[I]", "type": "number"},
                       {"key": "z", "label": "z 值", "type": "number"},
                       {"key": "p", "label": "P 值", "type": "number"},
                       {"key": "p_adj", "label": "FDR", "type": "number"},
                       {"key": "mean", "label": "平均表达", "type": "number"},
                       {"key": "var", "label": "方差", "type": "number"}],
                      frame, note=f"按 Moran's I 降序；JSON 只显示前 500 个，完整结果见 gene_autocorr_scan.csv")],
        figures=figures,
        warnings=warn,
        notes=["方法：全基因向量化 Moran's I（稀疏权重矩阵乘法）+ 正态近似 + BH FDR，内置实现。"],
        data_out=[data_file("autocorr_csv", frame_path, "全基因扫描表")],
    )
    payload["store"] = {"gene_autocorr_scan": {"top_genes": [str(g) for g in frame["gene"].head(500)],
                                               "n_genes_scanned": n_scan,
                                               "neighbor_radius": float(radius)}}
    return payload


# ===========================================================================
# PLOTS：可复用绘图入口（registry 合并成 kind → fn）
# ===========================================================================

def _plot_gene_list(ctx: dict, params: dict, default_n: int) -> tuple[list[int], list[str], dict]:
    """PLOTS 里解析基因参数（支持 gene / genes / 复用 morans_svg 的结果）。"""
    view = _view(ctx)
    warn = list(view["warnings"])
    requested = _p_list(params.get("genes")) or _p_list(params.get("gene"))
    if not requested:
        stored = _stored_svg_genes(ctx)
        if stored:
            requested = stored[:default_n]
            warn.append(f"未指定基因，复用最近一次 morans_svg 的 top {len(requested)} 基因")
    if requested:
        idx, names, pick_warn = _resolve_gene_set({**view, "warnings": warn}, {"genes": requested},
                                                  n_top=default_n, min_expr=0.0, default_n=default_n)
        if pick_warn:
            warn.append(pick_warn)
        return idx, names, {"view": view, "warn": warn}
    idx, names, pick_warn = _resolve_gene_set(view, {}, n_top=default_n, min_expr=0.0, default_n=default_n)
    if pick_warn:
        warn.append(pick_warn)
    warn.append(f"未指定基因，按表达方差取前 {len(idx)} 个基因绘图")
    return idx, names, {"view": view, "warn": warn}


def _stored_svg_genes(ctx: dict) -> list[str]:
    """从 ctx["results"] 读取最近一次 SVG 分析的基因列表（容错两种 store 形态）。"""
    results = (ctx or {}).get("results") or {}
    if not isinstance(results, dict):
        return []
    for key in ("morans_svg", "gene_autocorr_scan", "svg_permutation", "gearys_svg"):
        store = results.get(key)
        if not isinstance(store, dict):
            continue
        nested = store.get(key) if isinstance(store.get(key), dict) else store
        genes = nested.get("genes") or nested.get("top_genes")
        if isinstance(genes, (list, tuple)) and genes:
            return [str(g) for g in genes]
    return []


def plot_svg_map(ctx: dict, params: dict, outdir: str) -> dict:
    """``svg_map``：Top 空间高变基因表达图（默认复用 morans_svg 结果或按方差取前 4 个）。"""
    ensure_dir(outdir)
    params = dict(params or {})
    n_show = max(1, min(8, _p_int(params.get("n_top"), 4)))
    idx, names, extra = _plot_gene_list(ctx, params, n_show)
    return _render_svg_map(extra["view"], idx[:n_show], names[:n_show], outdir, params)


def plot_moran_plot(ctx: dict, params: dict, outdir: str) -> dict:
    """``moran_plot``：Moran's I 排序条形图（无结果时现场计算候选基因）。"""
    ensure_dir(outdir)
    params = dict(params or {})
    view = _view(ctx)
    warn = list(view["warnings"])
    n_top = max(1, min(200, _p_int(params.get("n_top"), 30)))
    idx, names, pick_warn = _resolve_gene_set(view, params, n_top=n_top, min_expr=0.0, default_n=n_top)
    if pick_warn:
        warn.append(pick_warn)
    W, radius, _info = _spatial_weights(view, params)
    stats = _weight_stats(W, view["n"])
    frame, _null = _moran_frame(view, idx, names, W, radius, stats, 0, np.random.default_rng(DEFAULT_SEED), warn)
    return _render_moran_plot(frame, outdir)


def plot_ripley(ctx: dict, params: dict, outdir: str) -> dict:
    """``ripley``：Ripley's K/L 曲线（默认 20 个半径档、19 次 CSR）。"""
    ensure_dir(outdir)
    params = dict(params or {})
    core = _ripley_core(ctx, params)
    return _render_ripley({"radii": core["radii"], "clusters": core["items"]}, outdir)


def plot_hotspot_map(ctx: dict, params: dict, outdir: str) -> dict:
    """``hotspot_map``：Getis-Ord Gi* 热点图（未指定基因时自动取方差最大基因）。"""
    ensure_dir(outdir)
    params = dict(params or {})
    core = _hotspot_core(ctx, params)
    return _render_hotspot_map(core["view"], core["cats"], core["values"], core["gene"], outdir, params)


def plot_gradient_plot(ctx: dict, params: dict, outdir: str) -> dict:
    """``gradient_plot``：空间梯度分箱曲线 + 线性趋势。"""
    ensure_dir(outdir)
    params = dict(params or {})
    core = _gradient_core(ctx, params)
    return _render_gradient_plot(core["bins"], core["gene"], core["axis_label"],
                                 core["slope"], core["p"], outdir)


# ===========================================================================
# 注册表
# ===========================================================================

_SVG_COMMON_PARAMS = [
    {"key": "genes", "label": "指定基因", "type": "string", "default": "",
     "hint": "逗号或空格分隔；留空 = 按表达方差自动挑选"},
    {"key": "n_top", "label": "候选基因数", "type": "number", "default": 20, "min": 1, "max": 2000, "step": 1,
     "hint": "未指定基因时按方差取前 n 个"},
    {"key": "min_expr", "label": "最低平均表达", "type": "number", "default": 0.1, "min": 0.0,
     "hint": "自动挑基因时过滤平均表达过低的基因"},
    {"key": "neighbor_radius", "label": "邻接半径", "type": "number", "default": None,
     "hint": "留空 = 用数据集估计的 neighbor_radius"},
]

_HOTSPOT_PARAMS = [
    {"key": "gene", "label": "基因", "type": "string", "default": "",
     "hint": "留空 = 自动选择方差最大的基因"},
    {"key": "k", "label": "kNN 邻居数", "type": "number", "default": 8, "min": 1, "max": 50, "step": 1},
    {"key": "z_thresh", "label": "热点 z 阈值", "type": "number", "default": 1.96, "min": 0.5, "max": 10,
     "hint": "|z| 超过该值判为热点/冷点（1.96 ≈ p<0.05）"},
]

_GRADIENT_PARAMS = [
    {"key": "gene", "label": "基因", "type": "string", "default": "",
     "hint": "留空 = 自动选择方差最大的基因"},
    {"key": "axis", "label": "投影轴", "type": "select", "default": "xy",
     "options": [{"value": "x", "label": "x 轴"}, {"value": "y", "label": "y 轴"},
                 {"value": "xy", "label": "任意角度投影"}]},
    {"key": "angle", "label": "投影角度", "type": "number", "default": 0, "min": -360, "max": 360,
     "hint": "axis=xy 时生效（度）；0=沿 x 轴"},
    {"key": "bins", "label": "分箱数", "type": "number", "default": 20, "min": 2, "max": 200, "step": 1},
]


def run_spatial_stat(ctx: dict, params: dict, outdir: str) -> dict:
    """兼容入口：按 ``params["analysis"]`` 分发（registry 一般不使用）。"""
    mapping = {item["id"]: item["run"] for item in ANALYSES}
    aid = _p_str((params or {}).get("analysis"))
    if aid not in mapping:
        raise bad_request(f"未知的空间统计分分析：{aid}")
    return mapping[aid](ctx, params, outdir)


ANALYSES = [
    {
        "id": "morans_svg",
        "title": "空间自相关（Moran's I）",
        "group": ANALYSIS_GROUP,
        "method": "Moran's I（Moran 1950）+ 正态近似与标签置换；行标准化半径邻接权重；内置实现",
        "needs": ["cluster"],
        "params": _SVG_COMMON_PARAMS + [
            {"key": "permutations", "label": "置换次数", "type": "number", "default": 100, "min": 0, "max": 5000,
             "step": 1, "hint": "0 = 只做解析 p"},
        ],
        "run": run_morans_svg,
    },
    {
        "id": "gearys_svg",
        "title": "空间自相关（Geary's C）",
        "group": ANALYSIS_GROUP,
        "method": "Geary's C（Geary 1954）+ 标签置换经验 p；行标准化半径邻接权重；内置实现",
        "needs": ["cluster"],
        "params": _SVG_COMMON_PARAMS + [
            {"key": "permutations", "label": "置换次数", "type": "number", "default": 99, "min": 0, "max": 5000,
             "step": 1},
        ],
        "run": run_gearys_svg,
    },
    {
        "id": "svg_permutation",
        "title": "空间置换检验",
        "group": ANALYSIS_GROUP,
        "method": "标签置换经验 p=(k+1)/(n+1)（SPARK-X 思想，Sun 2021）；内置实现",
        "needs": ["cluster"],
        "params": _SVG_COMMON_PARAMS + [
            {"key": "permutations", "label": "置换次数", "type": "number", "default": 200, "min": 1, "max": 5000,
             "step": 1},
        ],
        "run": run_svg_permutation,
    },
    {
        "id": "cross_moran",
        "title": "双变量空间自相关",
        "group": ANALYSIS_GROUP,
        "method": "双变量 Moran's I（Wartenberg 1985）+ 标签置换；内置实现",
        "needs": ["cluster"],
        "params": [
            {"key": "gene_a", "label": "基因 A", "type": "string", "default": ""},
            {"key": "gene_b", "label": "基因 B", "type": "string", "default": ""},
            {"key": "permutations", "label": "置换次数", "type": "number", "default": 199, "min": 0, "max": 5000,
             "step": 1},
            {"key": "neighbor_radius", "label": "邻接半径", "type": "number", "default": None},
        ],
        "run": run_cross_moran,
    },
    {
        "id": "hotspot",
        "title": "热点检测（Getis-Ord Gi*）",
        "group": ANALYSIS_GROUP,
        "method": "Getis-Ord Gi*（Getis & Ord 1992）；kNN 邻域二值权重（含自身）；内置实现",
        "needs": ["cluster"],
        "params": _HOTSPOT_PARAMS,
        "run": run_hotspot,
    },
    {
        "id": "spatial_gradient",
        "title": "空间梯度",
        "group": ANALYSIS_GROUP,
        "method": "沿轴投影等频分箱 + 均值±SE + 线性趋势检验（内置实现；组织梯度常用做法）",
        "needs": ["cluster"],
        "params": _GRADIENT_PARAMS,
        "run": run_spatial_gradient,
    },
    {
        "id": "spatial_corr_pairs",
        "title": "空间共表达基因对",
        "group": ANALYSIS_GROUP,
        "method": "Pearson 相关 + 邻域乘积一致性（内置实现；行标准化半径邻接）",
        "needs": ["cluster"],
        "params": [
            {"key": "genes", "label": "指定基因", "type": "string", "default": "",
             "hint": "逗号或空格分隔；留空 = 按方差取前 max_genes 个"},
            {"key": "top_pairs", "label": "返回基因对数", "type": "number", "default": 30, "min": 1, "max": 500},
            {"key": "max_genes", "label": "最多基因数", "type": "number", "default": 60, "min": 2, "max": 400},
            {"key": "min_expr", "label": "最低平均表达", "type": "number", "default": 0.1, "min": 0.0},
            {"key": "neighbor_radius", "label": "邻接半径", "type": "number", "default": None},
        ],
        "run": run_spatial_corr_pairs,
    },
    {
        "id": "ripley_k",
        "title": "Ripley's K/L 函数",
        "group": ANALYSIS_GROUP,
        "method": "Ripley's K/L（Ripley 1977）+ 完全空间随机(CSR) 随机点包络；内置实现",
        "needs": ["cluster"],
        "params": [
            {"key": "cluster_key", "label": "聚类列", "type": "string", "default": "clusters"},
            {"key": "radii", "label": "半径档数", "type": "number", "default": 20, "min": 4, "max": 200, "step": 1},
            {"key": "n_sim", "label": "CSR 模拟次数", "type": "number", "default": 19, "min": 0, "max": 999, "step": 1},
        ],
        "run": run_ripley_k,
    },
    {
        "id": "gene_autocorr_scan",
        "title": "全基因空间自相关扫描",
        "group": ANALYSIS_GROUP,
        "method": "全基因向量化 Moran's I（稀疏矩阵乘法）+ 正态近似 + BH FDR；内置实现",
        "needs": ["normalize"],
        "params": [
            {"key": "max_genes", "label": "最多扫描基因数", "type": "number", "default": 2000,
             "min": 1, "max": 60000, "step": 1, "hint": "超过则按表达方差取前 N 个"},
            {"key": "min_expr", "label": "最低平均表达", "type": "number", "default": 0.1, "min": 0.0},
            {"key": "neighbor_radius", "label": "邻接半径", "type": "number", "default": None},
        ],
        "run": run_gene_autocorr_scan,
    },
]

STEPS: list[dict] = []

PLOTS = {
    "svg_map": plot_svg_map,
    "moran_plot": plot_moran_plot,
    "ripley": plot_ripley,
    "hotspot_map": plot_hotspot_map,
    "gradient_plot": plot_gradient_plot,
}
