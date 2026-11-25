#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-spatial — 空间域与邻域模块（sp-adv）。

步骤
----
* ``svg``     空间高变基因：Moran's I（Moran 1950）/ Geary's C（Geary 1954），行标准化权重 + 正态近似或置换经验 p。
* ``domains`` 空间域：PCA 前 15 维 + 空间邻接图融合 → k-means 初始化 → 马尔可夫随机场（MRF）迭代平滑
  （BayesSpace 2021 思想，内置实现）；``method='graph'`` 时在融合图上跑 Louvain（networkx）。

分析（group=邻域与域）
--------------------
nhood_enrichment / co_occurrence / domains_mrf / domain_markers / domain_composition /
boundary_score / distance_to_domain / niche_cluster / tissue_diversity

分析（group=区域差异）
--------------------
de_domain / de_regions / de_groups / gene_domain_box

图（PLOTS）
----------
svg_rank_plot / domain_map / nhood_enrich / cooccurrence / domain_composition_bar /
boundary_plot / distance_plot / niche_map / diversity_map / de_volcano / de_heatmap /
gene_domain_box / domain_dotplot

约定
----
* 只依赖标准库 + numpy/pandas/scipy/sklearn/networkx/matplotlib + ``sp_common``（不 import 其它 sp_* 模块，避免循环）。
* 所有方法为**内置实现**，不调用 squidpy / BayesSpace / CellChat 等官方包；method 列如实标注依据的文献/产品思想。
* 缺前置 → ``prerequisite`` 中文；参数非法 → ``bad_request``；不吞异常返回假结果。
* 确定性：``np.random.default_rng(0)``；图内文字一律英文；matplotlib 只用 Agg（sp_common 已设定）。
"""

from __future__ import annotations

import math
import os
import tempfile

import numpy as np
import pandas as pd

import sp_common as sc

__all__ = [
    "STEPS",
    "ANALYSES",
    "PLOTS",
    "run_svg",
    "run_domains",
    "compute_domains",
    "get_domain_labels",
]

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

METHOD_SVG = ("moran", "geary")
METHOD_DOMAIN = ("mrf", "graph")

# 空间域列名候选（写回/读取时按顺序尝试）
DOMAIN_COLUMNS = ("domains", "domain", "spatial_domain", "domain_id")
GROUP_COLUMNS = ("slice", "group", "batch", "sample", "region")

MAX_PLOT_POINTS = 20000
MAX_MARKER_GENES_FOR_FIG = 60


# ---------------------------------------------------------------------------
# 参数归一化（UI 可能把数字/布尔传成字符串）
# ---------------------------------------------------------------------------

def _p_str(value, default: str = "") -> str:
    if value is None:
        return default
    text = str(value).strip()
    return text if text else default


def _p_int(value, default: int) -> int:
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return int(default)


def _p_float(value, default: float) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    return float(number) if math.isfinite(number) else float(default)


def _p_bool(value, default: bool = False) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in ("1", "true", "yes", "y", "on", "是"):
        return True
    if text in ("0", "false", "no", "n", "off", "否", ""):
        return False
    return bool(default)


def _p_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set, np.ndarray, pd.Index)):
        return [str(v).strip() for v in list(value) if str(v).strip()]
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in text.replace(";", ",").replace("|", ",").split(",") if item.strip()]


def _outdir(ctx, outdir) -> str:
    path = outdir or (ctx or {}).get("outdir")
    if not path:
        path = tempfile.mkdtemp(prefix="sp-domain-")
    return sc.ensure_dir(str(path))


def _rng(ctx=None) -> np.random.Generator:
    rng = (ctx or {}).get("rng")
    if isinstance(rng, np.random.Generator):
        return rng
    return np.random.default_rng(0)


# ---------------------------------------------------------------------------
# ctx 读取
# ---------------------------------------------------------------------------

def _matrix(ctx) -> tuple[np.ndarray, list[str], list[str]]:
    """返回 (X[n_spots, n_genes], genes, spots)。"""
    X = ctx.get("X")
    if X is None:
        raise sc.prerequisite("当前数据集没有表达矩阵，请先导入数据")
    X = np.asarray(X)
    if X.ndim != 2 or X.size == 0:
        raise sc.prerequisite("当前表达矩阵为空，无法做空间域/邻域分析（请先导入数据并完成标准化）")
    genes = [str(g) for g in (ctx.get("genes") or [])]
    spots = [str(s) for s in (ctx.get("spots") or [])]
    if X.shape[1] != len(genes) or X.shape[0] != len(spots):
        raise sc.EngineError(
            "internal",
            f"表达矩阵维度（{X.shape[0]}×{X.shape[1]}）与基因/spot 数（{len(spots)}×{len(genes)}）不一致",
        )
    return X, genes, spots


def _coords(ctx) -> np.ndarray:
    coords = ctx.get("coords")
    if coords is None:
        raise sc.prerequisite("该数据集没有空间坐标，无法做空间域/邻域分析（请导入带坐标的空间数据）")
    arr = np.asarray(coords, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2 or arr.shape[0] == 0:
        raise sc.prerequisite("空间坐标格式不正确（需要 n_spots × 2 的 x/y 坐标）")
    n = len(ctx.get("spots") or [])
    if n and arr.shape[0] != n:
        raise sc.EngineError("internal", f"坐标行数（{arr.shape[0]}）与 spot 数（{n}）不一致")
    return arr[:, :2]


def _labels_from_store(obj) -> list[str] | None:
    """从 dataset.results 里的对象里取出一列标签（兼容多种存储形状）。"""
    if obj is None:
        return None
    if isinstance(obj, dict):
        for key in ("labels", "values", "domain", "domains", "clusters"):
            value = obj.get(key)
            if isinstance(value, (list, tuple, np.ndarray)) and len(value):
                return [str(v) for v in list(value)]
        table_obj = obj.get("table")
        if isinstance(table_obj, pd.DataFrame):
            for col in ("domain", "label", "domain_id", "labels"):
                if col in table_obj.columns:
                    return [str(v) for v in table_obj[col].tolist()]
        if isinstance(table_obj, list) and table_obj and isinstance(table_obj[0], dict):
            for col in ("domain", "label", "domain_id", "labels"):
                if col in table_obj[0]:
                    return [str(row.get(col)) for row in table_obj]
    if isinstance(obj, pd.DataFrame):
        for col in ("domain", "label", "domain_id", "labels"):
            if col in obj.columns:
                return [str(v) for v in obj[col].tolist()]
    if isinstance(obj, (list, tuple, np.ndarray)) and len(obj):
        first = obj[0]
        if isinstance(first, (str, int, float, np.integer, np.floating)):
            return [str(v) for v in list(obj)]
        if isinstance(first, dict):
            for col in ("domain", "label", "domain_id", "labels"):
                if col in first:
                    return [str(row.get(col)) for row in obj]
    return None


def _column_labels(ctx, key: str, n: int) -> list[str] | None:
    """从 spot 级 sample_info 里取一列（按 ctx["spots"] 顺序对齐）。"""
    info = ctx.get("sample_info")
    if not isinstance(info, pd.DataFrame) or info.empty or not key or key not in info.columns:
        return None
    spots = [str(s) for s in (ctx.get("spots") or [])]
    if len(info.index) != n:
        return None
    try:
        if spots and all(s in info.index for s in spots):
            series = info.loc[spots, key]
        else:
            series = info[key]
    except Exception:  # noqa: BLE001
        return None
    if isinstance(series, pd.DataFrame):
        series = series.iloc[:, 0]
    if len(series) != n:
        return None
    out = []
    for value in series.tolist():
        if value is None or (isinstance(value, float) and math.isnan(value)) or pd.isna(value):
            out.append("NA")
        else:
            out.append(str(value).strip())
    if all(v in ("NA", "") for v in out):
        return None
    return out


def find_labels(ctx, keys, n: int | None = None) -> tuple[list[str] | None, str]:
    """按候选列名找 spot 级标签，返回 (labels, 命中的列名)。找不到返回 (None, "")。"""
    ctx = ctx or {}
    if n is None:
        n = len(ctx.get("spots") or [])
    results = ctx.get("results") if isinstance(ctx.get("results"), dict) else {}
    for key in keys:
        key = _p_str(key)
        if not key:
            continue
        labels = _column_labels(ctx, key, n)
        if labels is not None and len(labels) == n:
            return labels, key
        direct = ctx.get(key)
        if direct is not None:
            try:
                arr = np.asarray(direct, dtype=object).ravel()
            except Exception:  # noqa: BLE001
                arr = None
            if arr is not None and arr.size == n:
                return [str(v) for v in arr.tolist()], key
        labels = _labels_from_store(results.get(key))
        if labels is not None and len(labels) == n:
            return labels, key
        for alias in ("domains", "domain"):
            if key in ("clusters", "cluster", "cluster_key"):
                continue
            labels = _labels_from_store(results.get(alias))
            if labels is not None and len(labels) == n and alias == key:
                return labels, alias
    return None, ""


def _cluster_labels(ctx, key: str = "", n: int | None = None) -> tuple[list[str], str]:
    """读取聚类标签（缺前置 → prerequisite）。"""
    if n is None:
        n = len(ctx.get("spots") or [])
    key = _p_str(key) or _p_str(ctx.get("cluster_key")) or "clusters"
    labels, hit = find_labels(ctx, [key, "clusters", "cluster", "louvain", "leiden"])
    if labels is None and key in ("clusters", "cluster", "cluster_key"):
        direct = ctx.get("cluster")
        if direct is not None:
            try:
                arr = np.asarray(direct, dtype=object).ravel()
            except Exception:  # noqa: BLE001
                arr = None
            if arr is not None and arr.size == n:
                return [str(v) for v in arr.tolist()], _p_str(ctx.get("cluster_key"), "clusters")
    if labels is None:
        raise sc.prerequisite("还没有聚类结果，请先运行「聚类」(cluster) 步骤")
    return labels, hit or key


def _domain_labels(ctx, n: int | None = None) -> tuple[list[str] | None, str]:
    """读取空间域标签；没有就返回 (None, "")。"""
    if n is None:
        n = len(ctx.get("spots") or [])
    return find_labels(ctx, list(DOMAIN_COLUMNS), n)


def _match_label(value, labels: list[str]) -> str | None:
    """把参数里的 domain 值匹配到实际标签（支持 "1" / 1 / 1.0 互认）。"""
    text = _p_str(value)
    if not text:
        return None
    if text in labels:
        return text
    low = text.lower()
    for label in labels:
        if label.lower() == low:
            return label
    try:
        number = float(text)
    except ValueError:
        return None
    for label in labels:
        try:
            if float(label) == number:
                return label
        except ValueError:
            continue
    return None


def _require_label(value, labels: list[str], what: str = "域") -> str:
    hit = _match_label(value, labels)
    if hit is not None:
        return hit
    if not _p_str(value):
        uniq = list(dict.fromkeys(labels))
        raise sc.prerequisite(f"请指定{what}（可选：{'、'.join(uniq[:12])}）")
    uniq = list(dict.fromkeys(labels))
    raise sc.prerequisite(f"{what}「{value}」不存在（可选：{'、'.join(uniq[:12])}）")


def _group_values(ctx, n: int, group_col: str = "") -> tuple[list[str] | None, str]:
    """spot 级分组标签（用于 domain_composition / de_groups / deconv_roi）。"""
    key = _p_str(group_col)
    candidates = [key] if key else []
    ctx_group = _p_str(ctx.get("group_key"))
    if ctx_group:
        candidates.append(ctx_group)
    candidates.extend(GROUP_COLUMNS)
    for cand in candidates:
        if not cand:
            continue
        labels = _column_labels(ctx, cand, n)
        if labels is not None and len(set(labels)) >= 1:
            return labels, cand
    return None, ""


# ---------------------------------------------------------------------------
# 空间邻接 / 邻居
# ---------------------------------------------------------------------------

def _blocks(ctx) -> tuple[list[np.ndarray], list[str], str | None, list[str]]:
    """物理块（切片/样本）标签：``(blocks, names, key, warnings)``。

    ``sp_common.spatial_blocks`` 优先取 slide/library/slice/sample/specimen/patient 列；
    没有这类列时返回单个全量块（单切片数据，行为与旧版一致）。
    """
    return sc.spatial_blocks(ctx)


def _blocked(ctx) -> bool:
    """是否多切片/多样本（物理块 > 1）。"""
    try:
        blocks, _names, _key, _warns = _blocks(ctx)
    except Exception:  # noqa: BLE001
        return False
    return len(blocks) > 1


def _knn(coords: np.ndarray, k: int, blocks=None) -> tuple[np.ndarray, np.ndarray]:
    """按坐标的 k 近邻（不含自身），返回 (indices[n,k], distances[n,k])。用 KD-tree 控制内存。

    多切片（``blocks`` > 1）时邻居只在块内取（``sp_common.blocked_coord_knn``），
    绝不把「不同切片的相似坐标」当邻居。
    """
    from scipy.spatial import cKDTree

    arr = np.asarray(coords, dtype=float)
    n = arr.shape[0]
    if n < 2:
        return np.zeros((n, 0), dtype=int), np.zeros((n, 0))
    if blocks is not None and len(blocks) > 1:
        idx, dist, _k_eff = sc.blocked_coord_knn(arr, blocks, int(k))
        return np.asarray(idx, dtype=int), np.asarray(dist, dtype=float)
    k = int(max(1, min(int(k), n - 1)))
    tree = cKDTree(arr)
    dist, idx = tree.query(arr, k=k + 1)
    dist = np.atleast_2d(dist)
    idx = np.atleast_2d(idx)
    if idx.shape[1] < k + 1:
        k = idx.shape[1] - 1
    return idx[:, 1:].astype(int), dist[:, 1:]


def _sparse_row_norm(idx: np.ndarray, n: int, data=None):
    """把邻居索引构造成行标准化的稀疏平均矩阵。"""
    from scipy import sparse

    if idx.size == 0:
        return sparse.csr_matrix((n, n))
    rows = np.repeat(np.arange(n), idx.shape[1])
    cols = idx.ravel().astype(int)
    vals = np.ones(rows.size, dtype=float) if data is None else np.asarray(data, dtype=float).ravel()
    mat = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
    # 对称化（邻接关系天然对称，kNN 会因为边界不对称）
    mat = mat.maximum(mat.T)
    row_sum = np.asarray(mat.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    return sparse.diags(1.0 / row_sum) @ mat


def _spatial_graph(ctx, k: int = 6, blocks=None):
    """空间平滑用的行标准化稀疏邻接矩阵（对称 kNN，多切片时块内取邻居）。"""
    coords = _coords(ctx)
    if blocks is None:
        blocks = _blocks(ctx)[0]
    idx, _ = _knn(coords, k, blocks=blocks)
    return _sparse_row_norm(idx, coords.shape[0])


def _radius_weights(ctx):
    """Moran's I / Geary's C 用的行标准化空间权重。

    多切片（物理块 > 1）：按切片分块建边（``sp_common.blocked_spatial_weights``），半径按块估计
    （``blocked_radius_estimate``）—— 合并坐标里「同坐标不同切片」会把最近邻距离带偏成 0/fallback。
    单块：n ≤ 12000 时用 sp_common.spatial_weights（半径图，与 import 估计的 neighbor_radius 一致）；
    更大规模用 KD-tree 半径查询构造等价稀疏矩阵，避免 n×n 稠密距离矩阵爆内存。
    """
    from scipy import sparse

    coords = _coords(ctx)
    n = coords.shape[0]
    blocks, _names, _key, _warns = _blocks(ctx)
    if len(blocks) > 1:
        radius = float(sc.blocked_radius_estimate(coords, blocks))
        return sc.blocked_spatial_weights(coords, blocks, radius)
    radius = ctx.get("neighbor_radius")
    try:
        radius = float(radius) if radius is not None else 0.0
    except (TypeError, ValueError):
        radius = 0.0
    if n <= 12000:
        return sc.spatial_weights(coords, radius if radius > 0 else None)
    from scipy.spatial import cKDTree

    r = radius if radius > 0 else max(float(sc.estimate_radius(coords)), 1e-9)
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r, output_type="ndarray")
    if pairs.size == 0:
        # 半径太小：退化为 k 近邻
        k = min(6, n - 1)
        idx, _ = _knn(coords, k)
        return _sparse_row_norm(idx, n)
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    mat = sparse.csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n, n))
    row_sum = np.asarray(mat.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    return sparse.diags(1.0 / row_sum) @ mat


# ---------------------------------------------------------------------------
# 空间自相关统计（Moran's I / Geary's C，内置实现）
# ---------------------------------------------------------------------------

def _weight_moments(W) -> tuple[float, float, float]:
    """S0/S1/S2（Cliff & Ord 1981 的权重矩）。"""
    from scipy import sparse

    w = sparse.coo_matrix(W)
    data = w.data.astype(float)
    s0 = float(data.sum())
    sym = (W + W.T).tocsr()
    s1 = float(0.5 * (sym.multiply(sym)).sum())
    row = np.asarray(W.sum(axis=1)).ravel()
    col = np.asarray(W.sum(axis=0)).ravel()
    s2 = float(((row + col) ** 2).sum())
    return s0, s1, s2


def _moran_geary_block(block: np.ndarray, W, method: str, row_sum: np.ndarray):
    """一个基因块的统计量。block 为 (n, g) 已中心化；返回 (stat, den)。"""
    Wb = W @ block
    num = np.asarray((block * Wb).sum(axis=0)).ravel()
    den = np.asarray((block * block).sum(axis=0)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        if method == "moran":
            stat = num / den
        else:
            Wb2 = np.asarray((W @ (block * block)).sum(axis=0)).ravel()
            sq = np.asarray((block * block).sum(axis=0)).ravel()
            num_c = (row_sum * sq) + Wb2 - 2.0 * num
            stat = num_c / den
    stat = np.where(np.isfinite(stat), stat, np.nan)
    return stat


def _svg_scores(X: np.ndarray, W, method: str, chunk: int = 256,
                blocks=None) -> tuple[np.ndarray, np.ndarray, dict]:
    """全基因 Moran's I / Geary's C 与正态近似的 z 值。

    多切片：每个基因先做块内中心化再做空间自相关（去掉切片间均值差，避免伪信号）。
    """
    n, m = X.shape
    s0, s1, s2 = _weight_moments(W)
    row_sum = np.asarray(W.sum(axis=1)).ravel()
    stats = np.full(m, np.nan)
    zvals = np.full(m, np.nan)
    moments: dict = {"S0": s0, "S1": s1, "S2": s2}
    if n < 3 or s0 <= 0:
        return stats, zvals, moments
    if method == "moran":
        expect = -1.0 / (n - 1)
        var = (n * n * s1 - n * s2 + 3.0 * s0 * s0) / ((n * n - 1) * s0 * s0) - expect * expect
    else:
        expect = 1.0
        var = ((2.0 * s1 + s2) * (n - 1) - 4.0 * s0 * s0) / (2.0 * (n + 1) * s0 * s0)
    moments["expect"] = float(expect)
    moments["var"] = float(var) if np.isfinite(var) and var > 0 else float("nan")
    scale = (n - 1) / (2.0 * s0) if method == "geary" else 1.0
    for start in range(0, m, chunk):
        stop = min(m, start + chunk)
        block = np.asarray(X[:, start:stop], dtype=np.float64)
        if blocks is not None and len(blocks) > 1:
            block = sc.center_within_blocks(block, blocks)
        else:
            mean = block.mean(axis=0, keepdims=True)
            block = block - mean
        stat = _moran_geary_block(block, W, method, row_sum)
        stats[start:stop] = stat * scale
    if np.isfinite(moments["var"]) and moments["var"] > 0:
        zvals = (stats - expect) / math.sqrt(moments["var"])
    return stats, zvals, moments


def _normal_p(zvals: np.ndarray) -> np.ndarray:
    """双侧正态 p 值（scipy 不可用时用 erfc 回退）。"""
    arr = np.asarray(zvals, dtype=float)
    try:
        from scipy import stats

        p = 2.0 * stats.norm.sf(np.abs(arr))
    except Exception:  # noqa: BLE001
        p = np.array([math.erfc(abs(v) / math.sqrt(2.0)) for v in np.nan_to_num(arr, nan=0.0)])
    p = np.where(np.isfinite(p), p, 1.0)
    return np.clip(p, 0.0, 1.0)


def _permutation_p(X, W, gene_idx: np.ndarray, method: str, permutations: int, rng,
                   blocks=None) -> np.ndarray:
    """标签置换经验 p 值（Moran 单侧 greater；Geary 单侧 less）。

    多切片：只在块内打乱标签（``sp_common.within_block_permutation``），零分布不交换切片；
    并把观测/置换值都做块内中心化。
    """
    n = X.shape[0]
    idx = np.asarray(gene_idx, dtype=int)
    if idx.size == 0 or permutations <= 0:
        return np.full(idx.size, np.nan)
    permutations = int(max(1, min(permutations, 2000)))
    multi = blocks is not None and len(blocks) > 1
    block0 = np.asarray(X[:, idx], dtype=np.float64)
    if multi:
        block0 = sc.center_within_blocks(block0, blocks)
    else:
        block0 = block0 - block0.mean(axis=0, keepdims=True)
    s0, _, _ = _weight_moments(W)
    row_sum = np.asarray(W.sum(axis=1)).ravel()
    scale = (n - 1) / (2.0 * s0) if method == "geary" and s0 > 0 else 1.0
    observed = _moran_geary_block(block0, W, method, row_sum) * scale
    hits = np.zeros(idx.size)
    valid_perms = np.zeros(idx.size)
    for _ in range(permutations):
        if multi:
            perm = sc.within_block_permutation(blocks, n, rng)
        else:
            perm = rng.permutation(n)
        block = block0[perm]
        # 打乱后仍用同一空间权重；行和需与打乱后的矩阵一致（对称图不变）
        stat = _moran_geary_block(block, W, method, row_sum) * scale
        finite = np.isfinite(stat) & np.isfinite(observed)
        valid_perms += np.isfinite(stat).astype(float)
        if method == "moran":
            hits += np.where(finite, stat >= observed, False).astype(float)
        else:
            hits += np.where(finite, stat <= observed, False).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(valid_perms > 0, (1.0 + hits) / (1.0 + valid_perms), np.nan)
    p = np.where(np.isfinite(observed), p, np.nan)
    return np.clip(p, 0.0, 1.0)


def _svg_table(ctx, params, outdir) -> dict:
    """svg 步骤/分析的核心：返回 {frame, store_records, moments, method, n_top, warnings, notes}。"""
    X, genes, spots = _matrix(ctx)
    coords = _coords(ctx)
    method = _p_str(params.get("method"), "moran").lower()
    if method not in METHOD_SVG:
        raise sc.bad_request(f"不支持的 SVG 方法：{method}（可选 moran / geary）")
    n_top = _p_int(params.get("n_top"), 30)
    if n_top <= 0:
        raise sc.prerequisite("n_top 必须为正整数（要输出的空间高变基因数量）")
    permutations = max(0, min(_p_int(params.get("permutations"), 0), 2000))
    wanted = _p_list(params.get("genes"))
    min_expr = max(0.0, _p_float(params.get("min_expr"), 0.0))
    warnings: list[str] = []
    notes: list[str] = []

    n_spots = X.shape[0]
    if n_spots > 60000:
        # 空间自相关仍可用全量；这里只提示内存
        warnings.append(f"spot 数较多（{n_spots}），Moran's I 逐基因块计算以控制内存")

    blocks, block_names, block_key, block_warns = _blocks(ctx)
    if len(blocks) > 1:
        warnings.append(
            f"多切片（{len(blocks)} 个，列「{block_key}」）：Moran/Geary 先做切片内中心化，"
            "置换只在切片内进行；解析正态近似只作参考")
    warnings.extend(block_warns)
    W = _radius_weights(ctx)
    Xd = X.astype(np.float32, copy=False)
    stats, zvals, moments = _svg_scores(Xd, W, method, blocks=blocks)
    pvals = _normal_p(zvals)
    p_adj = sc.bh_fdr(pvals)

    # 低表达基因不参与排名（默认 min_expr=0 时不过滤，保持与契约一致）
    mean_expr = np.asarray(Xd.mean(axis=0)).ravel()
    keep = np.isfinite(stats)
    if min_expr > 0:
        keep &= mean_expr >= min_expr
    stat_key = "moran_i" if method == "moran" else "geary_c"

    frame = pd.DataFrame(
        {
            "gene": genes,
            "statistic": stats,
            stat_key: stats,
            "z": zvals,
            "p_value": pvals,
            "p_adj": p_adj,
            "mean_expr": mean_expr,
        }
    )
    frame = frame.loc[keep].copy()
    if frame.empty:
        raise sc.prerequisite("没有可用于空间自相关计算的基因（表达矩阵可能为空或全为常数）")
    # 排序：显著优先，然后按统计量（Moran 越大越聚集；Geary 越小越聚集）
    if method == "moran":
        frame = frame.sort_values(["p_adj", "statistic"], ascending=[True, False], na_position="last")
    else:
        frame = frame.sort_values(["p_adj", "statistic"], ascending=[True, True], na_position="last")
    frame = frame.reset_index(drop=True)

    # 置换经验 p（只对候选 top 基因做，控制耗时；随后统一 BH）
    if permutations > 0 and len(frame) > 0:
        cand_n = int(min(len(frame), max(n_top * 3, 100), 400))
        cand_genes = frame["gene"].head(cand_n).tolist()
        cand_idx = np.array([genes.index(g) for g in cand_genes if g in genes], dtype=int)
        p_perm = _permutation_p(Xd, W, cand_idx, method, permutations, _rng(ctx), blocks=blocks)
        pos = frame.index[frame["gene"].isin([genes[i] for i in cand_idx])].to_numpy()
        frame.loc[pos, "p_value"] = p_perm[: len(pos)]
        frame["p_adj"] = sc.bh_fdr(frame["p_value"].to_numpy())
        notes.append(
            f"对前 {len(cand_idx)} 个候选基因做了 {permutations} 次 spot 置换经验 p 值"
            f"（Moran 单侧 greater / Geary 单侧 less），其余基因用正态近似"
        )

    wanted_idx: list[str] = []
    if wanted:
        gene_set = set(genes)
        wanted_idx = [g for g in wanted if g in gene_set]
        missing = [g for g in wanted if g not in gene_set]
        if missing:
            warnings.append(f"指定基因中有 {len(missing)} 个不在数据集里（如 {'、'.join(missing[:5])}）")
    top = frame.head(max(1, min(n_top, len(frame))))
    store_records = [
        {
            "gene": str(row["gene"]),
            "statistic": float(row["statistic"]),
            stat_key: float(row["statistic"]),
            "z": None if not np.isfinite(row["z"]) else float(row["z"]),
            "p": float(row["p_value"]),
            "p_adj": float(row["p_adj"]),
            "mean_expr": float(row["mean_expr"]),
        }
        for _, row in frame.head(2000).iterrows()
    ]
    return {
        "frame": frame,
        "top": top,
        "store": store_records,
        "method": method,
        "stat_key": stat_key,
        "moments": moments,
        "n_top": n_top,
        "permutations": permutations,
        "wanted": wanted_idx,
        "notes": notes,
        "warnings": warnings,
        "genes": genes,
    }


# ---------------------------------------------------------------------------
# 空间域：PCA + MRF / 图 Louvain（内置实现）
# ---------------------------------------------------------------------------

def _pca_features(ctx, X: np.ndarray, n_components: int = 15) -> tuple[np.ndarray, int, list[str]]:
    """取 PCA 前若干维（优先 ctx.obsm['pca']），标准化并截尾。"""
    notes: list[str] = []
    obsm = ctx.get("obsm") if isinstance(ctx.get("obsm"), dict) else {}
    pca = obsm.get("pca")
    n_comp = int(max(2, min(int(n_components), 50)))
    if pca is not None:
        arr = np.asarray(pca, dtype=float)
        if arr.ndim == 2 and arr.shape[0] == X.shape[0] and arr.shape[1] >= 2:
            feat = arr[:, :n_comp]
            notes.append(f"使用 PCA 前 {feat.shape[1]} 维（来自 obsm['pca']）")
            return _scale_features(feat), feat.shape[1], notes
    from sklearn.decomposition import PCA

    n_spots = X.shape[0]
    n_use = int(min(n_comp + 5, X.shape[1], max(2, n_spots - 1)))
    mat = np.asarray(X, dtype=np.float32)
    # 去掉零方差基因
    var = mat.var(axis=0)
    keep = np.isfinite(var) & (var > 1e-12)
    if not keep.any():
        raise sc.prerequisite("表达矩阵没有变异基因，无法计算空间域")
    mat = mat[:, keep]
    n_use = int(max(2, min(n_use, mat.shape[1])))
    try:
        model = PCA(n_components=n_use, svd_solver="randomized", random_state=0)
        feat = model.fit_transform(mat)
    except Exception as exc:  # noqa: BLE001
        raise sc.EngineError("internal", f"PCA 计算失败，无法做空间域：{exc}") from exc
    notes.append(f"现场计算 PCA 前 {feat.shape[1]} 维（来自表达矩阵，{int(keep.sum())} 个变异基因）")
    return _scale_features(feat), feat.shape[1], notes


def _scale_features(feat: np.ndarray) -> np.ndarray:
    arr = np.asarray(feat, dtype=float)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    mean = arr.mean(axis=0, keepdims=True)
    std = arr.std(axis=0, keepdims=True)
    std = np.where(std > 0, std, 1.0)
    out = (arr - mean) / std
    return np.clip(out, -10.0, 10.0)


def _fused_graph(ctx, feat: np.ndarray, n_neighbors: int = 10):
    """表达 kNN 图与空间 kNN 图的融合（相加后行标准化）。

    物理空间 kNN 在切片/样本块内取邻居（表达 kNN 允许跨块：表达相似就是相似）。
    """
    from scipy import sparse

    coords = _coords(ctx)
    n = coords.shape[0]
    blocks, _names, _key, _warns = _blocks(ctx)
    sp_idx, _ = _knn(coords, 6, blocks=blocks)
    g_sp = _sparse_row_norm(sp_idx, n)
    k = int(max(2, min(n_neighbors, n - 1)))
    if feat.shape[1] >= 2:
        from sklearn.neighbors import NearestNeighbors

        nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(feat)
        _, idx = nn.kneighbors(feat)
        g_ex = _sparse_row_norm(idx[:, 1:].astype(int), n)
    else:
        g_ex = g_sp
    fused = g_sp + g_ex
    row = np.asarray(fused.sum(axis=1)).ravel()
    row[row == 0] = 1.0
    return sparse.diags(1.0 / row) @ fused


def _kmeans(feat: np.ndarray, k: int, seed: int = 0):
    from sklearn.cluster import KMeans

    k = int(max(2, min(k, feat.shape[0])))
    model = KMeans(n_clusters=k, n_init=10, random_state=seed)
    labels = model.fit_predict(feat)
    return labels, np.asarray(model.cluster_centers_, dtype=float)


def _mrf_domains(feat: np.ndarray, G, k: int, smooth: float, iterations: int):
    """MRF 迭代平滑（BayesSpace 2021 思想）。

    label_i ← argmax( (1-smooth)·sim_expr(i,c) + smooth·mean_{j∈N(i)} 1[label_j = c] )
    sim_expr 用「到域中心的负平方距离 softmax」表示（表达相似度），中心随标签更新。
    """
    n = feat.shape[0]
    labels, centers = _kmeans(feat, k, seed=0)
    smooth = float(np.clip(smooth, 0.0, 1.0))
    onehot = np.zeros((n, centers.shape[0]))
    n_iter = 0
    for step in range(max(1, int(iterations))):
        n_iter = step + 1
        # 表达相似度（softmax of negative squared distance）
        dist2 = ((feat[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
        scale = np.median(np.sqrt(np.maximum(dist2.min(axis=1), 0.0))) or 1.0
        sim = np.exp(-dist2 / (2.0 * scale * scale))
        sim = sim / np.maximum(sim.sum(axis=1, keepdims=True), 1e-12)
        # 邻居标签频率
        onehot[:] = 0.0
        onehot[np.arange(n), labels] = 1.0
        neigh = np.asarray(G @ onehot)
        score = (1.0 - smooth) * sim + smooth * neigh
        new_labels = np.argmax(score, axis=1)
        # 更新中心（空域保留旧中心）
        for c in range(centers.shape[0]):
            mask = new_labels == c
            if mask.any():
                centers[c] = feat[mask].mean(axis=0)
        if np.array_equal(new_labels, labels):
            labels = new_labels
            break
        labels = new_labels
    return labels, n_iter


def _louvain_domains(feat: np.ndarray, G) -> tuple[np.ndarray, int]:
    """融合图上的 Louvain 社区发现（networkx 内置实现）。"""
    import networkx as nx

    try:
        from scipy import sparse

        coo = sparse.triu(G, k=1).tocoo()
        edges = [(int(i), int(j), float(w)) for i, j, w in zip(coo.row, coo.col, coo.data) if w > 0]
    except Exception as exc:  # noqa: BLE001
        raise sc.EngineError("internal", f"构造融合图失败：{exc}") from exc
    graph = nx.Graph()
    graph.add_nodes_from(range(feat.shape[0]))
    graph.add_weighted_edges_from(edges)
    if graph.number_of_edges() == 0:
        raise sc.prerequisite("空间邻接图为空（坐标可能重复或只有 1 个 spot），无法做图聚类")
    communities = nx.community.louvain_communities(graph, seed=0)
    communities = sorted(communities, key=len, reverse=True)
    labels = np.zeros(feat.shape[0], dtype=int)
    for idx, members in enumerate(communities):
        for node in members:
            labels[int(node)] = idx
    return labels, len(communities)


def compute_domains(ctx, params=None, outdir=None, *, write_back: bool = False) -> dict:
    """计算空间域（不写回）；供 ``domains`` 步骤与分析共用。

    返回 {labels, n_domains, method, sizes, n_iter, notes, warnings, pca_dims}。
    """
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    X, genes, spots = _matrix(ctx)
    coords = _coords(ctx)
    n_spots = X.shape[0]
    if n_spots < 6:
        raise sc.prerequisite(f"只有 {n_spots} 个 spot，无法划分空间域（至少需要 6 个）")
    method = _p_str(params.get("method"), "mrf").lower()
    if method not in METHOD_DOMAIN:
        raise sc.bad_request(f"不支持的空间域方法：{method}（可选 mrf / graph）")
    n_domains = _p_int(params.get("n_domains"), 6)
    n_domains = int(max(2, min(n_domains, max(2, n_spots // 5))))
    smooth = float(np.clip(_p_float(params.get("smooth"), 0.5), 0.0, 1.0))
    iterations = int(max(1, min(_p_int(params.get("iterations"), 20), 200)))
    notes: list[str] = []
    warnings: list[str] = []

    feat, pca_dims, pca_notes = _pca_features(ctx, X, 15)
    notes.extend(pca_notes)

    if method == "graph":
        if n_spots > 20000:
            warnings.append(
                f"spot 数（{n_spots}）较大，Louvain 计算成本高；已退回 MRF 平滑（method=graph 仅适用于 ≤2 万 spot）"
            )
            method = "mrf"
        else:
            G = _fused_graph(ctx, feat, n_neighbors=10)
            labels, n_comm = _louvain_domains(feat, G)
            notes.append(f"融合图（空间 6NN + 表达 10NN）上的 Louvain 社区发现，得到 {n_comm} 个域")
            if n_comm != n_domains:
                notes.append(f"Louvain 自动决定域数（{n_comm}），参数 n_domains={n_domains} 未强制生效")
            n_iter = 0
            n_domains_eff = n_comm
    if method == "mrf":
        G = _fused_graph(ctx, feat, n_neighbors=10)
        raw_labels, n_iter = _mrf_domains(feat, G, n_domains, smooth, iterations)
        n_domains_eff = int(len(np.unique(raw_labels)))
        labels = raw_labels
        notes.append(
            f"PCA {pca_dims} 维 + 空间/表达融合图，k-means 初始化 {n_domains} 个域后做 MRF 平滑 "
            f"{n_iter} 轮（smooth={smooth:g}，收敛{'是' if n_iter < iterations else '否（达到迭代上限）'}）"
        )

    # 域编号按规模降序重排（domain 1 = 最大域），便于跨样本比较
    uniq, counts = np.unique(labels, return_counts=True)
    order = np.argsort(-counts)
    remap = {int(uniq[i]): rank + 1 for rank, i in enumerate(order)}
    labels_num = np.array([remap[int(v)] for v in labels], dtype=int)
    sizes = {str(rank + 1): int(counts[i]) for rank, i in enumerate(order)}
    labels_str = [str(int(v)) for v in labels_num]
    return {
        "labels": labels_str,
        "labels_num": labels_num,
        "n_domains": int(len(uniq)),
        "n_domains_requested": n_domains,
        "method": method,
        "sizes": sizes,
        "n_iter": n_iter,
        "pca_dims": pca_dims,
        "notes": notes,
        "warnings": warnings,
        "outdir": outdir,
    }


def get_domain_labels(ctx, params=None, *, allow_compute: bool = True) -> tuple[list[str], str, dict]:
    """取空间域标签：优先用数据集里已有结果，没有则按默认参数现场计算。

    返回 (labels, 来源说明, 计算信息 dict)；allow_compute=False 时缺结果抛 prerequisite。
    """
    n = len(ctx.get("spots") or [])
    labels, hit = _domain_labels(ctx, n)
    if labels is not None:
        return labels, f"数据集列「{hit}」", {}
    if not allow_compute:
        raise sc.prerequisite("还没有空间域结果，请先运行「空间域」(domains) 步骤")
    computed = compute_domains(ctx, params or {}, (ctx or {}).get("outdir"), write_back=False)
    note = (
        f"未找到空间域结果，已现场计算（method={computed['method']}，"
        f"n_domains={computed['n_domains']}）"
    )
    return computed["labels"], note, computed


# ---------------------------------------------------------------------------
# 差异 / 标记基因（Wilcoxon 秩和 + BH，scipy）
# ---------------------------------------------------------------------------

def _wilcoxon_frame(X: np.ndarray, genes: list[str], case_mask: np.ndarray, ctrl_mask: np.ndarray) -> pd.DataFrame:
    """两组 Wilcoxon 秩和检验（向量化，scipy mannwhitneyu axis=0）+ BH FDR + logFC。"""
    n_case = int(np.sum(case_mask))
    n_ctrl = int(np.sum(ctrl_mask))
    if n_case < 3 or n_ctrl < 3:
        raise sc.prerequisite(
            f"两组 spot 数不足（{n_case} vs {n_ctrl}），至少各需 3 个才能做 Wilcoxon 检验"
        )
    A = np.asarray(X[case_mask], dtype=float)
    B = np.asarray(X[ctrl_mask], dtype=float)
    try:
        from scipy import stats

        _, pvals = stats.mannwhitneyu(A, B, alternative="two-sided", axis=0, method="asymptotic")
    except Exception as exc:  # noqa: BLE001
        raise sc.EngineError("internal", f"Wilcoxon 检验失败：{exc}") from exc
    pvals = np.asarray(pvals, dtype=float).ravel()
    pvals = np.where(np.isfinite(pvals), pvals, 1.0)
    pvals = np.clip(pvals, 0.0, 1.0)
    mean_case = A.mean(axis=0)
    mean_ctrl = B.mean(axis=0)
    frame = pd.DataFrame(
        {
            "gene": genes,
            "logFC": mean_case - mean_ctrl,
            "mean_case": mean_case,
            "mean_ctrl": mean_ctrl,
            "pct_case": (A > 0).mean(axis=0) * 100.0,
            "pct_ctrl": (B > 0).mean(axis=0) * 100.0,
            "p_value": pvals,
            "p_adj": sc.bh_fdr(pvals),
        }
    )
    frame["regulation"] = np.where(
        (frame["p_adj"] <= 0.05) & (frame["logFC"] > 0),
        "up",
        np.where((frame["p_adj"] <= 0.05) & (frame["logFC"] < 0), "down", "ns"),
    )
    return frame


def _marker_frame(X: np.ndarray, genes: list[str], labels: np.ndarray, target: str) -> pd.DataFrame:
    """目标域 vs 其余（Wilcoxon + BH + logFC）。"""
    mask = labels == target
    return _wilcoxon_frame(X, genes, mask, ~mask)


def _de_columns(stat_label: str = "Wilcoxon") -> list[dict]:
    return [
        {"key": "gene", "label": "基因", "type": "string"},
        {"key": "logFC", "label": "logFC（均值差）", "type": "number"},
        {"key": "mean_case", "label": "组内均值", "type": "number"},
        {"key": "mean_ctrl", "label": "组外均值", "type": "number"},
        {"key": "pct_case", "label": "组内表达比例(%)", "type": "number"},
        {"key": "pct_ctrl", "label": "组外表达比例(%)", "type": "number"},
        {"key": "p_value", "label": f"{stat_label} P", "type": "number"},
        {"key": "p_adj", "label": "BH FDR", "type": "number"},
        {"key": "regulation", "label": "方向", "type": "string"},
    ]


def _de_head(frame: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """按显著性挑前 top_n 行（每个方向各一半）。"""
    top_n = int(max(1, min(top_n, 2000)))
    sig = frame[frame["regulation"] != "ns"]
    base = sig if not sig.empty else frame
    ordered = base.assign(_abs=base["logFC"].abs()).sort_values(["p_adj", "_abs"], ascending=[True, False])
    return ordered.drop(columns=["_abs"]).head(top_n)


# ---------------------------------------------------------------------------
# 通用绘图工具
# ---------------------------------------------------------------------------

def _subsample_idx(n: int, max_points: int, rng) -> np.ndarray:
    """确定性下采样索引（用于大数据的空间散点）。"""
    limit = int(max(100, max_points))
    if n <= limit:
        return np.arange(n)
    return np.sort(rng.choice(n, size=limit, replace=False))


def _point_size(params: dict) -> float:
    return float(max(1.0, _p_float((params or {}).get("point_size"), 8.0)))


def _max_points(params: dict) -> int:
    return int(max(500, _p_int((params or {}).get("max_points"), MAX_PLOT_POINTS)))


def _want_image(params: dict) -> bool:
    return _p_bool((params or {}).get("image"), False)


def _spatial_ax(ax, ctx, params, values=None, labels=None, legend_title=None):
    rng = _rng(ctx)
    coords = _coords(ctx)
    n = coords.shape[0]
    take = _subsample_idx(n, _max_points(params), rng)
    if labels is not None:
        labs = [labels[i] for i in take]
    else:
        labs = None
    vals = None if values is None else np.asarray(values, dtype=float)[take]
    return sc.scatter_spatial(
        ax,
        ctx,
        values=vals,
        labels=labs,
        point_size=_point_size(params),
        legend_title=legend_title,
        image=_want_image(params),
    ), take


def _fig_od(ctx, params, outdir, name: str, title: str) -> str:
    return _outdir(ctx, outdir)


def _fig_save(fig, ctx, params, outdir, name: str, title: str) -> dict:
    dpi = int(max(80, min(_p_int((params or {}).get("dpi"), 140), 300)))
    outdir = _outdir(ctx, outdir)
    return sc.save_fig(fig, outdir, name, title=title, dpi=dpi)


# ---------------------------------------------------------------------------
# svg 步骤 / 相关图
# ---------------------------------------------------------------------------

def _svg_rank_figure(ctx, params, outdir, *, computed: dict | None = None) -> dict:
    """Top SVG 排序条形图 + 最显著基因的空间表达图（英文图注）。"""
    if computed is None:
        computed = _svg_table(ctx, params or {}, outdir)
    frame = computed["top"]
    method = computed["method"]
    stat_label = "Moran's I" if method == "moran" else "Geary's C"
    genes = [str(g) for g in frame["gene"].tolist()][:20]
    values = frame["statistic"].to_numpy(dtype=float)[:20]
    pvals = frame["p_adj"].to_numpy(dtype=float)[:20]
    fig, axes = __import__("matplotlib.pyplot", fromlist=["pyplot"]).subplots(
        1, 2, figsize=(11.6, 5.0), gridspec_kw={"width_ratios": [1.45, 1.0]}
    )
    ax = axes[0]
    ypos = np.arange(len(genes))[::-1]
    colors = ["#c0392b" if (np.isfinite(p) and p <= 0.05) else "#7f8c9b" for p in pvals]
    ax.barh(ypos, values, color=colors, height=0.72)
    ax.set_yticks(ypos)
    ax.set_yticklabels(genes, fontsize=7.5)
    ax.set_xlabel(stat_label)
    ax.set_title(f"Top {len(genes)} spatially variable genes ({method})", fontsize=10)
    for y, value, p in zip(ypos, values, pvals):
        text = f"{value:.3f}" + ("*" if np.isfinite(p) and p <= 0.05 else "")
        ax.text(value, y, "  " + text, va="center", fontsize=6.2, color="#2b2f36")
    ax.grid(axis="y", visible=False)
    # 右图：最显著基因的空间表达
    ax2 = axes[1]
    X, all_genes, spots = _matrix(ctx)
    gene = genes[0] if genes else ""
    if gene and gene in all_genes:
        col = X[:, all_genes.index(gene)].astype(float)
        _spatial_ax(ax2, ctx, params or {}, values=col, legend_title=f"{gene} expression")
        ax2.set_title(f"Spatial expression: {gene}", fontsize=10)
    else:
        ax2.axis("off")
    return _fig_save(fig, ctx, params or {}, outdir, "svg_rank_plot", "Top spatially variable genes")


def _svg_store(ctx, computed: dict) -> None:
    """把 svg 结果写进 ctx['results']，方便同一进程内后续分析直接复用。"""
    results = ctx.get("results")
    if isinstance(results, dict):
        results["svg"] = computed["store"]


def run_svg(ctx, params=None, outdir=None) -> dict:
    """步骤 ``svg``：空间高变基因（Moran's I / Geary's C，内置实现）。"""
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    # 优先走 sp_spatial 的 squidpy 实现（squidpy.gr.spatial_autocorr）；不可用时回退内置。
    try:
        import sp_spatial  # 延迟导入，避免模块级循环依赖

        mode = str(params.get("method") or "moran").lower()
        runner = sp_spatial.run_gearys_svg if mode.startswith("geary") else sp_spatial.run_morans_svg
        payload = runner(ctx, params, outdir)
        summary = payload.get("summary") or {}
        if summary.get("method") and "内置实现" not in str(summary.get("method")):
            summary["title"] = "空间高变基因"
            payload["summary"] = summary
            return payload
    except Exception:  # noqa: BLE001 - 委托失败就回退内置
        pass
    computed = _svg_table(ctx, params, outdir)
    frame = computed["frame"]
    top = computed["top"]
    method = computed["method"]
    stat_key = computed["stat_key"]
    stat_label = "Moran's I" if method == "moran" else "Geary's C"

    full_path = sc.export_csv(frame, outdir, "svg_all", index=False)
    top_path = sc.export_csv(top, outdir, "svg_top", index=False)
    columns = [
        {"key": "gene", "label": "基因", "type": "string"},
        {"key": "statistic", "label": stat_label, "type": "number"},
        {"key": "z", "label": "z 值", "type": "number"},
        {"key": "p_value", "label": "P 值", "type": "number"},
        {"key": "p_adj", "label": "BH FDR", "type": "number"},
        {"key": "mean_expr", "label": "平均表达", "type": "number"},
    ]
    table_obj = sc.table(
        "svg",
        "空间高变基因（按显著性排序）",
        columns,
        top,
        note=f"{stat_label}：共 {len(frame)} 个基因，完整表见 svg_all.csv；{_svg_note(method)}",
    )
    figures: list[dict] = []
    warnings = list(computed["warnings"])
    try:
        figures.append(_svg_rank_figure(ctx, params, outdir, computed=computed))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"SVG 排序图绘制失败：{exc}")

    n_sig = int((frame["p_adj"] <= 0.05).sum())
    bullets = [
        f"{stat_label}（行标准化空间权重，半径 {_radius_text(ctx)}）：显著基因 {n_sig} / {len(frame)}（FDR<0.05）",
        "最显著基因：" + "、".join(
            f"{row['gene']}({row['statistic']:.3f})" for _, row in top.head(5).iterrows()
        ),
    ]
    metrics = {
        "method": method,
        "n_genes_tested": int(len(frame)),
        "n_sig": n_sig,
        "n_top": int(len(top)),
        "permutations": computed["permutations"],
        "expect": float(computed["moments"].get("expect", float("nan"))) if math.isfinite(
            float(computed["moments"].get("expect", float("nan")))
        ) else None,
        "top_gene": str(top["gene"].iloc[0]) if len(top) else "",
        "top_statistic": float(top["statistic"].iloc[0]) if len(top) else None,
    }
    notes = list(computed["notes"]) + [
        _svg_note(method),
        "空间权重来自 import 时估计的邻接半径，行标准化后计算（内置实现，不调用 squidpy）",
    ]
    out = sc.result(
        title="空间高变基因",
        bullets=bullets,
        metrics=metrics,
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=notes,
        data_out=[
            sc.data_file("svg_csv", full_path, "全部基因空间自相关"),
            sc.data_file("svg_top_csv", top_path, "Top SVG"),
        ],
    )
    out["store"] = {"svg": computed["store"]}
    _svg_store(ctx, computed)
    return out


def _svg_note(method: str) -> str:
    if method == "moran":
        return "方法：Moran's I（Moran 1950）+ 行标准化空间权重 + 正态近似 p（可选置换经验 p）；内置实现"
    return "方法：Geary's C（Geary 1954）+ 行标准化空间权重 + 正态近似 p（C<1 表示正空间自相关）；内置实现"


def _radius_text(ctx) -> str:
    """邻接半径的文字说明（多切片时按块估计，与合并坐标估计可能不同）。"""
    blocks, _names, _key, _warns = _blocks(ctx)
    try:
        if len(blocks) > 1:
            return f"按 {len(blocks)} 个切片分块估计 {sc.blocked_radius_estimate(_coords(ctx), blocks):.1f}"
        radius = float(ctx.get("neighbor_radius") or 0.0)
    except (TypeError, ValueError):
        radius = 0.0
    return f"{radius:.1f}" if radius > 0 else "自动估计"


def run_domains(ctx, params=None, outdir=None) -> dict:
    """步骤 ``domains``：空间域（MRF 平滑 / 融合图 Louvain），写回 ``domains`` 列。"""
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    computed = compute_domains(ctx, params, outdir, write_back=True)
    labels = computed["labels"]
    _method_text = ("MRF 空间域（BayesSpace 2021 思想：PCA + 空间邻接图 + 马尔可夫随机场迭代；"
                    "内置实现，squidpy 无空间域 API）")
    X, genes, spots = _matrix(ctx)
    labels_num = np.asarray(labels, dtype=object)
    uniq = list(dict.fromkeys([str(v) for v in labels]))
    sizes = computed["sizes"]

    size_frame = pd.DataFrame(
        [
            {
                "domain": str(d),
                "n_spots": int(np.sum(labels_num == str(d))),
                "pct": round(float(np.mean(labels_num == str(d)) * 100.0), 2),
            }
            for d in uniq
        ]
    )
    size_frame = size_frame.sort_values("n_spots", ascending=False).reset_index(drop=True)
    table_obj = sc.table(
        "domains",
        "空间域大小",
        [
            {"key": "domain", "label": "域", "type": "string"},
            {"key": "n_spots", "label": "spot 数", "type": "number"},
            {"key": "pct", "label": "占比(%)", "type": "number"},
        ],
        size_frame,
        note=f"共 {len(uniq)} 个域；编号按规模降序（域 1 = 最大）",
    )
    figures: list[dict] = []
    warnings = list(computed["warnings"])
    try:
        figures.append(_domain_map_figure(ctx, params, outdir, labels=labels))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"空间域图绘制失败：{exc}")
    csv_path = sc.export_csv(
        pd.DataFrame({"spot": [str(s) for s in spots], "domains": labels}), outdir, "domains", index=False
    )
    bullets = [
        f"划分出 {len(uniq)} 个空间域（方法 {computed['method']}，"
        f"PCA {computed['pca_dims']} 维 + 空间邻接图）",
        "域大小：" + "、".join(f"域{d}({sizes.get(str(d), 0)})" for d in uniq[:8]),
    ]
    metrics = {
        "method": computed["method"],
        "n_domains": int(len(uniq)),
        "n_domains_requested": int(computed["n_domains_requested"]),
        "iterations_used": int(computed["n_iter"]),
        "pca_dims": int(computed["pca_dims"]),
        "sizes": {str(k): int(v) for k, v in sizes.items()},
    }
    notes = list(computed["notes"]) + [
        "方法：PCA + 空间邻接图 + k-means 初始化 + 马尔可夫随机场（MRF）迭代平滑"
        "（BayesSpace 2021 思想；内置实现，不调用 BayesSpace/squidpy）",
        "域编号按大小降序重排，便于跨切片比较；结果写回 spot 级列 domains",
    ]
    try:
        summary_obj = {"title": "空间域", "bullets": bullets, "method": _method_text}
    except NameError:  # pragma: no cover - 防御
        summary_obj = None
    out = sc.result(
        summary=summary_obj or None,
        title="空间域" if summary_obj is None else None,
        bullets=None if summary_obj is not None else bullets,
        metrics=metrics,
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=notes,
        data_out=[sc.data_file("domains_csv", csv_path, "spot × 空间域")],
    )
    out["cluster_store"] = {"key": "domains", "labels": labels}
    out["store"] = {
        "domains": {
            "key": "domains",
            "labels": labels,
            "n_domains": int(len(uniq)),
            "method": computed["method"],
            "sizes": {str(k): int(v) for k, v in sizes.items()},
        }
    }
    results = ctx.get("results")
    if isinstance(results, dict):
        results["domains"] = out["store"]["domains"]
    return out


# ---------------------------------------------------------------------------
# 图：空间域 / 邻域 / 共现 / 边界 / 距离 / niche / 多样性 / 差异
# ---------------------------------------------------------------------------

def _domain_map_figure(ctx, params, outdir, *, labels: list[str] | None = None, title: str | None = None) -> dict:
    if labels is None:
        labels, src, _info = get_domain_labels(ctx, params, allow_compute=True)
        suffix = f" ({src})" if src.startswith("数据集") else ""
    else:
        suffix = ""
    fig, ax = sc.new_fig(6.6, 5.6)
    _spatial_ax(ax, ctx, params, labels=labels, legend_title="Domain")
    counts = pd.Series(labels).value_counts()
    ax.set_title(
        str(title or "Spatial domains") + f" (n={len(counts)})" + suffix,
        fontsize=10,
    )
    return _fig_save(fig, ctx, params, outdir, "domain_map", "Spatial domains")


def _cluster_pairs(labels: np.ndarray, k_cat: int) -> np.ndarray:
    """邻接计数矩阵（对称化）。"""
    return labels


def _nhood_enrichment_core(labels: list[str], coords: np.ndarray, n_perms: int, n_neighbors: int, rng,
                           blocks=None):
    """邻域富集：邻接计数 + 置换 z-score（squidpy nhood_enrichment 思想，内置实现）。

    多切片时邻居只在切片内取、标签置换也限制在切片内（不跨片交换），否则邻域计数被切片结构污染。
    """
    labels_arr = np.asarray([str(v) for v in labels], dtype=object)
    uniq = list(dict.fromkeys(labels_arr.tolist()))
    index = {name: i for i, name in enumerate(uniq)}
    codes = np.array([index[v] for v in labels_arr], dtype=int)
    k_cat = len(uniq)
    idx, _ = _knn(coords, n_neighbors, blocks=blocks)
    n = coords.shape[0]

    def _counts(code_vec: np.ndarray) -> np.ndarray:
        flat = (code_vec[:, None] * k_cat + code_vec[idx]).ravel()
        mat = np.bincount(flat, minlength=k_cat * k_cat).reshape(k_cat, k_cat).astype(float)
        return mat + mat.T

    observed = _counts(codes)
    n_perms = int(max(1, min(n_perms, 1000)))
    null = np.zeros((n_perms, k_cat, k_cat))
    for t in range(n_perms):
        null[t] = _counts(rng.permutation(codes))
    mean = null.mean(axis=0)
    std = null.std(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(std > 0, (observed - mean) / std, 0.0)
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    return {
        "uniq": uniq,
        "observed": observed,
        "null_mean": mean,
        "null_std": std,
        "z": z,
        "n_perms": n_perms,
        "n_neighbors": int(idx.shape[1]),
        "n_spots": n,
    }


def _nhood_table(res: dict) -> pd.DataFrame:
    uniq = res["uniq"]
    rows = []
    z = res["z"]
    obs = res["observed"]
    for i in range(len(uniq)):
        for j in range(i, len(uniq)):
            rows.append(
                {
                    "cluster_a": str(uniq[i]),
                    "cluster_b": str(uniq[j]),
                    "observed": float(obs[i, j]),
                    "null_mean": round(float(res["null_mean"][i, j]), 3),
                    "z": round(float(z[i, j]), 3),
                }
            )
    frame = pd.DataFrame(rows)
    return frame.sort_values("z", ascending=False).reset_index(drop=True)


def _nhood_figure(ctx, params, outdir, *, res: dict | None = None) -> dict:
    if res is None:
        labels, key = _cluster_labels(ctx, _p_str((params or {}).get("cluster_key")))
        coords = _coords(ctx)
        res = _nhood_enrichment_core(
            labels,
            coords,
            n_perms=_p_int((params or {}).get("n_perms"), 50),
            n_neighbors=_p_int((params or {}).get("n_neighbors"), 6),
            rng=_rng(ctx),
            blocks=_blocks(ctx)[0],
        )
    uniq = [str(u) for u in res["uniq"]]
    z = np.clip(res["z"], -8.0, 8.0)
    fig, axes = __import__("matplotlib.pyplot", fromlist=["pyplot"]).subplots(
        1, 2, figsize=(11.4, 4.9), gridspec_kw={"width_ratios": [1.0, 1.0]}
    )
    ax = axes[0]
    im = ax.imshow(z, cmap="RdBu_r", vmin=-8, vmax=8, aspect="auto")
    ax.set_xticks(range(len(uniq)))
    ax.set_xticklabels(uniq, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(uniq)))
    ax.set_yticklabels(uniq, fontsize=8)
    ax.set_title(f"Neighbourhood enrichment z-score (n={res['n_perms']} perms)", fontsize=10)
    ax.grid(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("z", fontsize=8)
    ax2 = axes[1]
    obs = res["observed"]
    im2 = ax2.imshow(np.log1p(obs), cmap="YlGnBu", aspect="auto")
    ax2.set_xticks(range(len(uniq)))
    ax2.set_xticklabels(uniq, rotation=45, ha="right", fontsize=8)
    ax2.set_yticks(range(len(uniq)))
    ax2.set_yticklabels(uniq, fontsize=8)
    ax2.set_title("Observed adjacency counts (log1p)", fontsize=10)
    ax2.grid(False)
    cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.02)
    cbar2.set_label("log1p(count)", fontsize=8)
    return _fig_save(fig, ctx, params or {}, outdir, "nhood_enrich", "Neighbourhood enrichment")


def _co_occurrence_core(labels: list[str], coords: np.ndarray, n_rings: int, interval=None,
                        blocks=None) -> dict:
    """环形分箱条件共现概率 P(label_b | label_a at ring r)（squidpy co_occurrence 思想，内置实现）。"""
    labels_arr = np.asarray([str(v) for v in labels], dtype=object)
    uniq = list(dict.fromkeys(labels_arr.tolist()))
    index = {name: i for i, name in enumerate(uniq)}
    codes = np.array([index[v] for v in labels_arr], dtype=int)
    k_cat = len(uniq)
    n_rings = int(max(1, min(n_rings, 20)))
    n = coords.shape[0]
    # 需要足够的邻居覆盖多个环
    k = int(min(n - 1, max(12, n_rings * 6)))
    if k < 2:
        raise sc.prerequisite("spot 数太少，无法计算共现（至少需要 3 个 spot）")
    idx, dist = _knn(coords, k, blocks=blocks)
    dmin = float(np.min(dist)) if dist.size else 0.0
    if interval is not None:
        try:
            lo, hi = float(interval[0]), float(interval[1])
        except (TypeError, IndexError, ValueError):
            raise sc.bad_request("interval 需要 [起始距离, 结束距离] 两个数字")
        if not (np.isfinite(lo) and np.isfinite(hi) and hi > lo):
            raise sc.bad_request("interval 必须是递增的两个有限数字（hi > lo）")
    else:
        # 用 5%–95% 分位避开极端值，保证环宽稳定
        finite = dist[np.isfinite(dist)]
        lo = float(np.percentile(finite, 5)) if finite.size else 0.0
        hi = float(np.percentile(finite, 95)) if finite.size else 1.0
        if not hi > lo:
            hi = lo + max(1e-6, (lo or 1.0) * 0.1)
    edges = np.linspace(lo, hi, n_rings + 1)
    ring_of = np.digitize(dist, edges) - 1
    valid = (ring_of >= 0) & (ring_of < n_rings) & np.isfinite(dist)
    # 双向计数（(i,j) 与 (j,i) 都计入 a→b），分母为 a 类在所有环上的邻居总数
    counts = np.zeros((k_cat, k_cat, n_rings))
    neigh_codes = codes[idx]
    for r in range(n_rings):
        mask = valid & (ring_of == r)
        if not mask.any():
            continue
        a = np.broadcast_to(codes[:, None], mask.shape)[mask]
        b = neigh_codes[mask]
        np.add.at(counts, (a, b, np.full(a.size, r)), 1.0)
        np.add.at(counts, (b, a, np.full(a.size, r)), 1.0)
    totals = counts.sum(axis=1)  # (a, r) 的邻居总数
    with np.errstate(divide="ignore", invalid="ignore"):
        prob = np.where(totals[:, None, :] > 0, counts / np.maximum(totals[:, None, :], 1e-12), 0.0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return {
        "uniq": uniq,
        "prob": np.clip(prob, 0.0, 1.0),
        "counts": counts,
        "totals": totals,
        "edges": edges,
        "centers": centers,
        "n_rings": n_rings,
        "n_neighbors": int(idx.shape[1]),
        "n_spots": n,
    }


def _cooccurrence_table(res: dict) -> pd.DataFrame:
    uniq = res["uniq"]
    rows = []
    for i, a in enumerate(uniq):
        for j, b in enumerate(uniq):
            for r in range(res["n_rings"]):
                rows.append(
                    {
                        "cluster_a": str(a),
                        "cluster_b": str(b),
                        "ring": r + 1,
                        "distance": round(float(res["centers"][r]), 4),
                        "prob": round(float(res["prob"][i, j, r]), 5),
                        "n_pairs": int(res["counts"][i, j, r]),
                    }
                )
    return pd.DataFrame(rows)


def _cooccurrence_figure(ctx, params, outdir, *, res: dict | None = None) -> dict:
    import matplotlib.pyplot as plt

    if res is None:
        labels, _key = _cluster_labels(ctx, _p_str((params or {}).get("cluster_key")))
        interval = (params or {}).get("interval")
        res = _co_occurrence_core(
            labels,
            _coords(ctx),
            n_rings=_p_int((params or {}).get("n_rings"), 5),
            interval=interval if interval not in (None, "", []) else None,
            blocks=_blocks(ctx)[0],
        )
    uniq = [str(u) for u in res["uniq"]]
    if len(uniq) > 6:
        keep = np.argsort(-res["totals"].sum(axis=1))[:6]
        keep = sorted(keep.tolist())
        uniq_show = [uniq[i] for i in keep]
    else:
        keep = list(range(len(uniq)))
        uniq_show = uniq
    n_panel = max(1, len(keep))
    ncol = min(3, n_panel)
    nrow = int(math.ceil(n_panel / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.3 * nrow), squeeze=False)
    colors = sc.palette(len(uniq))
    centers = res["centers"]
    for pos, i in enumerate(keep):
        ax = axes[pos // ncol][pos % ncol]
        for j, b in enumerate(uniq):
            ax.plot(centers, res["prob"][i, j, :], marker="o", ms=3.2, lw=1.2,
                    color=colors[j % len(colors)], label=str(b))
        ax.set_title(f"Anchor: {uniq[i]}", fontsize=9)
        ax.set_xlabel("Distance", fontsize=8)
        ax.set_ylabel("P(neighbour | anchor)", fontsize=8)
        ax.set_ylim(0.0, 1.02)
        ax.tick_params(labelsize=7)
    for extra in range(n_panel, nrow * ncol):
        axes[extra // ncol][extra % ncol].axis("off")
    try:
        handles, labels_ = axes[0][0].get_legend_handles_labels()
        fig.legend(handles, labels_, loc="upper center", ncol=min(6, len(uniq)), fontsize=7,
                   frameon=False, bbox_to_anchor=(0.5, 1.02), title="Neighbour cluster",
                   title_fontsize=7)
    except Exception:  # noqa: BLE001
        pass
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _fig_save(fig, ctx, params or {}, outdir, "cooccurrence", "Co-occurrence")


def _composition_table(labels: list[str], clusters: list[str]) -> pd.DataFrame:
    """每个 niche/域里各 cluster 的平均组成比例。"""
    lab = np.asarray(labels, dtype=object)
    clu = np.asarray([str(c) for c in clusters], dtype=object)
    clu_uniq = list(dict.fromkeys(clu.tolist()))
    rows = []
    for name in list(dict.fromkeys(lab.tolist())):
        mask = lab == name
        if not mask.any():
            continue
        row = {"label": str(name), "n_spots": int(mask.sum())}
        for c in clu_uniq:
            row[str(c)] = round(float(np.mean(clu[mask] == c)), 4)
        rows.append(row)
    return pd.DataFrame(rows)


def _neighborhood_composition(labels: list[str], coords: np.ndarray, k: int, blocks=None):
    """每个 spot 的邻域类别组成向量。"""
    labels_arr = np.asarray([str(v) for v in labels], dtype=object)
    uniq = list(dict.fromkeys(labels_arr.tolist()))
    index = {name: i for i, name in enumerate(uniq)}
    codes = np.array([index[v] for v in labels_arr], dtype=int)
    idx, _ = _knn(coords, k, blocks=blocks)
    k_eff = idx.shape[1]
    onehot = np.zeros((len(labels_arr), len(uniq)))
    onehot[np.arange(len(labels_arr)), codes] = 1.0
    comp = onehot[idx].mean(axis=1) if k_eff else onehot
    return comp, uniq


def _diversity_core(labels: list[str], coords: np.ndarray, k: int, blocks=None) -> dict:
    comp, uniq = _neighborhood_composition(labels, coords, k, blocks=blocks)
    p = np.clip(comp, 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        shannon = -np.sum(np.where(p > 0, p * np.log(np.maximum(p, 1e-12)), 0.0), axis=1)
    simpson = 1.0 - np.sum(p * p, axis=1)
    labels_arr = np.asarray([str(v) for v in labels], dtype=object)
    return {
        "comp": comp,
        "uniq": uniq,
        "shannon": shannon,
        "simpson": simpson,
        "labels": labels_arr,
        "k": int(_knn(coords, k, blocks=blocks)[0].shape[1]),
    }


def _boundary_core(ctx, params) -> dict:
    """边界强度：跨域邻接比例 + 环形表达梯度。"""
    labels, src, _info = get_domain_labels(ctx, params)
    coords = _coords(ctx)
    uniq = list(dict.fromkeys(labels))
    da = _require_label((params or {}).get("domain_a", 1), uniq, "域 A")
    db = _require_label((params or {}).get("domain_b", 2), uniq, "域 B")
    if da == db:
        raise sc.bad_request("domain_a 与 domain_b 不能相同")
    n_rings = int(max(1, min(_p_int((params or {}).get("n_rings"), 3), 10)))
    labels_arr = np.asarray(labels, dtype=object)
    blocks = _blocks(ctx)[0]
    idx, dist = _knn(coords, _p_int((params or {}).get("n_neighbors"), 6), blocks=blocks)
    mask_a = labels_arr == da
    mask_b = labels_arr == db
    if mask_a.sum() < 3 or mask_b.sum() < 3:
        raise sc.prerequisite(f"域 {da}（{int(mask_a.sum())} spot）或域 {db}（{int(mask_b.sum())} spot）太小，无法评估边界")
    # 跨域邻接比例
    ia = mask_a[:, None]
    ib = np.zeros_like(idx, dtype=bool)
    for col in range(idx.shape[1]):
        ib[:, col] = labels_arr[idx[:, col]] == db
    aa = int(np.sum(ia & ~ib))
    ab = int(np.sum(ia & ib))
    ba = int(np.sum(~ia & ib))
    bb_mask = (labels_arr == db)[:, None]
    bb = int(np.sum(bb_mask & (np.array([[labels_arr[j] == db for j in row] for row in idx]) if False else False)))
    # 用统一方式统计 BB 边
    neigh_labels = np.empty_like(idx, dtype=object)
    for col in range(idx.shape[1]):
        neigh_labels[:, col] = labels_arr[idx[:, col]]
    a_to_a = int(np.sum((labels_arr[:, None] == da) & (neigh_labels == da)))
    a_to_b = int(np.sum((labels_arr[:, None] == da) & (neigh_labels == db)))
    b_to_b = int(np.sum((labels_arr[:, None] == db) & (neigh_labels == db)))
    total = a_to_a + a_to_b + b_to_b
    crossing = float(a_to_b / total) if total > 0 else 0.0

    # 边界 spot：A 中与 B 相邻的 spot + B 中与 A 相邻的 spot
    is_ab_neighbor = (labels_arr[:, None] == da) & (neigh_labels == db)
    is_ba_neighbor = (labels_arr[:, None] == db) & (neigh_labels == da)
    boundary_mask = is_ab_neighbor.any(axis=1) | is_ba_neighbor.any(axis=1)
    boundary_coords = coords[boundary_mask]
    if boundary_coords.shape[0] == 0:
        raise sc.prerequisite(f"域 {da} 与域 {db} 在空间上不接触（没有跨域邻接 spot），无法评估边界")

    from scipy.spatial import cKDTree

    tree = cKDTree(boundary_coords)
    dist_to_boundary, _ = tree.query(coords, k=1)

    # 表达梯度：两域各自 top 标记基因的平均 z-score
    X, genes, _spots = _matrix(ctx)
    frame_ab = _wilcoxon_frame(X, genes, mask_a, mask_b)
    top_a = frame_ab.sort_values("logFC", ascending=False)["gene"].head(20).tolist()
    top_b = frame_ab.sort_values("logFC", ascending=True)["gene"].head(20).tolist()
    mean = np.asarray(X.mean(axis=0)).ravel()
    std = np.asarray(X.std(axis=0)).ravel()
    std = np.where(std > 0, std, 1.0)
    zmat = (np.asarray(X, dtype=float) - mean) / std
    gidx = {g: i for i, g in enumerate(genes)}
    score_a = zmat[:, [gidx[g] for g in top_a if g in gidx]].mean(axis=1)
    score_b = zmat[:, [gidx[g] for g in top_b if g in gidx]].mean(axis=1)

    width = float(np.median(dist_to_boundary))
    nn_med = float(np.median(dist)) if dist.size else 1.0
    width = max(width, nn_med, 1e-6)
    ring = np.clip((dist_to_boundary / width).astype(int), 0, n_rings - 1)
    rows = []
    for side, mask, score in ((da, mask_a, score_a), (db, mask_b, score_b)):
        for r in range(n_rings):
            sel = (ring == r) & mask
            if not sel.any():
                continue
            rows.append(
                {
                    "side": str(side),
                    "ring": r + 1,
                    "n_spots": int(sel.sum()),
                    "mean_distance": round(float(np.mean(dist_to_boundary[sel])), 4),
                    "mean_marker_score": round(float(np.mean(score[sel])), 4),
                }
            )
    return {
        "labels": labels,
        "domain_a": da,
        "domain_b": db,
        "crossing": crossing,
        "a_to_a": a_to_a,
        "a_to_b": a_to_b,
        "b_to_b": b_to_b,
        "boundary_mask": boundary_mask,
        "dist_to_boundary": dist_to_boundary,
        "gradient": pd.DataFrame(rows),
        "top_a": top_a[:10],
        "top_b": top_b[:10],
        "n_rings": n_rings,
        "ring_width": width,
        "top_a_genes": top_a,
        "top_b_genes": top_b,
        "marker_frame": frame_ab,
    }


def _distance_core(ctx, params) -> dict:
    """到目标域的欧氏距离场 + 分箱表达趋势。"""
    labels, src, _info = get_domain_labels(ctx, params)
    coords = _coords(ctx)
    uniq = list(dict.fromkeys(labels))
    target = _require_label((params or {}).get("domain", 1), uniq, "目标域")
    labels_arr = np.asarray(labels, dtype=object)
    mask = labels_arr == target
    if mask.sum() < 1:
        raise sc.prerequisite(f"域 {target} 没有 spot，无法计算距离场")
    X, genes, _spots = _matrix(ctx)
    bins = int(max(2, min(_p_int((params or {}).get("bins"), 10), 50)))
    wanted = _p_list((params or {}).get("genes"))
    gidx = {g: i for i, g in enumerate(genes)}
    picked = [g for g in wanted if g in gidx]
    notes: list[str] = []
    if not picked:
        frame = _wilcoxon_frame(X, genes, mask, ~mask)
        picked = frame.sort_values("logFC", ascending=False)["gene"].head(20).tolist()
        notes.append(f"未指定 genes，自动取域 {target} 的前 20 个标记基因（Wilcoxon + BH）")
    from scipy.spatial import cKDTree

    tree = cKDTree(coords[mask])
    dist_field, _ = tree.query(coords, k=1)
    order = np.argsort(dist_field, kind="stable")
    chunks = np.array_split(order, bins)
    rows = []
    for b, chunk in enumerate(chunks):
        if chunk.size == 0:
            continue
        for gene in picked:
            col = np.asarray(X[chunk, gidx[gene]], dtype=float)
            rows.append(
                {
                    "gene": gene,
                    "bin": b + 1,
                    "mean_distance": round(float(np.mean(dist_field[chunk])), 4),
                    "mean_expr": round(float(np.mean(col)), 4),
                    "n_spots": int(chunk.size),
                }
            )
    return {
        "target": target,
        "labels": labels,
        "dist_field": dist_field,
        "mask": mask,
        "genes": picked,
        "frame": pd.DataFrame(rows),
        "bins": bins,
        "notes": notes,
        "X": X,
        "all_genes": genes,
    }


def _marker_genes_for_domain(ctx, params, domain_value, top_n: int = 20) -> tuple[pd.DataFrame, str, list[str]]:
    labels, src, _info = get_domain_labels(ctx, params)
    uniq = list(dict.fromkeys(labels))
    target = _require_label(domain_value, uniq, "域")
    X, genes, _spots = _matrix(ctx)
    frame = _marker_frame(X, genes, np.asarray(labels, dtype=object), target)
    return frame, target, labels


# ---------------------------------------------------------------------------
# 分析实现（邻域与域）
# ---------------------------------------------------------------------------

def _analysis_nhood_enrichment(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    labels, key = _cluster_labels(ctx, _p_str(params.get("cluster_key")))
    coords = _coords(ctx)
    n_perms = int(max(1, min(_p_int(params.get("n_perms"), 100), 1000)))
    n_neighbors = int(max(1, _p_int(params.get("n_neighbors"), 6)))
    res = _nhood_enrichment_core(labels, coords, n_perms, n_neighbors, _rng(ctx),
                                 blocks=_blocks(ctx)[0])
    frame = _nhood_table(res)
    csv_path = sc.export_csv(frame, outdir, "nhood_enrichment", index=False)
    table_obj = sc.table(
        "nhood_enrichment",
        "邻域富集（类别对 z-score）",
        [
            {"key": "cluster_a", "label": "类别 A", "type": "string"},
            {"key": "cluster_b", "label": "类别 B", "type": "string"},
            {"key": "observed", "label": "观测邻接数", "type": "number"},
            {"key": "null_mean", "label": "置换均值", "type": "number"},
            {"key": "z", "label": "z-score", "type": "number"},
        ],
        frame,
        note=f"{n_perms} 次标签置换；z>0 表示该类别对在空间上比随机更常相邻",
    )
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_nhood_figure(ctx, params, outdir, res=res))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"邻域富集热图绘制失败：{exc}")
    top = frame.head(5)
    bullets = [
        f"按「{key}」的 {len(res['uniq'])} 个类别 × {res['n_neighbors']} 近邻统计邻接；{n_perms} 次置换得到 z-score",
        "最富集的类别对：" + "、".join(f"{r['cluster_a']}-{r['cluster_b']}(z={r['z']:+.2f})" for _, r in top.iterrows()),
    ]
    return sc.result(
        title="邻域富集",
        bullets=bullets,
        metrics={
            "cluster_key": key,
            "n_categories": int(len(res["uniq"])),
            "n_perms": n_perms,
            "n_neighbors": res["n_neighbors"],
            "max_z": float(np.max(res["z"])) if res["z"].size else None,
            "categories": [str(u) for u in res["uniq"]],
        },
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=[
            "方法：类别对邻接计数 + 标签置换 z-score（squidpy nhood_enrichment 思想；Schürch 2020 邻域分析；内置实现）",
            "置换保持类别大小不变、打乱空间位置，衡量「同类相邻」是否超出随机预期",
        ],
        data_out=[sc.data_file("nhood_csv", csv_path, "邻域富集表")],
    )


def _analysis_co_occurrence(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    labels, key = _cluster_labels(ctx, _p_str(params.get("cluster_key")))
    interval = params.get("interval")
    n_rings = int(max(1, min(_p_int(params.get("n_rings"), 5), 20)))
    res = _co_occurrence_core(
        labels, _coords(ctx), n_rings=n_rings,
        interval=interval if interval not in (None, "", []) else None,
        blocks=_blocks(ctx)[0],
    )
    frame = _cooccurrence_table(res)
    csv_path = sc.export_csv(frame, outdir, "co_occurrence", index=False)
    table_obj = sc.table(
        "co_occurrence",
        "共现概率（按距离环）",
        [
            {"key": "cluster_a", "label": "锚定类别", "type": "string"},
            {"key": "cluster_b", "label": "邻居类别", "type": "string"},
            {"key": "ring", "label": "环", "type": "number"},
            {"key": "distance", "label": "距离", "type": "number"},
            {"key": "prob", "label": "条件概率", "type": "number"},
            {"key": "n_pairs", "label": "邻居对数", "type": "number"},
        ],
        frame,
        note="P(邻居类别 | 锚定类别, 距离环)；距离环由 5%–95% 邻居距离等分",
    )
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_cooccurrence_figure(ctx, params, outdir, res=res))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"共现曲线绘制失败：{exc}")
    lo = float(res["edges"][0])
    hi = float(res["edges"][-1])
    return sc.result(
        title="共现分析",
        bullets=[
            f"按「{key}」的 {len(res['uniq'])} 个类别、{n_rings} 个距离环（{lo:.2f}–{hi:.2f}）统计条件共现概率",
            f"共 {len(frame)} 条「锚定类别 × 邻居类别 × 环」记录；完整表见 co_occurrence.csv",
        ],
        metrics={
            "cluster_key": key,
            "n_rings": n_rings,
            "distance_range": [round(lo, 4), round(hi, 4)],
            "n_categories": int(len(res["uniq"])),
            "max_prob": float(np.max(res["prob"])) if res["prob"].size else None,
        },
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=[
            "方法：环形分箱的条件共现概率（squidpy co_occurrence 思想；内置实现）",
            "环数为距离分箱数；概率按锚定类别的邻居总数归一化，取值 ∈[0,1]",
        ],
        data_out=[sc.data_file("cooccurrence_csv", csv_path, "共现概率表")],
    )


def _analysis_domains_mrf(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    X, genes, spots = _matrix(ctx)
    computed = compute_domains(ctx, params, outdir, write_back=False)
    labels = computed["labels"]
    clusters, ckey = _cluster_labels(ctx, _p_str(params.get("cluster_key")))
    try:
        from sklearn.metrics import adjusted_rand_score

        ari = float(adjusted_rand_score([str(c) for c in clusters], labels))
    except Exception:  # noqa: BLE001
        ari = float("nan")
    frame = pd.DataFrame(
        [
            {"domain": str(d), "n_spots": int(sizes), "pct": round(sizes / len(labels) * 100.0, 2)}
            for d, sizes in computed["sizes"].items()
        ]
    ).sort_values("n_spots", ascending=False)
    crosstab = pd.crosstab(
        pd.Series([str(c) for c in clusters], name="cluster"),
        pd.Series(labels, name="domain"),
    ).reset_index()
    tables = [
        sc.table(
            "domains_mrf",
            "空间域大小",
            [
                {"key": "domain", "label": "域", "type": "string"},
                {"key": "n_spots", "label": "spot 数", "type": "number"},
                {"key": "pct", "label": "占比(%)", "type": "number"},
            ],
            frame,
            note=f"与「{ckey}」聚类的一致性 ARI = {ari:.3f}（分析不写回数据集）",
        ),
        sc.table("domains_crosstab", "聚类 × 空间域 交叉表", [], crosstab),
    ]
    figures: list[dict] = []
    warnings = list(computed["warnings"])
    try:
        figures.append(_domain_map_figure(ctx, params, outdir, labels=labels,
                                          title="Spatial domains (MRF smooth)"))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"空间域图绘制失败：{exc}")
    csv_path = sc.export_csv(
        pd.DataFrame({"spot": [str(s) for s in spots], "domain": labels}), outdir, "domains_mrf", index=False
    )
    return sc.result(
        title="空间域（MRF 平滑）",
        bullets=[
            f"MRF 平滑得到 {computed['n_domains']} 个域（method={computed['method']}，"
            f"PCA {computed['pca_dims']} 维，迭代 {computed['n_iter']} 轮）",
            f"与现有聚类「{ckey}」的调整兰德指数 ARI = {ari:.3f}"
            + ("（一致性较好）" if ari >= 0.5 else "（两者划分不同，可能表达域与聚类粒度不同）"),
        ],
        metrics={
            "method": computed["method"],
            "n_domains": int(computed["n_domains"]),
            "iterations_used": int(computed["n_iter"]),
            "pca_dims": int(computed["pca_dims"]),
            "ari_vs_clusters": None if not np.isfinite(ari) else round(ari, 4),
            "sizes": computed["sizes"],
        },
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=list(computed["notes"]) + ["本分析只返回结果，不修改数据集的 domains 列（如需写回请运行「空间域」步骤）"],
        data_out=[sc.data_file("domains_csv", csv_path, "spot × 空间域（分析）")],
    )


def _analysis_domain_markers(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    top_n = int(max(1, min(_p_int(params.get("top_n"), 20), 500)))
    frame, target, labels = _marker_genes_for_domain(ctx, params, params.get("domain", 1))
    head = _de_head(frame, top_n)
    csv_path = sc.export_csv(frame, outdir, "domain_markers", index=False)
    table_obj = sc.table(
        "domain_markers",
        f"域 {target} vs 其余：标记基因",
        _de_columns(),
        head,
        note=f"Wilcoxon 秩和 + BH FDR；完整表见 domain_markers.csv（p_adj<0.05 视为显著）",
    )
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_de_heatmap_figure(ctx, params, outdir, frame=frame, labels=labels,
                                          group_title=f"Domain {target}", fig_name="de_heatmap"))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"标记基因热图绘制失败：{exc}")
    try:
        figures.append(_dotplot_figure(ctx, params, outdir, frame=frame, labels=labels,
                                       focus=target, fig_name="domain_dotplot"))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"标记基因点图绘制失败：{exc}")
    top_genes = "、".join(f"{row['gene']}({row['logFC']:+.2f})" for _, row in head.head(5).iterrows())
    return sc.result(
        title="区域标记基因",
        bullets=[
            f"域 {target}（{int(np.sum(np.asarray(labels, dtype=object) == target))} spot）vs 其余，"
            f"Wilcoxon + BH：显著 {int((frame['p_adj'] <= 0.05).sum())} / {len(frame)} 个基因",
            f"最显著的域 {target} 富集基因：{top_genes}",
        ],
        metrics={
            "domain": target,
            "n_sig_up": int(((frame["p_adj"] <= 0.05) & (frame["logFC"] > 0)).sum()),
            "n_sig_down": int(((frame["p_adj"] <= 0.05) & (frame["logFC"] < 0)).sum()),
            "n_genes": int(len(frame)),
            "top_n": top_n,
        },
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=[
            "方法：Wilcoxon 秩和检验（scipy mannwhitneyu，双侧）+ BH FDR；logFC 为两组均值差（内置实现）",
            "域标签来自数据集 domains 列；没有时按默认参数现场计算",
        ],
        data_out=[sc.data_file("domain_markers_csv", csv_path, "域标记基因全表")],
    )


def _analysis_domain_composition(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    labels, src, _info = get_domain_labels(ctx, params)
    n = len(labels)
    group_col = _p_str(params.get("group_col"))
    groups, gkey = _group_values(ctx, n, group_col)
    warnings: list[str] = []
    if groups is None:
        groups = ["all"] * n
        gkey = f"{group_col or '无分组'}"
        warnings.append("数据集没有切片/分组列，只能输出整体域组成（可用 group_col 指定分组列）")
    rows = []
    for g in list(dict.fromkeys(groups)):
        sel = np.asarray(groups, dtype=object) == g
        if not sel.any():
            continue
        for d in list(dict.fromkeys(labels)):
            m = sel & (np.asarray(labels, dtype=object) == d)
            rows.append(
                {
                    "group": str(g),
                    "domain": str(d),
                    "n_spots": int(m.sum()),
                    "pct": round(float(m.sum() / max(1, sel.sum()) * 100.0), 2),
                }
            )
    frame = pd.DataFrame(rows)
    csv_path = sc.export_csv(frame, outdir, "domain_composition", index=False)
    table_obj = sc.table(
        "domain_composition",
        f"各组在各域的 spot 比例（分组列：{gkey}）",
        [
            {"key": "group", "label": "分组", "type": "string"},
            {"key": "domain", "label": "域", "type": "string"},
            {"key": "n_spots", "label": "spot 数", "type": "number"},
            {"key": "pct", "label": "组内占比(%)", "type": "number"},
        ],
        frame,
        note="pct = 该组内落在该域的 spot 百分比（每个组内合计 100%）",
    )
    figures: list[dict] = []
    try:
        figures.append(_composition_bar_figure(ctx, params, outdir, frame=frame, gkey=gkey))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"域组成堆叠图绘制失败：{exc}")
    dominant = (
        frame.sort_values("pct", ascending=False).groupby("group").head(1)
        if not frame.empty
        else frame
    )
    bullets = [
        f"按「{gkey}」统计 {frame['group'].nunique()} 个分组 × {frame['domain'].nunique()} 个域的 spot 比例",
        "各组占比最高的域：" + "；".join(
            f"{row['group']}→域{row['domain']}({row['pct']:.1f}%)" for _, row in dominant.iterrows()
        ),
    ]
    return sc.result(
        title="区域组成",
        bullets=bullets,
        metrics={
            "group_col": gkey,
            "n_groups": int(frame["group"].nunique()) if not frame.empty else 0,
            "n_domains": int(frame["domain"].nunique()) if not frame.empty else 0,
        },
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=[
            "方法：按分组/切片统计各域 spot 比例（组成分析，内置实现）",
            "域标签来自数据集 domains 列；没有时按默认参数现场计算",
        ],
        data_out=[sc.data_file("domain_composition_csv", csv_path, "域组成表")],
    )


def _analysis_boundary_score(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    res = _boundary_core(ctx, params)
    gradient = res["gradient"]
    summary = pd.DataFrame(
        [
            {"metric": "cross_domain_pairs", "value": res["a_to_b"], "label": "跨域邻接对（A→B）"},
            {"metric": "within_a_pairs", "value": res["a_to_a"], "label": "域 A 内部邻接对"},
            {"metric": "within_b_pairs", "value": res["b_to_b"], "label": "域 B 内部邻接对"},
            {"metric": "crossing_ratio", "value": round(res["crossing"], 4), "label": "跨域邻接比例"},
            {"metric": "n_boundary_spots", "value": int(np.sum(res["boundary_mask"])), "label": "边界 spot 数"},
            {"metric": "n_rings", "value": res["n_rings"], "label": "环形梯度层数"},
        ]
    )
    tables = [
        sc.table(
            "boundary_summary",
            f"域 {res['domain_a']} / 域 {res['domain_b']} 边界统计",
            [
                {"key": "label", "label": "指标", "type": "string"},
                {"key": "value", "label": "数值", "type": "number"},
            ],
            summary,
            note="跨域邻接比例 = A↔B 邻接对 / (A内部 + B内部 + A↔B)",
        ),
        sc.table(
            "boundary_gradient",
            "边界环形表达梯度",
            [
                {"key": "side", "label": "域", "type": "string"},
                {"key": "ring", "label": "环", "type": "number"},
                {"key": "n_spots", "label": "spot 数", "type": "number"},
                {"key": "mean_distance", "label": "平均到边界距离", "type": "number"},
                {"key": "mean_marker_score", "label": "标记基因平均 z", "type": "number"},
            ],
            gradient,
            note="环 1 = 离边界最近；标记基因取两域差异最显著的前 20 个",
        ),
    ]
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_boundary_figure(ctx, params, outdir, res=res))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"边界图绘制失败：{exc}")
    g_csv = sc.export_csv(gradient, outdir, "boundary_gradient", index=False)
    metrics = {
        "domain_a": res["domain_a"],
        "domain_b": res["domain_b"],
        "crossing_ratio": round(res["crossing"], 4),
        "n_boundary_spots": int(np.sum(res["boundary_mask"])),
        "n_rings": res["n_rings"],
        "ring_width": round(res["ring_width"], 4),
    }
    return sc.result(
        title="边界强度",
        bullets=[
            f"域 {res['domain_a']} 与域 {res['domain_b']}：跨域邻接比例 {res['crossing'] * 100:.1f}%"
            f"（A→B {res['a_to_b']} 对，A 内部 {res['a_to_a']}，B 内部 {res['b_to_b']}）",
            f"边界 spot {int(np.sum(res['boundary_mask']))} 个；{res['n_rings']} 层环（环宽 ≈ {res['ring_width']:.1f}）",
        ],
        metrics=metrics,
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=[
            "方法：邻接对跨域比例 + 到边界的环形标记基因 z-score 梯度（内置实现）",
            "标记基因由域 A vs 域 B 的 Wilcoxon 检验给出；环宽取边界 spot 到全体 spot 的中位距离",
        ],
        data_out=[sc.data_file("boundary_csv", g_csv, "边界梯度表")],
    )


def _analysis_distance_to_domain(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    res = _distance_core(ctx, params)
    frame = res["frame"]
    csv_path = sc.export_csv(frame, outdir, "distance_to_domain", index=False)
    table_obj = sc.table(
        "distance_to_domain",
        f"到域 {res['target']} 的距离 × 表达",
        [
            {"key": "gene", "label": "基因", "type": "string"},
            {"key": "bin", "label": "距离分箱", "type": "number"},
            {"key": "mean_distance", "label": "平均距离", "type": "number"},
            {"key": "mean_expr", "label": "平均表达", "type": "number"},
            {"key": "n_spots", "label": "spot 数", "type": "number"},
        ],
        frame,
        note=f"{res['bins']} 个等频距离分箱；距离取到最近域 {res['target']} spot 的欧氏距离",
    )
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_distance_figure(ctx, params, outdir, res=res))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"距离趋势图绘制失败：{exc}")
    # 距离-表达相关性（用于摘要）
    corr_rows = []
    for gene in res["genes"][:20]:
        sub = frame[frame["gene"] == gene]
        if len(sub) < 3:
            continue
        d = sub["mean_distance"].to_numpy(dtype=float)
        e = sub["mean_expr"].to_numpy(dtype=float)
        if np.std(d) > 0 and np.std(e) > 0:
            corr_rows.append({"gene": gene, "pearson_r": float(np.corrcoef(d, e)[0, 1])})
    bullets = [
        f"目标域 {res['target']}：{int(res['mask'].sum())} 个 spot；距离范围 0–{res['dist_field'].max():.1f}",
        f"分析基因 {len(res['genes'])} 个（{', '.join(res['genes'][:6])}…），按 {res['bins']} 个距离分箱统计平均表达",
    ]
    if corr_rows:
        best = max(corr_rows, key=lambda r: abs(r["pearson_r"]))
        bullets.append(f"距离-表达相关最强：{best['gene']}（r={best['pearson_r']:+.2f}）")
    metrics = {
        "domain": res["target"],
        "n_target_spots": int(res["mask"].sum()),
        "bins": int(res["bins"]),
        "max_distance": round(float(np.max(res["dist_field"])), 4),
        "correlations": {r["gene"]: round(r["pearson_r"], 4) for r in corr_rows},
    }
    return sc.result(
        title="到区域距离",
        bullets=bullets,
        metrics=metrics,
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=list(res["notes"])
        + [
            "方法：欧氏距离场（scipy cKDTree 最近目标域 spot）+ 等频距离分箱的平均表达（内置实现）",
            "距离为 0 表示 spot 本身属于目标域",
        ],
        data_out=[sc.data_file("distance_csv", csv_path, "距离分箱表达表")],
    )


def _analysis_niche_cluster(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    labels, key = _cluster_labels(ctx, _p_str(params.get("cluster_key")))
    coords = _coords(ctx)
    k = int(max(2, min(_p_int(params.get("k"), 8), 30)))
    n_niches = int(max(2, min(_p_int(params.get("n_niches"), 5), 20)))
    comp, uniq = _neighborhood_composition(labels, coords, k, blocks=_blocks(ctx)[0])
    if len(uniq) < 2:
        raise sc.prerequisite(f"聚类类别只有 {len(uniq)} 个，无法构建邻域组成向量（至少需要 2 类）")
    n_niches = min(n_niches, max(2, comp.shape[0] // 10))
    from sklearn.cluster import KMeans

    model = KMeans(n_clusters=n_niches, n_init=10, random_state=0)
    niche = model.fit_predict(comp)
    niche_labels = [str(int(v) + 1) for v in niche]
    # 按规模降序重编号
    uniq_n, counts = np.unique(niche, return_counts=True)
    order = np.argsort(-counts)
    remap = {int(uniq_n[i]): str(rank + 1) for rank, i in enumerate(order)}
    niche_labels = [remap[int(v)] for v in niche]

    comp_frame = _composition_table(niche_labels, labels)
    size_frame = pd.DataFrame(
        [
            {"niche": str(n), "n_spots": int(v)}
            for n, v in pd.Series(niche_labels).value_counts().sort_index().items()
        ]
    )
    tables = [
        sc.table(
            "niche_sizes",
            "niche 大小",
            [
                {"key": "niche", "label": "niche", "type": "string"},
                {"key": "n_spots", "label": "spot 数", "type": "number"},
            ],
            size_frame,
            note=f"{k} 近邻组成向量 → k-means 聚成 {n_niches} 个 niche（按规模降序编号）",
        ),
        sc.table(
            "niche_composition",
            "niche 的邻域类别组成",
            [{"key": "label", "label": "niche", "type": "string"},
             {"key": "n_spots", "label": "spot 数", "type": "number"}]
            + [{"key": str(c), "label": f"类别 {c}", "type": "number"} for c in uniq],
            comp_frame,
            note="每行为该 niche 的平均邻域类别比例（行内合计≈1）",
        ),
    ]
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_niche_figure(ctx, params, outdir, niche_labels=niche_labels))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"niche 空间图绘制失败：{exc}")
    csv_path = sc.export_csv(
        pd.DataFrame({"spot": [str(s) for s in (ctx.get("spots") or [])], "niche": niche_labels}),
        outdir,
        "niche",
        index=False,
    )
    biggest = size_frame.sort_values("n_spots", ascending=False).iloc[0]
    bullets = [
        f"按「{key}」的 {k} 近邻组成向量聚出 {n_niches} 个 niche；最大 niche {biggest['niche']}"
        f"（{int(biggest['n_spots'])} spot）",
        "各 niche 主导类别："
        + "；".join(
            f"niche{row['label']}→{max(uniq, key=lambda c: float(row.get(str(c), 0.0)))}"
            for _, row in comp_frame.iterrows()
        ),
    ]
    out = sc.result(
        title="细胞邻域（niche）",
        bullets=bullets,
        metrics={
            "cluster_key": key,
            "k_neighbors": int(k),
            "n_niches": int(n_niches),
            "n_categories": int(len(uniq)),
            "sizes": {str(r["niche"]): int(r["n_spots"]) for _, r in size_frame.iterrows()},
        },
        tables=tables,
        figures=figures,
        warnings=warnings,
        notes=[
            "方法：k 近邻类别组成向量 + k-means（Schürch 2020 / Bhate 2022 邻域思想；内置实现）",
            "niche 标签已写回 spot 级列 niche，可在后续分析/导出中使用",
        ],
        data_out=[sc.data_file("niche_csv", csv_path, "spot × niche")],
        sample_columns=[{"key": "niche", "kind": "categorical"}],
        sample_values={"niche": niche_labels},
    )
    return out


def _analysis_tissue_diversity(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    labels, key = _cluster_labels(ctx, _p_str(params.get("cluster_key")))
    coords = _coords(ctx)
    k = int(max(2, min(_p_int(params.get("k"), 6), 30)))
    res = _diversity_core(labels, coords, k, blocks=_blocks(ctx)[0])
    shannon = res["shannon"]
    simpson = res["simpson"]
    labels_arr = res["labels"]
    rows = []
    for name in list(dict.fromkeys(labels_arr.tolist())):
        mask = labels_arr == name
        if not mask.any():
            continue
        rows.append(
            {
                "cluster": str(name),
                "n_spots": int(mask.sum()),
                "shannon_mean": round(float(np.mean(shannon[mask])), 4),
                "shannon_median": round(float(np.median(shannon[mask])), 4),
                "simpson_mean": round(float(np.mean(simpson[mask])), 4),
            }
        )
    frame = pd.DataFrame(rows).sort_values("shannon_mean", ascending=False)
    frame.loc[len(frame)] = {
        "cluster": "ALL",
        "n_spots": int(labels_arr.size),
        "shannon_mean": round(float(np.mean(shannon)), 4),
        "shannon_median": round(float(np.median(shannon)), 4),
        "simpson_mean": round(float(np.mean(simpson)), 4),
    }
    table_obj = sc.table(
        "tissue_diversity",
        f"邻域多样性（{k} 近邻）",
        [
            {"key": "cluster", "label": "类别", "type": "string"},
            {"key": "n_spots", "label": "spot 数", "type": "number"},
            {"key": "shannon_mean", "label": "Shannon 均值", "type": "number"},
            {"key": "shannon_median", "label": "Shannon 中位数", "type": "number"},
            {"key": "simpson_mean", "label": "Simpson 均值", "type": "number"},
        ],
        frame,
        note="Shannon = -Σp·ln p；Simpson = 1-Σp²（p 为邻域类别比例）；ALL 行为全体 spot",
    )
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_diversity_figure(ctx, params, outdir, shannon=shannon, simpson=simpson))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"多样性空间图绘制失败：{exc}")
    csv_path = sc.export_csv(frame, outdir, "tissue_diversity", index=False)
    top = frame.iloc[0]
    return sc.result(
        title="邻域多样性",
        bullets=[
            f"按「{key}」的 {k} 近邻组成计算 Shannon/Simpson 多样性；全体 spot 平均 Shannon = {np.mean(shannon):.3f}",
            f"邻域最多样的类别：{top['cluster']}（Shannon {top['shannon_mean']:.3f}）",
        ],
        metrics={
            "cluster_key": key,
            "k_neighbors": int(res["k"]),
            "n_categories": int(len(res["uniq"])),
            "shannon_mean": round(float(np.mean(shannon)), 4),
            "simpson_mean": round(float(np.mean(simpson)), 4),
            "shannon_range": [round(float(np.min(shannon)), 4), round(float(np.max(shannon)), 4)],
        },
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=[
            "方法：邻域类别组成向量的 Shannon / Simpson 多样性（内置实现）",
            "多样性高 = 该 spot 周围混合了更多类别（过渡带/免疫浸润带）；低 = 均一区域",
        ],
        data_out=[sc.data_file("diversity_csv", csv_path, "多样性表")],
    )


# ---------------------------------------------------------------------------
# 分析实现（区域差异）
# ---------------------------------------------------------------------------

def _de_result(ctx, params, outdir, *, frame, title, case_label, ctrl_label, extra_figures=None,
               warnings=None, fig_name="de_volcano", metrics_extra=None, bullet_head=None,
               notes_extra=None, table_note=None, notes_base=None) -> dict:
    top_n = int(max(1, min(_p_int(params.get("top_n"), 200), 2000)))
    head = _de_head(frame, top_n)
    csv_path = sc.export_csv(frame, outdir, f"{fig_name}_all", index=False)
    note = table_note or f"{case_label} vs {ctrl_label}；Wilcoxon 秩和 + BH FDR"
    table_obj = sc.table(
        fig_name,
        f"{title}（前 {min(top_n, len(head))} 行）",
        _de_columns(),
        head,
        note=f"{note}；完整表见 {fig_name}_all.csv",
    )
    figures: list[dict] = list(extra_figures or [])
    warns = list(warnings or [])
    labels = None
    try:
        figures.append(_de_volcano_figure(ctx, params, outdir, frame=frame, title=title, fig_name=fig_name))
    except Exception as exc:  # noqa: BLE001
        warns.append(f"火山图绘制失败：{exc}")
    n_up = int(((frame["p_adj"] <= 0.05) & (frame["logFC"] > 0)).sum())
    n_down = int(((frame["p_adj"] <= 0.05) & (frame["logFC"] < 0)).sum())
    n_sig = n_up + n_down
    bullets = [
        bullet_head or (
            f"{case_label} vs {ctrl_label}（Wilcoxon + BH）：显著 {n_sig} / {len(frame)} 个基因"
            f"（上调 {n_up} / 下调 {n_down}）"),
        "最显著基因：" + "、".join(
            f"{row['gene']}({row['logFC']:+.2f})" for _, row in head.head(5).iterrows()
        ),
    ]
    metrics = {
        "case": case_label,
        "control": ctrl_label,
        "n_genes": int(len(frame)),
        "n_sig": n_sig,
        "n_up": n_up,
        "n_down": n_down,
        "top_n": top_n,
    }
    metrics.update(metrics_extra or {})
    if notes_base is not None:
        notes = list(notes_base) + list(notes_extra or [])
    else:
        notes = [
            "方法：Wilcoxon 秩和检验（scipy mannwhitneyu，双侧）+ BH FDR；logFC = 两组均值差（内置实现）",
            "pct_case / pct_ctrl 为表达 >0 的 spot 百分比",
        ] + list(notes_extra or [])
    return sc.result(
        title=title,
        bullets=bullets,
        metrics=metrics,
        tables=[table_obj],
        figures=figures,
        warnings=warns,
        notes=notes,
        data_out=[sc.data_file(f"{fig_name}_csv", csv_path, "差异全表")],
    )


def _analysis_de_domain(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    frame, target, labels = _marker_genes_for_domain(ctx, params, params.get("domain", 1))
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_de_heatmap_figure(ctx, params, outdir, frame=frame, labels=labels,
                                          group_title=f"Domain {target}", fig_name="de_heatmap"))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"差异热图绘制失败：{exc}")
    return _de_result(
        ctx, params, outdir, frame=frame, title=f"域 {target} vs 其余 差异",
        case_label=f"域 {target}", ctrl_label="其余", extra_figures=figures, warnings=warnings,
    )


def _analysis_de_regions(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    labels, src, _info = get_domain_labels(ctx, params)
    uniq = list(dict.fromkeys(labels))
    da = _require_label(params.get("domain_a", 1), uniq, "域 A")
    db = _require_label(params.get("domain_b", 2), uniq, "域 B")
    if da == db:
        raise sc.bad_request("domain_a 与 domain_b 不能相同")
    X, genes, _spots = _matrix(ctx)
    labels_arr = np.asarray(labels, dtype=object)
    frame = _wilcoxon_frame(X, genes, labels_arr == da, labels_arr == db)
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_de_heatmap_figure(ctx, params, outdir, frame=frame, labels=labels,
                                          group_title=f"Domain {da} vs {db}", fig_name="de_heatmap",
                                          keep_labels=[da, db]))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"差异热图绘制失败：{exc}")
    return _de_result(
        ctx, params, outdir, frame=frame, title=f"域 {da} vs 域 {db} 差异",
        case_label=f"域 {da}", ctrl_label=f"域 {db}", extra_figures=figures, warnings=warnings,
    )


#: 区域列候选（de_groups 的「供体/切片 × 区域」单位里，区域优先用空间域，其次聚类）
_REGION_COLUMNS = ("domains", "domain", "spatial_domain", "domain_id", "region",
                   "clusters", "cluster")


def _region_labels(ctx, n: int, preferred: str = ""):
    """区域标签（domains 优先，其次 clusters）：返回 (labels|None, 命中的列名)。"""
    frame = ctx.get("sample_info")
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        return None, ""
    for cand in ([preferred] if preferred else []) + list(_REGION_COLUMNS):
        if not cand or cand not in frame.columns:
            continue
        labels = _column_labels(ctx, str(cand), n)
        if labels is None:
            continue
        valid = [v for v in labels if v not in ("NA", "")]
        if len(valid) >= 1:
            return labels, str(cand)
    return None, ""


def _unit_table(ctx, n: int, group_col: str = ""):
    """构造 供体/切片（物理块） × 区域 的观测单位表（不做统计，只分组）。

    返回 ``dict``：units=[{block, region, group, n_spots, indices}]、block_key、region_key、
    group_key、warnings；无分组时返回 None。
    """
    groups, gkey = _group_values(ctx, n, group_col)
    if groups is None:
        return None
    blocks, block_names, block_key, block_warns = _blocks(ctx)
    block_id = np.zeros(n, dtype=int)
    for i, idx in enumerate(blocks):
        block_id[np.asarray(idx, dtype=int)] = i
    region_labels, region_key = _region_labels(ctx, n, _p_str((ctx or {}).get("region_col")))
    if region_labels is None:
        region_labels = ["all"] * n
        region_key = "all"
    grouped: dict[tuple[int, str], list[int]] = {}
    for i in range(n):
        grouped.setdefault((int(block_id[i]), str(region_labels[i])), []).append(i)
    units: list[dict] = []
    warnings: list[str] = list(block_warns or [])
    mixed: list[str] = []
    for (bi, region), idx in sorted(grouped.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        idx_arr = np.asarray(idx, dtype=int)
        unit_groups = [str(groups[i]) for i in idx_arr]
        uniq_g = list(dict.fromkeys(unit_groups))
        if len(uniq_g) > 1:
            mixed.append(f"{block_names[bi]}/{region}")
        units.append({"block": str(block_names[bi]), "region": str(region),
                      "group": uniq_g[0] if uniq_g else "NA",
                      "n_spots": int(idx_arr.size), "indices": idx_arr})
    if mixed:
        warnings.append(
            "以下 切片×区域 单位里混有多个分组取值（切片列与分组列可能交叉）："
            + "、".join(mixed[:6]) + "；单位级检验按该单位首个 spot 的分组归类，结果需谨慎")
    return {"units": units, "block_key": block_key, "region_key": region_key,
            "group_key": gkey, "warnings": warnings}


def _unit_expression(ctx, units: list[dict], genes: list[str],
                     min_unit_spots: int = 5) -> tuple[np.ndarray, list[dict], list[str]]:
    """把每个 供体/切片×区域 单位聚合成一行表达（原始 counts 求和 → log1p CPM）。

    * 优先用 ``ctx["counts"]``（原始 counts）：单位内按基因求和 → log1p(1e6·count/库大小)，
      这就是 pseudobulk（单位级）表达，供体/切片是生物学重复；
    * 没有 counts 时退回对当前矩阵求和并写 warning（不是真正的 pseudobulk counts）。
    """
    X, _genes, _spots = _matrix(ctx)
    counts = ctx.get("counts")
    use_counts = counts is not None and np.asarray(counts).shape == X.shape
    source = np.asarray(counts if use_counts else X, dtype=np.float64)
    warnings: list[str] = []
    if not use_counts:
        warnings.append(
            "数据集没有原始 counts：单位聚合改用当前表达矩阵求和（不是真正的 pseudobulk counts），"
            "组间比较的量级需谨慎解读")
    kept: list[dict] = []
    rows: list[np.ndarray] = []
    dropped = 0
    for unit in units:
        idx = np.asarray(unit["indices"], dtype=int)
        if idx.size < int(max(1, min_unit_spots)):
            dropped += 1
            continue
        agg = np.nansum(source[idx], axis=0)
        lib = float(np.nansum(agg))
        if lib <= 0:
            dropped += 1
            continue
        rows.append(np.log1p(1e6 * agg / lib))
        kept.append(unit)
    if dropped:
        warnings.append(f"{dropped} 个单位 spot 数 < {min_unit_spots} 或总 counts 为 0，已从单位级检验中剔除")
    if not rows:
        raise sc.prerequisite(
            f"没有 spot 数 ≥ {min_unit_spots} 的 供体/切片×区域 单位：请降低 min_unit_spots，或检查切片/区域列")
    return np.vstack(rows), kept, warnings


def _unit_de_frame(U: np.ndarray, genes: list[str], case_mask: np.ndarray, ctrl_mask: np.ndarray,
                   *, test: bool = True) -> pd.DataFrame:
    """单位级差异表：每行一个基因。

    test=True：单位间 Wilcoxon 秩和（scipy mannwhitneyu，method=auto）+ BH；
    test=False：只给描述性汇总（均值/logFC/阳性比例），p 值与 FDR 留空（NaN）。
    """
    A = np.asarray(U[case_mask], dtype=float)
    B = np.asarray(U[ctrl_mask], dtype=float)
    n_genes = len(genes)
    mean_case = A.mean(axis=0) if A.size else np.full(n_genes, np.nan)
    mean_ctrl = B.mean(axis=0) if B.size else np.full(n_genes, np.nan)
    pct_case = (A > 0).mean(axis=0) * 100.0 if A.size else np.full(n_genes, np.nan)
    pct_ctrl = (B > 0).mean(axis=0) * 100.0 if B.size else np.full(n_genes, np.nan)
    if test and A.shape[0] >= 2 and B.shape[0] >= 2:
        from scipy import stats

        pvals = np.ones(n_genes, dtype=float)
        for j in range(n_genes):
            try:
                _stat, p = stats.mannwhitneyu(A[:, j], B[:, j], alternative="two-sided",
                                              method="auto")
                pvals[j] = float(p) if np.isfinite(p) else 1.0
            except Exception:  # noqa: BLE001 - 单基因退化时按不显著处理
                pvals[j] = 1.0
        pvals = np.clip(pvals, 0.0, 1.0)
        frame = pd.DataFrame({
            "gene": genes,
            "logFC": mean_case - mean_ctrl,
            "mean_case": mean_case,
            "mean_ctrl": mean_ctrl,
            "pct_case": pct_case,
            "pct_ctrl": pct_ctrl,
            "p_value": pvals,
            "p_adj": sc.bh_fdr(pvals),
        })
    else:
        frame = pd.DataFrame({
            "gene": genes,
            "logFC": mean_case - mean_ctrl,
            "mean_case": mean_case,
            "mean_ctrl": mean_ctrl,
            "pct_case": pct_case,
            "pct_ctrl": pct_ctrl,
            "p_value": np.nan,
            "p_adj": np.nan,
        })
    frame["regulation"] = np.where(
        (frame["p_adj"] <= 0.05) & (frame["logFC"] > 0), "up",
        np.where((frame["p_adj"] <= 0.05) & (frame["logFC"] < 0), "down", "ns"))
    return frame


def _analysis_de_groups(ctx, params=None, outdir=None):
    """按分组差异（单位级）：以 供体/切片 × 区域 为观测单位做 pseudobulk 比较。

    为什么不能再按 spot 做 Wilcoxon：同一张切片上的 spot 不是独立生物学重复，切片内 spot 数
    再多也只代表 1 个供体/切片 —— 这就是空间伪重复（pseudoreplication），会把 p 值压到 1e-9
    量级。正确做法：每个 供体/切片×区域 单位先把原始 counts 聚合成 pseudobulk（log1p CPM），
    单位之间再做 Wilcoxon + BH；每组只有 1 个供体/切片时只输出**描述性**（descriptive）汇总，
    **不做推断**（inferential）。
    """
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    X, genes, spots = _matrix(ctx)
    n = X.shape[0]
    group_col = _p_str(params.get("group_col"))
    groups, gkey = _group_values(ctx, n, group_col)
    if groups is None:
        raise sc.prerequisite(
            "没有 spot 级分组列（如切片/分组）可比较；请先在「样本」里导入分组表，或用 group_col 指定"
        )
    uniq = [g for g in dict.fromkeys(groups) if g not in ("NA", "")]
    if len(uniq) < 2:
        raise sc.prerequisite(f"分组列「{gkey}」只有 {len(uniq)} 个有效取值，无法做组间差异")
    case = _p_str(params.get("case"))
    control = _p_str(params.get("control"))
    if case and case not in uniq:
        raise sc.prerequisite(f"实验组「{case}」不在分组列「{gkey}」里（可选：{'、'.join(uniq)}）")
    if control and control not in uniq:
        raise sc.prerequisite(f"对照组「{control}」不在分组列「{gkey}」里（可选：{'、'.join(uniq)}）")
    if not case and not control:
        if len(uniq) > 2:
            raise sc.prerequisite(f"分组列「{gkey}」有 {len(uniq)} 个取值（{'、'.join(uniq)}），请指定 case 与 control")
        case, control = uniq[0], uniq[1]
    elif not case:
        case = [g for g in uniq if g != control][0]
    elif not control:
        control = [g for g in uniq if g != case][0]
    if case == control:
        raise sc.bad_request("case 与 control 相同，请选择不同的两组")

    min_unit_spots = int(max(2, min(_p_int(params.get("min_unit_spots"), 5), 1000)))
    unit_info = _unit_table(ctx, n, group_col)
    units = list(unit_info["units"])
    block_key = unit_info["block_key"] or "单块"
    region_key = unit_info["region_key"]
    warnings: list[str] = list(unit_info["warnings"])

    U, kept_units, agg_warnings = _unit_expression(ctx, units, genes,
                                                   min_unit_spots=min_unit_spots)
    warnings.extend(agg_warnings)
    # 单位矩阵的行顺序与 kept_units 一致，用布尔掩码取组
    unit_groups = np.asarray([u["group"] for u in kept_units], dtype=object)
    case_mask = unit_groups == case
    ctrl_mask = unit_groups == control
    n_case_units = int(case_mask.sum())
    n_ctrl_units = int(ctrl_mask.sum())
    # 推断的前提是「每组独立供体/切片 ≥ 2」；同一供体的多个区域不是独立生物学重复
    case_blocks = {u["block"] for u in kept_units if u["group"] == case}
    ctrl_blocks = {u["block"] for u in kept_units if u["group"] == control}
    n_case_blocks = len(case_blocks)
    n_ctrl_blocks = len(ctrl_blocks)
    inference = "unit_level" if (n_case_blocks >= 2 and n_ctrl_blocks >= 2) else "descriptive"
    unit_desc = f"{block_key}×{region_key}"
    notes_extra = [
        f"单位定义：{block_key} × {region_key}；每个单位把原始 counts 聚合成 pseudobulk（log1p CPM）"
        if ctx.get("counts") is not None else
        f"单位定义：{block_key} × {region_key}；每个单位聚合成一行表达（数据集没有 counts，用当前矩阵求和）",
        "空间数据里同一张切片的 spot 不是独立重复：按 spot 做检验会把 p 值压到 1e-9 量级（空间伪重复）；"
        "本分析按 供体/切片×区域 单位做检验",
    ]
    metrics_extra = {
        "inference": inference,
        "n_units": int(len(kept_units)),
        "n_case_units": n_case_units,
        "n_ctrl_units": n_ctrl_units,
        "n_case_blocks": n_case_blocks,
        "n_ctrl_blocks": n_ctrl_blocks,
        "unit": unit_desc,
        "n_units_dropped": int(len(units) - len(kept_units)),
        "min_unit_spots": min_unit_spots,
        "block_key": block_key,
        "region_key": region_key,
    }

    if inference == "descriptive":
        frame = _unit_de_frame(U, genes, case_mask, ctrl_mask, test=False)
        warnings.append(
            f"描述性（descriptive）：{case} 有 {n_case_blocks} 个独立供体/切片、{control} 有 {n_ctrl_blocks} 个"
            f"（共 {n_case_units}/{n_ctrl_units} 个单位，单位={unit_desc}）；每组独立供体/切片不足 2 个，"
            "属空间伪重复风险，不允许做推断（inferential），不报告 p 值")
        bullet = (f"{case} vs {control}：仅描述性汇总（独立供体/切片 {n_case_blocks} / {n_ctrl_blocks} 个，"
                  f"{n_case_units} / {n_ctrl_units} 个单位 {unit_desc}）；每组独立供体/切片 < 2，"
                  "不做推断、不给 p 值")
        desc_notes = [
            "方法：描述性（descriptive）汇总——单位级 pseudobulk 的组均值/占比；"
            "每组独立供体/切片 < 2，不做推断（inferential），不报告 p 值",
            "降级原因：每组只有 1 个供体/切片时，spot 级 p 值来自伪重复，不可作为推断证据",
        ]
        return _de_result(
            ctx, params, outdir, frame=frame, title=f"按「{gkey}」分组差异（{case} vs {control}，描述性）",
            case_label=case, ctrl_label=control, warnings=warnings,
            metrics_extra=metrics_extra, bullet_head=bullet,
            notes_extra=notes_extra, notes_base=desc_notes,
            table_note=f"{case} vs {control}；仅描述性汇总（n_units={len(kept_units)}，不报告推断 p 值）",
        )

    frame = _unit_de_frame(U, genes, case_mask, ctrl_mask, test=True)
    n_sig = int(((frame["p_adj"] <= 0.05) & (frame["logFC"] != 0)).sum())
    bullet = (f"{case} vs {control}（单位级 Wilcoxon + BH）：{n_case_units} vs {n_ctrl_units} 个单位"
              f"（{unit_desc}），显著 {n_sig} / {len(frame)} 个基因")
    unit_notes = [
        "方法：单位级 Wilcoxon 秩和检验（scipy mannwhitneyu，method=auto 小样本用精确分布）+ BH FDR；"
        "logFC = 单位 pseudobulk（log1p CPM）的组均值差；pct 为单位中表达 >0 的百分比",
        "空间数据里同一张切片的 spot 不是独立重复：按 spot 做检验会把 p 值压到 1e-9 量级（空间伪重复）；"
        "本分析按 供体/切片×区域 单位做检验",
    ]
    return _de_result(
        ctx, params, outdir, frame=frame, title=f"按「{gkey}」分组差异（{case} vs {control}，单位级）",
        case_label=case, ctrl_label=control, warnings=warnings,
        metrics_extra=metrics_extra, bullet_head=bullet, notes_extra=notes_extra,
        notes_base=unit_notes,
        table_note=f"{case} vs {control}；单位级（{unit_desc}）Wilcoxon + BH FDR",
    )


def _analysis_gene_domain_box(ctx, params=None, outdir=None):
    params = dict(params or {})
    outdir = _outdir(ctx, outdir)
    X, genes, _spots = _matrix(ctx)
    gene = _p_str(params.get("gene"))
    if gene and gene not in genes:
        raise sc.prerequisite(f"基因「{gene}」不在数据集里")
    domain_col = _p_str(params.get("domain_col"), "domains")
    key = domain_col if domain_col and domain_col not in DOMAIN_COLUMNS else ""
    if key:
        labels, hit = find_labels(ctx, [key])
        if labels is None:
            raise sc.prerequisite(f"数据集里没有列「{key}」，无法按该列分组（可用：domains / clusters）")
    else:
        labels, hit, _info = get_domain_labels(ctx, params)
    uniq = list(dict.fromkeys(labels))
    if len(uniq) < 2:
        raise sc.prerequisite("分组列只有 1 个取值，无法比较")
    if not gene:
        frame_m, target, _ = _marker_genes_for_domain(ctx, params, params.get("domain", 1))
        gene = str(_de_head(frame_m, 1)["gene"].iloc[0])
    labels_arr = np.asarray(labels, dtype=object)
    values = np.asarray(X[:, genes.index(gene)], dtype=float)
    rows = []
    groups_for_test = []
    for name in uniq:
        sel = labels_arr == name
        if not sel.any():
            continue
        vals = values[sel]
        groups_for_test.append(vals)
        rows.append(
            {
                "group": str(name),
                "n_spots": int(sel.sum()),
                "mean_expr": round(float(np.mean(vals)), 4),
                "median_expr": round(float(np.median(vals)), 4),
                "pct_expr": round(float(np.mean(vals > 0) * 100.0), 2),
            }
        )
    try:
        from scipy import stats

        kruskal = stats.kruskal(*groups_for_test) if len(groups_for_test) >= 2 else None
        stat = float(kruskal.statistic) if kruskal is not None else float("nan")
        pval = float(kruskal.pvalue) if kruskal is not None else float("nan")
    except Exception:  # noqa: BLE001
        stat, pval = float("nan"), float("nan")
    frame = pd.DataFrame(rows)
    table_obj = sc.table(
        "gene_domain_box",
        f"{gene} 在各组的表达",
        [
            {"key": "group", "label": "组", "type": "string"},
            {"key": "n_spots", "label": "spot 数", "type": "number"},
            {"key": "mean_expr", "label": "平均表达", "type": "number"},
            {"key": "median_expr", "label": "中位表达", "type": "number"},
            {"key": "pct_expr", "label": "表达比例(%)", "type": "number"},
        ],
        frame,
        note=f"Kruskal-Wallis H = {stat:.3f}，P = {pval:.3g}（分组列：{hit or domain_col}）",
    )
    figures: list[dict] = []
    warnings: list[str] = []
    try:
        figures.append(_gene_domain_box_figure(ctx, params, outdir, gene=gene, labels=labels,
                                               values=values, pval=pval, stat=stat))
    except Exception as exc:  # noqa: BLE001
        warnings.append(f"箱线图绘制失败：{exc}")
    return sc.result(
        title="区域表达箱线",
        bullets=[
            f"{gene} 在 {len(frame)} 个组的表达差异：Kruskal-Wallis H={stat:.3f}，P={pval:.3g}"
            + ("（显著）" if np.isfinite(pval) and pval < 0.05 else "（不显著）"),
            "各组均值：" + "、".join(f"{row['group']}={row['mean_expr']:.3f}" for _, row in frame.iterrows()),
        ],
        metrics={
            "gene": gene,
            "group_col": hit or domain_col,
            "kruskal_h": None if not np.isfinite(stat) else round(stat, 4),
            "p_value": None if not np.isfinite(pval) else float(pval),
            "n_groups": int(len(frame)),
        },
        tables=[table_obj],
        figures=figures,
        warnings=warnings,
        notes=["方法：分组箱线 + Kruskal-Wallis 检验（scipy）；内置实现"],
        data_out=[sc.data_file("gene_domain_box_csv", sc.export_csv(frame, outdir, "gene_domain_box", index=False),
                               "分组表达表")],
    )


# ---------------------------------------------------------------------------
# 图实现（供步骤/分析/PLOTS 共用）
# ---------------------------------------------------------------------------

def _de_volcano_figure(ctx, params, outdir, *, frame: pd.DataFrame, title: str = "", fig_name: str = "de_volcano") -> dict:
    import matplotlib.pyplot as plt

    padj_cut = 0.05
    lfc_cut = float(max(0.0, _p_float((params or {}).get("lfc"), 0.0)))
    top_n = int(max(1, min(_p_int((params or {}).get("top_n"), 20), 200)))
    plot = frame.copy()
    plot["neglog"] = -np.log10(np.clip(plot["p_adj"].to_numpy(dtype=float), 1e-300, 1.0))
    up = (plot["p_adj"] <= padj_cut) & (plot["logFC"] > lfc_cut)
    down = (plot["p_adj"] <= padj_cut) & (plot["logFC"] < -lfc_cut)
    other = ~(up | down)
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.scatter(plot.loc[other, "logFC"], plot.loc[other, "neglog"], s=9, color="#b9bfc8", alpha=0.55,
               label="ns", linewidths=0)
    ax.scatter(plot.loc[up, "logFC"], plot.loc[up, "neglog"], s=15, color="#c0392b", alpha=0.82,
               label=f"up ({int(up.sum())})", linewidths=0)
    ax.scatter(plot.loc[down, "logFC"], plot.loc[down, "neglog"], s=15, color="#2e6da4", alpha=0.82,
               label=f"down ({int(down.sum())})", linewidths=0)
    if lfc_cut > 0:
        for value in (lfc_cut, -lfc_cut):
            ax.axvline(value, color="#8a919c", linestyle="--", linewidth=0.8)
    ax.axhline(-math.log10(padj_cut), color="#8a919c", linestyle="--", linewidth=0.8)
    hits = plot.loc[up | down].nsmallest(top_n, "p_adj")
    if hits.empty:
        hits = plot.nsmallest(min(top_n, 8), "p_adj")
    for _, row in hits.iterrows():
        ax.annotate(str(row["gene"]), (float(row["logFC"]), float(row["neglog"])), fontsize=6.2,
                    alpha=0.9, xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("logFC (mean difference)")
    ax.set_ylabel("-log10(BH FDR)")
    ax.set_title(f"Volcano: {title}" if title else "Volcano plot", fontsize=10)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    return _fig_save(fig, ctx, params or {}, outdir, fig_name, "Volcano plot")


def _de_heatmap_figure(ctx, params, outdir, *, frame: pd.DataFrame, labels: list[str],
                       group_title: str = "", fig_name: str = "de_heatmap",
                       keep_labels: list[str] | None = None) -> dict:
    """Top 差异基因 × spot 的 z-score 热图（按组排序，最多每组 40 个 spot）。"""
    import matplotlib.pyplot as plt

    X, genes, _spots = _matrix(ctx)
    top_genes = _de_head(frame, int(max(3, min(_p_int((params or {}).get("top_n"), MAX_MARKER_GENES_FOR_FIG),
                                             MAX_MARKER_GENES_FOR_FIG))))["gene"].tolist()
    gidx = {g: i for i, g in enumerate(genes)}
    picked = [g for g in top_genes if g in gidx]
    if not picked:
        raise sc.prerequisite("没有可绘制的差异基因")
    labels_arr = np.asarray([str(v) for v in labels], dtype=object)
    uniq = [u for u in dict.fromkeys(labels_arr.tolist()) if (keep_labels is None or u in keep_labels)]
    if not uniq:
        uniq = list(dict.fromkeys(labels_arr.tolist()))
    rng = _rng(ctx)
    cols: list[int] = []
    for name in uniq:
        idx = np.flatnonzero(labels_arr == name)
        if idx.size > 40:
            idx = np.sort(rng.choice(idx, size=40, replace=False))
        cols.extend(idx.tolist())
    if not cols:
        raise sc.prerequisite("分组里没有 spot，无法画差异热图")
    mat = np.asarray(X[np.ix_(np.array(cols), [gidx[g] for g in picked])], dtype=float).T
    z = sc.zscore_rows(mat)
    fig, ax = plt.subplots(figsize=(max(7.0, len(cols) * 0.035 + 2.4), max(3.4, len(picked) * 0.11 + 1.4)))
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2, interpolation="nearest")
    ax.set_yticks(range(len(picked)))
    ax.set_yticklabels(picked, fontsize=5.6)
    step = max(1, len(cols) // 60)
    ticks = list(range(0, len(cols), step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([labels_arr[cols[i]] for i in ticks], rotation=90, fontsize=6)
    colors = sc.categorical_colors(uniq)
    for tick, pos in zip(ax.get_xticklabels(), ticks):
        tick.set_color(colors.get(labels_arr[cols[pos]], "#2b2f36"))
    ax.grid(False)
    ax.set_title(f"Top differential genes (row z-score) — {group_title}" if group_title
                 else "Top differential genes (row z-score)", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    return _fig_save(fig, ctx, params or {}, outdir, fig_name, "Differential gene heatmap")


def _dotplot_figure(ctx, params, outdir, *, frame: pd.DataFrame, labels: list[str],
                    focus: str, fig_name: str = "domain_dotplot") -> dict:
    """Top 标记基因 × 各域的点图（点大小=表达比例，颜色=平均 z-score）。"""
    import matplotlib.pyplot as plt

    X, genes, _spots = _matrix(ctx)
    top_n = int(max(3, min(_p_int((params or {}).get("top_n"), 20), 60)))
    cand = frame.assign(_abs=frame["logFC"].abs()).sort_values(["p_adj", "_abs"], ascending=[True, False])
    picked = [g for g in cand["gene"].head(top_n).tolist()]
    gidx = {g: i for i, g in enumerate(genes)}
    picked = [g for g in picked if g in gidx]
    if not picked:
        raise sc.prerequisite("没有可绘制的标记基因")
    labels_arr = np.asarray([str(v) for v in labels], dtype=object)
    uniq = list(dict.fromkeys(labels_arr.tolist()))
    mat = np.asarray(X[:, [gidx[g] for g in picked]], dtype=float)
    pidx = {gene: i for i, gene in enumerate(picked)}  # mat 只有 len(picked) 列，不能用全局 gidx 索引
    mean = mat.mean(axis=0)
    std = mat.std(axis=0)
    std = np.where(std > 0, std, 1.0)
    zmat = (mat - mean) / std
    fig, ax = plt.subplots(figsize=(max(6.6, len(picked) * 0.36 + 2.2), max(3.2, len(uniq) * 0.42 + 1.4)))
    xs, ys, sizes, colors, pts = [], [], [], [], []
    for yi, name in enumerate(uniq):
        mask = labels_arr == name
        if not mask.any():
            continue
        for xi, gene in enumerate(picked):
            col = mat[mask, pidx[gene]]
            xs.append(xi)
            ys.append(yi)
            sizes.append(float(np.mean(col > 0)) * 120.0 + 4.0)
            colors.append(float(np.mean(zmat[mask, pidx[gene]])))
    scatter = ax.scatter(xs, ys, s=sizes, c=colors, cmap="RdBu_r", vmin=-1.5, vmax=1.5,
                         edgecolors="#5a6068", linewidths=0.4)
    ax.set_xticks(range(len(picked)))
    ax.set_xticklabels(picked, rotation=90, fontsize=6.4)
    ax.set_yticks(range(len(uniq)))
    ax.set_yticklabels(uniq, fontsize=8)
    ax.set_xlabel("Top marker genes")
    ax.set_ylabel("Domain")
    ax.set_title(f"Domain marker expression (focus: {focus})", fontsize=10)
    ax.grid(True, linewidth=0.4)
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Mean z-score", fontsize=8)
    return _fig_save(fig, ctx, params or {}, outdir, fig_name, "Domain marker dot plot")


def _gene_domain_box_figure(ctx, params, outdir, *, gene: str, labels: list[str], values: np.ndarray,
                           pval: float, stat: float) -> dict:
    import matplotlib.pyplot as plt

    labels_arr = np.asarray([str(v) for v in labels], dtype=object)
    uniq = [u for u in dict.fromkeys(labels_arr.tolist())]
    data = [values[labels_arr == u] for u in uniq]
    fig, ax = plt.subplots(figsize=(max(5.6, len(uniq) * 0.85 + 2.0), 4.6))
    bp = ax.boxplot(data, patch_artist=True, widths=0.62, showfliers=False,
                    medianprops={"color": "#2b2f36", "linewidth": 1.2})
    colors = sc.palette(len(uniq))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
        patch.set_edgecolor("#5a6068")
    rng = _rng(ctx)
    for i, vals in enumerate(data):
        if vals.size == 0:
            continue
        take = vals if vals.size <= 200 else rng.choice(vals, size=200, replace=False)
        jitter = rng.uniform(-0.16, 0.16, size=take.size)
        ax.scatter(np.full(take.size, i + 1) + jitter, take, s=4, color="#4a4f57", alpha=0.35, linewidths=0)
    ax.set_xticks(range(1, len(uniq) + 1))
    ax.set_xticklabels(uniq, fontsize=9)
    ax.set_xlabel("Domain")
    ax.set_ylabel(f"{gene} expression")
    suffix = f"Kruskal-Wallis P = {pval:.3g}" if np.isfinite(pval) else ""
    ax.set_title(f"{gene} expression across domains" + (f" ({suffix})" if suffix else ""), fontsize=10)
    return _fig_save(fig, ctx, params or {}, outdir, "gene_domain_box", "Gene expression by domain")


def _composition_bar_figure(ctx, params, outdir, *, frame: pd.DataFrame, gkey: str) -> dict:
    import matplotlib.pyplot as plt

    if frame.empty:
        raise sc.prerequisite("没有可画的组成数据")
    pivot = frame.pivot_table(index="group", columns="domain", values="pct", aggfunc="sum", fill_value=0.0)
    domains = [str(c) for c in pivot.columns]
    groups = [str(i) for i in pivot.index]
    fig, ax = plt.subplots(figsize=(max(6.2, len(groups) * 0.72 + 2.2), 4.6))
    colors = sc.palette(len(domains))
    bottom = np.zeros(len(groups))
    for j, dom in enumerate(domains):
        vals = pivot[dom].to_numpy(dtype=float)
        ax.bar(groups, vals, bottom=bottom, color=colors[j % len(colors)], label=f"Domain {dom}",
               edgecolor="white", linewidth=0.6)
        bottom += vals
    ax.set_ylabel("Spot proportion (%)")
    ax.set_title(f"Domain composition by {gkey}", fontsize=10)
    ax.set_ylim(0, 100)
    ax.legend(fontsize=7, frameon=False, ncol=min(4, len(domains)), loc="upper center",
              bbox_to_anchor=(0.5, -0.16))
    ax.tick_params(axis="x", labelsize=8)
    return _fig_save(fig, ctx, params or {}, outdir, "domain_composition_bar", "Domain composition")


def _boundary_figure(ctx, params, outdir, *, res: dict) -> dict:
    import matplotlib.pyplot as plt

    labels = res["labels"]
    labels_arr = np.asarray(labels, dtype=object)
    coords = _coords(ctx)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0), gridspec_kw={"width_ratios": [1.05, 1.0]})
    ax = axes[0]
    other_mask = ~((labels_arr == res["domain_a"]) | (labels_arr == res["domain_b"]))
    ax.scatter(coords[other_mask, 0], coords[other_mask, 1], s=3.0, c="#d9dde3", linewidths=0,
               label="Other")
    ax.scatter(coords[labels_arr == res["domain_a"], 0], coords[labels_arr == res["domain_a"], 1],
               s=8.0, c=sc.PALETTE[0], linewidths=0, label=f"Domain {res['domain_a']}")
    ax.scatter(coords[labels_arr == res["domain_b"], 0], coords[labels_arr == res["domain_b"], 1],
               s=8.0, c=sc.PALETTE[2], linewidths=0, label=f"Domain {res['domain_b']}")
    bmask = np.asarray(res["boundary_mask"], dtype=bool)
    ax.scatter(coords[bmask, 0], coords[bmask, 1], s=16.0, facecolors="none", edgecolors="#111111",
               linewidths=0.7, label="Boundary")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(fontsize=7, frameon=False, loc="upper right", markerscale=1.4)
    ax.set_title(f"Boundary spots ({int(bmask.sum())}), crossing {res['crossing'] * 100:.1f}%", fontsize=10)
    ax2 = axes[1]
    g = res["gradient"]
    for side, color in ((res["domain_a"], sc.PALETTE[0]), (res["domain_b"], sc.PALETTE[2])):
        sub = g[g["side"] == str(side)].sort_values("ring")
        if sub.empty:
            continue
        ax2.plot(sub["mean_distance"], sub["mean_marker_score"], marker="o", ms=4, lw=1.4,
                 color=color, label=f"Domain {side} markers")
    ax2.set_xlabel("Distance to boundary")
    ax2.set_ylabel("Mean marker z-score")
    ax2.set_title(f"Ring gradient ({res['n_rings']} rings)", fontsize=10)
    ax2.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    return _fig_save(fig, ctx, params or {}, outdir, "boundary_plot", "Boundary strength")


def _distance_figure(ctx, params, outdir, *, res: dict) -> dict:
    import matplotlib.pyplot as plt

    coords = _coords(ctx)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0), gridspec_kw={"width_ratios": [1.05, 1.0]})
    ax = axes[0]
    sc.scatter_spatial(ax, ctx, values=res["dist_field"], point_size=_point_size(params or {}),
                       cmap="magma_r", legend_title=f"Distance to domain {res['target']}")
    ax.set_title(f"Distance field to domain {res['target']}", fontsize=10)
    ax2 = axes[1]
    frame = res["frame"]
    genes_show = res["genes"][:6]
    colors = sc.palette(len(genes_show))
    for i, gene in enumerate(genes_show):
        sub = frame[frame["gene"] == gene].sort_values("bin")
        if sub.empty:
            continue
        vals = sub["mean_expr"].to_numpy(dtype=float)
        vals = (vals - np.nanmean(vals)) / (np.nanstd(vals) or 1.0)
        ax2.plot(sub["mean_distance"], vals, marker="o", ms=3.4, lw=1.2, color=colors[i % len(colors)],
                 label=gene)
    ax2.set_xlabel(f"Distance to domain {res['target']}")
    ax2.set_ylabel("Mean expression (z per gene)")
    ax2.set_title(f"Expression vs distance ({res['bins']} bins)", fontsize=10)
    ax2.legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    return _fig_save(fig, ctx, params or {}, outdir, "distance_plot", "Distance to domain")


def _niche_figure(ctx, params, outdir, *, niche_labels: list[str]) -> dict:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    _spatial_ax(ax, ctx, params or {}, labels=niche_labels, legend_title="Niche")
    counts = pd.Series(niche_labels).value_counts()
    ax.set_title(f"Neighbourhood niches (n={len(counts)})", fontsize=10)
    return _fig_save(fig, ctx, params or {}, outdir, "niche_map", "Neighbourhood niches")


def _diversity_figure(ctx, params, outdir, *, shannon: np.ndarray, simpson: np.ndarray) -> dict:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0))
    sc.scatter_spatial(axes[0], ctx, values=shannon, point_size=_point_size(params or {}),
                       cmap="viridis", legend_title="Shannon")
    axes[0].set_title("Neighbourhood Shannon diversity", fontsize=10)
    sc.scatter_spatial(axes[1], ctx, values=simpson, point_size=_point_size(params or {}),
                       cmap="plasma", legend_title="Simpson")
    axes[1].set_title("Neighbourhood Simpson diversity", fontsize=10)
    fig.tight_layout()
    return _fig_save(fig, ctx, params or {}, outdir, "diversity_map", "Neighbourhood diversity")


# ---------------------------------------------------------------------------
# 数据获取辅助（PLOTS 用）
# ---------------------------------------------------------------------------

def _svg_records(ctx) -> list[dict]:
    results = ctx.get("results") if isinstance(ctx.get("results"), dict) else {}
    obj = results.get("svg")
    if isinstance(obj, list) and obj and isinstance(obj[0], dict):
        return obj
    if isinstance(obj, dict):
        for key in ("table", "rows", "records"):
            value = obj.get(key)
            if isinstance(value, list) and value and isinstance(value[0], dict):
                return value
    return []


def _de_frame_cached(ctx):
    results = ctx.get("results") if isinstance(ctx.get("results"), dict) else {}
    for key in ("de_domain", "de_regions", "de_groups", "de"):
        obj = results.get(key)
        frame = None
        if isinstance(obj, pd.DataFrame):
            frame = obj.copy()
        elif isinstance(obj, dict):
            table_obj = obj.get("table")
            if isinstance(table_obj, pd.DataFrame):
                frame = table_obj.copy()
            elif isinstance(table_obj, list) and table_obj and isinstance(table_obj[0], dict):
                frame = pd.DataFrame(table_obj)
        if frame is not None and not frame.empty and "gene" in frame.columns:
            return frame
    return None


def _resolve_de_frame(ctx, params) -> tuple[pd.DataFrame, str, str, list[str] | None]:
    """为 de_volcano / de_heatmap 选一个差异结果（优先复用缓存，否则现场计算）。"""
    params = dict(params or {})
    domain_a = params.get("domain_a")
    domain_b = params.get("domain_b")
    group_col = _p_str(params.get("group_col"))
    labels, src, _info = get_domain_labels(ctx, params)
    X, genes, _spots = _matrix(ctx)
    labels_arr = np.asarray(labels, dtype=object)
    if domain_a is not None and domain_b is not None:
        da = _require_label(domain_a, list(dict.fromkeys(labels)), "域 A")
        db = _require_label(domain_b, list(dict.fromkeys(labels)), "域 B")
        frame = _wilcoxon_frame(X, genes, labels_arr == da, labels_arr == db)
        return frame, f"Domain {da} vs {db}", f"domain {da} vs {db}", [da, db]
    if group_col:
        n = X.shape[0]
        groups, gkey = _group_values(ctx, n, group_col)
        if groups is not None:
            uniq = [g for g in dict.fromkeys(groups) if g not in ("NA", "")]
            if len(uniq) >= 2:
                case = _p_str(params.get("case")) or uniq[0]
                ctrl = _p_str(params.get("control")) or uniq[1]
                if case != ctrl:
                    garr = np.asarray(groups, dtype=object)
                    frame = _wilcoxon_frame(X, genes, garr == case, garr == ctrl)
                    return frame, f"{gkey}: {case} vs {ctrl}", gkey, None
    target = _require_label(params.get("domain", 1), list(dict.fromkeys(labels)), "域")
    frame = _wilcoxon_frame(X, genes, labels_arr == target, labels_arr != target)
    return frame, f"Domain {target} vs rest", f"domain {target}", None


# ---------------------------------------------------------------------------
# PLOTS（kind -> fn(ctx, params, outdir) -> figure）
# ---------------------------------------------------------------------------

def plot_svg_rank(ctx, params=None, outdir=None) -> dict:
    """Top SVG 排序图（优先用 svg 步骤缓存，否则现场计算）。"""
    params = dict(params or {})
    records = _svg_records(ctx)
    if records:
        frame = pd.DataFrame(records)
        if "statistic" not in frame.columns and "moran_i" in frame.columns:
            frame["statistic"] = frame["moran_i"]
        if "p_adj" not in frame.columns and "p" in frame.columns:
            frame["p_adj"] = frame["p"]
        n_top = int(max(1, min(_p_int(params.get("n_top"), 30), len(frame))))
        computed = {
            "top": frame.head(n_top),
            "method": str(params.get("method") or ("geary" if "geary_c" in frame.columns else "moran")),
        }
        return _svg_rank_figure(ctx, params, _outdir(ctx, outdir), computed=computed)
    step_params = {
        "n_top": _p_int(params.get("n_top"), 30),
        "method": _p_str(params.get("method"), "moran"),
        "permutations": _p_int(params.get("permutations"), 0),
        "genes": params.get("genes") or [],
    }
    return _svg_rank_figure(ctx, step_params, _outdir(ctx, outdir))


def plot_domain_map(ctx, params=None, outdir=None) -> dict:
    """空间域空间图（优先用数据集 domains 列）。"""
    return _domain_map_figure(ctx, params or {}, _outdir(ctx, outdir))


def plot_nhood_enrich(ctx, params=None, outdir=None) -> dict:
    """邻域富集 z-score 热图（现场计算，n_perms 默认 50）。"""
    params = dict(params or {})
    if params.get("n_perms") is None:
        params["n_perms"] = 50
    return _nhood_figure(ctx, params, _outdir(ctx, outdir))


def plot_cooccurrence(ctx, params=None, outdir=None) -> dict:
    """共现曲线（现场计算）。"""
    return _cooccurrence_figure(ctx, params or {}, _outdir(ctx, outdir))


def plot_domain_composition(ctx, params=None, outdir=None) -> dict:
    """域组成堆叠条形图。"""
    params = dict(params or {})
    labels, _src, _info = get_domain_labels(ctx, params)
    n = len(labels)
    groups, gkey = _group_values(ctx, n, _p_str(params.get("group_col")))
    if groups is None:
        groups = ["all"] * n
        gkey = "all spots"
    rows = []
    for g in list(dict.fromkeys(groups)):
        sel = np.asarray(groups, dtype=object) == g
        for d in list(dict.fromkeys(labels)):
            m = sel & (np.asarray(labels, dtype=object) == d)
            rows.append({"group": str(g), "domain": str(d), "n_spots": int(m.sum()),
                         "pct": round(float(m.sum() / max(1, sel.sum()) * 100.0), 2)})
    return _composition_bar_figure(ctx, params, _outdir(ctx, outdir), frame=pd.DataFrame(rows), gkey=gkey)


def plot_boundary(ctx, params=None, outdir=None) -> dict:
    """边界强度图（现场计算）。"""
    params = dict(params or {})
    res = _boundary_core(ctx, params)
    return _boundary_figure(ctx, params, _outdir(ctx, outdir), res=res)


def plot_distance(ctx, params=None, outdir=None) -> dict:
    """到目标域的距离场 + 分箱表达趋势（现场计算）。"""
    params = dict(params or {})
    res = _distance_core(ctx, params)
    return _distance_figure(ctx, params, _outdir(ctx, outdir), res=res)


def plot_niche_map(ctx, params=None, outdir=None) -> dict:
    """niche 空间图（优先用已有 niche 列，否则现场聚类）。"""
    params = dict(params or {})
    n = len(ctx.get("spots") or [])
    labels, hit = find_labels(ctx, ["niche", "niches"], n)
    if labels is None:
        labels_obj, _key = _cluster_labels(ctx, _p_str(params.get("cluster_key")))
        coords = _coords(ctx)
        k = int(max(2, min(_p_int(params.get("k"), 8), 30)))
        comp, uniq = _neighborhood_composition(labels_obj, coords, k, blocks=_blocks(ctx)[0])
        n_niches = int(max(2, min(_p_int(params.get("n_niches"), 5), 20)))
        n_niches = min(n_niches, max(2, comp.shape[0] // 10))
        from sklearn.cluster import KMeans

        niche = KMeans(n_clusters=n_niches, n_init=10, random_state=0).fit_predict(comp)
        labels = [str(int(v) + 1) for v in niche]
    return _niche_figure(ctx, params, _outdir(ctx, outdir), niche_labels=labels)


def plot_diversity(ctx, params=None, outdir=None) -> dict:
    """邻域多样性空间图（现场计算）。"""
    params = dict(params or {})
    labels, _key = _cluster_labels(ctx, _p_str(params.get("cluster_key")))
    res = _diversity_core(labels, _coords(ctx),
                          int(max(2, min(_p_int(params.get("k"), 6), 30))),
                          blocks=_blocks(ctx)[0])
    return _diversity_figure(ctx, params, _outdir(ctx, outdir), shannon=res["shannon"], simpson=res["simpson"])


def plot_de_volcano(ctx, params=None, outdir=None) -> dict:
    """区域差异火山图（有缓存优先，否则现场计算）。"""
    params = dict(params or {})
    frame = _de_frame_cached(ctx)
    if frame is None or params.get("domain_a") is not None or params.get("group_col"):
        frame, title, _tag, _keep = _resolve_de_frame(ctx, params)
    else:
        title = "cached DE"
    return _de_volcano_figure(ctx, params, _outdir(ctx, outdir), frame=frame, title=title,
                              fig_name="de_volcano")


def plot_de_heatmap(ctx, params=None, outdir=None) -> dict:
    """区域差异热图（有缓存优先，否则现场计算）。"""
    params = dict(params or {})
    labels, _src, _info = get_domain_labels(ctx, params)
    frame = _de_frame_cached(ctx)
    keep = None
    if frame is None or params.get("domain_a") is not None or params.get("group_col"):
        frame, title, _tag, keep = _resolve_de_frame(ctx, params)
    else:
        title = "cached DE"
    return _de_heatmap_figure(ctx, params, _outdir(ctx, outdir), frame=frame, labels=labels,
                              group_title=title, fig_name="de_heatmap", keep_labels=keep)


def plot_gene_domain_box(ctx, params=None, outdir=None) -> dict:
    """分域表达箱线图（gene 缺省时取域 1 的 top 标记）。"""
    params = dict(params or {})
    X, genes, _spots = _matrix(ctx)
    domain_col = _p_str(params.get("domain_col"), "domains")
    if domain_col and domain_col not in DOMAIN_COLUMNS:
        labels, hit = find_labels(ctx, [domain_col])
        if labels is None:
            raise sc.prerequisite(f"数据集里没有列「{domain_col}」（可用：domains / clusters）")
    else:
        labels, hit, _info = get_domain_labels(ctx, params)
    gene = _p_str(params.get("gene"))
    if gene and gene not in genes:
        raise sc.prerequisite(f"基因「{gene}」不在数据集里")
    if not gene:
        frame_m, target, _ = _marker_genes_for_domain(ctx, params, params.get("domain", 1))
        gene = str(_de_head(frame_m, 1)["gene"].iloc[0])
    labels_arr = np.asarray(labels, dtype=object)
    values = np.asarray(X[:, genes.index(gene)], dtype=float)
    groups = [values[labels_arr == u] for u in list(dict.fromkeys(labels_arr.tolist()))]
    pval, stat = float("nan"), float("nan")
    try:
        from scipy import stats

        if len(groups) >= 2:
            out = stats.kruskal(*groups)
            stat, pval = float(out.statistic), float(out.pvalue)
    except Exception:  # noqa: BLE001
        pass
    return _gene_domain_box_figure(ctx, params, _outdir(ctx, outdir), gene=gene, labels=labels,
                                   values=values, pval=pval, stat=stat)


def plot_domain_dotplot(ctx, params=None, outdir=None) -> dict:
    """域标记基因点图（现场计算域 1 或指定域的标记）。"""
    params = dict(params or {})
    frame, target, labels = _marker_genes_for_domain(ctx, params, params.get("domain", 1))
    return _dotplot_figure(ctx, params, _outdir(ctx, outdir), frame=frame, labels=labels, focus=target,
                           fig_name="domain_dotplot")


# ---------------------------------------------------------------------------
# 参数声明
# ---------------------------------------------------------------------------

_SVG_PARAMS = [
    {"key": "n_top", "label": "输出前 N 个 SVG", "type": "number", "default": 30},
    {"key": "method", "label": "统计量", "type": "select", "default": "moran",
     "options": [{"value": "moran", "label": "Moran's I"}, {"value": "geary", "label": "Geary's C"}]},
    {"key": "permutations", "label": "置换次数（0=正态近似）", "type": "number", "default": 0},
    {"key": "genes", "label": "只看指定基因（可空）", "type": "list", "default": []},
]

_DOMAIN_PARAMS = [
    {"key": "n_domains", "label": "域数量", "type": "number", "default": 6},
    {"key": "method", "label": "方法", "type": "select", "default": "mrf",
     "options": [{"value": "mrf", "label": "MRF 平滑"}, {"value": "graph", "label": "融合图 Louvain"}]},
    {"key": "smooth", "label": "平滑权重", "type": "number", "default": 0.5},
    {"key": "iterations", "label": "最大迭代", "type": "number", "default": 20},
]

_CLUSTER_PARAM = {"key": "cluster_key", "label": "聚类列", "type": "text", "default": "clusters"}


# ---------------------------------------------------------------------------
# STEPS / ANALYSES / PLOTS
# ---------------------------------------------------------------------------

STEPS = [
    {
        "id": "svg",
        "title": "空间高变基因",
        "icon": "svg",
        "requires": ["cluster"],
        "params": _SVG_PARAMS,
        "run": run_svg,
    },
    {
        "id": "domains",
        "title": "空间域",
        "icon": "domain",
        "requires": ["cluster"],
        "params": _DOMAIN_PARAMS,
        "run": run_domains,
    },
]

ANALYSES = [
    {
        "id": "nhood_enrichment",
        "title": "邻域富集",
        "group": "邻域与域",
        "method": "类别对邻接计数 + 标签置换 z-score（squidpy nhood_enrichment 思想；Schürch 2020）；内置实现",
        "needs": ["cluster"],
        "params": [
            _CLUSTER_PARAM,
            {"key": "n_perms", "label": "置换次数", "type": "number", "default": 100},
            {"key": "n_neighbors", "label": "空间近邻数", "type": "number", "default": 6},
        ],
        "run": _analysis_nhood_enrichment,
    },
    {
        "id": "co_occurrence",
        "title": "共现分析",
        "group": "邻域与域",
        "method": "环形分箱条件共现概率（squidpy co_occurrence 思想）；内置实现",
        "needs": ["cluster"],
        "params": [
            _CLUSTER_PARAM,
            {"key": "n_rings", "label": "距离环数", "type": "number", "default": 5},
            {"key": "interval", "label": "距离区间 [lo,hi]（可空）", "type": "list", "default": None},
        ],
        "run": _analysis_co_occurrence,
    },
    {
        "id": "domains_mrf",
        "title": "空间域（MRF 平滑）",
        "group": "邻域与域",
        "method": "PCA + 空间图 + MRF 迭代平滑（BayesSpace 2021 思想）；内置实现",
        "needs": ["cluster"],
        "params": _DOMAIN_PARAMS,
        "run": _analysis_domains_mrf,
    },
    {
        "id": "domain_markers",
        "title": "区域标记基因",
        "group": "邻域与域",
        "method": "Wilcoxon 秩和 + BH FDR（scipy）；内置实现",
        "needs": ["cluster"],
        "params": [
            {"key": "domain", "label": "目标域", "type": "text", "default": 1},
            {"key": "top_n", "label": "输出基因数", "type": "number", "default": 20},
        ],
        "run": _analysis_domain_markers,
    },
    {
        "id": "domain_composition",
        "title": "区域组成",
        "group": "邻域与域",
        "method": "各分组/切片在各域的 spot 比例（组成分析）；内置实现",
        "needs": ["cluster"],
        "params": [{"key": "group_col", "label": "分组列（可空=自动）", "type": "text", "default": None}],
        "run": _analysis_domain_composition,
    },
    {
        "id": "boundary_score",
        "title": "边界强度",
        "group": "邻域与域",
        "method": "跨域邻接比例 + 到边界的环形表达梯度；内置实现",
        "needs": ["cluster"],
        "params": [
            {"key": "domain_a", "label": "域 A", "type": "text", "default": 1},
            {"key": "domain_b", "label": "域 B", "type": "text", "default": 2},
            {"key": "n_rings", "label": "环数", "type": "number", "default": 3},
        ],
        "run": _analysis_boundary_score,
    },
    {
        "id": "distance_to_domain",
        "title": "到区域距离",
        "group": "邻域与域",
        "method": "欧氏距离场 + 等频分箱表达趋势；内置实现",
        "needs": ["cluster"],
        "params": [
            {"key": "domain", "label": "目标域", "type": "text", "default": 1},
            {"key": "genes", "label": "基因（空=自动取标记）", "type": "list", "default": []},
            {"key": "bins", "label": "距离分箱数", "type": "number", "default": 10},
        ],
        "run": _analysis_distance_to_domain,
    },
    {
        "id": "niche_cluster",
        "title": "细胞邻域（niche）",
        "group": "邻域与域",
        "method": "邻域组成向量 + k-means（Schürch 2020；Bhate 2022）；内置实现",
        "needs": ["cluster"],
        "params": [
            _CLUSTER_PARAM,
            {"key": "k", "label": "近邻数", "type": "number", "default": 8},
            {"key": "n_niches", "label": "niche 数量", "type": "number", "default": 5},
        ],
        "run": _analysis_niche_cluster,
    },
    {
        "id": "tissue_diversity",
        "title": "邻域多样性",
        "group": "邻域与域",
        "method": "邻域组成的 Shannon/Simpson 多样性；内置实现",
        "needs": ["cluster"],
        "params": [
            _CLUSTER_PARAM,
            {"key": "k", "label": "近邻数", "type": "number", "default": 6},
        ],
        "run": _analysis_tissue_diversity,
    },
    {
        "id": "de_domain",
        "title": "域 vs 其余 差异",
        "group": "区域差异",
        "method": "Wilcoxon 秩和 + BH FDR；内置实现",
        "needs": ["cluster"],
        "params": [
            {"key": "domain", "label": "目标域", "type": "text", "default": 1},
            {"key": "top_n", "label": "输出基因数", "type": "number", "default": 200},
        ],
        "run": _analysis_de_domain,
    },
    {
        "id": "de_regions",
        "title": "两区域差异",
        "group": "区域差异",
        "method": "Wilcoxon 秩和 + BH FDR + logFC；内置实现",
        "needs": ["cluster"],
        "params": [
            {"key": "domain_a", "label": "域 A", "type": "text", "default": 1},
            {"key": "domain_b", "label": "域 B", "type": "text", "default": 2},
            {"key": "top_n", "label": "输出基因数", "type": "number", "default": 200},
        ],
        "run": _analysis_de_regions,
    },
    {
        "id": "de_groups",
        "title": "按分组差异",
        "group": "区域差异",
        "method": "Wilcoxon 秩和 + BH FDR（spot 级分组）；内置实现",
        "needs": ["groups"],
        "params": [
            {"key": "group_col", "label": "分组列", "type": "text", "default": None},
            {"key": "case", "label": "实验组", "type": "text", "default": None},
            {"key": "control", "label": "对照组", "type": "text", "default": None},
            {"key": "top_n", "label": "输出基因数", "type": "number", "default": 200},
        ],
        "run": _analysis_de_groups,
    },
    {
        "id": "gene_domain_box",
        "title": "区域表达箱线",
        "group": "区域差异",
        "method": "分组箱线 + Kruskal-Wallis 检验（scipy）；内置实现",
        "needs": ["cluster"],
        "params": [
            {"key": "gene", "label": "基因（空=域 1 的 top 标记）", "type": "text", "default": ""},
            {"key": "domain_col", "label": "分组列", "type": "text", "default": "domains"},
        ],
        "run": _analysis_gene_domain_box,
    },
]

PLOTS = {
    "svg_rank_plot": plot_svg_rank,
    "domain_map": plot_domain_map,
    "nhood_enrich": plot_nhood_enrich,
    "cooccurrence": plot_cooccurrence,
    "domain_composition_bar": plot_domain_composition,
    "boundary_plot": plot_boundary,
    "distance_plot": plot_distance,
    "niche_map": plot_niche_map,
    "diversity_map": plot_diversity,
    "de_volcano": plot_de_volcano,
    "de_heatmap": plot_de_heatmap,
    "gene_domain_box": plot_gene_domain_box,
    "domain_dotplot": plot_domain_dotplot,
}
