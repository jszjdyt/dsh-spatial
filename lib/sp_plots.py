"""
dsh-spatial — 基础绘图模块（sp-base）

导出 `PLOTS`：`kind -> fn(ctx, params, outdir) -> save_fig 返回的 dict`
（registry 把所有模块的 PLOTS 合并成全局表，`plot` 命令按 kind 分发，DESIGN §5）。

本文件覆盖的 9 个 kind（DESIGN §5 白名单中属于基础图的部分）：
    qc_spatial     每 spot 的 counts / 检出基因数 / 线粒体占比空间图
    qc_violin      QC 指标小提琴图（可按分组）
    pca_variance   PCA 碎石图（解释方差）
    embedding      二维嵌入散点（t-SNE / UMAP / PCA）
    cluster_bar    各簇 spot 数柱状图
    spatial_scatter空间散点（按簇/域/分组/基因上色）
    gene_box       基因表达箱线图（按簇/分组）
    heatmap        基因 × (簇/spot) 热图，行 z-score + 层次聚类排序
    dotplot        簇 × 基因 平均表达/表达比例点图

通用参数（DESIGN §5）：
    color_by  "clusters|domains|group|gene"（或 sample_info 里的任意列名）
    gene      基因名（color_by=gene / gene_box 时使用）
    genes     基因列表（heatmap / dotplot / gene_box 批量）
    point_size=6, max_points=20000, dpi=140, image=false（叠加 H&E 底图）

约定：
  * 图上文字一律英文（DESIGN §2.5）；中文只出现在返回给 UI 的 title/warnings。
  * 只读 ctx：不修改 ctx 的 X/coords（需要子集时新建浅拷贝）。
  * 缺前置 → sp_common.prerequisite（中文原因）；参数非法 → bad_request。
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

import sp_common as C
from sp_common import (
    EngineError,
    bad_request,
    categorical_colors,
    new_fig,
    palette,
    prerequisite,
    save_fig,
    scatter_spatial,
    zscore_rows,
)

DEFAULT_POINT_SIZE = 6.0
DEFAULT_MAX_POINTS = 20000
DEFAULT_DPI = 140

QC_METRICS = ("counts", "genes", "pct_mt")
QC_LABELS = {"counts": "Total counts", "genes": "Genes detected", "pct_mt": "Mitochondrial %"}


# ---------------------------------------------------------------------------
# 通用工具
# ---------------------------------------------------------------------------

def _param(params, key, default=None):
    value = (params or {}).get(key)
    return default if value is None else value


def _as_float(params, key, default: float) -> float:
    try:
        return float(_param(params, key, default))
    except (TypeError, ValueError):
        return float(default)


def _as_int(params, key, default: int) -> int:
    try:
        return int(_param(params, key, default))
    except (TypeError, ValueError):
        return int(default)


def _as_bool(params, key, default: bool = False) -> bool:
    value = _param(params, key, default)
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


def _dense(matrix) -> np.ndarray:
    """scipy sparse / ndarray → ndarray（float32）；None → 空。"""
    if matrix is None:
        return np.zeros((0, 0), dtype=np.float32)
    if hasattr(matrix, "toarray"):
        return np.asarray(matrix.toarray(), dtype=np.float32)
    return np.asarray(matrix, dtype=np.float32)


def _X(ctx) -> np.ndarray:
    arr = _dense((ctx or {}).get("X"))
    if arr.ndim != 2 or arr.size == 0:
        raise prerequisite("表达矩阵为空，请先导入数据")
    return arr


def _coords(ctx) -> np.ndarray:
    coords = (ctx or {}).get("coords")
    if coords is None:
        raise prerequisite("该数据集没有空间坐标，无法画空间图")
    arr = np.asarray(coords, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2 or arr.shape[0] == 0:
        raise prerequisite("该数据集没有空间坐标，无法画空间图")
    return arr


def _genes(ctx) -> list[str]:
    return [str(g) for g in ((ctx or {}).get("genes") or [])]


def _gene_index(ctx) -> dict:
    return {str(gene): index for index, gene in enumerate(_genes(ctx))}


def _info(ctx) -> pd.DataFrame:
    info = (ctx or {}).get("sample_info")
    if isinstance(info, pd.DataFrame) and len(info) == len((ctx or {}).get("spots") or []):
        return info
    spots = [str(s) for s in ((ctx or {}).get("spots") or [])]
    return pd.DataFrame(index=pd.Index(spots, name="spot"))


def _n_spots(ctx) -> int:
    matrix = (ctx or {}).get("X")
    if matrix is None:
        return len((ctx or {}).get("spots") or [])
    return int(np.asarray(matrix).shape[0])


def _sort_labels(values) -> list[str]:
    """标签排序：能转数字的按数字升序，否则按字典序（保证图例顺序稳定）。"""
    uniq = list(dict.fromkeys([str(v) for v in values]))

    def key(value: str):
        try:
            return (0, float(value), "")
        except (TypeError, ValueError):
            return (1, 0.0, value)

    return sorted(uniq, key=key)


def _info_column(ctx, key) -> np.ndarray | None:
    """从 sample_info 取与 spot 等长的列；取不到返回 None。"""
    if not key:
        return None
    info = _info(ctx)
    if str(key) in info.columns and len(info) == len((ctx or {}).get("spots") or []):
        return np.asarray([str(v) for v in info[str(key)].tolist()], dtype=object)
    return None


def _cluster_labels(ctx) -> np.ndarray:
    """当前聚类标签（ctx.cluster → sample_info[cluster_key/clusters]）。"""
    labels = (ctx or {}).get("cluster")
    if labels is not None:
        arr = np.asarray(labels).ravel()
        if arr.size == len((ctx or {}).get("spots") or []):
            return np.asarray([str(v) for v in arr.tolist()], dtype=object)
    for key in ((ctx or {}).get("cluster_key"), "clusters", "domains"):
        found = _info_column(ctx, key)
        if found is not None:
            return found
    raise prerequisite("先跑聚类（cluster 步骤），再进行该图")


def _try_cluster_labels(ctx) -> np.ndarray | None:
    try:
        return _cluster_labels(ctx)
    except EngineError:
        return None


def _subsample(n: int, max_points: int, seed: int = 0) -> np.ndarray:
    """确定性下采样（超过 max_points 时随机取，seed 固定 → 同图可复现）。"""
    n = int(n)
    limit = int(max(1, min(int(max_points), n)))
    if n <= limit:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, size=limit, replace=False))


def _sub_ctx(ctx, index: np.ndarray) -> dict:
    """坐标/图像子集（只读：不修改原 ctx）。"""
    sub = dict(ctx or {})
    coords = (ctx or {}).get("coords")
    sub["coords"] = np.asarray(coords, dtype=float)[index]
    return sub


def _finalize(fig, outdir: str, name: str, title: str, dpi: int, params) -> dict:
    filename = str(_param(params, "_filename", name))
    return save_fig(fig, outdir, filename, title=title, dpi=int(dpi))


def _draw_categorical(ax, x, y, labels, point_size: float, *, title: str = "",
                      legend: bool = True, marker_scale: float = 2.0):
    """类别散点（颜色按标签首次出现顺序分配，标签排序后给顺序色）。"""
    labs = np.asarray([str(v) for v in labels], dtype=object)
    colors = categorical_colors(_sort_labels(labs))
    for label in _sort_labels(labs):
        mask = labs == label
        if not mask.any():
            continue
        ax.scatter(x[mask], y[mask], s=point_size, c=colors[label], linewidths=0,
                   alpha=0.9, label=label)
    if legend:
        n_label = len(colors)
        ax.legend(markerscale=marker_scale, fontsize=6.8 if n_label > 10 else 7.4,
                  frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  title=f"n={len(labs)}", title_fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title)
    return colors


def _draw_continuous(ax, x, y, values, point_size: float, *, title: str = "",
                     cmap: str = "viridis", label: str = ""):
    """连续数值散点 + 色带。"""
    vals = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(vals)
    if not finite.any():
        raise prerequisite("没有可绘制的数值（全为 NaN）")
    scatter = ax.scatter(x[finite], y[finite], s=point_size, c=vals[finite], cmap=cmap,
                         linewidths=0, alpha=0.9)
    try:
        bar = ax.figure.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
        bar.ax.tick_params(labelsize=7)
        if label:
            bar.set_label(label, fontsize=8)
    except Exception:  # noqa: BLE001
        pass
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title)
    return scatter


def _violin(ax, groups, positions, colors, *, ylabel: str = "", title: str = ""):
    """按组画小提琴（组内样本数少时加抖动散点）。"""
    data = [np.asarray(values, dtype=np.float64) for values in groups]
    data = [values[np.isfinite(values)] for values in data]
    if not any(values.size for values in data):
        raise prerequisite("没有可绘制的 QC 数值")
    try:
        parts = ax.violinplot(data, positions=positions, showmedians=True,
                              showextrema=False, widths=0.72)
        for body, color in zip(parts["bodies"], colors):
            body.set_facecolor(color)
            body.set_edgecolor("#5a6068")
            body.set_alpha(0.72)
        if "cmedians" in parts:
            parts["cmedians"].set_color("#2b2f36")
            parts["cmedians"].set_linewidth(1.1)
    except Exception:  # noqa: BLE001 - 极小样本时 violin 可能失败，退化为散点
        for position, values, color in zip(positions, data, colors):
            if values.size:
                ax.scatter(np.full(values.size, position), values, s=3, color=color, alpha=0.6)
    ax.set_xticks(list(positions))
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)


# ---------------------------------------------------------------------------
# QC 指标
# ---------------------------------------------------------------------------

def _qc_frame(ctx) -> tuple[pd.DataFrame, str]:
    """每 spot 的 counts / 检出基因数 / 线粒体占比。

    线粒体基因前缀：`MT-`（人）/ `mt-`（鼠），大小写不敏感；同时接受 `MT.` 记法。
    有原始 counts 时用 counts（X_kind=counts），否则退回当前表达矩阵并标注 value kind。
    """
    genes = _genes(ctx)
    counts = (ctx or {}).get("counts")
    if counts is not None:
        matrix = _dense(counts)
        kind = "counts"
    else:
        matrix = _dense((ctx or {}).get("X"))
        kind = str((ctx or {}).get("X_kind") or "unknown")
    if matrix.ndim != 2 or matrix.size == 0:
        raise prerequisite("表达矩阵为空，请先导入数据")
    if len(genes) != matrix.shape[1]:
        raise EngineError("internal", "gene 列表与表达矩阵列数不一致")
    row_sum = np.asarray(matrix.sum(axis=1), dtype=np.float64).ravel()
    detect = np.asarray((matrix > 0).sum(axis=1), dtype=np.float64).ravel()
    mito = np.array([str(g).upper().startswith(("MT-", "MT.")) for g in genes], dtype=bool)
    if mito.any():
        mito_sum = np.asarray(matrix[:, mito].sum(axis=1), dtype=np.float64).ravel()
    else:
        mito_sum = np.zeros(matrix.shape[0], dtype=np.float64)
    safe = np.maximum(row_sum, 1e-12)
    pct_mt = np.where(row_sum > 0, 100.0 * mito_sum / safe, 0.0)
    frame = pd.DataFrame({"counts": row_sum, "genes": detect, "pct_mt": pct_mt})
    return frame, kind


def _metric_list(params, default=QC_METRICS) -> list[str]:
    raw = _param(params, "metrics", list(default))
    if isinstance(raw, str):
        raw = [item.strip() for item in raw.replace(",", " ").split() if item.strip()]
    if not isinstance(raw, (list, tuple)):
        raise bad_request("metrics 需要是数组（counts/genes/pct_mt）")
    metrics = [str(item) for item in raw if str(item) in QC_METRICS]
    if not metrics:
        raise bad_request("metrics 只能取 counts / genes / pct_mt")
    return metrics


# ---------------------------------------------------------------------------
# 1. qc_spatial
# ---------------------------------------------------------------------------

def plot_qc_spatial(ctx, params, outdir):
    """空间 QC 分布图。

    每个指标一个空间面板（坐标 y 轴翻转、等比例，见 sp_common.scatter_spatial）；
    默认 counts / genes / pct_mt 三联图，可传 metrics=["counts"] 只要一张。
    """
    frame, kind = _qc_frame(ctx)
    metrics = _metric_list(params)
    coords = _coords(ctx)
    point_size = _as_float(params, "point_size", DEFAULT_POINT_SIZE)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    index = _subsample(frame.shape[0], _as_int(params, "max_points", DEFAULT_MAX_POINTS))
    sub = _sub_ctx(ctx, index)
    fig, axes = new_fig(4.1 * len(metrics) + 0.5, 4.3, ncols=len(metrics))
    axes = np.atleast_1d(axes).ravel()
    for ax, metric in zip(axes, metrics):
        values = frame[metric].to_numpy(dtype=np.float64)[index]
        label = QC_LABELS[metric] if kind == "counts" else f"{QC_LABELS[metric]} ({kind})"
        scatter_spatial(ax, sub, values=values, point_size=point_size,
                        cmap="magma_r" if metric == "pct_mt" else "viridis",
                        legend_title=label, image=_as_bool(params, "image", False))
        ax.set_title(label, fontsize=10)
    return _finalize(fig, outdir, "qc_spatial", "Spatial QC", dpi, params)


# ---------------------------------------------------------------------------
# 2. qc_violin
# ---------------------------------------------------------------------------

def plot_qc_violin(ctx, params, outdir):
    """QC 指标按分组的小提琴图（无分组时单组）；阈值线来自步骤参数。"""
    frame, kind = _qc_frame(ctx)
    metrics = _metric_list(params)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    group_col = _param(params, "group_col", None)
    labels = _info_column(ctx, group_col) if group_col else None
    thresholds = {
        "counts": _param(params, "min_counts", None),
        "genes": _param(params, "min_genes", None),
        "pct_mt": _param(params, "max_pct_mt", None),
    }
    fig, axes = new_fig(3.6 * len(metrics) + 0.5, 4.2, ncols=len(metrics))
    axes = np.atleast_1d(axes).ravel()
    for ax, metric in zip(axes, metrics):
        values = frame[metric].to_numpy(dtype=np.float64)
        if labels is not None:
            order = _sort_labels(labels)
            groups = [values[np.asarray(labels) == group] for group in order]
            positions = np.arange(1, len(order) + 1)
            colors = palette(len(order))
            _violin(ax, groups, positions, colors, ylabel=QC_LABELS[metric],
                    title=f"{QC_LABELS[metric]} by {group_col}")
            ax.set_xticklabels(order, rotation=45 if len(order) > 6 else 0, fontsize=7)
        else:
            _violin(ax, [values], [1], [palette(1)[0]], ylabel=QC_LABELS[metric],
                    title=QC_LABELS[metric])
            ax.set_xticklabels(["all spots"], fontsize=7)
        threshold = thresholds.get(metric)
        if threshold is not None:
            try:
                ax.axhline(float(threshold), color="#e15759", linestyle="--", linewidth=1.0,
                           label=f"threshold {float(threshold):g}")
                ax.legend(fontsize=7, frameon=False, loc="best")
            except (TypeError, ValueError):
                pass
    if kind != "counts":
        fig.suptitle(f"QC metrics (values: {kind}, not raw counts)", fontsize=10)
    return _finalize(fig, outdir, "qc_violin", "QC distributions", dpi, params)


# ---------------------------------------------------------------------------
# 3. pca_variance
# ---------------------------------------------------------------------------

def _pca_variance(ctx):
    """(variance ratios, source note)。优先取 pca 步骤缓存的真实解释方差。"""
    results = (ctx or {}).get("results") or {}
    cached = results.get("pca") if isinstance(results.get("pca"), dict) else {}
    ratios = cached.get("variance_ratio") if isinstance(cached, dict) else None
    if ratios:
        arr = np.asarray(ratios, dtype=np.float64)
        if arr.size:
            return arr, "explained variance (pca step)"
    obsm = (ctx or {}).get("obsm") or {}
    for key in ("pca", "X_pca"):
        if key in obsm:
            matrix = np.asarray(obsm[key], dtype=np.float64)
            if matrix.ndim == 2 and matrix.shape[1]:
                var = matrix.var(axis=0, ddof=1)
                total = float(var.sum())
                if total > 0:
                    return var / total, "relative variance of stored PCs (approximate)"
            raise prerequisite("PCA 坐标缺失解释方差信息，请重跑 pca 步骤")
    raise prerequisite("先跑 PCA（pca 步骤）")


def plot_pca_variance(ctx, params, outdir):
    """PCA 碎石图：单主成分解释方差柱 + 累计解释方差折线。"""
    ratios, source = _pca_variance(ctx)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    n_show = int(min(ratios.size, _as_int(params, "n_pcs", ratios.size)))
    ratios = ratios[:n_show]
    positions = np.arange(1, n_show + 1)
    fig, ax = new_fig(7.0, 4.0)
    ax.bar(positions, ratios * 100.0, color="#4e79a7", edgecolor="#2c4a66", linewidth=0.4,
           alpha=0.9, label="individual")
    ax.set_xlabel("principal component")
    ax.set_ylabel("explained variance (%)")
    ax2 = ax.twinx()
    ax2.plot(positions, np.cumsum(ratios) * 100.0, color="#e15759", linewidth=1.2,
             marker="o", markersize=3, label="cumulative")
    ax2.set_ylabel("cumulative (%)")
    ax2.set_ylim(0, 105)
    ax.set_xticks(positions[:: max(1, n_show // 15)])
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles + handles2, labels + labels2, fontsize=8, frameon=False, loc="center right")
    ax.set_title(f"PCA variance ({source}; PC1 = {ratios[0] * 100:.1f}%)")
    return _finalize(fig, outdir, "pca_variance", "PCA scree plot", dpi, params)


# ---------------------------------------------------------------------------
# 4. embedding
# ---------------------------------------------------------------------------

def _embedding(ctx) -> tuple[np.ndarray, str]:
    obsm = (ctx or {}).get("obsm") or {}
    for key in ("tsne", "umap", "X_tsne", "X_umap", "pca", "X_pca", "X_emb"):
        if key in obsm:
            matrix = np.asarray(obsm[key], dtype=np.float64)
            if matrix.ndim == 2 and matrix.shape[1] >= 2:
                name = key.replace("X_", "")
                return matrix[:, :2], name
    raise prerequisite("先跑 embed（二维嵌入）或 pca 步骤")


def plot_embedding(ctx, params, outdir):
    """二维嵌入散点（t-SNE / UMAP / PCA），按簇/分组/基因上色。"""
    matrix, layout = _embedding(ctx)
    resolved = _resolve_color(ctx, params)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    point_size = _as_float(params, "point_size", DEFAULT_POINT_SIZE * 0.8)
    index = _subsample(matrix.shape[0], _as_int(params, "max_points", DEFAULT_MAX_POINTS))
    x = matrix[index, 0]
    y = matrix[index, 1]
    fig, ax = new_fig(6.4, 5.2)
    if resolved["kind"] == "num":
        _draw_continuous(ax, x, y, resolved["values"][index], point_size,
                         title=f"Embedding ({layout}) - {resolved['title']}",
                         cmap="viridis", label=resolved["title"])
    else:
        _draw_categorical(ax, x, y, np.asarray(resolved["labels"], dtype=object)[index],
                          point_size, title=f"Embedding ({layout}) - {resolved['title']}")
    ax.set_xlabel("dimension 1")
    ax.set_ylabel("dimension 2")
    return _finalize(fig, outdir, "embedding", "Embedding", dpi, params)


# ---------------------------------------------------------------------------
# 5. cluster_bar
# ---------------------------------------------------------------------------

def plot_cluster_bar(ctx, params, outdir):
    """各簇 spot 数柱状图（按大小降序）。"""
    labels = _cluster_labels(ctx)
    order = _sort_labels(labels)
    counts = np.array([int((labels == label).sum()) for label in order], dtype=int)
    order = [label for label, count in sorted(zip(order, counts), key=lambda item: (-item[1], item[0]))]
    counts = np.array([int((labels == label).sum()) for label in order], dtype=int)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    colors = categorical_colors(order)
    fig, ax = new_fig(max(4.6, 0.62 * len(order) + 2.2), 3.9)
    positions = np.arange(len(order))
    ax.bar(positions, counts, color=[colors[label] for label in order],
           edgecolor="#3c4048", linewidth=0.5)
    for position, count in zip(positions, counts):
        ax.text(position, count, str(int(count)), ha="center", va="bottom", fontsize=7)
    ax.set_xticks(positions)
    ax.set_xticklabels(order, rotation=0 if len(order) <= 12 else 45, fontsize=8)
    ax.set_xlabel("cluster")
    ax.set_ylabel("spots")
    ax.set_title(f"Spots per cluster (n={int(counts.sum())})")
    ax.margins(y=0.14)
    return _finalize(fig, outdir, "cluster_bar", "Cluster sizes", dpi, params)


# ---------------------------------------------------------------------------
# 6. spatial_scatter
# ---------------------------------------------------------------------------

def plot_spatial_scatter(ctx, params, outdir):
    """空间散点：可按簇/域/分组/任意 sample_info 列，或按基因表达连续上色。"""
    coords = _coords(ctx)
    resolved = _resolve_color(ctx, params)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    point_size = _as_float(params, "point_size", DEFAULT_POINT_SIZE)
    index = _subsample(coords.shape[0], _as_int(params, "max_points", DEFAULT_MAX_POINTS))
    sub = _sub_ctx(ctx, index)
    fig, ax = new_fig(6.0, 5.4)
    image = _as_bool(params, "image", False)
    if resolved["kind"] == "num":
        scatter_spatial(ax, sub, values=np.asarray(resolved["values"], dtype=np.float64)[index],
                        point_size=point_size, cmap="viridis", legend_title=resolved["title"],
                        image=image)
        ax.set_title(f"Spatial expression - {resolved['title']}", fontsize=10)
    else:
        scatter_spatial(ax, sub, labels=np.asarray(resolved["labels"], dtype=object)[index],
                        point_size=point_size, legend_title=resolved["title"], image=image)
        ax.set_title(f"Spatial map - {resolved['title']}", fontsize=10)
    return _finalize(fig, outdir, "spatial_scatter", "Spatial scatter", dpi, params)


# ---------------------------------------------------------------------------
# color_by 解析
# ---------------------------------------------------------------------------

def _resolve_color(ctx, params) -> dict:
    """把 color_by/gene 解析成 {'kind':'cat|num','labels'|'values','title'}。

    支持 clusters / domains / group / gene，以及 sample_info 里的任意列名；
    默认 clusters，无聚类时给 prerequisite("先跑聚类…")。
    """
    color_by = str(_param(params, "color_by", "") or "").strip()
    lowered = color_by.lower()
    gene = _param(params, "gene", None)
    if isinstance(gene, (list, tuple)) and gene:
        gene = gene[0]
    gene = str(gene).strip() if gene else ""

    if gene and (not color_by or lowered in ("gene", "expression")):
        index = _gene_index(ctx)
        if gene not in index:
            raise bad_request(f"基因不存在于表达矩阵：{gene}")
        matrix = _X(ctx)
        values = np.asarray(matrix[:, index[gene]], dtype=np.float64).ravel()
        return {"kind": "num", "values": values, "title": f"{gene} expression", "gene": gene}

    if lowered in ("", "clusters", "cluster", "clustering"):
        labels = _cluster_labels(ctx)
        return {"kind": "cat", "labels": labels, "title": "clusters"}

    if lowered in ("domains", "domain"):
        labels = _info_column(ctx, "domains")
        if labels is None:
            results = (ctx or {}).get("results") or {}
            domains = results.get("domains") if isinstance(results.get("domains"), dict) else None
            if isinstance(domains, dict) and domains.get("labels") is not None:
                labels = np.asarray([str(v) for v in domains["labels"]], dtype=object)
        if labels is None:
            raise prerequisite("先跑空间域（domains 步骤），或改用 color_by=clusters")
        return {"kind": "cat", "labels": labels, "title": "domains"}

    if lowered in ("group", "groups", "group_key", "slice", "slices"):
        key = str((ctx or {}).get("group_key") or "group")
        column = key if lowered in ("group", "groups", "group_key") else color_by
        labels = _info_column(ctx, column) or _info_column(ctx, "group")
        if labels is None:
            raise prerequisite("样本信息里没有分组列（group），请先在样本表里设置分组")
        return {"kind": "cat", "labels": labels, "title": column}

    labels = _info_column(ctx, color_by)
    if labels is not None:
        return {"kind": "cat", "labels": labels, "title": color_by}
    index = _gene_index(ctx)
    if color_by in index:
        matrix = _X(ctx)
        values = np.asarray(matrix[:, index[color_by]], dtype=np.float64).ravel()
        return {"kind": "num", "values": values, "title": f"{color_by} expression", "gene": color_by}
    raise bad_request(
        f"未知着色方式：{color_by}（可用：clusters/domains/group/gene，或 sample_info 列名）")


# ---------------------------------------------------------------------------
# 基因选择
# ---------------------------------------------------------------------------

def _parse_gene_list(raw) -> list[str]:
    if raw is None:
        return []
    if isinstance(raw, str):
        items = raw.replace(",", " ").replace(";", " ").replace("\n", " ").split()
        return [item.strip() for item in items if item.strip()]
    if isinstance(raw, (list, tuple)):
        return [str(item).strip() for item in raw if str(item).strip()]
    return []


def _variable_genes(ctx, limit: int) -> list[str]:
    """兜底基因列表：先用 hvg 结果，其次按方差取 top。"""
    genes = _genes(ctx)
    results = (ctx or {}).get("results") or {}
    hvg = results.get("hvg")
    if isinstance(hvg, (list, tuple)) and hvg and isinstance(hvg[0], (int, np.integer)):
        picked = [genes[int(i)] for i in hvg if 0 <= int(i) < len(genes)]
        if picked:
            return picked[:limit]
    cached = results.get("hvg") if isinstance(results.get("hvg"), dict) else {}
    if isinstance(cached, dict):
        names = cached.get("genes") or cached.get("hvg_genes")
        if isinstance(names, (list, tuple)) and names:
            return [str(g) for g in names][:limit]
    matrix = _X(ctx)
    n_genes = matrix.shape[1]
    step = max(1, int(2_000_000 / max(matrix.shape[0], 1)))
    variances = np.zeros(n_genes, dtype=np.float64)
    for start in range(0, n_genes, step):
        chunk = np.asarray(matrix[:, start:start + step], dtype=np.float64)
        variances[start:start + step] = chunk.var(axis=0)
    order = np.argsort(-variances, kind="stable")
    return [genes[i] for i in order[:limit]]


def _pick_genes(ctx, params, limit: int) -> tuple[list[str], np.ndarray]:
    """params.genes/gene → 有效基因列表 + 列索引；为空时退回高变基因/高方差基因。"""
    raw = _param(params, "genes", None)
    if raw is None:
        raw = _param(params, "gene", None)
    requested = _parse_gene_list(raw)
    index = _gene_index(ctx)
    genes = [gene for gene in requested if gene in index]
    missing = [gene for gene in requested if gene not in index]
    if requested and not genes:
        raise bad_request(f"基因都不在表达矩阵中：{', '.join(missing[:8])}")
    if not genes:
        genes = _variable_genes(ctx, limit)[:limit]
    if not genes:
        raise prerequisite("没有可用于绘图的基因")
    return genes[:limit], np.array([index[gene] for gene in genes], dtype=np.int64)


# ---------------------------------------------------------------------------
# 7. gene_box
# ---------------------------------------------------------------------------

def _box_groups(ctx, params) -> tuple[list[str], np.ndarray | None, str]:
    """箱线图分组：优先 clusters，其次 group_col / group_key，没有就单组。"""
    group_by = str(_param(params, "group_by", "") or "").strip()
    lowered = group_by.lower()
    try:
        labels = _cluster_labels(ctx)
        if not group_by or lowered in ("clusters", "cluster"):
            return _sort_labels(labels), np.asarray(labels, dtype=object), "cluster"
    except EngineError:
        labels = None
    column = group_by or str((ctx or {}).get("group_key") or "group")
    found = _info_column(ctx, column)
    if found is not None:
        return _sort_labels(found), found, column
    if labels is not None:
        return _sort_labels(labels), np.asarray(labels, dtype=object), "cluster"
    return ["all spots"], None, "all spots"


def plot_gene_box(ctx, params, outdir):
    """基因表达箱线图（按簇/分组；最多 4 个基因，每个一子图）。"""
    genes, columns = _pick_genes(ctx, params, 4)
    matrix = _X(ctx)
    order, labels, group_name = _box_groups(ctx, params)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    colors = palette(len(order))
    fig, axes = new_fig(4.2 * len(genes) + 0.4, 4.1, ncols=len(genes))
    axes = np.atleast_1d(axes).ravel()
    for ax, gene, column in zip(axes, genes, columns):
        values = np.asarray(matrix[:, column], dtype=np.float64)
        if labels is None:
            groups = [values]
        else:
            groups = [values[np.asarray(labels) == group] for group in order]
        positions = np.arange(1, len(groups) + 1)
        boxes = ax.boxplot(groups, positions=positions, patch_artist=True, widths=0.62,
                           showfliers=False,
                           medianprops={"color": "#2b2f36", "linewidth": 1.1})
        for patch, color in zip(boxes["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.78)
            patch.set_edgecolor("#8a919c")
        ax.set_xticks(positions)
        ax.set_xticklabels(order, rotation=45 if len(order) > 6 else 0, fontsize=7)
        ax.set_ylabel("expression")
        ax.set_title(f"{gene} by {group_name}")
    return _finalize(fig, outdir, "gene_box", "Gene expression boxplot", dpi, params)


# ---------------------------------------------------------------------------
# 8. heatmap
# ---------------------------------------------------------------------------

def _hier_order(matrix: np.ndarray) -> np.ndarray:
    """按行/列轮廓做 average-linkage 层次聚类排序（失败时保持原序）。"""
    if matrix.shape[0] < 3:
        return np.arange(matrix.shape[0])
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import pdist

        if matrix.shape[1] < 3:
            return np.arange(matrix.shape[0])
        return np.asarray(leaves_list(linkage(pdist(matrix, metric="euclidean"),
                                              method="average")), dtype=np.int64)
    except Exception:  # noqa: BLE001
        return np.arange(matrix.shape[0])


def plot_gene_spatial(ctx, params, outdir):
    """单基因（或最多 4 个基因）的空间表达图：连续色带 + 可选 H&E 底图。"""
    genes, columns = _pick_genes(ctx, params, 4)
    matrix = _X(ctx)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    point_size = _as_float(params, "point_size", 6.0)
    image_on = _as_bool(params, "image", False)
    coords = _coords(ctx)
    n = len(genes)
    fig, axes = new_fig(4.4 * n + 0.4, 4.2, ncols=n)
    axes = np.atleast_1d(axes).ravel()
    for index, (gene, column) in enumerate(zip(genes, columns)):
        ax = axes[index]
        if image_on:
            try:
                add_image_background(ax, ctx)
            except Exception:  # noqa: BLE001
                pass
        values = np.asarray(matrix[:, column], dtype=float)
        finite = np.isfinite(values)
        if not finite.any():
            raise prerequisite(f"基因 {gene} 全为 NaN，无法画空间图")
        scatter = ax.scatter(coords[finite, 0], coords[finite, 1], c=values[finite], s=point_size,
                             cmap="viridis", linewidths=0, vmin=float(np.nanmin(values)),
                             vmax=float(np.nanmax(values)))
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(gene, fontsize=10)
        bar = fig.colorbar(scatter, ax=ax, fraction=0.04, pad=0.02)
        bar.ax.tick_params(labelsize=7)
        bar.set_label("expression", fontsize=8)
    return _finalize(fig, outdir, "gene_spatial", "Gene expression in space", dpi, params)


def plot_heatmap(ctx, params, outdir):
    """基因表达热图：行 z-score；cluster=true（默认）时按簇均值，否则按 spot（可下采样）。

    * 行（基因）按 z-score 轮廓层次聚类重排，列（簇或 spot）同法排序。
    * 色带 RdBu_r，vmin/vmax=±2.5，单位 = 行 z-score。
    """
    genes, columns = _pick_genes(ctx, params, int(_param(params, "max_genes", 60)))
    matrix = _X(ctx)
    gene_matrix = np.asarray(matrix[:, columns], dtype=np.float64).T  # (genes, spots)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    cluster_mode = _as_bool(params, "cluster", True)
    labels = _try_cluster_labels(ctx) if cluster_mode else None
    if labels is not None:
        order = _sort_labels(labels)
        rows = np.vstack([gene_matrix[:, np.asarray(labels) == label].mean(axis=1) for label in order])
        z = zscore_rows(rows.T)  # (genes, clusters)：每行（基因）跨簇 z-score
        z = np.clip(z, -2.5, 2.5)
        row_order = _hier_order(z)
        col_order = np.arange(len(order))
        columns_labels = order
        title = "Gene expression by cluster (row z-score)"
        xlabel = "cluster"
    else:
        index = _subsample(gene_matrix.shape[1], min(_as_int(params, "max_points", 800), 800))
        rows = gene_matrix[:, index]
        z = np.clip(zscore_rows(rows), -2.5, 2.5)
        row_order = _hier_order(z)
        col_order = _hier_order(z.T)
        columns_labels = None
        title = "Gene expression heatmap (row z-score, hierarchical order)"
        xlabel = f"spots (n={len(index)})"
    fig, ax = new_fig(max(5.4, 0.34 * len(genes) * 0.6 + 0.26 * len(columns_labels or []) + 3.0),
                      max(3.4, 0.26 * len(genes) + 1.8))
    image = ax.imshow(z[np.ix_(row_order, col_order)], aspect="auto", cmap="RdBu_r",
                      vmin=-2.5, vmax=2.5)
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels([genes[i] for i in row_order], fontsize=7)
    if columns_labels is not None:
        ax.set_xticks(np.arange(len(columns_labels)))
        ax.set_xticklabels([columns_labels[i] for i in col_order], fontsize=8,
                           rotation=0 if len(columns_labels) <= 12 else 45)
    else:
        ax.set_xticks([])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("gene")
    ax.set_title(title)
    ax.grid(False)
    bar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    bar.set_label("row z-score", fontsize=8)
    return _finalize(fig, outdir, "heatmap", "Expression heatmap", dpi, params)


# ---------------------------------------------------------------------------
# 9. dotplot
# ---------------------------------------------------------------------------

def _group_profile(matrix: np.ndarray, codes: np.ndarray, n_groups: int,
                   columns: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(组均值, 组内表达比例) 两个 n_groups × n_genes 矩阵（只算指定基因列）。"""
    means = np.zeros((n_groups, columns.size), dtype=np.float64)
    pct = np.zeros((n_groups, columns.size), dtype=np.float64)
    for code in range(n_groups):
        mask = codes == code
        if not mask.any():
            continue
        sub = np.asarray(matrix[np.ix_(mask, columns)], dtype=np.float64)
        means[code] = sub.mean(axis=0)
        pct[code] = (sub > 0).mean(axis=0)
    return means, pct


def _cluster_codes(labels: np.ndarray) -> tuple[list[str], np.ndarray]:
    order = _sort_labels(labels)
    mapping = {label: index for index, label in enumerate(order)}
    codes = np.array([mapping[str(value)] for value in labels], dtype=np.int64)
    return order, codes


def _auto_marker_genes(ctx, order, codes, top_per_cluster: int, limit: int) -> list[str]:
    """无 genes 参数时：用各簇「簇均值 - 全局均值」挑 top 基因（内置实现）。"""
    matrix = _X(ctx)
    genes = _genes(ctx)
    n_genes = matrix.shape[1]
    global_mean = np.asarray(matrix.mean(axis=0), dtype=np.float64).ravel()
    step = max(1, int(2_000_000 / max(matrix.shape[0], 1)))
    best: list[tuple[int, float, int]] = []
    for code in range(len(order)):
        mask = codes == code
        if not mask.any():
            continue
        scores = np.full(n_genes, -np.inf, dtype=np.float64)
        for start in range(0, n_genes, step):
            stop = min(start + step, n_genes)
            chunk = np.asarray(matrix[np.ix_(mask, np.arange(start, stop))], dtype=np.float64)
            scores[start:stop] = chunk.mean(axis=0) - global_mean[start:stop]
        scores[global_mean <= 0] = -np.inf
        take = np.argsort(-scores, kind="stable")[:top_per_cluster]
        for index in take:
            if np.isfinite(scores[index]) and scores[index] > 0:
                best.append((code, float(scores[index]), int(index)))
    picked: list[str] = []
    for _code, _score, index in sorted(best, key=lambda item: (item[0], -item[1])):
        gene = genes[index]
        if gene not in picked:
            picked.append(gene)
        if len(picked) >= limit:
            break
    return picked


def plot_dotplot(ctx, params, outdir):
    """簇 × 基因点图：点大小 = 该簇表达细胞比例，颜色 = 基因内 z-score 的平均表达。"""
    labels = _cluster_labels(ctx)
    order, codes = _cluster_codes(labels)
    top_per_cluster = int(_param(params, "top_n", 3))
    limit = int(_param(params, "max_genes", 30))
    requested = _parse_gene_list(_param(params, "genes", None)) or _parse_gene_list(
        _param(params, "gene", None))
    index = _gene_index(ctx)
    genes = [gene for gene in requested if gene in index]
    if not genes:
        genes = _auto_marker_genes(ctx, order, codes, top_per_cluster, limit)
    if not genes:
        raise prerequisite("没有可用于点图的基因（可显式传 genes）")
    genes = genes[:limit]
    columns = np.array([index[gene] for gene in genes], dtype=np.int64)
    matrix = _X(ctx)
    means, pct = _group_profile(matrix, codes, len(order), columns)
    scaled = zscore_rows(means.T).T  # 每个基因内跨簇 z-score
    scaled = np.clip(scaled, -2.5, 2.5)
    dpi = _as_int(params, "dpi", DEFAULT_DPI)
    fig, ax = new_fig(max(5.8, 0.30 * len(genes) + 2.6), max(3.0, 0.36 * len(order) + 1.8))
    xs, ys, sizes, colors = [], [], [], []
    for row, label in enumerate(order):
        for col in range(len(genes)):
            xs.append(col)
            ys.append(row)
            sizes.append(4.0 + 34.0 * float(pct[row, col]))
            colors.append(float(scaled[row, col]))
    scatter = ax.scatter(xs, ys, s=sizes, c=colors, cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                         edgecolor="#4a4f57", linewidth=0.35)
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    ax.set_xlabel("gene")
    ax.set_ylabel("cluster")
    ax.set_title("Dot plot (size = fraction expressing, color = scaled mean)")
    ax.grid(True, color="#f0f1f4", linewidth=0.5)
    bar = fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.02)
    bar.set_label("scaled mean expression", fontsize=8)
    return _finalize(fig, outdir, "dotplot", "Dot plot", dpi, params)


# ---------------------------------------------------------------------------
# PLOTS 注册表
# ---------------------------------------------------------------------------

PLOTS = {
    "qc_spatial": plot_qc_spatial,
    "qc_violin": plot_qc_violin,
    "pca_variance": plot_pca_variance,
    "embedding": plot_embedding,
    "cluster_bar": plot_cluster_bar,
    "spatial_scatter": plot_spatial_scatter,
    "gene_box": plot_gene_box,
    "gene_spatial": plot_gene_spatial,
    "heatmap": plot_heatmap,
    "dotplot": plot_dotplot,
}
