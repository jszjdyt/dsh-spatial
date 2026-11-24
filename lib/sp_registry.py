"""
dsh-spatial — 模块注册表与调度层

职责（DESIGN.md §2.3 / §2.5）：
  * 动态导入全部步骤/分析/绘图模块（sp_prep / sp_plots / sp_spatial / sp_domain /
    sp_deconv / sp_lr）；**某个模块文件缺失或 import 失败只跳过并记 modules[name].error，
    绝不拖垮整体**（增量交付可用）；
  * 汇总 STEPS / ANALYSES 目录（含 needs 可用性判断）与 **PLOTS**（每模块可导出
    `PLOTS = {kind: fn}`，合并成全局 kind→fn 表）；
  * 构造 ctx（§2.3）、执行步骤与分析、消费模块返回的写回键
    （replace / store / cluster_store / coords_store / sample_columns+sample_values）；
  * 数据集构造（sp_io.load_source → sp_state.Dataset）、预览、样本/切片表、
    自检 audit、建议 suggest、预设流程 recipe、导出 export；
  * 单细胞参考（h5ad）探测 / 元信息 / 挂载。

本文件只被 sp_engine.py 调用；步骤/分析模块只依赖 sp_common，不 import 本文件（避免循环）。
显式依赖 sp_state / sp_io / sp_geo（由 sp-core 交付）：只在真正用到时 import，
缺失时抛出中文可读的 prerequisite，插件本身照常启动。
"""

from __future__ import annotations

import importlib
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import sp_common
from sp_common import EngineError, bad_request, ensure_dir, jsonable, prerequisite, result, table

PLUGIN = "dsh-spatial"

# 模块注册顺序：先基础流程/图，再空间统计、域、反卷积、配体受体
MODULES = ["sp_prep", "sp_plots", "sp_spatial", "sp_domain", "sp_deconv", "sp_lr"]

# ---------------------------------------------------------------------------
# 目录常量（DESIGN.md §3 步骤顺序 / §5 图白名单 / §6 格式）
# ---------------------------------------------------------------------------

# 12 个流程步骤：UI 左栏固定顺序（模块缺失时仍显示为「未就绪」）
STEP_ORDER = [
    {"id": "qc", "title": "空间质控", "requires": [], "params": [
        {"key": "min_counts", "label": "最少 counts", "type": "number", "default": None},
        {"key": "min_genes", "label": "最少基因数", "type": "number", "default": None},
        {"key": "max_pct_mt", "label": "线粒体比例上限(%)", "type": "number", "default": 20},
        {"key": "group_col", "label": "分组列", "type": "text", "default": None},
    ]},
    {"id": "filter", "title": "过滤", "requires": ["qc"], "params": [
        {"key": "min_counts", "label": "最少 counts", "type": "number", "default": 500},
        {"key": "min_genes", "label": "最少基因数", "type": "number", "default": 200},
        {"key": "max_pct_mt", "label": "线粒体比例上限(%)", "type": "number", "default": 20},
        {"key": "min_cells", "label": "基因至少在几个 spot 中表达", "type": "number", "default": 3},
    ]},
    {"id": "normalize", "title": "标准化", "requires": ["filter"], "params": [
        {"key": "target_sum", "label": "目标总量", "type": "number", "default": 10000},
        {"key": "log", "label": "log1p 转换", "type": "boolean", "default": True},
        {"key": "method", "label": "方法", "type": "select", "default": "lognorm",
         "options": ["lognorm", "cpm", "none"]},
    ]},
    {"id": "hvg", "title": "高变基因", "requires": ["normalize"], "params": [
        {"key": "n_top_genes", "label": "高变基因数", "type": "number", "default": 2000},
        {"key": "flavor", "label": "算法", "type": "select", "default": "seurat",
         "options": ["seurat", "vst", "dispersion"]},
    ]},
    {"id": "pca", "title": "PCA", "requires": ["hvg"], "params": [
        {"key": "n_pcs", "label": "主成分数", "type": "number", "default": 30},
        {"key": "scale", "label": "先做 z-score", "type": "boolean", "default": True},
    ]},
    {"id": "neighbors", "title": "邻接图", "requires": ["pca"], "params": [
        {"key": "n_neighbors", "label": "近邻数", "type": "number", "default": 15},
        {"key": "metric", "label": "距离度量", "type": "select", "default": "cosine",
         "options": ["cosine", "euclidean"]},
        {"key": "spatial_weight", "label": "空间权重", "type": "number", "default": 0.3},
        {"key": "spatial_method", "label": "空间邻接", "type": "select", "default": "grid",
         "options": ["grid", "knn", "radius"]},
    ]},
    {"id": "cluster", "title": "聚类", "requires": ["neighbors"], "params": [
        {"key": "resolution", "label": "分辨率", "type": "number", "default": 1.0},
        {"key": "method", "label": "算法", "type": "select", "default": "louvain",
         "options": ["louvain", "leiden_like", "kmeans"]},
    ]},
    {"id": "embed", "title": "二维嵌入", "requires": ["pca"], "params": [
        {"key": "layout", "label": "布局", "type": "select", "default": "tsne", "options": ["tsne", "umap_like", "pca"]},
        {"key": "perplexity", "label": "perplexity", "type": "number", "default": 30},
    ]},
    {"id": "svg", "title": "空间高变基因", "requires": ["cluster"], "params": [
        {"key": "n_top", "label": "返回基因数", "type": "number", "default": 30},
        {"key": "method", "label": "统计量", "type": "select", "default": "moran", "options": ["moran", "geary"]},
        {"key": "permutations", "label": "置换次数", "type": "number", "default": 0},
    ]},
    {"id": "domains", "title": "空间域", "requires": ["cluster"], "params": [
        {"key": "n_domains", "label": "域数量", "type": "number", "default": 6},
        {"key": "method", "label": "方法", "type": "select", "default": "mrf", "options": ["mrf", "kmeans"]},
        {"key": "smooth", "label": "平滑强度", "type": "number", "default": 0.5},
        {"key": "iterations", "label": "迭代次数", "type": "number", "default": 20},
    ]},
    {"id": "deconv", "title": "反卷积", "requires": ["normalize"], "params": [
        {"key": "ref", "label": "参考 h5ad 路径", "type": "text", "default": ""},
        {"key": "method", "label": "方法", "type": "select", "default": "nnls", "options": ["nnls", "dwls"]},
        {"key": "cell_type_col", "label": "细胞类型列", "type": "text", "default": ""},
        {"key": "min_genes", "label": "最少共用基因", "type": "number", "default": 50},
    ]},
    {"id": "annotate", "title": "区域注释", "requires": ["domains"], "params": [
        {"key": "mapping", "label": "域→注释映射", "type": "object", "default": {}},
        {"key": "key", "label": "写回列名", "type": "text", "default": "domain_anno"},
    ]},
]
STEP_INDEX = {item["id"]: i for i, item in enumerate(STEP_ORDER)}
STEP_IDS = [item["id"] for item in STEP_ORDER]

# 图类型白名单（DESIGN.md §5；每个 kind 只能有一个 owner 模块）
PLOT_ORDER = [
    "qc_spatial", "qc_violin", "spatial_scatter", "embedding", "cluster_bar", "gene_spatial",
    "svg_map", "moran_plot", "domain_map", "deconv_spatial", "deconv_bar", "nhood_enrich",
    "cooccurrence", "ripley", "pca_variance", "heatmap", "dotplot", "gene_box",
    "de_volcano", "de_heatmap",
]

# 分析分组顺序（DESIGN.md §4）
GROUP_ORDER = ["空间统计", "邻域与域", "反卷积", "空间通讯", "区域差异", "质控与整合", "其它"]

# sp_io 尚未落盘时的格式清单（DESIGN.md §6），保证 /formats 与 ping 不为空
FORMAT_FALLBACK = [
    {"id": "visium", "label": "10x Visium outs/ 目录", "hint": "含 spatial/tissue_positions*.csv 与 filtered_feature_bc_matrix*"},
    {"id": "visium_h5", "label": "Visium h5 + 坐标", "hint": "单文件 filtered_feature_bc_matrix.h5，需同时给坐标 csv 或同目录 spatial/"},
    {"id": "h5ad", "label": "h5ad（obsm['spatial']）", "hint": "按 obsm['spatial'] 或 obs 的 x/y、pxl_row/pxl_col 取坐标"},
    {"id": "stereo", "label": "Stereo-seq BIN", "hint": "*_raw_bin_count.tsv.gz（geneID,x,y,count 长表）"},
    {"id": "slideseq", "label": "Slide-seq puck", "hint": "*_puck_*.txt / *_bead*.txt + beacon 坐标"},
    {"id": "matrix_text", "label": "表达矩阵 + 坐标表", "hint": "csv/tsv 矩阵，另给 coords（spot,x,y）"},
    {"id": "xlsx", "label": "Excel 矩阵 + 坐标", "hint": "同 matrix_text，可用 sheet/gene_col 指定"},
    {"id": "geo_series", "label": "GEO series matrix", "hint": "含 !series_matrix_table_begin；坐标缺失时只能当表达矩阵"},
    {"id": "archive", "label": "压缩包", "hint": "zip/tar/tar.gz，自动识别 Visium/BIN/矩阵；多候选时给出候选列表"},
    {"id": "rds", "label": "RDS（R 桥）", "hint": "Seurat/SingleCellExperiment → 矩阵 + obsm；R 不可用时给 unsupported"},
    {"id": "rdata", "label": "RData（R 桥）", "hint": "同上"},
]

EXPORT_KINDS = [
    "expr_csv", "svg_csv", "domain_csv", "deconv_csv", "de_csv", "sample_info_csv",
    "summary_json", "protocol", "methods", "report", "script",
]

_RUNTIME: dict = {"modules": {}, "errors": {}, "plots": {}, "plot_errors": {}, "loaded_at": 0.0}


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------

def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


def _require(name: str, hint: str):
    """按需 import 依赖模块（sp_state / sp_io / sp_geo），缺失时给中文可读错误。"""
    try:
        return importlib.import_module(name)
    except Exception as error:  # noqa: BLE001
        raise EngineError(
            "prerequisite",
            f"{hint}尚未就绪：无法导入 {name}.py（{type(error).__name__}: {error}）",
        )


def _as_array(value, dtype=np.float32):
    if value is None:
        return None
    try:
        arr = np.asarray(value, dtype=dtype)
    except Exception:  # noqa: BLE001
        return None
    return arr if arr.ndim >= 1 else None


def _n_spots(dataset) -> int:
    X = getattr(dataset, "X", None)
    if hasattr(X, "shape") and getattr(X, "ndim", 0) == 2:
        return int(X.shape[0])
    return len(getattr(dataset, "spots", None) or [])


def _n_genes(dataset) -> int:
    X = getattr(dataset, "X", None)
    if hasattr(X, "shape") and getattr(X, "ndim", 0) == 2:
        return int(X.shape[1])
    return len(getattr(dataset, "genes", None) or [])


def _spots(dataset) -> list[str]:
    raw = getattr(dataset, "spots", None)
    if raw is None or len(raw) == 0:
        raw = [f"spot{i + 1}" for i in range(_n_spots(dataset))]
    return [str(x) for x in raw]


def _genes(dataset) -> list[str]:
    raw = getattr(dataset, "genes", None)
    if raw is None or len(raw) == 0:
        raw = [f"GENE{i + 1}" for i in range(_n_genes(dataset))]
    return [str(x) for x in raw]


def _coords(dataset):
    coords = _as_array(getattr(dataset, "coords", None))
    if coords is None or coords.ndim != 2 or coords.shape[1] < 2:
        return None
    return coords[:, :2]


def _labels(dataset):
    """当前聚类标签（dataset.labels() 或 sample_info[cluster_key]）。"""
    labels = None
    fn = getattr(dataset, "labels", None)
    if callable(fn):
        try:
            labels = fn()
        except Exception:  # noqa: BLE001
            labels = None
    if labels is None:
        key = getattr(dataset, "cluster_key", None)
        frame = _sample_frame(dataset)
        if key and key in frame.columns:
            labels = frame[key].tolist()
    if labels is None:
        return None
    labels = [str(v) for v in labels]
    if len(labels) != _n_spots(dataset):
        return None
    return labels


def _sample_frame(dataset) -> pd.DataFrame:
    """返回与 spots 对齐的 spot 级表（index=spot id）。"""
    frame = getattr(dataset, "sample_info", None)
    spots = _spots(dataset)
    if not isinstance(frame, pd.DataFrame):
        frame = pd.DataFrame(index=pd.Index(spots, name="spot"))
    if len(frame) != len(spots) or [str(x) for x in frame.index] != spots:
        aligned = pd.DataFrame(index=pd.Index(spots, name=getattr(frame.index, "name", "spot") or "spot"))
        for col in list(frame.columns):
            try:
                aligned[col] = frame[col].reindex(spots).tolist()
            except Exception:  # noqa: BLE001
                continue
        frame = aligned
    return frame


def _column_kind(series: pd.Series) -> str:
    try:
        dtype = str(series.dtype)
        if dtype in ("object", "category", "bool", "string"):
            return "categorical"
        if int(series.nunique(dropna=True)) <= 12 and dtype.startswith("int"):
            return "categorical"
        return "numeric"
    except Exception:  # noqa: BLE001
        return "categorical"


def _sample_columns(dataset, limit: int = 40) -> list[dict]:
    frame = _sample_frame(dataset)
    out = []
    for col in list(frame.columns)[:limit]:
        try:
            out.append({"key": str(col), "kind": _column_kind(frame[col]),
                        "n_unique": int(frame[col].nunique(dropna=True))})
        except Exception:  # noqa: BLE001
            out.append({"key": str(col), "kind": "categorical", "n_unique": None})
    return out


def _guessed_group_key(dataset):
    frame = _sample_frame(dataset)
    for cand in ("group", "condition", "treatment", "disease", "region", "batch", "slice"):
        if cand in frame.columns:
            return cand
    return None


def _groups_of(frame: pd.DataFrame, key) -> list[str] | None:
    if not key or key not in frame.columns:
        return None
    try:
        values = [str(v) for v in frame[key].dropna().unique().tolist()]
    except Exception:  # noqa: BLE001
        return None
    return values[:50]


def _summary(dataset, history_limit: int = 50) -> dict:
    """dataset.summary()（§2.4）；签名/实现异常时给最小可用摘要，绝不抛错。"""
    fn = getattr(dataset, "summary", None)
    if callable(fn):
        for kwargs in ({}, {"history_limit": history_limit}):
            try:
                value = fn(**kwargs)
                if isinstance(value, dict):
                    return jsonable(value)
            except TypeError:
                continue
            except Exception:  # noqa: BLE001
                break
    return _fallback_summary(dataset)


def _fallback_summary(dataset) -> dict:
    coords = _coords(dataset)
    geometry = {
        "mode": getattr(dataset, "spatial_mode", "generic") or "generic",
        "coords_kind": getattr(dataset, "coords_kind", "generic") or "generic",
        "neighbor_radius": getattr(dataset, "neighbor_radius", None),
    }
    if coords is not None:
        geometry["x_range"] = [float(coords[:, 0].min()), float(coords[:, 0].max())]
        geometry["y_range"] = [float(coords[:, 1].min()), float(coords[:, 1].max())]
    image = getattr(dataset, "image", None) or {}
    geometry["image"] = image.get("path") if isinstance(image, dict) else image
    return jsonable({
        "id": getattr(dataset, "id", "ds1"),
        "name": getattr(dataset, "name", "dataset"),
        "path": getattr(dataset, "path", None),
        "source_format": getattr(dataset, "source_format", "unknown"),
        "n_spots": _n_spots(dataset),
        "n_genes": _n_genes(dataset),
        "species": getattr(dataset, "species", "unknown"),
        "geometry": geometry,
        "sample_columns": _sample_columns(dataset),
        "steps": list(getattr(dataset, "steps", None) or []),
        "current": getattr(dataset, "current", None),
        "cluster_keys": _cluster_keys(dataset),
        "ref": getattr(dataset, "ref", None),
        "warnings": list(getattr(dataset, "warnings", None) or [])[:20],
        "notes": list(getattr(dataset, "notes", None) or [])[:20],
        "history": list(getattr(dataset, "history", None) or [])[-history_limit:],
        "analyses_run": list(getattr(dataset, "analyses_run", None) or [])[-20:],
    })


def _cluster_keys(dataset) -> list[str]:
    frame = _sample_frame(dataset)
    keys = []
    current = getattr(dataset, "cluster_key", None)
    if current:
        keys.append(str(current))
    for cand in ("clusters", "cluster", "domains", "domain", "domain_anno", "cell_type", "annotation", "niche"):
        if cand in frame.columns and cand not in keys:
            keys.append(cand)
    return keys


# ---------------------------------------------------------------------------
# 模块加载与目录
# ---------------------------------------------------------------------------

def _norm_params(params) -> list[dict]:
    out = []
    for item in params or []:
        if not isinstance(item, dict):
            continue
        param = dict(item)
        param.setdefault("key", "")
        param.setdefault("label", param["key"])
        param.setdefault("type", "text")
        param.setdefault("default", None)
        out.append(jsonable(param))
    return out


def load_modules(force: bool = False) -> dict:
    """导入所有模块，收集 STEPS / ANALYSES / PLOTS；单模块失败不影响其它模块。"""
    if not force and _RUNTIME["modules"]:
        return _RUNTIME
    modules: dict[str, dict] = {}
    errors: dict[str, str] = {}
    plots: dict[str, dict] = {}
    plot_errors: dict[str, str] = {}
    for name in MODULES:
        # 不跳过上次失败的模块：文件可能刚落盘，再试一次（失败的 import 不会被缓存）
        try:
            mod = importlib.import_module(name)
        except Exception as error:  # noqa: BLE001
            modules[name] = {"steps": [], "analyses": [], "plots": [], "ok": False,
                             "error": f"{type(error).__name__}: {error}"}
            errors[name] = modules[name]["error"]
            continue
        steps = []
        for item in getattr(mod, "STEPS", []) or []:
            if not isinstance(item, dict) or "id" not in item or "run" not in item:
                continue
            steps.append({
                "id": str(item["id"]),
                "title": str(item.get("title") or item["id"]),
                "icon": str(item.get("icon") or item["id"]),
                "requires": [str(x) for x in (item.get("requires") or [])],
                "params": _norm_params(item.get("params") or []),
                "run": item["run"],
                "module": name,
            })
        analyses = []
        for item in getattr(mod, "ANALYSES", []) or []:
            if not isinstance(item, dict) or "id" not in item or "run" not in item:
                continue
            analyses.append({
                "id": str(item["id"]),
                "title": str(item.get("title") or item["id"]),
                "group": str(item.get("group") or "其它"),
                "method": str(item.get("method") or ""),
                "needs": [str(x) for x in (item.get("needs") or [])],
                "params": _norm_params(item.get("params") or []),
                "run": item["run"],
                "module": name,
            })
        module_plots = []
        raw_plots = getattr(mod, "PLOTS", {}) or {}
        if isinstance(raw_plots, dict):
            for kind, fn in raw_plots.items():
                kind = str(kind)
                if not callable(fn):
                    continue
                module_plots.append(kind)
                if kind in plots:
                    plot_errors[kind] = f"重复注册：{plots[kind]['module']} 与 {name}"
                    continue
                plots[kind] = {"kind": kind, "fn": fn, "module": name}
        modules[name] = {"steps": steps, "analyses": analyses, "plots": module_plots,
                         "ok": True, "error": None}
    _RUNTIME["modules"] = modules
    _RUNTIME["errors"] = errors
    _RUNTIME["plots"] = plots
    _RUNTIME["plot_errors"] = plot_errors
    _RUNTIME["loaded_at"] = time.time()
    return _RUNTIME


def module_status() -> dict:
    load_modules()
    return {
        name: {"ok": bool(info.get("ok")), "error": info.get("error"),
               "steps": [s["id"] for s in info.get("steps", [])],
               "analyses": [a["id"] for a in info.get("analyses", [])],
               "plots": list(info.get("plots", []))}
        for name, info in _RUNTIME["modules"].items()
    }


def all_steps() -> dict[str, dict]:
    found: dict[str, dict] = {}
    for info in load_modules()["modules"].values():
        for step in info.get("steps", []):
            found[step["id"]] = step
    return found


def all_analyses() -> dict[str, dict]:
    found: dict[str, dict] = {}
    for info in load_modules()["modules"].values():
        for analysis in info.get("analyses", []):
            found.setdefault(analysis["id"], analysis)
    return found


def all_plots() -> dict[str, dict]:
    load_modules()
    return dict(_RUNTIME["plots"])


def plot_kinds() -> list[str]:
    """已实现 + 契约白名单（保持 DESIGN §5 顺序，便于 UI 稳定渲染）。"""
    available = set(all_plots().keys())
    ordered = [kind for kind in PLOT_ORDER if kind in available]
    extra = sorted(available - set(PLOT_ORDER))
    return ordered + extra


def steps_catalog() -> list[dict]:
    """UI 左栏流程卡：固定 12 项 + 模块提供的参数与就绪状态。"""
    implemented = all_steps()
    out = []
    for base in STEP_ORDER:
        item = dict(base)
        impl = implemented.get(base["id"])
        if impl:
            item["params"] = impl["params"] or item["params"]
            item["requires"] = impl["requires"] or item["requires"]
            item["available"] = True
            item["module"] = impl["module"]
        else:
            item["available"] = False
            item["module"] = None
            item["reason"] = "该步骤的模块尚未就绪"
        out.append(jsonable(item))
    return out


def analyses_catalog(dataset) -> dict:
    found = all_analyses()
    items = []
    groups: list[str] = []
    for analysis in found.values():
        ok, reason = needs_ok(dataset, analysis["needs"])
        if analysis["group"] not in groups:
            groups.append(analysis["group"])
        items.append({
            "id": analysis["id"],
            "title": analysis["title"],
            "group": analysis["group"],
            "method": analysis["method"],
            "needs": analysis["needs"],
            "params": analysis["params"],
            "available": bool(ok),
            "reason": None if ok else reason,
            "native": True,
            "module": analysis["module"],
        })
    groups.sort(key=lambda g: (GROUP_ORDER.index(g) if g in GROUP_ORDER else len(GROUP_ORDER), g))
    items.sort(key=lambda it: (groups.index(it["group"]) if it["group"] in groups else 99, it["title"]))
    return {
        "analyses": jsonable(items),
        "groups": groups,
        "diagnostics": _diagnostics(dataset),
        "columns": _sample_columns(dataset),
    }


def _diagnostics(dataset) -> dict:
    coords = _coords(dataset)
    return jsonable({
        "n_spots": _n_spots(dataset),
        "n_genes": _n_genes(dataset),
        "value_kind": getattr(dataset, "X_kind", None) or getattr(dataset, "value_kind", "unknown"),
        "has_coords": coords is not None,
        "spatial_mode": getattr(dataset, "spatial_mode", None),
        "neighbor_radius": getattr(dataset, "neighbor_radius", None),
        "cluster_keys": _cluster_keys(dataset),
        "group_key": getattr(dataset, "group_key", None),
        "groups": getattr(dataset, "groups", None),
        "has_ref": bool(getattr(dataset, "ref", None)),
        "steps": list(getattr(dataset, "steps", None) or []),
    })


class _EmptyDataset:
    """ping 时计算分析可用性用（没有数据集 → 所有 needs 都不满足）。"""

    id = ""
    name = ""
    spots = []
    genes = []
    X = None
    coords = None
    sample_info = None
    results: dict = {}
    steps: list = []
    group_key = None
    groups = None
    cluster_key = None
    ref = None
    image = None
    spatial_mode = None

    def summary(self, **kwargs):
        return {}

    def labels(self):
        return None


def needs_ok(dataset, needs) -> tuple[bool, str | None]:
    """DESIGN.md §2.5 的语义前置判断：needs ∈ coords/x/y | normalize | cluster |
    domains | groups | ref | image（外加 hvg/pca/neighbors/filter/qc 等步骤名）。"""
    normalized = [str(n).strip() for n in (needs or []) if str(n).strip()]
    if not normalized:
        return True, None
    steps = set(getattr(dataset, "steps", None) or [])
    for need in normalized:
        low = need.lower()
        if low in ("coords", "coordinate", "x/y", "xy", "spatial", "coords_xy"):
            if _coords(dataset) is None:
                return False, "该数据集没有空间坐标，请先导入带坐标的数据（Visium/Stereo-seq/h5ad 的 obsm['spatial']）"
            continue
        if low == "normalize":
            kind = str(getattr(dataset, "X_kind", None) or getattr(dataset, "value_kind", "") or "")
            if "normalize" not in steps and kind not in ("lognorm", "lognormal", "cpm", "tpm", "fpkm"):
                return False, "需要先跑「标准化」步骤"
            continue
        if low == "cluster":
            if _labels(dataset) is None:
                return False, "需要先跑「聚类」步骤（得到 cluster 标签）"
            continue
        if low == "domains":
            frame = _sample_frame(dataset)
            if "domains" not in frame.columns and getattr(dataset, "cluster_key", None) != "domains":
                return False, "需要先跑「空间域」步骤（得到 domains 标签）"
            continue
        if low == "groups":
            if not getattr(dataset, "group_key", None):
                return False, "需要先指定分组列（样本/切片表的 group 列）"
            continue
        if low == "ref":
            if not getattr(dataset, "ref", None):
                return False, "需要先指定单细胞参考 h5ad（ref_attach 或 deconv 参数 ref）"
            continue
        if low == "image":
            image = getattr(dataset, "image", None) or {}
            path = image.get("path") if isinstance(image, dict) else image
            if not path or not os.path.exists(str(path)):
                return False, "该数据集没有组织图像（H&E），无法叠加底图"
            continue
        if low in ("hvg", "pca", "neighbors", "filter", "qc", "svg", "deconv", "embed"):
            if low not in steps:
                titles = {item["id"]: item["title"] for item in STEP_ORDER}
                return False, f"需要先跑「{titles.get(low, low)}」步骤"
            continue
        # 未知 need：不阻塞（模块自己会在 run 里抛 prerequisite）
        continue
    return True, None


# ---------------------------------------------------------------------------
# ctx（DESIGN.md §2.3）
# ---------------------------------------------------------------------------

class LazyRef(dict):
    """ctx["ref"]：X/genes/obs 第一次被取用时才读 h5ad（参考可能很大）。"""

    def __init__(self, path, cell_type_col=None, cell_types=None):
        super().__init__(path=str(path), cell_type_col=cell_type_col,
                         cell_types=list(cell_types or []), X=None, genes=None)
        self._loaded = False

    def _ensure(self):
        if self._loaded:
            return
        self._loaded = True
        path = str(self.get("path") or "")
        try:
            import anndata

            adata = anndata.read_h5ad(path)
            X = adata.X
            if hasattr(X, "toarray"):
                X = X.toarray()
            self["X"] = np.asarray(X, dtype=np.float32)
            self["genes"] = [str(g) for g in adata.var_names]
            self["obs"] = adata.obs
            if not self["cell_types"]:
                col = self.get("cell_type_col")
                if not col or str(col) not in adata.obs.columns:
                    col = next((c for c in ("cell_type", "celltype", "cell_type_l1", "annotation",
                                            "clusters", "leiden", "seurat_clusters", "louvain")
                                if c in adata.obs.columns), None)
                if col:
                    self["cell_type_col"] = str(col)
                    self["cell_types"] = [str(v) for v in
                                          list(dict.fromkeys(adata.obs[col].astype(str).tolist()))][:100]
        except Exception as error:  # noqa: BLE001
            self["error"] = f"{type(error).__name__}: {error}"

    def __getitem__(self, key):
        if key in ("X", "genes", "obs"):
            self._ensure()
        return super().__getitem__(key)

    def get(self, key, default=None):
        if key in ("X", "genes", "obs"):
            self._ensure()
        return super().get(key, default)

    def materialize(self) -> "LazyRef":
        self._ensure()
        return self


def _ref_payload(dataset) -> dict | None:
    ref = getattr(dataset, "ref", None)
    if not ref:
        return None
    if isinstance(ref, str):
        return {"path": ref, "cell_type_col": None, "cell_types": []}
    if isinstance(ref, dict):
        path = ref.get("path")
        if not path:
            return None
        return {"path": str(path), "cell_type_col": ref.get("cell_type_col"),
                "cell_types": list(ref.get("cell_types") or [])}
    return None


def build_ctx(dataset, outdir: str, engine_info: dict | None = None) -> dict:
    ensure_dir(outdir)
    ref_meta = _ref_payload(dataset)
    coords = _coords(dataset)
    image = getattr(dataset, "image", None)
    if isinstance(image, str):
        image = {"path": image, "extent": None, "scale": None}
    if not isinstance(image, dict):
        image = {"path": None, "extent": None, "scale": None}
    runtime_workdir = getattr(dataset, "_workdir", None) or os.path.dirname(
        os.path.dirname(outdir.rstrip(os.sep))) or outdir
    return {
        "dataset_id": getattr(dataset, "id", "ds1"),
        "name": getattr(dataset, "name", "dataset"),
        "genes": list(_genes(dataset)),
        "spots": list(_spots(dataset)),
        "X": getattr(dataset, "X", None),
        "X_kind": str(getattr(dataset, "X_kind", None) or getattr(dataset, "value_kind", "unknown") or "unknown"),
        "counts": getattr(dataset, "counts", None),
        "coords": coords,
        "coords_kind": str(getattr(dataset, "coords_kind", "generic") or "generic"),
        "image": image,
        "spatial_mode": str(getattr(dataset, "spatial_mode", "generic") or "generic"),
        "neighbor_radius": getattr(dataset, "neighbor_radius", None),
        "cluster_key": getattr(dataset, "cluster_key", None),
        "cluster": _labels(dataset),
        "sample_info": _sample_frame(dataset),
        "group_key": getattr(dataset, "group_key", None),
        "groups": getattr(dataset, "groups", None),
        "species": getattr(dataset, "species", "unknown") or "unknown",
        "obsm": getattr(dataset, "obsm", None) if isinstance(getattr(dataset, "obsm", None), dict) else {},
        "results": getattr(dataset, "results", None) if isinstance(getattr(dataset, "results", None), dict) else {},
        "ref": LazyRef(**ref_meta) if ref_meta else None,
        "rng": np.random.default_rng(0),
        "workdir": runtime_workdir,
        "outdir": outdir,
        "engine": engine_info or {},
        "notes": list(getattr(dataset, "notes", None) or []),
        "warnings": list(getattr(dataset, "warnings", None) or []),
    }


def _materialize_ref(ctx) -> None:
    ref = ctx.get("ref")
    if isinstance(ref, LazyRef):
        try:
            ref.materialize()
        except Exception:  # noqa: BLE001 - 参考读不了时让模块自己报 prerequisite
            pass


# ---------------------------------------------------------------------------
# 写回（§2.5 的 replace / store / cluster_store / coords_store / sample_columns）
# ---------------------------------------------------------------------------

_ATTR_ALIASES = {"value_kind": "X_kind", "kind": "X_kind", "samples": "spots"}


#: replace 里 sp_state.apply_replace 不认识的键（由 registry 直接写属性）
_REPLACE_EXTRA = ("spatial_mode", "neighbor_radius", "cluster_key", "species", "image",
                  "ref", "name", "path", "source_format", "group_key", "obsm")


def _apply_replace(dataset, replace: dict) -> None:
    """步骤 replace：替换矩阵/基因/spot/counts/坐标/value_kind（优先用 sp_state 的实现）。"""
    if not isinstance(replace, dict) or not replace:
        return
    fn = getattr(dataset, "apply_replace", None)
    if callable(fn):
        try:
            fn(replace)
        except EngineError:
            raise
        except Exception as error:  # noqa: BLE001
            raise EngineError("internal", f"写回 replace 失败：{type(error).__name__}: {error}")
        for key in _REPLACE_EXTRA:
            if key in replace and replace[key] is not None:
                try:
                    setattr(dataset, key, replace[key])
                except Exception:  # noqa: BLE001 - groups 这类只读属性跳过
                    continue
        return
    # 兼容实现（sp_state 未提供 apply_replace 时）
    for key, value in replace.items():
        attr = _ATTR_ALIASES.get(str(key), str(key))
        try:
            if attr == "X":
                dataset.X = np.asarray(value, dtype=np.float32)
            elif attr == "counts":
                dataset.counts = None if value is None else np.asarray(value, dtype=np.float32)
            elif attr == "coords":
                dataset.coords = None if value is None else np.asarray(value, dtype=np.float32)
            elif attr == "genes":
                dataset.genes = [str(g) for g in value]
            elif attr == "spots":
                dataset.spots = [str(s) for s in value]
                dataset.sample_info = _sample_frame(dataset)
            elif attr == "sample_info":
                if isinstance(value, pd.DataFrame):
                    dataset.sample_info = value
            elif attr != "groups":
                setattr(dataset, attr, value)
        except Exception:  # noqa: BLE001 - 单个写回键失败不能拖垮整步
            continue


def _apply_sample_columns(dataset, columns, values) -> tuple[list[str], list[str]]:
    warnings: list[str] = []
    written: list[str] = []
    n = _n_spots(dataset)
    if not isinstance(values, dict):
        return written, warnings
    accepted: dict[str, list] = {}
    for col in columns or []:
        key = str(col.get("key") if isinstance(col, dict) else col)
        if not key:
            continue
        vals = values.get(key)
        if vals is None:
            continue
        try:
            vals = list(vals.tolist() if isinstance(vals, np.ndarray) else vals)
        except Exception:  # noqa: BLE001
            continue
        if len(vals) != n:
            warnings.append(f"模块返回的 {key} 列长度（{len(vals)}）与 spot 数（{n}）不一致，已跳过写回")
            continue
        accepted[key] = vals
    if not accepted:
        return written, warnings
    fn = getattr(dataset, "apply_sample_columns", None)
    if callable(fn):
        try:
            fn([{"key": key} for key in accepted], accepted)
            written = list(accepted)
            return written, warnings
        except Exception:  # noqa: BLE001 - 回退到直接写
            written = []
    for key, vals in accepted.items():
        if _set_spot_column(dataset, key, vals):
            written.append(key)
    return written, warnings


def _set_spot_column(dataset, key: str, values, kind: str = "categorical") -> bool:
    if not key or key == "__noop__":
        return False
    spots = _spots(dataset)
    frame = _sample_frame(dataset)
    try:
        values = [jsonable(v) for v in list(values)]
    except Exception:  # noqa: BLE001
        return False
    if len(values) != len(spots):
        return False
    try:
        frame[key] = pd.Series(values, index=frame.index, dtype="object")
    except Exception:  # noqa: BLE001
        return False
    dataset.sample_info = frame
    return True


def _apply_cluster_store(dataset, payload) -> str | None:
    if not isinstance(payload, dict):
        return None
    key = str(payload.get("key") or payload.get("cluster_key") or getattr(dataset, "cluster_key", None) or "clusters")
    labels = payload.get("labels")
    if labels is None:
        for cand in ("values", "clusters", "cluster", "labels"):
            if isinstance(payload.get(cand), (list, tuple, np.ndarray)):
                labels = payload[cand]
                break
    if labels is None:
        # 形如 {"clusters": [...], "domains": [...]}
        for cand_key, cand_value in payload.items():
            if isinstance(cand_value, (list, tuple, np.ndarray)):
                key, labels = str(cand_key), cand_value
                break
    if labels is None:
        return None
    if _set_spot_column(dataset, key, [str(v) for v in labels]):
        dataset.cluster_key = key
        return key
    return None


def _apply_coords_store(dataset, payload) -> str | None:
    if not isinstance(payload, dict):
        return None
    obsm = getattr(dataset, "obsm", None)
    if not isinstance(obsm, dict):
        obsm = {}
    written = None
    items = {}
    if isinstance(payload.get("coords"), (list, np.ndarray)):
        items[str(payload.get("key") or payload.get("kind") or "embedding")] = payload["coords"]
    else:
        for key, value in payload.items():
            if isinstance(value, (list, np.ndarray)) and str(key) not in ("key", "kind"):
                items[str(key)] = value
    for key, value in items.items():
        try:
            arr = np.asarray(value, dtype=np.float32)
        except Exception:  # noqa: BLE001
            continue
        if arr.ndim != 2 or arr.shape[0] != _n_spots(dataset):
            continue
        obsm[key] = arr
        written = key
    dataset.obsm = obsm
    return written


def _normalize_store(value):
    """把模块返回的缓存值统一成可长期保存的结构（DataFrame 留内存，list[dict] → DataFrame）。"""
    if not isinstance(value, dict):
        return value
    out = dict(value)
    for key in ("table", "tables", "frame", "data", "full"):
        frame = out.get(key)
        if isinstance(frame, pd.DataFrame):
            continue
        if isinstance(frame, list) and frame and isinstance(frame[0], dict):
            try:
                out[key] = pd.DataFrame(frame)
            except Exception:  # noqa: BLE001
                continue
    return out


def _apply_result(dataset, payload) -> dict:
    """消费模块返回的副作用键，返回给前端/模型的干净信封。"""
    if not isinstance(payload, dict):
        return jsonable(result(title="完成"))
    outgoing = dict(payload)

    replace = outgoing.pop("replace", None) or outgoing.pop("_replace", None)
    if isinstance(replace, dict) and replace:
        _apply_replace(dataset, replace)

    store = {}
    for key in ("store", "result_store", "_store"):
        value = outgoing.pop(key, None)
        if isinstance(value, dict):
            store.update(value)
    if store:
        results = getattr(dataset, "results", None)
        if not isinstance(results, dict):
            results = {}
        for key, value in store.items():
            try:
                results[str(key)] = _normalize_store(value)
            except Exception:  # noqa: BLE001
                continue
        dataset.results = results

    cluster_key = _apply_cluster_store(dataset, outgoing.pop("cluster_store", None))
    coords_key = _apply_coords_store(dataset, outgoing.pop("coords_store", None))

    warnings: list[str] = []
    columns = outgoing.get("sample_columns")
    values = outgoing.get("sample_values")
    if columns and values:
        written, notes = _apply_sample_columns(dataset, columns, values)
        warnings.extend(notes)
        if written and cluster_key is None and any(
                str(w).lower() in ("clusters", "cluster", "domains") for w in written):
            dataset.cluster_key = written[0]

    summary = outgoing.get("summary")
    if not isinstance(summary, dict):
        summary = {"title": "", "bullets": [], "metrics": {}}
    summary.setdefault("title", "")
    summary.setdefault("bullets", [])
    summary.setdefault("metrics", {})
    outgoing["summary"] = summary

    module_warnings = outgoing.get("warnings")
    if not isinstance(module_warnings, list):
        module_warnings = []
    outgoing["warnings"] = [str(w) for w in module_warnings] + warnings

    for key in ("tables", "figures", "notes", "data_out"):
        if not isinstance(outgoing.get(key), list):
            outgoing[key] = []

    if cluster_key:
        outgoing.setdefault("cluster_key", cluster_key)
    if coords_key:
        outgoing.setdefault("coords_key", coords_key)
    # 写回数据集的 warnings/notes（UI 右侧提示）
    try:
        dataset_notes = getattr(dataset, "notes", None)
        if not isinstance(dataset_notes, list):
            dataset_notes = []
        for note in outgoing["notes"]:
            if str(note) not in dataset_notes:
                dataset_notes.append(str(note))
        dataset.notes = dataset_notes[-60:]
        dataset_warnings = getattr(dataset, "warnings", None)
        if not isinstance(dataset_warnings, list):
            dataset_warnings = []
        for warn in outgoing["warnings"]:
            if str(warn) not in dataset_warnings:
                dataset_warnings.append(str(warn))
        dataset.warnings = dataset_warnings[-60:]
    except Exception:  # noqa: BLE001
        pass
    return jsonable(outgoing)


# ---------------------------------------------------------------------------
# 执行记录（DESIGN.md §2.4 history / analyses_run）
# ---------------------------------------------------------------------------

def _record_step(dataset, step: str, params: dict, ok: bool, ms: int,
                 summary: str = "", error: str | None = None) -> None:
    history = getattr(dataset, "history", None)
    if not isinstance(history, list):
        history = []
    entry = {"at": _now(), "step": str(step), "params": jsonable(params or {}),
             "ok": bool(ok), "ms": int(ms), "summary": str(summary or "")}
    if error:
        entry["error"] = str(error)
    history.append(entry)
    dataset.history = history[-50:]


def _record_analysis(dataset, analysis_id: str, ok: bool, ms: int,
                     summary: str = "", error: str | None = None, params: dict | None = None) -> None:
    runs = getattr(dataset, "analyses_run", None)
    if not isinstance(runs, list):
        runs = []
    entry = {"at": _now(), "id": str(analysis_id), "ok": bool(ok), "ms": int(ms),
             "summary": str(summary or "")}
    if params:
        entry["params"] = jsonable(params)
    if error:
        entry["error"] = str(error)
    runs.append(entry)
    dataset.analyses_run = runs[-50:]


def _mark_step(dataset, step: str) -> list[str]:
    """标记步骤完成；重跑某步会清空其下游步骤（history 保留）。返回被清掉的下游步骤。"""
    order = STEP_INDEX.get(step)
    old = [s for s in (getattr(dataset, "steps", None) or []) if s in STEP_INDEX]
    if order is None:
        steps = list(old)
        if step not in steps:
            steps.append(step)
        dataset.steps = steps
        dataset.current = step
        return []
    steps = [s for s in old if STEP_INDEX[s] < order]
    if step in old:
        # 重跑：本步及下游都要重算
        steps = [s for s in old if STEP_INDEX[s] < order]
    if step not in steps:
        steps.append(step)
    steps.sort(key=lambda s: STEP_INDEX[s])
    dropped = [s for s in old if s not in steps]
    dataset.steps = steps
    dataset.current = step
    if dropped:
        results = getattr(dataset, "results", None)
        if isinstance(results, dict):
            for name in dropped:
                results.pop(name, None)
            dataset.results = results
        notes = getattr(dataset, "notes", None)
        if not isinstance(notes, list):
            notes = []
        notes.append(f"重跑「{step}」已清空下游步骤：{'、'.join(dropped)}")
        dataset.notes = notes[-60:]
    return dropped


# ---------------------------------------------------------------------------
# 数据集构造 / 导入
# ---------------------------------------------------------------------------

DATASET_DEFAULTS = {
    "id": "ds1", "name": "dataset", "path": None, "source_format": "unknown",
    "X": None, "counts": None, "X_kind": "unknown", "genes": [], "spots": [],
    "coords": None, "coords_kind": "generic", "spatial_mode": "generic",
    "neighbor_radius": None, "image": None, "sample_info": None,
    "group_key": None, "groups": None, "species": "unknown", "obsm": {},
    "cluster_key": None, "results": {}, "ref": None, "steps": [], "current": None,
    "history": [], "analyses_run": [], "notes": [], "warnings": [],
}


def ensure_attrs(dataset) -> None:
    """补齐 sp_state.Dataset 上缺失（或为 None 的核心字段），避免 registry/模块到处判空。"""
    for key, default in DATASET_DEFAULTS.items():
        if not hasattr(dataset, key):
            try:
                setattr(dataset, key, default() if callable(default) else
                        (dict(default) if isinstance(default, dict) else
                         (list(default) if isinstance(default, list) else default)))
            except Exception:  # noqa: BLE001
                continue
    X = getattr(dataset, "X", None)
    if X is not None and not isinstance(X, np.ndarray):
        try:
            dataset.X = np.asarray(X, dtype=np.float32)
        except Exception:  # noqa: BLE001
            dataset.X = None
    coords = getattr(dataset, "coords", None)
    if coords is not None and not isinstance(coords, np.ndarray):
        try:
            dataset.coords = np.asarray(coords, dtype=np.float32)
        except Exception:  # noqa: BLE001
            dataset.coords = None
    if not isinstance(getattr(dataset, "obsm", None), dict):
        dataset.obsm = {}
    if not isinstance(getattr(dataset, "results", None), dict):
        dataset.results = {}
    for key in ("steps", "notes", "warnings", "history", "analyses_run", "genes", "spots"):
        if not isinstance(getattr(dataset, key, None), list):
            setattr(dataset, key, [])
    frame = getattr(dataset, "sample_info", None)
    if not isinstance(frame, pd.DataFrame):
        spots = _spots(dataset)
        dataset.sample_info = pd.DataFrame(index=pd.Index(spots, name="spot"))
    if getattr(dataset, "group_key", None) is None:
        guess = _guessed_group_key(dataset)
        if guess:
            dataset.group_key = guess
            dataset.groups = _groups_of(_sample_frame(dataset), guess)
    elif getattr(dataset, "groups", None) is None:
        dataset.groups = _groups_of(_sample_frame(dataset), dataset.group_key)


def _infer_geometry(coords, meta: dict | None, source_format: str) -> dict:
    """估计空间邻接半径与空间模式（DESIGN.md §6 末条）。"""
    meta = meta or {}
    mode = str(meta.get("spatial_mode") or "").strip().lower()
    coords_kind = str(meta.get("coords_kind") or "").strip().lower()
    if not mode:
        fmt = str(source_format or "").lower()
        mode = {"visium": "visium", "visium_h5": "visium", "stereo": "stereo",
                "slideseq": "slideseq", "h5ad": "generic", "matrix_text": "generic"}.get(fmt, "generic")
    if not coords_kind:
        coords_kind = {"visium": "visium_hires", "visium_h5": "visium_hires",
                       "stereo": "grid", "slideseq": "pixel"}.get(mode, meta.get("coords_kind") or "generic")
    radius = meta.get("neighbor_radius")
    if radius in (None, 0, "") and coords is not None and len(coords) >= 3:
        # 用 cKDTree 求最近邻中位数 ×1.2（O(n log n)，200k spot 也不吃内存）
        try:
            from scipy.spatial import cKDTree

            tree = cKDTree(np.asarray(coords, dtype=float))
            dist, _ = tree.query(np.asarray(coords, dtype=float), k=2)
            nn = dist[:, 1] if dist.ndim == 2 else np.asarray(dist)
            nn = nn[np.isfinite(nn) & (nn > 0)]
            if nn.size:
                radius = float(np.median(nn) * 1.2)
        except Exception:  # noqa: BLE001
            radius = None
    try:
        radius = float(radius) if radius is not None else None
    except (TypeError, ValueError):
        radius = None
    return {"spatial_mode": mode or "generic", "coords_kind": coords_kind or "generic",
            "neighbor_radius": radius}


def make_dataset(state, loaded: dict, name: str | None = None, source: str | None = None):
    """sp_io.load_source 的结果 → sp_state.Dataset（§2.4 全部字段）。"""
    loaded = loaded or {}
    meta = loaded.get("meta") or {}
    X = _as_array(loaded.get("X"))
    genes = [str(g) for g in (loaded.get("genes") or [])]
    spots = [str(s) for s in (loaded.get("spots") or [])]
    if X is None:
        raise EngineError("internal", "导入结果缺少表达矩阵 X")
    if X.ndim != 2:
        raise EngineError("internal", f"表达矩阵维度不对（{X.shape}），应为 spot × 基因")
    if len(spots) != X.shape[0] or len(genes) != X.shape[1]:
        # 长度对不上时以矩阵为准，避免后续所有分析错位
        spots = spots if len(spots) == X.shape[0] else [f"spot{i + 1}" for i in range(X.shape[0])]
        genes = genes if len(genes) == X.shape[1] else [f"GENE{i + 1}" for i in range(X.shape[1])]
    coords = _as_array(loaded.get("coords"))
    if coords is not None and (coords.ndim != 2 or coords.shape != (X.shape[0], 2)):
        coords = None
    value_kind = str(loaded.get("X_kind") or loaded.get("value_kind") or meta.get("value_kind")
                     or meta.get("X_kind") or "unknown")
    counts = loaded.get("counts")
    if counts is None:
        counts = meta.get("counts")
    counts = _as_array(counts)
    if counts is None and value_kind == "counts":
        counts = X
    sample_info = loaded.get("sample_info")
    if not isinstance(sample_info, pd.DataFrame):
        sample_info = pd.DataFrame(index=pd.Index(spots, name="spot"))
    elif len(sample_info) != len(spots) or [str(x) for x in sample_info.index] != spots:
        frame = pd.DataFrame(index=pd.Index(spots, name="spot"))
        for col in list(sample_info.columns):
            try:
                frame[col] = sample_info[col].reindex(spots).tolist()
            except Exception:  # noqa: BLE001
                continue
        sample_info = frame
    source_format = str(loaded.get("format") or meta.get("format") or "unknown")
    geometry = _infer_geometry(coords, {**meta, **({"coords_kind": meta.get("coords_kind")} if meta else {})},
                               source_format)
    image = loaded.get("image") or meta.get("image") or None
    if isinstance(image, str):
        image = {"path": image, "extent": None, "scale": meta.get("image_scale")}
    if not isinstance(image, dict):
        image = None
    ds_name = name or loaded.get("name") or meta.get("name") or (
        os.path.basename(str(source).rstrip(os.sep)) if source else "dataset")
    payload = {
        "name": str(ds_name),
        "path": str(source) if source else (str(meta.get("path")) if meta.get("path") else None),
        "source_format": source_format,
        "X": X, "genes": genes, "spots": spots,
        "counts": counts, "X_kind": value_kind,
        "coords": coords, "coords_kind": geometry["coords_kind"],
        "spatial_mode": geometry["spatial_mode"],
        "neighbor_radius": geometry["neighbor_radius"],
        "image": image, "sample_info": sample_info,
        "species": str(meta.get("species") or loaded.get("species") or "unknown"),
        "obsm": {}, "results": {}, "ref": None,
        "steps": [], "current": None, "history": [], "analyses_run": [],
        "notes": [str(n) for n in (loaded.get("notes") or [])],
        "warnings": [str(w) for w in (loaded.get("warnings") or [])],
    }
    dataset = None
    try:
        dataset = state.new(**payload)
    except TypeError:
        # sp_state.State.new 只接受核心参数时，退回逐字段赋值
        core = {key: payload[key] for key in ("name", "X", "genes", "spots") if key in payload}
        try:
            dataset = state.new(**core)
        except TypeError as error:
            raise EngineError("internal", f"sp_state.State.new 签名不兼容：{error}")
        for key, value in payload.items():
            if key in core:
                continue
            if getattr(dataset, key, None) in (None, [], {}) or not hasattr(dataset, key):
                try:
                    setattr(dataset, key, value)
                except Exception:  # noqa: BLE001
                    continue
    except Exception as error:  # noqa: BLE001
        raise EngineError("internal", f"创建数据集失败：{type(error).__name__}: {error}")
    ensure_attrs(dataset)
    # sp_state 可能给了默认 id/name/geometry：只补空值，不覆盖它算好的结果
    for key in ("name", "path", "source_format", "coords_kind", "spatial_mode", "neighbor_radius", "species"):
        if getattr(dataset, key, None) in (None, "", "unknown") and payload.get(key) not in (None, "", "unknown"):
            try:
                setattr(dataset, key, payload[key])
            except Exception:  # noqa: BLE001
                continue
    if getattr(dataset, "neighbor_radius", None) in (None, 0):
        dataset.neighbor_radius = geometry["neighbor_radius"]
    return dataset


def load_dataset(state, workdir: str, path: str, options: dict | None = None, name: str | None = None):
    """路径/URL → Dataset（复用同路径已打开的数据集）。"""
    sp_io = _require("sp_io", "数据导入层（sp_io.py）")
    raw = str(path or "")
    if not raw:
        raise bad_request("请提供 path（本地文件/目录）或 url")
    if not raw.startswith(("http://", "https://", "ftp://")):
        absolute = os.path.abspath(raw)
        try:
            existing = state.find_by_path(absolute)
        except Exception:  # noqa: BLE001
            existing = None
        if existing is not None:
            return existing
        if not os.path.exists(absolute):
            raise EngineError("not_found", f"路径不存在：{absolute}")
        raw = absolute
    opts = {k: v for k, v in (options or {}).items() if v is not None}
    loaded = sp_io.load_source(raw, workdir=workdir, options=opts)
    return make_dataset(state, loaded, name=name, source=raw)


def inspect_source(workdir: str, args: dict) -> dict:
    sp_io = _require("sp_io", "数据导入层（sp_io.py）")
    path = args.get("path") or args.get("url")
    if not path:
        raise bad_request("请提供 path 或 url")
    raw = str(path)
    if not raw.startswith(("http://", "https://", "ftp://")):
        raw = os.path.abspath(raw)
    options = {
        "data_format": args.get("data_format") or args.get("format") or "auto",
        "member": args.get("member"),
        "object": args.get("object"),
        "delimiter": args.get("delimiter"),
        "sheet": args.get("sheet"),
        "gene_col": args.get("gene_col"),
        "coords": args.get("coords"),
        "image": args.get("image"),
        "visium_scale": args.get("visium_scale"),
    }
    options = {k: v for k, v in options.items() if v is not None}
    info = sp_io.inspect_source(raw, workdir=workdir, options=options)
    return jsonable(info)


def formats_catalog() -> dict:
    try:
        sp_io = importlib.import_module("sp_io")
        items = list(getattr(sp_io, "FORMAT_INFO", []) or [])
        ids = list(getattr(sp_io, "SUPPORTED_FORMATS", []) or [])
        if items:
            return {"formats": jsonable(items), "ids": ids}
    except Exception:  # noqa: BLE001 - sp_io 未就绪时给静态清单
        pass
    return {"formats": jsonable(FORMAT_FALLBACK), "ids": [item["id"] for item in FORMAT_FALLBACK]}


def make_demo(state, name: str | None = None):
    sp_state = _require("sp_state", "状态层（sp_state.py）")
    dataset = sp_state.make_demo(state, name=name or "demo-spatial")
    ensure_attrs(dataset)
    return dataset


# ---------------------------------------------------------------------------
# 执行：步骤 / 分析 / 图
# ---------------------------------------------------------------------------

def _step_requires(step_id: str, entry: dict | None) -> list[str]:
    if entry and entry.get("requires"):
        return list(entry["requires"])
    for base in STEP_ORDER:
        if base["id"] == step_id:
            return list(base["requires"])
    return []


def _step_title(step_id: str) -> str:
    entry = all_steps().get(step_id)
    if entry:
        return entry["title"]
    for base in STEP_ORDER:
        if base["id"] == step_id:
            return base["title"]
    return step_id


def run_step(dataset, step_id: str, params: dict, workdir: str, engine_info: dict | None = None) -> dict:
    step_id = str(step_id or "")
    implemented = all_steps()
    entry = implemented.get(step_id)
    if entry is None:
        raise EngineError("unsupported", f"步骤「{step_id}」尚未实现（对应模块未就绪）")
    requires = _step_requires(step_id, entry)
    steps_done = set(getattr(dataset, "steps", None) or [])
    missing = [r for r in requires if r not in steps_done]
    if missing:
        raise prerequisite(f"需要先完成：{'、'.join(_step_title(m) for m in missing)}")
    ensure_attrs(dataset)
    outdir = ensure_dir(os.path.join(workdir, "analyses", str(getattr(dataset, "id", "ds1")), f"step-{step_id}"))
    ctx = build_ctx(dataset, outdir, engine_info)
    if step_id == "deconv":
        _materialize_ref(ctx)
    call_params = dict(params or {})
    if getattr(dataset, "group_key", None) and "group_col" not in call_params:
        call_params["group_col"] = dataset.group_key
    started = time.time()
    try:
        payload = entry["run"](ctx, call_params, outdir)
    except EngineError as error:
        _record_step(dataset, step_id, call_params, False, int((time.time() - started) * 1000),
                     summary=error.message, error=error.message)
        raise
    except Exception as error:  # noqa: BLE001 - 未预期异常转成 internal，保留 trace 到 stderr
        _record_step(dataset, step_id, call_params, False, int((time.time() - started) * 1000),
                     summary=f"{type(error).__name__}: {error}", error=f"{type(error).__name__}: {error}")
        raise
    elapsed = int((time.time() - started) * 1000)
    clean = _apply_result(dataset, payload)
    dropped = _mark_step(dataset, step_id)
    bullets = (clean.get("summary") or {}).get("bullets") or []
    _record_step(dataset, step_id, call_params, True, elapsed,
                 summary=str(bullets[0]) if bullets else "")
    out = dict(clean)
    out.update({
        "step": step_id,
        "step_title": _step_title(step_id),
        "steps": list(getattr(dataset, "steps", None) or []),
        "current": getattr(dataset, "current", None),
        "dropped_steps": dropped,
        "dataset": _summary(dataset),
        "dataset_notes": list(getattr(dataset, "notes", None) or [])[-20:],
        "dataset_warnings": list(getattr(dataset, "warnings", None) or [])[-20:],
        "ms": elapsed,
    })
    return jsonable(out)


def run_analysis(dataset, analysis_id: str, params: dict, workdir: str, engine_info: dict | None = None) -> dict:
    analysis_id = str(analysis_id or "")
    analyses = all_analyses()
    entry = analyses.get(analysis_id)
    if entry is None:
        raise EngineError("not_found", f"分析「{analysis_id}」不存在（模块可能未就绪）")
    ok, reason = needs_ok(dataset, entry["needs"])
    if not ok:
        raise prerequisite(reason or "缺少前置数据")
    ensure_attrs(dataset)
    outdir = ensure_dir(os.path.join(workdir, "analyses", str(getattr(dataset, "id", "ds1")), f"an-{analysis_id}"))
    ctx = build_ctx(dataset, outdir, engine_info)
    if "ref" in entry["needs"]:
        _materialize_ref(ctx)
    call_params = dict(params or {})
    if getattr(dataset, "group_key", None) and "group_col" not in call_params:
        call_params["group_col"] = dataset.group_key
    started = time.time()
    try:
        payload = entry["run"](ctx, call_params, outdir)
    except EngineError as error:
        _record_analysis(dataset, analysis_id, False, int((time.time() - started) * 1000),
                         summary=error.message, error=error.message, params=call_params)
        raise
    except Exception as error:  # noqa: BLE001
        _record_analysis(dataset, analysis_id, False, int((time.time() - started) * 1000),
                         summary=f"{type(error).__name__}: {error}",
                         error=f"{type(error).__name__}: {error}", params=call_params)
        raise
    elapsed = int((time.time() - started) * 1000)
    clean = _apply_result(dataset, payload)
    bullets = (clean.get("summary") or {}).get("bullets") or []
    _record_analysis(dataset, analysis_id, True, elapsed,
                     summary=str(bullets[0]) if bullets else "", params=call_params)
    out = dict(clean)
    out.update({
        "analysis_id": analysis_id,
        "title": entry["title"],
        "group": entry["group"],
        "method": entry["method"],
        "dataset": _summary(dataset),
        "ms": elapsed,
    })
    return jsonable(out)


def run_plot(dataset, kind: str, params: dict, workdir: str, engine_info: dict | None = None) -> dict:
    kind = str(kind or "")
    plots = all_plots()
    entry = plots.get(kind)
    if entry is None:
        implemented = plot_kinds()
        hint = f"已实现：{'、'.join(implemented)}" if implemented else "绘图模块尚未就绪"
        raise EngineError("unsupported", f"图类型「{kind}」尚未实现（{hint}）")
    ensure_attrs(dataset)
    outdir = ensure_dir(os.path.join(workdir, "analyses", str(getattr(dataset, "id", "ds1")), "plots"))
    ctx = build_ctx(dataset, outdir, engine_info)
    payload = entry["fn"](ctx, dict(params or {}), outdir)
    figure = None
    if isinstance(payload, dict):
        if isinstance(payload.get("figures"), list) and payload["figures"]:
            figure = payload["figures"][0]
        elif payload.get("file"):
            figure = payload
    if not isinstance(figure, dict) or not figure.get("file"):
        raise EngineError("internal", f"图「{kind}」没有返回 PNG 路径")
    path = str(figure["file"])
    if not os.path.exists(path):
        raise EngineError("internal", f"图「{kind}」的 PNG 未生成：{path}")
    width = figure.get("width")
    height = figure.get("height")
    if not width or not height:
        try:
            width, height = sp_common.png_size(path)
        except Exception:  # noqa: BLE001
            width = height = None
    out = {
        "file": path,
        "id": str(figure.get("id") or sp_common.slug(os.path.splitext(os.path.basename(path))[0])),
        "title": str(figure.get("title") or kind),
        "width": width,
        "height": height,
        "bytes": os.path.getsize(path),
        "kind": kind,
        "module": entry["module"],
    }
    if isinstance(payload, dict) and isinstance(payload.get("summary"), dict):
        out["summary"] = payload["summary"]
    return jsonable(out)


# ---------------------------------------------------------------------------
# 预览 / 样本表
# ---------------------------------------------------------------------------

def preview(dataset) -> dict:
    """DESIGN.md §2.2 preview。"""
    ensure_attrs(dataset)
    coords = _coords(dataset)
    coords_meta = None
    if coords is not None and len(coords):
        width = float(coords[:, 0].max() - coords[:, 0].min())
        height = float(coords[:, 1].max() - coords[:, 1].min())
        coords_meta = {
            "x_range": [float(coords[:, 0].min()), float(coords[:, 0].max())],
            "y_range": [float(coords[:, 1].min()), float(coords[:, 1].max())],
            "size": [width, height],
        }
    image = getattr(dataset, "image", None)
    image_path = image.get("path") if isinstance(image, dict) else image
    counts = getattr(dataset, "counts", None)
    return jsonable({
        "genes": _genes(dataset)[:50],
        "spots": _spots(dataset)[:10],
        "n_spots": _n_spots(dataset),
        "n_genes": _n_genes(dataset),
        "coords_meta": coords_meta,
        "columns": _sample_columns(dataset),
        "steps": list(getattr(dataset, "steps", None) or []),
        "current": getattr(dataset, "current", None),
        "cluster_keys": _cluster_keys(dataset),
        "spatial_mode": getattr(dataset, "spatial_mode", None),
        "coords_kind": getattr(dataset, "coords_kind", None),
        "neighbor_radius": getattr(dataset, "neighbor_radius", None),
        "image": image_path,
        "X_kind": getattr(dataset, "X_kind", None) or getattr(dataset, "value_kind", None),
        "has_counts": counts is not None,
        "species": getattr(dataset, "species", "unknown"),
        "group_key": getattr(dataset, "group_key", None),
        "groups": getattr(dataset, "groups", None),
        "ref": _ref_payload(dataset),
        "warnings": list(getattr(dataset, "warnings", None) or [])[-10:],
    })


def _sample_rows(dataset, limit: int = 500) -> list[dict]:
    frame = _sample_frame(dataset).copy()
    frame.insert(0, "spot", frame.index.astype(str))
    return jsonable(frame.head(limit).reset_index(drop=True).to_dict(orient="records"))


def _samples_payload(dataset, **extra) -> dict:
    payload = {
        "dataset": _summary(dataset),
        "sample_info": _sample_rows(dataset),
        "groups": getattr(dataset, "groups", None),
        "group_column": getattr(dataset, "group_key", None),
        "total": _n_spots(dataset),
        "columns": _sample_columns(dataset),
    }
    payload.update(extra)
    return jsonable(payload)


def _read_mapping_table(path: str) -> pd.DataFrame:
    lower = str(path).lower()
    if lower.endswith((".xlsx", ".xls")):
        frame = pd.read_excel(path)
    else:
        sep = "\t" if lower.endswith((".tsv", ".txt")) else ","
        try:
            frame = pd.read_csv(path, sep=sep)
        except Exception:  # noqa: BLE001
            frame = pd.read_csv(path, sep=None, engine="python")
    frame.columns = [str(c) for c in frame.columns]
    return frame


def handle_samples(dataset, args: dict) -> dict:
    """DESIGN.md §2.2 samples：spot 级样本/切片表的增删改与分组。"""
    ensure_attrs(dataset)
    op = str(args.get("op") or "list")
    frame = _sample_frame(dataset)
    spots = _spots(dataset)

    if op == "list":
        return _samples_payload(dataset)

    if op == "update":
        entries = args.get("samples") or []
        updated = 0
        for entry in entries if isinstance(entries, list) else []:
            if not isinstance(entry, dict):
                continue
            ident = str(entry.get("id") or entry.get("sample") or entry.get("spot") or "")
            if not ident:
                continue
            mask = None
            if ident in frame.index:
                frame.loc[ident, "display_name"] = str(entry.get("name") or ident)
                if entry.get("group"):
                    frame.loc[ident, "group"] = str(entry["group"])
                updated += 1
                continue
            for col in ("sample", "slice", "sample_id", "group"):
                if col in frame.columns:
                    mask = frame[col].astype(str) == ident
                    if bool(mask.any()):
                        if entry.get("group"):
                            frame.loc[mask, "group"] = str(entry["group"])
                        if entry.get("name"):
                            frame.loc[mask, "sample"] = str(entry["name"])
                        updated += int(mask.sum())
                        break
        dataset.sample_info = frame
        groups = args.get("groups") or []
        if groups:
            notes = [n for n in (getattr(dataset, "notes", None) or []) if not str(n).startswith("分组：")]
            notes.append("分组：" + "/".join(str(g) for g in groups))
            dataset.notes = notes[-60:]
        if "group" in frame.columns:
            dataset.group_key = dataset.group_key or "group"
            dataset.groups = _groups_of(frame, dataset.group_key)
        return _samples_payload(dataset, updated=updated)

    if op == "auto_group":
        pattern = str(args.get("group_regex") or args.get("pattern") or "")
        if not pattern:
            raise bad_request("请提供 group_regex 正则（第 1 个捕获组=大组）")
        try:
            regex = re.compile(pattern)
        except re.error as error:
            raise bad_request(f"正则不合法：{error}")
        column = args.get("column")
        source = frame[column].astype(str).tolist() if column and str(column) in frame.columns else list(spots)
        values = []
        for name in source:
            match = regex.search(str(name))
            values.append(match.group(1) if match and match.groups() else "other")
        frame["group"] = pd.Series(values, index=frame.index, dtype="object")
        dataset.sample_info = frame
        dataset.group_key = "group"
        dataset.groups = _groups_of(frame, "group")
        return _samples_payload(dataset, matched=len(values))

    if op in ("import_map", "join_obs", "join"):
        path = str(args.get("path") or "")
        if not path:
            raise bad_request("请提供 mapping/obs 表路径（path）")
        path = os.path.abspath(path)
        if not os.path.exists(path):
            raise EngineError("not_found", f"文件不存在：{path}")
        table_df = _read_mapping_table(path)
        key = args.get("key") or table_df.columns[0]
        if str(key) not in table_df.columns:
            raise bad_request(f"表里没有列「{key}」（现有列：{list(table_df.columns)}）")
        value_columns = [c for c in table_df.columns if c != str(key)]
        mapping = {}
        for _, row in table_df.iterrows():
            try:
                mapping[str(row[str(key)])] = row[value_columns].to_dict() if len(value_columns) > 1 else row[value_columns[0]]
            except Exception:  # noqa: BLE001
                continue
        matched = 0
        if op == "import_map":
            key_col = args.get("column")
            if key_col and str(key_col) in frame.columns:
                source = frame[str(key_col)].astype(str).tolist()
            elif str(key) in frame.columns:
                source = frame[str(key)].astype(str).tolist()
            else:
                source = list(spots)
            groups = []
            for name in source:
                value = mapping.get(str(name))
                if isinstance(value, dict):
                    value = value.get("group", next(iter(value.values()), "other"))
                if value is None:
                    groups.append("other")
                else:
                    groups.append(str(value))
                    matched += 1
            frame["group"] = pd.Series(groups, index=frame.index, dtype="object")
            dataset.sample_info = frame
            dataset.group_key = "group"
            dataset.groups = _groups_of(frame, "group")
            return _samples_payload(dataset, matched=matched, total=len(groups))
        # join：按 spot 名左右连接任意 obs 列
        lower = {str(k).strip().lower(): k for k in mapping}
        for spot in list(frame.index):
            value = mapping.get(str(spot)) or mapping.get(lower.get(str(spot).strip().lower()))
            if value is None:
                continue
            matched += 1
            if isinstance(value, dict):
                for col, val in value.items():
                    if str(col) not in frame.columns:
                        frame[str(col)] = None
                    frame.loc[spot, str(col)] = val
        dataset.sample_info = frame
        if getattr(dataset, "group_key", None) is None and "group" in frame.columns:
            dataset.group_key = "group"
            dataset.groups = _groups_of(frame, "group")
        return _samples_payload(dataset, matched=matched, total=len(frame))

    if op == "derive":
        name = str(args.get("name") or args.get("samples_op_name") or "derived")
        pattern = str(args.get("group_regex") or args.get("source") or "")
        column = args.get("column")
        regex = None
        if pattern:
            try:
                regex = re.compile(pattern)
            except re.error as error:
                raise bad_request(f"正则不合法：{error}")
        values = []
        for spot in spots:
            if regex is not None:
                match = regex.search(str(spot))
                values.append(match.group(1) if match and match.groups() else "")
            elif column and str(column) in frame.columns:
                values.append(frame.loc[spot, str(column)])
            else:
                values.append("")
        frame[name] = pd.Series(values, index=frame.index, dtype="object")
        dataset.sample_info = frame
        return _samples_payload(dataset, derived=name)

    if op == "set_group":
        column = str(args.get("column") or "")
        if not column:
            raise bad_request("请提供分组列名（column）")
        if column not in frame.columns:
            raise EngineError("not_found", f"列「{column}」不存在（现有列：{list(frame.columns)}）")
        dataset.group_key = column
        dataset.groups = _groups_of(frame, column)
        notes = [n for n in (getattr(dataset, "notes", None) or []) if not str(n).startswith("分组列：")]
        notes.append(f"分组列：{column}")
        dataset.notes = notes[-60:]
        return _samples_payload(dataset)

    if op == "remove":
        names = {str(x) for x in (args.get("samples") or [])}
        if not names:
            raise bad_request("请提供要移除的 spot 名称列表（samples）")
        keep_mask = np.array([str(s) not in names for s in spots], dtype=bool)
        if not keep_mask.any():
            raise bad_request("不能移除全部 spot")
        dataset.X = np.asarray(dataset.X, dtype=np.float32)[keep_mask]
        if getattr(dataset, "counts", None) is not None:
            dataset.counts = np.asarray(dataset.counts, dtype=np.float32)[keep_mask]
        coords = _coords(dataset)
        if coords is not None:
            dataset.coords = coords[keep_mask]
        dataset.spots = [s for s, keep in zip(spots, keep_mask) if keep]
        dataset.sample_info = frame.loc[keep_mask]
        dropped = [s for s in (getattr(dataset, "steps", None) or []) if s not in ("qc", "filter")]
        dataset.steps = [s for s in (getattr(dataset, "steps", None) or []) if s in ("qc", "filter")]
        results = getattr(dataset, "results", None)
        if isinstance(results, dict):
            for name in dropped:
                results.pop(name, None)
            dataset.results = results
        return _samples_payload(dataset, removed=len(names), kept=int(keep_mask.sum()))

    raise bad_request(f"不支持的样本操作：{op}（可选 list/update/auto_group/import_map/join_obs/derive/set_group/remove）")


# ---------------------------------------------------------------------------
# 单细胞参考（h5ad）
# ---------------------------------------------------------------------------

def _h5ad_dims(path: str) -> tuple[int | None, int | None]:
    try:
        import h5py

        with h5py.File(path, "r") as handle:
            X = handle.get("X")
            if X is not None and hasattr(X, "shape") and len(getattr(X, "shape", ())) == 2:
                return int(X.shape[0]), int(X.shape[1])
            if X is not None and hasattr(X, "get"):
                shape = X.get("shape")
                if shape is not None:
                    values = list(shape[()])
                    if len(values) == 2:
                        return int(values[0]), int(values[1])
            obs = handle.get("obs")
            var = handle.get("var")
            if obs is not None and var is not None:
                return int(obs.attrs.get("_index", obs.__len__()) if hasattr(obs, "__len__") else 0), None
    except Exception:  # noqa: BLE001
        pass
    try:
        import anndata

        adata = anndata.read_h5ad(path, backed="r")
        dims = (int(adata.n_obs), int(adata.n_vars))
        try:
            adata.file.close()
        except Exception:  # noqa: BLE001
            pass
        return dims
    except Exception:  # noqa: BLE001
        return None, None


def refs_local(workdir: str) -> dict:
    """探测本机可用的 h5ad 参考（DESIGN §2.2 refs 无 path 分支）。"""
    import glob

    roots = []
    roots.append(os.path.join(workdir, "*.h5ad"))
    roots.append(os.path.join(workdir, "sc-ref", "*.h5ad"))
    roots.append(os.path.join(HERE, "..", "test", "fixtures", "*.h5ad"))
    roots.append(os.path.join(os.getcwd(), "*.h5ad"))
    roots.append(os.path.join(os.path.dirname(workdir.rstrip(os.sep)), "*.h5ad"))
    roots.append(os.path.join(HERE, "..", "..", "dsh-singlecell", "exports", "*.h5ad"))
    found: dict[str, dict] = {}
    for pattern in roots:
        try:
            items = glob.glob(pattern)
        except Exception:  # noqa: BLE001
            continue
        for item in items:
            absolute = os.path.abspath(item)
            if absolute in found:
                continue
            try:
                size = os.path.getsize(absolute)
            except OSError:
                continue
            n_obs, n_vars = (None, None)
            if size < 2_000_000_000:
                n_obs, n_vars = _h5ad_dims(absolute)
            found[absolute] = {"path": absolute, "name": os.path.basename(absolute),
                               "bytes": size, "n_obs": n_obs, "n_vars": n_vars}
    try:
        import anndata  # noqa: F401

        available = True
        reason = None
    except Exception as error:  # noqa: BLE001
        available = False
        reason = f"anndata 不可用：{error}"
    return {
        "refs": sorted(found.values(), key=lambda item: item["path"]),
        "available": available,
        "reason": reason,
        "search_dirs": [os.path.dirname(p) for p in roots],
    }


def ref_info(path: str) -> dict:
    """h5ad 参考元信息（refs 带 path 分支 / ref_set 用）。"""
    absolute = os.path.abspath(str(path or ""))
    if not absolute or not os.path.exists(absolute):
        raise EngineError("not_found", f"参考文件不存在：{absolute}")
    info = {"path": absolute, "name": os.path.basename(absolute),
            "bytes": os.path.getsize(absolute), "n_obs": None, "n_vars": None,
            "columns": [], "cell_types": [], "cell_type_col": None, "error": None}
    try:
        import anndata

        adata = anndata.read_h5ad(absolute, backed="r")
        info["n_obs"] = int(adata.n_obs)
        info["n_vars"] = int(adata.n_vars)
        columns = []
        for col in list(adata.obs.columns):
            try:
                dtype = str(adata.obs[col].dtype)
                kind = "categorical" if dtype in ("category", "object", "bool", "string") else "numeric"
                columns.append({"key": str(col), "kind": kind,
                                "n_unique": int(adata.obs[col].nunique(dropna=True))})
            except Exception:  # noqa: BLE001
                continue
        info["columns"] = columns
        for cand in ("cell_type", "celltype", "cell_type_l1", "cell_type_l2", "annotation",
                     "clusters", "leiden", "seurat_clusters", "louvain"):
            if cand in adata.obs.columns:
                try:
                    info["cell_types"] = [str(v) for v in
                                          list(dict.fromkeys(adata.obs[cand].astype(str).tolist()))][:100]
                    info["cell_type_col"] = cand
                except Exception:  # noqa: BLE001
                    pass
                break
        try:
            adata.file.close()
        except Exception:  # noqa: BLE001
            pass
    except Exception as error:  # noqa: BLE001
        info["error"] = f"{type(error).__name__}: {error}"
        n_obs, n_vars = _h5ad_dims(absolute)
        info["n_obs"], info["n_vars"] = n_obs, n_vars
    return jsonable(info)


def ref_set(dataset, ref_path: str | None = None, cell_type_col: str | None = None) -> dict:
    ensure_attrs(dataset)
    current = _ref_payload(dataset)
    path = str(ref_path or (current or {}).get("path") or "")
    if not path:
        raise bad_request("请提供 ref_path（单细胞 h5ad 绝对路径）")
    info = ref_info(path)
    cell_types = list(info.get("cell_types") or [])
    column = str(cell_type_col or info.get("cell_type_col") or "") or None
    if column and column not in [c["key"] for c in info.get("columns") or []]:
        raise EngineError("not_found", f"参考 h5ad 里没有 obs 列「{column}」")
    dataset.ref = {"path": info["path"], "cell_type_col": column, "cell_types": cell_types}
    notes = [n for n in (getattr(dataset, "notes", None) or []) if not str(n).startswith("单细胞参考：")]
    notes.append(f"单细胞参考：{os.path.basename(info['path'])}（{info.get('n_obs') or '?'} 细胞 × "
                 f"{info.get('n_vars') or '?'} 基因，类型列 {column or '未指定'}）")
    dataset.notes = notes[-60:]
    return {"dataset": _summary(dataset), "ref": jsonable(dataset.ref), "info": info}


# ---------------------------------------------------------------------------
# 自检 / 建议 / 预设流程
# ---------------------------------------------------------------------------

def audit(dataset) -> dict:
    ensure_attrs(dataset)
    steps = set(getattr(dataset, "steps", None) or [])
    checks: list[dict] = []

    def add(cid, level, title, detail, hint=""):
        checks.append({"id": cid, "level": level, "title": title, "detail": detail, "hint": hint})

    coords = _coords(dataset)
    n_spots = _n_spots(dataset)
    # 坐标完整性
    if coords is None:
        add("coords", "fail", "缺少空间坐标", "没有坐标就无法做任何空间分析", "导入 Visium outs/、h5ad（obsm['spatial']）或矩阵+坐标表")
    else:
        finite = np.isfinite(coords).all()
        if finite:
            add("coords", "ok", "空间坐标完整", f"{coords.shape[0]} 个 spot，坐标范围 "
                                f"x[{coords[:, 0].min():.1f},{coords[:, 0].max():.1f}] "
                                f"y[{coords[:, 1].min():.1f},{coords[:, 1].max():.1f}]")
        else:
            add("coords", "fail", "坐标含 NaN/Inf", "有 spot 的坐标缺失", "检查坐标表或改用 Visium 标准目录")
    # 邻接半径 / 空间模式
    radius = getattr(dataset, "neighbor_radius", None)
    if coords is not None:
        if radius is None:
            add("geometry", "warn", "还没有空间邻接半径", "空间统计会用自动估计值", "导入时重新估计或手动设置 neighbor_radius")
        else:
            add("geometry", "ok", "空间邻接半径已估计",
                f"{float(radius):.2f}（模式 {getattr(dataset, 'spatial_mode', 'generic')}）")
    # 单细胞参考
    if getattr(dataset, "ref", None):
        add("ref", "ok", "已挂单细胞参考", f"{_ref_payload(dataset).get('path')}")
    else:
        add("ref", "warn", "没有单细胞参考", "反卷积需要参考 h5ad", "用 ref_attach 挂载 h5ad（可取自 dsh-singlecell）")
    # 流程
    add("qc", "ok" if "qc" in steps else "warn", "空间质控",
        "已完成" if "qc" in steps else "还没做质控", "" if "qc" in steps else "运行「空间质控」")
    kind = str(getattr(dataset, "X_kind", None) or getattr(dataset, "value_kind", "") or "")
    if "normalize" in steps or kind in ("lognorm", "lognormal", "cpm", "tpm"):
        add("normalize", "ok", "已标准化", f"当前值类型：{kind or 'lognorm'}")
    else:
        add("normalize", "fail", "矩阵尚未标准化", f"当前值类型：{kind or 'unknown'}，空间统计会有偏",
            "运行「标准化」")
    has_cluster = _labels(dataset) is not None
    add("cluster", "ok" if has_cluster else "warn", "空间聚类",
        "已完成" if has_cluster else "还没聚类", "" if has_cluster else "运行「聚类」")
    svg_done = "svg" in steps or any(str(k).startswith("svg") for k in (getattr(dataset, "results", {}) or {}))
    add("svg", "ok" if svg_done else "warn", "空间高变基因",
        "已完成" if svg_done else "还没做 SVG（空间自相关）",
        "" if svg_done else "运行「空间高变基因」或分析 morans_svg")
    frame = _sample_frame(dataset)
    domains_done = "domains" in steps or "domains" in frame.columns
    add("domains", "ok" if domains_done else "warn", "空间域",
        "已完成" if domains_done else "还没做空间域", "" if domains_done else "运行「空间域」")
    deconv_done = "deconv" in steps or "deconv" in (getattr(dataset, "results", {}) or {})
    add("deconv", "ok" if deconv_done else "warn", "反卷积",
        "已完成" if deconv_done else ("缺少参考 h5ad" if not getattr(dataset, "ref", None) else "还没做反卷积"),
        "" if deconv_done else "运行「反卷积」或分析 deconv_nnls")
    # 多切片 / 分组
    group_key = getattr(dataset, "group_key", None)
    if group_key:
        groups = getattr(dataset, "groups", None) or _groups_of(frame, group_key)
        add("groups", "ok", "有 spot 级分组", f"分组列 {group_key}：{'/'.join(str(g) for g in (groups or [])[:8])}")
    elif "slice" in frame.columns:
        add("groups", "warn", "有切片列但未设为分组", "多切片/多组比较需要 group_key", "在样本表用 set_group 指定切片列")
    else:
        add("groups", "warn", "没有 spot 级分组", "分组比较类分析（de_groups）不可用", "导入样本表或设置分组列")
    # spot 数
    if n_spots < 200:
        add("spots", "warn", "spot 数偏少", f"当前 {n_spots} 个 spot，统计效能与空间邻域分析受限")
    else:
        add("spots", "ok", "spot 数可用", f"当前 {n_spots} 个 spot")
    # 图像
    image = getattr(dataset, "image", None) or {}
    image_path = image.get("path") if isinstance(image, dict) else image
    if image_path and os.path.exists(str(image_path)):
        add("image", "ok", "有组织图像", str(image_path))
    else:
        add("image", "warn", "没有组织图像", "空间图无法叠加 H&E 底图（不影响分析）", "导入 Visium 的 tissue_hires_image.png")
    fails = sum(1 for c in checks if c["level"] == "fail")
    warns = sum(1 for c in checks if c["level"] == "warn")
    score = max(0, min(100, 100 - 30 * fails - 10 * warns))
    level = "fail" if fails else ("warn" if warns else "ok")
    return {"score": int(score), "level": level, "checks": checks}


def suggest(dataset) -> dict:
    ensure_attrs(dataset)
    steps = set(getattr(dataset, "steps", None) or [])
    frame = _sample_frame(dataset)
    suggestions: list[dict] = []

    def add(sid, title, detail, action="run", step=None, params=None):
        item = {"id": sid, "title": title, "detail": detail, "action": action}
        if step:
            item["step"] = step
        item["params"] = params or {}
        suggestions.append(item)

    if not getattr(dataset, "path", None) and not getattr(dataset, "spots", None):
        add("demo", "载入示例数据", "先载入内置的 2400 spot 合成组织，立刻把流程走一遍", "demo")
        return {"suggestions": suggestions}
    if _coords(dataset) is None:
        add("open", "导入带坐标的数据", "空间工作台需要坐标：Visium outs/ 目录、Stereo-seq BIN 或 h5ad(obsm['spatial'])", "open")
        return {"suggestions": suggestions}
    if "qc" not in steps:
        add("run:qc", "空间质控", "先看每个 spot 的 counts、检出基因数与线粒体比例", step="qc",
            params={"group_col": getattr(dataset, "group_key", None)})
    if "qc" in steps and "filter" not in steps:
        add("run:filter", "过滤低质量 spot/基因", "按 counts/基因数/线粒体比例过滤，再过滤低表达基因", step="filter")
    if "filter" in steps and "normalize" not in steps:
        add("run:normalize", "标准化", "按 spot 总量归一 + log1p，空间统计与聚类都依赖它", step="normalize")
    if "normalize" in steps and "hvg" not in steps:
        add("run:hvg", "高变基因", "挑 2000 个高变基因给 PCA/聚类降噪", step="hvg")
    if "hvg" in steps and "pca" not in steps:
        add("run:pca", "PCA", "对高变基因做主成分分析", step="pca")
    if "pca" in steps and "neighbors" not in steps:
        add("run:neighbors", "邻接图", "表达近邻 + 空间近邻加权，供聚类使用", step="neighbors")
    if "neighbors" in steps and "cluster" not in steps:
        add("run:cluster", "空间聚类", "Louvain 聚类得到 cluster 标签", step="cluster")
    if "pca" in steps and "embed" not in steps:
        add("run:embed", "二维嵌入", "t-SNE/PCA 看聚类结构", step="embed")
    if _labels(dataset) is not None and "svg" not in steps:
        add("run:svg", "空间高变基因", "Moran's I 找空间上成片变化的基因", step="svg", params={"n_top": 30})
    if _labels(dataset) is not None and "domains" not in steps:
        add("run:domains", "空间域", "MRF 平滑把 spot 划分成若干空间域", step="domains", params={"n_domains": 6})
    if "normalize" in steps and "deconv" not in steps:
        if getattr(dataset, "ref", None):
            add("run:deconv", "反卷积", "用参考单细胞数据估计每个 spot 的细胞类型比例", step="deconv")
        else:
            add("refs", "挂载单细胞参考", "反卷积需要参考 h5ad：可上传或从 dsh-singlecell 拉取", "refs")
    if _labels(dataset) is not None:
        add("analysis:morans_svg", "空间自相关（Moran's I）", "对全部高变基因做空间自相关并置换检验", "analysis",
            params={"n_top": 20, "permutations": 100})
        add("analysis:nhood_enrichment", "邻域富集", "看哪些细胞类型/聚类在空间上互为邻居", "analysis")
        add("analysis:de_regions", "两区域差异", "比较两个空间域之间的差异基因", "analysis")
    if "deconv" in steps:
        add("analysis:deconv_roi", "区域细胞组成差异", "比较不同空间域/分组的细胞类型比例", "analysis")
    if steps:
        add("export:report", "导出分析报告", "把 Methods/结果/图整理成 Markdown 报告", "export", params={"what": "report"})
    seen = set()
    unique = []
    for item in suggestions:
        if item["id"] in seen:
            continue
        seen.add(item["id"])
        unique.append(item)
    return {"suggestions": unique[:8]}


RECIPES = [
    {
        "id": "standard",
        "title": "标准空间流程",
        "detail": "质控→过滤→标准化→高变基因→PCA→邻接图→聚类→嵌入→SVG",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {}},
            {"step": "normalize", "params": {"target_sum": 10000, "log": True}},
            {"step": "hvg", "params": {"n_top_genes": 2000}},
            {"step": "pca", "params": {"n_pcs": 30}},
            {"step": "neighbors", "params": {"n_neighbors": 15, "spatial_weight": 0.3}},
            {"step": "cluster", "params": {"resolution": 1.0}},
            {"step": "embed", "params": {"layout": "tsne"}},
            {"step": "svg", "params": {"n_top": 30}},
        ],
    },
    {
        "id": "spatial",
        "title": "空间统计流程",
        "detail": "标准流程 + 空间域（MRF）",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {}},
            {"step": "normalize", "params": {}},
            {"step": "hvg", "params": {}},
            {"step": "pca", "params": {}},
            {"step": "neighbors", "params": {}},
            {"step": "cluster", "params": {}},
            {"step": "svg", "params": {"n_top": 30}},
            {"step": "domains", "params": {"n_domains": 6}},
        ],
    },
    {
        "id": "domain",
        "title": "空间域划分",
        "detail": "标准化→高变基因→PCA→聚类→空间域",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {}},
            {"step": "normalize", "params": {}},
            {"step": "hvg", "params": {}},
            {"step": "pca", "params": {}},
            {"step": "neighbors", "params": {}},
            {"step": "cluster", "params": {}},
            {"step": "domains", "params": {"n_domains": 6, "smooth": 0.5}},
        ],
    },
    {
        "id": "deconv",
        "title": "反卷积流程",
        "detail": "质控→过滤→标准化→反卷积（需先挂参考 h5ad）",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {}},
            {"step": "normalize", "params": {}},
            {"step": "deconv", "params": {"method": "nnls"}},
        ],
    },
    {
        "id": "quick",
        "title": "快速探索",
        "detail": "质控→过滤→标准化→高变基因→PCA→邻接图→聚类",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {}},
            {"step": "normalize", "params": {}},
            {"step": "hvg", "params": {}},
            {"step": "pca", "params": {}},
            {"step": "neighbors", "params": {}},
            {"step": "cluster", "params": {}},
        ],
    },
]


def recipe(recipe_id: str | None = None):
    if not recipe_id:
        return {"recipes": [jsonable(item) for item in RECIPES]}
    for item in RECIPES:
        if item["id"] == recipe_id:
            return jsonable(item)
    raise EngineError("not_found", f"预设流程不存在：{recipe_id}")


# ---------------------------------------------------------------------------
# 导出（DESIGN.md §8 export：10 种）
# ---------------------------------------------------------------------------

def _store_frame(stored) -> pd.DataFrame | None:
    if isinstance(stored, pd.DataFrame):
        return stored
    if isinstance(stored, dict):
        for key in ("table", "full", "frame", "data", "de_table", "svg_table", "values", "df"):
            value = stored.get(key)
            if isinstance(value, pd.DataFrame):
                return value
            if isinstance(value, list) and value and isinstance(value[0], dict):
                try:
                    return pd.DataFrame(value)
                except Exception:  # noqa: BLE001
                    continue
        try:
            return pd.DataFrame([jsonable(stored)])
        except Exception:  # noqa: BLE001
            return None
    if isinstance(stored, list) and stored and isinstance(stored[0], dict):
        try:
            return pd.DataFrame(stored)
        except Exception:  # noqa: BLE001
            return None
    return None


def _find_result(dataset, keys, what: str) -> tuple[str, object]:
    results = getattr(dataset, "results", None)
    results = results if isinstance(results, dict) else {}
    for key in keys:
        if key in results and results[key]:
            return key, results[key]
    for key in results:
        for cand in keys:
            if str(cand) in str(key):
                return str(key), results[key]
    raise prerequisite(f"还没有 {what} 结果，先运行对应步骤/分析")


def _domain_frame(dataset) -> pd.DataFrame:
    frame = _sample_frame(dataset)
    cols = [c for c in frame.columns if str(c) in ("domains", "domain", "domain_anno", "niche")]
    if not cols:
        cols = [c for c in frame.columns if str(c).lower().startswith("domain")]
    if not cols:
        raise prerequisite("还没有空间域结果，先运行「空间域」步骤")
    out = frame[cols].copy()
    out.insert(0, "spot", out.index.astype(str))
    return out.reset_index(drop=True)


def _engine_line(engine_info: dict | None) -> str:
    engine = engine_info or {}
    packages = engine.get("packages") or {}
    pieces = []
    if engine.get("python"):
        pieces.append(f"Python {engine.get('python')}")
    for name in ("numpy", "pandas", "scipy", "sklearn", "matplotlib", "anndata"):
        version = packages.get(name)
        if version:
            pieces.append(f"{name} {version}")
    return "，".join(pieces) or "内置 Python 引擎"


def _fmt_value(value, digits: int = 4) -> str:
    try:
        if isinstance(value, float):
            if value != value:
                return "NA"
            return f"{value:.{digits}g}"
        return str(value)
    except Exception:  # noqa: BLE001
        return str(value)


def _methods_markdown(dataset, engine_info: dict | None = None) -> str:
    lines = ["# 方法（Methods）", "", "## 数据来源", ""]
    lines.append(f"- 数据集：{getattr(dataset, 'name', 'dataset')}")
    if getattr(dataset, "path", None):
        lines.append(f"- 来源路径：`{getattr(dataset, 'path')}`")
    lines.append(f"- 格式：{getattr(dataset, 'source_format', 'unknown')}")
    lines.append(f"- 规模：{_n_spots(dataset)} 个 spot × {_n_genes(dataset)} 个基因")
    lines.append(f"- 物种：{getattr(dataset, 'species', 'unknown')}")
    coords = _coords(dataset)
    if coords is not None:
        lines.append(f"- 空间模式：{getattr(dataset, 'spatial_mode', 'generic')}"
                     f"（坐标类型 {getattr(dataset, 'coords_kind', 'generic')}，"
                     f"邻接半径 {_fmt_value(getattr(dataset, 'neighbor_radius', None), 4)}）")
    lines += ["", "## 计算环境", "", f"- {_engine_line(engine_info)}", "- 所有步骤由 dsh-spatial 内置 Python 实现"
              "（numpy/scipy/sklearn 自行实现，不依赖 scanpy/squidpy）。"]
    lines += ["", "## 步骤与方法", ""]
    history = [h for h in (getattr(dataset, "history", None) or []) if h.get("ok")]
    if not history:
        lines.append("- 还没有成功执行的步骤。")
    for entry in history:
        step = str(entry.get("step") or "")
        lines.append(f"### {_step_title(step)}（{step}）")
        params = entry.get("params") or {}
        if params:
            pairs = "，".join(f"{k}={_fmt_value(v)}" for k, v in params.items() if v not in (None, "", {}))
            lines.append(f"- 参数：{pairs or '默认'}")
        if entry.get("summary"):
            lines.append(f"- 结果：{entry['summary']}")
        lines.append(f"- 耗时：{entry.get('ms', 0)} ms")
        lines.append("")
    analyses = all_analyses()
    runs = [r for r in (getattr(dataset, "analyses_run", None) or []) if r.get("ok")]
    if runs:
        lines += ["## 高级分析", ""]
        for run in runs:
            entry = analyses.get(str(run.get("id")))
            title = entry["title"] if entry else run.get("id")
            method = entry["method"] if entry else ""
            lines.append(f"### {title}（{run.get('id')}）")
            if method:
                lines.append(f"- 方法：{method}")
            if run.get("params"):
                pairs = "，".join(f"{k}={_fmt_value(v)}" for k, v in run["params"].items()
                                  if v not in (None, "", {}, []))
                if pairs:
                    lines.append(f"- 参数：{pairs}")
            if run.get("summary"):
                lines.append(f"- 结果：{run['summary']}")
            lines.append("")
    if getattr(dataset, "ref", None):
        ref = _ref_payload(dataset) or {}
        lines += ["## 反卷积参考", "",
                  f"- 参考 h5ad：`{ref.get('path')}`",
                  f"- 细胞类型列：{ref.get('cell_type_col') or '自动识别'}",
                  f"- 细胞类型数：{len(ref.get('cell_types') or [])}", ""]
    lines += ["## 统计说明", "",
              "- 空间邻接图：坐标最近邻中位数 × 1.2 作为半径，半径图行标准化后用于空间自相关。",
              "- 多重检验：Benjamini-Hochberg FDR。",
              "- 随机性：所有随机过程固定 seed=0（np.random.default_rng(0)），同数据同参数结果可复现。", ""]
    return "\n".join(lines)


def _report_markdown(dataset, engine_info: dict | None = None) -> str:
    lines = ["# 空间转录组分析报告", ""]
    lines += ["## 数据集", "",
              f"- 名称：{getattr(dataset, 'name', 'dataset')}",
              f"- 来源：`{getattr(dataset, 'path', '')}`（{getattr(dataset, 'source_format', 'unknown')}）",
              f"- 规模：{_n_spots(dataset)} spot × {_n_genes(dataset)} 基因；物种 {getattr(dataset, 'species', 'unknown')}",
              ""]
    steps = list(getattr(dataset, "steps", None) or [])
    lines += ["## 已完成步骤", ""]
    if steps:
        for step in steps:
            lines.append(f"- {_step_title(step)}（{step}）")
    else:
        lines.append("- 无")
    lines.append("")
    history = [h for h in (getattr(dataset, "history", None) or []) if h.get("ok")]
    if history:
        lines += ["## 执行记录", ""]
        for entry in history[-20:]:
            summary = entry.get("summary") or ""
            lines.append(f"- {entry.get('at', '')} `{entry.get('step')}` {summary}")
        lines.append("")
    results = getattr(dataset, "results", None) or {}
    if isinstance(results, dict) and results:
        lines += ["## 关键结果", ""]
        for key, value in results.items():
            frame = _store_frame(value)
            if isinstance(frame, pd.DataFrame) and len(frame):
                lines.append(f"### {key}（前 {min(5, len(frame))} 行 / 共 {len(frame)} 行）")
                cols = [str(c) for c in list(frame.columns)[:6]]
                lines.append("| " + " | ".join(cols) + " |")
                lines.append("|" + "|".join(["---"] * len(cols)) + "|")
                for _, row in frame.head(5).iterrows():
                    lines.append("| " + " | ".join(_fmt_value(row[c], 4) for c in cols) + " |")
                lines.append("")
    runs = [r for r in (getattr(dataset, "analyses_run", None) or []) if r.get("ok")]
    if runs:
        analyses = all_analyses()
        lines += ["## 高级分析", ""]
        for run in runs:
            entry = analyses.get(str(run.get("id")))
            lines.append(f"- {entry['title'] if entry else run.get('id')}：{run.get('summary') or '完成'}")
        lines.append("")
    warnings = list(getattr(dataset, "warnings", None) or [])
    if warnings:
        lines += ["## 警告", ""]
        lines += [f"- {w}" for w in warnings[-20:]]
        lines.append("")
    notes = list(getattr(dataset, "notes", None) or [])
    if notes:
        lines += ["## 备注", ""]
        lines += [f"- {n}" for n in notes[-20:]]
        lines.append("")
    lines += ["## 复现", "",
              f"- {_engine_line(engine_info)}",
              "- 用 export what=script 导出可复现脚本，逐步重放本报告中的步骤与参数。", ""]
    return "\n".join(lines)


def _script_python(dataset, engine_info: dict | None = None, filename: str = "replay.py") -> str:
    source = getattr(dataset, "path", None)
    data_format = getattr(dataset, "source_format", "auto")
    history = [h for h in (getattr(dataset, "history", None) or []) if h.get("ok")]
    runs = [r for r in (getattr(dataset, "analyses_run", None) or []) if r.get("ok")]
    steps_repr = ",\n    ".join(
        f"({entry.get('step')!r}, {json.dumps(jsonable(entry.get('params') or {}), ensure_ascii=False)})"
        for entry in history)
    analyses_repr = ",\n    ".join(
        f"({run.get('id')!r}, {json.dumps(jsonable(run.get('params') or {}), ensure_ascii=False)})"
        for run in runs)
    engine_json = json.dumps(jsonable(engine_info or {}), ensure_ascii=False, indent=2)
    return f'''#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""由 dsh-spatial 生成的复现脚本（数据集：{getattr(dataset, "name", "dataset")}）。

用法：
    python3 {os.path.basename(filename)}

脚本用与工作台完全相同的内置实现（sp_io / sp_state / sp_registry）重放已执行过的
步骤与高级分析；不依赖 scanpy/squidpy，也不需要联网。生成时间 {_now()}。
"""

import json
import os
import sys

LIB = {HERE!r}
WORKDIR = {os.environ.get("DS_SP_WORKDIR", os.path.join(os.path.expanduser("~"), "dsh-spatial-work"))!r}
SOURCE = {str(source or "")!r}
DATA_FORMAT = {str(data_format)!r}
STEPS = [
    {steps_repr}
]
ANALYSES = [
    {analyses_repr}
]
ENGINE = {engine_json}


def main():
    if not SOURCE or not os.path.exists(SOURCE):
        raise SystemExit("原始数据路径不存在：%s（请先重新下载/导入）" % SOURCE)
    sys.path.insert(0, LIB)
    import sp_registry  # noqa: E402  与工作台同一套实现
    import sp_state  # noqa: E402

    state = sp_state.State()
    options = {{"data_format": DATA_FORMAT}}
    dataset = sp_registry.load_dataset(state, WORKDIR, SOURCE, options=options,
                                       name={getattr(dataset, "name", "dataset")!r})
    for step, params in STEPS:
        out = sp_registry.run_step(dataset, step, params, WORKDIR, ENGINE)
        print("[step]", step, "ok", (out.get("summary") or {{}}).get("title", ""))
    for analysis_id, params in ANALYSES:
        out = sp_registry.run_analysis(dataset, analysis_id, params, WORKDIR, ENGINE)
        print("[analysis]", analysis_id, "ok", out.get("title", ""))
    summary = sp_registry._summary(dataset)
    print(json.dumps({{"dataset": summary.get("name"), "steps": summary.get("steps")}},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
'''


def export_dataset(dataset, what: str, out_path: str | None, workdir: str,
                   engine_info: dict | None = None) -> dict:
    ensure_attrs(dataset)
    what = str(what or "summary_json")
    if what not in EXPORT_KINDS:
        raise bad_request(f"不支持的导出类型：{what}（可选：{'、'.join(EXPORT_KINDS)}）")
    export_dir = ensure_dir(os.path.join(workdir, "exports"))
    stamp = time.strftime("%Y%m%d-%H%M%S")
    base = f"{getattr(dataset, 'id', 'ds1')}-{what}-{stamp}"

    def target(ext: str) -> str:
        if out_path:
            directory = os.path.dirname(os.path.abspath(out_path))
            if directory:
                ensure_dir(directory)
            return os.path.abspath(out_path)
        return os.path.join(export_dir, f"{base}.{ext}")

    def done(path: str, **extra) -> dict:
        payload = {"file": path, "bytes": os.path.getsize(path), "what": what}
        payload.update(extra)
        return jsonable(payload)

    if what == "expr_csv":
        path = target("csv")
        frame = pd.DataFrame(np.asarray(dataset.X, dtype=np.float32),
                             index=pd.Index(_spots(dataset), name="spot"), columns=_genes(dataset))
        frame.to_csv(path, encoding="utf-8")
        return done(path)

    if what == "de_csv":
        key, stored = _find_result(dataset, ["de", "de_table", "de_domain", "de_regions", "de_groups"],
                                   "差异分析（先跑 de_domain / de_regions / de_groups 之一）")
        frame = _store_frame(stored)
        if frame is None or not len(frame):
            raise prerequisite(f"还没有可导出的差异表：请先跑「区域差异」里的分析（缓存键 {key or '空'}）")
        path = target("csv")
        frame.to_csv(path, index=False, encoding="utf-8")
        return done(path, rows=int(len(frame)))

    if what == "sample_info_csv":
        path = target("csv")
        frame = _sample_frame(dataset).copy()
        frame.insert(0, "spot", frame.index.astype(str))
        frame.to_csv(path, index=False, encoding="utf-8")
        return done(path)

    if what == "svg_csv":
        key, stored = _find_result(dataset, ["svg", "svg_table", "moran", "morans"], "空间高变基因（SVG）")
        frame = _store_frame(stored)
        if frame is None or not len(frame):
            raise prerequisite("SVG 结果为空，先运行「空间高变基因」或分析 morans_svg")
        path = target("csv")
        frame.to_csv(path, index=False, encoding="utf-8")
        return done(path, source=key, rows=int(len(frame)))

    if what == "domain_csv":
        frame = _domain_frame(dataset)
        path = target("csv")
        frame.to_csv(path, index=False, encoding="utf-8")
        return done(path, rows=int(len(frame)))

    if what == "deconv_csv":
        key, stored = _find_result(dataset, ["deconv", "deconv_nnls", "deconv_dwls", "proportions"], "反卷积")
        frame = _store_frame(stored)
        if frame is None or not len(frame):
            raise prerequisite("反卷积结果为空，先运行「反卷积」或分析 deconv_nnls")
        path = target("csv")
        index = bool(frame.index.name or str(frame.index.dtype).startswith("object"))
        frame.to_csv(path, index=index, encoding="utf-8")
        return done(path, source=key, rows=int(len(frame)))

    if what == "summary_json":
        path = target("json")
        payload = {
            "plugin": PLUGIN,
            "dataset": _summary(dataset),
            "steps": list(getattr(dataset, "steps", None) or []),
            "history": list(getattr(dataset, "history", None) or []),
            "analyses_run": list(getattr(dataset, "analyses_run", None) or []),
            "engine": jsonable(engine_info or {}),
            "results_available": sorted((getattr(dataset, "results", None) or {}).keys())
            if isinstance(getattr(dataset, "results", None), dict) else [],
            "ref": _ref_payload(dataset),
            "exported_at": _now(),
        }
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(jsonable(payload), handle, ensure_ascii=False, indent=2)
        return done(path)

    if what == "protocol":
        path = target("json")
        history = list(getattr(dataset, "history", None) or [])
        payload = {
            "plugin": PLUGIN,
            "dataset": _summary(dataset),
            "source": {"path": getattr(dataset, "path", None),
                       "format": getattr(dataset, "source_format", "unknown")},
            "geometry": {
                "spatial_mode": getattr(dataset, "spatial_mode", None),
                "coords_kind": getattr(dataset, "coords_kind", None),
                "neighbor_radius": getattr(dataset, "neighbor_radius", None),
            },
            "engine": jsonable(engine_info or {}),
            "steps": [{"step": entry.get("step"), "params": entry.get("params"),
                       "ok": entry.get("ok"), "ms": entry.get("ms")} for entry in history],
            "analyses": [{"id": run.get("id"), "params": run.get("params"), "ok": run.get("ok"),
                          "ms": run.get("ms")} for run in (getattr(dataset, "analyses_run", None) or [])],
            "ref": _ref_payload(dataset),
            "seed": 0,
            "exported_at": _now(),
        }
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(jsonable(payload), handle, ensure_ascii=False, indent=2)
        return done(path)

    if what == "methods":
        path = target("md")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(_methods_markdown(dataset, engine_info))
        return done(path)

    if what == "report":
        path = target("md")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(_report_markdown(dataset, engine_info))
        return done(path)

    if what == "script":
        path = target("py")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(_script_python(dataset, engine_info, filename=os.path.basename(path)))
        return done(path)

    raise bad_request(f"不支持的导出类型：{what}")


# 兼容 engine 里对下划线函数的调用（测试/脚本复用）
summary = _summary
