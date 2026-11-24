"""
dsh-spatial — 数据集状态层（DESIGN §2.3 / §2.4 / §9）

职责：
  * Dataset：一个空间数据集的全部内存状态
      - 表达矩阵 X / counts（(n_spots, n_genes) float32）
      - 空间坐标 coords（(n_spots, 2)，图像像素坐标，Visium 已按 scalefactors 缩放）
      - spot 级样本信息 sample_info（index=spots，可含 slice/group 等列）
      - 步骤状态 steps/current、history、analyses_run、results、obsm、ref
  * State：多数据集注册表（new/get/close/reset/find_by_path/summaries）
  * estimate_geometry()：由坐标估计空间几何（mode/coords_kind/x_range/y_range/neighbor_radius）
  * make_demo()：DESIGN §9 的内置示例（2400 spots × 1000 基因 × 5 空间域 × 2 切片 × 2 分组）

约定：
  * X / counts 一律 float32、形状 (n_spots, n_genes)；
  * coords 一律 float32、形状 (n_spots, 2)，x = 图像列方向、y = 图像行方向（y 向下）；
  * 本文件只依赖标准库 + numpy/pandas/scipy/sklearn + sp_common（禁止 import 其它 sp_* 模块）。
"""

from __future__ import annotations

import math
import os
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from sp_common import EngineError, jsonable

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

#: 步骤顺序（重跑某步时清空其下游；DESIGN §3）
STEP_ORDER = [
    "qc", "filter", "normalize", "hvg", "pca", "neighbors", "cluster", "embed",
    "svg", "domains", "deconv", "annotate",
]

#: 坐标语义（DESIGN §2.3 coords_kind）
COORDS_KINDS = ("visium_hires", "visium_lowres", "pixel", "grid", "generic")

#: 空间模式（DESIGN §2.3 spatial_mode）
SPATIAL_MODES = ("visium", "grid", "generic", "stereo", "slideseq")

#: 默认聚类列
DEFAULT_CLUSTER_KEY = "clusters"

#: history / analyses_run 上限（DESIGN §2.4：history 最多 50 条，返回 UI 截断到 20）
MAX_HISTORY = 50
MAX_ANALYSES_RUN = 100
SUMMARY_HISTORY_LIMIT = 20

#: demo 规模（DESIGN §9）
DEMO_N_SPOTS = 2400
DEMO_N_GENES = 1000
DEMO_N_DOMAINS = 5
DEMO_MARKERS_PER_DOMAIN = 40
DEMO_SEED = 0

#: demo 里被识别为"分组列"的列名优先级
_GROUP_COLUMN_PREFERENCE = (
    "group", "Group", "GROUP", "condition", "Condition", "tissue", "Tissue",
    "slice", "Slice", "disease", "Disease", "type", "Type", "category",
    "phenotype", "region", "sample", "patient",
)

#: 这些列是"聚类结果"而不是实验分组，不能当 group_key
_CLUSTER_LIKE_COLUMNS = {"clusters", "domains", "niches", "niche", "cluster", "domain", "leiden", "louvain"}


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------

def now_iso() -> str:
    """UTC ISO8601（毫秒）时间戳，供 history/analyses_run 使用。"""
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def _dedupe(names, prefix: str = "spot") -> list[str]:
    """去空/去引号/同名加后缀（保持顺序），保证 spots/genes 唯一。"""
    seen: dict[str, int] = {}
    out: list[str] = []
    for index, raw in enumerate(names or []):
        name = str(raw).strip().strip('"').strip("'")
        if not name:
            name = f"{prefix}{index + 1}"
        if name in seen:
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 0
        out.append(name)
    return out


def column_kind(series: pd.Series) -> str:
    """判断 spot 级样本信息列的种类：categorical / numeric（足够 UI 区分即可）。"""
    try:
        if series.dropna().empty:
            return "categorical"
        if pd.api.types.is_numeric_dtype(series):
            if series.nunique(dropna=True) > 12:
                return "numeric"
            return "categorical"
        return "categorical"
    except Exception:  # noqa: BLE001
        return "categorical"


def _float(value):
    """None / NaN / Inf 安全转 float。"""
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(out):
        return None
    return out


def _median_nn_distance(coords) -> float | None:
    """中位最近邻距离（cKDTree，O(n log n)，大样本也不爆内存）。"""
    arr = np.asarray(coords, dtype=float)
    n = int(arr.shape[0])
    if n < 2:
        return None
    if n == 2:
        dist = float(np.linalg.norm(arr[0] - arr[1]))
        return dist if dist > 0 else None
    try:
        from scipy.spatial import cKDTree

        tree = cKDTree(arr)
        dists, _ = tree.query(arr, k=2)
        nn = np.asarray(dists, dtype=float)[:, 1]
    except Exception:  # noqa: BLE001 - scipy 不可用时退回 numpy（小数据）
        if n > 4000:
            return None
        d = np.sqrt(((arr[:, None, :] - arr[None, :, :]) ** 2).sum(axis=2))
        np.fill_diagonal(d, np.inf)
        nn = np.min(d, axis=1)
    nn = nn[np.isfinite(nn) & (nn > 0)]
    if nn.size == 0:
        return None
    return float(np.median(nn))


def _mode_from_kind(kind: str) -> str | None:
    """coords_kind → spatial_mode 的默认映射（pixel 可能是 Stereo/Slide-seq，按 generic 处理）。"""
    key = str(kind or "").strip().lower()
    if key in ("visium_hires", "visium_lowres", "visium"):
        return "visium"
    if key in ("grid", "array", "hex"):
        return "grid"
    if key in ("pixel", "generic", "image", "px"):
        return "generic"
    if key in ("stereo", "bin"):
        return "stereo"
    if key in ("slideseq", "bead"):
        return "slideseq"
    return None


def _looks_regular_lattice(coords) -> bool:
    """用 1-7 近邻距离比例判断坐标是否落在规则格点（六边形/方形）。"""
    arr = np.asarray(coords, dtype=float)
    n = int(arr.shape[0])
    if n < 12:
        return False
    try:
        from scipy.spatial import cKDTree

        tree = cKDTree(arr)
        dists, _ = tree.query(arr, k=min(8, n))
        dists = np.asarray(dists, dtype=float)[:, 1:]  # 去掉自身
    except Exception:  # noqa: BLE001
        return False
    d1 = dists[:, 0]
    ok = (d1 > 0) & np.isfinite(d1)
    if ok.sum() < 8:
        return False
    d1 = d1[ok]
    ratio = dists[ok] / np.where(d1 > 0, d1, 1.0)[:, None]
    d2 = float(np.median(ratio[:, 1]))
    d4 = float(np.median(ratio[:, 3]))
    # 六边形：第 2~6 近邻几乎等距（≤1.08×d1）；方形：第 2~4 近邻 ≈√2×d1
    hex_like = d2 <= 1.08
    square_like = abs(d2 - math.sqrt(2.0)) <= 0.15 and abs(d4 - math.sqrt(2.0)) <= 0.15
    return bool(hex_like or square_like)


def _coord_ranges(arr: np.ndarray) -> tuple[list[float] | None, list[float] | None]:
    if arr is None or arr.size == 0:
        return None, None
    x = arr[:, 0]
    y = arr[:, 1]
    if not np.isfinite(x).any() or not np.isfinite(y).any():
        return None, None
    return [float(np.nanmin(x)), float(np.nanmax(x))], [float(np.nanmin(y)), float(np.nanmax(y))]


def estimate_geometry(coords, coords_kind=None) -> dict:
    """由坐标估计空间几何（DESIGN §6 / §2.4 geometry）。

    返回 {"mode","coords_kind","x_range","y_range","neighbor_radius"}：
      * mode：visium|grid|generic|stereo|slideseq（coords_kind 优先，其次按格点规律猜测）；
      * neighbor_radius：Visium 取 spot 中心距（≈中位最近邻距离，留 2% 容差），
        grid/generic 取 1.2×中位最近邻距离（邻居图能覆盖 4~6 个相邻 spot）。
    方法依据：10x Visium 阵列间距 100 µm（Visium 用户指南）；半径邻接常用 1.2×NN 经验值。
    """
    kind = str(coords_kind or "").strip().lower()
    if coords is None:
        return {"mode": _mode_from_kind(kind) or "generic", "coords_kind": kind or None,
                "x_range": None, "y_range": None, "neighbor_radius": None}
    arr = np.asarray(coords, dtype=float)
    if arr.ndim != 2 or arr.shape[0] == 0:
        raise EngineError("bad_request", "坐标必须是二维数组（n_spots × 2）")
    if arr.shape[1] < 2:
        raise EngineError("bad_request", f"坐标需要 x/y 两列，当前只有 {arr.shape[1]} 列")
    arr = arr[:, :2]
    x_range, y_range = _coord_ranges(arr)
    spacing = _median_nn_distance(arr)
    mode = _mode_from_kind(kind)
    if mode is None:
        mode = "grid" if _looks_regular_lattice(arr) else "generic"
    if spacing is None:
        radius = None
    elif mode == "visium":
        radius = spacing * 1.02
    else:
        radius = spacing * 1.2
    if radius is not None:
        radius = float(max(radius, 1e-6))
    return {
        "mode": mode,
        "coords_kind": kind or None,
        "x_range": x_range,
        "y_range": y_range,
        "neighbor_radius": radius,
        "spacing": float(spacing) if spacing is not None else None,
    }


def build_sample_info(spots, info=None, extra=None) -> pd.DataFrame:
    """把导入层给的 spot 级信息整理成规整 DataFrame（index=spots，按 spots 顺序对齐）。"""
    index = pd.Index([str(s) for s in spots], name="spot")
    frame = pd.DataFrame(index=index)
    if isinstance(info, pd.DataFrame) and len(info):
        src = info.copy()
        src.columns = [str(c) for c in src.columns]
        src.index = [str(i) for i in src.index]
        try:
            aligned = src.reindex(index)
            if aligned.isna().all(axis=None) and len(src) == len(index):
                aligned = src.copy()
                aligned.index = index
            frame = frame.join(aligned)
        except Exception:  # noqa: BLE001
            pass
    if extra:
        for key, values in extra.items():
            if key in frame.columns:
                continue
            if isinstance(values, (list, tuple, np.ndarray)) and len(values) == len(index):
                frame[str(key)] = list(values)
    frame.index = index
    return frame


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

class Dataset:
    """一个空间转录组数据集：表达矩阵 + 坐标 + spot 信息 + 流程状态。

    属性（DESIGN §2.3）：
      id, name, path, source_format, X, counts, X_kind, genes, spots, coords, coords_kind,
      spatial_mode, neighbor_radius, image, sample_info, group_key, groups, species, obsm,
      cluster_key, results, ref, steps, current, history, analyses_run, notes, warnings
    """

    def __init__(
        self,
        ds_id: str | None = None,
        *,
        name: str = "spatial",
        X,
        genes=None,
        spots=None,
        counts=None,
        coords=None,
        coords_kind: str | None = None,
        spatial_mode: str | None = None,
        neighbor_radius: float | None = None,
        image=None,
        sample_info: pd.DataFrame | None = None,
        group_key: str | None = None,
        species: str = "unknown",
        X_kind: str | None = None,
        path: str | None = None,
        source_format: str = "unknown",
        obsm=None,
        ref=None,
        notes=None,
        warnings=None,
        cluster_key: str = DEFAULT_CLUSTER_KEY,
        **extra,
    ):
        arr = np.asarray(X, dtype=np.float32)
        if arr.ndim != 2:
            raise EngineError("bad_request", "表达矩阵必须是二维（spot × 基因）")
        n_spots, n_genes = int(arr.shape[0]), int(arr.shape[1])
        if n_spots == 0 or n_genes == 0:
            raise EngineError("bad_request", "表达矩阵为空（0 个 spot 或 0 个基因）")
        self.id = str(ds_id or "")
        self.name = str(name or "spatial")
        self.path = None if path is None else str(path)
        self.source_format = str(source_format or "unknown")
        self.X = arr
        self.counts = None if counts is None else np.asarray(counts, dtype=np.float32)
        if self.counts is not None and self.counts.shape != arr.shape:
            raise EngineError(
                "bad_request",
                f"counts 形状 {self.counts.shape} 与矩阵形状 {arr.shape} 不一致",
            )
        self.genes = _dedupe(genes if genes is not None else [f"GENE{i + 1}" for i in range(n_genes)], "gene")
        if len(self.genes) != n_genes:
            raise EngineError("bad_request", f"基因数（{len(self.genes)}）与矩阵列数（{n_genes}）不一致")
        self.spots = _dedupe(spots if spots is not None else [f"spot{i + 1}" for i in range(n_spots)], "spot")
        if len(self.spots) != n_spots:
            raise EngineError("bad_request", f"spot 数（{len(self.spots)}）与矩阵行数（{n_spots}）不一致")
        self.X_kind = self._resolve_value_kind(X_kind)

        # 坐标
        if coords is None:
            self.coords = None
        else:
            carr = np.asarray(coords, dtype=np.float32)
            if carr.ndim != 2 or carr.shape[1] < 2:
                raise EngineError("bad_request", "坐标必须是二维数组（n_spots × 2）")
            if carr.shape[0] != n_spots:
                raise EngineError("bad_request", f"坐标行数（{carr.shape[0]}）与 spot 数（{n_spots}）不一致")
            self.coords = carr[:, :2]
        self.coords_kind = str(coords_kind or ("generic" if self.coords is not None else "")).strip() or None

        geometry = estimate_geometry(self.coords, self.coords_kind)
        mode = str(spatial_mode or geometry.get("mode") or "generic").strip().lower()
        if mode not in SPATIAL_MODES:
            mode = "generic"
        self.spatial_mode = mode
        geometry["mode"] = mode
        radius = _float(neighbor_radius)
        if radius is None:
            radius = _float(geometry.get("neighbor_radius"))
        self.neighbor_radius = radius
        geometry["neighbor_radius"] = radius
        self.geometry = geometry

        self.image = self._normalize_image(image)
        self.species = str(species or "unknown").strip().lower() or "unknown"
        self.obsm: dict = dict(obsm or {})
        if self.coords is not None:
            self.obsm.setdefault("spatial", self.coords)

        self.sample_info = build_sample_info(self.spots, sample_info)
        self.notes = [str(item) for item in (notes or [])]
        self.warnings = [str(item) for item in (warnings or [])]
        self.results: dict = {}
        self.ref = dict(ref) if isinstance(ref, dict) else None
        self.cluster_key = str(cluster_key or DEFAULT_CLUSTER_KEY)
        self.cluster_keys: list[str] = []
        self.steps: list[str] = []
        self.current: str | None = None
        self.history: list[dict] = []
        self.analyses_run: list[dict] = []
        self.created_at = now_iso()
        self.group_key: str | None = None
        self.group_levels: list[str] = []
        if group_key:
            self.set_group(group_key)
        else:
            guess = self._guess_group()
            if guess:
                self.group_key = guess
        if extra:
            # 允许上层传额外元数据；不认识的键只在 notes 里记一句，避免静默丢参数
            known = ", ".join(sorted(extra))[:200]
            self.notes.append(f"忽略未识别的数据集参数：{known}")

    # -- 构造期辅助 -------------------------------------------------------
    def _resolve_value_kind(self, X_kind: str | None) -> str:
        kind = str(X_kind or "").strip().lower()
        if kind in ("lognorm", "log", "lognormalized", "normalized"):
            return "lognorm"
        if kind in ("counts", "count", "raw"):
            return "counts"
        if kind in ("unknown", "", "none"):
            pass
        if kind:
            return kind
        source = self.counts if self.counts is not None else self.X
        sample = np.asarray(source).ravel()
        if sample.size == 0:
            return "unknown"
        if sample.size > 200_000:
            sample = sample[:: max(1, sample.size // 200_000)]
        finite = sample[np.isfinite(sample)]
        if finite.size and np.all(finite >= 0) and np.allclose(finite, np.round(finite), atol=1e-4):
            return "counts"
        return "unknown"

    @staticmethod
    def _normalize_image(image) -> dict:
        """image 支持 None / 路径字符串 / {"path","extent","scale"}。"""
        if not image:
            return {"path": None, "extent": None, "scale": None}
        if isinstance(image, str):
            return {"path": str(image), "extent": None, "scale": None}
        if isinstance(image, dict):
            extent = image.get("extent")
            if extent is not None:
                try:
                    extent = [float(v) for v in extent]
                except (TypeError, ValueError):
                    extent = None
                if extent is not None and len(extent) != 4:
                    extent = None
            return {
                "path": image.get("path"),
                "extent": extent,
                "scale": _float(image.get("scale")),
            }
        return {"path": None, "extent": None, "scale": None}

    # -- 基本属性 ---------------------------------------------------------
    @property
    def n_spots(self) -> int:
        return int(self.X.shape[0])

    @property
    def n_genes(self) -> int:
        return int(self.X.shape[1])

    @property
    def has_counts(self) -> bool:
        return self.counts is not None

    @property
    def groups(self) -> list[str] | None:
        """当前分组列的取值列表（无分组→None）。"""
        if not self.group_key or self.group_key not in self.sample_info.columns:
            return None
        values = self.sample_info[self.group_key].astype(str)
        out = [value for value in dict.fromkeys(values.tolist()) if value not in ("nan", "None", "")]
        return out or None

    def group_labels(self) -> np.ndarray:
        """每个 spot 的分组标签（无分组时全为 "all"）。"""
        if not self.group_key or self.group_key not in self.sample_info.columns:
            return np.array(["all"] * self.n_spots)
        return self.sample_info[self.group_key].astype(str).values

    def _guess_group(self) -> str | None:
        candidates = [col for col in self.sample_info.columns
                      if str(col).lower() not in _CLUSTER_LIKE_COLUMNS]
        for prefer in _GROUP_COLUMN_PREFERENCE:
            if prefer in candidates:
                try:
                    n_unique = int(self.sample_info[prefer].nunique(dropna=True))
                except Exception:  # noqa: BLE001
                    continue
                if 1 < n_unique <= 20:
                    return prefer
        for col in candidates:
            try:
                series = self.sample_info[col]
                n_unique = int(series.nunique(dropna=True))
            except Exception:  # noqa: BLE001
                continue
            if 1 < n_unique <= 10 and not pd.api.types.is_numeric_dtype(series):
                return str(col)
        return None

    def set_group(self, key: str) -> None:
        """把 sample_info 的某列设为分组列（spot 级）。"""
        key = str(key or "").strip()
        if not key or key not in self.sample_info.columns:
            raise EngineError("bad_request", f"样本信息里没有列「{key}」；可用列：{list(self.sample_info.columns)[:12]}")
        self.group_key = key
        values = self.sample_info[key].astype(str)
        self.group_levels = [v for v in dict.fromkeys(values.tolist()) if v not in ("nan", "None", "")]

    def sample_columns(self) -> list[dict]:
        """spot 级信息列清单（DESIGN §2.4 sample_columns）。"""
        out = []
        for col in self.sample_info.columns:
            series = self.sample_info[col]
            try:
                n_unique = int(series.nunique(dropna=True))
            except Exception:  # noqa: BLE001
                n_unique = 0
            out.append({"key": str(col), "kind": column_kind(series), "n_unique": n_unique})
        return out

    def labels(self):
        """当前聚类标签 np.ndarray（字符串）或 None（还没聚类）。"""
        key = self.cluster_key or DEFAULT_CLUSTER_KEY
        if key not in self.sample_info.columns:
            return None
        try:
            values = self.sample_info[key].astype(str).values
        except Exception:  # noqa: BLE001
            return None
        if values.size != self.n_spots:
            return None
        return values

    def all_cluster_keys(self) -> list[str]:
        """所有可用的聚类/域列名（DESIGN §2.4 cluster_keys）。"""
        out = list(self.cluster_keys)
        for col in self.sample_info.columns:
            key = str(col)
            if key in _CLUSTER_LIKE_COLUMNS and key not in out:
                out.append(key)
        return out

    # -- 摘要 / 预览 ------------------------------------------------------
    def summary(self, history_limit: int = SUMMARY_HISTORY_LIMIT) -> dict:
        """DESIGN §2.4 dataset 摘要（键名逐字一致；history 返回截断到 20 条）。"""
        image = self.image or {}
        out = {
            "id": self.id,
            "name": self.name,
            "path": self.path,
            "source_format": self.source_format,
            "n_spots": self.n_spots,
            "n_genes": self.n_genes,
            "species": self.species,
            "geometry": {
                "mode": self.spatial_mode,
                "coords_kind": self.coords_kind,
                "x_range": self.geometry.get("x_range"),
                "y_range": self.geometry.get("y_range"),
                "image": image.get("path"),
                "scale": image.get("scale"),
                "neighbor_radius": self.neighbor_radius,
            },
            "sample_columns": self.sample_columns(),
            "steps": list(self.steps),
            "current": self.current,
            "cluster_keys": self.all_cluster_keys(),
            "ref": self._ref_summary(),
            "warnings": list(self.warnings),
            "notes": list(self.notes),
            "history": list(self.history[-max(0, int(history_limit)):]),
            "analyses_run": list(self.analyses_run[-max(0, int(history_limit)):]),
            "group_key": self.group_key,
            "groups": self.groups,
            "has_counts": self.has_counts,
            "X_kind": self.X_kind,
            "created_at": self.created_at,
        }
        return jsonable(out)

    def _ref_summary(self) -> dict | None:
        if not isinstance(self.ref, dict) or not self.ref:
            return None
        return jsonable({
            "path": self.ref.get("path"),
            "cell_type_col": self.ref.get("cell_type_col"),
            "cell_types": list(self.ref.get("cell_types") or []),
        })

    def preview(self) -> dict:
        """DESIGN §2.2 preview（键名逐字一致；不返回整张矩阵）。"""
        x_range = self.geometry.get("x_range")
        y_range = self.geometry.get("y_range")
        out = {
            "genes": list(self.genes[:50]),
            "spots": list(self.spots[:10]),
            "n_spots": self.n_spots,
            "n_genes": self.n_genes,
            "coords_meta": {
                "x_range": x_range,
                "y_range": y_range,
                "size": [self.n_spots, 2] if self.coords is not None else [self.n_spots, 0],
            },
            "columns": self.sample_columns(),
            "steps": list(self.steps),
            "current": self.current,
            "cluster_keys": self.all_cluster_keys(),
            "spatial_mode": self.spatial_mode,
            "image": (self.image or {}).get("path"),
        }
        return jsonable(out)

    @property
    def coords_meta(self) -> dict:
        """坐标元信息（preview/geom 之外的地方也常用）。"""
        return {
            "x_range": self.geometry.get("x_range"),
            "y_range": self.geometry.get("y_range"),
            "coords_kind": self.coords_kind,
            "neighbor_radius": self.neighbor_radius,
        }

    # -- 历史 -------------------------------------------------------------
    def add_history(self, step: str, params: dict, ok: bool, ms: int, summary: str = "") -> None:
        """记录一次步骤执行（DESIGN §2.4 history 条目）。"""
        self.history.append({
            "at": now_iso(),
            "step": str(step),
            "params": jsonable(params or {}),
            "ok": bool(ok),
            "ms": int(ms),
            "summary": str(summary or ""),
        })
        if len(self.history) > MAX_HISTORY:
            self.history = self.history[-MAX_HISTORY:]

    def add_analysis_run(self, analysis_id: str, ok: bool, ms: int, summary: str = "") -> None:
        """记录一次分析执行（DESIGN §2.4 analyses_run 条目）。"""
        self.analyses_run.append({
            "at": now_iso(),
            "id": str(analysis_id),
            "ok": bool(ok),
            "ms": int(ms),
            "summary": str(summary or ""),
        })
        if len(self.analyses_run) > MAX_ANALYSES_RUN:
            self.analyses_run = self.analyses_run[-MAX_ANALYSES_RUN:]

    def mark_step(self, step: str) -> None:
        """标记某步完成：清空其下游 steps 与 results 缓存（history 保留）。"""
        step = str(step)
        if step in STEP_ORDER:
            index = STEP_ORDER.index(step)
            downstream = set(STEP_ORDER[index + 1:])
            self.steps = [item for item in self.steps if item not in downstream and item != step]
            self.steps.append(step)
            for key in list(self.results.keys()):
                if key in downstream:
                    self.results.pop(key, None)
        else:
            if step not in self.steps:
                self.steps.append(step)
        self.current = step
        self.reorder_steps()

    def reorder_steps(self) -> None:
        """按 DESIGN §3 的顺序排列 steps（未知步骤排最后）。"""
        self.steps = ([item for item in STEP_ORDER if item in self.steps]
                      + [item for item in self.steps if item not in STEP_ORDER])

    # -- 写回（registry 消费步骤/分析返回值时调用）------------------------
    def store_result(self, key: str, value) -> None:
        """写入 dataset.results[key]（跨分析缓存）。"""
        if not key:
            return
        self.results[str(key)] = value

    def apply_replace(self, replace: dict | None) -> None:
        """步骤 replace：替换矩阵/基因/spot/counts/坐标/value_kind（保持形状一致）。"""
        if not replace:
            return
        if "X" in replace:
            arr = np.asarray(replace["X"], dtype=np.float32)
            if arr.ndim != 2:
                raise EngineError("internal", "replace.X 必须是二维（spot × 基因）")
            if "genes" in replace:
                genes = [str(g) for g in replace["genes"]]
                if len(genes) != arr.shape[1]:
                    raise EngineError("internal", "replace.genes 与矩阵列数不一致")
                self.genes = _dedupe(genes, "gene")
            if "spots" in replace:
                spots = [str(s) for s in replace["spots"]]
                if len(spots) != arr.shape[0]:
                    raise EngineError("internal", "replace.spots 与矩阵行数不一致")
                self.spots = _dedupe(spots, "spot")
                self.sample_info = self.sample_info.reindex(pd.Index(self.spots, name="spot"))
                self.sample_info.index = pd.Index(self.spots, name="spot")
                if self.coords is not None and self.coords.shape[0] != arr.shape[0]:
                    self.coords = None  # 行数变了但没给新坐标：坐标失效（调用方应同时 give coords）
                    self.geometry = estimate_geometry(None, self.coords_kind)
            self.X = arr
        if "counts" in replace:
            counts = replace["counts"]
            if counts is None:
                self.counts = None
            else:
                carr = np.asarray(counts, dtype=np.float32)
                if carr.shape != self.X.shape:
                    raise EngineError("internal", "replace.counts 与矩阵形状不一致")
                self.counts = carr
        if "coords" in replace and replace["coords"] is not None:
            carr = np.asarray(replace["coords"], dtype=np.float32)
            if carr.ndim != 2 or carr.shape[0] != self.n_spots:
                raise EngineError("internal", "replace.coords 行数与 spot 数不一致")
            self.coords = carr[:, :2]
            self.obsm["spatial"] = self.coords
            self.geometry = estimate_geometry(self.coords, replace.get("coords_kind") or self.coords_kind)
            self.neighbor_radius = _float(self.geometry.get("neighbor_radius")) or self.neighbor_radius
        if "coords_kind" in replace and replace["coords_kind"]:
            self.coords_kind = str(replace["coords_kind"])
            self.geometry = estimate_geometry(self.coords, self.coords_kind)
        if "value_kind" in replace and replace["value_kind"]:
            self.X_kind = self._resolve_value_kind(replace["value_kind"])
        if "X_kind" in replace and replace["X_kind"]:
            self.X_kind = self._resolve_value_kind(replace["X_kind"])
        if "sample_info" in replace and isinstance(replace["sample_info"], pd.DataFrame):
            self.sample_info = build_sample_info(self.spots, replace["sample_info"])
        if "genes" in replace and "X" not in replace:
            self.genes = _dedupe([str(g) for g in replace["genes"]], "gene")
        if "spots" in replace and "X" not in replace:
            self.spots = _dedupe([str(s) for s in replace["spots"]], "spot")

    def apply_sample_columns(self, columns, values) -> None:
        """分析/步骤返回 spot 级列（sample_columns + sample_values）写回 sample_info。"""
        if not columns or not isinstance(values, dict):
            return
        if isinstance(columns, dict):
            columns = [columns]
        for item in columns:
            key = str(item.get("key") if isinstance(item, dict) else item)
            if not key or key not in values:
                continue
            vals = values[key]
            if vals is None or len(vals) != self.n_spots:
                continue
            try:
                self.sample_info[key] = list(vals)
            except Exception:  # noqa: BLE001
                continue
            if key.lower() in _CLUSTER_LIKE_COLUMNS and key not in self.cluster_keys:
                self.cluster_keys.append(key)

    def apply_cluster_store(self, cluster_store: dict | None) -> None:
        """写回聚类标签：{"key":"clusters","labels":[...]}。"""
        if not cluster_store:
            return
        key = str(cluster_store.get("key") or DEFAULT_CLUSTER_KEY)
        labels = cluster_store.get("labels")
        if labels is None or len(labels) != self.n_spots:
            return
        self.sample_info[key] = [str(v) for v in labels]
        if key not in self.cluster_keys:
            self.cluster_keys.append(key)
        self.cluster_key = key

    def apply_coords_store(self, coords_store: dict | None) -> None:
        """写回派生二维坐标：{"key":"tsne"|"umap","values":[[x,y],...]}。"""
        if not coords_store:
            return
        key = str(coords_store.get("key") or "")
        values = coords_store.get("values")
        if not key or values is None:
            return
        arr = np.asarray(values, dtype=np.float32)
        if arr.ndim != 2 or arr.shape[0] != self.n_spots:
            return
        self.obsm[key] = arr

    def set_ref(self, ref: dict | None) -> None:
        """设置/清除反卷积参考（ref={path,cell_type_col,cell_types,X,genes}）。"""
        if ref is None:
            self.ref = None
            return
        out = dict(ref)
        if out.get("cell_types") is not None:
            out["cell_types"] = [str(item) for item in out["cell_types"]]
        self.ref = out

    # -- 杂项 -------------------------------------------------------------
    def describe(self) -> str:
        """一句话描述（日志/提示用）。"""
        return f"{self.id} {self.name}（{self.n_spots} spots × {self.n_genes} genes，{self.spatial_mode}）"

    def __repr__(self) -> str:  # pragma: no cover - 仅调试
        return f"<Dataset {self.describe()}>"


# ---------------------------------------------------------------------------
# State
# ---------------------------------------------------------------------------

class State:
    """多数据集注册表（DESIGN §2.2 list/select/close/reset）。"""

    def __init__(self):
        self.datasets: dict[str, Dataset] = {}
        self.active: str | None = None
        self.counter = 0

    def new(self, **kwargs) -> Dataset:
        """新建数据集：自动分配 id（ds1/ds2/...），无 active 时设为 active。"""
        self.counter += 1
        ds_id = str(kwargs.pop("id", None) or kwargs.pop("ds_id", None) or f"ds{self.counter}")
        if ds_id in self.datasets:
            ds_id = f"ds{self.counter}"
            while ds_id in self.datasets:
                self.counter += 1
                ds_id = f"ds{self.counter}"
        dataset = Dataset(ds_id, **kwargs)
        dataset.id = ds_id
        self.datasets[ds_id] = dataset
        if not self.active or self.active not in self.datasets:
            self.active = ds_id
        return dataset

    def add(self, dataset: Dataset) -> Dataset:
        """把已构造好的 Dataset 注册进来（id 冲突时重新编号）。"""
        if not dataset.id or dataset.id in self.datasets:
            self.counter += 1
            new_id = f"ds{self.counter}"
            while new_id in self.datasets:
                self.counter += 1
                new_id = f"ds{self.counter}"
            dataset.id = new_id
        self.datasets[dataset.id] = dataset
        if not self.active or self.active not in self.datasets:
            self.active = dataset.id
        return dataset

    def get(self, dataset_id: str | None = None) -> Dataset:
        """取数据集：None → active；不存在 → EngineError('not_found', 中文)。"""
        target = str(dataset_id).strip() if dataset_id else self.active
        if not target:
            raise EngineError("not_found", "还没有打开任何数据集：请先导入数据，或载入内置示例（demo）")
        dataset = self.datasets.get(target)
        if dataset is None:
            raise EngineError("not_found", f"数据集 {target} 不存在（可能已关闭）；"
                                          f"当前可用：{', '.join(self.datasets) or '无'}")
        return dataset

    def select(self, dataset_id: str | None = None) -> Dataset:
        """切换 active 数据集（等价 get + 设置 active）。"""
        dataset = self.get(dataset_id)
        self.active = dataset.id
        return dataset

    def summaries(self) -> list[dict]:
        """所有数据集摘要（DESIGN §2.2 list）。"""
        return [dataset.summary() for dataset in self.datasets.values()]

    def close(self, dataset_id: str) -> None:
        """关闭并移除数据集；关的是 active 时自动切到下一个。"""
        target = str(dataset_id or "").strip()
        if not target:
            target = self.active or ""
        if target not in self.datasets:
            raise EngineError("not_found", f"数据集 {target or '(空)'} 不存在")
        self.datasets.pop(target, None)
        if self.active == target:
            self.active = next(iter(self.datasets), None)

    def find_by_path(self, path) -> Dataset | None:
        """按导入路径查重（同一路径不重复打开）。"""
        if not path:
            return None
        try:
            target = os.path.abspath(os.path.expanduser(str(path)))
        except Exception:  # noqa: BLE001
            return None
        for dataset in self.datasets.values():
            if not dataset.path:
                continue
            try:
                if os.path.abspath(str(dataset.path)) == target:
                    return dataset
            except Exception:  # noqa: BLE001
                continue
        return None

    def reset(self) -> None:
        """清空全部数据集。"""
        self.datasets.clear()
        self.active = None

    def __len__(self) -> int:  # pragma: no cover - 便于日志
        return len(self.datasets)


# ---------------------------------------------------------------------------
# demo 数据（DESIGN §9）
# ---------------------------------------------------------------------------

#: demo 各空间域的标记基因（真实符号；不足 40 个时用 "EPI24" 这类编号补齐）
DOMAIN_MARKERS = [
    ("epithelial", (
        "EPCAM KRT8 KRT18 KRT19 KRT7 KRT17 MUC1 CDH1 LAMA3 LAMB3 ITGA6 ITGB4 S100A14 "
        "SERPINB5 CLDN4 CLDN7 TACSTD2 EMP2 EHF ELF3 GRHL2 OVOL1 IRF6 CLDN3"
    ).split()),
    ("proliferative", (
        "MKI67 TOP2A CCNB1 CCNB2 CDK1 AURKA AURKB PLK1 BUB1 BUB1B TYMS RRM2 PCNA "
        "MCM2 MCM3 MCM4 MCM5 MCM6 CCNE1 CCNE2 CDK2 E2F1 E2F3 NDC80"
    ).split()),
    ("stromal", (
        "COL1A1 COL1A2 COL3A1 COL5A2 ACTA2 TAGLN FN1 VIM S100A4 POSTN LUM DCN FBLN1 "
        "SPARC MMP2 TIMP1 TGFBI PCOLCE FSTL1 CTHRC1 IGFBP7 NNMT APOD CFD"
    ).split()),
    ("immune", (
        "CD3D CD3E CD3G CD2 CD8A CD8B GZMA GZMB NKG7 GNLY PRF1 IL7R CCL5 CD52 "
        "TRBC1 CTLA4 LAG3 PDCD1 TIGIT FOXP3 MS4A1 CD79A IGHM CD19"
    ).split()),
    ("vascular", (
        "PECAM1 VWF CDH5 KDR ENG ESAM RAMP2 PLVAP EMCN TIE1 CLDN5 VEGFA HIF1A LDHA "
        "PGK1 HK2 SLC2A1 ANGPT2 ESM1 NOS3 ADGRL4 EPAS1"
    ).split()),
]

#: demo 组织轮廓的域种子（归一化 uv 坐标）：5 个种子在单位圆内大致等分，得到 5 块面积相近的域
_DEMO_DOMAIN_SEEDS = np.array([
    [0.00, 0.45],
    [-0.43, 0.14],
    [-0.26, -0.36],
    [0.26, -0.36],
    [0.43, 0.14],
], dtype=float)

_DEMO_BASES = "ACGT"


def _hex_lattice(spacing: float, n_rows: int, n_cols: int) -> np.ndarray:
    """生成六边形格点（row 交错 0.5 spacing），返回 (n,2)。"""
    points = []
    dy = spacing * math.sqrt(3.0) / 2.0
    for row in range(n_rows):
        offset = 0.5 * spacing if row % 2 else 0.0
        for col in range(n_cols):
            points.append((col * spacing + offset, row * dy))
    return np.asarray(points, dtype=float)


def _demo_slice(rng, n_target: int, spacing: float, rotation: float = 0.0,
                offset=(0.0, 0.0)) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """造一个不规则组织轮廓的切片：返回 (coords[n,2], domain_labels[n], u_normalized[n])。"""
    n_rows = max(8, int(round(math.sqrt(n_target * 1.55 / 1.15))))
    n_cols = max(8, int(math.ceil(n_target * 1.55 / n_rows)))
    pts = _hex_lattice(spacing, n_rows, n_cols)
    xmin, ymin = pts[:, 0].min(), pts[:, 1].min()
    xmax, ymax = pts[:, 0].max(), pts[:, 1].max()
    cx, cy = (xmin + xmax) / 2.0, (ymin + ymax) / 2.0
    a = max((xmax - xmin) / 2.0, spacing)
    b = max((ymax - ymin) / 2.0, spacing)
    u = (pts[:, 0] - cx) / a
    v = (pts[:, 1] - cy) / b
    theta = np.arctan2(v, u)
    # 不规则轮廓：椭圆 + 低频径向扰动（相位由 rng 决定：两个切片轮廓不同但整体确定）
    phase = float(rng.uniform(0.0, 2.0 * math.pi))
    wobble = (0.11 * np.sin(3.0 * theta + phase * 0.7) + 0.07 * np.cos(5.0 * theta + 2.10)
              + 0.04 * np.sin(7.0 * theta + phase))
    field = 1.0 - (u ** 2 + v ** 2) + wobble
    # 取 field 最大的 n_target 个点（并列时按原始顺序稳定排序）
    order = np.argsort(-(field + 1e-9 * np.arange(field.size)), kind="stable")
    keep = np.sort(order[:n_target])
    coords = pts[keep].copy()
    u_keep, v_keep = u[keep], v[keep]
    if rotation:
        cos_r, sin_r = math.cos(rotation), math.sin(rotation)
        u_keep, v_keep = (u_keep * cos_r - v_keep * sin_r,
                          u_keep * sin_r + v_keep * cos_r)
    # 域：5 个种子的最近邻（块状/分层分布）
    d2 = ((u_keep[:, None] - _DEMO_DOMAIN_SEEDS[None, :, 0]) ** 2
          + (v_keep[:, None] - _DEMO_DOMAIN_SEEDS[None, :, 1]) ** 2)
    domains = np.argmin(d2, axis=1).astype(np.int16)
    coords[:, 0] += offset[0]
    coords[:, 1] += offset[1]
    return coords, domains, u_keep


def _demo_gene_names() -> tuple[list[str], list[list[int]], np.ndarray]:
    """1000 个基因名 + 每域 40 个标记基因下标 + 每域标记的 log2FC。"""
    rng = np.random.default_rng(DEMO_SEED + 7)
    genes: list[str] = []
    seen: set[str] = set()
    for domain, (tag, markers) in enumerate(DOMAIN_MARKERS):
        need = DEMO_MARKERS_PER_DOMAIN
        names = list(dict.fromkeys(markers))[:need]
        index = 1
        while len(names) < need:
            name = f"{tag[:3].upper()}{index:02d}"
            index += 1
            if name not in names:
                names.append(name)
        for name in names:
            if name in seen:
                continue
            genes.append(name)
            seen.add(name)
    filler = 1
    while len(genes) < DEMO_N_GENES:
        name = f"GENE{filler:04d}"
        filler += 1
        if name in seen:
            continue
        genes.append(name)
        seen.add(name)
    genes = genes[:DEMO_N_GENES]
    lookup = {name: idx for idx, name in enumerate(genes)}
    marker_idx: list[list[int]] = []
    for domain, (_, markers) in enumerate(DOMAIN_MARKERS):
        need = DEMO_MARKERS_PER_DOMAIN
        names = list(dict.fromkeys(markers))[:need]
        index = 1
        tag = DOMAIN_MARKERS[domain][0]
        while len(names) < need:
            name = f"{tag[:3].upper()}{index:02d}"
            index += 1
            if name not in names:
                names.append(name)
        marker_idx.append([lookup[name] for name in names if name in lookup])
    logfc = [rng.uniform(1.2, 2.5, len(items)) for items in marker_idx]
    return genes, marker_idx, logfc


def _demo_barcodes(rng, n: int, suffix: str) -> list[str]:
    """生成 10x 风格的 spot barcode（16 个 ACGT + "-<slice>"）。"""
    out: list[str] = []
    seen: set[str] = set()
    while len(out) < n:
        code = "".join(rng.choice(list(_DEMO_BASES), size=16))
        if code in seen:
            continue
        seen.add(code)
        out.append(f"{code}-{suffix}")
    return out


def make_demo(state: State, name: str | None = None) -> Dataset:
    """内置示例数据（DESIGN §9）。

    合成一块不规则组织：2400 spots（六边形格点）、5 个块状空间域、
    1000 个基因（每域 40 个标记基因，log2FC 1.2–2.5）、2 个切片、2 个分组。
    确定性：seed=0；生成时间 < 3 秒（2400×1000 的 Poisson 抽样约 0.1 s）。
    """
    rng = np.random.default_rng(DEMO_SEED)
    spacing = 100.0
    per_slice = DEMO_N_SPOTS // 2
    genomes, marker_idx, marker_logfc = _demo_gene_names()
    n_genes = len(genomes)

    coords_all: list[np.ndarray] = []
    domain_all: list[np.ndarray] = []
    u_all: list[np.ndarray] = []
    matrix_parts: list[np.ndarray] = []
    barcodes_all: list[str] = []
    slice_all: list[str] = []
    group_all: list[str] = []

    # 两个切片：slice1=tumor（无旋转），slice2=adjacent（旋转 20° + 右侧平移），
    # group 与 slice 一一对应（2 组），这样 de_groups / 多切片诊断都可用。
    slice_plan = [
        ("slice1", "tumor", 0.0, (0.0, 0.0), 0),
        ("slice2", "adjacent", 0.35, (3000.0, 420.0), 1),
    ]
    # 批次基因（slice2 整体上移），让多切片/分组比较有信号
    batch_idx = np.arange(300, 400)
    gradient_idx = np.arange(200, 300)

    for suffix, (slice_name, group_name, rotation, offset, batch_sign) in enumerate(slice_plan, start=1):
        coords, domains, u = _demo_slice(rng, per_slice, spacing, rotation=rotation, offset=offset)
        lib = rng.lognormal(mean=0.0, sigma=0.16, size=coords.shape[0])
        mu = np.exp(rng.normal(1.55, 0.72, size=(coords.shape[0], n_genes))) * lib[:, None]
        # 每域 40 个标记基因（log2FC 1.2–2.5）
        for domain in range(DEMO_N_DOMAINS):
            rows = np.where(domains == domain)[0]
            if rows.size == 0:
                continue
            ids = marker_idx[domain]
            mu[np.ix_(rows, ids)] *= (2.0 ** marker_logfc[domain])[None, :]
        # 100 个基因沿组织切向（u）单调升高：供 spatial_gradient 分析使用
        u_norm = (u - np.nanmin(u)) / max(1e-9, float(np.nanmax(u) - np.nanmin(u)))
        mu[:, gradient_idx] *= np.exp(0.9 * (u_norm[:, None] - 0.5))
        # 切片批次效应
        if batch_sign:
            mu[:, batch_idx] *= 1.55
        counts = rng.poisson(mu).astype(np.float32)
        matrix_parts.append(counts)
        coords_all.append(coords.astype(np.float32))
        domain_all.append(domains)
        u_all.append(u)
        barcodes_all.extend(_demo_barcodes(rng, coords.shape[0], str(suffix)))
        slice_all.extend([slice_name] * coords.shape[0])
        group_all.extend([group_name] * coords.shape[0])

    counts = np.concatenate(matrix_parts, axis=0).astype(np.float32)
    coords = np.concatenate(coords_all, axis=0).astype(np.float32)
    domains = np.concatenate(domain_all, axis=0)
    sample_info = pd.DataFrame(
        {
            "slice": slice_all,
            "group": group_all,
            "domain_true": [f"D{int(v) + 1}" for v in domains],
        },
        index=pd.Index(barcodes_all, name="spot"),
    )
    notes = [
        "内置示例数据：2400 spots × 1000 基因，2 个切片（slice1/slice2）× 2 个分组（tumor/adjacent），"
        "5 个块状空间域，每域 40 个标记基因（log2FC 1.2–2.5）",
        "坐标为合成六边形格点（间距 100，像素量纲），无 H&E 底图",
    ]
    dataset = state.new(
        name=name or "demo-spatial",
        X=counts,
        counts=counts,
        genes=genomes,
        spots=barcodes_all,
        coords=coords,
        coords_kind="grid",
        spatial_mode="grid",
        sample_info=sample_info,
        group_key="group",
        species="human",
        X_kind="counts",
        path=None,
        source_format="demo",
        notes=notes,
    )
    # demo 真值（供测试/分析对照，不参与协议）
    dataset.results["demo_truth"] = {
        "domain": [int(v) for v in domains],
        "slice": list(slice_all),
        "group": list(group_all),
        "marker_genes": {str(d): [genomes[i] for i in marker_idx[d]] for d in range(DEMO_N_DOMAINS)},
        "marker_logfc": {str(d): [float(v) for v in marker_logfc[d]] for d in range(DEMO_N_DOMAINS)},
        "gradient_genes": [genomes[i] for i in gradient_idx],
        "batch_genes": [genomes[i] for i in batch_idx],
    }
    return dataset


__all__ = [
    "Dataset", "State", "estimate_geometry", "build_sample_info", "make_demo",
    "column_kind", "now_iso", "STEP_ORDER", "COORDS_KINDS", "SPATIAL_MODES",
    "DEMO_N_SPOTS", "DEMO_N_GENES", "DEMO_N_DOMAINS", "DEMO_MARKERS_PER_DOMAIN",
]
