"""
dsh-spatial — 全格式导入层（DESIGN.md §6）

对外接口（sp_registry 逐字调用）：
  * SUPPORTED_FORMATS: list[str]
  * FORMAT_INFO: list[dict]                      # [{"id","label","hint"}]
  * load_source(path_or_url, workdir=None, options=None) -> dict
  * inspect_source(path_or_url, workdir=None, options=None) -> dict

load_source 返回（DESIGN §6）：
  {"format", "X", "genes", "spots", "coords", "sample_info", "image",
   "notes", "warnings", "meta", "X_kind", "counts"?}

支持的空间平台：
  visium       10x Visium `outs/` 目录（h5/mtx + spatial/tissue_positions*.csv + scalefactors + 图）
  visium_h5    单个 `filtered_feature_bc_matrix.h5`（需配 coords= 坐标 csv）
  h5ad         AnnData（obsm['spatial'] / obs 的 pxl_row|pxl_col|x|y）
  stereo       Stereo-seq BIN 长表（geneID,x,y,count）
  slideseq     Slide-seq puck 矩阵（需配 coords=）
  matrix_text  csv/tsv/txt(.gz) 表达矩阵（需配 coords=）
  xlsx         Excel（需配 coords=）
  geo_series   GEO series matrix（表段基因×样本；无坐标 → 只有表达）
  archive      zip/tar/tar.gz（解包后自动找 Visium/BIN/矩阵）
  rds/rdata    R 对象（经 Rscript 桥接；R 不可用给 unsupported）
  coords_csv   —— 不是数据格式，而是 `options["coords"]` 的坐标表格式

约定：
  * 只依赖标准库 + numpy/pandas(+h5py/anndata/openpyxl/PIL 若已装)，禁止新增第三方依赖。
  * 任何失败都给中文可读原因；不确定（多个候选）时抛 EngineError("ambiguous", candidates=[...])。
"""

from __future__ import annotations

import gzip
import io
import json
import os
import re
import shutil
import tarfile
import tempfile
import uuid
import zipfile

import numpy as np
import pandas as pd

from sp_common import EngineError, bad_request, ensure_dir, prerequisite

# ---------------------------------------------------------------------------
# 格式清单
# ---------------------------------------------------------------------------

SUPPORTED_FORMATS = [
    "auto", "visium", "visium_h5", "h5ad", "stereo", "slideseq",
    "matrix_text", "xlsx", "geo_series", "archive", "rds", "rdata",
]

FORMAT_INFO = [
    {"id": "auto", "label": "自动识别", "hint": "按路径/文件内容猜；目录优先当 Visium outs/"},
    {"id": "visium", "label": "10x Visium outs/", "hint": "含 spatial/ 与 filtered_feature_bc_matrix 的目录"},
    {"id": "visium_h5", "label": "Visium .h5", "hint": "单文件矩阵，需另给 coords 坐标表"},
    {"id": "h5ad", "label": "AnnData .h5ad", "hint": "obsm['spatial'] 或 obs 的 x/y、pxl_row/pxl_col"},
    {"id": "stereo", "label": "Stereo-seq BIN", "hint": "*_raw_bin_count.tsv.gz（geneID,x,y,count）"},
    {"id": "slideseq", "label": "Slide-seq puck", "hint": "puck 矩阵 + coords 坐标表"},
    {"id": "matrix_text", "label": "表达矩阵 csv/tsv", "hint": "需配 coords（spot,x,y）"},
    {"id": "xlsx", "label": "Excel 矩阵", "hint": "需配 coords；可用 sheet/gene_col"},
    {"id": "geo_series", "label": "GEO series matrix", "hint": "GSE 下载的系列矩阵；无坐标时只能当表达矩阵"},
    {"id": "archive", "label": "压缩包 zip/tar", "hint": "解包后自动找 Visium/BIN/矩阵"},
    {"id": "rds", "label": "R .rds", "hint": "Seurat/SCE/矩阵（需要 Rscript）"},
    {"id": "rdata", "label": "R .RData/.rda", "hint": "同上，可用 object 指定对象名"},
]

_GENE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._\-]{0,24}$")
_MT_RE = re.compile(r"^(MT-|mt-|Mt-)")

VISIUM_POSITION_NAMES = (
    "tissue_positions_list.csv", "tissue_positions.csv", "tissue_positions_list.txt",
)
VISIUM_MATRIX_H5 = ("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5")


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------

def _is_url(value: str) -> bool:
    return str(value or "").startswith(("http://", "https://", "ftp://"))


def _read_text(path: str, limit: int | None = None) -> str:
    """读文本（自动 gzip）。limit=None 读全部。"""
    try:
        if path.lower().endswith((".gz", ".gzip")):
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                return handle.read() if limit is None else handle.read(limit)
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            return handle.read() if limit is None else handle.read(limit)
    except OSError as error:
        raise EngineError("bad_request", f"读取文件失败：{os.path.basename(path)}（{error}）")


def _sniff_delimiter(sample: str) -> str:
    first = ""
    for line in sample.splitlines():
        if line.strip():
            first = line
            break
    if not first:
        return "\t"
    counts = {sep: first.count(sep) for sep in ("\t", ",", ";", "|")}
    best = max(counts, key=counts.get)
    if counts[best] == 0:
        return r"\s+"
    return best


def _read_table(path: str, options: dict | None = None, nrows: int | None = None) -> pd.DataFrame:
    """读 csv/tsv/txt（含 gz）。分隔符自动嗅探，也可 options['delimiter'] 指定。"""
    options = options or {}
    delimiter = options.get("delimiter")
    if delimiter in ("tab", "\\t"):
        delimiter = "\t"
    if not delimiter:
        head = _read_text(path, limit=256 * 1024)
        delimiter = _sniff_delimiter(head)
    try:
        return pd.read_csv(path, sep=delimiter, index_col=0, nrows=nrows,
                           engine="python", compression="infer", on_bad_lines="skip")
    except Exception as error:  # noqa: BLE001
        raise EngineError("bad_request", f"解析表格失败：{os.path.basename(path)}（{error}）")


def _numeric_frame(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """只保留数值列；返回 (数值帧, 丢弃的非数值行数)。"""
    out = frame.copy()
    out.columns = [str(col).strip() for col in out.columns]
    converted = out.apply(pd.to_numeric, errors="coerce")
    bad_rows = int(converted.isna().all(axis=1).sum())
    converted = converted.loc[~converted.isna().all(axis=1)]
    keep = [col for col in converted.columns if converted[col].notna().any()]
    dropped_cols = len(converted.columns) - len(keep)
    return converted[keep], bad_rows + dropped_cols


def _looks_like_genes(labels) -> float:
    labels = [str(item) for item in labels if str(item).strip()]
    if not labels:
        return 0.0
    hits = sum(1 for item in labels if _GENE_RE.match(item))
    return hits / len(labels)


def _looks_like_spots(labels) -> float:
    labels = [str(item) for item in labels if str(item).strip()]
    if not labels:
        return 0.0
    hits = sum(1 for item in labels if re.search(r"[-_:]", item) or re.fullmatch(r"[ACGTN]{12,}", item))
    return hits / len(labels)


def _value_kind(matrix: np.ndarray) -> str:
    arr = np.asarray(matrix, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return "unknown"
    if finite.min() >= 0 and np.allclose(finite, np.round(finite)):
        return "counts"
    if finite.max() < 60 and finite.min() >= 0:
        return "lognorm"
    return "unknown"


def _dedupe(names, prefix: str) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for index, name in enumerate(names):
        base = str(name).strip() or f"{prefix}{index + 1}"
        if base in seen:
            seen[base] += 1
            out.append(f"{base}-{seen[base]}")
        else:
            seen[base] = 0
            out.append(base)
    return out


# ---------------------------------------------------------------------------
# 坐标表
# ---------------------------------------------------------------------------

_COORD_PATTERNS = [
    ("barcode", "x", "y"),
    ("spot", "x", "y"),
    ("cell", "x", "y"),
    ("id", "x", "y"),
]


def read_coords(path: str, options: dict | None = None) -> tuple[list[str], np.ndarray, dict]:
    """读坐标表：返回 (labels, coords[n,2], meta)。

    支持列名：barcode/spot/cell/id + x/y（或 pxl_col/pxl_row、col/row、array_col/array_row）。
    """
    options = options or {}
    delimiter = options.get("delimiter")
    if delimiter in ("tab", "\\t"):
        delimiter = "\t"
    if not delimiter:
        delimiter = _sniff_delimiter(_read_text(path, limit=256 * 1024))
    try:
        # 坐标表必须保留 barcode 列（不能当索引吞掉）
        frame = pd.read_csv(path, sep=delimiter, engine="python", compression="infer",
                            on_bad_lines="skip")
    except Exception as error:  # noqa: BLE001
        raise EngineError("bad_request", f"解析坐标表失败：{os.path.basename(path)}（{error}）")
    if "Unnamed: 0" in frame.columns:
        frame = frame.rename(columns={"Unnamed: 0": "barcode"})
    frame.columns = [str(col).strip().lower() for col in frame.columns]
    label_col = None
    for candidate in ("barcode", "spot", "cell", "id", "barcodes", "sample"):
        if candidate in frame.columns:
            label_col = candidate
            break
    if label_col is None:
        label_col = frame.columns[0]
    pairs = [
        ("x", "y"), ("pxl_col_in_fullres", "pxl_row_in_fullres"),
        ("pxl_col", "pxl_row"), ("array_col", "array_row"), ("col", "row"),
        ("x_coord", "y_coord"), ("imagecol", "imagerow"),
    ]
    xcol = ycol = None
    for cx, cy in pairs:
        if cx in frame.columns and cy in frame.columns:
            xcol, ycol = cx, cy
            break
    if xcol is None:
        # 退而求其次：两个数值列
        numeric = [col for col in frame.columns if col != label_col and pd.api.types.is_numeric_dtype(frame[col])]
        if len(numeric) >= 2:
            xcol, ycol = numeric[0], numeric[1]
    if xcol is None or ycol is None:
        raise prerequisite(
            f"坐标表 {os.path.basename(path)} 里找不到 x/y 列（可用列：{', '.join(frame.columns[:12])}）")
    labels = [str(item) for item in frame[label_col].tolist()]
    coords = np.asarray(frame[[xcol, ycol]].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32))
    keep = np.isfinite(coords).all(axis=1)
    labels = [label for label, ok in zip(labels, keep) if ok]
    coords = coords[keep]
    if not len(labels):
        raise prerequisite(f"坐标表 {os.path.basename(path)} 没有有效行")
    meta = {"coords_kind": "pixel" if "pxl" in xcol or xcol in ("imagecol",) else "generic",
            "x_col": xcol, "y_col": ycol, "label_col": label_col}
    return labels, coords, meta


# ---------------------------------------------------------------------------
# 10x h5 / mtx
# ---------------------------------------------------------------------------

def _read_10x_h5(path: str) -> tuple[np.ndarray, list[str], list[str]]:
    """读 10x `.h5`（v3：/matrix/data|indices|indptr|shape + features/name + barcodes）。"""
    try:
        import h5py
    except Exception as error:  # noqa: BLE001
        raise EngineError("unsupported", f"需要 h5py 读取 10x .h5（{error}）")
    try:
        from scipy import sparse
    except Exception as error:  # noqa: BLE001
        raise EngineError("internal", f"需要 scipy 读取稀疏矩阵（{error}）")
    with h5py.File(path, "r") as handle:
        group = handle["matrix"] if "matrix" in handle else handle
        data = np.asarray(group["data"][:])
        indices = np.asarray(group["indices"][:])
        indptr = np.asarray(group["indptr"][:])
        shape = tuple(int(x) for x in np.asarray(group["shape"][:]).ravel()[:2])
        matrix = sparse.csc_matrix((data, indices, indptr), shape=shape)

        def _names(*keys):
            for key in keys:
                node = group
                for part in key.split("/"):
                    if part not in node:
                        node = None
                        break
                    node = node[part]
                if node is None:
                    continue
                try:
                    raw = node[:]
                except Exception:  # noqa: BLE001
                    continue
                out = []
                for item in raw:
                    out.append(item.decode("utf-8") if isinstance(item, bytes) else str(item))
                if out:
                    return out
            return []

        genes = _names("features/name", "features/id", "genes", "gene_names")
        spots = _names("barcodes", "spot_names")
    n_genes, n_spots = shape
    if not genes:
        genes = [f"GENE{i + 1}" for i in range(n_genes)]
    if not spots:
        spots = [f"SPOT{i + 1}" for i in range(n_spots)]
    dense = np.asarray(matrix.todense(), dtype=np.float32).T  # → spots × genes
    return dense, genes, spots


def _read_10x_mtx_dir(directory: str) -> tuple[np.ndarray, list[str], list[str]]:
    """读 10x `filtered_feature_bc_matrix/`（matrix.mtx(.gz) + features.tsv(.gz) + barcodes.tsv(.gz)）。"""
    try:
        from scipy.io import mmread
        from scipy import sparse
    except Exception as error:  # noqa: BLE001
        raise EngineError("internal", f"需要 scipy 读取 mtx（{error}）")
    matrix_path = None
    for name in ("matrix.mtx.gz", "matrix.mtx"):
        candidate = os.path.join(directory, name)
        if os.path.exists(candidate):
            matrix_path = candidate
            break
    if not matrix_path:
        raise EngineError("not_found", f"{directory} 里没有 matrix.mtx")
    matrix = mmread(matrix_path)
    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    matrix = matrix.tocsc()
    features = None
    for name in ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"):
        candidate = os.path.join(directory, name)
        if os.path.exists(candidate):
            features = _read_text(candidate).splitlines()
            break
    barcodes = None
    for name in ("barcodes.tsv.gz", "barcodes.tsv"):
        candidate = os.path.join(directory, name)
        if os.path.exists(candidate):
            barcodes = _read_text(candidate).splitlines()
            break
    genes = [line.split("\t")[1] if line.count("\t") >= 1 else line.split("\t")[0]
             for line in (features or []) if line.strip()]
    spots = [line.split("\t")[0] for line in (barcodes or []) if line.strip()]
    n_genes, n_spots = matrix.shape
    dense = np.asarray(matrix.todense(), dtype=np.float32).T
    if not genes:
        genes = [f"GENE{i + 1}" for i in range(n_genes)]
    if not spots:
        spots = [f"SPOT{i + 1}" for i in range(n_spots)]
    return dense, genes, spots


# ---------------------------------------------------------------------------
# Visium
# ---------------------------------------------------------------------------

def _find_visium_matrix(root: str) -> tuple[str, str]:
    """在 Visium 目录里找矩阵：返回 (kind, path)。kind ∈ h5|mtx。"""
    for name in VISIUM_MATRIX_H5:
        candidate = os.path.join(root, name)
        if os.path.exists(candidate):
            return "h5", candidate
    for name in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix",
                 "filtered_gene_bc_matrices", "raw_gene_bc_matrices"):
        candidate = os.path.join(root, name)
        if os.path.isdir(candidate) and (
                os.path.exists(os.path.join(candidate, "matrix.mtx.gz"))
                or os.path.exists(os.path.join(candidate, "matrix.mtx"))):
            return "mtx", candidate
        # v1 布局：filtered_gene_bc_matrices/<genome>/matrix.mtx
        if os.path.isdir(candidate):
            for sub in sorted(os.listdir(candidate)):
                subdir = os.path.join(candidate, sub)
                if os.path.isdir(subdir) and (
                        os.path.exists(os.path.join(subdir, "matrix.mtx.gz"))
                        or os.path.exists(os.path.join(subdir, "matrix.mtx"))):
                    return "mtx", subdir
    # 上级 outs/ 目录情况：root 可能是 spatial/ 本身
    parent = os.path.dirname(root.rstrip(os.sep))
    if parent and parent != root:
        try:
            return _find_visium_matrix(parent)
        except EngineError:
            pass
    raise EngineError("not_found",
                      f"{root} 里没找到 filtered_feature_bc_matrix.h5 / filtered_feature_bc_matrix/")


def _find_spatial_dir(root: str) -> str | None:
    candidate = os.path.join(root, "spatial")
    if os.path.isdir(candidate):
        return candidate
    if os.path.basename(root.rstrip(os.sep)) == "spatial":
        return root
    parent = os.path.dirname(root.rstrip(os.sep))
    candidate = os.path.join(parent, "spatial")
    if os.path.isdir(candidate):
        return candidate
    return None


#: Visium 图像 → （scalefactor 键, 坐标语义）
_VISIUM_IMAGE_PLAN = (
    ("tissue_hires_image.png", "tissue_hires_scalef", "hires", "visium_hires"),
    ("tissue_lowres_image.png", "tissue_lowres_scalef", "lowres", "visium_lowres"),
    ("tissue_image.png", None, "custom", "pixel"),
)


def _visium_positions_file(spatial_dir: str) -> tuple[str, str]:
    """找 tissue_positions 文件并判定 Visium 输出版本：返回 (path, version)。

    * v2（Space Ranger ≥ 2.0）：``tissue_positions.csv``（带表头）；
    * v1（Space Ranger 1.x）：``tissue_positions_list.csv``（无表头，6 列固定顺序）。
    """
    for name, version in (("tissue_positions.csv", "v2"),):
        candidate = os.path.join(spatial_dir, name)
        if os.path.exists(candidate):
            return candidate, version
    for name in ("tissue_positions_list.csv", "tissue_positions_list.txt"):
        candidate = os.path.join(spatial_dir, name)
        if os.path.exists(candidate):
            head = _read_text(candidate, limit=4096)
            first = next((line for line in head.splitlines() if line.strip()), "")
            # 少数导出工具会在 v1 文件上也加表头：有表头就按 v2 语义解析
            version = "v2" if "barcode" in first.lower() else "v1"
            return candidate, version
    raise EngineError("not_found", f"{spatial_dir} 里没有 tissue_positions*.csv")
    # pragma: no cover - 上面必然 raise 或 return


def _read_visium_positions(spatial_dir: str, scale_key: str | None = "tissue_hires_scalef",
                           scale_override: float | None = None) -> tuple[list[str], np.ndarray, dict]:
    """读 tissue_positions*.csv 并校验 Visium 输出版本/列。

    坐标 = ``pxl_col/row_in_fullres × 实际使用的那张图的 scalefactor``：
    scalefactor 由调用方按「真正会叠加的组织图」传入（hires → tissue_hires_scalef，
    lowres → tissue_lowres_scalef），**不得**在只有低分辨率图时误乘 hires 比例。

    非法输出（缺 pxl_row/col 列、v1 列数不足、坐标为全 NaN）→ ``bad_request`` 可读中文。
    """
    path, version = _visium_positions_file(spatial_dir)
    head = _read_text(path, limit=4096)
    first = next((line for line in head.splitlines() if line.strip()), "")
    has_header = "barcode" in first.lower()
    if has_header:
        try:
            frame = pd.read_csv(path)
        except Exception as error:  # noqa: BLE001
            raise EngineError("bad_request", f"读取 Visium 坐标表失败（{os.path.basename(path)}）：{error}")
        frame.columns = [str(col).strip() for col in frame.columns]
        required = ("barcode", "pxl_row_in_fullres", "pxl_col_in_fullres")
        missing = [col for col in required if col not in frame.columns]
        if missing:
            raise EngineError(
                "bad_request",
                f"Visium 输出版本异常（{os.path.basename(path)}，判定为 {version}）：缺列 "
                + "、".join(missing)
                + "。请确认给的是 Space Ranger outs/ 目录（spatial/tissue_positions*.csv 需含 "
                  "pxl_row_in_fullres/pxl_col_in_fullres）")
    else:
        columns = [col for col in re.split(r"[\t,]", first.strip()) if col != ""]
        if len(columns) < 6:
            raise EngineError(
                "bad_request",
                f"Visium 输出版本异常（{os.path.basename(path)}，判定为 {version}）："
                f"无表头的坐标表应有 6 列（barcode,in_tissue,array_row,array_col,"
                f"pxl_row_in_fullres,pxl_col_in_fullres），实际 {len(columns)} 列")
        frame = pd.read_csv(path, header=None,
                            names=["barcode", "in_tissue", "array_row", "array_col",
                                   "pxl_row_in_fullres", "pxl_col_in_fullres"])
    frame["barcode"] = frame["barcode"].astype(str)
    if "in_tissue" in frame.columns:
        in_tissue = pd.to_numeric(frame["in_tissue"], errors="coerce").fillna(0)
        frame = frame.loc[in_tissue > 0]
    # scalefactor：显式覆盖 > 实际使用图像的 scalefactor 键 > 不缩放（1.0，并写 warning）
    warnings: list[str] = []
    factors: dict = {}
    scalefactors_path = os.path.join(spatial_dir, "scalefactors_json.json")
    if os.path.exists(scalefactors_path):
        try:
            with open(scalefactors_path, "r", encoding="utf-8") as handle:
                factors = dict(json.load(handle))
        except Exception:  # noqa: BLE001
            factors = {}
            warnings.append("scalefactors_json.json 解析失败：坐标按 fullres 像素使用（未缩放）")
    scale = 1.0
    applied_key = None
    scale_applied = False
    if scale_override is not None:
        try:
            candidate = float(scale_override)
        except (TypeError, ValueError):
            candidate = 0.0
        if np.isfinite(candidate) and candidate > 0:
            scale = candidate
            scale_applied = True
            applied_key = "visium_scale(参数覆盖)"
    elif scale_key and factors:
        try:
            candidate = float(factors.get(scale_key) or 0)
        except (TypeError, ValueError):
            candidate = 0.0
        if np.isfinite(candidate) and candidate > 0:
            scale = candidate
            scale_applied = True
            applied_key = scale_key
        else:
            warnings.append(
                f"scalefactors_json.json 里缺有效的 {scale_key}：坐标按 fullres 像素使用（未缩放）")
    elif not factors:
        warnings.append(
            "没有找到可用的 scalefactors_json.json：坐标按 fullres 像素使用（未缩放）")
    x = pd.to_numeric(frame.get("pxl_col_in_fullres"), errors="coerce").to_numpy(dtype=np.float32) * scale
    y = pd.to_numeric(frame.get("pxl_row_in_fullres"), errors="coerce").to_numpy(dtype=np.float32) * scale
    keep = np.isfinite(x) & np.isfinite(y)
    if not keep.any():
        raise EngineError("bad_request",
                          f"Visium 坐标表里的 pxl_row/pxl_col_in_fullres 全是非数值（{os.path.basename(path)}）")
    frame = frame.loc[keep]
    coords = np.column_stack([x[keep], y[keep]]).astype(np.float32)
    if scale_applied and applied_key == "tissue_lowres_scalef":
        coords_kind = "visium_lowres"
    elif scale_applied and applied_key == "tissue_hires_scalef":
        coords_kind = "visium_hires"
    elif scale_applied:
        coords_kind = "visium_hires" if applied_key == "visium_scale(参数覆盖)" else "pixel"
    else:
        coords_kind = "pixel"
    meta = {"coords_kind": coords_kind, "scale": float(scale), "scale_key": applied_key,
            "scale_applied": bool(scale_applied), "visium_version": version,
            "positions_file": os.path.basename(path),
            "warnings": warnings,
            "array_row": frame.get("array_row", pd.Series(dtype=float)).tolist(),
            "array_col": frame.get("array_col", pd.Series(dtype=float)).tolist()}
    return [str(item) for item in frame["barcode"].tolist()], coords, meta


def _find_tissue_image(spatial_dir: str) -> tuple[str | None, str | None, float | None, str | None]:
    """按实际会叠加的组织图返回 ``(path, scale_key, scale, kind)``。

    注意：只有低分辨率图时返回 lowres 的 scalefactor —— 坐标必须按低分辨率图缩放，
    不能误用 hires 比例（否则坐标被放大 5×）。
    """
    scalefactors_path = os.path.join(spatial_dir, "scalefactors_json.json")
    factors: dict = {}
    if os.path.exists(scalefactors_path):
        try:
            with open(scalefactors_path, "r", encoding="utf-8") as handle:
                factors = dict(json.load(handle))
        except Exception:  # noqa: BLE001
            factors = {}
    for name, scale_key, kind, _coords_kind in _VISIUM_IMAGE_PLAN:
        candidate = os.path.join(spatial_dir, name)
        if not os.path.exists(candidate):
            continue
        scale = None
        if scale_key:
            try:
                scale = float(factors.get(scale_key) or 0) or None
            except (TypeError, ValueError):
                scale = None
        return candidate, scale_key, scale, kind
    return None, None, None, None


def _image_extent(path: str | None, coords: np.ndarray) -> list[float] | None:
    if not path or coords.size == 0:
        return None
    try:
        from PIL import Image

        with Image.open(path) as handle:
            width, height = handle.size
        return [0.0, float(width), float(height), 0.0]
    except Exception:  # noqa: BLE001
        return None


def load_visium(root: str, options: dict | None = None) -> dict:
    """读 10x Visium `outs/` 目录（或 spatial/、或含矩阵的目录）。"""
    options = options or {}
    matrix_kind, matrix_path = _find_visium_matrix(root)
    if matrix_kind == "h5":
        dense, genes, spots = _read_10x_h5(matrix_path)
    else:
        dense, genes, spots = _read_10x_mtx_dir(matrix_path)
    spatial_dir = _find_spatial_dir(root) or _find_spatial_dir(os.path.dirname(matrix_path))
    notes: list[str] = []
    warnings: list[str] = []
    coords = None
    coords_kind = "visium_hires"
    image = None
    visium_version = None
    position_meta: dict = {}
    if spatial_dir:
        # 先确定「实际会叠加的组织图」→ 用那张图的 scalefactor 缩放坐标
        image_path, image_scale_key, image_scale, image_kind = _find_tissue_image(spatial_dir)
        scale_key = image_scale_key or "tissue_hires_scalef"
        scale_override = options.get("visium_scale")
        labels, raw_coords, meta = _read_visium_positions(
            spatial_dir, scale_key=scale_key, scale_override=scale_override)
        coords_kind = meta.get("coords_kind", "visium_hires")
        visium_version = meta.get("visium_version")
        position_meta = dict(meta)
        notes.append(
            f"Visium 输出版本 {visium_version}（{meta.get('positions_file')}）；"
            f"坐标缩放键 {meta.get('scale_key') or '无（fullres）'} = {float(meta.get('scale') or 1.0):.4f}")
        warnings.extend(meta.get("warnings") or [])
        position = {str(label): index for index, label in enumerate(labels)}
        coords = np.full((len(spots), 2), np.nan, dtype=np.float32)
        matched = 0
        for index, spot in enumerate(spots):
            pos = position.get(str(spot))
            if pos is not None:
                coords[index] = raw_coords[pos]
                matched += 1
        if matched == 0:
            coords = None
            warnings.append("矩阵 barcode 与 tissue_positions 的 barcode 全都不匹配（可能不是同一份 outs）")
        else:
            if matched < len(spots):
                warnings.append(f"{len(spots) - matched} 个 barcode 没有坐标（已被坐标分析排除，表达矩阵仍保留）")
            notes.append(f"Visium 坐标：{matched}/{len(spots)} 个 spot 匹配 tissue_positions")
            if meta.get("scale_applied") and abs(float(meta.get("scale") or 1.0) - 1.0) > 1e-9:
                notes.append(f"坐标已按 {meta.get('scale_key')}={float(meta['scale']):.4f} 缩放")
        if image_path:
            image = {"path": image_path,
                     "extent": _image_extent(image_path, raw_coords if raw_coords is not None else np.zeros((1, 2))),
                     "scale": image_scale, "kind": image_kind, "scale_key": image_scale_key}
            notes.append(f"已识别组织图：{os.path.basename(image_path)}（{image_kind}）")
    else:
        warnings.append("没有找到 spatial/ 目录：只能当表达矩阵用（空间分析需要坐标）")
    sample_info = pd.DataFrame(index=pd.Index(spots, name="spot"))
    if spatial_dir:
        try:
            positions_path, _version = _visium_positions_file(spatial_dir)
            head = _read_text(positions_path, limit=1024)
            if "barcode" in head.lower():
                frame_meta = pd.read_csv(positions_path)
                frame_meta.columns = [str(col).strip() for col in frame_meta.columns]
            else:
                frame_meta = pd.read_csv(positions_path, header=None,
                                         names=["barcode", "in_tissue", "array_row", "array_col",
                                                "pxl_row_in_fullres", "pxl_col_in_fullres"])
            frame_meta["barcode"] = frame_meta["barcode"].astype(str)
            frame_meta = frame_meta.set_index("barcode")
            keep_cols = [col for col in ("in_tissue", "array_row", "array_col") if col in frame_meta.columns]
            if keep_cols:
                sample_info = sample_info.join(frame_meta[keep_cols], how="left")
                sample_info.index = pd.Index(spots, name="spot")
        except Exception:  # noqa: BLE001
            pass
    if "in_tissue" in sample_info.columns:
        values = pd.to_numeric(sample_info["in_tissue"], errors="coerce").fillna(0).astype(int)
        sample_info["in_tissue"] = values
        if coords is not None:
            missing = ~np.isfinite(coords).all(axis=1)
            if missing.any():
                keep = ~missing
                warnings.append(f"剔除 {int(missing.sum())} 个无坐标 spot（保持矩阵与坐标一致）")
                dense = dense[keep]
                spots = [spot for spot, ok in zip(spots, keep) if ok]
                coords = coords[keep]
                sample_info = sample_info.loc[spots]
    return {
        "format": "visium",
        "X": dense.astype(np.float32),
        "genes": genes,
        "spots": spots,
        "coords": coords,
        "coords_kind": coords_kind,
        "X_kind": _value_kind(dense),
        "counts": dense.astype(np.float32) if _value_kind(dense) == "counts" else None,
        "sample_info": sample_info,
        "image": image,
        "notes": notes,
        "warnings": warnings,
        "meta": {"format": "visium", "coords_kind": coords_kind,
                 "spatial_dir": spatial_dir, "image": image,
                 "visium_version": visium_version,
                 "scale": position_meta.get("scale"),
                 "scale_key": position_meta.get("scale_key"),
                 "scale_applied": position_meta.get("scale_applied"),
                 "positions_file": position_meta.get("positions_file")},
    }


# ---------------------------------------------------------------------------
# h5ad
# ---------------------------------------------------------------------------

def load_h5ad(path: str, options: dict | None = None) -> dict:
    options = options or {}
    try:
        import anndata
    except Exception as error:  # noqa: BLE001
        raise EngineError("unsupported", f"需要 anndata 读取 .h5ad（{error}）")
    try:
        adata = anndata.read_h5ad(path)
    except Exception as error:  # noqa: BLE001
        raise EngineError("bad_request", f"读取 h5ad 失败：{error}")
    coords = None
    coords_kind = "generic"
    coords_source = None
    try:
        if "spatial" in getattr(adata, "obsm", {}):
            coords = np.asarray(adata.obsm["spatial"], dtype=np.float32)[:, :2]
            coords_kind = "pixel"
            coords_source = "obsm['spatial']"
        else:
            obs = adata.obs
            lowered = {str(col).lower(): str(col) for col in obs.columns}
            pairs = [("pxl_col_in_fullres", "pxl_row_in_fullres"), ("pxl_col", "pxl_row"),
                     ("x", "y"), ("x_coord", "y_coord"), ("imagecol", "imagerow"),
                     ("spatial_1", "spatial_2")]
            for cx, cy in pairs:
                if cx in lowered and cy in lowered:
                    coords = np.column_stack([
                        pd.to_numeric(obs[lowered[cx]], errors="coerce").to_numpy(dtype=np.float32),
                        pd.to_numeric(obs[lowered[cy]], errors="coerce").to_numpy(dtype=np.float32),
                    ])
                    coords_kind = "pixel" if "pxl" in cx else "generic"
                    coords_source = f"obs({lowered[cx]},{lowered[cy]})"
                    break
    except Exception:  # noqa: BLE001
        coords = None
    matrix = adata.X
    try:
        if hasattr(matrix, "toarray"):
            matrix = matrix.toarray()
        dense = np.asarray(matrix, dtype=np.float32)
    except Exception as error:  # noqa: BLE001
        raise EngineError("bad_request", f"h5ad 的 X 不是可读矩阵：{error}")
    genes = [str(item) for item in adata.var_names]
    spots = [str(item) for item in adata.obs_names]
    sample_info = adata.obs.copy()
    sample_info.index = pd.Index(spots, name="spot")
    notes = [f"h5ad：{len(spots)} spot × {len(genes)} 基因"]
    warnings: list[str] = []
    if coords is None:
        warnings.append("h5ad 没有 obsm['spatial'] 或 x/y 列：只能当表达矩阵用")
    else:
        if coords.shape[0] != len(spots):
            warnings.append("坐标行数与细胞数不一致，已忽略坐标")
            coords = None
    image = None
    scalefactors_note = None
    try:
        uns = dict(adata.uns) if hasattr(adata, "uns") else {}
        for key in ("spatial", "image", "images"):
            value = uns.get(key)
            if isinstance(value, dict) and value.get("path"):
                image = {"path": str(value["path"]), "extent": value.get("extent"), "scale": value.get("scale")}
                break
        # Visium 风格 h5ad：uns["spatial"][library_id] 里带 scalefactors（坐标通常是 fullres/hires 像素）
        spatial_uns = uns.get("spatial")
        if isinstance(spatial_uns, dict) and spatial_uns:
            scale_keys = []
            for lib, entry in spatial_uns.items():
                if isinstance(entry, dict) and isinstance(entry.get("scalefactors"), dict):
                    scale_keys = sorted(entry["scalefactors"].keys())
                    break
            if scale_keys:
                scalefactors_note = (f"uns['spatial'] 里带 scalefactors（{', '.join(scale_keys[:6])}）："
                                     "坐标按文件里 obsm['spatial'] 的原始像素使用，未二次缩放")
    except Exception:  # noqa: BLE001
        image = None
    if scalefactors_note:
        notes.append(scalefactors_note)

    # ---- 矩阵标度 & counts 层（P0：X=lognorm + layers["counts"] 时不能丢掉 counts）----
    kind = _value_kind(dense)
    finite_vals = dense[np.isfinite(dense)]
    if kind == "counts":
        matrix_scale = "counts"
    elif kind == "lognorm":
        matrix_scale = "lognorm"
    elif finite_vals.size and float(finite_vals.min()) < 0:
        matrix_scale = "scaled"
    else:
        matrix_scale = "unknown"
    counts = None
    counts_layer = None
    counts_source = None
    layers = getattr(adata, "layers", None) or {}
    for name in ("counts", "count", "raw_counts", "umi_counts", "counts_raw"):
        if name not in layers:
            continue
        try:
            layer = layers[name]
            if hasattr(layer, "toarray"):
                layer = layer.toarray()
            candidate = np.asarray(layer, dtype=np.float32)
        except Exception:  # noqa: BLE001
            continue
        if candidate.shape != dense.shape or candidate.size == 0:
            warnings.append(f"h5ad 的 layers[{name}] 形状 {candidate.shape} 与 X {dense.shape} 不一致，已忽略")
            continue
        counts = candidate
        counts_layer = name
        counts_source = f'layers["{name}"]'
        break
    if counts is None and kind == "counts":
        counts = dense
        counts_source = "X"
    if counts is not None:
        finite_counts = counts[np.isfinite(counts)]
        if finite_counts.size and (finite_counts.min() < 0
                                   or not np.allclose(finite_counts, np.round(finite_counts), atol=1e-6)):
            warnings.append(
                f"{counts_source} 不是非负整数 counts（min={float(finite_counts.min()):.3g}）："
                "RCTD 等需要原始 counts 的方法会在使用前再校验")
        notes.append(f"保留原始 counts：{counts_source}（{counts.shape[0]} × {counts.shape[1]}）")
    elif matrix_scale in ("lognorm", "scaled"):
        warnings.append(
            "h5ad 的 X 已 log/标准化，但没有 layers['counts']（或 raw/counts 层）："
            "RCTD 等需要原始 counts 的方法缺前置；如需 counts 请在导出 h5ad 时保留 counts 层")
    if coords is None:
        coords_source = None
    return {
        "format": "h5ad",
        "X": dense, "genes": genes, "spots": spots, "coords": coords,
        "X_kind": kind, "counts": counts,
        "sample_info": sample_info, "image": image,
        "notes": notes, "warnings": warnings,
        "meta": {"format": "h5ad", "coords_kind": coords_kind,
                 "coords_source": coords_source,
                 "matrix_scale": matrix_scale,
                 "counts_layer": counts_layer,
                 "counts_source": counts_source,
                 "x_kind": kind,
                 "has_scalefactors": bool(scalefactors_note)},
    }


# ---------------------------------------------------------------------------
# Stereo-seq BIN 长表
# ---------------------------------------------------------------------------

def load_stereo(path: str, options: dict | None = None) -> dict:
    options = options or {}
    head = _read_text(path, limit=8192)
    lines = [line for line in head.splitlines() if line.strip()]
    if not lines:
        raise EngineError("bad_request", f"{os.path.basename(path)} 是空文件")
    header = re.split(r"[\t,]", lines[0].lower())
    if len(header) < 4 or not any("gene" in item for item in header):
        raise EngineError("bad_request",
                          "Stereo BIN 表需要 4 列（geneID,x,y,count）并带表头；请用 data_format=matrix_text 并配 coords")
    try:
        frame = pd.read_csv(path, sep="\t" if "\t" in lines[0] else ",")
    except Exception as error:  # noqa: BLE001
        raise EngineError("bad_request", f"解析 Stereo BIN 表失败：{error}")
    cols = {str(col).strip().lower(): str(col) for col in frame.columns}
    gene_col = next((cols[key] for key in cols if "gene" in key), None)
    x_col = next((cols[key] for key in cols if key in ("x", "x_coord", "col")), None)
    y_col = next((cols[key] for key in cols if key in ("y", "y_coord", "row")), None)
    count_col = next((cols[key] for key in cols if "count" in key), None)
    if not all([gene_col, x_col, y_col, count_col]):
        raise EngineError("bad_request", f"Stereo BIN 表列不认识：{list(frame.columns)[:8]}")
    frame["_spot"] = frame[x_col].astype(str) + "x" + frame[y_col].astype(str)
    frame["_gene"] = frame[gene_col].astype(str)
    frame["_count"] = pd.to_numeric(frame[count_col], errors="coerce").fillna(0)
    wide = frame.pivot_table(index="_spot", columns="_gene", values="_count",
                             aggfunc="sum", fill_value=0)
    spots = [str(item) for item in wide.index]
    genes = [str(item) for item in wide.columns]
    dense = wide.to_numpy(dtype=np.float32)
    coords = []
    for spot in spots:
        try:
            x_text, y_text = spot.rsplit("x", 1)
            coords.append([float(x_text), float(y_text)])
        except Exception:  # noqa: BLE001
            coords.append([np.nan, np.nan])
    coords_arr = np.asarray(coords, dtype=np.float32)
    notes = [f"Stereo-seq BIN：{len(spots)} 个 bin × {len(genes)} 基因"]
    return {
        "format": "stereo", "X": dense, "genes": genes, "spots": spots,
        "coords": coords_arr if np.isfinite(coords_arr).all() else None,
        "X_kind": "counts", "counts": dense,
        "sample_info": pd.DataFrame(index=pd.Index(spots, name="spot")),
        "image": None, "notes": notes, "warnings": [],
        "meta": {"format": "stereo", "coords_kind": "grid"},
    }


# ---------------------------------------------------------------------------
# 通用矩阵（+ 坐标）
# ---------------------------------------------------------------------------

def _orient_matrix(frame: pd.DataFrame, options: dict) -> tuple[pd.DataFrame, str]:
    """判定方向，返回 (spots×genes 的数值帧, 说明)。"""
    transpose = options.get("transpose")
    if transpose is not None:
        if bool(transpose):
            return frame.T, "按 transpose=true 转置成 spot×gene"
        return frame, "按 transpose=false 保留原方向（行=spot）"
    n_rows, n_cols = frame.shape
    row_score = _looks_like_genes(frame.index)
    col_score = _looks_like_genes(frame.columns)
    row_spot = _looks_like_spots(frame.index)
    col_spot = _looks_like_spots(frame.columns)
    if n_rows > n_cols * 2.5 and row_score >= col_score:
        return frame.T, "按行名像基因、行数远多于列数判断为 基因×样本，已转置"
    if n_cols > n_rows * 2.5 and col_score >= row_score:
        return frame, "按列名像基因、列数远多于行数判断为 样本×基因，保持原方向"
    if col_spot > row_spot and row_score >= 0.5:
        return frame.T, "按行名像基因、列名像样本判断并转置"
    if row_score >= 0.6 and col_score < 0.5:
        return frame.T, "按行名像基因判断并转置"
    if col_score >= 0.6 and row_score < 0.5:
        return frame, "按列名像基因判断，保持原方向"
    return frame, "方向不明确，按 样本×基因 处理（如不对请指定 transpose）"


def _attach_coords(result: dict, options: dict, workdir: str | None) -> dict:
    """给矩阵补坐标：options['coords'] 指定坐标表。"""
    coords_path = str(options.get("coords") or "").strip()
    if not coords_path:
        result["warnings"].append("没有提供坐标表（coords）：只能当表达矩阵用，空间分析不可用")
        return result
    if not _is_url(coords_path):
        coords_path = os.path.abspath(os.path.expanduser(coords_path))
    labels, coords, meta = read_coords(coords_path, options)
    spots = result["spots"]
    position = {str(label): index for index, label in enumerate(labels)}
    aligned = np.full((len(spots), 2), np.nan, dtype=np.float32)
    matched = 0
    for index, spot in enumerate(spots):
        pos = position.get(str(spot))
        if pos is not None:
            aligned[index] = coords[pos]
            matched += 1
    if matched == 0:
        result["warnings"].append("坐标表与表达矩阵的 spot 名完全不匹配（检查 barcode 列）")
        return result
    if matched < len(spots):
        result["warnings"].append(f"{len(spots) - matched} 个 spot 没有坐标")
    result["coords"] = aligned
    result["meta"]["coords_kind"] = meta.get("coords_kind", "generic")
    result["notes"].append(f"坐标：{matched}/{len(spots)} 个 spot 匹配")
    return result


def load_matrix_text(path: str, options: dict | None = None, workdir: str | None = None) -> dict:
    options = options or {}
    frame = _read_table(path, options)
    numeric, dropped = _numeric_frame(frame)
    if numeric.shape[0] < 2 or numeric.shape[1] < 2:
        raise EngineError("bad_request",
                          f"矩阵太小或不是数值表：{os.path.basename(path)}（{numeric.shape[0]} 行 × {numeric.shape[1]} 列）")
    oriented, note = _orient_matrix(numeric, options)
    dense = np.asarray(oriented.to_numpy(dtype=np.float32))
    kind = _value_kind(dense)
    result = {
        "format": "matrix_text",
        "X": dense,
        "genes": _dedupe([str(item) for item in oriented.columns.tolist()], "GENE"),
        "spots": _dedupe([str(item) for item in oriented.index.tolist()], "spot"),
        "coords": None,
        "X_kind": kind,
        "counts": dense if kind == "counts" else None,
        "sample_info": pd.DataFrame(index=pd.Index(_dedupe([str(i) for i in oriented.index.tolist()], "spot"), name="spot")),
        "image": None,
        "notes": [note, f"矩阵尺寸：{dense.shape[0]} spot × {dense.shape[1]} 基因"],
        "warnings": [f"丢弃 {dropped} 个非数值行/列"] if dropped else [],
        "meta": {"format": "matrix_text", "coords_kind": "generic"},
    }
    return _attach_coords(result, options, workdir)


def load_xlsx(path: str, options: dict | None = None, workdir: str | None = None) -> dict:
    options = options or {}
    try:
        sheets = pd.read_excel(path, sheet_name=None, index_col=0)
    except Exception as error:  # noqa: BLE001
        raise EngineError("bad_request", f"读取 Excel 失败：{os.path.basename(path)}（{error}）")
    wanted = options.get("sheet")
    if wanted and wanted in sheets:
        frame = sheets[wanted]
    else:
        usable = {name: sheet for name, sheet in sheets.items() if sheet.shape[0] > 2 and sheet.shape[1] > 2}
        if not usable:
            raise EngineError("bad_request", f"{os.path.basename(path)} 里没有可用的表达表")
        if len(usable) > 1 and not wanted:
            raise EngineError("ambiguous",
                              f"Excel 有多个候选 sheet：{', '.join(list(usable)[:6])}；请用 sheet 指定",
                              candidates=[{"name": name, "shape": list(sheet.shape)} for name, sheet in list(usable.items())[:12]])
        frame = list(usable.values())[0]
    numeric, dropped = _numeric_frame(frame)
    oriented, note = _orient_matrix(numeric, options)
    dense = np.asarray(oriented.to_numpy(dtype=np.float32))
    kind = _value_kind(dense)
    spots = _dedupe([str(item) for item in oriented.index.tolist()], "spot")
    result = {
        "format": "xlsx", "X": dense,
        "genes": _dedupe([str(item) for item in oriented.columns.tolist()], "GENE"),
        "spots": spots, "coords": None, "X_kind": kind,
        "counts": dense if kind == "counts" else None,
        "sample_info": pd.DataFrame(index=pd.Index(spots, name="spot")),
        "image": None,
        "notes": [note],
        "warnings": [f"丢弃 {dropped} 个非数值行/列"] if dropped else [],
        "meta": {"format": "xlsx", "coords_kind": "generic"},
    }
    return _attach_coords(result, options, workdir)


def load_geo_series(path: str, options: dict | None = None, workdir: str | None = None) -> dict:
    """GEO series matrix：`!series_matrix_table_begin` 之后是 基因×样本 表。"""
    text = _read_text(path)
    lines = text.splitlines()
    start = end = None
    for index, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            start = index + 1
        elif line.startswith("!series_matrix_table_end"):
            end = index
            break
    if start is None:
        raise EngineError("bad_request", "不是 GEO series matrix（找不到 !series_matrix_table_begin）")
    body = "\n".join(lines[start:end if end is not None else len(lines)])
    frame = pd.read_csv(io.StringIO(body), sep="\t", index_col=0)
    frame = frame.apply(pd.to_numeric, errors="coerce")
    frame = frame.loc[~frame.isna().all(axis=1)]
    samples = [str(col).strip('"') for col in frame.columns]
    genes = [str(item).strip('"') for item in frame.index]
    dense = np.asarray(frame.to_numpy(dtype=np.float32).T)  # → 样本 × 基因
    if dense.size == 0 or dense.shape[0] == 0 or dense.shape[1] == 0:
        raise prerequisite(
            f"series matrix 里没有表达数据表（只有 !Series_/!Sample_ 头部）：{os.path.basename(path)}。"
            "这种 GSE 的数据在 suppl/ 的补充文件里（如 *RAW.tar、*counts.csv），请下那些文件再导入")
    chars: dict[str, list[str]] = {}
    for line in lines[:start - 1]:
        if line.startswith("!Sample_characteristics_ch1"):
            parts = [part.strip().strip('"') for part in line.split("\t")[1:]]
            key = parts[0].split(":")[0].strip() if parts and ":" in parts[0] else f"char{len(chars) + 1}"
            values = [part.split(":", 1)[1].strip() if ":" in part else part for part in parts]
            chars[key or f"char{len(chars) + 1}"] = values
        elif line.startswith("!Sample_title"):
            chars["title"] = [part.strip().strip('"') for part in line.split("\t")[1:]]
    sample_info = pd.DataFrame(index=pd.Index(samples, name="spot"))
    for key, values in chars.items():
        if len(values) == len(samples):
            sample_info[key] = values
    kind = _value_kind(dense)
    result = {
        "format": "geo_series", "X": dense, "genes": _dedupe(genes, "GENE"),
        "spots": samples, "coords": None, "X_kind": kind,
        "counts": dense if kind == "counts" else None,
        "sample_info": sample_info, "image": None,
        "notes": [f"GEO 系列矩阵：{len(samples)} 样本 × {len(genes)} 基因"],
        "warnings": ["GEO series matrix 通常没有空间坐标：空间分析需要另给坐标表（coords）"],
        "meta": {"format": "geo_series", "coords_kind": "generic"},
    }
    return _attach_coords(result, options, workdir)


# ---------------------------------------------------------------------------
# 压缩包
# ---------------------------------------------------------------------------

def _safe_extract(archive_path: str, dest: str, limit_bytes: int = 40 * 1024 ** 3) -> None:
    ensure_dir(dest)
    base = os.path.realpath(dest)
    if zipfile.is_zipfile(archive_path):
        with zipfile.ZipFile(archive_path) as handle:
            total = 0
            for info in handle.infolist():
                total += int(info.file_size)
                if total > limit_bytes:
                    raise EngineError("unsupported", "压缩包解压后超过 40GB 上限")
                target = os.path.realpath(os.path.join(dest, info.filename))
                if not target.startswith(base + os.sep) and target != base:
                    raise EngineError("bad_request", f"压缩包含越权路径：{info.filename}")
                if info.is_dir():
                    continue
                os.makedirs(os.path.dirname(target), exist_ok=True)
                with handle.open(info) as src, open(target, "wb") as out:
                    shutil.copyfileobj(src, out)
        return
    try:
        with tarfile.open(archive_path) as handle:
            total = 0
            for member in handle.getmembers():
                if member.issym() or member.islnk():
                    continue
                total += int(member.size)
                if total > limit_bytes:
                    raise EngineError("unsupported", "压缩包解压后超过 40GB 上限")
                target = os.path.realpath(os.path.join(dest, member.name))
                if not target.startswith(base + os.sep) and target != base:
                    raise EngineError("bad_request", f"压缩包含越权路径：{member.name}")
            handle.extractall(dest)
    except tarfile.TarError as error:
        raise EngineError("bad_request", f"不认识的压缩包格式：{error}")


_MATRIX_SUFFIX = (".tsv", ".csv", ".txt", ".tsv.gz", ".csv.gz", ".txt.gz", ".xlsx", ".xls", ".h5ad", ".h5")


def _find_candidates(root: str) -> list[dict]:
    """在解包目录里找候选数据源（Visium 目录 / BIN 表 / 矩阵 / h5ad）。"""
    candidates: list[dict] = []
    for current, dirs, files in os.walk(root):
        dirs[:] = [name for name in dirs if not name.startswith(".")]
        lowered = {name.lower() for name in files}
        subdirs = {name.lower() for name in dirs}
        has_spatial_dir = "spatial" in subdirs or "spatial" in lowered
        has_matrix = bool(set(lowered) & set(VISIUM_MATRIX_H5)) or "filtered_feature_bc_matrix" in subdirs \
            or "filtered_feature_bc_matrix" in lowered
        if has_spatial_dir and has_matrix:
            candidates.append({"path": current, "format": "visium", "name": os.path.basename(current) or "visium"})
            dirs[:] = []
            continue
        for name in files:
            lower = name.lower()
            full = os.path.join(current, name)
            if "bin_count" in lower or "bin1" in lower and "raw_bin" in lower:
                candidates.append({"path": full, "format": "stereo", "name": name})
            elif lower.endswith(".h5ad"):
                candidates.append({"path": full, "format": "h5ad", "name": name})
            elif lower.endswith(_MATRIX_SUFFIX):
                candidates.append({"path": full, "format": "matrix_text", "name": name})
    return candidates


def load_archive(path: str, options: dict | None = None, workdir: str | None = None) -> dict:
    options = options or {}
    base = workdir or tempfile.gettempdir()
    dest = os.path.join(ensure_dir(os.path.join(base, "unpack")), f"arc-{uuid.uuid4().hex[:10]}")
    try:
        _safe_extract(path, dest)
    except EngineError:
        shutil.rmtree(dest, ignore_errors=True)
        raise
    candidates = _find_candidates(dest)
    member = str(options.get("member") or "").strip()
    if member:
        chosen = os.path.join(dest, member)
        if not os.path.exists(chosen):
            raise EngineError("not_found", f"压缩包里没有 {member}", candidates=candidates[:24])
        fmt = str(options.get("data_format") or "auto")
        if fmt in ("auto", "archive"):
            fmt = "visium" if os.path.isdir(chosen) else ("h5ad" if chosen.endswith(".h5ad") else "matrix_text")
        return _dispatch_format(chosen, fmt, options, workdir, extract_dir=dest)
    if not candidates:
        raise EngineError("not_found", "压缩包里没有找到 Visium 目录 / BIN 表 / 表达矩阵 / h5ad")
    if len(candidates) > 1:
        visium = [item for item in candidates if item["format"] == "visium"]
        if len(visium) == 1:
            candidates = visium
        else:
            raise EngineError("ambiguous",
                              f"压缩包里有 {len(candidates)} 个候选数据源，请用 member 指定",
                              candidates=candidates[:24])
    chosen = candidates[0]
    result = _dispatch_format(chosen["path"], chosen["format"], options, workdir, extract_dir=dest)
    result["notes"].append(f"压缩包已解到 {dest}（可用 member 精确指定）")
    return result


# ---------------------------------------------------------------------------
# R 桥
# ---------------------------------------------------------------------------

_R_READER = r"""
args <- commandArgs(trailingOnly = TRUE)
in_path <- args[1]; out_dir <- args[2]
obj_name <- if (length(args) >= 3 && nzchar(args[3])) args[3] else ""
read_obj <- function(p, name) {
  if (grepl("\\.rds$", p, ignore.case=TRUE)) return(readRDS(p))
  env <- new.env(); load(p, envir=env)
  if (nzchar(name)) return(get(name, envir=env))
  objs <- ls(env); if (!length(objs)) stop("RData 里没有对象")
  get(objs[1], envir=env)
}
obj <- read_obj(in_path, obj_name)
get_mtx <- function(o) {
  if (inherits(o, "Seurat")) {
    m <- tryCatch(SeuratObject::GetAssayData(o, layer="counts"), error=function(e) NULL)
    if (is.null(m)) m <- tryCatch(o@assays[[1]]@counts, error=function(e) NULL)
    if (is.null(m)) m <- get_mtx(o@assays[[1]]@data)
    return(as.matrix(m))
  }
  if (inherits(o, "SingleCellExperiment")) return(as.matrix(SummarizedExperiment::assay(o)))
  if (inherits(o, "ExpressionSet")) return(exprs(o))
  if (is.matrix(o)) return(o)
  if (is.data.frame(o)) return(as.matrix(o))
  m <- tryCatch(as.matrix(o), error=function(e) NULL)
  if (is.null(m)) stop("无法把对象转成矩阵")
  m
}
m <- get_mtx(obj)
write.csv(m, file.path(out_dir, "matrix.csv"), quote=FALSE)
if (inherits(obj, "Seurat")) {
  try(write.csv(obj@meta.data, file.path(out_dir, "obs.csv"), quote=FALSE), silent=TRUE)
  emb <- tryCatch(SeuratObject::Embeddings(obj, "spatial"), error=function(e) NULL)
  if (!is.null(emb)) write.csv(emb, file.path(out_dir, "spatial.csv"), quote=FALSE)
}
cat("OK\n")
"""


def load_r_object(path: str, options: dict | None = None, workdir: str | None = None) -> dict:
    options = options or {}
    rscript = None
    for candidate in (
        "/Library/Frameworks/R.framework/Resources/bin/Rscript",
        "/opt/homebrew/bin/Rscript", "/usr/local/bin/Rscript", "Rscript",
    ):
        if os.path.sep not in candidate or os.path.exists(candidate):
            rscript = candidate
            break
    if not rscript:
        raise EngineError("unsupported", "读取 rds/RData 需要 Rscript（未找到）；可先在 R 里导出成 csv/tsv")
    import subprocess

    base = workdir or tempfile.gettempdir()
    outdir = os.path.join(ensure_dir(base), "r-bridge", uuid.uuid4().hex[:10])
    ensure_dir(outdir)
    script_path = os.path.join(outdir, "read.R")
    with open(script_path, "w", encoding="utf-8") as handle:
        handle.write(_R_READER)
    try:
        proc = subprocess.run([rscript, "--vanilla", script_path, os.path.abspath(path), outdir,
                               str(options.get("object") or "")],
                              capture_output=True, text=True, timeout=600)
    except Exception as error:  # noqa: BLE001
        raise EngineError("internal", f"调用 Rscript 失败：{error}")
    matrix_path = os.path.join(outdir, "matrix.csv")
    if proc.returncode != 0 or not os.path.exists(matrix_path):
        detail = (proc.stderr or proc.stdout or "").strip()[-800:]
        raise EngineError("unsupported", f"R 读取失败（可能需要 SeuratObject/SummarizedExperiment）：{detail}")
    frame = pd.read_csv(matrix_path, index_col=0)
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.loc[~numeric.isna().all(axis=1)]
    oriented, note = _orient_matrix(numeric, options)
    dense = np.asarray(oriented.to_numpy(dtype=np.float32))
    kind = _value_kind(dense)
    spots = _dedupe([str(item) for item in oriented.index.tolist()], "spot")
    sample_info = pd.DataFrame(index=pd.Index(spots, name="spot"))
    obs_path = os.path.join(outdir, "obs.csv")
    if os.path.exists(obs_path):
        try:
            obs = pd.read_csv(obs_path, index_col=0)
            obs.index = [str(item) for item in obs.index]
            sample_info = obs.reindex(spots)
            sample_info.index = pd.Index(spots, name="spot")
        except Exception:  # noqa: BLE001
            pass
    coords = None
    spatial_path = os.path.join(outdir, "spatial.csv")
    if os.path.exists(spatial_path):
        try:
            emb = pd.read_csv(spatial_path, index_col=0)
            emb.index = [str(item) for item in emb.index]
            values = emb.reindex(spots)[[emb.columns[0], emb.columns[1]]].to_numpy(dtype=np.float32)
            if np.isfinite(values).any():
                coords = values
        except Exception:  # noqa: BLE001
            coords = None
    result = {
        "format": "rds" if path.lower().endswith(".rds") else "rdata",
        "X": dense,
        "genes": _dedupe([str(item) for item in oriented.columns.tolist()], "GENE"),
        "spots": spots, "coords": coords, "X_kind": kind,
        "counts": dense if kind == "counts" else None,
        "sample_info": sample_info, "image": None,
        "notes": [note, "经 Rscript 读取（R 桥）"],
        "warnings": [] if coords is not None else ["R 对象里没有空间坐标（仅有表达矩阵）"],
        "meta": {"format": "rds", "coords_kind": "pixel" if coords is not None else "generic"},
    }
    return result


# ---------------------------------------------------------------------------
# 格式识别与分发
# ---------------------------------------------------------------------------

def detect_format(path: str) -> str:
    if _is_url(path):
        lower = path.lower().split("?")[0]
        for suffix, fmt in ((".h5ad", "h5ad"), (".h5", "visium_h5"), (".rds", "rds"),
                            (".rdata", "rdata"), (".rda", "rdata"), (".xlsx", "xlsx"),
                            (".xls", "xlsx"), (".zip", "archive"), (".tar.gz", "archive"),
                            (".tgz", "archive"), (".tar", "archive")):
            if lower.endswith(suffix):
                return fmt
        return "matrix_text"
    if os.path.isdir(path):
        try:
            _find_visium_matrix(path)
            return "visium"
        except EngineError:
            pass
        raise EngineError("bad_request",
                          f"目录 {os.path.basename(path)} 不是 Visium outs/（缺 filtered_feature_bc_matrix）；"
                          "如果是 mtx 三件套目录请改用 zip 或指定 data_format=visium")
    lower = path.lower()
    if lower.endswith(".h5ad"):
        return "h5ad"
    if lower.endswith((".h5", ".hdf5")):
        return "visium_h5"
    if lower.endswith((".rds",)):
        return "rds"
    if lower.endswith((".rdata", ".rda")):
        return "rdata"
    if lower.endswith((".xlsx", ".xls")):
        return "xlsx"
    if lower.endswith((".zip", ".tar", ".tar.gz", ".tgz", ".tar.bz2")):
        return "archive"
    head = _read_text(path, limit=64 * 1024)
    if "!series_matrix_table_begin" in head:
        return "geo_series"
    first = next((line for line in head.splitlines() if line.strip()), "")
    header = re.split(r"[\t,]", first.lower())
    if len(header) >= 4 and any("gene" in item for item in header) and "x" in header and "y" in header:
        return "stereo"
    return "matrix_text"


def _dispatch_format(path: str, data_format: str, options: dict, workdir: str | None,
                     extract_dir: str | None = None) -> dict:
    fmt = str(data_format or "auto").strip().lower()
    if fmt == "auto":
        fmt = detect_format(path)
    if fmt == "visium":
        return load_visium(path, options)
    if fmt == "visium_h5":
        dense, genes, spots = _read_10x_h5(path)
        kind = _value_kind(dense)
        result = {
            "format": "visium_h5", "X": dense, "genes": genes, "spots": spots, "coords": None,
            "X_kind": kind, "counts": dense if kind == "counts" else None,
            "sample_info": pd.DataFrame(index=pd.Index(spots, name="spot")),
            "image": None, "notes": ["10x .h5 表达矩阵"], "warnings": [],
            "meta": {"format": "visium_h5", "coords_kind": "generic"},
        }
        return _attach_coords(result, options, workdir)
    if fmt == "h5ad":
        return load_h5ad(path, options)
    if fmt == "stereo":
        return load_stereo(path, options)
    if fmt == "slideseq":
        return load_matrix_text(path, options, workdir)
    if fmt == "matrix_text":
        return load_matrix_text(path, options, workdir)
    if fmt == "xlsx":
        return load_xlsx(path, options, workdir)
    if fmt == "geo_series":
        return load_geo_series(path, options, workdir)
    if fmt == "archive":
        return load_archive(path, options, workdir)
    if fmt in ("rds", "rdata"):
        return load_r_object(path, options, workdir)
    raise bad_request(f"未知的 data_format：{data_format}（可用：{', '.join(SUPPORTED_FORMATS)}）")


def _download_url(url: str, workdir: str | None) -> str:
    """下载远端的空间数据到 workdir/downloads 并返回本地路径。"""
    import ssl
    import urllib.parse
    import urllib.request

    base = workdir or tempfile.gettempdir()
    dest_dir = ensure_dir(os.path.join(base, "downloads"))
    name = os.path.basename(urllib.parse.unquote(urllib.parse.urlparse(url).path)) or f"spatial-{uuid.uuid4().hex[:8]}"
    target = os.path.join(dest_dir, f"{uuid.uuid4().hex[:6]}-{name}")
    try:
        import certifi

        context = ssl.create_default_context(cafile=certifi.where())
    except Exception:  # noqa: BLE001
        context = ssl.create_default_context()
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "dsh-spatial/1.0"}),
                                    timeout=300, context=context) as response, open(target, "wb") as handle:
            shutil.copyfileobj(response, handle)
    except Exception as error:  # noqa: BLE001
        raise EngineError("unsupported", f"下载失败：{error}")
    return target


def load_source(path: str, workdir: str | None = None, options: dict | None = None) -> dict:
    """主入口：路径/URL → 统一的结果 dict（DESIGN §6）。"""
    options = dict(options or {})
    data_format = str(options.get("data_format") or options.get("format") or "auto").strip().lower()
    raw = str(path or "").strip()
    if not raw:
        raise bad_request("请提供 path（本地文件/目录）或 url")
    if _is_url(raw):
        raw = _download_url(raw, workdir)
        options["name"] = options.get("name") or os.path.basename(raw)
    else:
        raw = os.path.abspath(os.path.expanduser(raw))
        if not os.path.exists(raw):
            raise EngineError("not_found", f"路径不存在：{raw}")
    result = _dispatch_format(raw, data_format, options, workdir)
    result.setdefault("notes", [])
    result.setdefault("warnings", [])
    result.setdefault("meta", {})
    result["meta"].setdefault("format", result.get("format"))
    X = np.asarray(result.get("X"))
    if X.ndim != 2:
        raise EngineError("internal", "导入结果不是二维矩阵")
    if result.get("coords") is not None:
        coords = np.asarray(result["coords"], dtype=np.float32)
        if coords.shape != (X.shape[0], 2):
            result["coords"] = None
            result["warnings"].append("坐标形状与矩阵不一致，已忽略坐标")
    return result


# ---------------------------------------------------------------------------
# inspect：只探测，不加载整表
# ---------------------------------------------------------------------------

def inspect_source(path: str, workdir: str | None = None, options: dict | None = None) -> dict:
    options = dict(options or {})
    data_format = str(options.get("data_format") or options.get("format") or "auto").strip().lower()
    raw = str(path or "").strip()
    if not raw:
        raise bad_request("请提供 path 或 url")
    remote = _is_url(raw)
    if not remote:
        raw = os.path.abspath(os.path.expanduser(raw))
        if not os.path.exists(raw):
            raise EngineError("not_found", f"路径不存在：{raw}")
    fmt = data_format if data_format != "auto" else ("visium" if (not remote and os.path.isdir(raw)) else "auto")
    if remote and data_format == "auto":
        return {"format": detect_format(raw), "size": None,
                "notes": ["远端文件：先用 open 下载后再探测更准确"]}
    if fmt == "visium":
        matrix_kind, matrix_path = _find_visium_matrix(raw)
        spatial_dir = _find_spatial_dir(raw)
        n_spots = n_genes = None
        try:
            if matrix_kind == "h5":
                dense, genes, spots = _read_10x_h5(matrix_path)
                n_spots, n_genes = dense.shape[0], dense.shape[1]
            else:
                n_spots, n_genes = None, None
        except Exception:  # noqa: BLE001
            pass
        has_image = bool(spatial_dir and _find_tissue_image(spatial_dir)[0])
        return {
            "format": "visium", "n_spots": n_spots, "n_genes": n_genes,
            "geometry": {"mode": "visium", "spatial_dir": spatial_dir, "has_image": has_image},
            "notes": [f"矩阵：{os.path.basename(matrix_path)}", "spatial/ 目录：" + ("有" if spatial_dir else "无")],
            "candidates": [],
            "size": _dir_size(raw) if os.path.isdir(raw) else os.path.getsize(raw),
        }
    if fmt == "auto":
        fmt = detect_format(raw)
    size = os.path.getsize(raw) if os.path.isfile(raw) else _dir_size(raw)
    info: dict = {"format": fmt, "size": size, "notes": [], "candidates": []}
    try:
        if fmt in ("matrix_text", "stereo"):
            head = _read_text(raw, limit=64 * 1024)
            lines = head.splitlines()
            first = next((line for line in lines if line.strip()), "")
            sep = "\t" if "\t" in first else ("," if "," in first else None)
            n_cols = len(re.split(r"[\t,]", first)) if sep else len(first.split())
            info["notes"].append(f"表头列数约 {n_cols}；显示前 {len(lines)} 行样本")
            if fmt == "matrix_text":
                frame = pd.read_csv(raw, sep=sep or r"\s+", index_col=0, nrows=200,
                                    engine="python", on_bad_lines="skip") if sep else None
                if frame is not None:
                    info["n_samples"] = int(frame.shape[0])
                    info["n_features"] = int(frame.shape[1])
                    info["notes"].append("行数按前 200 行估计；实际以 open 为准")
        elif fmt == "h5ad":
            try:
                import anndata

                adata = anndata.read_h5ad(raw, backed="r")
                info["n_spots"] = int(adata.n_obs)
                info["n_genes"] = int(adata.n_vars)
                info["has_coords"] = "spatial" in getattr(adata, "obsm", {})
                info["notes"].append(f"h5ad：{adata.n_obs} × {adata.n_vars}；"
                                     f"spatial 坐标{'有' if info['has_coords'] else '无'}")
                try:
                    adata.file.close()
                except Exception:  # noqa: BLE001
                    pass
            except Exception as error:  # noqa: BLE001
                info["notes"].append(f"读取 h5ad 头失败：{error}")
        elif fmt == "xlsx":
            try:
                sheets = pd.read_excel(raw, sheet_name=None, nrows=5)
                info["candidates"] = [{"name": name, "shape_note": "前 5 行"} for name in list(sheets)[:12]]
                if sheets:
                    first_sheet = list(sheets.values())[0]
                    info["notes"].append(f"sheet 数 {len(sheets)}；首个 sheet 列数 {first_sheet.shape[1]}")
            except Exception as error:  # noqa: BLE001
                info["notes"].append(f"读取 Excel 头失败：{error}")
        elif fmt == "archive":
            names = []
            if zipfile.is_zipfile(raw):
                with zipfile.ZipFile(raw) as handle:
                    names = [item.filename for item in handle.infolist()][:200]
            else:
                try:
                    with tarfile.open(raw) as handle:
                        names = [member.name for member in handle.getmembers()][:200]
                except tarfile.TarError:
                    names = []
            info["candidates"] = [{"name": name} for name in names[:30]]
            info["notes"].append(f"压缩包内 {len(names)} 个条目（显示前 30）")
        elif fmt == "geo_series":
            text = _read_text(raw, limit=256 * 1024)
            info["notes"].append("GEO series matrix（表段为 基因×样本）")
            info["has_coords"] = False
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        info["notes"].append(f"探测时遇到问题（不影响 open）：{error}")
    return info


def _dir_size(path: str) -> int:
    total = 0
    for current, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += os.path.getsize(os.path.join(current, name))
            except OSError:
                continue
    return total


__all__ = ["SUPPORTED_FORMATS", "FORMAT_INFO", "load_source", "inspect_source",
           "detect_format", "read_coords", "load_visium", "load_h5ad", "load_stereo",
           "load_matrix_text", "load_xlsx", "load_geo_series", "load_archive", "load_r_object"]
