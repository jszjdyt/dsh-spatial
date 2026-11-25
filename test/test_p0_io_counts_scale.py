#!/usr/bin/env python3
"""P0 回归：H5AD 的 counts 层 / 矩阵标度元数据 与 Visium scalefactor / 输出版本。

修复前反例（probe_before.py 实测）：
  * load_h5ad（X=log1p(counts) + layers["counts"]）：counts=None、X_kind=lognorm、
    meta 里没有 matrix_scale / counts_layer（RCTD 因此拿不到参考原始 counts）
  * Visium 只有低分辨率图（tissue_lowres_scalef=1.0、tissue_hires_scalef=5.0）：
    坐标 [10,10]/[20,10] 被错误乘上 hires 5.0 → [50,50]/[100,50]

修复后要求：
  * h5ad 保留 layers["counts"]，X 仍是 lognorm，meta 记录 matrix_scale/counts_layer/坐标系
  * Visium 按「实际使用的那张图」取 scalefactor：低分辨率图 ×lowres，绝不 ×hires
  * Visium 输出版本（v1 tissue_positions_list.csv / v2 tissue_positions.csv）被识别与校验

运行：MPLBACKEND=Agg MPLCONFIGDIR=<tmp> python3 test/test_p0_io_counts_scale.py
"""

from __future__ import annotations

import gzip
import json
import os
import sys
import tempfile
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(os.path.dirname(HERE), "lib")
if LIB not in sys.path:
    sys.path.insert(0, LIB)

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="sp-mpl-"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import sp_common as sc  # noqa: E402
import sp_io  # noqa: E402

WORK = tempfile.mkdtemp(prefix="sp-p0-io-")


# ---------------------------------------------------------------------------
# h5ad fixtures
# ---------------------------------------------------------------------------

def make_h5ad(name: str, X: np.ndarray, counts_layer=None, coords=None) -> str:
    import anndata as ad

    obs = pd.DataFrame(index=[f"c{i}" for i in range(X.shape[0])])
    adata = ad.AnnData(X=np.asarray(X, dtype=np.float32), obs=obs,
                       var=pd.DataFrame(index=[f"G{j + 1}" for j in range(X.shape[1])]))
    if counts_layer is not None:
        adata.layers["counts"] = np.asarray(counts_layer, dtype=np.int32)
    if coords is not None:
        adata.obsm["spatial"] = np.asarray(coords, dtype=np.float32)
    path = os.path.join(WORK, name)
    adata.write_h5ad(path)
    return path


def check_h5ad_counts_layer() -> list[str]:
    failures: list[str] = []
    raw = np.array([[3, 0, 7], [0, 4, 2], [5, 1, 0], [2, 8, 3]], dtype=np.int32)
    path = make_h5ad("lognorm_counts.h5ad", np.log1p(raw.astype(np.float32)), counts_layer=raw)
    loaded = sp_io.load_h5ad(path, {})
    meta = loaded.get("meta") or {}
    print(f"  .. X_kind={loaded['X_kind']} counts is None={loaded['counts'] is None}")
    print(f"  .. meta={ {k: meta.get(k) for k in ('matrix_scale', 'counts_layer', 'coords_kind', 'coords_source')} }")
    if loaded["counts"] is None:
        failures.append("layers['counts'] 没有被保留（修复前 counts=None，RCTD 拿不到原始 counts）")
    else:
        if not np.array_equal(np.asarray(loaded["counts"], dtype=np.int32), raw):
            failures.append("counts 层数值与写入的原始 counts 不一致")
    if loaded["X_kind"] != "lognorm":
        failures.append(f"X=log1p(counts) 时 X_kind 应为 lognorm，实际 {loaded['X_kind']!r}")
    if not np.allclose(np.asarray(loaded["X"]), np.log1p(raw), atol=1e-6):
        failures.append("X 被改写（应保留原始 lognorm 矩阵，不得换成 counts）")
    if meta.get("matrix_scale") != "lognorm":
        failures.append(f"meta.matrix_scale 应为 lognorm，实际 {meta.get('matrix_scale')!r}")
    if str(meta.get("counts_layer") or "") != "counts":
        failures.append(f"meta.counts_layer 应为 counts，实际 {meta.get('counts_layer')!r}")

    # X 本身是 counts、没有 counts 层
    path2 = make_h5ad("counts_only.h5ad", raw)
    loaded2 = sp_io.load_h5ad(path2, {})
    meta2 = loaded2.get("meta") or {}
    print(f"  .. 纯 counts：X_kind={loaded2['X_kind']} meta.matrix_scale={meta2.get('matrix_scale')}")
    if loaded2["X_kind"] != "counts":
        failures.append(f"整数矩阵 X_kind 应为 counts，实际 {loaded2['X_kind']!r}")
    if loaded2["counts"] is None:
        failures.append("X 本身是 counts 时也应填 counts")
    if meta2.get("matrix_scale") != "counts":
        failures.append(f"纯 counts 的 meta.matrix_scale 应为 counts，实际 {meta2.get('matrix_scale')!r}")

    # X 是 scaled（含负值）+ counts 层：matrix_scale 记录 scaled，counts 仍保留
    scaled = (raw - raw.mean()) / (raw.std() + 1e-9)
    path3 = make_h5ad("scaled_counts.h5ad", scaled, counts_layer=raw)
    loaded3 = sp_io.load_h5ad(path3, {})
    meta3 = loaded3.get("meta") or {}
    print(f"  .. scaled+counts：X_kind={loaded3['X_kind']} matrix_scale={meta3.get('matrix_scale')} "
          f"counts_layer={meta3.get('counts_layer')}")
    if loaded3["counts"] is None:
        failures.append("scaled X + counts 层时 counts 丢失")
    if meta3.get("matrix_scale") != "scaled":
        failures.append(f"含负值 X 的 matrix_scale 应为 scaled，实际 {meta3.get('matrix_scale')!r}")

    # 带 obsm['spatial'] 时记录坐标系来源
    path4 = make_h5ad("spatial.h5ad", raw, coords=np.array([[1, 2], [3, 4], [5, 6], [7, 8]]))
    loaded4 = sp_io.load_h5ad(path4, {})
    meta4 = loaded4.get("meta") or {}
    print(f"  .. obsm spatial：coords_kind={meta4.get('coords_kind')} 来源={meta4.get('coords_source')}")
    if meta4.get("coords_source") != "obsm['spatial']":
        failures.append(f"坐标系来源应记录 obsm['spatial']，实际 {meta4.get('coords_source')!r}")
    return failures


# ---------------------------------------------------------------------------
# Visium fixtures
# ---------------------------------------------------------------------------

def make_visium(name: str, *, hires: bool = True, lowres: bool = True,
                positions: str = "v2", scalefactors: dict | None = None,
                invalid_positions: bool = False) -> str:
    """造一个最小 Visium outs/ 目录（mtx 三件套 + spatial/）。"""
    root = os.path.join(WORK, name)
    mtx_dir = os.path.join(root, "filtered_feature_bc_matrix")
    spatial = os.path.join(root, "spatial")
    os.makedirs(mtx_dir, exist_ok=True)
    os.makedirs(spatial, exist_ok=True)
    from scipy import sparse
    from scipy.io import mmwrite

    counts = np.array([[3, 0], [0, 4], [5, 1]], dtype=np.int32)  # genes × spots
    plain_mtx = os.path.join(mtx_dir, "matrix.mtx")
    mmwrite(plain_mtx, sparse.csr_matrix(counts))
    with open(plain_mtx, "rb") as src, gzip.open(plain_mtx + ".gz", "wb") as dst:
        dst.write(src.read())
    os.remove(plain_mtx)
    with gzip.open(os.path.join(mtx_dir, "barcodes.tsv.gz"), "wt") as handle:
        handle.write("b1\nb2\n")
    with gzip.open(os.path.join(mtx_dir, "features.tsv.gz"), "wt") as handle:
        handle.write("G1\tG1\tGene Expression\nG2\tG2\tGene Expression\nG3\tG3\tGene Expression\n")
    scale = scalefactors if scalefactors is not None else {
        "tissue_hires_scalef": 5.0, "tissue_lowres_scalef": 1.0, "spot_diameter_fullres": 10.0}
    with open(os.path.join(spatial, "scalefactors_json.json"), "w", encoding="utf-8") as handle:
        json.dump(scale, handle)
    if invalid_positions:
        frame = pd.DataFrame({"barcode": ["b1", "b2"], "in_tissue": [1, 1],
                              "array_row": [0, 0], "array_col": [0, 1]})
    else:
        frame = pd.DataFrame({"barcode": ["b1", "b2"], "in_tissue": [1, 1],
                              "array_row": [0, 0], "array_col": [0, 1],
                              "pxl_row_in_fullres": [10, 10], "pxl_col_in_fullres": [10, 20]})
    if positions == "v1":
        frame.to_csv(os.path.join(spatial, "tissue_positions_list.csv"), index=False, header=False)
    else:
        frame.to_csv(os.path.join(spatial, "tissue_positions.csv"), index=False)
    from PIL import Image

    if hires:
        Image.fromarray(np.zeros((50, 50, 3), dtype=np.uint8)).save(
            os.path.join(spatial, "tissue_hires_image.png"))
    if lowres:
        Image.fromarray(np.zeros((10, 10, 3), dtype=np.uint8)).save(
            os.path.join(spatial, "tissue_lowres_image.png"))
    return root


def _visium_info(loaded: dict) -> dict:
    meta = loaded.get("meta") or {}
    return {"coords": loaded.get("coords"), "coords_kind": loaded.get("coords_kind"),
            "image": loaded.get("image"), "meta": meta,
            "notes": loaded.get("notes"), "warnings": loaded.get("warnings")}


def check_visium_lowres_only() -> list[str]:
    failures: list[str] = []
    root = make_visium("visium_lowres", hires=False, lowres=True)
    loaded = sp_io.load_visium(root, {})
    info = _visium_info(loaded)
    coords = np.asarray(info["coords"], dtype=float)
    print(f"  .. coords={coords.tolist()} coords_kind={info['coords_kind']}")
    print(f"  .. image={ {k: info['image'].get(k) for k in ('kind', 'scale')} }")
    print(f"  .. meta.scale={info['meta'].get('scale')} 用于={info['meta'].get('scale_key')} "
          f"version={info['meta'].get('visium_version')}")
    expected = np.array([[10.0, 10.0], [20.0, 10.0]])
    if coords.shape != expected.shape or not np.allclose(coords, expected):
        failures.append(f"只有低分辨率图时坐标应为 {expected.tolist()}（×lowres 1.0），"
                        f"实际 {coords.tolist()}（修复前 ×hires 5.0 → [[50,50],[100,50]]）")
    if info["coords_kind"] != "visium_lowres":
        failures.append(f"coords_kind 应为 visium_lowres，实际 {info['coords_kind']!r}")
    image = info["image"] or {}
    if image.get("kind") != "lowres":
        failures.append(f"image.kind 应为 lowres，实际 {image.get('kind')!r}")
    if image.get("scale") is None or abs(float(image["scale"]) - 1.0) > 1e-9:
        failures.append(f"image.scale 应为 lowres=1.0，实际 {image.get('scale')!r}")
    if info["meta"].get("scale_key") != "tissue_lowres_scalef":
        failures.append(f"meta.scale_key 应为 tissue_lowres_scalef，实际 {info['meta'].get('scale_key')!r}")
    return failures


def check_visium_hires_preferred() -> list[str]:
    failures: list[str] = []
    # 两张图都在：按契约用 hires
    root = make_visium("visium_both", hires=True, lowres=True)
    loaded = sp_io.load_visium(root, {})
    info = _visium_info(loaded)
    coords = np.asarray(info["coords"], dtype=float)
    print(f"  .. 两图都在 coords={coords.tolist()} kind={info['coords_kind']} "
          f"scale_key={info['meta'].get('scale_key')}")
    if not np.allclose(coords, np.array([[50.0, 50.0], [100.0, 50.0]])):
        failures.append(f"两图都在时应按契约用 hires（×5）→ [[50,50],[100,50]]，实际 {coords.tolist()}")
    if info["meta"].get("scale_key") != "tissue_hires_scalef":
        failures.append(f"两图都在时 scale_key 应为 tissue_hires_scalef，实际 {info['meta'].get('scale_key')}")
    # 只有 hires 图
    root = make_visium("visium_hires", hires=True, lowres=False)
    loaded = sp_io.load_visium(root, {})
    info = _visium_info(loaded)
    coords = np.asarray(info["coords"], dtype=float)
    print(f"  .. 仅 hires coords={coords.tolist()} kind={info['coords_kind']}")
    if not np.allclose(coords, np.array([[50.0, 50.0], [100.0, 50.0]])):
        failures.append(f"仅 hires 图时坐标应 ×5，实际 {coords.tolist()}")
    return failures


def check_visium_version_validation() -> list[str]:
    failures: list[str] = []
    # v1：tissue_positions_list.csv（无表头）
    root = make_visium("visium_v1", hires=True, lowres=False, positions="v1")
    loaded = sp_io.load_visium(root, {})
    info = _visium_info(loaded)
    print(f"  .. v1 version={info['meta'].get('visium_version')} coords={np.asarray(info['coords']).tolist()}")
    if info["meta"].get("visium_version") != "v1":
        failures.append(f"v1 版本识别错误：{info['meta'].get('visium_version')!r}")
    if not np.allclose(np.asarray(info["coords"], dtype=float), np.array([[50.0, 50.0], [100.0, 50.0]])):
        failures.append("v1 坐标缩放错误")
    # v2：tissue_positions.csv（有表头）
    root = make_visium("visium_v2", hires=True, lowres=False, positions="v2")
    loaded = sp_io.load_visium(root, {})
    info = _visium_info(loaded)
    print(f"  .. v2 version={info['meta'].get('visium_version')}")
    if info["meta"].get("visium_version") != "v2":
        failures.append(f"v2 版本识别错误：{info['meta'].get('visium_version')!r}")
    # 无效的 positions（缺 pxl 列）→ 可读的 bad_request，不崩
    root = make_visium("visium_bad", hires=True, lowres=False, invalid_positions=True)
    try:
        sp_io.load_visium(root, {})
        failures.append("缺 pxl_row/col 的 tissue_positions 没有被校验（应报 bad_request）")
    except sc.EngineError as error:
        print(f"  .. 坏 positions → {error.code}: {error.message[:90]}")
        if error.code != "bad_request":
            failures.append(f"坏 positions 的错误码应为 bad_request，实际 {error.code!r}")
    return failures


def check_visium_scale_override() -> list[str]:
    failures: list[str] = []
    root = make_visium("visium_override", hires=True, lowres=False)
    loaded = sp_io.load_visium(root, {"visium_scale": 2.0})
    coords = np.asarray(loaded.get("coords"), dtype=float)
    print(f"  .. visium_scale=2.0 → coords={coords.tolist()} meta.scale={loaded['meta'].get('scale')}")
    if not np.allclose(coords, np.array([[20.0, 20.0], [40.0, 20.0]])):
        failures.append(f"visium_scale=2.0 应覆盖 scalefactor → [[20,20],[40,20]]，实际 {coords.tolist()}")
    return failures


def check_visium_missing_scalefactor() -> list[str]:
    failures: list[str] = []
    root = make_visium("visium_noscale", hires=True, lowres=False,
                       scalefactors={"spot_diameter_fullres": 10.0})
    loaded = sp_io.load_visium(root, {})
    coords = np.asarray(loaded.get("coords"), dtype=float)
    warnings = " ".join(loaded.get("warnings") or [])
    print(f"  .. 缺 scalefactor coords={coords.tolist()} warnings={warnings[:120]}")
    if not np.allclose(coords, np.array([[10.0, 10.0], [20.0, 10.0]])):
        failures.append(f"缺 scalefactor 时不得乘任何比例 → 原坐标，实际 {coords.tolist()}")
    if "scalefactor" not in warnings.lower() and "比例" not in warnings:
        failures.append(f"缺 scalefactor 时应有可读 warning，实际 {warnings!r}")
    return failures


def main() -> int:
    print("=" * 96)
    print("P0-2 I/O：H5AD counts 层 / 矩阵标度元数据 / Visium scalefactor 与输出版本")
    print("=" * 96)
    failures: list[str] = []
    for fn in (check_h5ad_counts_layer, check_visium_lowres_only, check_visium_hires_preferred,
               check_visium_version_validation, check_visium_scale_override,
               check_visium_missing_scalefactor):
        print(f"\n[{fn.__name__}]")
        try:
            failures.extend(fn())
        except Exception:  # noqa: BLE001
            traceback.print_exc()
            failures.append(f"{fn.__name__} 抛异常")
    print("\n" + "=" * 96)
    if failures:
        for item in failures:
            print(f"FAIL  {item}")
        print(f"\nFAIL  P0-2 I/O（{len(failures)} 项）")
        return 1
    print("PASS  P0-2 I/O：h5ad counts/matrix_scale/坐标系 + Visium 按实际图取 scalefactor")
    return 0


if __name__ == "__main__":
    sys.exit(main())
