#!/usr/bin/env python3
"""P0 回归：反卷积的 NA 处理（RCTD 未拟合/拒绝的 spot）与区域组成的单位级模型。

修复前反例（sp_deconv._deconvolve_rctd 代码路径，修复前行为）：
  * RCTD 对齐后缺失的 spot（create.RCTD 按 UMI_min/counts_MIN 过滤掉）→
    「用其余 spot 的平均比例填充」，并写进比例表参与下游汇总
  * 全 0 权重行 → 「比例记为均匀分布」fill 成 1/K
  → 覆盖率被掩盖：下游 deconv_roi / 图 / 表都把未拟合 spot 当成真实观测

修复后要求：
  * 未拟合/被拒绝的 spot 保留 NA + fit_status + spot_class，不做均值/均匀分布填补
  * 报告反卷积覆盖率（n_fitted / n_spots）
  * deconv_roi 改为 供体/切片 × 区域 单位级模型；每组只有 1 个供体/切片时只给描述性

运行：MPLBACKEND=Agg MPLCONFIGDIR=<tmp> python3 test/test_p0_deconv_na.py
"""

from __future__ import annotations

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
import sp_deconv  # noqa: E402


def check_rctd_na_retention() -> list[str]:
    """_rctd_align_weights：缺失/拒绝/全 0 行必须保留 NA + fit_status + spot_class。"""
    failures: list[str] = []
    if not hasattr(sp_deconv, "_rctd_align_weights"):
        return ["sp_deconv 没有 _rctd_align_weights（RCTD 后处理仍是均值/均匀分布填补）"]
    spots = ["s0", "s1", "s2", "s3", "s4"]
    weights = pd.DataFrame(
        [[0.4, 0.4, 0.2], [np.nan, np.nan, np.nan], [0.0, 0.0, 0.0], [0.1, 0.1, 0.8]],
        index=["s0", "s1", "s2", "s3"], columns=["T", "B", "M"])
    results_df = pd.DataFrame({"spot_class": ["singlet", "reject", "doublet", "singlet"]},
                              index=["s0", "s1", "s2", "s3"])
    out = sp_deconv._rctd_align_weights(weights, spots, results_df)
    props = np.asarray(out["proportions"], dtype=float)
    print(f"  .. proportions=\n{np.round(props, 3)}")
    print(f"  .. fit_status={out['fit_status']}")
    print(f"  .. spot_class={out['spot_class']}")
    print(f"  .. coverage={out['coverage']}")
    if props.shape != (5, 3):
        failures.append(f"比例矩阵形状应为 (5,3)，实际 {props.shape}")
        return failures
    # s0/s3 正常：行和归一
    if not np.isclose(np.nansum(props[0]), 1.0) or not np.isclose(np.nansum(props[3]), 1.0):
        failures.append(f"已拟合 spot 的行和应为 1，实际 {np.nansum(props[0])!r}/{np.nansum(props[3])!r}")
    if not np.allclose(np.nan_to_num(props[3]), np.array([0.1, 0.1, 0.8]), atol=1e-6):
        failures.append(f"已拟合 spot 的比例被改动：{props[3]}")
    # 未拟合/拒绝/全 0/缺失 spot：必须整行 NA，不得均值或均匀分布填补
    for row, why in ((1, "NaN 权重"), (2, "全 0 权重"), (4, "缺失 spot")):
        if not np.isnan(props[row]).all():
            failures.append(f"{why} 的 spot 行未被保留为 NA：{props[row]}（修复前用均值/均匀分布填补）")
    expected_status = ["ok", "rejected", "zero_weights", "ok", "missing"]
    if list(out["fit_status"]) != expected_status:
        failures.append(f"fit_status 应为 {expected_status}，实际 {out['fit_status']}")
    expected_class = ["singlet", "reject", "doublet", "singlet", "NA"]
    if list(out["spot_class"]) != expected_class:
        failures.append(f"spot_class 应为 {expected_class}，实际 {out['spot_class']}")
    coverage = out["coverage"] or {}
    print(f"  .. coverage keys={sorted(coverage)}")
    if int(coverage.get("n_spots") or 0) != 5:
        failures.append(f"coverage.n_spots 应为 5，实际 {coverage.get('n_spots')}")
    if int(coverage.get("n_fitted") or 0) != 2:
        failures.append(f"coverage.n_fitted 应为 2（s0/s3），实际 {coverage.get('n_fitted')}")
    if coverage.get("fraction_fitted") is None or abs(float(coverage["fraction_fitted"]) - 0.4) > 1e-9:
        failures.append(f"覆盖率应为 0.4，实际 {coverage.get('fraction_fitted')}")
    warn_text = " ".join(out.get("warnings") or [])
    if "覆盖率" not in warn_text and "coverage" not in warn_text.lower():
        failures.append(f"warnings 应报告反卷积覆盖率，实际 {warn_text!r}")
    return failures


def make_ctx(slices=("S1", "S2", "S3", "S4"),
             slice_groups=("tumor", "tumor", "adjacent", "adjacent"),
             n_per_slice: int = 25, nan_spots: int = 0) -> dict:
    spacing = 10.0
    cols = int(np.ceil(np.sqrt(n_per_slice)))
    rows = int(np.ceil(n_per_slice / cols))
    base = np.array([[x, y] for y in range(rows) for x in range(cols)],
                    dtype=float)[:n_per_slice] * spacing
    coords, slice_ids, group_ids = [], [], []
    for i, name in enumerate(slices):
        coords.append(base + np.array([i * 200.0, 0.0]))
        slice_ids += [name] * n_per_slice
        group_ids += [slice_groups[i]] * n_per_slice
    coords = np.vstack(coords)
    n = coords.shape[0]
    spots = [f"{slice_ids[i]}_{i}" for i in range(n)]
    # 2 个区域：片内按 x 切成左右两半
    region_x = np.tile(base[:, 0], len(slices))
    region = np.where(region_x < spacing * cols / 2.0, "D1", "D2")
    info = pd.DataFrame({"slice": slice_ids, "group": group_ids, "domains": region,
                         "clusters": ["C1" if i % 2 == 0 else "C2" for i in range(n)]},
                        index=pd.Index(spots, name="spot"))
    types = ["T", "B", "M"]
    props = np.zeros((n, len(types)), dtype=np.float32)
    signal = {"tumor": (0.55, 0.25, 0.20), "adjacent": (0.20, 0.60, 0.20)}
    rng = np.random.default_rng(7)
    for i in range(n):
        base_prop = np.array(signal[group_ids[i]], dtype=float)
        noise = rng.normal(0, 0.02, size=len(types))
        val = np.clip(base_prop + noise, 1e-6, None)
        props[i] = (val / val.sum()).astype(np.float32)
    if nan_spots:
        props[np.arange(0, n, max(1, n // nan_spots))[:nan_spots]] = np.nan
    genes = ["G1", "G2", "G3"]
    X = rng.poisson(5.0, size=(n, len(genes))).astype(np.float32)
    return {"genes": genes, "spots": spots, "X": X, "X_kind": "counts", "counts": X.copy(),
            "coords": coords, "coords_kind": "grid", "spatial_mode": "grid",
            "neighbor_radius": None, "image": {}, "sample_info": info,
            "group_key": "group", "groups": ["tumor", "adjacent"], "species": "human",
            "obsm": {}, "results": {"deconv": {"cell_types": types, "proportions": props,
                                               "method": "test-rctd"}},
            "ref": None, "rng": np.random.default_rng(0),
            "workdir": tempfile.mkdtemp(prefix="sp-p0-"),
            "outdir": tempfile.mkdtemp(prefix="sp-p0-"), "notes": []}


def _rows(payload: dict) -> list[dict]:
    tables = payload.get("tables") or []
    return list(tables[0].get("rows") or []) if tables else []


def check_deconv_roi_unit_level() -> list[str]:
    failures: list[str] = []
    ctx = make_ctx()
    payload = sp_deconv.analysis_deconv_roi(ctx, {"group_col": "group"}, ctx["outdir"])
    metrics = payload.get("summary", {}).get("metrics") or {}
    rows = _rows(payload)
    print(f"  .. metrics={metrics}")
    print(f"  .. rows={rows}")
    print(f"  .. warnings={payload.get('warnings')}")
    inference = str(metrics.get("inference") or "")
    if inference != "unit_level":
        failures.append(f"2 片/组时 deconv_roi 应为单位级模型（inference=unit_level），实际 {inference!r}")
    if int(metrics.get("n_units") or 0) != 8:
        failures.append(f"单位数应为 4 片 × 2 区域 = 8，实际 {metrics.get('n_units')}")
    if int(metrics.get("n_case_units") or 0) != 4 or int(metrics.get("n_ctrl_units") or 0) != 4:
        failures.append(f"单位数应 4/4，实际 {metrics.get('n_case_units')}/{metrics.get('n_ctrl_units')}")
    # 单位级 4v4：p 不应是 spot 级量级
    finite_p = [row["p"] for row in rows if row.get("p") is not None]
    if finite_p and min(finite_p) < 1e-3:
        failures.append(f"deconv_roi 仍是 spot 级 p={min(finite_p):.3g}（伪重复）")
    by_type = {row["cell_type"]: row for row in rows}
    if "T" in by_type:
        means = {k: v for k, v in by_type["T"].items() if str(k).startswith("mean_")}
        print(f"  .. T 各均值: {means}")
        tumor_mean = by_type["T"].get("mean_tumor")
        adj_mean = by_type["T"].get("mean_adjacent")
        if tumor_mean is None or adj_mean is None or not (tumor_mean > adj_mean):
            failures.append(f"T 在 tumor 组应更高，实际 tumor={tumor_mean} adjacent={adj_mean}")
    return failures


def check_deconv_roi_descriptive_guard() -> list[str]:
    failures: list[str] = []
    ctx = make_ctx(slices=("S1", "S2"), slice_groups=("tumor", "adjacent"))
    payload = sp_deconv.analysis_deconv_roi(ctx, {"group_col": "group"}, ctx["outdir"])
    metrics = payload.get("summary", {}).get("metrics") or {}
    rows = _rows(payload)
    print(f"  .. 每组 1 片 metrics={metrics}")
    print(f"  .. rows={rows}")
    if str(metrics.get("inference") or "") != "descriptive":
        failures.append(f"每组 1 片时应降级为描述性，实际 inference={metrics.get('inference')!r}")
    bad = [row for row in rows if row.get("p") is not None or row.get("p_adj") is not None]
    if bad:
        failures.append(f"每组 1 片时仍给推断 p 值：{bad[0]}")
    bullets = " ".join((payload.get("summary") or {}).get("bullets") or [])
    if "描述" not in bullets and "descriptive" not in bullets.lower():
        failures.append(f"summary 未说明描述性降级：{bullets[:160]}")
    return failures


def check_deconv_roi_nan_coverage() -> list[str]:
    failures: list[str] = []
    ctx = make_ctx(nan_spots=8)
    props = np.asarray(ctx["results"]["deconv"]["proportions"], dtype=float)
    n_nan = int(np.isnan(props).any(axis=1).sum())
    print(f"  .. 注入 NaN spot 数={n_nan}")
    payload = sp_deconv.analysis_deconv_roi(ctx, {"group_col": "group"}, ctx["outdir"])
    metrics = payload.get("summary", {}).get("metrics") or {}
    bullets = " ".join((payload.get("summary") or {}).get("bullets") or [])
    print(f"  .. metrics={metrics} bullets={bullets[:200]}")
    coverage = metrics.get("coverage") if isinstance(metrics.get("coverage"), dict) else {}
    if not coverage:
        failures.append(f"deconv_roi 未报告覆盖率：metrics={metrics}")
    else:
        used = int(coverage.get("n_spots_used") or 0)
        if used != props.shape[0] - n_nan:
            failures.append(f"覆盖率里使用的 spot 数应为 {props.shape[0] - n_nan}，实际 {used}")
    if "覆盖" not in bullets and "coverage" not in bullets.lower():
        failures.append("summary 未说明 NA spot 的覆盖率")
    return failures


def check_composition_figure_with_nan() -> list[str]:
    """图/汇总必须容忍 NaN（nanmean），不得因 NaN 崩溃或静默变 0。"""
    failures: list[str] = []
    ctx = make_ctx(nan_spots=8)
    props = np.asarray(ctx["results"]["deconv"]["proportions"], dtype=float)
    types = ctx["results"]["deconv"]["cell_types"]
    try:
        figure = sp_deconv._fig_composition(ctx, {}, ctx["outdir"], props, types,
                                            name="p0-nan-bar", title="NaN tolerance")
        print(f"  .. figure={figure.get('file')}")
        if not figure or not os.path.exists(figure.get("file", "")):
            failures.append("含 NaN 的比例组成图没有生成")
    except Exception as error:  # noqa: BLE001
        failures.append(f"含 NaN 的比例组成图崩溃：{type(error).__name__}: {error}")
    return failures


def main() -> int:
    print("=" * 96)
    print("P0-2 反卷积：NA 保留/覆盖率 + deconv_roi 供体×区域单位级模型")
    print("=" * 96)
    failures: list[str] = []
    for fn in (check_rctd_na_retention, check_deconv_roi_unit_level,
               check_deconv_roi_descriptive_guard, check_deconv_roi_nan_coverage,
               check_composition_figure_with_nan):
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
        print(f"\nFAIL  P0-2 反卷积 NA（{len(failures)} 项）")
        return 1
    print("PASS  P0-2 反卷积：未拟合 spot 保留 NA/fit_status/spot_class；deconv_roi 单位级")
    return 0


if __name__ == "__main__":
    sys.exit(main())
