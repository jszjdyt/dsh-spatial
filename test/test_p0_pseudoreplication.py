#!/usr/bin/env python3
"""P0 回归：空间伪重复（把 spot 当独立样本）与块内置换。

修复前反例（probe_before.py 实测，2 片 × 25 spot、每组 1 片）：
  * sp_domain._analysis_de_groups：6 个基因里 2 个显著，最小 p_adj = 3.896e-09
    —— 每组只有 1 片，spot 被当成独立重复，属于伪重复（pseudoreplication）
  * 切片常数基因（片 A=10、片 B=0）：Moran's I = -0.0204，块间置换 p = 0.0196（假显著）
  * sp_domain._svg_scores 用合并坐标估半径 → fallback 100，把跨片均值差当成空间信号

修复后要求：
  * 每组只有 1 片/供体 → 只出描述性汇总（descriptive），不做推断（无 p 值）
  * 每组 ≥2 片 → 按 供体/切片 × 区域 单位聚合原始 counts，单位级检验（n_units 正确）
  * Moran's I/Geary 先做块内中心化、置换限制在块内；切片常数基因不得判为显著空间自相关

运行：MPLBACKEND=Agg MPLCONFIGDIR=<tmp> python3 test/test_p0_pseudoreplication.py
"""

from __future__ import annotations

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


def make_ctx(n_per_slice: int = 25, slices=("S1", "S2"), slice_groups=("tumor", "adjacent"),
             offset: float = 0.0, n_genes: int = 6, seed: int = 3,
             domains=("D1", "D2")) -> dict:
    """合成多切片数据：每片一块 n×n 规则格点。

    offset=0 时各片坐标完全相同（复现跨片伪邻接）；
    slice_groups 长度必须等于切片数（一片一组 或 两片一组）。
    """
    spacing = 10.0
    cols = int(np.ceil(np.sqrt(n_per_slice)))
    rows = int(np.ceil(n_per_slice / cols))
    base = np.array([[x, y] for y in range(rows) for x in range(cols)],
                    dtype=float)[:n_per_slice] * spacing
    coords, slice_ids, group_ids = [], [], []
    for i, slice_name in enumerate(slices):
        coords.append(base + np.array([i * offset, 0.0]))
        slice_ids += [slice_name] * n_per_slice
        group_ids += [slice_groups[i]] * n_per_slice
    coords = np.vstack(coords)
    n = coords.shape[0]
    rng = np.random.default_rng(seed)
    X = rng.poisson(5.0, size=(n, n_genes)).astype(np.float32)
    # 组效应基因（G1 高、G2 低）：只在组水平，每个切片内部同分布
    uniq_groups = list(dict.fromkeys(slice_groups))
    group_codes = np.array([uniq_groups.index(g) for g in group_ids])
    X[:, 0] += np.where(group_codes == 0, 40.0, 0.0)
    X[:, 1] += np.where(group_codes == 1, 40.0, 0.0)
    # 区域标签（域）：片内按 x 坐标切成两半，两片一致
    region = np.where((np.tile(base[:, 0], len(slices)) < spacing * cols / 2.0), domains[0], domains[1])
    genes = [f"G{j + 1}" for j in range(n_genes)]
    spots = [f"{slice_ids[i]}_{i}" for i in range(n)]
    info = pd.DataFrame({"slice": slice_ids, "group": group_ids,
                         "domains": region,
                         "clusters": ["C1" if i % 2 == 0 else "C2" for i in range(n)]},
                        index=pd.Index(spots, name="spot"))
    return {"genes": genes, "spots": spots, "X": X, "X_kind": "counts", "counts": X.copy(),
            "coords": coords, "coords_kind": "grid", "spatial_mode": "grid",
            "neighbor_radius": None, "image": {}, "sample_info": info,
            "group_key": "group", "groups": list(dict.fromkeys(slice_groups)),
            "species": "human", "obsm": {}, "results": {}, "ref": None,
            "rng": np.random.default_rng(0),
            "workdir": tempfile.mkdtemp(prefix="sp-p0-"),
            "outdir": tempfile.mkdtemp(prefix="sp-p0-"), "notes": []}


def _table_rows(payload: dict) -> list[dict]:
    tables = payload.get("tables") or []
    return list(tables[0].get("rows") or []) if tables else []


def check_one_slide_per_group_descriptive() -> list[str]:
    """每组 1 片 → 只给描述性，不得做推断。"""
    import sp_domain

    failures: list[str] = []
    ctx = make_ctx(slices=("S1", "S2"), slice_groups=("tumor", "adjacent"))
    payload = sp_domain._analysis_de_groups(ctx, {"group_col": "group"}, ctx["outdir"])
    metrics = payload.get("summary", {}).get("metrics") or {}
    rows = _table_rows(payload)
    print(f"  .. metrics: {json.dumps(metrics, ensure_ascii=False, default=str)}")
    print(f"  .. 前 2 行: {rows[:2]}")
    print(f"  .. warnings: {payload.get('warnings')}")

    inference = str(metrics.get("inference") or "")
    if inference != "descriptive":
        failures.append(f"每组 1 片时 metrics.inference 应为 descriptive，实际 {inference!r}")
    bad_p = [row for row in rows if row.get("p_value") is not None or row.get("p_adj") is not None]
    if bad_p:
        failures.append(f"每组 1 片时仍给出推断 p 值：{bad_p[0]}")
    bullets = " ".join((payload.get("summary") or {}).get("bullets") or [])
    if "描述性" not in bullets and "descriptive" not in bullets.lower():
        failures.append(f"summary 未说明降级为描述性：{bullets[:160]}")
    return failures


def check_two_slides_per_group_unit_level() -> list[str]:
    """每片一组伪重复修复 vs 两片一组单位级检验。"""
    import sp_domain

    failures: list[str] = []
    # A) 每组 1 片：spot 级 Wilcoxon 原本会给 3.9e-09 级别的假显著；现在必须描述性
    ctx_bad = make_ctx(slices=("S1", "S2"), slice_groups=("tumor", "adjacent"))
    payload_bad = sp_domain._analysis_de_groups(ctx_bad, {"group_col": "group"}, ctx_bad["outdir"])
    rows_bad = _table_rows(payload_bad)
    min_adj = min((row.get("p_adj") for row in rows_bad
                   if row.get("p_adj") is not None), default=None)
    print(f"  .. 每组 1 片 min p_adj: {min_adj}（修复前 3.896e-09）")
    if min_adj is not None and min_adj < 1e-3:
        failures.append(f"每组 1 片仍报出伪重复 p_adj={min_adj:.3g}（修复前 3.90e-09）")

    # B) 每组 2 片（S1,S2=tumor；S3,S4=adjacent），切片间坐标不重叠
    ctx = make_ctx(n_per_slice=25, slices=("S1", "S2", "S3", "S4"),
                   slice_groups=("tumor", "tumor", "adjacent", "adjacent"), offset=200.0)
    payload = sp_domain._analysis_de_groups(ctx, {"group_col": "group"}, ctx["outdir"])
    metrics = payload.get("summary", {}).get("metrics") or {}
    rows = _table_rows(payload)
    print(f"  .. 2 片/组 metrics: {json.dumps(metrics, ensure_ascii=False, default=str)}")
    print(f"  .. 2 片/组 G1/G2 行: {[r for r in rows if r.get('gene') in ('G1', 'G2')]}")
    inference = str(metrics.get("inference") or "")
    if inference != "unit_level":
        failures.append(f"每组 2 片时应按单位级检验（inference=unit_level），实际 {inference!r}")
    n_units = int(metrics.get("n_units") or 0)
    # 4 片 × 2 个区域 = 8 个单位
    if n_units != 8:
        failures.append(f"单位数应为 4 片 × 2 区域 = 8，实际 {n_units}")
    if int(metrics.get("n_case_units") or 0) != 4 or int(metrics.get("n_ctrl_units") or 0) != 4:
        failures.append(f"每组单位数应为 4/4，实际 {metrics.get('n_case_units')}/{metrics.get('n_ctrl_units')}")
    # 4 vs 4 单位时 Mann-Whitney 精确双侧最小 p = 2/70 ≈ 0.0286；不允许出现 spot 级（1e-18）量级的 p
    finite_p = [row["p_value"] for row in rows if row.get("p_value") is not None]
    if finite_p and min(finite_p) < 0.02:
        failures.append(f"单位级检验出现 spot 级量级的 p={min(finite_p):.3g}（伪重复未修，修复前 5.7e-18）")
    # 组效应基因方向必须正确
    by_gene = {row["gene"]: row for row in rows}
    if "G1" in by_gene and not (by_gene["G1"].get("logFC") or 0) > 0:
        failures.append(f"组效应基因 G1 的 logFC 方向错误：{by_gene['G1']}")
    return failures


def check_within_block_moran() -> list[str]:
    """块内中心化 Moran + 块内置换：切片常数基因不得判为显著。"""
    import sp_spatial
    import sp_domain

    failures: list[str] = []
    ctx = make_ctx(slices=("S1", "S2"), slice_groups=("tumor", "adjacent"))
    n_per = 25
    X = np.asarray(ctx["X"], dtype=np.float32).copy()
    # G1：切片常数（A=10, B=0）—— 合并坐标下会被误判成空间自相关
    X[:, 0] = np.where(np.arange(X.shape[0]) < n_per, 10.0, 0.0)
    # G2：片内 iid 噪声（无空间结构）—— 块内置换 p 不应显著
    rng = np.random.default_rng(11)
    X[:, 1] = rng.normal(5.0, 2.0, size=X.shape[0]).astype(np.float32)
    ctx = dict(ctx, X=X, counts=X.copy())

    view = sp_spatial._view(ctx)
    W, radius, _info = sp_spatial._spatial_weights(view, {})
    stats = sp_spatial._weight_stats(W, view["n"])
    X_sub = sp_spatial._submatrix(view, [0, 1])
    observed = sp_spatial._moran_vector(X_sub, W, view["n"], stats["S0"], blocks=view["blocks"])
    pvals, _null = sp_spatial._moran_permutation(
        X_sub, W, view["n"], stats["S0"], observed, 99, np.random.default_rng(0),
        blocks=view["blocks"])
    print(f"  .. 半径={radius:.3g} 块内 Moran I={np.round(observed, 4)} p_perm={np.round(pvals, 4)}")
    if np.isfinite(observed[0]):
        failures.append(f"切片常数基因的块内 Moran 应为缺失（块内无方差），实际 I={observed[0]:.4f}")
    if np.isfinite(pvals[0]) and pvals[0] < 0.05:
        failures.append(f"切片常数基因块内置换 p={pvals[0]:.4f}（修复前 0.0196，假显著）")
    if np.isfinite(pvals[1]) and pvals[1] < 0.01:
        failures.append(f"片内 iid 噪声的块内置换 p={pvals[1]:.4f} 过小（零分布错误）")

    # sp_domain 的 SVG 路径同样必须块内中心化
    Wd = sp_domain._radius_weights(ctx)
    stats_d, zvals, moments = sp_domain._svg_scores(np.asarray(ctx["X"], dtype=np.float32), Wd,
                                                    "moran", blocks=sp_domain._blocks(ctx)[0])
    pdom = sp_domain._permutation_p(np.asarray(ctx["X"], dtype=np.float32), Wd, np.array([0, 1]),
                                    "moran", 99, np.random.default_rng(0),
                                    blocks=sp_domain._blocks(ctx)[0])
    print(f"  .. sp_domain 块内 stats={np.round(stats_d[:2], 4)} p_perm={np.round(pdom, 4)}")
    if np.isfinite(pdom[0]) and pdom[0] < 0.05:
        failures.append(f"sp_domain 切片常数基因块内置换 p={pdom[0]:.4f}（修复前 0.0196）")
    return failures


def main() -> int:
    print("=" * 96)
    print("P0-1 空间伪重复：单位级检验 + 块内 Moran/置换")
    print("=" * 96)
    failures: list[str] = []
    for fn in (check_one_slide_per_group_descriptive, check_two_slides_per_group_unit_level,
               check_within_block_moran):
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
        print(f"\nFAIL  P0-1 空间伪重复（{len(failures)} 项）")
        return 1
    print("PASS  P0-1 空间伪重复：每组 1 片只出描述性；2 片/组按单位级检验；块内置换")
    return 0


if __name__ == "__main__":
    sys.exit(main())
