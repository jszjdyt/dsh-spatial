#!/usr/bin/env python3
"""P0 回归：同坐标双切片的空间邻接必须按切片（物理块）分别构造。

修复前反例（临时脚本 probe_before.py 实测）：
  * sp_spatial._spatial_weights：半径退化为 fallback 100，2450 条边中 1250 条跨片（51.0%）
  * sp_domain._radius_weights：同上 1250/2450（51.0%）
  * sp_lr._neighbour_matrix（6 近邻）：300 个邻居位中 168 个跨片（56%）
  * sp_common.radius_graph 无法按切片分组，4 点两片测试里 6 条边有 4 条跨片

修复后要求：
  * 块半径估计 = 12.0（5×5 格点、间距 10，1.2×最近邻），不再是 fallback 100
  * 物理图跨片边 = 0；sp_lr 邻居跨片 = 0

运行：MPLBACKEND=Agg MPLCONFIGDIR=<tmp> python3 test/test_p0_slice_adjacency.py
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


def make_ctx(n_per_slice: int = 25, n_slices: int = 2, offset: float = 0.0,
             groups=("tumor", "adjacent"), n_genes: int = 6, seed: int = 3) -> dict:
    """两块（或多块）坐标完全相同/平移的合成数据；group 与切片一一对应。"""
    spacing = 10.0
    cols = int(np.ceil(np.sqrt(n_per_slice)))
    rows = int(np.ceil(n_per_slice / cols))
    base = np.array([[x, y] for y in range(rows) for x in range(cols)],
                    dtype=float)[:n_per_slice] * spacing
    coords, slice_ids, group_ids = [], [], []
    for i in range(n_slices):
        coords.append(base + np.array([i * offset, 0.0]))
        slice_ids += [f"S{i + 1}"] * n_per_slice
        group_ids += [groups[i % len(groups)]] * n_per_slice
    coords = np.vstack(coords)
    n = coords.shape[0]
    rng = np.random.default_rng(seed)
    X = np.zeros((n, n_genes), dtype=np.float32)
    for g in range(n_genes):
        X[:, g] = rng.poisson(5.0, size=n)
    for i, gid in enumerate(group_ids):
        X[i, 0 if gid == groups[0] else 1] += 40.0
    genes = [f"G{j + 1}" for j in range(n_genes)]
    spots = [f"{slice_ids[i]}_{i}" for i in range(n)]
    info = pd.DataFrame({"slice": slice_ids, "group": group_ids,
                         "clusters": ["C1" if i % 2 == 0 else "C2" for i in range(n)]},
                        index=pd.Index(spots, name="spot"))
    return {"genes": genes, "spots": spots, "X": X, "X_kind": "counts", "counts": X.copy(),
            "coords": coords, "coords_kind": "pixel", "spatial_mode": "grid",
            "neighbor_radius": None, "image": {}, "sample_info": info,
            "group_key": "group", "groups": list(groups), "species": "human",
            "obsm": {}, "results": {}, "ref": None, "rng": np.random.default_rng(0),
            "workdir": tempfile.mkdtemp(prefix="sp-p0-"),
            "outdir": tempfile.mkdtemp(prefix="sp-p0-"), "notes": []}


def _report(label: str, value) -> None:
    print(f"  .. {label}: {value}")


def check_common_blocks() -> list[str]:
    """sp_common 的块工具 + radius_graph 的按块分组入口（lead 验收脚本同款断言）。"""
    failures: list[str] = []
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 0.0], [1.0, 0.0]])
    slides = ["A", "A", "B", "B"]
    base_graph = sc.radius_graph(coords, 10.0)
    base = sc.block_edge_counts(base_graph, sc.block_indices(slides)[0])
    _report("radius_graph(base) 边", base)
    grouped = None
    for kwargs in ({"slide": slides}, {"group": slides}, {"groups": slides},
                   {"slides": slides}, {"block": slides}, {"slide_ids": slides}):
        try:
            grouped = sc.radius_graph(coords, 10.0, **kwargs)
            break
        except TypeError:
            continue
    if grouped is None:
        failures.append(f"sp_common.radius_graph 不支持按切片分组（base 跨片 {base['n_cross']}/{base['n_edges']}）")
    else:
        counts = sc.block_edge_counts(grouped, sc.block_indices(slides)[0])
        _report("radius_graph(按块) 边", counts)
        if counts["n_cross"] != 0:
            failures.append(f"radius_graph 分组后仍有跨片边 {counts['n_cross']}/{counts['n_edges']}")

    # 块半径：25 点 5×5 格点、间距 10 → 12.0（不是 fallback 100）
    ctx = make_ctx()
    blocks, names, key, warns = sc.spatial_blocks(ctx)
    _report("blocks", f"{len(blocks)} 块 / key={key} / 名称={names}")
    if len(blocks) != 2 or key != "slice":
        failures.append(f"spatial_blocks 识别错误：blocks={len(blocks)} key={key}")
    radius = sc.blocked_radius_estimate(ctx["coords"], blocks)
    _report("blocked_radius_estimate", radius)
    if not np.isfinite(radius) or abs(radius - 12.0) > 1e-6:
        failures.append(f"块半径估计应为 12.0，实际 {radius!r}（修复前 fallback=100）")
    plain = sc.estimate_radius(ctx["coords"])
    _report("estimate_radius(合并坐标，修复前口径)", plain)
    return failures


def check_spatial_weights() -> list[str]:
    import sp_spatial

    failures: list[str] = []
    ctx = make_ctx()
    blocks, _names, _key, _warns = sc.spatial_blocks(ctx)
    view = sp_spatial._view(ctx)
    W, radius, info = sp_spatial._spatial_weights(view, {})
    counts = sc.block_edge_counts(W, view["blocks"])
    _report("spatial_weights radius", radius)
    _report("spatial_weights 边", counts)
    _report("spatial_weights info", {k: v for k, v in info.items() if k not in ("blocks",)})
    if not np.isfinite(radius) or abs(radius - 12.0) > 1e-6:
        failures.append(f"sp_spatial 半径应为 12.0，实际 {radius!r}（修复前 100.0）")
    if counts["n_cross"] != 0:
        failures.append(
            f"sp_spatial 跨片边 {counts['n_cross']}/{counts['n_edges']}"
            f"（{counts['cross_fraction'] * 100:.1f}%，修复前 1250/2450=51.0%）")
    return failures


def check_domain_weights() -> list[str]:
    import sp_domain

    failures: list[str] = []
    ctx = make_ctx()
    _blocks, _names, _key, _warns = sc.spatial_blocks(ctx)
    W = sp_domain._radius_weights(ctx)
    counts = sc.block_edge_counts(W, sc.block_indices(ctx["sample_info"]["slice"].tolist())[0])
    _report("sp_domain._radius_weights 边", counts)
    if counts["n_edges"] == 0:
        failures.append("sp_domain._radius_weights 没有任何边")
    if counts["n_cross"] != 0:
        failures.append(
            f"sp_domain 跨片边 {counts['n_cross']}/{counts['n_edges']}"
            f"（{counts['cross_fraction'] * 100:.1f}%，修复前 1250/2450=51.0%）")
    return failures


def check_lr_neighbours() -> list[str]:
    import sp_lr

    failures: list[str] = []
    ctx = make_ctx()
    n_per = 25
    neigh = sp_lr._neighbour_matrix(ctx, 6)
    left = np.arange(neigh.shape[0])[:, None] < n_per
    right = neigh < n_per
    cross = int(np.sum(left != right))
    _report("sp_lr 邻居矩阵", f"shape={neigh.shape} 跨片={cross}/{neigh.size}")
    if cross != 0:
        failures.append(f"sp_lr 邻居跨片 {cross}/{neigh.size}"
                        f"（{cross / neigh.size * 100:.1f}%，修复前 168/300=56%）")
    return failures


def check_prep_physical_graph() -> list[str]:
    import sp_prep

    failures: list[str] = []
    ctx = make_ctx()
    if not hasattr(sp_prep, "_physical_graph"):
        failures.append("sp_prep 没有 _physical_graph（物理图未按块构造）")
        return failures
    graph, label, info = sp_prep._physical_graph(ctx, 15)
    counts = sc.block_edge_counts(graph, info.get("blocks") or [])
    _report("sp_prep 物理图", f"label={label} info={ {k: v for k, v in info.items() if k != 'blocks'} }")
    if counts["n_cross"] != 0:
        failures.append(f"sp_prep 物理图跨片边 {counts['n_cross']}/{counts['n_edges']}")
    # 融合图里物理部分按块（表达 kNN 允许跨块），所以只验说明文字与物理图本身
    fused, weight, desc = sp_prep._fused_graph(ctx, {"n_neighbors": 15, "spatial_weight": 0.3})
    _report("sp_prep 融合图说明", desc)
    if "分块" not in desc and "block" not in desc.lower():
        failures.append(f"融合图说明没有标注物理图分块：{desc}")
    return failures


def main() -> int:
    print("=" * 96)
    print("P0-1 空间伪重复：同坐标双切片不得跨片建边")
    print("=" * 96)
    failures: list[str] = []
    for fn in (check_common_blocks, check_spatial_weights, check_domain_weights,
               check_lr_neighbours, check_prep_physical_graph):
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
    print("PASS  P0-1 空间伪重复：物理图跨片边 0，块半径 12.0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
