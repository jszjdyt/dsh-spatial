#!/usr/bin/env python3
"""dsh-spatial 集成回归测试：demo → 12 步 → 全部图 → 全部分析 → 导出。

运行：MPLBACKEND=Agg MPLCONFIGDIR=<tmp> python3 test/py_pipeline_test.py
不依赖网络、不依赖 npm；失败抛 AssertionError 并以非 0 退出。
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

import sp_registry as R  # noqa: E402
import sp_state  # noqa: E402

STEPS = ["qc", "filter", "normalize", "hvg", "pca", "neighbors", "cluster", "embed",
         "svg", "domains"]

#: 需要显式参数的图/分析（其余用默认参数）
ANALYSIS_PARAMS: dict = {}
SKIP_PLOTS = {"deconv_bar", "deconv_spatial", "lr_map"}  # 需要单细胞参考或指定基因对
SKIP_ANALYSES = {
    "deconv_nnls", "deconv_dwls", "deconv_spatial_map", "deconv_compare",
    "deconv_roi", "deconv_marker_corr",  # 需要单细胞参考，见下方 deconv 专项
    "lr_map",                            # 需要指定配体/受体
    "cross_moran",                       # 需要两个基因
}


def _ok(label: str, value=None) -> None:
    print(f"  OK  {label}" + (f" — {value}" if value is not None else ""))


def main() -> int:
    workdir = tempfile.mkdtemp(prefix="sp-pipeline-")
    state = sp_state.State()
    dataset = sp_state.make_demo(state, name="demo-spatial")
    assert dataset.X.shape == (2400, 1000), f"demo 形状异常：{dataset.X.shape}"
    _ok("demo", f"{dataset.X.shape[0]} spot × {dataset.X.shape[1]} 基因")

    for step in STEPS:
        payload = R.run_step(dataset, step, {}, workdir, {})
        summary = payload.get("summary") or {}
        assert summary.get("bullets"), f"步骤 {step} 没有 summary.bullets"
        _ok(f"step {step}", summary.get("title"))

    assert "clusters" in dataset.sample_info.columns, "聚类没有写回 clusters 列"
    from sklearn.metrics import adjusted_rand_score

    if "domain_true" in dataset.sample_info.columns:
        ari = adjusted_rand_score(dataset.sample_info["domain_true"].astype(str),
                                  dataset.sample_info["clusters"].astype(str))
        assert ari >= 0.4, f"聚类与注入域的 ARI 太低：{ari:.3f}"
        _ok("ARI(clusters, domain_true)", f"{ari:.3f}")

    # 全部图
    kinds = sorted(R.all_plots().keys())
    failed = []
    for kind in kinds:
        if kind in SKIP_PLOTS:
            continue
        try:
            figure = R.run_plot(dataset, kind, {}, workdir, {})
            assert os.path.exists(figure["file"]) and os.path.getsize(figure["file"]) > 1000
        except Exception as error:  # noqa: BLE001
            failed.append((kind, f"{type(error).__name__}: {error}"))
    assert not failed, f"图失败：{failed}"
    _ok("plots", f"{len(kinds) - len(SKIP_PLOTS)} 种可出图（跳过需外部输入的 {len(SKIP_PLOTS)} 种）")

    # 全部分析（跳过需要外部输入的）
    analyses = R.all_analyses()
    failed = []
    for aid in sorted(analyses):
        if aid in SKIP_ANALYSES:
            continue
        params = dict(ANALYSIS_PARAMS.get(aid) or {})
        if aid == "de_domain":
            params.setdefault("domain", sorted(set(dataset.sample_info["domains"]))[0])
        try:
            payload = R.run_analysis(dataset, aid, params, workdir, {})
            assert (payload.get("summary") or {}).get("bullets"), f"{aid} 没有 summary"
        except Exception as error:  # noqa: BLE001
            failed.append((aid, f"{type(error).__name__}: {error}"))
    assert not failed, f"分析失败：{failed}"
    _ok("analyses", f"{len(analyses) - len(SKIP_ANALYSES)} 项可跑（跳过需外部输入的 {len(SKIP_ANALYSES)} 项）")

    # 反卷积专项：用 demo 的簇均值造一个小参考
    try:
        import anndata as ad
        import pandas as pd

        labels = dataset.sample_info["clusters"].astype(str).to_numpy()
        uniq = sorted(set(labels))[:4]
        X = np.asarray(dataset.X, dtype=np.float32)
        rng = np.random.default_rng(0)
        rows, obs = [], []
        for label in uniq:
            base = X[labels == label].mean(axis=0)
            for _ in range(60):
                rows.append(np.maximum(base + rng.normal(0, base.std() * 0.3 + 0.1, size=base.shape), 0))
                obs.append(label)
        ref_X = np.log1p(np.asarray(rows, dtype=np.float32))
        adata = ad.AnnData(X=ref_X,
                           obs=pd.DataFrame({"cell_type": obs},
                                            index=[f"c{i}" for i in range(len(obs))]),
                           var=pd.DataFrame(index=list(dataset.genes)))
        ref_path = os.path.join(workdir, "ref4.h5ad")
        adata.write_h5ad(ref_path)
        payload = R.run_step(dataset, "deconv", {"ref": ref_path}, workdir, {})
        assert (payload.get("summary") or {}).get("bullets")
        props = (payload.get("store") or {})
        _ok("step deconv", "参考 4 类细胞")
        for aid in ("deconv_nnls", "deconv_dwls", "deconv_spatial_map", "deconv_compare"):
            out = R.run_analysis(dataset, aid, {"ref": ref_path}, workdir, {})
            assert (out.get("summary") or {}).get("bullets")
        _ok("deconv analyses", "4 项可跑")
    except ImportError:
        print("  SKIP deconv（没有 anndata）")

    # 导出
    for what in ("svg_csv", "summary_json", "methods"):
        try:
            out = R.export_dataset(dataset, what, None, workdir) if hasattr(R, "export_dataset") \
                else R.export(dataset, what, None, workdir)
            assert os.path.exists(out["file"]), f"{what} 导出文件不存在"
        except AttributeError:
            # registry 的导出函数名可能不同
            fn = getattr(R, "export", None) or getattr(R, "export_dataset", None)
            assert fn is not None
            out = fn(dataset, what, None, workdir)
        _ok(f"export {what}", os.path.basename(out["file"]))

    score = R.audit(dataset)
    assert 0 <= score["score"] <= 100
    _ok("audit", f"{score['score']}/100")
    print("\nPASS  dsh-spatial pipeline")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print(f"\nFAIL  {error}")
        sys.exit(1)
    except Exception:  # noqa: BLE001
        traceback.print_exc()
        print("\nFAIL  未预期异常")
        sys.exit(1)
