"""
dsh-spatial — 空间配体-受体分析（DESIGN.md §4 group=空间通讯）

方法（文献）：
  * 配体/受体在簇内的表达比例 × 空间邻接权重，构成「簇对通讯强度」
    （CellPhoneDB 2020 / CellChat 2021 的共表达打分思路；squidpy `ligrec` 的簇对思路）。
  * 置换检验：打乱 spot 的簇标签后重算，得到零分布 → 经验 p 值。
  * 内置 LR 库（人类 symbol，≥80 对）来自 CellPhoneDB / CellChat 公开配对的常用子集；
    只做"共表达 + 邻近"的筛选，不声称等同于官方数据库的完整统计检验。
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import sp_common as sc

# ---------------------------------------------------------------------------
# 内置配体-受体库（ligand, receptor, pathway）
# ---------------------------------------------------------------------------

LR_PAIRS = [
    ("TGFB1", "TGFBR1", "TGFb"), ("TGFB1", "TGFBR2", "TGFb"), ("TGFB2", "TGFBR1", "TGFb"),
    ("IL6", "IL6R", "IL6"), ("IL6", "IL6ST", "IL6"), ("IL10", "IL10RA", "IL10"),
    ("IL1B", "IL1R1", "IL1"), ("IL1A", "IL1R1", "IL1"), ("IL4", "IL4R", "IL4"),
    ("IL13", "IL13RA1", "IL13"), ("IL2", "IL2RA", "IL2"), ("IL7", "IL7R", "IL7"),
    ("IL15", "IL2RB", "IL15"), ("IL17A", "IL17RA", "IL17"), ("IL18", "IL18R1", "IL18"),
    ("IL33", "IL1RL1", "IL33"), ("IL12A", "IL12RB1", "IL12"), ("IL23A", "IL23R", "IL23"),
    ("IFNG", "IFNGR1", "IFNg"), ("IFNG", "IFNGR2", "IFNg"), ("IFNB1", "IFNAR1", "IFNb"),
    ("TNF", "TNFRSF1A", "TNF"), ("TNF", "TNFRSF1B", "TNF"), ("LTA", "TNFRSF1A", "TNF"),
    ("CXCL12", "CXCR4", "CXCL"), ("CXCL12", "ACKR3", "CXCL"), ("CXCL9", "CXCR3", "CXCL"),
    ("CXCL10", "CXCR3", "CXCL"), ("CXCL13", "CXCR5", "CXCL"), ("CXCL8", "CXCR1", "CXCL"),
    ("CXCL8", "CXCR2", "CXCL"), ("CXCL16", "CXCR6", "CXCL"), ("CXCL1", "CXCR2", "CXCL"),
    ("CCL2", "CCR2", "CCL"), ("CCL3", "CCR1", "CCL"), ("CCL4", "CCR5", "CCL"),
    ("CCL5", "CCR5", "CCL"), ("CCL5", "CCR1", "CCL"), ("CCL17", "CCR4", "CCL"),
    ("CCL19", "CCR7", "CCL"), ("CCL21", "CCR7", "CCL"), ("CCL20", "CCR6", "CCL"),
    ("CCL22", "CCR4", "CCL"), ("CCL28", "CCR10", "CCL"),
    ("VEGFA", "KDR", "VEGF"), ("VEGFA", "FLT1", "VEGF"), ("VEGFB", "FLT1", "VEGF"),
    ("VEGFC", "FLT4", "VEGF"), ("PGF", "FLT1", "VEGF"),
    ("EGF", "EGFR", "EGF"), ("TGFA", "EGFR", "EGF"), ("AREG", "EGFR", "EGF"),
    ("EREG", "EGFR", "EGF"), ("HBEGF", "EGFR", "EGF"), ("NRG1", "ERBB3", "NRG"),
    ("HGF", "MET", "HGF"), ("IGF1", "IGF1R", "IGF"), ("IGF2", "IGF1R", "IGF"),
    ("FGF1", "FGFR1", "FGF"), ("FGF2", "FGFR1", "FGF"), ("FGF7", "FGFR2", "FGF"),
    ("FGF10", "FGFR2", "FGF"), ("PDGFA", "PDGFRA", "PDGF"), ("PDGFB", "PDGFRB", "PDGF"),
    ("PDGFC", "PDGFRA", "PDGF"), ("CSF1", "CSF1R", "CSF"), ("CSF2", "CSF2RA", "CSF"),
    ("KITLG", "KIT", "KIT"), ("FLT3LG", "FLT3", "FLT3"),
    ("WNT5A", "FZD5", "WNT"), ("WNT5A", "ROR2", "WNT"), ("WNT2", "FZD1", "WNT"),
    ("WNT3A", "FZD1", "WNT"), ("DKK1", "LRP6", "WNT"),
    ("JAG1", "NOTCH1", "Notch"), ("JAG1", "NOTCH2", "Notch"), ("JAG2", "NOTCH1", "Notch"),
    ("DLL1", "NOTCH1", "Notch"), ("DLL4", "NOTCH1", "Notch"),
    ("SHH", "PTCH1", "Hedgehog"), ("IHH", "PTCH1", "Hedgehog"),
    ("BMP2", "BMPR1A", "BMP"), ("BMP4", "BMPR1A", "BMP"), ("BMP7", "BMPR2", "BMP"),
    ("LAMA1", "ITGB1", "ECM"), ("LAMB1", "ITGB1", "ECM"), ("FN1", "ITGB1", "ECM"),
    ("FN1", "ITGAV", "ECM"), ("COL1A1", "ITGB1", "ECM"), ("COL4A1", "ITGB1", "ECM"),
    ("THBS1", "CD47", "ECM"), ("SPP1", "CD44", "ECM"), ("SPP1", "ITGAV", "ECM"),
    ("MIF", "CD74", "MIF"), ("MIF", "CXCR4", "MIF"), ("HLA-DRA", "CD4", "MHC"),
    ("CD80", "CD28", "Costim"), ("CD86", "CD28", "Costim"), ("CD274", "PDCD1", "Checkpoint"),
    ("PDCD1LG2", "PDCD1", "Checkpoint"), ("LGALS9", "HAVCR2", "Checkpoint"),
    ("CD47", "SIRPA", "Checkpoint"), ("APP", "CD74", "APP"),
    ("ANXA1", "FPR1", "ANXA"), ("SAA1", "FPR2", "SAA"), ("PECAM1", "PECAM1", "Adhesion"),
    ("CDH1", "CDH1", "Adhesion"), ("ICAM1", "ITGAL", "Adhesion"), ("VCAM1", "ITGA4", "Adhesion"),
    ("ADM", "CALCRL", "ADM"), ("EDN1", "EDNRA", "EDN"), ("AGT", "AGTR1", "RAS"),
    ("PENK", "OPRM1", "Opioid"), ("TAC1", "TACR1", "Tachykinin"), ("NPY", "NPY1R", "NPY"),
    ("APOE", "LRP1", "Lipid"), ("APOE", "TREM2", "Lipid"), ("LPL", "LRP1", "Lipid"),
    ("SELENOP", "LRP8", "Selenoprotein"), ("TTR", "RBP4", "Transport"),
    ("GCG", "GCGR", "Glucagon"), ("INS", "INSR", "Insulin"), ("LEP", "LEPR", "Leptin"),
    ("ADIPOQ", "ADIPOR1", "Adiponectin"), ("RETN", "CAP1", "Resistin"),
    ("KRT17", "ITGB4", "Keratin"), ("KRT8", "ITGB4", "Keratin"),
]


def _ligand_receptor_arrays(genes: list[str]):
    index = {gene: i for i, gene in enumerate(genes)}
    picked = []
    for ligand, receptor, pathway in LR_PAIRS:
        if ligand in index and receptor in index:
            picked.append((ligand, receptor, pathway, index[ligand], index[receptor]))
    return picked


def _p(params, key, default=None):
    value = (params or {}).get(key, default)
    return default if value in (None, "", "null") else value


def _p_int(params, key, default):
    try:
        return int(float(_p(params, key, default)))
    except (TypeError, ValueError):
        return int(default)


def _p_float(params, key, default):
    try:
        return float(_p(params, key, default))
    except (TypeError, ValueError):
        return float(default)


def _X(ctx):
    X = ctx.get("X")
    if X is None:
        raise sc.prerequisite("数据集没有表达矩阵")
    return np.asarray(X, dtype=np.float32)


def _genes(ctx):
    return [str(item) for item in (ctx.get("genes") or [])]


def _spots(ctx):
    return [str(item) for item in (ctx.get("spots") or [])]


def _info(ctx) -> pd.DataFrame:
    frame = ctx.get("sample_info")
    return frame if isinstance(frame, pd.DataFrame) else pd.DataFrame(index=pd.Index(_spots(ctx), name="spot"))


def _outdir(ctx, outdir) -> str:
    return sc.ensure_dir(outdir or ctx.get("outdir") or ".")


def _labels(ctx, key: str) -> np.ndarray:
    frame = _info(ctx)
    for candidate in ([key] if key else []) + ["clusters", "domains", ctx.get("cluster_key") or ""]:
        if candidate and candidate in frame.columns:
            return np.asarray([str(v) for v in frame[candidate].tolist()], dtype=object)
    labels = ctx.get("cluster")
    if labels is None:
        raise sc.prerequisite("需要先聚类（或指定 cluster_key）：空间通讯按簇对计算")
    return np.asarray([str(v) for v in labels], dtype=object)


# ---------------------------------------------------------------------------
# 打分
# ---------------------------------------------------------------------------


def _neighbour_matrix(ctx, k: int) -> np.ndarray:
    coords = ctx.get("coords")
    if coords is None:
        raise sc.prerequisite("空间通讯需要空间坐标")
    coords = np.asarray(coords, dtype=float)[:, :2]
    return sc.coord_knn(coords, k)[0]


def _cluster_adjacency(labels: np.ndarray, neighbours: np.ndarray, uniq: list[str]) -> np.ndarray:
    """W[a,b] = 簇 a 的邻居里属于簇 b 的比例。"""
    position = {name: i for i, name in enumerate(uniq)}
    codes = np.array([position[str(item)] for item in labels], dtype=int)
    neigh_codes = codes[neighbours]
    k = len(uniq)
    mat = np.zeros((k, k), dtype=float)
    for a in range(k):
        rows = codes == a
        if not rows.any():
            continue
        for b in range(k):
            mat[a, b] = float(np.mean(neigh_codes[rows] == b))
    return mat


def _fractions(X: np.ndarray, codes: np.ndarray, n_clusters: int, gene_idx: np.ndarray) -> np.ndarray:
    """F[cluster, gene] = 该簇内表达该基因的 spot 比例。"""
    expr = X[:, gene_idx] > 0
    out = np.zeros((n_clusters, gene_idx.size), dtype=float)
    for a in range(n_clusters):
        rows = codes == a
        if rows.any():
            out[a] = expr[rows].mean(axis=0)
    return out


def _pair_scores(F: np.ndarray, W: np.ndarray, l_idx: np.ndarray, r_idx: np.ndarray) -> np.ndarray:
    """返回 S[cluster_a, cluster_b, pair] = sqrt(F_a,L * F_b,R) * W_ab。"""
    left = np.sqrt(np.clip(F[:, l_idx], 0, None))     # (K, P)
    right = np.sqrt(np.clip(F[:, r_idx], 0, None))    # (K, P)
    return left[:, None, :] * right[None, :, :] * W[:, :, None]


def _analyse(ctx, params) -> dict:
    X = _X(ctx)
    genes = _genes(ctx)
    key = str(_p(params, "cluster_key", "") or "clusters")
    labels = _labels(ctx, key)
    uniq = list(dict.fromkeys([str(v) for v in labels.tolist()]))
    if len(uniq) < 2:
        raise sc.prerequisite("只有 1 个簇：空间通讯至少需要 2 个簇/区域")
    k_neighbors = _p_int(params, "n_neighbors", 6)
    min_frac = _p_float(params, "min_frac", 0.1)
    n_perms = _p_int(params, "permutations", 100)
    top_n = _p_int(params, "top_n", 30)
    pairs = _ligand_receptor_arrays(genes)
    if not pairs:
        raise sc.prerequisite("内置 LR 库里没有任何一对同时出现在本数据集的基因中（基因命名可能是 Ensembl ID）")
    neighbours = _neighbour_matrix(ctx, k_neighbors)
    position = {name: i for i, name in enumerate(uniq)}
    codes = np.array([position[str(item)] for item in labels], dtype=int)
    l_idx = np.array([p[3] for p in pairs], dtype=int)
    r_idx = np.array([p[4] for p in pairs], dtype=int)
    F = _fractions(X, codes, len(uniq), np.concatenate([l_idx, r_idx]))
    F = F / np.maximum(F.max(), 1e-9)  # 归一化到 0-1，避免表达量级主导
    F_l, F_r = F[:, :len(pairs)], F[:, len(pairs):]
    W = _cluster_adjacency(labels, neighbours, uniq)
    pair_index = np.arange(len(pairs))
    observed = _pair_scores(F, W, pair_index, pair_index + len(pairs))
    # 零分布：打乱簇标签后重算
    rng = np.random.default_rng(0)
    null_values = []
    for _ in range(max(0, n_perms)):
        perm = rng.permutation(codes)
        F_perm_all = _fractions(X, perm, len(uniq), np.concatenate([l_idx, r_idx]))
        F_perm_all = F_perm_all / np.maximum(F_perm_all.max(), 1e-9)
        W_perm = _cluster_adjacency(labels, neighbours, uniq) if False else _cluster_adjacency_permuted(perm, neighbours, len(uniq))
        scores = _pair_scores(F_perm_all, W_perm, np.arange(len(pairs)), np.arange(len(pairs)) + len(pairs))
        null_values.append(scores.ravel())
    if null_values:
        null = np.concatenate(null_values)
        mu, sd = float(np.mean(null)), float(np.std(null) + 1e-12)
    else:
        mu, sd = 0.0, 1.0
    rows = []
    for p_index, (ligand, receptor, pathway, _li, _ri) in enumerate(pairs):
        for a, name_a in enumerate(uniq):
            for b, name_b in enumerate(uniq):
                if a == b:
                    continue
                score = float(observed[a, b, p_index])
                if score <= 0:
                    continue
                p_lig = float(F_l[a, p_index])
                p_rec = float(F_r[b, p_index])
                if p_lig < min_frac or p_rec < min_frac:
                    continue
                z = (score - mu) / sd
                p_val = float(min(1.0, max(0.0, 0.5 * (1.0 + _erf(-z / np.sqrt(2.0))))))
                rows.append({
                    "ligand": ligand, "receptor": receptor, "pathway": pathway,
                    "cluster_a": name_a, "cluster_b": name_b,
                    "score": score, "z": z, "p": p_val,
                    "pct_ligand": p_lig * 100, "pct_receptor": p_rec * 100,
                    "adjacency": float(W[a, b]),
                })
    if not rows:
        raise sc.prerequisite("没有满足 min_frac 的配体-受体簇对：请降低 min_frac 或换用更细的聚类")
    frame = pd.DataFrame(rows).sort_values(["z", "score"], ascending=False)
    try:
        frame["p_adj"] = sc.bh_fdr(frame["p"].to_numpy())
    except Exception:  # noqa: BLE001
        frame["p_adj"] = frame["p"]
    frame = frame.head(max(10, top_n * 10)).reset_index(drop=True)
    return {"frame": frame, "uniq": uniq, "observed": observed, "pairs": pairs,
            "W": W, "labels": labels, "n_perms": n_perms}


def _cluster_adjacency_permuted(perm: np.ndarray, neighbours: np.ndarray, k: int) -> np.ndarray:
    neigh_codes = perm[neighbours]
    mat = np.zeros((k, k), dtype=float)
    for a in range(k):
        rows = perm == a
        if not rows.any():
            continue
        for b in range(k):
            mat[a, b] = float(np.mean(neigh_codes[rows] == b))
    return mat


def _erf(x):
    from math import erf

    return np.vectorize(erf)(x)


# ---------------------------------------------------------------------------
# 分析
# ---------------------------------------------------------------------------


def analysis_ligrec(ctx, params, outdir):
    out = _analyse(ctx, params)
    frame = out["frame"]
    uniq = out["uniq"]
    outdir = _outdir(ctx, outdir)
    top_n = _p_int(params, "top_n", 30)
    # 热图：cluster_a × cluster_b，值 = 该簇对的最高 z
    mat = np.zeros((len(uniq), len(uniq)), dtype=float)
    for _, row in frame.iterrows():
        a = uniq.index(row["cluster_a"]); b = uniq.index(row["cluster_b"])
        mat[a, b] = max(mat[a, b], float(row["z"]))
    fig, ax = sc.new_fig(6.4, 5.2)
    image = ax.imshow(mat, cmap="viridis")
    ax.set_xticks(np.arange(len(uniq))); ax.set_xticklabels(uniq, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(uniq))); ax.set_yticklabels(uniq, fontsize=8)
    ax.set_xlabel("receiver cluster (receptor)", fontsize=9)
    ax.set_ylabel("sender cluster (ligand)", fontsize=9)
    ax.set_title("Spatial ligand-receptor communication (max z)", fontsize=10)
    bar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    bar.set_label("z-score vs permutation", fontsize=8)
    figure = sc.save_fig(fig, outdir, "lr-heatmap", title=None)
    rows = frame.head(200).to_dict(orient="records")
    significant = int((frame["p_adj"] < 0.05).sum()) if "p_adj" in frame.columns else 0
    bullets = [
        f"{len(out['pairs'])} 对配体-受体 × {len(uniq)} 个簇，得到 {len(frame)} 条候选通讯（置换 {out['n_perms']} 次）",
        "最强信号：" + "、".join(f"{row['ligand']}→{row['receptor']}({row['cluster_a']}→{row['cluster_b']}, z={row['z']:.1f})"
                                 for _, row in frame.head(3).iterrows()),
        f"p_adj<0.05 的候选：{significant} 条",
    ]
    return sc.result(
        summary={"title": "空间配体-受体", "bullets": bullets,
                 "metrics": {"n_pairs": len(out["pairs"]), "n_clusters": len(uniq),
                             "n_calls": int(len(frame)), "n_significant": significant}},
        tables=[sc.table("ligrec", "配体-受体通讯（按 z 排序）",
                         columns=[{"key": "ligand", "label": "配体", "type": "string"},
                                  {"key": "receptor", "label": "受体", "type": "string"},
                                  {"key": "pathway", "label": "通路", "type": "string"},
                                  {"key": "cluster_a", "label": "发送簇", "type": "string"},
                                  {"key": "cluster_b", "label": "接收簇", "type": "string"},
                                  {"key": "score", "label": "得分", "type": "number"},
                                  {"key": "z", "label": "z", "type": "number"},
                                  {"key": "p", "label": "p", "type": "number"},
                                  {"key": "p_adj", "label": "p_adj", "type": "number"},
                                  {"key": "pct_ligand", "label": "配体阳性%", "type": "number"},
                                  {"key": "pct_receptor", "label": "受体阳性%", "type": "number"},
                                  {"key": "adjacency", "label": "空间相邻度", "type": "number"}],
                         rows=rows)],
        figures=[figure],
        data_out=[{"kind": "ligrec_csv", "path": sc.export_csv(frame, outdir, "ligrec.csv", index=False),
                   "label": "配体-受体通讯表"}],
        store={"ligrec": {"frame": frame, "n_pairs": len(out["pairs"])}},
    )


def analysis_lr_map(ctx, params, outdir):
    ligand = str(_p(params, "ligand", "") or "")
    receptor = str(_p(params, "receptor", "") or "")
    genes = _genes(ctx)
    if not ligand or not receptor:
        raise sc.prerequisite("请指定 ligand 与 receptor（如 TGFB1 / TGFBR1）")
    if ligand not in genes or receptor not in genes:
        raise sc.prerequisite(f"基因不在数据集里：{ligand} / {receptor}")
    X = _X(ctx)
    left = np.asarray(X[:, genes.index(ligand)], dtype=float)
    right = np.asarray(X[:, genes.index(receptor)], dtype=float)
    product = left * right
    outdir = _outdir(ctx, outdir)
    values = product
    fig, axes = sc.new_fig(11.0, 4.0)
    fig.set_size_inches(11.0, 4.0)
    fig.clf()
    ax1 = fig.add_subplot(1, 3, 1)
    sc.scatter_spatial(ax1, ctx, values=left, point_size=6.0, legend_title=ligand)
    ax1.set_title(f"Ligand: {ligand}", fontsize=9)
    ax2 = fig.add_subplot(1, 3, 2)
    sc.scatter_spatial(ax2, ctx, values=right, point_size=6.0, legend_title=receptor)
    ax2.set_title(f"Receptor: {receptor}", fontsize=9)
    ax3 = fig.add_subplot(1, 3, 3)
    sc.scatter_spatial(ax3, ctx, values=product, point_size=6.0, legend_title="L x R")
    ax3.set_title("Co-expression", fontsize=9)
    figure = sc.save_fig(fig, outdir, "lr-map", title=None)
    positive = int((product > 0).sum())
    bullets = [f"{ligand} × {receptor}：{positive} 个 spot（{positive / max(len(product), 1) * 100:.1f}%）同时表达两者",
               "共表达位置说明两个细胞群在空间上相邻，可结合 ligrec 表一起看"]
    return sc.result(summary={"title": "通讯强度空间图", "bullets": bullets,
                              "metrics": {"ligand": ligand, "receptor": receptor, "coexp_spots": positive}},
                     figures=[figure])


# ---------------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------------


def plot_lr_heatmap(ctx, params=None, outdir=None):
    params = dict(params or {})
    out = _analyse(ctx, params)
    frame, uniq = out["frame"], out["uniq"]
    mat = np.zeros((len(uniq), len(uniq)), dtype=float)
    for _, row in frame.iterrows():
        mat[uniq.index(row["cluster_a"]), uniq.index(row["cluster_b"])] = max(
            mat[uniq.index(row["cluster_a"]), uniq.index(row["cluster_b"])], float(row["z"]))
    outdir = _outdir(ctx, outdir)
    fig, ax = sc.new_fig(6.4, 5.2)
    image = ax.imshow(mat, cmap="viridis")
    ax.set_xticks(np.arange(len(uniq))); ax.set_xticklabels(uniq, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(uniq))); ax.set_yticklabels(uniq, fontsize=8)
    ax.set_xlabel("receiver cluster", fontsize=9); ax.set_ylabel("sender cluster", fontsize=9)
    ax.set_title("Spatial LR communication (max z)", fontsize=10)
    bar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02); bar.set_label("z", fontsize=8)
    return sc.save_fig(fig, outdir, "lr-heatmap", title=None)


def plot_lr_map(ctx, params=None, outdir=None):
    params = dict(params or {})
    ligand = str(_p(params, "ligand", "") or "")
    receptor = str(_p(params, "receptor", "") or "")
    genes = _genes(ctx)
    if not ligand or not receptor or ligand not in genes or receptor not in genes:
        raise sc.prerequisite("plot lr_map 需要 ligand 与 receptor（都在数据集里）")
    X = _X(ctx)
    product = np.asarray(X[:, genes.index(ligand)], dtype=float) * np.asarray(
        X[:, genes.index(receptor)], dtype=float)
    outdir = _outdir(ctx, outdir)
    fig, ax = sc.new_fig(6.2, 5.0)
    sc.scatter_spatial(ax, ctx, values=product, point_size=6.0, legend_title="L x R")
    ax.set_title(f"{ligand} x {receptor}", fontsize=10)
    return sc.save_fig(fig, outdir, "lr-map", title=None)


ANALYSES = [
    {
        "id": "ligrec_spatial", "title": "空间配体-受体", "group": "空间通讯",
        "method": "LR 共表达比例 + 空间邻接 + 标签置换（CellPhoneDB 2020 / CellChat 2021 思路；squidpy ligrec 语义）",
        "needs": ["cluster"],
        "params": [
            {"key": "cluster_key", "label": "簇列", "type": "text", "default": "clusters"},
            {"key": "n_neighbors", "label": "邻居数", "type": "number", "default": 6},
            {"key": "min_frac", "label": "最低表达比例", "type": "number", "default": 0.1},
            {"key": "permutations", "label": "置换次数", "type": "number", "default": 100},
            {"key": "top_n", "label": "返回条数", "type": "number", "default": 30},
        ],
        "run": analysis_ligrec,
    },
    {
        "id": "lr_map", "title": "通讯强度空间图", "group": "空间通讯",
        "method": "配体×受体共表达映射（内置实现）",
        "needs": ["cluster"],
        "params": [
            {"key": "ligand", "label": "配体", "type": "text", "default": ""},
            {"key": "receptor", "label": "受体", "type": "text", "default": ""},
        ],
        "run": analysis_lr_map,
    },
]

PLOTS = {
    "lr_heatmap": plot_lr_heatmap,
    "lr_map": plot_lr_map,
}

__all__ = ["STEPS", "ANALYSES", "PLOTS", "LR_PAIRS"]
