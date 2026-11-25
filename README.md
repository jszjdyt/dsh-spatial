# dsh-spatial — 空间转录组工作台

DeepSeek Harness 的**空间转录组可视化工作台**插件（Host 工具 `spatial` + 整页 UI），
和 `dsh-geobulk`（bulk）、`dsh-singlecell`（单细胞）同形态，风格与交互对齐。

一句话：**搜 GEO（中文可搜）→ 下 Visium/Stereo-seq/h5ad → 空间质控/聚类 → 空间高变基因/空间域/邻域/反卷积 → 出图出表**，页面上点着用，也可以在对话里让 Agent 一句话跑完。

## 功能

**数据**
- GEO 检索（E-utilities；中文关键词用内置 245 条词典反查成英文，结果太少自动放宽；`spatial_score` 与平台识别；勾选文件下载到指定目录，`.part` 原子写入）
- 导入：10x Visium `outs/`（h5 或 mtx + `spatial/` 坐标 + `scalefactors` + H&E 图）、`visium_h5`、`h5ad`（`obsm['spatial']`）、Stereo-seq BIN、Slide-seq、表达矩阵 + 坐标表、xlsx、GEO series matrix、zip/tar（自动找 Visium/BIN/矩阵）、rds/RData（R 桥）

**常规流程**（12 步，左栏从上到下）
`qc` 空间质控 → `filter` 过滤 → `normalize` 标准化 → `hvg` 高变基因 → `pca` → `neighbors` 邻接图（表达图 + 空间图融合）→ `cluster` Louvain 聚类 → `embed` t-SNE → `svg` 空间高变基因（Moran's I / Geary's C）→ `domains` 空间域（MRF 平滑）→ `deconv` 反卷积（NNLS / DWLS）→ `annotate` 区域注释

**高级分析**（34 项，按组）
- 空间统计：Moran's I、Geary's C、置换检验、双变量空间自相关、Getis-Ord Gi* 热点、空间梯度、空间共表达基因对、Ripley's K/L、全基因自相关扫描
- 邻域与域：邻域富集、共现、MRF 空间域、域标记基因、域组成、边界强度、到域距离、细胞邻域（niche）、邻域多样性
- 反卷积：NNLS、DWLS、细胞类型空间分布、方法一致性、区域组成差异、比例-标记基因一致性
- 空间通讯：配体-受体（内置 ≥100 对 LR + 空间邻接 + 置换）、通讯强度空间图
- 区域差异：域 vs 其余、两域差异、按分组差异、区域表达箱线
- 质控与整合：质控汇总、空间离群点、多切片对齐诊断、空间轨迹（扩散伪时序）

**出图/导出**：31 种图（空间散点/表达图/域图/比例图/邻域热图/Ripley/火山/热图…）；导出表达矩阵、SVG 表、域表、反卷积比例、样本表、汇总 JSON、Methods、Markdown 报告、可复现脚本。

**Agent 工具** `spatial`：`status/search/detail/download/open/demo/inspect/preview/run/plot/analyses/analysis/export/audit/suggest/recipe/samples/refs/ref_attach/formats/reset/close`。

## 安装

```bash
# 1) 依赖（本机 python3 已具备，无需安装）
python3 -c "import numpy,pandas,scipy,sklearn,matplotlib,h5py,anndata,PIL,certifi"

# 2) 安装插件（Full access 或按提示批准）
#    设置 → 插件 → 添加本地插件 → 选 plugins/dsh-spatial
# 或让 Agent 调 plugin_manager: action=install_bundle target=<plugins/dsh-spatial 绝对路径>
```

装好后重载页面：侧边栏出现「空间工作台」，对话里多一个 `spatial` 工具。

## 首次使用

1. 侧边栏打开「空间工作台」→ 点「载入示例数据」（2400 spot × 1000 基因 × 5 个空间域 × 2 切片）。
2. 左栏从上到下点「运行」：质控 → 过滤 → 标准化 → 高变基因 → PCA → 邻接图 → 聚类 → 二维嵌入 → 空间高变基因 → 空间域。
3. 想用真实数据：顶栏「GEO 检索」搜 `乳腺 空间转录组 Visium`（或 `spatial transcriptomics breast cancer`）→ 详情 → 选下载目录 → 下载 → 「导入并开始清洗」；也可以直接填 Visium `outs/` 目录路径。
4. 反卷积：先去「反卷积」tab 选参考 h5ad（可点「从单细胞工作台拉取」），再跑 `deconv`。

## 说明与边界

**计算实现（2026-11 更新：已装真包，默认走主流实现）**

专用环境 `~/.dsh/runtimes/bio312`（Python 3.12，用 uv 建）：`scanpy 1.12.6`、`squidpy 1.8.4`、`anndata 0.13.4`、`leidenalg 0.12`、`igraph 1.0`、`umap-learn`、`scikit-misc`、`harmonypy`。
引擎会自动优先这个解释器；每个步骤/分析的 `method` 字段写的是**真实调用**，例如：

- `qc` → `scanpy.pp.calculate_qc_metrics(qc_vars=[mt])`
- `normalize` → `scanpy.pp.normalize_total` + `scanpy.pp.log1p`
- `hvg` → `scanpy.pp.highly_variable_genes(flavor=seurat_v3)`
- `pca` → `scanpy.pp.scale` + `scanpy.tl.pca`
- `neighbors` → `scanpy.pp.neighbors` + `squidpy.gr.spatial_neighbors`
- `cluster` → `scanpy.tl.leiden(flavor="igraph")`
- `embed` → `scanpy.tl.tsne` / `scanpy.tl.umap`
- `morans_svg` / `gearys_svg` → `squidpy.gr.spatial_autocorr(mode="moran"/"geary")`
- `nhood_enrichment` → `squidpy.gr.nhood_enrichment`；`co_occurrence` → `squidpy.gr.co_occurrence`；`ripley_k` → `squidpy.gr.ripley`
- `ligrec_spatial` → `squidpy.gr.ligrec`（OmniPath；不可达时用内置 LR 库并注明）

**仍不是主流包的地方（squidpy 没有对应 API，已在 method 写明"内置实现"）**：
Getis-Ord Gi* 热点、沿轴空间梯度、双基因空间共表达对、双变量 Moran、全基因扫描的独立置换表、空间离群点/切片对齐诊断/空间轨迹；
空间域用 MRF 平滑（BayesSpace/GraphST 思想，R/TF 实现未装）、反卷积在 RCTD 不可用时用 NNLS/DWLS（Tsoucas 2019 语义）。
**scanpy/squidpy 不在时**（宿主自带 python3）会自动回退到等价内置实现，method 会标"回退内置实现"。

**GEO 检索行为（2026-11 修订）**
- 选「空间平台」（Visium/Stereo-seq/Slide-seq/Xenium/MERFISH/全部）时，平台词会**并入检索式**（如 `胃癌 AND (Visium OR "spatial transcriptomics")`），从源头只回该平台记录；事后过滤也是**严格的**——不命中就返回空并提示，不再回退成一堆无关数据。
- 不选平台时按**空间相关度排序**（high → mid → low），非空间记录仍在但排后面；`共 N 条` 显示的是过滤后的数量。
- NCBI 限流时（HTTP 403/429）会给可读提示并重试；如需高频检索，设置环境变量 `DSH_NCBI_API_KEY`（免费申请）可把配额提到 10 次/秒，`DSH_NCBI_EMAIL` 可选。

其它：
- 反卷积需要单细胞参考（h5ad，带细胞类型列；RCTD 需要 counts）。
- GEO 检索依赖 NCBI 可达；失败给中文原因（重试 ≤2 次）。
- 图内文字用英文 + 中文回落字体（Helvetica/PingFang SC）。
