# dsh-spatial 内部契约（空间转录组工作台）

> 本文件是并行开发者的**唯一事实来源**。实现与契约冲突时以契约为准；确需改协议，先改本文件再改代码。
> 参考实现（只读，不要改）：`../dsh-singlecell/`（同形态工作台）、`../dsh-geobulk/`（同形态 + GEO + jobs）。
> 文献与产品基线见 §4；每一步/分析都要能说清"用了什么方法、依据哪篇文献/哪个产品"。

## 0. 目标

在 DeepSeek Harness Web GUI 里提供一个**可视化、可操作、Agent 可驱动**的空间转录组工作台：

1. **GEO 检索**：中文/英文关键词检索空间转录组数据集（NCBI E-utilities），中文化字段，筛选（物种/是否有矩阵/空间平台/样本数），看图选目录下载。
2. **导入**：10x Visium（`outs/` 目录、`.h5`、坐标 csv、H&E 图）、Stereo-seq BIN、Slide-seq、`h5ad`（`obsm['spatial']`）、表达矩阵 + 坐标表、GEO 补充文件（zip/tar/系列矩阵）。
3. **常规流程**：空间质控 → 过滤 → 标准化 → 高变基因 → PCA → 邻接图 → 聚类 → 二维嵌入 → 空间高变基因（SVG）。
4. **空间高级分析**：空间域（MRF）、邻域富集、共现、Ripley's K、热点（Getis-Ord）、空间梯度、区域标记/区域差异、边界强度、到区域距离、细胞邻域（niche）、反卷积（NNLS/DWLS，参考可取自 dsh-singlecell）、空间配体-受体、空间轨迹、多切片对齐诊断。
5. **出图与数据**：空间散点/表达图/域图/比例图/曲线图；导出矩阵、SVG 表、域表、反卷积比例、汇总 JSON、Methods、报告、可复现脚本。
6. UI 风格对齐 singlecell/geobulk：清新淡雅 + 高科技感，浅色、留白、细边框、微渐变、等宽数字。

**硬性要求**：Python 未就绪时 UI 给可读原因而不是白屏；单步失败不影响整页；重计算有超时与进度；无网络时 GEO 给可读降级。

---

## 1. 目录与写权限（严格隔离）

```
plugins/dsh-spatial/
  DESIGN.md                      # Lead（本文件）
  README.md / INSTALL.md         # Lead
  package.json / cordis.patch.yml / icon.svg / locale/{zh,en}.json   # Lead
  index.js                       # Lead  Host 入口（/spatial 前缀路由 + spatial 工具 + prompt section）
  lib/impl.mjs                   # Lead  唯一操作实现（HTTP + 工具共用）
  lib/engine.mjs                 # Lead  Python worker 管理（探测/spawn/队列/超时/重启）
  lib/sp_engine.py               # Lead  worker 主循环 + 命令分发（协议 §2）
  lib/sp_registry.py             # Lead  模块汇总 + ctx 构造 + 步骤/分析调度 + 结果整形
  lib/sp_common.py               # sp-core  EngineError / table / save_fig / jsonable / 文本工具
  lib/sp_state.py                # sp-core  Dataset、坐标容器、样本表、demo、history
  lib/sp_io.py                   # sp-core  全格式导入（Visium/Stereo/Slide-seq/h5ad/矩阵+坐标）
  lib/sp_geo.py                  # sp-core  GEO 检索/详情/下载（中文化）
  lib/sp_prep.py                 # sp-base  步骤 qc/filter/normalize/hvg/pca/neighbors/cluster/embed/annotate
  lib/sp_plots.py                # sp-base  基础图（qc_spatial/qc_violin/pca/embedding/cluster_bar/heatmap/dotplot/gene_box/pca_variance）
  lib/sp_spatial.py              # sp-stat  空间统计（Moran/Geary/置换/双变量/热点/梯度/共表达对/Ripley/全基因扫描）
  lib/sp_domain.py               # sp-adv   步骤 svg/domains + 域/邻域/边界/距离/niche/多样性/区域差异
  lib/sp_deconv.py               # sp-adv   步骤 deconv + 反卷积分析
  lib/sp_lr.py                   # sp-adv   空间配体-受体
  client.js                      # sp-ui    浏览器 UI
  test/                          # 各自作用域
  examples/ figures/             # Lead
```

- 只写自己作用域的文件；改别人文件一律禁止，需要改就发消息。
- **写之前先 `read`**。
- 文件间只允许 `import sp_common`（以及 A1 交付的 `sp_state` 的**纯数据类**）；模块之间不要互相 import，避免循环。

---

## 2. Python worker 协议

### 2.1 启动与传输（与 dsh-geobulk 完全一致）

- Host 以 `python3 -u lib/sp_engine.py --workdir <abs_dir>` 启动；`cwd = lib/`。
- stdin 一行一个 JSON 请求 `{"id":"..","cmd":"..","args":{..}}`；stdout 一行一个 JSON 响应：
  成功 `{"id":..,"ok":true,"result":{..}}`；失败 `{"id":..,"ok":false,"error":{"code":"..","message":"中文","trace":"可选","candidates":[...]?}}`；
  进度 `{"id":..,"event":"progress","message":"..","pct":0-100}`。**stdout 只允许协议行**，日志一律 stderr。
- `code` ∈ `bad_request|not_found|prerequisite|unsupported|ambiguous|internal|timeout|no_python|engine_exited|busy|disposed`。
- 单条命令异常只影响该条；`shutdown` 或 stdin EOF 退出。

### 2.2 命令表

| cmd | args | 返回 |
|---|---|---|
| `ping` | `{}` | `{python,platform,mode:"builtin",workdir,packages:{numpy,pandas,scipy,sklearn,matplotlib,h5py,anndata,PIL,networkx},capabilities:{steps:[],analyses:[],plots:[],formats:[]},geo:true,ref:{available}}` |
| `open` | `{path?,url?,name?,data_format?,member?,object?,transpose?,delimiter?,gene_col?,coords?,sample_sheet?,image?,visium_scale?,demo?}` | `{dataset:<摘要>, preview:{...}}` |
| `inspect` | 同 open 的探测子集 | `{format,n_spots?,n_genes?,geometry?,candidates:[],notes:[],size}` |
| `demo` | `{name?}` | `{dataset, preview}`（等同 `open {demo:true}`） |
| `list` | `{}` | `{datasets:[摘要],active}` |
| `select`/`close`/`reset` | `{dataset?}` | `{active,datasets}` / `{closed,datasets}` / `{reset:true}` |
| `preview` | `{dataset}` | `{genes:[前50],spots:[前10],n_spots,n_genes,coords_meta:{x_range,y_range,size},columns:[{key,kind,n_unique}],steps,current,cluster_keys:[],spatial_mode,image}` |
| `samples` | `{dataset,op:"list|update|auto_group|import_map|join_obs|derive|set_group|remove",...}` | `{dataset, sample_info:[{id,name,group,kind,n_spots}], groups, group_column, mapping?, matched?, total?}` |
| `run` | `{dataset,step,params}` | `{step,summary,steps,current,dataset:<摘要>,tables?,figures?,dataset_notes?}` |
| `plot` | `{dataset,kind,params}` | `{file,width,height,bytes,kind}` |
| `analyses` | `{dataset}` | `{analyses:[{id,title,group,method,needs,params,available,reason,native}],groups:[],diagnostics:{},columns:[]}` |
| `analysis` | `{dataset,id,params}` | 结果信封 + `{analysis_id,ms}` |
| `export` | `{dataset,what,path?}` | `{file,bytes,what}` |
| `audit` | `{dataset}` | `{score,level,checks:[{id,level,title,detail,hint}]}` |
| `suggest` | `{dataset?}` | `{suggestions:[{id,title,detail,action,step?,params?}]}` |
| `recipe` | `{id?}` | 无 id：`{recipes:[{id,title,detail,steps:[{step,params}]}]}` |
| `formats` | `{}` | `{formats:[{id,label,hint}],ids:[]}` |
| `steps` | `{}` | `{steps:[目录]}` |
| `geo_search` | §7.1 | §7.1 |
| `geo_detail` | `{gse}` | §7.2 |
| `geo_download` | `{gse,files,dest_dir,overwrite?}` | §7.3 |
| `refs` | `{path?}` | 无 path：`{refs:[{path,name,bytes,n_obs?,n_vars?}],available}`；有 path：`{path,n_obs,n_vars,columns:[],cell_types:[]}` |
| `ref_set` | `{dataset,ref_path?,cell_type_col?}` | `{dataset,ref:{path,cell_type_col,cell_types}}` |
| `shutdown` | `{}` | `{bye:true}` |

### 2.3 `ctx`（传给所有步骤/分析模块）

```python
ctx = {
  "dataset_id": "ds1", "name": "GSE...",
  "genes": ["A1BG", ...],                     # list[str] 长度 = X.shape[1]
  "spots": ["AAACAACGAATAGTTC-1", ...],       # list[str] 长度 = X.shape[0]
  "X": np.ndarray float32,                    # (n_spots, n_genes)，已标准化则是 lognorm
  "X_kind": "counts|lognorm|unknown",
  "counts": np.ndarray|None,
  "coords": np.ndarray float32,               # (n_spots, 2)，图像像素坐标（Visium 已按 scalefactors 缩放到 hires）
  "coords_kind": "visium_hires|visium_lowres|pixel|grid|generic",
  "image": {"path": str|None, "extent": [x0,x1,y0,y1]|None, "scale": float|None},
  "spatial_mode": "visium|grid|generic|stereo|slideseq",
  "neighbor_radius": 150.0,                   # 空间邻接半径（像素/单位），由 import 时估计并写进 dataset
  "cluster_key": "clusters",                  # 当前聚类列
  "cluster": np.ndarray|None,                 # 当前聚类标签（与 spots 等长）
  "sample_info": pd.DataFrame,                # index=spots 或样本级 DataFrame（见下）
  "group_key": "group"|None, "groups": ["a","b"]|None,
  "species": "human|mouse|unknown",
  "obsm": {},                                 # 派生坐标 {"pca","tsne","umap","spatial"}
  "results": {},                              # 跨分析缓存 {"svg":..,"domains":..,"deconv":..,"de":..}
  "ref": {"path","cell_type_col","cell_types","X","genes"} | None,   # 反卷积参考（懒加载）
  "rng": np.random.Generator,                 # seed=0
  "workdir": str, "outdir": str,              # outdir=<workdir>/analyses/<ds>/<analysis>
  "engine": {"python":"..","mode":"builtin"},
  "notes": [], "warnings": [],
}
```

- **只读**：模块不得原地改 ctx 的 X/coords；需要改就复制。写回由 registry 消费返回键完成（§2.5）。
- 样本级信息：空间数据集没有"bulk 样本"概念，`sample_info` 可以是 spot 级（有 `slice`/`group` 列）或数据集级 1 行。凡按分组比较的分析，用 **spot 级 `group_key`**（如切片分组），没有就要求 prereq。
- 大数据上限：n_spots ≤ 200000、n_genes ≤ 60000；超过 60000 spots 时聚类/嵌入/绘图内部下采样并在 warnings 说明（SVG 仍全量）。

### 2.4 `dataset` 摘要

```json
{"id":"ds1","name":"Visium-human-breast","path":"/abs/outs","source_format":"visium",
 "n_spots":3800,"n_genes":33000,"species":"human",
 "geometry":{"mode":"visium","coords_kind":"visium_hires","x_range":[0,1200],"y_range":[0,1400],
             "image":"/abs/tissue_hires_image.png","scale":0.15,"neighbor_radius":150},
 "sample_columns":[{"key":"group","kind":"categorical","n_unique":2},{"key":"slice","kind":"categorical"}],
 "steps":["qc","filter","normalize"],"current":"normalize","cluster_keys":["clusters","domains"],
 "ref":{"path":"...","cell_type_col":"cell_type","cell_types":["T cell"]} ,
 "warnings":[],"notes":[],"history":[{"at":"ISO","step":"qc","params":{},"ok":true,"ms":812,"summary":".."}],
 "analyses_run":[{"at":"ISO","id":"morans_svg","ok":true,"ms":120}]}
```
`history` 最多 50 条，返回 UI 时截断到 20。

### 2.5 步骤 / 分析模块契约

模块文件导出 `STEPS` / `ANALYSES`（结构同 dsh-geobulk）：

```python
STEPS = [{
  "id":"qc","title":"空间质控","icon":"qc","requires":[],
  "params":[{"key":"min_counts","label":"最少 counts","type":"number","default":None}, ...],
  "run": callable,          # run(ctx, params, outdir) -> result 信封
}]
ANALYSES = [{
  "id":"morans_svg","title":"空间自相关（Moran's I）","group":"空间统计",
  "method":"Moran's I + 置换检验（Moran 1950；squidpy spatial_autocorr）",
  "needs":["cluster"],      # 语义前置：coords/x/y | normalize | cluster | domains | groups | ref | image
  "params":[...], "run": callable,
}]
```


**图的归属（重要）**：每个模块除 `STEPS`/`ANALYSES` 外还可以导出
```python
PLOTS = {"svg_map": fn, "domain_map": fn}   # kind -> fn(ctx, params, outdir) -> figure dict
```
registry 会把所有模块的 `PLOTS` 合并成全局 `kind → fn` 表，`plot` 命令按 kind 分发（每个 kind 只能有一个 owner，
归属按 DESIGN.md §5 的分组：基础图在 `<prefix>_prep.py`/`<prefix>_plots.py`，其余各自模块）。

**result 信封（固定键名）**：

```python
{"summary":{"title":"...","bullets":["..."],"metrics":{}},
 "tables":[{"id":"svg","title":"空间高变基因","columns":[{"key":"gene","label":"基因","type":"string"},
            {"key":"moran_i","label":"Moran's I","type":"number"}],"rows":[...500 行以内...],"note":"可选"}],
 "figures":[{"id":"svg_map","title":"Top SVG","file":"/abs/x.png","width":760,"height":520}],
 "sample_columns":[{"key":"domains","kind":"categorical"}],       # 可选：写回 spot 级列
 "sample_values":{"domains":[...]} ,                              # 与 ctx.spots 等长
 "data_out":[{"kind":"svg_csv","path":"/abs/svg.csv","label":"SVG 表"}],
 "warnings":[],"notes":[]}
```

- 表格完整写文件到 `data_out`；返回 JSON 的 `rows` 每表 ≤ 500 行（超出 `"truncated":true`）。
- 图：PNG，白底，dpi≈140，**图上文字一律英文**（避免缺 CJK 字体），存 `outdir`，返回绝对路径。
- 缺前置 → `EngineError("prerequisite","中文原因")`；不得编造结果。
- 确定性：同数据同参数同结果（`rng = default_rng(0)`）。

registry 额外消费的顶层键（**不发给前端**）：`"replace"`（`{"X":..,"genes":..,"spots":..,"counts":..,"value_kind":..}`）、`"store"`（写入 `dataset.results[id]`）、`"cluster_store"`（`{"key":"clusters","labels":[...]}`）、`"coords_store"`。

### 2.6 `index.js` / `impl.mjs` / 路由（Lead 固定，模块无关）

- 前缀 `/spatial`，只接受回环地址；`GET /figure?id=<basename>.png`（`[A-Za-z0-9._-]+\.png`）。
- 所有 op 函数在 impl.mjs；HTTP 路由与 Agent 工具调同一函数（一个操作两个调用方）。
- 长任务（GEO 下载、反卷积、聚类）走 `POST` 立即返回 `{job}`，前端轮询 `GET /jobs?id=`；工具内部等待完成（带超时）。

---

## 3. 步骤清单（顺序即左栏卡片顺序）

| step | 标题 | requires | 参数（默认） | 输出 |
|---|---|---|---|---|
| `qc` | 空间质控 | — | `min_counts=null, min_genes=null, max_pct_mt=20, group_col=null` | spot 级 counts/genes/pct_mt 表；`qc_spatial` + `qc_violin` |
| `filter` | 过滤 | qc | `min_counts=500, min_genes=200, max_pct_mt=20, min_cells=3` | spot/基因数变化；notes 写明阈值 |
| `normalize` | 标准化 | filter | `target_sum=10000, log=true, method=lognorm` | `X_kind=lognorm`；分布图（可选） |
| `hvg` | 高变基因 | normalize | `n_top_genes=2000, flavor=seurat` | HVG 表/图 |
| `pca` | PCA | hvg | `n_pcs=30, scale=true` | 碎石图；`obsm["pca"]` |
| `neighbors` | 邻接图 | pca | `n_neighbors=15, metric=cosine, spatial_weight=0.3, spatial_method=grid` | 图统计（连通分量/平均度） |
| `cluster` | 聚类 | neighbors | `resolution=1.0, method=louvain` | `clusters` 列 + `cluster_bar` + `spatial_scatter` |
| `embed` | 二维嵌入 | pca | `layout=tsne, perplexity=30` | `embedding` 图；`obsm["tsne"]` |
| `svg` | 空间高变基因 | cluster | `n_top=30, method=moran, permutations=0` | SVG 表 + `svg_map` + `moran_plot` |
| `domains` | 空间域 | cluster | `n_domains=6, method=mrf, smooth=0.5, iterations=20` | `domains` 列 + `domain_map` + 域组成表 |
| `deconv` | 反卷积 | normalize | `ref="", method=nnls, cell_type_col="", min_genes=50` | 比例表 + `deconv_bar` + `deconv_spatial` |
| `annotate` | 区域注释 | domains | `mapping={}, key="domain_anno"` | 写回注释列 |

- 重跑某步清空其下游（`history` 保留）。
- `svg` 在 `n_top≤0` 或坐标缺失时给 `prerequisite`。

---

## 4. 分析清单（≥30 项；`needs` 只写真实依赖）

> 方法列必须诚实：我们自己实现的一律写"内置实现（依据 X）"，不要写成官方包。

**空间统计（group=空间统计）**
| id | 标题 | method | needs | 主要参数 |
|---|---|---|---|---|
| `morans_svg` | 空间自相关（Moran's I） | Moran's I + 正态近似/置换（Moran 1950；squidpy `spatial_autocorr`） | cluster | `n_top=20, genes=[], permutations=100, min_expr=0.1` |
| `gearys_svg` | 空间自相关（Geary's C） | Geary's C（Geary 1954） | cluster | `n_top=20, genes=[]` |
| `svg_permutation` | 空间置换检验 | 标签置换经验 p（SPARK-X 思想，Sun 2021） | cluster | `genes=[], permutations=200, n_top=20` |
| `cross_moran` | 双变量空间自相关 | 双变量 Moran's I（Wartenberg 1985） | cluster | `gene_a, gene_b` |
| `hotspot` | 热点检测（Getis-Ord Gi*） | Gi* 统计量（Getis & Ord 1992） | cluster | `gene, k=8, z_thresh=1.96` |
| `spatial_gradient` | 空间梯度 | 沿轴分箱 + 趋势/斜率（tissue-gradient 常用做法） | cluster | `gene, axis=xy, angle=0, bins=20` |
| `spatial_corr_pairs` | 空间共表达基因对 | Pearson/Spearman + 邻域一致性 | cluster | `genes=[], top_pairs=30` |
| `ripley_k` | Ripley's K/L 函数 | Ripley 1977；完全空间随机(CSR) 包络 | cluster | `cluster_key=clusters, radii=20, n_sim=19` |
| `gene_autocorr_scan` | 全基因空间自相关扫描 | 向量化 Moran's I（对角化近似） | normalize | `max_genes=2000, min_expr=0.1` |

**邻域与域（group=邻域与域）**
| id | 标题 | method | needs | 主要参数 |
|---|---|---|---|---|
| `nhood_enrichment` | 邻域富集 | 类别对计数 + 置换 z-score（squidpy `nhood_enrichment`；Schürch 2020 邻域思想） | cluster | `cluster_key=clusters, n_perms=100, n_neighbors=6` |
| `co_occurrence` | 共现分析 | 环形分箱条件概率（squidpy `co_occurrence`） | cluster | `cluster_key=clusters, n_rings=5, interval=null` |
| `domains_mrf` | 空间域（MRF 平滑） | PCA + kNN + 马尔可夫随机场迭代（BayesSpace 2021 思想） | cluster | `n_domains=6, smooth=0.5, iterations=20` |
| `domain_markers` | 区域标记基因 | Wilcoxon 秩和 + BH | cluster | `domain=1, top_n=20` |
| `domain_composition` | 区域组成 | 各分组/切片在各域的 spot 比例 | cluster | `group_col=null` |
| `boundary_score` | 边界强度 | 邻接对中跨域比例 + 环形梯度 | cluster | `domain_a=1, domain_b=2, n_rings=3` |
| `distance_to_domain` | 到区域距离 | 欧氏距离场 + 分箱表达 | cluster | `domain=1, genes=[], bins=10` |
| `niche_cluster` | 细胞邻域（niche） | 邻域组成向量 + k-means（Schürch 2020；Bhate 2022） | cluster | `k=8, n_niches=5` |
| `tissue_diversity` | 邻域多样性 | Shannon/Simpson 多样性（邻域组成） | cluster | `k=6` |

**反卷积（group=反卷积；needs 含 ref）**
| id | 标题 | method | needs | 参数 |
|---|---|---|---|---|
| `deconv_nnls` | 反卷积（NNLS） | 非负最小二乘（Tsoucas 2019 DWLS 基础；SPOTlight 2021） | ref | `ref, cell_type_col, min_genes=50` |
| `deconv_dwls` | 反卷积（加权 NNLS） | 阻尼加权最小二乘（Tsoucas 2019） | ref | `ref, cell_type_col, weight=expression` |
| `deconv_spatial_map` | 细胞类型空间分布 | 比例映射到坐标 | ref | `cell_type, point_size=8` |
| `deconv_compare` | 反卷积方法一致性 | Pearson/Spearman 相关矩阵 | ref | — |
| `deconv_roi` | 区域细胞组成差异 | 域/分组内比例 + Kruskal-Wallis | ref, cluster | `group_col=null, cell_type=null` |
| `deconv_marker_corr` | 比例-标记基因一致性 | 比例与标志基因表达相关 | ref | `cell_type, genes=[]` |

**空间通讯（group=空间通讯）**
| id | 标题 | method | needs | 参数 |
|---|---|---|---|---|
| `ligrec_spatial` | 空间配体-受体 | LR 共表达 + 空间邻近（CellPhoneDB 2020；CellChat 2021；squidpy `ligrec`） | cluster | `cluster_key=clusters, n_neighbors=6, min_frac=0.1, top_n=30` |
| `lr_map` | 通讯强度空间图 | 邻近簇间 LR 乘积映射 | cluster | `ligand, receptor` |

**区域差异（group=区域差异）**
| id | 标题 | method | needs | 参数 |
|---|---|---|---|---|
| `de_domain` | 域 vs 其余 差异 | Wilcoxon + BH | cluster | `domain=1, top_n=200` |
| `de_regions` | 两区域差异 | Wilcoxon + BH + logFC | cluster | `domain_a, domain_b, top_n=200` |
| `de_groups` | 按分组差异 | Wilcoxon/t + BH（需 spot 级分组） | groups | `group_col, case, control, top_n=200` |
| `gene_domain_box` | 区域表达箱线 | 分域箱线 + Kruskal | cluster | `gene, domain_col` |

**质控与整合（group=质控与整合）**
| id | 标题 | method | needs | 参数 |
|---|---|---|---|---|
| `qc_summary` | 质控汇总 | 汇总表 + 空间 QC 图 | coords | — |
| `outlier_spots` | 空间离群点 | 邻域表达 z-score（kNN 平滑残差） | cluster | `k=8, z=3` |
| `slice_align` | 多切片对齐诊断 | 各切片 PCA + 重心/分位数对齐诊断 | normalize | `slice_col` |
| `spatial_trajectory` | 空间轨迹 | 扩散伪时序（Haghverdi 2016 DPT 思想） | cluster | `start_gene="", genes=[]` |

**验收要求**：以上 id 全部出现在 `/spatial/analyses`；`needs` 缺前置时返回 `prerequisite` 且理由中文可读。

---

## 5. 图（`kind` 白名单，与 impl.mjs 一致）

`qc_spatial` `qc_violin` `spatial_scatter` `embedding` `cluster_bar` `gene_spatial`
`svg_map` `moran_plot` `domain_map` `deconv_spatial` `deconv_bar` `nhood_enrich`
`cooccurrence` `ripley` `pca_variance` `heatmap` `dotplot` `gene_box`
`de_volcano` `de_heatmap`

- 所有图都接受 `{color_by?:"clusters|domains|group|cell_type|gene", gene?, point_size=6, max_points=20000, dpi=140}`。
- 空间图：坐标 y 轴翻转（图像坐标原点在左上），上色用 `sp_common.PALETTE`，无边框、等比、留白。
- 允许叠加 H&E 底图（有 `ctx.image.path` 时 `params.image=true`）。

## 6. 数据格式（`sp_io.py`）

| data_format | 识别 | 规则 |
|---|---|---|
| `visium` | 目录含 `spatial/tissue_positions*.csv` + `filtered_feature_bc_matrix*` | 读 10x h5 或 mtx；坐标取 hires（有 scalefactors 时 ×`tissue_hires_scalef`）；读 `tissue_hires_image.png` |
| `visium_h5` | `filtered_feature_bc_matrix.h5` 单文件 | 要求同时给 `coords`（坐标 csv）或同目录 `spatial/` |
| `h5ad` | `.h5ad` | `obsm['spatial']` 或 `obs` 的 `x/y`/`pxl_row/pxl_col`；`var_names`=基因；保留 obs 列与 image |
| `stereo` | `*_raw_bin_count.tsv.gz` 或含 `geneID,x,y,count` | 长表 → 宽表 + 坐标；`bin1/bin20` 自动识别 |
| `slideseq` | `*_puck_*.txt` / `*_bead*.txt` | 矩阵 + beacon 坐标表（row/col × 分辨率） |
| `matrix_text` | `.csv/.tsv/.txt(.gz)` | 需配 `coords`（含 spot,x,y 或 barcode,x,y）；否则 `prerequisite` |
| `xlsx` | `.xlsx/.xls` | 同 matrix_text；`sheet`/`gene_col` 可指定 |
| `geo_series` | 含 `!series_matrix_table_begin` | 表段基因×样本；空间坐标缺失 → 只当表达矩阵，提示需坐标 |
| `archive` | `.zip/.tar/.tar.gz` | 解压到 `<workdir>/unpack/<uuid>`（防穿越），自动找 Visium/BIN/矩阵；多候选 → `ambiguous`+candidates |
| `rds`/`rdata` | R 桥 | 可读则读（Seurat/SCE→矩阵+obsm），不可用 → `unsupported` |
| `coords_csv` | — | 作为 `coords` 参数使用：列 `barcode/spot,x,y`（大小写不敏感） |

- 返回 `{"format","X","genes","spots","coords","sample_info","image","notes","warnings","meta"}`。
- 空间邻接半径估计：Visium=中心距（默认 100×scale 后 ~`hires` 100px），grid=1.2×中位最近邻距离，generic=1.2×中位最近邻距离；写 `dataset.geometry.neighbor_radius`。

## 7. GEO（`sp_geo.py`，只依赖标准库；可整体参考 `../dsh-geobulk/lib/bg_geo.py`）

### 7.1 `geo_search`
`{term, species:"human|mouse|all", study_type:"all|high throughput sequencing|array", min_samples, has_matrix:bool, platform:"all|visium|stereo|slideseq|10x", limit=20, page=0, sort:"relevance|date|samples"}`
→ `{total,page,items:[{gse,title,title_zh,summary,species,n_samples,platform,study_type,release,has_matrix,spatial_score:"high|mid|low",spatial_platform:["Visium"],tags:[]}],query?,query_zh?,query_en?,translated?,relaxed?,used_term?}`

- 中文查询用内置词典反查成英文（词典在 `sp_geo.py` 的 `ZH_DICT`，≥180 条，含空间技术词：Visium/Stereo-seq/Slide-seq/MERFISH/CODEX/空间转录组/切片/组织/肿瘤微环境等）。
- 少结果时逐级放宽（去泛词 → 前 3 词 → 前 2 词），返回 `relaxed=true`。
- `spatial_score`/`spatial_platform` 由标题/摘要关键词判定（"spatial transcriptomics"/"Visium"/"Stereo-seq"/"Slide-seq"/"MERFISH"/"spatial"）。
- 无结果 → `EngineError("not_found", 可读中文)`；网络失败重试 ≤2 次后 → `EngineError("unsupported","NCBI 暂时不可达：...")`。

### 7.2 `geo_detail`
`{gse}` → `{gse,title,title_zh,summary,design,species,n_samples,platforms,pubmed,release,contacts,has_matrix,has_clinical,spatial_platform,preview_files:[{name,size,kind}],files:[{name,size,url,kind:"matrix|suppl|raw"}],suggested,url}`
- 同时解析 `<GSE>/matrix/` 与 `<GSE>/suppl/` 目录；`suppl` 里 Visium 类通常是 `<GSE>_RAW.tar`（含 `spatial/` 图与坐标）。

### 7.3 `geo_download`
`{gse, files:[], dest_dir, overwrite=false}` → `{dir,gse,files:[{name,path,bytes,ok,skipped?}],ms,notes:[],warnings?}`
- `dest_dir` 必须绝对路径且不是 `/`；单文件 ≤8GB；`.part` + 完整性校验 + 原子替换；不做断点续传。
- 下载完由 Host 调 `open`（用户点"导入并清洗"）。

## 8. Agent 工具 `spatial`

单一工具，`action` 枚举 + 扁平参数（**参数名不能叫 `format`/`pattern`**）：

```
action: status | search | detail | download | open | demo | inspect | preview |
        run | plot | analyses | analysis | export | audit | suggest | recipe |
        samples | refs | ref_attach | formats | reset | close
```
关键参数：`term, gse, dest_dir, path, url, data_format, dataset, step, params, analysis, kind, what, op, ref_path, cell_type_col, out_path, recipe, actor`。
返回给模型：`summary.bullets` + 表前 20 行 + 图绝对路径 + 建议下一步。工具描述必须写清何时用 `search`/`open`/`run`/`analysis`。

## 9. demo 数据（`sp_state.py`）

- 合成组织：**2400 spots**，六边形网格排布成不规则组织轮廓，**5 个空间域**（分层分布，模拟肿瘤/间质/免疫浸润带），**1000 基因**（每域 40 个标记基因 logFC 1.2–2.5），几何 `neighbor_radius` 自动估计；`slice` 列 2 个切片（构成"多切片"场景），`group` 列 2 组。
- 生成 < 3 秒；`test/make_demo.py` 复用（`import sp_state; sp_state.make_demo()`）。
- 另置一个 400 细胞 × 800 基因 × 6 细胞类型的小 h5ad（`test/fixtures/demo_sc.h5ad`，由 `test/make_fixtures.py` 生成）供 `deconv` 演示；没有时 `deconv` 给 `prerequisite`。

## 10. 环境事实（禁止越界）

- 本机 `python3` = `/Library/Frameworks/Python.framework/Versions/3.14/bin/python3`（3.14.2）：
  numpy 2.5.1 / pandas 2.3.3 / scipy 1.17.1 / sklearn 1.9.0 / matplotlib 3.11.0 / anndata 0.12.19 / h5py 3.16 / networkx 3.6.1 / statsmodels 0.14.6 / openpyxl 3.1.5 / matplotlib only Agg / **无 scanpy、无 squidpy、无 umap-learn、无 igraph/leidenalg、无 pyarrow**。
  → 聚类自己实现（kNN + Louvain，用 networkx.community 或自写两阶段），嵌入用 sklearn TSNE / PCA；**禁止 pip/conda 安装**。
- R 4.4.2 可用（Seurat/Matrix/limma/preprocessCore/ComplexHeatmap/igraph 等），但**空间模块不依赖 R**（可选：`rds` 读取时用 Rscript 桥）。
- `MPLBACKEND=Agg` + `MPLCONFIGDIR=<workdir>/mpl-cache` 必须设置（`~/.matplotlib` 不可写）。
- 网络：NCBI E-utilities/FTP 可达（urllib 需用 `certifi` CA，见 bg_geo）；其它外网按不可达处理。
- 禁止：pip/conda 安装、改 profile、pnpm、写 `$DSH_HOME`、改 `../dsh-singlecell` 与 `../dsh-geobulk`。
- 提交前自检：`python3 -m py_compile lib/*.py`；`node --check client.js`；`node --input-type=module -e "import('./lib/impl.mjs')"`。

## 11. 测试要求（各 owner 自己跑通再交）

- **sp-core**：`test/py_io_test.py` + `test/make_fixtures.py`：造 Visium 目录（h5+坐标+scalefactors）、BIN 长表、h5ad、矩阵+坐标、zip；断言 `genes/spots/X/coords` 形状与 truth 一致；`inspect` 不加载全表；坏文件给可读错误不崩。`test/py_geo_test.py`：mock `_urlopen` 返回固定 JSON/HTML，断言 search/detail/download 解析与重试、中文翻译、放宽策略。
- **sp-base**：`test/py_prep_test.py`：demo 上 qc→filter→normalize→hvg→pca→neighbors→cluster→embed 全通；`X_kind` 演进正确；聚类在合成域上 ARI ≥0.6；缺前置给 `prerequisite`。
- **sp-stat**：`test/py_spatial_test.py`：注入的 40×5 域标记基因 Moran's I 显著高于随机基因；置换 p 值 ∈[0,1]；Ripley K 单调；Gi* 热点落在注入的高表达块。
- **sp-adv**：`test/py_adv_test.py`：MRF 域与注入域 ARI ≥0.5；反卷积用 `demo_sc.h5ad` 比例和 ≈1、与注入真值相关 >0.5；`ligrec` 返回 ≥5 对且 p 值合法。
- **sp-ui**：`test/client-logic.test.mjs`：纯函数（格式化/坐标变换/搜索 payload/API 客户端/建议指令/任务轮询）+ 静态约束。
- 测试不依赖 npm 包（node 用内置 `node:test`；python 用已装包）。

## 12. 交付验收（Lead）

1. `plugin_manager install_bundle` 成功；`GET /spatial/health` 200 且 `engine.ready=true`。
2. demo 全流程：qc→…→svg→domains→deconv 每步成功且有 summary/图。
3. 真实 GEO 检索一次（关键词"spatial transcriptomics breast cancer"）+ 真实下载一个 suppl 文件成功（网络受限时用 mock 证明链路）。
4. 导入 10x Visium `outs/` 目录（可用测试 fixture 造）→ 出图正确。
5. UI：侧边栏出现「空间工作台」，整页可打开，主要 tab 不白屏，GEO/流程卡/分析可交互。
6. 反卷积能用 dsh-singlecell 的数据集（`ref_attach` 拉取）。
```
