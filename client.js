/**
 * dsh-spatial — 浏览器端整页 UI（Client 半边）
 *
 * 形态照 dsh-geobulk：window.__ModuleLoader__.load({ id:'dsh-spatial', factory(require){...} })，
 * 只 require('react')，不碰任何 @deepseek-ai Client 包（宿主版本一变就会炸）。
 *
 * 数据面：全部走同源相对路径 /spatial/*（见 DESIGN.md §2.6 / §7）。
 *
 * 文件分两层：
 *   1) 纯函数层（本段之前到 SP_TEST）：STEPS 元数据 / 数值与时间格式化 / 分组取色 /
 *      空间坐标→像素变换 / GEO 搜索 payload / 建议→中文指令 / 可注入 fetch 的 API 客户端 /
 *      异步任务轮询状态机 / 服务端响应归一化。不碰 DOM、不 require react；
 *      node 环境（typeof window === 'undefined'）挂到 globalThis.__SP_TEST__
 *      供 test/client-logic.test.mjs 直接 import。
 *   2) 浏览器层（factory 内）：React 组件 + 样式 + 槽注册。
 *
 * 约定：类名一律 sp- 前缀；颜色只用 --dsw-alias-* 主题变量（见 CSS 段），
 * 数据色（分组/簇/连续值）用 JS 常量 PALETTE / CONTINUOUS_STOPS。
 */
;(function () {
  'use strict'

  var NS = 'dsh-spatial'
  var API = '/spatial'
  var POLL_MS = 1500
  var POLL_TIMEOUT_MS = 30 * 60 * 1000
  var DEFAULT_API_TIMEOUT_MS = 60000
  var HEAVY_API_TIMEOUT_MS = 600000
  var MAX_HISTORY_ROWS = 20
  var MAX_TABLE_RENDER_ROWS = 200
  var DEFAULT_GEO_DIR = '~/.dsh/spatial/GEO'

  /** 分组/簇/细胞类型调色板（数据色，不参与主题）。 */
  var PALETTE = [
    '#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f', '#edc948',
    '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac', '#1f77b4', '#ff7f0e',
    '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
    '#bcbd22', '#17becf',
  ]

  /** 连续值色带（数据色）：低→高，用于基因表达/比例映射。 */
  var CONTINUOUS_STOPS = ['#3b4cc0', '#6a9fdb', '#9ec6e8', '#e8d4a8', '#f0a35e', '#c0392b']

  /** 自检等级样式（由主题色派生）。 */
  var AUDIT_LEVEL_STYLE = {
    ok: { tone: 'ok', labelKey: 'auditOk' },
    warn: { tone: 'warn', labelKey: 'auditWarn' },
    fail: { tone: 'fail', labelKey: 'auditFail' },
  }

  /** 空间图 kind 白名单（与 DESIGN.md §5 逐字一致）。 */
  var PLOT_KINDS = [
    'qc_spatial', 'qc_violin', 'spatial_scatter', 'embedding', 'cluster_bar', 'gene_spatial',
    'svg_map', 'moran_plot', 'domain_map', 'deconv_spatial', 'deconv_bar', 'nhood_enrich',
    'cooccurrence', 'ripley', 'pca_variance', 'heatmap', 'dotplot', 'gene_box',
    'de_volcano', 'de_heatmap',
  ]

  /** 记录 tab 的导出按钮组（what 与 DESIGN.md §2.2 export / §5 一致）。 */
  var EXPORT_KINDS = [
    { id: 'expr_csv', labelKey: 'exportExprCsv', hint: 'CSV' },
    { id: 'de_csv', labelKey: 'exportDeCsv', hint: 'CSV' },
    { id: 'svg_csv', labelKey: 'exportSvgCsv', hint: 'CSV' },
    { id: 'domain_csv', labelKey: 'exportDomainCsv', hint: 'CSV' },
    { id: 'deconv_csv', labelKey: 'exportDeconvCsv', hint: 'CSV' },
    { id: 'sample_info_csv', labelKey: 'exportSampleCsv', hint: 'CSV' },
    { id: 'summary_json', labelKey: 'exportSummary', hint: 'JSON' },
    { id: 'protocol', labelKey: 'exportProtocol', hint: 'JSON' },
    { id: 'methods', labelKey: 'exportMethods', hint: 'Markdown' },
    { id: 'report', labelKey: 'exportReport', hint: 'Markdown' },
    { id: 'script', labelKey: 'exportScript', hint: 'Python' },
  ]

  /** 记录 tab 上方快速导出（只放最常用的三个）。 */
  var QUICK_EXPORTS = ['expr_csv', 'summary_json', 'report']

  /** 空间图着色维度（DESIGN.md §5 color_by 白名单）。 */
  var COLOR_BY = [
    { value: 'clusters', labelKey: 'colorClusters' },
    { value: 'domains', labelKey: 'colorDomains' },
    { value: 'group', labelKey: 'colorGroup' },
    { value: 'cell_type', labelKey: 'colorCellType' },
    { value: 'gene', labelKey: 'colorGene' },
  ]

  /** GEO 空间平台筛选（DESIGN.md §7.1 platform 取值）。 */
  var GEO_PLATFORMS = [
    { value: 'all', labelKey: 'geoPlatformAll' },
    { value: 'visium', label: 'Visium' },
    { value: 'stereo', label: 'Stereo-seq' },
    { value: 'slideseq', label: 'Slide-seq' },
    { value: 'merfish', label: 'MERFISH' },
    { value: '10x', label: '10x Genomics' },
  ]

  /** 空间数据格式（服务端 /formats 不可用时的静态兜底，与 DESIGN.md §6 一致）。 */
  var DATA_FORMATS = [
    'auto', 'visium', 'visium_h5', 'h5ad', 'stereo', 'slideseq',
    'matrix_text', 'xlsx', 'geo_series', 'archive', 'rds', 'rdata',
  ]

  // ==========================================================================
  // 流程步骤：顺序即左栏卡片顺序；requires 与 DESIGN.md §3 逐字一致。
  // params 用通用渲染（number/boolean/select/column/text/json）。
  // ==========================================================================
  var STEPS = [
    {
      id: 'qc',
      title: '空间质控',
      icon: 'qc',
      requires: [],
      params: [
        { key: 'min_counts', label: '最少 counts', type: 'number', default: null, min: 0, step: 100 },
        { key: 'min_genes', label: '最少基因数', type: 'number', default: null, min: 0, step: 50 },
        { key: 'max_pct_mt', label: '线粒体比例上限（%）', type: 'number', default: 20, min: 0, max: 100, step: 1 },
        { key: 'group_col', label: '分组列', type: 'column', default: null },
      ],
    },
    {
      id: 'filter',
      title: '过滤',
      icon: 'filter',
      requires: ['qc'],
      params: [
        { key: 'min_counts', label: '最少 counts', type: 'number', default: 500, min: 0, step: 100 },
        { key: 'min_genes', label: '最少基因数', type: 'number', default: 200, min: 0, step: 50 },
        { key: 'max_pct_mt', label: '线粒体比例上限（%）', type: 'number', default: 20, min: 0, max: 100, step: 1 },
        { key: 'min_cells', label: '基因最少 spot 数', type: 'number', default: 3, min: 0, step: 1 },
      ],
    },
    {
      id: 'normalize',
      title: '标准化',
      icon: 'normalize',
      requires: ['filter'],
      params: [
        { key: 'target_sum', label: '目标总量', type: 'number', default: 10000, min: 100, step: 1000 },
        { key: 'log', label: 'log1p 变换', type: 'boolean', default: true },
        {
          key: 'method', label: '方法', type: 'select', default: 'lognorm',
          options: [{ value: 'lognorm', label: 'LogNormalize' }, { value: 'cpm', label: 'CPM' }, { value: 'none', label: '不标准化' }],
        },
      ],
    },
    {
      id: 'hvg',
      title: '高变基因',
      icon: 'hvg',
      requires: ['normalize'],
      params: [
        { key: 'n_top_genes', label: '高变基因数', type: 'number', default: 2000, min: 50, step: 100 },
        {
          key: 'flavor', label: '算法', type: 'select', default: 'seurat',
          options: [{ value: 'seurat', label: 'Seurat' }, { value: 'cell_ranger', label: 'Cell Ranger' }, { value: 'vst', label: 'VST' }],
        },
      ],
    },
    {
      id: 'pca',
      title: 'PCA',
      icon: 'pca',
      requires: ['hvg'],
      params: [
        { key: 'n_pcs', label: '主成分数', type: 'number', default: 30, min: 2, max: 300, step: 5 },
        { key: 'scale', label: '标准化基因', type: 'boolean', default: true },
      ],
    },
    {
      id: 'neighbors',
      title: '邻接图',
      icon: 'neighbors',
      requires: ['pca'],
      params: [
        { key: 'n_neighbors', label: '近邻数', type: 'number', default: 15, min: 3, max: 200, step: 1 },
        {
          key: 'metric', label: '距离度量', type: 'select', default: 'cosine',
          options: [{ value: 'cosine', label: 'cosine' }, { value: 'euclidean', label: 'euclidean' }, { value: 'correlation', label: 'correlation' }],
        },
        { key: 'spatial_weight', label: '空间权重', type: 'number', default: 0.3, min: 0, max: 1, step: 0.05 },
        {
          key: 'spatial_method', label: '空间邻接', type: 'select', default: 'grid',
          options: [{ value: 'grid', label: '网格（Visium）' }, { value: 'knn', label: 'kNN' }, { value: 'radius', label: '半径' }],
        },
      ],
    },
    {
      id: 'cluster',
      title: '聚类',
      icon: 'cluster',
      requires: ['neighbors'],
      params: [
        { key: 'resolution', label: '分辨率', type: 'number', default: 1, min: 0.1, max: 5, step: 0.1 },
        {
          key: 'method', label: '算法', type: 'select', default: 'louvain',
          options: [{ value: 'louvain', label: 'Louvain' }, { value: 'kmeans', label: 'k-means' }, { value: 'gmm', label: 'GMM' }],
        },
      ],
    },
    {
      id: 'embed',
      title: '二维嵌入',
      icon: 'embed',
      requires: ['pca'],
      params: [
        {
          key: 'layout', label: '布局', type: 'select', default: 'tsne',
          options: [{ value: 'tsne', label: 't-SNE' }, { value: 'umap', label: 'UMAP' }, { value: 'pca', label: 'PCA' }],
        },
        { key: 'perplexity', label: 'perplexity', type: 'number', default: 30, min: 5, max: 100, step: 1 },
      ],
    },
    {
      id: 'svg',
      title: '空间高变基因',
      icon: 'svg',
      requires: ['cluster'],
      params: [
        { key: 'n_top', label: '返回基因数', type: 'number', default: 30, min: 1, max: 500, step: 5 },
        {
          key: 'method', label: '统计量', type: 'select', default: 'moran',
          options: [{ value: 'moran', label: 'Moran\'s I' }, { value: 'geary', label: 'Geary\'s C' }],
        },
        { key: 'permutations', label: '置换次数', type: 'number', default: 0, min: 0, max: 1000, step: 50 },
      ],
    },
    {
      id: 'domains',
      title: '空间域',
      icon: 'domains',
      requires: ['cluster'],
      params: [
        { key: 'n_domains', label: '域数量', type: 'number', default: 6, min: 2, max: 20, step: 1 },
        {
          key: 'method', label: '方法', type: 'select', default: 'mrf',
          options: [{ value: 'mrf', label: 'MRF 平滑' }, { value: 'kmeans', label: 'k-means' }],
        },
        { key: 'smooth', label: '平滑强度', type: 'number', default: 0.5, min: 0, max: 1, step: 0.05 },
        { key: 'iterations', label: '迭代次数', type: 'number', default: 20, min: 1, max: 200, step: 5 },
      ],
    },
    {
      id: 'deconv',
      title: '反卷积',
      icon: 'deconv',
      requires: ['normalize'],
      params: [
        { key: 'ref', label: '参考 h5ad 路径', type: 'text', default: '', placeholder: '/absolute/path/ref.h5ad' },
        {
          key: 'method', label: '方法', type: 'select', default: 'nnls',
          options: [{ value: 'nnls', label: 'NNLS（内置）' }, { value: 'dwls', label: 'DWLS（内置）' }, { value: 'rctd', label: 'RCTD（spacexr，需 R）' }],
        },
        { key: 'cell_type_col', label: '细胞类型列', type: 'text', default: '', placeholder: 'cell_type' },
        { key: 'min_genes', label: '最少共享基因', type: 'number', default: 50, min: 10, step: 10 },
      ],
    },
    {
      id: 'annotate',
      title: '区域注释',
      icon: 'annotate',
      requires: ['domains'],
      params: [
        { key: 'mapping', label: '域→注释映射（JSON）', type: 'json', default: '', placeholder: '{"1":"肿瘤上皮","2":"间质"}' },
        { key: 'key', label: '写回列名', type: 'text', default: 'domain_anno' },
      ],
    },
  ]

  var STEP_BY_ID = {}
  for (var stepIndex = 0; stepIndex < STEPS.length; stepIndex++) STEP_BY_ID[STEPS[stepIndex].id] = STEPS[stepIndex]

  function stepById(id) {
    var key = typeof id === 'string' ? id : ''
    return STEP_BY_ID[key] || null
  }

  /** 中栏 tab 顺序（概览 / 空间图 / 表达 / 差异 / 空间统计 / 域与邻域 / 反卷积 / 通讯 / 记录 / 日志）。 */
  var TAB_KEYS = ['overview', 'spatial', 'expr', 'de', 'stats', 'domains', 'deconv', 'comm', 'records', 'log']

  var TAB_META = {
    overview: { titleKey: 'tabOverview' },
    spatial: { titleKey: 'tabSpatial' },
    expr: { titleKey: 'tabExpr' },
    de: { titleKey: 'tabDe' },
    stats: { titleKey: 'tabStats' },
    domains: { titleKey: 'tabDomains' },
    deconv: { titleKey: 'tabDeconv' },
    comm: { titleKey: 'tabComm' },
    records: { titleKey: 'tabRecords' },
    log: { titleKey: 'tabLog' },
  }

  /** 步骤跑完自动切到对应 tab（没有映射的留在原地）。 */
  var TAB_FOR_STEP = {
    qc: 'overview',
    filter: 'overview',
    normalize: 'overview',
    hvg: 'overview',
    pca: 'overview',
    neighbors: 'overview',
    cluster: 'spatial',
    embed: 'overview',
    svg: 'stats',
    domains: 'domains',
    deconv: 'deconv',
    annotate: 'domains',
  }

  /**
   * 分析分组 → tab。服务端 GROUP 文案取自 DESIGN.md §4（空间统计/邻域与域/反卷积/空间通讯/区域差异/质控与整合）；
   * 未映射的分组统一落到「空间统计」tab，保证 34 项分析全部可达。
   */
  var ANALYSIS_TAB = {
    '空间统计': 'stats',
    '邻域与域': 'domains',
    '反卷积': 'deconv',
    '空间通讯': 'comm',
    '区域差异': 'de',
    '质控与整合': 'stats',
  }

  var ANALYSIS_TAB_TITLES = {
    stats: 'tabStats',
    domains: 'tabDomains',
    deconv: 'tabDeconv',
    comm: 'tabComm',
    de: 'tabDe',
  }

  // ==========================================================================
  // 文案字典（zh 完整；en 与 zh 同 key 集）。apply() 里 locale.register 注册；
  // locale/zh.json 与 locale/en.json 由同一份词表生成，保证三处一致。
  // 参数 label 等来自服务端 schema 的文案属于数据，不走字典。
  // ==========================================================================
  var DICT_ZH = {
    nav: '空间工作台',
    appTitle: '空间转录组工作台',
    appSubtitle: 'GEO 检索 · 坐标导入 · 空间聚类 · SVG · 空间域 · 反卷积',
    refresh: '刷新',
    refreshHint: '重新读取工作台状态',
    retry: '重试',
    loading: '加载中…',
    loadingData: '正在读取工作台…',
    noData: '尚未载入数据集',
    errorTitle: '出错了',
    crashTitle: '这个面板渲染失败',
    crashHint: '其它面板仍可使用，可以点击重试恢复。',
    close: '关闭',
    cancel: '取消',
    copy: '复制',
    copied: '已复制',
    copyFailed: '复制失败（浏览器不允许写入剪贴板）',
    search: '搜索',
    all: '全部',
    none: '无',
    yes: '是',
    no: '否',
    unknown: '未知',
    dataset: '数据集',
    spots: 'spots',
    samples: '样本',
    genes: '基因',
    slices: '切片',
    clusters: '聚类',
    domainsLabel: '空间域',
    species: '物种',
    spatialMode: '空间模式',
    engine: '引擎',
    warnings: '警告',
    notes: '说明',
    current: '当前',
    completed: '已完成',
    pending: '待运行',
    running: '运行中',
    failed: '失败',
    plot: '出图',
    plotting: '出图中…',
    noFigure: '暂无图',
    figureFail: '图片加载失败（文件可能已被清理）',
    exportBtn: '导出',
    exported: '已导出',
    exportFail: '导出失败',
    copyPath: '复制路径',
    pathLabel: '文件路径',
    tableRows: '共 {v} 行',
    tableTruncated: '仅展示前 {v} 行',
    emptyTable: '空表',
    total: '总计',
    selectAll: '全选',
    clearSelection: '清空',
    of: '/',
    selectDataset: '切换数据集',
    closeDataset: '关闭数据集',
    workdir: '工作目录',
    engineReady: '引擎就绪',
    engineBooting: '引擎启动中',
    engineBad: '引擎不可用',
    engineGuideTitle: 'Python 引擎尚未就绪',
    engineGuideBody: '工作台依赖本机 Python 环境（numpy / pandas / scipy / sklearn / matplotlib / anndata / h5py / networkx）。整页仍可浏览，但导入、计算与出图需要引擎先就绪。',
    engineGuideSteps: '排查顺序：① 查看 engine.log 的报错；② 确认插件设置里的 python3 路径可执行；③ 重启插件或宿主后点「重试」。',
    pythonLabel: 'Python',
    engineInfo: '引擎信息',
    statusOk: '正常',
    statusWarn: '注意',
    statusFail: '异常',
    pipelineTitle: '分析流程',
    pipelineHint: '按顺序运行；重跑某步会清空其下游结果',
    datasetsTitle: '工作台',
    datasetNoDataset: '还没有数据集',
    datasetOpenHint: '搜 GEO、导入 Visium outs/ 目录或 h5ad，或载入内置示例立刻试一遍。',
    datasetStepsDone: '已完成 {v} 步',
    datasetCurrentStep: '当前 {v}',
    datasetSpots: 'spots {v}',
    datasetGenes: '基因 {v}',
    datasetSlicesCount: '切片 {v}',
    datasetClustersCount: '聚类 {v}',
    datasetRadius: '邻域半径 {v}',
    datasetRef: '参考 {v}',
    stepStatusDone: '已完成',
    stepStatusRunning: '运行中',
    stepStatusReady: '可运行',
    stepStatusBlocked: '待前置',
    stepStatusFailed: '失败',
    stepRun: '运行',
    stepRerun: '重跑',
    stepRunning: '运行中…',
    stepParams: '参数',
    stepNeedPrereq: '需要先完成：',
    stepSummaryEmpty: '—',
    runStepDone: '步骤完成',
    analysisCatalog: '分析目录',
    analysisHint: '从 /spatial/analyses 动态读取的全部分析',
    analysisTitle: '空间分析目录',
    analysisSubtitle: '按分组浏览；参数表单由服务端 schema 生成',
    analysisClose: '收起',
    analysisRun: '运行分析',
    analysisRunning: '运行中…',
    analysisNoParams: '该分析没有可调参数',
    analysisNeedDataset: '请先载入数据集',
    analysisUnavailable: '不可用',
    analysisNeeds: '前置',
    analysisEmpty: '这个分组下还没有可用分析（引擎未就绪或缺少前置步骤）',
    analysisGroupOther: '其它',
    analysisSelectHint: '左侧选一个分析，右侧填参数后运行',
    analysisResult: '分析结果',
    analysisMethod: '方法',
    analysisRetry: '重试',
    tabOverview: '概览',
    tabSpatial: '空间图',
    tabExpr: '表达',
    tabDe: '差异',
    tabStats: '空间统计',
    tabDomains: '域与邻域',
    tabDeconv: '反卷积',
    tabComm: '通讯',
    tabRecords: '记录',
    tabLog: '日志',
    overviewTitle: '数据集摘要',
    overviewHint: '先看数据形状与空间布局，再决定走哪条分析路线',
    metricSpots: 'Spots',
    metricGenes: '基因数',
    metricSpatialMode: '空间模式',
    metricSlices: '切片数',
    metricClusters: '当前聚类',
    metricDomains: '当前域',
    metricEngine: '引擎',
    clusterMap: '空间聚类图',
    pcaScree: 'PCA 碎石图',
    qcMap: '空间质控图',
    overviewClustersEmpty: '先运行聚类步骤，聚类图上会按 cluster 上色',
    miniMapTitle: '坐标缩略图',
    miniMapHint: '来自预览接口的坐标抽样；完整空间图见「空间图」tab',
    geneQuickTitle: '快速看基因',
    geneQuickHint: '输入基因名（回车搜索），出空间表达图与分组箱线',
    geneInputPlaceholder: '如 CDH1, VIM, MKI67',
    nextSteps: '建议下一步',
    noSuggestions: '暂无建议（引擎未就绪或已验证到终点）',
    spatialTitle: '空间图',
    spatialHint: '按坐标画 spot，可选聚类/域/分组/基因着色，可叠加 H&E 底图',
    spatialColorBy: '着色',
    spatialGeneSearch: '基因',
    spatialGenePh: '输入基因名，如 COL1A1',
    spatialOverlay: '叠加底图',
    spatialPointSize: '点大小',
    spatialMaxPoints: '最多点',
    spatialPlot: '画空间图',
    spatialGenePlot: '画基因空间图',
    spatialNeedGene: '请先输入基因名',
    spatialNeedCluster: '还没有聚类结果，可先在左栏运行「聚类」',
    spatialOverlayMissing: '该数据集没有底图，勾选无效',
    colorClusters: '聚类',
    colorDomains: '空间域',
    colorGroup: '分组',
    colorCellType: '细胞类型',
    colorGene: '基因表达',
    exprTitle: '表达',
    exprHint: '检索基因，看空间表达、分组箱线与热图',
    exprGeneSearch: '搜基因',
    exprAddGene: '加入',
    exprGeneNotInData: '数据里没有匹配的基因',
    exprSelected: '已选基因',
    exprEmpty: '还没选基因',
    exprPlotMap: '空间表达图',
    exprPlotBox: '分组箱线',
    exprPlotHeatmap: '表达热图',
    exprExport: '导出表达矩阵',
    exprExportHint: '表达矩阵以 CSV 导出（基因 × spot）',
    exprTableFallback: '热图与箱线按当前分组/聚类计算；服务端自行裁剪最多 200 行。',
    deTitle: '区域差异',
    deHint: '域 vs 其余 / 两域之间 / 按分组比较，Wilcoxon + BH',
    deEmpty: '还没有差异结果，选一个分析运行',
    deExport: '导出差异表',
    deVolcano: '火山图',
    deHeatmap: '差异热图',
    deQuickPlots: '快捷出图',
    deNeedCluster: '需要先有聚类/域结果',
    statsTitle: '空间统计',
    statsHint: 'Moran/Geary/置换/热点/梯度/共表达/Ripley/全基因扫描',
    statsEmpty: '还没有空间统计结果',
    domainsTitle: '域与邻域',
    domainsHint: 'MRF 空间域、区域标记、邻域富集、共现、边界与距离',
    domainsEmpty: '还没有空间域结果，先运行「空间域」步骤或 domains_mrf',
    deconvTitle: '反卷积',
    deconvHint: '用单细胞参考估计每个 spot 的细胞类型比例（NNLS/DWLS）',
    deconvEmpty: '还没有反卷积结果',
    deconvRefTitle: '参考数据',
    deconvRefHint: '选择参考 h5ad：本机路径，或直接从单细胞工作台拉取当前数据集',
    deconvRefPath: '参考路径',
    deconvRefPathPh: '/absolute/path/ref.h5ad',
    deconvRefAttach: '挂载参考',
    deconvRefAttaching: '挂载中…',
    deconvRefAttached: '已挂载参考',
    deconvRefNone: '未挂载参考',
    deconvRefCells: '细胞类型',
    deconvScList: '单细胞工作台数据集',
    deconvScPull: '拉取',
    deconvScUnavailable: '单细胞工作台不可用（未检测到 /singlecell）',
    deconvLocalRefs: '本机参考',
    deconvNoRefs: '没有找到可用的 h5ad 参考',
    deconvCellTypeCol: '细胞类型列',
    deconvRefSet: '设为反卷积参考',
    deconvRefDone: '参考已就绪，可以运行反卷积分析',
    deconvNeedRef: '反卷积分析需要先挂载参考（参考 h5ad）',
    commTitle: '空间通讯',
    commHint: '配体-受体共表达 + 空间邻近（CellPhoneDB / CellChat 思想）',
    commEmpty: '还没有通讯结果',
    recordsTitle: '记录与导出',
    recordsHint: '历史步骤、可复现脚本与 0–100 自检',
    recordsHistory: '步骤历史',
    recordsNoHistory: '还没有运行记录',
    recordsAudit: '确定性自检',
    recordsRunAudit: '运行自检',
    recordsAuditScore: '自检得分',
    recordsAuditEmpty: '点「运行自检」检查流程完整性、坐标与质控情况',
    recordsExport: '导出',
    contractTitle: '统计口径',
    contractHint: '口径 = 数据集现状 + 最近一次分析的设计/检验/校正；随状态刷新，导出的 protocol 同源。',
    contractOk: '可引用',
    contractWarn: '待注意',
    contractNoResult: '还没有分析记录',
    contractDesign: '设计',
    contractInference: '最近分析',
    contractCaveats: '引用前必须处理',
    contractGroups: '分组',
    contractScale: '尺度',
    contractSpots: 'spot / 基因',
    contractSamples: '切片 / 供体',
    recordsAnalysis: '分析目录',
    recordsQuickExports: '快速导出',
    recordsExporting: '导出中…',
    auditHint: '得分不是模型评分，只是流程完整性检查；有问题按提示补齐',
    auditOk: '通过',
    auditWarn: '注意',
    auditFail: '异常',
    exportExprCsv: '表达矩阵',
    exportDeCsv: '差异表',
    exportSvgCsv: 'SVG 表',
    exportDomainCsv: '域表',
    exportDeconvCsv: '反卷积比例',
    exportSampleCsv: '样本/spot 表',
    exportSummary: '汇总 JSON',
    exportProtocol: '分析协议',
    exportMethods: 'Methods',
    exportReport: '报告',
    exportScript: '可复现脚本',
    logTitle: '活动日志',
    logRefresh: '刷新日志',
    logEmpty: '还没有活动记录',
    at: '时间',
    actor: '来源',
    action: '动作',
    duration: '耗时',
    resultColumn: '结果',
    activityAgent: 'Agent',
    activityUser: '用户',
    activityDataset: '数据集',
    agentTitle: 'Agent',
    agentPanel: 'Agent',
    agentHint: '也可以在对话里直接让 Agent 操作这个工作台。',
    agentStatus: '状态',
    agentNoDataset: '载入数据集后这里显示当前状态',
    agentSuggestions: '建议下一步',
    agentActivity: '活动',
    agentNoActivity: '暂无活动',
    agentCopy: '复制指令',
    agentCoordsHint: '坐标与邻接半径由导入时估计',
    suggestionFallback: '继续分析',
    geoTitle: 'GEO 检索',
    geoSubtitle: 'NCBI E-utilities 检索空间转录组数据集；中文关键词会自动翻成英文',
    geoTermPh: '如：乳腺 空间转录组 Visium',
    geoSearchBtn: '检索',
    geoSearching: '检索中…',
    geoHint: '中文可直接搜（内置词典翻译）；结果少时会自动放宽关键词。',
    geoSpecies: '物种',
    geoAll: '全部',
    geoHuman: '人',
    geoMouse: '小鼠',
    geoStudyType: '研究类型',
    geoStudyAll: '全部',
    geoStudyArray: '芯片',
    geoStudySeq: '测序',
    geoMinSamples: '最少样本',
    geoHasMatrix: '只看有矩阵',
    geoPlatform: '空间平台',
    geoPlatformAll: '全部平台',
    geoResults: '检索结果',
    geoTotal: '共 {v} 条',
    geoNoResults: '没有结果，换个关键词或放宽筛选',
    geoLoadingDetail: '读取详情中…',
    geoErrorSearch: '检索失败',
    geoErrorDetail: '详情读取失败',
    geoDetailTitle: '数据集详情',
    geoDetailBtn: '点左侧结果查看详情',
    geoSummary: '摘要',
    geoDesign: '实验设计',
    geoFiles: '文件列表',
    geoNoFiles: '没有可下载文件',
    geoSelectAll: '全选',
    geoClear: '清空',
    geoDownloadDir: '下载目录',
    geoDownloadDirHint: '绝对路径；默认 ~/.dsh/spatial/GEO',
    geoPickDir: '选择目录',
    geoPicking: '选择中…',
    geoPickFailed: '目录选择不可用，请手动填写路径',
    geoPickTitle: '选择 GEO 下载目录',
    geoDownloadBtn: '下载',
    geoDownloading: '下载中…',
    geoJobDone: '下载完成',
    geoJobFailed: '下载失败',
    geoImportBtn: '导入并开始清洗',
    geoImporting: '导入中…',
    geoImportDone: '已导入并完成质控',
    geoOpenAfter: '下载并导入',
    geoNeedFile: '请至少勾选一个文件',
    geoRelease: '发布时间',
    geoPubmed: 'PubMed',
    geoMatrix: '有矩阵',
    geoScore: '空间相关度',
    geoScoreHigh: '高',
    geoScoreMid: '中',
    geoScoreLow: '低',
    geoInvalidDir: '下载目录必须是绝对路径',
    jobTimeout: '任务超时，引擎可能仍在计算，请稍后在记录里查看',
    importTitle: '导入数据',
    importSubtitle: '支持 Visium outs/、h5ad、Stereo-seq BIN、Slide-seq、矩阵+坐标、压缩包',
    importSourceLocal: '本机路径',
    importSourceUrl: 'URL',
    importChooseFile: '选择本地文件（浏览器上传）',
    importUploadHint: '上传适合小矩阵；大文件请用本机路径，避免浏览器拷贝',
    importAuto: '自动识别',
    importFormat: '数据格式',
    importCoords: '坐标表路径',
    importCoordsPh: '/absolute/path/coords.csv（列 barcode,x,y）',
    importImage: '底图路径（可选）',
    importImagePh: '/absolute/path/tissue_hires_image.png',
    importAdvanced: '高级选项',
    importSheet: 'Excel sheet',
    importGeneCol: '基因列名',
    importMember: '压缩包内路径',
    importObject: 'RData 对象名',
    importTranspose: '矩阵转置（基因×spot → spot×基因）',
    importSampleSheet: '样本信息表（可选）',
    importHint: '常见错误：矩阵缺坐标 → 会提示需要 coords；压缩包多候选 → 会列出候选文件让你点选。',
    importError: '导入失败',
    importNoPath: '请填写路径或 URL',
    importCandidates: '发现多个候选，点选一个再试',
    importInspect: '只探测格式',
    importInspectBusy: '探测中…',
    importInspectTitle: '格式探测',
    importStart: '导入',
    importBusy: '导入中…',
    importPathPh: '/absolute/path/outs 或 /path/file.h5ad',
    importUrlPh: 'https://…（直链）',
    importPick: '选择文件/目录',
    importPickHint: '原生选择框需 Host 支持；不可用时手动输入绝对路径',
    demo: '载入示例数据',
    demoLoad: '载入示例',
    searchGeo: 'GEO 检索',
    importData: '导入数据',
    geneNotFound: '没有匹配的基因',
    genePickerHint: '输入后回车选中第一个匹配项',
    copyInstruction: '复制指令',
    noSelection: '不选择',
    tableNoteEmpty: '—',
    empty: '—',
    confirm: '确定',
    ok: '好',
    apply: '应用',
    stepParamsTitle: '步骤参数',
    stepParamsHint: '运行时会带上这里的参数',
    analysisParamsTitle: '分析参数',
    clusterKeysTitle: '聚类列',
    currentDataset: '当前数据集',
    refreshPreview: '刷新预览',
    rowsLabel: '行',
  }

  var DICT_EN = {
    nav: 'Spatial',
    appTitle: 'Spatial Transcriptomics Workbench',
    appSubtitle: 'GEO search · Coordinate import · Spatial clustering · SVG · Domains · Deconvolution',
    refresh: 'Refresh',
    refreshHint: 'Re-read workbench state',
    retry: 'Retry',
    loading: 'Loading…',
    loadingData: 'Reading workbench…',
    noData: 'No dataset loaded',
    errorTitle: 'Something went wrong',
    crashTitle: 'This panel failed to render',
    crashHint: 'Other panels still work; retry to recover.',
    close: 'Close',
    cancel: 'Cancel',
    copy: 'Copy',
    copied: 'Copied',
    copyFailed: 'Copy failed (clipboard not allowed)',
    search: 'Search',
    all: 'All',
    none: 'None',
    yes: 'Yes',
    no: 'No',
    unknown: 'Unknown',
    dataset: 'Dataset',
    spots: 'spots',
    samples: 'Samples',
    genes: 'genes',
    slices: 'Slices',
    clusters: 'Clusters',
    domainsLabel: 'Domains',
    species: 'Species',
    spatialMode: 'Spatial mode',
    engine: 'Engine',
    warnings: 'Warnings',
    notes: 'Notes',
    current: 'Current',
    completed: 'Done',
    pending: 'Pending',
    running: 'Running',
    failed: 'Failed',
    plot: 'Plot',
    plotting: 'Plotting…',
    noFigure: 'No figure yet',
    figureFail: 'Image failed to load (file may be cleaned)',
    exportBtn: 'Export',
    exported: 'Exported',
    exportFail: 'Export failed',
    copyPath: 'Copy path',
    pathLabel: 'File path',
    tableRows: '{v} rows',
    tableTruncated: 'showing first {v}',
    emptyTable: 'Empty table',
    total: 'Total',
    selectAll: 'Select all',
    clearSelection: 'Clear',
    of: '/',
    selectDataset: 'Switch dataset',
    closeDataset: 'Close dataset',
    workdir: 'Workdir',
    engineReady: 'Engine ready',
    engineBooting: 'Engine starting',
    engineBad: 'Engine unavailable',
    engineGuideTitle: 'Python engine is not ready',
    engineGuideBody: 'The workbench needs a local Python environment (numpy / pandas / scipy / sklearn / matplotlib / anndata / h5py / networkx). You can browse the page, but import, compute and plotting need the engine.',
    engineGuideSteps: 'Check: 1) engine.log for the traceback; 2) python3 path in plugin settings; 3) restart the plugin or host, then retry.',
    pythonLabel: 'Python',
    engineInfo: 'Engine info',
    statusOk: 'OK',
    statusWarn: 'Warn',
    statusFail: 'Error',
    pipelineTitle: 'Pipeline',
    pipelineHint: 'Run in order; re-running a step clears its downstream results',
    datasetsTitle: 'Workbench',
    datasetNoDataset: 'No dataset yet',
    datasetOpenHint: 'Search GEO, import a Visium outs/ folder or h5ad, or load the built-in demo.',
    datasetStepsDone: '{v} steps done',
    datasetCurrentStep: 'current {v}',
    datasetSpots: 'spots {v}',
    datasetGenes: 'genes {v}',
    datasetSlicesCount: 'slices {v}',
    datasetClustersCount: 'clusters {v}',
    datasetRadius: 'radius {v}',
    datasetRef: 'ref {v}',
    stepStatusDone: 'Done',
    stepStatusRunning: 'Running',
    stepStatusReady: 'Ready',
    stepStatusBlocked: 'Blocked',
    stepStatusFailed: 'Failed',
    stepRun: 'Run',
    stepRerun: 'Re-run',
    stepRunning: 'Running…',
    stepParams: 'Params',
    stepNeedPrereq: 'Needs: ',
    stepSummaryEmpty: '—',
    runStepDone: 'Step finished',
    analysisCatalog: 'Analyses',
    analysisHint: 'Every analysis is read from /spatial/analyses',
    analysisTitle: 'Spatial analyses',
    analysisSubtitle: 'Browse by group; param form is generated from the server schema',
    analysisClose: 'Collapse',
    analysisRun: 'Run analysis',
    analysisRunning: 'Running…',
    analysisNoParams: 'This analysis has no tunable params',
    analysisNeedDataset: 'Load a dataset first',
    analysisUnavailable: 'Unavailable',
    analysisNeeds: 'Needs',
    analysisEmpty: 'No analysis in this group (engine not ready or prerequisites missing)',
    analysisGroupOther: 'Other',
    analysisSelectHint: 'Pick an analysis on the left, tune params, then run',
    analysisResult: 'Result',
    analysisMethod: 'Method',
    analysisRetry: 'Retry',
    tabOverview: 'Overview',
    tabSpatial: 'Spatial map',
    tabExpr: 'Expression',
    tabDe: 'DE',
    tabStats: 'Spatial stats',
    tabDomains: 'Domains',
    tabDeconv: 'Deconvolution',
    tabComm: 'Communication',
    tabRecords: 'Records',
    tabLog: 'Log',
    overviewTitle: 'Dataset summary',
    overviewHint: 'Check data shape and spatial layout first, then pick an analysis route',
    metricSpots: 'Spots',
    metricGenes: 'Genes',
    metricSpatialMode: 'Spatial mode',
    metricSlices: 'Slices',
    metricClusters: 'Clusters',
    metricDomains: 'Domains',
    metricEngine: 'Engine',
    clusterMap: 'Spatial clusters',
    pcaScree: 'PCA scree',
    qcMap: 'Spatial QC',
    overviewClustersEmpty: 'Run the cluster step first to colour by cluster',
    miniMapTitle: 'Coordinate preview',
    miniMapHint: 'Sampled coordinates from the preview endpoint; full map is under “Spatial map”',
    geneQuickTitle: 'Quick gene view',
    geneQuickHint: 'Type a gene (Enter to search) for a spatial map and grouped boxplot',
    geneInputPlaceholder: 'e.g. CDH1, VIM, MKI67',
    nextSteps: 'Next steps',
    noSuggestions: 'No suggestions (engine not ready, or at the end of validated steps)',
    spatialTitle: 'Spatial map',
    spatialHint: 'Plot spots by coordinate; colour by cluster/domain/group/gene; overlay H&E',
    spatialColorBy: 'Colour by',
    spatialGeneSearch: 'Gene',
    spatialGenePh: 'Gene symbol, e.g. COL1A1',
    spatialOverlay: 'Overlay image',
    spatialPointSize: 'Point size',
    spatialMaxPoints: 'Max points',
    spatialPlot: 'Draw map',
    spatialGenePlot: 'Draw gene map',
    spatialNeedGene: 'Enter a gene first',
    spatialNeedCluster: 'No cluster result yet — run the “Cluster” step in the left panel',
    spatialOverlayMissing: 'This dataset has no image; the toggle has no effect',
    colorClusters: 'Cluster',
    colorDomains: 'Domain',
    colorGroup: 'Group',
    colorCellType: 'Cell type',
    colorGene: 'Gene expression',
    exprTitle: 'Expression',
    exprHint: 'Search genes; get spatial expression, grouped boxplot and heatmap',
    exprGeneSearch: 'Search gene',
    exprAddGene: 'Add',
    exprGeneNotInData: 'No matching gene in this dataset',
    exprSelected: 'Selected genes',
    exprEmpty: 'No gene selected yet',
    exprPlotMap: 'Spatial map',
    exprPlotBox: 'Grouped boxplot',
    exprPlotHeatmap: 'Expression heatmap',
    exprExport: 'Export matrix',
    exprExportHint: 'The expression matrix is exported as CSV (genes x spots)',
    exprTableFallback: 'Heatmap/boxplot use the current grouping; the server caps the table at 200 rows.',
    deTitle: 'Regional DE',
    deHint: 'Domain vs rest / between domains / by group, Wilcoxon + BH',
    deEmpty: 'No DE result yet — pick an analysis and run it',
    deExport: 'Export DE table',
    deVolcano: 'Volcano',
    deHeatmap: 'DE heatmap',
    deQuickPlots: 'Quick plots',
    deNeedCluster: 'Requires cluster/domain results',
    statsTitle: 'Spatial statistics',
    statsHint: 'Moran/Geary/permutation/hotspot/gradient/co-expression/Ripley/genome scan',
    statsEmpty: 'No spatial statistics yet',
    domainsTitle: 'Domains and neighborhoods',
    domainsHint: 'MRF domains, markers, neighborhood enrichment, co-occurrence, boundary, distance',
    domainsEmpty: 'No domain result yet — run the “Domains” step or domains_mrf',
    deconvTitle: 'Deconvolution',
    deconvHint: 'Estimate cell-type proportions per spot from a single-cell reference (NNLS/DWLS)',
    deconvEmpty: 'No deconvolution result yet',
    deconvRefTitle: 'Reference',
    deconvRefHint: 'Pick a reference h5ad: a local path, or pull the current dataset from the single-cell workbench',
    deconvRefPath: 'Reference path',
    deconvRefPathPh: '/absolute/path/ref.h5ad',
    deconvRefAttach: 'Attach reference',
    deconvRefAttaching: 'Attaching…',
    deconvRefAttached: 'Reference attached',
    deconvRefNone: 'No reference attached',
    deconvRefCells: 'Cell types',
    deconvScList: 'Single-cell datasets',
    deconvScPull: 'Pull',
    deconvScUnavailable: 'Single-cell workbench unavailable (/singlecell not detected)',
    deconvLocalRefs: 'Local references',
    deconvNoRefs: 'No usable h5ad reference found',
    deconvCellTypeCol: 'Cell-type column',
    deconvRefSet: 'Use as deconvolution reference',
    deconvRefDone: 'Reference ready — you can run deconvolution now',
    deconvNeedRef: 'Deconvolution needs a reference h5ad first',
    commTitle: 'Cell communication',
    commHint: 'Ligand-receptor co-expression + spatial proximity (CellPhoneDB / CellChat style)',
    commEmpty: 'No communication result yet',
    recordsTitle: 'Records and export',
    recordsHint: 'Step history, reproducible script and a 0–100 self-check',
    recordsHistory: 'Step history',
    recordsNoHistory: 'No runs yet',
    recordsAudit: 'Deterministic self-check',
    recordsRunAudit: 'Run check',
    recordsAuditScore: 'Score',
    recordsAuditEmpty: 'Run the check to inspect pipeline completeness, coordinates and QC',
    recordsExport: 'Export',
    contractTitle: 'Statistical contract',
    contractHint: 'Contract = dataset state + design/test/correction of the latest analysis; refreshes with state and matches the exported protocol.',
    contractOk: 'Citable',
    contractWarn: 'Needs attention',
    contractNoResult: 'No analysis yet',
    contractDesign: 'Design',
    contractInference: 'Latest analysis',
    contractCaveats: 'Fix before citing',
    contractGroups: 'Groups',
    contractScale: 'Scale',
    contractSpots: 'Spots / genes',
    contractSamples: 'Slices / donors',
    recordsAnalysis: 'Analyses',
    recordsQuickExports: 'Quick export',
    recordsExporting: 'Exporting…',
    auditHint: 'The score is not a model score; it only checks pipeline completeness',
    auditOk: 'Pass',
    auditWarn: 'Warn',
    auditFail: 'Fail',
    exportExprCsv: 'Expression matrix',
    exportDeCsv: 'DE table',
    exportSvgCsv: 'SVG table',
    exportDomainCsv: 'Domain table',
    exportDeconvCsv: 'Deconvolution proportions',
    exportSampleCsv: 'Sample/spot table',
    exportSummary: 'Summary JSON',
    exportProtocol: 'Analysis protocol',
    exportMethods: 'Methods',
    exportReport: 'Report',
    exportScript: 'Reproducible script',
    logTitle: 'Activity log',
    logRefresh: 'Refresh log',
    logEmpty: 'No activity yet',
    at: 'Time',
    actor: 'Actor',
    action: 'Action',
    duration: 'Duration',
    resultColumn: 'Result',
    activityAgent: 'Agent',
    activityUser: 'User',
    activityDataset: 'Dataset',
    agentTitle: 'Agent',
    agentPanel: 'Agent',
    agentHint: 'You can also drive this workbench from the chat with the Agent.',
    agentStatus: 'Status',
    agentNoDataset: 'Load a dataset to see its status here',
    agentSuggestions: 'Next steps',
    agentActivity: 'Activity',
    agentNoActivity: 'No activity yet',
    agentCopy: 'Copy prompt',
    agentCoordsHint: 'Coordinates and neighbor radius are estimated on import',
    suggestionFallback: 'Continue analysis',
    geoTitle: 'GEO search',
    geoSubtitle: 'Search spatial transcriptomics datasets via NCBI E-utilities; Chinese terms are translated automatically',
    geoTermPh: 'e.g. breast spatial transcriptomics Visium',
    geoSearchBtn: 'Search',
    geoSearching: 'Searching…',
    geoHint: 'Chinese queries work (built-in dictionary); keywords are relaxed automatically when few results.',
    geoSpecies: 'Species',
    geoAll: 'All',
    geoHuman: 'Human',
    geoMouse: 'Mouse',
    geoStudyType: 'Study type',
    geoStudyAll: 'All',
    geoStudyArray: 'Array',
    geoStudySeq: 'Sequencing',
    geoMinSamples: 'Min samples',
    geoHasMatrix: 'Has matrix',
    geoPlatform: 'Spatial platform',
    geoPlatformAll: 'All platforms',
    geoResults: 'Results',
    geoTotal: '{v} hits',
    geoNoResults: 'No hits; try another keyword or relax the filters',
    geoLoadingDetail: 'Loading details…',
    geoErrorSearch: 'Search failed',
    geoErrorDetail: 'Failed to load details',
    geoDetailTitle: 'Dataset details',
    geoDetailBtn: 'Pick a result on the left to see details',
    geoSummary: 'Summary',
    geoDesign: 'Design',
    geoFiles: 'Files',
    geoNoFiles: 'No downloadable file',
    geoSelectAll: 'Select all',
    geoClear: 'Clear',
    geoDownloadDir: 'Download dir',
    geoDownloadDirHint: 'Absolute path; default ~/.dsh/spatial/GEO',
    geoPickDir: 'Choose dir',
    geoPicking: 'Choosing…',
    geoPickFailed: 'Directory picker unavailable; type the path manually',
    geoPickTitle: 'Choose GEO download directory',
    geoDownloadBtn: 'Download',
    geoDownloading: 'Downloading…',
    geoJobDone: 'Download finished',
    geoJobFailed: 'Download failed',
    geoImportBtn: 'Import and start QC',
    geoImporting: 'Importing…',
    geoImportDone: 'Imported and QC finished',
    geoOpenAfter: 'Download & import',
    geoNeedFile: 'Select at least one file',
    geoRelease: 'Released',
    geoPubmed: 'PubMed',
    geoMatrix: 'matrix',
    geoScore: 'Spatial relevance',
    geoScoreHigh: 'high',
    geoScoreMid: 'mid',
    geoScoreLow: 'low',
    geoInvalidDir: 'Download dir must be an absolute path',
    jobTimeout: 'Job timed out; the engine may still be computing — check Records later',
    importTitle: 'Import data',
    importSubtitle: 'Visium outs/, h5ad, Stereo-seq BIN, Slide-seq, matrix + coordinates, archives',
    importSourceLocal: 'Local path',
    importSourceUrl: 'URL',
    importChooseFile: 'Choose a local file (browser upload)',
    importUploadHint: 'Upload suits small matrices; use a local path for big files',
    importAuto: 'Auto detect',
    importFormat: 'Data format',
    importCoords: 'Coordinates path',
    importCoordsPh: '/absolute/path/coords.csv (columns barcode,x,y)',
    importImage: 'Image path (optional)',
    importImagePh: '/absolute/path/tissue_hires_image.png',
    importAdvanced: 'Advanced',
    importSheet: 'Excel sheet',
    importGeneCol: 'Gene column',
    importMember: 'Path inside archive',
    importObject: 'RData object',
    importTranspose: 'Transpose (genes x spots to spots x genes)',
    importSampleSheet: 'Sample sheet (optional)',
    importHint: 'Common errors: missing coordinates (server asks for coords); archives with several candidates (pick one).',
    importError: 'Import failed',
    importNoPath: 'Fill in a path or URL',
    importCandidates: 'Several candidates found — pick one and retry',
    importInspect: 'Inspect only',
    importInspectBusy: 'Inspecting…',
    importInspectTitle: 'Format probe',
    importStart: 'Import',
    importBusy: 'Importing…',
    importPathPh: '/absolute/path/outs or /path/file.h5ad',
    importUrlPh: 'https://… (direct link)',
    importPick: 'Choose file/folder',
    importPickHint: 'Native picker needs Host support; type an absolute path otherwise',
    demo: 'Load demo',
    demoLoad: 'Load demo',
    searchGeo: 'GEO search',
    importData: 'Import',
    geneNotFound: 'No matching gene',
    genePickerHint: 'Press Enter to pick the first match',
    copyInstruction: 'Copy prompt',
    noSelection: 'None',
    tableNoteEmpty: '—',
    empty: '—',
    confirm: 'OK',
    ok: 'OK',
    apply: 'Apply',
    stepParamsTitle: 'Step params',
    stepParamsHint: 'These params are sent when you run the step',
    analysisParamsTitle: 'Analysis params',
    clusterKeysTitle: 'Cluster columns',
    currentDataset: 'Current dataset',
    refreshPreview: 'Refresh preview',
    rowsLabel: 'rows',
  }
  // ==========================================================================
  // 纯函数层：格式化 / 取色 / 坐标变换 / payload / 参数 / 归一化 / API / 轮询
  // 不碰 DOM、不 require；可在 node 里直接单测（见文件尾 __SP_TEST__）。
  // ==========================================================================

  function isNil(value) {
    return value === null || value === undefined
  }

  function clamp(value, min, max) {
    var number = Number(value)
    if (!isFinite(number)) return min
    if (number < min) return min
    if (number > max) return max
    return number
  }

  function asNumber(value) {
    if (isNil(value) || value === '') return null
    var number = Number(value)
    return isFinite(number) ? number : null
  }

  function plainText(value, max) {
    var text = isNil(value) ? '' : String(value).replace(/\s+/g, ' ').trim()
    if (typeof max === 'number' && max > 0 && text.length > max) return text.slice(0, max)
    return text
  }

  function clampInt(value, min, max, fallback) {
    var number = asNumber(value)
    if (number === null) return fallback
    return Math.round(clamp(number, min, max))
  }

  function trimFixed(number, digits) {
    return Number(number).toFixed(clamp(digits, 0, 6))
  }

  /** 1234 → 1,234；25000 → 25K；1500000 → 1.5M；空值 → —。 */
  function fmtNumber(value, digits) {
    var number = asNumber(value)
    if (number === null) return '—'
    var fixed = isFinite(Number(digits)) ? clamp(digits, 0, 6) : 2
    var abs = Math.abs(number)
    var sign = number < 0 ? '-' : ''
    if (abs >= 1e9) return sign + trimFixed(abs / 1e9, fixed > 0 ? fixed : 1).replace(/\.0+$/, '') + 'B'
    if (abs >= 1e6) return sign + trimFixed(abs / 1e6, fixed > 0 ? fixed : 1).replace(/\.0+$/, '') + 'M'
    if (abs >= 1e4) return sign + trimFixed(abs / 1e3, fixed > 0 ? fixed : 1).replace(/\.0+$/, '') + 'K'
    if (abs >= 1000) {
      var parts = trimFixed(abs, fixed).split('.')
      var head = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, ',')
      return sign + head + (parts[1] && parts[1].length ? '.' + parts[1] : '')
    }
    var fixedText = trimFixed(abs, fixed)
    return sign + (fixedText.indexOf('.') >= 0 ? fixedText.replace(/0+$/, '').replace(/\.$/, '') : fixedText)
  }

  /** 0.5 → 50%；37.5 → 37.5%；空值 → —。 */
  function fmtPct(value, digits) {
    var number = asNumber(value)
    if (number === null) return '—'
    var fixed = isFinite(Number(digits)) ? clamp(digits, 0, 4) : 1
    var percent = Math.abs(number) <= 1 ? number * 100 : number
    var text = trimFixed(percent, fixed)
    if (text.indexOf('.') >= 0) text = text.replace(/0+$/, '').replace(/\.$/, '')
    return text + '%'
  }

  /** P 值：0.05 → 0.050；5e-4 → 5.00e-4；0 → 0。 */
  function fmtP(value) {
    var number = asNumber(value)
    if (number === null) return '—'
    if (number === 0) return '0'
    if (number < 1e-3) {
      var text = Number(number).toExponential(2)
      return text.replace('e', 'e')
    }
    if (number < 0.01) return Number(number).toFixed(4)
    return Number(number).toFixed(3)
  }

  /** 字节 → 1024 进制短单位。 */
  function formatBytes(bytes) {
    var number = asNumber(bytes)
    if (number === null || number < 0) return '—'
    var units = ['B', 'KB', 'MB', 'GB', 'TB']
    var exponent = Math.min(Math.floor(Math.log(number) / Math.log(1024)), units.length - 1)
    if (number <= 0) return '0 B'
    var value = number / Math.pow(1024, exponent)
    return (Math.round(value * 10) / 10).toFixed(1).replace(/\.0$/, '') + ' ' + units[exponent]
  }

  function formatDuration(ms) {
    var number = asNumber(ms)
    if (number === null || number < 0) return '—'
    if (number < 1000) return Math.round(number) + ' ms'
    if (number < 60000) return (Math.round(number / 100) / 10).toFixed(1).replace(/\.0$/, '') + ' s'
    var totalSec = Math.round(number / 1000)
    var minutes = Math.floor(totalSec / 60)
    var seconds = totalSec % 60
    if (minutes < 60) return minutes + 'm ' + (seconds < 10 ? '0' + seconds : String(seconds)) + 's'
    var hours = Math.floor(minutes / 60)
    return hours + 'h ' + (minutes % 60) + 'm'
  }

  function formatClock(value) {
    if (isNil(value) || value === '') return '—'
    var date = value instanceof Date ? value : new Date(value)
    if (isNaN(date.getTime())) {
      var text = String(value)
      return text.length > 16 ? text.slice(5, 16).replace('T', ' ') : text
    }
    function pad(number) { return number < 10 ? '0' + number : String(number) }
    return pad(date.getMonth() + 1) + '-' + pad(date.getDate()) + ' ' + pad(date.getHours()) + ':' + pad(date.getMinutes())
  }

  function formatError(error) {
    if (isNil(error)) return ''
    if (typeof error === 'string') return error
    if (error && typeof error.message === 'string' && error.message) return error.message
    return String(error)
  }

  function hexToRgb(hex) {
    var text = String(hex || '').replace('#', '')
    if (text.length === 3) text = text.charAt(0) + text.charAt(0) + text.charAt(1) + text.charAt(1) + text.charAt(2) + text.charAt(2)
    if (!/^[0-9a-fA-F]{6}$/.test(text)) return { r: 128, g: 128, b: 128 }
    return {
      r: parseInt(text.slice(0, 2), 16),
      g: parseInt(text.slice(2, 4), 16),
      b: parseInt(text.slice(4, 6), 16),
    }
  }

  /** 同名同色：数字当调色板下标，字符串取稳定哈希。 */
  function colorForGroup(value) {
    var index = null
    if (typeof value === 'number' && isFinite(value)) {
      index = Math.abs(Math.round(value)) % PALETTE.length
      return PALETTE[index]
    }
    var text = isNil(value) ? '' : String(value)
    if (!text) return PALETTE[PALETTE.length - 1]
    var hash = 0
    for (var cursor = 0; cursor < text.length; cursor++) {
      hash = ((hash << 5) - hash + text.charCodeAt(cursor)) | 0
    }
    index = Math.abs(hash) % PALETTE.length
    return PALETTE[index]
  }

  /** 连续值 → 色带颜色（低→高）；min/max 相同或非法时落到中位色。 */
  function continuousColor(value, min, max) {
    var number = asNumber(value)
    if (number === null) return PALETTE[PALETTE.length - 1]
    var low = asNumber(min)
    var high = asNumber(max)
    if (low === null || high === null || !(high > low)) return CONTINUOUS_STOPS[Math.floor(CONTINUOUS_STOPS.length / 2)]
    var ratio = clamp((number - low) / (high - low), 0, 1)
    var span = CONTINUOUS_STOPS.length - 1
    var position = ratio * span
    var leftIndex = Math.floor(position)
    var rightIndex = Math.min(leftIndex + 1, span)
    var mix = position - leftIndex
    var left = hexToRgb(CONTINUOUS_STOPS[leftIndex])
    var right = hexToRgb(CONTINUOUS_STOPS[rightIndex])
    var channel = function (a, b) {
      var mixed = Math.round(a + (b - a) * mix)
      var hex = mixed.toString(16)
      return hex.length === 1 ? '0' + hex : hex
    }
    return '#' + channel(left.r, right.r) + channel(left.g, right.g) + channel(left.b, right.b)
  }

  function groupBy(list, keyFn) {
    var source = Array.isArray(list) ? list : []
    var pick = typeof keyFn === 'function'
      ? keyFn
      : function (item) { return item && typeof item === 'object' ? item[keyFn] : item }
    var out = {}
    for (var index = 0; index < source.length; index++) {
      var key = String(pick(source[index], index))
      if (!out[key]) out[key] = []
      out[key].push(source[index])
    }
    return out
  }

  // --------------------------------------------------------------------------
  // 空间坐标 → 像素：把 (x, y) 等比映射进 width × height 的画布/SVG。
  // 语义：
  //   - y 轴与图像坐标一致（原点左上、y 向下增），直接线性映射，不翻转；
  //     需要数学坐标（y 向上增）时传 options.flipY = true。
  //   - options.uniform !== false 时保持等比（默认），多余空间居中留白。
  //   - bbox 省略或非法时用点的包围盒。
  // 返回 [{ x, y, index }]，无效点被跳过（不会打断整批）。
  // --------------------------------------------------------------------------

  function normalizeExtent(bbox) {
    var raw = bbox && typeof bbox === 'object' ? bbox : null
    if (!raw) return null
    var xRange = Array.isArray(raw.x) ? raw.x : (Array.isArray(raw.x_range) ? raw.x_range : null)
    var yRange = Array.isArray(raw.y) ? raw.y : (Array.isArray(raw.y_range) ? raw.y_range : null)
    var x0 = xRange ? asNumber(xRange[0]) : null
    var x1 = xRange ? asNumber(xRange[1]) : null
    var y0 = yRange ? asNumber(yRange[0]) : null
    var y1 = yRange ? asNumber(yRange[1]) : null
    if (x0 === null || x1 === null || y0 === null || y1 === null) return null
    return {
      x: [Math.min(x0, x1), Math.max(x0, x1)],
      y: [Math.min(y0, y1), Math.max(y0, y1)],
    }
  }

  function pointXY(point) {
    if (!point) return null
    if (Array.isArray(point)) {
      var ax = asNumber(point[0])
      var ay = asNumber(point[1])
      if (ax === null || ay === null) return null
      return { x: ax, y: ay }
    }
    if (typeof point === 'object') {
      var x = asNumber(point.x !== undefined ? point.x : point[0])
      var y = asNumber(point.y !== undefined ? point.y : point[1])
      if (x === null || y === null) return null
      return { x: x, y: y }
    }
    return null
  }

  /** 点集包围盒 → { x: [min, max], y: [min, max] }；无有效点返回 null。 */
  function extentOf(points) {
    var source = Array.isArray(points) ? points : []
    var minX = null
    var maxX = null
    var minY = null
    var maxY = null
    for (var index = 0; index < source.length; index++) {
      var point = pointXY(source[index])
      if (!point) continue
      if (minX === null || point.x < minX) minX = point.x
      if (maxX === null || point.x > maxX) maxX = point.x
      if (minY === null || point.y < minY) minY = point.y
      if (maxY === null || point.y > maxY) maxY = point.y
    }
    if (minX === null) return null
    return { x: [minX, maxX], y: [minY, maxY] }
  }

  function spatialScreen(points, bbox, width, height, options) {
    var opts = options && typeof options === 'object' ? options : {}
    var boxWidth = asNumber(width)
    var boxHeight = asNumber(height)
    if (boxWidth === null || boxWidth <= 0) boxWidth = 640
    if (boxHeight === null || boxHeight <= 0) boxHeight = 480
    var padding = asNumber(opts.padding)
    if (padding === null || padding < 0) padding = 8
    var source = Array.isArray(points) ? points : []
    var valid = []
    for (var index = 0; index < source.length; index++) {
      var point = pointXY(source[index])
      if (point) valid.push({ x: point.x, y: point.y, index: index })
    }
    if (!valid.length) return []
    var extent = normalizeExtent(bbox) || extentOf(valid)
    var minX = extent.x[0]
    var maxX = extent.x[1]
    var minY = extent.y[0]
    var maxY = extent.y[1]
    var spanX = maxX - minX
    var spanY = maxY - minY
    if (!(spanX > 0)) spanX = 1
    if (!(spanY > 0)) spanY = 1
    var innerWidth = Math.max(boxWidth - padding * 2, 1)
    var innerHeight = Math.max(boxHeight - padding * 2, 1)
    var uniform = opts.uniform !== false
    var scaleX = innerWidth / spanX
    var scaleY = innerHeight / spanY
    var scale = uniform ? Math.min(scaleX, scaleY) : 1
    var drawWidth = uniform ? spanX * scale : innerWidth
    var drawHeight = uniform ? spanY * scale : innerHeight
    var offsetX = padding + (innerWidth - drawWidth) / 2
    var offsetY = padding + (innerHeight - drawHeight) / 2
    var result = []
    for (var cursor = 0; cursor < valid.length; cursor++) {
      var item = valid[cursor]
      var sx = offsetX + (item.x - minX) * (uniform ? scale : scaleX)
      var sy = offsetY + (item.y - minY) * (uniform ? scale : scaleY)
      if (opts.flipY === true) sy = boxHeight - sy
      result.push({ x: Math.round(sx * 100) / 100, y: Math.round(sy * 100) / 100, index: item.index })
    }
    return result
  }

  // --------------------------------------------------------------------------
  // GEO 搜索 payload：只带非空/非默认的可选字段（服务端自己补默认值）。
  // --------------------------------------------------------------------------

  function buildSearchPayload(filters) {
    var raw = filters && typeof filters === 'object' ? filters : {}
    var payload = {
      term: plainText(raw.term, 300),
      limit: clampInt(raw.limit, 1, 50, 20),
      page: clampInt(raw.page, 0, 1000, 0),
      sort: 'relevance',
    }
    var sort = plainText(raw.sort)
    if (sort === 'date' || sort === 'samples' || sort === 'relevance') payload.sort = sort
    var species = plainText(raw.species)
    if (species && species !== 'all') payload.species = species
    var studyType = plainText(raw.study_type)
    if (studyType && studyType !== 'all') payload.study_type = studyType
    var minSamples = asNumber(raw.min_samples)
    if (minSamples !== null && minSamples > 0) payload.min_samples = Math.round(minSamples)
    if (raw.has_matrix === true) payload.has_matrix = true
    var platform = plainText(raw.platform)
    if (platform && platform !== 'all') payload.platform = platform
    return payload
  }

  /** 导入参数：空值一律剔除；data_format=auto 不发。 */
  function buildImportPayload(options) {
    var raw = options && typeof options === 'object' ? options : {}
    var out = {}
    var path = plainText(raw.path)
    if (path) out.path = path
    var url = plainText(raw.url)
    if (url) out.url = url
    var name = plainText(raw.name)
    if (name) out.name = name
    var format = plainText(raw.data_format)
    if (format && format !== 'auto') out.data_format = format
    var coords = plainText(raw.coords)
    if (coords) out.coords = coords
    var image = plainText(raw.image)
    if (image) out.image = image
    var member = plainText(raw.member)
    if (member) out.member = member
    var object = plainText(raw.object)
    if (object) out.object = object
    var sheet = plainText(raw.sheet)
    if (sheet) out.sheet = sheet
    var geneCol = plainText(raw.gene_col)
    if (geneCol) out.gene_col = geneCol
    var sampleSheet = plainText(raw.sample_sheet)
    if (sampleSheet) out.sample_sheet = sampleSheet
    var delimiter = plainText(raw.delimiter)
    if (delimiter) out.delimiter = delimiter
    var visiumScale = asNumber(raw.visium_scale)
    if (visiumScale !== null && visiumScale > 0) out.visium_scale = visiumScale
    if (raw.transpose === true) out.transpose = true
    return out
  }

  /** 建议 → 可直接粘给 Agent 的中文指令。 */
  function suggestionInstruction(suggestion) {
    var raw = suggestion && typeof suggestion === 'object' ? suggestion : {}
    var title = plainText(raw.title, 120) || '继续分析'
    var detail = plainText(raw.detail, 400)
    var step = plainText(raw.step, 40)
    var meta = step && STEP_BY_ID[step] ? STEP_BY_ID[step] : null
    var text = '请在工作台执行「' + title + '」'
    if (meta) text += '（步骤：' + meta.title + '，step=' + step + '）'
    else if (step) text += '（step=' + step + '）'
    if (detail) text += '：' + detail
    return text + '。'
  }

  // --------------------------------------------------------------------------
  // 步骤参数：默认值 / 收集（按 descriptor.type 收敛）/ 依赖判定
  // --------------------------------------------------------------------------

  function parseTokenList(value) {
    if (Array.isArray(value)) {
      var out = []
      for (var index = 0; index < value.length; index++) {
        var text = plainText(value[index], 200)
        if (text) out.push(text)
      }
      return out
    }
    var raw = plainText(value, 4000)
    if (!raw) return []
    return raw.split(/[,，;；\s|]+/).map(function (item) { return item.trim() }).filter(Boolean)
  }

  function cloneDefault(value) {
    if (Array.isArray(value)) return value.slice()
    if (value && typeof value === 'object') {
      var out = {}
      var keys = Object.keys(value)
      for (var index = 0; index < keys.length; index++) out[keys[index]] = value[keys[index]]
      return out
    }
    return value
  }

  /** 步骤默认参数（数组/对象复制，避免多次运行共享同一引用）。 */
  function defaultParams(step) {
    var meta = typeof step === 'string' ? STEP_BY_ID[step] : step
    var out = {}
    if (!meta || !Array.isArray(meta.params)) return out
    for (var index = 0; index < meta.params.length; index++) {
      var descriptor = meta.params[index]
      if (!descriptor || !descriptor.key) continue
      out[descriptor.key] = cloneDefault(descriptor.default)
    }
    return out
  }

  function parseJsonParam(text) {
    var raw = plainText(text, 8000)
    if (!raw) return null
    try {
      var parsed = JSON.parse(raw)
      if (parsed && typeof parsed === 'object') return parsed
      return null
    } catch (ignored) {
      return null
    }
  }

  /** 按 descriptor 类型收敛表单值；空值统一成 null（数组类型成 []）。 */
  function collectStepParams(step, values) {
    var meta = typeof step === 'string' ? STEP_BY_ID[step] : step
    var merged = defaultParams(meta)
    var source = values && typeof values === 'object' ? values : {}
    var out = {}
    if (!meta || !Array.isArray(meta.params)) return out
    for (var index = 0; index < meta.params.length; index++) {
      var descriptor = meta.params[index]
      if (!descriptor || !descriptor.key) continue
      var value = source[descriptor.key]
      if (value === undefined) value = merged[descriptor.key]
      var type = descriptor.type || 'text'
      if (type === 'number') {
        out[descriptor.key] = asNumber(value)
        continue
      }
      if (type === 'boolean') {
        out[descriptor.key] = value === true
        continue
      }
      if (type === 'multiselect' || type === 'genes' || type === 'columns' || type === 'array') {
        out[descriptor.key] = parseTokenList(value)
        continue
      }
      if (type === 'json') {
        out[descriptor.key] = parseJsonParam(value)
        continue
      }
      if (type === 'select' || type === 'column' || type === 'file' || type === 'text' || type === 'string') {
        var text = plainText(value, 2000)
        out[descriptor.key] = text || null
        continue
      }
      out[descriptor.key] = value === undefined ? null : value
    }
    return out
  }

  function formatParamValue(descriptor, value) {
    var type = descriptor && descriptor.type ? descriptor.type : 'text'
    if (type === 'multiselect') {
      var labels = {}
      var options = descriptor && Array.isArray(descriptor.options) ? descriptor.options : []
      for (var index = 0; index < options.length; index++) labels[String(options[index].value)] = options[index].label
      return parseTokenList(value).map(function (item) { return labels[item] || item }).join(' / ')
    }
    if (Array.isArray(value)) return value.join(', ')
    if (isNil(value) || value === '') return '—'
    if (typeof value === 'boolean') return value ? '是' : '否'
    return String(value)
  }

  function doneStepSet(dataset) {
    var set = {}
    var steps = dataset && Array.isArray(dataset.steps) ? dataset.steps : []
    for (var index = 0; index < steps.length; index++) set[String(steps[index])] = true
    return set
  }

  /** 前置是否全部完成（requires 与 DESIGN.md §3 一致）。 */
  function isStepReady(step, dataset) {
    var meta = typeof step === 'string' ? STEP_BY_ID[step] : step
    if (!meta || !dataset) return false
    var done = doneStepSet(dataset)
    var requires = Array.isArray(meta.requires) ? meta.requires : []
    for (var index = 0; index < requires.length; index++) {
      if (!done[requires[index]]) return false
    }
    return true
  }

  /** 左栏状态：running > failed > done > ready > blocked。 */
  function stepStatusOf(step, dataset, busyStepId, failedStepId) {
    var meta = typeof step === 'string' ? STEP_BY_ID[step] : step
    if (!meta) return 'blocked'
    if (busyStepId === meta.id) return 'running'
    if (failedStepId === meta.id) return 'failed'
    if (dataset && doneStepSet(dataset)[meta.id]) return 'done'
    if (dataset && isStepReady(meta, dataset)) return 'ready'
    return 'blocked'
  }

  function stepStatusKey(status) {
    if (status === 'done') return 'stepStatusDone'
    if (status === 'running') return 'stepStatusRunning'
    if (status === 'ready') return 'stepStatusReady'
    if (status === 'failed') return 'stepStatusFailed'
    return 'stepStatusBlocked'
  }

  // --------------------------------------------------------------------------
  // 图 URL / 结果归一化（DESIGN.md §2.5 result 信封）
  // --------------------------------------------------------------------------

  function figureSrc(raw) {
    if (typeof raw === 'string' && raw) {
      if (raw.indexOf('data:') === 0 || raw.charAt(0) === '/') return raw
      var nameOnly = raw.split('/').pop()
      if (/^[A-Za-z0-9._-]+\.png$/.test(nameOnly)) return API + '/figure?id=' + encodeURIComponent(nameOnly)
      return ''
    }
    if (!raw || typeof raw !== 'object') return ''
    if (typeof raw.url === 'string' && raw.url) return raw.url
    if (typeof raw.file === 'string' && raw.file) {
      var name = raw.file.split('/').pop()
      if (/^[A-Za-z0-9._-]+\.png$/.test(name)) return API + '/figure?id=' + encodeURIComponent(name)
    }
    return ''
  }

  function normalizeFigure(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.figure && typeof payload.figure === 'object' ? payload.figure : payload)
      : {}
    var file = typeof raw.file === 'string' ? raw.file : ''
    var url = figureSrc({ url: typeof raw.url === 'string' ? raw.url : '', file: file })
    if (!url) return null
    return {
      id: isNil(raw.id) ? '' : String(raw.id),
      title: plainText(raw.title || raw.kind || '', 160),
      url: url,
      file: file,
      width: asNumber(raw.width) || 0,
      height: asNumber(raw.height) || 0,
      kind: typeof raw.kind === 'string' ? raw.kind : '',
      bytes: asNumber(raw.bytes),
    }
  }

  function normalizeTable(table, index) {
    var raw = table && typeof table === 'object' ? table : {}
    var columns = []
    if (Array.isArray(raw.columns)) {
      for (var columnIndex = 0; columnIndex < raw.columns.length; columnIndex++) {
        var column = raw.columns[columnIndex] && typeof raw.columns[columnIndex] === 'object' ? raw.columns[columnIndex] : {}
        var key = plainText(column.key || column.name || ('' + columnIndex))
        if (!key) continue
        columns.push({
          key: key,
          label: plainText(column.label || key),
          type: column.type === 'number' || column.numeric === true ? 'number' : 'string',
        })
      }
    }
    var rows = Array.isArray(raw.rows) ? raw.rows : []
    if (!columns.length && rows.length && rows[0] && typeof rows[0] === 'object') {
      var firstKeys = Object.keys(rows[0])
      for (var keyIndex = 0; keyIndex < firstKeys.length; keyIndex++) {
        columns.push({
          key: firstKeys[keyIndex],
          label: firstKeys[keyIndex],
          type: typeof rows[0][firstKeys[keyIndex]] === 'number' ? 'number' : 'string',
        })
      }
    }
    return {
      id: plainText(raw.id || ('table-' + index)),
      title: plainText(raw.title || ''),
      columns: columns,
      rows: rows,
      note: plainText(raw.note || '', 300),
      truncated: raw.truncated === true,
    }
  }

  function normalizeStringList(value) {
    var list = Array.isArray(value) ? value : (isNil(value) || value === '' ? [] : [value])
    var out = []
    for (var index = 0; index < list.length; index++) {
      var text = plainText(list[index], 500)
      if (text) out.push(text)
    }
    return out
  }

  /** 步骤/分析/导出共用的结果信封归一化。 */
  function normalizeResult(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    if (raw.result && typeof raw.result === 'object' && !raw.summary) {
      raw = Object.assign({}, raw.result, {
        analysis_id: raw.analysis_id || raw.result.analysis_id,
        ms: raw.ms || raw.result.ms,
        dataset: raw.dataset || raw.result.dataset,
      })
    }
    var summary = raw.summary && typeof raw.summary === 'object' ? raw.summary : {}
    var bullets = []
    if (Array.isArray(summary.bullets)) {
      for (var index = 0; index < summary.bullets.length; index++) {
        var bullet = plainText(summary.bullets[index], 500)
        if (bullet) bullets.push(bullet)
      }
    }
    var metrics = {}
    if (summary.metrics && typeof summary.metrics === 'object') {
      var metricKeys = Object.keys(summary.metrics)
      for (var metricIndex = 0; metricIndex < metricKeys.length; metricIndex++) {
        metrics[metricKeys[metricIndex]] = summary.metrics[metricKeys[metricIndex]]
      }
    }
    var tables = []
    if (Array.isArray(raw.tables)) {
      for (var tableIndex = 0; tableIndex < raw.tables.length; tableIndex++) {
        tables.push(normalizeTable(raw.tables[tableIndex], tableIndex))
      }
    }
    var figures = []
    if (Array.isArray(raw.figures)) {
      for (var figureIndex = 0; figureIndex < raw.figures.length; figureIndex++) {
        var figure = normalizeFigure(raw.figures[figureIndex])
        if (figure) figures.push(figure)
      }
    }
    var dataOut = []
    if (Array.isArray(raw.data_out)) {
      for (var outIndex = 0; outIndex < raw.data_out.length; outIndex++) {
        var item = raw.data_out[outIndex] && typeof raw.data_out[outIndex] === 'object' ? raw.data_out[outIndex] : {}
        if (!item.path) continue
        dataOut.push({ kind: plainText(item.kind || ''), path: String(item.path), label: plainText(item.label || item.kind || '') })
      }
    }
    return {
      title: plainText(summary.title || raw.title || '', 200),
      summary: { title: plainText(summary.title || '', 200), bullets: bullets, metrics: metrics },
      bullets: bullets,
      metrics: metrics,
      tables: tables,
      figures: figures,
      data_out: dataOut,
      warnings: normalizeStringList(raw.warnings),
      notes: normalizeStringList(raw.notes),
      analysis_id: plainText(raw.analysis_id || ''),
      step: plainText(raw.step || ''),
      ms: asNumber(raw.ms),
      raw: raw,
    }
  }

  /** dataset 摘要归一化：服务端缺字段（老版本/引擎忙）时也能渲染。 */
  /**
   * 工作台口径契约 → UI 模型（缺失一律 null/空，不编默认值）。
   * 数据来源：Host `lib/impl.mjs` 的 buildWorkbenchContract（引擎真实字段）。
   */
  function normalizeStatsContract(payload) {
    var raw = payload && typeof payload === 'object' ? payload : null
    if (!raw) return null
    function objectOrEmpty(value) {
      return value && typeof value === 'object' && !Array.isArray(value) ? value : {}
    }
    function textList(value) {
      return Array.isArray(value) ? value.map(function (entry) { return plainText(entry, 120) }).filter(Boolean) : []
    }
    function boolOrNull(value) {
      if (value === undefined || value === null) return null
      return value === true
    }
    var dataset = objectOrEmpty(raw.dataset)
    var design = objectOrEmpty(raw.design)
    var inference = raw.inference && typeof raw.inference === 'object' ? raw.inference : null
    var identity = objectOrEmpty(raw.identity)
    var nPerGroup = objectOrEmpty(design.n_per_group)
    return {
      schema: plainText(raw.schema || ''),
      plugin: plainText(raw.plugin || ''),
      pluginVersion: plainText(raw.plugin_version || ''),
      okToCite: raw.ok_to_cite === true,
      dataset: {
        id: plainText(dataset.id || ''),
        name: plainText(dataset.name || ''),
        nSamples: asNumber(dataset.n_samples),
        nGenes: asNumber(dataset.n_genes),
        nSpots: asNumber(dataset.n_spots),
        species: plainText(dataset.species || ''),
        valueKind: plainText(dataset.value_kind || ''),
        dataKind: plainText(dataset.data_kind || ''),
        spatialMode: plainText(dataset.spatial_mode || ''),
        coordsKind: plainText(dataset.coords_kind || ''),
        hasCounts: boolOrNull(dataset.has_counts),
        hasCoords: boolOrNull(dataset.has_coords),
        fingerprint: plainText(dataset.import_fingerprint || ''),
      },
      design: {
        groupKey: plainText(design.group_key || ''),
        groups: textList(design.groups),
        nPerGroup: nPerGroup,
        donorKey: plainText(design.donor_key || ''),
        nDonors: asNumber(design.n_donors),
        sampleKey: plainText(design.sample_key || ''),
        sampleColumns: textList(design.sample_columns),
        clinicalColumns: textList(design.clinical_columns),
      },
      inference: inference ? {
        lastStep: plainText(inference.last_step || ''),
        ok: inference.ok !== false,
        at: plainText(inference.at || ''),
        ms: asNumber(inference.ms),
        summary: plainText(inference.summary || ''),
        method: objectOrEmpty(inference.method),
        designParams: objectOrEmpty(inference.design),
        fdr: objectOrEmpty(inference.fdr),
      } : null,
      recentAnalyses: Array.isArray(raw.recent_analyses) ? raw.recent_analyses.map(function (entry) {
        var item = entry && typeof entry === 'object' ? entry : {}
        return { analysis: plainText(item.analysis || ''), at: plainText(item.at || ''), ok: item.ok !== false }
      }).filter(function (item) { return Boolean(item.analysis) }) : [],
      identity: {
        fingerprint: plainText(identity.fingerprint || ''),
        createdAt: plainText(identity.created_at || ''),
        nStepsDone: asNumber(identity.n_steps_done),
        nHistory: asNumber(identity.n_history),
      },
      caveats: Array.isArray(raw.caveats) ? raw.caveats.map(function (entry) {
        var item = entry && typeof entry === 'object' ? entry : {}
        return { code: plainText(item.code || ''), level: plainText(item.level || 'info'), text: plainText(item.text || '') }
      }).filter(function (item) { return Boolean(item.text) }) : [],
    }
  }

  function normalizeDataset(payload) {
    var raw = payload && typeof payload === 'object'
      ? (payload.dataset && typeof payload.dataset === 'object' ? payload.dataset : payload)
      : {}
    var id = plainText(raw.id || raw.dataset || raw.dataset_id || '')
    var steps = Array.isArray(raw.steps) ? raw.steps.map(function (item) { return plainText(item) }).filter(Boolean) : []
    var sampleColumns = []
    if (Array.isArray(raw.sample_columns)) {
      for (var columnIndex = 0; columnIndex < raw.sample_columns.length; columnIndex++) {
        var column = raw.sample_columns[columnIndex] && typeof raw.sample_columns[columnIndex] === 'object' ? raw.sample_columns[columnIndex] : {}
        var key = plainText(column.key || column.name || '')
        if (!key) continue
        sampleColumns.push({
          key: key,
          kind: plainText(column.kind || 'categorical'),
          n_unique: asNumber(column.n_unique),
        })
      }
    }
    var geometryRaw = raw.geometry && typeof raw.geometry === 'object' ? raw.geometry : {}
    var xRange = Array.isArray(geometryRaw.x_range) ? geometryRaw.x_range : null
    var yRange = Array.isArray(geometryRaw.y_range) ? geometryRaw.y_range : null
    var clusterKeys = []
    if (Array.isArray(raw.cluster_keys)) {
      for (var keyIndex = 0; keyIndex < raw.cluster_keys.length; keyIndex++) {
        var clusterKey = plainText(raw.cluster_keys[keyIndex])
        if (clusterKey) clusterKeys.push(clusterKey)
      }
    }
    var ref = null
    if (raw.ref && typeof raw.ref === 'object' && plainText(raw.ref.path)) {
      ref = {
        path: plainText(raw.ref.path),
        cell_type_col: plainText(raw.ref.cell_type_col || ''),
        cell_types: Array.isArray(raw.ref.cell_types) ? raw.ref.cell_types.map(function (item) { return plainText(item) }).filter(Boolean) : [],
      }
    } else if (typeof raw.ref === 'string' && raw.ref) {
      ref = { path: raw.ref, cell_type_col: '', cell_types: [] }
    }
    var history = []
    if (Array.isArray(raw.history)) {
      var entries = raw.history.slice(-MAX_HISTORY_ROWS)
      for (var historyIndex = 0; historyIndex < entries.length; historyIndex++) {
        var entry = entries[historyIndex] && typeof entries[historyIndex] === 'object' ? entries[historyIndex] : {}
        history.push({
          at: plainText(entry.at || entry.time || ''),
          step: plainText(entry.step || entry.action || ''),
          ok: entry.ok !== false && entry.status !== 'failed',
          ms: asNumber(entry.ms),
          summary: plainText(entry.summary || entry.message || '', 240),
          actor: plainText(entry.actor || ''),
        })
      }
    }
    var runs = []
    if (Array.isArray(raw.analyses_run)) {
      for (var runIndex = 0; runIndex < raw.analyses_run.length; runIndex++) {
        var run = raw.analyses_run[runIndex] && typeof raw.analyses_run[runIndex] === 'object' ? raw.analyses_run[runIndex] : {}
        runs.push({ at: plainText(run.at || ''), id: plainText(run.id || ''), ok: run.ok !== false, ms: asNumber(run.ms) })
      }
    }
    return {
      id: id,
      name: plainText(raw.name || id || ''),
      path: typeof raw.path === 'string' ? raw.path : '',
      source_format: plainText(raw.source_format || raw.format || ''),
      n_spots: asNumber(raw.n_spots !== undefined ? raw.n_spots : raw.n_samples),
      n_genes: asNumber(raw.n_genes),
      species: plainText(raw.species || ''),
      geometry: {
        mode: plainText(geometryRaw.mode || raw.spatial_mode || ''),
        coords_kind: plainText(geometryRaw.coords_kind || ''),
        x_range: xRange && xRange.length >= 2 ? [asNumber(xRange[0]), asNumber(xRange[1])] : null,
        y_range: yRange && yRange.length >= 2 ? [asNumber(yRange[0]), asNumber(yRange[1])] : null,
        image: typeof geometryRaw.image === 'string' ? geometryRaw.image : '',
        scale: asNumber(geometryRaw.scale),
        neighbor_radius: asNumber(geometryRaw.neighbor_radius),
      },
      spatial_mode: plainText(raw.spatial_mode || (raw.geometry && raw.geometry.mode) || ''),
      sample_columns: sampleColumns,
      steps: steps,
      current: plainText(raw.current || ''),
      cluster_keys: clusterKeys,
      ref: ref,
      warnings: normalizeStringList(raw.warnings),
      notes: normalizeStringList(raw.notes),
      history: history,
      analyses_run: runs,
      stats_contract: normalizeStatsContract(raw.stats_contract),
    }
  }

  function normalizePreview(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var coordsMeta = raw.coords_meta && typeof raw.coords_meta === 'object' ? raw.coords_meta : {}
    var columns = []
    if (Array.isArray(raw.sample_columns)) {
      for (var index = 0; index < raw.sample_columns.length; index++) {
        var item = raw.sample_columns[index] && typeof raw.sample_columns[index] === 'object' ? raw.sample_columns[index] : {}
        var key = plainText(item.key || item.name || '')
        if (!key) continue
        columns.push({ key: key, kind: plainText(item.kind || 'categorical'), n_unique: asNumber(item.n_unique) })
      }
    }
    var coords = []
    var sourceCoords = Array.isArray(raw.coords) ? raw.coords : (Array.isArray(raw.coords_sample) ? raw.coords_sample : [])
    for (var coordIndex = 0; coordIndex < sourceCoords.length; coordIndex++) {
      var point = pointXY(sourceCoords[coordIndex])
      if (point) coords.push([point.x, point.y])
    }
    return {
      genes: Array.isArray(raw.genes) ? raw.genes.map(function (gene) { return plainText(gene) }).filter(Boolean) : [],
      spots: Array.isArray(raw.spots) ? raw.spots.map(function (spot) { return plainText(spot) }).filter(Boolean) : [],
      n_spots: asNumber(raw.n_spots !== undefined ? raw.n_spots : raw.n_samples),
      n_genes: asNumber(raw.n_genes),
      columns: columns,
      coords_meta: {
        x_range: Array.isArray(coordsMeta.x_range) ? coordsMeta.x_range : null,
        y_range: Array.isArray(coordsMeta.y_range) ? coordsMeta.y_range : null,
        size: asNumber(coordsMeta.size),
      },
      coords: coords,
      steps: Array.isArray(raw.steps) ? raw.steps.map(function (step) { return plainText(step) }) : [],
      current: plainText(raw.current || ''),
      cluster_keys: Array.isArray(raw.cluster_keys) ? raw.cluster_keys.map(function (key) { return plainText(key) }) : [],
      spatial_mode: plainText(raw.spatial_mode || ''),
      image: typeof raw.image === 'string' ? raw.image : '',
    }
  }

  function normalizeHealth(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var engine = raw.engine && typeof raw.engine === 'object' ? raw.engine : raw
    var capabilities = engine.capabilities && typeof engine.capabilities === 'object' ? engine.capabilities : {}
    var packages = engine.packages && typeof engine.packages === 'object' ? engine.packages : {}
    var ready = engine.ready === true
    if (isNil(engine.ready) && (typeof engine.python === 'string' || typeof engine.mode === 'string')) ready = true
    var steps = Array.isArray(capabilities.steps)
      ? capabilities.steps.map(function (item) { return typeof item === 'string' ? item : plainText(item && item.id) }).filter(Boolean)
      : []
    return {
      ready: ready,
      mode: plainText(engine.mode || ''),
      python: plainText(engine.python || ''),
      platform: plainText(engine.platform || ''),
      version: plainText(engine.version || ''),
      workdir: plainText(raw.workdir || engine.workdir || ''),
      packages: packages,
      capabilities: {
        steps: steps,
        analyses: Array.isArray(capabilities.analyses) ? capabilities.analyses : [],
        plots: Array.isArray(capabilities.plots) ? capabilities.plots : [],
        formats: Array.isArray(capabilities.formats) ? capabilities.formats : [],
      },
      geo: raw.geo !== false,
      error: plainText(raw.error || engine.error || ''),
      lastProgress: engine.lastProgress && typeof engine.lastProgress === 'object' ? engine.lastProgress : null,
    }
  }

  /** 引擎徽章：未就绪时把可读原因带出来（顶栏显示）。 */
  function engineBadge(engine) {
    var health = engine && (engine.capabilities || engine.ready !== undefined) ? engine : normalizeHealth(engine)
    var detail = []
    if (health.python) detail.push('Python ' + health.python)
    if (health.mode) detail.push(health.mode)
    if (!health.ready) {
      var reason = health.error || '引擎尚未启动'
      return {
        ready: false,
        tone: health.error ? 'bad' : 'warn',
        labelKey: health.error ? 'engineBad' : 'engineBooting',
        detail: detail.join(' · '),
        reason: reason,
      }
    }
    return { ready: true, tone: 'ok', labelKey: 'engineReady', detail: detail.join(' · '), reason: '' }
  }

  function normalizeSuggestion(item, index) {
    var raw = item && typeof item === 'object' ? item : { title: String(item || '') }
    return {
      id: plainText(raw.id || ('suggestion-' + index)),
      title: plainText(raw.title || raw.detail || ''),
      detail: plainText(raw.detail || ''),
      action: plainText(raw.action || ''),
      step: plainText(raw.step || ''),
      params: raw.params && typeof raw.params === 'object' ? raw.params : null,
    }
  }

  function normalizeSuggestions(list) {
    if (!Array.isArray(list)) return []
    var out = []
    for (var index = 0; index < list.length; index++) out.push(normalizeSuggestion(list[index], index))
    return out
  }

  function normalizeActivityEntry(entry, index) {
    var raw = entry && typeof entry === 'object' ? entry : {}
    var ok = raw.ok !== false && raw.status !== 'failed' && raw.status !== 'error'
    var action = plainText(raw.action || raw.cmd || raw.step || raw.kind || raw.op || '')
    var text = plainText(raw.summary || raw.message || raw.detail || raw.error || (typeof entry === 'string' ? entry : ''), 300)
    return {
      id: plainText(raw.id || ('activity-' + index)),
      at: plainText(raw.at || raw.time || raw.ts || ''),
      actor: plainText(raw.actor || 'user'),
      action: action,
      text: text,
      ok: ok,
      ms: asNumber(raw.ms),
      dataset: plainText(raw.dataset || ''),
      step: plainText(raw.step || ''),
    }
  }

  function normalizeState(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var datasets = []
    if (Array.isArray(raw.datasets)) {
      for (var index = 0; index < raw.datasets.length; index++) datasets.push(normalizeDataset(raw.datasets[index]))
    }
    var activity = []
    if (Array.isArray(raw.activity)) {
      for (var activityIndex = 0; activityIndex < raw.activity.length; activityIndex++) {
        activity.push(normalizeActivityEntry(raw.activity[activityIndex], activityIndex))
      }
    }
    var active = plainText(raw.active || raw.active_id || '')
    if (!active && datasets.length) active = datasets[0].id
    return {
      active: active,
      datasets: datasets,
      engine: normalizeHealth(raw.engine ? { engine: raw.engine } : raw),
      suggestions: normalizeSuggestions(raw.suggestions),
      recipes: Array.isArray(raw.recipes) ? raw.recipes : [],
      steps: Array.isArray(raw.steps) ? raw.steps : [],
      activity: activity,
      workdir: plainText(raw.workdir || ''),
    }
  }

  /** /spatial/analyses 目录归一化：分组 + 参数 schema + 可用性。 */
  function normalizeAnalyses(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var list = Array.isArray(raw.analyses) ? raw.analyses : []
    var out = []
    for (var index = 0; index < list.length; index++) {
      var item = list[index] && typeof list[index] === 'object' ? list[index] : {}
      var params = []
      if (Array.isArray(item.params)) {
        for (var paramIndex = 0; paramIndex < item.params.length; paramIndex++) {
          var descriptor = item.params[paramIndex] && typeof item.params[paramIndex] === 'object' ? item.params[paramIndex] : {}
          if (!descriptor.key) continue
          params.push({
            key: plainText(descriptor.key),
            label: plainText(descriptor.label || descriptor.key),
            type: plainText(descriptor.type || 'string'),
            default: descriptor.default,
            min: asNumber(descriptor.min),
            max: asNumber(descriptor.max),
            step: asNumber(descriptor.step),
            placeholder: plainText(descriptor.placeholder || ''),
            options: Array.isArray(descriptor.options)
              ? descriptor.options.map(function (option) {
                  if (option && typeof option === 'object') return { value: plainText(option.value), label: plainText(option.label || option.value) }
                  return { value: plainText(option), label: plainText(option) }
                })
              : null,
          })
        }
      }
      out.push({
        id: plainText(item.id || ''),
        title: plainText(item.title || item.id || ''),
        group: plainText(item.group || '其它'),
        method: plainText(item.method || ''),
        needs: Array.isArray(item.needs) ? item.needs.map(function (need) { return plainText(need) }) : [],
        params: params,
        available: item.available !== false,
        reason: plainText(item.reason || ''),
        native: item.native !== false,
      })
    }
    var groups = Array.isArray(raw.groups) ? raw.groups.map(function (item) { return plainText(item) }).filter(Boolean) : []
    if (!groups.length) {
      var seen = {}
      for (var groupIndex = 0; groupIndex < out.length; groupIndex++) {
        if (!seen[out[groupIndex].group]) { seen[out[groupIndex].group] = true; groups.push(out[groupIndex].group) }
      }
    }
    return {
      analyses: out,
      groups: groups,
      diagnostics: raw.diagnostics && typeof raw.diagnostics === 'object' ? raw.diagnostics : {},
      columns: Array.isArray(raw.columns) ? raw.columns : [],
    }
  }

  /** audit 归一化：score 0–100 取整，level 三态。 */
  function normalizeAudit(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var checks = []
    var list = Array.isArray(raw.checks) ? raw.checks : []
    for (var index = 0; index < list.length; index++) {
      var item = list[index] && typeof list[index] === 'object' ? list[index] : {}
      checks.push({
        id: plainText(item.id || ('check-' + index)),
        level: item.level === 'fail' ? 'fail' : (item.level === 'warn' ? 'warn' : 'ok'),
        title: plainText(item.title || item.id || ''),
        detail: plainText(item.detail || ''),
        hint: plainText(item.hint || ''),
      })
    }
    var score = asNumber(raw.score)
    return {
      score: score === null ? 0 : Math.round(clamp(score, 0, 100)),
      level: raw.level === 'fail' ? 'fail' : (raw.level === 'warn' ? 'warn' : 'ok'),
      checks: checks,
    }
  }

  function normalizeJob(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    var job = raw.job && typeof raw.job === 'object' ? raw.job : raw
    var pct = asNumber(job.pct)
    return {
      id: plainText(job.id || ''),
      kind: plainText(job.kind || ''),
      status: plainText(job.status || ''),
      pct: pct === null ? 0 : clamp(pct, 0, 100),
      message: plainText(job.message || job.detail || ''),
      result: job.result && typeof job.result === 'object' ? job.result : null,
      error: plainText(job.error || ''),
    }
  }

  /**
   * /jobs 响应状态机：running | done | failed | timeout
   * timeout 由调用方传 elapsedMs/timeoutMs 判定；done 优先级最高（超时判定不覆盖已完成）。
   */
  function jobPolling(response, options) {
    var job = normalizeJob(response)
    var opts = options && typeof options === 'object' ? options : {}
    var elapsedMs = asNumber(opts.elapsedMs)
    if (elapsedMs === null || elapsedMs < 0) elapsedMs = 0
    var timeoutMs = asNumber(opts.timeoutMs)
    if (timeoutMs === null || timeoutMs <= 0) timeoutMs = POLL_TIMEOUT_MS
    var done = job.status === 'done' || job.status === 'success' || job.status === 'completed'
    var failed = job.status === 'failed' || job.status === 'error'
    var timedOut = !done && !failed && elapsedMs >= timeoutMs
    var status = done ? 'done' : (failed ? 'failed' : (timedOut ? 'timeout' : 'running'))
    return {
      id: job.id,
      kind: job.kind,
      status: status,
      pct: job.pct,
      message: job.message,
      result: job.result,
      error: job.error,
      done: done,
      failed: failed,
      timedOut: timedOut,
      running: status === 'running',
    }
  }

  function jobDone(response) {
    return jobPolling(response).done
  }

  function pickJobId(payload) {
    var raw = payload && typeof payload === 'object' ? payload : {}
    if (raw.job && typeof raw.job === 'object' && raw.job.id) return plainText(raw.job.id)
    if (raw.job_id) return plainText(raw.job_id)
    if (raw.id) return plainText(raw.id)
    return ''
  }

  // ==========================================================================
  // API 客户端：全部走同源 /spatial/*；fetch 可注入（单测用假 fetch）。
  // 每个方法只描述「一个操作」，HTTP 路由与 Agent 工具共用同一套实现。
  // ==========================================================================

  function createApi(fetchImpl, baseUrl) {
    var base = typeof baseUrl === 'string' && baseUrl ? baseUrl.replace(/\/$/, '') : API

    function resolveFetch() {
      if (typeof fetchImpl === 'function') return fetchImpl
      if (typeof globalThis.fetch === 'function') return globalThis.fetch.bind(globalThis)
      return null
    }

    function apiError(message, data, status) {
      var error = new Error(message)
      error.__spApiError = true
      error.code = ''
      if (data && typeof data === 'object') {
        error.data = data
        if (typeof data.code === 'string' && data.code) error.code = data.code
        if (Array.isArray(data.candidates)) error.candidates = data.candidates
        if (data.reason) error.reason = String(data.reason)
      }
      if (isFinite(Number(status)) && Number(status) > 0) error.status = Number(status)
      return error
    }

    function pickError(data) {
      if (!data || typeof data !== 'object') return ''
      if (typeof data.error === 'string' && data.error) return data.error
      if (typeof data.message === 'string' && data.message) return data.message
      return ''
    }

    function parseResponse(response) {
      if (!response) return Promise.reject(apiError('Host 返回了空响应'))
      var ok = response.ok !== false
      var status = Number(response.status) || 0
      var textPromise = typeof response.text === 'function'
        ? Promise.resolve().then(function () { return response.text() }).catch(function () { return '' })
        : Promise.resolve('')
      return textPromise.then(function (text) {
        var data = null
        if (text) {
          try { data = JSON.parse(text) } catch (ignored) { data = null }
        }
        if (!ok) throw apiError(pickError(data) || ('HTTP ' + status), data, status)
        if (data && data.ok === false) throw apiError(pickError(data) || '请求失败', data, status)
        return data && typeof data === 'object' ? data : {}
      })
    }

    function request(path, settings, extra) {
      var doFetch = resolveFetch()
      if (!doFetch) return Promise.reject(apiError('当前环境不支持 fetch'))
      var options = Object.assign({}, settings || {})
      var extraOptions = extra && typeof extra === 'object' ? extra : {}
      var timeoutMs = isFinite(Number(extraOptions.timeoutMs)) && Number(extraOptions.timeoutMs) > 0
        ? Number(extraOptions.timeoutMs)
        : DEFAULT_API_TIMEOUT_MS
      var controller = null
      var timer = null
      if (extraOptions.signal) {
        options.signal = extraOptions.signal
      } else if (typeof AbortController === 'function' && timeoutMs > 0) {
        controller = new AbortController()
        options.signal = controller.signal
        timer = setTimeout(function () {
          try { controller.abort() } catch (ignored) { /* ignore */ }
        }, timeoutMs)
      }
      function clear() {
        if (timer) { clearTimeout(timer); timer = null }
      }
      var promise = Promise.resolve()
        .then(function () { return doFetch(base + path, options) })
        .then(parseResponse)
      return promise.then(function (value) {
        clear()
        return value
      }, function (error) {
        clear()
        if (error && error.__spApiError) throw error
        var message = formatError(error)
        if (error && (error.name === 'AbortError' || /aborted/i.test(message))) {
          if (extraOptions.signal) throw error
          throw apiError('请求超时（' + Math.round(timeoutMs / 1000) + ' 秒）：引擎可能还在计算，请稍后刷新重试。')
        }
        throw apiError('网络请求失败：' + (message || '无法连接 Host'))
      })
    }

    function jsonPost(path, body, extra) {
      return request(path, {
        method: 'POST',
        headers: { 'content-type': 'application/json' },
        body: JSON.stringify(body || {}),
      }, extra)
    }

    function withActor(body) {
      return Object.assign({ actor: 'user' }, body && typeof body === 'object' ? body : {})
    }

    /** 主路由不存在时（老 Host）回落到兼容路由，只在该 code 下重试一次。 */
    function requestWithFallback(primaryPath, fallbackPath, settings, extra) {
      return request(primaryPath, settings, extra).catch(function (error) {
        var notFound = error && (error.status === 404 || error.code === 'not_found')
        if (!fallbackPath || !notFound) throw error
        return request(fallbackPath, settings, extra)
      })
    }

    return {
      base: base,
      request: request,
      getHealth: function (extra) { return request('/health', { method: 'GET' }, extra) },
      getState: function (extra) { return request('/state', { method: 'GET' }, extra) },
      getActivity: function (extra) { return request('/activity', { method: 'GET' }, extra) },
      getFormats: function (extra) { return request('/formats', { method: 'GET' }, extra) },
      getAnalyses: function (dataset, extra) {
        var query = dataset ? '?dataset=' + encodeURIComponent(plainText(dataset)) : ''
        return request('/analyses' + query, { method: 'GET' }, extra)
      },
      getJob: function (id, extra) {
        return request('/jobs?id=' + encodeURIComponent(plainText(id)), { method: 'GET' }, extra)
      },
      geoSearch: function (filters, extra) {
        return jsonPost('/geo/search', buildSearchPayload(filters), Object.assign({ timeoutMs: 120000 }, extra || {}))
      },
      geoDetail: function (gse, extra) {
        return jsonPost('/geo/detail', withActor({ gse: plainText(gse).toUpperCase() }), Object.assign({ timeoutMs: 180000 }, extra || {}))
      },
      geoDownload: function (options, extra) {
        var body = {
          gse: plainText(options && options.gse).toUpperCase(),
          files: Array.isArray(options && options.files) ? options.files.slice() : [],
          dest_dir: plainText(options && options.dest_dir) || DEFAULT_GEO_DIR,
          overwrite: options && options.overwrite === true,
        }
        if (options && options.open_after === true) body.open_after = true
        return jsonPost('/geo/download', withActor(body), Object.assign({ timeoutMs: 90000 }, extra || {}))
      },
      pick: function (kind, title, extra) {
        return jsonPost('/pick', withActor({
          kind: plainText(kind) || 'directory',
          title: plainText(title),
        }), Object.assign({ timeoutMs: 120000 }, extra || {}))
      },
      pickDir: function (title, extra) { return this.pick('directory', title, extra) },
      openDataset: function (options, extra) {
        return jsonPost('/dataset/open', withActor(buildImportPayload(options)),
          Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      demo: function (extra) { return jsonPost('/dataset/demo', withActor({}), extra) },
      selectDataset: function (dataset, extra) {
        return jsonPost('/dataset/select', withActor({ dataset: plainText(dataset) }), extra)
      },
      closeDataset: function (dataset, extra) {
        return jsonPost('/dataset/close', withActor({ dataset: plainText(dataset) }), extra)
      },
      inspectSource: function (options, extra) {
        return jsonPost('/inspect', withActor(buildImportPayload(options)),
          Object.assign({ timeoutMs: 180000 }, extra || {}))
      },
      preview: function (dataset, extra) {
        return jsonPost('/preview', withActor({ dataset: plainText(dataset) }), extra)
      },
      samples: function (payload, extra) {
        return jsonPost('/samples', withActor(payload && typeof payload === 'object' ? payload : {}), extra)
      },
      runStep: function (dataset, step, params, extra) {
        return jsonPost('/run', withActor({
          dataset: plainText(dataset),
          step: plainText(step),
          params: params && typeof params === 'object' ? params : {},
        }), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      plot: function (dataset, kind, params, extra) {
        return jsonPost('/plot', withActor({
          dataset: plainText(dataset),
          kind: plainText(kind),
          params: params && typeof params === 'object' ? params : {},
        }), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      analyses: function (dataset, extra) {
        var query = dataset ? '?dataset=' + encodeURIComponent(plainText(dataset)) : ''
        return request('/analyses' + query, { method: 'GET' }, extra)
      },
      analysis: function (dataset, id, params, extra) {
        return jsonPost('/analysis', withActor({
          dataset: plainText(dataset),
          id: plainText(id),
          params: params && typeof params === 'object' ? params : {},
        }), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS * 2 }, extra || {}))
      },
      exportData: function (dataset, what, path, extra) {
        var body = { dataset: plainText(dataset), what: plainText(what) }
        if (path) body.path = plainText(path)
        return jsonPost('/export', withActor(body), Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      audit: function (dataset, extra) {
        return jsonPost('/audit', withActor({ dataset: plainText(dataset) }), extra)
      },
      suggest: function (dataset, extra) {
        return jsonPost('/suggest', withActor({ dataset: plainText(dataset) }), extra)
      },
      refs: function (extra) {
        return requestWithFallback('/refs', '/sc/refs', { method: 'GET' }, extra)
      },
      refAttach: function (payload, extra) {
        var body = withActor({
          dataset: plainText(payload && payload.dataset),
          ref_path: plainText(payload && payload.ref_path),
          sc_dataset: plainText(payload && payload.sc_dataset),
          cell_type_col: plainText(payload && payload.cell_type_col),
        })
        if (!body.sc_dataset) delete body.sc_dataset
        return requestWithFallback('/ref/attach', '/sc/attach', {
          method: 'POST',
          headers: { 'content-type': 'application/json' },
          body: JSON.stringify(body),
        }, Object.assign({ timeoutMs: HEAVY_API_TIMEOUT_MS }, extra || {}))
      },
      refSet: function (dataset, refPath, cellTypeCol, extra) {
        return jsonPost('/ref/set', withActor({
          dataset: plainText(dataset),
          ref_path: plainText(refPath),
          cell_type_col: plainText(cellTypeCol),
        }), extra)
      },
      reset: function (extra) { return jsonPost('/reset', withActor({}), extra) },
      restartEngine: function (extra) { return jsonPost('/engine/restart', withActor({}), extra) },
    }
  }

  // --------------------------------------------------------------------------
  // 展示小工具：基因检索 / 表格排序 / 图缓存键 / 取词
  // --------------------------------------------------------------------------

  function matchGenes(allGenes, query, limit) {
    var source = Array.isArray(allGenes) ? allGenes : []
    var needle = plainText(query).toUpperCase()
    var max = clampInt(limit, 1, 500, 50)
    var prefix = []
    var contains = []
    for (var index = 0; index < source.length; index++) {
      var gene = String(source[index] || '')
      var upper = gene.toUpperCase()
      if (!needle || upper === needle) { prefix.push(gene); continue }
      if (upper.indexOf(needle) === 0) prefix.push(gene)
      else if (upper.indexOf(needle) > 0) contains.push(gene)
      if (prefix.length >= max) break
    }
    return prefix.concat(contains).slice(0, max)
  }

  function tableColumns(table) {
    return table && Array.isArray(table.columns) ? table.columns : []
  }

  function sortRows(rows, key, direction, columns) {
    var list = Array.isArray(rows) ? rows.slice() : []
    var numeric = false
    if (Array.isArray(columns)) {
      for (var index = 0; index < columns.length; index++) {
        if (columns[index] && columns[index].key === key && columns[index].type === 'number') numeric = true
      }
    }
    var factor = direction === 'desc' ? -1 : 1
    return list.sort(function (a, b) {
      var left = a && typeof a === 'object' ? a[key] : null
      var right = b && typeof b === 'object' ? b[key] : null
      if (isNil(left) && isNil(right)) return 0
      if (isNil(left)) return 1
      if (isNil(right)) return -1
      if (numeric) {
        var leftNumber = asNumber(left)
        var rightNumber = asNumber(right)
        if (leftNumber !== null && rightNumber !== null) return (leftNumber - rightNumber) * factor
      }
      return String(left).localeCompare(String(right)) * factor
    })
  }

  function formatCell(value, numeric) {
    if (isNil(value) || value === '') return '—'
    if (typeof value === 'number') {
      if (!numeric) return String(value)
      var rounded = Math.abs(value) >= 1000 ? value : Number(value.toPrecision(4))
      return String(rounded)
    }
    if (typeof value === 'boolean') return value ? '是' : '否'
    if (Array.isArray(value)) return value.join(', ')
    return String(value)
  }

  function formatTableRows(table, limit) {
    var columns = tableColumns(table)
    var rows = table && Array.isArray(table.rows) ? table.rows : []
    var max = clampInt(limit, 1, 5000, MAX_TABLE_RENDER_ROWS)
    var shown = rows.slice(0, max)
    var out = []
    for (var index = 0; index < shown.length; index++) {
      var row = shown[index] && typeof shown[index] === 'object' ? shown[index] : {}
      var cells = []
      for (var columnIndex = 0; columnIndex < columns.length; columnIndex++) {
        var column = columns[columnIndex]
        cells.push({
          key: column.key,
          text: formatCell(row[column.key], column.type === 'number'),
          numeric: column.type === 'number',
        })
      }
      out.push({ key: 'row-' + index, cells: cells })
    }
    return {
      columns: columns,
      rows: out,
      total: rows.length,
      shown: out.length,
      truncated: rows.length > out.length || (table && table.truncated === true),
    }
  }

  /** PlotBlock 缓存键：kind + 稳定序列化参数（键排序、跳过空值）。 */
  function plotKey(kind, params) {
    var name = plainText(kind)
    var source = params && typeof params === 'object' && !Array.isArray(params) ? params : {}
    var keys = Object.keys(source).sort()
    var parts = []
    for (var index = 0; index < keys.length; index++) {
      var key = keys[index]
      var value = source[key]
      if (isNil(value) || value === '' || value === false) continue
      if (Array.isArray(value)) {
        var joined = value.map(function (item) { return plainText(item) }).filter(Boolean).join(',')
        if (joined) parts.push(key + '=' + joined)
        continue
      }
      if (typeof value === 'object') {
        try { parts.push(key + '=' + JSON.stringify(value)) } catch (ignored) { /* 不可序列化就跳过 */ }
        continue
      }
      parts.push(key + '=' + String(value))
    }
    return name + (parts.length ? '|' + parts.join('&') : '')
  }

  /**
   * 宿主 locale 服务 → 取词函数。宿主不可用时静默回落 DICT_ZH，
   * 保证任何环境（含单测）下 T() 都不会返回裸 key 或抛异常。
   */
  function createTranslator(locale) {
    var bound = null
    try {
      if (locale && typeof locale.bind === 'function') bound = locale.bind(NS)
    } catch (ignored) {
      bound = null
    }
    return function translate(key) {
      if (typeof bound === 'function') {
        try {
          var value = bound(key)
          if (typeof value === 'string' && value && value !== key) return value
        } catch (ignored) { /* 宿主 locale 不可用时静默回落 */ }
      }
      var fallback = DICT_ZH[key]
      return typeof fallback === 'string' ? fallback : String(key)
    }
  }

  function analysisGroupTab(group) {
    var name = plainText(group)
    return ANALYSIS_TAB[name] || 'stats'
  }

  /** 数据集当前聚类列（cluster_keys 首项，缺省 clusters）。 */
  function currentClusterKey(dataset) {
    var keys = dataset && Array.isArray(dataset.cluster_keys) ? dataset.cluster_keys : []
    if (keys.length) return keys[0]
    return 'clusters'
  }

  /** 有 domains 列时返回 true（决定空间图是否提供「空间域」着色）。 */
  function hasDomainKey(dataset) {
    var keys = dataset && Array.isArray(dataset.cluster_keys) ? dataset.cluster_keys : []
    for (var index = 0; index < keys.length; index++) {
      if (/domain/i.test(keys[index])) return true
    }
    var steps = dataset && Array.isArray(dataset.steps) ? dataset.steps : []
    for (var stepIndex = 0; stepIndex < steps.length; stepIndex++) {
      if (steps[stepIndex] === 'domains') return true
    }
    return false
  }

  // ==========================================================================
  // 纯函数测试出口。node 环境（typeof window === 'undefined'）挂到
  // globalThis.__SP_TEST__；浏览器里也会挂进 factory 返回的 __test__，便于联调。
  // ==========================================================================
  var SP_TEST = {
    NS: NS,
    API: API,
    POLL_MS: POLL_MS,
    POLL_TIMEOUT_MS: POLL_TIMEOUT_MS,
    DEFAULT_GEO_DIR: DEFAULT_GEO_DIR,
    MAX_HISTORY_ROWS: MAX_HISTORY_ROWS,
    DEFAULT_API_TIMEOUT_MS: DEFAULT_API_TIMEOUT_MS,
    STEPS: STEPS,
    STEP_BY_ID: STEP_BY_ID,
    TAB_KEYS: TAB_KEYS,
    TAB_META: TAB_META,
    TAB_FOR_STEP: TAB_FOR_STEP,
    ANALYSIS_TAB: ANALYSIS_TAB,
    PALETTE: PALETTE,
    CONTINUOUS_STOPS: CONTINUOUS_STOPS,
    PLOT_KINDS: PLOT_KINDS,
    EXPORT_KINDS: EXPORT_KINDS,
    QUICK_EXPORTS: QUICK_EXPORTS,
    COLOR_BY: COLOR_BY,
    GEO_PLATFORMS: GEO_PLATFORMS,
    DATA_FORMATS: DATA_FORMATS,
    AUDIT_LEVEL_STYLE: AUDIT_LEVEL_STYLE,
    DICT_ZH: DICT_ZH,
    DICT_EN: DICT_EN,
    clamp: clamp,
    asNumber: asNumber,
    plainText: plainText,
    fmtNumber: fmtNumber,
    fmtPct: fmtPct,
    fmtP: fmtP,
    fmtBytes: formatBytes,
    formatBytes: formatBytes,
    formatDuration: formatDuration,
    formatClock: formatClock,
    formatError: formatError,
    hexToRgb: hexToRgb,
    colorForGroup: colorForGroup,
    continuousColor: continuousColor,
    groupBy: groupBy,
    extentOf: extentOf,
    normalizeExtent: normalizeExtent,
    spatialScreen: spatialScreen,
    buildSearchPayload: buildSearchPayload,
    buildImportPayload: buildImportPayload,
    suggestionInstruction: suggestionInstruction,
    defaultParams: defaultParams,
    collectStepParams: collectStepParams,
    parseTokenList: parseTokenList,
    formatParamValue: formatParamValue,
    doneStepSet: doneStepSet,
    isStepReady: isStepReady,
    stepStatusOf: stepStatusOf,
    stepStatusKey: stepStatusKey,
    stepById: stepById,
    figureSrc: figureSrc,
    normalizeFigure: normalizeFigure,
    normalizeTable: normalizeTable,
    normalizeResult: normalizeResult,
    normalizeDataset: normalizeDataset,
    normalizeStatsContract: normalizeStatsContract,
    normalizePreview: normalizePreview,
    normalizeHealth: normalizeHealth,
    normalizeState: normalizeState,
    normalizeSuggestions: normalizeSuggestions,
    normalizeActivityEntry: normalizeActivityEntry,
    normalizeAnalyses: normalizeAnalyses,
    normalizeAudit: normalizeAudit,
    normalizeJob: normalizeJob,
    normalizeStringList: normalizeStringList,
    engineBadge: engineBadge,
    jobPolling: jobPolling,
    jobDone: jobDone,
    pickJobId: pickJobId,
    createApi: createApi,
    apiClient: createApi,
    createTranslator: createTranslator,
    matchGenes: matchGenes,
    sortRows: sortRows,
    formatTableRows: formatTableRows,
    tableColumns: tableColumns,
    formatCell: formatCell,
    plotKey: plotKey,
    analysisGroupTab: analysisGroupTab,
    currentClusterKey: currentClusterKey,
    hasDomainKey: hasDomainKey,
  }
  // ==========================================================================
  // 浏览器层：React 组件 + 样式 + 槽注册。只 require('react')。
  // ==========================================================================

  function factory(require) {
    var React = require('react')
    var useState = React.useState
    var useEffect = React.useEffect
    var useRef = React.useRef
    var useMemo = React.useMemo

    /** React.createElement 简写（不写 JSX，避免宿主构建链依赖）。 */
    function h(type, props) {
      var children = Array.prototype.slice.call(arguments, 2)
      return React.createElement.apply(React, [type, props].concat(children))
    }

    /** t()：apply() 里用宿主 locale 服务重新绑定；默认回落中文词典。 */
    var T = function (key) {
      var fallback = DICT_ZH[key]
      return typeof fallback === 'string' ? fallback : String(key)
    }

    function warnPlugin(message, error) {
      try { console.warn('[dsh-spatial] ' + message, error || '') } catch (ignored) { /* ignore */ }
    }

    /** 模板取词：字典里用 `{v}` 占位，避免依赖宿主 t() 的插值实现。 */
    function tf(translate, key, value) {
      var raw = typeof translate === 'function' ? translate(key) : ''
      var template = typeof raw === 'string' && raw && raw !== key ? raw : (DICT_ZH[key] || '')
      return String(template).replace(/\{v\}/g, isNil(value) ? '' : String(value))
    }

    // ==========================================================================
    // 样式：一个 CSS 字符串，由 <style> 渲染在组件内部（生命周期跟随组件）。
    // 颜色只用 --dsw-alias-* 主题变量，派生色统一走 color-mix，深浅色自适应。
    // 类名一律 sp- 前缀；数据色（PALETTE / CONTINUOUS_STOPS）走内联样式。
    // ==========================================================================
    var CSS = [
      '.sp-root{display:flex;flex-direction:column;height:100%;min-height:72vh;overflow:hidden;position:relative;isolation:isolate;color:var(--dsw-alias-label-primary);background:var(--dsw-alias-bg-base);border-radius:12px;font-size:13px;line-height:1.5;font-variant-numeric:tabular-nums}',
      '.sp-root *,.sp-root *::before,.sp-root *::after{box-sizing:border-box}',
      '.sp-root button,.sp-root input,.sp-root select,.sp-root textarea{font-family:inherit;font-variant-numeric:tabular-nums}',
      '.sp-root ::selection{background:color-mix(in srgb,var(--dsw-alias-brand-primary) 24%,transparent)}',

      /* 顶栏 + 1px 渐变线 */
      '.sp-top{position:relative;display:flex;align-items:center;gap:10px;flex-wrap:wrap;flex:none;min-height:52px;padding:9px 14px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sp-top::after{content:"";position:absolute;left:0;right:0;bottom:-1px;height:1px;background:linear-gradient(90deg,transparent,color-mix(in srgb,var(--dsw-alias-brand-primary) 55%,transparent) 26%,color-mix(in srgb,var(--dsw-alias-state-success-primary) 45%,transparent) 58%,transparent);pointer-events:none}',
      '.sp-brand{display:flex;align-items:center;gap:9px;min-width:0}',
      '.sp-brand-mark{display:inline-flex;align-items:center;justify-content:center;flex:none;width:28px;height:28px;border-radius:9px;border:1px solid color-mix(in srgb,var(--dsw-alias-brand-primary) 34%,var(--dsw-alias-border-l1));background:linear-gradient(135deg,color-mix(in srgb,var(--dsw-alias-brand-primary) 16%,var(--dsw-alias-bg-layer-1)),var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary)}',
      '.sp-title-block{display:flex;flex-direction:column;gap:1px;min-width:0}',
      '.sp-title{font-size:14px;font-weight:600;letter-spacing:.01em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-subtitle{font-size:11px;color:var(--dsw-alias-label-secondary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-badge{display:inline-flex;align-items:center;gap:6px;max-width:340px;height:22px;padding:0 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font-size:11px;line-height:1;cursor:default;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-badge-dot{width:6px;height:6px;border-radius:50%;flex:none;background:currentColor;box-shadow:0 0 0 2px color-mix(in srgb,currentColor 18%,transparent)}',
      '.sp-badge-ok{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 40%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 9%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-success-primary)}',
      '.sp-badge-warn{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 42%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 9%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-warn-primary)}',
      '.sp-badge-bad{border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 42%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 9%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-error-primary)}',
      '.sp-badge-reason{max-width:200px;overflow:hidden;text-overflow:ellipsis;color:var(--dsw-alias-state-error-primary);font-size:10.5px}',
      '.sp-contract{display:flex;flex-direction:column;gap:5px}',
      '.sp-contract-hint{font-size:11px;color:var(--dsw-alias-label-secondary);line-height:1.5}',
      '.sp-contract-row{display:flex;gap:8px;font-size:12px;align-items:baseline}',
      '.sp-contract-k{flex:0 0 88px;color:var(--dsw-alias-label-secondary)}',
      '.sp-contract-v{flex:1 1 auto;color:var(--dsw-alias-label-primary);word-break:break-word}',
      '.sp-contract-warn{margin-top:4px;padding:7px 9px;border-radius:8px;font-size:11.5px;line-height:1.5;background:var(--dsw-alias-bg-layer-2);border:1px solid var(--dsw-alias-border-l1)}',
      '.sp-contract-warn-title{font-weight:600;margin-bottom:2px;color:var(--dsw-alias-state-warn-primary)}',
      '.sp-contract-warn-item{color:var(--dsw-alias-label-primary)}',
      '.sp-contract-info{margin-top:4px;font-size:11px;line-height:1.5;color:var(--dsw-alias-label-secondary)}',
      '.sp-chips{display:flex;align-items:center;gap:5px;flex-wrap:wrap;min-width:0;flex:1 1 200px}',
      '.sp-chip{display:inline-flex;align-items:center;gap:5px;max-width:210px;height:23px;padding:0 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font:inherit;font-size:11.5px;cursor:pointer;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;transition:border-color .14s ease,background .14s ease,color .14s ease}',
      '.sp-chip:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);color:var(--dsw-alias-label-primary)}',
      '.sp-chip-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 45%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 9%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary)}',
      '.sp-chip-label{border:none;background:transparent;color:inherit;font:inherit;font-size:11.5px;cursor:pointer;padding:0;max-width:170px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sp-chip-x{display:inline-flex;align-items:center;justify-content:center;width:14px;height:14px;padding:0;border:none;border-radius:50%;background:transparent;font-size:10px;color:var(--dsw-alias-label-tertiary);cursor:pointer}',
      '.sp-chip-x:hover{background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 18%,transparent);color:var(--dsw-alias-state-error-primary)}',
      '.sp-top-actions{display:flex;align-items:center;gap:6px;flex-wrap:wrap;margin-left:auto;justify-content:flex-end}',

      /* 按钮 / 输入 */
      '.sp-btn{display:inline-flex;align-items:center;justify-content:center;gap:5px;height:28px;padding:0 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-primary);font:inherit;font-size:12px;line-height:1;cursor:pointer;white-space:nowrap;transition:border-color .14s ease,background .14s ease,opacity .14s ease,box-shadow .14s ease}',
      '.sp-btn:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-overlay)}',
      '.sp-btn:active:not(:disabled){transform:translateY(.5px)}',
      '.sp-btn:disabled{opacity:.45;cursor:default}',
      '.sp-btn-primary{background:var(--dsw-alias-brand-primary);border-color:var(--dsw-alias-brand-primary);color:var(--dsw-alias-bg-base);box-shadow:0 1px 2px color-mix(in srgb,var(--dsw-alias-brand-primary) 32%,transparent)}',
      '.sp-btn-primary:hover:not(:disabled){opacity:.9;background:var(--dsw-alias-brand-primary);border-color:var(--dsw-alias-brand-primary)}',
      '.sp-btn-ghost{background:transparent;border-color:transparent;color:var(--dsw-alias-label-secondary)}',
      '.sp-btn-ghost:hover:not(:disabled){background:var(--dsw-alias-bg-layer-2);border-color:var(--dsw-alias-border-l1);color:var(--dsw-alias-label-primary)}',
      '.sp-btn-sm{height:24px;padding:0 8px;font-size:11.5px;border-radius:6px}',
      '.sp-btn-danger{color:var(--dsw-alias-state-error-primary);border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 36%,var(--dsw-alias-border-l1))}',
      '.sp-input{height:28px;padding:0 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base);color:var(--dsw-alias-label-primary);font:inherit;font-size:12px;outline:none;min-width:0;transition:border-color .14s ease,box-shadow .14s ease}',
      '.sp-input:focus{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 48%,var(--dsw-alias-border-l1));box-shadow:0 0 0 2px color-mix(in srgb,var(--dsw-alias-brand-primary) 12%,transparent)}',
      '.sp-input:disabled{opacity:.55}',
      '.sp-input::placeholder{color:var(--dsw-alias-label-tertiary)}',
      '.sp-select{cursor:pointer;padding-right:5px}',
      '.sp-textarea{min-height:54px;padding:6px 9px;resize:vertical;line-height:1.45}',
      '.sp-check{display:inline-flex;align-items:center;gap:6px;font-size:12px;color:var(--dsw-alias-label-secondary);cursor:pointer;user-select:none}',
      '.sp-check input{accent-color:var(--dsw-alias-brand-primary);width:14px;height:14px;cursor:pointer}',
      '.sp-form{display:flex;flex-direction:column;gap:7px}',
      '.sp-field{display:flex;align-items:center;justify-content:space-between;gap:8px;min-width:0}',
      '.sp-field-col{flex-direction:column;align-items:stretch;gap:4px}',
      '.sp-field-label{flex:0 1 auto;min-width:0;font-size:11.5px;color:var(--dsw-alias-label-secondary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-field .sp-input{flex:1 1 auto;max-width:152px;height:25px;font-size:11.5px}',
      '.sp-field-col .sp-input{max-width:none;width:100%}',
      '.sp-multi{display:flex;flex-wrap:wrap;gap:4px}',
      '.sp-multi-opt{display:inline-flex;align-items:center;gap:4px;height:22px;padding:0 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:6px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font:inherit;font-size:11px;cursor:pointer;transition:border-color .14s ease,background .14s ease,color .14s ease}',
      '.sp-multi-opt:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);color:var(--dsw-alias-label-primary)}',
      '.sp-multi-opt-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 46%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 10%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary)}',

      /* 三栏骨架 */
      '.sp-body{display:flex;flex:1 1 auto;min-height:0;overflow:hidden}',
      '.sp-left{display:flex;flex-direction:column;flex:0 0 268px;width:268px;min-height:0;border-right:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sp-center{display:flex;flex-direction:column;flex:1 1 auto;min-width:0;min-height:0;background:var(--dsw-alias-bg-base)}',
      '.sp-right{display:flex;flex-direction:column;flex:0 0 300px;width:300px;min-height:0;overflow-y:auto;border-left:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sp-panel-head{display:flex;align-items:center;justify-content:space-between;gap:8px;flex:none;padding:9px 12px;border-bottom:1px solid var(--dsw-alias-border-l1);font-size:11px;font-weight:600;letter-spacing:.04em;color:var(--dsw-alias-label-secondary)}',
      '.sp-panel-note{font-weight:400;font-size:10.5px;letter-spacing:0;color:var(--dsw-alias-label-tertiary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-panel-body{display:flex;flex-direction:column;gap:9px;padding:10px 12px 16px;overflow-y:auto;min-height:0}',
      '.sp-banner{display:flex;align-items:center;justify-content:space-between;gap:10px;flex:none;padding:6px 14px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-state-warn-primary) 32%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 7%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-label-primary);font-size:11.5px;animation:sp-fade .16s ease}',
      '.sp-banner-bad{border-bottom-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 34%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 6%,var(--dsw-alias-bg-layer-1))}',
      '.sp-banner-text{flex:1 1 auto;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',

      /* 左栏：数据集卡 + 12 张流程卡 */
      '.sp-dataset-card{display:flex;flex-direction:column;gap:5px;padding:10px 11px;border:1px solid var(--dsw-alias-border-l1);border-radius:12px;background:linear-gradient(160deg,color-mix(in srgb,var(--dsw-alias-brand-primary) 5%,var(--dsw-alias-bg-base)),var(--dsw-alias-bg-base) 58%);box-shadow:0 1px 2px color-mix(in srgb,var(--dsw-alias-label-primary) 5%,transparent)}',
      '.sp-dataset-name{font-size:12.5px;font-weight:600;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-dataset-meta{display:flex;flex-wrap:wrap;gap:4px 8px;font-size:11px;color:var(--dsw-alias-label-tertiary)}',
      '.sp-dataset-meta b{font-weight:600;color:var(--dsw-alias-label-secondary)}',
      '.sp-steps{display:flex;flex-direction:column;gap:8px;padding:10px 12px 16px;overflow-y:auto;min-height:0}',
      '.sp-step{display:flex;flex-direction:column;gap:5px;padding:8px 10px;border:1px solid var(--dsw-alias-border-l1);border-radius:12px;background:var(--dsw-alias-bg-base);box-shadow:0 1px 2px color-mix(in srgb,var(--dsw-alias-label-primary) 4%,transparent);transition:border-color .16s ease,background .16s ease,opacity .16s ease}',
      '.sp-step-done{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 28%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 4%,var(--dsw-alias-bg-base))}',
      '.sp-step-running{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 48%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 5%,var(--dsw-alias-bg-base))}',
      '.sp-step-blocked{opacity:.8}',
      '.sp-step-head{display:flex;align-items:center;gap:7px;min-width:0}',
      '.sp-step-dot{position:relative;display:inline-flex;align-items:center;justify-content:center;flex:none;width:20px;height:20px;border-radius:50%;border:1px solid var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font-size:10.5px;font-weight:600}',
      '.sp-step-dot-done{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 55%,transparent);background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 14%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-state-success-primary)}',
      '.sp-step-dot-running{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 58%,transparent);background:color-mix(in srgb,var(--dsw-alias-brand-primary) 14%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary);animation:sp-pulse 1.6s ease-in-out infinite}',
      '.sp-step-dot-pending{border-color:var(--dsw-alias-border-l2);color:var(--dsw-alias-label-tertiary)}',
      '.sp-step-dot-blocked{border-color:var(--dsw-alias-border-l1);color:var(--dsw-alias-label-tertiary)}',
      '.sp-step-title{flex:1 1 auto;min-width:0;font-size:12.5px;font-weight:500;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-step-status{flex:none;font-size:10.5px;color:var(--dsw-alias-label-tertiary);white-space:nowrap}',
      '.sp-status-done{color:var(--dsw-alias-state-success-primary)}',
      '.sp-status-running{color:var(--dsw-alias-brand-primary)}',
      '.sp-status-failed{color:var(--dsw-alias-state-error-primary)}',
      '.sp-status-ready{color:var(--dsw-alias-label-secondary)}',
      '.sp-step-caret{flex:none;display:inline-flex;align-items:center;justify-content:center;width:18px;height:18px;padding:0;border:none;border-radius:6px;background:transparent;color:var(--dsw-alias-label-tertiary);cursor:pointer;transition:background .14s ease,color .14s ease}',
      '.sp-step-caret:hover{background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-primary)}',
      '.sp-step-params{display:flex;flex-direction:column;gap:6px;padding:7px 1px 2px;border-top:1px dashed var(--dsw-alias-border-l1);margin-top:2px}',
      '.sp-step-foot{display:flex;align-items:center;gap:6px;flex-wrap:wrap}',
      '.sp-step-summary{flex:1 1 auto;min-width:0;font-size:11px;color:var(--dsw-alias-label-tertiary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',

      /* 中栏 tabs */
      '.sp-tabs{display:flex;align-items:center;gap:2px;flex:none;padding:6px 10px 0;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1);overflow-x:auto;scrollbar-width:thin}',
      '.sp-tab{display:inline-flex;align-items:center;justify-content:center;gap:5px;height:30px;padding:0 11px;border:none;border-bottom:2px solid transparent;border-radius:7px 7px 0 0;background:transparent;color:var(--dsw-alias-label-secondary);font:inherit;font-size:12px;cursor:pointer;white-space:nowrap;transition:color .14s ease,background .14s ease,border-color .14s ease}',
      '.sp-tab:hover{color:var(--dsw-alias-label-primary);background:var(--dsw-alias-bg-layer-2)}',
      '.sp-tab-on{color:var(--dsw-alias-brand-primary);border-bottom-color:var(--dsw-alias-brand-primary);font-weight:600}',
      '.sp-pane{display:flex;flex:1 1 auto;flex-direction:column;gap:11px;min-height:0;padding:12px 14px 18px;overflow-y:auto}',
      '.sp-toolbar{display:flex;align-items:center;gap:8px;flex-wrap:wrap;flex:none;min-height:28px}',
      '.sp-toolbar-label{font-size:11px;color:var(--dsw-alias-label-tertiary);white-space:nowrap}',
      '.sp-spacer{flex:1 1 auto}',

      /* 卡片 / 指标 / 区块 */
      '.sp-card{display:flex;flex-direction:column;gap:8px;padding:11px 13px;border:1px solid var(--dsw-alias-border-l1);border-radius:12px;background:var(--dsw-alias-bg-layer-1);box-shadow:0 1px 2px color-mix(in srgb,var(--dsw-alias-label-primary) 4%,transparent)}',
      '.sp-card-title{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:11.5px;font-weight:600;color:var(--dsw-alias-label-secondary);letter-spacing:.02em}',
      '.sp-metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(104px,1fr));gap:8px}',
      '.sp-metric{display:flex;flex-direction:column;gap:1px;padding:8px 11px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.sp-metric-k{font-size:10.5px;color:var(--dsw-alias-label-tertiary);letter-spacing:.03em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-metric-v{font-size:18px;font-weight:300;line-height:1.25;color:var(--dsw-alias-label-primary);font-variant-numeric:tabular-nums;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-grid-2{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:11px;align-items:start}',
      '.sp-section-title{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:11px;font-weight:600;letter-spacing:.04em;color:var(--dsw-alias-label-secondary);margin-top:2px}',
      '.sp-muted{font-size:11.5px;color:var(--dsw-alias-label-secondary);line-height:1.5}',
      '.sp-inline-row{display:flex;align-items:center;gap:6px;flex-wrap:wrap;min-width:0}',
      '.sp-inline-err{padding:5px 8px;border-radius:6px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-error-primary) 40%,transparent);background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 8%,transparent);color:var(--dsw-alias-state-error-primary);font-size:11.5px;word-break:break-word}',
      '.sp-warn-text{color:var(--dsw-alias-state-warn-primary);font-size:11.5px}',
      '.sp-ok-text{color:var(--dsw-alias-state-success-primary);font-size:11.5px}',
      '.sp-divider{flex:none;height:1px;margin:2px 0;background:var(--dsw-alias-border-l1)}',
      '.sp-path{font-size:10.5px;color:var(--dsw-alias-label-tertiary);word-break:break-all;font-family:ui-monospace,SFMono-Regular,Menlo,monospace}',
      '.sp-tag{display:inline-flex;align-items:center;height:18px;padding:0 7px;border:1px solid var(--dsw-alias-border-l1);border-radius:999px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font-size:10.5px;white-space:nowrap}',
      '.sp-tag-ok{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 40%,var(--dsw-alias-border-l1));color:var(--dsw-alias-state-success-primary)}',
      '.sp-tag-warn{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 40%,var(--dsw-alias-border-l1));color:var(--dsw-alias-state-warn-primary)}',

      /* 图 */
      '.sp-figure{display:flex;flex-direction:column;gap:6px;min-width:0}',
      '.sp-figure-title{font-size:11px;color:var(--dsw-alias-label-secondary);font-weight:500}',
      '.sp-figure img{display:block;max-width:100%;height:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-2)}',
      '.sp-figure-meta{font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sp-figure-empty{display:flex;align-items:center;justify-content:center;min-height:120px;border:1px dashed var(--dsw-alias-border-l2);border-radius:10px;background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-tertiary);font-size:11.5px;padding:12px;text-align:center}',

      /* 表格 */
      '.sp-table-wrap{max-height:420px;overflow:auto;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.sp-table{width:100%;border-collapse:separate;border-spacing:0;font-size:11.5px}',
      '.sp-table th{position:sticky;top:0;z-index:1;padding:6px 10px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-secondary);font-weight:600;text-align:left;white-space:nowrap}',
      '.sp-table th.sp-th-sort{cursor:pointer;user-select:none}',
      '.sp-table th.sp-th-sort:hover{color:var(--dsw-alias-label-primary)}',
      '.sp-table td{padding:4px 10px;border-bottom:1px solid color-mix(in srgb,var(--dsw-alias-border-l1) 55%,transparent);color:var(--dsw-alias-label-primary);white-space:nowrap;max-width:340px;overflow:hidden;text-overflow:ellipsis}',
      '.sp-table tbody tr:hover td{background:color-mix(in srgb,var(--dsw-alias-brand-primary) 4%,transparent)}',
      '.sp-td-num{text-align:right;font-variant-numeric:tabular-nums}',
      '.sp-table-note{padding:5px 2px 0;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sp-empty-table{padding:18px 12px;text-align:center;color:var(--dsw-alias-label-tertiary);font-size:11.5px}',

      /* 结果视图 */
      '.sp-result{display:flex;flex-direction:column;gap:10px}',
      '.sp-bullets{margin:0;padding-left:16px;display:flex;flex-direction:column;gap:3px;font-size:12px;color:var(--dsw-alias-label-primary)}',
      '.sp-bullets li::marker{color:var(--dsw-alias-brand-primary)}',
      '.sp-result-block{display:flex;flex-direction:column;gap:6px}',
      '.sp-result-block-title{font-size:11px;font-weight:600;color:var(--dsw-alias-label-secondary);letter-spacing:.03em}',
      '.sp-warn-list{display:flex;flex-direction:column;gap:3px;font-size:11.5px;color:var(--dsw-alias-state-warn-primary)}',
      '.sp-note-list{display:flex;flex-direction:column;gap:3px;font-size:11.5px;color:var(--dsw-alias-label-secondary)}',
      '.sp-file-chip{display:inline-flex;align-items:center;gap:6px;max-width:100%;padding:4px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-2);font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sp-file-chip .sp-path{flex:1 1 auto;min-width:0}',

      /* 引导 / 空态 / 崩溃 */
      '.sp-guide{display:flex;flex-direction:column;align-items:flex-start;gap:9px;justify-content:center;padding:22px 20px;border:1.5px dashed var(--dsw-alias-border-l2);border-radius:12px;background:var(--dsw-alias-bg-layer-1);text-align:left}',
      '.sp-guide-title{font-size:13.5px;font-weight:600}',
      '.sp-guide-body{font-size:12px;color:var(--dsw-alias-label-secondary);max-width:640px;line-height:1.55}',
      '.sp-guide-hint{font-size:11.5px;color:var(--dsw-alias-label-tertiary);max-width:660px;line-height:1.5}',
      '.sp-notice{display:flex;align-items:center;gap:8px;padding:9px 12px;border:1px dashed var(--dsw-alias-border-l2);border-radius:10px;background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-secondary);font-size:11.5px}',
      '.sp-crash{display:flex;flex-direction:column;align-items:flex-start;gap:8px;margin:14px;padding:16px;border:1px solid color-mix(in srgb,var(--dsw-alias-state-error-primary) 40%,var(--dsw-alias-border-l1));border-radius:12px;background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 5%,var(--dsw-alias-bg-layer-1))}',
      '.sp-crash-title{font-size:13px;font-weight:600;color:var(--dsw-alias-state-error-primary)}',
      '.sp-crash-hint{font-size:11.5px;color:var(--dsw-alias-label-secondary)}',
      '.sp-crash-detail{font-size:10.5px;color:var(--dsw-alias-label-tertiary);font-family:ui-monospace,SFMono-Regular,Menlo,monospace;word-break:break-all}',
      '.sp-crash-mini{display:inline-flex;align-items:center;justify-content:center;width:16px;height:16px;border-radius:50%;background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 16%,transparent);color:var(--dsw-alias-state-error-primary);font-size:11px;font-weight:700}',
      '.sp-loading{display:inline-flex;align-items:center;gap:7px;font-size:12px;color:var(--dsw-alias-label-secondary)}',
      '.sp-spinner{display:inline-block;width:13px;height:13px;border-radius:50%;border:1.6px solid color-mix(in srgb,var(--dsw-alias-brand-primary) 30%,transparent);border-top-color:var(--dsw-alias-brand-primary);animation:sp-spin .8s linear infinite;flex:none}',
      '@keyframes sp-spin{to{transform:rotate(360deg)}}',
      '@keyframes sp-pulse{0%,100%{box-shadow:0 0 0 0 color-mix(in srgb,var(--dsw-alias-brand-primary) 30%,transparent)}50%{box-shadow:0 0 0 4px color-mix(in srgb,var(--dsw-alias-brand-primary) 6%,transparent)}}',
      '@keyframes sp-fade{from{opacity:0;transform:translateY(2px)}to{opacity:1;transform:none}}',

      /* 弹层 */
      '.sp-modal-mask{position:fixed;inset:0;z-index:60;display:flex;align-items:center;justify-content:center;padding:24px;background:color-mix(in srgb,var(--dsw-alias-label-primary) 34%,transparent);backdrop-filter:blur(3px);animation:sp-fade .14s ease}',
      '.sp-modal{display:flex;flex-direction:column;width:min(920px,94vw);max-height:88vh;border:1px solid var(--dsw-alias-border-l1);border-radius:14px;background:var(--dsw-alias-bg-base);box-shadow:0 24px 80px color-mix(in srgb,var(--dsw-alias-label-primary) 24%,transparent);overflow:hidden}',
      '.sp-modal-wide{width:min(1100px,95vw)}',
      '.sp-modal-head{display:flex;align-items:flex-start;justify-content:space-between;gap:10px;flex:none;padding:13px 16px;border-bottom:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sp-modal-title{font-size:14px;font-weight:600}',
      '.sp-modal-subtitle{font-size:11.5px;color:var(--dsw-alias-label-secondary);margin-top:1px}',
      '.sp-modal-body{display:flex;flex-direction:column;gap:11px;padding:13px 16px 18px;overflow-y:auto;min-height:0}',
      '.sp-modal-foot{display:flex;align-items:center;justify-content:flex-end;gap:8px;flex-wrap:wrap;flex:none;padding:10px 16px;border-top:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-1)}',
      '.sp-advanced{margin-top:1px;padding:9px 12px;border:1px dashed var(--dsw-alias-border-l2);border-radius:10px;background:var(--dsw-alias-bg-layer-1)}',
      '.sp-advanced summary{cursor:pointer;font-size:11.5px;font-weight:600;color:var(--dsw-alias-label-secondary)}',
      '.sp-advanced[open] summary{margin-bottom:7px;color:var(--dsw-alias-label-primary)}',
      '.sp-form-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:8px 10px}',
      '.sp-form-grid .sp-field{flex-direction:column;align-items:stretch;gap:3px}',
      '.sp-form-grid .sp-input{max-width:none;width:100%}',

      /* GEO 检索 */
      '.sp-geo-search{display:flex;align-items:center;gap:7px;flex-wrap:wrap}',
      '.sp-geo-term{flex:1 1 300px;min-width:220px;height:32px;font-size:12.5px}',
      '.sp-geo-filters{display:flex;align-items:center;gap:10px;flex-wrap:wrap;font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sp-geo-layout{display:grid;grid-template-columns:minmax(260px,340px) 1fr;gap:12px;min-height:0;flex:1 1 auto}',
      '.sp-geo-list{display:flex;flex-direction:column;gap:7px;overflow-y:auto;max-height:50vh;padding-right:2px}',
      '.sp-geo-card{display:flex;flex-direction:column;gap:4px;padding:9px 11px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-1);cursor:pointer;text-align:left;font:inherit;transition:border-color .14s ease,background .14s ease,box-shadow .14s ease}',
      '.sp-geo-card:hover{border-color:var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-layer-2)}',
      '.sp-geo-card-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 50%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 6%,var(--dsw-alias-bg-layer-1));box-shadow:0 0 0 1px color-mix(in srgb,var(--dsw-alias-brand-primary) 24%,transparent)}',
      '.sp-geo-gse{font-size:12px;font-weight:600;color:var(--dsw-alias-brand-primary)}',
      '.sp-geo-title{font-size:12px;color:var(--dsw-alias-label-primary);line-height:1.4}',
      '.sp-geo-title-en{font-size:10.5px;color:var(--dsw-alias-label-tertiary);line-height:1.4;display:-webkit-box;-webkit-line-clamp:2;-webkit-box-orient:vertical;overflow:hidden}',
      '.sp-geo-meta{display:flex;flex-wrap:wrap;align-items:center;gap:4px 7px;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',
      '.sp-geo-tag{padding:0 6px;height:17px;display:inline-flex;align-items:center;border-radius:999px;border:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-2);font-size:10px;color:var(--dsw-alias-label-secondary)}',
      '.sp-geo-tag-ok{border-color:color-mix(in srgb,var(--dsw-alias-state-success-primary) 40%,var(--dsw-alias-border-l1));color:var(--dsw-alias-state-success-primary)}',
      '.sp-geo-tag-high{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 44%,var(--dsw-alias-border-l1));color:var(--dsw-alias-brand-primary)}',
      '.sp-geo-tag-mid{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 40%,var(--dsw-alias-border-l1));color:var(--dsw-alias-state-warn-primary)}',
      '.sp-geo-tag-low{color:var(--dsw-alias-label-tertiary)}',
      '.sp-geo-detail{display:flex;flex-direction:column;gap:9px;min-height:0;overflow-y:auto;padding:1px 2px 4px;max-height:54vh}',
      '.sp-geo-detail-title{font-size:13px;font-weight:600;line-height:1.4}',
      '.sp-geo-kv{display:grid;grid-template-columns:84px 1fr;gap:3px 8px;font-size:11.5px}',
      '.sp-geo-kv dt{color:var(--dsw-alias-label-tertiary)}',
      '.sp-geo-kv dd{margin:0;color:var(--dsw-alias-label-primary);word-break:break-word}',
      '.sp-geo-file{display:flex;align-items:center;gap:7px;padding:5px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base);font-size:11.5px}',
      '.sp-geo-file-name{flex:1 1 auto;min-width:0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-progress{display:flex;flex-direction:column;gap:4px;min-width:180px;flex:1 1 auto}',
      '.sp-progress-bar{position:relative;height:7px;border-radius:999px;border:1px solid var(--dsw-alias-border-l1);background:var(--dsw-alias-bg-layer-2);overflow:hidden}',
      '.sp-progress-fill{position:absolute;left:0;top:0;bottom:0;border-radius:999px;background:linear-gradient(90deg,color-mix(in srgb,var(--dsw-alias-brand-primary) 70%,transparent),var(--dsw-alias-brand-primary));transition:width .2s ease}',
      '.sp-progress-text{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:10.5px;color:var(--dsw-alias-label-tertiary)}',

      /* 分析目录 / 运行器 */
      '.sp-analysis-wrap{display:grid;grid-template-columns:minmax(240px,300px) 1fr;gap:12px;align-items:start;min-height:0}',
      '.sp-analysis-list{display:flex;flex-direction:column;gap:8px;max-height:64vh;overflow-y:auto;padding-right:2px}',
      '.sp-analysis-group{display:flex;flex-direction:column;gap:5px}',
      '.sp-analysis-group-title{font-size:10.5px;font-weight:600;letter-spacing:.05em;color:var(--dsw-alias-label-tertiary)}',
      '.sp-analysis-item{display:flex;flex-direction:column;gap:2px;padding:7px 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:9px;background:var(--dsw-alias-bg-layer-1);cursor:pointer;text-align:left;font:inherit;transition:border-color .14s ease,background .14s ease}',
      '.sp-analysis-item:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-layer-2)}',
      '.sp-analysis-item:disabled{opacity:.55;cursor:not-allowed}',
      '.sp-analysis-item-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 48%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 6%,var(--dsw-alias-bg-layer-1))}',
      '.sp-analysis-main{display:flex;align-items:center;justify-content:space-between;gap:8px;min-width:0}',
      '.sp-analysis-name{font-size:12px;font-weight:500;color:var(--dsw-alias-label-primary);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',
      '.sp-analysis-method{font-size:10.5px;color:var(--dsw-alias-label-tertiary);line-height:1.4;word-break:break-word}',
      '.sp-analysis-reason{font-size:10.5px;color:var(--dsw-alias-state-warn-primary);line-height:1.4}',
      '.sp-analysis-runner{display:flex;flex-direction:column;gap:10px;min-width:0}',
      '.sp-analysis-head{display:flex;flex-direction:column;gap:2px;min-width:0}',
      '.sp-analysis-title{font-size:13.5px;font-weight:600;line-height:1.35}',

      /* 反卷积参考面板 */
      '.sp-ref-panel{display:flex;flex-direction:column;gap:8px}',
      '.sp-ref-list{display:flex;flex-direction:column;gap:4px;max-height:180px;overflow-y:auto}',
      '.sp-ref-item{display:flex;align-items:center;justify-content:space-between;gap:8px;padding:5px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-base);font-size:11px}',
      '.sp-ref-item-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 46%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 6%,var(--dsw-alias-bg-base))}',
      '.sp-ref-path{flex:1 1 auto;min-width:0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}',

      /* 空间图控制条 */
      '.sp-color-row{display:flex;align-items:center;gap:4px;flex-wrap:wrap}',
      '.sp-color-opt{display:inline-flex;align-items:center;gap:5px;height:24px;padding:0 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:6px;background:var(--dsw-alias-bg-layer-2);color:var(--dsw-alias-label-secondary);font:inherit;font-size:11px;cursor:pointer;transition:border-color .14s ease,background .14s ease,color .14s ease}',
      '.sp-color-opt:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);color:var(--dsw-alias-label-primary)}',
      '.sp-color-opt-on{border-color:color-mix(in srgb,var(--dsw-alias-brand-primary) 46%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-brand-primary) 10%,var(--dsw-alias-bg-layer-1));color:var(--dsw-alias-brand-primary)}',
      '.sp-color-dot{width:7px;height:7px;border-radius:50%;flex:none;background:currentColor}',
      '.sp-gene-row{display:flex;align-items:center;gap:6px;flex-wrap:wrap}',
      '.sp-minimap{display:flex;align-items:center;justify-content:center;padding:6px;border:1px solid var(--dsw-alias-border-l1);border-radius:10px;background:var(--dsw-alias-bg-layer-2)}',
      '.sp-minimap svg{display:block;width:100%;height:auto;max-height:220px}',

      /* 导出按钮组 */
      '.sp-export-grid{display:flex;flex-wrap:wrap;gap:6px}',
      '.sp-export-item{display:flex;flex-direction:column;gap:2px;align-items:flex-start;min-width:118px;padding:7px 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:9px;background:var(--dsw-alias-bg-layer-1);color:var(--dsw-alias-label-primary);font:inherit;font-size:11.5px;cursor:pointer;text-align:left;transition:border-color .14s ease,background .14s ease}',
      '.sp-export-item:hover:not(:disabled){border-color:var(--dsw-alias-border-l2);background:var(--dsw-alias-bg-layer-2)}',
      '.sp-export-item:disabled{opacity:.5;cursor:default}',
      '.sp-export-hint{font-size:10px;color:var(--dsw-alias-label-tertiary)}',

      /* 右栏 Agent 面板 */
      '.sp-agent-section{display:flex;flex-direction:column;gap:6px;padding:10px 12px;border-bottom:1px solid var(--dsw-alias-border-l1)}',
      '.sp-agent-title{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:10.5px;font-weight:600;letter-spacing:.05em;color:var(--dsw-alias-label-tertiary)}',
      '.sp-suggestion{display:flex;flex-direction:column;gap:4px;padding:8px 9px;border:1px solid var(--dsw-alias-border-l1);border-radius:9px;background:var(--dsw-alias-bg-layer-2);animation:sp-fade .16s ease}',
      '.sp-suggestion-title{font-size:11.5px;font-weight:600;color:var(--dsw-alias-label-primary);line-height:1.4}',
      '.sp-suggestion-detail{font-size:10.5px;color:var(--dsw-alias-label-tertiary);line-height:1.45}',
      '.sp-activity{display:flex;flex-direction:column;gap:6px}',
      '.sp-activity-item{display:flex;gap:7px;font-size:11px;line-height:1.4}',
      '.sp-activity-dot{position:relative;flex:none;width:6px;height:6px;margin-top:4px;border-radius:50%;background:var(--dsw-alias-state-success-primary);box-shadow:0 0 0 2px color-mix(in srgb,var(--dsw-alias-state-success-primary) 16%,transparent)}',
      '.sp-activity-dot-bad{background:var(--dsw-alias-state-error-primary);box-shadow:0 0 0 2px color-mix(in srgb,var(--dsw-alias-state-error-primary) 16%,transparent)}',
      '.sp-activity-main{display:flex;flex-direction:column;gap:1px;min-width:0}',
      '.sp-activity-head{display:flex;align-items:center;gap:5px;flex-wrap:wrap;color:var(--dsw-alias-label-secondary)}',
      '.sp-activity-actor{font-size:9.5px;padding:0 5px;height:15px;display:inline-flex;align-items:center;border-radius:999px;border:1px solid var(--dsw-alias-border-l1);color:var(--dsw-alias-label-tertiary)}',
      '.sp-activity-text{color:var(--dsw-alias-label-tertiary);font-size:10.5px;word-break:break-word}',
      '.sp-stat-row{display:flex;align-items:center;justify-content:space-between;gap:8px;font-size:11px;color:var(--dsw-alias-label-secondary)}',
      '.sp-stat-row b{font-weight:600;color:var(--dsw-alias-label-primary);font-variant-numeric:tabular-nums;max-width:60%;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}',
      '.sp-history{display:flex;flex-direction:column;gap:4px}',
      '.sp-history-item{display:flex;align-items:center;gap:7px;padding:5px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-1);font-size:11px}',
      '.sp-history-at{flex:none;font-size:10.5px;color:var(--dsw-alias-label-tertiary);font-variant-numeric:tabular-nums}',
      '.sp-history-step{flex:none;font-weight:600;color:var(--dsw-alias-label-secondary)}',
      '.sp-history-summary{flex:1 1 auto;min-width:0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;color:var(--dsw-alias-label-tertiary)}',
      '.sp-history-ms{flex:none;color:var(--dsw-alias-label-tertiary);font-size:10.5px}',
      '.sp-audit-item{display:flex;align-items:flex-start;gap:7px;padding:6px 8px;border:1px solid var(--dsw-alias-border-l1);border-radius:8px;background:var(--dsw-alias-bg-layer-1);font-size:11px}',
      '.sp-audit-item-warn{border-color:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 36%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 5%,var(--dsw-alias-bg-layer-1))}',
      '.sp-audit-item-fail{border-color:color-mix(in srgb,var(--dsw-alias-state-error-primary) 36%,var(--dsw-alias-border-l1));background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 5%,var(--dsw-alias-bg-layer-1))}',
      '.sp-audit-mark{flex:none;width:16px;height:16px;display:inline-flex;align-items:center;justify-content:center;border-radius:50%;font-size:9.5px;font-weight:700;background:color-mix(in srgb,var(--dsw-alias-state-success-primary) 16%,transparent);color:var(--dsw-alias-state-success-primary)}',
      '.sp-audit-mark-warn{background:color-mix(in srgb,var(--dsw-alias-state-warn-primary) 16%,transparent);color:var(--dsw-alias-state-warn-primary)}',
      '.sp-audit-mark-fail{background:color-mix(in srgb,var(--dsw-alias-state-error-primary) 16%,transparent);color:var(--dsw-alias-state-error-primary)}',
      '.sp-audit-main{display:flex;flex-direction:column;gap:1px;min-width:0}',
      '.sp-audit-title{font-weight:600;color:var(--dsw-alias-label-primary)}',
      '.sp-audit-detail{color:var(--dsw-alias-label-tertiary);font-size:10.5px;line-height:1.45}',
      '.sp-score{display:flex;align-items:baseline;gap:6px}',
      '.sp-score-v{font-size:30px;font-weight:250;line-height:1;color:var(--dsw-alias-label-primary)}',
      '.sp-score-label{font-size:11px;color:var(--dsw-alias-label-tertiary)}',

      /* 响应式 */
      '@media (max-width:1180px){.sp-left{flex-basis:240px;width:240px}.sp-right{flex-basis:262px;width:262px}}',
      '@media (max-width:980px){',
      '  .sp-left{flex-basis:216px;width:216px}.sp-right{display:none}',
      '  .sp-geo-layout{grid-template-columns:1fr}',
      '  .sp-analysis-wrap{grid-template-columns:1fr}',
      '  .sp-field .sp-input{max-width:120px}',
      '}',
      '@media (prefers-reduced-motion:reduce){.sp-spinner{animation:none}.sp-step-dot-running{animation:none}}',
    ].join('')

    // ==========================================================================
    // 通用小组件：加载 / 错误 / 图标 / 按钮 / 错误边界 / 侧边栏图标
    // ==========================================================================

    function Spinner() {
      return h('span', { className: 'sp-spinner', 'aria-hidden': 'true' })
    }

    function LoadingLine(props) {
      return h('div', { className: 'sp-loading' }, h(Spinner), h('span', null, (props && props.text) || T('loading')))
    }

    function ErrorBanner(props) {
      return h('div', { className: 'sp-crash', role: 'alert' },
        h('div', { className: 'sp-crash-title' }, (props && props.title) || T('errorTitle')),
        props && props.message ? h('div', { className: 'sp-crash-detail' }, String(props.message)) : null,
        props && props.onRetry
          ? h('button', { type: 'button', className: 'sp-btn sp-btn-sm', onClick: props.onRetry },
              (props && props.retryLabel) || T('retry'))
          : null)
    }

    function InlineError(props) {
      if (!props || !props.message) return null
      return h('div', { className: 'sp-inline-err', role: 'alert' }, String(props.message))
    }

    function Notice(props) {
      return h('div', { className: 'sp-notice' }, (props && props.text) || '')
    }

    /** 内联 SVG 图标（不用 emoji；stroke 用 currentColor 跟随主题）。 */
    function Icon(props) {
      var name = (props && props.name) || 'dot'
      var size = (props && props.size) || 14
      var paths = {
        refresh: 'M20 11a8 8 0 1 0-2.3 5.7M20 4v7h-7',
        search: 'M10.5 17a6.5 6.5 0 1 0 0-13 6.5 6.5 0 0 0 0 13zM15.2 15.2 21 21',
        importIcon: 'M12 3v12m0 0 4-4m-4 4-4-4M4 17v3h16v-3',
        download: 'M12 3v12m0 0 4-4m-4 4-4-4M4 20h16',
        close: 'M6 6l12 12M18 6 6 18',
        play: 'M7 5l12 7-12 7z',
        chart: 'M4 20V10M10 20V4M16 20v-7M22 20H2',
        chevron: 'M9 6l6 6-6 6',
        chevronDown: 'M6 9l6 6 6-6',
        copy: 'M9 9h11v11H9zM5 15V4h11',
        check: 'M4 12.5l5 5L20 6.5',
        folder: 'M3 6h6l2 2h10v11H3z',
        plus: 'M12 5v14M5 12h14',
        flask: 'M9 3h6M10 3v6l-5 9a2 2 0 0 0 1.7 3h10.6a2 2 0 0 0 1.7-3l-5-9V3',
        shield: 'M12 3l8 3v6c0 5-3.5 8-8 9-4.5-1-8-4-8-9V6z',
        clock: 'M12 21a9 9 0 1 0 0-18 9 9 0 0 0 0 18zM12 7v5l3.2 2',
        database: 'M12 3c4.4 0 8 1.3 8 3s-3.6 3-8 3-8-1.3-8-3 3.6-3 8-3zM4 6v6c0 1.7 3.6 3 8 3s8-1.3 8-3V6M4 12v6c0 1.7 3.6 3 8 3s8-1.3 8-3v-6',
        grid: 'M4 4h7v7H4zM13 4h7v7h-7zM4 13h7v7H4zM13 13h7v7h-7z',
        dot: 'M12 12h.01',
      }
      return h('svg', {
        width: size,
        height: size,
        viewBox: '0 0 24 24',
        fill: name === 'play' ? 'currentColor' : 'none',
        stroke: 'currentColor',
        strokeWidth: props && props.strokeWidth ? props.strokeWidth : 1.7,
        strokeLinecap: 'round',
        strokeLinejoin: 'round',
        'aria-hidden': 'true',
        focusable: 'false',
        style: { flex: 'none', display: 'inline-block' },
      }, h('path', { d: paths[name] || paths.dot }))
    }

    function Btn(props) {
      var children = Array.prototype.slice.call(arguments, 1)
      var extra = props || {}
      // 两种调用方式都支持：Btn(props, ...children) 直呼；h(Btn, props, ...children) 经 React。
      if (!children.length && extra.children !== undefined && extra.children !== null && extra.children !== false) {
        children = Array.isArray(extra.children) ? extra.children : [extra.children]
      }
      var className = 'sp-btn' +
        (extra.variant === 'primary' ? ' sp-btn-primary' : '') +
        (extra.variant === 'ghost' ? ' sp-btn-ghost' : '') +
        (extra.variant === 'danger' ? ' sp-btn-danger' : '') +
        (extra.small ? ' sp-btn-sm' : '') +
        (extra.className ? ' ' + extra.className : '')
      return h('button', {
        key: extra.key,
        type: 'button',
        className: className,
        title: extra.title || undefined,
        disabled: extra.disabled === true,
        onClick: extra.onClick,
      }, extra.icon ? h(Icon, { name: extra.icon, size: extra.small ? 12 : 14 }) : null, children.filter(Boolean))
    }

    /** 渲染崩溃兜底：每个主要区块都包一层，避免整页白屏。 */
    class PluginErrorBoundary extends React.Component {
      constructor(props) {
        super(props)
        this.state = { error: null }
        this.handleReset = this.handleReset.bind(this)
      }

      static getDerivedStateFromError(error) {
        return { error: error || new Error('unknown render error') }
      }

      componentDidCatch(error, info) {
        warnPlugin('UI crashed: ' + formatError(error), info)
      }

      handleReset() {
        this.setState({ error: null })
      }

      render() {
        if (this.state.error) {
          if (this.props && this.props.mini) {
            return h('span', { className: 'sp-crash-mini', title: T('crashTitle') }, '!')
          }
          return h('div', { className: 'sp-crash', role: 'alert' },
            h('div', { className: 'sp-crash-title' }, T('crashTitle')),
            h('div', { className: 'sp-crash-hint' }, T('crashHint')),
            h('div', { className: 'sp-crash-detail' }, formatError(this.state.error)),
            h('button', { type: 'button', className: 'sp-btn sp-btn-primary sp-btn-sm', onClick: this.handleReset }, T('retry')))
        }
        return this.props ? this.props.children : null
      }
    }

    // ==========================================================================
    // 侧边栏图标（纯 SVG：组织切片轮廓 + 网格点邻域）
    // ==========================================================================

    function SidebarGlyph(props) {
      var size = (props && props.size) || 18
      return h('svg', {
        width: size,
        height: size,
        viewBox: '0 0 24 24',
        fill: 'none',
        stroke: 'currentColor',
        strokeWidth: 1.6,
        strokeLinecap: 'round',
        strokeLinejoin: 'round',
        'aria-hidden': 'true',
        focusable: 'false',
      },
        h('path', { d: 'M12 3.4c3.2 0 5.6 1.4 5.6 3.2v5.3c0 1.3 2.1 2 2.1 4.2 0 2.6-3.4 4.5-7.7 4.5S4.3 18.7 4.3 16.1c0-2.2 2.1-2.9 2.1-4.2V6.6c0-1.8 2.4-3.2 5.6-3.2z' }),
        h('path', { d: 'M8.6 9.2h6.8M8.6 12.6h6.8M8.6 16h6.8' }),
        h('circle', { cx: 9.4, cy: 9.2, r: 1.05, fill: 'currentColor', stroke: 'none' }),
        h('circle', { cx: 14.6, cy: 12.6, r: 1.05, fill: 'currentColor', stroke: 'none' }),
        h('circle', { cx: 9.4, cy: 16, r: 1.05, fill: 'currentColor', stroke: 'none' }))
    }

    function SidebarIcon(props) {
      return h(PluginErrorBoundary, { mini: true }, h(SidebarGlyph, { size: props && props.size }))
    }

    // ==========================================================================
    // 通用展示组件：指标卡 / 卡片 / 图 / 表 / 结果信封
    // ==========================================================================

    function MetricCard(props) {
      return h('div', { className: 'sp-metric', title: props.hint || undefined },
        h('div', { className: 'sp-metric-k' }, props.label),
        h('div', { className: 'sp-metric-v' }, props.value, props.unit ? h('small', null, props.unit) : null))
    }

    function SectionCard(props) {
      return h('div', { className: 'sp-card' },
        props.title || props.extra
          ? h('div', { className: 'sp-card-title' },
              h('span', null, props.title || ''),
              props.extra ? h('span', null, props.extra) : null)
          : null,
        props.children)
    }

    /** 单张图：加载 / 失败 / 空态 / 图片（点击在新标签打开原图）。 */
    function FigureBlock(props) {
      var brokenState = useState(false)
      var broken = brokenState[0]
      var setBroken = brokenState[1]
      var figure = props.figure || null
      var url = figure && figure.url ? figure.url : ''
      useEffect(function () {
        setBroken(false)
      }, [url])
      if (props.loading) {
        return h('div', { className: 'sp-figure' },
          props.title ? h('div', { className: 'sp-figure-title' }, props.title) : null,
          h('div', { className: 'sp-figure-empty' }, h(LoadingLine, { text: props.loadingText || T('plotting') })))
      }
      if (props.error) {
        return h('div', { className: 'sp-figure' },
          props.title ? h('div', { className: 'sp-figure-title' }, props.title) : null,
          h('div', { className: 'sp-inline-err' }, props.error),
          props.onReload ? h('div', null, Btn({ small: true, icon: 'refresh', onClick: props.onReload }, T('retry'))) : null)
      }
      if (!url) {
        return h('div', { className: 'sp-figure' },
          props.title ? h('div', { className: 'sp-figure-title' }, props.title) : null,
          h('div', { className: 'sp-figure-empty' }, props.emptyText || T('noFigure')))
      }
      return h('div', { className: 'sp-figure' },
        props.title ? h('div', { className: 'sp-figure-title' }, props.title) : null,
        broken
          ? h('div', { className: 'sp-figure-empty' }, T('figureFail'))
          : h('a', {
              href: url,
              target: '_blank',
              rel: 'noopener noreferrer',
              title: T('plot'),
              style: { display: 'block' },
            }, h('img', {
              src: url,
              alt: figure.title || props.title || 'figure',
              loading: 'lazy',
              onError: function () { setBroken(true) },
            })),
        figure.file && figure.file !== url ? h('div', { className: 'sp-figure-meta', title: figure.file }, figure.file) : null,
        props.onReload && !broken ? h('div', null, Btn({ small: true, variant: 'ghost', icon: 'refresh', onClick: props.onReload }, T('retry'))) : null)
    }

    /** 自动加载 + 缓存展示一张服务端图（固定参数用；动态参数由各 tab 自己管）。 */
    function PlotBlock(props) {
      var key = plotKey(props.kind, props.params)
      var entry = (props.plots && props.plots[key]) || null
      var loadPlot = props.loadPlot
      var enabled = props.enabled !== false
      useEffect(function () {
        if (!enabled || typeof loadPlot !== 'function') return
        if (entry && (entry.loading || entry.figure || entry.error)) return
        loadPlot(props.kind, props.params)
      }, [enabled, key, entry ? (entry.loading ? 2 : (entry.figure ? 1 : 3)) : 0])
      var reload = function () { if (typeof loadPlot === 'function') loadPlot(props.kind, props.params, true) }
      return h(FigureBlock, {
        title: props.title,
        figure: entry && entry.figure ? entry.figure : null,
        loading: !!(entry && entry.loading),
        error: entry && entry.error ? entry.error : '',
        emptyText: props.emptyText || T('noFigure'),
        onReload: enabled ? reload : null,
      })
    }

    /** 服务端表格：列头可排序，最多渲染 MAX_TABLE_RENDER_ROWS 行。 */
    function DataTable(props) {
      var table = props.table || { columns: [], rows: [] }
      var sortState = useState({ key: '', dir: 'asc' })
      var sort = sortState[0]
      var setSort = sortState[1]
      var limit = isFinite(Number(props.maxRows)) ? Number(props.maxRows) : MAX_TABLE_RENDER_ROWS
      var rows = Array.isArray(table.rows) ? table.rows : []
      var view
      if (props.sortable !== false) {
        var ordered = sort.key ? sortRows(rows, sort.key, sort.dir, table.columns) : rows
        view = formatTableRows({ columns: table.columns, rows: ordered, truncated: table.truncated }, limit)
      } else {
        view = formatTableRows(table, limit)
      }
      function clickHeader(column) {
        if (props.sortable === false) return
        setSort(function (previous) {
          if (previous.key === column.key) return { key: column.key, dir: previous.dir === 'asc' ? 'desc' : 'asc' }
          return { key: column.key, dir: 'asc' }
        })
      }
      if (!view.columns.length) {
        return h('div', { className: 'sp-empty-table' }, props.emptyText || T('emptyTable'))
      }
      return h('div', null,
        h('div', { className: 'sp-table-wrap', style: props.maxHeight ? { maxHeight: props.maxHeight } : undefined },
          h('table', { className: 'sp-table' },
            h('thead', null,
              h('tr', null, view.columns.map(function (column) {
                var active = sort.key === column.key
                return h('th', {
                  key: column.key,
                  className: (props.sortable === false ? '' : 'sp-th-sort') + (column.type === 'number' ? ' sp-td-num' : ''),
                  onClick: function () { clickHeader(column) },
                }, column.label + (active ? (sort.dir === 'asc' ? ' ↑' : ' ↓') : ''))
              }))),
            h('tbody', null, view.rows.map(function (row) {
              return h('tr', { key: row.key }, row.cells.map(function (cell) {
                return h('td', { key: cell.key, className: cell.numeric ? 'sp-td-num' : '', title: cell.text }, cell.text)
              }))
            })))),
        h('div', { className: 'sp-table-note' },
          tf(T, 'tableRows', view.total) +
          (view.truncated ? ' · ' + tf(T, 'tableTruncated', view.shown) : '') +
          (table.note ? ' · ' + table.note : '')))
    }

    /**
     * 通用结果信封渲染（DESIGN §2.5）：bullets / metrics / tables / figures /
     * warnings / notes / data_out。标题优先用服务端文案，缺省回落字典词条。
     */
    function ResultView(props) {
      var result = props.result || null
      if (!result) {
        return props.emptyText ? h('div', { className: 'sp-notice' }, props.emptyText) : null
      }
      var blocks = []
      if (typeof props.onExport === 'function') {
        blocks.push(h('div', { key: 'export', className: 'sp-toolbar' },
          Btn({ small: true, icon: 'download', disabled: props.exportBusy === true, onClick: props.onExport },
            props.exportBusy ? T('recordsExporting') : (props.exportLabel || T('exportBtn'))),
          props.exportFile ? h('span', { className: 'sp-file-chip' }, h('span', { className: 'sp-path' }, props.exportFile)) : null))
      }
      var bullets = Array.isArray(result.bullets) ? result.bullets : []
      if (bullets.length) {
        blocks.push(h('div', { key: 'bullets', className: 'sp-card' },
          h('div', { className: 'sp-card-title' }, h('span', null, result.title || T('analysisResult'))),
          h('ul', { className: 'sp-bullets' }, bullets.map(function (bullet, index) {
            return h('li', { key: 'b-' + index }, bullet)
          }))))
      }
      var metrics = result.metrics && typeof result.metrics === 'object' ? result.metrics : {}
      var metricKeys = Object.keys(metrics)
      if (metricKeys.length) {
        blocks.push(h('div', { key: 'metrics', className: 'sp-metrics' }, metricKeys.map(function (key) {
          var value = metrics[key]
          return h(MetricCard, {
            key: key,
            label: key,
            value: typeof value === 'number' ? fmtNumber(value) : formatCell(value, false),
          })
        })))
      }
      if (result.tables && result.tables.length) {
        blocks.push(h('div', { key: 'tables', className: 'sp-card' },
          h('div', { className: 'sp-card-title' },
            h('span', null, T('resultColumn')),
            h('span', null, String(result.tables.length))),
          result.tables.map(function (table, index) {
            return h('div', { key: 'table-' + (table.id || index), className: 'sp-result-block' },
              h('div', { className: 'sp-figure-title' }, table.title || tf(T, 'tableRows', (table.rows || []).length)),
              h(DataTable, { table: table, maxRows: props.maxTableRows || MAX_TABLE_RENDER_ROWS }))
          })))
      }
      if (result.figures && result.figures.length) {
        blocks.push(h('div', { key: 'figures', className: 'sp-grid-2' }, result.figures.map(function (figure, index) {
          return h(FigureBlock, {
            key: 'figure-' + (figure.id || index),
            title: figure.title || T('plot'),
            figure: figure,
          })
        })))
      }
      if (result.warnings && result.warnings.length) {
        blocks.push(h('div', { key: 'warnings', className: 'sp-card' },
          h('div', { className: 'sp-card-title' }, h('span', null, T('warnings')), h('span', null, String(result.warnings.length))),
          h('div', { className: 'sp-warn-list' }, result.warnings.map(function (warning, index) {
            return h('div', { key: 'w-' + index }, '· ' + warning)
          }))))
      }
      if (result.notes && result.notes.length) {
        blocks.push(h('div', { key: 'notes', className: 'sp-card' },
          h('div', { className: 'sp-card-title' }, h('span', null, T('notes')), h('span', null, String(result.notes.length))),
          h('div', { className: 'sp-note-list' }, result.notes.map(function (note, index) {
            return h('div', { key: 'n-' + index }, '· ' + note)
          }))))
      }
      if (result.data_out && result.data_out.length) {
        blocks.push(h('div', { key: 'files', className: 'sp-card' },
          h('div', { className: 'sp-card-title' }, h('span', null, T('pathLabel'))),
          h('div', { className: 'sp-inline-row' }, result.data_out.map(function (item, index) {
            return h('span', { key: 'f-' + index, className: 'sp-file-chip', title: item.path },
              h(Icon, { name: 'download', size: 11 }),
              h('span', null, item.label || item.kind),
              h('span', { className: 'sp-path' }, item.path))
          }))))
      }
      if (!blocks.length) {
        return h('div', { className: 'sp-notice' }, T('stepSummaryEmpty'))
      }
      return h('div', { className: 'sp-result' }, blocks)
    }

    // ==========================================================================
    // 参数表单：按 schema 通用渲染
    //   number / boolean / select / column / multiselect / genes / columns /
    //   text / string / file / json / object
    // ==========================================================================

    function ParamControl(props) {
      var descriptor = props.descriptor || {}
      var value = props.value
      var disabled = props.disabled === true
      var onChange = props.onChange || function () {}
      var type = descriptor.type || 'text'
      if (type === 'number') {
        return h('input', {
          className: 'sp-input',
          type: 'number',
          min: isNil(descriptor.min) ? undefined : descriptor.min,
          max: isNil(descriptor.max) ? undefined : descriptor.max,
          step: isNil(descriptor.step) ? 'any' : descriptor.step,
          value: isNil(value) ? '' : String(value),
          disabled: disabled,
          onChange: function (event) { onChange(event.target.value === '' ? null : Number(event.target.value)) },
        })
      }
      if (type === 'boolean') {
        return h('label', { className: 'sp-check' },
          h('input', {
            type: 'checkbox',
            checked: value === true,
            disabled: disabled,
            onChange: function (event) { onChange(event.target.checked) },
          }),
          h('span', null, value === true ? T('yes') : T('no')))
      }
      if (type === 'select') {
        var options = Array.isArray(descriptor.options) ? descriptor.options : []
        var current = isNil(value) ? '' : String(value)
        var known = false
        for (var optionIndex = 0; optionIndex < options.length; optionIndex++) {
          if (String(options[optionIndex].value) === current) known = true
        }
        var list = options.slice()
        if (current && !known) list.unshift({ value: current, label: current })
        return h('select', {
          className: 'sp-input sp-select',
          value: current,
          disabled: disabled,
          onChange: function (event) { onChange(event.target.value) },
        }, list.map(function (option) {
          return h('option', { key: String(option.value), value: String(option.value) }, option.label)
        }))
      }
      if (type === 'column') {
        var columns = Array.isArray(props.columns) ? props.columns : []
        var columnOptions = [{ value: '', label: T('noSelection') }]
        for (var columnIndex = 0; columnIndex < columns.length; columnIndex++) {
          var column = columns[columnIndex]
          columnOptions.push({
            value: column.key,
            label: column.key + (column.kind && column.kind !== 'categorical' ? ' · ' + column.kind : ''),
          })
        }
        var columnValue = isNil(value) ? '' : String(value)
        var hasColumn = false
        for (var checkIndex = 0; checkIndex < columnOptions.length; checkIndex++) {
          if (columnOptions[checkIndex].value === columnValue) hasColumn = true
        }
        if (columnValue && !hasColumn) columnOptions.push({ value: columnValue, label: columnValue })
        return h('select', {
          className: 'sp-input sp-select',
          value: columnValue,
          disabled: disabled,
          onChange: function (event) { onChange(event.target.value || null) },
        }, columnOptions.map(function (option) {
          return h('option', { key: option.value || 'none', value: option.value }, option.label)
        }))
      }
      if (type === 'multiselect') {
        var optionList = Array.isArray(descriptor.options) ? descriptor.options : []
        var selected = Array.isArray(value) ? value : parseTokenList(value)
        return h('div', { className: 'sp-multi' }, optionList.map(function (option) {
          var on = selected.indexOf(option.value) >= 0
          return h('button', {
            key: String(option.value),
            type: 'button',
            disabled: disabled,
            className: 'sp-multi-opt' + (on ? ' sp-multi-opt-on' : ''),
            onClick: function () {
              var next = selected.slice()
              var index = next.indexOf(option.value)
              if (index >= 0) next.splice(index, 1)
              else next.push(option.value)
              onChange(next)
            },
          }, (on ? '✓ ' : '') + (option.label || option.value))
        }))
      }
      if (type === 'json' || type === 'object') {
        var jsonText = ''
        if (value && typeof value === 'object') {
          try { jsonText = JSON.stringify(value) } catch (ignored) { jsonText = '' }
        } else if (!isNil(value)) {
          jsonText = String(value)
        }
        return h('textarea', {
          className: 'sp-input sp-textarea',
          defaultValue: jsonText,
          disabled: disabled,
          spellCheck: false,
          placeholder: descriptor.placeholder || '{"1":"region A"}',
          onChange: function (event) { onChange(parseJsonParam(event.target.value)) },
        })
      }
      var textValue = isNil(value) ? '' : (Array.isArray(value) ? value.join(', ') : String(value))
      var isList = type === 'genes' || type === 'columns' || type === 'array'
      if (isList) {
        return h('textarea', {
          className: 'sp-input sp-textarea',
          value: textValue,
          disabled: disabled,
          spellCheck: false,
          placeholder: descriptor.placeholder || (type === 'genes' ? 'CDH1, VIM, MKI67' : ''),
          onChange: function (event) { onChange(parseTokenList(event.target.value)) },
        })
      }
      return h('input', {
        className: 'sp-input',
        type: 'text',
        value: textValue,
        disabled: disabled,
        spellCheck: false,
        placeholder: descriptor.placeholder || (type === 'file' ? '/absolute/path' : ''),
        onChange: function (event) { onChange(event.target.value) },
      })
    }

    function ParamForm(props) {
      var params = Array.isArray(props.params) ? props.params : []
      var values = props.values && typeof props.values === 'object' ? props.values : {}
      var onChange = props.onChange || function () {}
      var columns = Array.isArray(props.columns) ? props.columns : []
      if (!params.length) return h('div', { className: 'sp-muted' }, T('analysisNoParams'))
      return h('div', { className: 'sp-form' }, params.map(function (descriptor) {
        if (!descriptor || !descriptor.key) return null
        var label = descriptor.label || descriptor.key
        var control = h(ParamControl, {
          key: 'control-' + descriptor.key,
          descriptor: descriptor,
          columns: columns,
          value: values[descriptor.key],
          disabled: props.disabled,
          onChange: function (next) {
            var patch = {}
            patch[descriptor.key] = next
            onChange(patch)
          },
        })
        var full = descriptor.type === 'multiselect' || descriptor.type === 'genes' ||
          descriptor.type === 'columns' || descriptor.type === 'array' ||
          descriptor.type === 'json' || descriptor.type === 'object' || descriptor.type === 'boolean'
        return h('div', { key: descriptor.key, className: 'sp-field' + (full ? ' sp-field-col' : '') },
          h('span', { className: 'sp-field-label', title: label }, label),
          control)
      }))
    }

    function Modal(props) {
      function onMask(event) {
        if (event.target === event.currentTarget && !props.busy && props.onClose) props.onClose()
      }
      return h('div', { className: 'sp-modal-mask', onMouseDown: onMask },
        h('div', { className: 'sp-modal' + (props.wide ? ' sp-modal-wide' : '') },
          h('div', { className: 'sp-modal-head' },
            h('div', null,
              h('div', { className: 'sp-modal-title' }, props.title || ''),
              props.subtitle ? h('div', { className: 'sp-modal-subtitle' }, props.subtitle) : null),
            h('button', {
              type: 'button',
              className: 'sp-btn sp-btn-sm sp-btn-ghost',
              disabled: props.busy === true,
              title: T('close'),
              onClick: props.onClose,
            }, h(Icon, { name: 'close', size: 13 }))),
          props.children,
          props.footer ? h('div', { className: 'sp-modal-foot' }, props.footer) : null))
    }

    function GuideBlock(props) {
      return h('div', { className: 'sp-guide' },
        h('div', { className: 'sp-guide-title' }, props.title || ''),
        props.body ? h('div', { className: 'sp-guide-body' }, props.body) : null,
        props.hint ? h('div', { className: 'sp-guide-hint' }, props.hint) : null,
        props.actions ? h('div', { className: 'sp-inline-row' }, props.actions) : null)
    }

    /** 弹层里的字段行。 */
    function FormField(props) {
      return h('label', { className: 'sp-field' + (props.full ? ' sp-field-col' : '') },
        h('span', { className: 'sp-field-label' }, props.label),
        props.children)
    }

    // ==========================================================================
    // 顶栏：品牌 + 引擎徽章 + 数据集 chips + 四个动作
    // ==========================================================================

    function TopBar(props) {
      var badge = engineBadge(props.health)
      var datasets = Array.isArray(props.datasets) ? props.datasets : []
      var tone = props.health && props.health.ready ? 'ok' : (props.health && props.health.error ? 'bad' : 'warn')
      return h('header', { className: 'sp-top' },
        h('div', { className: 'sp-brand' },
          h('span', { className: 'sp-brand-mark' }, h(SidebarGlyph, { size: 17 })),
          h('div', { className: 'sp-title-block' },
            h('div', { className: 'sp-title', title: T('appTitle') }, T('appTitle')),
            h('div', { className: 'sp-subtitle', title: T('appSubtitle') }, T('appSubtitle')))),
        h('span', { className: 'sp-badge sp-badge-' + tone, title: badge.detail || badge.reason || '' },
          h('span', { className: 'sp-badge-dot' }),
          T(badge.labelKey) + (badge.detail ? ' · ' + badge.detail : '')),
        !badge.ready && badge.reason
          ? h('span', { className: 'sp-badge-reason', title: badge.reason }, badge.reason)
          : null,
        h('div', { className: 'sp-chips' }, datasets.map(function (dataset) {
          var on = dataset.id === props.activeId
          var title = (dataset.path || dataset.name || dataset.id) +
            (dataset.n_spots ? ' · ' + fmtNumber(dataset.n_spots) + ' ' + T('spots') : '')
          return h('span', { key: dataset.id, className: 'sp-chip' + (on ? ' sp-chip-on' : ''), title: title },
            h('button', {
              type: 'button',
              className: 'sp-chip-label',
              disabled: props.busy === true,
              title: on ? T('current') : T('selectDataset'),
              onClick: function () { if (!on && props.onSelect) props.onSelect(dataset.id) },
            }, dataset.name || dataset.id),
            on && props.onCloseDataset
              ? h('button', {
                  type: 'button',
                  className: 'sp-chip-x',
                  title: T('closeDataset'),
                  disabled: props.busy === true,
                  onClick: function (event) {
                    event.stopPropagation()
                    props.onCloseDataset(dataset.id)
                  },
                }, '×')
              : null)
        })),
        h('div', { className: 'sp-top-actions' },
          Btn({ small: true, icon: 'search', onClick: props.onGeo, disabled: props.busy === true }, T('searchGeo')),
          Btn({ small: true, icon: 'importIcon', onClick: props.onImport, disabled: props.busy === true }, T('importData')),
          Btn({ small: true, icon: 'flask', onClick: props.onDemo, disabled: props.busy === true, title: T('demoLoad') }, T('demo')),
          Btn({ small: true, icon: 'refresh', onClick: props.onRefresh, disabled: props.refreshing === true, title: T('refreshHint') }, T('refresh'))))
    }

    // ==========================================================================
    // 左栏：数据集摘要卡 + 12 张流程卡
    // ==========================================================================

    function columnSummary(dataset, keys) {
      var columns = dataset && Array.isArray(dataset.sample_columns) ? dataset.sample_columns : []
      for (var index = 0; index < columns.length; index++) {
        for (var keyIndex = 0; keyIndex < keys.length; keyIndex++) {
          if (columns[index].key === keys[keyIndex]) return columns[index]
        }
      }
      return null
    }

    function DatasetCard(props) {
      var dataset = props.dataset
      if (!dataset) {
        return h('div', { className: 'sp-dataset-card' },
          h('div', { className: 'sp-dataset-name' }, T('datasetNoDataset')),
          h('div', { className: 'sp-dataset-meta' }, T('datasetOpenHint')))
      }
      var geometry = dataset.geometry || {}
      var sliceColumn = columnSummary(dataset, ['slice', 'slices'])
      var clusterKeys = Array.isArray(dataset.cluster_keys) ? dataset.cluster_keys : []
      var meta = []
      if (!isNil(dataset.n_spots)) meta.push(h('span', { key: 's' }, T('spots') + ' ', h('b', null, fmtNumber(dataset.n_spots))))
      if (!isNil(dataset.n_genes)) meta.push(h('span', { key: 'g' }, T('genes') + ' ', h('b', null, fmtNumber(dataset.n_genes))))
      if (geometry.mode) meta.push(h('span', { key: 'm' }, h('b', null, geometry.mode)))
      if (sliceColumn) meta.push(h('span', { key: 'slice' }, tf(T, 'datasetSlicesCount', sliceColumn.n_unique || '—')))
      if (clusterKeys.length) meta.push(h('span', { key: 'c' }, tf(T, 'datasetClustersCount', clusterKeys.length)))
      if (!isNil(geometry.neighbor_radius)) meta.push(h('span', { key: 'r' }, tf(T, 'datasetRadius', fmtNumber(geometry.neighbor_radius))))
      if (dataset.ref && dataset.ref.path) {
        meta.push(h('span', { key: 'ref' }, tf(T, 'datasetRef', (dataset.ref.cell_types || []).length || plainText(dataset.ref.path.split('/').pop(), 40))))
      }
      return h('div', { className: 'sp-dataset-card' },
        h('div', { className: 'sp-dataset-name', title: dataset.path || dataset.name }, dataset.name || dataset.id),
        h('div', { className: 'sp-dataset-meta' }, meta),
        h('div', { className: 'sp-dataset-meta' },
          h('span', null, tf(T, 'datasetStepsDone', dataset.steps.length)),
          dataset.current ? h('span', null, tf(T, 'datasetCurrentStep', (stepById(dataset.current) || { title: dataset.current }).title)) : null),
        dataset.geometry && dataset.geometry.image
          ? h('div', { className: 'sp-dataset-meta' }, h('span', { className: 'sp-tag sp-tag-ok' }, T('spatialOverlay')))
          : null,
        dataset.warnings.length
          ? h('div', { className: 'sp-warn-text', title: dataset.warnings.join('\n') }, T('warnings') + ' · ' + dataset.warnings.length)
          : null)
    }

    function StepCard(props) {
      var step = props.step
      var status = props.status
      var expandedState = useState(status === 'running' || status === 'failed')
      var expanded = expandedState[0]
      var setExpanded = expandedState[1]
      var statusKey = stepStatusKey(status)
      var dotClass = status === 'done'
        ? ' sp-step-dot-done'
        : (status === 'running' ? ' sp-step-dot-running' : (status === 'blocked' ? ' sp-step-dot-blocked' : ' sp-step-dot-pending'))
      var cardClass = 'sp-step' +
        (status === 'done' ? ' sp-step-done' : (status === 'running' ? ' sp-step-running' : (status === 'blocked' ? ' sp-step-blocked' : '')))
      var blockedRequires = []
      var done = doneStepSet(props.dataset)
      var requires = Array.isArray(step.requires) ? step.requires : []
      for (var index = 0; index < requires.length; index++) {
        if (!done[requires[index]]) {
          var required = stepById(requires[index])
          blockedRequires.push(required ? required.title : requires[index])
        }
      }
      var disabled = props.busy === true || !props.dataset || status === 'blocked'
      var summary = props.summary || ''
      var fallback = props.busy ? T('stepRunning')
        : (status === 'blocked' ? T('stepNeedPrereq') + blockedRequires.join('、') : T('stepSummaryEmpty'))
      return h('div', { className: cardClass, 'data-step': step.id },
        h('div', { className: 'sp-step-head' },
          h('span', { className: 'sp-step-dot' + dotClass, title: T(statusKey) }, String(props.index + 1)),
          h('span', { className: 'sp-step-title', title: step.title }, step.title),
          h('span', { className: 'sp-step-status sp-status-' + status }, T(statusKey)),
          h('button', {
            type: 'button',
            className: 'sp-step-caret',
            title: expanded ? T('analysisClose') : T('stepParams'),
            'aria-expanded': expanded ? 'true' : 'false',
            onClick: function () { setExpanded(!expanded) },
          }, h(Icon, { name: expanded ? 'chevronDown' : 'chevron', size: 12 }))),
        expanded
          ? h('div', { className: 'sp-step-params' },
              h(ParamForm, {
                params: step.params,
                values: props.values || {},
                columns: props.columns || [],
                disabled: props.busy === true,
                onChange: props.onParamChange,
              }))
          : null,
        h('div', { className: 'sp-step-foot' },
          Btn({
            small: true,
            variant: status === 'done' ? undefined : 'primary',
            icon: status === 'done' ? 'refresh' : 'play',
            disabled: disabled,
            title: status === 'blocked' ? T('stepNeedPrereq') + blockedRequires.join('、') : T('stepRun'),
            onClick: function () { if (props.onRun) props.onRun(step.id) },
          }, props.busy ? T('stepRunning') : (status === 'done' ? T('stepRerun') : T('stepRun'))),
          h('span', { className: 'sp-step-summary', title: summary || fallback }, summary || fallback)),
        props.error ? h(InlineError, { message: props.error }) : null)
    }

    function PipelinePanel(props) {
      return h('aside', { className: 'sp-left' },
        h('div', { className: 'sp-panel-head' },
          h('span', null, T('datasetsTitle')),
          h('span', { className: 'sp-panel-note' }, T('pipelineHint'))),
        h('div', { className: 'sp-panel-body', style: { paddingBottom: '8px' } },
          h(DatasetCard, { dataset: props.dataset })),
        h('div', { className: 'sp-steps' }, STEPS.map(function (step, index) {
          return h(StepCard, {
            key: step.id,
            step: step,
            index: index,
            status: props.statuses[step.id],
            dataset: props.dataset,
            values: props.values[step.id],
            columns: props.columns,
            summary: props.summaries[step.id],
            error: props.errors[step.id],
            busy: props.busyStep === step.id,
            onParamChange: function (patch) { if (props.onParamChange) props.onParamChange(step.id, patch) },
            onRun: props.onRun,
          })
        })))
    }

    // ==========================================================================
    // 右栏：Agent 面板（状态摘要 / 建议 + 复制指令 / 活动时间线）
    // ==========================================================================

    function AgentPanel(props) {
      var dataset = props.dataset
      var suggestions = Array.isArray(props.suggestions) ? props.suggestions : []
      var activity = Array.isArray(props.activity) ? props.activity : []
      var sections = []
      if (!dataset) {
        sections.push(h('div', { key: 'empty', className: 'sp-agent-section' },
          h('div', { className: 'sp-agent-title' }, T('agentStatus')),
          h('div', { className: 'sp-muted' }, T('agentNoDataset'))))
      } else {
        var currentStep = stepById(dataset.current)
        sections.push(h('div', { key: 'status', className: 'sp-agent-section' },
          h('div', { className: 'sp-agent-title' }, T('agentStatus')),
          h('div', { className: 'sp-stat-row' },
            h('span', null, dataset.name || dataset.id),
            h('b', null, fmtNumber(dataset.n_spots) + ' × ' + fmtNumber(dataset.n_genes))),
          h('div', { className: 'sp-stat-row' }, h('span', null, T('completed')), h('b', null, dataset.steps.length + ' / ' + STEPS.length)),
          h('div', { className: 'sp-stat-row' }, h('span', null, T('current')), h('b', null, currentStep ? currentStep.title : T('none'))),
          h('div', { className: 'sp-stat-row' }, h('span', null, T('spatialMode')), h('b', null, dataset.geometry.mode || T('unknown'))),
          h('div', { className: 'sp-stat-row' },
            h('span', null, T('engine')),
            h('b', null, props.health && props.health.ready ? T('engineReady') : T('engineBooting'))),
          h('div', { className: 'sp-muted' }, T('agentCoordsHint'))))
      }
      sections.push(h('div', { key: 'suggestions', className: 'sp-agent-section' },
        h('div', { className: 'sp-agent-title' }, T('agentSuggestions')),
        suggestions.length
          ? suggestions.map(function (suggestion) {
              var stepMeta = suggestion.step ? stepById(suggestion.step) : null
              var copied = props.copiedKey === ('suggestion-' + suggestion.id)
              return h('div', { key: suggestion.id, className: 'sp-suggestion' },
                h('div', { className: 'sp-suggestion-title' }, suggestion.title || T('suggestionFallback')),
                suggestion.detail ? h('div', { className: 'sp-suggestion-detail' }, suggestion.detail) : null,
                h('div', { className: 'sp-inline-row' },
                  stepMeta ? h('span', { className: 'sp-tag' }, stepMeta.title) : null,
                  Btn({
                    small: true,
                    variant: 'ghost',
                    icon: copied ? 'check' : 'copy',
                    onClick: function () { if (props.onCopy) props.onCopy(suggestionInstruction(suggestion), 'suggestion-' + suggestion.id) },
                  }, copied ? T('copied') : T('agentCopy'))))
            })
          : h('div', { className: 'sp-muted' }, T('noSuggestions'))))
      sections.push(h('div', { key: 'activity', className: 'sp-agent-section', style: { borderBottom: 'none' } },
        h('div', { className: 'sp-agent-title' }, T('agentActivity')),
        activity.length
          ? h('div', { className: 'sp-activity' }, activity.slice(0, 10).map(function (entry) {
              return h('div', { key: entry.id, className: 'sp-activity-item' },
                h('span', { className: 'sp-activity-dot' + (entry.ok ? '' : ' sp-activity-dot-bad') }),
                h('div', { className: 'sp-activity-main' },
                  h('div', { className: 'sp-activity-head' },
                    h('span', { className: 'sp-activity-actor' }, entry.actor === 'agent' ? T('activityAgent') : T('activityUser')),
                    h('span', null, entry.action || T('action')),
                    h('span', { className: 'sp-history-at' }, formatClock(entry.at))),
                  entry.text ? h('div', { className: 'sp-activity-text' }, entry.text) : null))
            }))
          : h('div', { className: 'sp-muted' }, T('agentNoActivity'))))
      return h('aside', { className: 'sp-right' },
        h('div', { className: 'sp-panel-head' },
          h('span', null, T('agentTitle')),
          h('span', { className: 'sp-panel-note' }, T('agentPanel'))),
        sections)
    }

    // ==========================================================================
    // 中栏 Tab 框架
    // ==========================================================================

    function TabsBar(props) {
      return h('nav', { className: 'sp-tabs', role: 'tablist' }, TAB_KEYS.map(function (key) {
        return h('button', {
          key: key,
          type: 'button',
          role: 'tab',
          'aria-selected': props.active === key ? 'true' : 'false',
          className: 'sp-tab' + (props.active === key ? ' sp-tab-on' : ''),
          onClick: function () { if (props.onChange) props.onChange(key) },
        }, T(TAB_META[key].titleKey))
      }))
    }

    function TabGuide(props) {
      return h(GuideBlock, {
        title: props.title,
        body: props.body,
        hint: props.hint,
        actions: props.action
          ? Btn({ small: true, variant: 'primary', icon: 'play', onClick: props.action.onClick, disabled: props.action.disabled }, props.action.label)
          : null,
      })
    }

    /** 坐标缩略图（preview.coords 可能为空：引擎只回采样坐标）。 */
    function MiniMap(props) {
      var points = Array.isArray(props.points) ? props.points : []
      if (!points.length) {
        return h('div', { className: 'sp-figure-empty' }, T('miniMapHint'))
      }
      var box = spatialScreen(points, props.bbox || null, 320, 220, { padding: 10 })
      return h('div', { className: 'sp-minimap' },
        h('svg', {
          viewBox: '0 0 320 220',
          width: '100%',
          height: '220',
          role: 'img',
          'aria-label': T('miniMapTitle'),
        }, box.map(function (point, index) {
          return h('circle', {
            key: 'p-' + index,
            cx: point.x,
            cy: point.y,
            r: 1.8,
            fill: colorForGroup(index % PALETTE.length),
            opacity: 0.85,
          })
        })))
    }

    function overviewColumn(dataset, keys) {
      var column = columnSummary(dataset, keys)
      return column && !isNil(column.n_unique) ? column.n_unique : null
    }

    function OverviewTab(props) {
      var dataset = props.dataset
      if (!dataset) {
        return h('div', { className: 'sp-pane' },
          h(TabGuide, {
            title: T('datasetNoDataset'),
            body: T('datasetOpenHint'),
            actions: null,
            action: { label: T('importData'), onClick: props.onOpenImport },
          }))
      }
      var preview = props.preview || normalizePreview({})
      var health = props.health || {}
      var geometry = dataset.geometry || {}
      var engineInfo = health.ready
        ? (health.mode || T('engineReady')) + (health.python ? ' · ' + T('pythonLabel') + ' ' + health.python : '')
        : T('engineBooting')
      var metrics = [
        { label: T('metricSpots'), value: fmtNumber(dataset.n_spots) },
        { label: T('metricGenes'), value: fmtNumber(dataset.n_genes) },
        { label: T('metricSpatialMode'), value: geometry.mode || T('unknown') },
        { label: T('metricSlices'), value: fmtNumber(overviewColumn(dataset, ['slice', 'slices'])) },
        { label: T('metricClusters'), value: fmtNumber((dataset.cluster_keys || []).length) },
        { label: T('metricDomains'), value: hasDomainKey(dataset) ? T('yes') : T('no') },
        { label: T('metricEngine'), value: engineInfo },
      ]
      var suggestions = Array.isArray(props.suggestions) ? props.suggestions : []
      var clusterReady = (dataset.cluster_keys || []).length > 0
      return h('div', { className: 'sp-pane' },
        h('div', { className: 'sp-toolbar' },
          h('span', { className: 'sp-section-title' }, T('overviewTitle')),
          h('span', { className: 'sp-muted' }, T('overviewHint')),
          h('span', { className: 'sp-spacer' }),
          props.busyStep ? h('span', { className: 'sp-loading' }, h(Spinner), h('span', null, props.busyStep)) : null,
          Btn({ small: true, icon: 'grid', onClick: props.onOpenAnalysis, disabled: props.busy === true }, T('analysisCatalog'))),
        h('div', { className: 'sp-metrics' }, metrics.map(function (metric) {
          return h(MetricCard, { key: metric.label, label: metric.label, value: metric.value })
        })),
        h('div', { className: 'sp-grid-2' },
          h(SectionCard, { title: T('clusterMap'), extra: clusterReady ? null : T('overviewClustersEmpty') },
            h(PlotBlock, {
              plots: props.plots,
              loadPlot: props.loadPlot,
              kind: 'spatial_scatter',
              params: { color_by: 'clusters' },
              enabled: !!dataset && clusterReady,
              emptyText: T('overviewClustersEmpty'),
            })),
          h(SectionCard, { title: T('pcaScree') },
            h(PlotBlock, {
              plots: props.plots,
              loadPlot: props.loadPlot,
              kind: 'pca_variance',
              params: {},
              enabled: !!dataset,
            }))),
        h('div', { className: 'sp-grid-2' },
          h(SectionCard, { title: T('qcMap') },
            h(PlotBlock, {
              plots: props.plots,
              loadPlot: props.loadPlot,
              kind: 'qc_spatial',
              params: {},
              enabled: !!dataset,
            })),
          h(SectionCard, { title: T('miniMapTitle'), extra: T('miniMapHint') },
            h(MiniMap, {
              points: preview.coords,
              bbox: geometry.x_range && geometry.y_range ? { x: geometry.x_range, y: geometry.y_range } : null,
            }))),
        h(SectionCard, { title: T('nextSteps'), extra: T('agentSuggestions') },
          suggestions.length
            ? h('div', { className: 'sp-form' }, suggestions.slice(0, 6).map(function (suggestion) {
                var stepMeta = suggestion.step ? stepById(suggestion.step) : null
                return h('div', { key: suggestion.id, className: 'sp-suggestion' },
                  h('div', { className: 'sp-suggestion-title' }, suggestion.title || T('suggestionFallback')),
                  suggestion.detail ? h('div', { className: 'sp-suggestion-detail' }, suggestion.detail) : null,
                  h('div', { className: 'sp-inline-row' },
                    stepMeta
                      ? Btn({
                          small: true,
                          variant: 'ghost',
                          icon: 'play',
                          disabled: props.stepBusy === stepMeta.id,
                          onClick: function () { if (props.onRunStep) props.onRunStep(stepMeta.id) },
                        }, stepMeta.title)
                      : null,
                    Btn({
                      small: true,
                      variant: 'ghost',
                      icon: 'copy',
                      onClick: function () { if (props.onCopy) props.onCopy(suggestionInstruction(suggestion), 'suggestion-' + suggestion.id) },
                    }, props.copiedKey === ('suggestion-' + suggestion.id) ? T('copied') : T('agentCopy'))))
              }))
            : h('div', { className: 'sp-muted' }, T('noSuggestions'))),
        dataset.warnings.length
          ? h(SectionCard, { title: T('warnings'), extra: String(dataset.warnings.length) },
              h('div', { className: 'sp-warn-list' }, dataset.warnings.map(function (warning, index) {
                return h('div', { key: 'w-' + index }, '· ' + warning)
              })))
          : null,
        dataset.notes.length
          ? h(SectionCard, { title: T('notes'), extra: String(dataset.notes.length) },
              h('div', { className: 'sp-note-list' }, dataset.notes.slice(0, 8).map(function (note, index) {
                return h('div', { key: 'n-' + index }, '· ' + note)
              })))
          : null)
    }

    // ==========================================================================
    // 中栏：空间图 / 表达 / 差异 / 空间统计 / 域与邻域 / 反卷积 / 通讯
    // ==========================================================================

    /** 单图请求：本地 busy/error/figure 三态（参数可变的图用这个）。 */
    function usePlotState() {
      var state = useState({ busy: false, error: '', figure: null })
      var entry = state[0]
      var setEntry = state[1]
      var run = function (request) {
        setEntry({ busy: true, error: '', figure: null })
        Promise.resolve()
          .then(request)
          .then(function (figure) {
            setEntry({ busy: false, error: figure ? '' : T('noFigure'), figure: figure || null })
            return figure
          })
          .catch(function (error) {
            setEntry({ busy: false, error: formatError(error), figure: null })
            return null
          })
      }
      return [entry, run]
    }

    function colorByOptions(dataset) {
      var options = []
      for (var index = 0; index < COLOR_BY.length; index++) {
        var item = COLOR_BY[index]
        if (item.value === 'domains' && !hasDomainKey(dataset)) continue
        options.push({ value: item.value, label: T(item.labelKey) })
      }
      return options
    }

    function SpatialMapTab(props) {
      var dataset = props.dataset
      var colorState = useState('clusters')
      var colorBy = colorState[0]
      var setColorBy = colorState[1]
      var geneState = useState('')
      var gene = geneState[0]
      var setGene = geneState[1]
      var sizeState = useState(6)
      var pointSize = sizeState[0]
      var setPointSize = sizeState[1]
      var maxState = useState(20000)
      var maxPoints = maxState[0]
      var setMaxPoints = maxState[1]
      var overlayState = useState(false)
      var overlay = overlayState[0]
      var setOverlay = overlayState[1]
      var plotState = usePlotState()
      var entry = plotState[0]
      var runPlot = plotState[1]
      var geometry = (dataset && dataset.geometry) || {}
      var hasImage = !!geometry.image
      var geneMode = colorBy === 'gene'
      function draw(withGene) {
        var params = {
          color_by: geneMode ? 'gene' : colorBy,
          point_size: asNumber(pointSize) || 6,
          max_points: asNumber(maxPoints) || 20000,
        }
        if (hasImage && overlay) params.image = true
        if (geneMode) {
          var target = plainText(withGene || gene)
          if (!target) return
          params.gene = target
          params.genes = [target]
        }
        runPlot(function () { return props.plotOnce('spatial_scatter', params) })
      }
      useEffect(function () {
        if (!dataset || !dataset.id) return
        if (!(dataset.cluster_keys || []).length) return
        runPlot(function () { return props.plotOnce('spatial_scatter', { color_by: 'clusters', point_size: 6, max_points: 20000 }) })
      }, [dataset ? dataset.id : ''])
      if (!dataset) {
        return h('div', { className: 'sp-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      }
      var options = colorByOptions(dataset)
      return h('div', { className: 'sp-pane' },
        h('div', { className: 'sp-toolbar' },
          h('span', { className: 'sp-section-title' }, T('spatialTitle')),
          h('span', { className: 'sp-muted' }, T('spatialHint')),
          h('span', { className: 'sp-spacer' }),
          h('span', { className: 'sp-toolbar-label' }, T('spatialColorBy')),
          h('div', { className: 'sp-color-row' }, options.map(function (option, index) {
            var on = colorBy === option.value
            return h('button', {
              key: option.value,
              type: 'button',
              className: 'sp-color-opt' + (on ? ' sp-color-opt-on' : ''),
              onClick: function () { setColorBy(option.value) },
            },
              h('span', { className: 'sp-color-dot', style: { background: PALETTE[index % PALETTE.length] } }),
              option.label)
          }))),
        h('div', { className: 'sp-toolbar' },
          h('span', { className: 'sp-toolbar-label' }, T('spatialGeneSearch')),
          h('input', {
            className: 'sp-input',
            value: gene,
            spellCheck: false,
            placeholder: T('spatialGenePh'),
            style: { minWidth: '160px' },
            onChange: function (event) { setGene(event.target.value) },
            onKeyDown: function (event) {
              if (event.key !== 'Enter') return
              event.preventDefault()
              if (!plainText(gene)) return
              draw(gene)
            },
          }),
          h('span', { className: 'sp-toolbar-label' }, T('spatialPointSize')),
          h('input', {
            className: 'sp-input',
            type: 'number',
            min: 1,
            max: 40,
            step: 1,
            value: String(pointSize),
            style: { width: '68px' },
            onChange: function (event) { setPointSize(event.target.value) },
          }),
          h('span', { className: 'sp-toolbar-label' }, T('spatialMaxPoints')),
          h('input', {
            className: 'sp-input',
            type: 'number',
            min: 100,
            max: 200000,
            step: 1000,
            value: String(maxPoints),
            style: { width: '90px' },
            onChange: function (event) { setMaxPoints(event.target.value) },
          }),
          h('label', { className: 'sp-check', title: hasImage ? T('spatialOverlay') : T('spatialOverlayMissing') },
            h('input', {
              type: 'checkbox',
              checked: overlay && hasImage,
              disabled: !hasImage,
              onChange: function (event) { setOverlay(event.target.checked) },
            }),
            h('span', null, T('spatialOverlay')))),
        !hasImage ? h('div', { className: 'sp-muted' }, T('spatialOverlayMissing')) : null,
        h('div', { className: 'sp-toolbar' },
          Btn({
            small: true,
            variant: 'primary',
            icon: 'chart',
            disabled: entry.busy,
            onClick: function () { draw() },
          }, entry.busy ? T('plotting') : T('spatialPlot')),
          Btn({
            small: true,
            icon: 'chart',
            disabled: entry.busy || !plainText(gene),
            title: plainText(gene) ? T('spatialGenePlot') : T('spatialNeedGene'),
            onClick: function () { draw(gene) },
          }, T('spatialGenePlot')),
          !(dataset.cluster_keys || []).length
            ? h('span', { className: 'sp-warn-text' }, T('spatialNeedCluster'))
            : null),
        h(FigureBlock, {
          title: T('spatialTitle'),
          figure: entry.figure,
          loading: entry.busy,
          error: entry.error,
          emptyText: T('noFigure'),
          onReload: function () { draw() },
        }))
    }

    function GenePicker(props) {
      var queryState = useState('')
      var query = queryState[0]
      var setQuery = queryState[1]
      var allGenes = Array.isArray(props.allGenes) ? props.allGenes : []
      var selected = Array.isArray(props.value) ? props.value : []
      var matches = useMemo(function () { return matchGenes(allGenes, query, 8) }, [allGenes.length, query])
      function addGene(gene) {
        var name = plainText(gene)
        if (!name || selected.indexOf(name) >= 0) return
        if (props.onChange) props.onChange(selected.concat([name]))
        setQuery('')
      }
      return h('div', { className: 'sp-form' },
        h('div', { className: 'sp-gene-row' },
          h('input', {
            className: 'sp-input',
            style: { flex: '1 1 180px' },
            value: query,
            spellCheck: false,
            placeholder: props.placeholder || T('exprGeneSearch'),
            onChange: function (event) { setQuery(event.target.value) },
            onKeyDown: function (event) {
              if (event.key !== 'Enter') return
              event.preventDefault()
              addGene(matches.length ? matches[0] : query)
            },
          }),
          Btn({ small: true, icon: 'plus', onClick: function () { addGene(matches.length ? matches[0] : query) } }, T('exprAddGene'))),
        query
          ? h('div', { className: 'sp-multi' }, matches.length
              ? matches.map(function (name) {
                  return h('button', {
                    key: name,
                    type: 'button',
                    className: 'sp-multi-opt',
                    title: T('genePickerHint'),
                    onClick: function () { addGene(name) },
                  }, name)
                })
              : h('span', { className: 'sp-muted' }, T('exprGeneNotInData')))
          : null,
        selected.length
          ? h('div', { className: 'sp-multi' }, selected.map(function (name) {
              return h('button', {
                key: name,
                type: 'button',
                className: 'sp-multi-opt sp-multi-opt-on',
                title: T('clearSelection'),
                onClick: function () { if (props.onChange) props.onChange(selected.filter(function (item) { return item !== name })) },
              }, name + ' ×')
            }))
          : h('div', { className: 'sp-muted' }, T('exprEmpty')))
    }

    function ExprTab(props) {
      var dataset = props.dataset
      var preview = props.preview || normalizePreview({})
      var genesState = useState([])
      var genes = genesState[0]
      var setGenes = genesState[1]
      var mapState = usePlotState()
      var boxState = usePlotState()
      var heatState = usePlotState()
      if (!dataset) {
        return h('div', { className: 'sp-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      }
      var firstGene = genes.length ? genes[0] : ''
      function plotMap() {
        if (!firstGene) return
        mapState[1](function () { return props.plotOnce('gene_spatial', { color_by: 'gene', gene: firstGene, genes: [firstGene] }) })
      }
      function plotBox() {
        if (!genes.length) return
        boxState[1](function () { return props.plotOnce('gene_box', { genes: genes }) })
      }
      function plotHeatmap() {
        if (!genes.length) return
        heatState[1](function () { return props.plotOnce('heatmap', { genes: genes }) })
      }
      var exported = (props.exports && props.exports.expr_csv) || {}
      return h('div', { className: 'sp-pane' },
        h('div', { className: 'sp-toolbar' },
          h('span', { className: 'sp-section-title' }, T('exprTitle')),
          h('span', { className: 'sp-muted' }, T('exprHint')),
          h('span', { className: 'sp-spacer' }),
          Btn({
            small: true,
            icon: 'download',
            disabled: exported.busy === true,
            onClick: function () { if (props.runExport) props.runExport('expr_csv') },
          }, exported.busy ? T('recordsExporting') : T('exprExport'))),
        exported.file ? h('span', { className: 'sp-file-chip' }, h(Icon, { name: 'check', size: 11 }), h('span', { className: 'sp-path' }, exported.file)) : null,
        exported.error ? h(InlineError, { message: exported.error }) : null,
        h(SectionCard, { title: T('exprSelected'), extra: T('exprHint') },
          h(GenePicker, { allGenes: preview.genes || [], value: genes, onChange: setGenes })),
        h('div', { className: 'sp-toolbar' },
          Btn({ small: true, variant: 'primary', icon: 'chart', disabled: mapState[0].busy || !firstGene, title: T('exprPlotMap'), onClick: plotMap },
            mapState[0].busy ? T('plotting') : T('exprPlotMap')),
          Btn({ small: true, variant: 'primary', icon: 'chart', disabled: boxState[0].busy || !genes.length, title: T('exprPlotBox'), onClick: plotBox },
            boxState[0].busy ? T('plotting') : T('exprPlotBox')),
          Btn({ small: true, variant: 'primary', icon: 'chart', disabled: heatState[0].busy || !genes.length, title: T('exprPlotHeatmap'), onClick: plotHeatmap },
            heatState[0].busy ? T('plotting') : T('exprPlotHeatmap'))),
        h('div', { className: 'sp-grid-2' },
          h(SectionCard, { title: T('exprPlotMap'), extra: firstGene || null },
            h(FigureBlock, {
              figure: mapState[0].figure,
              loading: mapState[0].busy,
              error: mapState[0].error,
              emptyText: T('exprEmpty'),
              onReload: plotMap,
            })),
          h(SectionCard, { title: T('exprPlotBox') },
            h(FigureBlock, {
              figure: boxState[0].figure,
              loading: boxState[0].busy,
              error: boxState[0].error,
              emptyText: T('exprEmpty'),
              onReload: plotBox,
            }))),
        h(SectionCard, { title: T('exprPlotHeatmap') },
          h(FigureBlock, {
            figure: heatState[0].figure,
            loading: heatState[0].busy,
            error: heatState[0].error,
            emptyText: T('exprEmpty'),
            onReload: plotHeatmap,
          })),
        h(Notice, { text: T('exprTableFallback') }),
        h('div', { className: 'sp-muted' }, T('exprExportHint')))
    }

    // --------------------------------------------------------------------------
    // 通用分析 tab：按服务端 /analyses 目录过滤 group，参数表单 + 运行 + 结果
    //   空间统计 / 域与邻域 / 反卷积 / 通讯 / 区域差异 全部复用这一个实现。
    // --------------------------------------------------------------------------

    /** 分析参数默认值（对象/数组复制，避免多次运行共享引用）。 */
    function analysisDefaults(descriptors) {
      var out = {}
      var list = Array.isArray(descriptors) ? descriptors : []
      for (var index = 0; index < list.length; index++) {
        var descriptor = list[index]
        if (!descriptor || !descriptor.key) continue
        out[descriptor.key] = cloneDefault(isNil(descriptor.default) ? null : descriptor.default)
      }
      return out
    }

    function groupAnalyses(catalog, tab) {
      var list = catalog && Array.isArray(catalog.analyses) ? catalog.analyses : []
      var groups = catalog && Array.isArray(catalog.groups) ? catalog.groups : []
      var grouped = {}
      for (var index = 0; index < list.length; index++) {
        var item = list[index]
        if (analysisGroupTab(item.group) !== tab) continue
        if (!grouped[item.group]) grouped[item.group] = []
        grouped[item.group].push(item)
      }
      var order = []
      for (var groupIndex = 0; groupIndex < groups.length; groupIndex++) {
        if (grouped[groups[groupIndex]]) order.push(groups[groupIndex])
      }
      var keys = Object.keys(grouped)
      for (var keyIndex = 0; keyIndex < keys.length; keyIndex++) {
        if (order.indexOf(keys[keyIndex]) < 0) order.push(keys[keyIndex])
      }
      return { groups: order, items: grouped }
    }

    function AnalysisTab(props) {
      var dataset = props.dataset
      var catalog = props.catalog || { loading: false, error: '', data: normalizeAnalyses({}) }
      var filtered = groupAnalyses(catalog.data, props.tab)
      var allItems = []
      for (var groupIndex = 0; groupIndex < filtered.groups.length; groupIndex++) {
        var items = filtered.items[filtered.groups[groupIndex]] || []
        for (var itemIndex = 0; itemIndex < items.length; itemIndex++) allItems.push(items[itemIndex])
      }
      var selectState = useState('')
      var selectedId = selectState[0]
      var setSelectedId = selectState[1]
      var valuesState = useState({})
      var values = valuesState[0]
      var setValues = valuesState[1]
      var meta = null
      for (var index = 0; index < allItems.length; index++) {
        if (allItems[index].id === selectedId) meta = allItems[index]
      }
      if (!meta && allItems.length) meta = allItems[0]
      useEffect(function () {
        if (!meta) return
        if (selectedId !== meta.id) {
          setSelectedId(meta.id)
          setValues(analysisDefaults(meta.params))
        }
      }, [dataset ? dataset.id : '', catalog.data, meta ? meta.id : ''])
      var run = (props.runs && meta && props.runs[meta.id]) || null
      function patch(next) {
        setValues(function (previous) { return Object.assign({}, previous, next || {}) })
      }
      function runAnalysis() {
        if (!meta || !dataset || props.onRunAnalysis === undefined) return
        props.onRunAnalysis(meta.id, values || {})
      }
      var groupNodes = filtered.groups.map(function (group) {
        var groupItems = filtered.items[group] || []
        if (!groupItems.length) return null
        return h('div', { key: 'group-' + group, className: 'sp-analysis-group' },
          h('div', { className: 'sp-analysis-group-title' }, group),
          groupItems.map(function (item) {
            var disabled = item.available === false
            var on = meta && meta.id === item.id
            return h('button', {
              key: item.id,
              type: 'button',
              className: 'sp-analysis-item' + (on ? ' sp-analysis-item-on' : ''),
              disabled: disabled,
              title: disabled ? (item.reason || T('analysisUnavailable')) : (item.method || item.title),
              onClick: function () {
                setSelectedId(item.id)
                setValues(analysisDefaults(item.params))
              },
            },
              h('div', { className: 'sp-analysis-main' },
                h('div', { className: 'sp-analysis-name' }, item.title || item.id),
                h('span', { className: 'sp-tag' + (disabled ? '' : ' sp-tag-ok') },
                  disabled ? T('analysisUnavailable') : T('statusOk'))),
              item.method ? h('div', { className: 'sp-analysis-method' }, item.method) : null,
              disabled && item.reason ? h('div', { className: 'sp-analysis-reason' }, item.reason) : null)
          }))
      })
      var runner = null
      if (!dataset) {
        runner = h(Notice, { text: T('analysisNeedDataset') })
      } else if (!allItems.length) {
        runner = h(Notice, { text: props.emptyText || T('analysisEmpty') })
      } else if (meta) {
        runner = h('div', { className: 'sp-analysis-runner' },
          h('div', { className: 'sp-analysis-head' },
            h('div', { className: 'sp-analysis-title' }, meta.title || meta.id),
            meta.method ? h('div', { className: 'sp-analysis-method' }, T('analysisMethod') + '：' + meta.method) : null),
          meta.needs && meta.needs.length
            ? h('div', { className: 'sp-inline-row' },
                h('span', { className: 'sp-toolbar-label' }, T('analysisNeeds')),
                meta.needs.map(function (need) {
                  var stepMeta = stepById(need)
                  return h('span', { key: 'need-' + need, className: 'sp-tag' }, stepMeta ? stepMeta.title : need)
                }))
            : null,
          meta.available === false && meta.reason ? h('div', { className: 'sp-warn-text' }, meta.reason) : null,
          h(ParamForm, {
            params: meta.params,
            values: values,
            columns: props.columns || [],
            disabled: run && run.busy === true,
            onChange: patch,
          }),
          h('div', { className: 'sp-inline-row' },
            Btn({
              small: true,
              variant: 'primary',
              icon: 'play',
              disabled: !dataset || (run && run.busy === true) || meta.available === false,
              onClick: runAnalysis,
            }, run && run.busy ? T('analysisRunning') : T('analysisRun')),
            props.extraActions || null),
          run && run.error ? h(InlineError, { message: run.error }) : null,
          run && run.result
            ? h(ResultView, {
                result: run.result,
                onExport: props.tabExport && props.runExport
                  ? function () { props.runExport(props.tabExport.what) }
                  : null,
                exportBusy: props.tabExport ? !!(props.exports && props.exports[props.tabExport.what] && props.exports[props.tabExport.what].busy) : false,
                exportLabel: props.tabExport ? T(props.tabExport.labelKey) : '',
                exportFile: props.tabExport && props.exports && props.exports[props.tabExport.what]
                  ? props.exports[props.tabExport.what].file || ''
                  : '',
              })
            : null)
      }
      return h('div', { className: 'sp-pane' },
        h('div', { className: 'sp-toolbar' },
          h('span', { className: 'sp-section-title' }, props.title),
          h('span', { className: 'sp-muted' }, props.hint),
          h('span', { className: 'sp-spacer' }),
          props.toolbarExtra || null),
        catalog.loading ? h(LoadingLine, { text: T('loading') }) : null,
        catalog.error ? h(InlineError, { message: catalog.error }) : null,
        h('div', { className: 'sp-analysis-wrap' },
          h('div', { className: 'sp-analysis-list' }, groupNodes),
          runner))
    }

    /** 结果导出按钮（各分析 tab 顶栏共用；what 与 DESIGN §2.2 export 一致）。 */
    function exportButton(props, what, label) {
      var state = (props && props.exports && props.exports[what]) || {}
      return Btn({
        small: true,
        icon: 'download',
        disabled: state.busy === true || !props.dataset,
        title: state.file || label,
        onClick: function () { if (props.runExport) props.runExport(what) },
      }, state.busy ? T('recordsExporting') : label)
    }

    function DeTab(props) {
      return h(AnalysisTab, Object.assign({}, props, {
        tab: 'de',
        tabExport: { what: 'de_csv', labelKey: 'deExport' },
        title: T('deTitle'),
        hint: T('deHint'),
        emptyText: T('deEmpty'),
        toolbarExtra: exportButton(props, 'de_csv', T('deExport')),
      }))
    }

    function SvgTab(props) {
      return h(AnalysisTab, Object.assign({}, props, {
        tab: 'stats',
        tabExport: { what: 'svg_csv', labelKey: 'exportSvgCsv' },
        title: T('statsTitle'),
        hint: T('statsHint'),
        emptyText: T('statsEmpty'),
        toolbarExtra: exportButton(props, 'svg_csv', T('exportSvgCsv')),
      }))
    }

    function DomainTab(props) {
      return h(AnalysisTab, Object.assign({}, props, {
        tab: 'domains',
        tabExport: { what: 'domain_csv', labelKey: 'exportDomainCsv' },
        title: T('domainsTitle'),
        hint: T('domainsHint'),
        emptyText: T('domainsEmpty'),
        toolbarExtra: exportButton(props, 'domain_csv', T('exportDomainCsv')),
      }))
    }

    function CommTab(props) {
      return h(AnalysisTab, Object.assign({}, props, {
        tab: 'comm',
        title: T('commTitle'),
        hint: T('commHint'),
        emptyText: T('commEmpty'),
      }))
    }

    /** 参考数据面板（反卷积用）：本机路径 / 单细胞工作台拉取 / 已挂载参考。 */
    function RefPanel(props) {
      var dataset = props.dataset
      var refsState = useState({ loading: true, error: '', available: false, datasets: [], local: [] })
      var refs = refsState[0]
      var setRefs = refsState[1]
      var pathState = useState('')
      var refPath = pathState[0]
      var setRefPath = pathState[1]
      var columnState = useState('cell_type')
      var cellTypeCol = columnState[0]
      var setCellTypeCol = columnState[1]
      var scState = useState('')
      var scDataset = scState[0]
      var setScDataset = scState[1]
      var attachState = useState({ busy: false, error: '', done: '' })
      var attach = attachState[0]
      var setAttach = attachState[1]
      var datasetId = dataset ? dataset.id : ''
      useEffect(function () {
        var cancelled = false
        setRefs(function (previous) { return Object.assign({}, previous, { loading: true, error: '' }) })
        props.api.refs()
          .then(function (payload) {
            if (cancelled) return
            var raw = payload && typeof payload === 'object' ? payload : {}
            var list = Array.isArray(raw.local) ? raw.local : []
            var local = []
            for (var index = 0; index < list.length; index++) {
              var item = list[index] && typeof list[index] === 'object' ? list[index] : {}
              var path = plainText(item.path || '')
              if (!path) continue
              local.push({
                path: path,
                name: plainText(item.name || path.split('/').pop()),
                bytes: asNumber(item.bytes),
                n_obs: asNumber(item.n_obs),
                n_vars: asNumber(item.n_vars),
              })
            }
            setRefs({
              loading: false,
              error: '',
              available: raw.available === true,
              datasets: Array.isArray(raw.datasets) ? raw.datasets : [],
              local: local,
            })
          })
          .catch(function (error) {
            if (cancelled) return
            setRefs({ loading: false, error: formatError(error), available: false, datasets: [], local: [] })
          })
        return function () { cancelled = true }
      }, [datasetId])
      function attachRef(payload) {
        if (!datasetId) return
        setAttach({ busy: true, error: '', done: '' })
        props.api.refAttach(Object.assign({ dataset: datasetId }, payload))
          .then(function () {
            setAttach({ busy: false, error: '', done: T('deconvRefDone') })
            if (props.onAttached) props.onAttached()
          })
          .catch(function (error) { setAttach({ busy: false, error: formatError(error), done: '' }) })
      }
      var current = dataset && dataset.ref ? dataset.ref : null
      return h(SectionCard, { title: T('deconvRefTitle'), extra: T('deconvRefHint') },
        h('div', { className: 'sp-ref-panel' },
          current
            ? h('div', { className: 'sp-inline-row' },
                h('span', { className: 'sp-tag sp-tag-ok' }, T('deconvRefAttached')),
                h('span', { className: 'sp-path', title: current.path }, current.path),
                current.cell_type_col ? h('span', { className: 'sp-tag' }, current.cell_type_col) : null,
                current.cell_types && current.cell_types.length
                  ? h('span', { className: 'sp-tag' }, tf(T, 'deconvRefCells', current.cell_types.length))
                  : null)
            : h('div', { className: 'sp-muted' }, T('deconvRefNone')),
          h('div', { className: 'sp-inline-row' },
            h('span', { className: 'sp-field-label' }, T('deconvRefPath')),
            h('input', {
              className: 'sp-input',
              style: { flex: '1 1 240px' },
              value: refPath,
              spellCheck: false,
              placeholder: T('deconvRefPathPh'),
              onChange: function (event) { setRefPath(event.target.value) },
            }),
            h('span', { className: 'sp-field-label' }, T('deconvCellTypeCol')),
            h('input', {
              className: 'sp-input',
              style: { width: '110px' },
              value: cellTypeCol,
              spellCheck: false,
              onChange: function (event) { setCellTypeCol(event.target.value) },
            }),
            Btn({
              small: true,
              variant: 'primary',
              icon: 'database',
              disabled: attach.busy || !datasetId || !plainText(refPath),
              onClick: function () {
                attachRef({ ref_path: plainText(refPath), cell_type_col: plainText(cellTypeCol) })
              },
            }, attach.busy ? T('deconvRefAttaching') : T('deconvRefAttach'))),
          refs.loading ? h(LoadingLine, { text: T('loading') }) : null,
          refs.error ? h(InlineError, { message: refs.error }) : null,
          refs.datasets.length
            ? h('div', { className: 'sp-inline-row' },
                h('span', { className: 'sp-field-label' }, T('deconvScList')),
                h('select', {
                  className: 'sp-input sp-select',
                  style: { flex: '1 1 200px' },
                  value: scDataset,
                  disabled: attach.busy || !datasetId,
                  onChange: function (event) { setScDataset(event.target.value) },
                }, [
                  h('option', { key: 'none', value: '' }, T('noSelection')),
                ].concat(refs.datasets.map(function (item) {
                  var id = plainText(item.id || '')
                  return h('option', { key: 'sc-' + id, value: id },
                    (item.name || id) + (item.n_obs ? ' · ' + fmtNumber(item.n_obs) + ' ' + T('spots') : ''))
                }))),
                Btn({
                  small: true,
                  icon: 'database',
                  disabled: attach.busy || !datasetId || !scDataset,
                  onClick: function () { attachRef({ sc_dataset: scDataset, cell_type_col: plainText(cellTypeCol) }) },
                }, T('deconvScPull')))
            : h('div', { className: 'sp-muted' }, T('deconvScUnavailable')),
          refs.local.length
            ? h('div', { className: 'sp-ref-list' }, refs.local.slice(0, 12).map(function (item) {
                var on = refPath === item.path
                return h('button', {
                  key: item.path,
                  type: 'button',
                  className: 'sp-ref-item' + (on ? ' sp-ref-item-on' : ''),
                  title: item.path,
                  onClick: function () {
                    setRefPath(item.path)
                    if (dataset) attachRef({ ref_path: item.path, cell_type_col: plainText(cellTypeCol) })
                  },
                },
                  h('span', { className: 'sp-ref-path' }, item.name),
                  h('span', { className: 'sp-tag' }, isNil(item.bytes) ? T('deconvLocalRefs') : formatBytes(item.bytes)))
              }))
            : h('div', { className: 'sp-muted' }, T('deconvNoRefs')),
          attach.error ? h(InlineError, { message: attach.error }) : null,
          attach.done ? h('div', { className: 'sp-ok-text' }, attach.done) : null))
    }

    function DeconvTab(props) {
      var dataset = props.dataset
      if (!dataset) {
        return h('div', { className: 'sp-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      }
      return h('div', { className: 'sp-pane' },
        h(RefPanel, { api: props.api, dataset: dataset, onAttached: props.onAttached }),
        dataset.ref && dataset.ref.path
          ? null
          : h(Notice, { text: T('deconvNeedRef') }),
        h(AnalysisTab, Object.assign({}, props, {
          tab: 'deconv',
        tabExport: { what: 'deconv_csv', labelKey: 'exportDeconvCsv' },
          title: T('deconvTitle'),
          hint: T('deconvHint'),
          emptyText: T('deconvEmpty'),
          toolbarExtra: exportButton(props, 'deconv_csv', T('exportDeconvCsv')),
        })))
    }

    // ==========================================================================
    // 中栏：记录（历史 + 自检 + 导出）/ 日志
    // ==========================================================================

    /**
     * 统计口径卡：数据集现状（spot/基因/切片·尺度·counts/坐标）+ 最近分析的设计/检验/校正
     * + 引用前必须处理的项。数据来自 Host 的 buildWorkbenchContract（取不到的字段为 null）。
     */
    function StatsContractCard(props) {
      var contract = props.contract || null
      function kv(label, value) {
        if (!value) return null
        return h('div', { className: 'sp-contract-row' },
          h('span', { className: 'sp-contract-k' }, label),
          h('span', { className: 'sp-contract-v' }, value))
      }
      function joinCounts(map) {
        var keys = Object.keys(map || {})
        if (!keys.length) return ''
        return keys.map(function (key) {
          var count = map[key]
          return key + ' ' + (count === null || count === undefined ? '?' : fmtNumber(count, 0))
        }).join(' / ')
      }
      var badge = contract
        ? h('span', { className: 'sp-badge ' + (contract.okToCite ? 'sp-badge-ok' : 'sp-badge-warn') },
            h('span', { className: 'sp-badge-dot' }),
            contract.okToCite ? T('contractOk') : T('contractWarn'))
        : null
      if (!contract) {
        return h(SectionCard, { title: T('contractTitle'), extra: badge },
          h('div', { className: 'sp-muted' }, T('contractHint')))
      }
      var dataset = contract.dataset || {}
      var design = contract.design || {}
      var inference = contract.inference
      var warns = (contract.caveats || []).filter(function (item) { return item.level === 'warn' })
      var infos = (contract.caveats || []).filter(function (item) { return item.level !== 'warn' })
      var spotsText = [
        dataset.nSpots !== null ? fmtNumber(dataset.nSpots, 0) : null,
        dataset.nGenes !== null ? fmtNumber(dataset.nGenes, 0) : null,
      ].filter(function (item) { return item !== null }).join(' / ')
      var sampleParts = []
      if (dataset.nSamples !== null) {
        sampleParts.push(fmtNumber(dataset.nSamples, 0) + (design.sampleKey ? '（' + design.sampleKey + '）' : ''))
      }
      if (design.nDonors !== null) {
        sampleParts.push(fmtNumber(design.nDonors, 0) + (design.donorKey ? '（' + design.donorKey + '）' : ''))
      }
      var samplesText = sampleParts.join(' · ')
      var groupText = [
        joinCounts(design.nPerGroup),
        design.groupKey ? '（' + design.groupKey + '）' : '',
      ].filter(Boolean).join(' ')
      var scaleText = [
        dataset.valueKind || '',
        dataset.species || '',
        dataset.hasCounts === true ? 'counts=yes' : (dataset.hasCounts === false ? 'counts=no' : ''),
        dataset.hasCoords === true ? 'coords=yes' : (dataset.hasCoords === false ? 'coords=no' : ''),
        dataset.spatialMode || '',
      ].filter(Boolean).join(' · ')
      var methodText = inference
        ? Object.keys(inference.method || {}).map(function (key) { return key + '=' + inference.method[key] }).join(' · ')
        : ''
      var designText = inference
        ? Object.keys(inference.designParams || {}).map(function (key) { return key + '=' + inference.designParams[key] }).join(' · ')
        : ''
      var fdrText = inference
        ? Object.keys(inference.fdr || {}).map(function (key) { return key + '=' + inference.fdr[key] }).join(' · ')
        : ''
      var inferenceText = inference
        ? [inference.lastStep, inference.ok ? T('completed') : T('failed'),
            inference.ms !== null ? formatDuration(inference.ms) : '',
            methodText, fdrText].filter(Boolean).join(' · ')
        : T('contractNoResult')
      return h(SectionCard, { title: T('contractTitle'), extra: badge },
        h('div', { className: 'sp-contract' },
          h('div', { className: 'sp-contract-hint' }, T('contractHint')),
          kv(T('contractSpots'), spotsText || null),
          kv(T('contractSamples'), samplesText),
          kv(T('contractGroups'), groupText || null),
          kv(T('contractScale'), scaleText || null),
          kv(T('contractInference'), inferenceText),
          kv(T('contractDesign'), designText || null),
          warns.length
            ? h('div', { className: 'sp-contract-warn' },
              h('div', { className: 'sp-contract-warn-title' }, T('contractCaveats') + '（' + warns.length + '）'),
              warns.map(function (item, index) {
                return h('div', { key: 'warn-' + index, className: 'sp-contract-warn-item' }, '· ' + item.text)
              }))
            : null,
          infos.length
            ? h('div', { className: 'sp-contract-info' },
              infos.map(function (item, index) { return h('div', { key: 'info-' + index }, '· ' + item.text) }))
            : null))
    }

    function RecordsTab(props) {
      var dataset = props.dataset
      var auditState = useState({ busy: false, error: '', audit: null })
      var audit = auditState[0]
      var setAudit = auditState[1]
      if (!dataset) {
        return h('div', { className: 'sp-pane' }, h(TabGuide, { title: T('datasetNoDataset'), body: T('datasetOpenHint') }))
      }
      function runAudit() {
        setAudit({ busy: true, error: '', audit: null })
        props.api.audit(dataset.id)
          .then(function (payload) { setAudit({ busy: false, error: '', audit: normalizeAudit(payload) }) })
          .catch(function (error) { setAudit({ busy: false, error: formatError(error), audit: null }) })
      }
      var history = (dataset.history || []).slice().reverse()
      var analysisRuns = (dataset.analyses_run || []).slice().reverse()
      return h('div', { className: 'sp-pane' },
        h('div', { className: 'sp-toolbar' },
          h('span', { className: 'sp-section-title' }, T('recordsTitle')),
          h('span', { className: 'sp-muted' }, T('recordsHint')),
          h('span', { className: 'sp-spacer' }),
          Btn({ small: true, icon: 'grid', onClick: props.onOpenAnalysis }, T('recordsAnalysis'))),
        h('div', { className: 'sp-grid-2' },
          h(SectionCard, { title: T('recordsHistory'), extra: tf(T, 'tableRows', history.length) },
            history.length
              ? h('div', { className: 'sp-history' }, history.map(function (entry, index) {
                  return h('div', { key: 'h-' + index, className: 'sp-history-item' },
                    h('span', { className: 'sp-history-at' }, formatClock(entry.at)),
                    h('span', { className: 'sp-history-step' }, (stepById(entry.step) || {}).title || entry.step || '—'),
                    h('span', { className: 'sp-history-summary', title: entry.summary || '' }, entry.summary || ''),
                    h('span', { className: 'sp-history-ms' }, entry.ms ? formatDuration(entry.ms) : ''),
                    entry.ok ? h('span', { className: 'sp-ok-text' }, T('completed')) : h('span', { className: 'sp-warn-text' }, T('failed')))
                }))
              : h('div', { className: 'sp-muted' }, T('recordsNoHistory'))),
          h(SectionCard, {
            title: T('recordsAudit'),
            extra: Btn({ small: true, icon: 'shield', disabled: audit.busy, onClick: runAudit },
              audit.busy ? T('running') : T('recordsRunAudit')),
          },
            audit.error ? h(InlineError, { message: audit.error }) : null,
            audit.audit
              ? h('div', { className: 'sp-form' },
                  h('div', { className: 'sp-score' },
                    h('span', { className: 'sp-score-v' }, String(audit.audit.score)),
                    h('span', { className: 'sp-score-label' },
                      T('recordsAuditScore') + ' · ' +
                      T(AUDIT_LEVEL_STYLE[audit.audit.level] ? AUDIT_LEVEL_STYLE[audit.audit.level].labelKey : 'auditOk'))),
                  h('div', { className: 'sp-muted' }, T('auditHint')),
                  audit.audit.checks.map(function (check) {
                    return h('div', {
                      key: check.id,
                      className: 'sp-audit-item' +
                        (check.level === 'warn' ? ' sp-audit-item-warn' : (check.level === 'fail' ? ' sp-audit-item-fail' : '')),
                    },
                      h('span', {
                        className: 'sp-audit-mark' +
                          (check.level === 'warn' ? ' sp-audit-mark-warn' : (check.level === 'fail' ? ' sp-audit-mark-fail' : '')),
                      }, check.level === 'ok' ? '✓' : (check.level === 'warn' ? '!' : '×')),
                      h('div', { className: 'sp-audit-main' },
                        h('span', { className: 'sp-audit-title' }, check.title || check.id),
                        check.detail ? h('span', { className: 'sp-audit-detail' }, check.detail) : null,
                        check.hint ? h('span', { className: 'sp-audit-detail' }, check.hint) : null))
                  }))
              : h('div', { className: 'sp-muted' }, T('recordsAuditEmpty')))),
        h(SectionCard, { title: T('recordsHistory'), extra: tf(T, 'tableRows', analysisRuns.length) },
          analysisRuns.length
            ? h('div', { className: 'sp-history' }, analysisRuns.map(function (entry, index) {
                return h('div', { key: 'a-' + index, className: 'sp-history-item' },
                  h('span', { className: 'sp-history-at' }, formatClock(entry.at)),
                  h('span', { className: 'sp-history-step' }, entry.id || '—'),
                  h('span', { className: 'sp-history-ms' }, entry.ms ? formatDuration(entry.ms) : ''),
                  entry.ok ? h('span', { className: 'sp-ok-text' }, T('completed')) : h('span', { className: 'sp-warn-text' }, T('failed')))
              }))
            : h('div', { className: 'sp-muted' }, T('recordsNoHistory'))),
        h(SectionCard, { title: T('recordsExport'), extra: T('exportBtn') },
          h('div', { className: 'sp-export-grid' }, EXPORT_KINDS.map(function (item) {
            var state = (props.exports && props.exports[item.id]) || {}
            return h('button', {
              key: item.id,
              type: 'button',
              className: 'sp-export-item',
              disabled: state.busy === true,
              title: state.file || T(item.labelKey),
              onClick: function () { if (props.runExport) props.runExport(item.id) },
            },
              h('span', null, state.busy ? T('recordsExporting') : T(item.labelKey)),
              h('span', { className: 'sp-export-hint' }, state.file ? state.file : item.hint),
              state.error ? h('span', { className: 'sp-warn-text' }, state.error) : null)
          }))),
        h(StatsContractCard, { contract: dataset ? dataset.stats_contract : null }),
        h(Notice, { text: T('auditHint') }))
    }

    function LogTab(props) {
      var state = useState({ loading: true, error: '', entries: [] })
      var log = state[0]
      var setLog = state[1]
      var mountedRef = useRef(false)
      function load() {
        setLog(function (previous) { return Object.assign({}, previous, { loading: true, error: '' }) })
        props.api.getActivity()
          .then(function (payload) {
            var list = Array.isArray(payload && payload.entries)
              ? payload.entries
              : (Array.isArray(payload && payload.activity) ? payload.activity : [])
            setLog({
              loading: false,
              error: '',
              entries: list.map(function (entry, index) { return normalizeActivityEntry(entry, index) }),
            })
          })
          .catch(function (error) { setLog({ loading: false, error: formatError(error), entries: [] }) })
      }
      useEffect(function () {
        if (mountedRef.current) return
        mountedRef.current = true
        load()
      }, [])
      return h('div', { className: 'sp-pane' },
        h('div', { className: 'sp-toolbar' },
          h('span', { className: 'sp-section-title' }, T('logTitle')),
          h('span', { className: 'sp-spacer' }),
          Btn({ small: true, icon: 'refresh', onClick: load, disabled: log.loading }, T('logRefresh'))),
        log.error ? h(InlineError, { message: log.error }) : null,
        log.loading
          ? h(LoadingLine, { text: T('loading') })
          : (log.entries.length
              ? h('div', { className: 'sp-table-wrap' },
                  h('table', { className: 'sp-table' },
                    h('thead', null, h('tr', null,
                      h('th', null, T('at')),
                      h('th', null, T('actor')),
                      h('th', null, T('action')),
                      h('th', null, T('dataset')),
                      h('th', null, T('duration')),
                      h('th', null, T('resultColumn')))),
                    h('tbody', null, log.entries.map(function (entry) {
                      return h('tr', { key: entry.id },
                        h('td', null, formatClock(entry.at)),
                        h('td', null, entry.actor === 'agent' ? T('activityAgent') : (entry.actor === 'system' ? T('engine') : T('activityUser'))),
                        h('td', null, entry.action || '—'),
                        h('td', null, entry.dataset || '—'),
                        h('td', { className: 'sp-td-num' }, entry.ms ? formatDuration(entry.ms) : '—'),
                        h('td', null, h('span', { className: entry.ok ? 'sp-ok-text' : 'sp-warn-text' }, entry.ok ? T('completed') : T('failed'))))
                    }))))
              : h(Notice, { text: T('logEmpty') })))
    }

    // ==========================================================================
    // GEO 检索弹层：搜索 → 详情 → 指定目录下载（轮询）→ 导入并质控
    // ==========================================================================

    function normalizeGeoItem(item) {
      var raw = item && typeof item === 'object' ? item : {}
      return {
        gse: plainText(raw.gse || raw.id || ''),
        title: plainText(raw.title || ''),
        titleZh: plainText(raw.title_zh || ''),
        summary: plainText(raw.summary || '', 600),
        species: plainText(raw.species || ''),
        nSamples: asNumber(raw.n_samples),
        platform: plainText(raw.platform || ''),
        platforms: Array.isArray(raw.platforms) ? raw.platforms.map(function (name) { return plainText(name) }) : [],
        studyType: plainText(raw.study_type || ''),
        release: plainText(raw.release || ''),
        hasMatrix: raw.has_matrix === true,
        spatialScore: plainText(raw.spatial_score || ''),
        spatialPlatform: Array.isArray(raw.spatial_platform)
          ? raw.spatial_platform.map(function (name) { return plainText(name) }).filter(Boolean)
          : [],
        tags: Array.isArray(raw.tags) ? raw.tags.map(function (tag) { return plainText(tag) }).filter(Boolean) : [],
      }
    }

    function normalizeGeoDetail(payload) {
      var raw = payload && typeof payload === 'object' ? payload : {}
      var files = []
      var source = Array.isArray(raw.files) ? raw.files : []
      for (var index = 0; index < source.length; index++) {
        var file = source[index] && typeof source[index] === 'object' ? source[index] : {}
        var name = plainText(file.name || '')
        if (!name) continue
        files.push({
          name: name,
          size: asNumber(file.size),
          url: plainText(file.url || ''),
          kind: plainText(file.kind || 'other'),
        })
      }
      var clinical = raw.clinical && typeof raw.clinical === 'object' ? raw.clinical : {}
      var fields = []
      var fieldSource = Array.isArray(clinical.fields) ? clinical.fields : []
      for (var fieldIndex = 0; fieldIndex < fieldSource.length; fieldIndex++) {
        var field = fieldSource[fieldIndex] && typeof fieldSource[fieldIndex] === 'object' ? fieldSource[fieldIndex] : {}
        fields.push({
          key: plainText(field.key || field.name || ''),
          n: asNumber(field.n),
          example: plainText(field.example || ''),
        })
      }
      var suggested = plainText(raw.suggested || '')
      if (!suggested) {
        var matrix = files.filter(function (item) { return item.kind === 'series_matrix' || /series_matrix/i.test(item.name) })[0]
        if (matrix) suggested = matrix.name
      }
      return {
        gse: plainText(raw.gse || ''),
        title: plainText(raw.title || ''),
        summary: plainText(raw.summary || '', 1200),
        design: plainText(raw.design || '', 600),
        species: plainText(raw.species || ''),
        nSamples: asNumber(raw.n_samples),
        platforms: Array.isArray(raw.platforms) ? raw.platforms.map(function (name) { return plainText(name) }) : [],
        pubmed: plainText(raw.pubmed || ''),
        release: plainText(raw.release || ''),
        visiumHint: plainText(raw.visium_hint || ''),
        clinical: { fields: fields },
        files: files,
        suggested: suggested,
      }
    }

    function GeoCard(props) {
      var item = props.item
      var scoreKey = item.spatialScore === 'high' ? 'geoScoreHigh' : (item.spatialScore === 'mid' ? 'geoScoreMid' : 'geoScoreLow')
      var scoreClass = item.spatialScore === 'high' ? ' sp-geo-tag-high' : (item.spatialScore === 'mid' ? ' sp-geo-tag-mid' : ' sp-geo-tag-low')
      return h('button', {
        type: 'button',
        className: 'sp-geo-card' + (props.active ? ' sp-geo-card-on' : ''),
        onClick: function () { if (props.onSelect) props.onSelect(item) },
      },
        h('div', { className: 'sp-inline-row' },
          h('span', { className: 'sp-geo-gse' }, item.gse),
          item.spatialScore ? h('span', { className: 'sp-geo-tag' + scoreClass }, T('geoScore') + ' ' + T(scoreKey)) : null,
          item.hasMatrix ? h('span', { className: 'sp-geo-tag sp-geo-tag-ok' }, T('geoMatrix')) : null),
        h('div', { className: 'sp-geo-title' }, item.titleZh || item.title),
        item.titleZh && item.title ? h('div', { className: 'sp-geo-title-en' }, item.title) : null,
        h('div', { className: 'sp-geo-meta' },
          !isNil(item.nSamples) ? h('span', null, item.nSamples + ' ' + T('samples')) : null,
          item.platform ? h('span', null, item.platform) : null,
          item.release ? h('span', null, item.release.slice(0, 4)) : null,
          item.species ? h('span', null, item.species) : null),
        item.spatialPlatform.length || item.tags.length
          ? h('div', { className: 'sp-geo-meta' },
              item.spatialPlatform.slice(0, 3).map(function (name) {
                return h('span', { key: 'p-' + name, className: 'sp-geo-tag sp-geo-tag-high' }, name)
              }),
              item.tags.slice(0, 4).map(function (tag) {
                return h('span', { key: 't-' + tag, className: 'sp-geo-tag' }, tag)
              }))
          : null)
    }

    function GeoDialog(props) {
      var filtersState = useState({
        term: '', species: 'all', study_type: 'all', min_samples: '', has_matrix: false, platform: 'all', limit: 20, page: 0,
      })
      var filters = filtersState[0]
      var setFilters = filtersState[1]
      var searchState = useState({ busy: false, error: '', total: 0, items: [] })
      var search = searchState[0]
      var setSearch = searchState[1]
      var detailState = useState({ gse: '', busy: false, error: '', detail: null })
      var detail = detailState[0]
      var setDetail = detailState[1]
      var selectedState = useState({})
      var selectedFiles = selectedState[0]
      var setSelectedFiles = selectedState[1]
      var dirState = useState(DEFAULT_GEO_DIR)
      var destDir = dirState[0]
      var setDestDir = dirState[1]
      var pickState = useState({ busy: false, error: '' })
      var pick = pickState[0]
      var setPick = pickState[1]
      var downloadState = useState({ busy: false, pct: 0, message: '', status: '', error: '', result: null })
      var download = downloadState[0]
      var setDownload = downloadState[1]
      var importState = useState({ busy: false, error: '' })
      var importInfo = importState[0]
      var setImportInfo = importState[1]
      var timerRef = useRef(null)
      var startRef = useRef(0)
      useEffect(function () {
        return function () {
          if (timerRef.current) { clearInterval(timerRef.current); timerRef.current = null }
        }
      }, [])
      function patchFilters(next) {
        setFilters(function (previous) { return Object.assign({}, previous, next || {}) })
      }
      function runSearch(event) {
        if (event && event.preventDefault) event.preventDefault()
        if (!plainText(filters.term)) return
        setSearch({ busy: true, error: '', notes: [], total: 0, items: [] })
        setDetail({ gse: '', busy: false, error: '', detail: null })
        props.api.geoSearch(filters)
          .then(function (response) {
            var items = (Array.isArray(response && response.items) ? response.items : []).map(normalizeGeoItem)
            setSearch({
              busy: false,
              error: '',
              notes: normalizeStringList(response && response.notes),
              total: asNumber(response && response.total) === null ? items.length : asNumber(response.total),
              items: items,
            })
          })
          .catch(function (error) { setSearch({ busy: false, error: formatError(error), notes: [], total: 0, items: [] }) })
      }
      function openDetail(item) {
        setDetail({ gse: item.gse, busy: true, error: '', detail: null })
        setSelectedFiles({})
        setDownload({ busy: false, pct: 0, message: '', status: '', error: '', result: null })
        props.api.geoDetail(item.gse)
          .then(function (payload) {
            var normalized = normalizeGeoDetail(payload)
            var defaults = {}
            if (normalized.suggested) defaults[normalized.suggested] = true
            setSelectedFiles(defaults)
            setDetail({ gse: item.gse, busy: false, error: '', detail: normalized })
          })
          .catch(function (error) { setDetail({ gse: item.gse, busy: false, error: formatError(error), detail: null }) })
      }
      function pollJob(jobId) {
        startRef.current = Date.now()
        if (timerRef.current) clearInterval(timerRef.current)
        timerRef.current = setInterval(function () {
          props.api.getJob(jobId)
            .then(function (payload) {
              var state = jobPolling(payload, { elapsedMs: Date.now() - startRef.current })
              if (state.running) {
                setDownload(function (previous) {
                  return Object.assign({}, previous, { pct: state.pct, message: state.message })
                })
                return
              }
              if (timerRef.current) { clearInterval(timerRef.current); timerRef.current = null }
              if (state.done) {
                setDownload({
                  busy: false, pct: 100, message: state.message || T('geoJobDone'),
                  status: 'done', error: '', result: state.result,
                })
                if (state.result && state.result.opened && props.onImported) props.onImported(state.result.opened)
              } else {
                setDownload({
                  busy: false, pct: state.pct, message: state.message,
                  status: state.failed ? 'failed' : 'timeout',
                  error: state.error || (state.failed ? T('geoJobFailed') : T('jobTimeout')),
                  result: null,
                })
              }
            })
            .catch(function (error) {
              if (timerRef.current) { clearInterval(timerRef.current); timerRef.current = null }
              setDownload({ busy: false, pct: 0, message: '', status: 'failed', error: formatError(error), result: null })
            })
        }, POLL_MS)
      }
      function startDownload(openAfter) {
        if (!detail.detail) return
        var files = Object.keys(selectedFiles).filter(function (name) { return selectedFiles[name] })
        if (!files.length) {
          setDownload({ busy: false, pct: 0, message: '', status: 'failed', error: T('geoNeedFile'), result: null })
          return
        }
        var dir = plainText(destDir)
        if (!dir || dir.charAt(0) !== '/') {
          setDownload({ busy: false, pct: 0, message: '', status: 'failed', error: T('geoInvalidDir'), result: null })
          return
        }
        setDownload({ busy: true, pct: 0, message: '', status: 'running', error: '', result: null })
        props.api.geoDownload({ gse: detail.gse, files: files, dest_dir: dir, open_after: openAfter === true })
          .then(function (payload) {
            var job = payload && payload.job && typeof payload.job === 'object' ? payload.job : payload
            if (job && job.id) { pollJob(job.id); return }
            if (job && (job.status === 'done' || job.status === 'completed')) {
              setDownload({ busy: false, pct: 100, message: T('geoJobDone'), status: 'done', error: '', result: job.result || null })
              return
            }
            setDownload({ busy: false, pct: 0, message: '', status: 'failed', error: T('geoJobFailed'), result: null })
          })
          .catch(function (error) {
            setDownload({ busy: false, pct: 0, message: '', status: 'failed', error: formatError(error), result: null })
          })
      }
      function pickDir() {
        setPick({ busy: true, error: '' })
        props.api.pickDir(T('geoPickTitle'))
          .then(function (payload) {
            var path = plainText(payload && (payload.path || payload.dir || payload.directory))
            if (path) setDestDir(path)
            setPick({ busy: false, error: path ? '' : T('geoPickFailed') })
          })
          .catch(function () { setPick({ busy: false, error: T('geoPickFailed') }) })
      }
      function importedFileFrom(result) {
        var files = result && Array.isArray(result.files) ? result.files : []
        var candidates = files.filter(function (file) { return file && file.ok !== false && file.path })
        var matrix = candidates.filter(function (file) { return /series_matrix/i.test(String(file.name || file.path)) })[0]
        return matrix || candidates[0] || null
      }
      function importAndQc() {
        var target = importedFileFrom(download.result)
        if (!target) {
          setImportInfo({ busy: false, error: T('geoNeedFile') })
          return
        }
        setImportInfo({ busy: true, error: '' })
        var options = { path: target.path, name: detail.gse || '' }
        if (/series_matrix/i.test(String(target.name || target.path))) options.data_format = 'geo_series'
        props.api.openDataset(options)
          .then(function (payload) {
            var dataset = payload && payload.dataset && typeof payload.dataset === 'object' ? payload.dataset : payload
            var datasetId = plainText(dataset && dataset.id ? dataset.id : (payload && payload.active))
            if (!datasetId) throw new Error(T('importError'))
            return props.api.runStep(datasetId, 'qc', defaultParams('qc')).then(function () {
              setImportInfo({ busy: false, error: '' })
              setDownload(function (previous) { return Object.assign({}, previous, { message: T('geoImportDone') }) })
              if (props.onImported) props.onImported(datasetId)
            })
          })
          .catch(function (error) { setImportInfo({ busy: false, error: formatError(error) }) })
      }
      var items = search.items
      var activeDetail = detail.detail
      var selectedCount = Object.keys(selectedFiles).filter(function (name) { return selectedFiles[name] }).length
      var busy = download.busy || importInfo.busy
      return h(Modal, {
        title: T('geoTitle'),
        subtitle: T('geoSubtitle'),
        wide: true,
        busy: busy,
        onClose: props.onClose,
      },
        h('div', { className: 'sp-modal-body' },
          h('form', { className: 'sp-geo-search', onSubmit: runSearch },
            h('input', {
              className: 'sp-input sp-geo-term',
              value: filters.term,
              placeholder: T('geoTermPh'),
              spellCheck: false,
              onChange: function (event) { patchFilters({ term: event.target.value }) },
            }),
            Btn({ variant: 'primary', icon: 'search', disabled: search.busy, onClick: runSearch },
              search.busy ? T('geoSearching') : T('geoSearchBtn'))),
          h('div', { className: 'sp-geo-filters' },
            h('label', { className: 'sp-field' },
              h('span', { className: 'sp-field-label' }, T('geoSpecies')),
              h('select', {
                className: 'sp-input sp-select',
                value: filters.species,
                onChange: function (event) { patchFilters({ species: event.target.value }) },
              },
                h('option', { value: 'all' }, T('geoAll')),
                h('option', { value: 'human' }, T('geoHuman')),
                h('option', { value: 'mouse' }, T('geoMouse')))),
            h('label', { className: 'sp-field' },
              h('span', { className: 'sp-field-label' }, T('geoStudyType')),
              h('select', {
                className: 'sp-input sp-select',
                value: filters.study_type,
                onChange: function (event) { patchFilters({ study_type: event.target.value }) },
              },
                h('option', { value: 'all' }, T('geoStudyAll')),
                h('option', { value: 'expression profiling by array' }, T('geoStudyArray')),
                h('option', { value: 'high throughput sequencing' }, T('geoStudySeq')))),
            h('label', { className: 'sp-field' },
              h('span', { className: 'sp-field-label' }, T('geoPlatform')),
              h('select', {
                className: 'sp-input sp-select',
                value: filters.platform,
                onChange: function (event) { patchFilters({ platform: event.target.value }) },
              }, GEO_PLATFORMS.map(function (option) {
                return h('option', { key: option.value, value: option.value },
                  option.labelKey ? T(option.labelKey) : option.label)
              }))),
            h('label', { className: 'sp-field' },
              h('span', { className: 'sp-field-label' }, T('geoMinSamples')),
              h('input', {
                className: 'sp-input',
                type: 'number',
                min: 0,
                style: { width: '76px' },
                value: filters.min_samples,
                onChange: function (event) { patchFilters({ min_samples: event.target.value }) },
              })),
            h('label', { className: 'sp-check' },
              h('input', {
                type: 'checkbox',
                checked: filters.has_matrix,
                onChange: function (event) { patchFilters({ has_matrix: event.target.checked }) },
              }),
              h('span', null, T('geoHasMatrix')))),
          h('div', { className: 'sp-muted' }, T('geoHint')),
          search.error ? h(InlineError, { message: T('geoErrorSearch') + '：' + search.error }) : null,
          (search.notes && search.notes.length)
            ? h('div', { className: 'sp-muted', style: { marginTop: '6px' } },
                search.notes.map(function (note, index) { return h('div', { key: 'gn' + index }, '· ' + note) }))
            : null,
          search.busy
            ? h(LoadingLine, { text: T('geoSearching') })
            : h('div', { className: 'sp-geo-layout' },
                h('div', null,
                  h('div', { className: 'sp-section-title' },
                    h('span', null, T('geoResults')),
                    h('span', { className: 'sp-muted' }, items.length ? tf(T, 'geoTotal', search.total) : '')),
                  items.length
                    ? h('div', { className: 'sp-geo-list' }, items.map(function (item) {
                        return h(GeoCard, {
                          key: item.gse,
                          item: item,
                          active: detail.gse === item.gse,
                          onSelect: openDetail,
                        })
                      }))
                    : h('div', { className: 'sp-muted' }, T('geoNoResults'))),
                h('div', { className: 'sp-geo-detail' },
                  detail.busy ? h(LoadingLine, { text: T('geoLoadingDetail') }) : null,
                  detail.error ? h(InlineError, { message: T('geoErrorDetail') + '：' + detail.error }) : null,
                  activeDetail
                    ? h(React.Fragment, null,
                        h('div', { className: 'sp-geo-detail-title' }, activeDetail.title || detail.gse),
                        h('dl', { className: 'sp-geo-kv' },
                          h('dt', null, T('dataset')), h('dd', null, activeDetail.gse),
                          h('dt', null, T('geoPlatform')), h('dd', null, (activeDetail.platforms || []).join(', ') || '—'),
                          h('dt', null, T('samples')), h('dd', null, isNil(activeDetail.nSamples) ? '—' : String(activeDetail.nSamples)),
                          h('dt', null, T('geoRelease')), h('dd', null, activeDetail.release || '—'),
                          h('dt', null, T('geoPubmed')), h('dd', null, activeDetail.pubmed || '—')),
                        activeDetail.summary
                          ? h('div', null,
                              h('div', { className: 'sp-field-label' }, T('geoSummary')),
                              h('div', { className: 'sp-muted' }, activeDetail.summary))
                          : null,
                        activeDetail.design
                          ? h('div', null,
                              h('div', { className: 'sp-field-label' }, T('geoDesign')),
                              h('div', { className: 'sp-muted' }, activeDetail.design))
                          : null,
                        activeDetail.visiumHint
                          ? h('div', { className: 'sp-warn-text' }, activeDetail.visiumHint)
                          : null,
                        activeDetail.clinical.fields.length
                          ? h('div', null,
                              h('div', { className: 'sp-field-label' }, T('geoFiles') + ' · ' + activeDetail.clinical.fields.length),
                              h('div', { className: 'sp-multi' }, activeDetail.clinical.fields.slice(0, 12).map(function (field) {
                                return h('span', { key: field.key, className: 'sp-geo-tag' }, field.key + (field.n ? ' · ' + field.n : ''))
                              })))
                          : null,
                        h('div', { className: 'sp-field-label' }, T('geoFiles') + ' · ' + activeDetail.files.length),
                        activeDetail.files.length
                          ? h('div', { className: 'sp-form' }, activeDetail.files.map(function (file) {
                              var on = selectedFiles[file.name] === true
                              return h('label', { key: file.name, className: 'sp-geo-file' },
                                h('input', {
                                  type: 'checkbox',
                                  checked: on,
                                  onChange: function (event) {
                                    var next = Object.assign({}, selectedFiles)
                                    next[file.name] = event.target.checked
                                    setSelectedFiles(next)
                                  },
                                }),
                                h('span', { className: 'sp-geo-file-name', title: file.name }, file.name),
                                h('span', { className: 'sp-muted' }, file.kind || ''),
                                h('span', { className: 'sp-muted' }, isNil(file.size) ? '' : formatBytes(file.size)))
                            }))
                          : h('div', { className: 'sp-muted' }, T('geoNoFiles')),
                        h('div', { className: 'sp-inline-row' },
                          Btn({
                            small: true,
                            icon: 'check',
                            onClick: function () {
                              var all = {}
                              activeDetail.files.forEach(function (file) { all[file.name] = true })
                              setSelectedFiles(all)
                            },
                          }, T('geoSelectAll')),
                          Btn({ small: true, icon: 'close', onClick: function () { setSelectedFiles({}) } }, T('geoClear')),
                          h('span', { className: 'sp-muted' }, selectedCount + ' / ' + activeDetail.files.length)))
                    : h('div', { className: 'sp-muted' }, T('geoDetailTitle') + ' · ' + T('geoDetailBtn')))),
          h('div', { className: 'sp-divider' }),
          h('div', { className: 'sp-inline-row' },
            h('span', { className: 'sp-field-label' }, T('geoDownloadDir')),
            h('input', {
              className: 'sp-input',
              style: { flex: '1 1 260px' },
              value: destDir,
              spellCheck: false,
              placeholder: DEFAULT_GEO_DIR,
              title: T('geoDownloadDirHint'),
              onChange: function (event) { setDestDir(event.target.value) },
            }),
            Btn({ small: true, icon: 'folder', disabled: pick.busy, onClick: pickDir },
              pick.busy ? T('geoPicking') : T('geoPickDir')),
            h('span', { className: 'sp-muted' }, T('geoDownloadDirHint'))),
          pick.error ? h('InlineError', { message: pick.error }) : null,
          download.busy || download.status
            ? h('div', { className: 'sp-progress' },
                h('div', { className: 'sp-progress-bar' },
                  h('div', { className: 'sp-progress-fill', style: { width: clamp(download.pct, 0, 100) + '%' } })),
                h('div', { className: 'sp-progress-text' },
                  h('span', null, download.message || (download.busy
                    ? T('geoDownloading')
                    : (download.status === 'done' ? T('geoJobDone') : (download.status === 'failed' ? T('geoJobFailed') : '')))),
                  h('span', null, Math.round(clamp(download.pct, 0, 100)) + '%')))
            : null,
          download.error ? h(InlineError, { message: download.error }) : null,
          importInfo.error ? h(InlineError, { message: importInfo.error }) : null),
        h('div', { className: 'sp-modal-foot' },
          Btn({
            variant: 'primary',
            icon: 'download',
            disabled: !activeDetail || download.busy || selectedCount === 0 || !plainText(destDir),
            onClick: function () { startDownload(false) },
          }, download.busy ? T('geoDownloading') : T('geoDownloadBtn')),
          Btn({
            icon: 'importIcon',
            title: T('geoOpenAfter'),
            disabled: !activeDetail || download.busy || selectedCount === 0 || !plainText(destDir),
            onClick: function () { startDownload(true) },
          }, T('geoOpenAfter')),
          download.result
            ? Btn({ icon: 'importIcon', disabled: importInfo.busy, onClick: importAndQc },
                importInfo.busy ? T('geoImporting') : T('geoImportBtn'))
            : null,
          Btn({ icon: 'close', onClick: props.onClose, disabled: busy }, T('close'))))
    }

    // ==========================================================================
    // 导入弹层：本机路径 / URL / 上传 + 格式 + 坐标/底图 + 高级项 + 候选点选
    // ==========================================================================

    function formatOptions(payload) {
      var raw = payload && typeof payload === 'object' ? payload : {}
      var list = Array.isArray(raw) ? raw : (Array.isArray(raw.formats) ? raw.formats : [])
      var out = []
      for (var index = 0; index < list.length; index++) {
        var item = list[index]
        if (typeof item === 'string') {
          out.push({ value: item, label: item })
          continue
        }
        var format = item && typeof item === 'object' ? item : {}
        var value = plainText(format.id || format.value || format.name || '')
        if (!value) continue
        out.push({ value: value, label: plainText(format.label || format.title || format.id || format.name || value) })
      }
      if (!out.length) {
        for (var fallbackIndex = 0; fallbackIndex < DATA_FORMATS.length; fallbackIndex++) {
          out.push({ value: DATA_FORMATS[fallbackIndex], label: DATA_FORMATS[fallbackIndex] })
        }
      }
      return out
    }

    function ImportDialog(props) {
      var sourceState = useState('local')
      var source = sourceState[0]
      var setSource = sourceState[1]
      var formState = useState({
        path: '', url: '', name: '', data_format: 'auto', coords: '', image: '',
        sheet: '', gene_col: '', member: '', object: '', transpose: false, sample_sheet: '',
      })
      var form = formState[0]
      var setForm = formState[1]
      var fileState = useState(null)
      var file = fileState[0]
      var setFile = fileState[1]
      var busyState = useState({ import: false, inspect: false, error: '', inspectError: '', result: null, candidates: [] })
      var busy = busyState[0]
      var setBusy = busyState[1]
      function patchForm(next) {
        setForm(function (previous) { return Object.assign({}, previous, next || {}) })
      }
      var formats = formatOptions(props.formats)
      function payloadOf() {
        var payload = Object.assign({}, form)
        if (source === 'url') payload.path = ''
        if (source === 'local') payload.url = ''
        return payload
      }
      function handleAmbiguous(error) {
        var candidates = Array.isArray(error && error.candidates) ? error.candidates : []
        var names = []
        for (var index = 0; index < candidates.length; index++) {
          var candidate = candidates[index]
          var name = typeof candidate === 'string'
            ? candidate
            : plainText(candidate && (candidate.name || candidate.path || candidate.sheet || candidate.id))
          if (name) names.push(name)
        }
        setBusy(function (previous) {
          return Object.assign({}, previous, { error: formatError(error), candidates: names })
        })
      }
      function importDone(payload) {
        setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: '', result: payload }) })
        if (props.onImported) props.onImported(payload)
      }
      function submit() {
        if (source === 'upload') {
          if (!file) {
            setBusy(function (previous) { return Object.assign({}, previous, { error: T('importNoPath'), candidates: [] }) })
            return
          }
          setBusy(function (previous) { return Object.assign({}, previous, { import: true, error: '', candidates: [] }) })
          Promise.resolve(typeof file.arrayBuffer === 'function' ? file.arrayBuffer() : file)
            .then(function (buffer) {
              return props.api.request('/dataset/upload', {
                method: 'POST',
                headers: { 'x-dsh-filename': file.name || 'upload.bin', 'content-type': 'application/octet-stream' },
                body: buffer,
              }, { timeoutMs: 600000 })
            })
            .then(importDone)
            .catch(function (error) {
              if (error && error.code === 'ambiguous') { handleAmbiguous(error); return }
              setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: formatError(error), candidates: [] }) })
            })
          return
        }
        var payload = payloadOf()
        if (!plainText(payload.path) && !plainText(payload.url)) {
          setBusy(function (previous) { return Object.assign({}, previous, { error: T('importNoPath'), candidates: [] }) })
          return
        }
        setBusy(function (previous) { return Object.assign({}, previous, { import: true, error: '', candidates: [] }) })
        props.api.openDataset(payload)
          .then(importDone)
          .catch(function (error) {
            if (error && error.code === 'ambiguous') { handleAmbiguous(error); return }
            setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: formatError(error), candidates: [] }) })
          })
      }
      function inspect() {
        var payload = payloadOf()
        if (!plainText(payload.path) && !plainText(payload.url)) {
          setBusy(function (previous) { return Object.assign({}, previous, { inspectError: T('importNoPath') }) })
          return
        }
        setBusy(function (previous) { return Object.assign({}, previous, { inspect: true, inspectError: '', result: null }) })
        props.api.inspectSource(payload)
          .then(function (response) {
            setBusy(function (previous) { return Object.assign({}, previous, { inspect: false, inspectError: '', result: response }) })
          })
          .catch(function (error) {
            setBusy(function (previous) { return Object.assign({}, previous, { inspect: false, inspectError: formatError(error) }) })
          })
      }
      function loadDemo() {
        setBusy(function (previous) { return Object.assign({}, previous, { import: true, error: '', candidates: [] }) })
        props.api.demo()
          .then(importDone)
          .catch(function (error) {
            setBusy(function (previous) { return Object.assign({}, previous, { import: false, error: formatError(error), candidates: [] }) })
          })
      }
      function pickCandidate(candidate) {
        if (/\.xlsx?$/i.test(candidate)) patchForm({ sheet: candidate })
        else patchForm({ member: candidate })
        setBusy(function (previous) { return Object.assign({}, previous, { candidates: [] }) })
      }
      var sourceTabs = [
        { id: 'local', labelKey: 'importSourceLocal' },
        { id: 'url', labelKey: 'importSourceUrl' },
        { id: 'upload', labelKey: 'importChooseFile' },
      ]
      var inspected = busy.result && typeof busy.result === 'object' ? busy.result : null
      return h(Modal, {
        title: T('importTitle'),
        subtitle: T('importSubtitle'),
        busy: busy.import,
        onClose: props.onClose,
      },
        h('div', { className: 'sp-modal-body' },
          h('div', { className: 'sp-inline-row' }, sourceTabs.map(function (item) {
            return h('button', {
              key: item.id,
              type: 'button',
              className: 'sp-multi-opt' + (source === item.id ? ' sp-multi-opt-on' : ''),
              onClick: function () { setSource(item.id) },
            }, T(item.labelKey))
          })),
          source === 'upload'
            ? h('div', { className: 'sp-field sp-field-col' },
                h('span', { className: 'sp-field-label' }, T('importChooseFile')),
                h('input', {
                  type: 'file',
                  className: 'sp-input',
                  style: { height: 'auto', padding: '5px 8px' },
                  onChange: function (event) { setFile(event.target.files && event.target.files[0] ? event.target.files[0] : null) },
                }),
                h('div', { className: 'sp-muted' }, T('importUploadHint')))
            : h('div', { className: 'sp-field sp-field-col' },
                h('span', { className: 'sp-field-label' }, source === 'url' ? T('importSourceUrl') : T('importSourceLocal')),
                h('input', {
                  className: 'sp-input',
                  value: source === 'url' ? form.url : form.path,
                  spellCheck: false,
                  placeholder: source === 'url' ? T('importUrlPh') : T('importPathPh'),
                  onChange: function (event) {
                    patchForm(source === 'url' ? { url: event.target.value } : { path: event.target.value })
                  },
                })),
          h('div', { className: 'sp-form-grid' },
            h(FormField, { label: T('importFormat') },
              h('select', {
                className: 'sp-input sp-select',
                value: form.data_format,
                onChange: function (event) { patchForm({ data_format: event.target.value }) },
              }, h('option', { value: 'auto' }, T('importAuto')), formats.map(function (item) {
                return h('option', { key: item.value, value: item.value }, item.label)
              }))),
            h(FormField, { label: T('importCoords') },
              h('input', {
                className: 'sp-input',
                value: form.coords,
                spellCheck: false,
                placeholder: T('importCoordsPh'),
                onChange: function (event) { patchForm({ coords: event.target.value }) },
              })),
            h(FormField, { label: T('importImage') },
              h('input', {
                className: 'sp-input',
                value: form.image,
                spellCheck: false,
                placeholder: T('importImagePh'),
                onChange: function (event) { patchForm({ image: event.target.value }) },
              }))),
          h('div', { className: 'sp-muted' }, T('importPickHint')),
          h('details', { className: 'sp-advanced' },
            h('summary', null, T('importAdvanced')),
            h('div', { className: 'sp-form-grid', style: { marginTop: '8px' } },
              h(FormField, { label: T('importSheet') },
                h('input', {
                  className: 'sp-input', value: form.sheet, spellCheck: false,
                  onChange: function (event) { patchForm({ sheet: event.target.value }) },
                })),
              h(FormField, { label: T('importGeneCol') },
                h('input', {
                  className: 'sp-input', value: form.gene_col, spellCheck: false,
                  onChange: function (event) { patchForm({ gene_col: event.target.value }) },
                })),
              h(FormField, { label: T('importMember') },
                h('input', {
                  className: 'sp-input', value: form.member, spellCheck: false,
                  onChange: function (event) { patchForm({ member: event.target.value }) },
                })),
              h(FormField, { label: T('importObject') },
                h('input', {
                  className: 'sp-input', value: form.object, spellCheck: false,
                  onChange: function (event) { patchForm({ object: event.target.value }) },
                })),
              h(FormField, { label: T('importSampleSheet') },
                h('input', {
                  className: 'sp-input', value: form.sample_sheet, spellCheck: false, placeholder: '/absolute/path/samples.csv',
                  onChange: function (event) { patchForm({ sample_sheet: event.target.value }) },
                })),
              h('label', { className: 'sp-check' },
                h('input', {
                  type: 'checkbox',
                  checked: form.transpose,
                  onChange: function (event) { patchForm({ transpose: event.target.checked }) },
                }),
                h('span', null, T('importTranspose'))))),
          h('div', { className: 'sp-muted' }, T('importHint')),
          busy.error ? h(InlineError, { message: T('importError') + '：' + busy.error }) : null,
          busy.candidates.length
            ? h('div', { className: 'sp-card' },
                h('div', { className: 'sp-card-title' }, T('importCandidates')),
                h('div', { className: 'sp-multi' }, busy.candidates.map(function (candidate) {
                  return h('button', {
                    key: candidate,
                    type: 'button',
                    className: 'sp-multi-opt',
                    onClick: function () { pickCandidate(candidate) },
                  }, candidate)
                })))
            : null,
          busy.inspectError ? h(InlineError, { message: busy.inspectError }) : null,
          inspected && (inspected.format || inspected.notes)
            ? h('div', { className: 'sp-card' },
                h('div', { className: 'sp-card-title' }, T('importInspectTitle')),
                h('div', { className: 'sp-inline-row' },
                  inspected.format ? h('span', { className: 'sp-tag sp-tag-ok' }, String(inspected.format)) : null,
                  !isNil(inspected.n_spots) ? h('span', { className: 'sp-tag' }, fmtNumber(inspected.n_spots) + ' ' + T('spots')) : null,
                  !isNil(inspected.n_genes) ? h('span', { className: 'sp-tag' }, fmtNumber(inspected.n_genes) + ' ' + T('genes')) : null,
                  !isNil(inspected.size) ? h('span', { className: 'sp-tag' }, formatBytes(inspected.size)) : null),
                normalizeStringList(inspected.notes).length
                  ? h('div', { className: 'sp-note-list' }, normalizeStringList(inspected.notes).map(function (note, index) {
                      return h('div', { key: 'n-' + index }, '· ' + note)
                    }))
                  : null)
            : null),
        h('div', { className: 'sp-modal-foot' },
          Btn({ small: true, icon: 'flask', disabled: busy.import || busy.inspect, title: T('demoLoad'), onClick: loadDemo }, T('demoLoad')),
          Btn({ small: true, icon: 'search', disabled: busy.import || busy.inspect, onClick: inspect },
            busy.inspect ? T('importInspectBusy') : T('importInspect')),
          Btn({
            variant: 'primary',
            icon: 'importIcon',
            disabled: busy.import || (source === 'upload' ? !file : (!plainText(form.path) && !plainText(form.url))),
            onClick: submit,
          }, busy.import ? T('importBusy') : T('importStart')),
          Btn({ icon: 'close', onClick: props.onClose, disabled: busy.import }, T('close'))))
    }

    // ==========================================================================
    // 分析目录抽屉：/analyses 按 group 分组 + 参数表单 + 运行 + 结果
    // ==========================================================================

    function AnalysisDrawer(props) {
      var catalog = props.catalog || { loading: false, error: '', data: normalizeAnalyses({}) }
      var dataset = props.dataset
      var data = catalog.data
      var selectState = useState('')
      var selectedId = selectState[0]
      var setSelectedId = selectState[1]
      var valuesState = useState({})
      var values = valuesState[0]
      var setValues = valuesState[1]
      var meta = null
      var list = Array.isArray(data.analyses) ? data.analyses : []
      for (var index = 0; index < list.length; index++) {
        if (list[index].id === selectedId) meta = list[index]
      }
      function select(item) {
        if (!item || item.available === false) return
        setSelectedId(item.id)
        setValues(analysisDefaults(item.params))
      }
      function patch(next) {
        setValues(function (previous) { return Object.assign({}, previous, next || {}) })
      }
      var grouped = groupBy(list, 'group')
      var groupNames = Array.isArray(data.groups) ? data.groups.slice() : []
      var extraGroups = Object.keys(grouped)
      for (var groupIndex = 0; groupIndex < extraGroups.length; groupIndex++) {
        if (groupNames.indexOf(extraGroups[groupIndex]) < 0) groupNames.push(extraGroups[groupIndex])
      }
      var run = (props.runs && meta && props.runs[meta.id]) || null
      var groupNodes = groupNames.map(function (group) {
        var items = grouped[group] || []
        if (!items.length) return null
        return h('div', { key: 'group-' + group, className: 'sp-analysis-group' },
          h('div', { className: 'sp-analysis-group-title' }, group),
          items.map(function (item) {
            var disabled = item.available === false
            return h('button', {
              key: item.id,
              type: 'button',
              className: 'sp-analysis-item' + (meta && meta.id === item.id ? ' sp-analysis-item-on' : ''),
              disabled: disabled,
              title: disabled ? (item.reason || T('analysisUnavailable')) : (item.method || item.title),
              onClick: function () { select(item) },
            },
              h('div', { className: 'sp-analysis-main' },
                h('div', { className: 'sp-analysis-name' }, item.title || item.id),
                h('span', { className: 'sp-tag' + (disabled ? '' : ' sp-tag-ok') },
                  disabled ? T('analysisUnavailable') : T('statusOk'))),
              item.method ? h('div', { className: 'sp-analysis-method' }, item.method) : null,
              disabled && item.reason ? h('div', { className: 'sp-analysis-reason' }, item.reason) : null)
          }))
      })
      var runner = null
      if (!dataset) {
        runner = h(Notice, { text: T('analysisNeedDataset') })
      } else if (!list.length) {
        runner = h(Notice, { text: catalog.loading ? T('loading') : T('analysisEmpty') })
      } else if (!meta) {
        runner = h(Notice, { text: T('analysisSelectHint') })
      } else {
        runner = h('div', { className: 'sp-analysis-runner' },
          h('div', { className: 'sp-analysis-head' },
            h('div', { className: 'sp-analysis-title' }, meta.title || meta.id),
            meta.method ? h('div', { className: 'sp-analysis-method' }, T('analysisMethod') + '：' + meta.method) : null),
          meta.needs && meta.needs.length
            ? h('div', { className: 'sp-inline-row' },
                h('span', { className: 'sp-toolbar-label' }, T('analysisNeeds')),
                meta.needs.map(function (need) {
                  var stepMeta = stepById(need)
                  return h('span', { key: 'need-' + need, className: 'sp-tag' }, stepMeta ? stepMeta.title : need)
                }))
            : null,
          meta.available === false && meta.reason ? h('div', { className: 'sp-warn-text' }, meta.reason) : null,
          h(ParamForm, {
            params: meta.params,
            values: values,
            columns: props.columns || [],
            disabled: run && run.busy === true,
            onChange: patch,
          }),
          h('div', { className: 'sp-inline-row' },
            Btn({
              small: true,
              variant: 'primary',
              icon: 'play',
              disabled: !dataset || (run && run.busy === true) || meta.available === false,
              onClick: function () { if (props.onRun) props.onRun(meta.id, values || {}) },
            }, run && run.busy ? T('analysisRunning') : T('analysisRun'))),
          run && run.error ? h(InlineError, { message: run.error }) : null,
          run && run.result
            ? h(ResultView, {
                result: run.result,
                onExport: props.tabExport && props.runExport
                  ? function () { props.runExport(props.tabExport.what) }
                  : null,
                exportBusy: props.tabExport ? !!(props.exports && props.exports[props.tabExport.what] && props.exports[props.tabExport.what].busy) : false,
                exportLabel: props.tabExport ? T(props.tabExport.labelKey) : '',
                exportFile: props.tabExport && props.exports && props.exports[props.tabExport.what]
                  ? props.exports[props.tabExport.what].file || ''
                  : '',
              })
            : null)
      }
      return h(Modal, {
        title: T('analysisTitle'),
        subtitle: T('analysisSubtitle'),
        wide: true,
        onClose: props.onClose,
      },
        h('div', { className: 'sp-modal-body' },
          catalog.loading ? h(LoadingLine, { text: T('loading') }) : null,
          catalog.error ? h(InlineError, { message: catalog.error }) : null,
          h('div', { className: 'sp-analysis-wrap' },
            h('div', { className: 'sp-analysis-list' }, groupNodes),
            runner)),
        h('div', { className: 'sp-modal-foot' },
          Btn({ icon: 'close', onClick: props.onClose, disabled: run && run.busy === true }, T('close'))))
    }

    // ==========================================================================
    // 整页工作台：状态容器 + 三栏布局 + 两个弹层 + 分析抽屉
    //   · 左栏状态由 STEPS + dataset.steps 派生（stepStatusOf）
    //   · 步骤/分析结果只保存在本会话（切数据集清空；后端 history/steps 仍保留）
    //   · 所有故障单点包 PluginErrorBoundary，不黑整页
    // ==========================================================================

    function pickDataset(datasets, id) {
      var list = Array.isArray(datasets) ? datasets : []
      for (var index = 0; index < list.length; index++) {
        if (list[index].id === id) return list[index]
      }
      return null
    }

    function datasetIdOf(payload) {
      if (!payload) return ''
      if (typeof payload === 'string') return plainText(payload)
      if (typeof payload !== 'object') return ''
      if (payload.dataset && typeof payload.dataset === 'object') return plainText(payload.dataset.id || '')
      return plainText(payload.id || payload.active || '')
    }

    /** 步骤卡上的结果摘要：优先 bullets，其次 metrics，最后图表计数。 */
    function resultSummaryText(result) {
      if (!result) return ''
      var bullets = Array.isArray(result.bullets) ? result.bullets : []
      if (bullets.length) return bullets.slice(0, 2).join('；')
      var metrics = result.metrics && typeof result.metrics === 'object' ? result.metrics : {}
      var keys = Object.keys(metrics)
      if (keys.length) {
        var parts = []
        for (var index = 0; index < keys.length && parts.length < 4; index++) {
          var value = metrics[keys[index]]
          parts.push(keys[index] + ' ' + (typeof value === 'number' ? fmtNumber(value) : formatCell(value, false)))
        }
        return parts.join(' · ')
      }
      var counts = []
      if (result.tables && result.tables.length) counts.push(result.tables.length + ' ' + T('resultColumn'))
      if (result.figures && result.figures.length) counts.push(result.figures.length + ' ' + T('plot'))
      return counts.join(' · ')
    }

    function WorkbenchBody() {
      var api = useMemo(function () { return createApi() }, [])
      var bootState = useState({ loading: true, error: '', health: normalizeHealth({}), formats: null })
      var boot = bootState[0]
      var setBoot = bootState[1]
      var bootNonceState = useState(0)
      var bootNonce = bootNonceState[0]
      var setBootNonce = bootNonceState[1]
      var stateState = useState(function () { return normalizeState({}) })
      var snapshot = stateState[0]
      var setSnapshot = stateState[1]
      var previewState = useState(function () { return normalizePreview({}) })
      var preview = previewState[0]
      var setPreview = previewState[1]
      var valuesState = useState({})
      var stepValues = valuesState[0]
      var setStepValues = valuesState[1]
      var resultsState = useState({})
      var stepResults = resultsState[0]
      var setStepResults = resultsState[1]
      var summariesState = useState({})
      var stepSummaries = summariesState[0]
      var setStepSummaries = summariesState[1]
      var errorsState = useState({})
      var stepErrors = errorsState[0]
      var setStepErrors = errorsState[1]
      var busyStepState = useState('')
      var stepBusy = busyStepState[0]
      var setStepBusy = busyStepState[1]
      var failedState = useState('')
      var failedStep = failedState[0]
      var setFailedStep = failedState[1]
      var plotsState = useState({})
      var plots = plotsState[0]
      var setPlots = plotsState[1]
      var exportsState = useState({})
      var exportsMap = exportsState[0]
      var setExports = exportsState[1]
      var catalogState = useState({ loading: false, error: '', data: normalizeAnalyses({}) })
      var catalog = catalogState[0]
      var setCatalog = catalogState[1]
      var runsState = useState({})
      var runs = runsState[0]
      var setRuns = runsState[1]
      var tabState = useState('overview')
      var tab = tabState[0]
      var setTab = tabState[1]
      var dialogState = useState({ geo: false, import: false, analysis: false })
      var dialogs = dialogState[0]
      var setDialogs = dialogState[1]
      var copiedState = useState('')
      var copiedKey = copiedState[0]
      var setCopiedKey = copiedState[1]
      var refreshingState = useState(false)
      var refreshing = refreshingState[0]
      var setRefreshing = refreshingState[1]
      var busyState = useState(false)
      var busy = busyState[0]
      var setBusy = busyState[1]
      var bannerState = useState({ kind: '', text: '' })
      var banner = bannerState[0]
      var setBanner = bannerState[1]
      var activeRef = useRef('')
      var catalogRef = useRef('')
      var copyTimerRef = useRef(null)

      // 首屏：health + state + formats 并发拉取；bootNonce 用于「重试」。
      useEffect(function () {
        var cancelled = false
        setBoot(function (previous) { return Object.assign({}, previous, { loading: true, error: '' }) })
        Promise.all([
          api.getHealth().catch(function () { return null }),
          api.getState(),
          api.getFormats().catch(function () { return null }),
        ]).then(function (list) {
          if (cancelled) return
          var next = normalizeState(list[1])
          var health = normalizeHealth(list[0] || {})
          if (!health.ready && next.engine.ready) health = next.engine
          setSnapshot(next)
          setBoot({ loading: false, error: '', health: health, formats: list[2] })
        }).catch(function (error) {
          if (cancelled) return
          setBoot({ loading: false, error: formatError(error), health: normalizeHealth({}), formats: null })
        })
        return function () { cancelled = true }
      }, [bootNonce, api])

      // 切换数据集：清掉本会话的图/结果/分析结果，重新拉 preview 与 analyses 目录。
      useEffect(function () {
        if (activeRef.current === snapshot.active) return undefined
        activeRef.current = snapshot.active
        catalogRef.current = snapshot.active
        setStepResults({})
        setStepSummaries({})
        setStepErrors({})
        setStepBusy('')
        setFailedStep('')
        setExports({})
        setPlots({})
        setStepValues({})
        setRuns({})
        setCopiedKey('')
        if (!snapshot.active) {
          setPreview(normalizePreview({}))
          setCatalog({ loading: false, error: '', data: normalizeAnalyses({}) })
          return undefined
        }
        var cancelled = false
        api.preview(snapshot.active)
          .then(function (payload) { if (!cancelled) setPreview(normalizePreview(payload)) })
          .catch(function () { if (!cancelled) setPreview(normalizePreview({})) })
        setCatalog(function (previous) { return Object.assign({}, previous, { loading: true, error: '' }) })
        api.analyses(snapshot.active)
          .then(function (payload) {
            if (cancelled) return
            setCatalog({ loading: false, error: '', data: normalizeAnalyses(payload) })
          })
          .catch(function (error) {
            if (cancelled) return
            setCatalog({ loading: false, error: formatError(error), data: normalizeAnalyses({}) })
          })
        return function () { cancelled = true }
      }, [snapshot.active, api])

      useEffect(function () {
        return function () {
          if (copyTimerRef.current) clearTimeout(copyTimerRef.current)
        }
      }, [])

      var dataset = pickDataset(snapshot.datasets, snapshot.active)
      var health = snapshot.engine && snapshot.engine.ready ? snapshot.engine : boot.health
      var columns = dataset && dataset.sample_columns && dataset.sample_columns.length
        ? dataset.sample_columns
        : (preview.columns && preview.columns.length ? preview.columns : [])

      var statuses = {}
      var paramValues = {}
      for (var stepsLoop = 0; stepsLoop < STEPS.length; stepsLoop++) {
        var stepMeta = STEPS[stepsLoop]
        statuses[stepMeta.id] = stepStatusOf(stepMeta, dataset, stepBusy, failedStep)
        paramValues[stepMeta.id] = Object.assign({}, defaultParams(stepMeta), stepValues[stepMeta.id] || {})
      }

      function refreshState(options) {
        var opts = options && typeof options === 'object' ? options : {}
        if (!opts.silent) setRefreshing(true)
        return api.getState().then(function (payload) {
          var next = normalizeState(payload)
          setSnapshot(next)
          setBoot(function (previous) {
            var merged = previous.health
            if (next.engine.ready && !merged.ready) merged = next.engine
            return Object.assign({}, previous, { health: merged })
          })
          if (!opts.silent) setRefreshing(false)
          return next
        }).catch(function (error) {
          if (!opts.silent) {
            setRefreshing(false)
            setBanner({ kind: 'bad', text: formatError(error) })
          }
          return null
        })
      }

      function refreshPreview(targetId) {
        var id = plainText(targetId) || snapshot.active
        if (!id) {
          setPreview(normalizePreview({}))
          return Promise.resolve(null)
        }
        return api.preview(id).then(function (payload) {
          var normalized = normalizePreview(payload)
          setPreview(normalized)
          return normalized
        }).catch(function () { return null })
      }

      function refreshCatalog() {
        var id = snapshot.active
        if (!id) return Promise.resolve(null)
        setCatalog(function (previous) { return Object.assign({}, previous, { loading: true, error: '' }) })
        return api.analyses(id).then(function (payload) {
          var data = normalizeAnalyses(payload)
          setCatalog({ loading: false, error: '', data: data })
          return data
        }).catch(function (error) {
          setCatalog({ loading: false, error: formatError(error), data: normalizeAnalyses({}) })
          return null
        })
      }

      function loadPlot(kind, params, force) {
        var id = snapshot.active
        if (!id) return
        var safeParams = params && typeof params === 'object' ? params : {}
        var key = plotKey(kind, safeParams)
        var entry = plots[key]
        if (!force && entry && (entry.loading || entry.figure)) return
        setPlots(function (previous) {
          var next = Object.assign({}, previous)
          next[key] = { loading: true, error: '', figure: null }
          return next
        })
        api.plot(id, kind, safeParams).then(function (payload) {
          var figure = normalizeFigure(payload && payload.figure ? payload.figure : payload)
          setPlots(function (previous) {
            var next = Object.assign({}, previous)
            next[key] = { loading: false, error: figure ? '' : T('noFigure'), figure: figure }
            return next
          })
        }).catch(function (error) {
          setPlots(function (previous) {
            var next = Object.assign({}, previous)
            next[key] = { loading: false, error: formatError(error), figure: null }
            return next
          })
        })
      }

      function plotOnce(kind, params) {
        var id = snapshot.active
        if (!id) return Promise.reject(new Error(T('analysisNeedDataset')))
        return api.plot(id, kind, params && typeof params === 'object' ? params : {}).then(function (payload) {
          var figure = normalizeFigure(payload && payload.figure ? payload.figure : payload)
          if (!figure) throw new Error(T('noFigure'))
          return figure
        })
      }

      function runExport(what) {
        var id = snapshot.active
        if (!id) return
        setExports(function (previous) {
          var next = Object.assign({}, previous)
          next[what] = Object.assign({}, previous[what] || {}, { busy: true, error: '' })
          return next
        })
        api.exportData(id, what).then(function (payload) {
          var file = plainText(payload && (payload.file || payload.path)) ||
            plainText(payload && payload.result && payload.result.file)
          setExports(function (previous) {
            var next = Object.assign({}, previous)
            next[what] = { busy: false, file: file, error: '' }
            return next
          })
          refreshState({ silent: true })
        }).catch(function (error) {
          setExports(function (previous) {
            var next = Object.assign({}, previous)
            next[what] = { busy: false, file: '', error: formatError(error) }
            return next
          })
        })
      }

      function runAnalysis(analysisId, params) {
        var id = snapshot.active
        if (!id) return
        setRuns(function (previous) {
          var next = Object.assign({}, previous)
          var previousRun = previous[analysisId] || {}
          next[analysisId] = { busy: true, error: '', result: previousRun.result || null }
          return next
        })
        api.analysis(id, analysisId, params && typeof params === 'object' ? params : {}).then(function (payload) {
          var result = normalizeResult(payload)
          setRuns(function (previous) {
            var next = Object.assign({}, previous)
            next[analysisId] = { busy: false, error: '', result: result }
            return next
          })
          refreshState({ silent: true })
        }).catch(function (error) {
          setRuns(function (previous) {
            var next = Object.assign({}, previous)
            next[analysisId] = { busy: false, error: formatError(error), result: null }
            return next
          })
        })
      }

      function handleRunStep(stepId, params) {
        var step = stepById(stepId)
        var id = snapshot.active
        if (!step || !id) return Promise.resolve(null)
        var requestParams = params && typeof params === 'object' && !Array.isArray(params)
          ? params
          : collectStepParams(step, stepValues[stepId] || {})
        setStepBusy(stepId)
        setFailedStep('')
        setBanner({ kind: '', text: '' })
        setStepErrors(function (previous) {
          var next = Object.assign({}, previous)
          delete next[stepId]
          return next
        })
        return api.runStep(id, stepId, requestParams).then(function (payload) {
          var result = normalizeResult(payload)
          setStepResults(function (previous) {
            var next = Object.assign({}, previous)
            next[stepId] = result
            return next
          })
          setStepSummaries(function (previous) {
            var next = Object.assign({}, previous)
            next[stepId] = resultSummaryText(result)
            return next
          })
          setStepValues(function (previous) {
            var next = Object.assign({}, previous)
            next[stepId] = Object.assign({}, requestParams)
            return next
          })
          setStepBusy('')
          setPlots({})
          setBanner({ kind: '', text: tf(T, 'runStepDone', step.title) })
          refreshCatalog()
          if (TAB_FOR_STEP[stepId]) setTab(TAB_FOR_STEP[stepId])
          refreshState({ silent: true }).then(function () {
            refreshPreview()
          })
          return result
        }).catch(function (error) {
          setStepBusy('')
          setFailedStep(stepId)
          setStepErrors(function (previous) {
            var next = Object.assign({}, previous)
            next[stepId] = formatError(error)
            return next
          })
          return null
        })
      }

      function handleParamChange(stepId, patch) {
        var step = stepById(stepId)
        setStepValues(function (previous) {
          var next = Object.assign({}, previous)
          next[stepId] = Object.assign({}, step ? defaultParams(step) : {}, previous[stepId] || {}, patch || {})
          return next
        })
      }

      function handleSelectDataset(id) {
        if (!id || id === snapshot.active || busy) return
        setBusy(true)
        setBanner({ kind: '', text: '' })
        api.selectDataset(id).then(function () {
          return refreshState({ silent: true })
        }).then(function () {
          setBusy(false)
        }).catch(function (error) {
          setBusy(false)
          setBanner({ kind: 'bad', text: formatError(error) })
        })
      }

      function handleCloseDataset(id) {
        if (!id || busy) return
        setBusy(true)
        api.closeDataset(id).then(function () {
          return refreshState({ silent: true })
        }).then(function () {
          setBusy(false)
        }).catch(function (error) {
          setBusy(false)
          setBanner({ kind: 'bad', text: formatError(error) })
        })
      }

      function handleImported(payload) {
        var targetId = datasetIdOf(payload)
        setDialogs({ geo: false, import: false, analysis: false })
        setTab('overview')
        setBusy(true)
        refreshState({ silent: true }).then(function (next) {
          var active = next ? next.active : ''
          if (targetId && active && targetId !== active) {
            return api.selectDataset(targetId).then(function () { return refreshState({ silent: true }) })
          }
          return null
        }).then(function () {
          setBusy(false)
          refreshPreview(targetId)
          refreshCatalog()
        }).catch(function (error) {
          setBusy(false)
          setBanner({ kind: 'bad', text: formatError(error) })
        })
      }

      function loadDemo() {
        if (busy) return
        setBusy(true)
        setBanner({ kind: '', text: '' })
        api.demo()
          .then(function (payload) { handleImported(payload) })
          .catch(function (error) {
            setBusy(false)
            setBanner({ kind: 'bad', text: formatError(error) })
          })
      }

      function copyText(text, key) {
        var value = plainText(text, 4000)
        if (!value) return
        var markCopied = function () {
          setCopiedKey(key || 'copy')
          if (copyTimerRef.current) clearTimeout(copyTimerRef.current)
          copyTimerRef.current = setTimeout(function () { setCopiedKey('') }, 1800)
        }
        var fallback = function () {
          try {
            if (typeof document === 'undefined' || !document.body) return false
            var area = document.createElement('textarea')
            area.value = value
            area.setAttribute('readonly', 'readonly')
            area.style.position = 'fixed'
            area.style.top = '-1000px'
            area.style.opacity = '0'
            document.body.appendChild(area)
            area.select()
            var ok = typeof document.execCommand === 'function' ? document.execCommand('copy') : false
            document.body.removeChild(area)
            return ok === true
          } catch (ignored) {
            return false
          }
        }
        var clipboard = null
        try { clipboard = typeof navigator === 'undefined' ? null : navigator.clipboard } catch (ignored) { clipboard = null }
        if (clipboard && typeof clipboard.writeText === 'function') {
          Promise.resolve(clipboard.writeText(value)).then(markCopied).catch(function () {
            if (fallback()) markCopied()
            else setBanner({ kind: 'bad', text: T('copyFailed') })
          })
          return
        }
        if (fallback()) markCopied()
        else setBanner({ kind: 'bad', text: T('copyFailed') })
      }

      function openDialog(name) {
        setDialogs(function (previous) {
          var next = Object.assign({}, previous)
          next[name] = true
          return next
        })
      }

      function closeDialog(name) {
        setDialogs(function (previous) {
          var next = Object.assign({}, previous)
          next[name] = false
          return next
        })
      }

      var tabProps = {
        dataset: dataset,
        preview: preview,
        columns: columns,
        api: api,
        catalog: catalog,
        runs: runs,
        plots: plots,
        loadPlot: loadPlot,
        stepResults: stepResults,
        stepBusy: stepBusy,
        exports: exportsMap,
        runExport: runExport,
        runAnalysis: runAnalysis,
        onRunAnalysis: runAnalysis,
        plotOnce: plotOnce,
        runStep: handleRunStep,
        onRunStep: handleRunStep,
        onAttached: function () { refreshState({ silent: true }); refreshCatalog() },
        onOpenAnalysis: function () { openDialog('analysis') },
        busy: busy,
      }

      function renderTab() {
        if (tab === 'spatial') return h(SpatialMapTab, tabProps)
        if (tab === 'expr') return h(ExprTab, tabProps)
        if (tab === 'de') return h(DeTab, tabProps)
        if (tab === 'stats') return h(SvgTab, tabProps)
        if (tab === 'domains') return h(DomainTab, tabProps)
        if (tab === 'deconv') return h(DeconvTab, tabProps)
        if (tab === 'comm') return h(CommTab, tabProps)
        if (tab === 'records') return h(RecordsTab, tabProps)
        if (tab === 'log') return h(LogTab, { api: api })
        return h(OverviewTab, Object.assign({}, tabProps, {
          health: health,
          suggestions: snapshot.suggestions,
          busyStep: stepBusy ? (stepById(stepBusy) || { title: stepBusy }).title : '',
          copiedKey: copiedKey,
          onCopy: copyText,
          onOpenImport: function () { openDialog('import') },
        }))
      }

      var bannerNode = banner.text
        ? h('div', { className: 'sp-banner' + (banner.kind === 'bad' ? ' sp-banner-bad' : '') },
            h('span', { className: 'sp-banner-text', title: banner.text }, banner.text),
            Btn({
              small: true,
              variant: 'ghost',
              icon: 'close',
              onClick: function () { setBanner({ kind: '', text: '' }) },
            }, T('close')))
        : null

      var engineNode = !boot.loading && !boot.error && health && !health.ready
        ? h('div', { className: 'sp-banner' },
            h('span', { className: 'sp-banner-text', title: T('engineGuideBody') + ' ' + T('engineGuideSteps') },
              T('engineGuideTitle') + ' · ' + (health.error || T('engineBooting'))))
        : null

      var centerNode = boot.loading
        ? h('div', { className: 'sp-pane' }, h(LoadingLine, { text: T('loadingData') }))
        : (boot.error
            ? h('div', { className: 'sp-pane' },
                h(ErrorBanner, {
                  title: T('errorTitle'),
                  message: boot.error,
                  retryLabel: T('retry'),
                  onRetry: function () { setBootNonce(function (previous) { return previous + 1 }) },
                }))
            : h(PluginErrorBoundary, { key: 'pane-' + tab }, renderTab()))

      return h('div', { className: 'sp-root' },
        h('style', { 'data-plugin': NS }, CSS),
        h(PluginErrorBoundary, null,
          h(TopBar, {
            health: health,
            datasets: snapshot.datasets,
            activeId: snapshot.active,
            busy: busy || refreshing,
            refreshing: refreshing,
            onSelect: handleSelectDataset,
            onCloseDataset: handleCloseDataset,
            onGeo: function () { openDialog('geo') },
            onImport: function () { openDialog('import') },
            onDemo: loadDemo,
            onRefresh: function () {
              refreshState({ silent: false }).then(function (next) {
                if (next) { refreshPreview(); refreshCatalog() }
              })
            },
          })),
        bannerNode,
        h('div', { className: 'sp-body' },
          h(PluginErrorBoundary, null,
            h(PipelinePanel, {
              dataset: dataset,
              statuses: statuses,
              values: paramValues,
              columns: columns,
              summaries: stepSummaries,
              errors: stepErrors,
              busyStep: stepBusy,
              onParamChange: handleParamChange,
              onRun: handleRunStep,
            })),
          h('div', { className: 'sp-center' },
            engineNode,
            h(PluginErrorBoundary, null, h(TabsBar, { active: tab, onChange: setTab })),
            centerNode),
          h(PluginErrorBoundary, null,
            h(AgentPanel, {
              dataset: dataset,
              health: health,
              suggestions: snapshot.suggestions,
              activity: snapshot.activity,
              copiedKey: copiedKey,
              onCopy: copyText,
            }))),
        dialogs.geo
          ? h(PluginErrorBoundary, { key: 'dialog-geo' },
              h(GeoDialog, {
                api: api,
                onImported: handleImported,
                onClose: function () { closeDialog('geo') },
              }))
          : null,
        dialogs.import
          ? h(PluginErrorBoundary, { key: 'dialog-import' },
              h(ImportDialog, {
                api: api,
                formats: boot.formats,
                onImported: handleImported,
                onClose: function () { closeDialog('import') },
              }))
          : null,
        dialogs.analysis
          ? h(PluginErrorBoundary, { key: 'dialog-analysis' },
              h(AnalysisDrawer, {
                dataset: dataset,
                columns: columns,
                catalog: catalog,
                runs: runs,
                onRun: runAnalysis,
                onClose: function () { closeDialog('analysis') },
              }))
          : null)
    }

    /** 整页入口：整体再包一层错误边界（单点异常不黑整页）。 */
    function WorkbenchPage() {
      return h(PluginErrorBoundary, null, h(WorkbenchBody))
    }

    // ==========================================================================
    // 插件入口：locale 字典 + 两个槽（整页 main + 侧边栏入口 sidebar.panellist）
    // ==========================================================================

    function apply(ctx) {
      var locale = ctx && ctx.locale
      if (locale && typeof locale.register === 'function') {
        var registerDictionaries = function () {
          try {
            return locale.register(NS, { zh: DICT_ZH, en: DICT_EN })
          } catch (error) {
            warnPlugin('failed to register locale dictionaries', error)
            return undefined
          }
        }
        try {
          if (typeof ctx.effect === 'function') ctx.effect(registerDictionaries, NS + ': locale')
          else registerDictionaries()
        } catch (error) {
          warnPlugin('failed to install locale effect', error)
        }
        T = createTranslator(locale)
      }

      var slots = ctx && ctx.slots
      if (!slots && ctx && typeof ctx.get === 'function') {
        try { slots = ctx.get('slots') } catch (error) { slots = null }
      }
      if (!slots || typeof slots.register !== 'function' || typeof slots.inject !== 'function') {
        warnPlugin('slots service unavailable; UI not registered')
        return
      }

      var disposers = []
      function keep(handle) {
        if (typeof handle === 'function') disposers.push(handle)
        return handle
      }

      try {
        if (typeof ctx.effect === 'function') {
          ctx.effect(function () {
            return function () {
              while (disposers.length) {
                var dispose = disposers.pop()
                try { dispose() } catch (ignored) { /* 单个 disposer 失败不影响其它 */ }
              }
            }
          }, NS + ': slot disposers')
        }
      } catch (error) {
        warnPlugin('failed to install slot cleanup', error)
      }

      try {
        keep(slots.inject('main', function () {
          return keep(slots.register({ name: 'main', key: 'spatial', locale: NS }, WorkbenchPage))
        }))
      } catch (error) {
        warnPlugin('failed to register main slot', error)
      }

      try {
        keep(slots.inject('sidebar.panellist', function () {
          return keep(slots.register({
            name: 'sidebar.panellist',
            id: 'spatial',
            order: 5,
            label: function () { return T('nav') },
            locale: NS,
          }, SidebarIcon))
        }))
      } catch (error) {
        warnPlugin('failed to register sidebar.panellist slot', error)
      }
    }

    // __test__ 只给自测/联调用（宿主不消费）：纯函数 + 组件引用。
    return {
      name: NS,
      inject: ['slots','locale'],
      apply: apply,
      __test__: Object.assign({}, SP_TEST, {
        components: {
          WorkbenchPage: WorkbenchPage,
          WorkbenchBody: WorkbenchBody,
          PipelinePanel: PipelinePanel,
          StepCard: StepCard,
          DatasetCard: DatasetCard,
          TopBar: TopBar,
          AgentPanel: AgentPanel,
          TabsBar: TabsBar,
          OverviewTab: OverviewTab,
          SpatialMapTab: SpatialMapTab,
          ExprTab: ExprTab,
          DeTab: DeTab,
          SvgTab: SvgTab,
          DomainTab: DomainTab,
          DeconvTab: DeconvTab,
          CommTab: CommTab,
          AnalysisTab: AnalysisTab,
          RefPanel: RefPanel,
          RecordsTab: RecordsTab,
          LogTab: LogTab,
          GeoDialog: GeoDialog,
          ImportDialog: ImportDialog,
          AnalysisDrawer: AnalysisDrawer,
          ResultView: ResultView,
          DataTable: DataTable,
          ParamForm: ParamForm,
          ParamControl: ParamControl,
          FigureBlock: FigureBlock,
          PlotBlock: PlotBlock,
          MiniMap: MiniMap,
          Modal: Modal,
          MetricCard: MetricCard,
          SectionCard: SectionCard,
          FormField: FormField,
          GuideBlock: GuideBlock,
          GenePicker: GenePicker,
          InlineError: InlineError,
          Notice: Notice,
          Spinner: Spinner,
          LoadingLine: LoadingLine,
          ErrorBanner: ErrorBanner,
          Icon: Icon,
          Btn: Btn,
          PluginErrorBoundary: PluginErrorBoundary,
          SidebarIcon: SidebarIcon,
          SidebarGlyph: SidebarGlyph,
        },
      }),
    }
  }

  // ==========================================================================
  // 顶层守卫：node 里只暴露纯函数给单测；浏览器里注册 Client 模块。
  // ==========================================================================

  if (typeof window === 'undefined') {
    globalThis.__SP_TEST__ = SP_TEST
  } else if (window.__ModuleLoader__ && typeof window.__ModuleLoader__.load === 'function') {
    window.__ModuleLoader__.load({ id: NS, factory: factory })
  }
})()
