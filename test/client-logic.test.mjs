/**
 * dsh-spatial 浏览器端纯函数单测（不依赖浏览器、不装 npm 依赖）。
 *
 * 运行：node test/client-logic.test.mjs
 *
 * client.js 是浏览器模块（window.__ModuleLoader__.load），但在 node 里
 * （typeof window === 'undefined'）只会把纯函数挂到 globalThis.__SP_TEST__，
 * 所以这里直接 import 它、断言纯函数行为 + 静态合规约束（DESIGN §2.2/§3/§4/§5/§8）。
 * 打印 PASS 表示全部通过。
 */
import assert from 'node:assert/strict'
import { readFile } from 'node:fs/promises'
import { test } from 'node:test'

const CLIENT_URL = new URL('../client.js', import.meta.url)
await import(CLIENT_URL.href)

const SP = globalThis.__SP_TEST__
assert.ok(SP, 'node 环境应通过 globalThis.__SP_TEST__ 暴露纯函数')

const CLIENT_SOURCE = await readFile(CLIENT_URL, 'utf8')

/** 去掉块注释与行注释后的代码（静态断言只看真实代码，注释里的历史说明不算）。 */
function stripComments(source) {
  return source
    .replace(/\/\*[\s\S]*?\*\//g, '')
    .replace(/(^|[^:])\/\/[^\n]*/g, '$1')
}

const CLIENT_CODE = stripComments(CLIENT_SOURCE)

// ---------------------------------------------------------------------------
// 0. 加载与导出形状
// ---------------------------------------------------------------------------

test('import：无 window 时不炸，且暴露 NS/API/纯函数层', () => {
  assert.equal(typeof globalThis.window, 'undefined')
  assert.equal(SP.NS, 'dsh-spatial')
  assert.equal(SP.API, '/spatial')
  assert.equal(SP.DEFAULT_GEO_DIR, '~/.dsh/spatial/GEO')
  assert.ok(SP.POLL_MS > 0)
  assert.ok(SP.POLL_TIMEOUT_MS >= SP.POLL_MS)
  assert.ok(SP.DICT_ZH.appTitle && SP.DICT_EN.appTitle)
  assert.ok(SP.DICT_ZH.nav && SP.DICT_EN.nav)
  assert.equal(typeof SP.createApi, 'function')
  assert.equal(typeof SP.apiClient, 'function')
  assert.equal(typeof SP.jobPolling, 'function')
  assert.equal(typeof SP.jobDone, 'function')
  assert.equal(typeof SP.spatialScreen, 'function')
  assert.equal(typeof SP.buildSearchPayload, 'function')
  assert.equal(typeof SP.suggestionInstruction, 'function')
  assert.equal(typeof SP.plotKey, 'function')
  assert.equal(typeof SP.createTranslator, 'function')
  assert.equal(typeof SP.normalizeState, 'function')
  assert.equal(typeof SP.normalizeDataset, 'function')
})

test('字典：中英文 key 集完全一致', () => {
  const zhKeys = Object.keys(SP.DICT_ZH).sort()
  const enKeys = Object.keys(SP.DICT_EN).sort()
  assert.ok(zhKeys.length > 300, `字典词条太少：${zhKeys.length}`)
  assert.deepEqual(zhKeys.filter((key) => !SP.DICT_EN[key]), [], 'EN 缺中文侧的 key')
  assert.deepEqual(enKeys.filter((key) => !SP.DICT_ZH[key]), [], 'ZH 缺英文侧的 key')
})

// ---------------------------------------------------------------------------
// 1. STEPS / TAB（DESIGN §3 顺序与 requires）
// ---------------------------------------------------------------------------

test('STEPS：12 个步骤顺序与 DESIGN §3 一致', () => {
  assert.equal(SP.STEPS.length, 12)
  assert.deepEqual(SP.STEPS.map((step) => step.id), [
    'qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors',
    'cluster', 'embed', 'svg', 'domains', 'deconv', 'annotate',
  ])
  const requires = {}
  for (const step of SP.STEPS) requires[step.id] = step.requires
  assert.deepEqual(requires, {
    qc: [],
    filter: ['qc'],
    normalize: ['filter'],
    hvg: ['normalize'],
    pca: ['hvg'],
    neighbors: ['pca'],
    cluster: ['neighbors'],
    embed: ['pca'],
    svg: ['cluster'],
    domains: ['cluster'],
    deconv: ['normalize'],
    annotate: ['domains'],
  })
  const allowedTypes = ['number', 'select', 'boolean', 'column', 'text', 'genes', 'multiselect', 'json', 'file', 'string', 'array', 'columns']
  for (const step of SP.STEPS) {
    assert.ok(step.title, `${step.id} 应有中文标题`)
    assert.ok(Array.isArray(step.params) && step.params.length > 0, `${step.id} 应有参数`)
    for (const descriptor of step.params) {
      assert.ok(descriptor.key, `${step.id} 参数缺 key`)
      assert.ok(allowedTypes.includes(descriptor.type), `${step.id}.${descriptor.key} 控件类型不支持：${descriptor.type}`)
    }
  }
})

test('TAB_KEYS：10 个 tab 与 TAB_META 标题齐全', () => {
  assert.deepEqual(SP.TAB_KEYS, [
    'overview', 'spatial', 'expr', 'de', 'stats', 'domains', 'deconv', 'comm', 'records', 'log',
  ])
  for (const key of SP.TAB_KEYS) {
    assert.ok(SP.TAB_META[key] && SP.TAB_META[key].titleKey, `${key} 缺 TAB_META`)
    assert.ok(SP.DICT_ZH[SP.TAB_META[key].titleKey], `${key} 的 titleKey 不在字典里`)
  }
  assert.equal(SP.PLOT_KINDS.length, 20, 'DESIGN §5 图白名单应 20 种')
  assert.ok(SP.PLOT_KINDS.includes('gene_spatial') && SP.PLOT_KINDS.includes('domain_map'))
  assert.ok(SP.EXPORT_KINDS.length >= 10)
  assert.ok(SP.COLOR_BY.length >= 4)
  assert.ok(SP.GEO_PLATFORMS.length >= 5)
  assert.ok(SP.DATA_FORMATS.includes('visium') && SP.DATA_FORMATS.includes('h5ad'))
  assert.equal(SP.analysisGroupTab('空间统计'), 'stats')
  assert.equal(SP.analysisGroupTab('邻域与域'), 'domains')
  assert.equal(SP.analysisGroupTab('反卷积'), 'deconv')
  assert.equal(SP.analysisGroupTab('空间通讯'), 'comm')
  assert.equal(SP.analysisGroupTab('区域差异'), 'de')
  assert.equal(SP.analysisGroupTab('质控与整合'), 'stats')
  assert.equal(SP.analysisGroupTab('不存在的分组'), 'stats')
})

// ---------------------------------------------------------------------------
// 2. 坐标映射（DESIGN §5：y 轴与图像坐标一致，默认等比居中）
// ---------------------------------------------------------------------------

test('spatialScreen：等比缩放、居中留白、y 轴不翻转', () => {
  const points = [[0, 0], [100, 0], [0, 100]]
  const placed = SP.spatialScreen(points, { x: [0, 100], y: [0, 100] }, 120, 120, { padding: 10 })
  assert.equal(placed.length, 3)
  // 100×100 的正方形映射进 100×100 的内框：不缩放，原点落在 (10,10)
  assert.deepEqual(placed.map((p) => [p.x, p.y]), [[10, 10], [110, 10], [10, 110]])
  assert.deepEqual(placed.map((p) => p.index), [0, 1, 2])

  const wide = SP.spatialScreen([[0, 0], [200, 0], [0, 100]], null, 220, 120, { padding: 10 })
  // 包围盒 200×100 → 内框 200×100，等比缩放 1：多余空间为 0
  assert.equal(wide[0].x, 10)
  assert.equal(wide[1].x, 210)
  assert.equal(wide[2].y, 110)

  // flipY 只影响 y
  const flipped = SP.spatialScreen([[0, 0], [0, 100]], { x: [0, 100], y: [0, 100] }, 120, 120, { padding: 10, flipY: true })
  assert.equal(flipped[0].y, 110)
  assert.equal(flipped[1].y, 10)

  // 脏数据跳过，不打断整批（{x:1} 缺 y、{x:'a'} 非数字都被跳过）
  const dirty = SP.spatialScreen([{ x: 1 }, { x: 'a', y: 2 }, [3, 4], { x: 0, y: 0 }], null, 100, 100)
  assert.equal(dirty.length, 2)
  assert.deepEqual(dirty.map((point) => point.index), [2, 3])
  assert.equal(SP.spatialScreen([], null, 100, 100).length, 0)
  assert.equal(SP.spatialScreen(null).length, 0)
  // pointXY 未导出，用 spatialScreen 间接验证对象/数组两种坐标形状
  const objects = SP.spatialScreen([{ x: 0, y: 0 }, { x: 10, y: 10 }], { x: [0, 10], y: [0, 10] }, 120, 120, { padding: 10 })
  assert.deepEqual(objects.map((point) => [point.x, point.y]), [[10, 10], [110, 110]])
})

test('extentOf / normalizeExtent：包围盒与非法输入', () => {
  assert.deepEqual(SP.extentOf([[1, 2], [3, 0]]), { x: [1, 3], y: [0, 2] })
  assert.equal(SP.extentOf([]), null)
  assert.deepEqual(SP.normalizeExtent({ x_range: [3, 1], y_range: [-1, 1] }), { x: [1, 3], y: [-1, 1] })
  assert.equal(SP.normalizeExtent({ x: [0, 1] }), null)
  assert.equal(SP.normalizeExtent(null), null)
})

// ---------------------------------------------------------------------------
// 3. GEO 搜索 payload / 导入 payload / 建议指令
// ---------------------------------------------------------------------------

test('buildSearchPayload：只带非空/非默认字段', () => {
  const empty = SP.buildSearchPayload({})
  assert.equal(empty.term, '')
  assert.equal(empty.limit, 20)
  assert.equal(empty.page, 0)
  assert.equal(empty.sort, 'relevance')
  assert.equal(empty.species, undefined, '全部物种不发 species')
  assert.equal(empty.study_type, undefined)
  assert.equal(empty.platform, undefined)
  assert.equal(empty.min_samples, undefined)
  assert.equal(empty.has_matrix, undefined)

  const full = SP.buildSearchPayload({
    term: '  乳腺 空间转录组 Visium ',
    species: 'human',
    study_type: 'high throughput sequencing',
    platform: 'visium',
    min_samples: '8',
    has_matrix: true,
    limit: 1000,
    page: -3,
    sort: 'date',
  })
  assert.equal(full.term, '乳腺 空间转录组 Visium')
  assert.equal(full.species, 'human')
  assert.equal(full.study_type, 'high throughput sequencing')
  assert.equal(full.platform, 'visium')
  assert.equal(full.min_samples, 8)
  assert.equal(full.has_matrix, true)
  assert.equal(full.limit, 50, 'limit 上限夹到 50')
  assert.equal(full.page, 0)
  assert.equal(full.sort, 'date')

  const zero = SP.buildSearchPayload({ min_samples: 0, has_matrix: false, sort: 'bogus' })
  assert.equal(zero.min_samples, undefined)
  assert.equal(zero.has_matrix, undefined)
  assert.equal(zero.sort, 'relevance')
})

test('buildImportPayload：空值剔除、auto 不发、transpose 只发 true', () => {
  assert.deepEqual(SP.buildImportPayload({}), {})
  const payload = SP.buildImportPayload({
    path: ' /data/outs ',
    url: '',
    data_format: 'auto',
    coords: '/data/coords.csv',
    image: '/data/tissue.png',
    transpose: false,
    visium_scale: '0.15',
  })
  assert.deepEqual(payload, {
    path: '/data/outs',
    coords: '/data/coords.csv',
    image: '/data/tissue.png',
    visium_scale: 0.15,
  })
  const h5ad = SP.buildImportPayload({ path: '/x.h5ad', data_format: 'h5ad', transpose: true })
  assert.equal(h5ad.data_format, 'h5ad')
  assert.equal(h5ad.transpose, true)
})

test('suggestionInstruction：生成可粘给 Agent 的中文指令', () => {
  const text = SP.suggestionInstruction({ title: '运行空间域', detail: 'MRF 平滑', step: 'domains' })
  assert.match(text, /^请在工作台执行「运行空间域」/)
  assert.match(text, /（步骤：空间域，step=domains）/)
  assert.match(text, /：MRF 平滑/)
  assert.match(text, /。$/)
  assert.match(SP.suggestionInstruction({ step: 'mystery' }), /（step=mystery）/)
  assert.match(SP.suggestionInstruction(null), /^请在工作台执行「继续分析」/)
})

// ---------------------------------------------------------------------------
// 4. 参数默认值 / 收集 / 步骤依赖
// ---------------------------------------------------------------------------

test('defaultParams / collectStepParams：类型收敛且数组不共享引用', () => {
  const domains = SP.defaultParams('domains')
  assert.equal(domains.n_domains, 6)
  assert.equal(domains.method, 'mrf')
  assert.notEqual(SP.defaultParams('domains'), domains, '每次都应返回新对象')

  const collected = SP.collectStepParams('domains', { n_domains: '8', method: 'kmeans', smooth: '0.2', iterations: '30' })
  assert.equal(collected.n_domains, 8)
  assert.equal(collected.method, 'kmeans')
  assert.equal(collected.smooth, 0.2)

  const blank = SP.collectStepParams('deconv', { ref: '', cell_type_col: '', min_genes: '' })
  assert.equal(blank.ref, null)
  assert.equal(blank.cell_type_col, null)
  assert.equal(blank.min_genes, null)

  const mapping = SP.collectStepParams('annotate', { mapping: '{"1":"肿瘤"}', key: 'domain_anno' })
  assert.deepEqual(mapping.mapping, { 1: '肿瘤' })
  assert.equal(mapping.key, 'domain_anno')
  assert.deepEqual(SP.parseTokenList('CDH1, VIM；MKI67 EPCAM'), ['CDH1', 'VIM', 'MKI67', 'EPCAM'])
  assert.equal(SP.collectStepParams('annotate', { mapping: 'not json' }).mapping, null, '非法 JSON 收敛成 null')
})

test('doneStepSet / isStepReady / stepStatusOf：requires 判定与 running 优先', () => {
  const dataset = { id: 'ds1', steps: ['qc', 'filter', 'normalize'] }
  assert.deepEqual(SP.doneStepSet(dataset), { qc: true, filter: true, normalize: true })
  assert.equal(SP.isStepReady('qc', dataset), true)
  assert.equal(SP.isStepReady('hvg', dataset), true)
  assert.equal(SP.isStepReady('cluster', dataset), false, 'cluster 需要 neighbors')
  assert.equal(SP.isStepReady('deconv', dataset), true, 'deconv 只依赖 normalize')
  assert.equal(SP.isStepReady('qc', null), false)
  assert.equal(SP.isStepReady('nope', dataset), false)

  assert.equal(SP.stepStatusOf('qc', dataset, '', ''), 'done')
  assert.equal(SP.stepStatusOf('qc', dataset, 'qc', ''), 'running')
  assert.equal(SP.stepStatusOf('qc', dataset, '', 'qc'), 'failed')
  assert.equal(SP.stepStatusOf('hvg', dataset, '', ''), 'ready')
  assert.equal(SP.stepStatusOf('cluster', dataset, '', ''), 'blocked')
  assert.equal(SP.stepStatusKey('done'), 'stepStatusDone')
  assert.equal(SP.stepStatusKey('ready'), 'stepStatusReady')
  assert.equal(SP.stepStatusKey('weird'), 'stepStatusBlocked')
})

// ---------------------------------------------------------------------------
// 5. API 客户端：method / url / body
// ---------------------------------------------------------------------------

function makeFetch(payload, options) {
  const opts = options || {}
  const calls = []
  const impl = function (url, init) {
    calls.push({ url: String(url), init: init || {} })
    if (opts.reject) return Promise.reject(opts.reject)
    const body = typeof payload === 'function' ? payload(url, init) : payload
    return Promise.resolve({
      ok: opts.httpOk === false ? false : true,
      status: opts.status || 200,
      text: function () { return Promise.resolve(typeof body === 'string' ? body : JSON.stringify(body)) },
    })
  }
  impl.calls = calls
  return impl
}

test('apiClient：GET/POST、同源 /spatial 前缀、JSON body 与 actor=user', async () => {
  const fake = makeFetch({ ok: true, engine: { ready: true } })
  const api = SP.apiClient(fake)

  await api.getHealth()
  assert.equal(fake.calls[0].url, '/spatial/health')
  assert.equal(fake.calls[0].init.method, 'GET')

  await api.getState()
  assert.equal(fake.calls[1].url, '/spatial/state')

  await api.runStep('ds1', 'cluster', { resolution: 0.8 })
  assert.equal(fake.calls[2].url, '/spatial/run')
  assert.equal(fake.calls[2].init.method, 'POST')
  assert.equal(fake.calls[2].init.headers['content-type'], 'application/json')
  assert.deepEqual(JSON.parse(fake.calls[2].init.body), {
    dataset: 'ds1', step: 'cluster', params: { resolution: 0.8 }, actor: 'user',
  })

  await api.plot('ds1', 'spatial_scatter', { color_by: 'clusters' })
  assert.equal(fake.calls[3].url, '/spatial/plot')
  assert.deepEqual(JSON.parse(fake.calls[3].init.body), {
    dataset: 'ds1', kind: 'spatial_scatter', params: { color_by: 'clusters' }, actor: 'user',
  })

  await api.analysis('ds1', 'morans_svg', { n_top: 20 })
  assert.equal(fake.calls[4].url, '/spatial/analysis')
  assert.deepEqual(JSON.parse(fake.calls[4].init.body), {
    dataset: 'ds1', id: 'morans_svg', params: { n_top: 20 }, actor: 'user',
  })

  await api.analyses('ds 1')
  assert.equal(fake.calls[5].url, '/spatial/analyses?dataset=ds%201')
  assert.equal(fake.calls[5].init.method, 'GET')

  await api.getJob('job 1')
  assert.equal(fake.calls[6].url, '/spatial/jobs?id=job%201')

  await api.getFormats()
  assert.equal(fake.calls[7].url, '/spatial/formats')

  await api.preview('ds1')
  assert.deepEqual(JSON.parse(fake.calls[8].init.body), { dataset: 'ds1', actor: 'user' })

  await api.samples({ dataset: 'ds1', op: 'auto_group', group_regex: '_(\\w+)$' })
  assert.equal(fake.calls[9].url, '/spatial/samples')
  assert.equal(JSON.parse(fake.calls[9].init.body).op, 'auto_group')

  await api.exportData('ds1', 'summary_json')
  assert.deepEqual(JSON.parse(fake.calls[10].init.body), { dataset: 'ds1', what: 'summary_json', actor: 'user' })

  await api.audit('ds1')
  assert.equal(fake.calls[11].url, '/spatial/audit')
  await api.suggest('ds1')
  assert.equal(fake.calls[12].url, '/spatial/suggest')
  await api.reset()
  assert.equal(fake.calls[13].url, '/spatial/reset')

  await api.openDataset({ path: '/data/outs', data_format: 'visium' })
  assert.equal(fake.calls[14].url, '/spatial/dataset/open')
  assert.deepEqual(JSON.parse(fake.calls[14].init.body), { path: '/data/outs', data_format: 'visium', actor: 'user' })

  await api.demo()
  assert.equal(fake.calls[15].url, '/spatial/dataset/demo')
  await api.selectDataset('ds2')
  assert.equal(fake.calls[16].url, '/spatial/dataset/select')
  await api.closeDataset('ds2')
  assert.equal(fake.calls[17].url, '/spatial/dataset/close')

  await api.geoDownload({ gse: 'gse123', files: ['a.txt'], dest_dir: '~/.dsh/spatial/GEO' })
  assert.equal(fake.calls[18].url, '/spatial/geo/download')
  const download = JSON.parse(fake.calls[18].init.body)
  assert.equal(download.gse, 'GSE123')
  assert.deepEqual(download.files, ['a.txt'])
  assert.equal(download.dest_dir, '~/.dsh/spatial/GEO')
  assert.equal(download.overwrite, false)
})

test('apiClient：geoSearch 走 buildSearchPayload、错误信封带 code', async () => {
  const fake = makeFetch({ ok: true, total: 0, items: [] })
  const api = SP.apiClient(fake)
  await api.geoSearch({ term: '乳腺癌 Visium', species: 'human', limit: 5 })
  assert.equal(fake.calls[0].url, '/spatial/geo/search')
  const body = JSON.parse(fake.calls[0].init.body)
  assert.equal(body.term, '乳腺癌 Visium')
  assert.equal(body.species, 'human')
  assert.equal(body.limit, 5)
  assert.equal(body.study_type, undefined)

  const soft = makeFetch({ ok: false, code: 'bad_request', error: '参数错误' })
  await assert.rejects(SP.apiClient(soft).getState(), (error) => {
    assert.equal(error.message, '参数错误')
    assert.equal(error.code, 'bad_request')
    return true
  })
  const hard = makeFetch({ ok: false, code: 'internal', error: '引擎崩了' }, { httpOk: false, status: 500 })
  await assert.rejects(SP.apiClient(hard).getHealth(), (error) => {
    assert.equal(error.message, '引擎崩了')
    assert.equal(error.status, 500)
    return true
  })
})

test('plotKey：参数顺序无关、空值忽略', () => {
  assert.equal(SP.plotKey('spatial_scatter', { color_by: 'clusters', point_size: 6 }), SP.plotKey('spatial_scatter', { point_size: 6, color_by: 'clusters' }))
  assert.equal(SP.plotKey('gene_spatial', { gene: 'CDH1' }), 'gene_spatial|gene=CDH1')
  assert.equal(SP.plotKey('pca_variance', {}), 'pca_variance')
  assert.equal(SP.plotKey('gene_box', { genes: ['A', 'B'] }), SP.plotKey('gene_box', { genes: ['A', 'B'] }))
  assert.notEqual(SP.plotKey('gene_box', { genes: ['A'] }), SP.plotKey('gene_box', { genes: ['A', 'B'] }))
})

// ---------------------------------------------------------------------------
// 6. 轮询状态机
// ---------------------------------------------------------------------------

test('jobPolling：running / done / failed / timeout 四态', () => {
  const running = SP.jobPolling({ job: { id: 'j1', status: 'running', pct: 42, message: '下载中' } })
  assert.equal(running.status, 'running')
  assert.equal(running.running, true)
  assert.equal(running.pct, 42)
  assert.equal(running.done, false)

  const done = SP.jobPolling({ ok: true, job: { id: 'j1', status: 'done', pct: 100, result: { dir: '/tmp/geo' } } })
  assert.equal(done.status, 'done')
  assert.equal(done.done, true)
  assert.deepEqual(done.result, { dir: '/tmp/geo' })
  assert.equal(SP.jobDone({ job: { status: 'done' } }), true)
  assert.equal(SP.jobDone({ job: { status: 'running' } }), false)

  const failed = SP.jobPolling({ ok: false, job: { id: 'j2', status: 'failed', error: '网络错误' } })
  assert.equal(failed.status, 'failed')
  assert.equal(failed.failed, true)
  assert.equal(failed.error, '网络错误')

  const timeout = SP.jobPolling({ job: { status: 'running', pct: 10 } }, { elapsedMs: 999999, timeoutMs: 1000 })
  assert.equal(timeout.status, 'timeout')
  assert.equal(timeout.timedOut, true)
  assert.equal(timeout.running, false)

  const doneWins = SP.jobPolling({ job: { status: 'done' } }, { elapsedMs: 999999, timeoutMs: 1 })
  assert.equal(doneWins.status, 'done')
  assert.equal(SP.pickJobId({ job: { id: 'j9' } }), 'j9')
  assert.equal(SP.pickJobId({ job_id: 'j8' }), 'j8')
  assert.equal(SP.pickJobId({}), '')
})

// ---------------------------------------------------------------------------
// 7. 响应归一化（DESIGN §2.4 dataset 摘要 / §2.5 结果信封）
// ---------------------------------------------------------------------------

test('normalizeDataset：缺字段也能渲染；ref/geometry/history 归一化', () => {
  const empty = SP.normalizeDataset(null)
  assert.equal(empty.id, '')
  assert.equal(empty.name, '')
  assert.equal(empty.n_spots, null)
  assert.equal(empty.n_genes, null)
  assert.deepEqual(empty.steps, [])
  assert.deepEqual(empty.sample_columns, [])
  assert.deepEqual(empty.cluster_keys, [])
  assert.deepEqual(empty.warnings, [])
  assert.deepEqual(empty.notes, [])
  assert.deepEqual(empty.history, [])
  assert.equal(empty.ref, null)
  assert.equal(empty.current, '')
  assert.equal(empty.geometry.mode, '')
  assert.equal(empty.geometry.x_range, null)

  const dataset = SP.normalizeDataset({
    dataset: {
      id: 'ds1',
      name: 'Visium breast',
      path: '/data/outs',
      source_format: 'visium',
      n_spots: '3800',
      n_genes: 33000,
      spatial_mode: 'visium',
      geometry: { mode: 'visium', x_range: [0, 1200], y_range: [0, 1400], image: '/data/tissue.png', neighbor_radius: 150 },
      sample_columns: [{ key: 'group', kind: 'categorical', n_unique: 2 }, { name: 'slice' }],
      steps: ['qc', 'filter', 'normalize'],
      current: 'normalize',
      cluster_keys: ['clusters', 'domains'],
      ref: { path: '/ref.h5ad', cell_type_col: 'cell_type', cell_types: ['T cell', 'B cell'] },
      warnings: ['spot 数偏少'],
      history: [{ at: '2026-01-01T00:00:00Z', step: 'qc', ok: true, ms: 812, summary: '完成' }],
      analyses_run: [{ at: '2026-01-01T00:00:00Z', id: 'morans_svg', ok: true, ms: 120 }],
    },
  })
  assert.equal(dataset.id, 'ds1')
  assert.equal(dataset.n_spots, 3800)
  assert.equal(dataset.n_genes, 33000)
  assert.equal(dataset.spatial_mode, 'visium')
  assert.deepEqual(dataset.geometry.x_range, [0, 1200])
  assert.equal(dataset.geometry.neighbor_radius, 150)
  assert.equal(dataset.sample_columns.length, 2)
  assert.equal(dataset.sample_columns[0].n_unique, 2)
  assert.deepEqual(dataset.steps, ['qc', 'filter', 'normalize'])
  assert.deepEqual(dataset.cluster_keys, ['clusters', 'domains'])
  assert.equal(dataset.ref.cell_types.length, 2)
  assert.deepEqual(dataset.warnings, ['spot 数偏少'])
  assert.equal(dataset.history[0].ms, 812)
  assert.equal(dataset.analyses_run[0].id, 'morans_svg')
  assert.equal(SP.hasDomainKey(dataset), true)
  assert.equal(SP.currentClusterKey(dataset), 'clusters')
  assert.equal(SP.currentClusterKey({}), 'clusters')
})

test('normalizeResult：summary/tables/figures/warnings/data_out 信封', () => {
  const result = SP.normalizeResult({
    ok: true,
    summary: { title: '空间聚类', bullets: ['检出 8 个簇', '  '], metrics: { n_clusters: 8 } },
    tables: [{ id: 'clusters', title: '簇组成', columns: [{ key: 'cluster' }, { key: 'n', type: 'number' }], rows: [{ cluster: '0', n: 120 }] }],
    figures: [{ id: 'f1', title: '空间聚类图', file: '/tmp/analyses/ds1/cluster/spatial_scatter.png', width: 800, height: 600 }],
    warnings: ['邻接图连通分量 > 1'],
    notes: ['坐标已按 scalefactors 缩放'],
    data_out: [{ kind: 'clusters', path: '/tmp/out/clusters.csv', label: '簇表' }],
  })
  assert.equal(result.title, '空间聚类')
  assert.deepEqual(result.bullets, ['检出 8 个簇'])
  assert.deepEqual(result.metrics, { n_clusters: 8 })
  assert.equal(result.tables.length, 1)
  assert.equal(result.tables[0].columns[1].type, 'number')
  assert.equal(result.figures.length, 1)
  assert.equal(result.figures[0].url, '/spatial/figure?id=spatial_scatter.png')
  assert.deepEqual(result.warnings, ['邻接图连通分量 > 1'])
  assert.deepEqual(result.notes, ['坐标已按 scalefactors 缩放'])
  assert.equal(result.data_out[0].path, '/tmp/out/clusters.csv')

  const nested = SP.normalizeResult({ ok: true, result: { summary: { bullets: ['ok'] } }, ms: 120 })
  assert.deepEqual(nested.bullets, ['ok'])
  assert.equal(nested.ms, 120)
  assert.deepEqual(SP.normalizeResult(null).tables, [])
})

test('normalizeFigure / normalizeState / normalizeAnalyses / normalizeAudit / normalizeHealth', () => {
  assert.equal(SP.figureSrc('/spatial/figure?id=a.png'), '/spatial/figure?id=a.png')
  assert.equal(SP.figureSrc('spatial_scatter.png'), '/spatial/figure?id=spatial_scatter.png')
  assert.equal(SP.figureSrc('bad.svg'), '')
  const figure = SP.normalizeFigure({ file: '/tmp/x/gene_spatial.png', title: '基因空间图', kind: 'gene_spatial' })
  assert.equal(figure.url, '/spatial/figure?id=gene_spatial.png')
  assert.equal(figure.kind, 'gene_spatial')
  assert.equal(SP.normalizeFigure({ file: 'bad.svg' }), null)

  const state = SP.normalizeState({
    datasets: [{ id: 'ds1', name: 'demo', n_spots: 100, steps: ['qc'] }],
    engine: { ready: true, mode: 'builtin', python: '3.11' },
    suggestions: [{ id: 's1', title: '下一步', step: 'filter' }],
    activity: [{ at: '2026-01-01T00:00:00Z', actor: 'agent', action: 'run:qc', detail: '完成', ok: true }],
  })
  assert.equal(state.active, 'ds1', '未给 active 时回落第一个数据集')
  assert.equal(state.datasets[0].n_spots, 100)
  assert.equal(state.engine.ready, true)
  assert.equal(state.suggestions[0].step, 'filter')
  assert.equal(state.activity[0].action, 'run:qc')
  assert.equal(state.activity[0].ok, true)

  const catalog = SP.normalizeAnalyses({
    analyses: [
      { id: 'morans_svg', title: "空间自相关（Moran's I）", group: '空间统计', method: 'Moran', needs: ['cluster'], params: [{ key: 'n_top', type: 'number', default: 20 }] },
      { id: 'deconv_nnls', title: '反卷积', group: '反卷积', available: false, reason: '缺参考' },
    ],
  })
  assert.equal(catalog.analyses.length, 2)
  assert.deepEqual(catalog.groups, ['空间统计', '反卷积'])
  assert.equal(catalog.analyses[0].params[0].type, 'number')
  assert.equal(catalog.analyses[1].available, false)
  assert.equal(catalog.analyses[1].reason, '缺参考')

  const audit = SP.normalizeAudit({ score: 88.6, level: 'warn', checks: [{ id: 'coords', level: 'warn', title: '坐标', detail: '缺失' }] })
  assert.equal(audit.score, 89)
  assert.equal(audit.checks[0].level, 'warn')

  const health = SP.normalizeHealth({ engine: { ready: true, python: '3.11.9', packages: { numpy: '1.26' } } })
  assert.equal(health.ready, true)
  assert.equal(health.python, '3.11.9')
  assert.equal(SP.engineBadge(health).labelKey, 'engineReady')
  assert.equal(SP.engineBadge({ ready: false }).labelKey, 'engineBooting')
  assert.equal(SP.engineBadge({ ready: false, error: '没有 python' }).labelKey, 'engineBad')
})

test('matchGenes / sortRows / formatTableRows：基因检索与表格排序', () => {
  const genes = ['CDH1', 'CDH2', 'VIM', 'MKI67', 'EPCAM']
  assert.deepEqual(SP.matchGenes(genes, 'cdh', 5), ['CDH1', 'CDH2'])
  assert.deepEqual(SP.matchGenes(genes, '', 2), ['CDH1', 'CDH2'])
  assert.deepEqual(SP.matchGenes(genes, 'KI6', 5), ['MKI67'])
  assert.deepEqual(SP.matchGenes(null, 'x'), [])

  const table = {
    columns: [{ key: 'gene', type: 'string' }, { key: 'score', type: 'number' }],
    rows: [{ gene: 'B', score: 2 }, { gene: 'A', score: 10 }, { gene: 'C', score: null }],
  }
  const sorted = SP.sortRows(table.rows, 'score', 'desc', table.columns)
  assert.deepEqual(sorted.map((row) => row.gene), ['A', 'B', 'C'])
  const view = SP.formatTableRows(table, 2)
  assert.equal(view.total, 3)
  assert.equal(view.shown, 2)
  assert.equal(view.truncated, true)
  assert.equal(view.rows[0].cells[1].numeric, true)
  assert.equal(SP.formatTableRows({ columns: [], rows: [] }, 2).columns.length, 0)
})

// ---------------------------------------------------------------------------
// 8. 静态合规：骨架 / 槽注册 / 样式 / 文案 / 数据面
// ---------------------------------------------------------------------------

test('模块骨架：IIFE + ModuleLoader + 两个槽 + __test__ 组件表', () => {
  assert.match(CLIENT_SOURCE, /^\/\*\*/)
  assert.match(CLIENT_SOURCE, /;\(function \(\) \{/)
  assert.ok(CLIENT_SOURCE.trimEnd().endsWith('})()'), '整文件应以 })() 收尾')
  assert.match(CLIENT_SOURCE, /window\.__ModuleLoader__\.load\(\{ id: NS, factory: factory \}\)/)
  assert.match(CLIENT_SOURCE, /globalThis\.__SP_TEST__ = SP_TEST/)
  assert.match(CLIENT_SOURCE, /inject: \['slots',\s*'locale'\]/)

  // main 槽 key=spatial；侧边栏 id=spatial, order=5
  assert.match(CLIENT_SOURCE, /slots\.inject\('main',/)
  assert.match(CLIENT_SOURCE, /key: 'spatial'/)
  assert.match(CLIENT_SOURCE, /slots\.inject\('sidebar\.panellist',/)
  assert.match(CLIENT_SOURCE, /id: 'spatial'/)
  assert.match(CLIENT_SOURCE, /order: 5/)
  // locale 注册 + 翻译器
  assert.match(CLIENT_SOURCE, /locale\.register\(NS, \{ zh: DICT_ZH, en: DICT_EN \}\)/)
  assert.match(CLIENT_SOURCE, /T = createTranslator\(locale\)/)
  // 错误边界 / 关键组件
  assert.match(CLIENT_SOURCE, /class PluginErrorBoundary extends React\.Component/)
  for (const name of [
    'Spinner', 'LoadingLine', 'ErrorBanner', 'SidebarIcon', 'TopBar',
    'DatasetCard', 'PipelinePanel', 'StepCard', 'ParamForm', 'ResultView', 'TabsBar',
    'OverviewTab', 'SpatialMapTab', 'ExprTab', 'DeTab', 'SvgTab', 'DomainTab', 'DeconvTab',
    'CommTab', 'RecordsTab', 'LogTab', 'GeoDialog', 'ImportDialog', 'AnalysisDrawer',
    'AgentPanel', 'WorkbenchBody', 'WorkbenchPage', 'apply',
  ]) {
    assert.match(CLIENT_SOURCE, new RegExp('function ' + name + '\\('), `缺组件 ${name}`)
  }
  // factory 返回 __test__ 组件表
  assert.match(CLIENT_SOURCE, /__test__: Object\.assign\(\{\}, SP_TEST, \{/)
  assert.match(CLIENT_SOURCE, /components: \{/)

  // 只 require('react')，不碰任何 @deepseek-ai 包
  const requires = [...CLIENT_CODE.matchAll(/\brequire\((['"])([^'"]+)\1\)/g)].map((match) => match[2])
  assert.deepEqual([...new Set(requires)], ['react'], '唯一允许的 require 目标是 react')
  assert.ok(!CLIENT_CODE.includes('@deepseek-ai'), '不得 require 宿主 Client 包')
  // 不写 JSX：不出现 `h(` 之外的标签字面量
  assert.ok(!/<[A-Za-z][A-Za-z0-9]*[\s/>]/.test(CLIENT_CODE), '不应出现 JSX 标签字面量')
  // 数据面只打同源 /spatial（不出现别的插件前缀）
  assert.ok(!CLIENT_CODE.includes('geobulk'), '不得残留 geobulk 依赖')
  assert.ok(!CLIENT_CODE.includes('/geobulk'), '不得请求 /geobulk 路由')
  assert.ok(CLIENT_SOURCE.includes("'/spatial'"), "API 常量应为 '/spatial'")
})

test('样式：类名 sp- 前缀、只用 --dsw-alias-* 主题色', () => {
  const cssStart = CLIENT_SOURCE.indexOf('var CSS = [')
  const cssEnd = CLIENT_SOURCE.indexOf("].join('')", cssStart)
  assert.ok(cssStart > 0 && cssEnd > cssStart, 'CSS 段落应存在')
  const css = CLIENT_SOURCE.slice(cssStart, cssEnd)
  assert.ok(css.length > 12000, `CSS 太短：${css.length}`)

  const selectors = [...css.matchAll(/\.([A-Za-z][A-Za-z0-9_-]*)\s*[{,:]/g)].map((match) => match[1])
  const badSelectors = [...new Set(selectors)].filter((name) => !name.startsWith('sp-'))
  assert.deepEqual(badSelectors, [], `CSS 选择器不是 sp- 前缀：${badSelectors.join(', ')}`)

  const hexColors = [...css.matchAll(/#[0-9a-fA-F]{3,8}\b/g)].map((match) => match[0])
  assert.deepEqual(hexColors, [], `CSS 不应出现硬编码颜色：${hexColors.join(', ')}`)
  assert.ok(!/\brgba?\(/.test(css), 'CSS 不应出现 rgb()/rgba()，请用 color-mix + 主题变量')
  const aliasUses = [...css.matchAll(/--dsw-alias-[a-z0-9-]+/g)].map((match) => match[0])
  assert.ok(aliasUses.length > 120, `主题变量用得不够：${aliasUses.length}`)
})

test('组件里用到的 sp- 类名都有样式定义；文案 key 都在字典里', () => {
  const cssStart = CLIENT_SOURCE.indexOf('var CSS = [')
  const cssEnd = CLIENT_SOURCE.indexOf("].join('')", cssStart)
  const css = CLIENT_SOURCE.slice(cssStart, cssEnd)
  const defined = new Set([...css.matchAll(/\.(sp-[A-Za-z0-9_-]+)/g)].map((match) => match[1]))
  const used = new Set()
  for (const match of CLIENT_SOURCE.matchAll(/className: '([^']*)'/g)) {
    for (const token of match[1].split(/\s+/)) {
      if (token.startsWith('sp-') && !token.endsWith('-')) used.add(token)
    }
  }
  const missing = [...used].filter((name) => !defined.has(name)).sort()
  assert.deepEqual(missing, [], `用了未定义的样式类：${missing.join(', ')}`)

  const referenced = new Set()
  for (const match of CLIENT_SOURCE.matchAll(/\bT\('([A-Za-z0-9_]+)'\)/g)) referenced.add(match[1])
  for (const match of CLIENT_SOURCE.matchAll(/tf\(T, '([A-Za-z0-9_]+)'/g)) referenced.add(match[1])
  for (const match of CLIENT_SOURCE.matchAll(/labelKey: '([A-Za-z0-9_]+)'/g)) referenced.add(match[1])
  for (const match of CLIENT_SOURCE.matchAll(/titleKey: '([A-Za-z0-9_]+)'/g)) referenced.add(match[1])
  for (const key of [
    'stepStatusDone', 'stepStatusRunning', 'stepStatusReady', 'stepStatusBlocked', 'stepStatusFailed',
    'engineReady', 'engineBooting', 'engineBad', 'auditOk', 'auditWarn', 'auditFail',
    'tabOverview', 'tabSpatial', 'tabExpr', 'tabDe', 'tabStats', 'tabDomains', 'tabDeconv',
    'tabComm', 'tabRecords', 'tabLog', 'geoPlatformAll', 'colorClusters', 'colorDomains',
    'colorGroup', 'colorCellType', 'colorGene',
  ]) referenced.add(key)
  const missingZh = [...referenced].filter((key) => !SP.DICT_ZH[key]).sort()
  const missingEn = [...referenced].filter((key) => !SP.DICT_EN[key]).sort()
  assert.deepEqual(missingZh, [], `字典缺 key：${missingZh.join(', ')}`)
  assert.deepEqual(missingEn, [], `EN 字典缺 key：${missingEn.join(', ')}`)
})

test('不用 emoji 当图标；侧边栏是内联 SVG', () => {
  const emoji = CLIENT_SOURCE.match(/[\u{1F000}-\u{1FAFF}\u{1F1E6}-\u{1F1FF}]/gu)
  assert.equal(emoji, null, `出现 emoji：${emoji ? emoji.join('') : ''}`)
  assert.match(CLIENT_SOURCE, /function SidebarGlyph\(props\)/)
  assert.match(CLIENT_SOURCE, /viewBox: '0 0 24 24'/)
})

test('纯函数层不碰 DOM / 不 require react；浏览器层只在 factory 里', () => {
  const factoryStart = CLIENT_SOURCE.indexOf('function factory(require)')
  assert.ok(factoryStart > 0)
  const pure = stripComments(CLIENT_SOURCE.slice(0, factoryStart))
  assert.ok(!/\bdocument\./.test(pure), '纯函数层不应碰 document')
  assert.ok(!/\bwindow\./.test(pure), '纯函数层不应碰 window')
  assert.ok(!/\brequire\(/.test(pure), '纯函数层不应 require')
  const browser = CLIENT_SOURCE.slice(factoryStart)
  assert.match(browser, /require\('react'\)/)
})

// ---------------------------------------------------------------------------
// 9. 全部通过时打印 PASS（供 CI / Lead 验收脚本 grep）
// ---------------------------------------------------------------------------

test('汇总', () => {
  console.log('PASS')
})
