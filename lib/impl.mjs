/**
 * dsh-spatial — Host 业务实现（唯一操作实现）
 *
 * 一条操作只写一次：浏览器 UI 走 HTTP 路由，Agent 走 spatial 工具，
 * 两边入口都调本文件里的 op* 函数（见 cordis-plugin-dev 的「一个操作两个调用方」）。
 *
 * 本文件按需用带 rev 的 URL 动态 import（见 index.js 的 loadImpl），
 * 因此不要在 lib 模块之间写静态 import；需要 engine.mjs 时用 loadEngineMod()。
 */

import { createReadStream, createWriteStream, existsSync, mkdirSync, readdirSync, statSync, unlinkSync } from 'node:fs'
import path from 'node:path'
import { randomUUID } from 'node:crypto'
import { pipeline } from 'node:stream/promises'

const PREFIX = '/spatial'
const MAX_JSON_BODY = 4 * 1024 * 1024
const MAX_UPLOAD_BYTES = 8 * 1024 * 1024 * 1024
const ACTIVITY_LIMIT = 200
const SLOW_STEP_TIMEOUT = 15 * 60 * 1000
const ANALYSIS_TIMEOUT = 20 * 60 * 1000
const DOWNLOAD_TIMEOUT = 45 * 60 * 1000
const REF_TIMEOUT = 10 * 60 * 1000

/** 12 个流程步骤（顺序即 UI 左栏顺序，与 DESIGN §3、sp_registry.STEP_ORDER 一致）。 */
export const STEPS = [
  { id: 'qc', title: '空间质控', requires: [], params: { min_counts: null, min_genes: null, max_pct_mt: 20, group_col: null } },
  { id: 'filter', title: '过滤', requires: ['qc'], params: { min_counts: 500, min_genes: 200, max_pct_mt: 20, min_cells: 3 } },
  { id: 'normalize', title: '标准化', requires: ['filter'], params: { target_sum: 10000, log: true, method: 'lognorm' } },
  { id: 'hvg', title: '高变基因', requires: ['normalize'], params: { n_top_genes: 2000, flavor: 'seurat' } },
  { id: 'pca', title: 'PCA', requires: ['hvg'], params: { n_pcs: 30, scale: true } },
  { id: 'neighbors', title: '邻接图', requires: ['pca'], params: { n_neighbors: 15, metric: 'cosine', spatial_weight: 0.3, spatial_method: 'grid' } },
  { id: 'cluster', title: '聚类', requires: ['neighbors'], params: { resolution: 1.0, method: 'louvain' } },
  { id: 'embed', title: '二维嵌入', requires: ['pca'], params: { layout: 'tsne', perplexity: 30 } },
  { id: 'svg', title: '空间高变基因', requires: ['cluster'], params: { n_top: 30, method: 'moran', permutations: 0 } },
  { id: 'domains', title: '空间域', requires: ['cluster'], params: { n_domains: 6, method: 'mrf', smooth: 0.5, iterations: 20 } },
  { id: 'deconv', title: '反卷积', requires: ['normalize'], params: { ref: '', method: 'nnls', cell_type_col: '', min_genes: 50 } },
  { id: 'annotate', title: '区域注释', requires: ['domains'], params: { mapping: {}, key: 'domain_anno' } },
]
export const STEP_IDS = STEPS.map((s) => s.id)
/** 20 种图（DESIGN §5、sp_registry.PLOT_ORDER）。 */
export const PLOT_KINDS = [
  'qc_spatial', 'qc_violin', 'spatial_scatter', 'embedding', 'cluster_bar', 'gene_spatial',
  'svg_map', 'moran_plot', 'domain_map', 'deconv_spatial', 'deconv_bar', 'nhood_enrich',
  'cooccurrence', 'ripley', 'pca_variance', 'heatmap', 'dotplot', 'gene_box',
  'de_volcano', 'de_heatmap',
]
/** 导出类型（与 sp_registry.EXPORT_KINDS 一致）。 */
export const EXPORT_KINDS = [
  'expr_csv', 'svg_csv', 'domain_csv', 'deconv_csv', 'de_csv', 'sample_info_csv',
  'summary_json', 'protocol', 'methods', 'report', 'script',
]
/** 数据格式（与 sp_io.SUPPORTED_FORMATS 一致，DESIGN §6）。 */
export const DATA_FORMATS = [
  'auto', 'visium', 'visium_h5', 'h5ad', 'stereo', 'slideseq',
  'matrix_text', 'xlsx', 'geo_series', 'archive', 'rds', 'rdata',
]

// ---------------------------------------------------------------------------
// HTTP 小工具
// ---------------------------------------------------------------------------

const LOOPBACK_HOSTS = new Set(['127.0.0.1', 'localhost', '[::1]', '::1', '0.0.0.0'])

export function hostIsLoopback(req) {
  const raw = String((req && req.headers && req.headers.host) || '')
  if (!raw) return true
  const name = raw.startsWith('[') ? raw.slice(0, raw.indexOf(']') + 1) : raw.split(':')[0]
  return LOOPBACK_HOSTS.has(name.toLowerCase())
}

function sendJson(res, status, value) {
  const body = Buffer.from(JSON.stringify(value), 'utf8')
  res.writeHead(status, {
    'content-type': 'application/json; charset=utf-8',
    'cache-control': 'no-store',
    'content-length': String(body.length),
  })
  res.end(body)
}

function opError(code, message, status = 400, details = null) {
  const body = { ok: false, code, error: message }
  if (details && typeof details === 'object') Object.assign(body, details)
  return { status, body }
}

function toFailure(error) {
  const code = (error && error.code) || 'internal'
  const message = String((error && error.message) || error)
  const details = error && error.details && typeof error.details === 'object' ? error.details : null
  if (code === 'prerequisite' || code === 'bad_request' || code === 'not_found' || code === 'unsupported' || code === 'ambiguous') {
    return opError(code, message, 400, details)
  }
  if (code === 'timeout') return opError(code, message, 504, details)
  if (code === 'no_python' || code === 'no_r' || code === 'engine_exited' || code === 'internal' || code === 'disposed') {
    return opError(code, message, 500, details)
  }
  if (code === 'busy') return opError(code, message, 429, details)
  return opError(code, message, 500, details)
}

function readBody(req, limit) {
  return new Promise((resolve, reject) => {
    const chunks = []
    let size = 0
    req.on('data', (chunk) => {
      size += chunk.length
      if (size > limit) {
        reject(Object.assign(new Error(`请求体超过上限（${Math.round(limit / 1024)} KB）`), { code: 'bad_request' }))
        req.destroy()
        return
      }
      chunks.push(chunk)
    })
    req.on('end', () => resolve(Buffer.concat(chunks)))
    req.on('error', (error) => reject(Object.assign(error, { code: 'bad_request' })))
  })
}

async function readJson(req) {
  const buffer = await readBody(req, MAX_JSON_BODY)
  if (!buffer.length) return {}
  try {
    const value = JSON.parse(buffer.toString('utf8'))
    return value && typeof value === 'object' ? value : {}
  } catch (error) {
    throw Object.assign(new Error(`请求 JSON 解析失败：${error.message}`), { code: 'bad_request' })
  }
}

function sanitizeName(name) {
  const base = String(name || 'upload')
  return path.basename(base).replace(/[^\w.\-\u4e00-\u9fa5]+/g, '_').slice(0, 160) || 'upload'
}

function nowIso() {
  return new Date().toISOString()
}

async function loadEngineMod() {
  return import(new URL(`./engine.mjs?rev=${Date.now()}-${Math.random().toString(36).slice(2)}`, import.meta.url).href)
}

// ---------------------------------------------------------------------------
// 引擎
// ---------------------------------------------------------------------------

async function ensureEngine(runtime, { allowStart = true } = {}) {
  if (runtime.engine) {
    if (runtime.engine.running) return runtime.engine
    if (!allowStart) return runtime.engine
  } else {
    const mod = await loadEngineMod()
    const config = runtime.config && typeof runtime.config === 'object' ? runtime.config : {}
    runtime.engine = mod.createEngine({
      python: typeof config.python === 'string' && config.python.trim() ? config.python.trim() : undefined,
      workdir: typeof config.workdir === 'string' && config.workdir.trim() ? config.workdir.trim() : undefined,
      timeout: Number(config.timeoutMs) > 0 ? Number(config.timeoutMs) : undefined,
      config,
    })
    if (typeof runtime.ctx?.effect === 'function') {
      const engine = runtime.engine
      runtime.ctx.effect(() => () => { engine.dispose().catch(() => {}) }, 'dsh-spatial: engine')
    }
  }
  if (!allowStart) return runtime.engine
  const now = Date.now()
  if (runtime.engineFail && now - runtime.engineFail.at < 60_000) throw runtime.engineFail.error
  try {
    await runtime.engine.start()
    runtime.engineFail = null
  } catch (error) {
    runtime.engineFail = { at: now, error }
    throw error
  }
  return runtime.engine
}

export async function engineStatus(runtime) {
  let engine = runtime.engine
  if (!engine) {
    try {
      await ensureEngine(runtime, { allowStart: false })
      engine = runtime.engine
    } catch (error) {
      return { ready: false, mode: 'none', error: String(error && error.message ? error.message : error), python: null, packages: {}, ref: null }
    }
  }
  if (!engine) return { ready: false, mode: 'none', error: 'Python 引擎未初始化', python: null, packages: {}, ref: null }
  const status = typeof engine.status === 'function' ? engine.status() : { ready: false }
  if (!status.ready && !status.error && runtime.engineFail) {
    status.error = String(runtime.engineFail.error && runtime.engineFail.error.message ? runtime.engineFail.error.message : runtime.engineFail.error)
  }
  return status
}

/** 引擎请求（带长超时选项）。 */
async function engineRequest(runtime, cmd, args, timeout) {
  const engine = await ensureEngine(runtime)
  return engine.request(cmd, args, timeout ? { timeout } : {})
}

/** ping 结果缓存（30s）：健康检查/状态轮询不必每次都跑完整探测。 */
async function enginePing(runtime) {
  const now = Date.now()
  if (runtime.pingInfo && now - runtime.pingAt < 30_000) return runtime.pingInfo
  try {
    const info = await engineRequest(runtime, 'ping', {}, 90_000)
    runtime.pingInfo = info
    runtime.pingAt = now
    return info
  } catch (error) {
    if (runtime.pingInfo) return runtime.pingInfo
    throw error
  }
}

// ---------------------------------------------------------------------------
// runtime 状态：活动日志 / 任务
// ---------------------------------------------------------------------------

export function logActivity(runtime, entry) {
  if (!runtime.activity) runtime.activity = []
  const item = {
    at: nowIso(),
    actor: entry.actor || 'system',
    action: String(entry.action || ''),
    dataset: entry.dataset || null,
    detail: entry.detail == null ? '' : String(entry.detail).slice(0, 400),
    ok: entry.ok !== false,
    ms: Number.isFinite(entry.ms) ? entry.ms : null,
  }
  runtime.activity.push(item)
  if (runtime.activity.length > ACTIVITY_LIMIT) runtime.activity = runtime.activity.slice(-ACTIVITY_LIMIT)
  return item
}

export function activityOf(runtime, since = 0) {
  const list = runtime.activity || []
  return { entries: list.slice(Math.max(0, Number(since) || 0)) }
}

function jobs(runtime) {
  if (!runtime.jobs) runtime.jobs = new Map()
  return runtime.jobs
}

function newJob(runtime, kind, meta = {}) {
  const id = randomUUID().slice(0, 10)
  const job = {
    id,
    kind,
    status: 'running',
    pct: 0,
    message: '准备中…',
    startedAt: nowIso(),
    finishedAt: null,
    result: null,
    error: null,
    code: null,
    ...meta,
  }
  jobs(runtime).set(id, job)
  return job
}

function finishJob(job, { ok, result, error, code }) {
  job.status = ok ? 'done' : 'failed'
  job.finishedAt = nowIso()
  job.pct = ok ? 100 : job.pct
  if (ok) job.result = result
  else {
    job.error = String(error || '未知错误')
    job.code = code || 'internal'
  }
  return job
}

function jobView(runtime, id) {
  const job = jobs(runtime).get(String(id || ''))
  if (!job) throw Object.assign(new Error('任务不存在或已过期'), { code: 'not_found' })
  // 下载中的进度从引擎 lastProgress 同步
  if (job.status === 'running' && runtime.engine && runtime.engine.progress) {
    const progress = runtime.engine.progress
    if (progress && progress.message) job.message = progress.message
    if (progress && Number.isFinite(progress.pct)) job.pct = Math.max(0, Math.min(99, progress.pct))
  }
  return { job }
}

// ---------------------------------------------------------------------------
// 基础操作
// ---------------------------------------------------------------------------

export async function opHealth(runtime) {
  let engine = await engineStatus(runtime)
  // 健康检查同时负责「首次拉起 worker」：GET /health 必须能反映真实就绪状态（DESIGN §12）。
  if (!engine.ready && !engine.error) {
    try {
      await ensureEngine(runtime)
    } catch (error) {
      engine.error = String(error && error.message ? error.message : error)
      engine.mode = 'none'
    }
    engine = await engineStatus(runtime)
  }
  let ping = null
  if (engine.ready) {
    try { ping = await enginePing(runtime) } catch { ping = null }
  }
  return {
    ok: true,
    engine: {
      ready: Boolean(engine.ready),
      mode: (ping && ping.mode) || engine.mode || (engine.ready ? 'builtin' : 'none'),
      python: (ping && ping.python) || engine.python || null,
      version: engine.version || null,
      packages: (ping && ping.packages) || engine.packages || {},
      modules: (ping && ping.modules) || null,
      ref: (ping && ping.ref) || null,
      capabilities: (ping && ping.capabilities) || null,
      error: engine.error || null,
      lastProgress: engine.lastProgress || null,
    },
    workdir: runtime.engine ? runtime.engine.workdir : null,
    geo: true,
  }
}

export async function opState(runtime) {
  const engine = await engineStatus(runtime)
  let datasets = []
  let active = null
  let suggestions = []
  let recipes = []
  if (engine.ready) {
    try {
      const list = await engineRequest(runtime, 'list', {})
      datasets = list.datasets || []
      active = list.active || null
      if (active) {
        const sug = await engineRequest(runtime, 'suggest', { dataset: active })
        suggestions = sug.suggestions || []
      }
      const rec = await engineRequest(runtime, 'recipe', {})
      recipes = rec.recipes || []
    } catch { /* 引擎忙/未就绪时保持空 */ }
  }
  let refs = { available: false, origin: null, refs: [], error: null }
  try {
    refs = await refsPayload(runtime)
  } catch { /* ignore */ }
  let ping = null
  if (engine.ready) {
    try { ping = await enginePing(runtime) } catch { ping = null }
  }
  return {
    ok: true,
    active,
    datasets,
    engine: {
      ready: Boolean(engine.ready),
      mode: (ping && ping.mode) || engine.mode || (engine.ready ? 'builtin' : 'none'),
      python: (ping && ping.python) || engine.python || null,
      version: engine.version || null,
      packages: (ping && ping.packages) || engine.packages || {},
      modules: (ping && ping.modules) || null,
      ref: (ping && ping.ref) || null,
      capabilities: (ping && ping.capabilities) || null,
      error: engine.error || null,
    },
    steps: (ping && ping.capabilities && ping.capabilities.steps) || STEPS.map((s) => ({ ...s, available: false })),
    recipes,
    suggestions,
    activity: (runtime.activity || []).slice(-50),
    refs,
    workdir: runtime.engine ? runtime.engine.workdir : null,
  }
}

// ---------------------------------------------------------------------------
// GEO
// ---------------------------------------------------------------------------

export async function opGeoSearch(runtime, args = {}) {
  const started = Date.now()
  const engine = await ensureEngine(runtime)
  const payload = {
    term: String(args.term || ''),
    species: args.species || 'all',
    study_type: args.study_type || 'all',
    min_samples: Number(args.min_samples) || 0,
    platform: args.platform || 'all',
    limit: Number(args.limit) || 20,
    page: Number(args.page) || 0,
    sort: args.sort || 'relevance',
  }
  if (args.has_matrix !== undefined && args.has_matrix !== null) payload.has_matrix = Boolean(args.has_matrix)
  const result = await engine.request('geo_search', payload, { timeout: 120000 })
  logActivity(runtime, {
    actor: args.actor || 'user',
    action: 'geo:search',
    detail: `${args.term || ''} → ${result.total || 0} 条`,
    ok: true,
    ms: Date.now() - started,
  })
  return { ok: true, ...result }
}

export async function opGeoDetail(runtime, args = {}) {
  const started = Date.now()
  const engine = await ensureEngine(runtime)
  const result = await engine.request('geo_detail', { gse: String(args.gse || '') }, { timeout: 180000 })
  logActivity(runtime, {
    actor: args.actor || 'user',
    action: 'geo:detail',
    detail: String(args.gse || ''),
    ok: true,
    ms: Date.now() - started,
  })
  return { ok: true, ...result }
}

/** 异步下载：立即返回 job，前端轮询 /jobs。 */
export async function opGeoDownload(runtime, args = {}) {
  const gse = String(args.gse || '')
  if (!gse) throw Object.assign(new Error('请提供 GSE 编号'), { code: 'bad_request' })
  const destDir = String(args.dest_dir || '').trim()
  const job = newJob(runtime, 'geo_download', { gse, destDir })
  const actor = args.actor || 'user'
  ;(async () => {
    const started = Date.now()
    try {
      const engine = await ensureEngine(runtime)
      job.message = '连接 NCBI…'
      const result = await engine.request('geo_download', {
        gse,
        files: Array.isArray(args.files) ? args.files : [],
        dest_dir: destDir || undefined,
        overwrite: Boolean(args.overwrite),
      }, { timeout: DOWNLOAD_TIMEOUT })
      job.pct = 100
      job.message = '下载完成'
      let opened = null
      if (args.open_after && result && Array.isArray(result.files)) {
        const first = result.files.find((f) => f.ok !== false)
        if (first) {
          try {
            opened = await opOpen(runtime, { path: first.path, name: `${gse}`, actor })
            await engine.request('run', { dataset: opened.dataset.id, step: 'qc', params: {} }, { timeout: SLOW_STEP_TIMEOUT })
            job.message = '已下载并导入，质控完成'
          } catch (error) {
            job.message = `下载完成，但自动导入失败：${error && error.message ? error.message : error}`
          }
        }
      }
      finishJob(job, { ok: true, result: { ...result, opened } })
      logActivity(runtime, {
        actor,
        action: 'geo:download',
        detail: `${gse} → ${result.dir || destDir}（${(result.files || []).length} 个文件）`,
        ok: true,
        ms: Date.now() - started,
      })
    } catch (error) {
      finishJob(job, { ok: false, error: error && error.message ? error.message : String(error), code: error && error.code })
      logActivity(runtime, {
        actor,
        action: 'geo:download',
        detail: `${gse}：${error && error.message ? error.message : error}`,
        ok: false,
        ms: Date.now() - started,
      })
    }
  })()
  return { ok: true, job: jobView(runtime, job.id).job }
}

export function opJob(runtime, id) {
  return { ok: true, ...jobView(runtime, id) }
}

export function opJobs(runtime) {
  return { ok: true, jobs: Array.from(jobs(runtime).values()).slice(-30) }
}

// ---------------------------------------------------------------------------
// 数据集
// ---------------------------------------------------------------------------

export async function opOpen(runtime, args = {}) {
  const started = Date.now()
  const engine = await ensureEngine(runtime)
  const payload = {
    path: args.path || undefined,
    url: args.url || undefined,
    name: args.name || undefined,
    data_format: args.data_format || args.format || 'auto',
    member: args.member || undefined,
    object: args.object || undefined,
    transpose: args.transpose === undefined ? undefined : Boolean(args.transpose),
    delimiter: args.delimiter || undefined,
    sheet: args.sheet || undefined,
    gene_col: args.gene_col || undefined,
    coords: args.coords || undefined,
    sample_sheet: args.sample_sheet || undefined,
    image: args.image || undefined,
    visium_scale: args.visium_scale === undefined ? undefined : Number(args.visium_scale),
    species: args.species || 'auto',
    demo: Boolean(args.demo),
  }
  if (!Number.isFinite(payload.visium_scale)) payload.visium_scale = undefined
  for (const key of Object.keys(payload)) if (payload[key] === undefined) delete payload[key]
  const result = await engine.request('open', payload, { timeout: 20 * 60 * 1000 })
  const dataset = result.dataset
  logActivity(runtime, {
    actor: args.actor || 'user',
    action: 'dataset:open',
    dataset: dataset.id,
    detail: `${dataset.name}（${dataset.n_spots} spot × ${dataset.n_genes} 基因，格式 ${dataset.source_format || 'auto'}）`,
    ok: true,
    ms: Date.now() - started,
  })
  return { ok: true, ...result }
}

export async function opDemo(runtime, args = {}) {
  const started = Date.now()
  const engine = await ensureEngine(runtime)
  const result = await engine.request('open', { demo: true, name: args.name || 'demo-spatial' }, { timeout: 120000 })
  logActivity(runtime, {
    actor: args.actor || 'user',
    action: 'dataset:demo',
    dataset: result.dataset.id,
    detail: `${result.dataset.n_spots} spot × ${result.dataset.n_genes} 基因（合成空间组织）`,
    ok: true,
    ms: Date.now() - started,
  })
  return { ok: true, ...result }
}

export async function opUpload(runtime, { filename, buffer, actor } = {}) {
  const started = Date.now()
  const engine = await ensureEngine(runtime)
  const safe = sanitizeName(filename)
  const dest = path.join(engine.dirs.uploads, `${Date.now().toString(36)}-${randomUUID().slice(0, 8)}-${safe}`)
  if (buffer) {
    await new Promise((resolve, reject) => {
      const stream = createWriteStream(dest)
      stream.on('error', reject)
      stream.on('finish', resolve)
      stream.end(buffer)
    })
  }
  try {
    return await opOpen(runtime, { path: dest, name: safe, actor: actor || 'user' })
  } catch (error) {
    try { unlinkSync(dest) } catch { /* ignore */ }
    throw error
  } finally {
    logActivity(runtime, {
      actor: actor || 'user',
      action: 'dataset:upload',
      detail: safe,
      ok: true,
      ms: Date.now() - started,
    })
  }
}

export async function opSelect(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  const list = await engine.request('list', {})
  const target = String(args.dataset || '')
  if (!list.datasets.some((item) => item.id === target)) {
    throw Object.assign(new Error(`数据集 ${target} 不存在`), { code: 'not_found' })
  }
  await engine.request('select', { dataset: target })
  runtime.activeDataset = target
  return { ok: true, active: target, datasets: list.datasets }
}

export async function opClose(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  const result = await engine.request('close', { dataset: args.dataset || undefined })
  logActivity(runtime, { actor: args.actor || 'user', action: 'dataset:close', dataset: args.dataset || null, ok: true })
  return { ok: true, ...result }
}

export async function opPreview(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  const result = await engine.request('preview', { dataset: args.dataset || undefined })
  return { ok: true, ...result }
}

export async function opInspect(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  const result = await engine.request('inspect', {
    path: args.path || undefined,
    url: args.url || undefined,
    data_format: args.data_format || 'auto',
    member: args.member || undefined,
    object: args.object || undefined,
    delimiter: args.delimiter || undefined,
    sheet: args.sheet || undefined,
    gene_col: args.gene_col || undefined,
    coords: args.coords || undefined,
    image: args.image || undefined,
    visium_scale: args.visium_scale === undefined ? undefined : Number(args.visium_scale),
  }, { timeout: 5 * 60 * 1000 })
  return { ok: true, ...result }
}

export async function opFormats(runtime) {
  let formats = []
  try {
    const engine = await ensureEngine(runtime)
    const result = await engine.request('formats', {})
    formats = result.formats || []
  } catch { /* 引擎未就绪时给静态清单 */ }
  if (!formats.length) {
    formats = DATA_FORMATS.map((id) => ({ id, label: id, hint: '' }))
  }
  return { ok: true, formats }
}

export async function opSamples(runtime, args = {}) {
  const started = Date.now()
  const engine = await ensureEngine(runtime)
  const payload = { dataset: args.dataset || undefined, op: args.op || 'list' }
  for (const key of ['samples', 'groups', 'column', 'group_regex', 'samples_op_name', 'name', 'source', 'path', 'key']) {
    if (args[key] !== undefined) payload[key] = args[key]
  }
  const result = await engine.request('samples', payload, { timeout: 5 * 60 * 1000 })
  logActivity(runtime, {
    actor: args.actor || 'user',
    action: `samples:${payload.op}`,
    dataset: args.dataset || null,
    detail: payload.op === 'list' ? '' : `${(result.sample_info || []).length} 个 spot`,
    ok: true,
    ms: Date.now() - started,
  })
  return { ok: true, ...result }
}

export async function opRun(runtime, args = {}) {
  const started = Date.now()
  const step = String(args.step || '')
  if (!STEP_IDS.includes(step)) throw Object.assign(new Error(`未知步骤：${step}`), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const result = await engine.request('run', {
    dataset: args.dataset || undefined,
    step,
    params: args.params && typeof args.params === 'object' ? args.params : {},
  }, { timeout: SLOW_STEP_TIMEOUT })
  logActivity(runtime, {
    actor: args.actor || 'agent',
    action: `run:${step}`,
    dataset: result.dataset ? result.dataset.id : args.dataset || null,
    detail: ((result.summary && result.summary.bullets) || []).slice(0, 2).join('；'),
    ok: true,
    ms: Date.now() - started,
  })
  return { ok: true, ...result }
}

export async function opPlot(runtime, args = {}) {
  const started = Date.now()
  const kind = String(args.kind || '')
  if (!PLOT_KINDS.includes(kind)) throw Object.assign(new Error(`未知图类型：${kind}`), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const result = await engine.request('plot', {
    dataset: args.dataset || undefined,
    kind,
    params: args.params && typeof args.params === 'object' ? args.params : {},
  }, { timeout: 5 * 60 * 1000 })
  const figure = { id: result.id, url: `${PREFIX}/figure?id=${encodeURIComponent(result.id + '.png')}`, file: result.file, width: result.width, height: result.height, kind }
  figure.url = `${PREFIX}/figure?id=${encodeURIComponent(path.basename(result.file))}`
  logActivity(runtime, { actor: args.actor || 'user', action: `plot:${kind}`, dataset: args.dataset || null, ok: true, ms: Date.now() - started })
  return { ok: true, figure, ...result }
}

export async function opAnalyses(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  const result = await engine.request('analyses', { dataset: args.dataset || undefined })
  return { ok: true, ...result }
}

export async function opAnalysis(runtime, args = {}) {
  const started = Date.now()
  const id = String(args.analysis || args.id || '')
  if (!id) throw Object.assign(new Error('请提供分析 id'), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const result = await engine.request('analysis', {
    dataset: args.dataset || undefined,
    id,
    params: args.params && typeof args.params === 'object' ? args.params : {},
  }, { timeout: ANALYSIS_TIMEOUT })
  logActivity(runtime, {
    actor: args.actor || 'agent',
    action: `analysis:${id}`,
    dataset: result.dataset ? result.dataset.id : args.dataset || null,
    detail: ((result.summary && result.summary.bullets) || []).slice(0, 2).join('；'),
    ok: true,
    ms: Date.now() - started,
  })
  return { ok: true, ...result }
}

export async function opExport(runtime, args = {}) {
  const started = Date.now()
  const what = String(args.what || 'summary_json')
  if (!EXPORT_KINDS.includes(what)) throw Object.assign(new Error(`未知导出类型：${what}`), { code: 'bad_request' })
  const engine = await ensureEngine(runtime)
  const result = await engine.request('export', {
    dataset: args.dataset || undefined,
    what,
    path: args.out_path || args.path || undefined,
  }, { timeout: 5 * 60 * 1000 })
  logActivity(runtime, {
    actor: args.actor || 'user',
    action: `export:${what}`,
    dataset: args.dataset || null,
    detail: result.file,
    ok: true,
    ms: Date.now() - started,
  })
  return { ok: true, ...result }
}

export async function opAudit(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  const result = await engine.request('audit', { dataset: args.dataset || undefined })
  return { ok: true, ...result }
}

export async function opSuggest(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  const result = await engine.request('suggest', { dataset: args.dataset || undefined })
  return { ok: true, ...result }
}

export async function opRecipe(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  if (!args.recipe && !args.id) {
    const result = await engine.request('recipe', {})
    return { ok: true, ...result }
  }
  const result = await engine.request('recipe', { id: args.recipe || args.id })
  return { ok: true, recipe: result }
}

/** 运行预设流程：顺序执行，每步一条活动日志（便于失败恢复）。 */
export async function opRecipeRun(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  const recipeId = String(args.recipe || args.id || 'standard')
  const definition = await engine.request('recipe', { id: recipeId })
  const steps = definition.steps || []
  const results = []
  for (const item of steps) {
    try {
      const result = await engine.request('run', {
        dataset: args.dataset || undefined,
        step: item.step,
        params: item.params || {},
      }, { timeout: SLOW_STEP_TIMEOUT })
      results.push({ step: item.step, ok: true, summary: result.summary })
      logActivity(runtime, {
        actor: args.actor || 'user',
        action: `run:${item.step}`,
        dataset: result.dataset ? result.dataset.id : null,
        detail: ((result.summary && result.summary.bullets) || [])[0] || '',
        ok: true,
        ms: result.ms,
      })
    } catch (error) {
      results.push({ step: item.step, ok: false, error: error && error.message ? error.message : String(error) })
      break
    }
  }
  return { ok: true, recipe: recipeId, results, done: results.filter((r) => r.ok).length, total: steps.length }
}

/** 开发/排障用：重启 Python worker（重载 lib/*.py 改动）。 */
export async function opEngineRestart(runtime) {
  const engine = runtime.engine
  runtime.engine = null
  runtime.pingInfo = null
  runtime.pingAt = 0
  runtime.engineFail = null
  if (engine && typeof engine.dispose === 'function') {
    try { await engine.dispose() } catch { /* ignore */ }
  }
  const fresh = await ensureEngine(runtime)
  const ping = await fresh.request('ping', {}, { timeout: 90_000 })
  return { ok: true, engine: { ready: true, python: ping.python, ref: ping.ref || null, modules: ping.modules } }
}

export async function opReset(runtime, args = {}) {
  const engine = await ensureEngine(runtime)
  await engine.request('reset', {}).catch(() => {})
  runtime.activity = []
  jobs(runtime).clear()
  logActivity(runtime, { actor: args.actor || 'user', action: 'reset', ok: true })
  return { ok: true, datasets: [] }
}

export function opFigure(runtime, id) {
  const name = sanitizeName(String(id || ''))
  if (!/^[A-Za-z0-9._-]+\.png$/.test(name)) {
    throw Object.assign(new Error('非法图片 id'), { code: 'bad_request' })
  }
  const engine = runtime.engine
  if (!engine) throw Object.assign(new Error('引擎未就绪'), { code: 'no_python' })
  const candidates = [
    path.join(engine.workdir || engine.dirs?.root || '', 'analyses'),
  ]
  // 图可能落在 analyses/<dataset>/<step|an|plots>/ 下：从表格缓存里找不到时按文件名递归找一层
  const found = findFigure(engine, name)
  if (found) return { file: found }
  throw Object.assign(new Error('图片不存在（可能已被清理，请重新生成）'), { code: 'not_found' })
}

function findFigure(engine, name) {
  if (!engine) return null
  const root = engine.workdir || (engine.dirs && engine.dirs.root)
  if (!root) return null
  const stack = [path.join(root, 'analyses'), path.join(root, 'figures'), path.join(root, 'exports')]
  const seen = new Set()
  while (stack.length) {
    const dir = stack.pop()
    if (!dir || seen.has(dir)) continue
    seen.add(dir)
    let entries
    try { entries = readdirSync(dir, { withFileTypes: true }) } catch { continue }
    for (const entry of entries) {
      const full = path.join(dir, entry.name)
      if (entry.isDirectory()) stack.push(full)
      else if (entry.name === name) return full
    }
  }
  return null
}

// ---------------------------------------------------------------------------
// 单细胞参考（反卷积用 h5ad，可取自 dsh-singlecell 工作台）
// ---------------------------------------------------------------------------

async function singlecellOrigin(runtime) {
  try {
    const ctx = runtime.ctx
    const webServer = ctx && typeof ctx.get === 'function' ? ctx.get('webServer') : ctx && ctx.webServer
    if (webServer && typeof webServer.getOrigin === 'function') {
      const origin = webServer.getOrigin()
      if (origin) return String(origin).replace(/\/$/, '')
    }
  } catch { /* ignore */ }
  return null
}

async function fetchJson(url, options = {}, timeoutMs = 15000) {
  const controller = new AbortController()
  const timer = setTimeout(() => controller.abort(), timeoutMs)
  try {
    const response = await fetch(url, { ...options, signal: controller.signal })
    const text = await response.text()
    let parsed = null
    try { parsed = text ? JSON.parse(text) : null } catch { parsed = null }
    return { status: response.status, body: parsed, text }
  } finally {
    clearTimeout(timer)
  }
}

/** 探测可用参考：本机 h5ad（引擎扫）+ dsh-singlecell 工作台的数据集列表。 */
async function refsPayload(runtime) {
  const origin = await singlecellOrigin(runtime)
  const result = { available: false, origin, refs: [], datasets: [], error: null, reason: null }
  // 本地 h5ad（引擎 refs 无 path 分支）
  try {
    const engine = await ensureEngine(runtime, { allowStart: false })
    if (engine && engine.running) {
      const local = await engine.request('refs', {}, { timeout: 60000 })
      result.refs = local.refs || []
      result.available = Boolean(local.available)
      result.reason = local.reason || null
    }
  } catch { /* ignore */ }
  if (!origin) {
    if (!result.available) result.error = '未找到宿主 origin'
    return result
  }
  try {
    const response = await fetchJson(`${origin}/singlecell/state`, {}, 8000)
    if (response.status === 200 && response.body && response.body.ok) {
      result.available = true
      result.datasets = (response.body.datasets || []).map((item) => ({
        id: item.id,
        name: item.name,
        n_obs: item.n_obs,
        n_vars: item.n_vars,
        steps: item.steps || [],
      }))
    } else {
      result.error = '单细胞插件未响应（仍可用本地 h5ad）'
    }
  } catch (error) {
    result.error = String(error && error.message ? error.message : error)
  }
  return result
}

/** GET /refs：无 path 列出可用参考；带 path 读该 h5ad 的元信息（细胞类型列、obs 列）。 */
export async function opRefs(runtime, args = {}) {
  if (args.path) {
    const engine = await ensureEngine(runtime)
    const info = await engine.request('refs', { path: String(args.path) }, { timeout: 180000 })
    return { ok: true, info }
  }
  const refs = await refsPayload(runtime)
  return { ok: true, ...refs }
}

/** POST /ref/attach：把单细胞参考挂到当前数据集（本地 h5ad，或从 dsh-singlecell 导出一份）。 */
export async function opRefAttach(runtime, args = {}) {
  const started = Date.now()
  const engine = await ensureEngine(runtime)
  let refPath = String(args.ref_path || '').trim()
  let source = refPath ? 'local' : ''
  if (!refPath && args.sc_dataset) {
    const origin = await singlecellOrigin(runtime)
    if (!origin) throw Object.assign(new Error('单细胞工作台不可用（未找到宿主 origin）'), { code: 'unsupported' })
    const target = path.join(engine.workdir, 'sc-ref', `${sanitizeName(args.sc_dataset)}.h5ad`)
    mkdirSync(path.dirname(target), { recursive: true })
    const response = await fetchJson(`${origin}/singlecell/export`, {
      method: 'POST',
      headers: { 'content-type': 'application/json' },
      body: JSON.stringify({ dataset: args.sc_dataset, what: 'h5ad', path: target }),
    }, REF_TIMEOUT)
    if (response.status !== 200 || !response.body || !response.body.ok) {
      throw Object.assign(new Error(`单细胞工作台导出失败：${(response.body && (response.body.error || response.body.code)) || response.status}`), { code: 'unsupported' })
    }
    refPath = response.body.file || target
    source = 'singlecell'
  }
  if (!refPath) throw Object.assign(new Error('请提供 ref_path（参考 h5ad 绝对路径）或 sc_dataset（dsh-singlecell 的数据集 id）'), { code: 'bad_request' })
  const payload = { dataset: args.dataset || undefined, ref_path: refPath }
  if (args.cell_type_col) payload.cell_type_col = String(args.cell_type_col)
  const result = await engine.request('ref_set', payload, { timeout: 60000 })
  const info = result.info || null
  logActivity(runtime, {
    actor: args.actor || 'user',
    action: 'ref:attach',
    dataset: args.dataset || null,
    detail: `${source}: ${refPath}`,
    ok: true,
    ms: Date.now() - started,
  })
  return { ok: true, ...result, info, source }
}

// ---------------------------------------------------------------------------
// Agent 工具
// ---------------------------------------------------------------------------

/** 工具输出必须是 lossless JSON：去掉 undefined、把 NaN/Infinity 转 null。 */
function sanitizeValue(value, depth = 0) {
  if (depth > 16) return null
  if (value === undefined) return undefined
  if (typeof value === 'number') return Number.isFinite(value) ? value : null
  if (typeof value === 'bigint') return Number(value)
  if (typeof value === 'function' || typeof value === 'symbol') return null
  if (Array.isArray(value)) {
    return value.map((item) => {
      const next = sanitizeValue(item, depth + 1)
      return next === undefined ? null : next
    })
  }
  if (value && typeof value === 'object') {
    const out = {}
    for (const [key, item] of Object.entries(value)) {
      const next = sanitizeValue(item, depth + 1)
      if (next !== undefined) out[key] = next
    }
    return out
  }
  return value
}

function compactResult(result) {
  if (!result || typeof result !== 'object') return result
  const out = { ...result }
  if (Array.isArray(out.tables)) {
    out.tables = out.tables.map((table) => ({
      id: table.id,
      title: table.title,
      columns: table.columns,
      total: table.total,
      rows: Array.isArray(table.rows) ? table.rows.slice(0, 20) : [],
      truncated: table.rows && table.rows.length > 20 ? true : table.truncated,
    }))
  }
  if (out.dataset && typeof out.dataset === 'object') {
    out.dataset = {
      id: out.dataset.id,
      name: out.dataset.name,
      source_format: out.dataset.source_format,
      n_spots: out.dataset.n_spots,
      n_genes: out.dataset.n_genes,
      species: out.dataset.species,
      geometry: out.dataset.geometry,
      steps: out.dataset.steps,
      current: out.dataset.current,
      cluster_keys: out.dataset.cluster_keys,
      sample_columns: out.dataset.sample_columns,
    }
  }
  if (out.preview && typeof out.preview === 'object') {
    out.preview = {
      genes: (out.preview.genes || []).slice(0, 12),
      spots: (out.preview.spots || []).slice(0, 5),
      columns: out.preview.columns,
      coords_meta: out.preview.coords_meta,
      spatial_mode: out.preview.spatial_mode,
      cluster_keys: out.preview.cluster_keys,
    }
  }
  return out
}

function bulletText(result) {
  const bullets = (result && result.summary && result.summary.bullets) || []
  return bullets.length ? bullets.join('；') : ''
}

function figureLines(result) {
  const figures = (result && result.figures) || []
  return figures.map((fig) => fig.file).filter(Boolean)
}

export async function runTool(runtime, args = {}) {
  const action = String(args.action || 'status')
  const dataset = args.dataset || undefined
  const actor = 'agent'
  switch (action) {
    case 'status': {
      const state = await opState(runtime)
      runtime.lastState = { datasets: state.datasets, active: state.active }
      return {
        message: state.datasets.length
          ? `空间工作台就绪：${state.datasets.length} 个数据集，当前 ${state.active || '—'}（Python ${state.engine.python || '?'}）。` +
            ((state.suggestions || [])[0] ? `建议下一步：${state.suggestions[0].title}` : '')
          : `空间工作台已就绪但没有数据集（Python ${state.engine.python || '?'}）。可以 action=search 搜 GEO 空间转录组数据，或 action=open 导入 Visium outs/ 目录 / h5ad，或 action=demo 载入合成组织示例。`,
        ...state,
      }
    }
    case 'search': {
      const result = await opGeoSearch(runtime, { ...args, actor })
      const items = result.items || []
      return {
        message: `${result.total || 0} 条命中，当前返回 ${items.length} 条：` +
          items.slice(0, 8).map((item) => `${item.gse}（${item.n_samples || '?'} 样本${item.has_matrix ? '，有矩阵' : ''}${Array.isArray(item.spatial_platform) && item.spatial_platform.length ? `，${item.spatial_platform.slice(0, 2).join('/')}` : ''}）`).join('；'),
        ...result,
      }
    }
    case 'detail': {
      const result = await opGeoDetail(runtime, { ...args, actor })
      return {
        message: `${result.gse}：${result.n_samples || '?'} 样本，空间平台 ${result.spatial_platform || '—'}，${(result.files || []).length} 个补充文件，建议下载 ${result.suggested || '—'}`,
        ...result,
      }
    }
    case 'download': {
      const result = await opGeoDownload(runtime, { ...args, actor })
      const job = result.job
      // 工具调用等待任务完成（最多 20 分钟）
      const finished = await waitJob(runtime, job.id, 20 * 60 * 1000)
      if (finished.status === 'failed') throw Object.assign(new Error(finished.error || '下载失败'), { code: finished.code || 'internal' })
      return { message: `下载完成 → ${(finished.result && finished.result.dir) || args.dest_dir}`, ...finished }
    }
    case 'open': {
      const result = await opOpen(runtime, { ...args, actor })
      const geometry = (result.dataset && result.dataset.geometry) || {}
      const valueKind = (result.preview && result.preview.X_kind) || 'unknown'
      return {
        message: `已导入 ${result.dataset.name}：${result.dataset.n_spots} spot × ${result.dataset.n_genes} 基因（${result.dataset.source_format || geometry.mode || 'auto'}，空间模式 ${geometry.mode || '—'}，坐标 ${geometry.coords_kind || '—'}，值类型 ${valueKind}）。建议下一步：run step=qc。`,
        ...compactResult(result),
      }
    }
    case 'demo': {
      const result = await opDemo(runtime, { actor })
      return { message: `已载入内置示例：${result.dataset.n_spots} spot × ${result.dataset.n_genes} 基因（合成组织，5 个空间域 / 2 切片 / 2 分组）。建议下一步：run step=qc。`, ...compactResult(result) }
    }
    case 'inspect': {
      const result = await opInspect(runtime, args)
      if (!result.ok) throw Object.assign(new Error(result.error), { code: result.code })
      return { message: `格式 ${result.format}${result.n_spots ? `，${result.n_spots} spot × ${result.n_genes} 基因` : ''}${result.geometry ? `，几何 ${result.geometry}` : ''}${result.candidates ? `，${result.candidates.length} 个候选` : ''}`, ...result }
    }
    case 'preview': {
      const result = await opPreview(runtime, { dataset })
      return {
        message: `${result.n_spots} spot × ${result.n_genes} 基因；空间模式 ${result.spatial_mode || '—'}（邻接半径 ${result.neighbor_radius == null ? '未估计' : Math.round(result.neighbor_radius)}）；分组列 ${result.group_key || '未设置'}；值类型 ${result.X_kind || 'unknown'}；已跑步骤 ${(result.steps || []).join('/') || '无'}。`,
        ...result,
      }
    }
    case 'samples': {
      const result = await opSamples(runtime, { ...args, dataset, actor })
      return { message: `切片/样本表已更新（${(result.sample_info || []).length} 个 spot${result.group_column ? `，分组列 ${result.group_column}` : ''}）`, ...result }
    }
    case 'run': {
      const result = await opRun(runtime, { ...args, dataset, actor })
      const figs = figureLines(result)
      return {
        message: `${result.summary.title}：${bulletText(result)}${figs.length ? `\n图：${figs.join('、')}` : ''}`,
        ...compactResult(result),
      }
    }
    case 'plot': {
      const result = await opPlot(runtime, { ...args, dataset, actor })
      return { message: `已出图：${result.file}`, ...result }
    }
    case 'analyses': {
      const result = await opAnalyses(runtime, { dataset })
      const available = (result.analyses || []).filter((item) => item.available)
      return {
        message: `${(result.analyses || []).length} 个分析（可用 ${available.length} 个），分组：${(result.groups || []).join('、')}。可用：${available.slice(0, 12).map((a) => a.id).join('、')}`,
        ...result,
      }
    }
    case 'analysis': {
      const result = await opAnalysis(runtime, { ...args, dataset, actor })
      const figs = figureLines(result)
      return {
        message: `${result.title || result.analysis_id}（${result.method || ''}）：${bulletText(result)}${figs.length ? `\n图：${figs.join('、')}` : ''}`,
        ...compactResult(result),
      }
    }
    case 'export': {
      const result = await opExport(runtime, { ...args, dataset, actor })
      return { message: `已导出 ${result.what} → ${result.file}（${Math.round((result.bytes || 0) / 1024)} KB）`, ...result }
    }
    case 'audit': {
      const result = await opAudit(runtime, { dataset })
      const bad = (result.checks || []).filter((c) => c.level !== 'ok')
      return {
        message: `自检得分 ${result.score}/100（${result.level}）` + (bad.length ? `；${bad.length} 项待改进：${bad.slice(0, 6).map((c) => c.title).join('、')}` : ''),
        ...result,
      }
    }
    case 'suggest': {
      const result = await opSuggest(runtime, { dataset })
      return { message: (result.suggestions || []).map((s) => s.title).join('；') || '暂无建议', ...result }
    }
    case 'recipe': {
      // 没给具体 id：只列预设流程清单（不执行）；给了 id：真正跑 opRecipeRun（与 UI 的 /recipe/run 同一实现）。
      if (!args.recipe && !args.id) {
        const result = await opRecipe(runtime, args)
        return { message: (result.recipes || []).map((r) => `${r.id}（${r.title}）`).join('；'), ...result }
      }
      const result = await opRecipeRun(runtime, { ...args, dataset, actor })
      const lines = (result.results || []).map((item) => {
        const bullets = ((item.summary && item.summary.bullets) || []).slice(0, 2).join('；')
        const detail = bullets || item.error || ''
        return `${item.step} ${item.ok ? '✓' : '✗'}${detail ? `：${detail}` : ''}`
      })
      let next = ''
      try {
        const suggestion = await opSuggest(runtime, { dataset })
        const first = (suggestion.suggestions || [])[0]
        if (first) next = `建议下一步：${first.title}`
      } catch { /* 没数据集/引擎忙时不给建议 */ }
      const interrupted = result.done < result.total
      return {
        message: `预设流程 ${result.recipe} 已完成 ${result.done}/${result.total} 步：\n${lines.join('\n')}` +
          (interrupted ? '\n有步骤失败，已停在失败处；按上面的原因修好后重试未完成步骤。' : '') +
          (next ? `\n${next}` : ''),
        ...result,
      }
    }
    case 'refs': {
      const result = await opRefs(runtime, args)
      if (result.info) {
        const info = result.info
        return {
          message: `${path.basename(String(info.path || ''))}：${info.n_obs || '?'} 细胞 × ${info.n_vars || '?'} 基因；细胞类型列 ${info.cell_type_col || '未识别'}${Array.isArray(info.cell_types) && info.cell_types.length ? `（${info.cell_types.slice(0, 6).join('/')}）` : ''}`,
          ...result,
        }
      }
      const names = (result.refs || []).map((item) => item.path || item.name)
      const datasets = (result.datasets || []).map((item) => `ds:${item.id}（${item.name}，${item.n_obs || '?'} 细胞）`)
      const all = [...names.slice(0, 5), ...datasets.slice(0, 5)]
      return { message: all.length ? `可用单细胞参考：${all.join('、')}` : '没有找到可用的单细胞参考（可上传 h5ad，或先打开 dsh-singlecell 数据集后用 action=ref_attach + sc_dataset 拉取）', ...result }
    }
    case 'ref_attach': {
      const result = await opRefAttach(runtime, { ...args, dataset, actor })
      const ref = result.ref || {}
      return { message: `已挂载单细胞参考（${result.source}）：${ref.cell_types ? `${ref.cell_types.length} 个细胞类型，` : ''}类型列 ${ref.cell_type_col || '未指定'}；现在可以跑 run step=deconv 或分析 deconv_nnls。`, ...result }
    }
    case 'formats': {
      const result = await opFormats(runtime)
      return { message: `支持 ${(result.formats || []).length} 种空间数据格式：${(result.formats || []).map((f) => f.id).join('、')}`, ...result }
    }
    case 'close': {
      const result = await opClose(runtime, { ...args, dataset, actor })
      return { message: result.closed ? `已关闭数据集 ${args.dataset || ''}。` : '没有可关闭的数据集。', ...result }
    }
    case 'reset': {
      const result = await opReset(runtime, { actor })
      return { message: '空间工作台已清空。', ...result }
    }
    default:
      throw Object.assign(new Error(`未知 action "${action}"，可用：status/search/detail/download/open/demo/inspect/preview/run/plot/analyses/analysis/export/audit/suggest/recipe/samples/refs/ref_attach/formats/reset/close`), { code: 'bad_request' })
  }
}

async function waitJob(runtime, id, timeoutMs) {
  const started = Date.now()
  for (;;) {
    const view = jobView(runtime, id).job
    if (view.status !== 'running') return view
    if (Date.now() - started > timeoutMs) {
      throw Object.assign(new Error('下载任务超时，请稍后用 /jobs 查询'), { code: 'timeout' })
    }
    await new Promise((resolve) => setTimeout(resolve, 800))
  }
}

function buildToolDefinition(hostRuntime) {
  const def = {
    name: 'spatial',
    description:
      '空间转录组分析工作台（本机 Host 插件，带可视化页面与 Python 计算引擎）。用于：检索并下载 GEO 空间转录组数据（Visium / Stereo-seq / Slide-seq / MERFISH 等）、导入 Visium outs/ 目录 / h5ad / Stereo-seq BIN / 表达矩阵+坐标表、跑空间流程（空间质控/过滤/标准化/高变基因/PCA/邻接图/聚类/二维嵌入/空间高变基因 SVG/空间域/反卷积/区域注释）、做空间高级分析（Moran\'s I / Geary\'s C、邻域富集、共现、Ripley\'s K、Getis-Ord 热点、空间梯度、MRF 空间域、边界强度、到区域距离、细胞邻域 niche、NNLS/DWLS 反卷积、空间配体-受体、区域差异、多切片对齐、空间轨迹）与出图/导出。' +
      'action: status=查看工作台状态与建议下一步; search=搜索 GEO 空间转录组数据集（term 关键词，中文或英文；可带 platform/species/study_type/min_samples/has_matrix/limit/page/sort）; detail=GSE 详情与可下载文件列表; download=下载到指定目录（gse+dest_dir；下载完再用 action=open 导入）; ' +
      'open=导入本地/URL 数据（path 或 url；auto 识别 visium/visium_h5/h5ad/stereo/slideseq/matrix_text/xlsx/geo_series/archive/rds/rdata；Visium 直接给 outs/ 目录，h5ad 用 obsm[\'spatial\']，矩阵另配坐标表）; demo=载入内置合成组织示例（2400 spot / 5 空间域 / 2 切片）; inspect=只探测格式、维度与几何; preview=查看基因、spot、坐标范围与坐标列; ' +
      'samples=切片/样本表管理（op=list/update/auto_group/import_map/join_obs/derive/set_group/remove；按切片或分组比较用 op=set_group）; ' +
      'run=按顺序跑流程步骤（qc→filter→normalize→hvg→pca→neighbors→cluster→embed→svg→domains→deconv→annotate，step 必填，params 可选）; ' +
      'analyses=列出全部空间分析及其可用性与缺什么前置; analysis=运行空间分析（analysis 传 id，params 传参数）; plot=出图（kind=qc_spatial/qc_violin/spatial_scatter/embedding/cluster_bar/gene_spatial/svg_map/moran_plot/domain_map/deconv_spatial/deconv_bar/nhood_enrich/cooccurrence/ripley/pca_variance/heatmap/dotplot/gene_box/de_volcano/de_heatmap）; ' +
      'export=导出（what=expr_csv/svg_csv/domain_csv/deconv_csv/sample_info_csv/summary_json/protocol/methods/report/script，可用 out_path）; ' +
      'audit=确定性自检（坐标/邻接半径/标准化/域/反卷积/图像逐项打分）; suggest=建议下一步; recipe=预设流程：省略 id 只列清单（standard/spatial/domain/deconv/quick），给 recipe 或 id 则真正执行该流程并返回每步摘要（与 UI 的 recipe/run 同一实现）; refs=列出可用单细胞参考或查看某个 h5ad 的细胞类型列; ref_attach=挂载反卷积参考 h5ad（ref_path，或从 dsh-singlecell 工作台拉取 sc_dataset）; formats=支持的空间数据格式; close=关闭数据集; reset=清空。' +
      '用户提到空间转录组/空间转录组学/spatial transcriptomics、Visium/Stereo-seq/Slide-seq/10x Visium、空间高变基因/SVG/Moran\'s I、空间域/空间聚类、邻域/共现/Ripley、反卷积/细胞类型比例、配体-受体、GEO 空间数据时优先用本工具驱动这个工作台（与页面按钮同一套实现），不要自己写脚本重跑。' +
      '给结论前建议先 audit；不确定状态先 status。图的绝对路径可直接用 read_image 查看。',
    parameters: {
      action: {
        type: 'string',
        enum: ['status', 'search', 'detail', 'download', 'open', 'demo', 'inspect', 'preview', 'run', 'plot', 'analyses', 'analysis', 'export', 'audit', 'suggest', 'recipe', 'samples', 'refs', 'ref_attach', 'formats', 'reset', 'close'],
        description: '要执行的操作',
      },
      term: { type: 'string', description: 'GEO 检索关键词（中文或英文，如「空间转录组 乳腺癌 Visium」）' },
      gse: { type: 'string', description: 'GSE 编号（如 GSE111672）' },
      dest_dir: { type: 'string', description: '下载目录绝对路径（用户可指定；默认 workdir/downloads/GEO）' },
      path: { type: 'string', description: '本地文件/目录绝对路径（action=open/inspect；Visium 给 outs/ 目录；action=refs 时给参考 h5ad）' },
      url: { type: 'string', description: 'http(s)/ftp 直链（GEO FTP 等）' },
      // 参数名不能叫 `format`（JSON Schema 保留字），用 data_format
      data_format: { type: 'string', enum: DATA_FORMATS, description: '强制数据格式（默认 auto）' },
      dataset: { type: 'string', description: '数据集 id（如 ds1）；省略则用当前活动数据集' },
      step: { type: 'string', enum: STEP_IDS, description: '流程步骤（action=run 必填）' },
      params: { type: 'object', description: '步骤/分析参数（各步骤默认值见 steps 目录；如 qc: min_counts/min_genes/max_pct_mt/group_col；filter: min_counts/min_genes/min_cells；normalize: target_sum/log/method；cluster: resolution/method；svg: n_top/method/permutations；domains: n_domains/smooth/iterations；deconv: ref/method/cell_type_col/min_genes；反卷积可视化的 color_by/gene/point_size）', additionalProperties: true },
      analysis: { type: 'string', description: 'action=analysis 的分析 id（如 morans_svg/nhood_enrichment/domains_mrf/deconv_nnls/ligrec_spatial/de_domain）' },
      kind: { type: 'string', enum: PLOT_KINDS, description: 'action=plot 的图类型' },
      what: { type: 'string', enum: EXPORT_KINDS, description: 'action=export 的导出类型' },
      op: { type: 'string', description: 'action=samples 的子操作：list/update/auto_group/import_map/join_obs/derive/set_group/remove' },
      ref_path: { type: 'string', description: 'action=refs/ref_attach 的参考 h5ad 绝对路径' },
      sc_dataset: { type: 'string', description: 'action=ref_attach 从 dsh-singlecell 工作台拉取的数据集 id（如 ds1）' },
      cell_type_col: { type: 'string', description: 'action=ref_attach 参考 h5ad 里表示细胞类型的 obs 列名' },
      out_path: { type: 'string', description: 'action=export 的目标文件绝对路径（可选）' },
      recipe: { type: 'string', description: 'action=recipe 的预设 id：standard / spatial / domain / deconv / quick（给了 id 会直接执行该流程并返回每步摘要）' },
      source: { type: 'string', description: 'action=samples op=derive 的取值来源（正则或 obs 列名）' },
      platform: { type: 'string', description: 'GEO 空间平台筛选：all / visium / stereo / slideseq / 10x（默认 all）' },
      species: { type: 'string', description: '物种筛选：human / mouse / all（默认 all）' },
      study_type: { type: 'string', description: '研究类型：high throughput sequencing / expression profiling by array / all' },
      min_samples: { type: 'number', description: '最少样本数（默认 0）' },
      has_matrix: { type: 'boolean', description: '只看有表达矩阵数据的系列' },
      limit: { type: 'number', description: '每页条数（默认 20）' },
      page: { type: 'number', description: '页码（从 0 开始）' },
      sort: { type: 'string', description: '排序：relevance / date / samples（默认 relevance）' },
      overwrite: { type: 'boolean', description: 'action=download 时覆盖已存在的同名文件' },
    },
    output: {
      schema: { type: 'object' },
      render: (_args, value) => {
        const text = value && typeof value.message === 'string' ? value.message : JSON.stringify(value)
        return [{ type: 'text', text: String(text) }]
      },
    },
    execute: async (args) => {
      const fresh = await import(new URL(`./impl.mjs?rev=${Date.now()}-${Math.random().toString(36).slice(2)}`, import.meta.url).href)
      return sanitizeValue(await fresh.runTool(hostRuntime, args || {}))
    },
  }
  // 防呆：参数名撞 JSON Schema 保留关键字会让 API 判整个 tools 数组非法
  const reservedKeywords = new Set(['type', 'enum', 'const', 'multipleOf', 'maximum', 'exclusiveMaximum', 'minimum', 'exclusiveMinimum', 'maxLength', 'minLength', 'pattern', 'maxItems', 'minItems', 'uniqueItems', 'maxContains', 'minContains', 'maxProperties', 'minProperties', 'required', 'dependentRequired', 'properties', 'patternProperties', 'additionalProperties', 'propertyNames', 'items', 'prefixItems', 'additionalItems', 'contains', 'title', 'description', 'default', 'deprecated', 'readOnly', 'writeOnly', 'examples', 'format', 'contentEncoding', 'contentMediaType', 'contentSchema', 'if', 'then', 'else', 'allOf', 'anyOf', 'oneOf', 'not'])
  const clash = Object.keys(def.parameters).filter((key) => reservedKeywords.has(key))
  if (clash.length) throw new Error(`[dsh-spatial] 参数名 ${clash.join('、')} 是 JSON Schema 保留关键字，请改名`)
  return def
}

const TOOL_SCHEMA_VERSION = 3

export function ensureToolRegistered(runtime) {
  if (runtime.toolRegistered && runtime.toolSchemaVersion === TOOL_SCHEMA_VERSION) return runtime.toolRegistered
  try {
    const ctx = runtime.ctx
    const tools = ctx && typeof ctx.get === 'function' ? ctx.get('tools') : ctx && ctx.tools
    if (!tools || typeof tools.register !== 'function') {
      runtime.toolRegistered = false
      return false
    }
    if (runtime.toolRegistered) {
      if (typeof runtime.toolRegistered !== 'function') return runtime.toolRegistered
      try { runtime.toolRegistered() } catch (error) {
        try { console.warn('[dsh-spatial] old tool dispose failed:', error && error.message ? error.message : error) } catch { /* ignore */ }
      }
      runtime.toolRegistered = false
    }
    const dispose = tools.register(buildToolDefinition(runtime))
    runtime.toolRegistered = typeof dispose === 'function' ? dispose : true
    runtime.toolSchemaVersion = TOOL_SCHEMA_VERSION
    return true
  } catch (error) {
    runtime.toolRegistered = false
    try { console.warn('[dsh-spatial] agent tool unavailable:', error && error.message ? error.message : error) } catch { /* ignore */ }
    return false
  }
}

// ---------------------------------------------------------------------------
// Prompt section
// ---------------------------------------------------------------------------

export function renderPromptSection(runtime) {
  const state = runtime.lastState
  if (!state || !state.datasets || state.datasets.length === 0) return ''
  const active = state.datasets.find((item) => item.id === state.active)
  const lines = ['空间转录组分析工作台（Host 插件 dsh-spatial）已在使用中。']
  if (active) {
    const done = (active.steps || []).length ? active.steps.join('、') : '无'
    const mode = (active.geometry && active.geometry.mode) || '未识别'
    lines.push(`当前数据集：${active.name || active.id}（${active.n_spots} spot × ${active.n_genes} 基因，空间模式 ${mode}），已完成步骤：${done}。`)
  }
  lines.push(
    '用户提到空间转录组/空间转录组学、GEO 空间数据、Visium/Stereo-seq/Slide-seq 导入、空间质控/聚类、空间高变基因（SVG）/Moran\'s I、空间域、邻域/共现/Ripley、反卷积/细胞类型比例、配体-受体、区域差异时，' +
      '优先调用 spatial 工具操作这个工作台（与页面按钮同一套实现），不要另写脚本重跑；图可用 read_image 查看返回的 PNG 路径。',
  )
  return lines.join('\n')
}

export function ensurePromptSection(runtime) {
  if (runtime.promptRegistered) return runtime.promptRegistered
  try {
    const ctx = runtime.ctx
    const systemPrompt = ctx && typeof ctx.get === 'function' ? ctx.get('systemPrompt') : ctx && ctx.systemPrompt
    if (!systemPrompt || typeof systemPrompt.section !== 'function') {
      runtime.promptRegistered = false
      return false
    }
    const dispose = systemPrompt.section({
      name: 'dsh-spatial-workbench',
      order: 770,
      text: () => renderPromptSection(runtime),
    })
    runtime.promptRegistered = typeof dispose === 'function' ? dispose : true
    return true
  } catch (error) {
    runtime.promptRegistered = false
    try { console.warn('[dsh-spatial] prompt section unavailable:', error && error.message ? error.message : error) } catch { /* ignore */ }
    return false
  }
}

// ---------------------------------------------------------------------------
// 上传
// ---------------------------------------------------------------------------

async function handleUpload(runtime, req, res) {
  const filename = sanitizeName(req.headers['x-dsh-filename'] || 'upload.csv')
  const engine = await ensureEngine(runtime)
  const dest = path.join(engine.dirs.uploads, `${Date.now().toString(36)}-${randomUUID().slice(0, 8)}-${filename}`)
  let size = 0
  await pipeline(
    req,
    async function* count(source) {
      for await (const chunk of source) {
        size += chunk.length
        if (size > MAX_UPLOAD_BYTES) throw Object.assign(new Error('文件超过 8 GiB 上限'), { code: 'bad_request' })
        yield chunk
      }
    },
    createWriteStream(dest),
  )
  const result = await opOpen(runtime, { path: dest, name: filename, actor: 'user' })
  return { ok: true, ...result }
}

// ---------------------------------------------------------------------------
// 路由
// ---------------------------------------------------------------------------

async function serveFigure(runtime, url, res) {
  const id = url.searchParams.get('id')
  const entry = opFigure(runtime, id)
  let stat
  try { stat = statSync(entry.file) } catch {
    throw Object.assign(new Error('图片文件已被清理，请重新生成'), { code: 'not_found' })
  }
  res.writeHead(200, {
    'content-type': 'image/png',
    'cache-control': 'private, max-age=3600',
    'content-length': String(stat.size),
  })
  await pipeline(createReadStream(entry.file), res)
}

export async function handle(runtime, req, res) {
  ensureToolRegistered(runtime)
  ensurePromptSection(runtime)
  if (!hostIsLoopback(req)) {
    return sendJson(res, 403, { ok: false, code: 'forbidden', error: '拒绝：只接受本机（回环地址）请求' })
  }
  const url = new URL(req.url || '/', 'http://127.0.0.1')
  if (!url.pathname.startsWith(PREFIX)) {
    return sendJson(res, 404, { ok: false, code: 'not_found', error: '未知路由' })
  }
  const route = url.pathname.slice(PREFIX.length) || '/'
  const method = (req.method || 'GET').toUpperCase()
  const key = `${method} ${route}`
  try {
    switch (key) {
      case 'GET /health':
        return sendJson(res, 200, await opHealth(runtime))
      case 'GET /state': {
        const state = await opState(runtime)
        runtime.lastState = { datasets: state.datasets, active: state.active }
        return sendJson(res, 200, state)
      }
      case 'GET /formats':
        return sendJson(res, 200, await opFormats(runtime))
      case 'GET /analyses':
        return sendJson(res, 200, await opAnalyses(runtime, { dataset: url.searchParams.get('dataset') || undefined }))
      case 'GET /activity':
        return sendJson(res, 200, activityOf(runtime, Number(url.searchParams.get('since') || 0)))
      case 'GET /jobs':
        return sendJson(res, 200, url.searchParams.get('id') ? opJob(runtime, url.searchParams.get('id')) : opJobs(runtime))
      case 'GET /figure':
        return await serveFigure(runtime, url, res)
      case 'GET /refs': {
        const pathQuery = url.searchParams.get('path') || undefined
        return sendJson(res, 200, await opRefs(runtime, pathQuery ? { path: pathQuery } : {}))
      }
      case 'POST /geo/search':
        return sendJson(res, 200, await opGeoSearch(runtime, await readJson(req)))
      case 'POST /geo/detail':
        return sendJson(res, 200, await opGeoDetail(runtime, await readJson(req)))
      case 'POST /geo/download':
        return sendJson(res, 200, await opGeoDownload(runtime, await readJson(req)))
      case 'POST /dataset/open':
        return sendJson(res, 200, await opOpen(runtime, await readJson(req)))
      case 'POST /dataset/demo':
        return sendJson(res, 200, await opDemo(runtime, await readJson(req)))
      case 'POST /dataset/upload':
        return sendJson(res, 200, await handleUpload(runtime, req, res))
      case 'POST /dataset/select':
        return sendJson(res, 200, await opSelect(runtime, await readJson(req)))
      case 'POST /dataset/close':
        return sendJson(res, 200, await opClose(runtime, await readJson(req)))
      case 'POST /inspect':
        return sendJson(res, 200, await opInspect(runtime, await readJson(req)))
      case 'POST /preview':
        return sendJson(res, 200, await opPreview(runtime, await readJson(req)))
      case 'POST /samples':
        return sendJson(res, 200, await opSamples(runtime, await readJson(req)))
      case 'POST /run':
        return sendJson(res, 200, await opRun(runtime, await readJson(req)))
      case 'POST /plot':
        return sendJson(res, 200, await opPlot(runtime, await readJson(req)))
      case 'POST /analysis':
        return sendJson(res, 200, await opAnalysis(runtime, await readJson(req)))
      case 'POST /ref/attach':
        return sendJson(res, 200, await opRefAttach(runtime, await readJson(req)))
      case 'POST /export':
        return sendJson(res, 200, await opExport(runtime, await readJson(req)))
      case 'POST /audit':
        return sendJson(res, 200, await opAudit(runtime, await readJson(req)))
      case 'POST /suggest':
        return sendJson(res, 200, await opSuggest(runtime, await readJson(req)))
      case 'POST /recipe':
        return sendJson(res, 200, await opRecipe(runtime, await readJson(req)))
      case 'POST /recipe/run':
        return sendJson(res, 200, await opRecipeRun(runtime, await readJson(req)))
      case 'POST /reset':
        return sendJson(res, 200, await opReset(runtime, await readJson(req)))
      case 'POST /engine/restart':
        return sendJson(res, 200, await opEngineRestart(runtime))
      default:
        return sendJson(res, 404, { ok: false, code: 'not_found', error: `未知路由 ${key}` })
    }
  } catch (error) {
    const failure = toFailure(error)
    if (failure && failure.body) return sendJson(res, failure.status, failure.body)
    return sendJson(res, 500, { ok: false, code: 'internal', error: String(error && error.message ? error.message : error) })
  }
}
