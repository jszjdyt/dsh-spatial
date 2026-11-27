/**
 * 工作台口径契约测试（node:test）
 *
 * 口径卡的职责是「把这次分析用了什么设计/检验/校正」说清楚，且**取不到就写 null**。
 * 这里用空间引擎 `sp_state.Dataset.summary()` 的真实形状做夹具钉死：
 * 齐全 → ok_to_cite、缺字段 → warn caveat、版本常量与 package.json 不漂移。
 */
import test from 'node:test'
import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'

import { buildWorkbenchContract, PLUGIN_VERSION } from '../lib/impl.mjs'

const PACKAGE_JSON = JSON.parse(readFileSync(new URL('../package.json', import.meta.url), 'utf8'))

/** 夹具：sp_state.Dataset.summary() 的真实键（+ 这次为口径卡新增的键）。 */
function baseDataset(overrides = {}) {
  return {
    id: 'ds1',
    name: 'demo-spatial',
    path: null,
    source_format: 'visium',
    n_spots: 2400,
    n_genes: 1000,
    species: 'human',
    geometry: {
      mode: 'visium',
      coords_kind: 'visium_hires',
      x_range: [0, 6400],
      y_range: [0, 6400],
      image: null,
      scale: 1,
      neighbor_radius: 120,
    },
    sample_columns: [
      { key: 'sample', kind: 'categorical', n_unique: 4 },
      { key: 'group', kind: 'categorical', n_unique: 2 },
      { key: 'batch', kind: 'categorical', n_unique: 2 },
      { key: 'donor', kind: 'categorical', n_unique: 4 },
    ],
    steps: ['qc', 'filter', 'normalize', 'hvg', 'pca', 'neighbors', 'cluster'],
    current: 'cluster',
    cluster_keys: ['clusters'],
    ref: { path: '/refs/pbmc.h5ad', cell_type_col: 'cell_type', cell_types: ['T cell', 'B cell'] },
    warnings: [],
    notes: [],
    history: [
      { at: '2026-11-27T10:00:00.000Z', step: 'normalize', params: { method: 'lognorm', target_sum: 10000 }, ok: true, ms: 300, summary: 'lognorm' },
      { at: '2026-11-27T10:05:00.000Z', step: 'cluster', params: { method: 'leiden', resolution: 0.6, n_neighbors: 15 }, ok: true, ms: 1200, summary: '8 clusters' },
    ],
    analyses_run: [
      { at: '2026-11-27T10:02:00.000Z', id: 'morans_svg', ok: true, ms: 900, summary: '30 SVG' },
      { at: '2026-11-27T10:03:00.000Z', id: 'domains_mrf', ok: true, ms: 4000, summary: '6 domains' },
    ],
    group_key: 'group',
    groups: ['tumor', 'adjacent'],
    has_counts: true,
    X_kind: 'counts',
    created_at: '2026-11-27T09:00:00.000Z',
    n_history: 2,
    import_fingerprint: '0f1e2d3c4b5a6978',
    slice_key: 'sample',
    n_slices: 4,
    group_counts: { tumor: 2, adjacent: 2 },
    donor_key: 'donor',
    n_donors: 4,
    value_kind: 'counts',
    has_coords: true,
    ...overrides,
  }
}

test('口径 1) 齐全：切片/供体 n / 最近步骤的设计与校正 / 导入指纹都可回查', () => {
  const contract = buildWorkbenchContract(baseDataset())
  assert.equal(contract.schema, 'dsh-workbench.stats_contract/1')
  assert.equal(contract.plugin, 'dsh-spatial')
  assert.equal(contract.plugin_version, PLUGIN_VERSION)
  // 版本常量必须与 package.json 一致（CI：node --test 即校验，防漂移）
  assert.equal(PLUGIN_VERSION, PACKAGE_JSON.version)
  // spot 数是矩阵行数，切片/取样单元数才是生物学 n（不拿 spot 冒充重复）
  assert.equal(contract.dataset.n_spots, 2400)
  assert.equal(contract.dataset.n_genes, 1000)
  assert.equal(contract.dataset.n_samples, 4)
  assert.equal(contract.dataset.n_spots !== contract.dataset.n_samples, true)
  assert.equal(contract.dataset.value_kind, 'counts')
  assert.equal(contract.dataset.data_kind, 'spatial')
  assert.equal(contract.dataset.spatial_mode, 'visium')
  assert.equal(contract.dataset.coords_kind, 'visium_hires')
  assert.equal(contract.dataset.has_counts, true)
  assert.equal(contract.dataset.has_coords, true)
  assert.deepEqual(contract.design.groups, ['tumor', 'adjacent'])
  assert.deepEqual(contract.design.n_per_group, { tumor: 2, adjacent: 2 })
  assert.equal(contract.design.sample_key, 'sample')
  assert.equal(contract.design.donor_key, 'donor')
  assert.equal(contract.design.n_donors, 4)
  assert.equal(contract.design.batches, true)
  // 最近事件按 at 比较 = cluster（10:05），参数原样搬运
  assert.equal(contract.inference.last_step, 'cluster')
  assert.equal(contract.inference.method.method, 'leiden')
  assert.equal(contract.inference.design.resolution, 0.6)
  assert.equal(contract.inference.ok, true)
  assert.equal(contract.inference.ms, 1200)
  assert.deepEqual(contract.recent_analyses.map((item) => item.analysis), ['morans_svg', 'domains_mrf'])
  assert.equal(contract.identity.fingerprint, '0f1e2d3c4b5a6978')
  assert.equal(contract.identity.n_steps_done, 7)
  assert.deepEqual(contract.caveats.filter((item) => item.level === 'warn'), [])
  assert.equal(contract.ok_to_cite, true)
  assert.doesNotThrow(() => JSON.stringify(contract), '契约必须可 JSON 序列化')
})

test('口径 2) 缺分组 / 无切片列 / 尺度未知 / 无 counts / 无坐标 / 最近失败 → 全部进 warn caveats', () => {
  const contract = buildWorkbenchContract(baseDataset({
    group_key: null,
    groups: [],
    group_counts: {},
    slice_key: null,
    n_slices: null,
    donor_key: null,
    n_donors: null,
    value_kind: 'unknown',
    X_kind: 'unknown',
    has_counts: false,
    has_coords: false,
    history: [
      { at: '2026-11-27T10:00:00.000Z', step: 'neighbors', params: { n_neighbors: 15 }, ok: false, ms: 500, summary: 'pca 缺失' },
    ],
    analyses_run: [],
    sample_columns: [],
  }))
  const codes = contract.caveats.map((item) => item.code)
  assert.ok(codes.includes('no_groups'))
  assert.ok(codes.includes('unit_unknown'))
  assert.ok(codes.includes('scale_unknown'))
  assert.ok(codes.includes('counts_missing'))
  assert.ok(codes.includes('coords_missing'))
  assert.ok(codes.includes('last_step_failed'))
  assert.ok(codes.includes('no_sample_info'))
  assert.equal(contract.ok_to_cite, false)

  // 有分组但拿不到每组切片数 → group_n_unknown
  const noCounts = buildWorkbenchContract(baseDataset({ group_counts: {} }))
  assert.ok(noCounts.caveats.some((item) => item.code === 'group_n_unknown' && item.level === 'warn'))
  assert.equal(noCounts.design.n_per_group.tumor, null)
})

test('口径 3) 每组切片 n<2 必须警告；无历史时只给 info（不误报成错误）', () => {
  const lowN = buildWorkbenchContract(baseDataset({ group_counts: { tumor: 1, adjacent: 2 } }))
  assert.ok(lowN.caveats.some((item) => item.code === 'min_group_n_below_2'))
  assert.equal(lowN.ok_to_cite, false)

  const noHistory = buildWorkbenchContract(baseDataset({ history: [], analyses_run: [], n_history: 0 }))
  assert.equal(noHistory.inference, null)
  assert.ok(noHistory.caveats.some((item) => item.code === 'no_history' && item.level === 'info'))
  assert.ok(noHistory.caveats.some((item) => item.code === 'no_analyses' && item.level === 'info'))
  assert.ok(!noHistory.caveats.some((item) => item.code === 'no_history' && item.level === 'warn'))
})

test('口径 4) 空数据集不抛错，字段为 null/空而不是编造', () => {
  const contract = buildWorkbenchContract({})
  assert.equal(contract.dataset.n_samples, null)
  assert.equal(contract.dataset.n_spots, null)
  assert.equal(contract.dataset.value_kind, 'unknown')
  assert.equal(contract.dataset.has_counts, null)
  assert.equal(contract.dataset.has_coords, null)
  assert.deepEqual(contract.design.groups, [])
  assert.equal(contract.inference, null)
  assert.equal(contract.identity.fingerprint, '')
  assert.ok(contract.caveats.some((item) => item.code === 'no_groups'))
  assert.equal(buildWorkbenchContract(null).ok_to_cite, false)
})
