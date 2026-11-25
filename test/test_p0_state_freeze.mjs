/**
 * P0-2 / P0-3 回归测试：单步失败不返回 ok:true、切 dataset 不漂移、
 * 同参数重复提交幂等、不同参数产生新 run_id、已提交产物不可变可回查。
 *
 * 用真实 lib/engine.mjs（队列 + run store + withTask）驱动 impl.mjs 的真实 op* 实现，
 * Python 侧换成只讲 JSON-lines 的假 worker，因此不依赖 numpy/pandas。
 *
 * 运行：node test/test_p0_state_freeze.mjs
 */

import { mkdirSync, mkdtempSync, existsSync, readFileSync, readdirSync, rmSync, writeFileSync } from 'node:fs'
import { tmpdir } from 'node:os'
import path from 'node:path'
import { createEngine } from '../lib/engine.mjs'
import * as impl from '../lib/impl.mjs'

const PLUGIN = 'spatial'
const PLOT_KIND = 'de_volcano'
const TMP = mkdtempSync(path.join(tmpdir(), `dsh-p0-state-${PLUGIN}-`))
const workerPath = path.join(TMP, 'fake_worker.py')

const FAKE_WORKER = `#!/usr/bin/env python3
import json, os, sys, time

def arg(name, default=None):
    try:
        return sys.argv[sys.argv.index(name) + 1]
    except Exception:
        return default

WORKDIR = arg('--workdir', '.')

def note(cmd, args):
    try:
        with open(os.path.join(WORKDIR, 'calls.jsonl'), 'a', encoding='utf-8') as handle:
            handle.write(json.dumps({'cmd': cmd, 'args': args}, ensure_ascii=False) + '\\n')
    except Exception:
        pass

def reply(msg_id, result=None, error=None):
    payload = {'id': msg_id, 'ok': error is None}
    if error is None:
        payload['result'] = result if result is not None else {}
    else:
        payload['error'] = {'code': 'internal', 'message': error}
    sys.stdout.write(json.dumps(payload) + '\\n')
    sys.stdout.flush()

def dataset(summary_id, name, n):
    return {'id': summary_id, 'name': name, 'n_samples': n, 'n_genes': 100, 'steps': [], 'path': None,
            'groups': ['case', 'control']}

DATASETS = [dataset('ds1', 'A', 4), dataset('ds2', 'B', 6)]
ACTIVE = 'ds1'
RECIPES = {
    'standard': [{'step': 'qc'}, {'step': 'de'}],
    'bad': [{'step': 'qc'}, {'step': 'boom'}, {'step': 'enrich'}],
    'plain': [{'step': 'qc'}],
}

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    try:
        msg = json.loads(line)
    except Exception:
        continue
    cmd = msg.get('cmd')
    args = msg.get('args') or {}
    note(cmd, args)
    if cmd == 'ping':
        reply(msg.get('id'), {'mode': 'builtin', 'packages': {}, 'layouts': ['pca'], 'platform': 'test'})
    elif cmd == 'list':
        reply(msg.get('id'), {'datasets': DATASETS, 'active': ACTIVE})
    elif cmd == 'state':
        reply(msg.get('id'), {'datasets': DATASETS, 'active': ACTIVE})
    elif cmd == 'select':
        ACTIVE = str(args.get('dataset') or ACTIVE)
        reply(msg.get('id'), {'active': ACTIVE})
    elif cmd == 'recipe':
        rid = str(args.get('id') or '')
        if rid:
            reply(msg.get('id'), {'id': rid, 'title': rid, 'steps': RECIPES.get(rid, RECIPES['standard'])})
        else:
            reply(msg.get('id'), {'recipes': [{'id': k, 'title': k, 'steps': v} for k, v in RECIPES.items()]})
    elif cmd == 'run':
        step = str(args.get('step') or '')
        if step == 'boom':
            reply(msg.get('id'), error='合成失败：de 步骤失败')
            continue
        # 模拟「用户在两步之间切了活动数据集」：qc 跑完后把 active 改成 ds2
        try:
            if step == 'qc' and os.path.exists(os.path.join(WORKDIR, 'switch-dataset')):
                ACTIVE = 'ds2'
        except Exception:
            pass
        reply(msg.get('id'), {'dataset': {'id': args.get('dataset') or ACTIVE},
                              'summary': {'title': step, 'bullets': [step + ' ok']},
                              'ms': 3})
    elif cmd == 'plot':
        kind = str(args.get('kind') or 'plot')
        name = 'f_%s_%s.png' % (kind, str(int(time.time() * 1000))[-6:])
        target = os.path.join(WORKDIR, 'figures', name)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with open(target, 'wb') as handle:
            handle.write(b'\\x89PNG\\r\\n\\x1a\\n fake')
        reply(msg.get('id'), {'id': name, 'file': target, 'filename': name, 'width': 100, 'height': 80,
                              'dataset': {'id': args.get('dataset') or ACTIVE}})
    elif cmd == 'export':
        what = str(args.get('what') or 'file')
        target = os.path.join(WORKDIR, 'exports', what + '.csv')
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with open(target, 'w', encoding='utf-8') as handle:
            handle.write('gene,logfc\\nA,1\\n')
        reply(msg.get('id'), {'what': what, 'file': target, 'bytes': 16})
    elif cmd == 'geo_download':
        dest = str(args.get('dest_dir') or os.path.join(WORKDIR, 'downloads'))
        os.makedirs(dest, exist_ok=True)
        files = []
        for name in (args.get('files') or ['series_matrix.txt']):
            target = os.path.join(dest, name)
            with open(target, 'w', encoding='utf-8') as handle:
                handle.write('ID_REF' + chr(9) + 'GSM1' + chr(10) + '1' + chr(9) + '2' + chr(10))
            files.append({'name': name, 'path': target, 'bytes': 17, 'ok': True})
        reply(msg.get('id'), {'dir': dest, 'gse': args.get('gse'), 'files': files, 'ms': 1})
    elif cmd == 'shutdown':
        reply(msg.get('id'), {'bye': True})
        break
    else:
        reply(msg.get('id'), {'echo': cmd})
`

writeFileSync(workerPath, FAKE_WORKER, { mode: 0o755 })

let passed = 0
let failed = 0
const failures = []

function check(name, condition, detail = '') {
  if (condition) {
    passed += 1
    console.log(`  PASS ${name}${detail ? `  [${detail}]` : ''}`)
  } else {
    failed += 1
    failures.push(`${name}${detail ? `：${detail}` : ''}`)
    console.log(`  FAIL ${name}${detail ? `  [${detail}]` : ''}`)
  }
}

function readCalls(workdir) {
  const file = path.join(workdir, 'calls.jsonl')
  if (!existsSync(file)) return []
  return readFileSync(file, 'utf8').split('\n').filter(Boolean).map((line) => {
    try { return JSON.parse(line) } catch { return { cmd: '?', args: {} } }
  })
}

function stagingLeftovers(workdir) {
  try { return readdirSync(path.join(workdir, 'staging')) } catch { return [] }
}

async function waitJob(runtime, jobId, timeoutMs = 5000) {
  const started = Date.now()
  for (;;) {
    const view = impl.opJob(runtime, jobId).job
    if (view.status !== 'running') return view
    if (Date.now() - started > timeoutMs) throw new Error(`job ${jobId} 超时`)
    await new Promise((resolve) => setTimeout(resolve, 30))
  }
}

async function main() {
  const workdir = path.join(TMP, 'work')
  const engine = createEngine({
    workdir,
    engineScript: workerPath,
    python: process.env.DSH_TEST_PYTHON || 'python3',
  })
  const runtime = {
    engine,
    config: {},
    ctx: null,
    activity: [],
    jobs: new Map(),
  }

  check('0 引擎与 run store 就绪', (await engine.start()).mode === 'builtin' && Boolean(engine.store))

  // -------------------------------------------------------------------------
  // 1. 单步失败 → 整体 ok:false（旧 bug：ok:true）
  // -------------------------------------------------------------------------
  const bad = await impl.opRecipeRun(runtime, { recipe: 'bad' })
  check('1a 单步失败整体 ok:false（旧 bug：ok:true）', bad.ok === false, `ok=${bad.ok}`)
  check('1b status=partial 且失败步骤保留原因', bad.status === 'partial' &&
    bad.results.some((item) => item.step === 'boom' && item.status === 'failed' && /合成失败/.test(item.error || '')),
    JSON.stringify(bad.results))
  check('1c 失败后的步骤标记 skipped', bad.results.some((item) => item.step === 'enrich' && item.status === 'skipped'),
    JSON.stringify(bad.results.map((item) => `${item.step}:${item.status}`)))
  check('1d 失败返回带 run_id 与明确 code', Boolean(bad.run_id) && bad.code === 'step_failed',
    `run=${bad.run_id} code=${bad.code}`)
  const badRun = engine.getRun(bad.run_id)
  check('1e 失败 run 有独立记录与逐步状态', Boolean(badRun) && badRun.steps.length === 3 &&
    badRun.steps.map((s) => `${s.step}:${s.status}`).join(',') === 'qc:ok,boom:failed,enrich:skipped',
    badRun ? JSON.stringify(badRun.steps.map((s) => `${s.step}:${s.status}`)) : 'no-run')
  check('1f 失败 run 不以成功提交且没有产物', Boolean(badRun) && badRun.status !== 'committed' && badRun.artifacts.length === 0,
    badRun ? `status=${badRun.status} artifacts=${badRun.artifacts.length}` : 'no-run')

  // -------------------------------------------------------------------------
  // 2. 切 dataset 不漂移：qc 执行期间 active 被改成 ds2，后续步骤仍用 ds1
  // -------------------------------------------------------------------------
  writeFileSync(path.join(workdir, 'switch-dataset'), '1')
  const callsBefore = readCalls(workdir).length
  const drift = await impl.opRecipeRun(runtime, { recipe: 'standard', dataset: 'ds1' })
  const runCalls = readCalls(workdir).slice(callsBefore).filter((item) => item.cmd === 'run')
  check('2a 任务开始时冻结数据集 id', drift.ok !== false && drift.frozen && drift.frozen.id === 'ds1',
    drift.frozen ? `frozen=${drift.frozen.id}` : 'no-frozen')
  check('2b 每一步都显式带冻结的 dataset（不读活动数据集）',
    runCalls.length === 2 && runCalls.every((item) => item.args && item.args.dataset === 'ds1'),
    JSON.stringify(runCalls.map((item) => item.args && item.args.dataset)))
  const standardRun = engine.listRuns({ op: 'recipe:standard' })[0]
  check('2c run 记录里绑定的是 ds1 的版本', Boolean(standardRun) && standardRun.dataset &&
    standardRun.dataset.id === 'ds1' && /^[0-9a-f]{64}$/.test(standardRun.dataset.revision || ''),
    standardRun && standardRun.dataset ? JSON.stringify(standardRun.dataset) : 'no-run')

  // -------------------------------------------------------------------------
  // 3. 同参数重复提交幂等
  // -------------------------------------------------------------------------
  // 幂等要显式指定同一个数据集：不指定时按「任务开始时的活动数据集」冻结，
  // 而上面的用例已经把活动数据集切成了 ds2（那属于不同数据版本 → 本就该是新 run）。
  const runCountBefore = readCalls(workdir).filter((item) => item.cmd === 'run').length
  const again = await impl.opRecipeRun(runtime, { recipe: 'standard', dataset: 'ds1' })
  const runCountAfter = readCalls(workdir).filter((item) => item.cmd === 'run').length
  check('3a 同参数重复提交复用同一个 run_id', again.run_id === drift.run_id && again.run_reused === true,
    `first=${drift.run_id} second=${again.run_id} reused=${again.run_reused}`)
  check('3b 幂等复用没有再执行步骤', runCountAfter === runCountBefore,
    `before=${runCountBefore} after=${runCountAfter}`)

  // -------------------------------------------------------------------------
  // 4. 不同参数 → 新 run_id
  // -------------------------------------------------------------------------
  const variant = await impl.opRecipeRun(runtime, { recipe: 'standard', dataset: 'ds1', params: { de: { method: 'deseq2' } } })
  check('4a 不同参数产生新 run_id', variant.run_id !== drift.run_id,
    `standard=${drift.run_id} variant=${variant.run_id}`)

  // -------------------------------------------------------------------------
  // 5. 图：写入不可变 run 目录、可回查、同参数复用、不同参数新 run
  // -------------------------------------------------------------------------
  const plot1 = await impl.opPlot(runtime, { kind: PLOT_KIND, params: { lfc: 1 } })
  check('5a 出图返回 run_id 与图片路径', Boolean(plot1.run_id) && Boolean(plot1.figure && plot1.figure.file),
    `run=${plot1.run_id}`)
  const plotRun = engine.getRun(plot1.run_id)
  const figureArtifact = plotRun && plotRun.artifacts.find((item) => item.kind === 'figure')
  check('5b 图被提交进 runs/<run_id>/figures', Boolean(figureArtifact) && figureArtifact.file.startsWith(plotRun.dir) &&
    existsSync(figureArtifact.file), figureArtifact ? figureArtifact.file : 'no-artifact')
  check('5c 产物带 sha256 校验', Boolean(figureArtifact) && /^[0-9a-f]{64}$/.test(figureArtifact.sha256 || ''),
    figureArtifact && figureArtifact.sha256 ? figureArtifact.sha256.slice(0, 12) : 'no-hash')
  check('5d 按文件名可回查已提交产物', (() => {
    const found = engine.resolveArtifact(figureArtifact.name, 'figure')
    return Boolean(found) && found.file === figureArtifact.file
  })(), figureArtifact ? figureArtifact.name : 'n/a')
  const plot2 = await impl.opPlot(runtime, { kind: PLOT_KIND, params: { lfc: 1 } })
  check('5e 同参出图幂等（复用同一 run）', plot2.run_id === plot1.run_id && plot2.run_reused === true,
    `first=${plot1.run_id} second=${plot2.run_id}`)
  const plot3 = await impl.opPlot(runtime, { kind: PLOT_KIND, params: { lfc: 2 } })
  check('5f 不同参数出图产生新 run_id', plot3.run_id !== plot1.run_id,
    `a=${plot1.run_id} b=${plot3.run_id}`)

  // -------------------------------------------------------------------------
  // 6. 导出：linked artifact 登记不复制
  // -------------------------------------------------------------------------
  const exported = await impl.opExport(runtime, { what: 'de_csv' })
  const exportRun = engine.getRun(exported.run_id)
  const exportArtifact = exportRun && exportRun.artifacts.find((item) => item.kind === 'file')
  check('6a 导出登记进 run 记录且有哈希', Boolean(exportArtifact) && /^[0-9a-f]{64}$/.test(exportArtifact.sha256 || ''),
    exportArtifact ? exportArtifact.relPath || 'link' : 'no-artifact')

  // -------------------------------------------------------------------------
  // 7. 失败/非法请求不污染已提交结果，staging 不残留
  // -------------------------------------------------------------------------
  const committedBefore = engine.listRuns({ limit: 100 }).filter((item) => item.status === 'committed').length
  let threw = false
  try {
    await impl.opRecipeRun(runtime, { recipe: 'standard', dataset: 'ds-never' })
  } catch (error) {
    threw = error && error.code === 'not_found'
  }
  const committedAfter = engine.listRuns({ limit: 100 }).filter((item) => item.status === 'committed').length
  check('7a 非法数据集明确失败', threw === true)
  check('7b 失败请求不新增已提交 run', committedAfter === committedBefore,
    `before=${committedBefore} after=${committedAfter}`)
  check('7c staging 无残留目录', stagingLeftovers(workdir).length === 0,
    JSON.stringify(stagingLeftovers(workdir)))
  check('7d 旧已提交产物仍在且可读', existsSync(figureArtifact.file))

  // -------------------------------------------------------------------------
  // 8. 下载：路径已存在 → 按哈希复用（不再访问 NCBI）；overwrite 才重新下载
  // -------------------------------------------------------------------------
  const destDir = path.join(TMP, 'dl')
  mkdirSync(destDir, { recursive: true })
  const localFile = path.join(destDir, 'series_matrix.txt')
  writeFileSync(localFile, 'ID_REF\tGSM1\n1\t2\n')
  engine.acquireFile(localFile, { name: 'series_matrix.txt' })
  const dlCallsBefore = readCalls(workdir).filter((item) => item.cmd === 'geo_download').length
  const dlJob = await impl.opGeoDownload(runtime, { gse: 'GSE000', dest_dir: destDir, files: ['series_matrix.txt'] })
  const dlDone = await waitJob(runtime, dlJob.job.id)
  const dlCallsAfter = readCalls(workdir).filter((item) => item.cmd === 'geo_download').length
  check('8a 已存在文件按哈希复用', dlDone.status === 'done' &&
    (dlDone.result.files || []).every((item) => item.reused === true) && Boolean((dlDone.result.files || [])[0].sha256),
    JSON.stringify(dlDone.result && dlDone.result.files))
  check('8b 复用时不触发真正的下载', dlCallsAfter === dlCallsBefore,
    `before=${dlCallsBefore} after=${dlCallsAfter}`)
  const dlJob2 = await impl.opGeoDownload(runtime, { gse: 'GSE000', dest_dir: destDir, files: ['series_matrix.txt'], overwrite: true })
  const dlDone2 = await waitJob(runtime, dlJob2.job.id)
  const dlCallsAfter2 = readCalls(workdir).filter((item) => item.cmd === 'geo_download').length
  check('8c overwrite=true 才重新下载', dlCallsAfter2 === dlCallsBefore + 1 && dlDone2.status === 'done',
    `calls=${dlCallsAfter2 - dlCallsBefore}`)

  await engine.dispose()
  try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }

  console.log(`\n== P0-2/P0-3 状态回归（${PLUGIN}）：${passed} 通过 / ${failed} 失败 ==`)
  if (failed) {
    for (const item of failures) console.log(`  ✗ ${item}`)
    process.exit(1)
  }
}

main().catch((error) => {
  console.error('测试本身异常：', error)
  try { rmSync(TMP, { recursive: true, force: true }) } catch { /* ignore */ }
  process.exit(1)
})
