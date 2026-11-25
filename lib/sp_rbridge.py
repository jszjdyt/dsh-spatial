"""
dsh-spatial — 通用 Rscript 桥（spacexr/RCTD 等官方 R 实现的调用通道）

用途：把 Python 侧数据打包成 bundle，调
`Rscript --vanilla lib/<script>.R --in <bundle> --out <out> --params <json>`，
再把 R 写出的 result.json / tables/*.csv 读回来。

约定（与 dsh-metabolomics/lib/mb_rbridge.py 同一套设计；两个插件各自独立，这里不 import 它）：
  * bundle 目录（本文件按 RCTD 的输入形态裁剪，见 build_rctd_bundle）：
      ref_counts.mtx | ref_counts.tsv.gz    参考 counts（gene × cell）
      ref_types.tsv                         cell, cell_type 两列
      spot_counts.mtx | spot_counts.tsv.gz  空间 counts（gene × spot）
      coords.tsv                            spot, x, y 三列
      genes.txt                             mtx 时的行名
      params.json                           R 侧参数
  * R 脚本必须写 <out>/result.json，表写 <out>/tables/<name>.csv（首列=行名，Python 侧按
    index_col=0 读回；spot/细胞名不会因为 R 的 check.names 被改写）
  * 退出码：0=成功；3=缺前置（转 prerequisite）；4=运行期程序错误（转 runtime，**不得**据此
    改换方法）；5=缺依赖包（转 unsupported，中文原因）；其它=失败（stderr 尾部作消息）
  * stdout 一律只进返回值，不进协议：桥里没有任何 print()，R 的 stdout/stderr 都被 capture
  * 超时默认 300s、上限 600s（RCTD full 模式实测：2400 spot × 4 类型 ≈ 35s，留足余量）

本文件只依赖标准库 + numpy/pandas + sp_common。禁止新增第三方依赖。
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import threading
import time

import numpy as np
import pandas as pd

import sp_common as sc

#: Rscript 候选（顺序即优先级；可用环境变量指定）
R_CANDIDATES = [
    os.environ.get("DSH_SPATIAL_RSCRIPT") or "",
    os.environ.get("DS_SP_RSCRIPT") or "",
    "/Library/Frameworks/R.framework/Resources/bin/Rscript",
    "/usr/local/bin/Rscript",
    "/opt/homebrew/bin/Rscript",
    "/usr/bin/Rscript",
]

#: 探测时检查的 R 包（spacexr = RCTD 的官方实现，本桥唯一必需的包）
R_PACKAGES = ["spacexr"]

#: 安装提示（与 INSTALL.md 一致）
INSTALL_HINT = (
    "按 INSTALL.md 安装："
    "Rscript -e 'remotes::install_github(\"dmcable/spacexr\")'"
    "（Bioconductor 主站不可达时用 GitHub），装完 Rscript -e 'requireNamespace(\"spacexr\")' 验证"
)

#: 单次 R 调用允许的最长等待（上限 600s）
MAX_TIMEOUT = 600
DEFAULT_TIMEOUT = 300
CACHE_SECONDS = 300

#: 元素数超过该阈值时用 MatrixMarket（大矩阵写 gz TSV 又慢又占内存）
MTX_ELEMENT_THRESHOLD = 5_000_000

_LOCK = threading.Lock()
_STATE: dict = {
    "checked_at": 0.0,
    "available": False,
    "script": None,
    "version": None,
    "packages": {},
    "spacexr_version": None,
    "error": None,
}


class RUnavailable(sc.EngineError):
    """R / spacexr 不可用：调用方可以据此给 unsupported（其它错误不得回退方法）。"""

    def __init__(self, message: str, **details):
        super().__init__("unsupported", message, **details)


# ---------------------------------------------------------------------------
# 探测
# ---------------------------------------------------------------------------


def candidate_scripts() -> list[str]:
    out: list[str] = []
    for cand in R_CANDIDATES:
        if cand and os.path.exists(cand):
            out.append(cand)
    which = shutil.which("Rscript")
    if which:
        out.append(which)
    seen: set[str] = set()
    uniq: list[str] = []
    for item in out:
        if item in seen:
            continue
        seen.add(item)
        uniq.append(item)
    return uniq


def _run(cmd: list[str], timeout: int) -> subprocess.CompletedProcess:
    # R 的 stdout/stderr 可能不是 UTF-8（非 UTF-8 locale 下的中文消息）：errors="replace"
    # 保证桥本身不会因为解码失败而崩；真实错误信息以 R 写出的 error.json 为准。
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                          env=dict(os.environ), encoding="utf-8", errors="replace")


def probe(force: bool = False) -> dict:
    """探测 Rscript 与 spacexr，结果缓存 5 分钟。永不抛异常。"""
    now = time.time()
    with _LOCK:
        if not force and _STATE["checked_at"] and now - _STATE["checked_at"] < CACHE_SECONDS:
            return dict(_STATE)
        state = {
            "checked_at": now,
            "available": False,
            "script": None,
            "version": None,
            "packages": {},
            "spacexr_version": None,
            "error": None,
        }
        for script in candidate_scripts():
            try:
                res = _run([script, "--version"], timeout=20)
                text = (res.stdout or "") + (res.stderr or "")
                if "Rscript" not in text and "R version" not in text:
                    continue
                state["script"] = script
                state["available"] = True
                raw = text.strip()
                match = re.search(r"version\s+([0-9]+\.[0-9]+\.[0-9]+)", raw)
                state["version"] = match.group(1) if match else (raw.splitlines()[0] if raw else "R")
                break
            except Exception as error:  # noqa: BLE001 - 探测失败只记录
                state["error"] = str(error)
                continue
        if state["available"]:
            code = (
                "cat(paste(vapply(c(%s), function(p) requireNamespace(p, quietly=TRUE), logical(1)),"
                " collapse=','))" % (",".join('"%s"' % p for p in R_PACKAGES))
            )
            try:
                res = _run([state["script"], "--vanilla", "-e", code], timeout=120)
                flags = (res.stdout or "").strip().split(",")
                if len(flags) == len(R_PACKAGES):
                    state["packages"] = {p: flags[i].strip() == "TRUE"
                                         for i, p in enumerate(R_PACKAGES)}
                else:
                    state["packages"] = {}
            except Exception as error:  # noqa: BLE001
                state["packages"] = {}
                state["error"] = str(error)
            if (state["packages"] or {}).get("spacexr"):
                try:
                    res = _run([state["script"], "--vanilla", "-e",
                                'cat(as.character(utils::packageVersion("spacexr")))'], timeout=60)
                    text = (res.stdout or "").strip()
                    state["spacexr_version"] = text.splitlines()[-1].strip() if text else None
                except Exception:  # noqa: BLE001
                    state["spacexr_version"] = None
        _STATE.update(state)
        return dict(_STATE)


def r_available() -> bool:
    return bool(probe().get("available"))


def has_pkg(name: str) -> bool:
    return bool((probe().get("packages") or {}).get(name))


def rctd_available() -> bool:
    """RCTD 是否可用（Rscript + spacexr）。"""
    return r_available() and has_pkg("spacexr")


def spacexr_version() -> str | None:
    return probe().get("spacexr_version")


def r_info() -> dict:
    state = probe()
    return {
        "available": bool(state.get("available")),
        "version": state.get("version"),
        "script": state.get("script"),
        "packages": state.get("packages") or {},
        "spacexr_version": state.get("spacexr_version"),
        "error": state.get("error"),
    }


def unavailable_reason() -> str:
    """R / spacexr 不可用的中文原因（用于 unsupported 消息与 method 说明）。"""
    state = probe()
    if not state.get("available"):
        return ("本机没有可用的 Rscript：RCTD（spacexr）不可用（已尝试 "
                + "、".join(candidate_scripts() or ["PATH"]) + "）")
    if not (state.get("packages") or {}).get("spacexr"):
        return "R 已找到（" + str(state.get("script") or "Rscript") + "）但 spacexr 包未安装"
    return ""


# ---------------------------------------------------------------------------
# bundle 写出
# ---------------------------------------------------------------------------


def _as_counts(mat) -> np.ndarray:
    """把计数矩阵转成「整数就用 int64、否则 float64」的数组（写文件前的最小规整）。"""
    arr = np.asarray(mat)
    if arr.dtype == object:
        arr = arr.astype(np.float64)
    arr = np.asarray(arr, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        return arr
    if np.all(arr >= 0) and np.all(np.abs(arr - np.rint(arr)) < 1e-9):
        return np.rint(arr).astype(np.int64)
    return arr


def _write_matrix(bundle_dir: str, base: str, mat, row_names, col_names, fmt: str) -> str:
    """写 gene × sample 矩阵：fmt=mtx → MatrixMarket + genes.txt；fmt=tsv → .tsv.gz。"""
    arr = np.asarray(mat)
    if len(row_names) != arr.shape[0] or len(col_names) != arr.shape[1]:
        raise sc.EngineError("internal", f"run_r: {base} 矩阵形状与行/列名数量不一致")
    if not all(str(name).strip() for name in list(row_names)[:200]):
        raise sc.EngineError("internal", f"run_r: {base} 的行名里有空值")
    if fmt == "mtx":
        from scipy import sparse
        from scipy.io import mmwrite

        counts = _as_counts(arr)
        field = "integer" if np.issubdtype(counts.dtype, np.integer) else "real"
        mmwrite(os.path.join(bundle_dir, f"{base}.mtx"), sparse.coo_matrix(counts), field=field)
        with open(os.path.join(bundle_dir, "genes.txt"), "w", encoding="utf-8") as handle:
            handle.write("\n".join(str(g) for g in row_names) + "\n")
        return os.path.join(bundle_dir, f"{base}.mtx")
    frame = pd.DataFrame(_as_counts(arr), index=[str(v) for v in row_names],
                         columns=[str(v) for v in col_names])
    path = os.path.join(bundle_dir, f"{base}.tsv.gz")
    frame.to_csv(path, sep="\t", index=True, header=True, encoding="utf-8",
                 compression="gzip")
    return path


def build_rctd_bundle(payload: dict, bundle_dir: str) -> str:
    """把 RCTD 输入写成 R 侧 bundle（ref_counts / ref_types / spot_counts / coords / params）。

    payload 键：
      ref_counts   (n_cells, n_genes) 参考原始 counts
      ref_cells    [cell]             参考细胞名
      ref_types    [cell_type]        每个参考细胞的类型
      spot_counts  (n_spots, n_genes) 空间 counts
      spot_names   [spot]             空间 spot 名
      genes        [gene]             两者共用的基因名（顺序即矩阵列序）
      coords       (n_spots, 2) 或 None
      params       dict               mode/max_cores/min_genes/cell_type_col ...
      fmt          "auto"|"mtx"|"tsv"  (默认 auto：大矩阵用 mtx)
    """
    os.makedirs(bundle_dir, exist_ok=True)
    genes = [str(g) for g in (payload.get("genes") or [])]
    ref_cells = [str(c) for c in (payload.get("ref_cells") or [])]
    ref_types = [str(t) for t in (payload.get("ref_types") or [])]
    spots = [str(s) for s in (payload.get("spot_names") or [])]
    ref_counts = payload.get("ref_counts")
    spot_counts = payload.get("spot_counts")
    if ref_counts is None or spot_counts is None:
        raise sc.EngineError("internal", "build_rctd_bundle: 缺少 ref_counts / spot_counts")
    ref_arr = np.asarray(ref_counts, dtype=np.float64)
    spot_arr = np.asarray(spot_counts, dtype=np.float64)
    if ref_arr.ndim != 2 or spot_arr.ndim != 2:
        raise sc.EngineError("internal", "build_rctd_bundle: 计数矩阵必须是二维")
    if len(genes) != ref_arr.shape[1] or len(genes) != spot_arr.shape[1]:
        raise sc.EngineError("internal", "build_rctd_bundle: genes 数量与计数矩阵列数不一致")
    if len(ref_cells) != ref_arr.shape[0] or len(ref_types) != ref_arr.shape[0]:
        raise sc.EngineError("internal", "build_rctd_bundle: 参考细胞名/类型数量与矩阵行数不一致")
    if len(spots) != spot_arr.shape[0]:
        raise sc.EngineError("internal", "build_rctd_bundle: spot 数量与空间矩阵行数不一致")
    if not spots or not ref_cells:
        raise sc.EngineError("internal", "build_rctd_bundle: 参考或空间矩阵为空")
    fmt = str(payload.get("fmt") or "auto").lower()
    if fmt not in ("auto", "mtx", "tsv"):
        fmt = "auto"
    if fmt == "auto":
        fmt = "mtx" if max(ref_arr.size, spot_arr.size) > MTX_ELEMENT_THRESHOLD else "tsv"
    # ref_counts / spot_counts 在 Python 侧是 sample × gene，R 侧要 gene × sample
    _write_matrix(bundle_dir, "ref_counts", ref_arr.T, genes, ref_cells, fmt)
    _write_matrix(bundle_dir, "spot_counts", spot_arr.T, genes, spots, fmt)
    pd.DataFrame({"cell": ref_cells, "cell_type": ref_types}).to_csv(
        os.path.join(bundle_dir, "ref_types.tsv"), sep="\t", index=False, encoding="utf-8")
    coords = payload.get("coords")
    if coords is None:
        coords = np.zeros((len(spots), 2), dtype=np.float64)
    coord_arr = np.asarray(coords, dtype=np.float64)
    if coord_arr.ndim != 2 or coord_arr.shape[0] != len(spots) or coord_arr.shape[1] < 2:
        raise sc.EngineError("internal", "build_rctd_bundle: coords 形状与 spot 数量不一致")
    pd.DataFrame({"spot": spots, "x": coord_arr[:, 0], "y": coord_arr[:, 1]}).to_csv(
        os.path.join(bundle_dir, "coords.tsv"), sep="\t", index=False, encoding="utf-8")
    params = payload.get("params") or {}
    with open(os.path.join(bundle_dir, "params.json"), "w", encoding="utf-8") as handle:
        json.dump(sc.jsonable(params), handle, ensure_ascii=False, indent=2)
    return bundle_dir


# ---------------------------------------------------------------------------
# 结果读写
# ---------------------------------------------------------------------------


def _read_tables(tables_dir: str) -> dict:
    """读 tables/*.csv（首列是行名）；解析失败的表跳过但不丢整个结果。"""
    tables: dict[str, pd.DataFrame] = {}
    if not os.path.isdir(tables_dir):
        return tables
    for name in sorted(os.listdir(tables_dir)):
        if not name.endswith((".csv", ".tsv")):
            continue
        path = os.path.join(tables_dir, name)
        key = os.path.splitext(name)[0]
        try:
            if os.path.getsize(path) == 0:
                continue
            if name.endswith(".tsv"):
                tables[key] = pd.read_csv(path, sep="\t", index_col=0)
            else:
                tables[key] = pd.read_csv(path, index_col=0)
        except Exception:  # noqa: BLE001 - 单表读回失败不影响其它表
            continue
    return tables


def _read_result(outdir: str) -> tuple[dict, dict]:
    result_path = os.path.join(outdir, "result.json")
    payload: dict = {}
    if os.path.exists(result_path):
        try:
            with open(result_path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
        except Exception as error:  # noqa: BLE001
            raise sc.EngineError("internal", f"R 输出 result.json 解析失败：{error}")
    if not isinstance(payload, dict):
        payload = {}
    tables = _read_tables(os.path.join(outdir, "tables"))
    payload.setdefault("warnings", [])
    return payload, tables


def _read_error(outdir: str) -> dict | None:
    path = os.path.join(outdir, "error.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        return payload if isinstance(payload, dict) else None
    except Exception:  # noqa: BLE001
        return None


# ---------------------------------------------------------------------------
# 调用
# ---------------------------------------------------------------------------


def run_r(script_name: str, payload: dict, outdir: str,
          timeout: int = DEFAULT_TIMEOUT) -> dict:
    """调 R 脚本并读回结果。

    返回 {"ok":True,"result":dict,"tables":{name:DataFrame},"warnings":[],
          "ms":int,"stdout":str,"stderr":str,"r_version":str,"spacexr_version":str,"script":str}
    失败抛 EngineError：unsupported（缺 R / 缺包）、prerequisite（数据不满足）、
    runtime（程序错误，调用方不得静默换方法）、timeout、internal。
    """
    state = probe()
    if not state.get("available"):
        raise RUnavailable(unavailable_reason() + "；" + INSTALL_HINT)
    if not (state.get("packages") or {}).get("spacexr"):
        raise RUnavailable(
            "R 包 spacexr（RCTD）未安装：无法调用 spacexr::create.RCTD / run.RCTD 官方实现；"
            + INSTALL_HINT)
    lib_dir = os.path.dirname(os.path.abspath(__file__))
    script = os.path.join(lib_dir, script_name)
    if not os.path.exists(script):
        raise sc.EngineError("unsupported", f"R 脚本不存在：{script_name}")
    timeout = int(max(1, min(int(timeout or DEFAULT_TIMEOUT), MAX_TIMEOUT)))
    outdir = os.path.abspath(outdir)
    bundle = os.path.join(outdir, "r_bundle")
    r_out = os.path.join(outdir, "r_out")
    # 每次调用前清掉上一轮的 bundle / 输出，避免残留表被当成新结果读回
    for path in (bundle, r_out):
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
        os.makedirs(path, exist_ok=True)
    build_rctd_bundle(payload, bundle)
    params_json = json.dumps(sc.jsonable(payload.get("params") or {}), ensure_ascii=False)
    cmd = [state["script"], "--vanilla", script, "--in", bundle, "--out", r_out,
           "--params", params_json]
    t0 = time.time()
    try:
        res = _run(cmd, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise sc.EngineError(
            "timeout",
            f"R 计算超时（>{timeout}s）：{script_name}（RCTD 在 spot 很多时较慢；"
            "可先用 method=nnls 拿到快速结果，或减少 spot 数）")
    except Exception as error:  # noqa: BLE001
        raise sc.EngineError("internal", f"R 调用失败：{error}")
    elapsed = int((time.time() - t0) * 1000)
    stderr_text = (res.stderr or "")
    stdout_text = (res.stdout or "")
    stderr_tail = stderr_text[-2000:]
    error_json = _read_error(r_out)
    if res.returncode != 0:
        message = ""
        if isinstance(error_json, dict) and error_json.get("message"):
            message = str(error_json["message"])
        reason = message or stderr_tail.strip() or stdout_text.strip() or f"Rscript 退出码 {res.returncode}"
        if res.returncode == 5 or "MISSING_PACKAGE:" in stderr_text:
            raise RUnavailable(f"R 依赖包缺失：{reason[-800:]}")
        if res.returncode == 3:
            raise sc.EngineError("prerequisite", reason[-800:])
        if res.returncode == 4:
            # 运行期程序错误：调用方不得据此改换反卷积方法（GEO-02 同一原则）
            raise sc.EngineError("runtime", f"R 计算内部错误：{reason[-800:]}")
        raise sc.EngineError("internal", f"R 计算失败（退出码 {res.returncode}）：{reason[-800:]}")
    parsed, tables = _read_result(r_out)
    return {
        "ok": True,
        "result": parsed,
        "tables": tables,
        "warnings": list(parsed.get("warnings") or []),
        "ms": elapsed,
        "stdout": stdout_text[-2000:],
        "stderr": stderr_tail,
        "r_version": state.get("version"),
        "spacexr_version": state.get("spacexr_version"),
        "script": script,
    }


def try_run_r(script_name: str, payload: dict, outdir: str,
              timeout: int = DEFAULT_TIMEOUT) -> dict | None:
    """run_r 的软版本：R / spacexr 不可用时返回 None。

    注意：prerequisite / runtime / timeout **必须上抛** —— 数据不满足要让用户看到真实原因，
    程序错误绝不能被静默换成另一种算法。
    """
    try:
        return run_r(script_name, payload, outdir, timeout=timeout)
    except RUnavailable:
        return None


__all__ = [
    "MAX_TIMEOUT",
    "DEFAULT_TIMEOUT",
    "MTX_ELEMENT_THRESHOLD",
    "INSTALL_HINT",
    "RUnavailable",
    "candidate_scripts",
    "probe",
    "r_available",
    "r_info",
    "has_pkg",
    "rctd_available",
    "spacexr_version",
    "unavailable_reason",
    "build_rctd_bundle",
    "run_r",
    "try_run_r",
]
