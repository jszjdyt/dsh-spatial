#!/usr/bin/env python3
"""dsh-spatial — Host/Python 协议冒烟测试（只依赖标准库）。

验证 DESIGN §2.2 命令表与 sp_registry 目录的连通性：
  1. ping           → ok=true，capabilities.steps 含 12 个空间流程步骤（顺序正确）
  2. formats        → 空间数据格式清单非空
  3. analyses       → 空间分析目录 ≥ 25 项
  4. demo           → 内置合成组织示例（dataset.n_spots > 0）
  5. bogus          → 未知命令返回 error.code = bad_request

用法：
  cd plugins/dsh-spatial && python3 test/host_smoke_test.py

不写工作目录到仓库内：工作目录用 tempfile.mkdtemp()，测试结束自动清理。
"""

from __future__ import annotations

import json
import os
import queue
import shutil
import subprocess
import sys
import tempfile
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                     # plugins/dsh-spatial
ENGINE = os.path.join(ROOT, "lib", "sp_engine.py")
PYTHON = sys.executable or "python3"

EXPECTED_STEPS = [
    "qc", "filter", "normalize", "hvg", "pca", "neighbors",
    "cluster", "embed", "svg", "domains", "deconv", "annotate",
]
MIN_ANALYSES = 25
DEFAULT_TIMEOUT = 180.0
DEMO_TIMEOUT = 600.0


class Worker:
    """sp_engine.py 子进程 + 一行一 JSON 的请求/响应封装。"""

    def __init__(self, workdir: str):
        self.workdir = workdir
        self.proc = subprocess.Popen(
            [PYTHON, "-u", ENGINE, "--workdir", workdir],
            cwd=ROOT,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        self.lines: "queue.Queue[str | None]" = queue.Queue()
        threading.Thread(target=self._pump, name="sp-engine-stdout", daemon=True).start()

    def _pump(self) -> None:
        assert self.proc.stdout is not None
        for line in self.proc.stdout:
            self.lines.put(line)
        self.lines.put(None)

    def call(self, request_id: int, cmd: str, args=None, timeout: float = DEFAULT_TIMEOUT) -> dict:
        payload = {"id": request_id, "cmd": cmd}
        if args is not None:
            payload["args"] = args
        assert self.proc.stdin is not None
        self.proc.stdin.write(json.dumps(payload, ensure_ascii=False) + "\n")
        self.proc.stdin.flush()
        deadline = time.time() + timeout
        while True:
            remaining = deadline - time.time()
            if remaining <= 0:
                raise TimeoutError(f"等待 {cmd} 响应超时（{timeout}s）")
            try:
                line = self.lines.get(timeout=min(remaining, 5.0))
            except queue.Empty:
                if self.proc.poll() is not None and self.lines.empty():
                    raise RuntimeError(f"引擎提前退出（code={self.proc.returncode}）：{self.stderr_tail()}")
                continue
            if line is None:
                raise RuntimeError(f"引擎 stdout 已关闭（code={self.proc.poll()}）：{self.stderr_tail()}")
            line = line.strip()
            if not line:
                continue
            try:
                message = json.loads(line)
            except json.JSONDecodeError:
                raise AssertionError(f"stdout 出现非协议行：{line[:200]}")
            if message.get("event") == "progress":
                continue
            if message.get("id") == request_id:
                return message

    def stderr_tail(self, limit: int = 2000) -> str:
        assert self.proc.stderr is not None
        try:
            return self.proc.stderr.read()[-limit:]
        except Exception:  # noqa: BLE001
            return ""

    def close(self) -> None:
        try:
            if self.proc.poll() is None:
                assert self.proc.stdin is not None
                self.proc.stdin.write(json.dumps({"id": 99, "cmd": "shutdown"}) + "\n")
                self.proc.stdin.flush()
                self.proc.wait(timeout=15)
        except Exception:  # noqa: BLE001
            pass
        finally:
            if self.proc.poll() is None:
                self.proc.kill()
            for stream in (self.proc.stdin, self.proc.stdout, self.proc.stderr):
                try:
                    if stream is not None:
                        stream.close()
                except Exception:  # noqa: BLE001
                    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def unwrap(response: dict, cmd: str) -> dict:
    check(response.get("ok") is True, f"{cmd} 失败：{json.dumps(response.get('error'), ensure_ascii=False)}")
    result = response.get("result")
    check(isinstance(result, dict), f"{cmd} 的 result 不是对象：{type(result).__name__}")
    return result


def main() -> int:
    workdir = tempfile.mkdtemp(prefix="dsh-spatial-smoke-")
    worker = Worker(workdir)
    summary: list[str] = []
    try:
        # 1) ping：12 个空间流程步骤，且顺序与 DESIGN §3 一致
        ping = unwrap(worker.call(1, "ping"), "ping")
        capabilities = ping.get("capabilities") or {}
        raw_steps = capabilities.get("steps") or []
        step_ids = [item["id"] if isinstance(item, dict) else str(item) for item in raw_steps]
        check(step_ids == EXPECTED_STEPS,
              f"ping.capabilities.steps 与 DESIGN §3 不一致：{step_ids}")
        check(ping.get("geo") is True, "ping.geo 应为 true")
        check(bool(ping.get("python")), "ping.python 为空")
        summary.append(f"ping: python {ping.get('python')} / steps {len(step_ids)} / "
                       f"analyses {len(capabilities.get('analyses') or [])} / "
                       f"plots {len(capabilities.get('plots') or [])}")

        # 2) formats：空间数据格式清单非空（visium/h5ad/stereo…）
        formats = unwrap(worker.call(2, "formats"), "formats")
        format_items = formats.get("formats") or []
        check(isinstance(format_items, list) and format_items, "formats 清单为空")
        format_ids = [item.get("id") if isinstance(item, dict) else str(item) for item in format_items]
        for required in ("auto", "visium", "h5ad", "stereo", "matrix_text", "archive"):
            check(required in format_ids, f"formats 缺少空间格式 {required}（实际：{format_ids}）")
        summary.append(f"formats: {len(format_items)} 个（{','.join(format_ids[:6])}…）")

        # 3) analyses：空间分析目录 ≥ 25 项
        #    注意：sp_engine.cmd_analyses 目前要求「已有 active 数据集」（sp_state.State.get(None)
        #    会抛 not_found），所以没有数据集时先记下并等 demo 之后重试，避免把引擎行为当测试失败。
        analyses_resp = worker.call(3, "analyses")
        need_dataset_for_analyses = False
        if analyses_resp.get("ok") is not True:
            error = analyses_resp.get("error") or {}
            check(error.get("code") == "not_found",
                  f"analyses 既没成功也不是 not_found：{json.dumps(error, ensure_ascii=False)}")
            need_dataset_for_analyses = True

        # 4) demo：内置合成组织（2400 spot / 5 空间域）能生成且 n_spots > 0
        demo = unwrap(worker.call(4, "demo", timeout=DEMO_TIMEOUT), "demo")
        dataset = demo.get("dataset") or {}
        n_spots = dataset.get("n_spots")
        check(isinstance(n_spots, (int, float)) and n_spots > 0, f"demo.dataset.n_spots 非法：{n_spots}")
        check(bool(dataset.get("id")), "demo.dataset.id 为空")
        preview = demo.get("preview") or {}
        check(bool(preview.get("genes")), "demo.preview.genes 为空")
        check(preview.get("coords_meta") is not None, "demo 数据应带空间坐标（coords_meta 为空）")
        summary.append(f"demo: {dataset.get('name')} {n_spots} spot × {dataset.get('n_genes')} 基因")

        if need_dataset_for_analyses:
            analyses_resp = worker.call(30, "analyses")
            summary.append("analyses: 无数据集时返回 not_found（引擎现状），demo 之后重试")
        analyses = unwrap(analyses_resp, "analyses")
        analysis_items = analyses.get("analyses") or []
        check(len(analysis_items) >= MIN_ANALYSES,
              f"analyses 只有 {len(analysis_items)} 项（要求 ≥ {MIN_ANALYSES}）")
        groups = analyses.get("groups") or []
        check(bool(groups), "analyses.groups 为空")
        summary.append(f"analyses: {len(analysis_items)} 项 / 分组 {len(groups)}")

        # 5) bogus：未知命令必须返回 bad_request，而不是崩掉 worker
        bogus = worker.call(5, "bogus")
        check(bogus.get("ok") is False, "未知命令应返回 ok=false")
        error = bogus.get("error") or {}
        check(error.get("code") == "bad_request",
              f"未知命令的 error.code 应为 bad_request，实际 {error.get('code')!r}：{error.get('message')!r}")
        check(bool(error.get("message")), "未知命令缺少可读中文 message")
        summary.append(f"bogus: code={error.get('code')} msg={error.get('message')}")

        # 6) 出错后 worker 仍然存活（单条错误不影响进程）
        ping2 = unwrap(worker.call(6, "ping"), "ping(again)")
        check(bool(ping2.get("python")), "bogus 之后 ping 不再返回 python（worker 可能已崩溃）")
        summary.append("worker: 错误命令后仍可继续 ping")
    finally:
        worker.close()
        shutil.rmtree(workdir, ignore_errors=True)

    for line in summary:
        print(f"  - {line}")
    print("PASS")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
    except Exception as error:  # noqa: BLE001
        print(f"FAIL({type(error).__name__}): {error}", file=sys.stderr)
        sys.exit(2)
