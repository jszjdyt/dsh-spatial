#!/usr/bin/env python3
"""
dsh-spatial — Python worker（长驻子进程）

协议（DESIGN.md §2.1，与 dsh-geobulk / dsh-singlecell 完全一致）：
  stdin  : 一行一个 JSON 请求 {"id":"..","cmd":"..","args":{..}}
  stdout : 一行一个 JSON 响应 {"id":"..","ok":true,"result":{..}} /
           {"id":"..","ok":false,"error":{"code":"..","message":"..","trace":"可选"}}
           可选 {"id":"..","event":"progress","message":"..","pct":0-100}
  stdout 只允许协议行；日志一律 stderr。

启动：python3 -u lib/sp_engine.py --workdir <abs_dir>

命令（DESIGN.md §2.2）：
  ping/open/demo/list/select/close/reset/preview/samples/run/plot/analyses/analysis/
  export/audit/suggest/recipe/formats/steps/inspect/geo_search/geo_detail/geo_download/
  refs/ref_set/shutdown

依赖分层：
  * 本文件 + sp_registry + sp_common 是 Host 层（Lead 维护）；
  * sp_state / sp_io / sp_geo 由 sp-core 交付，sp_prep / sp_plots / sp_spatial /
    sp_domain / sp_deconv / sp_lr 由各模块 owner 交付；
  * 任何一个缺失都只在真正用到时给出中文可读的 prerequisite，worker 本身照常 ping。
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)


def _prepare_env(workdir: str) -> None:
    """matplotlib 只能用 Agg，且 ~/.matplotlib 不可写（DESIGN §10）。"""
    os.environ["MPLBACKEND"] = "Agg"
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(workdir, "mpl-cache"))
    os.environ["DS_SP_WORKDIR"] = workdir
    try:
        os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
    except OSError:
        os.environ["MPLCONFIGDIR"] = "/tmp"
    for sub in ("figures", "exports", "uploads", "downloads", "analyses", "unpack",
                "sc-ref", "refs", "mpl-cache"):
        try:
            os.makedirs(os.path.join(workdir, sub), exist_ok=True)
        except OSError:
            pass


class Worker:
    def __init__(self, workdir: str):
        self.workdir = workdir
        self.started = time.time()
        self.state = None
        self._engine_info_cache = None

    # -- 基础设施 ---------------------------------------------------------
    def load_state(self):
        if self.state is None:
            import sp_state

            self.state = sp_state.State()
        return self.state

    def engine_info(self) -> dict:
        if self._engine_info_cache is not None:
            return self._engine_info_cache
        import platform

        packages: dict = {}
        for name in ("numpy", "pandas", "scipy", "sklearn", "matplotlib", "h5py",
                     "anndata", "statsmodels", "openpyxl", "networkx", "PIL",
                     "scanpy", "squidpy", "umap", "leidenalg", "igraph", "skmisc", "harmonypy"):
            try:
                module = __import__(name)
                packages[name] = getattr(module, "__version__", "?")
            except Exception:  # noqa: BLE001
                packages[name] = None
        info = {
            "python": platform.python_version(),
            "platform": sys.platform,
            "mode": "builtin",
            "workdir": self.workdir,
            "packages": packages,
            "started_at": self.started,
        }
        self._engine_info_cache = info
        return info

    def _ref_status(self) -> dict:
        """单细胞参考能力（DESIGN §2.2 ping 的 ref:{available}）。"""
        try:
            import anndata  # noqa: F401

            return {"available": True, "reason": None}
        except Exception as error:  # noqa: BLE001
            return {"available": False, "reason": f"anndata 不可用：{error}"}

    # -- 基础命令 ---------------------------------------------------------
    def cmd_ping(self, args):
        import sp_registry

        info = dict(self.engine_info())
        info["modules"] = sp_registry.module_status()
        empty = sp_registry._EmptyDataset()
        catalog = sp_registry.analyses_catalog(empty)
        formats = sp_registry.formats_catalog()
        info["capabilities"] = {
            "steps": sp_registry.steps_catalog(),
            "analyses": [item["id"] for item in catalog["analyses"]],
            "plots": sp_registry.plot_kinds(),
            "formats": list(formats.get("ids") or []),
            "recipes": [item["id"] for item in sp_registry.RECIPES],
        }
        info["geo"] = True
        info["ref"] = self._ref_status()
        return info

    def cmd_formats(self, args):
        import sp_registry

        return sp_registry.formats_catalog()

    def cmd_steps(self, args):
        import sp_registry

        return {"steps": sp_registry.steps_catalog(), "modules": sp_registry.module_status()}

    def cmd_open(self, args):
        import sp_registry

        state = self.load_state()
        if args.get("demo"):
            dataset = sp_registry.make_demo(state, name=args.get("name") or "demo-spatial")
            try:
                state.select(dataset.id)   # 新建/载入的数据集应当成为当前数据集
            except Exception:  # noqa: BLE001
                state.active = dataset.id
            return {"dataset": sp_registry._summary(dataset), "preview": sp_registry.preview(dataset)}
        path = args.get("path") or args.get("url")
        options = {
            "data_format": args.get("data_format") or args.get("format") or "auto",
            "member": args.get("member"),
            "object": args.get("object"),
            "transpose": args.get("transpose"),
            "delimiter": args.get("delimiter"),
            "sheet": args.get("sheet"),
            "gene_col": args.get("gene_col"),
            "coords": args.get("coords"),
            "sample_sheet": args.get("sample_sheet"),
            "image": args.get("image"),
            "visium_scale": args.get("visium_scale"),
            "species": args.get("species") or "auto",
        }
        dataset = sp_registry.load_dataset(state, self.workdir, path, options=options,
                                           name=args.get("name"))
        try:
            state.select(dataset.id)   # 新导入的数据集应当成为当前数据集
        except Exception:  # noqa: BLE001
            state.active = dataset.id
        payload = {"dataset": sp_registry._summary(dataset), "preview": sp_registry.preview(dataset)}
        if getattr(dataset, "_reused", False):
            payload["reused"] = True
        return payload

    def cmd_demo(self, args):
        return self.cmd_open({"demo": True, "name": args.get("name")})

    def cmd_inspect(self, args):
        import sp_registry

        return sp_registry.inspect_source(self.workdir, args or {})

    def cmd_list(self, args):
        state = self.load_state()
        return {"datasets": state.summaries(), "active": state.active}

    def cmd_select(self, args):
        state = self.load_state()
        target = str(args.get("dataset") or "")
        if target not in getattr(state, "datasets", {}):
            raise _err("not_found", f"数据集 {target} 不存在")
        state.active = target
        import sp_registry

        return {"active": state.active, "dataset": sp_registry._summary(state.datasets[target])}

    def cmd_close(self, args):
        state = self.load_state()
        target = str(args.get("dataset") or state.active or "")
        if not target or target not in getattr(state, "datasets", {}):
            return {"closed": False, "datasets": state.summaries(), "active": state.active}
        state.close(target)
        return {"closed": True, "datasets": state.summaries(), "active": state.active}

    def cmd_reset(self, args):
        state = self.load_state()
        state.reset()
        return {"reset": True, "datasets": [], "active": None}

    # -- 数据集操作 -------------------------------------------------------
    def cmd_preview(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        return sp_registry.preview(dataset)

    def cmd_samples(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        payload = dict(args or {})
        payload.pop("dataset", None)
        return sp_registry.handle_samples(dataset, payload)

    def cmd_run(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        step = str(args.get("step") or "")
        if not step:
            raise _err("bad_request", "请提供 step")
        return sp_registry.run_step(dataset, step, args.get("params") or {}, self.workdir,
                                    self.engine_info())

    def cmd_plot(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        kind = str(args.get("kind") or "")
        if not kind:
            raise _err("bad_request", "请提供图类型 kind")
        return sp_registry.run_plot(dataset, kind, args.get("params") or {}, self.workdir,
                                    self.engine_info())

    def cmd_analyses(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        return sp_registry.analyses_catalog(dataset)

    def cmd_analysis(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        analysis_id = str(args.get("id") or "")
        if not analysis_id:
            raise _err("bad_request", "请提供分析 id")
        return sp_registry.run_analysis(dataset, analysis_id, args.get("params") or {},
                                        self.workdir, self.engine_info())

    def cmd_export(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        return sp_registry.export_dataset(dataset, str(args.get("what") or "summary_json"),
                                          args.get("path"), self.workdir, self.engine_info())

    def cmd_audit(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        return sp_registry.audit(dataset)

    def cmd_suggest(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        return sp_registry.suggest(dataset)

    def cmd_recipe(self, args):
        import sp_registry

        return sp_registry.recipe(str(args.get("id")) if args.get("id") else None)

    # -- 单细胞参考 -------------------------------------------------------
    def cmd_refs(self, args):
        import sp_registry

        path = args.get("path")
        if path:
            return sp_registry.ref_info(str(path))
        return sp_registry.refs_local(self.workdir)

    def cmd_ref_set(self, args):
        import sp_registry

        dataset = self.load_state().get(args.get("dataset"))
        return sp_registry.ref_set(dataset, args.get("ref_path"), args.get("cell_type_col"))

    # -- GEO --------------------------------------------------------------
    def _geo(self):
        try:
            import sp_geo

            return sp_geo
        except Exception as error:  # noqa: BLE001
            raise _err("prerequisite",
                       f"GEO 模块（sp_geo.py）尚未就绪：{type(error).__name__}: {error}")

    def cmd_geo_search(self, args):
        geo = self._geo()
        return geo.geo_search(
            term=str(args.get("term") or ""),
            species=str(args.get("species") or "all"),
            study_type=str(args.get("study_type") or "all"),
            min_samples=int(args.get("min_samples") or 0),
            has_matrix=bool(args.get("has_matrix")) if args.get("has_matrix") is not None else None,
            platform=str(args.get("platform") or "all"),
            limit=int(args.get("limit") or 20),
            page=int(args.get("page") or 0),
            sort=str(args.get("sort") or "relevance"),
        )

    def cmd_geo_detail(self, args):
        geo = self._geo()
        return geo.geo_detail(str(args.get("gse") or ""))

    def cmd_geo_download(self, args):
        geo = self._geo()
        request_id = args.get("__request_id__")

        def progress(message, pct=None):
            if not request_id:
                return
            payload = {"id": request_id, "event": "progress", "message": str(message)}
            if pct is not None:
                payload["pct"] = float(pct)
            _emit(payload)

        gse = str(args.get("gse") or "")
        if not gse:
            raise _err("bad_request", "请提供 GSE 编号")
        files = args.get("files") or []
        dest = str(args.get("dest_dir") or os.path.join(self.workdir, "downloads", "GEO"))
        return geo.geo_download(gse, files, dest, overwrite=bool(args.get("overwrite")),
                                progress=progress)

    def cmd_shutdown(self, args):
        return {"bye": True}

    # -- 分发 -------------------------------------------------------------
    def dispatch(self, cmd: str, args: dict):
        handler = getattr(self, f"cmd_{cmd}", None)
        if handler is None:
            raise _err("bad_request", f"未知命令：{cmd}")
        return handler(args or {})


class _EngineFailure(Exception):
    def __init__(self, code: str, message: str, **details):
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details


def _err(code: str, message: str, **details):
    return _EngineFailure(code, message, **details)


def _emit(payload: dict) -> None:
    sys.stdout.write(json.dumps(payload, ensure_ascii=False) + "\n")
    sys.stdout.flush()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", default=os.path.join(os.path.dirname(HERE), "work"))
    parsed, _ = parser.parse_known_args()
    workdir = os.path.abspath(parsed.workdir)
    _prepare_env(workdir)

    import sp_common

    worker = Worker(workdir)
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        request_id = None
        try:
            payload = json.loads(line)
            request_id = payload.get("id")
            cmd = str(payload.get("cmd") or "")
            args = payload.get("args") or {}
            if cmd == "geo_download":
                args = dict(args)
                args["__request_id__"] = request_id
            result = worker.dispatch(cmd, args)
            _emit({"id": request_id, "ok": True, "result": sp_common.jsonable(result)})
            if cmd == "shutdown":
                return 0
        except _EngineFailure as error:
            _emit({"id": request_id, "ok": False,
                   "error": {"code": error.code, "message": error.message,
                             **sp_common.jsonable(error.details)}})
        except sp_common.EngineError as error:
            _emit({"id": request_id, "ok": False, "error": error.to_dict()})
        except Exception as error:  # noqa: BLE001
            trace = traceback.format_exc()
            sys.stderr.write(trace + "\n")
            sys.stderr.flush()
            _emit({"id": request_id, "ok": False,
                   "error": {"code": "internal", "message": f"{type(error).__name__}: {error}",
                             "trace": trace[-2000:]}})
    return 0


if __name__ == "__main__":
    sys.exit(main())
