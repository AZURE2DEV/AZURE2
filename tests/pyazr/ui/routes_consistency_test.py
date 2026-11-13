#!/usr/bin/env python3
"""The route table of the web UI and its handlers agree (docs/upgrade/web_ui.md section 4.2):

  1. every row names a callable that exists: a method of pyazr.Project, an operation of the adapter
     pyazr.ui.engine.Engine, a method of pyazr.ui.server.Server, and every pyazr call a server row says it serves;
  2. every row answering 501 for a missing pyazr call names one that does not exist (when it appears, this fails:
     serve the route);
  3. every public name of pyazr.Project is served by a row or listed, with its reason, in PROJECT_NOT_SERVED;
  4. every public operation of the adapter is a route or a call the server makes itself; every public method of
     today's pyazr.azure2 is called by the adapter or listed, with its reason, in engine.AZURE2_NOT_SERVED; no other
     file of pyazr/ui names the class azure2;
  5. the app serves exactly the rows (method and path), and the OpenAPI document lists every one of them;
  6. ``import pyazr`` imports neither fastapi nor uvicorn nor pyazr.ui;
  7. every job row runs a recipe or names the call it waits for;
  8. no file of pyazr/ui names or branches on a data kind, record kind or parameter kind of the engine (a new one
     needs no edit of the UI);
  9. every test of tests/pyazr/ui is registered in tests/pyazr/CMakeLists.txt.

Run: python3 tests/pyazr/ui/routes_consistency_test.py
"""
import glob
import importlib
import inspect
import os
import re
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import ROOT, copy_project, make_client  # noqa: E402

UI = os.path.join(ROOT, "pyazr", "ui")


def resolve(dotted):
    """The object named ``pyazr.a.b`` / ``engine.Engine.op`` / ``server.Server.m``; AttributeError if none."""
    head, *rest = dotted.split(".")
    obj = importlib.import_module("pyazr" if head == "pyazr" else f"pyazr.ui.{head}")
    for part in rest:
        obj = getattr(obj, part)
    return obj


def test_every_row_names_an_existing_callable():
    from pyazr.ui.routes import ROUTES
    for r in ROUTES:
        if r.handler in ("project", "session", "server"):
            assert callable(resolve(r.call)), r
            for served in r.serves:
                resolve(served)
        elif r.handler == "missing":
            with pytest.raises((AttributeError, ImportError)):
                resolve(r.call)
        else:
            assert r.handler == "later" and "step" in r.call, r


def test_every_public_name_of_project_is_served_or_listed():
    import pyazr
    from pyazr.ui.routes import PROJECT_NOT_SERVED, ROUTES
    served = {r.call.rsplit(".", 1)[1] for r in ROUTES if r.call.startswith("pyazr.Project.")}
    served |= {s.rsplit(".", 1)[1] for r in ROUTES for s in r.serves if s.startswith("pyazr.Project.")}
    names = {n for n in dir(pyazr.Project) if not n.startswith("_")}
    names = {"add_<kind>" if n.startswith("add_") and n != "add_level" else n for n in names}
    left = sorted(names - served - set(PROJECT_NOT_SERVED))
    assert not left, f"public names of pyazr.Project with no route and no reason: {left}"
    stale = sorted(set(PROJECT_NOT_SERVED) - names)
    assert not stale, f"PROJECT_NOT_SERVED names what Project does not have: {stale}"


def test_adapter_and_azure2_covered():
    import pyazr
    from pyazr.ui.engine import AZURE2_NOT_SERVED, Engine
    from pyazr.ui.routes import ROUTES
    rows = {r.call.rsplit(".", 1)[1] for r in ROUTES if r.handler == "session"}
    with open(os.path.join(UI, "server.py")) as f:
        server_source = f.read()
    for op in [n for n, _ in inspect.getmembers(Engine, inspect.isfunction) if not n.startswith("_")]:
        assert op in rows or f'"{op}"' in server_source, f"Engine.{op}: neither a route nor called by the server"
    with open(os.path.join(UI, "engine.py")) as f:
        adapter = f.read()
    listed = {n for names in AZURE2_NOT_SERVED.values() for n in names.split()}
    public = {n for n in dir(pyazr.azure2) if not n.startswith("_") and n.upper() != n}
    left = sorted(n for n in public - listed if not re.search(rf"\.{n}\b", adapter))
    assert not left, f"public names of pyazr.azure2 neither called by the adapter nor listed: {left}"
    assert not sorted(listed - public), f"AZURE2_NOT_SERVED names what azure2 does not have: {sorted(listed - public)}"
    for path in glob.glob(os.path.join(UI, "*.py")):
        if os.path.basename(path) != "engine.py":
            with open(path) as f:
                assert not re.search(r"\bazure2\(|pyazr\.azure2\b|import azure2", f.read()), path


def test_app_serves_the_rows_and_openapi_lists_them(tmp_path):
    from pyazr.ui.routes import ROUTES
    with make_client(copy_project("T1", tmp_path)) as client:
        served = {(m, r.path) for r in client.app.routes if r.path.startswith("/api")
                  for m in getattr(r, "methods", ())}
        spec = client.get("/openapi.json").json()
    table = {(r.method, r.path.replace(":path", "")) for r in ROUTES}
    assert {(m, p.replace(":path", "")) for m, p in served} == table
    listed = {(m.upper(), p) for p, ops in spec["paths"].items() for m in ops if p.startswith("/api")}
    assert listed == table, f"not in OpenAPI: {sorted(table - listed)}"
    assert len(ROUTES) == len(table), "two rows with the same method and path"


def test_import_pyazr_does_not_import_the_server():
    code = "import sys, pyazr; print(sorted(m for m in ('fastapi', 'uvicorn', 'pyazr.ui') if m in sys.modules))"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=ROOT,
                         env=dict(os.environ, PYTHONPATH=ROOT)).stdout.strip().splitlines()
    assert out[-1] == "[]", out


def test_job_rows():
    from pyazr.ui.jobs import JOBS
    for row in JOBS:
        assert (row.run is None) != (row.missing is None), row
    assert len({(r.kind, r.method) for r in JOBS}) == len(JOBS)


# Words of the engine's tables that occur in pyazr/ui as something else: argument names of Project.add_level and
# Project.set_channel ("energy", "width") and a column of the program's chiSquared.out ("norm").
NOT_KINDS = {"energy", "width", "norm"}


def test_no_engine_kind_is_spelled_in_the_ui():
    """A new data kind, record kind or parameter kind of the engine needs no edit of pyazr/ui: no file of it names a
    quantity, a record kind or a parameter kind of the engine's tables, or branches on one."""
    import pyazr
    from pyazr.parameters import parameter_kinds
    schema = pyazr.schema()
    words = (set(schema["records"]) | set(schema["quantities"]) | {k["key"] for k in parameter_kinds()}) - NOT_KINDS
    for path in glob.glob(os.path.join(UI, "*.py")):
        with open(path) as f:
            source = f.read()
        spelled = sorted(w for w in words if re.search(rf"[\"']{re.escape(w)}[\"']", source))
        assert not spelled, f"{path} names kinds of the engine: {spelled}"
        assert not re.search(r"\b(kind|quantity|record_kind)\s*(==|!=|\bin\b)\s*[\"'(\[{]", source), path


def test_ui_tests_registered():
    with open(os.path.join(ROOT, "tests", "pyazr", "CMakeLists.txt")) as f:
        cmake = f.read()
    for path in glob.glob(os.path.join(uitest.HERE, "*_test.py")):
        assert f"ui/{os.path.basename(path)}" in cmake, f"{path} is not registered in tests/pyazr/CMakeLists.txt"


if __name__ == "__main__":
    uitest.run(__file__)
