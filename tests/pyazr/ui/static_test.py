"""The page's files and how the server serves them (docs/upgrade/web_ui.md sections 7 to 10):

- the page, its modules and the vendored plotly.js are served by the static rows (pyazr.ui.routes.STATIC_ROUTES), with
  the token exchanged for a cookie on the first load and a Content-Security-Policy on every answer; nothing outside
  pyazr/ui/static is served;
- every file the page loads (index.html's scripts and style sheet, every module import) exists and is tracked;
- the vendored plotly.js is the bundle whose hash was checked when it was fetched, with its licence beside it;
- forms.js has a widget for every value type of the engine's schema (by name, in its one widget table);
- no module of the page names a record kind, quantity or parameter kind of the engine outside the layout hints;
- the page's own sources keep the repository's measures (lines of at most 120 columns, files of about 500 lines);
  the vendored directory is the stated exception (minified third-party code).
"""
import glob
import hashlib
import os
import re
import subprocess

import uitest
from uitest import AUTH, TOKEN, copy_project, make_client

STATIC = os.path.join(uitest.ROOT, "pyazr", "ui", "static")
PLOTLY = "vendor/plotly-cartesian-4.1.2.min.js"
PLOTLY_SHA256 = "aaf12ff1d4671f7e93f6df4291651ab4ca6b76c0ea5f1576ec1030e3c2819a20"
MAX_COLUMNS, MAX_LINES = 120, 500


def page_sources():
    """The page's own files (relative to static/), the vendored directory left out."""
    out = []
    for path in glob.glob(os.path.join(STATIC, "**", "*"), recursive=True):
        rel = os.path.relpath(path, STATIC)
        if os.path.isfile(path) and not rel.startswith("vendor" + os.sep):
            out.append(rel)
    return sorted(out)


def test_page_served_with_token_cookie_and_csp(tmp_path):
    with make_client(copy_project("T1", tmp_path)) as client:
        bare = client.get("/", headers={"Authorization": ""})
        assert bare.status_code == 401
        first = client.get(f"/?token={TOKEN}", headers={"Authorization": ""}, follow_redirects=False)
        assert first.status_code == 303 and first.headers["location"] == "/"
        assert "azure2_ui_token=" in first.headers["set-cookie"] and "HttpOnly" in first.headers["set-cookie"]
        client.cookies.set("azure2_ui_token", TOKEN)
        page = client.get("/", headers={"Authorization": ""})
        assert page.status_code == 200 and page.headers["content-type"].startswith("text/html")
        csp = page.headers["content-security-policy"]
        assert "script-src 'self'" in csp and "unsafe-eval" not in csp and "default-src 'none'" in csp
        module = client.get("/static/app.js", headers={"Authorization": ""})
        assert module.status_code == 200 and module.headers["content-type"].startswith("text/javascript")
        assert client.get("/static/" + PLOTLY).status_code == 200
        for outside in ("/static/../routes.py", "/static/%2e%2e/server.py", "/static/nothing.js"):
            assert client.get(outside).status_code == 404, outside
        assert client.get("/api/project").headers["content-security-policy"] == csp


def loaded_files():
    """Every file the page loads: index.html's scripts and style sheets, then every module import, followed."""
    with open(os.path.join(STATIC, "index.html")) as f:
        html = f.read()
    todo = [re.sub(r"^static/", "", p) for p in re.findall(r'(?:src|href)="([^"]+)"', html)]
    seen = set()
    while todo:
        rel = os.path.normpath(todo.pop())
        if rel in seen:
            continue
        seen.add(rel)
        if rel.endswith(".js") and not rel.startswith("vendor"):
            with open(os.path.join(STATIC, rel)) as f:
                source = f.read()
            for spec in re.findall(r'(?:^|\s)(?:import|export)\s[^;]*?from\s+"([^"]+)"|import\(\s*"([^"]+)"\s*\)',
                                   source, re.M):
                target = spec[0] or spec[1]
                todo.append(os.path.join(os.path.dirname(rel), target))
    return sorted(seen)


def test_every_file_the_page_loads_exists_and_is_tracked():
    files = loaded_files()
    assert "app.js" in files and PLOTLY in files and "tabs/data.js" in files and "selftest.js" in files, files
    missing = [f for f in files if not os.path.isfile(os.path.join(STATIC, f))]
    assert not missing, f"the page loads files that do not exist: {missing}"
    tracked = subprocess.run(["git", "ls-files", "pyazr/ui/static"], cwd=uitest.ROOT, capture_output=True, text=True)
    if tracked.returncode == 0 and tracked.stdout:
        names = {os.path.relpath(os.path.join(uitest.ROOT, p), STATIC) for p in tracked.stdout.split()}
        untracked = [f for f in files if f not in names]
        assert not untracked, f"the page loads files git does not track: {untracked}"
    unused = [f for f in page_sources() if f.endswith(".js") and f not in files]
    assert not unused, f"modules no page loads: {unused}"


def test_vendored_plotly_is_the_checked_bundle_with_its_licence():
    with open(os.path.join(STATIC, PLOTLY), "rb") as f:
        assert hashlib.sha256(f.read()).hexdigest() == PLOTLY_SHA256
    with open(os.path.join(STATIC, "vendor", "plotly-LICENSE.txt")) as f:
        assert "MIT License" in f.read()
    with open(os.path.join(uitest.ROOT, "README.md")) as f:
        assert "plotly.js" in f.read(), "README.md lists the vendored third-party code"


def widget_names():
    with open(os.path.join(STATIC, "forms.js")) as f:
        source = f.read()
    block = source.split("export const WIDGETS = {", 1)[1].split("\n};", 1)[0]
    return set(re.findall(r"^  (\w+):", block, re.M))


def test_widget_table_has_every_value_type_of_the_schema():
    import pyazr
    types = set(pyazr.schema()["value_types"])
    names = widget_names()
    assert not types - names, f"value types of the engine without a widget in forms.js: {sorted(types - names)}"
    assert not names - types, f"widgets of value types the engine does not have: {sorted(names - types)}"


def test_no_engine_kind_is_spelled_in_the_page():
    import pyazr
    from pyazr.parameters import parameter_kinds
    from routes_consistency_test import NOT_KINDS
    schema = pyazr.schema()
    words = (set(schema["records"]) | set(schema["quantities"]) | {k["key"] for k in parameter_kinds()}) - NOT_KINDS
    for rel in page_sources():
        if not rel.endswith(".js") or rel == "layout.js":  # layout.js: the one table of layout hints per kind
            continue
        with open(os.path.join(STATIC, rel)) as f:
            source = f.read()
        spelled = sorted(w for w in words if re.search(rf"[\"'`]{re.escape(w)}[\"'`]", source))
        assert not spelled, f"{rel} names kinds of the engine: {spelled}"


def test_quantity_axes_are_ones_the_plot_row_knows():
    import pyazr
    for name, q in pyazr.schema()["quantities"].items():
        cols = q["columns"]
        assert cols["angle_axis"] in ("none", "required", "optional"), (name, cols)
        if "dataset" in q["record_kinds"] and cols["energy"]:
            assert q["engine_unit"], name


def test_page_sources_keep_the_measures():
    for rel in page_sources():
        with open(os.path.join(STATIC, rel), encoding="utf-8") as f:
            lines = f.read().split("\n")
        long = [i + 1 for i, line in enumerate(lines) if len(line) > MAX_COLUMNS]
        assert not long, f"{rel}: lines longer than {MAX_COLUMNS} columns: {long}"
        assert len(lines) <= MAX_LINES + 20, f"{rel}: {len(lines)} lines (about {MAX_LINES} at most)"


def test_file_head_route_stays_inside_the_root(tmp_path):
    with make_client(copy_project("T1", tmp_path)) as client:
        head = client.get("/api/files/T1.azr", params={"rows": 2}).json()
        assert len(head["lines"]) == 2 and head["more"] is True
        assert client.get("/api/files/..%2F..%2Fetc%2Fpasswd").status_code in (403, 404)
        assert client.get("/api/files/no_such_file.dat").status_code == 404


if __name__ == "__main__":
    uitest.run(__file__)
