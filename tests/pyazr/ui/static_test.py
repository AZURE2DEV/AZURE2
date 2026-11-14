"""The page's files and how the server serves them (docs/upgrade/web_ui.md sections 7 to 10):

- the page, its modules and the vendored plotly.js are served by the static rows (pyazr.ui.routes.STATIC_ROUTES), with
  the token exchanged for a cookie on the first load and a Content-Security-Policy on every answer; nothing outside
  pyazr/ui/static is served;
- every file the page loads (index.html's scripts and style sheet, every module import) exists and is tracked;
- the vendored plotly.js is the bundle whose hash was checked when it was fetched, with its licence beside it;
- forms.js has a widget for every value type of the engine's schema (by name, in its one widget table);
- no module of the page names a record kind, quantity or parameter kind of the engine outside the layout hints;
- the page's own sources keep the repository's measures (lines of at most 120 columns, files of about 500 lines);
  the vendored directory is the stated exception (minified third-party code);
- the page shows the original program's icon (the file of packaging/icons, served as a PNG, the tab's icon and the
  header's), its colour tokens give every text at least 4.5:1 on its backgrounds in both themes, and its About box
  cites what the documentation's front page cites.
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


def test_browser_without_token_gets_a_page_that_asks_for_it(tmp_path):
    """A navigation (GET accepting HTML, outside /api) without the token or with a wrong one: a page with a field
    for the token, not the JSON refusal; the field's form is the printed URL's door (GET /?token=)."""
    from pyazr.ui.server import cookie_name
    html = {"Authorization": "", "Accept": "text/html,application/xhtml+xml,*/*;q=0.8"}
    with make_client(copy_project("T1", tmp_path)) as client:
        for path, sentence in (("/", "opened without its token"), ("/?token=nope", "not this server's")):
            asked = client.get(path, headers=html, follow_redirects=False)
            assert asked.status_code == 401 and asked.headers["content-type"].startswith("text/html"), path
            assert sentence in asked.text and '<form method="get" action="/">' in asked.text, asked.text
            assert 'name="token"' in asked.text and TOKEN not in asked.text
            assert "form-action 'self'" in asked.headers["content-security-policy"]
            assert "script-src" not in asked.headers["content-security-policy"] and "<script" not in asked.text
            assert "set-cookie" not in asked.headers
        # a program (no HTML asked for) and the API keep the JSON refusal with its three parts
        for path, headers in (("/", {"Authorization": ""}), ("/api/project", html)):
            refused = client.get(path, headers=headers)
            assert refused.status_code == 401 and refused.json()["error"] == "Unauthorized", path
            assert "?token=<token>" in refused.json()["fix"]
        # what the form sends is the printed URL
        taken = client.get(f"/?token={TOKEN}", headers=html, follow_redirects=False)
        assert taken.status_code == 303 and taken.headers["location"] == "/"
    assert cookie_name("127.0.0.1:28701") == "azure2_ui_token_28701" and cookie_name("[::1]:9000").endswith("_9000")
    assert cookie_name("localhost") == cookie_name("[::1]") == cookie_name("") == "azure2_ui_token_80"


def test_two_servers_on_one_host_keep_their_own_cookies(tmp_path):
    """A browser has one cookie jar per host whatever the port: the cookie of the server on one port is not read as
    the token of the server on another (it was: opening a second project logged the first one out)."""
    with make_client(copy_project("T1", tmp_path)) as client:
        client.cookies.set("azure2_ui_token_9000", "anotherServersToken")
        client.cookies.set("azure2_ui_token", "theNameBeforeThePortWasPartOfIt")
        assert client.get("/", headers={"Authorization": ""}).status_code == 401
        client.cookies.set("azure2_ui_token_80", TOKEN)
        assert client.get("/", headers={"Authorization": ""}).status_code == 200
        assert client.get("/", headers={"Authorization": "", "Host": "127.0.0.1:9000"}).status_code == 401


def test_page_served_with_token_cookie_and_csp(tmp_path):
    with make_client(copy_project("T1", tmp_path)) as client:
        bare = client.get("/", headers={"Authorization": ""})
        assert bare.status_code == 401
        first = client.get(f"/?token={TOKEN}", headers={"Authorization": ""}, follow_redirects=False)
        assert first.status_code == 303 and first.headers["location"] == "/"
        # the cookie is named for the port (two servers on one machine do not overwrite each other's) and is lax:
        # a strict one is kept back when the printed URL is followed from another site, and the redirect then fails
        cookie = first.headers["set-cookie"]
        assert cookie.startswith(f"azure2_ui_token_80={TOKEN};") and "HttpOnly" in cookie, cookie
        assert "samesite=lax" in cookie.lower(), cookie
        client.cookies.set("azure2_ui_token_80", TOKEN)
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
        if rel.endswith(".png"):
            continue
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



def test_file_head_is_bounded(tmp_path):
    """The rows asked are 1 to 200 (refused outside with the range); a line is read to 4096 bytes and marked cut, so
    a binary file without line ends is not read into memory whole."""
    from pyazr.ui.server import FILE_HEAD_LINE_BYTES, FILE_HEAD_ROWS
    project = copy_project("T1", tmp_path)
    directory = os.path.dirname(project)
    with open(os.path.join(directory, "one_line.bin"), "wb") as f:
        f.write(b"x" * (10 * FILE_HEAD_LINE_BYTES))
    with open(os.path.join(directory, "lines.txt"), "wb") as f:
        f.write(b"short\n" + b"y" * (FILE_HEAD_LINE_BYTES + 5) + b"\nlast\n")
    with make_client(project) as client:
        for rows in (0, -1, FILE_HEAD_ROWS + 1):
            r = client.get("/api/files/T1.azr", params={"rows": rows})
            assert r.status_code == 422 and str(FILE_HEAD_ROWS) in r.json()["allowed"], (rows, r.json())
        assert len(client.get("/api/files/T1.azr", params={"rows": FILE_HEAD_ROWS}).json()["lines"]) > 20
        head = client.get("/api/files/one_line.bin").json()
        assert len(head["lines"]) == 1 and len(head["lines"][0]) == FILE_HEAD_LINE_BYTES and head["cut"] == [1]
        assert head["more"] is False
        head = client.get("/api/files/lines.txt", params={"rows": 3}).json()
        assert head["lines"][0] == "short" and len(head["lines"][1]) == FILE_HEAD_LINE_BYTES
        assert head["lines"][2] == "last" and head["cut"] == [2] and head["more"] is False

if __name__ == "__main__":
    uitest.run(__file__)


ICON = "azure2_icon.png"


def test_icon_is_the_original_programs(tmp_path):
    with open(os.path.join(STATIC, ICON), "rb") as a, open(os.path.join(uitest.ROOT, "packaging", "icons", ICON),
                                                           "rb") as b:
        assert a.read() == b.read(), "the page's icon is not packaging/icons/azure2_icon.png"
    assert ICON in loaded_files(), "index.html does not name the icon (the browser tab's icon)"
    with open(os.path.join(STATIC, "app.js")) as f:
        assert f'src: "static/{ICON}"' in f.read(), "the header does not show the icon"
    with make_client(copy_project("T1", tmp_path)) as client:
        r = client.get(f"/static/{ICON}", headers=AUTH)
        assert r.status_code == 200 and r.headers["content-type"] == "image/png", r.headers


def css_tokens(css):
    """{theme: {token: colour}} of style.css: the light tokens of the first :root block, the dark ones of the
    :root[data-theme="dark"] block; the dark block of the system's preference must be the same."""
    def block(selector):
        m = re.search(re.escape(selector) + r"\s*\{([^}]*)\}", css)
        assert m, f"style.css has no block {selector}"
        return dict(re.findall(r"--([\w-]+):\s*(#[0-9a-fA-F]{6})", m.group(1)))
    light, dark = block(":root"), block(':root[data-theme="dark"]')
    assert block(':root:not([data-theme="light"])') == dark, "the two dark blocks of style.css differ"
    return {"light": light, "dark": dark}


def contrast(a, b):
    """The WCAG 2 contrast ratio of two #rrggbb colours."""
    def luminance(c):
        rgb = [int(c[i:i + 2], 16) / 255 for i in (1, 3, 5)]
        lin = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb]
        return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    hi, lo = sorted((luminance(a), luminance(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


# Every text token on every background a text is drawn on (style.css), and the pairs of a box: (text, background).
TEXT_TOKENS = ["fg", "muted", "accent", "error", "warn", "ok"]
BACKGROUNDS = ["bg", "panel", "accent-soft", "row-hover"]
BOX_PAIRS = [("accent-fg", "accent"), ("fg", "error-bg"), ("muted", "error-bg"), ("error", "error-bg"),
             ("fg", "warn-bg"), ("muted", "warn-bg"), ("warn", "warn-bg")]


def test_text_tokens_have_contrast_in_both_themes():
    with open(os.path.join(STATIC, "style.css")) as f:
        themes = css_tokens(f.read())
    low = []
    for theme, t in themes.items():
        pairs = [(a, b) for a in TEXT_TOKENS for b in BACKGROUNDS] + BOX_PAIRS
        for a, b in pairs:
            ratio = contrast(t[a], t[b])
            if ratio < 4.5:
                low.append(f"{theme}: --{a} {t[a]} on --{b} {t[b]}: {ratio:.2f}")
    assert not low, "text below 4.5:1: " + "; ".join(low)


def test_contrast_rule_bites():
    """The check above fails on a pair that is too pale (the measure is not vacuous)."""
    assert contrast("#777777", "#ffffff") < 4.5 < contrast("#004293", "#ffffff")


def test_accent_is_the_icons_blue():
    with open(os.path.join(STATIC, "style.css")) as f:
        themes = css_tokens(f.read())
    assert themes["light"]["accent"] == "#004293" and themes["light"]["error"] == "#990000", themes["light"]


def test_about_cites_what_the_documentation_cites():
    with open(os.path.join(uitest.ROOT, "docs", "source", "index.rst")) as f:
        rst = f.read()
    m = re.search(r"should cite:\s*\n\s*\n(.*?)\n\s*\n", rst, re.S)
    assert m, "docs/source/index.rst: no citation found"
    cited = re.sub(r"\s+", " ", re.sub(r"\*+", "", m.group(1))).strip().rstrip(".")
    with open(os.path.join(STATIC, "app.js")) as f:
        app = f.read()
    m = re.search(r"citation: ((?:\"[^\"]*\"\s*\+?\s*)+),", app)
    assert m, "app.js: no ABOUT.citation"
    shown = "".join(re.findall(r'"([^"]*)"', m.group(1)))
    assert shown == cited, (shown, cited)


def test_files_route_lists_the_project_folder(tmp_path):
    """GET /api/files: the project folder's files (not output/, checks/, hidden ones), paths as a field names them,
    a numeric table with its columns; a symbolic link leading out of the root is not listed."""
    project = copy_project("T2", tmp_path)
    directory = os.path.dirname(project)
    os.symlink("/etc/hostname", os.path.join(directory, "data", "outside.dat"))
    with open(os.path.join(directory, "output", "x.dat"), "w") as f:
        f.write("1 2\n")
    with make_client(project) as client:
        answer = client.get("/api/files").json()
    files = {f["path"]: f for f in answer["files"]}
    assert files["T2.azr"]["kind"] == "project"
    assert files["data/s01_int_thin.dat"]["kind"] == "table" and files["data/s01_int_thin.dat"]["columns"] >= 3
    size = os.path.getsize(os.path.join(directory, "data", "s01_int_thin.dat"))
    assert files["data/s01_int_thin.dat"]["size"] == size
    assert "data/outside.dat" not in files and not any(p.startswith(("output/", "checks/", ".")) for p in files)
    assert answer["more"] is False


def test_field_preview_is_the_engines_reading_on_a_copy(tmp_path):
    """GET /api/project/records/{name}/preview: the description of the record with the field set, as
    Project.set then Project.describe give it on a copy; a file the engine cannot read is refused with its
    diagnostics; outside the root 403; the project is not changed."""
    from pyazr import Project
    project = copy_project("T2", tmp_path)
    with make_client(project) as client:
        before = client.get("/api/project/text").text
        got = client.get("/api/project/records/s01_int_thin/preview",
                         params={"field": "file", "value": "data/s04_thick.dat"}).json()
        direct = Project.open(project)
        direct.set("s01_int_thin", file="data/s04_thick.dat")
        assert got["description"] == direct.describe("s01_int_thin", as_json=True)
        refused = client.get("/api/project/records/s01_int_thin/preview", params={"field": "file", "value": "T2.azr"})
        assert refused.status_code == 422 and refused.json()["diagnostics"], refused.text
        outside = client.get("/api/project/records/s01_int_thin/preview",
                             params={"field": "file", "value": "../../x.dat"})
        assert outside.status_code == 403
        assert client.get("/api/project/text").text == before
