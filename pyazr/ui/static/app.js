// The page of AZURE2's web UI (docs/upgrade/web_ui.md sections 5 to 7; the navigation: docs/upgrade/
// web_navigation.md): the top bar (the project, its saved state, the menus, Save / Validate / Calculate always
// visible), the row of main tabs in the old setup program's order and a second row for a group's tabs, their hash
// routes (#/data/s02_diff090), the keys and their overlay, the status line (the engine's last message, the job that
// runs, the validator's verdict), the first visit's overview, and the self-test (#selftest).
//
// The page holds no state of the project.  After every change of revision (an edit here, Undo, Redo, another
// client's edit seen by the poll) it re-reads the project's state, asks the server for the evaluation of the
// current text (the header's χ² is the engine's number, as the server sends it) and draws the open tab again.

import { createApi, watchRevision } from "./api.js";
import { button, errorBox, h } from "./ui.js";
import { baselineOf, changeText, remember, snapshotOf } from "./calculate.js";
import { START, mainTabs, mainTitleOf, shownTabs as tabsOfDepth, tabByArrow, tabByKey, tabOfRoute } from "./layout.js";
import { aboutOverlay, askPath, docsOf, jobMark, keysOverlay, menuBar, menus, overviewDue, showOverview,
  validationMark } from "./shell.js";
import { formatNumber } from "./tables.js";
import * as projectTab from "./tabs/project.js";
import * as nucleusTab from "./tabs/nucleus.js";
import * as dataTab from "./tabs/data.js";
import * as fitTab from "./tabs/fit.js";
import * as resultsTab from "./tabs/results.js";
import * as physicsTab from "./tabs/physics.js";
import * as thmTab from "./tabs/thm.js";
import * as rateTab from "./tabs/rate.js";
import * as covarianceTab from "./tabs/covariance.js";
import * as endfTab from "./tabs/endf.js";
import * as priorTab from "./tabs/prior.js";

/** The tabs: id (the hash route's first word), title, the module that draws it.  Which of them the tab rows show,
 * under which names and in which order, is layout.js FIRST_VIEW. */
export const TABS = [
  { id: "project", title: "Project", module: projectTab },
  { id: "nucleus", title: "Nucleus", module: nucleusTab },
  { id: "data", title: "Data", module: dataTab },
  { id: "fit", title: "Fit", module: fitTab },
  { id: "results", title: "Results", module: resultsTab },
  { id: "physics", title: "Physics", module: physicsTab },
  { id: "thm", title: "THM", module: thmTab, shown: thmTab.shown },
  { id: "rate", title: "Rate", module: rateTab },
  { id: "covariance", title: "Covariance", module: covarianceTab },
  { id: "endf", title: "ENDF", module: endfTab },
  { id: "prior", title: "Prior", module: priorTab },
];

/** The keys of the page, shown by "?" (the overlay lists exactly this table). */
export const KEYS = [
  ["1 … 9, 0", "switch to the main tab of that number, as the tab row shows them (0: the tenth)"],
  ["[ ]", "the previous and the next main tab"],
  ["← → Home End", "in the tab row: move to the previous, next, first, last tab"],
  ["c", "Calculate: every data set with the model, the residuals and χ², and what changed since the last time"],
  ["p", "Validate: the engine validator's messages, each a link to its row"],
  ["e", "show every field and column (the expert's depth) or the usual ones"],
  ["Ctrl+Z", "undo the last edit (outside a text box)"],
  ["Ctrl+Shift+Z, Ctrl+Y", "redo"],
  ["Ctrl+S", "save the project to its file"],
  ["/", "filter the current table"],
  ["↑ ↓ (j k), Enter", "move through a table's rows, open the row"],
  ["Ctrl+Enter", "apply the form being edited"],
  ["Esc", "revert the form being edited; close a menu, a list or this overlay"],
  ["t", "switch the theme: system, light, dark"],
  ["?", "this list"],
];

/** The About box: the program, the reference every publication that uses its results cites (as the documentation's
 * front page gives it, docs/source/index.rst: tests/pyazr/ui/static_test.py compares the two) and the documentation. */
export const ABOUT = {
  citation: "R.E. Azuma, E. Uberseder, E.C. Simpson, C.R. Brune, H. Costantini, R.J. deBoer, J. Gorres, M. Heil, "
    + "P.J. Leblanc, C. Ugalde, and M. Wiescher, Physical Review C 81, 045805 (2010)",
  documentation: "https://azure2dev.github.io/AZURE2/",
};

/** The tabs (layout.js shownTabs over TABS; the same in both depths). */
export const shownTabs = (expert) => tabsOfDepth(expert, TABS);

/** A hash route: "#/data/s02_diff090?theme=dark" -> {tab, args: ["s02_diff090"], params: {theme: "dark"}}; no route:
 * the first tab of the work (layout.js START). */
export function parseRoute(hash) {
  const [path, query = ""] = String(hash || "").replace(/^#\/?/, "").split("?");
  const words = path.split("/").filter(Boolean).map(decodeURIComponent);
  const params = Object.fromEntries(new URLSearchParams(query));
  if (!words.length) return { ...parseRoute(START), params };
  return { tab: words[0], args: words.slice(1), params };
}

const ctx = {
  api: null, schema: null, kinds: null, describe: null, project: null, evaluation: null, evalError: null,
  route: parseRoute(""), statusLine: null, expert: false, calculations: [],
};

/** Show a message of the engine or the page in the status line. */
ctx.status = (text, kind = "info") => {
  if (!ctx.statusLine) return;
  ctx.statusLine.textContent = text;
  ctx.statusLine.className = `status ${kind}`;
  ctx.statusLine.title = text;
};

/** Show a refusal in the status line and, when `where` is given, as a box in that element. */
ctx.showError = (err, where) => {
  const view = err?.view
    || { title: "Error", what: String(err?.message || err), allowed: "", fix: "", diagnostics: [] };
  ctx.status(`${view.title}: ${view.what}`, "error");
  if (where) where.prepend(errorBox(view, { onLink: (d) => ctx.navigate(`#/data/${encodeURIComponent(d.record)}`) }));
  return view;
};

/** Run an edit (a function returning the request's promise), then re-read everything; a refusal is shown. */
ctx.edit = async (what, run, where) => {
  ctx.status(`${what}…`);
  ctx.editing = true;
  try {
    const answer = await run();
    ctx.editing = false;
    const warnings = answer?.warnings || [];
    await ctx.refresh();
    if (warnings.length) ctx.status(`${what}: done; ${warnings.join("; ")}`, "warn");
    else ctx.status(`${what}: done${chi2Text()}${changeSinceText()}`, "ok");
    return answer;
  } catch (e) {
    ctx.editing = false;
    ctx.showError(e, where);
    throw e;
  }
};

ctx.navigate = (hash) => { if (location.hash !== hash) location.hash = hash; else draw(); };

/** Re-read the project's state and its evaluation, then draw the header and the open tab again. */
ctx.refresh = async () => {
  ctx.pairs = null;
  ctx.project = await ctx.api.get("/api/project");
  ctx.api.rev = ctx.project.rev;
  drawHeader();
  await showTabs();
  ctx.validation.refresh();
  ctx.jobMark.poke();
  await evaluate();
  await draw();
};

/** The conditions a tab is shown under (layout.js FIRST_VIEW `when`): the project's own content as the server reads
 * it.  A link of a tab whose condition fails is hidden (its route still draws it). */
const CONDITIONS = { thm: (c) => thmTab.shown(c) };

async function showTabs() {
  await Promise.all(Object.entries(CONDITIONS).map(async ([when, holds]) => {
    let on = false;
    try { on = await holds(ctx); } catch { /* a failed read keeps the tab hidden */ }
    document.querySelectorAll(`a.tab[data-when="${when}"]`).forEach((link) => { link.hidden = !on; });
  }));
}

async function evaluate() {
  ctx.chi2Box.textContent = "χ² …";
  try {
    // a project the engine cannot evaluate (no data sets, a refused set-up) is shown, not a failure of the page
    ctx.evaluation = await ctx.api.request("POST", "/api/evaluate", {}, { quiet: true });
    ctx.evalError = null;
    const w = ctx.evaluation.warnings || [];
    if (w.length) ctx.status(`engine: ${w.map((x) => x.message || x).join("; ")}`, "warn");
    else if (!ctx.editing) ctx.status(`evaluated${chi2Text()}`, "ok");
  } catch (e) {
    ctx.evaluation = null;
    ctx.evalError = e.view || { what: String(e) };
    ctx.status(`no evaluation: ${ctx.evalError.what}`, "warn");
  }
  drawChi2();
}

/** "Calculate": the engine's evaluation of the current text, remembered for the next comparison, then the Plot tab
 * with every data set, its model and residuals, and the changes since the calculation before (calculate.js). */
ctx.calculate = async () => {
  ctx.status("calculating…");
  await evaluate();
  if (ctx.evaluation) ctx.calculations = remember(ctx.calculations, snapshotOf(ctx.evaluation, ctx.project.rev));
  ctx.navigate("#/results?calculated=1");
};

/** The change of the data χ² since the last calculation, for the status line ("" before the first one). */
function changeSinceText() {
  if (!ctx.evaluation || !ctx.calculations.length) return "";
  const before = baselineOf(ctx.calculations, snapshotOf(ctx.evaluation, ctx.project.rev));
  return before ? ` · Δχ² ${changeText(ctx.evaluation.chi2 - before.chi2)} since the last Calculate` : "";
}

/** The full χ² of the last evaluation for the status line (the header shows it to 8 digits). */
const chi2Text = () => (ctx.evaluation ? ` · χ² ${ctx.evaluation.chi2} (data), objective ${ctx.evaluation.objective}`
  : "");

function drawChi2() {
  const e = ctx.evaluation;
  if (!e) {
    ctx.chi2Box.textContent = "χ² —";
    ctx.chi2Box.title = ctx.evalError ? `no evaluation: ${ctx.evalError.what}` : "";
    return;
  }
  const n = e.totals?.["Total-N"];
  ctx.chi2Box.textContent = `χ² ${formatNumber(e.chi2, 8)}${n !== undefined ? ` · N ${n}` : ""}`;
  ctx.chi2Box.title = `data χ² ${e.chi2} (the engine's, at the working x); objective with priors ${e.objective}`;
}

function drawHeader() {
  const p = ctx.project;
  const name = p.path.split(/[\\/]/).pop();
  ctx.nameBox.textContent = name;
  ctx.nameBox.title = p.path;
  ctx.unsavedBox.textContent = p.unsaved ? "● unsaved" : "saved";
  ctx.unsavedBox.className = p.unsaved ? "unsaved" : "saved";
  ctx.unsavedBox.title = p.unsaved ? "edits not yet in the file: Save (Ctrl+S) writes them" : "the file holds the text";
  ctx.revBox.textContent = `rev ${p.rev}`;
  ctx.undoButton.disabled = !p.can_undo;
  ctx.redoButton.disabled = !p.can_redo;
  document.title = `${p.unsaved ? "● " : ""}${name} — AZURE2`;
}

const undo = () => ctx.edit("undo", () => ctx.api.post("/api/project/undo", {}));
const redo = () => ctx.edit("redo", () => ctx.api.post("/api/project/redo", {}));
const save = () => ctx.edit("save", () => ctx.api.post("/api/project/save", {}));

// -- theme ----------------------------------------------------------------------------------------------------------

const THEMES = ["system", "light", "dark"];
function storedTheme() { try { return localStorage.getItem("azure2-theme") || "system"; } catch { return "system"; } }
function applyTheme(theme) {
  if (theme === "system") document.documentElement.removeAttribute("data-theme");
  else document.documentElement.setAttribute("data-theme", theme);
  ctx.theme = theme;
}
function cycleTheme() {
  const next = THEMES[(THEMES.indexOf(storedTheme()) + 1) % THEMES.length];
  try { localStorage.setItem("azure2-theme", next); } catch { /* the theme lasts for this load only */ }
  applyTheme(next);
  draw();
}

// -- depth: the usual fields or every field --------------------------------------------------------------------

function storedDepth() { try { return localStorage.getItem("azure2-depth") || "first"; } catch { return "first"; } }
/** The depth of this view: a hash parameter (?depth=expert | first) for one view, else the browser's choice.  The
 * depth changes the fields and columns shown, never the tabs. */
function applyDepth(depth) {
  const expert = depth === "expert";
  document.documentElement.dataset.depth = expert ? "expert" : "first";
  ctx.expert = expert;
  if (!ctx.tabBar.childElementCount) drawTabBar();
}
function toggleDepth() {
  const next = ctx.expert ? "first" : "expert";
  try { localStorage.setItem("azure2-depth", next); } catch { /* the depth lasts for this load only */ }
  const r = parseRoute(location.hash);
  if (r.params.depth) {  // a view linked with its depth: the switch wins
    delete r.params.depth;
    const q = new URLSearchParams(r.params).toString();
    location.hash = `#/${[r.tab, ...r.args].map(encodeURIComponent).join("/")}${q ? `?${q}` : ""}`;
  }
  applyDepth(next);
  draw();
}

/** The two tab rows: the main tabs (a group's link opens its first member), and under a group its members. */
function drawTabBar() {
  ctx.tabBar.textContent = "";
  ctx.subBar.textContent = "";
  const link = (t, cls, extra = {}) => {
    const data = { tab: t.tab, view: t.hash, ...(t.when ? { when: t.when } : {}), ...extra };
    const a = h(`a.${cls}`, { href: t.hash, role: "tab", title: t.doc || t.title, data, tabindex: "-1" }, t.title);
    if (t.when) a.hidden = true;  // until showTabs has read the project
    return a;
  };
  mainTabs(shownTabs(ctx.expert)).forEach((m, i) => {
    const a = link(m, "tab", { main: m.title });
    a.prepend(h("span.tab-key", { "aria-hidden": "true" }, i < 9 ? String(i + 1) : "0"));
    if (m.group) a.append(h("span.group-mark", { "aria-hidden": "true" }, " ▾"));
    delete a.dataset.when;
    a.hidden = false;
    ctx.tabBar.appendChild(a);
    if (m.group) {
      ctx.subBar.appendChild(h("div.subgroup", { role: "tablist", "aria-label": m.title, data: { main: m.title },
        hidden: true }, m.members.map((v) => link(v, "tab.sub"))));
    }
  });
}

/** Mark the tab of the route in both rows (the main tab of a group's member too) and show that group's row. */
function markTabs(on) {
  const main = mainTitleOf(on);
  ctx.tabBar.querySelectorAll("a.tab").forEach((a) => {
    const active = a.dataset.main === main;
    a.classList.toggle("active", active);
    a.setAttribute("aria-selected", active ? "true" : "false");
    a.tabIndex = active ? 0 : -1;
  });
  ctx.subBar.querySelectorAll("div.subgroup").forEach((g) => { g.hidden = g.dataset.main !== main; });
  ctx.subBar.hidden = !on?.group;
  ctx.subBar.querySelectorAll("a.tab").forEach((a) => {
    const active = Boolean(on) && a.dataset.view === on.hash;
    a.classList.toggle("active", active);
    a.setAttribute("aria-selected", active ? "true" : "false");
    a.tabIndex = active ? 0 : -1;
  });
  // the main tab of a group opens the member last shown
  if (on?.group) {
    const a = ctx.tabBar.querySelector(`a.tab[data-main="${on.group}"]`);
    if (a) a.href = location.hash || on.hash;
  }
}

/** ← → Home End in a tab row move the focus along it (the ARIA tab pattern); Enter or a click opens the tab. */
function onTabRowKey(e) {
  const row = e.target.closest("div.subgroup") || e.currentTarget;
  const links = [...row.querySelectorAll("a.tab:not([hidden])")];
  const next = tabByArrow(links, links.indexOf(document.activeElement), e.key);
  if (next) { e.preventDefault(); next.focus(); }
}

// -- layout -----------------------------------------------------------------------------------------------------------

const openProject = () => askPath(ctx, "Open a project", "A .azr file of the project's folder (the page then serves "
  + "it; the open one stays as saved).", "other.azr", (path) => ctx.edit(`open ${path}`,
  () => ctx.api.post("/api/project/open", { path })));
const newProject = () => askPath(ctx, "New project from this one", "A new file of the folder with this project's "
  + "pairs and levels and no data (unsaved until saved).", "new.azr", (path) => ctx.edit(`new project ${path}`,
  () => ctx.api.post("/api/project/new", { path, nuclear_part: ctx.project.path })));
const saveAs = () => askPath(ctx, "Save as", "A file of the project's folder; the page then edits that file.",
  "copy.azr", (path) => ctx.edit(`save as ${path}`, () => ctx.api.post("/api/project/save", { path })));

function build() {
  const root = document.getElementById("app");
  root.textContent = "";
  ctx.nameBox = h("span.project-name");
  ctx.unsavedBox = h("span");
  ctx.revBox = h("span.rev");
  ctx.undoButton = button("Undo", "Back to the previous text (Ctrl+Z)", undo);
  ctx.redoButton = button("Redo", "Forward to the next text (Ctrl+Shift+Z)", redo);
  ctx.chi2Box = h("span.chi2", { "aria-live": "polite" }, "χ² …");
  ctx.busy = h("span.busy", { "aria-hidden": "true" });
  ctx.validation = validationMark(ctx);
  ctx.jobMarkEl = jobMark(ctx);
  ctx.jobMark = { poke: () => ctx.jobMarkEl.poke() };
  const calculate = button("Calculate", "Every data set with the model, the residuals and χ², and what changed since "
    + "the last calculation (c)", () => ctx.calculate().catch((e) => ctx.showError(e)), "primary calculate");
  const validate = button("Validate", "The engine validator's messages, each a link to its row (p)",
    () => ctx.validation.open());
  validate.dataset.action = "validate";
  const brand = h("button.brand.quiet", { type: "button", title: "AZURE2: about the program", onclick: toggleAbout },
    h("img", { src: "static/azure2_icon.png", alt: "", width: 22, height: 22 }), "AZURE2");
  ctx.menuBar = menuBar(menus(ctx, { undo, redo, save, saveAs, openProject, newProject, toggleDepth, cycleTheme,
    toggleKeys, toggleAbout, docs: ABOUT.documentation }));
  const header = h("header.top", {}, brand, ctx.menuBar, h("span.project", {}, ctx.nameBox, ctx.unsavedBox,
    ctx.revBox), h("span.spacer"), ctx.undoButton, ctx.redoButton,
  button("Save", "Write the project to its file (Ctrl+S)", save), validate, calculate, ctx.chi2Box, ctx.busy);
  ctx.tabBar = h("nav.tabs", { role: "tablist", "aria-label": "tabs" });
  ctx.subBar = h("nav.tabs.subtabs-row", { "aria-label": "the tabs of the group", hidden: true });
  for (const row of [ctx.tabBar, ctx.subBar]) row.addEventListener("keydown", onTabRowKey);
  ctx.main = h("main.content", { id: "main", tabindex: "-1" });
  ctx.statusLine = h("div.status", { role: "status", "aria-live": "polite" }, "");
  root.append(h("div.top-frame", {}, header, ctx.tabBar, ctx.subBar), ctx.main, ctx.validation.panel,
    h("footer.statusbar", {}, ctx.statusLine, h("span.spacer"), ctx.jobMarkEl, ctx.validation.mark));
}

let drawing = Promise.resolve();
/** Draw the tab of the current hash route (one draw at a time; the newest route wins). */
async function draw() {
  drawing = drawing.then(async () => {
    ctx.route = parseRoute(location.hash);
    if (ctx.route.params.theme) applyTheme(ctx.route.params.theme);
    applyDepth(ctx.route.params.depth || storedDepth());
    const tab = TABS.find((t) => t.id === ctx.route.tab) || TABS.find((t) => t.id === parseRoute(START).tab);
    if (!TABS.some((t) => t.id === ctx.route.tab)) ctx.route = parseRoute(START);
    ctx.section = null;  // the Data tab says which of its groups it drew (its tab is marked by it)
    const mark = () => {
      const on = tabOfRoute(shownTabs(ctx.expert), ctx.route, ctx.section);
      markTabs(on);
      return on;
    };
    ctx.view = mark();
    const panel = h("div.tab-panel", { role: "tabpanel", "aria-label": ctx.view?.title || tab.title });
    // a tab that waits for the engine (the Plot tab's curves) says so instead of an empty page
    const waiting = h("p.loading.tab-loading", {}, `${ctx.view?.title || tab.title}: asking the engine…`);
    ctx.main.textContent = "";
    ctx.main.append(waiting, panel);
    try {
      await tab.module.render(ctx, panel, ctx.route.args);
    } catch (e) {
      ctx.showError(e, panel);
      if (!e.view) throw e;
    } finally {
      waiting.remove();
    }
    ctx.view = mark();
    // a section the address names (#/project?show=history): brought into view
    const show = ctx.route.params.show && panel.querySelector(`[data-show="${ctx.route.params.show}"]`);
    if (show) {
      for (let d = show.closest("details"); d; d = d.parentElement?.closest("details")) d.open = true;
      if (show.tagName === "DETAILS") show.open = true;
      show.scrollIntoView({ block: "start" });
    }
  });
  return drawing;
}

// -- keys ---------------------------------------------------------------------------------------------------------

const toggleKeys = () => keysOverlay(KEYS);
const toggleAbout = () => aboutOverlay(ctx, ABOUT);

function inTextBox(e) {
  const t = e.target;
  return t && (t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.tagName === "SELECT" || t.isContentEditable);
}

function onKey(e) {
  const mod = e.ctrlKey || e.metaKey;
  if (mod && e.key.toLowerCase() === "s") { e.preventDefault(); save(); return; }
  if (e.key === "Escape" && document.getElementById("keys-overlay")) { toggleKeys(); return; }
  if (e.key === "Escape" && document.getElementById("about-overlay")) { toggleAbout(); return; }
  if (e.key === "Escape") ctx.menuBar?.close();
  if (e.key === "F1") { e.preventDefault(); window.open(docsOf(ABOUT.documentation, ctx.view), "_blank", "noopener"); }
  if (inTextBox(e)) return;
  const key = e.key.toLowerCase();
  if (mod && (key === "y" || (key === "z" && e.shiftKey))) { e.preventDefault(); redo(); }
  else if (mod && e.key.toLowerCase() === "z") { e.preventDefault(); undo(); }
  else if (mod || e.altKey) return;
  else if (/^[0-9]$/.test(e.key) || e.key === "[" || e.key === "]") {
    const link = tabByKey([...ctx.tabBar.querySelectorAll("a.tab:not([hidden])")], e.key);  // the main tabs
    if (link) ctx.navigate(link.dataset.view);
  } else if (e.key === "c") ctx.calculate().catch((err) => ctx.showError(err));
  else if (e.key === "p") ctx.validation.open();
  else if (e.key === "e") toggleDepth();
  else if (e.key === "?") toggleKeys();
  else if (e.key === "t") cycleTheme();
  else if (e.key === "/") {
    const f = ctx.main.querySelector("input.table-filter");
    if (f) { e.preventDefault(); f.focus(); }
  }
}

// -- start ----------------------------------------------------------------------------------------------------------

async function start() {
  build();
  applyTheme(storedTheme());
  applyDepth(parseRoute(location.hash).params.depth || storedDepth());
  ctx.api = createApi({
    // a newer revision seen outside an edit of this page: another client's edit (or a reload's first answer)
    onRevision: () => { if (!ctx.editing && ctx.project) ctx.refresh().catch((e) => ctx.showError(e)); },
    onStale: () => ctx.status("another client changed the project first: the page shows the current text", "warn"),
    onActivity: (n) => { ctx.busy.textContent = n ? "working…" : ""; },
  });
  // a job started, stopped or accepted on any tab: the status line's job mark looks at once
  const post = ctx.api.post;
  ctx.api.post = (path, body) => post(path, body).then((answer) => {
    if (/^\/api\/jobs(\/|$)/.test(path)) ctx.jobMark.poke();
    return answer;
  });
  // the self-test needs the API object before the first request
  if (location.hash.startsWith("#selftest")) {
    const { selftest } = await import("./selftest.js");
    selftest(ctx, TABS, draw);
  }
  try {
    [ctx.describe, ctx.schema] = await Promise.all([ctx.api.get("/api/describe"), ctx.api.get("/api/schema")]);
    ctx.project = await ctx.api.get("/api/project");
    ctx.api.rev = ctx.project.rev;
  } catch (e) {
    document.getElementById("app").prepend(errorBox(e.view || { title: "Error", what: String(e) }));
    throw e;
  }
  drawHeader();
  await showTabs();
  ctx.validation.refresh();
  ctx.jobMark.poke();
  const engine = ctx.describe.package?.engine || {};
  ctx.status(`AZURE2 ${engine.version || ""} (${engine.git || ""}) · pyazr ${ctx.describe.package?.pyazr || ""}`);
  window.addEventListener("hashchange", () => { draw(); });
  document.addEventListener("keydown", onKey);
  watchRevision(ctx.api, { isVisible: () => document.visibilityState === "visible" });
  if (location.hash.startsWith("#selftest")) return;
  const firstVisit = overviewDue(location.hash, globalThis.localStorage);
  // the page opened without a route (the printed URL): the address shows the first tab's, so Reload keeps it
  if (!location.hash || location.hash === "#") history.replaceState(null, "", START);
  await draw();  // the first paint does not wait for the evaluation (a session build: about a second)
  if (firstVisit) showOverview(ctx);
  await evaluate();
  await draw();  // the χ² per data set and the plots come with the evaluation
}

start().catch((e) => { console.error(e); });
