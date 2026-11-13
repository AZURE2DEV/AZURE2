// The page of AZURE2's web UI (docs/upgrade/web_ui.md sections 5 to 7): the header, the tabs and their hash routes
// (#/data/s02_diff090), the keys and their overlay, the status line, and the self-test (#selftest).
//
// The page holds no state of the project.  After every change of revision (an edit here, Undo, Redo, another
// client's edit seen by the poll) it re-reads the project's state, asks the server for the evaluation of the
// current text (the header's χ² is the engine's number, as the server sends it) and draws the open tab again.

import { createApi, watchRevision } from "./api.js";
import { button, errorBox, h } from "./ui.js";
import { formatNumber } from "./tables.js";
import * as projectTab from "./tabs/project.js";
import * as nucleusTab from "./tabs/nucleus.js";
import * as dataTab from "./tabs/data.js";

/** The tabs: id (the hash route's first word), title, key, the module that draws it or the step that brings it. */
export const TABS = [
  { id: "project", title: "Project", module: projectTab },
  { id: "nucleus", title: "Nucleus", module: nucleusTab },
  { id: "data", title: "Data", module: dataTab },
  { id: "fit", title: "Fit", later: "the fitting step: the parameter table, recipe fits as jobs with progress and "
    + "Stop, MINOS, the sampler. Until then a fit runs through the API (POST /api/jobs) or pyazr." },
  { id: "results", title: "Results", later: "the results step: plots of every data set and extrapolation, residuals, "
    + "the χ² table, the program's output files. The Data tab shows one data set's preview plot today." },
  { id: "physics", title: "Physics", later: "the physics step: one generated form and plot per physics function of "
    + "the engine (Coulomb functions, penetrabilities, stopping, the transfer library)." },
  { id: "thm", title: "THM", later: "the physics and THM step: THM model options, experiments, channels and "
    + "diagnostics. The <thm> section is edited through the API today (PATCH /api/project/thm)." },
];

/** The keys of the page, shown by "?" (the overlay lists exactly this table). */
export const KEYS = [
  ["1 … 7", "switch to the tab of that number"],
  ["Ctrl+Z", "undo the last edit (outside a text box)"],
  ["Ctrl+Shift+Z, Ctrl+Y", "redo"],
  ["Ctrl+S", "save the project to its file"],
  ["/", "filter the current table"],
  ["↑ ↓ (j k), Enter", "move through a table's rows, open the row"],
  ["Ctrl+Enter", "apply the form being edited"],
  ["Esc", "revert the form being edited; close this overlay"],
  ["t", "switch the theme: system, light, dark"],
  ["?", "this list"],
];

/** A hash route: "#/data/s02_diff090?theme=dark" -> {tab, args: ["s02_diff090"], params: {theme: "dark"}}. */
export function parseRoute(hash) {
  const [path, query = ""] = String(hash || "").replace(/^#\/?/, "").split("?");
  const words = path.split("/").filter(Boolean).map(decodeURIComponent);
  const params = Object.fromEntries(new URLSearchParams(query));
  return { tab: words[0] || "project", args: words.slice(1), params };
}

const ctx = {
  api: null, schema: null, kinds: null, describe: null, project: null, evaluation: null, evalError: null,
  route: parseRoute(""), statusLine: null,
};

/** Show a message of the engine or the page in the status line. */
ctx.status = (text, kind = "info") => {
  if (!ctx.statusLine) return;
  ctx.statusLine.textContent = text;
  ctx.statusLine.className = `status ${kind}`;
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
    else ctx.status(`${what}: done`, "ok");
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
  await evaluate();
  await draw();
};

async function evaluate() {
  ctx.chi2Box.textContent = "χ² …";
  try {
    // a project the engine cannot evaluate (no data sets, a refused set-up) is shown, not a failure of the page
    ctx.evaluation = await ctx.api.request("POST", "/api/evaluate", {}, { quiet: true });
    ctx.evalError = null;
    const w = ctx.evaluation.warnings || [];
    if (w.length) ctx.status(`engine: ${w.map((x) => x.message || x).join("; ")}`, "warn");
  } catch (e) {
    ctx.evaluation = null;
    ctx.evalError = e.view || { what: String(e) };
    ctx.status(`no evaluation: ${ctx.evalError.what}`, "warn");
  }
  drawChi2();
}

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
  if (ctx.themeButton) ctx.themeButton.textContent = `theme: ${theme}`;
}
function cycleTheme() {
  const next = THEMES[(THEMES.indexOf(storedTheme()) + 1) % THEMES.length];
  try { localStorage.setItem("azure2-theme", next); } catch { /* the theme lasts for this load only */ }
  applyTheme(next);
  draw();
}

// -- layout -----------------------------------------------------------------------------------------------------------

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
  ctx.themeButton = button("theme", "Switch the theme: system, light, dark (t)", cycleTheme, "quiet");
  const header = h("header.top", {}, h("span.brand", {}, "AZURE2"), ctx.nameBox, ctx.unsavedBox, ctx.revBox,
    ctx.undoButton, ctx.redoButton, button("Save", "Write the project to its file (Ctrl+S)", save), ctx.chi2Box,
    ctx.busy, h("span.spacer"), ctx.themeButton, button("?", "The keys of the page (?)", toggleKeys, "quiet"));
  ctx.tabBar = h("nav.tabs", { role: "tablist", "aria-label": "tabs" });
  TABS.forEach((t, i) => {
    const a = h("a.tab", { href: `#/${t.id}`, role: "tab", title: `${t.title} (${i + 1})`, data: { tab: t.id } },
      t.title, t.later ? h("span.later-mark", { title: "comes with a later step" }, " ·") : null);
    ctx.tabBar.appendChild(a);
  });
  ctx.main = h("main.content", { id: "main", tabindex: "-1" });
  ctx.statusLine = h("div.status", { role: "status", "aria-live": "polite" }, "");
  root.append(header, ctx.tabBar, ctx.main, h("footer.statusbar", {}, ctx.statusLine));
}

let drawing = Promise.resolve();
/** Draw the tab of the current hash route (one draw at a time; the newest route wins). */
async function draw() {
  drawing = drawing.then(async () => {
    ctx.route = parseRoute(location.hash);
    if (ctx.route.params.theme) applyTheme(ctx.route.params.theme);
    const tab = TABS.find((t) => t.id === ctx.route.tab) || TABS[0];
    ctx.tabBar.querySelectorAll("a.tab").forEach((a) => {
      const on = a.dataset.tab === tab.id;
      a.classList.toggle("active", on);
      a.setAttribute("aria-selected", on ? "true" : "false");
    });
    const panel = h("div.tab-panel", { role: "tabpanel", "aria-label": tab.title });
    ctx.main.textContent = "";
    ctx.main.appendChild(panel);
    if (tab.later) {
      panel.appendChild(h("section.panel.later", {}, h("h2", {}, tab.title),
        h("p", {}, `This tab comes with ${tab.later}`)));
      return;
    }
    try {
      await tab.module.render(ctx, panel, ctx.route.args);
    } catch (e) {
      ctx.showError(e, panel);
      if (!e.view) throw e;
    }
  });
  return drawing;
}

// -- keys ---------------------------------------------------------------------------------------------------------

function toggleKeys() {
  const old = document.getElementById("keys-overlay");
  if (old) { old.remove(); return; }
  const rows = KEYS.map(([k, what]) => h("tr", {}, h("th", {}, h("kbd", {}, k)), h("td", {}, what)));
  const box = h("div.overlay", { id: "keys-overlay", role: "dialog", "aria-label": "keys", tabindex: "-1" },
    h("div.overlay-box", {}, h("h2", {}, "Keys"), h("table.keys", {}, h("tbody", {}, rows)),
      button("Close", "Close (Esc)", toggleKeys)));
  box.addEventListener("click", (e) => { if (e.target === box) toggleKeys(); });
  document.body.appendChild(box);
  box.focus();
}

function inTextBox(e) {
  const t = e.target;
  return t && (t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.tagName === "SELECT" || t.isContentEditable);
}

function onKey(e) {
  const mod = e.ctrlKey || e.metaKey;
  if (mod && e.key.toLowerCase() === "s") { e.preventDefault(); save(); return; }
  if (e.key === "Escape" && document.getElementById("keys-overlay")) { toggleKeys(); return; }
  if (inTextBox(e)) return;
  const key = e.key.toLowerCase();
  if (mod && (key === "y" || (key === "z" && e.shiftKey))) { e.preventDefault(); redo(); }
  else if (mod && e.key.toLowerCase() === "z") { e.preventDefault(); undo(); }
  else if (mod || e.altKey) return;
  else if (/^[1-7]$/.test(e.key)) ctx.navigate(`#/${TABS[+e.key - 1].id}`);
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
  ctx.api = createApi({
    // a newer revision seen outside an edit of this page: another client's edit (or a reload's first answer)
    onRevision: () => { if (!ctx.editing && ctx.project) ctx.refresh().catch((e) => ctx.showError(e)); },
    onStale: () => ctx.status("another client changed the project first: the page shows the current text", "warn"),
    onActivity: (n) => { ctx.busy.textContent = n ? "working…" : ""; },
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
  const engine = ctx.describe.package?.engine || {};
  ctx.status(`AZURE2 ${engine.version || ""} (${engine.git || ""}) · pyazr ${ctx.describe.package?.pyazr || ""}`);
  window.addEventListener("hashchange", () => { draw(); });
  document.addEventListener("keydown", onKey);
  watchRevision(ctx.api, { isVisible: () => document.visibilityState === "visible" });
  if (location.hash.startsWith("#selftest")) return;
  await draw();  // the first paint does not wait for the evaluation (a session build: about a second)
  await evaluate();
  await draw();  // the χ² per data set and the plots come with the evaluation
}

start().catch((e) => { console.error(e); });
