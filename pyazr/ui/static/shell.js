// The frame around the tabs (docs/upgrade/web_navigation.md): the menu bar (File, Edit, Configure, View, Help: the
// old setup program's menus, each item a click, never a hover), the status line's two marks (the job that runs or
// ran last, with its progress and Stop; the engine validator's verdict, a click away from its list, each message a
// link to the row it is about), and the overview a first visit shows.  app.js builds the header and the tab rows and
// calls these; the rules (what a mark says, where a message leads) are plain functions node's tests call.

import { pathOf } from "./api.js";
import { isFinal, progressView } from "./fitting.js";
import { FIRST_VIEW, mainTabs } from "./layout.js";
import { formatNumber } from "./tables.js";
import { button, diagnosticLine, h } from "./ui.js";
import { recordRoute } from "./tabs/records.js";

// -- rules (no DOM) ---------------------------------------------------------------------------------------------------

/** The job the status line shows: the running one (the newest that runs), else the newest of all; null without jobs.
 * `jobs` as GET /api/jobs sends them. */
export function jobToMark(jobs) {
  const newest = (list) => list.slice().sort((a, b) => (b.started || 0) - (a.started || 0))[0] || null;
  return newest((jobs || []).filter((j) => !isFinal(j.state))) || newest(jobs || []);
}

/** The words of the job mark: {text, live, canStop, hash} ("fit engine_migrad: running · 152 calls · best 480.3"). */
export function jobMarkView(job) {
  if (!job) return null;
  const p = progressView(job);
  const parts = [`job ${job.id} ${p.title}: ${p.words}`];
  if (p.calls) parts.push(`${p.calls} calls`);
  if (p.best !== null && p.best !== undefined) parts.push(`best ${formatNumber(p.best, 8)}`);
  return { text: parts.join(" · "), live: !p.final, canStop: p.canStop, failed: job.state === "failed",
    hash: `#/fit/job/${encodeURIComponent(job.id)}` };
}

/** The verdict of the engine's validator (GET /api/project/validate) in words: {text, kind, errors, warnings}. */
export function validationView(problems) {
  const list = problems || [];
  const errors = list.filter((d) => (d.severity || "error") === "error").length;
  const warnings = list.length - errors;
  if (!list.length) return { text: "valid", kind: "ok", errors, warnings };
  const said = [errors ? `${errors} error${errors > 1 ? "s" : ""}` : "", warnings
    ? `${warnings} warning${warnings > 1 ? "s" : ""}` : ""].filter(Boolean).join(", ");
  return { text: said, kind: errors ? "error" : "warn", errors, warnings };
}

/** The route of the row a message of the validator is about: its record's (in the tab that shows the record's kind),
 * a level's (the record names a level label), else the Project tab's list. */
export function diagnosticRoute(d, records, schema, levels = []) {
  if (!d.record) return "#/project";
  const r = (records || []).find((x) => x.name === d.record);
  if (r) return recordRoute(schema, r);
  if ((levels || []).some((l) => l.label === d.record)) return `#/nucleus/level/${encodeURIComponent(d.record)}`;
  return "#/project";
}

/** The page of the documentation that describes a tab (the user guide's section of that title). */
export function docsOf(base, view) {
  const slug = String(view?.group || view?.title || "the-tabs").toLowerCase().replace(/[^a-z0-9]+/g, "-")
    .replace(/^-|-$/g, "");
  return `${base}user_guide/web_interface.html#${slug}`;
}

// -- the menu bar ----------------------------------------------------------------------------------------------------

/** The menu bar: `menus` [{title, items: [{label, key?, run, checked?, disabled?} | "-"]}]; a click opens a menu,
 * ↑ ↓ move in it, Enter runs, Esc closes; items are buttons (role menuitem), so the keyboard reaches each. */
export function menuBar(menus) {
  const bar = h("nav.menubar", { role: "menubar", "aria-label": "menu" });
  let open = null;
  const close = () => {
    if (!open) return;
    open.list.hidden = true;
    open.head.setAttribute("aria-expanded", "false");
    open = null;
  };
  for (const m of menus) {
    const head = h("button.menu-head", { type: "button", role: "menuitem", "aria-haspopup": "true",
      "aria-expanded": "false", data: { menu: m.title } }, m.title);
    const list = h("div.menu-list", { role: "menu", "aria-label": m.title, hidden: true });
    const wrap = h("div.menu", {}, head, list);
    const fill = () => {
      list.textContent = "";
      for (const item of (typeof m.items === "function" ? m.items() : m.items)) {
        if (item === "-") { list.appendChild(h("hr", { role: "separator" })); continue; }
        const b = h("button.menu-item", { type: "button", role: item.checked === undefined ? "menuitem"
          : "menuitemcheckbox", "aria-checked": item.checked === undefined ? null : String(Boolean(item.checked)),
        disabled: item.disabled ? true : null, title: item.doc || null, data: { item: item.label } },
        h("span", {}, `${item.checked ? "✓ " : ""}${item.label}`), item.key ? h("kbd", {}, item.key) : null);
        b.addEventListener("click", () => { close(); item.run(); });
        list.appendChild(b);
      }
    };
    const show = () => {
      close();
      fill();
      list.hidden = false;
      head.setAttribute("aria-expanded", "true");
      open = { head, list };
    };
    head.addEventListener("click", () => { if (open?.head === head) close(); else show(); });
    head.addEventListener("keydown", (e) => {
      if (e.key === "ArrowDown") { e.preventDefault(); show(); list.querySelector("button:not(:disabled)")?.focus(); }
    });
    list.addEventListener("keydown", (e) => {
      const items = [...list.querySelectorAll("button:not(:disabled)")];
      const at = items.indexOf(document.activeElement);
      if (e.key === "ArrowDown" || e.key === "ArrowUp") {
        e.preventDefault();
        items[(at + (e.key === "ArrowDown" ? 1 : -1) + items.length) % items.length]?.focus();
      } else if (e.key === "Escape") { e.stopPropagation(); close(); head.focus(); }
    });
    bar.appendChild(wrap);
  }
  document.addEventListener("click", (e) => { if (open && !bar.contains(e.target)) close(); });
  bar.close = close;
  return bar;
}

/** The menus (the old program's File, Configure and Documentation, with Edit and View); `actions` are the page's
 * own (app.js). */
export function menus(ctx, { undo, redo, save, saveAs, openProject, newProject, toggleDepth, cycleTheme,
  toggleKeys, toggleAbout, docs }) {
  const go = (hash) => () => ctx.navigate(hash);
  const doc = (url) => () => window.open(url, "_blank", "noopener");
  return [
    { title: "File", items: [
      { label: "Open…", doc: "Open another project of the folder (POST /api/project/open)", run: openProject },
      { label: "New project from this one…", doc: "A new file with this project's pairs and levels, no data "
        + "(POST /api/project/new)", run: newProject },
      "-",
      { label: "Save", key: "Ctrl+S", run: save },
      { label: "Save as…", doc: "Write the project to another file of the folder", run: saveAs },
      { label: "Write the output files", doc: "The program's output files (AZUREOut_*, chiSquared.out) into the "
        + "output directory", run: () => ctx.edit("write the program's output files",
        () => ctx.api.post("/api/outputs", {})).catch(() => {}) },
      "-",
      { label: "New project from an ENDF-6 file…", run: go("#/endf") },
      { label: "Export as ENDF-6…", run: go("#/endf/export") },
      "-",
      { label: "About AZURE2…", run: toggleAbout },
    ] },
    { title: "Edit", items: () => [
      { label: "Undo", key: "Ctrl+Z", disabled: !ctx.project?.can_undo, run: undo },
      { label: "Redo", key: "Ctrl+Y", disabled: !ctx.project?.can_redo, run: redo },
      "-",
      { label: "History…", run: go("#/project?show=history") },
      { label: "The project's text…", run: go("#/project?show=text") },
    ] },
    { title: "Configure", items: [
      { label: "Run settings…", doc: "Formalism, checks, directories, runtime options: the <config> entries",
        run: go("#/project?show=settings") },
      { label: "Session options…", run: go("#/project?show=options") },
      { label: "Prior on the parameters…", run: go("#/prior") },
      { label: "Engine and threads…", run: go("#/project?show=engine") },
    ] },
    { title: "View", items: () => [
      { label: "Every field and column (expert)", key: "e", checked: ctx.expert, run: toggleDepth },
      { label: `Theme: ${ctx.theme || "system"} (next)`, key: "t", run: cycleTheme },
      "-",
      { label: "The engine validator's messages", key: "p", run: () => ctx.validation.open() },
      { label: "The job shown in the status line", run: () => ctx.navigate(ctx.jobMarkEl?.querySelector("a")
        ?.getAttribute("href") || "#/fit") },
    ] },
    { title: "Help", items: () => [
      { label: "Overview of the tabs", run: () => showOverview(ctx) },
      { label: "Documentation for this tab", key: "F1", run: doc(docsOf(docs, ctx.view)) },
      { label: "Documentation", run: doc(docs) },
      { label: "Keys", key: "?", run: toggleKeys },
      "-",
      { label: "About AZURE2…", run: toggleAbout },
    ] },
  ];
}


// -- the status line's marks -------------------------------------------------------------------------------------

/** The job mark of the status line, kept up to date: GET /api/jobs every 2 s while a job runs, every 10 s otherwise
 * (a job another client started), and at once after `poke()`.  Stop sends POST /api/jobs/{id}/stop. */
export function jobMark(ctx) {
  const text = h("a.job-text", { href: "#/fit", title: "the job's panel (its progress, result, Accept)" });
  const stop = button("Stop", "Stop the job before its next objective call; its best point is kept", async () => {
    if (!shown) return;
    stop.disabled = true;
    try { await ctx.api.post(`${pathOf("/api/jobs", shown.id)}/stop`, {}); } catch (e) { ctx.showError(e); }
    poke();
  }, "danger");
  stop.dataset.action = "status-stop";
  const el = h("span.job-mark", { "aria-live": "polite", hidden: true }, text, stop);
  let shown = null, timer = null, last = null;
  async function poke() {
    clearTimeout(timer);
    let jobs = [];
    try {
      jobs = await ctx.api.request("GET", "/api/jobs", undefined, { quiet: true, background: true });
    } catch { /* the next poll tries again */ }
    shown = jobToMark(jobs);
    const v = jobMarkView(shown);
    el.hidden = !v;
    if (v) {
      text.textContent = v.text;
      text.href = v.hash;
      el.className = `job-mark ${v.live ? "live" : v.failed ? "failed" : "done"}`;
      stop.hidden = !v.live;
      stop.disabled = !v.canStop;
      // a job that ended since the last look: the project may have a new best point to accept
      if (last && last.id === shown.id && !isFinal(last.state) && isFinal(shown.state)) {
        ctx.status(`job ${shown.id}: ${v.text.split(": ").slice(1).join(": ")} (the Calculate and Fit tab shows it)`,
          shown.state === "failed" ? "error" : "ok");
      }
    }
    last = shown;
    if (!document.hidden || v?.live) timer = setTimeout(poke, v?.live ? 2000 : 10000);
    else timer = setTimeout(poke, 10000);
  }
  el.poke = poke;
  return el;
}

/** The validator's mark of the status line and its list: `refresh()` asks GET /api/project/validate (after every
 * change of revision), a click opens the list above the status line, each message a link to its row. */
export function validationMark(ctx) {
  const mark = h("button.validation", { type: "button", title: "The engine's validator: the list of its messages",
    "aria-expanded": "false" }, "…");
  mark.dataset.action = "validation";
  const panel = h("div.diag-panel", { role: "dialog", "aria-label": "problems", hidden: true });
  let problems = [];
  async function refresh() {
    try {
      problems = await ctx.api.request("GET", "/api/project/validate", undefined, { quiet: true, background: true });
    } catch (e) {
      problems = [{ code: "", severity: "error", text: e.view?.what || String(e) }];
    }
    const v = validationView(problems);
    mark.textContent = v.kind === "ok" ? "✓ valid" : `${v.kind === "error" ? "✕" : "⚠"} ${v.text}`;
    mark.className = `validation ${v.kind}`;
    mark.dataset.errors = String(v.errors);
    if (!panel.hidden) await fill();
    return v;
  }
  async function fill() {
    const [records, levels] = await Promise.all([ctx.api.get("/api/project/records", { quiet: true }).catch(() => []),
      ctx.api.get("/api/project/levels", { quiet: true }).catch(() => [])]);
    panel.textContent = "";
    panel.append(h("div.diag-head", {}, h("strong", {}, `The engine's validator: ${validationView(problems).text}`),
      h("span.spacer"), button("Close", "Close the list (Esc)", () => toggle(false), "quiet")));
    if (!problems.length) panel.appendChild(h("p.ok", {}, "The engine's validator reports nothing."));
    for (const d of problems) {
      const route = diagnosticRoute(d, records, ctx.schema, levels);
      const line = diagnosticLine({ ...d, record: "" });
      if (d.record) {
        line.prepend(h("a.diag-link", { href: route, title: "go to the row it is about" }, d.record), " ");
      }
      panel.appendChild(line);
    }
  }
  async function toggle(on = panel.hidden) {
    panel.hidden = !on;
    mark.setAttribute("aria-expanded", String(on));
    if (on) { await fill(); panel.querySelector("a, button")?.focus(); }
  }
  mark.addEventListener("click", () => toggle());
  panel.addEventListener("click", (e) => { if (e.target.closest("a.diag-link")) toggle(false); });
  panel.addEventListener("keydown", (e) => { if (e.key === "Escape") { e.stopPropagation(); toggle(false); } });
  return { mark, panel, refresh, open: async () => { await refresh(); await toggle(true); } };
}

// -- the overview ----------------------------------------------------------------------------------------------------

const SEEN = "azure2-overview-seen";
/** Whether the overview opens by itself: on the first visit of this browser, when the page is opened without a route
 * (the printed URL; a link to a place goes there at once). */
export function overviewDue(hash, storage) {
  if (hash && hash !== "#" && hash !== "#/") return false;
  try { return !storage?.getItem(SEEN); } catch { return false; }
}

/** The overview: the tabs in the order of the work, what each is for, where the menus and the status line are. */
export function showOverview(ctx, { onClose } = {}) {
  document.getElementById("overview")?.remove();
  const close = () => {
    try { localStorage.setItem(SEEN, "1"); } catch { /* shown again next time */ }
    box.remove();
    onClose?.();
  };
  const rows = mainTabs(FIRST_VIEW).map((m, i) => h("li", {}, h("a", { href: m.hash, onclick: () => close() },
    m.title), m.group ? ` (${m.members.map((x) => x.title).join(", ")})` : "", ": ",
  m.members.map((x) => x.doc).join("; "), h("kbd.key-of", {}, i < 9 ? String(i + 1) : "0")));
  const box = h("div.overlay", { id: "overview", role: "dialog", "aria-label": "overview", tabindex: "-1" },
    h("div.overlay-box.overview", {}, h("h2", {}, "AZURE2: the tabs in the order of the work"),
      h("ol.overview-tabs", {}, rows),
      h("p", {}, "The first six tabs are the old setup program's, in its order. Every tab is one click away in the "
        + "row at the top (or its number key); its table is on the left with Add, Edit and Remove at the table, the "
        + "selected row's detail beside it. The menus above hold the file (File), Undo and Redo (Edit), the run's "
        + "settings (Configure), the depth and the theme (View) and this overview (Help)."),
      h("p", {}, "The status line at the bottom says what the engine said last, the job that runs (with Stop) and "
        + "whether the project is valid: click it for the list, each message a link to its row. The address follows "
        + "the tab and the row, so Back, Forward and Reload keep your place."),
      h("div.form-actions", {}, button(`Start with ${FIRST_VIEW[0].title}`, "Close and open the first tab", () => {
        close();
        ctx.navigate(FIRST_VIEW[0].hash);
      }, "primary"), button("Close", "Close (Esc); Help › Overview opens it again", close))));
  box.addEventListener("click", (e) => { if (e.target === box) close(); });
  box.addEventListener("keydown", (e) => { if (e.key === "Escape") { e.stopPropagation(); close(); } });
  document.body.appendChild(box);
  box.focus();
  return box;
}

// -- small dialogs ---------------------------------------------------------------------------------------------------

/** A dialog with one text box (the folder's .azr files offered), OK and Cancel; `run(path)` refused keeps it open
 * with the refusal. */
export async function askPath(ctx, title, doc, placeholder, run) {
  document.getElementById("ask-overlay")?.remove();
  const input = h("input.mono", { type: "text", size: 40, placeholder, "aria-label": title, list: "ask-files" });
  const list = h("datalist", { id: "ask-files" });
  const msg = h("div.messages");
  const ok = async () => {
    msg.textContent = "";
    try { await run(input.value.trim()); box.remove(); } catch (e) { ctx.showError(e, msg); }
  };
  const box = h("div.overlay", { id: "ask-overlay", role: "dialog", "aria-label": title, tabindex: "-1" },
    h("div.overlay-box", {}, h("h2", {}, title), h("p.hint", {}, doc), input, list, msg,
      h("div.form-actions", {}, button("OK", "", ok, "primary"), button("Cancel", "Close (Esc)", () => box.remove()))));
  input.addEventListener("keydown", (e) => { if (e.key === "Enter") ok(); });
  box.addEventListener("keydown", (e) => { if (e.key === "Escape") { e.stopPropagation(); box.remove(); } });
  document.body.appendChild(box);
  input.focus();
  ctx.api.get("/api/files", { quiet: true }).then((r) => r.files.filter((f) => f.kind === "project")
    .forEach((f) => list.appendChild(h("option", { value: f.path })))).catch(() => {});
}

/** The overlay of the keys (app.js KEYS), opened and closed by "?". */
export function keysOverlay(keys) {
  const old = document.getElementById("keys-overlay");
  if (old) { old.remove(); return; }
  const rows = keys.map(([k, what]) => h("tr", {}, h("th", {}, h("kbd", {}, k)), h("td", {}, what)));
  const box = h("div.overlay", { id: "keys-overlay", role: "dialog", "aria-label": "keys", tabindex: "-1" },
    h("div.overlay-box", {}, h("h2", {}, "Keys"), h("table.keys", {}, h("tbody", {}, rows)),
      button("Close", "Close (Esc)", () => keysOverlay(keys))));
  box.addEventListener("click", (e) => { if (e.target === box) keysOverlay(keys); });
  document.body.appendChild(box);
  box.focus();
}

/** The About box: name and version of the program and of pyazr, the reference to cite, the documentation. */
export function aboutOverlay(ctx, about) {
  const old = document.getElementById("about-overlay");
  if (old) { old.remove(); return; }
  const pkg = ctx.describe?.package || {};
  const engine = pkg.engine || {};
  const box = h("div.overlay", { id: "about-overlay", role: "dialog", "aria-label": "about", tabindex: "-1" },
    h("div.overlay-box", {}, h("h2.brand", {}, h("img", { src: "static/azure2_icon.png", alt: "", width: 32,
      height: 32 }), "AZURE2"),
    h("p", {}, "A multi-level, multi-channel R-matrix code for nuclear reactions and scattering."),
    h("dl.about", {},
      h("dt", {}, "Version"), h("dd.mono", {}, `AZURE2 ${engine.version ?? ""} (${engine.git ?? ""}), pyazr `
        + `${pkg.pyazr ?? ""}`),
      h("dt", {}, "Cite"), h("dd", {}, about.citation),
      h("dt", {}, "Documentation"), h("dd", {}, h("a", { href: about.documentation, target: "_blank",
        rel: "noopener" }, about.documentation))),
    button("Close", "Close (Esc)", () => aboutOverlay(ctx, about))));
  box.addEventListener("click", (e) => { if (e.target === box) aboutOverlay(ctx, about); });
  document.body.appendChild(box);
  box.focus();
}
