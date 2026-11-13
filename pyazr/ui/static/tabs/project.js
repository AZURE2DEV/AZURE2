// The Project tab (docs/upgrade/web_ui.md section 6.1): the file, the <config> entries (edited through the config
// route), the session's options as the server reports them (read-only until they are <config> keys), the problems
// of the engine's validator with links to what they name, the history with Undo and Redo, the text (read-only, with
// "replace text"), and the engine's line.

import { formatNumber, renderTable } from "../tables.js";
import { button, diagnosticLine, h, section } from "../ui.js";
import { recordRoute } from "./records.js";

export async function render(ctx, el) {
  const [config, problems, history, records] = await Promise.all([ctx.api.get("/api/project/config"),
    ctx.api.get("/api/project/validate"), ctx.api.get("/api/project/history"), ctx.api.get("/api/project/records")]);
  const msg = h("div.messages");
  el.appendChild(msg);
  el.appendChild(fileSection(ctx, msg));
  el.appendChild(problemsSection(ctx, problems, records));
  const columns = h("div.two-columns");
  columns.append(configSection(ctx, config, msg), optionsSection(ctx));
  el.appendChild(columns);
  el.appendChild(historySection(ctx, history));
  el.appendChild(textSection(ctx, msg));
  el.appendChild(engineSection(ctx));
}

function fileSection(ctx, msg) {
  const p = ctx.project;
  const dl = h("dl.facts", {},
    h("dt", {}, "File"), h("dd.mono", {}, p.path),
    h("dt", {}, "Working copy"),
    h("dd.mono", {}, `${p.working_file} (the current text while it differs from the file)`),
    h("dt", {}, "Root"), h("dd.mono", {}, `${p.root} (files may be read and written inside it)`),
    h("dt", {}, "Revision"), h("dd", {}, `${p.rev}${p.unsaved ? ", unsaved edits" : ", as saved"}`));
  const s = section("File", dl);
  if (p.crash_copy) {
    s.appendChild(h("div.warn-box", {}, `Unsaved edits of an earlier server were found (${p.crash_copy.path}). `,
      button("Restore them", "Take the crash copy as an edit (Undo goes back)", () => ctx.edit("restore",
        () => ctx.api.post("/api/project/restore", {}), msg).catch(() => {}))));
  }
  return s;
}

/** The <config> entries: one text box per key; Apply sends the changed ones through PATCH /api/project/config. */
function configSection(ctx, config, msg) {
  const rows = Object.entries(config);
  const inputs = [];
  const tbody = h("tbody");
  for (const [k, v] of rows) {
    const input = h("input", { type: "text", value: String(v), size: 22, "aria-label": k });
    inputs.push([k, String(v), input]);
    tbody.appendChild(h("tr", {}, h("th", {}, k), h("td", {}, input)));
  }
  const newKey = h("input", { type: "text", placeholder: "key", size: 16, "aria-label": "new key" });
  const newValue = h("input", { type: "text", placeholder: "value", size: 22, "aria-label": "new value" });
  tbody.appendChild(h("tr", {}, h("th", {}, newKey), h("td", {}, newValue)));
  const apply = () => {
    const entries = Object.fromEntries(inputs.filter(([, v, i]) => i.value !== v)
      .map(([k, , i]) => [k, i.value.trim() === "" ? null : i.value.trim()]));
    if (newKey.value.trim()) entries[newKey.value.trim()] = newValue.value.trim();
    if (!Object.keys(entries).length) { ctx.status("<config>: nothing changed"); return; }
    ctx.edit("set <config>", () => ctx.api.patch("/api/project/config", { entries }), msg).catch(() => {});
  };
  const table = h("table.kv", {}, tbody);
  table.addEventListener("keydown", (e) => { if (e.key === "Enter") { e.preventDefault(); apply(); } });
  return section("<config>", h("p.hint", {}, "An empty value removes the key; Enter applies."), table,
    h("div.form-actions", {}, button("Apply", "set_config with the changed entries (Enter)", apply, "primary")));
}

/** The session's options as GET /api/describe reports them: keyword, value in use, its source, default, flag, doc. */
function optionsSection(ctx) {
  const opts = ctx.describe?.session_options || {};
  const rows = Object.entries(opts).map(([k, o]) => ({ name: k, value: o.value, source: o.source, default: o.default,
    flag: [o.flag, o.negation].filter(Boolean).join(" / "), doc: o.doc, type: o.type }));
  const table = renderTable(document, { rows, withFilter: false, label: "session options",
    columns: [{ key: "name", title: "option" }, { key: "value", title: "in use" }, { key: "source", title: "from" },
      { key: "default", title: "default" }, { key: "flag", title: "program flag", wrap: true },
      { key: "doc", title: "what it does", wrap: true }] });
  return section("Session options", h("p.hint", {}, "The options of the calculation are <config> keys of the "
    + "project under the same name (value 0 or 1): set one in the <config> table above; an empty value returns it to "
    + "the default."), table.element);
}

/** The problems of Project.validate(), each with a link to the record or level it names. */
function problemsSection(ctx, problems, records) {
  const s = section(`Problems (${problems.length})`);
  if (!problems.length) s.appendChild(h("p.ok", {}, "The engine's validator reports nothing."));
  const byName = Object.fromEntries(records.map((r) => [r.name, r]));
  for (const d of problems) {
    s.appendChild(diagnosticLine(d, (x) => {
      const r = byName[x.record];
      ctx.navigate(r ? recordRoute(ctx.schema, r) : `#/nucleus/level/${encodeURIComponent(x.record)}`);
    }));
  }
  return s;
}

function historySection(ctx, history) {
  const rows = history.slice().reverse().map((e) => ({ ...e, when: new Date(e.time * 1000).toLocaleTimeString(),
    now: e.current ? "◀ current" : "" }));
  const table = renderTable(document, { rows, rowKey: (r) => String(r.rev), withFilter: false, label: "history",
    columns: [{ key: "rev", title: "rev", numeric: true }, { key: "when", title: "time" },
      { key: "what", title: "edit" },
      { key: "now", title: "" }] });
  return section("History", h("div.form-actions", {},
    button("Undo", "Back to the previous text (Ctrl+Z)", () => ctx.edit("undo", () => ctx.api.post("/api/project/undo",
      {})).catch(() => {})),
    button("Redo", "Forward to the next text (Ctrl+Shift+Z)", () => ctx.edit("redo",
      () => ctx.api.post("/api/project/redo", {})).catch(() => {}))), table.element);
}

/** The project's text, read-only; "Replace text" opens it for an edit of the whole (validated by the engine). */
function textSection(ctx, msg) {
  const pre = h("pre.project-text");
  const details = h("details", {}, h("summary", {}, "Show the text"), pre);
  details.addEventListener("toggle", () => {
    if (details.open && !pre.textContent) ctx.api.text("/api/project/text").then((t) => { pre.textContent = t; });
  });
  const editor = h("div");
  const replace = button("Replace text…", "Edit the whole text; the engine validates it first", async () => {
    const text = await ctx.api.text("/api/project/text");
    const area = h("textarea.project-text", { rows: 24, spellcheck: "false", "aria-label": "the project's text" });
    area.value = text;
    editor.textContent = "";
    editor.append(area, h("div.form-actions", {},
      button("Replace", "PUT /api/project/text", () => ctx.edit("replace text",
        () => ctx.api.put("/api/project/text", { text: area.value }), msg).catch(() => {}), "primary"),
      button("Cancel", "Keep the text", () => { editor.textContent = ""; })));
    area.focus();
  });
  return section("Text", details, h("div.form-actions", {}, replace), editor);
}

/** The engine: versions, build options, threads, Coulomb table, stopping and mass data, the compute budget. */
function engineSection(ctx) {
  const pkg = ctx.describe?.package || {};
  const e = pkg.engine || {};
  const b = pkg.build || {};
  const ct = pkg.coulomb_table || {};
  const budget = ctx.describe?.budget || {};
  const dl = h("dl.facts", {},
    h("dt", {}, "Engine"), h("dd", {}, `AZURE2 ${e.version ?? ""} (${e.git ?? ""}), pyazr ${pkg.pyazr ?? ""}, `
      + `Python ${ctx.describe?.python ?? ""}`),
    h("dt", {}, "Build"), h("dd", {}, Object.entries(b).map(([k, v]) => `${k} ${v ? "on" : "off"}`).join(" · ")),
    h("dt", {}, "Threads"), h("dd", {}, `budget ${budget.threads ?? "?"} (job slots ${budget.job_slots ?? "?"}); `
      + `engine ${formatNumber(pkg.threads?.max_threads)}`),
    h("dt", {}, "Coulomb table"),
    h("dd", {}, ct.in_use === "table" ? `${ct.path} (${ct.identity})` : (ct.in_use ?? "")),
    h("dt", {}, "Stopping"), h("dd", {}, pkg.stopping?.implementation ?? ""),
    h("dt", {}, "Masses"), h("dd", {}, `${pkg.masses?.version ?? ""} ${pkg.masses?.file ?? ""}`));
  return section("Engine", dl);
}
