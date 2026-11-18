// The Project tab (docs/upgrade/web_ui.md section 6.1): the file, the <config> entries (edited through the config
// route), the session's options as the server reports them (read-only until they are <config> keys), the problems
// of the engine's validator with links to what they name, the history with Undo and Redo, the text (read-only, with
// "replace text"), and the engine's line.

import { formatNumber, renderTable } from "../tables.js";
import { button, diagnosticLine, h, section } from "../ui.js";
import { FIRST_VIEW } from "../layout.js";
import { recordRoute } from "./records.js";
import { settingsSection } from "./settings.js";

/**
 * The next step of a project, from what it has (the first view's guide): {text, hash, label} of the first thing
 * missing in the old program's order (pairs, levels, data, a calculation), or the way on when nothing is.  `counts`:
 * {pairs, levels, data} (data: the records of a kind with a quantity); `evalError`: the
 * engine's refusal to evaluate, if any.
 */
export function nextStep(counts, evalError, firstView = FIRST_VIEW) {
  const view = (title) => firstView.find((v) => v.title === title)?.hash;
  if (!counts.pairs) return { text: "Add the particle pairs: the entrance pair and every exit pair.",
    label: "Particle Pairs", hash: view("Particle Pairs") };
  if (!counts.levels) return { text: "Add the levels of the compound nucleus; their channels come from the pairs.",
    label: "Levels and Channels", hash: view("Levels and Channels") };
  if (!counts.data) return { text: "Add the data: a segment from data (a data file) or without data (a grid).",
    label: "Segments", hash: view("Segments") };
  if (evalError) return { text: `The engine does not calculate the project yet: ${evalError.what}`,
    label: "Problems", hash: "#/project?show=problems" };
  return { text: "Calculate shows every data set against the model; Fit adjusts the free parameters.",
    label: "Calculate", hash: "calculate" };
}

export async function render(ctx, el) {
  if (!ctx.kinds) ctx.kinds = await ctx.api.get("/api/kinds");
  const [config, problems, history, records] = await Promise.all([ctx.api.get("/api/project/config"),
    ctx.api.get("/api/project/validate"), ctx.api.get("/api/project/history"), ctx.api.get("/api/project/records")]);
  const msg = h("div.messages");
  el.appendChild(msg);
  el.appendChild(fileSection(ctx, msg));
  const send = (entries) => ctx.edit("set <config>", () => ctx.api.patch("/api/project/config", { entries }), msg)
    .catch(() => {});
  const settings = settingsSection(ctx, ctx.kinds.config?.fields || {}, config, send, configSection(ctx, config, msg));
  // the places the menus open (#/project?show=settings: Configure › Run settings, ...)
  const tag = (name, part) => { part.dataset.show = name; return part; };
  tag("settings", settings);
  const optionsSection = (c) => tag("options", optionsPart(c));
  const historySection = (c, x) => tag("history", historyPart(c, x));
  const textSection = (c, m) => tag("text", textPart(c, m));
  const engineSection = (c) => tag("engine", enginePart(c));
  const problemsSection = (c, p, r) => tag("problems", problemsPart(c, p, r));
  if (!ctx.expert) {
    // the first view: the next step, the problems, the run's settings; the rest one click away
    el.appendChild(await nextSection(ctx, records));
    el.appendChild(problemsSection(ctx, problems, records));
    el.appendChild(settings);
    el.appendChild(h("details.more", {}, h("summary", {}, "More: session options, history, the project's text, "
      + "the engine"), optionsSection(ctx), historySection(ctx, history), textSection(ctx, msg),
    engineSection(ctx)));
    return;
  }
  el.appendChild(problemsSection(ctx, problems, records));
  el.appendChild(settings);
  el.appendChild(optionsSection(ctx));
  el.appendChild(historySection(ctx, history));
  el.appendChild(textSection(ctx, msg));
  el.appendChild(engineSection(ctx));
}

/** The next step (nextStep) with the button that takes it. */
async function nextSection(ctx, records) {
  const [pairs, levels] = await Promise.all([ctx.api.get("/api/project/pairs"), ctx.api.get("/api/project/levels")]);
  // the data: records of a kind with a quantity (measured values, values on a grid, model data)
  const dataKinds = Object.entries(ctx.schema.records).filter(([, r]) => r.fields?.quantity).map(([k]) => k);
  const counts = { pairs: pairs.length, levels: levels.length,
    data: records.filter((r) => dataKinds.includes(r.kind)).length };
  const step = nextStep(counts, ctx.evalError);
  const go = step.hash === "calculate" ? () => ctx.calculate().catch((e) => ctx.showError(e))
    : step.hash ? () => ctx.navigate(step.hash) : null;
  return section("Next step", h("p", {}, step.text, " ", go ? button(step.label, step.text, go, "primary") : null),
    h("p.hint", {}, `This project: ${counts.pairs} pairs, ${counts.levels} levels, ${counts.data} segments (from `
      + "data, without data, model data). The tabs follow the old program's order; View › Every field and column "
      + "(the expert's depth) shows every field."));
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
  return h("div.config-text", {}, h("p.hint", {}, "Every <config> entry as text: an empty value removes the key, a "
    + "new key and value go in the last row; Enter applies."), table,
    h("div.form-actions", {}, button("Apply", "set_config with the changed entries (Enter)", apply, "primary")));
}

/** The session's options as GET /api/describe reports them: keyword, value in use, its source, default, flag, doc. */
function optionsPart(ctx) {
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
function problemsPart(ctx, problems, records) {
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

function historyPart(ctx, history) {
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
function textPart(ctx, msg) {
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
function enginePart(ctx) {
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
