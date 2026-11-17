// The THM tab (docs/upgrade/web_ui.md section 6.7), shown when the project has THM content: the model options of
// <thm> and its experiments as forms of the form generator (the keys are the engine's table, pyazr.kinds()["thm"],
// GET /api/thm/keys; a nucleus or global optical potential offered from pyazr's lists), the THM pairs and channels
// of the session (edited in the Nucleus tab), and the session's THM diagnostics, drawn as the Physics tab draws its
// functions (functionview.js).
//
// Every edit is a project edit through the THM routes (PATCH /api/project/thm, PUT and DELETE
// /api/project/thm/experiments/{name}): one history entry each, refused with the engine's words.  The tab holds no
// model state: only the diagnostics' form texts and view choices while the page lives (#/thm/<experiment>?fn=<row>).

import { pathOf } from "../api.js";
import { describeForm, tokenize } from "../forms.js";
import { renderForm } from "../formview.js";
import { computeRow } from "../functionview.js";
import { applyParams, newState, rowForm } from "../physicsform.js";
import { renderTable } from "../tables.js";
import { button, h, section } from "../ui.js";
import { loadRows } from "./physics.js";

/** The tab's id: the session's rows that name it (`tab`) are its diagnostics. */
export const ID = "thm";
const SOURCES = [{ rows: "/api/session/functions", call: "/api/session/functions", session: true }];
/** The key of the form's field of data sets that join an experiment (not a key of <thm>). */
export const JOIN = "@join";
const NEW = "+new";

/** Page state: the key schema, pyazr's lists and the rows read once per page, each diagnostic's state. */
const page = { keys: null, choices: null, rows: null, states: {}, project: null };

/** The suggestions of the key forms' free-text fields, by the grammar token that names them (forms.js
 * suggestionsOf): the nuclides and the global optical potentials the engine knows (GET /api/kinds). */
export function choicesOf(kinds) {
  return { "<nucleus>": Object.keys(kinds?.nuclides || {}),
    "<global set>": Object.keys(kinds?.optical_potentials || {}) };
}

/** Whether the project has THM content: a <thm> section, or a data set that names an experiment (its `thm`
 * field). */
export function hasThm(thm, records) {
  return thm !== null && thm !== undefined ? true : (records || []).some((r) => r.fields && "thm" in r.fields);
}

/** Whether the tab is listed (app.js TABS): the project's THM content. */
export async function shown(ctx) {
  const [thm, records] = await Promise.all([ctx.api.get("/api/project/thm"), ctx.api.get("/api/project/records")]);
  return hasThm(thm, records);
}

/**
 * The schema the forms are made from: the key table's record kind with, for every key the project has and the table
 * does not list by its name (spectatorEnergy[<pair key>]), a text field of its own; `extra`: more fields.
 */
export function keySchema(keys, kind, written, extra = {}) {
  const rec = keys[kind];
  const fields = { ...rec.fields };
  for (const k of Object.keys(written || {})) {
    if (!fields[k]) fields[k] = { type: "text", title: k, doc: "a key the key table does not list by this name (the "
      + "engine reads it)" };
  }
  Object.assign(fields, extra);
  return { records: { [kind]: { ...rec, fields } } };
}

/** The data sets of the project that may join an experiment (those with a `thm` field) and the members of each. */
export function membersOf(records) {
  const thm = (records || []).filter((r) => r.fields && "thm" in r.fields);
  const members = {};
  for (const r of thm) (members[String(r.fields.thm)] ||= []).push(r.name);
  return { names: thm.map((r) => r.name), members };
}

/** The request of an experiment form's changed fields: the keys to set (null removes one) and the data sets that
 * join. */
export function experimentEdit(changed) {
  const keys = { ...changed };
  const join = tokenize(keys[JOIN] ?? "");
  delete keys[JOIN];
  return { keys, datasets: join };
}

export async function render(ctx, el, args) {
  if (page.project !== ctx.project.path) Object.assign(page, { project: ctx.project.path, states: {}, rows: null });
  const [thm, records] = await Promise.all([ctx.api.get("/api/project/thm"), ctx.api.get("/api/project/records")]);
  if (!page.keys) page.keys = await ctx.api.get("/api/thm/keys");
  if (!page.choices) page.choices = choicesOf(await ctx.api.get("/api/kinds"));
  if (!page.rows) page.rows = await loadRows(ctx, SOURCES, ID);
  const msg = h("div.messages");
  el.appendChild(msg);
  const section_ = thm || { options: {}, experiments: {} };
  if (!hasThm(thm, records)) {
    el.appendChild(h("p.hint.no-thm", {}, "This project has no THM content (no <thm> section and no data set that "
      + "names an experiment), so the tab is not listed. Setting an option below adds the <thm> section."));
  }
  const groups = membersOf(records);
  el.appendChild(h("div.two-columns", {}, optionsSection(ctx, section_.options, msg),
    experimentsSection(ctx, section_.experiments, groups, args[0], msg)));
  const diagnostics = h("div");
  el.append(diagnostics, await channelsSection(ctx));
  diagnostics.appendChild(await diagnosticsSection(ctx, Object.keys(section_.experiments), args[0]));
}

/** The global options: one generated form; Apply sends the changed keys (PATCH /api/project/thm); "other key" adds
 * one the table does not list (spectatorEnergy[<pair key>]). */
function optionsSection(ctx, options, msg) {
  const schema = keySchema(page.keys, "<thm>", options);
  const form = describeForm(schema, "<thm>", options, { schema, choices: page.choices });
  const element = renderForm(document, form, { schema, choices: page.choices }, {
    onApply: (changed) => ctx.edit("set <thm> options", () => ctx.api.patch("/api/project/thm",
      { options: changed }), msg) });
  element.onrevert = () => ctx.refresh();
  const key = h("input", { type: "text", placeholder: "spectatorEnergy[2]", size: 22, "aria-label": "other key" });
  const value = h("input", { type: "text", placeholder: "value", size: 12, "aria-label": "its value" });
  const other = h("div.form-actions", {}, h("span.k", {}, "other key "), key, value,
    button("Set", "Set a key the form does not list (PATCH /api/project/thm)", () => {
      if (!key.value.trim()) return;
      ctx.edit(`set <thm> ${key.value.trim()}`, () => ctx.api.patch("/api/project/thm",
        { options: { [key.value.trim()]: value.value.trim() === "" ? null : value.value.trim() } }), msg)
        .catch(() => {});
    }));
  return section("Model options", h("p.hint", {}, page.keys["<thm>"].doc + " An empty value removes the key "
    + "(the engine's default applies)."), element, other);
}

/** The experiments: a table (name, data sets, keys), the form of the selected one with Remove, and a new one. */
function experimentsSection(ctx, experiments, groups, selected, msg) {
  const rows = Object.entries(experiments).map(([name, keys]) => ({ name,
    datasets: (groups.members[name] || []).join(" "),
    keys: Object.entries(keys).map(([k, v]) => `${k}=${v}`).join(" ") }));
  const table = renderTable(document, { rows, selected: selected || null, withFilter: false, label: "experiments",
    onSelect: (r) => ctx.navigate(`#/thm/${encodeURIComponent(r.name)}`),
    columns: [{ key: "name", title: "experiment" }, { key: "datasets", title: "data sets", wrap: true },
      { key: "keys", title: "keys as written", wrap: true }] });
  const s = section(`Experiments (${rows.length})`, h("p.hint", {}, page.keys.experiment.doc + " A data set joins "
    + "an experiment by its own thm field (the Data tab), or here when the experiment is added or edited."),
  table.element, h("div.form-actions", {}, button("+ experiment", "Add an experiment",
    () => ctx.navigate(`#/thm/${NEW}`))));
  const join = { [JOIN]: { type: "text", title: "data sets that join", doc: "THM data sets (names separated by "
    + "blanks) whose thm field is set to this experiment in the same edit", required: false } };
  const fctx = (schema) => ({ schema, pairLabels: [], spins: [], choices: page.choices });
  if (selected === NEW) {
    const schema = keySchema(page.keys, "experiment", {}, join);
    const form = renderForm(document, describeForm(schema, "experiment", {}, fctx(schema)), fctx(schema), {
      nameInput: true, applyLabel: "Add", onApply: (name, changed) => {
        const { keys, datasets } = experimentEdit(changed);
        return ctx.edit(`add experiment[${name}]`, () => ctx.api.put(pathOf("/api/project/thm/experiments", name),
          { keys, datasets }), msg).then(() => ctx.navigate(`#/thm/${encodeURIComponent(name)}`));
      } });
    form.onrevert = () => ctx.navigate("#/thm");
    s.append(h("h3", {}, "New experiment"), hintOfNames(groups), form);
  } else if (selected && experiments[selected]) {
    const written = experiments[selected];
    const schema = keySchema(page.keys, "experiment", written, join);
    const form = renderForm(document, describeForm(schema, "experiment", written, fctx(schema)), fctx(schema), {
      onApply: (changed) => {
        const { keys, datasets } = experimentEdit(changed);
        return ctx.edit(`set experiment[${selected}]`, () => ctx.api.put(pathOf("/api/project/thm/experiments",
          selected), { keys, datasets }), msg);
      },
      actions: [{ label: "Remove", title: "Remove the experiment (refused while a data set names it)",
        run: () => ctx.edit(`remove experiment[${selected}]`, () => ctx.api.del(pathOf(
          "/api/project/thm/experiments", selected), {}), msg).then(() => ctx.navigate("#/thm"), () => {}) }] });
    form.onrevert = () => ctx.refresh();
    s.append(h("h3", {}, `experiment[${selected}]`), h("p.hint", {}, `Its data sets: `
      + `${(groups.members[selected] || []).join(", ") || "none"}.`), hintOfNames(groups), form);
  } else if (selected) {
    s.appendChild(h("p.hint", {}, `No experiment ${selected} in the project.`));
  }
  return s;
}

function hintOfNames(groups) {
  return h("p.hint", {}, `THM data sets of the project: ${groups.names.join(", ") || "none"}.`);
}

/** The THM pairs of the session (binding energy B) and their channels (whether the width is an amplitude), as the
 * session reads them; both are edited in the Nucleus tab (a pair's B, a channel's amplitude flag). */
async function channelsSection(ctx) {
  const why = "As the session reads them. A pair's binding energy B is edited with the pair, a channel's amplitude "
    + "flag with its level, in the Nucleus tab (set_pair, set_channel): one edit each, Undo takes it back.";
  let session;
  try {
    session = await ctx.api.get("/api/session", { quiet: true });
  } catch (e) {
    const s = section("Channels", h("p.hint", {}, why));
    ctx.showError(e, s);
    return s;
  }
  const pairs = (session.pairs || []).filter((p) => p.binding_energy);
  const keys = new Set(pairs.map((p) => p.key));
  const channels = (session.parameters || []).filter((p) => "width_is_amplitude" in p && keys.has(p.pair_key))
    .map((p) => ({ name: p.name, level: p.level_label, pair: p.pair_key, L: p.L, S: p.S,
      amplitude: p.width_is_amplitude }));
  const pairTable = renderTable(document, { rows: pairs.map((p) => ({ ...p, name: String(p.key) })),
    onSelect: (r) => ctx.navigate(`#/nucleus/pair/${r.key}`),
    withFilter: false, label: "THM pairs", caption: "THM pairs (a nonzero binding energy)",
    columns: [{ key: "key", title: "pair", numeric: true }, { key: "label", title: "particles (J^π)" },
      { key: "binding_energy", title: "B", unit: "MeV", numeric: true, doc: "the THM binding energy of the pair's "
        + "composite projectile" }, { key: "radius", title: "radius", unit: "fm", numeric: true,
        value: (r) => r.channel_radius }] });
  const channelTable = renderTable(document, { rows: channels, label: "THM channels",
    onSelect: (r) => ctx.navigate(`#/nucleus/level/${encodeURIComponent(r.level)}`),
    caption: "their channels", columns: [{ key: "level", title: "level" }, { key: "pair", title: "pair",
      numeric: true }, { key: "L", title: "l", numeric: true }, { key: "S", title: "s", numeric: true },
    { key: "amplitude", title: "width is an amplitude" }, { key: "name", title: "parameter" }] });
  return section("Channels", h("p.hint", {}, why), pairs.length ? pairTable.element
    : h("p.hint", {}, "The session has no pair with a binding energy."), channels.length ? channelTable.element : null);
}

/** The session's THM diagnostics: one row chosen (?fn=), its generated form, its answer drawn by functionview.js. */
async function diagnosticsSection(ctx, experiments, selected) {
  const rows = page.rows;
  const s = section("Diagnostics");
  if (!rows.length) { s.appendChild(h("p.hint", {}, "The server lists no THM diagnostic.")); return s; }
  const row = rows.find((r) => r.name === ctx.route.params.fn) || rows[0];
  const fctx = { experiments };
  if (!page.states[row.name]) {
    const state = newState(row, fctx);
    if (selected && experiments.includes(selected) && "experiment" in state.texts) state.texts.experiment = selected;
    page.states[row.name] = state;
  }
  const state = page.states[row.name];
  applyParams(row, state, ctx.route.params, fctx);
  const at = (r) => `#/thm${selected ? `/${encodeURIComponent(selected)}` : ""}?fn=${encodeURIComponent(r.name)}`;
  s.appendChild(h("nav.fn-choices", { "aria-label": "diagnostics" }, ...rows.map((r) => h(
    `a.fn${r === row ? ".active" : ""}`, { href: at(r), title: r.summary, data: { fn: r.name } }, r.name))));
  s.appendChild(h("div.fn-head", {}, h("h3.mono", {}, `${row.name}(${row.arguments.map((a) => a.name).join(", ")})`),
    h("div", {}, row.summary), h("div.hint", { title: row.doc }, `${row.doc} Library: ${row.library}.`)));
  if (row.arguments.some((a) => a.kind === "experiment") && !experiments.length) {
    s.appendChild(h("p.hint", {}, "The project has no THM experiment: this diagnostic needs one."));
  }
  const out = h("div.physics-out");
  const form = renderForm(document, rowForm(row, state.axis, state.texts, fctx), fctx, {
    whole: true, applyLabel: "Compute", onApply: (texts) => {
      state.texts = texts;
      return computeRow(ctx, out, row, state, { where: "the session" });
    } });
  form.onrevert = () => { delete page.states[row.name]; ctx.navigate(at(row)); };
  s.appendChild(h("div.physics-main", {}, form, out));  // laid out as the Physics tab's form and answer
  await computeRow(ctx, out, row, state, { where: "the session" });
  return s;
}
