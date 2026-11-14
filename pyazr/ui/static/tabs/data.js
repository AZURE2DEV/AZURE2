// The Data tab (docs/upgrade/web_ui.md section 6.3): the records of every section of the schema (except those the
// Nucleus tab claims), grouped by section; the χ² per data set from the evaluation; one record's generated form with
// the engine's description; the first rows of its data file; a preview plot of the data and the model.

import { pathOf } from "../api.js";
import { plotOf } from "../plots.js";
import { renderTable } from "../tables.js";
import { button, drawPlot, h } from "../ui.js";
import { NUCLEUS_SECTIONS, newRecordPanel, recordPanel } from "./records.js";

/** The sections the Data tab shows, in the schema's order, with their record kinds. */
export function dataSections(schema) {
  return Object.entries(schema.sections).filter(([s]) => !NUCLEUS_SECTIONS.includes(s));
}

/** Up to three fields a record kind is summarised by in the list: its first required fields of a one-word type. */
export function summaryFields(schema, kinds) {
  const plain = new Set(["enum", "unit", "reaction_expr", "reference", "nucleus", "integer", "spin_parity", "pair"]);
  const keys = [];
  for (const k of kinds) {
    for (const [key, f] of Object.entries(schema.records[k]?.fields || {})) {
      if (f.required && plain.has(f.type) && !keys.includes(key) && keys.length < 3) keys.push(key);
    }
  }
  return keys;
}

export async function render(ctx, el, args) {
  const [records, problems] = await Promise.all([ctx.api.get("/api/project/records"),
    ctx.api.get("/api/project/validate")]);
  const sections = dataSections(ctx.schema);
  const selected = args[0] && args[0] !== "+new" ? records.find((r) => r.name === args[0]) : null;
  const sectionOf = (kind) => sections.find(([, kinds]) => kinds.includes(kind))?.[0];
  const current = ctx.route.params.section || (selected ? sectionOf(selected.kind) : null) || sections[0][0];
  const kinds = ctx.schema.sections[current] || [];

  // the group bar: one button per section with its number of records
  const bar = h("div.group-bar", { role: "toolbar", "aria-label": "sections" });
  for (const [s, ks] of sections) {
    const n = records.filter((r) => ks.includes(r.kind)).length;
    const a = h("a.group", { href: `#/data?section=${s}`, title: ks.join(", ") }, `${s} `, h("span.count", {}, n));
    if (s === current) a.classList.add("active");
    bar.appendChild(a);
  }
  bar.append(h("span.spacer"), button("+ new", `Add a record of ${kinds.join(", ")}`,
    () => ctx.navigate(`#/data/+new?section=${current}`)));
  el.appendChild(bar);

  // the list of the section's records, with χ² and points per data set from the evaluation
  const evaluated = Object.fromEntries((ctx.evaluation?.datasets || []).map((d) => [d.name, d]));
  const problemsOf = (name) => problems.filter((d) => d.record === name);
  const rows = records.filter((r) => kinds.includes(r.kind)).map((r) => {
    const row = { name: r.name, kind: r.kind, line: r.line, problems: problemsOf(r.name).length || "" };
    for (const k of summaryFields(ctx.schema, kinds)) row[k] = Array.isArray(r.fields[k]) ? r.fields[k].join("; ")
      : r.fields[k] ?? "";
    const e = evaluated[r.name];
    if (e) Object.assign(row, { chi2: e.chi2, n: e.n, norm: e.norm });
    return row;
  });
  const columns = [{ key: "name", title: "name" }];
  if (kinds.length > 1) columns.push({ key: "kind", title: "kind" });
  for (const k of summaryFields(ctx.schema, kinds)) columns.push({ key: k, title: k });
  if (rows.some((r) => r.chi2 !== undefined)) {
    columns.push({ key: "chi2", title: "χ²", numeric: true, digits: 6, doc: "the program's χ² of the data set" },
      { key: "n", title: "N", numeric: true, doc: "points used" },
      { key: "norm", title: "norm", numeric: true, digits: 6, doc: "the normalization in use" });
  }
  columns.push({ key: "problems", title: "⚠", numeric: true, doc: "messages of the engine's validator" },
    { key: "line", title: "line", numeric: true });
  const table = renderTable(document, { columns, rows, selected: selected?.name, label: `${current} records`,
    caption: `${current}: ${kinds.map((k) => ctx.schema.records[k]?.title || k).join(", ")}`,
    onSelect: (r) => ctx.navigate(`#/data/${encodeURIComponent(r.name)}`) });
  const list = h("section.panel.list", {}, table.element);
  if (ctx.evalError) list.appendChild(h("p.hint", {}, `No evaluation: ${ctx.evalError.what}`));
  el.appendChild(list);

  if (args[0] === "+new") {
    el.appendChild(await newRecordPanel(ctx, kinds, records,
      { onAdded: (r) => ctx.navigate(`#/data/${encodeURIComponent(r.name)}`) }));
    return;
  }
  if (!selected) {
    if (args[0]) el.appendChild(h("p.hint", {}, `No record named ${args[0]}.`));
    else el.appendChild(h("p.hint", {}, "Choose a record (click, or ↓ and Enter in the list) to see its form."));
    return;
  }
  el.appendChild(await recordPanel(ctx, selected, records, problems,
    { onGone: () => ctx.navigate(`#/data?section=${current}`) }));
  const extras = h("div.record-extras");
  el.appendChild(extras);
  const file = selected.fields.file;
  if (typeof file === "string" && file) extras.appendChild(fileRows(ctx, file));
  if (ctx.schema.records[selected.kind]?.fields?.covariance) extras.appendChild(button("Covariance…",
    "Build the covariance from parts", () => ctx.navigate(`#/covariance/${encodeURIComponent(selected.name)}`)));
  if (evaluated[selected.name]) extras.appendChild(await preview(ctx, selected.name));
}

/** The first rows of a data file, as the server reads them (the text, not a parse). */
function fileRows(ctx, file) {
  const box = h("section.panel.file-rows", {}, h("h3", {}, `First rows of ${file}`));
  const pre = h("pre.file-head", {}, "…");
  box.appendChild(pre);
  ctx.api.get(`${pathOf("/api/files", ...file.replace(/^"|"$/g, "").split("/"))}?rows=8`, { quiet: true })
    .then((r) => { pre.textContent = r.lines.join("\n") + (r.more ? "\n…" : ""); })
    .catch((e) => { pre.textContent = e.view ? e.view.what : String(e); });
  return box;
}

/** The preview plot of one data set at the working x. */
async function preview(ctx, name) {
  const box = h("section.panel.preview", {}, h("h3", {}, `Preview: ${name}, data and model at the working x`));
  const area = h("div.plot");
  box.appendChild(area);
  try {
    const data = await ctx.api.get(pathOf("/api/evaluation/data", name));
    const plot = plotOf(data, { schema: ctx.schema });
    if (plot) area.dataset.drawn = drawPlot(area, plot, name, { usePlotly: ctx.route.params.plot !== "table" });
    else area.textContent = "No plot row applies to this data set.";
  } catch (e) {
    ctx.showError(e, box);
  }
  return box;
}
