// The Data tab (docs/upgrade/web_ui.md section 6.3): the records of every section of the schema (except those the
// Nucleus tab claims), grouped by section; the χ² per data set from the evaluation; one record's generated form with
// the engine's description; the first rows of its data file; a preview plot of the data and the model.

import { pathOf } from "../api.js";
import { plotOf } from "../plots.js";
import { renderTable } from "../tables.js";
import { button, drawPlot, h } from "../ui.js";
import { FIRST_VIEW, SECTION_TITLES } from "../layout.js";
import { listColumns, listKinds, onlyColumns, usersOf } from "../lists.js";
import { NUCLEUS_SECTIONS, newRecordPanel, recordPanel } from "./records.js";

/** The sections the Data tab shows, in the schema's order, with their record kinds. */
export function dataSections(schema) {
  return Object.entries(schema.sections).filter(([s]) => !NUCLEUS_SECTIONS.includes(s));
}

/** The groups a view of the first view's Data tabs shows (layout.js FIRST_VIEW): its `sections`, or, for the view
 * without them, every group no other view names (a section the schema gains appears there). */
export function sectionsOfView(sections, view, firstView = FIRST_VIEW) {
  if (!view) return sections;
  if (view.sections) return sections.filter(([s]) => view.sections.includes(s));
  const named = firstView.filter((v) => v.tab === "data" && v.sections).flatMap((v) => v.sections);
  return sections.filter(([s]) => !named.includes(s));
}

/** What the engine makes of the records it describes, by name (GET /api/project/describe; none while the project
 * has an error: the lists then show the records' fields). */
async function describedRecords(ctx) {
  try {
    const d = await ctx.api.get("/api/project/describe", { quiet: true });
    return Object.fromEntries(Object.values(d).filter(Array.isArray).flat().filter((e) => e && e.name)
      .map((e) => [e.name, e]));
  } catch (e) {
    return {};
  }
}

export async function render(ctx, el, args) {
  const [records, problems] = await Promise.all([ctx.api.get("/api/project/records"),
    ctx.api.get("/api/project/validate")]);
  if (!ctx.kinds) ctx.kinds = await ctx.api.get("/api/kinds");
  const every = dataSections(ctx.schema);
  const selected = args[0] && args[0] !== "+new" ? records.find((r) => r.name === args[0]) : null;
  const sectionOf = (kind) => every.find(([, kinds]) => kinds.includes(kind))?.[0];
  const current = ctx.route.params.section || (selected ? sectionOf(selected.kind) : null) || every[0][0];
  ctx.section = current;  // the first view marks its tab by the group drawn
  const kinds = ctx.schema.sections[current] || [];
  // the groups of the tab the group drawn belongs to (Segments or Experimental Effects), in the old program's words
  const sections = sectionsOfView(every, firstViewOf(current));
  const title = (s) => SECTION_TITLES[s] || s;

  // the group bar: one button per section with its number of records
  const bar = h("div.group-bar", { role: "toolbar", "aria-label": "sections" });
  for (const [s, ks] of sections) {
    const n = records.filter((r) => ks.includes(r.kind)).length;
    const a = h("a.group", { href: `#/data?section=${s}`, title: `${s}: ${ks.join(", ")}` }, `${title(s)} `,
      h("span.count", {}, n));
    if (s === current) a.classList.add("active");
    bar.appendChild(a);
  }
  el.appendChild(bar);

  // one list per record kind of the group, in words (lists.js): the data sets first, then the records they name
  const evaluated = Object.fromEntries((ctx.evaluation?.datasets || []).map((d) => [d.name, d]));
  const problemsOf = (name) => problems.filter((d) => d.record === name);
  const [described, references] = await Promise.all([describedRecords(ctx), ctx.api.get(
    "/api/project/references", { quiet: true }).catch(() => [])]);
  const users = usersOf(references);
  // Add at the table (the old program's +); the selected record's Duplicate, Rename and Delete head its panel
  const add = button(`Add to ${title(current)}`, `Add a record of ${kinds.join(", ")}`,
    () => ctx.navigate(`#/data/+new?section=${current}`), "primary");
  add.dataset.action = "tool-add-record";
  const list = h(`section.panel.list${selected || args[0] === "+new" ? ".with-detail" : ""}`, {},
    h("div.table-tools", {}, add, h("span.hint", {}, selected ? `${selected.name} below: Duplicate, Rename and `
      + "Delete at its head" : "select a row: its form opens below the list")));
  const mine = records.filter((r) => kinds.includes(r.kind));
  for (const kind of listKinds(ctx.schema, kinds)) {
    const rows = mine.filter((r) => r.kind === kind).map((r) => ({ ...r, problems: problemsOf(r.name).length || "" }));
    if (!rows.length) continue;
    const all = [...listColumns({ schema: ctx.schema, kind, rows, described,
      descriptionFields: ctx.kinds?.description?.fields || {}, users, evaluated }),
    { key: "problems", title: "⚠", numeric: true, doc: "messages of the engine's validator" },
    { key: "line", title: "line", numeric: true, doc: "the record's first line in the project file" }];
    // the first view: the old program's columns of a list the engine describes (the server's table)
    const describedKind = rows.some((r) => described[r.name]);
    const columns = !ctx.expert && describedKind ? onlyColumns(all, ctx.describe?.first_view?.description) : all;
    const rec = ctx.schema.records[kind] || {};
    const table = renderTable(document, { columns, rows, selected: selected?.name, label: `${rec.title || kind}`,
      withFilter: rows.length > 8, onSelect: (r) => ctx.navigate(`#/data/${encodeURIComponent(r.name)}`) });
    list.appendChild(h("div.kind-list", { data: { kind } }, h("h3", { title: `record kind ${kind}` },
      `${rec.title || kind}`, h("span.count", {}, ` ${rows.length}`)), h("p.hint.kind-doc", {}, rec.doc || ""),
    table.element));
  }
  if (!mine.length) {
    // an empty group says what belongs there and offers the first one
    const said = kinds.map((k) => `${ctx.schema.records[k]?.title || k}: `
      + `${firstSentence(ctx.schema.records[k]?.doc)}`);
    list.appendChild(h("p.hint.empty", {}, `No ${title(current)} yet. ${said.join(" ")} `, button(`+ new`,
      `Add the first record of ${kinds.join(", ")}`,
      () => ctx.navigate(`#/data/+new?section=${current}`), "primary")));
  }
  if (ctx.evalError) list.appendChild(h("p.hint", {}, `No evaluation: ${ctx.evalError.what}`));
  if (!(ctx.route.params.preview && selected)) el.appendChild(list);  // a linked preview shows the record alone

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

/** The first sentence of the engine's description of a kind. */
const firstSentence = (doc) => (String(doc || "").match(/^.*?[.;](\s|$)/)?.[0] || String(doc || "")).trim();

/** The first view's tab a group of the Data tab belongs to (a view naming it, else the view without sections). */
function firstViewOf(section) {
  const views = FIRST_VIEW.filter((v) => v.tab === "data");
  return views.find((v) => v.sections?.includes(section)) || views.find((v) => !v.sections) || null;
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
