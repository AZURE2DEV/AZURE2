// The Physics tab (docs/upgrade/web_ui.md section 6.6): the engine's functions on a range, one generated form and
// one plot per row of the function tables (the physics table, the session's functions), with the status of every
// element and the pyazr call that makes the same arrays.  The rules (form, request, plots, status line, CSV) are
// physicsform.js; this module draws them.  Every number drawn comes from a route; the tab holds no model state,
// only the form's texts and view choices per function while the page lives (hash parameters set them:
// #/physics/<fn>?v.l=0 1&logy=1).  Computing a row and drawing its answer is functionview.js, shared with the THM tab.

import { renderForm } from "../formview.js";
import { computeRow } from "../functionview.js";
import {
  axisChoices, groupedRows, newState, rowForm, startTexts, applyParams, MAX_ELEMENTS, SPACING,
} from "../physicsform.js";
import { h, section } from "../ui.js";

/** The function tables the tab lists: where the rows are read and where a row's request goes.  A new table is a
 * row here (a server route answering rows of the physics table's shape). */
export const SOURCES = [
  { rows: "/api/physics", call: "/api/physics", doc: "the engine's physics table (pyazr.physics)" },
  { rows: "/api/stopping", call: "/api/stopping", doc: "the stopping library (SRIM-2013; cite SRIM)" },
  { rows: "/api/session/functions", call: "/api/session/functions", session: true,
    doc: "the session's functions (need the project to build)" },
];

/** Page state: the rows read once per page, the state of each function, the search text. */
const page = { rows: null, states: {}, search: "", project: null };

/** The rows of every table of SOURCES, each with its source; a row that names a tab (the THM diagnostics) is drawn by
 * that tab, not here. */
export async function loadRows(ctx, sources = SOURCES, tab = null) {
  const lists = await Promise.all(sources.map((s) => ctx.api.get(s.rows).then(
    (rows) => Object.values(rows).filter((r) => (r.tab || null) === tab).map((r) => ({ ...r, source: s })),
    (e) => { ctx.status(`${s.rows}: ${e.view?.what || e}`, "warn"); return []; })));
  return lists.flat();
}

export async function render(ctx, el, args) {
  if (!page.rows) page.rows = await loadRows(ctx);
  if (page.project !== ctx.project.path) { page.project = ctx.project.path; page.states = {}; }
  let pairs = [];
  try { pairs = await ctx.api.get("/api/project/pairs", { quiet: true }); } catch { /* no pairs: shown in the form */ }
  const fctx = { pairs: Array.isArray(pairs) ? pairs : pairs.result || [] };
  const row = page.rows.find((r) => r.name === args[0]) || page.rows[0];
  if (!row) { el.appendChild(section("Physics", h("p.hint", {}, "The server lists no function."))); return; }
  if (!page.states[row.name]) page.states[row.name] = newState(row, fctx);
  const state = page.states[row.name];
  applyParams(row, state, ctx.route.params, fctx);
  const panel = h("div.physics-main");
  el.appendChild(h("div.physics", {}, listPanel(ctx, row), panel));
  await drawFunction(ctx, panel, row, state, fctx);
}

/** The list of functions, grouped as the tables group them, with a search box. */
function listPanel(ctx, current) {
  const list = h("nav.physics-list", { "aria-label": "functions" });
  const box = h("input.table-filter", { type: "search", placeholder: "search (/)", "aria-label": "search functions",
    value: page.search });
  const fill = () => {
    list.querySelectorAll(".fn-group").forEach((g) => g.remove());
    for (const g of groupedRows(page.rows, page.search)) {
      list.appendChild(h("div.fn-group", {}, h("div.group-title", {}, g.group),
        ...g.rows.map((r) => h(`a.fn${r === current ? ".active" : ""}`, { href: `#/physics/${encodeURIComponent(
          r.name)}`, title: r.summary, data: { fn: r.name } }, r.name))));
    }
  };
  box.addEventListener("input", () => { page.search = box.value; fill(); });
  list.appendChild(box);
  fill();
  return list;
}

/** One function: its title, the axis choice, the generated form, then the plot of its answer (the promise of its
 * computation). */
function drawFunction(ctx, el, row, state, fctx) {
  el.textContent = "";
  const signature = `${row.name}(${row.arguments.map((a) => a.name).join(", ")})`;
  el.appendChild(h("div.fn-head", {}, h("h2.mono", {}, signature), h("div", {}, row.summary),
    h("div.hint", { title: row.doc }, `${row.doc} Library: ${row.library}.`)));
  const axes = axisChoices(row);
  const axisBox = h("select", { "aria-label": "plot against" },
    ...axes.map((a) => h("option", { value: a.name }, a.name)));
  axisBox.value = state.axis;
  axisBox.addEventListener("change", () => {
    state.axis = axisBox.value;
    state.texts = startTexts(row, state.axis, fctx);
    state.logx = axes.find((a) => a.name === state.axis)?.sweep?.spacing === "log";
    drawFunction(ctx, el, row, state, fctx);
  });
  el.appendChild(h("label.choice", {}, "plot against ", axisBox));
  if (row.source.session && !fctx.pairs.length) {
    el.appendChild(h("p.hint", {}, "The project has no particle pairs: this function needs one."));
  }
  const out = h("div.physics-out");
  const form = renderForm(document, rowForm(row, state.axis, state.texts, fctx), fctx, {
    whole: true, applyLabel: "Plot", onApply: (texts) => {
      state.texts = texts;
      state.logx = texts[SPACING] === "log";  // the axis as spaced; the box changes it after
      return compute(ctx, out, row, state);
    },
  });
  form.onrevert = () => {
    state.texts = startTexts(row, state.axis, fctx);
    drawFunction(ctx, el, row, state, fctx);
  };
  el.append(form, h("div.hint", {}, `At most ${MAX_ELEMENTS} elements per request (points × lines); several values `
    + "of an argument (separated by blanks) give one line each."), out);
  return compute(ctx, out, row, state);
}

/** Send the request of the form's texts, then draw its answer (or the server's refusal). */
function compute(ctx, out, row, state) {
  return computeRow(ctx, out, row, state, row.source.session
    ? { where: "the session (the free functions above need no project)",
      refusalNote: "The functions of the engine's tables above need no project and keep working." }
    : {});
}
