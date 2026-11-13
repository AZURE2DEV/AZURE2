// The Physics tab (docs/upgrade/web_ui.md section 6.6): the engine's functions on a range, one generated form and
// one plot per row of the function tables (the physics table, the session's functions), with the status of every
// element and the pyazr call that makes the same arrays.  The rules (form, request, plots, status line, CSV) are
// physicsform.js; this module draws them.  Every number drawn comes from a route; the tab holds no model state,
// only the form's texts and view choices per function while the page lives (hash parameters set them:
// #/physics/<fn>?v.l=0 1&logy=1).

import { pathOf } from "../api.js";
import { renderForm } from "../formview.js";
import {
  answerCsv, axisChoices, groupedRows, newState, physicsPlots, requestOf, rowForm, statusLine, startTexts,
  applyParams, MAX_ELEMENTS, SPACING, unitGroups,
} from "../physicsform.js";
import { button, download, drawPlot, h, section } from "../ui.js";

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

async function loadRows(ctx) {
  const lists = await Promise.all(SOURCES.map((s) => ctx.api.get(s.rows).then(
    (rows) => Object.values(rows).map((r) => ({ ...r, source: s })),
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
async function compute(ctx, out, row, state) {
  const request = requestOf(row, state.axis, state.texts);
  out.textContent = "";
  if (request.error) { out.appendChild(h("div.error-box", { role: "alert" }, request.error)); return; }
  const where = row.source.session ? "the session (the free functions above need no project)" : "the server";
  out.appendChild(h("p.hint", {}, `computing in ${where}…`));
  let answer;
  try {
    answer = await ctx.api.request("POST", pathOf(row.source.call, row.name), { arguments: request.arguments },
      { quiet: true });
  } catch (e) {
    out.textContent = "";
    ctx.showError(e, out);
    if (row.source.session) out.appendChild(h("p.hint", {}, "The functions of the engine's tables above need no "
      + "project and keep working."));
    return;
  }
  state.answer = { answer, request };
  drawAnswer(ctx, out, row, state);
}

/** The status line, the field and log choices, the plots, the exports and the Python text of an answer. */
function drawAnswer(ctx, out, row, state) {
  const { answer, request } = state.answer;
  out.textContent = "";
  const st = statusLine(answer);
  out.appendChild(h("div.physics-status", { role: "status", data: { computed: st.computed, elements: st.elements } },
    st.text));
  const redraw = () => drawAnswer(ctx, out, row, state);
  const flag = (label, key) => {
    const input = h("input", { type: "checkbox", "aria-label": label });
    input.checked = !!state[key];
    input.addEventListener("change", () => { state[key] = input.checked; redraw(); });
    return h("label.choice", {}, input, ` ${label}`);
  };
  const fields = row.returns.map((r) => {
    const input = h("input", { type: "checkbox", "aria-label": `show ${r.name}` });
    input.checked = state.shown.has(r.name);
    input.addEventListener("change", () => {
      if (input.checked) state.shown.add(r.name); else state.shown.delete(r.name);
      redraw();
    });
    return h("label.choice", { title: `${r.doc}${r.unit ? ` (${r.unit})` : ""}` }, input, ` ${r.name}`);
  });
  out.appendChild(h("div.choices", {}, h("span.k", {}, "fields "), ...fields, flag("log x", "logx"),
    flag("log y", "logy")));
  // the citation a row owes (SRIM's for the stopping rows) goes with every plot, the CSV and the Python text
  const citation = answer.citation || row.citation || "";
  const csv = () => answerCsv(row, answer, citation ? [`cite: ${citation}`] : []);
  const file = `${row.name}`;
  for (const plot of physicsPlots(row, answer, request, state)) {
    const area = h("div.plot.physics-plot");
    out.appendChild(area);
    area.dataset.drawn = drawPlot(area, plot, file, { usePlotly: !state.table, csv });
    if (citation) area.appendChild(h("p.citation", { data: { citation } }, `Cite: ${citation}`));
  }
  if (!unitGroups(row, state.shown).length) out.appendChild(h("p.hint", {}, "No field is shown."));
  const code = h("pre.python.mono", { "aria-label": "the pyazr call" }, answer.python);
  const copy = button("Copy as Python", "The pyazr call that makes the same arrays, to the clipboard", async () => {
    try {
      await navigator.clipboard.writeText(answer.python);
      ctx.status("copied the pyazr call", "ok");
    } catch {
      getSelection().selectAllChildren(code);
      ctx.status("the clipboard is not available: the call is selected, copy it with Ctrl+C", "warn");
    }
  });
  out.append(h("div.plot-bar", {}, copy, button("Python file", "The pyazr call as a .py file",
    () => download(`${file}.py`, answer.python, "text/x-python"))), code);
}
