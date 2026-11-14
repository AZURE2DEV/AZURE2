// A function of a table (the physics table, the stopping table, the session's functions) computed and drawn: the
// request of a form's texts sent to the row's route, then the answer: its status line, its facts (the call's scalars
// with units), its tables (an answer that is not a function of an axis), the field and log choices, one plot per unit
// of the fields shown, the CSV of the arrays as returned, and the pyazr call that makes the same arrays.  The Physics
// tab and the THM tab draw their rows with it; the rules (requests, plots, CSV) are physicsform.js.

import { pathOf } from "./api.js";
import { answerCsv, physicsPlots, requestOf, statusLine, unitGroups } from "./physicsform.js";
import { formatNumber, renderTable } from "./tables.js";
import { button, download, drawPlot, h } from "./ui.js";

/** Send the request of the state's texts to the row's route (`row.source.call`), then draw the answer in `out` (or
 * the server's refusal).  `where`: the words of the status line while it computes. */
export async function computeRow(ctx, out, row, state, { where = "the server", refusalNote = "" } = {}) {
  const request = requestOf(row, state.axis, state.texts);
  out.textContent = "";
  if (request.error) { out.appendChild(h("div.error-box", { role: "alert" }, request.error)); return null; }
  out.appendChild(h("p.hint", {}, `computing in ${where}…`));
  let answer;
  try {
    answer = await ctx.api.request("POST", pathOf(row.source.call, row.name), { arguments: request.arguments },
      { quiet: true });
  } catch (e) {
    out.textContent = "";
    ctx.showError(e, out);
    if (refusalNote) out.appendChild(h("p.hint", {}, refusalNote));
    return null;
  }
  state.answer = { answer, request };
  drawAnswer(ctx, out, row, state);
  return answer;
}

/** A value of a fact or a table cell: a number as JavaScript writes it (every digit), a text as it is. */
const shown = (v) => (typeof v === "number" ? formatNumber(v) : (v === null || v === undefined ? "" : String(v)));

/** The facts of an answer (name, value, unit) as a definition list. */
function factsList(facts) {
  return h("dl.facts.function-facts", {}, ...facts.flatMap((f) => [
    h("dt", { title: f.doc || "" }, f.unit ? `${f.name} (${f.unit})` : f.name),
    h("dd.mono", { data: { fact: f.name } }, shown(f.value))]));
}

/** One table of an answer: {title, columns: [{key, title, unit}], rows}. */
function answerTable(t) {
  const columns = t.columns.map((c) => ({ key: c.key, title: c.title, unit: c.unit,
    numeric: t.rows.some((r) => typeof r[c.key] === "number"), format: (v) => shown(v) }));
  const table = renderTable(document, { columns, rows: t.rows.map((r, i) => ({ ...r, _k: i })), rowKey: (r) => r._k,
    caption: t.title, withFilter: false, label: t.title });
  table.element.dataset.table = t.title;
  return table.element;
}

/** The status line, facts, tables, field and log choices, plots, exports and the Python text of an answer. */
export function drawAnswer(ctx, out, row, state) {
  const { answer, request } = state.answer;
  out.textContent = "";
  const plotted = Boolean(row.axis);
  if (plotted) {
    const st = statusLine(answer);
    out.appendChild(h("div.physics-status", { role: "status", data: { computed: st.computed,
      elements: st.elements } }, st.text));
  }
  if (answer.facts?.length) out.appendChild(factsList(answer.facts));
  for (const t of answer.tables || []) out.appendChild(answerTable(t));
  if (plotted) drawPlots(ctx, out, row, state, answer, request);
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
    () => download(`${row.name}.py`, answer.python, "text/x-python"))), code);
}

function drawPlots(ctx, out, row, state, answer, request) {
  const redraw = () => drawAnswer(ctx, out, row, state);
  const flag = (label, key) => {
    const input = h("input", { type: "checkbox", "aria-label": label });
    input.checked = !!state[key];
    input.addEventListener("change", () => { state[key] = input.checked; redraw(); });
    return h("label.choice", {}, input, ` ${label}`);
  };
  const fields = row.returns.filter((r) => answer.returns?.[r.name] !== undefined).map((r) => {
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
  for (const plot of physicsPlots(row, answer, request, state)) {
    const area = h("div.plot.physics-plot");
    out.appendChild(area);
    area.dataset.drawn = drawPlot(area, plot, row.name, { usePlotly: !state.table, csv });
    if (citation) area.appendChild(h("p.citation", { data: { citation } }, `Cite: ${citation}`));
  }
  if (!unitGroups(row, state.shown, answer).length) out.appendChild(h("p.hint", {}, "No field is shown."));
}
