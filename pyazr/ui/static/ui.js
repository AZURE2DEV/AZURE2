// Small helpers of the page's tabs: building elements, showing a refusal of the server (what / allowed / fix and
// the engine's diagnostics), drawing a plot with plotly.js or as a table of numbers, and saving a CSV file.

import { csvOf, plotlyFigure, tableOf } from "./plots.js";
import { renderTable } from "./tables.js";

/** An element: h("div.cls", {attr: v, onclick: f}, child, "text", ...). */
export function h(spec, attrs = {}, ...children) {
  const [tag, ...classes] = spec.split(".");
  const el = document.createElement(tag || "div");
  if (classes.length) el.className = classes.join(" ");
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === undefined || v === null || v === false) continue;
    if (k.startsWith("on") && typeof v === "function") el.addEventListener(k.slice(2), v);
    else if (k === "text") el.textContent = v;
    else if (k === "data") Object.assign(el.dataset, v);  // data-* attributes
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const c of children.flat()) {
    if (c === null || c === undefined || c === false) continue;
    el.appendChild(typeof c === "string" || typeof c === "number" ? document.createTextNode(String(c)) : c);
  }
  return el;
}

/** A button with its keyboard-reachable action and a title that names its key, if any. */
export function button(label, title, run, cls = "") {
  return h(`button${cls ? "." + cls : ""}`, { type: "button", title, onclick: run }, label);
}

/** A section with a heading. */
export function section(title, ...children) {
  return h("section.panel", {}, h("h2", {}, title), ...children);
}

/** The view of a refusal: what was wrong, what is allowed, how to fix it, and the engine's diagnostics. */
export function errorBox(view, { onLink } = {}) {
  const box = h("div.error-box", { role: "alert" });
  box.appendChild(h("div.what", {}, h("strong", {}, `${view.title}${view.status ? ` (${view.status})` : ""}: `),
    view.what));
  if (view.allowed) box.appendChild(h("div.allowed", {}, h("span.k", {}, "Allowed: "), view.allowed));
  if (view.fix) box.appendChild(h("div.fix", {}, h("span.k", {}, "Fix: "), view.fix));
  if (view.known && view.known.length) box.appendChild(h("div.allowed", {}, h("span.k", {}, "Known: "),
    view.known.slice(0, 20).join(", ")));
  for (const d of view.diagnostics || []) box.appendChild(diagnosticLine(d, onLink));
  return box;
}

/** One diagnostic of the engine, with a link to the record or level it names. */
export function diagnosticLine(d, onLink) {
  const where = [d.record, d.field].filter(Boolean).join(".");
  const line = h(`div.diag.${d.severity || "error"}`, {}, h("span.code", {}, d.code), " ",
    where ? h("span.where", {}, where) : null, where ? ": " : "", d.text || "",
    d.line ? h("span.lineno", {}, ` (line ${d.line})`) : null);
  if (onLink && d.record) {
    line.appendChild(h("a", { href: `#/data/${encodeURIComponent(d.record)}`, onclick: (e) => {
      e.preventDefault(); onLink(d);
    } }, " → go to"));
  }
  return line;
}

/** Offer `text` as a file download named `name`. */
export function download(name, text, type = "text/csv") {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = h("a", { href: url, download: name });
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/** The page's theme colours, read from its CSS custom properties (plots follow the theme). */
export function themeColours() {
  const s = getComputedStyle(document.documentElement);
  const v = (n) => s.getPropertyValue(n).trim();
  return { fg: v("--fg"), bg: v("--panel"), grid: v("--grid"), font: v("--font-ui") };
}

/**
 * Draw a plot description into `el`: with plotly.js where the page loaded it, as a table of its numbers otherwise
 * (section 7); "Export data (CSV)" in both cases.  Returns "plotly" or "table".
 */
export function drawPlot(el, plot, fileName, { usePlotly = true } = {}) {
  el.textContent = "";
  const bar = h("div.plot-bar", {}, button("Export data (CSV)", "The numbers drawn, with their units, as a CSV file",
    () => download(`${fileName}.csv`, csvOf(plot))));
  if (plot.note) bar.appendChild(h("span.hint", {}, plot.note));
  const Plotly = globalThis.Plotly;
  if (usePlotly && Plotly && typeof Plotly.newPlot === "function") {
    const area = h("div.plot-area", { role: "img",
      "aria-label": `${plot.title}; ${plot.y.title} against ${plot.x.title}` });
    el.append(area, bar);
    const fig = plotlyFigure(plot, themeColours());
    Plotly.newPlot(area, fig.data, fig.layout, fig.config);
    return "plotly";
  }
  const { columns, rows } = tableOf(plot);
  const cols = columns.map((c, i) => ({ key: String(i), title: c, numeric: true }));
  const why = Plotly ? "the table view was asked for" : "plotly.js is not loaded";
  const table = renderTable(document, { columns: cols, rows: rows.map((r, i) => ({ ...r, _k: i })),
    rowKey: (r) => r._k, caption: `${plot.title} (${why}: the numbers instead)`, withFilter: false });
  el.append(h("div.plot-fallback", {}, table.element), bar);
  return "table";
}
