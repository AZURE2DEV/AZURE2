// The plot table of the page (docs/upgrade/web_ui.md sections 7 and 10): each row turns an answer of the server into
// a plot description (axes with their quantity and unit, traces) that plotly.js draws, or that the page shows as a
// table of numbers where plotly.js is absent.  Every plot can be exported as CSV with the numbers drawn and their
// units.  The axes, titles and units come from the engine's tables (the quantity row of pyazr.schema(), the
// dimensions' engine units); the page converts nothing: the numbers drawn are the server's.
//
// One row today: a data set against energy or angle.  A plot unlike "value against energy or angle" adds a row here
// (and a step of the smoke test).  No DOM in the table itself (tests/pyazr/ui/js/plots.test.mjs).

/** Okabe and Ito's colour-blind-safe palette; after eight series the markers change as well (MARKERS). */
export const PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#7F7F7F"];
export const MARKERS = ["circle", "square", "diamond", "triangle-up", "triangle-down", "x", "cross", "star"];

/** The style of series i: a colour of the palette and, past its length, another marker. */
export function seriesStyle(i) {
  return { color: PALETTE[i % PALETTE.length], marker: MARKERS[Math.floor(i / PALETTE.length) % MARKERS.length] };
}

const distinct = (values) => [...new Set(values)];
/** The most groups (traces of data) a data set is split into; more distinct values draw as one trace. */
export const MAX_GROUPS = 16;

/** The axis a data set is drawn against: angle when the energy is one value and the angle varies, else energy. */
export function axisOf(data, quantityRow, schema) {
  const dims = schema?.dimensions || {};
  const angleAxis = quantityRow?.columns?.angle_axis ?? "none";
  const energies = distinct(data.e_cm || []);
  const angles = distinct(data.angle_cm || []);
  if (angleAxis !== "none" && energies.length === 1 && angles.length > 1) {
    return { key: "angle_cm", title: "θ_cm", unit: dims.angle?.engine_unit ?? "", group: "e_cm",
      groupTitle: "E_cm", groupUnit: dims.energy?.engine_unit ?? "" };
  }
  // several angles, each with several energies (an excitation function per angle): one trace per angle; angles that
  // differ from point to point (a lab angle seen in the c.m. frame) make no groups
  const byAngle = angleAxis !== "none" && angles.length > 1 && angles.length <= MAX_GROUPS
    && angles.length * 3 <= (data.value || []).length;
  return { key: "e_cm", title: "E_cm", unit: dims.energy?.engine_unit ?? "", group: byAngle ? "angle_cm" : null,
    groupTitle: "θ_cm", groupUnit: dims.angle?.engine_unit ?? "" };
}

/** Indices of the points grouped by the values of `key` (all in one group when key is null), in the data's order. */
export function groupsOf(data, key) {
  const n = (data.value || []).length;
  if (!key) return [{ value: null, index: [...Array(n).keys()] }];
  const out = new Map();
  (data[key] || []).forEach((v, i) => { if (!out.has(v)) out.set(v, []); out.get(v).push(i); });
  return [...out.entries()].map(([value, index]) => ({ value, index }));
}

const pick = (arr, idx) => idx.map((i) => (arr || [])[i]);

/** The plot of one data set at the working x: points with errors (data) and a line (the model). */
function datasetPlot(data, { schema }) {
  const q = schema?.quantities?.[data.quantity] || {};
  const axis = axisOf(data, q, schema);
  const groups = groupsOf(data, axis.group);
  const unit = data.engine_unit || q.engine_unit || "";
  const traces = [];
  groups.forEach((g, gi) => {
    // one group: data in the text colour, the model in the palette's first colour; several: a colour per group
    const style = groups.length === 1 ? { color: null, marker: MARKERS[0], model: PALETTE[0] }
      : { ...seriesStyle(gi), model: seriesStyle(gi).color };
    const suffix = g.value === null ? "" : ` ${axis.groupTitle} = ${g.value} ${axis.groupUnit}`;
    const order = [...g.index].sort((a, b) => data[axis.key][a] - data[axis.key][b]);
    traces.push({ role: "data", name: `data${suffix}`, x: pick(data[axis.key], order), y: pick(data.value, order),
      error: pick(data.error, order), color: style.color, marker: style.marker });
    if (data.model) {
      traces.push({ role: "model", name: `model at the data points${suffix}`, x: pick(data[axis.key], order),
        y: pick(data.model, order), color: style.model });
    }
  });
  return {
    title: `${data.name}: ${q.title || data.quantity}${data.observable ? ` (${data.observable})` : ""}`,
    x: { title: axis.title, unit: axis.unit },
    y: { title: q.title || data.quantity, unit },
    traces,
    note: q.columns?.order ? "one trace per order needs the order of each point, which the data route does not send"
      : "",
  };
}

/** The plot table: (name, what it draws, when it applies, the builder). */
export const PLOTS = [
  { name: "against_energy_or_angle",
    doc: "a data set against energy or angle: data with errors and the model at the working x",
    applies: (d) => d && Array.isArray(d.value) && (Array.isArray(d.e_cm) || Array.isArray(d.angle_cm)),
    build: datasetPlot },
];

/** The plot description of an answer: the first row of PLOTS that applies, or null. */
export function plotOf(data, ctx) {
  const row = PLOTS.find((r) => r.applies(data, ctx));
  return row ? { row: row.name, ...row.build(data, ctx) } : null;
}

const head = (title, unit) => (unit ? `${title} (${unit})` : title);

/** The columns of a plot's table and CSV: x, then per trace its y (and the data's errors). */
export function tableOf(plot) {
  const columns = [head(plot.x.title, plot.x.unit)];
  const rows = [];
  for (const t of plot.traces) {
    const cols = t.role === "data" ? [head(t.name, plot.y.unit), head(`${t.name} error`, plot.y.unit)]
      : [head(t.name, plot.y.unit)];
    const at = columns.length;
    columns.push(...cols);
    t.x.forEach((x, i) => {
      const row = new Array(at).fill("");
      row[0] = x;
      row.push(t.y[i]);
      if (t.role === "data") row.push(t.error[i]);
      rows.push(row);
    });
  }
  const width = columns.length;
  return { columns, rows: rows.map((r) => r.concat(new Array(width - r.length).fill(""))) };
}

/** The CSV text of a plot: a header line with units, then one line per point of each trace; numbers as sent. */
export function csvOf(plot) {
  const quote = (s) => (/[",\n]/.test(String(s)) ? `"${String(s).replace(/"/g, '""')}"` : String(s));
  const { columns, rows } = tableOf(plot);
  return [columns.map(quote).join(","), ...rows.map((r) => r.map((v) => quote(v ?? "")).join(","))].join("\n") + "\n";
}

/** plotly.js traces and layout of a plot description, coloured from the page's theme (`theme`: CSS colours). */
export function plotlyFigure(plot, theme = {}) {
  const fg = theme.fg || "#222", grid = theme.grid || "#ddd", bg = theme.bg || "#fff";
  const data = plot.traces.map((t) => (t.role === "data" ? {
    type: "scatter", mode: "markers", name: t.name, x: t.x, y: t.y,
    error_y: { type: "data", array: t.error, visible: true, color: t.color || fg, thickness: 1, width: 2 },
    marker: { color: t.color || fg, symbol: t.marker, size: 6 },
  } : { type: "scatter", mode: "lines", name: t.name, x: t.x, y: t.y, line: { color: t.color || fg, width: 2 } }));
  const axis = (a) => ({ title: { text: head(a.title, a.unit) }, color: fg, gridcolor: grid, zeroline: false,
    ticks: "outside", mirror: true, showline: true, linecolor: fg, exponentformat: "power" });
  return {
    data,
    layout: { title: { text: plot.title, font: { size: 13 }, x: 0, xref: "paper", xanchor: "left" },
      paper_bgcolor: bg, plot_bgcolor: bg,
      font: { color: fg, family: theme.font || "sans-serif", size: 12 }, xaxis: axis(plot.x), yaxis: axis(plot.y),
      margin: { l: 80, r: 20, t: 36, b: 50 }, showlegend: true,
      legend: { x: 1, xanchor: "right", y: 1, yanchor: "top", bgcolor: "rgba(0,0,0,0)", font: { size: 11 } } },
    config: { displaylogo: false, responsive: true, modeBarButtonsToRemove: ["select2d", "lasso2d"] },
  };
}
