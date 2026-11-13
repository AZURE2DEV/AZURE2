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

// -- the Results tab's plots ------------------------------------------------------------------------------------

/** The x choices of the Results tab: the payload's array, its title and the dimension of its unit. */
export const X_CHOICES = [
  { key: "e_cm", title: "E_cm", dimension: "energy", doc: "c.m. energy of each point (after its energy shift)" },
  { key: "e_x", title: "E_x", dimension: "energy", doc: "excitation energy of the compound nucleus" },
  { key: "angle_cm", title: "θ_cm", dimension: "angle", doc: "c.m. angle of each point" },
];

/** The y choices: the datum times its normalization with its error and the model, in the engine's unit; or the S
 * factor where the server sends it (the engine's conversion of each point). */
export const Y_CHOICES = [
  { key: "measured", title: "as measured",
    of: (d) => ({ value: d.value_norm, error: d.error_norm, model: d.model, unit: d.engine_unit }) },
  { key: "s_factor", title: "S factor",
    of: (d) => (d.s_factor ? { ...d.s_factor, title: "S factor" } : null) },
];

/** The residual's axis title: the engine's standardized residual, model m minus the datum d times the norm n, over
 * the error σ times n (the server's `residual` description says it in full). */
export const RESIDUAL_TITLE = "(m − n·d)/(n·σ)";

/** The groups of one set's points against x: one per value of the other coordinate where it repeats on a few values
 * (an excitation function per angle, an angular distribution per energy), else one group. */
export function groupsAgainst(d, xKey, index) {
  const other = xKey === "angle_cm" ? "e_cm" : "angle_cm";
  const values = distinct(index.map((i) => d[other]?.[i]));
  if (values.length < 2 || values.length > MAX_GROUPS || values.length * 3 > index.length) {
    return [{ value: null, key: other, index }];
  }
  return values.map((v) => ({ value: v, key: other, index: index.filter((i) => d[other][i] === v) }));
}

// the y title: the converted quantity's own (the S factor), else the titles of the quantities drawn
const yTitle = (sets, y) => sets[0]?.y?.title || distinct(sets.map((s) => s.title)).join(", ") || y.title;

/**
 * The plots of the Results tab: one per unit of the chosen y (a cross section and a differential cross section are
 * not drawn on one axis).  `sets`: [{d (the payload's entry), color, marker, extrapolation, title (the quantity's),
 * index (the points to draw, all by default)}]; `view`: {x, y, logx, logy, model, residuals, levels};
 * `levels`: [{label, energy}] (MeV, drawn where x is the excitation energy).  No number is computed: the points are
 * the payload's, sorted by x within a group.
 */
export function resultsPlots(sets, view, { schema, levels = [], title = null } = {}) {
  const xc = X_CHOICES.find((c) => c.key === view.x) || X_CHOICES[0];
  const yc = Y_CHOICES.find((c) => c.key === view.y) || Y_CHOICES[0];
  const xUnit = schema?.dimensions?.[xc.dimension]?.engine_unit ?? "";
  const byUnit = new Map();
  for (const s of sets) {
    const y = yc.of(s.d);
    if (!y) continue;
    if (!byUnit.has(y.unit)) byUnit.set(y.unit, []);
    byUnit.get(y.unit).push({ ...s, y });
  }
  return [...byUnit.entries()].map(([unit, group]) => {
    const traces = [];
    let hasResiduals = false;
    for (const s of group) {
      const d = s.d;
      const all = s.index || [...(d[xc.key] || []).keys()];
      for (const g of groupsAgainst(d, xc.key, all)) {
        const order = [...g.index].sort((a, b) => d[xc.key][a] - d[xc.key][b]);
        const other = X_CHOICES.find((c) => c.key === g.key);
        const suffix = g.value === null ? "" : ` (${other.title} = ${g.value} ${
          schema?.dimensions?.[other.dimension]?.engine_unit ?? ""})`;
        const x = pick(d[xc.key], order);
        const style = { color: s.color, marker: s.marker, group: d.name };
        if (!s.extrapolation && s.y.value) {
          traces.push({ role: "data", name: `${d.name}${suffix}`, x, y: pick(s.y.value, order),
            error: pick(s.y.error, order), ...style });
        }
        if (view.model !== false && s.y.model) {
          traces.push({ role: "model", name: s.extrapolation ? `${d.name}${suffix}: model on its grid`
            : `${d.name}${suffix}: model at the points`, x, y: pick(s.y.model, order), dash: s.extrapolation
            ? "dash" : "solid", ...style });
        }
        if (view.residuals !== false && !s.extrapolation && d.residual) {
          hasResiduals = true;
          traces.push({ role: "residual", name: `${d.name}${suffix}: residual`, x, y: pick(d.residual, order),
            ...style });
        }
      }
    }
    if (view.levels && xc.key === "e_x" && levels.length) {
      traces.push({ role: "levels", name: "levels (project text)", x: levels.map((l) => l.energy),
        labels: levels.map((l) => l.label) });
    }
    return {
      title: title || `${yTitle(group, yc)} against ${xc.title}`,
      x: { title: xc.title, unit: xUnit, log: !!view.logx },
      y: { title: yTitle(group, yc), unit, log: !!view.logy },
      residual: hasResiduals ? { title: RESIDUAL_TITLE, unit: "" } : null,
      traces,
      note: notesOf(traces, view),
    };
  });
}

/** What the plot cannot show as asked: values a log axis drops (counted, not changed). */
export function notesOf(traces, view) {
  const notes = [];
  const count = (role, axis) => traces.filter((t) => role.includes(t.role))
    .reduce((n, t) => n + t[axis].filter((v) => !(typeof v === "number" && v > 0)).length, 0);
  if (view.logy) {
    const n = count(["data", "model"], "y");
    if (n) notes.push(`${n} values ≤ 0 (or not numbers) are not drawn on the log y axis`);
  }
  if (view.logx) {
    const n = count(["data", "model", "residual"], "x");
    if (n) notes.push(`${n} points at x ≤ 0 are not drawn on the log x axis`);
  }
  const points = traces.filter((t) => t.role === "data").reduce((n, t) => n + t.x.length, 0);
  if (points > LARGE) notes.push(`${points} points, all drawn (no decimation): select fewer data sets to speed it up`);
  return notes.join("; ");
}

/** Above this many data points a plot says that it draws them all. */
export const LARGE = 5000;

/** The plot table: (name, what it draws, when it applies, the builder). */
export const PLOTS = [
  { name: "against_energy_or_angle",
    doc: "a data set against energy or angle: data with errors and the model at the working x",
    applies: (d) => d && Array.isArray(d.value) && (Array.isArray(d.e_cm) || Array.isArray(d.angle_cm)),
    build: datasetPlot },
  { name: "results_against_x",
    doc: "the Results tab: the chosen data sets (data with errors, the model at the points), extrapolations (the "
      + "model on their grids), residuals under the plot on the same x axis, level markers; one plot per y unit",
    applies: (d) => d && Array.isArray(d.sets) && d.view,
    build: (d, ctx) => resultsPlots(d.sets, d.view, ctx)[0] || null },
  { name: "angular_distribution",
    doc: "the Results tab: the points of the chosen data sets at one c.m. energy against the c.m. angle, with the "
      + "model at the points",
    applies: (d) => d && Array.isArray(d.sets) && d.energy !== undefined,
    build: (d, ctx) => resultsPlots(d.sets, { x: "angle_cm", y: "measured", model: true, residuals: true },
      { ...ctx, title: `angular distribution at E_cm = ${d.energy} MeV` })[0] || null },
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
    const unit = t.role === "residual" ? plot.residual?.unit : plot.y.unit;
    const cols = t.role === "data" ? [head(t.name, unit), head(`${t.name} error`, unit)]
      : [t.role === "levels" ? t.name : head(t.name, unit)];
    const at = columns.length;
    columns.push(...cols);
    t.x.forEach((x, i) => {
      const row = new Array(at).fill("");
      row[0] = x;
      row.push(t.role === "levels" ? t.labels[i] : t.y[i]);
      if (t.role === "data") row.push(t.error[i]);
      rows.push(row);
    });
  }
  const width = columns.length;
  return { columns, rows: rows.map((r) => r.concat(new Array(width - r.length).fill(""))) };
}

/** A number of the CSV: 17 significant digits (the double as the engine computed it); other values as sent. */
export const csvNumber = (v) => (typeof v === "number" && Number.isFinite(v) ? v.toPrecision(17) : v ?? "");

/** The CSV text of a plot: a header line with units, then one line per point of each trace, every number with its
 * 17 significant digits. */
export function csvOf(plot) {
  const quote = (s) => (/[",\n]/.test(String(s)) ? `"${String(s).replace(/"/g, '""')}"` : String(s));
  const { columns, rows } = tableOf(plot);
  return [columns.map(quote).join(","), ...rows.map((r) => r.map((v) => quote(csvNumber(v))).join(","))].join("\n")
    + "\n";
}

/** plotly.js traces and layout of a plot description, coloured from the page's theme (`theme`: CSS colours).  A
 * plot with `residual` gets a second y axis under the first on the same x axis (one zoom for both); level markers are
 * vertical lines on a hidden axis from 0 to 1 over the plot's height. */
export function plotlyFigure(plot, theme = {}) {
  const fg = theme.fg || "#222", grid = theme.grid || "#ddd", bg = theme.bg || "#fff";
  const muted = theme.muted || "#888";
  const trace = (t) => {
    const common = { type: "scatter", name: t.name, x: t.x, y: t.y, legendgroup: t.group };
    if (t.role === "data") {
      return { ...common, mode: "markers",
        error_y: { type: "data", array: t.error, visible: true, color: t.color || fg, thickness: 1, width: 2 },
        marker: { color: t.color || fg, symbol: t.marker, size: 6, line: { color: fg, width: 0.5 } } };
    }
    if (t.role === "residual") {
      return { ...common, mode: "markers", yaxis: "y2", marker: { color: t.color || fg, symbol: t.marker, size: 5 } };
    }
    if (t.role === "levels") {
      return { type: "scatter", mode: "lines", name: t.name, yaxis: "y3", hoverinfo: "text",
        x: t.x.flatMap((e) => [e, e, null]), y: t.x.flatMap(() => [0, 1, null]),
        text: t.labels.flatMap((l) => [l, l, null]), line: { color: muted, width: 1, dash: "dot" } };
    }
    return { ...common, mode: "lines", line: { color: t.color || fg, width: 2, dash: t.dash || "solid" } };
  };
  const axis = (a) => ({ title: { text: head(a.title, a.unit) }, color: fg, gridcolor: grid, zeroline: false,
    ticks: "outside", mirror: true, showline: true, linecolor: fg, exponentformat: "power",
    type: a.log ? "log" : "linear" });
  const many = plot.traces.length > 4;
  const layout = { title: { text: plot.title, font: { size: 13 }, x: 0, xref: "paper", xanchor: "left" },
    paper_bgcolor: bg, plot_bgcolor: bg,
    font: { color: fg, family: theme.font || "sans-serif", size: 12 }, xaxis: axis(plot.x), yaxis: axis(plot.y),
    margin: { l: 80, r: 20, t: 36, b: 50 }, showlegend: true,
    legend: many ? { x: 1.02, xanchor: "left", y: 1, yanchor: "top", font: { size: 11 } }
      : { x: 1, xanchor: "right", y: 1, yanchor: "top", bgcolor: "rgba(0,0,0,0)", font: { size: 11 } } };
  if (plot.residual) {
    layout.yaxis.domain = [0.32, 1];
    layout.yaxis2 = { ...axis({ ...plot.residual, log: false }), domain: [0, 0.24], anchor: "x", zeroline: true,
      zerolinecolor: muted };
    layout.xaxis.anchor = "y2";
  }
  if (plot.traces.some((t) => t.role === "levels")) {
    layout.yaxis3 = { overlaying: "y", range: [0, 1], visible: false, fixedrange: true };
  }
  return { data: plot.traces.map(trace), layout,
    config: { displaylogo: false, responsive: true, modeBarButtonsToRemove: ["select2d", "lasso2d"] } };
}
