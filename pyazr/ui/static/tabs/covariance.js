// The covariance builder (#/covariance/<data set>), opened by "Covariance…" on a data set of the Data tab: the parts
// of the data set's covariance as rows (kind and its fields, from pyazr.covariance.PART_KINDS as the server sends
// them), Build (POST /api/covariance/build: pyazr builds and checks the matrix), the total error per point beside
// the file's, the correlation on the fitted rows as a heat map, the engine's refusal with the parts it names, the
// notes, and "Write and use" (POST /api/covariance/write: the file and the record's field in one history entry, so
// one Undo takes back both).  The parts are kept in the file the builder writes (its #pyazr.covariance line), and the
// builder reopens with them; the page computes nothing of the matrix.

import { pathOf } from "../api.js";
import { renderTable } from "../tables.js";
import { button, drawPlot, h, section, themeColours } from "../ui.js";

/** The tab is reached from a data set; it has no link of its own in the header. */
export const shown = () => false;

export async function render(ctx, el, args) {
  const name = args[0];
  if (!name) {
    el.appendChild(h("p.hint", {}, "Open the builder from a data set: Data tab, a data set, \"Covariance…\"."));
    return;
  }
  const views = ctx.covarianceViews || (ctx.covarianceViews = {});
  let opened;
  try {
    opened = await ctx.api.get(pathOf("/api/covariance", name));
  } catch (e) {
    ctx.showError(e, el);
    return;
  }
  const view = views[name] || (views[name] = { parts: opened.parts.map((p) => ({ ...p })),
    datasets: opened.datasets.join(", "), path: opened.proposed, summary: null });
  const kinds = opened.kinds;

  el.appendChild(section(`Covariance of ${name}`,
    h("p", {}, "The data set's points enter χ² as vᵀC⁻¹v (v = model − data); C is the sum of the parts below. ",
      h("a", { href: `#/data/${encodeURIComponent(name)}` }, "→ back to the data set")),
    h("p.hint", {}, opened.covariance
      ? `The record names ${opened.covariance}${opened.hand_written ? " (written by hand: it has no parts; writing "
        + "here replaces it)" : "; the builder opened with its parts"}.`
      : "The record names no covariance file yet.")));

  const partsBox = h("div.cov-parts");
  const drawParts = () => {
    partsBox.textContent = "";
    view.parts.forEach((part, i) => partsBox.appendChild(partRow(kinds, view, part, i, drawParts)));
    const add = h("select", { "aria-label": "kind of the new part" }, Object.keys(kinds).map((k) =>
      h("option", { value: k, title: kinds[k].doc }, k)));
    partsBox.appendChild(h("div.cov-add", {}, add, button("+ part", "Add a part of the kind chosen", () => {
      view.parts.push({ kind: add.value });
      drawParts();
    })));
  };
  drawParts();
  const datasets = h("input.mono", { type: "text", size: 40, value: view.datasets, "aria-label": "data sets",
    title: "the data sets this covariance joins (they name one file); separated by commas" });
  datasets.addEventListener("change", () => { view.datasets = datasets.value; });
  const out = h("div.cov-result");
  const list = () => view.datasets.split(/[\s,]+/).filter(Boolean);
  const build = button("Build", "Build and check the matrix (pyazr.covariance.build)", async () => {
    out.textContent = "…";
    try {
      view.summary = await ctx.api.post("/api/covariance/build", { datasets: list(), parts: view.parts });
      out.textContent = "";
      drawSummary(ctx, out, view.summary);
    } catch (e) {
      out.textContent = "";
      ctx.showError(e, out);
    }
  });
  const path = h("input.mono", { type: "text", size: 32, value: view.path, "aria-label": "covariance file",
    title: "the file written, relative to the project's directory" });
  path.addEventListener("change", () => { view.path = path.value.trim(); });
  const write = button("Write and use", "Write the file and set the data sets' covariance field (one Undo takes "
    + "back both)", () => ctx.edit(`covariance ${view.path}`, () => ctx.api.post("/api/covariance/write",
    { datasets: list(), parts: view.parts, path: view.path, rev: ctx.project.rev }), out).catch(() => {}));
  el.appendChild(section("Parts", partsBox, h("div.cov-actions", {}, h("label", {}, "data sets ", datasets), " ",
    build, " ", h("label", {}, "file ", path), " ", write)));
  el.appendChild(out);
  if (view.summary) drawSummary(ctx, out, view.summary);
}

/** One part as a row: its kind, a field per field of the kind (with unit and doc), Remove. */
function partRow(kinds, view, part, i, redraw) {
  const kind = kinds[part.kind];
  const row = h("div.cov-part", {}, h("span.k", { title: kind?.doc || "" }, `${i + 1}. ${part.kind}`), " ");
  for (const f of kind?.fields || []) {
    const value = part[f.name] ?? f.default ?? "";
    let input;
    if (f.type.startsWith("choice:")) {
      input = h("select", { "aria-label": f.name, title: f.doc }, f.type.slice(7).split(",").map((c) =>
        h("option", { value: c }, c)));
      input.value = value;
    } else {
      input = h("input", { type: "text", size: f.type === "path" ? 24 : 8, value, "aria-label": f.name,
        title: f.doc, placeholder: f.default === null ? "required" : "" });
    }
    input.addEventListener("change", () => {
      const v = input.value.trim();
      if (v === "") delete part[f.name]; else part[f.name] = v;
    });
    row.append(h("label", { title: f.doc }, `${f.name}${f.unit ? ` (${f.unit})` : ""} `, input), " ");
  }
  row.appendChild(button("Remove", "Remove this part", () => { view.parts.splice(i, 1); redraw(); }));
  row.appendChild(h("div.hint", {}, kind?.doc || ""));
  return row;
}

/** The built matrix: the refusal or the checks, the notes, the errors per point, the correlation. */
function drawSummary(ctx, el, s) {
  if (s.refusal) {
    el.appendChild(h("div.error-box", { role: "alert" }, h("strong", {}, "The engine would refuse it: "), s.refusal,
      s.refusal_parts.length ? h("div", {}, `Parts concerned: ${s.refusal_parts.map((p) => `${p + 1} (${s.parts[p]
        .kind})`).join(", ")}`) : null));
  } else {
    const c = s.checks;
    el.appendChild(h("div.ok-box", {}, `Positive definite on the ${c.fitted_rows} fitted rows`
      + (c.smallest_eigenvalue !== undefined ? ` (eigenvalues ${c.smallest_eigenvalue.toPrecision(3)} to `
        + `${c.largest_eigenvalue.toPrecision(3)} ${s.unit}²)` : ` (${c.dense || ""})`) + "."));
  }
  if (s.notes.length) el.appendChild(section("What the engine does with it", h("ul", {}, s.notes.map((n) =>
    h("li", {}, n)))));
  const r = s.rows;
  const fitted = r.energy.map((_, i) => i).filter((i) => r.fitted[i]);
  const plotArea = h("div.plot");
  el.appendChild(section("Errors per point", plotArea));
  drawPlot(plotArea, {
    title: "the total error sqrt(C_ii) (data standing in for the model) and the file's error column",
    x: { title: "E_lab", unit: "MeV" }, y: { title: "error", unit: s.unit, log: true },
    traces: [
      { role: "line", name: "total sqrt(C_ii)", x: fitted.map((i) => r.energy[i]), y: fitted.map((i) =>
        r.total_error[i]) },
      { role: "line", name: "error column", dash: "dot", x: fitted.map((i) => r.energy[i]), y: fitted.map((i) =>
        r.error[i]) },
    ],
  }, "covariance_errors", { usePlotly: ctx.route.params.plot !== "table" });
  const table = renderTable(document, {
    columns: [{ key: "set", title: "data set" }, { key: "record", title: "row", numeric: true },
      { key: "energy", title: "E_lab (MeV)", numeric: true, digits: 6 },
      { key: "value", title: `value (${s.unit})`, numeric: true, digits: 6 },
      { key: "error", title: "error column", numeric: true, digits: 4 },
      { key: "total", title: "sqrt(C_ii)", numeric: true, digits: 4 }],
    rows: fitted.map((i) => ({ k: i, set: r.dataset[i], record: r.record[i], energy: r.energy[i],
      value: r.value[i], error: r.error[i], total: r.total_error[i] })),
    rowKey: (row) => row.k, caption: "The fitted rows", withFilter: false });
  el.appendChild(section("Rows", table.element));
  const heat = h("div.plot-area.cov-heat", { role: "img", "aria-label": "correlation on the fitted rows" });
  el.appendChild(section("Correlation on the fitted rows", heat));
  const Plotly = globalThis.Plotly;
  if (!s.correlation) heat.textContent = "Not formed (too many rows).";
  else if (!Plotly || ctx.route.params.plot === "table") heat.textContent = "plotly.js is not loaded: no heat map.";
  else {
    const t = themeColours();
    // the fitted rows in file order on both axes (several data sets or angles can share an energy)
    const k = s.correlation.rows.map((_, j) => j + 1);
    const label = s.correlation.rows.map((i) => `${r.dataset[i]} row ${r.record[i]}, E_lab ${r.energy[i]} MeV`);
    Plotly.newPlot(heat, [{ type: "heatmap", z: s.correlation.matrix, x: k, y: k, zmin: -1, zmax: 1,
      colorscale: "RdBu", reversescale: true, colorbar: { title: { text: "ρ" } },
      text: label.map((a) => label.map((b) => `${a}<br>${b}`)), hoverinfo: "text+z" }],
    { paper_bgcolor: t.bg, plot_bgcolor: t.bg, font: { color: t.fg }, margin: { l: 70, r: 20, t: 20, b: 50 },
      xaxis: { title: { text: "fitted row" } }, yaxis: { title: { text: "fitted row" }, scaleanchor: "x" } },
    { displaylogo: false, responsive: true });
  }
}
