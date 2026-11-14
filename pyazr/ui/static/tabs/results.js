// The Results tab (docs/upgrade/web_ui.md section 6.5): what a physicist looks at after every edit and every fit.
// The data sets and extrapolations to show (grouped by section as the Data tab groups them), x and y choices, the
// plots (data with errors, the model, residuals under the plot on the same x axis, level markers), angular
// distributions at a chosen energy, the χ² table as the program splits it, the model data, the program's output
// files, and the places of the reaction rate and the bands.
//
// Every number comes from GET /api/evaluation (the Results payload, the engine's numbers at the working x); the page
// selects, sorts and draws them and computes none.  The view's choices are also hash parameters, so a view can be
// linked: #/results?x=e_x&y=s_factor&logy=1&levels=1&show=cap_g0,cap_g1&energy=1.2

import { pathOf } from "../api.js";
import { X_CHOICES, Y_CHOICES, resultsPlots, seriesStyle } from "../plots.js";
import { formatNumber, renderTable } from "../tables.js";
import { button, drawPlot, h, section } from "../ui.js";

/** The view's choices (what is shown, not the project): kept while the page lives, overridden by hash parameters. */
const view = { x: "e_cm", y: "measured", logx: false, logy: false, model: true, residuals: true, levels: false,
  together: false, band: true, shown: null, project: null, energy: null };

const FLAGS = ["logx", "logy", "model", "residuals", "levels", "together", "band"];

/** The hash parameters of the route over the view (x, y, the flags as 0/1, show=a,b, energy). */
export function applyParams(v, params) {
  if (X_CHOICES.some((c) => c.key === params.x)) v.x = params.x;
  if (Y_CHOICES.some((c) => c.key === params.y)) v.y = params.y;
  for (const f of FLAGS) if (params[f] !== undefined) v[f] = params[f] === "1";
  if (params.show !== undefined) v.shown = new Set(params.show.split(",").filter(Boolean));
  if (params.energy !== undefined) v.energy = Number(params.energy);
  return v;
}

/**
 * The series of the payload: data sets, then extrapolations, each with its colour and marker (the same in every plot
 * and in the χ² table), its section (the Data tab's group of its record kind) and its quantity's title.
 */
export function seriesOf(evaluation, records, schema) {
  const kindOf = Object.fromEntries(records.map((r) => [r.name, r.kind]));
  const sections = Object.entries(schema.sections || {});
  const sectionOf = (name) => sections.find(([, kinds]) => kinds.includes(kindOf[name]))?.[0] ?? "";
  const list = [...(evaluation.datasets || []).map((d) => ({ d, extrapolation: false })),
    ...(evaluation.extrapolations || []).map((d) => ({ d, extrapolation: true }))];
  return list.map((s, i) => ({ ...s, ...seriesStyle(i), section: sectionOf(s.d.name),
    title: schema.quantities?.[s.d.quantity]?.title || s.d.quantity }));
}

/** The newest band job that finished, among the project's jobs (GET /api/jobs), or null. */
export function newestBandJob(jobs) {
  return jobs.filter((j) => j.kind === "bands" && j.state === "finished")
    .sort((a, b) => (b.started || 0) - (a.started || 0))[0] || null;
}

/** The series with the band of a band job's result where it has one of the same name and points (`d.band`: lower,
 * upper and the S factor's, as the server sent them); the names of the series it left without one. */
export function withBands(series, result) {
  const byName = Object.fromEntries((result?.bands?.sets || []).map((b) => [b.name, b]));
  const missing = [];
  const out = series.map((s) => {
    const b = byName[s.d.name];
    if (!b || b.lower.length !== (s.d.e_cm || []).length) { missing.push(s.d.name); return s; }
    return { ...s, d: { ...s.d, band: b } };
  });
  return { series: out, missing };
}

/** The c.m. energies at which a data set has points at three or more angles: {energy: [{series, index}]}. */
export function angularEnergies(series) {
  const out = new Map();
  for (const s of series) {
    if (s.extrapolation) continue;
    const at = new Map();
    (s.d.e_cm || []).forEach((e, i) => { if (!at.has(e)) at.set(e, []); at.get(e).push(i); });
    for (const [e, index] of at) {
      if (new Set(index.map((i) => s.d.angle_cm[i])).size < 3) continue;
      if (!out.has(e)) out.set(e, []);
      out.get(e).push({ ...s, index });
    }
  }
  return new Map([...out.entries()].sort((a, b) => a[0] - b[0]));
}

export async function render(ctx, el) {
  let evaluation;
  try {
    evaluation = await ctx.api.get("/api/evaluation", { quiet: true });
  } catch (e) {
    el.appendChild(section("Results", h("p.hint", {}, `No evaluation: ${e.view?.what || e}`)));
    return;
  }
  const [records, levels, jobs] = await Promise.all([ctx.api.get("/api/project/records"),
    ctx.api.get("/api/project/levels"), ctx.api.get("/api/jobs")]);
  const newest = newestBandJob(jobs);
  const bandJob = newest ? await ctx.api.get(pathOf("/api/jobs", newest.id)) : null;
  const banded = withBands(seriesOf(evaluation, records, ctx.schema), bandJob?.result);
  const series = banded.series;
  if (view.project !== ctx.project.path || !view.shown) {
    view.project = ctx.project.path;
    view.shown = new Set(series.map((s) => s.d.name));
  }
  applyParams(view, ctx.route.params);
  const recordOf = Object.fromEntries(records.map((r) => [r.name, r]));
  const marks = levels.filter((l) => l.active !== false).map((l) => ({ label: l.label, energy: l.energy }));

  el.appendChild(headPanel(evaluation, bandJob, banded.missing, ctx));
  const plots = h("div.results-plots");
  const angular = h("div.results-angular");
  const redraw = () => {
    drawPlots(ctx, plots, series, marks);
    drawAngular(ctx, angular, series);
  };
  el.append(choicePanel(ctx, series, redraw), plots, angular);
  redraw();
  el.append(chi2Panel(evaluation, series, recordOf), modelDataPanel(evaluation), outputsPanel(ctx),
    laterPanel());
}

/** Where the numbers come from: the working x; the χ² totals as the program writes them; the band drawn. */
function headPanel(e, bandJob, missing, ctx) {
  const totals = Object.entries(e.totals || {}).map(([k, v]) => `${k} ${v}`).join(" · ");
  const band = bandJob ? h("div.hint.band-source", { data: { job: bandJob.id } }, `Band: job ${bandJob.id}, `
    + `${bandJob.result.message}; on the text of rev ${bandJob.text_rev}${bandJob.text_rev === ctx.project.rev ? ""
      : ` (the project is at rev ${ctx.project.rev}: the band belongs to the fit of that text)`}; the band is the `
    + "model's around the band job's start point (a fit's or a sampler's best point), the line the model at the "
    + "working x: accept that point into the project to see them together"
    + (missing.length ? `; none for ${missing.join(", ")}` : ""))
    : h("div.hint", {}, "No band yet: a band job of the Fit tab (from a fit's covariance or a sampler's chain) "
      + "draws its band here.");
  return h("section.panel.results-head", {}, band,
    h("div", {}, h("strong", {}, "At: the working x"), h("span.hint", {}, " (as the Fit tab's table shows it; a "
      + "job's best point is seen here once accepted into the project; a choice of point here is not built yet)")),
    h("div.mono", { title: "the Total- line of the program's chiSquared.out at this x" }, totals),
    h("div.hint", {}, `objective (what a fit minimises: data, normalizations, shifts and every prior) ${e.objective}`));
}

/** The selector of the series (grouped by section, all / none per group) and the x, y and display choices. */
function choicePanel(ctx, series, redraw) {
  const box = section("Shown");
  const groups = [...new Set(series.map((s) => s.section))];
  for (const g of groups) {
    const members = series.filter((s) => s.section === g);
    const boxes = members.map((s) => {
      const input = h("input", { type: "checkbox", "aria-label": `show ${s.d.name}` });
      input.checked = view.shown.has(s.d.name);
      input.addEventListener("change", () => {
        if (input.checked) view.shown.add(s.d.name); else view.shown.delete(s.d.name);
        redraw();
      });
      return h("label.series", { title: `${s.title} (${s.d.engine_unit})` }, input,
        h("span.swatch", { style: `color: ${s.color}` }, "●"), ` ${s.d.name}`);
    });
    const setAll = (on) => () => {
      for (const s of members) { if (on) view.shown.add(s.d.name); else view.shown.delete(s.d.name); }
      boxes.forEach((b) => { b.querySelector("input").checked = on; });
      redraw();
    };
    box.appendChild(h("div.series-group", {}, h("span.group-title", {}, g || "other"),
      button("all", `show every record of ${g}`, setAll(true), "quiet"),
      button("none", `show none of ${g}`, setAll(false), "quiet"), ...boxes));
  }
  const select = (label, choices, key) => {
    const s = h("select", { "aria-label": label });
    for (const c of choices) s.appendChild(h("option", { value: c.key, title: c.doc || "" }, c.title));
    s.value = view[key];
    s.addEventListener("change", () => { view[key] = s.value; redraw(); });
    return h("label.choice", {}, `${label} `, s);
  };
  const flag = (label, key, title) => {
    const input = h("input", { type: "checkbox", "aria-label": label });
    input.checked = !!view[key];
    input.addEventListener("change", () => { view[key] = input.checked; redraw(); });
    return h("label.choice", { title }, input, ` ${label}`);
  };
  box.appendChild(h("div.choices", {},
    select("x", X_CHOICES, "x"), flag("log x", "logx"), select("y", Y_CHOICES, "y"), flag("log y", "logy"),
    flag("model", "model", "the model at the data points, and extrapolations on their grids"),
    flag("residuals", "residuals", "the engine's standardized residuals under the plot, on the same x axis"),
    flag("levels", "levels", "the levels of the project text at their excitation energies (x = E_x)"),
    flag("together", "together", "the series shown in one plot per y unit, instead of one plot per series"),
    flag("band", "band", "the band of the newest band job (the Fit tab's), shaded around the model")));
  box.appendChild(h("p.hint", {}, "Not offered: the lab energy (the engine's lab energy of each point is not in "
    + "pyazr's arrays) and the cross section of a data set measured as an S factor (the program's conversion is "
    + "not bound): docs/upgrade/web_ui.md section 12. The S factor is offered where the engine unit is a cross "
    + "section and the entrance pair is charged."));
  return box;
}

const shownSeries = (series) => series.filter((s) => view.shown.has(s.d.name));

/** One plot per unit of the chosen y, each with its residuals, legend, CSV, SVG and PNG. */
function drawPlots(ctx, el, series, marks) {
  el.textContent = "";
  const chosen = shownSeries(series);
  const options = { schema: ctx.schema, levels: marks };
  const plots = view.together ? resultsPlots(chosen, view, options)
    : chosen.flatMap((s) => resultsPlots([s], view, { ...options, title: `${s.d.name}: ${s.title}` }));
  const panel = section("Plots");
  el.appendChild(panel);
  if (!plots.length) {
    panel.appendChild(h("p.hint", {}, chosen.length ? "None of the series shown has the chosen y (the S factor is "
      + "offered for cross sections with a charged entrance)." : "No series chosen: tick a data set or an "
      + "extrapolation above."));
  }
  for (const plot of plots) {
    const area = h("div.plot.results-plot");
    panel.appendChild(area);
    const file = `${plot.title.split(":")[0]}_${view.y}_${view.x}`.replace(/[^A-Za-z0-9_.-]+/g, "_");
    const height = (plot.residual ? 140 : 0) + (view.together ? 420 : 360);
    area.dataset.drawn = drawPlot(area, plot, file, { usePlotly: ctx.route.params.plot !== "table", height });
  }
  const orders = chosen.filter((s) => ctx.schema.quantities?.[s.d.quantity]?.columns?.order);
  panel.appendChild(h("p.hint", {}, "The model line of a data set joins the engine's values at the data points: a "
    + "resonance narrower than the spacing of the points is not drawn. An extrapolation is drawn on its grid (an "
    + "adaptive grid is refined by the engine where the model changes). The model on a dense grid between the data "
    + "points waits for pyazr's Session.evaluate (POST /api/evaluate/grid, docs/upgrade/web_ui.md section 12)."
    + (orders.length ? ` Legendre coefficients (${orders.map((s) => s.d.name).join(", ")}) are drawn as one series: `
      + "the order of each point is not in the data arrays (finding 14)." : "")));
}

/** The angular distribution at a chosen c.m. energy: the shown data sets with points at three or more angles. */
function drawAngular(ctx, el, series) {
  el.textContent = "";
  const energies = angularEnergies(shownSeries(series));
  const panel = section("Angular distributions");
  el.appendChild(panel);
  if (!energies.size) {
    panel.appendChild(h("p.hint", {}, "No data set shown has points at three or more angles at one energy."));
    return;
  }
  if (!energies.has(view.energy)) view.energy = [...energies.keys()][0];
  const s = h("select", { "aria-label": "c.m. energy of the angular distribution" });
  for (const [e, sets] of energies) {
    s.appendChild(h("option", { value: String(e) }, `${formatNumber(e, 7)} MeV (${sets.map((x) => x.d.name)
      .join(", ")})`));
  }
  s.value = String(view.energy);
  s.addEventListener("change", () => { view.energy = Number(s.value); drawAngular(ctx, el, series); });
  panel.appendChild(h("label.choice", {}, "E_cm ", s));
  const plot = resultsPlots(energies.get(view.energy), { x: "angle_cm", y: "measured", model: true,
    residuals: true }, { schema: ctx.schema, title: `angular distribution at E_cm = ${view.energy} MeV` })[0];
  if (!plot) return;
  const area = h("div.plot.results-plot");
  panel.appendChild(area);
  area.dataset.drawn = drawPlot(area, plot, `angular_${view.energy}`,
    { usePlotly: ctx.route.params.plot !== "table", height: 460 });
}

/** The χ² per data set as the program splits it (chiSquared.out), with the sum of the residuals' squares beside it. */
function chi2Panel(e, series, recordOf) {
  const colour = Object.fromEntries(series.map((s) => [s.d.name, s.color]));
  const rows = (e.datasets || []).map((d) => ({ name: d.name, chi2: d.chi2, n: d.n, norm: d.norm,
    norm_chi2: d.norm_chi2, residual_chi2: d.residual_chi2,
    note: recordOf[d.name]?.fields?.covariance ? `data covariance ${recordOf[d.name].fields.covariance}: residuals `
      + "whitened; their squares split the group's χ² differently" : "" }));
  const table = renderTable(document, { rows, withFilter: false, label: "chi-squared per data set",
    caption: "χ² per data set: the program's chiSquared.out at the working x (its six digits)",
    columns: [
      { key: "colour", title: "", render: (r) => h("span.swatch", { style: `color: ${colour[r.name]}` }, "●") },
      { key: "name", title: "data set" },
      { key: "chi2", title: "χ²", numeric: true, doc: "the program's χ² of the data set" },
      { key: "n", title: "N", numeric: true, doc: "points" },
      { key: "norm", title: "norm", numeric: true, doc: "the normalization in use" },
      { key: "norm_chi2", title: "norm χ²", numeric: true, doc: "the normalization's prior term" },
      { key: "residual_chi2", title: "Σ r²", numeric: true, digits: 6,
        doc: "the sum of the squares of the residuals drawn (pyazr's segment_chi2)" },
      { key: "note", title: "note", wrap: true }] });
  const terms = Object.entries(e.terms || {}).filter(([, v]) => v !== 0)
    .map(([k, v]) => h("li", {}, `${k}: ${v}`));
  return section("χ² (the program's split)", table.element,
    h("p.mono", {}, Object.entries(e.totals || {}).map(([k, v]) => `${k} ${v}`).join("   ")),
    h("div.hint", {}, "The objective's non-zero terms (full precision):"), h("ul.terms", {}, terms),
    h("p.hint", {}, `Residual r = ${e.residual?.definition || ""}: ${e.residual?.doc || ""}`));
}

/** The model data (ωγ, widths, ANCs, scattering lengths, ...) against their measured values. */
function modelDataPanel(e) {
  const rows = e.model_data || [];
  if (!rows.length) return section("Model data", h("p.hint", {}, "The project has no model data."));
  const table = renderTable(document, { rows, withFilter: false, label: "model data",
    caption: "model data at the working x: the engine's value of each datum against the value written",
    columns: [{ key: "name", title: "name" }, { key: "quantity", title: "quantity" },
      { key: "model", title: "model", numeric: true, digits: 8 }, { key: "value", title: "value", numeric: true },
      { key: "error", title: "error", numeric: true }, { key: "unit", title: "unit" },
      { key: "residual", title: "(model − value)/error", numeric: true, digits: 5 },
      { key: "chi2", title: "χ²", numeric: true, digits: 6 }, { key: "source", title: "source", wrap: true }] });
  return section("Model data", table.element);
}

/** The program's output files: written at the working x into the project's output directory, listed, first rows. */
function outputsPanel(ctx) {
  const list = h("div.output-files", {}, "…");
  const head = h("pre.file-head", { hidden: true });
  const show = async () => {
    const answer = await ctx.api.get("/api/outputs");
    list.textContent = "";
    if (!answer.files.length) list.appendChild(h("p.hint", {}, `${answer.directory}/ holds no files yet.`));
    for (const f of answer.files) {
      const open = async () => {
        const r = await ctx.api.get(`${pathOf("/api/files", ...answer.directory.split("/"), f.name)}?rows=20`);
        head.hidden = false;
        head.textContent = `${answer.directory}/${f.name}\n${r.lines.join("\n")}${r.more ? "\n…" : ""}`;
      };
      list.appendChild(h("div.output-file", {}, button(f.name, `the first rows of ${f.name}`, open, "quiet"),
        h("span.hint", {}, ` ${f.size} bytes, ${new Date(f.mtime * 1000).toLocaleString()}`)));
    }
  };
  const write = () => ctx.edit("write the program's output files", () => ctx.api.post("/api/outputs", {}));
  show().catch((e) => ctx.showError(e, list));
  return section("Output files", h("p", {}, button("Write the program's output files", "AZUREOut_*, "
    + "chiSquared.out, parameters.out, ... at the working x into the project's output directory", write),
  h("span.hint", {}, " as the program writes them at the end of a run (pyazr's write_output_files)")), list, head);
}

/** What waits for a later step, said once. */
function laterPanel() {
  return section("Reaction rate", h("p.hint", {}, "The reaction rate waits for pyazr's rate call "
    + "(Session.reaction_rate, docs/upgrade/web_ui.md finding 6) and a job kind of the Fit tab."));
}
