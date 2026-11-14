// The Rate tab (#/rate, #/rate/<job id>): the thermonuclear reaction rate N_A<σv> of the project, computed by the
// engine as the program's rate command computes it (pyazr.rate), as a job (kind rate, method engine: POST /api/jobs):
// the reaction (an entrance and an exit pair), the temperature grid (the named grids of GET /api/rate/grids, or the
// program's min / max / step), the parameters (the working x, the project's values, or a finished fit or sampler
// job's best point, whose covariance or chain gives the band), a comparison table of the project folder; the job's
// progress; the rate against T9 on log-log axes with the band shaded and the comparison dashed, the ratios to the
// median below (the job's numbers), the Gamow window per temperature (pyazr.rate.gamow_window) and the table as CSV
// with its header (pyazr.rate.write_rate_table).  The page computes no rate, ratio or window.

import { pathOf } from "../api.js";
import { isFinal, progressView } from "../fitting.js";
import { renderTable } from "../tables.js";
import { button, download, drawPlot, h, section } from "../ui.js";

export async function render(ctx, el, args) {
  const view = ctx.rateView || (ctx.rateView = { entrance: null, exit: null, grid: "table", t_min: "", t_max: "",
    t_step: "", start: "x", band: "none", comparison: "", draws: "100", starlib_label: "azr2", starlib_flag: "",
    reaclib_label: "azr2", reaclib_flag: "n", reaclib_tolerance: "0.05" });
  const [pairs, grids, jobs] = await Promise.all([ctx.api.get("/api/project/pairs"), ctx.api.get("/api/rate/grids"),
    ctx.api.get("/api/jobs")]);
  const particles = pairs.filter((p) => p.kind === 0);
  if (view.entrance === null && particles.length) view.entrance = particles[0].key;
  if (view.exit === null) view.exit = pairs.find((p) => p.key !== view.entrance)?.key ?? -1;
  const rateJobs = jobs.filter((j) => j.kind === "rate");
  const sources = jobs.filter((j) => (j.kind === "fit" || j.kind === "sample") && isFinal(j.state) && j.best?.x);

  const pairName = (p) => `${p.key}: ${p.label || ""}`;
  const entrance = h("select", { "aria-label": "entrance pair" }, particles.map((p) =>
    h("option", { value: p.key }, pairName(p))));
  entrance.value = String(view.entrance);
  const exit = h("select", { "aria-label": "exit pair" }, [h("option", { value: -1 }, "-1: total radiative capture"),
    ...pairs.map((p) => h("option", { value: p.key }, pairName(p)))]);
  exit.value = String(view.exit);
  const reaction = h("span.hint");
  const showReaction = async () => {
    view.entrance = +entrance.value;
    view.exit = +exit.value;
    try {
      reaction.textContent = await ctx.api.get(`/api/project/reaction?entrance=${view.entrance}&exit=${view.exit}`,
        { quiet: true });
    } catch (e) {
      reaction.textContent = e.view ? e.view.what : String(e);
    }
  };
  entrance.addEventListener("change", showReaction);
  exit.addEventListener("change", showReaction);
  showReaction();

  const grid = h("select", { "aria-label": "temperature grid" }, [...Object.entries(grids.grids).map(([k, g]) =>
    h("option", { value: k, title: g.doc }, `${k}: ${g.doc}`)), h("option", { value: "custom" },
    "custom: min to max in steps (GK)")]);
  grid.value = view.grid;
  const num = (key, label) => {
    const input = h("input.num", { type: "text", size: 7, value: view[key], "aria-label": label, placeholder: "GK" });
    input.addEventListener("change", () => { view[key] = input.value.trim(); });
    return h("label", {}, `${label} `, input);
  };
  const custom = h("span.rate-custom", {}, num("t_min", "min"), " ", num("t_max", "max"), " ", num("t_step", "step"));
  const showCustom = () => { view.grid = grid.value; custom.hidden = grid.value !== "custom"; };
  grid.addEventListener("change", showCustom);
  showCustom();

  const start = h("select", { "aria-label": "parameters" }, [h("option", { value: "x" }, "the working x"),
    h("option", { value: "text" }, "the values of the project text"),
    ...sources.map((j) => h("option", { value: j.id }, `the best point of ${j.kind} job ${j.id} (${j.method})`))]);
  start.value = sources.some((j) => j.id === view.start) || ["x", "text"].includes(view.start) ? view.start : "x";
  const band = h("select", { "aria-label": "band" }, [h("option", { value: "none" }, "no band"),
    h("option", { value: "covariance" }, "from the fit's covariance (linear in ln rate)"),
    h("option", { value: "draws" }, "from the sampler's draws (16, 50, 84 %)")]);
  band.value = view.band;
  start.addEventListener("change", () => { view.start = start.value; });
  band.addEventListener("change", () => { view.band = band.value; });
  const draws = h("input.num", { type: "text", size: 6, value: view.draws, "aria-label": "draws" });
  draws.addEventListener("change", () => { view.draws = draws.value.trim(); });
  const comparison = h("input.mono", { type: "text", size: 28, value: view.comparison, "aria-label": "comparison",
    placeholder: "a file of the project folder", title: "two or three columns: T9 (GK), N_A<σv> (cm³/mol/s), "
      + "optionally its uncertainty" });
  comparison.addEventListener("change", () => { view.comparison = comparison.value.trim(); });

  // the libraries' labels: a text box for the label, the format's allowed values for a flag
  const text = (key, label, size) => {
    const input = h("input.mono", { type: "text", size, value: view[key], "aria-label": label });
    input.addEventListener("change", () => { view[key] = input.value.trim(); });
    return h("label", {}, `${label} `, input);
  };
  const choice = (key, label, values) => {
    const select = h("select", { "aria-label": label }, values.map(([v, words]) => h("option", { value: v }, words)));
    select.value = view[key];
    select.addEventListener("change", () => { view[key] = select.value; });
    return h("label", {}, `${label} `, select);
  };
  const libraries = h("div", {}, "STARLIB ", text("starlib_label", "label", 4), " ",
    choice("starlib_flag", "flag", [["", "none"], ["w", "w: weak"], ["g", "g: γ-ray transition"],
      ["v", "v: detailed balance"]]), " · REACLIB ", text("reaclib_label", "set label", 4), " ",
    choice("reaclib_flag", "flag", [["n", "n: non-resonant"], ["r", "r: resonant"], ["w", "w: weak"]]), " ",
    text("reaclib_tolerance", "tolerance", 5));
  const calculate = button("Calculate", "Start the rate job (POST /api/jobs, kind rate)", async () => {
    const settings = { entrance: view.entrance, exit: view.exit, grid: view.grid, band: view.band,
      comparison: view.comparison, draws: view.draws, starlib_label: view.starlib_label,
      starlib_flag: view.starlib_flag, reaclib_label: view.reaclib_label, reaclib_flag: view.reaclib_flag,
      reaclib_tolerance: view.reaclib_tolerance };
    if (view.grid === "custom") Object.assign(settings, { t_min: view.t_min, t_max: view.t_max, t_step: view.t_step });
    const body = { kind: "rate", method: "engine", settings };
    if (view.start === "text") body.x0 = "text";
    else if (view.start !== "x") body.x0 = { job: view.start };
    try {
      const started = await ctx.api.post("/api/jobs", body);
      ctx.status(`rate job ${started.id} started`, "ok");
      ctx.navigate(`#/rate/${started.id}`);
    } catch (e) {
      ctx.showError(e, el);
    }
  });
  const options = h("div", {},
    h("div", {}, h("label", {}, "parameters ", start), " ", h("label", {}, "band ", band), " ",
      h("label", {}, "draws ", draws)),
    libraries,
    h("div", {}, h("label", {}, "compare with ", comparison), ctx.expert ? [" ", calculate] : null));
  el.appendChild(section("Reaction rate",
    h("p.hint", {}, "N_A<σv> integrated by the engine over the entrance pair's c.m. energy from 10⁻⁷ to 50 MeV "
      + "(σ on a grid with a cluster of points around each level, a spline of ln σ against ln E), as the program's "
      + "rate command; a job in its own process."),
    h("div.rate-form", {},
      h("div", {}, h("label", {}, "entrance ", entrance), " ", h("label", {}, "exit ", exit), " ", reaction),
      h("div", {}, h("label", {}, "temperatures ", grid), " ", custom, ctx.expert ? null : [" ", calculate]),
      // the first view: what the task needs above; the parameters, the band, the comparison and the libraries'
      // labels one click away (the expert's depth shows them as they are)
      ctx.expert ? options : h("details.more", {}, h("summary", {}, "More: parameters, uncertainty band, "
        + "comparison, library labels"), options))));

  const shownId = args[0] || rateJobs[rateJobs.length - 1]?.id;
  if (!shownId) {
    el.appendChild(h("p.hint", {}, "No rate yet: choose the reaction and the temperatures, and press Calculate."));
    return;
  }
  let job;
  try {
    job = await ctx.api.get(pathOf("/api/jobs", shownId));
  } catch (e) {
    ctx.showError(e, el);
    return;
  }
  const panel = h("div.rate-job");
  el.appendChild(panel);
  drawJob(ctx, panel, job);
  if (!isFinal(job.state)) follow(ctx, panel, job);
  if (rateJobs.length > 1) {
    el.appendChild(section("Rate jobs", h("ul", {}, rateJobs.slice().reverse().map((j) => h("li", {},
      h("a", { href: `#/rate/${j.id}` }, `${j.id}`), ` ${j.state}`, j.result?.message ? ` — ${j.result.message}`
        : "")))));
  }
}

/** Follow a running rate job with one long poll in flight; at its end draw it again. */
async function follow(ctx, panel, rec) {
  let r = rec;
  while (panel.isConnected && !isFinal(r.state)) {
    try {
      r = await ctx.api.get(`${pathOf("/api/jobs", r.id)}?wait=1`, { background: true });
    } catch (e) {
      ctx.showError(e, panel);
      return;
    }
    if (panel.isConnected) drawJob(ctx, panel, r);
  }
}

/** A rate job: its progress, or its refusal, or its result (plots, Gamow window, table). */
function drawJob(ctx, panel, job) {
  panel.textContent = "";
  const p = progressView(job);
  if (!p.final) {
    panel.appendChild(section(`Rate job ${job.id}`, h("p", {}, `${p.words}: ${p.calls} rate(s) computed`),
      button("Stop", "Stop between two rates of the band", () => ctx.api.post(pathOf("/api/jobs", job.id, "stop"))
        .catch((e) => ctx.showError(e, panel)))));
    return;
  }
  const res = job.result || {};
  if (job.state === "failed" || !res.rate) {
    panel.appendChild(section(`Rate job ${job.id}`, h("div.error-box", { role: "alert" },
      `${job.state}: ${res.message || "no result"}`)));
    return;
  }
  const rate = res.rate;
  const t = rate.T9;
  const traces = [];
  if (res.band) traces.push({ role: "band", name: "1σ band", x: t, y: res.band.high, lower: res.band.low });
  traces.push({ role: "line", name: `N_A<σv> (${rate.reaction})`, x: t, y: rate.rate });
  if (res.comparison) traces.push({ role: "line", name: `comparison: ${res.comparison.path.split("/").pop()}`,
    dash: "dash", x: t, y: res.comparison.at_T9 });
  const top = h("div.plot");
  const bottom = h("div.plot");
  panel.appendChild(section(`Rate job ${job.id}: ${rate.reaction}`, h("p.hint", {}, `Parameters: ${rate.source}. `
    + (res.band ? `Band: ${res.band.how}.` : "No band.")), top, bottom,
    button("Download the table (CSV)", "The table with its header: project, parameters, units",
      () => download(`rate_${job.id}.csv`, res.table)), h("span.hint", {}, ` also written to ${res.csv}`),
    libraryExports(job, res.formats || {})));
  const usePlotly = ctx.route.params.plot !== "table";
  drawPlot(top, { title: `N_A<σv> of ${rate.reaction}`, x: { title: "T9", unit: "GK", log: true },
    y: { title: "N_A<σv>", unit: rate.units.rate, log: true }, traces }, `rate_${job.id}`, { usePlotly });
  const ratio = [];
  if (res.band) ratio.push({ role: "band", name: "1σ band / median", x: t, y: res.ratios.high, lower: res.ratios.low });
  ratio.push({ role: "line", name: "rate / median", x: t, y: res.ratios.rate });
  if (res.ratios.comparison) ratio.push({ role: "line", name: "comparison / median", dash: "dash", x: t,
    y: res.ratios.comparison });
  drawPlot(bottom, { title: "ratio to the median", x: { title: "T9", unit: "GK", log: true },
    y: { title: "ratio", unit: "", log: false }, traces: ratio }, `rate_ratio_${job.id}`, { usePlotly, height: 220 });
  const g = res.gamow;
  const rows = t.map((T, i) => ({ k: i, T9: T, rate: rate.rate[i], E0: g.E0[i], low: g.low[i], high: g.high[i],
    kT: g.kT[i] }));
  const table = renderTable(document, {
    columns: [{ key: "T9", title: "T9 (GK)", numeric: true, digits: 4 },
      { key: "rate", title: `N_A<σv> (${rate.units.rate})`, numeric: true, digits: 6 },
      { key: "kT", title: "kT (MeV)", numeric: true, digits: 4 },
      { key: "E0", title: "Gamow peak E0 (MeV)", numeric: true, digits: 4 },
      { key: "low", title: "window from (MeV)", numeric: true, digits: 4 },
      { key: "high", title: "to (MeV)", numeric: true, digits: 4 }],
    rows, rowKey: (r) => r.k, withFilter: false,
    caption: "The Gamow window E0 ∓ Δ/2 (c.m.) per temperature: the data in that range of the entrance pair's "
      + "energy constrain the rate there (a narrow resonance inside it dominates instead); none for a neutral "
      + "entrance" });
  panel.appendChild(section("Gamow window", table.element));
}

/** The rate in the libraries' formats (pyazr.rate_formats, made by the job): a download each, with what it holds,
 * or the refusal's words (STARLIB needs a band and its 60 temperatures; REACLIB the range 0.01-10 GK). */
function libraryExports(job, formats) {
  const sp = formats.species;
  const box = h("div.rate-formats", {}, sp ? h("p.hint", {}, `${sp.reactants.join(" + ")} → ${sp.products.join(" + ")}`
    + `, chapter ${sp.chapter}, Q = ${sp.Q} MeV (from the project's pairs)`) : null);
  for (const [name, file] of [["starlib", "STARLIB"], ["reaclib", "REACLIB"]]) {
    const f = formats[name];
    if (!f) continue;
    if (f.refused) {
      box.appendChild(h("div.warn-box", {}, h("strong", {}, `${file}: `), f.what || f.refused,
        f.allowed ? h("div", {}, `Allowed: ${f.allowed}`) : null, f.fix ? h("div", {}, `Fix: ${f.fix}`) : null));
      continue;
    }
    box.appendChild(h("div", {}, button(`Download ${file}`, `The rate as a ${file} file`, () =>
      download(`rate_${job.id}_${name}.dat`, f.text, "text/plain")), " ", h("span.hint", {},
      `${f.summary || f.uncertainty || ""}; written to ${f.path}`)));
  }
  return box;
}
