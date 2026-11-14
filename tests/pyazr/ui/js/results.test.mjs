// The Results tab's plot rows (plots.js: resultsPlots, groupsAgainst, notesOf, the figure's residual axis) and the
// payload's mapping (tabs/results.js: seriesOf, angularEnergies, applyParams): which numbers of the payload become
// which trace, unchanged; one plot per y unit; residuals under the plot on the same x axis; level markers only on the
// excitation-energy axis; an S factor only where the server sent one; what a log axis drops is said.
import { test } from "node:test";
import assert from "node:assert/strict";
import { mod } from "./common.mjs";

const p = await mod("plots.js");
const r = await mod("tabs/results.js");
const schema = {
  dimensions: { energy: { engine_unit: "MeV" }, angle: { engine_unit: "deg" } },
  sections: { data: ["k_set"], extrapolations: ["k_grid"] },
  quantities: { q_s: { title: "Cross section", columns: { angle_axis: "none" } },
    q_d: { title: "Differential cross section", columns: { angle_axis: "required" } },
    q_l: { title: "Legendre coefficients", columns: { angle_axis: "none", order: true } } },
};
const set = (name, quantity, unit, extra = {}) => ({ name, quantity, engine_unit: unit, e_cm: [0.3, 0.1, 0.2],
  e_x: [7.3, 7.1, 7.2], angle_cm: [0, 0, 0], value_norm: [3, 1, 2], error_norm: [0.3, 0.1, 0.2],
  model: [2.9, 1.1, 2.1], residual: [-1, 1, 0.5], ...extra });
const payload = {
  datasets: [set("a", "q_s", "b", { s_factor: { model: [29, 11, 21], value: [30, 10, 20], error: [3, 1, 2],
    unit: "MeV*b" } }), set("b", "q_d", "b/sr")],
  extrapolations: [{ name: "g", quantity: "q_s", engine_unit: "b", e_cm: [0.05, 0.4], e_x: [7.05, 7.4],
    angle_cm: [0, 0], model: [0.5, 4] }],
};
const records = [{ name: "a", kind: "k_set" }, { name: "b", kind: "k_set" }, { name: "g", kind: "k_grid" }];
const all = { x: "e_cm", y: "measured", model: true, residuals: true, levels: false };

test("series: data sets then extrapolations, a colour each, grouped by the Data tab's sections", () => {
  const s = r.seriesOf(payload, records, schema);
  assert.deepEqual(s.map((x) => [x.d.name, x.section, x.extrapolation]),
    [["a", "data", false], ["b", "data", false], ["g", "extrapolations", true]]);
  assert.deepEqual(s.map((x) => x.color), p.PALETTE.slice(0, 3));
  assert.equal(s[1].title, "Differential cross section");
});

test("one plot per y unit; data, model and residual traces carry the payload's numbers sorted by x", () => {
  const plots = p.resultsPlots(r.seriesOf(payload, records, schema), all, { schema });
  assert.deepEqual(plots.map((x) => x.y.unit), ["b", "b/sr"]);
  const [first] = plots;
  assert.equal(first.x.unit, "MeV");
  assert.equal(first.residual.title, p.RESIDUAL_TITLE);
  const byRole = (role) => first.traces.filter((t) => t.role === role);
  assert.deepEqual(byRole("data")[0].x, [0.1, 0.2, 0.3]);
  assert.deepEqual(byRole("data")[0].y, [1, 2, 3]);
  assert.deepEqual(byRole("data")[0].error, [0.1, 0.2, 0.3]);
  assert.deepEqual(byRole("model")[0].y, [1.1, 2.1, 2.9]);
  assert.deepEqual(byRole("residual")[0].y, [1, 0.5, -1]);
  const grid = byRole("model").find((t) => t.name.startsWith("g"));
  assert.equal(grid.dash, "dash");
  assert.deepEqual(grid.y, [0.5, 4]);
  assert.ok(!first.traces.some((t) => t.role === "data" && t.name.startsWith("g")), "no data of an extrapolation");
  // every trace of a set has the set's colour (the χ² table uses the same)
  assert.ok(first.traces.filter((t) => t.group === "a").every((t) => t.color === p.PALETTE[0]));
});

test("the S factor: only the series the server sent it for, its own numbers and unit", () => {
  const plots = p.resultsPlots(r.seriesOf(payload, records, schema), { ...all, y: "s_factor" }, { schema });
  assert.equal(plots.length, 1);
  assert.equal(plots[0].y.unit, "MeV*b");
  assert.equal(plots[0].y.title, "S factor");
  assert.deepEqual(plots[0].traces.find((t) => t.role === "data").y, [10, 20, 30]);
  assert.deepEqual(plots[0].traces.find((t) => t.role === "model").y, [11, 21, 29]);
});

test("level markers on the excitation-energy axis only, in the legend, and in the figure on a hidden axis", () => {
  const levels = [{ label: "1/2+ 7.2", energy: 7.2 }];
  const s = r.seriesOf(payload, records, schema);
  assert.ok(!p.resultsPlots(s, { ...all, levels: true }, { schema, levels })[0].traces.some((t) => t.role
    === "levels"));
  const plot = p.resultsPlots(s, { ...all, x: "e_x", levels: true }, { schema, levels })[0];
  const marks = plot.traces.find((t) => t.role === "levels");
  assert.deepEqual(marks.x, [7.2]);
  const fig = p.plotlyFigure(plot);
  const drawn = fig.data.find((t) => t.yaxis === "y3");
  assert.deepEqual(drawn.x, [7.2, 7.2, null]);
  assert.equal(fig.layout.yaxis3.overlaying, "y");
  // residuals on the second axis under the first, one x axis for both (one zoom)
  assert.ok(fig.data.filter((t) => t.yaxis === "y2").length >= 1);
  assert.equal(fig.layout.xaxis.anchor, "y2");
  assert.ok(fig.layout.yaxis.domain[0] > fig.layout.yaxis2.domain[1]);
  assert.ok(fig.data.every((t) => t.name), "every trace has a legend entry");
});

test("log axes: the values a log axis drops are counted in the note, the numbers are not changed", () => {
  const neg = { datasets: [set("a", "q_s", "b", { value_norm: [3, -1, 0] })], extrapolations: [] };
  const plot = p.resultsPlots(r.seriesOf(neg, records, schema), { ...all, logy: true }, { schema })[0];
  assert.match(plot.note, /2 values ≤ 0/);
  assert.equal(plot.y.log, true);
  assert.equal(p.plotlyFigure(plot).layout.yaxis.type, "log");
  assert.deepEqual(plot.traces[0].y, [-1, 0, 3]);
});

test("excitation functions at a few angles against energy, angular distributions against angle", () => {
  const n = 12;
  const d = { ...set("b", "q_d", "b/sr"), e_cm: [...Array(n).keys()].map((i) => i % 6),
    angle_cm: [...Array(n).keys()].map((i) => (i < 6 ? 90 : 150)), value_norm: Array(n).fill(1),
    error_norm: Array(n).fill(0.1), model: Array(n).fill(1), residual: Array(n).fill(0) };
  assert.equal(p.groupsAgainst(d, "e_cm", [...Array(n).keys()]).length, 2);
  const at = r.angularEnergies(r.seriesOf({ datasets: [{ ...d, angle_cm: [...Array(n).keys()].map((i) => 10 * i),
    e_cm: [...Array(n).keys()].map((i) => (i < 4 ? 1.5 : 2.5 + i)) }], extrapolations: [] }, records, schema));
  assert.deepEqual([...at.keys()], [1.5]);
  assert.deepEqual(at.get(1.5)[0].index, [0, 1, 2, 3]);
});

test("the view's hash parameters: x, y, flags, the series shown, the energy; unknown values ignored", () => {
  const v = r.applyParams({ x: "e_cm", y: "measured", logy: false, shown: new Set(["a"]) },
    { x: "e_x", y: "nonsense", logy: "1", show: "a,g", energy: "1.5" });
  assert.equal(v.x, "e_x");
  assert.equal(v.y, "measured");
  assert.equal(v.logy, true);
  assert.deepEqual([...v.shown], ["a", "g"]);
  assert.equal(v.energy, 1.5);
});

test("the plot table has the Results rows, each with its doc", () => {
  const names = p.PLOTS.map((x) => x.name);
  for (const n of ["against_energy_or_angle", "results_against_x", "angular_distribution"]) assert.ok(names
    .includes(n), n);
  const sets = r.seriesOf(payload, records, schema);
  assert.equal(p.plotOf({ sets, view: all }, { schema }).row, "results_against_x");
});

test("bands: the newest finished band job's band on the series of the same name and points, sorted by x", () => {
  const jobs = [{ id: "old", kind: "bands", state: "finished", started: 1 },
    { id: "run", kind: "bands", state: "running", started: 3 }, { id: "new", kind: "bands", state: "finished",
      started: 2 }, { id: "fit", kind: "fit", state: "converged", started: 4 }];
  assert.equal(r.newestBandJob(jobs).id, "new");
  assert.equal(r.newestBandJob([]), null);
  const result = { bands: { how: "linear", sets: [
    { name: "a", lower: [2.8, 0.9, 1.9], upper: [3.0, 1.3, 2.3],
      s_factor: { lower: [28, 9, 19], upper: [30, 13, 23] } },
    { name: "b", lower: [1], upper: [2] }, { name: "g", lower: [0.4, 3.5], upper: [0.6, 4.5] }] } };
  const { series, missing } = r.withBands(r.seriesOf(payload, records, schema), result);
  assert.deepEqual(missing, ["b"], "a band of other points is not drawn");
  const [plot] = p.resultsPlots(series, all, { schema });
  const bands = plot.traces.filter((t) => t.role === "band");
  assert.deepEqual(bands.map((t) => t.group), ["a", "g"]);
  assert.deepEqual(bands[0].x, [0.1, 0.2, 0.3]);
  assert.deepEqual(bands[0].lower, [0.9, 1.9, 2.8]);
  assert.deepEqual(bands[0].y, [1.3, 2.3, 3.0]);
  assert.ok(!p.resultsPlots(series, { ...all, band: false }, { schema })[0].traces.some((t) => t.role === "band"));
  const [sf] = p.resultsPlots(series, { ...all, y: "s_factor" }, { schema });
  assert.deepEqual(sf.traces.find((t) => t.role === "band").y, [13, 23, 30]);
  const fig = p.plotlyFigure(plot);
  const region = fig.data.find((t) => t.fill === "toself");
  assert.deepEqual(region.x, [0.1, 0.2, 0.3, 0.3, 0.2, 0.1]);
  assert.deepEqual(region.y, [1.3, 2.3, 3.0, 2.8, 1.9, 0.9]);
  const csv = p.csvOf(plot).split("\n");
  assert.ok(csv[0].includes("a: band lower (b)") && csv[0].includes("a: band upper (b)"), csv[0]);
});

test("a band's legend names its job and the point it is drawn around; the job whose best point the text is", () => {
  assert.equal(r.bandLabel(null), "band");
  assert.equal(r.bandLabel({ id: "b", result: { source: "f" } }), "band of job b, around the best point of job f");
  const jobs = [{ id: "f", accepted: { rev: 4, time: 1 } }, { id: "g", accepted: { rev: 6, time: 2 } }, { id: "h" }];
  assert.equal(r.acceptedJob(jobs, 6).id, "g");
  assert.equal(r.acceptedJob(jobs, 7), null);
  const s = r.withBands([{ d: { name: "a", e_cm: [1, 2] } }], { bands: { sets: [{ name: "a", lower: [0, 0],
    upper: [1, 1] }] } }, "band of job b");
  assert.equal(s.series[0].d.band.label, "band of job b");
});

test("the model curve: one line per branch along its variable, the server's numbers; the points stay points", () => {
  const curve = { variable: "energy", e_cm: [0.1, 0.15, 0.2, 0.1, 0.2], e_x: [7.1, 7.15, 7.2, 7.1, 7.2],
    angle_cm: [90, 90, 90, 150, 150], model: [1, 5, 2, 3, 4], s_factor: { model: [10, 50, 20, 30, 40] },
    branches: [{ start: 0, stop: 3, fixed: "angle_lab", value: 90 }, { start: 3, stop: 5, fixed: "angle_lab",
      value: 150 }] };
  const s = r.seriesOf({ datasets: [set("a", "q_s", "b", { curve, s_factor: { model: [29, 11, 21],
    value: [30, 10, 20], error: [3, 1, 2], unit: "MeV*b" } })] }, records, schema);
  const [plot] = p.resultsPlots(s, all, { schema });
  const models = plot.traces.filter((t) => t.role === "model");
  assert.deepEqual(models.map((t) => t.name), ["a (θ_lab = 90 deg): model, 3 points",
    "a (θ_lab = 150 deg): model, 2 points"]);
  assert.deepEqual(models.map((t) => [t.x, t.y]), [[[0.1, 0.15, 0.2], [1, 5, 2]], [[0.1, 0.2], [3, 4]]]);
  assert.deepEqual(plot.traces.find((t) => t.role === "data").x, [0.1, 0.2, 0.3]);
  assert.deepEqual(plot.traces.find((t) => t.role === "residual").y, [1, 0.5, -1]);
  const sf = p.resultsPlots(s, { ...all, y: "s_factor" }, { schema })[0].traces.filter((t) => t.role === "model");
  assert.deepEqual(sf[0].y, [10, 50, 20]);
  // against the angle an energy curve is not drawn: the line through the points is
  const byAngle = p.resultsPlots(s, { ...all, x: "angle_cm" }, { schema })[0].traces.filter((t) => t.role === "model");
  assert.deepEqual(byAngle.map((t) => t.name), ["a: model at the points"]);
  // a refused curve: the line through the points
  const refused = r.seriesOf({ datasets: [set("a", "q_s", "b", { curve: { refused: { what: "w" } } })] }, records,
    schema);
  assert.deepEqual(p.resultsPlots(refused, all, { schema })[0].traces.filter((t) => t.role === "model")
    .map((t) => t.name), ["a: model at the points"]);
});
