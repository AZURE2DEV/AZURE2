// The ENDF tab's pure parts (tabs/endf.js): the plots of a data set of the curves' answer (which numbers of the
// answer become which trace, unchanged; a colour per file in the order shown; the ratios only for files that gave
// values; a file's refusal draws nothing), the rows of a check's findings, the count line and the routes.
import { test } from "node:test";
import assert from "node:assert/strict";
import { mod } from "./common.mjs";

const e = await mod("tabs/endf.js");
const p = await mod("plots.js");
const schema = {
  dimensions: { energy: { engine_unit: "MeV" }, angle: { engine_unit: "deg" } },
  quantities: { q_d: { title: "Differential cross section", columns: { angle_axis: "required" } } },
};
const entry = {
  name: "s", quantity: "q_d", engine_unit: "b/sr", e_cm: [0.3, 0.1, 0.2], e_x: [7.3, 7.1, 7.2], angle_cm: [90, 90, 90],
  value_norm: [3, 1, 2], error_norm: [0.3, 0.1, 0.2], model: [2.9, 1.1, 2.1],
  curves: {
    "a.endf": { values: [6, null, 4], ratio_to_data: [2, null, 2], ratio_to_data_error: [0.2, null, 0.2],
      ratio_to_model: [6 / 2.9, null, 4 / 2.1], coverage: [true, false, true] },
    "b.endf": { refused: { what: "the quantity 'q_d'", allowed: "x", fix: "y" } },
  },
};

test("a data set's plots: data, model and each file's values sorted by x; ratios below, a colour per file", () => {
  const { main, ratio } = e.evaluationPlots(entry, ["b.endf", "a.endf"], schema);
  assert.deepEqual(main.x, { title: "E_cm", unit: "MeV" });
  assert.equal(main.y.unit, "b/sr");
  assert.deepEqual(main.traces.map((t) => [t.role, t.name]), [["data", "data × norm"], ["model", "model"],
    ["line", "a.endf"]]);
  assert.deepEqual(main.traces[0].x, [0.1, 0.2, 0.3]);
  assert.deepEqual(main.traces[0].y, [1, 2, 3]);
  assert.deepEqual(main.traces[0].error, [0.1, 0.2, 0.3]);
  assert.deepEqual(main.traces[2].y, [null, 4, 6]);
  assert.equal(main.traces[2].color, e.fileColour(1));
  assert.equal(main.traces[1].color, p.PALETTE[0]);
  assert.deepEqual(ratio.traces.map((t) => t.name), ["a.endf / data", "a.endf / model"]);
  assert.deepEqual(ratio.traces[0].y, [null, 2, 2]);
  assert.deepEqual(ratio.traces[0].error, [null, 0.2, 0.2]);
  assert.deepEqual(ratio.traces[1].y, [null, 4 / 2.1, 6 / 2.9]);
});

test("a data set no file gives values for: no ratio plot (the refusals are drawn instead)", () => {
  const { ratio } = e.evaluationPlots(entry, ["b.endf"], schema);
  assert.equal(ratio, null);
});

test("an angular distribution (one energy) is drawn against the angle", () => {
  const ang = { ...entry, e_cm: [1, 1, 1], angle_cm: [30, 90, 60] };
  const { main } = e.evaluationPlots(ang, ["a.endf"], schema);
  assert.equal(main.x.title, "θ_cm");
  assert.deepEqual(main.traces[0].x, [30, 60, 90]);
  assert.deepEqual(main.traces[2].y, [6, 4, null]);
});

test("findings: the module's own (line, page) and a code's (record, meaning) as rows; the count line", () => {
  const rows = e.findingRows([{ severity: "error", line: 2, mat: 9, mf: 1, mt: 451, message: "m", page: "41",
    record: 0, meaning: "" }, { severity: "warning", line: 0, mat: 9, mf: 2, mt: 151, message: "AP", page: "",
    record: 5, meaning: "not a radius" }]);
  assert.deepEqual(rows.map((r) => [r.line, r.page, r.record, r.meaning]), [[2, "p. 41", null, ""],
    [null, "", 5, "not a radius"]]);
  assert.equal(e.countText({ error: 1, warning: 0, note: 2 }), "1 error, 0 warnings, 2 notes");
});

test("routes: the files shown and the file checked", () => {
  assert.equal(e.routeOf({ files: ["a.endf", "d/b.endf"] }), "#/endf?files=a.endf,d%2Fb.endf");
  assert.equal(e.routeOf({ file: "d/b.endf", job: "j1" }), "#/endf/check?file=d%2Fb.endf&job=j1");
  assert.equal(e.routeOf(), "#/endf");
  assert.deepEqual(e.filesOf("a.endf, d/b.endf,"), ["a.endf", "d/b.endf"]);
});

test("with the model curve: the model and each file on the curve's grid, the ratios at the points", () => {
  const withCurve = { ...entry, curve: { variable: "energy", e_cm: [0.1, 0.2, 0.3], angle_cm: [90, 90, 90],
    model: [1, 2, 3], branches: [{ start: 0, stop: 3, fixed: "angle_cm", value: 90 }] },
  curves: { ...entry.curves, "a.endf": { ...entry.curves["a.endf"], on_curve: { values: [5, 6, 7] } } } };
  const { main, ratio } = e.evaluationPlots(withCurve, ["a.endf"], schema);
  assert.deepEqual(main.traces.map((t) => [t.role, t.name]), [["model", "s (θ_cm = 90 deg): model, 3 points"],
    ["line", "a.endf"], ["data", "data × norm"]]);
  assert.deepEqual(main.traces[1].y, [5, 6, 7]);
  assert.deepEqual(ratio.traces[0].y, [null, 2, 2]);
});

test("a finding of no section (MF 0, MT 0) reads 'whole file'", () => {
  const rows = e.findingRows([{ severity: "note", line: 0, mat: 9, mf: 0, mt: 0, message: "no sequence numbers",
    page: "", record: 0 }, { severity: "error", line: 3, mat: 9, mf: 1, mt: 451, message: "m", page: "41" }]);
  assert.deepEqual(rows.map((r) => [r.mf, r.mt]), [["whole file", ""], [1, 451]]);
});

test("a new project's report: the pairs and the levels as rows, the values with their units", async () => {
  const proj = await mod("tabs/endfproject.js");
  const { pairs, levels } = proj.reportRows({ pairs: [{ key: 1, mt: 2, q_eV: 0 }], levels: [{ jpi: "3/2-",
    er_eV: 1.5e6, formal_MeV: 2.1, energy_MeV: 2.0, values: [1.5, 0.2], kinds: ["width", "anc"],
    status: "converged", background: true }] });
  assert.deepEqual(pairs, [{ key: 1, mt: 2, q_eV: 0 }]);
  assert.equal(levels[0].values, "1.5 eV, 0.2 fm^-1/2");
  assert.equal(levels[0].background, "above EH: background");
});
