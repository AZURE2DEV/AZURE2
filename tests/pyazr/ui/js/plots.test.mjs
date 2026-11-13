// plots.js: the plot table's row for a data set (axis from the quantity row, units from the engine's dimensions,
// one trace of data with errors and one of the model, grouped by angle only where angles repeat); the CSV and the
// table hold exactly the numbers sent; the palette is Okabe and Ito's.
import { test } from "node:test";
import assert from "node:assert/strict";
import { mod } from "./common.mjs";

const p = await mod("plots.js");
const schema = {
  dimensions: { energy: { engine_unit: "MeV" }, angle: { engine_unit: "deg" } },
  quantities: { q_e: { title: "Cross section", engine_unit: "b", columns: { angle_axis: "none", order: false } },
    q_a: { title: "Differential cross section", engine_unit: "b/sr", columns: { angle_axis: "required" } } },
};
const base = { name: "s", value: [3, 1, 2], error: [0.3, 0.1, 0.2], model: [2.9, 1.1, 2.1] };

test("against energy: units from the dimensions, points sorted by energy, numbers as sent", () => {
  const plot = p.plotOf({ ...base, quantity: "q_e", engine_unit: "b", e_cm: [0.3, 0.1, 0.2], angle_cm: [0, 0, 0] },
    { schema });
  assert.equal(plot.x.unit, "MeV");
  assert.equal(plot.y.unit, "b");
  assert.equal(plot.traces.length, 2);
  assert.deepEqual(plot.traces[0].x, [0.1, 0.2, 0.3]);
  assert.deepEqual(plot.traces[0].y, [1, 2, 3]);
  assert.deepEqual(plot.traces[0].error, [0.1, 0.2, 0.3]);
  assert.deepEqual(plot.traces[1].y, [1.1, 2.1, 2.9]);
});

test("an angular distribution (one energy, several angles) is drawn against angle", () => {
  const plot = p.plotOf({ ...base, quantity: "q_a", e_cm: [1, 1, 1], angle_cm: [30, 60, 90] }, { schema });
  assert.equal(plot.x.title, "θ_cm");
  assert.equal(plot.x.unit, "deg");
});

test("excitation functions at a few angles: one trace pair per angle; angles that differ per point: one", () => {
  const n = 12;
  const d = { name: "s", quantity: "q_a", value: Array(n).fill(1), error: Array(n).fill(0.1), model: Array(n).fill(1),
    e_cm: [...Array(n).keys()].map((i) => i % 6), angle_cm: [...Array(n).keys()].map((i) => (i < 6 ? 90 : 150)) };
  assert.equal(p.plotOf(d, { schema }).traces.length, 4);
  const spread = { ...d, angle_cm: [...Array(n).keys()].map((i) => 91 + i * 0.01) };
  assert.equal(p.plotOf(spread, { schema }).traces.length, 2);
});

test("CSV: a header with units, one line per point and trace, the numbers as sent", () => {
  const plot = p.plotOf({ ...base, quantity: "q_e", e_cm: [0.1, 0.2, 0.3], value: [1e-7, 2, 3] }, { schema });
  const lines = p.csvOf(plot).trim().split("\n");
  assert.equal(lines[0], "E_cm (MeV),data (b),data error (b),model at the data points (b)");
  assert.equal(lines[1], "0.1,1e-7,0.3,");
  assert.equal(lines.length, 1 + 3 + 3);
  assert.equal(lines[4], "0.1,,,2.9");
});

test("palette: eight colour-blind-safe colours, then the marker changes", () => {
  assert.equal(p.PALETTE.length, 8);
  assert.notEqual(p.seriesStyle(0).marker, p.seriesStyle(8).marker);
  assert.equal(p.seriesStyle(0).color, p.seriesStyle(8).color);
  assert.ok(p.PLOTS.every((r) => r.name && r.doc && typeof r.build === "function"));
});
