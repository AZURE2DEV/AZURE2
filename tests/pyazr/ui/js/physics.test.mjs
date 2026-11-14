// The Physics tab's rules (physicsform.js): every row of the engine's function tables becomes a working form from its
// fields alone (a row injected here too: the extensibility claim), several values become broadcast dimensions and
// lines, elements whose status is not 0 are gaps (null) in the plot, the status line counts every status with its
// words, the CSV holds exactly the answer's arrays, the bound is the server's.
import { readFileSync } from "node:fs";
import { test } from "node:test";
import assert from "node:assert/strict";
import { mod } from "./common.mjs";

const p = await mod("physicsform.js");
const forms = await mod("forms.js");
const FUNCTIONS = process.env.AZURE2_UI_FUNCTIONS
  ? JSON.parse(readFileSync(process.env.AZURE2_UI_FUNCTIONS, "utf8")) : null;
const ctx = { pairs: [{ key: 1, label: "p+12C" }, { key: 2, label: "g+13N" }], experiments: ["LC09", "E"] };

/** A row of a kind no table has yet, written here: it must get a working form and plot with no code. */
const INJECTED = {
  name: "made_up_function", group: "a new group", summary: "a row injected by the test", doc: "", library: "none",
  axis: "x",
  arguments: [
    { name: "k", kind: "integer", unit: "", doc: "a whole number", domain: "integer >= 0", default: null, example: 2,
      sweep: { lower: 0, upper: 4, spacing: "linear" }, choices: [] },
    { name: "x", kind: "real", unit: "MeV", doc: "an energy", domain: "finite > 0", default: null, example: 1,
      sweep: { lower: 0.1, upper: 10, spacing: "log" }, choices: [] },
    { name: "how", kind: "choice", unit: "", doc: "a method", domain: "one of a, b", default: "a", example: null,
      sweep: null, choices: ["a", "b"] },
    { name: "on", kind: "flag", unit: "", doc: "a flag", domain: "0 or 1", default: false, example: null,
      sweep: null, choices: [] },
  ],
  returns: [{ name: "u", unit: "MeV", doc: "" }, { name: "v", unit: "", doc: "" }, { name: "w", unit: "MeV", doc: "" }],
};

const rows = () => [...(FUNCTIONS ? FUNCTIONS.tables.flat() : []), INJECTED];

function checkForm(row, axis) {
  const texts = p.startTexts(row, axis, ctx);
  const form = p.rowForm(row, axis, texts, ctx);
  const entries = form.lines.flat();
  const keys = entries.map((e) => e.key);
  for (const a of row.arguments) assert.ok(keys.includes(a.name), `${row.name}: no field for ${a.name}`);
  const integerAxis = row.arguments.find((a) => a.name === axis).kind === "integer";
  assert.equal(keys.includes(p.POINTS), !integerAxis, `${row.name}: points field`);
  assert.ok(keys.includes(p.SPACING), `${row.name}: spacing field`);
  for (const e of entries) {
    assert.ok(forms.WIDGETS[e.type], `${row.name}.${e.key}: value type ${e.type} has no widget`);
    if (e.required) assert.notEqual((texts[e.key] ?? "").trim(), "", `${row.name}.${e.key}: required and empty`);
    // the widget shows the start text: composing its controls gives the text back
    const w = e.widgets[0];
    assert.equal(w.compose(forms.valuesOf(w)), (texts[e.key] ?? "").trim(), `${row.name}.${e.key}`);
  }
  const req = p.requestOf(row, axis, texts);
  assert.equal(req.error, "");
  const ax = req.arguments[axis];
  assert.ok(typeof ax.from === "number" && typeof ax.to === "number" && ax.from < ax.to, `${row.name}: axis ${ax}`);
  if (!integerAxis) assert.equal(ax.points, p.DEFAULT_POINTS);
  for (const a of row.arguments) {
    // a required argument (a row's own `required`, else one without a default) is always sent
    if (a.required ?? (a.default === null || a.default === undefined)) {
      assert.ok(a.name in req.arguments, `${row.name}: ${a.name}`);
    }
  }
  assert.equal(req.lines.length, 1);
}

/** A row answered as a table (no axis): a field per argument, a request of its arguments alone. */
function checkTableForm(row) {
  const texts = p.startTexts(row, null, ctx);
  const keys = p.rowForm(row, null, texts, ctx).lines.flat().map((e) => e.key);
  assert.deepEqual(keys, row.arguments.map((a) => a.name), `${row.name}: one field per argument`);
  const req = p.requestOf(row, null, texts);
  assert.equal(req.error, "");
  assert.deepEqual(Object.keys(req.arguments), row.arguments.map((a) => a.name), row.name);
}

test("every row of the tables (and an injected one) gets a working form from its fields alone, on every axis", () => {
  if (!FUNCTIONS) return;  // run by node_test.py, which writes the engine's tables
  for (const row of rows()) {
    if (!row.axis) { checkTableForm(row); continue; }
    const axes = p.axisChoices(row);
    assert.equal(axes[0].name, row.axis, `${row.name}: the row's axis first`);
    for (const a of axes) checkForm(row, a.name);
  }
});

test("a widget-making row of ARGUMENT_KINDS for every kind of argument the tables have", () => {
  if (!FUNCTIONS) return;
  const kinds = new Set(rows().flatMap((r) => r.arguments.map((a) => a.kind)));
  for (const k of kinds) assert.ok(p.ARGUMENT_KINDS[k], `argument kind ${k} has no row in ARGUMENT_KINDS`);
});

test("a text argument (a stopping row's formula) is one value, sent as typed", () => {
  if (!FUNCTIONS) return;
  const row = FUNCTIONS.tables.flat().find((r) => r.arguments.some((a) => a.kind === "text"));
  assert.ok(row && row.citation, "a stopping row with its citation");
  const text = row.arguments.find((a) => a.kind === "text");
  const req = p.requestOf(row, row.axis, { ...p.startTexts(row, row.axis), [text.name]: " Ti N " });
  assert.equal(req.arguments[text.name], "Ti N");
  assert.equal(req.lines.length, 1);
});

test("the bound of a request is the server's", () => {
  if (!FUNCTIONS) return;
  assert.equal(p.MAX_ELEMENTS, FUNCTIONS.max_elements);
});

test("several values: a broadcast dimension each, in the row's order, the axis last; one line per combination", () => {
  const texts = { ...p.startTexts(INJECTED, "x"), k: "0 1 2", how: "a b" };
  const req = p.requestOf(INJECTED, "x", texts);
  assert.deepEqual(req.arguments.k, [[[0]], [[1]], [[2]]]);
  assert.deepEqual(req.arguments.how, [["a"], ["b"]]);
  assert.deepEqual(req.arguments.x, { from: 0.1, to: 10, spacing: "log", points: p.DEFAULT_POINTS });
  assert.equal(req.arguments.on, 0);
  assert.deepEqual(req.lines.map((l) => l.label), ["k = 0, how = a", "k = 0, how = b", "k = 1, how = a",
    "k = 1, how = b", "k = 2, how = a", "k = 2, how = b"]);
  const many = p.requestOf(INJECTED, "x", { ...texts, k: Array.from({ length: 13 }, (_, i) => i).join(" ") });
  assert.match(many.error, /26 lines/);
  // an integer axis has no points; a text that is not a number goes to the server as typed (it refuses it)
  const int = p.requestOf(INJECTED, "k", { ...p.startTexts(INJECTED, "k"), x: "one" });
  assert.deepEqual(int.arguments.k, { from: 0, to: 4, spacing: "integer" });
  assert.equal(int.arguments.x, "one");
});

// an answer of shape [2, 3] (two lines of k, three points of x): status 2 and -1 at two elements; the second holds
// a number (the session's functions send 0 where a channel is closed): a gap all the same
const answer = {
  shape: [2, 3], axis: "x",
  arguments: { k: [[0], [1]], x: [0.5, 1, 2], how: "a", on: false },
  returns: { u: [[1, 2, "nan"], [0, 0.25, 0.125]], v: [[0, 0, "nan"], ["nan", 3, 4]], w: [[5, 6, "nan"],
    ["nan", 7, 8]] },
  status: [[0, 0, 2], [-1, 0, 0]], status_words: { 0: "computed", 2: "x outside its domain (finite > 0)",
    "-1": "not converged" },
  units: { u: "MeV", v: "", w: "MeV" },
};
const request = { axis: "x", lines: [{ label: "k = 0" }, { label: "k = 1" }] };

test("elements not computed are gaps (null), never zeros; one plot per unit, a trace per field and line", () => {
  const plots = p.physicsPlots(INJECTED, answer, request, { shown: new Set(["u", "v", "w"]) });
  assert.deepEqual(plots.map((x) => x.y.unit), ["MeV", ""]);
  assert.deepEqual(plots[0].traces.map((t) => t.name), ["u, k = 0", "u, k = 1", "w, k = 0", "w, k = 1"]);
  assert.deepEqual(plots[0].traces[0].y, [1, 2, null]);
  assert.deepEqual(plots[0].traces[1].y, [null, 0.25, 0.125]);
  assert.deepEqual(plots[1].traces[1].y, [null, 3, 4]);
  assert.deepEqual(plots[0].traces[0].x, [0.5, 1, 2]);
  assert.equal(plots[0].x.unit, "MeV");
  assert.match(plots[0].note, /2 elements not computed are gaps/);
  // a field hidden: its unit's plot goes when none of its fields is shown
  const only = p.physicsPlots(INJECTED, answer, request, { shown: new Set(["v"]), logy: true });
  assert.equal(only.length, 1);
  assert.match(only[0].note, /2 values ≤ 0 are not drawn on the log y axis/);
});

test("the status line: computed of all, then each other status's words with its count", () => {
  const s = p.statusLine(answer);
  assert.equal(s.computed, 4);
  assert.equal(s.elements, 6);
  assert.equal(s.text, "4 of 6 elements computed; x outside its domain (finite > 0): 1; not converged: 1");
  assert.equal(p.statusLine({ status: 0, status_words: { 0: "computed" } }).text, "1 of 1 elements computed");
});

test("the CSV is the answer's arrays: one line per element, arguments broadcast, 17 digits, nan, status", () => {
  const lines = p.answerCsv(INJECTED, answer, ["a note"]).trim().split("\n");
  assert.equal(lines[0], "# a note");
  assert.equal(lines[1], "k,x (MeV),how,on,u (MeV),v,w (MeV),status");
  assert.equal(lines.length, 2 + 6);
  const cells = lines.slice(2).map((l) => l.split(","));
  assert.deepEqual(cells.map((c) => Number(c[0])), [0, 0, 0, 1, 1, 1]);
  assert.deepEqual(cells.map((c) => Number(c[1])), [0.5, 1, 2, 0.5, 1, 2]);
  assert.deepEqual(cells.map((c) => c[2]), ["a", "a", "a", "a", "a", "a"]);
  assert.deepEqual(cells.map((c) => c[4]), [1, 2, "nan", 0, 0.25, 0.125].map((v) => (typeof v === "number"
    ? v.toPrecision(17) : v)));
  assert.deepEqual(cells.map((c) => Number(c[7])), [0, 0, 2, -1, 0, 0]);
  // every number reads back as the double sent
  const u = [1, 2, null, 0, 0.25, 0.125];
  cells.forEach((c, i) => { if (u[i] !== null) assert.equal(Number(c[4]), u[i]); });
  assert.deepEqual(p.valueAt([[0], [1]], [1, 2]), 1);
  assert.deepEqual(p.valueAt([0.5, 1, 2], [1, 2]), 2);
});

test("hash parameters set the view: v.<argument>, axis, from, to, points, spacing, fields, logx, logy, table", () => {
  const state = p.newState(INJECTED);
  assert.equal(state.axis, "x");
  assert.equal(state.logx, true);
  p.applyParams(INJECTED, state, { "v.k": "0 3", from: "-1", points: "40", spacing: "linear", fields: "u,w",
    logy: "1", plot: "table", "v.nothing": "1" });
  assert.equal(state.texts.k, "0 3");
  assert.equal(state.texts.x, "-1 10");
  assert.equal(state.texts[p.POINTS], "40");
  assert.equal(state.texts[p.SPACING], "linear");
  assert.deepEqual([...state.shown], ["u", "w"]);
  assert.ok(state.logy && state.table && !("nothing" in state.texts));
  p.applyParams(INJECTED, state, { axis: "k" });
  assert.equal(state.axis, "k");
  assert.equal(state.texts[p.SPACING], "integer");
});

test("the list: rows grouped as the tables group them, in their order; the search filters", () => {
  const list = p.groupedRows([{ name: "a1", group: "A" }, { name: "b1", group: "B", summary: "Coulomb" },
    { name: "a2", group: "A" }]);
  assert.deepEqual(list.map((g) => [g.group, g.rows.map((r) => r.name)]), [["A", ["a1", "a2"]], ["B", ["b1"]]]);
  assert.deepEqual(p.groupedRows([{ name: "a1", group: "A" }, { name: "b1", group: "B", summary: "Coulomb" }],
    "coul").map((g) => g.group), ["B"]);
});

// A THM diagnostic's answer: fields over channels and levels have lines of their own (one list per line), a field on
// the energies alone has none; one status per energy, which every line of a field shares.
const thmRow = {
  name: "made_up_vertex", axis: "energies", returns: [{ name: "M2", unit: "", doc: "" },
    { name: "q", unit: "fm^-1", doc: "" }, { name: "absent", unit: "", doc: "" }],
  arguments: [{ name: "experiment", kind: "experiment", unit: "", doc: "", domain: "", default: null, example: null,
    sweep: null, choices: [] }, { name: "energies", kind: "real", unit: "MeV", doc: "", domain: "finite",
    default: null, example: 0.5, sweep: { lower: 0.05, upper: 2, spacing: "linear" }, choices: [] }],
};
const thmAnswer = {
  shape: [3], axis: "energies", arguments: { experiment: "E", energies: [0.1, 0.2, 0.3] },
  returns: { M2: [[1, 2, 3], [4, 5, 6]], q: [0.5, 0.25, 0.125] },
  lines: { M2: ["1/2+ level 1", "1/2+ level 2"] }, status: [0, -3, 0],
  status_words: { 0: "computed", "-3": "beyond the reach" },
};

test("an experiment argument: a choice of the project's experiments, the first at start, sent as one text", () => {
  const texts = p.startTexts(thmRow, "energies", ctx);
  assert.equal(texts.experiment, "LC09");
  const entry = p.rowForm(thmRow, "energies", texts, ctx).lines.flat().find((e) => e.key === "experiment");
  assert.equal(entry.type, "enum");
  assert.deepEqual(entry.field.values, ["LC09", "E"]);
  const req = p.requestOf(thmRow, "energies", { ...texts, experiment: " E " });
  assert.equal(req.arguments.experiment, "E");
  assert.equal(req.lines.length, 1);
});

test("a field with lines: a trace per line, the energy's status a gap on every line; absent fields go", () => {
  const req = p.requestOf(thmRow, "energies", p.startTexts(thmRow, "energies", ctx));
  const plots = p.physicsPlots(thmRow, thmAnswer, req, { shown: new Set(["M2", "q", "absent"]) });
  assert.deepEqual(plots.map((x) => x.y.unit), ["", "fm^-1"]);
  assert.deepEqual(plots[0].traces.map((t) => t.name), ["M2, 1/2+ level 1", "M2, 1/2+ level 2"]);
  assert.deepEqual(plots[0].traces.map((t) => t.y), [[1, null, 3], [4, null, 6]]);
  assert.deepEqual(plots[1].traces.map((t) => t.y), [[0.5, null, 0.125]]);
  assert.deepEqual(p.unitGroups(thmRow, null, thmAnswer).map((g) => g.fields), [["M2"], ["q"]]);
});

test("the CSV of an answer with lines: a line per energy, a column per field and line, the status", () => {
  const lines = p.answerCsv(thmRow, thmAnswer).trim().split("\n");
  assert.equal(lines[0], "energies (MeV),experiment,M2 [1/2+ level 1],M2 [1/2+ level 2],q (fm^-1),status");
  const n = (v) => (typeof v === "number" ? v.toPrecision(17) : v);  // the CSV's 17 significant digits
  assert.deepEqual(lines.slice(1).map((l) => l.split(",")), [[0.1, "E", 1, 4, 0.5, 0], [0.2, "E", 2, 5, 0.25, -3],
    [0.3, "E", 3, 6, 0.125, 0]].map((r) => r.map(n)));
});
