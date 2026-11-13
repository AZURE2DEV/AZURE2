// forms.js: the widget table covers every value type of the engine's schema; every record kind gets a form; the
// widgets give back the text they were given (split into boxes and joined again); the conditions hide fields.
import { test } from "node:test";
import assert from "node:assert/strict";
import { engineSchema, mod } from "./common.mjs";

const forms = await mod("forms.js");
const { LAYOUT } = await mod("layout.js");
const schema = engineSchema();
const skipNoSchema = schema ? false : "the engine's schema was not given (AZURE2_UI_SCHEMA)";

const norm = (t) => forms.tokenize(t).join(" ");
const roundTrip = (field, text, ctx = {}) => {
  const w = forms.widgetOf(field, text, { schema, ...ctx });
  return w.compose(forms.valuesOf(w));
};

test("every value type of the schema has a row of WIDGETS", { skip: skipNoSchema }, () => {
  const missing = Object.keys(schema.value_types).filter((t) => !(t in forms.WIDGETS));
  assert.deepEqual(missing, [], `value types without a widget: ${missing}`);
  const used = new Set();
  const walk = (f) => { used.add(f.type); for (const a of f.alternatives || []) for (const p of a.parts) walk(p); };
  for (const r of Object.values(schema.records)) for (const f of Object.values(r.fields)) walk(f);
  for (const t of used) assert.ok(t in forms.WIDGETS, `the type ${t} of a field has no widget`);
});

test("every record kind gets a form, with and without layout hints", { skip: skipNoSchema }, () => {
  for (const kind of Object.keys(schema.records)) {
    const form = forms.describeForm(schema, kind, {}, { records: [] });
    const keys = form.lines.flat().map((e) => e.key);
    assert.ok(keys.length > 0, kind);
    assert.equal(new Set(keys).size, keys.length, `${kind}: a field twice`);
    for (const k of keys) assert.ok(k in schema.records[kind].fields, `${kind}.${k}`);
    for (const e of form.lines.flat()) {
      // an absent repeatable field starts with no occurrence ("+" adds one): its empty widget must have controls
      const w = e.widgets[0] || forms.widgetOf(e.field, "", { schema, records: [] });
      assert.ok(w.controls.length > 0, `${kind}.${e.key}: no control`);
      assert.equal(e.widgets.length === 0, e.repeatable && e.texts.length === 0, `${kind}.${e.key}`);
    }
  }
  assert.ok(Object.keys(schema.records).some((k) => !LAYOUT[k]), "a kind without hints exists (the generic case)");
});

test("layout hints name only fields of the schema's kinds", { skip: skipNoSchema }, () => {
  for (const [kind, lines] of Object.entries(LAYOUT)) {
    assert.ok(kind in schema.records, `LAYOUT names the kind ${kind}, which the schema lacks`);
    for (const k of lines.flat()) assert.ok(k in schema.records[kind].fields, `LAYOUT ${kind}.${k}`);
  }
});

test("the widgets give back the text of T2's records", { skip: skipNoSchema }, () => {
  const R = schema.records;
  const records = [{ kind: "normalization", name: "D" }, { kind: "energy_scale", name: "G" }];
  const cases = [
    [R.dataset.fields.norm, "1 fit prior 1 5%"], [R.dataset.fields.norm, "D"], [R.dataset.fields.norm, "profiled"],
    [R.dataset.fields.shift, "offset 0 MeV fit prior 0 0.003"], [R.dataset.fields.shift, "G"],
    [R.dataset.fields.energy_window, "0.199999 1.400001 MeV"], [R.dataset.fields.quadrature,
      "points 20 width_multiplier 10 per_width 20"], [R.dataset.fields.columns, "E_lab[MeV] _ value error"],
    [R.energy_scale.fields.function, '"a0*(x-0.8)" MeV'], [R.energy_scale.fields.param, "a0 0 fit prior 0 0.01"],
    [R.gamma_decay.fields.branching, "g+15O[6.7931] 0.23037345432917689 fit"],
    [R.gamma_decay.fields.total_width, "0.0385027 eV fit"], [R.target.fields.stopping, "evaluation"],
    [R.potential.fields.central, "woods_saxon depth 50 MeV radius 3 fm diffuseness 0.6 fm"],
  ];
  for (const [field, text] of cases) assert.equal(roundTrip(field, text, { records }), norm(text), text);
});

test("a text a widget cannot split is kept whole in one box", () => {
  const field = { type: "number_with_unit", grammar: "<number> <unit>", units: ["MeV"] };
  const w = forms.widgetOf(field, "1 2 3 4", {});
  assert.equal(w.controls.length, 1);
  assert.equal(w.compose(forms.valuesOf(w)), "1 2 3 4");
  const unknown = forms.widgetOf({ type: "a_type_the_table_lacks" }, "x y", {});
  assert.equal(unknown.compose(forms.valuesOf(unknown)), "x y");
});

test("switching a structured field's alternative starts its parts empty", { skip: skipNoSchema }, () => {
  const f = schema.records.target.fields.stopping;
  const w = forms.widgetOf(f, "evaluation", { schema });
  const values = { ...forms.valuesOf(w), alt: "table" };
  assert.equal(w.compose(values), "table");
});

test("conditions: field_present and field_is hide; others show with their text", () => {
  const present = { applies_when: { if: "field_present", field: "thm", text: "thm is given" } };
  assert.equal(forms.applies("w", present, {}, {}).show, false);
  assert.equal(forms.applies("w", present, { thm: "LC09" }, {}).show, true);
  const is = { applies_when: { if: "field_is", field: "kind", in: ["layer"], text: "" } };
  assert.equal(forms.applies("a", is, { kind: "infinite" }, {}).show, false);
  assert.equal(forms.applies("a", is, { kind: "layer" }, {}).show, true);
  const other = { applies_when: { if: "fold", text: "the record is folded" } };
  assert.deepEqual(forms.applies("q", other, {}, {}), { show: true, note: "the record is folded" });
  assert.equal(forms.applies("a", is, { kind: "infinite", a: "1" }, {}).show, true, "a field with a value shows");
});

test("changedFields sends only what changed; an emptied field is removed (null)", () => {
  const form = { lines: [[{ key: "a", texts: ["1 MeV"], repeatable: false }, { key: "b", texts: ["x"],
    repeatable: false }], [{ key: "r", texts: ["p 1", "q 2"], repeatable: true }]] };
  assert.deepEqual(forms.changedFields(form, { a: ["1 MeV"], b: [""], r: ["p 1"] }), { b: null, r: ["p 1"] });
  assert.deepEqual(forms.changedFields(form, { a: ["2 MeV"], b: ["x"], r: ["p 1", "q 2"] }), { a: "2 MeV" });
});
