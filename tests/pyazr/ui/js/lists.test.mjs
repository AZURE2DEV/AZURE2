// The lists of records in words (lists.js): one list per record kind, heads from the engine's schema and pyazr's
// tables, a data set's description, who uses a record, the pairs' columns.
import { test } from "node:test";
import assert from "node:assert/strict";
import { engineKinds, engineSchema, mod } from "./common.mjs";

const lists = await mod("lists.js");
const tables = await mod("tables.js");
const schema = engineSchema();
const kinds = engineKinds();
const skip = schema && kinds ? false : "the engine's schema or pyazr's tables were not given";

const heads = (columns) => columns.map((c) => tables.columnHead(c));
const cells = (columns, row) => Object.fromEntries(columns.map((c) => [tables.columnHead(c), tables.cellText(row, c)]));

test("a target's list: the schema's titles of the fields its records set, and who uses it", { skip }, () => {
  const rows = [{ name: "tgt", kind: "target", fields: { kind: "layer", active: "15N", composition: "Ta 2 N 5",
    doppler: "free_gas 300 K fixed" } },
    { name: "tgt2", kind: "target", fields: { kind: "infinite", active: "15N" } }];
  const users = lists.usersOf([{ record: "s01", field: "target", name: "tgt", epoch: "" },
    { record: "s04", field: "target", name: "tgt", epoch: "" }]);
  const columns = lists.listColumns({ schema, kind: "target", rows, users });
  const f = schema.records.target.fields;
  assert.deepEqual(heads(columns), ["Name", f.kind.title, f.active.title, f.composition.title, "Used by"]);
  assert.equal(cells(columns, rows[0])[f.active.title], "15N");  // the nucleus under its own head
  assert.equal(cells(columns, rows[0])["Used by"], "s01 (target), s04 (target)");
  assert.equal(cells(columns, rows[1])["Used by"], "");
  assert.ok(schema.records.target.fields.doppler.rare);  // a rare field set stays in the form, not the list
  for (const c of columns.slice(1, 4)) assert.match(c.doc, /\(\w+\)$/);  // the key stays as the hint
});

test("a beam and a detector get their own heads (no union of kinds)", { skip }, () => {
  const beam = lists.listColumns({ schema, kind: "beam", rows: [{ name: "b", fields: { spread: "2 keV lab" } }] });
  assert.deepEqual(heads(beam), ["Name", schema.records.beam.fields.spread.title]);
  const det = lists.listColumns({ schema, kind: "detector", rows: [{ name: "d", fields: { attenuation: "1 0.98" } }] });
  assert.deepEqual(heads(det), ["Name", schema.records.detector.fields.attenuation.title]);
});

test("a data set's list: the description's entries pyazr marks for the list, the norm in use and its status",
  { skip }, () => {
    const rows = [{ name: "s01", kind: "dataset", fields: {} }, { name: "s02", kind: "dataset", fields: {} }];
    const d = (name, norm) => ({ name, tree: "15N(p,a)12C", observable: "angle-integrated cross section", unit: "b",
      data_file: `data/${name}.dat`, kept_E_lab_MeV: [0.2, 1.4], rows_kept: 40, active: true,
      parameters: [{ name: norm, kind: "norm", status: "fit" }] });
    const described = { s01: d("s01", "s01.norm"), s02: d("s02", "D.value") };
    const evaluated = { s01: { chi2: 697.5, norm: 1 }, s02: { chi2: 99.2, norm: 1.0234 } };
    const columns = lists.listColumns({ schema, kind: "dataset", rows, described,
      descriptionFields: kinds.description.fields, evaluated });
    const t = kinds.description.fields;
    assert.deepEqual(heads(columns), ["Name", t.tree.title, t.observable.title, t.data_file.title,
      `${t.kept_E_lab_MeV.title} (${t.kept_E_lab_MeV.unit})`, t.rows_kept.title, t.active.title,
      schema.records.dataset.fields.norm.title, "χ²"]);
    const c = cells(columns, rows[1]);
    assert.equal(c[`${t.kept_E_lab_MeV.title} (${t.kept_E_lab_MeV.unit})`], "0.2 – 1.4");
    assert.equal(c[schema.records.dataset.fields.norm.title], "1.0234 · fitted · shared D");
    assert.equal(cells(columns, rows[0])[schema.records.dataset.fields.norm.title], "1 · fitted");
    assert.equal(c[t.active.title], "yes");
  });

test("the data sets' list comes before the records they name", { skip }, () => {
  assert.deepEqual(lists.listKinds(schema, schema.sections.data), ["dataset", "normalization", "energy_scale"]);
  assert.deepEqual(lists.listKinds(schema, ["normalization", "dataset"]), ["dataset", "normalization"]);
});

test("the pairs' list: the particles by name, heads in words with units, rare columns only where set", { skip },
  () => {
    const pair = { key: 1, label: "p+15N", light_name: "p", heavy_name: "15N", light_jpi: "1/2+", heavy_jpi: "1/2-",
      light_mass: 1.007825, heavy_mass: 15.000109, light_excitation: 0, heavy_excitation: 0,
      separation_energy: 12.1274, radius: 4.85, kind_name: "particle", binding_energy: 0 };
    const f = kinds.pairs.fields;
    let columns = lists.pairColumns(f, [pair]);
    const h = heads(columns);
    assert.ok(!h.some((x) => /kind|E_x2|_/.test(x)), h);
    assert.deepEqual(h.slice(0, 4), [f.key.title, f.light_name.title, f.light_jpi.title,
      `${f.light_mass.title} (u)`]);
    assert.ok(h.includes(`${f.radius.title} (fm)`) && h.includes(`${f.separation_energy.title} (MeV)`));
    assert.ok(!h.includes(`${f.binding_energy.title} (MeV)`));
    columns = lists.pairColumns(f, [pair, { ...pair, key: 2, binding_energy: 2.224566 }]);
    assert.ok(heads(columns).includes(`${f.binding_energy.title} (MeV)`));
    assert.ok(columns.every((c) => c.doc.endsWith(`(${c.key})`)));
  });
