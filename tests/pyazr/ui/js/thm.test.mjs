// What the THM tab (tabs/thm.js) computes itself: whether a project has THM content, the schema of a form (the key
// table with the keys it lacks as text fields; the names offered from pyazr's lists), the members of each experiment,
// and the request of an experiment form (the keys changed and the data sets that join).  Everything else it shows is
// the server's.
import { test } from "node:test";
import assert from "node:assert/strict";
import { mod } from "./common.mjs";

const t = await mod("tabs/thm.js");
const forms = await mod("forms.js");

const KEYS = {
  "<thm>": { title: "<thm> options", doc: "", fields: { vertex: { type: "enum", values: ["perlevel", "constant"],
    doc: "" } } },
  experiment: { title: "experiment", doc: "", fields: { Ebeam: { type: "number", title: "Ebeam (MeV)", doc: "" } } },
};
const records = [
  { name: "s1", kind: "dataset", fields: { thm: "LC09" } }, { name: "s2", kind: "dataset", fields: { thm: "LC09" } },
  { name: "s3", kind: "dataset", fields: { thm: "global" } }, { name: "s4", kind: "dataset", fields: {} },
];

test("THM content: a <thm> section, or a data set with a thm field", () => {
  assert.equal(t.hasThm(null, records.slice(3)), false);
  assert.equal(t.hasThm(null, records), true);
  assert.equal(t.hasThm({ options: {}, experiments: {} }, []), true);
  assert.equal(t.hasThm(undefined, undefined), false);
});

test("a key the table lacks becomes a text field; the table's own fields and extra ones stay", () => {
  const schema = t.keySchema(KEYS, "<thm>", { vertex: "constant", "spectatorEnergy[2]": "0.4" });
  const fields = schema.records["<thm>"].fields;
  assert.deepEqual(Object.keys(fields), ["vertex", "spectatorEnergy[2]"]);
  assert.equal(fields["spectatorEnergy[2]"].type, "text");
  assert.equal(KEYS["<thm>"].fields["spectatorEnergy[2]"], undefined, "the table itself is not changed");
  const form = forms.describeForm(schema, "<thm>", { vertex: "constant", "spectatorEnergy[2]": "0.4" }, { schema });
  assert.deepEqual(form.lines.flat().map((e) => [e.key, e.texts]), [["vertex", ["constant"]],
    ["spectatorEnergy[2]", ["0.4"]]]);
  const exp = t.keySchema(KEYS, "experiment", {}, { [t.JOIN]: { type: "text", title: "join" } });
  assert.deepEqual(Object.keys(exp.records.experiment.fields), ["Ebeam", t.JOIN]);
});

test("members: the data sets with a thm field, grouped by the experiment they name", () => {
  const m = t.membersOf(records);
  assert.deepEqual(m.names, ["s1", "s2", "s3"]);
  assert.deepEqual(m.members, { LC09: ["s1", "s2"], global: ["s3"] });
});

test("an experiment form's edit: the keys changed (null removes one) and the data sets that join", () => {
  assert.deepEqual(t.experimentEdit({ Ebeam: "60", ps: null, [t.JOIN]: "s3  s4" }),
    { keys: { Ebeam: "60", ps: null }, datasets: ["s3", "s4"] });
  assert.deepEqual(t.experimentEdit({ Ebeam: "60" }), { keys: { Ebeam: "60" }, datasets: [] });
});

test("a nucleus and a global optical potential are offered from pyazr's lists (GET /api/kinds), by the grammar", () => {
  const choices = t.choicesOf({ nuclides: { d: {}, "15N": {} }, optical_potentials: { kd03: {}, ancai06: {} } });
  const ctx = { choices };
  const beam = forms.WIDGETS.text({ type: "text", grammar: "<nucleus>" }, "15N", ctx).controls[0];
  assert.deepEqual(beam.list, ["d", "15N"]);
  const optical = forms.WIDGETS.text({ type: "text", grammar: "plane | <global set>[:extrapolate] | <parameters>" },
    "", ctx).controls[0];
  assert.deepEqual(optical.list, ["kd03", "ancai06"]);
  assert.deepEqual(forms.WIDGETS.text({ type: "text", grammar: "all | thmin-thmax" }, "", ctx).controls[0].list, []);
  assert.deepEqual(t.choicesOf(null), { "<nucleus>": [], "<global set>": [] });
});
