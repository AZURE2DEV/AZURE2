// What the Nucleus tab (tabs/nucleus.js) makes itself for its pair and level edits: the fields of a pair form (the
// route's arguments, from the server's route table), the body of add_pair (the partners as dicts, a photon without a
// light partner, the other arguments where typed), and the body of add_level (the channels checked among the allowed
// ones, or "allowed" when all are).  Every number is the engine's to check.
import { test } from "node:test";
import assert from "node:assert/strict";
import { engineKinds, mod } from "./common.mjs";

const n = await mod("tabs/nucleus.js");
const kinds = engineKinds();
const skipNoKinds = kinds ? false : "pyazr's tables were not given (AZURE2_UI_KINDS)";

const ROUTES = [
  { method: "PATCH", path: "/api/project/pairs/{key}", arguments: ["key:int", "radius", "binding_energy", "new_one"] },
  { method: "POST", path: "/api/project/pairs", arguments: ["kind", "light", "heavy", "separation_energy", "radius",
    "binding_energy", "key"] },
];

test("a pair form's fields are the route's arguments: a field the route gains is a field of the form", () => {
  assert.deepEqual(n.routeArguments(ROUTES, "PATCH", "/api/project/pairs/{key}"),
    ["radius", "binding_energy", "new_one"]);
  assert.deepEqual(n.routeArguments(ROUTES, "POST", "/api/project/pairs"),
    ["kind", "light", "heavy", "separation_energy", "radius", "binding_energy", "key"]);
  assert.deepEqual(n.routeArguments(ROUTES, "GET", "/nowhere"), []);
});

test("add_pair's body: partners as dicts, numbers as numbers, other texts as typed, empty ones left out",
  { skip: skipNoKinds }, () => {
  const table = kinds.pairs;
  const fields = n.routeArguments(ROUTES, "POST", "/api/project/pairs");
  const texts = { kind: "particle", "light.J": "0.5", "light.parity": "1", "light.Z": "1", "light.mass": "1.007825",
    "heavy.J": "0", "heavy.parity": "1", "heavy.Z": "6", "heavy.mass": "12", separation_energy: "1.94",
    radius: "4.2", binding_energy: "", key: "abc" };
  assert.deepEqual(n.newPairRequest(texts, fields, table), { kind: "particle",
    light: { J: 0.5, parity: 1, Z: 1, mass: 1.007825 }, heavy: { J: 0, parity: 1, Z: 6, mass: 12 },
    separation_energy: 1.94, radius: 4.2, key: "abc" });
  const photon = n.newPairRequest({ ...texts, kind: "photon" }, fields, table);
  assert.equal(photon.light, undefined, "a photon pair has no light partner");
  assert.deepEqual(photon.heavy, { J: 0, parity: 1, Z: 6, mass: 12 });
});

test("the pair form's words and units are pyazr's table; the parity a choice of its values", { skip: skipNoKinds },
  () => {
    const table = kinds.pairs;
    assert.deepEqual(n.partnerFields(table).map(([k]) => k), ["J", "parity", "Z", "mass", "g"]);
    assert.deepEqual(table.partner_fields.parity.values, [1, -1]);
    assert.equal(table.partner_fields.mass.unit, "u");
    assert.deepEqual(n.pairWords(table, "separation_energy"), { title: table.fields.separation_energy.title,
      unit: "MeV", doc: table.fields.separation_energy.doc });
    assert.equal(n.pairWords(table, "radius").unit, "fm");
    assert.deepEqual(n.pairWords(table, "a_new_argument"), { title: "a new argument", unit: "", doc: "" });
    assert.deepEqual(table.kinds, ["particle", "photon", "beta"]);
  });

test("add_level's body: the channels checked, width 0, fixed; every one checked: \"allowed\"", () => {
  const allowed = [{ pair: 1, l: 0, s: 0.5 }, { pair: 1, l: 2, s: 0.5 }, { pair: 2, l: 1, s: 0 }];
  const level = { J: 0.5, parity: "-1", energy: "8.5", fixed: true };
  assert.deepEqual(n.newLevelRequest(level, allowed, [true, false, true]), { J: 0.5, parity: -1, energy: 8.5,
    fixed: true, channels: [{ pair: 1, l: 0, s: 0.5, width: 0, fixed: true }, { pair: 2, l: 1, s: 0, width: 0,
      fixed: true }] });
  assert.equal(n.newLevelRequest(level, allowed, [true, true, true]).channels, "allowed");
  assert.deepEqual(n.newLevelRequest(level, [], []).channels, []);
});
