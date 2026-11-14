// The run's settings (tabs/settings.js): a row per entry of the engine's <config> table (pyazr.kinds()["config"]),
// in its order, with the value in use, whether the project sets it, its group, the rare rows folded unless set, and
// the search.
import { test } from "node:test";
import assert from "node:assert/strict";
import { mod } from "./common.mjs";

const s = await mod("tabs/settings.js");

const FIELDS = {
  amatrix: { type: "enum", values: ["true", "false"], positional: true, group: "Runtime options", title: "Formalism",
    doc: "true: A-matrix", rare: false },
  check_data: { type: "enum", values: ["none", "screen", "file"], positional: true, group: "Checks", title: "Data",
    doc: "where", rare: true },
  use_brune: { type: "enum", values: ["0", "1"], default: "1", group: "Runtime options", title: "Use Brune formalism",
    doc: "Brune's level matrix", rare: false, option: true },
  maxOrbitalMomentum: { type: "integer", min: 0, group: "Physics", title: "Maximum orbital momentum",
    doc: "the largest l", absent: "the largest l of any channel", rare: false },
  use_ec_table: { type: "enum", values: ["0", "1"], default: "1", group: "Numerics", title: "EC tables",
    doc: "tables", rare: true, option: true },
};

test("a row per entry, in the table's order, with the value in use and whether the project sets it", () => {
  const rows = s.settingRows(FIELDS, { amatrix: "true", check_data: "screen", use_ec_table: "0" });
  assert.deepEqual(rows.map((r) => r.name), Object.keys(FIELDS));
  const by = Object.fromEntries(rows.map((r) => [r.name, r]));
  assert.deepEqual([by.amatrix.set, by.check_data.set, by.use_brune.set, by.use_ec_table.set], [false, true, false,
    true], "a positional entry is set when it is not its first value; a key when present");
  assert.equal(by.use_brune.inUse, "1", "an option not set: its default");
  assert.equal(by.maxOrbitalMomentum.inUse, null, "a key without a default: nothing in use but the engine's rule");
  assert.deepEqual([by.check_data.rare, by.use_ec_table.rare], [false, false], "a rare row the project sets is shown");
  assert.equal(s.settingRows(FIELDS, {}).find((r) => r.name === "use_ec_table").rare, true);
  assert.deepEqual(s.groupsOf(rows), ["Runtime options", "Checks", "Physics", "Numerics"]);
});

test("the search matches the key, the title, the meaning and the group", () => {
  const shown = (q) => s.settingRows(FIELDS, {}, q).filter((r) => r.shown).map((r) => r.name);
  assert.deepEqual(shown("brune"), ["use_brune"]);
  assert.deepEqual(shown("orbital largest"), ["maxOrbitalMomentum"]);
  assert.deepEqual(shown("numerics"), ["use_ec_table"]);
  assert.equal(shown("").length, 5);
});

test("a switch is a key of 0 or 1; its value as the engine reads it", () => {
  assert.ok(s.isSwitch(FIELDS.use_brune) && !s.isSwitch(FIELDS.amatrix) && !s.isSwitch(FIELDS.check_data));
  assert.deepEqual([s.switchValue(true), s.switchValue(false)], ["1", "0"]);
});
