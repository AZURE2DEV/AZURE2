// tables.js: sorting (numbers as numbers, natural order of names, empty cells last, stable), filtering (words,
// column:value), number formatting (every digit kept unless digits are asked for) and column heads with units.
import { test } from "node:test";
import assert from "node:assert/strict";
import { mod } from "./common.mjs";

const t = await mod("tables.js");

test("sorting: numbers numerically, names naturally, empty last, stable", () => {
  const rows = [{ n: "s10", v: 2 }, { n: "s2", v: 10 }, { n: "s1", v: null }, { n: "s3", v: 2 }];
  assert.deepEqual(t.sortRows(rows, "v").map((r) => r.n), ["s10", "s3", "s2", "s1"]);
  assert.deepEqual(t.sortRows(rows, "v", true).map((r) => r.n), ["s2", "s10", "s3", "s1"], "empty last, descending");
  assert.deepEqual(t.sortRows(rows, "n").map((r) => r.n), ["s1", "s2", "s3", "s10"]);
  assert.deepEqual(rows.map((r) => r.n), ["s10", "s2", "s1", "s3"], "the input is left as it is");
  assert.ok(t.compareValues("inf", 1e300) > 0 && t.compareValues("-inf", -1e300) < 0);
});

test("filtering: every word, case-insensitive; column:value names a column", () => {
  const cols = [{ key: "name", title: "name" }, { key: "quantity", title: "quantity" }];
  const rows = [{ name: "s01_int_thin", quantity: "sigma" }, { name: "s02_diff090", quantity: "dsigma_domega" }];
  assert.deepEqual(t.filterRows(rows, "SIGMA", cols).length, 2);
  assert.deepEqual(t.filterRows(rows, "name:diff", cols).map((r) => r.name), ["s02_diff090"]);
  assert.deepEqual(t.filterRows(rows, "quantity:dsig s02", cols).map((r) => r.name), ["s02_diff090"]);
  assert.deepEqual(t.filterRows(rows, "s03", cols), []);
  assert.deepEqual(t.filterRows(rows, "  ", cols).length, 2);
});

test("numbers: every digit kept unless digits are asked for; the server's inf and nan pass", () => {
  assert.equal(t.formatNumber(10905.2313020616), "10905.2313020616");
  assert.equal(t.formatNumber(10905.2313020616, 6), "10905.2");
  assert.equal(t.formatNumber(1.2e-7, 3), "1.20e-7");
  assert.equal(t.formatNumber("inf"), "inf");
  assert.equal(t.formatNumber(null), "");
  assert.equal(t.cellText({ x: 0.1 + 0.2 }, { key: "x" }), String(0.1 + 0.2));
});

test("a column head carries its unit", () => {
  assert.equal(t.columnHead({ title: "E_λ", unit: "MeV" }), "E_λ (MeV)");
  assert.equal(t.columnHead({ title: "name" }), "name");
});
