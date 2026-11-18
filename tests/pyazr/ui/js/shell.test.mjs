// The frame around the tabs (pyazr/ui/static/shell.js, docs/upgrade/web_navigation.md): what the status line's job
// mark and validator mark say, where a validator message leads, which page of the documentation a tab opens, when
// the first visit's overview opens.
import { test } from "node:test";
import assert from "node:assert/strict";
import { mod } from "./common.mjs";

const shell = await mod("shell.js");

test("the job mark shows the running job, else the newest; its words carry the progress", () => {
  assert.equal(shell.jobToMark([]), null);
  assert.equal(shell.jobMarkView(null), null);
  const jobs = [
    { id: "1", kind: "fit", method: "engine_migrad", state: "converged", started: 10 },
    { id: "2", kind: "fit", method: "engine_lm", state: "running", started: 5,
      progress: { calls: 152, objective: 481.2, best: 480.25 } },
    { id: "3", kind: "bands", method: "linear", state: "finished", started: 20 },
  ];
  assert.equal(shell.jobToMark(jobs).id, "2", "a running job before a newer finished one");
  const v = shell.jobMarkView(shell.jobToMark(jobs));
  assert.ok(v.live && v.canStop, v);
  assert.match(v.text, /^job 2 fit engine_lm: .* · 152 calls · best 480\.25/);
  assert.equal(v.hash, "#/fit/job/2");
  const done = jobs.filter((j) => j.id !== "2");
  assert.equal(shell.jobToMark(done).id, "3", "without a running job, the newest");
  const d = shell.jobMarkView(shell.jobToMark(done));
  assert.ok(!d.live && !d.canStop && !d.failed, d);
  assert.ok(shell.jobMarkView({ id: "4", kind: "fit", method: "x", state: "failed" }).failed);
});

test("the validator mark counts errors and warnings in words", () => {
  assert.deepEqual(shell.validationView([]), { text: "valid", kind: "ok", errors: 0, warnings: 0 });
  const two = shell.validationView([{ severity: "error" }, { code: "E_X" }, { severity: "warning" }]);
  assert.deepEqual(two, { text: "2 errors, 1 warning", kind: "error", errors: 2, warnings: 1 });
  assert.equal(shell.validationView([{ severity: "warning" }]).kind, "warn");
});

test("a validator message leads to the row it is about", () => {
  const schema = { records: { dataset: { sections: ["data"] }, gamma_decay: { sections: ["gammaDecays"] } } };
  const records = [{ name: "s01", kind: "dataset" }, { name: "g1", kind: "gamma_decay" }];
  const levels = [{ label: "1/2+ 2.365" }];
  assert.equal(shell.diagnosticRoute({ record: "s01" }, records, schema, levels), "#/data/s01");
  assert.equal(shell.diagnosticRoute({ record: "g1" }, records, schema, levels), "#/nucleus/g1");
  assert.equal(shell.diagnosticRoute({ record: "1/2+ 2.365" }, records, schema, levels),
    "#/nucleus/level/1%2F2%2B%202.365");
  assert.equal(shell.diagnosticRoute({ record: "" }, records, schema, levels), "#/project");
  assert.equal(shell.diagnosticRoute({ record: "gone" }, records, schema, levels), "#/project");
});

test("the documentation of a tab is the user guide's section of its title", () => {
  const base = "https://example.org/docs/";
  assert.equal(shell.docsOf(base, { title: "Levels and Channels" }),
    `${base}user_guide/web_interface.html#levels-and-channels`);
  assert.equal(shell.docsOf(base, { title: "THM", group: "Physics and THM" }),
    `${base}user_guide/web_interface.html#physics-and-thm`);
  assert.equal(shell.docsOf(base, null), `${base}user_guide/web_interface.html#the-tabs`);
});

test("the overview opens on a first visit without a route, never on a link to a place", () => {
  const store = (v) => ({ getItem: () => v });
  assert.equal(shell.overviewDue("", store(null)), true);
  assert.equal(shell.overviewDue("#", store(null)), true);
  assert.equal(shell.overviewDue("", store("1")), false);
  assert.equal(shell.overviewDue("#/data/s01", store(null)), false);
  assert.equal(shell.overviewDue("", { getItem: () => { throw new Error("blocked"); } }), false);
});
