// api.js: an edit carries the revision the page last saw; a newer revision in an answer is reported; a 412 makes
// the page re-read (onStale) and takes the server's revision; a refusal's body becomes what / allowed / fix and the
// engine's diagnostics; failed requests are listed (the self-test reports them), expected ones (quiet) are not.
import { test } from "node:test";
import assert from "node:assert/strict";
import { mod } from "./common.mjs";

const api = await mod("api.js");

function fakeFetch(answers) {
  const sent = [];
  const fetchFn = async (url, init) => {
    sent.push({ url, init, body: init.body ? JSON.parse(init.body) : undefined });
    const [status, body] = answers.shift();
    return { ok: status < 400, status, text: async () => (typeof body === "string" ? body : JSON.stringify(body)) };
  };
  return { fetchFn, sent };
}

test("an edit carries the revision; an answer's newer revision is reported", async () => {
  const seen = [];
  const { fetchFn, sent } = fakeFetch([[200, { rev: 4 }], [200, { rev: 5 }], [200, { objective: 1 }]]);
  const c = api.createApi({ fetchFn, onRevision: (r) => seen.push(r) });
  await c.get("/api/project");
  await c.patch("/api/project/records/D", { fields: { value: "1 fixed" } });
  await c.post("/api/evaluate", {});
  assert.equal(sent[1].body.rev, 4);
  assert.equal(sent[1].init.method, "PATCH");
  assert.equal(sent[2].body.rev, undefined, "an evaluation is not an edit");
  assert.deepEqual(seen, [4, 5]);
  assert.equal(c.rev, 5);
});

test("a stale revision (412): onStale, the server's revision taken, the error thrown", async () => {
  let stale = null;
  const { fetchFn } = fakeFetch([[412, { error: "StaleRevision", message: "made on 3, at 7", rev: 7 }]]);
  const c = api.createApi({ fetchFn, onStale: (v) => { stale = v; } });
  c.rev = 3;
  await assert.rejects(c.post("/api/project/undo", {}), (e) => e.status === 412);
  assert.equal(stale.rev, 7);
  assert.equal(c.rev, 7);
});

test("a refusal: what, allowed, fix and the engine's diagnostics", async () => {
  const body = { error: "FormatError", message: "refused", what: "the edit is refused",
    diagnostics: [{ code: "E_UNKNOWN_REFERENCE", severity: "error", record: "artemov", field: "norm",
      text: "no normalization record is named 'abc'", allowed: "D, E", fix: "name a record" }] };
  const v = api.errorView(422, body);
  assert.equal(v.what, "the edit is refused");
  assert.equal(v.diagnostics[0].field, "norm");
  assert.equal(api.errorView(422, { ...body, allowed: "D, E", fix: "name a record" }).allowed, "D, E");
  assert.equal(api.errorView(0, null).title, "Network error");
  const { fetchFn } = fakeFetch([[422, body], [404, { message: "no such" }]]);
  const c = api.createApi({ fetchFn });
  await assert.rejects(c.patch("/api/project/records/x", {}), (e) => e.view.diagnostics.length === 1);
  await assert.rejects(c.get("/api/project/describe/x", { quiet: true }));
  assert.equal(c.failures.length, 1, "a quiet refusal is not listed among the failures");
});

test("paths: segments encoded (a level label holds a blank and may hold /)", () => {
  assert.equal(api.pathOf("/api/project/levels", "1- 12.4258"), "/api/project/levels/1-%2012.4258");
  assert.equal(api.isEdit("POST", "/api/project/undo"), true);
  assert.equal(api.isEdit("POST", "/api/jobs"), false);
  assert.equal(api.isEdit("GET", "/api/project"), false);
});
