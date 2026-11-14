// fitting.js and the Fit tab's tables: the parameter rows from the server's answers (working x for free parameters,
// a job's best point, errors, start and MINOS columns only when a job carries them); the job form from the server's
// job table (rows that wait for pyazr disabled with their reason, a widget per setting type); the progress reducer;
// the result's followers (MINOS after an iminuit fit); the objective-against-call plot row (log y); matrices as CSV.
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { mod } from "./common.mjs";

const f = await mod("fitting.js");
const plots = await mod("plots.js");
const fitjob = await mod("tabs/fitjob.js");
const TABLE = process.env.AZURE2_UI_JOBS ? JSON.parse(readFileSync(process.env.AZURE2_UI_JOBS, "utf8")) : null;

// rows shaped as GET /api/parameters answers them (T2)
const PARAMS = [
  { name: "energy_1", kind: "k_a", unit: "MeV", x: 12.4258, physical_value: 12.4258, free: true, bounds: null },
  { name: "width_1_1", kind: "k_b", unit: "MeV^1/2", x: null, physical_value: 547.9, free: false, bounds: null },
  { name: "tgt.areal", kind: "k_c", unit: "atoms/cm^2", x: 2e17, physical_value: 2e17, free: true, bounds: [0, null],
    prior_centre: 2e17, prior_sigma: 2e16 },
];

test("parameter rows: working x for free ones, nothing for fixed ones; job columns only with a job", () => {
  const rows = f.parameterRows(PARAMS, { x: { energy_1: 12.5, "tgt.areal": 2e17 } });
  assert.equal(rows[0].x, 12.5, "the working x, not the row's x0");
  assert.equal(rows[1].x, null);
  assert.equal(rows[1].physical, 547.9);
  assert.equal(rows[2].prior, "200000000000000000 ± 20000000000000000");
  assert.equal(rows[2].bounds, "[0, ]");
  const keys = f.shownColumns(rows).map((c) => c.key);
  assert.ok(!keys.includes("fit") && !keys.includes("error") && !keys.includes("start"), keys.join());
  const job = { best: { x: { energy_1: 12.43, "tgt.areal": 2.1e17 } },
    result: { errors: { energy_1: 0.001 }, start_x: { energy_1: 12.4258 }, minos: { energy_1: { lower: -0.002,
      upper: 0.003 } } } };
  const withJob = f.parameterRows(PARAMS, { x: {}, job });
  assert.equal(withJob[0].fit, 12.43);
  assert.equal(withJob[0].error, 0.001);
  assert.equal(withJob[0].minos_upper, 0.003);
  assert.equal(withJob[1].fit, null);
  const shown = f.shownColumns(withJob).map((c) => c.key);
  assert.ok(shown.includes("fit") && shown.includes("minos_lower") && !shown.includes("start"));
  assert.ok(f.shownColumns(withJob, { compare: true }).map((c) => c.key).includes("start"));
  assert.ok(f.PARAMETER_COLUMNS.every((c) => c.doc), "every column says where its words come from");
});

test("kinds in the order of the engine's table, with their docs; the filters", () => {
  const kinds = f.kindsOf(PARAMS, { k_c: { doc: "density" }, k_a: { doc: "energy" } });
  assert.deepEqual(kinds.map((k) => k.kind), ["k_c", "k_a", "k_b"]);
  assert.equal(kinds[0].doc, "density");
  const rows = f.parameterRows(PARAMS);
  assert.equal(f.filterParameters(rows, { freeOnly: true }).length, 2);
  assert.deepEqual(f.filterParameters(rows, { kind: "k_b" }).map((r) => r.name), ["width_1_1"]);
});

test("the job form from the server's job table: waiting rows disabled with their reason", { skip: !TABLE }, () => {
  const choices = f.jobChoices(TABLE);
  assert.deepEqual(choices.map((k) => k.kind), [...new Set(TABLE.map((r) => r.kind))]);
  for (const row of TABLE) {
    const m = choices.find((k) => k.kind === row.kind).methods.find((x) => x.method === row.method);
    assert.equal(m.available, row.available);
    if (!row.available) assert.ok(m.reason.includes(row.missing || row.needs), m.reason);
  }
  const fit = choices.find((k) => k.kind === "fit");
  assert.equal(f.defaultMethod(fit, "no_such").method, TABLE.find((r) => r.kind === "fit" && r.available).method);
  const waiting = TABLE.find((r) => !r.available);  // none where every package is installed and no call waits
  if (waiting) assert.equal(f.defaultMethod(fit, waiting.method)?.method === waiting.method, false);
  const types = new Set(TABLE.flatMap((r) => r.settings.map((s) => s.type)));
  for (const t of types) assert.ok(fitjob.SETTING_WIDGETS[t], `no widget for the setting type ${t}`);
  const m = fit.methods.find((x) => x.settings.length);
  assert.deepEqual(f.settingsOf(m, { [m.settings[0].name]: "" }), {}, "an empty box is the default");
  assert.deepEqual(f.settingsOf(m, { [m.settings[0].name]: "12" }), { [m.settings[0].name]: "12" });
});

test("start points: the working x, the text, a job's best point when it has one", () => {
  assert.deepEqual(f.startPoints(null).map((p) => p.x0), [undefined, "text"]);
  assert.deepEqual(f.startPoints({ id: "ab", best: { x: { a: 1 } } })[2].x0, { job: "ab" });
});

test("the progress reducer: numbers as sent, Stop and Kill by state", () => {
  const rec = { id: "j1", kind: "fit", method: "m", state: "running", build_s: 0.4,
    progress: { calls: 12, objective: 700.5, best: 690.1, elapsed: 3.2 }, trace: [[1, 900, 900], [12, 700.5, 690.1]],
    log: ["a", "b"] };
  const v = f.progressView(rec);
  assert.equal(v.calls, 12);
  assert.equal(v.best, 690.1);
  assert.ok(v.canStop && v.canKill && !v.final);
  assert.equal(v.trace.length, 2);
  const stopping = f.progressView({ ...rec, state: "stopping" });
  assert.ok(!stopping.canStop && stopping.canKill);
  const done = f.progressView({ ...rec, state: "stopped" });
  assert.ok(done.final && !done.canStop && !done.canKill);
  assert.equal(f.progressView({ id: "j", state: "building" }).calls, 0);
});

test("the result: the recipe's words, MINOS offered after an iminuit fit", { skip: !TABLE }, () => {
  const rec = { id: "j1", kind: "fit", method: "iminuit_migrad", state: "converged", best: { objective: 5, x: {} },
    result: { status: "converged", message: "MIGRAD: valid minimum", calls: 40, elapsed: 1.5, start_objective: 9,
      errors: {} } };
  const v = f.resultView(rec, TABLE);
  assert.equal(v.message, "MIGRAD: valid minimum");
  const minos = (view) => view.followers.filter((r) => r.kind === "minos").map((r) => `${r.kind}/${r.method}`);
  assert.deepEqual(minos(v), ["minos/iminuit"]);
  assert.equal(minos(f.resultView({ ...rec, method: "scipy_lbfgsb" }, TABLE)).length, 0);
});

test("the objective against call: its own plot row, log y, numbers as sent", () => {
  const rec = { id: "j1", kind: "fit", method: "m", trace: [[1, 900, 900], [2, 950, 900], [3, 700, 700]] };
  const plot = plots.plotOf(rec, {});
  assert.equal(plot.row, "objective_against_call");
  assert.equal(plot.y.log, true);
  assert.deepEqual(plot.traces[0].y, [900, 950, 700]);
  assert.deepEqual(plot.traces[1].y, [900, 900, 700]);
  assert.equal(plots.plotlyFigure(plot).layout.yaxis.type, "log");
  assert.equal(plots.csvOf(plot).split("\n")[0], "call,objective at each call,best objective so far");
});

test("matrices: rows by name, CSV with the numbers as sent; the job a visit shows", () => {
  const m = { names: ["a", "b"], matrix: [[1, -0.25], [-0.25, 1]] };
  assert.deepEqual(f.matrixRows(m)[1], { name: "b", a: -0.25, b: 1 });
  assert.equal(f.matrixCsv(m), ",a,b\na,1,-0.25\nb,-0.25,1\n");
  const jobs = [{ id: "old", state: "stopped", started: 1 }, { id: "new", state: "converged", started: 2 }];
  assert.equal(f.jobToShow(jobs).id, "new");
  assert.equal(f.jobToShow([...jobs, { id: "run", state: "running", started: 0 }]).id, "run");
  assert.equal(f.jobToShow(jobs, "old").id, "old");
});

test("the engine's fitters: facts as words, the engine's MINOS offered after its fit, errors said", { skip: !TABLE },
  () => {
    const engineFits = TABLE.filter((r) => r.kind === "fit" && r.method.startsWith("engine_"));
    assert.ok(engineFits.length >= 1 && engineFits.every((r) => r.available), "the engine's rows run");
    assert.equal(f.jobChoices(TABLE).find((k) => k.kind === "fit").methods[0].method, engineFits[0].method,
      "the default fit is the engine's");
    for (const row of engineFits) {
      const words = f.factsText(row.facts);
      assert.ok(words.includes(`covariance: ${row.facts.covariance ? "yes" : "no"}`), words);
      const rec = { id: "e", kind: "fit", method: row.method, state: "converged", best: { objective: 1, x: {} },
        result: { status: "converged", message: "the engine's: converged", errors: { a: 1 },
          errors_from: "HESSE after the engine's MIGRAD" } };
      const v = f.resultView(rec, TABLE);
      assert.equal(v.errorsFrom, "HESSE after the engine's MIGRAD");
      assert.deepEqual(v.followers.filter((r) => r.kind === "minos").map((r) => `${r.kind}/${r.method}`),
        row.facts["MINOS can follow"] ? [`minos/${row.method}`] : []);
    }
  });

test("the progress reducer: numbers beyond the four common ones shown by name", () => {
  const v = f.progressView({ id: "j", state: "running", progress: { calls: 3, objective: 2, best: 1, elapsed: 0.5,
    stage: "hesse", measure: null } });
  assert.deepEqual(v.more, [["stage", "hesse"], ["measure", null]]);
});

test("followers of any method of a kind, offered only with what they take; the sampler's summary rows",
  { skip: !TABLE }, () => {
    const fit = { id: "f", kind: "fit", method: "scipy_lbfgsb", state: "converged", best: { objective: 1, x: {} },
      result: { status: "converged" } };
    const names = (rec) => f.resultView(rec, TABLE).followers.map((r) => `${r.kind}/${r.method}`);
    assert.ok(names(fit).includes("sample/emcee"), "the sampler follows any fit");
    assert.ok(!names(fit).includes("bands/covariance"), "no band without a covariance");
    const withCov = { ...fit, result: { ...fit.result, covariance: { names: ["a"], matrix: [[1]] } } };
    assert.ok(names(withCov).includes("bands/covariance"));
    const sample = { id: "s", kind: "sample", method: "emcee", state: "finished", best: { objective: 1, x: {} },
      result: { status: "finished", chain: "output/chain_s.npy" } };
    assert.deepEqual(names(sample).filter((n) => n.startsWith("bands")), ["bands/draws"]);
    assert.equal(f.follows({ follows: ["fit", null] }, fit), true);
    // the sampler follows a fit and a sampler (its chain continues), whatever the sampler's state
    assert.equal(f.follows({ follows: [["fit", null], ["sample", "emcee"]] }, { kind: "sample", method: "emcee" }),
      true);
    const killed = { id: "k", kind: "sample", method: "emcee", state: "killed", best: { objective: 1, x: {} } };
    assert.ok(names(killed).includes("sample/emcee"), "a killed sampler is offered its continuation");
    assert.ok(!names(killed).includes("bands/draws"), "no band from a chain the job's result lacks");
    assert.equal(f.variedText({ varied: ["a"], free: ["a", "n"] }),
      "a (1 of 2 free parameters; held at the start point: n)");
    assert.equal(f.variedText({}), "");
    const lm = f.resultView({ ...fit, result: { status: "converged", notices: ["falling back to MIGRAD."],
      fitter: { asked: "lm", ran: "migrad" } } }, TABLE);
    assert.deepEqual([lm.notices, lm.fitter.ran], [["falling back to MIGRAD."], "migrad"]);
    assert.equal(f.follows({ follows: ["fit", "engine_lm"] }, fit), false);
    assert.equal(f.STATES.finished.final, true);
    const rows = f.chainRows({ quantiles: { a: { "2.5": 1, "16": 2, "50": 3, "84": 4, "97.5": 5 } },
      mean: { a: 3 }, errors: { a: 1 }, tau: { a: 20 }, effective_draws: { a: 100 }, rhat: { a: 1.01 } });
    assert.deepEqual(rows[0], { name: "a", mean: 3, sd: 1, "q2.5": 1, q16: 2, q50: 3, q84: 4, "q97.5": 5, tau: 20,
      n_eff: 100, rhat: 1.01 });
  });

test("MINOS validity: per side the flags the engine set, one word for iminuit's; the outcome and invalid sides", () => {
  const engine = { lower: -0.1, upper: 0.2, valid: false, sides: {
    lower: { valid: true, at_limit: false, call_limit: false, new_minimum: false },
    upper: { valid: false, at_limit: true, call_limit: false, new_minimum: false } } };
  assert.equal(f.minosValidity(engine), "− valid; + at_limit");
  const both = { ...engine, sides: { lower: engine.sides.lower, upper: engine.sides.lower } };
  assert.equal(f.minosValidity(both), "valid");
  assert.equal(f.minosValidity({ lower: -1, upper: 1, valid: false }), "not valid");
  assert.equal(f.minosValidity(null), null);
  const rows = f.parameterRows([{ name: "a", free: true, x: 1 }], { job: { result: { minos: { a: engine } } } });
  assert.equal(rows[0].minos_valid, "− valid; + at_limit");
  assert.ok(f.shownColumns(rows).some((c) => c.key === "minos_valid"));
  const v = f.resultView({ state: "ended", result: { outcome: "call_limit", minos: { a: engine } } });
  assert.equal(v.outcome, "call_limit");
  assert.deepEqual(v.minosInvalid, ["a"]);
});
