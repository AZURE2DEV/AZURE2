// What the Fit tab shows, without a document (docs/upgrade/web_ui.md sections 4.4 and 6.4): the rows of the
// parameter table from the server's answers, the job form from the server's job table, the view of a running job
// (the progress reducer) and of its result.  Every number is the server's, as sent; the page computes nothing.
// node's test runner loads this module (tests/pyazr/ui/js/fitting.test.mjs); tabs/fit.js and tabs/fitjob.js draw it.

/** The states of a job (pyazr/ui/jobs.py) in the words the page shows; `final`: the job has ended. */
export const STATES = {
  queued: { words: "waiting for a thread", final: false },
  building: { words: "building its session", final: false },
  running: { words: "running", final: false },
  stopping: { words: "stopping: it ends before its next objective call, with the best point", final: false },
  stopped: { words: "stopped on request: the best point so far is kept", final: true },
  converged: { words: "converged", final: true },
  ended: { words: "ended without convergence", final: true },
  finished: { words: "finished: it ran all its steps", final: true },
  failed: { words: "failed", final: true },
  killed: { words: "killed: the best point reported before the kill is kept", final: true },
  interrupted: { words: "interrupted: the server ended while it ran; its checkpoint is kept", final: true },
};

export const isFinal = (state) => Boolean(STATES[state]?.final);

/**
 * The columns of the parameter table and where their words come from (each column's `doc`, shown on its head).
 * `edit` marks the column whose cells of free parameters are editable (a value of the working x).
 */
export const PARAMETER_COLUMNS = [
  { key: "name", title: "name", doc: "the engine's name of the parameter (GET /api/parameters)" },
  { key: "kind", title: "kind", doc: "the engine's parameter kind; its description on hover (GET /api/kinds)" },
  { key: "free", title: "free", format: (v) => (v ? "free" : "fixed"),
    doc: "free: a value of the working x, edited here; fixed: a value of the project text, edited in its record "
      + "or level" },
  { key: "unit", title: "x unit", doc: "the unit of x in the fit basis, as the engine's row gives it" },
  { key: "x", title: "x", numeric: true, edit: true,
    doc: "the working x (fit basis) of a free parameter (GET /api/x); Enter sends a new value (PATCH /api/x)" },
  { key: "start", title: "start", numeric: true, compare: true,
    doc: "the job's start point (its result's start_x), shown with Compare with start" },
  { key: "fit", title: "fit", numeric: true, digits: 10,
    doc: "the best point of the job shown (GET /api/jobs/{id}, best.x); every digit sent on hover" },
  { key: "error", title: "error", numeric: true, digits: 4,
    doc: "the error of the job shown, in the fit basis (its result says where it comes from: HESSE, the engine's "
      + "(JᵀJ + priors)⁻¹)" },
  { key: "minos_lower", title: "MINOS −", numeric: true, digits: 4,
    doc: "the MINOS lower interval of the job shown (the engine's or iminuit's), in the fit basis" },
  { key: "minos_upper", title: "MINOS +", numeric: true, digits: 4,
    doc: "the MINOS upper interval of the job shown (the engine's or iminuit's), in the fit basis" },
  { key: "minos_valid", title: "MINOS validity", wrap: true,
    doc: "whether MINOS found each side's crossing: valid, or the flags the fitter set (at_limit: the scan reached "
      + "the parameter's limit; call_limit: it ran out of calls; new_minimum: it found a lower point, the fit was "
      + "not at its minimum), per side where the fitter gives them (the engine's), else one word (iminuit's)" },
  { key: "physical", title: "physical value", numeric: true, digits: 8,
    doc: "the session's physical value of the parameter at the text's start values (GET /api/parameters, "
      + "physical_value); its unit is not sent: a width is in the project's convention" },
  { key: "bounds", title: "bounds", doc: "the engine's bounds of the parameter (none: unbounded)" },
  { key: "prior", title: "prior", doc: "centre ± width of the parameter's own prior as the engine's row gives it; "
    + "the priors of normalizations and shifts are in their records (Data tab)" },
];

/** The kinds of the parameters present, in the order of the engine's table (GET /api/kinds), each with its doc. */
export function kindsOf(params, kindsTable = {}) {
  const present = new Set(params.map((p) => p.kind));
  const ordered = Object.keys(kindsTable).filter((k) => present.has(k));
  const rest = [...present].filter((k) => !ordered.includes(k)).sort();
  return [...ordered, ...rest].map((k) => ({ kind: k, doc: kindsTable[k]?.doc || "",
    n: params.filter((p) => p.kind === k).length, free: params.filter((p) => p.kind === k && p.free).length }));
}

const priorOf = (p) => (p.prior_centre !== undefined && p.prior_centre !== null
  ? `${p.prior_centre} ± ${p.prior_sigma}` : "");
const boundsOf = (b) => (Array.isArray(b) ? `[${b.map((v) => (v === null ? "" : v)).join(", ")}]` : "");

/**
 * The rows of the parameter table: the server's parameter rows (GET /api/parameters) with the working x (GET /api/x)
 * and, when a job is shown, its best point, start point, errors and MINOS intervals (GET /api/jobs/{id}).
 */
export function parameterRows(params, { x = {}, job = null } = {}) {
  const best = job?.best?.x || {};
  const result = job?.result || {};
  const errors = result.errors || {};
  const minos = result.minos || {};
  const start = result.start_x || {};
  return params.map((p) => ({
    name: p.name, kind: p.kind, unit: p.unit, free: p.free,
    x: p.free ? (x[p.name] ?? p.x) : null, physical: p.physical_value ?? null,
    bounds: boundsOf(p.bounds), prior: priorOf(p),
    start: start[p.name] ?? null, fit: best[p.name] ?? null, error: errors[p.name] ?? null,
    minos_lower: minos[p.name]?.lower ?? null, minos_upper: minos[p.name]?.upper ?? null,
    minos_valid: minosValidity(minos[p.name]),
  }));
}

/** The words of a MINOS interval's validity: per side the flags set (the engine's {valid, at_limit, call_limit,
 * new_minimum}), "valid" where only valid is set; one word where the fitter gives no sides (iminuit). */
export function minosValidity(m) {
  if (!m) return null;
  if (m.sides) {
    const side = (flags) => {
      const set = Object.entries(flags || {}).filter(([k, v]) => v && k !== "valid").map(([k]) => k);
      return set.length ? set.join(", ") : (flags?.valid ? "valid" : "not valid");
    };
    const lower = side(m.sides.lower);
    const upper = side(m.sides.upper);
    return lower === upper ? lower : `− ${lower}; + ${upper}`;
  }
  return m.valid === true ? "valid" : (m.valid === false ? "not valid" : null);
}

/** The rows a filter keeps: one kind ("" for all) and, with `freeOnly`, the free parameters only. */
export function filterParameters(rows, { kind = "", freeOnly = false } = {}) {
  return rows.filter((r) => (!kind || r.kind === kind) && (!freeOnly || r.free));
}

/** The columns to show: those with numbers in some row (a job's columns appear with a job), start when comparing. */
export function shownColumns(rows, { compare = false } = {}) {
  const jobColumns = new Set(["start", "fit", "error", "minos_lower", "minos_upper", "minos_valid"]);
  return PARAMETER_COLUMNS.filter((c) => {
    if (c.compare && !compare) return false;
    if (!jobColumns.has(c.key)) return true;
    return rows.some((r) => r[c.key] !== null && r[c.key] !== undefined);
  });
}

/**
 * The job form's choices from the server's job table (GET /api/jobs/kinds): the kinds in the table's order, each with
 * its methods; a method that cannot run says why (the pyazr call it waits for, or the package it needs).  Only the
 * rows of the tab `tab` (a row's `tab`: the tab whose form starts it; the rate and the ENDF jobs have their own).
 */
export function jobChoices(table, tab = "fit") {
  const kinds = [];
  for (const row of table) {
    if ((row.tab || "fit") !== tab) continue;
    let k = kinds.find((x) => x.kind === row.kind);
    if (!k) { k = { kind: row.kind, methods: [], firstView: Boolean(row.first_view) }; kinds.push(k); }
    const reason = row.available ? "" : (row.missing ? `waits for pyazr: ${row.missing}`
      : `needs the package ${row.needs}, not installed`);
    k.methods.push({ method: row.method, doc: row.doc, available: row.available, reason,
      settings: row.settings || [], follows: row.follows || null, takes: row.takes || null,
      facts: row.facts || {} });
  }
  kinds.forEach((k) => { k.available = k.methods.some((m) => m.available); });
  return kinds;
}

/** The method the form starts on: the chosen one when it can run, else the kind's first that can. */
export function defaultMethod(kindChoice, chosen) {
  const ok = kindChoice.methods.filter((m) => m.available);
  return ok.find((m) => m.method === chosen) || ok[0] || null;
}

/** The settings a request carries: the values typed (as text; the server converts by type), empty ones left out. */
export function settingsOf(method, typed = {}) {
  const out = {};
  for (const s of method.settings) {
    const v = typed[s.name];
    if (v === undefined || v === null || v === "") continue;
    out[s.name] = v;
  }
  return out;
}

/** The start points a job can take (the x0 of POST /api/jobs); `job`: a job with a best point, or null. */
export function startPoints(job) {
  const out = [
    { id: "x", title: "the working x (the header's χ²)", x0: undefined },
    { id: "text", title: "the values of the project text", x0: "text" },
  ];
  if (job?.best?.x) out.push({ id: `job:${job.id}`, title: `the best point of job ${job.id}`, x0: { job: job.id } });
  return out;
}

/** A row's facts (name -> value, GET /api/jobs/kinds) as words: "covariance: yes · MINOS can follow: no". */
export function factsText(facts = {}) {
  const word = (v) => (v === true ? "yes" : v === false ? "no" : String(v));
  return Object.entries(facts).map(([k, v]) => `${k}: ${word(v)}`).join(" · ");
}

/** The numbers of a job's progress every job sends; others (a fitter's stage, a sampler's acceptance) are shown by
 * their names, as sent. */
export const PROGRESS_KEYS = ["calls", "objective", "best", "elapsed"];

/** The view of a job's progress (the progress reducer): what the running panel shows, from the job's record. */
export function progressView(rec) {
  const p = rec.progress || {};
  const state = STATES[rec.state] || { words: rec.state, final: false };
  return {
    id: rec.id, title: `${rec.kind} ${rec.method}`, state: rec.state, words: state.words, final: state.final,
    calls: p.calls ?? 0, objective: p.objective ?? null, best: p.best ?? rec.best?.objective ?? null,
    elapsed: p.elapsed ?? null, build: rec.build_s ?? null,
    canStop: !state.final && rec.state !== "stopping", canKill: !state.final,
    more: Object.entries(p).filter(([k]) => !PROGRESS_KEYS.includes(k)),
    trace: Array.isArray(rec.trace) ? rec.trace : [], log: Array.isArray(rec.log) ? rec.log : [],
  };
}

/** Whether a row of the job table follows the job `rec`: its `follows` names the job's kind and its method (or any
 * method of that kind, null), or is a list of such pairs (the sampler follows a fit, and a sampler: its chain). */
export function follows(row, rec) {
  if (!Array.isArray(row.follows) || !rec) return false;
  const pairs = Array.isArray(row.follows[0]) ? row.follows : [row.follows];
  return pairs.some(([kind, method]) => kind === rec.kind && (method === null || method === rec.method));
}

/** The words of a band job's result about what it varied: "energy_5, width_5_1 (2 of 3 free parameters)". */
export function variedText(result) {
  const varied = result?.varied;
  if (!Array.isArray(varied)) return "";
  const free = Array.isArray(result.free) ? result.free.length : null;
  const left = Array.isArray(result.free) ? result.free.filter((n) => !varied.includes(n)) : [];
  return `${varied.join(", ") || "none"} (${varied.length}${free === null ? "" : ` of ${free}`} free parameters`
    + `${left.length ? `; held at the start point: ${left.join(", ")}` : ""})`;
}

/** The rows of a sampler's summary (its result: mean, errors, quantiles, tau, effective_draws, rhat by name). */
export function chainRows(result) {
  return Object.keys(result.quantiles || {}).map((name) => ({ name, mean: result.mean?.[name] ?? null,
    sd: result.errors?.[name] ?? null, ...Object.fromEntries(Object.entries(result.quantiles[name])
      .map(([q, v]) => [`q${q}`, v])), tau: result.tau?.[name] ?? null,
    n_eff: result.effective_draws?.[name] ?? null, rhat: result.rhat?.[name] ?? null }));
}

/** The followers of a job's result shown as buttons, and the kinds folded behind "More": in the first view
 * (`firstViewOnly`) only those of the kinds the job table marks for it (`first_view`, the rule of the job form),
 * else all. */
export function followersShown(followers, firstViewOnly) {
  const shown = firstViewOnly ? followers.filter((f) => f.first_view) : followers;
  return { shown, folded: [...new Set(followers.filter((f) => !shown.includes(f)).map((f) => f.kind))] };
}

/** The view of a finished job: the recipe's own words, its numbers, and what can be done with it. */
export function resultView(rec, table = []) {
  const r = rec.result || {};
  const followers = table.filter((row) => row.available && follows(row, rec) && (!row.takes || r[row.takes]));
  return {
    status: r.status || rec.state, words: STATES[rec.state]?.words || rec.state, message: r.message || r.error || "",
    calls: r.calls ?? null, elapsed: r.elapsed ?? null, startObjective: r.start_objective ?? null,
    bestObjective: rec.best?.objective ?? null, edm: r.edm ?? null, valid: r.valid ?? null,
    hasBest: Boolean(rec.best?.x), covariance: r.covariance || null, correlation: r.correlation || null,
    errorsFrom: r.errors ? (r.errors_from || "the recipe") : "", minos: r.minos || null, followers,
    notices: Array.isArray(r.notices) ? r.notices : [], fitter: r.fitter || null, outcome: r.outcome ?? null,
    minosInvalid: Object.entries(r.minos || {}).filter(([, m]) => m.valid === false).map(([n]) => n),
  };
}

/** A square matrix of the server ({names, matrix}) as table rows: one per name, a column per name. */
export function matrixRows(m) {
  return m.names.map((n, i) => Object.fromEntries([["name", n], ...m.names.map((c, j) => [c, m.matrix[i][j]])]));
}

/** The CSV of a matrix: a header of names, one line per name; numbers as sent. */
export function matrixCsv(m) {
  const quote = (s) => (/[",\n]/.test(String(s)) ? `"${String(s).replace(/"/g, '""')}"` : String(s));
  return [["", ...m.names].map(quote).join(","),
    ...m.names.map((n, i) => [quote(n), ...m.matrix[i].map(String)].join(","))].join("\n") + "\n";
}

/** The job a visit of the tab shows: the one named, else a running one, else the newest. */
/** The tab of a job (the row of the job table of its kind and method: `JobRow.tab`; "fit" where none says). */
export function tabOfJob(job, table) {
  const row = table.find((r) => r.kind === job.kind && r.method === job.method)
    || table.find((r) => r.kind === job.kind);
  return row?.tab || "fit";
}

/** The jobs a tab lists: those of its own rows of the job table (the Fit tab: fits, MINOS, the sampler, the bands;
 * the ENDF tab's and the Rate tab's jobs are theirs). */
export function jobsOfTab(jobs, table, tab) {
  return jobs.filter((j) => tabOfJob(j, table) === tab);
}

export function jobToShow(jobs, id) {
  if (id) return jobs.find((j) => j.id === id) || null;
  const running = jobs.find((j) => !isFinal(j.state));
  if (running) return running;
  return jobs.slice().sort((a, b) => (b.started || 0) - (a.started || 0))[0] || null;
}
