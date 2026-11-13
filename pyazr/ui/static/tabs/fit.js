// The Fit tab (docs/upgrade/web_ui.md section 6.4): the job form and the job shown (its progress while it runs, its
// result when it has ended), the table of every parameter, the covariance and correlation of a fit, and the list of
// the project's jobs.  Routes: #/fit (the running job, else the newest), #/fit/job/<id>.
//
// The page holds no state of the project.  A free parameter's value typed here changes the working x (PATCH /api/x;
// the header's χ² follows); a fixed one is a value of the project text, edited in its record or level.  A running
// job is followed with one long poll in flight (GET /api/jobs/{id}?wait=1: the answer comes when its state changes,
// or after a second with its progress); a reload finds it again in GET /api/jobs.

import { pathOf } from "../api.js";
import { filterParameters, isFinal, jobToShow, kindsOf, matrixCsv, matrixRows, parameterRows,
  shownColumns } from "../fitting.js";
import { formatNumber, renderTable } from "../tables.js";
import { button, download, h, section } from "../ui.js";
import { jobPanel, jobsList, jobForm } from "./fitjob.js";
import { recordRoute } from "./records.js";

export async function render(ctx, el, args) {
  const view = ctx.fitView || (ctx.fitView = { kind: "", freeOnly: false, filter: "", compare: false, form: {},
    evals: {}, written: {} });
  if (!ctx.kinds) ctx.kinds = await ctx.api.get("/api/kinds");
  const [table, jobs, records] = await Promise.all([ctx.api.get("/api/jobs/kinds"), ctx.api.get("/api/jobs"),
    ctx.api.get("/api/project/records")]);
  const listed = jobToShow(jobs, args[0] === "job" ? args[1] : null);
  const job = listed ? await ctx.api.get(pathOf("/api/jobs", listed.id)) : null;

  // the parameter rows need the session of the current text; refused only while a job holds the budget and the
  // session was never built (a reload after a restart)
  let params = null, x = {}, refused = null;
  try {
    const answers = await Promise.all([ctx.api.get("/api/parameters", { quiet: true }),
      ctx.api.get("/api/x", { quiet: true })]);
    [params, x] = [answers[0], answers[1].x];
  } catch (e) {
    refused = e.view || { what: String(e) };
  }
  if (job?.accepted?.at_best && !view.evals[job.id]) view.evals[job.id] = job.accepted.at_best;
  if (job && isFinal(job.state) && job.best?.x && job.text_rev === ctx.project.rev && !view.evals[job.id]) {
    view.evals[job.id] = await evaluateAt(ctx, job.best.x);
  }
  if (job?.accepted?.rev === ctx.project.rev && !view.written[job.id]) {
    // accepted at the current revision: the project as written, at its own values (x0)
    const xs = await ctx.api.get("/api/x");
    const e = await evaluateAt(ctx, xs.x0);
    if (e.chi2 !== undefined) view.written[job.id] = { ...e, rev: ctx.project.rev };
  }

  const left = h("div", {}, jobForm(ctx, table, job, view), jobsList(ctx, jobs, job));
  const right = h("div", {}, job ? jobPanel(ctx, job, table, view)
    : section("Job shown", h("p.hint", {}, "No job yet: choose a recipe and press Run. A job runs in its own process "
      + "on the text of the current revision; the page follows its progress, and Stop ends it before its next "
      + "objective call with the best point.")));
  el.append(h("div.two-columns.fit-top", {}, left, right));
  el.appendChild(params ? parametersSection(ctx, params, x, job, records, view)
    : section("Parameters", h("div.warn-box", {}, `The parameter table needs the session of the current text: `
      + `${refused?.what || ""}`)));
  const fitted = job?.result?.covariance;
  if (fitted) el.appendChild(matricesSection(job));
  if (job && !isFinal(job.state)) follow(ctx, job, right);
}

/** The engine's χ², objective and N at a point x ({name: value}); null when the engine refuses. */
async function evaluateAt(ctx, x) {
  try {
    const e = await ctx.api.post("/api/evaluate", { x });
    return { chi2: e.chi2, objective: e.objective, n: e.totals?.["Total-N"] ?? null, rev: e.rev };
  } catch (err) {
    return { refused: err.view?.what || String(err) };
  }
}

/** Follow a running job: one long poll at a time, the panel updated in place; at its end the page re-reads. */
async function follow(ctx, rec, holder) {
  const panel = holder.querySelector(".job-panel");
  if (!panel || panel.dataset.following) return;
  panel.dataset.following = "1";
  let r = rec;
  while (panel.isConnected && !isFinal(r.state)) {
    try {
      r = await ctx.api.get(`${pathOf("/api/jobs", r.id)}?wait=1`, { background: true });
    } catch (e) {
      ctx.showError(e);
      return;
    }
    if (panel.isConnected && panel.update) panel.update(r);
    // a pause the page's clock counts between two polls (headless Chrome's virtual time stands still while a
    // request is in flight, and advances only here)
    if (!isFinal(r.state)) await new Promise((done) => setTimeout(done, 200));
  }
  if (panel.isConnected) {
    ctx.status(`job ${r.id}: ${r.state}${r.result?.message ? ` (${r.result.message})` : ""}`,
      r.state === "failed" ? "error" : "ok");
    await ctx.refresh();
  }
}

/** The parameter table: kind filter (the kinds the server sent), free only, text filter, the job's columns. */
function parametersSection(ctx, params, x, job, records, view) {
  const kinds = kindsOf(params, ctx.kinds?.parameters || {});
  const kindSelect = h("select", { "aria-label": "parameter kind" },
    h("option", { value: "" }, `every kind (${params.length})`),
    kinds.map((k) => h("option", { value: k.kind, title: k.doc }, `${k.kind} (${k.n}, ${k.free} free)`)));
  kindSelect.value = view.kind;
  kindSelect.addEventListener("change", () => { view.kind = kindSelect.value; draw(); });
  const freeOnly = h("input", { type: "checkbox", "aria-label": "free parameters only" });
  freeOnly.checked = view.freeOnly;
  freeOnly.addEventListener("change", () => { view.freeOnly = freeOnly.checked; draw(); });
  const compare = h("input", { type: "checkbox", "aria-label": "compare with start", disabled: !job?.result?.start_x });
  compare.checked = view.compare && Boolean(job?.result?.start_x);
  compare.addEventListener("change", () => { view.compare = compare.checked; draw(); });
  const nFree = params.filter((p) => p.free).length;
  const kindDoc = h("div.hint.kind-doc");
  const bar = h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "kind"), kindSelect),
    h("label.ctl", {}, freeOnly, h("span.ctl-label", {}, "free only")),
    h("label.ctl", { title: "the job's start point beside its best point" }, compare,
      h("span.ctl-label", {}, "compare with start")),
    h("span.hint", {}, `${params.length} parameters, ${nFree} free`
      + (job?.best?.x ? ` · fit, error: job ${job.id}` : "")),
    h("span.spacer"),
    button("Reset x", "Set the working x back to the values of the project text (PATCH /api/x, reset)",
      () => ctx.edit("reset x to the text's values", () => ctx.api.patch("/api/x", { reset: true })).catch(() => {})));
  const holder = h("div.parameter-table");
  const box = section("Parameters", bar, kindDoc, holder, h("p.hint", {}, "x is the engine's fit vector: level "
    + "energies in MeV, widths as reduced-width amplitudes (MeV^1/2), nuisances in their own units. Each column "
    + "head says where its words come from (hover). Enter in an x box sets that value of the working x; Esc "
    + "puts it back."));
  const recordByName = Object.fromEntries(records.map((r) => [r.name, r]));
  const all = parameterRows(params, { x, job });

  function draw() {
    kindDoc.textContent = kinds.find((k) => k.kind === view.kind)?.doc || "";
    const rows = filterParameters(all, view);
    const columns = shownColumns(all, view).map((c) => decorate(c));
    const t = renderTable(document, { columns, rows, rowKey: (r) => r.name, label: "parameters",
      caption: "Every parameter of the session: free ones are the fit vector x", filter: view.filter });
    t.filterInput.addEventListener("input", () => { view.filter = t.filterInput.value; });
    holder.textContent = "";
    holder.appendChild(t.element);
  }

  function decorate(c) {
    const col = { ...c };
    if (c.edit) col.render = (r) => (r.free ? xEditor(ctx, r) : h("span.hint", {}, ""));
    if (c.key === "free") col.render = (r) => (r.free ? h("span", {}, "free") : fixedCell(ctx, r, recordByName));
    if (c.key === "kind") col.render = (r) => h("span", { title: ctx.kinds?.parameters?.[r.kind]?.doc || "" }, r.kind);
    return col;
  }
  draw();
  return box;
}

/** A fixed parameter: a value of the project text; a link to its record where the name begins with one. */
function fixedCell(ctx, r, recordByName) {
  const owner = recordByName[r.name.split(/[.[]/)[0]];
  const why = "fixed: a value of the project text (a project edit, one history entry); levels and channels are "
    + "edited in the Nucleus tab, records in the Data tab";
  if (!owner) return h("span", { title: why }, "fixed");
  return h("span", { title: why }, "fixed ", h("a", { href: recordRoute(ctx.schema, owner),
    title: `edit the record ${owner.name}` }, "→ record"));
}

/** The box of a free parameter's working x: Enter sends it (PATCH /api/x); the status line gives the χ² change. */
function xEditor(ctx, r) {
  const text = String(r.x);
  const input = h("input.num", { type: "text", value: text, size: 14, "aria-label": `x of ${r.name}` });
  input.addEventListener("keydown", (e) => {
    if (e.key === "Escape") { input.value = text; return; }
    if (e.key !== "Enter" || input.value.trim() === text) return;
    e.preventDefault();
    const before = ctx.evaluation?.chi2;
    ctx.edit(`x of ${r.name}`, () => ctx.api.patch("/api/x", { values: { [r.name]: Number(input.value) } }))
      .then(() => {
        const after = ctx.evaluation?.chi2;
        if (before !== undefined && after !== undefined) {
          ctx.status(`x of ${r.name} = ${input.value}: χ² ${formatNumber(before, 10)} → ${formatNumber(after, 10)}`,
            "ok");
        }
      }, () => {});
  });
  return input;
}

/** The covariance and the correlation of a fit (the recipe's, by name), each as a table and as CSV. */
function matricesSection(job) {
  const r = job.result;
  const one = (title, m, digits, name) => {
    const cols = [{ key: "name", title: "" }, ...m.names.map((n) => ({ key: n, title: n, numeric: true, digits,
      sortable: false }))];
    const t = renderTable(document, { columns: cols, rows: matrixRows(m), rowKey: (row) => row.name,
      withFilter: false, label: title, caption: `${title} of job ${job.id} (fit basis, order of x)` });
    return h("div.matrix", {}, t.element, button(`${title} (CSV)`, `The ${title} as sent, as a CSV file`,
      () => download(`${name}_${job.id}.csv`, matrixCsv(m))));
  };
  const details = h("details.matrices", {}, h("summary", {}, `Covariance and correlation (${r.covariance.names.length}`
    + ` free parameters, from ${r.errors_from || "the recipe"})`),
  one("correlation", r.correlation, 3, "correlation"), one("covariance", r.covariance, 4, "covariance"));
  return section("Covariance", details);
}
