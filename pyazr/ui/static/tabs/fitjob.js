// The job parts of the Fit tab (docs/upgrade/web_ui.md sections 4.4 and 6.4): the job form generated from the
// server's job table (kinds, methods, and each method's settings with their types and defaults), the panel of the
// job shown (progress, the objective against call, the last lines of its log, Stop and Kill; at its end the
// recipe's words, the χ² at its best point, Accept into project, Compare with start, the jobs that follow it), and the
// list of the project's jobs.  A row of the job table that waits for pyazr is shown with what it waits for.

import { pathOf } from "../api.js";
import { defaultMethod, factsText, isFinal, jobChoices, progressView, resultView, settingsOf,
  startPoints } from "../fitting.js";
import { plotOf, plotlyFigure } from "../plots.js";
import { formatNumber, renderTable } from "../tables.js";
import { button, drawPlot, h, section, themeColours } from "../ui.js";

/** One widget per type of a job's setting (pyazr/ui/jobs.py, SETTING_TYPES): (setting, typed text) -> {el, read}. */
export const SETTING_WIDGETS = {
  integer: (s, typed) => {
    const input = h("input.num", { type: "text", inputmode: "numeric", size: 10, value: typed ?? (s.default ?? ""),
      placeholder: s.default === null ? "default" : String(s.default), "aria-label": s.name });
    return { el: input, read: () => input.value.trim() };
  },
  boolean: (s, typed) => {
    const input = h("input", { type: "checkbox", "aria-label": s.name });
    input.checked = typed === undefined ? Boolean(s.default) : typed === "true";
    return { el: input, read: () => String(input.checked) };
  },
  names: (s, typed) => {
    const input = h("input.mono", { type: "text", size: 36, value: typed ?? "", placeholder: "every free parameter",
      "aria-label": s.name });
    return { el: input, read: () => input.value.trim() };
  },
};

/** The job form: kind, method (those that wait for pyazr disabled, with the reason), settings, start, threads, Run. */
export function jobForm(ctx, table, job, view) {
  const choices = jobChoices(table);
  const form = view.form;
  const kind = choices.find((k) => k.kind === form.kind && k.available) || choices.find((k) => k.available);
  const method = defaultMethod(kind, form.method);
  const kindSelect = h("select", { "aria-label": "job kind" }, choices.map((k) => h("option",
    { value: k.kind, disabled: !k.available, title: k.available ? "" : k.methods[0].reason },
    k.available ? k.kind : `${k.kind} (waits)`)));
  kindSelect.value = kind.kind;
  const methodSelect = h("select", { "aria-label": "method" }, kind.methods.map((m) => h("option",
    { value: m.method, disabled: !m.available, title: m.available ? m.doc : m.reason },
    m.available ? m.method : `${m.method} — ${m.reason}`)));
  methodSelect.value = method.method;
  const redraw = () => ctx.navigate(location.hash);
  kindSelect.addEventListener("change", () => { view.form = { kind: kindSelect.value }; redraw(); });
  methodSelect.addEventListener("change", () => { view.form = { kind: kind.kind, method: methodSelect.value };
    redraw(); });

  const typed = form.typed || {};
  const widgets = method.settings.map((s) => {
    const w = SETTING_WIDGETS[s.type] ? SETTING_WIDGETS[s.type](s, typed[s.name])
      : { el: h("span.warn-box", {}, `no widget for the setting type ${s.type}`), read: () => "" };
    w.el.addEventListener("change", () => { form.typed = { ...(form.typed || {}), [s.name]: w.read() }; });
    return [s, w];
  });
  const points = startPoints(job);
  const follows = method.follows && job && job.kind === method.follows[0] && job.method === method.follows[1];
  const startSelect = h("select", { "aria-label": "start point" }, points.map((p) => h("option", { value: p.id },
    p.title)));
  startSelect.value = form.start && points.some((p) => p.id === form.start) ? form.start
    : (follows && job.best?.x ? `job:${job.id}` : "x");
  startSelect.addEventListener("change", () => { form.start = startSelect.value; });
  const budget = ctx.describe?.budget || {};
  const run = button("Run", "Start the job (POST /api/jobs)", async () => {
    const values = Object.fromEntries(widgets.map(([s, w]) => [s.name, w.read()]));
    form.typed = values;
    const point = points.find((p) => p.id === startSelect.value);
    const body = { kind: kind.kind, method: method.method, settings: settingsOf(method, values) };
    if (point.x0 !== undefined) body.x0 = point.x0;
    try {
      const started = await ctx.api.post("/api/jobs", body);
      ctx.status(`job ${started.id} started: ${kind.kind} ${method.method}`, "ok");
      ctx.navigate(`#/fit/job/${started.id}`);
    } catch (e) {
      ctx.showError(e, box);
    }
  }, "primary");
  run.dataset.action = "run";
  const box = section("Job",
    h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "kind"), kindSelect),
      h("label.ctl", {}, h("span.ctl-label", {}, "method"), methodSelect)),
    h("p.hint", {}, method.doc), Object.keys(method.facts).length ? h("p.hint.job-facts", {}, factsText(method.facts))
      : null,
    h("div.inline-form.settings", {}, widgets.map(([s, w]) => h("label.ctl", { title: `${s.doc} (${s.allowed})` },
      h("span.ctl-label", {}, s.name.replace(/_/g, " ")), w.el))),
    h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "start"), startSelect),
      h("span.hint", { title: "each job runs with one thread of the server's budget (pyazr ui --threads)" },
        `threads: 1 of ${budget.threads ?? "?"} (${budget.job_slots ?? "?"} job slot(s))`), h("span.spacer"), run),
    jobTable(table));
  return box;
}

/** Every row of the server's job table with what it runs or what it waits for (the rows that wait are not runs). */
function jobTable(table) {
  const rows = table.map((r) => ({ key: `${r.kind}/${r.method}`, kind: r.kind, method: r.method, doc: r.doc,
    status: r.available ? "runs today" : (r.missing ? `waits for pyazr: ${r.missing}` : `needs ${r.needs}`),
    facts: factsText(r.facts), _class: r.available ? "" : "unavailable" }));
  const t = renderTable(document, { rows, rowKey: (r) => r.key, withFilter: false, label: "job table",
    caption: "The server's job table (GET /api/jobs/kinds)",
    columns: [{ key: "kind", title: "kind" }, { key: "method", title: "method" },
      { key: "facts", title: "gives", wrap: true }, { key: "status", title: "status", wrap: true }] });
  return h("details.job-table", {}, h("summary", {}, `job table: ${table.filter((r) => r.available).length} of `
    + `${table.length} rows run today`), t.element);
}

/** The panel of a job: progress while it runs (updated in place by `panel.update(record)`), its result at its end. */
export function jobPanel(ctx, rec, table, view) {
  const panel = h("section.panel.job-panel", { "aria-label": `job ${rec.id}` });
  const head = h("h2");
  const state = h("div.job-state", { role: "status" });
  const facts = h("dl.facts.job-progress");
  const stop = button("Stop", "Stop between objective calls: the job ends with its best point (POST .../stop)",
    () => ctx.api.post(pathOf("/api/jobs", rec.id, "stop"), {}).then((r) => panel.update(r), (e) => ctx.showError(e)));
  const kill = button("Kill", "End the job's process now: the best point reported so far is kept (POST .../kill)",
    () => ctx.api.post(pathOf("/api/jobs", rec.id, "kill"), {}).then((r) => panel.update(r), (e) => ctx.showError(e)),
    "danger");
  stop.dataset.action = "stop";
  kill.dataset.action = "kill";
  const plotBox = h("div.progress-plot");
  const log = h("pre.job-log", { "aria-label": "the last lines of the job's log" });
  panel.append(head, state, h("div.job-numbers", {}, facts, h("div.job-buttons", {}, stop, kill)), plotBox,
    h("div.hint", {}, "the job's log (the engine's notes), last lines:"), log);

  panel.update = (r) => {
    const v = progressView(r);
    head.textContent = `Job ${v.id} · ${v.title}`;
    state.textContent = v.words;
    state.dataset.state = v.state;
    state.className = `job-state ${v.final ? (v.state === "failed" ? "error" : "done") : "live"}`;
    const n = (x, d) => (x === null || x === undefined ? "—" : formatNumber(x, d));
    facts.textContent = "";
    facts.dataset.calls = String(v.calls);
    for (const [k, value, title] of [["calls", String(v.calls), "objective calls (with gradient) so far"],
      ["objective", n(v.objective, 10), "χ² + priors at the last call (the engine's)"],
      ["best", n(v.best, 10), "the lowest objective so far"],
      ["time", v.elapsed === null ? "—" : `${v.elapsed.toFixed(1)} s`, "since the recipe started"],
      ["session build", v.build === null ? "—" : `${v.build.toFixed(2)} s`, "the job's own session"],
      ...v.more.map(([k, x]) => [k.replace(/_/g, " "), typeof x === "number" ? n(x, 6) : (x ?? "—"),
        "as the job reports it"])]) {
      facts.append(h("dt", { title }, k), h("dd", {}, value));
    }
    stop.disabled = !v.canStop;
    kill.disabled = !v.canKill;
    drawTrace(plotBox, r);
    log.textContent = v.log.join("\n") || "(empty)";
  };
  panel.update(rec);
  if (isFinal(rec.state)) panel.appendChild(resultBlock(ctx, rec, table, view));
  return panel;
}

/** The objective against call (log y), redrawn in place while the job runs; the trace's table without plotly. */
function drawTrace(el, rec) {
  if (!rec.trace?.length) { el.textContent = ""; return; }
  const plot = plotOf(rec, {});
  const area = el.querySelector(".plot-area");
  if (area && globalThis.Plotly?.react) {
    const fig = plotlyFigure(plot, themeColours());
    globalThis.Plotly.react(area, fig.data, fig.layout, fig.config);
    return;
  }
  el.dataset.drawn = drawPlot(el, plot, `job_${rec.id}_trace`);
}

/** The result of an ended job: the recipe's words and numbers, the χ² at the best point, Accept, the followers. */
function resultBlock(ctx, rec, table, view) {
  const v = resultView(rec, table);
  const ev = view.evals[rec.id];
  const written = view.written[rec.id];
  const n = (x, d = 10) => (x === null || x === undefined ? "—" : formatNumber(x, d));
  const box = h("div.job-result", { data: { status: v.status } });
  box.append(h("h3", {}, "Result"), h("p.recipe-words", {}, h("strong", {}, v.words),
    v.message ? ` · in its own words: “${v.message}”` : ""));
  const facts = h("dl.facts", {}, h("dt", {}, "objective"), h("dd", {}, `${n(v.startObjective)} at the start → `
    + `${n(v.bestObjective)} at the best point`),
  h("dt", {}, "calls, time"), h("dd", {}, `${v.calls ?? "—"} calls, ${v.elapsed === null ? "—"
    : `${v.elapsed.toFixed(1)} s`}${v.edm !== null ? ` · EDM ${n(v.edm, 3)}` : ""}`),
  h("dt", {}, "errors"),
  h("dd", {}, v.errorsFrom || "none (this job gives none: see the job table's column \u201cgives\u201d)"));
  if (ev?.chi2 !== undefined) {
    facts.append(h("dt", {}, "at the best point"), h("dd", { data: { chi2Fit: String(ev.chi2) } },
      `χ² ${n(ev.chi2)} · N ${ev.n ?? "—"} (the engine's, on the text of rev ${rec.text_rev})`));
  } else if (ev?.refused) {
    facts.append(h("dt", {}, "at the best point"), h("dd", {}, `not evaluated: ${ev.refused}`));
  } else if (rec.text_rev !== ctx.project.rev) {
    facts.append(h("dt", {}, "at the best point"), h("dd", {}, `not evaluated: the job ran on rev ${rec.text_rev}, the `
      + `project is at rev ${ctx.project.rev}`));
  }
  box.appendChild(facts);
  if (written) box.appendChild(writtenBlock(ctx, ev, written));

  const actions = h("div.form-actions");
  if (v.hasBest) {
    const accept = button("Accept into project", "Write the best point into the project text: one history entry, "
      + "Undo takes it back (POST /api/jobs/{id}/accept)", () => accept_(ctx, rec), "primary");
    accept.dataset.action = "accept";
    actions.appendChild(accept);
  }
  for (const f of v.followers) {
    actions.appendChild(button(`${f.kind} (${f.method})…`, `${f.doc}: fills the job form, from this job's best point`,
      () => { view.form = { kind: f.kind, method: f.method, start: `job:${rec.id}` }; ctx.navigate(location.hash); }));
  }
  box.appendChild(actions);
  return box;
}

/** The fit's χ² beside the χ² of the project as written after Accept, both the engine's; a difference is said. */
function writtenBlock(ctx, ev, w) {
  // every digit as sent: the two numbers may differ in their last digits only
  const n = (x) => (x === null || x === undefined ? "—" : formatNumber(x));
  const box = h("div.written", { data: { chi2Written: String(w.chi2), rev: String(w.rev) } },
    h("dl.facts", {}, h("dt", {}, "the fit"), h("dd", {}, `χ² ${n(ev?.chi2)} · objective ${n(ev?.objective)}`),
      h("dt", {}, `as written (rev ${w.rev})`), h("dd", {}, `χ² ${n(w.chi2)} · objective ${n(w.objective)}`
        + ` · N ${w.n ?? "—"}`)));
  if (ev && ev.chi2 !== w.chi2) {
    box.appendChild(h("div.warn-box", {}, "The project as written gives another χ² than the fit at its best point; "
      + "both numbers are the engine's. Two known causes: the values written into the text carry a finite number of "
      + "digits (a difference in the last digits), and, an open point, a fitted target density, Doppler "
      + "temperature, stopping scale or kernel parameter is laid out on grids from its start value when the "
      + "project is read, so the written project is computed on other grids than the fit."));
  }
  if (ctx.project.rev !== w.rev) {
    box.appendChild(h("p.hint", {}, `The project is now at rev ${ctx.project.rev}: an Undo or another edit came after `
      + "the accepted one."));
  }
  return box;
}

/** Accept: the best point into the project text (one history entry); the tab, drawn again, evaluates the project as
 * written (the job's record names the revision its acceptance made). */
async function accept_(ctx, rec) {
  try {
    await ctx.edit(`accept job ${rec.id}`, () => ctx.api.post(pathOf("/api/jobs", rec.id, "accept"),
      { rev: ctx.api.rev }));
  } catch {
    return;
  }
  ctx.status(`accept job ${rec.id}: done; the χ² of the project as written is beside the fit's (Undo: Ctrl+Z)`, "ok");
}

/** The project's jobs, newest first; a click shows one. */
export function jobsList(ctx, jobs, shown) {
  const rows = jobs.slice().sort((a, b) => (b.started || 0) - (a.started || 0)).map((j) => ({
    // a running job's numbers are in its panel (followed); the list shows those of ended jobs
    id: j.id, what: `${j.kind} ${j.method}`, state: j.state, calls: isFinal(j.state) ? j.progress?.calls ?? "" : "",
    best: isFinal(j.state) ? j.progress?.best ?? "" : "",
    rev: j.text_rev, started: j.started ? new Date(j.started * 1000).toLocaleTimeString() : "" }));
  if (!rows.length) return h("div");
  const t = renderTable(document, { rows, rowKey: (r) => r.id, selected: shown?.id, withFilter: false, label: "jobs",
    caption: "Jobs of this project (GET /api/jobs); a click shows one",
    columns: [{ key: "id", title: "job" }, { key: "what", title: "recipe" }, { key: "state", title: "state" },
      { key: "calls", title: "calls", numeric: true }, { key: "best", title: "best", numeric: true, digits: 8,
        doc: "the best objective (χ² + priors)" }, { key: "started", title: "started" }],
    onSelect: (r) => ctx.navigate(`#/fit/job/${r.id}`) });
  return section("Jobs", t.element);
}
