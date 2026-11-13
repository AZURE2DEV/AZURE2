// The page's self-test (#selftest; docs/upgrade/web_ui.md section 9): visits every tab and the views that matter
// (a level, a data set's form and plot, the plot's table, a new record's form), records every script error,
// unhandled rejection, console error and failed request, and writes the result as JSON into
// <pre id="selftest-result" data-done="1">.  tests/pyazr/ui/smoke_test.py reads it from Chrome's dumped document.
// "#selftest?theme=dark" runs it in the dark theme; "#selftest?edit=1" also edits a width and undoes it.

import { LAYOUT } from "./layout.js";
import { NUCLEUS_SECTIONS } from "./tabs/records.js";

export function selftest(ctx, tabs, draw) {
  const errors = [];
  addEventListener("error", (e) => errors.push(`${e.message} (${e.filename}:${e.lineno})`));
  addEventListener("unhandledrejection", (e) => errors.push(`unhandled rejection: ${e.reason?.stack || e.reason}`));
  const consoleError = console.error.bind(console);
  console.error = (...a) => { errors.push(`console.error: ${a.map(String).join(" ")}`); consoleError(...a); };
  const extra = location.hash.includes("?") ? `?${location.hash.split("?")[1]}` : "";
  const join = (route) => (route.includes("?") ? `${route}&${extra.slice(1)}` : `${route}${extra}`);
  const steps = [];
  const edit = {};

  async function run() {
    while (!ctx.project || !ctx.schema) await new Promise((r) => setTimeout(r, 20));
    const t0 = performance.now();
    await ctx.refresh();
    const routes = ["#/project", "#/nucleus"];
    const levels = await ctx.api.get("/api/project/levels");
    if (levels.length) routes.push(`#/nucleus/level/${encodeURIComponent(levels[0].label)}`);
    const records = await ctx.api.get("/api/project/records");
    const nucleusKinds = NUCLEUS_SECTIONS.flatMap((s) => ctx.schema.sections[s] || []);
    const nuc = records.find((r) => nucleusKinds.includes(r.kind));
    if (nuc) routes.push(`#/nucleus/${encodeURIComponent(nuc.name)}`);
    routes.push("#/data");
    const evaluated = ctx.evaluation?.datasets?.[0]?.name;
    const at = `#/data/${encodeURIComponent(evaluated)}`;
    if (evaluated) routes.push(at, `${at}?plot=table`);
    // a record of a kind the layout hints do not name: its form is the schema's alone
    const other = records.find((r) => !nucleusKinds.includes(r.kind) && !LAYOUT[r.kind]);
    if (other) routes.push(`#/data/${encodeURIComponent(other.name)}`);
    routes.push("#/data/+new", ...tabs.filter((t) => t.later).map((t) => `#/${t.id}`));
    for (const route of routes) {
      const s = performance.now();
      location.hash = join(route);
      await new Promise((r) => setTimeout(r, 0));
      await draw();
      await new Promise((r) => setTimeout(r, 30));  // forms and descriptions that load after the draw
      const plot = document.querySelector(".plot[data-drawn]");
      steps.push({ route, ms: Math.round(performance.now() - s), plot: plot ? plot.dataset.drawn : null,
        forms: document.querySelectorAll("form.record-form").length,
        errorBoxes: [...document.querySelectorAll(".error-box")].map((b) => b.textContent.slice(0, 300)) });
    }
    if (/[?&]edit=1/.test(extra)) await editAndUndo(levels);
    return Math.round(performance.now() - t0);
  }

  /** "#selftest?edit=1": a free channel's width typed into its box (Enter), then the Undo button, as a person
   * would; records the header's χ² before, after and back (the engine's numbers) in `edit`. */
  async function editAndUndo(levels) {
    const level = levels.find((l) => l.channels.some((c) => !c.fixed && c.width !== 0));
    if (!level) { edit.skipped = "no level with a free, non-zero width"; return; }
    if (!ctx.evaluation) { edit.skipped = "the engine does not evaluate this project"; return; }
    const channel = level.channels.find((c) => !c.fixed && c.width !== 0);
    location.hash = join(`#/nucleus/level/${encodeURIComponent(level.label)}`);
    await new Promise((r) => setTimeout(r, 0));
    await draw();
    const until = async (test) => {
      for (let i = 0; i < 600 && !test(); i += 1) await new Promise((r) => setTimeout(r, 50));
      return test();
    };
    Object.assign(edit, { level: level.label, channel: channel.channel, width: channel.width_text,
      chi2_before: ctx.evaluation?.chi2, rev_before: ctx.project.rev });
    const input = document.querySelector(`input[aria-label="width of channel ${channel.channel}"]`);
    input.value = String(channel.width * 1.25);
    edit.new_width = input.value;
    input.dispatchEvent(new KeyboardEvent("keydown", { key: "Enter", bubbles: true }));
    edit.edited = await until(() => ctx.statusLine.textContent.includes("→"));
    Object.assign(edit, { chi2_after: ctx.evaluation?.chi2, rev_after: ctx.project.rev,
      status: ctx.statusLine.textContent });
    ctx.undoButton.click();
    edit.undone = await until(() => ctx.statusLine.textContent.startsWith("undo: done"));
    Object.assign(edit, { chi2_back: ctx.evaluation?.chi2, rev_back: ctx.project.rev });
  }

  run().then((ms) => finish(ms), (e) => { errors.push(`self-test: ${e?.stack || e}`); finish(null); });

  function finish(ms) {
    const out = document.createElement("pre");
    out.id = "selftest-result";
    out.textContent = JSON.stringify({ done: true, ms, plotly: typeof globalThis.Plotly === "object",
      errors, failed: ctx.api.failures, steps, edit, chi2: ctx.evaluation?.chi2 ?? null,
      evaluation_refused: ctx.evalError?.what ?? null }, null, 1);
    out.dataset.done = "1";
    out.style.display = "none";
    document.body.appendChild(out);
    document.title = `selftest: ${errors.length || ctx.api.failures.length ? "FAILED" : "done"}`;
  }
}
