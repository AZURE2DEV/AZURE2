// The page's self-test (#selftest; docs/upgrade/web_ui.md section 9): visits every tab and the views that matter
// (a level, a data set's form and plot, the plot's table, a new record's form), records every script error,
// unhandled rejection, console error and failed request, and writes the result as JSON into
// <pre id="selftest-result" data-done="1">.  tests/pyazr/ui/smoke_test.py reads it from Chrome's dumped document.
// "#selftest?theme=dark" runs it in the dark theme; "#selftest?edit=1" also edits a width and undoes it, and adds a
// pair with the Nucleus tab's form, removes it with its button and undoes both.  Where the project folder holds ENDF-6
// files, the ENDF tab draws them over the data and presses Check on the first.  The views are visited in the
// expert's depth (every tab and field), then each tab of the first view (FIRST_VIEW) and Calculate twice.

import { FIRST_VIEW, LAYOUT } from "./layout.js";
import { fitCheck } from "./selftestfit.js";
import { NUCLEUS_SECTIONS } from "./tabs/records.js";

/** The Physics tab's views the self-test visits: the first row as it opens, a Coulomb function with elements
 * outside the domain (rho <= 0) and four lines of l, its table view, the session's channel functions, and a
 * stopping power with SRIM's citation. */
export const PHYSICS_ROUTES = ["#/physics",
  "#/physics/coulomb_waves?from=-2&to=20&points=111&spacing=linear&v.l=0 1 2 3&fields=F,G",
  "#/physics/coulomb_waves?plot=table", "#/physics/channel_functions",
  "#/physics/stopping_power?logy=1"];

/** The THM tab's views, where the project has THM content: the tab as it opens (the experiments' profile), then
 * each THM diagnostic of the session's rows (a refusal of the engine, a diagnostic the experiment does not have, is
 * recorded with its step). */
async function thmRoutes(ctx, tabs) {
  const tab = tabs.find((t) => t.id === "thm");
  if (!tab?.shown || !(await tab.shown(ctx))) return ["#/thm"];
  const rows = await ctx.api.get("/api/session/functions");
  return ["#/thm", ...Object.values(rows).filter((r) => r.tab === "thm").map((r) => `#/thm?fn=${r.name}`)];
}

export function selftest(ctx, tabs, draw) {
  const errors = [];
  addEventListener("error", (e) => errors.push(`${e.message} (${e.filename}:${e.lineno})`));
  addEventListener("unhandledrejection", (e) => errors.push(`unhandled rejection: ${e.reason?.stack || e.reason}`));
  const consoleError = console.error.bind(console);
  console.error = (...a) => { errors.push(`console.error: ${a.map(String).join(" ")}`); consoleError(...a); };
  let extra = location.hash.includes("?") ? `?${location.hash.split("?")[1]}` : "";
  if (!/[?&]depth=/.test(extra)) extra = `${extra}${extra ? "&" : "?"}depth=expert`;
  const join = (route) => (route.includes("?") ? `${route}&${extra.slice(1)}` : `${route}${extra}`);
  const first = (route) => join(route).replace(/depth=expert/, "depth=first");
  const steps = [];
  const edit = {};
  const pair = {};
  const endf = {};
  const preview = {};
  let fit = {};

  async function run() {
    while (!ctx.project || !ctx.schema) await new Promise((r) => setTimeout(r, 20));
    const t0 = performance.now();
    await ctx.refresh();
    const routes = ["#/project", "#/nucleus", "#/nucleus/pair/+new", "#/nucleus/level/+new"];
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
    routes.push("#/data/+new", "#/fit", ...tabs.filter((t) => t.later).map((t) => `#/${t.id}`));
    routes.push(...await thmRoutes(ctx, tabs));
    routes.push("#/results", "#/results?plot=table", "#/results?x=e_x&levels=1&logy=1");
    routes.push(...PHYSICS_ROUTES);
    // the ENDF tab: as it opens, the project folder's first ENDF-6 file over the data, and its check screen
    const endfFiles = (await ctx.api.get("/api/endf/files")).files.map((f) => f.path);
    routes.push("#/endf");
    if (endfFiles.length) {
      routes.push(`#/endf?files=${endfFiles.map(encodeURIComponent).join(",")}`,
        `#/endf/check?file=${encodeURIComponent(endfFiles[0])}`, `#/endf/new?file=${encodeURIComponent(endfFiles[0])}`);
    }
    routes.push("#/endf/export");
    for (const route of routes) await visit(route, join(route));
    if (endfFiles.length) await endfCheck(endfFiles[0]);
    if (endfFiles.length) await endfConvert(endfFiles[0]);
    await firstView();
    if (/[?&]edit=1/.test(extra)) await editAndUndo(levels);
    if (/[?&]edit=1/.test(extra)) await pairAddRemove();
    if (/[?&]fit=1/.test(extra)) fit = await fitCheck(ctx, draw, join);
    return Math.round(performance.now() - t0);
  }

  /** The first view: each of its tabs as it opens, a data set's form (its rare fields folded), and Calculate twice
   * (the second shows what changed: nothing, the text being the same). */
  async function firstView() {
    for (const v of FIRST_VIEW) await visit(v.hash, first(v.hash));
    const evaluated = ctx.evaluation?.datasets?.[0]?.name;
    const at = `#/data/${encodeURIComponent(evaluated)}`;
    if (evaluated) {
      await visit(at, first(at));
      // a file typed into the data set's file field: the preview of its rows and the engine's reading, before Apply
      const box = document.querySelector('[data-field="file"] input');
      if (box) {
        box.value = box.value.replace(/^"|"$/g, "");
        box.dispatchEvent(new Event("change", { bubbles: true }));
        for (let i = 0; i < 100 && !document.querySelector(".file-read > *"); i++) {
          await new Promise((r) => setTimeout(r, 20));
        }
        preview.text = document.querySelector(".file-read")?.textContent || "";
        preview.head = document.querySelector(".file-preview pre.file-head")?.textContent || "";
        preview.offered = document.querySelector('[data-field="file"] datalist')?.options.length || 0;
      }
      // "new …" of the data set's target: the form of a new target opens, the data set unchanged
      const target = document.querySelector('[data-field="target"] select');
      const choice = target && [...target.options].find((o) => o.value.startsWith("+new:"));
      if (choice) {
        target.value = choice.value;
        target.dispatchEvent(new Event("change", { bubbles: true }));
        await new Promise((r) => setTimeout(r, 0));
        await draw();
        await new Promise((r) => setTimeout(r, 30));
        preview.newRoute = location.hash;
        preview.newKind = document.querySelector('[aria-label="kind of the new record"]')?.value || "";
        preview.rev = ctx.project.rev;
      }
    }
    for (let i = 0; i < 2; i++) {
      location.hash = first("#/project");
      await new Promise((r) => setTimeout(r, 0));
      await draw();
      await ctx.calculate();
      await new Promise((r) => setTimeout(r, 0));
      await draw();
      await record("Calculate", performance.now());
    }
  }

  /** The ENDF tab's Check pressed on a file: the job followed until its result is drawn; records the counts the
   * page shows (the own check's heading and one heading per code) in `endf`. */
  async function endfCheck(file) {
    const until = async (test) => {
      for (let i = 0; i < 600 && !test(); i += 1) await new Promise((r) => setTimeout(r, 50));
      return test();
    };
    const route = `#/endf/check?file=${encodeURIComponent(file)}`;
    location.hash = join(route);
    await new Promise((r) => setTimeout(r, 0));
    await draw();
    const before = document.querySelector(".endf-check-summary")?.dataset.job || null;
    document.querySelector("section.check button.primary").click();
    endf.checked = await until(() => {
      const shown = document.querySelector(".endf-check-summary");
      return shown && shown.dataset.job !== before;
    });
    Object.assign(endf, { file, headings: [...document.querySelectorAll(".endf-check-job h3")]
      .map((x) => x.textContent), meanings: [...document.querySelectorAll(".endf-check-job td")]
      .filter((td) => td.textContent.includes("p. ")).length });
    await record(`${route} (checked)`, performance.now());
  }

  /** The ENDF tab's Convert pressed on a file (a new project from its MF=2): the job followed until its report is
   * drawn; records its heading and the "Write the new project" button in `endf.converted`, `endf.report`. */
  async function endfConvert(file) {
    const route = `#/endf/new?file=${encodeURIComponent(file)}`;
    location.hash = join(route);
    await new Promise((r) => setTimeout(r, 0));
    await draw();
    const convert = [...document.querySelectorAll("button.primary")].find((b) => b.textContent === "Convert");
    if (!convert || convert.disabled) { endf.converted = false; return; }
    convert.click();
    for (let i = 0; i < 600 && !document.querySelector(".endf-import-job h3"); i += 1) {
      await new Promise((r) => setTimeout(r, 50));
    }
    endf.report = document.querySelector(".endf-import-job h3")?.textContent || "";
    endf.converted = Boolean(endf.report) && [...document.querySelectorAll(".endf-keep button")]
      .some((b) => b.textContent === "Write the new project");
    await record(`${route} (converted)`, performance.now());
  }

  /** One view: the route set, the page drawn, then its step recorded. */
  async function visit(route, hash) {
    const s = performance.now();
    location.hash = hash;
    await new Promise((r) => setTimeout(r, 0));
    await draw();
    await new Promise((r) => setTimeout(r, 30));  // forms and descriptions that load after the draw
    await record(route, s);
  }

  /** A tab link's title (its key number and group mark left out). */
  function titleOf(a) {
    if (!a) return "";
    return [...a.childNodes].filter((n) => n.nodeType === 3).map((n) => n.textContent).join("").trim();
  }
  /** The tab shown: a group's member where one is active, else the main tab. */
  function activeTab() {
    return titleOf(document.querySelector("nav.subtabs-row:not([hidden]) a.tab.active")
      || document.querySelector("nav.tabs:not(.subtabs-row) a.tab.active"));
  }

  async function record(route, s) {
    {
      const plot = document.querySelector(".plot[data-drawn]");
      const physics = document.querySelector(".physics-status");
      steps.push({ route, at: Math.round(s), ms: Math.round(performance.now() - s),
        plot: plot ? plot.dataset.drawn : null,
        forms: document.querySelectorAll("form.record-form").length,
        physics: physics ? { computed: +physics.dataset.computed, elements: +physics.dataset.elements,
          text: physics.textContent, python: document.querySelector("pre.python")?.textContent || "",
          citation: document.querySelector("p.citation")?.dataset.citation || "" } : null,
        tables: [...document.querySelectorAll(".physics-out table[data-table], .physics-out .table-wrap")].length,
        facts: document.querySelectorAll("dl.function-facts dd").length,
        tabListed: !document.querySelector('a.tab[data-tab="thm"]')?.hidden,
        // the main row's tabs and every group's tabs as listed (a group's own row shows under its main tab)
        tabsShown: [...document.querySelectorAll("nav.tabs a.tab:not([hidden])")].map((a) => titleOf(a)),
        mainTab: titleOf(document.querySelector("nav.tabs:not(.subtabs-row) a.tab.active")),
        jobKinds: [...document.querySelectorAll('[aria-label="job kind"] option')].map((o) => o.value),
        moreKinds: Boolean(document.querySelector('[data-action="more-kinds"]')),
        depth: document.documentElement.dataset.depth, tab: activeTab(),
        hash: location.hash, statusMarks: { job: document.querySelector(".job-mark:not([hidden]) a")?.textContent
          || "", validation: document.querySelector("button.validation")?.textContent || "" },
        folded: document.querySelectorAll("details.form-more").length,
        changes: document.querySelectorAll('[aria-label^="changes"]').length,
        errorBoxes: [...document.querySelectorAll(".error-box")].map((b) => b.textContent.slice(0, 300)) });
    }
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
    // Calculate (the header's button) after the edit: the change since the first view's Calculate, the data set
    // that moved most first
    document.querySelector("button.calculate").click();
    edit.calculated = Boolean(await until(() => document.querySelector('[aria-label^="changes"] tr')));
    Object.assign(edit, { calc_chi2: ctx.calculations.slice(-2).map((c) => c.chi2),
      calc_total: document.querySelector(".tab-panel section.panel strong")?.textContent || "",
      calc_first: document.querySelector('[aria-label^="changes"] tr td')?.textContent || "",
      calc_tab: activeTab() });
    ctx.undoButton.click();
    edit.undone = await until(() => ctx.statusLine.textContent.startsWith("undo: done"));
    Object.assign(edit, { chi2_back: ctx.evaluation?.chi2, rev_back: ctx.project.rev });
  }

  /** "#selftest?edit=1": a pair copied from the first one (its partners, S, radius) typed into the "+ pair" form
   * and added with its button, then removed with "Remove pair", then both undone with the Undo button; records the
   * pairs' count and the χ² at each step in `pair`. */
  async function pairAddRemove() {
    const until = async (test) => {
      for (let i = 0; i < 600 && !test(); i += 1) await new Promise((r) => setTimeout(r, 50));
      return test();
    };
    const visit = async (route) => {
      location.hash = join(route);
      await new Promise((r) => setTimeout(r, 0));
      await draw();
    };
    const pairs = await ctx.api.get("/api/project/pairs");
    const first = pairs[0];
    if (!first || first.kind !== 0) { pair.skipped = "no particle pair to copy"; return; }
    Object.assign(pair, { before: pairs.length, chi2_before: ctx.evaluation?.chi2, rev_before: ctx.project.rev });
    await visit("#/nucleus/pair/+new");
    const type = (label, value) => {
      const input = document.querySelector(`.new-pair [aria-label="${label}"]`);  // a box or a choice
      input.value = String(value);
      input.dispatchEvent(new Event(input.tagName === "SELECT" ? "change" : "input", { bubbles: true }));
    };
    // a copy of the first pair would make the data sets' reactions ambiguous: the heavy partner's charge differs
    type("light J", first.two_j_light / 2); type("light parity", first.light_parity); type("light Z", first.light[0]);
    type("light mass", first.light_mass); type("heavy J", first.two_j_heavy / 2);
    type("heavy parity", first.heavy_parity); type("heavy Z", first.heavy[0] + 1); type("heavy mass", first.heavy_mass);
    type("separation energy (MeV)", first.separation_energy + 5); type("radius (fm)", first.radius);
    document.querySelector('button[data-action="add-pair"]').click();
    pair.added = await until(() => ctx.statusLine.textContent.startsWith("add a pair: done"));
    pair.status_added = ctx.statusLine.textContent;
    const after = await ctx.api.get("/api/project/pairs");
    Object.assign(pair, { after_add: after.length, chi2_added: ctx.evaluation?.chi2 });
    const key = after.find((q) => !pairs.some((o) => o.key === q.key))?.key;
    await visit(`#/nucleus/pair/${key}`);
    document.querySelector('button[data-action="remove-pair"]').click();
    pair.removed = await until(() => ctx.statusLine.textContent.startsWith(`remove pair ${key}: done`));
    pair.after_remove = (await ctx.api.get("/api/project/pairs")).length;
    ctx.undoButton.click();
    await until(() => ctx.statusLine.textContent.startsWith("undo: done"));
    ctx.statusLine.textContent = "";
    ctx.undoButton.click();
    pair.undone = await until(() => ctx.statusLine.textContent.startsWith("undo: done"));
    Object.assign(pair, { after_undo: (await ctx.api.get("/api/project/pairs")).length,
      chi2_back: ctx.evaluation?.chi2, rev_back: ctx.project.rev });
  }

  run().then((ms) => finish(ms), (e) => { errors.push(`self-test: ${e?.stack || e}`); finish(null); });

  function finish(ms) {
    const out = document.createElement("pre");
    out.id = "selftest-result";
    out.textContent = JSON.stringify({ done: true, ms, plotly: typeof globalThis.Plotly === "object",
      errors, failed: ctx.api.failures, steps, edit, pair, endf, preview, fit, chi2: ctx.evaluation?.chi2 ?? null,
      evaluation_refused: ctx.evalError?.what ?? null }, null, 1);
    out.dataset.done = "1";
    out.style.display = "none";
    document.body.appendChild(out);
    document.title = `selftest: ${errors.length || ctx.api.failures.length ? "FAILED" : "done"}`;
  }
}
