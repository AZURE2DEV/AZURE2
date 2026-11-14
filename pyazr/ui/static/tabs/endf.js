// The ENDF tab (#/endf, #/endf/check): ENDF-6 evaluations beside the project (pyazr.endf), in two screens.
//
// Over the data (#/endf?files=a.endf,b.endf): the ENDF-6 files of the project folder (GET /api/endf/files) or a path
// typed inside the project's root; what each shown file holds (POST /api/endf/summary: Material.summary()); per data
// set the files' values at its points in its quantity and unit, drawn over the datum and the model, with the ratios
// to the datum and to the model below (POST /api/endf/curves: Material.curves(data, session=s); the ratios are the
// server's); where a file has nothing for a data set, the refusal's what, allowed and fix in its place.
//
// Check and correct (#/endf/check?file=a.endf): the checking codes found (GET /api/endf/checkers); Check, a job
// (endf_check: pyazr.endf.check): the module's own findings apart from each code's, with the meaning pyazr gives a
// known message, the counts and the raw report; Correct, a job (endf_standardize: pyazr.endf.standardize) on a copy,
// with what the program does to numbers said before the run, the candidate's difference per section, then "Keep as
// a new file" or "Replace" after a second question (POST /api/endf/accept).
//
// Which files are shown is the route's (?files=, ?file=): the page holds no state of the project (web_ui.md section
// 5); the tab remembers the last route per project for this page's life only.  The page computes no value, ratio
// or count: it draws the server's numbers.

import { pathOf } from "../api.js";
import { isFinal, progressView } from "../fitting.js";
import { PALETTE, Y_CHOICES, axisOf, curveTraces, groupsAgainst } from "../plots.js";
import { formatNumber, renderTable } from "../tables.js";
import { button, drawPlot, h, section } from "../ui.js";
import { renderExport, renderNew } from "./endfproject.js";

/** The colour of the i-th file shown (the palette's first colour is the model's). */
export const fileColour = (i) => PALETTE[(i + 1) % PALETTE.length];

/** The files of a route parameter (comma-separated, as parseRoute decoded it; a file name with a comma cannot be
 * shown by its route) and back. */
export const filesOf = (param) => String(param || "").split(",").map((s) => s.trim()).filter(Boolean);
export const filesParam = (files) => files.map(encodeURIComponent).join(",");

/** The route of a screen: over the data with `files`, or the check of `file` (with a job shown). */
export function routeOf({ files = null, file = null, job = null } = {}) {
  if (file !== null) return `#/endf/check?file=${encodeURIComponent(file)}${job ? `&job=${job}` : ""}`;
  return `#/endf${files && files.length ? `?files=${filesParam(files)}` : ""}`;
}

/** Open a screen of this tab, keeping the view's depth where the route gave one (?depth=). */
function go(ctx, hash) {
  const depth = ctx.route.params.depth;
  ctx.navigate(depth ? `${hash}${hash.includes("?") ? "&" : "?"}depth=${depth}` : hash);
}

/**
 * The plots of one data set of the curves' answer: the datum times its normalization with its error, the model, and
 * each file's values (a colour per file, in the order shown), against the data set's axis (energy, or angle at one
 * energy: plots.js axisOf), one group per value of the other coordinate where it repeats; the model and the files on
 * the grid of the data set's model curve where the answer has one (`entry.curve`, each file's `on_curve`), else at
 * the points; below, per file its ratio to the datum (with its error) and to the model, at the points.  `entry`: a
 * data set of the answer; `files`: the names in the order shown.  No number is computed: the points are the
 * server's, sorted by x within a group.
 */
export function evaluationPlots(entry, files, schema) {
  const q = schema?.quantities?.[entry.quantity] || {};
  const data = { e_cm: entry.e_cm, angle_cm: entry.angle_cm, value: entry.value_norm };
  const axis = axisOf(data, q, schema);
  const xs = entry[axis.key];
  const all = [...xs.keys()];
  const groups = axis.key === "angle_cm" ? [{ value: null, index: all }] : groupsAgainst(entry, axis.key, all);
  const pick = (arr, idx) => idx.map((i) => (arr || [])[i]);
  const main = [], ratio = [];
  const drawn = files.filter((f) => entry.curves[f] && !entry.curves[f].refused);
  // the model and each file on the grid of the data set's model curve, where the server sent one along this axis
  const curve = curveTraces(entry, axis.key, Y_CHOICES[0], { color: PALETTE[0] });
  if (curve) {
    main.push(...curve);
    const c = entry.curve;
    for (const f of drawn.filter((name) => entry.curves[name].on_curve?.values)) {
      const values = entry.curves[f].on_curve.values;
      for (const b of c.branches) {
        main.push({ role: "line", name: `${f}${curve.length > 1 ? ` (${b.value})` : ""}`,
          x: c[axis.key].slice(b.start, b.stop), y: values.slice(b.start, b.stop),
          color: fileColour(files.indexOf(f)), dash: "dash" });
      }
    }
  }
  for (const g of groups) {
    const order = [...g.index].sort((a, b) => xs[a] - xs[b]);
    const suffix = g.value === null ? "" : ` (${g.key === "angle_cm" ? "θ_cm" : "E_cm"} = ${g.value})`;
    const x = pick(xs, order);
    main.push({ role: "data", name: `data × norm${suffix}`, x, y: pick(entry.value_norm, order),
      error: pick(entry.error_norm, order), color: null, marker: "circle" });
    if (!curve) main.push({ role: "model", name: `model${suffix}`, x, y: pick(entry.model, order), color: PALETTE[0] });
    for (const f of drawn) {
      const c = entry.curves[f];
      const colour = fileColour(files.indexOf(f));
      if (!curve || !c.on_curve?.values) {
        main.push({ role: "line", name: `${f}${suffix}`, x, y: pick(c.values, order), color: colour, dash: "dash" });
      }
      ratio.push({ role: "data", name: `${f} / data${suffix}`, x, y: pick(c.ratio_to_data, order),
        error: pick(c.ratio_to_data_error, order), color: colour, marker: "circle" });
      ratio.push({ role: "line", name: `${f} / model${suffix}`, x, y: pick(c.ratio_to_model, order), color: colour });
    }
  }
  const title = `${entry.name}: ${q.title || entry.quantity}`;
  const xAxis = { title: axis.title, unit: axis.unit };
  return {
    main: { title, x: xAxis, y: { title: q.title || entry.quantity, unit: entry.engine_unit }, traces: main },
    ratio: drawn.length ? { title: `${entry.name}: the file over the data and over the model`, x: xAxis,
      y: { title: "ratio", unit: "" }, traces: ratio } : null,
  };
}

/** The rows of a check's findings: the module's own (line, MF, MT, message, page) or a code's (MF, MT, record,
 * text, meaning), in the report's order; a finding of no section (MF 0, MT 0: the file as a whole) reads "whole
 * file" in the MF column and nothing in the MT column. */
export function findingRows(findings) {
  return findings.map((f, i) => {
    const whole = !f.mf && !f.mt;
    return { k: i, severity: f.severity, line: f.line || null, mat: f.mat, mf: whole ? "whole file" : f.mf,
      mt: whole ? "" : f.mt, record: f.record || null, message: f.message, page: f.page ? `p. ${f.page}` : "",
      meaning: f.meaning || "" };
  });
}

/** A count line: "2 errors, 1 warning, 0 notes" from the server's counts. */
export function countText(counts) {
  const word = (n, w) => `${n} ${w}${n === 1 ? "" : "s"}`;
  return [word(counts.error, "error"), word(counts.warning, "warning"), word(counts.note, "note")].join(", ");
}

// -- the tab ------------------------------------------------------------------------------------------------------

export async function render(ctx, el, args) {
  const remembered = ctx.endfShown || (ctx.endfShown = {});
  const project = ctx.project?.path || "";
  const params = ctx.route.params;
  // the files of the route; the tab opened from its link shows those of its last route in this page's life
  const files = params.files !== undefined ? filesOf(params.files) : remembered[project] || [];
  remembered[project] = files;
  const listing = await ctx.api.get("/api/endf/files");
  const screen = ["check", "new", "export"].includes(args[0]) ? args[0] : "over";
  const first = params.file || files[0] || listing.files[0]?.path || "";
  const link = (name, href, title) => h("a", { href, class: screen === name ? "active" : null }, title);
  el.appendChild(h("nav.subtabs", { "aria-label": "ENDF screens" },
    link("over", routeOf({ files }), "Over the data"), " · ",
    link("check", routeOf({ file: first }), "Check and correct"), " · ",
    link("new", `#/endf/new${first ? `?file=${encodeURIComponent(first)}` : ""}`, "New project from a file"), " · ",
    link("export", "#/endf/export", "Export this project")));
  if (screen === "check") await renderCheck(ctx, el, listing, params);
  else if (screen === "new") await renderNew(ctx, el, listing, params, { showJob, go });
  else if (screen === "export") await renderExport(ctx, el, params, { showJob, routeOf, go });
  else await renderOver(ctx, el, listing, files);
}

/** The list of the project folder's ENDF-6 files with a box each (shown or not), and a typed path. */
function fileChooser(ctx, listing, files, onChange) {
  const rows = listing.files.map((f) => ({ ...f, shown: files.includes(f.path) }));
  for (const f of files) {
    if (!rows.some((r) => r.path === f)) rows.push({ path: f, size: null, mtime: null, shown: true });
  }
  const table = renderTable(document, {
    columns: [
      { key: "shown", title: "shown", sortable: false, render: (r) => {
        const box = h("input", { type: "checkbox", "aria-label": `show ${r.path}` });
        box.checked = r.shown;
        box.addEventListener("change", () => onChange(box.checked ? [...files, r.path]
          : files.filter((f) => f !== r.path)));
        return box;
      } },
      { key: "colour", title: "", sortable: false, render: (r) => (r.shown
        ? h("span.swatch", { style: `color: ${fileColour(files.indexOf(r.path))}` }, "●") : null) },
      { key: "path", title: "file" },
      { key: "size", title: "size", unit: "bytes", numeric: true },
      { key: "check", title: "", sortable: false, render: (r) => h("span", {}, h("a", { href: routeOf({
        file: r.path }) }, "Check…"), " ", h("a", { href: `#/endf/new?file=${encodeURIComponent(r.path)}`,
        title: "a new project from the file's resonance parameters" }, "New project…")) },
    ],
    rows, rowKey: (r) => r.path, withFilter: rows.length > 8, label: "ENDF-6 files",
    caption: `ENDF-6 files of the project folder (${listing.how})` });
  const typed = h("input.mono", { type: "text", size: 36, "aria-label": "a path inside the project's root",
    placeholder: "or a path inside the project's root" });
  const add = button("Show", "Show this file too (a path inside the project's root)", () => {
    const v = typed.value.trim();
    if (v && !files.includes(v)) onChange([...files, v]);
  });
  typed.addEventListener("keydown", (e) => { if (e.key === "Enter") { e.preventDefault(); add.click(); } });
  return h("div.endf-files", {}, listing.files.length ? table.element
    : h("p.hint", {}, `No ENDF-6 file in the project folder (${listing.how}).`), h("div.inline-form", {}, typed, " ",
    add));
}

async function renderOver(ctx, el, listing, files) {
  const show = (next) => go(ctx, routeOf({ files: next }));
  el.appendChild(section("Evaluations over the data",
    h("p.hint", {}, "An evaluated file's cross sections and distributions at every data set's points, in the data "
      + "set's quantity and unit (pyazr.endf: Material.curves), drawn over the data and the model at the working x. "
      + "Tick the files to show; each has a colour."), fileChooser(ctx, listing, files, show)));
  if (!files.length) {
    el.appendChild(h("p.hint", {}, "No file shown: tick one above, or type its path."));
    return;
  }
  const summaries = h("div.endf-summaries");
  el.appendChild(summaries);
  await Promise.all(files.map(async (f, i) => {
    const box = h("div.endf-summary");
    summaries.appendChild(box);
    try {
      box.appendChild(summaryView(await ctx.api.post("/api/endf/summary", { path: f }), fileColour(i)));
    } catch (e) {
      box.appendChild(h("h3", {}, h("span.swatch", { style: `color: ${fileColour(i)}` }, "● "), f));
      ctx.showError(e, box);
    }
  }));
  const over = h("div.endf-over");
  el.appendChild(section("Each data set", over));
  let answer;
  try {
    answer = await ctx.api.post("/api/endf/curves", { files });
  } catch (e) {
    ctx.showError(e, over);
    return;
  }
  const usePlotly = ctx.route.params.plot !== "table";
  over.appendChild(h("p.hint", {}, `Ratios: ${answer.ratios.ratio_to_data}; ${answer.ratios.ratio_to_model}. `
    + "Each curve says the sections it comes from and the points the file covers."));
  for (const entry of answer.datasets) over.appendChild(datasetView(ctx, entry, answer.files, files, usePlotly));
}

/** What a file holds (Material.summary()): the material, the library, the sections in an evaluator's words, the
 * header text. */
function summaryView(s, colour) {
  const sections = renderTable(document, {
    columns: [{ key: "mf", title: "MF", numeric: true }, { key: "mt", title: "MT", numeric: true },
      { key: "file", title: "file" }, { key: "reaction", title: "reaction" },
      { key: "lines", title: "lines", numeric: true }, { key: "what", title: "what it holds", wrap: true }],
    rows: s.sections, rowKey: (r) => `${r.mf}/${r.mt}`, withFilter: false, label: `sections of ${s.name}` });
  return h("details.endf-file", {},
    h("summary", {}, h("span.swatch", { style: `color: ${colour}` }, "● "), h("strong", {}, s.name),
      ` — MAT ${s.mat}, ${s.target} + ${s.projectile} (${s.sublibrary}); ${s.library}`),
    h("dl.endf-facts", {},
      h("dt", {}, "lab, date, authors"), h("dd", {}, `${s.lab}, ${s.date}, ${s.authors}`),
      h("dt", {}, "LRP"), h("dd", {}, s.lrp),
      h("dt", {}, "EMAX"), h("dd", {}, `${formatNumber(s.emax_mev, 6)} MeV`),
      h("dt", {}, "AWR, AWI"), h("dd", {}, `${s.awr}, ${s.awi} (neutron masses)`)),
    sections.element,
    h("details.more", {}, h("summary", {}, `the header text (${s.text.length} lines)`),
      h("pre.endf-text", {}, s.text.join("\n"))));
}

/** One data set: its plots with each file that gives values, and the refusals of the others in their place. */
function datasetView(ctx, entry, read, files, usePlotly) {
  const box = h("div.endf-dataset", { data: { name: entry.name } }, h("h3", {}, entry.name));
  const unread = read.filter((r) => r.refused);
  const refusals = files.filter((f) => entry.curves[f]?.refused);
  const plots = evaluationPlots(entry, files, ctx.schema);
  if (plots.ratio) {
    const top = h("div.plot");
    const bottom = h("div.plot");
    box.append(top, bottom);
    top.dataset.drawn = drawPlot(top, plots.main, `endf_${entry.name}`, { usePlotly });
    drawPlot(bottom, plots.ratio, `endf_ratio_${entry.name}`, { usePlotly, height: 220 });
    for (const f of files) {
      const c = entry.curves[f];
      if (c && !c.refused) {
        box.appendChild(h("p.hint", {}, h("span.swatch", { style: `color: ${fileColour(files.indexOf(f))}` }, "● "),
          `${f}: ${c.source}, ${c.covered} of ${c.values.length} points covered (${c.unit}, ${c.frame})`
          + (c.notes.length ? `; ${c.notes.join("; ")}` : "")));
      }
    }
  }
  for (const f of refusals) box.appendChild(refusalBox(f, entry.curves[f].refused, fileColour(files.indexOf(f))));
  for (const r of unread) box.appendChild(refusalBox(r.name, r.refused, fileColour(files.indexOf(r.name))));
  return box;
}

/** A file's refusal for a data set (or for every data set: the file is not read): what, allowed, fix. */
function refusalBox(file, parts, colour) {
  return h("div.warn-box.endf-refusal", {}, h("span.swatch", { style: `color: ${colour}` }, "● "),
    h("strong", {}, `${file}: `), parts.what || parts.message,
    parts.allowed ? h("div", {}, h("span.k", {}, "Allowed: "), parts.allowed) : null,
    parts.fix ? h("div", {}, h("span.k", {}, "Fix: "), parts.fix) : null);
}

// -- check and correct ------------------------------------------------------------------------------------------

async function renderCheck(ctx, el, listing, params) {
  const file = params.file || "";
  const [codes, jobs] = await Promise.all([ctx.api.get("/api/endf/checkers"), ctx.api.get("/api/jobs")]);
  const choose = h("select", { "aria-label": "ENDF-6 file" }, [h("option", { value: "" }, "choose a file"),
    ...listing.files.map((f) => h("option", { value: f.path }, f.path))]);
  if (file && !listing.files.some((f) => f.path === file)) choose.appendChild(h("option", { value: file }, file));
  choose.value = file;
  choose.addEventListener("change", () => go(ctx, routeOf({ file: choose.value })));
  const typed = h("input.mono", { type: "text", size: 36, "aria-label": "a path inside the project's root",
    placeholder: "or a path inside the project's root" });
  typed.addEventListener("keydown", (e) => {
    if (e.key === "Enter" && typed.value.trim()) go(ctx, routeOf({ file: typed.value.trim() }));
  });
  el.appendChild(section("Check an ENDF-6 file", h("div.inline-form", {}, h("label", {}, "file ", choose), " ",
    typed), codesView(codes)));
  if (!file) {
    el.appendChild(h("p.hint", {}, "Choose a file to check."));
    return;
  }
  const mine = (kind) => jobs.filter((j) => j.kind === kind && j.settings?.path === file)
    .sort((a, b) => (a.started || 0) - (b.started || 0));
  const shown = (kind) => mine(kind).find((j) => j.id === params.job) || mine(kind).slice(-1)[0] || null;

  // Check
  const found = codes.found.filter((c) => c.kind === "check");
  const boxes = found.map((c) => {
    const box = h("input", { type: "checkbox", "aria-label": `run ${c.name}` });
    box.checked = codes.check_default.includes(c.name);
    return [c.name, box];
  });
  const timeout = h("input.num", { type: "text", size: 5, value: "300", "aria-label": "time limit (s)" });
  const check = button("Check", "Start the check (a job: pyazr.endf.check)", () => startJob(ctx, el, file,
    "endf_check", "codes", { path: file, codes: boxes.filter(([, b]) => b.checked).map(([n]) => n).join(" "),
      timeout: timeout.value.trim() }), "primary");
  const checkPanel = h("div.endf-check-job");
  el.appendChild(section("Check", h("p.hint", {}, "pyazr.endf's own format check always runs (each finding with "
    + "the page of ENDF-102 that states the rule); then each code ticked, on a copy, in a directory of its own."),
  h("div.inline-form", {}, ...boxes.flatMap(([n, b]) => [h("label", {}, b, ` ${n}`), " "]),
    h("label", {}, "time limit ", timeout, " s"), " ", check), checkPanel));
  const checkJob = shown("endf_check");
  if (checkJob) await showJob(ctx, checkPanel, checkJob.id, drawCheck);

  // Correct
  const rewriters = codes.found.filter((c) => c.kind === "rewrite");
  const correctPanel = h("div.endf-correct-job");
  const box = section("Correct", correctPanel);
  el.appendChild(box);
  if (!rewriters.length) {
    box.insertBefore(h("p.hint", {}, `No program that rewrites a file is found (stanef, fixup, dictin, endf2c). `
      + codes.how), correctPanel);
  } else {
    const radios = rewriters.map((c, i) => {
      const r = h("input", { type: "radio", name: "endf-rewriter", value: c.name, "aria-label": c.name });
      r.checked = i === 0;
      return h("div.endf-rewriter", {}, h("label", {}, r, " ", h("strong", {}, c.name), ` — ${c.does}`),
        h("div.hint", {}, `What it does to numbers: ${c.numbers}`));
    });
    const mf2 = h("input", { type: "checkbox", "aria-label": "MF=2 too" });
    const run = button("Run on a copy", "Start the rewrite (a job: pyazr.endf.standardize); nothing is replaced",
      () => startJob(ctx, el, file, "endf_standardize", "program", { path: file,
        program: box.querySelector("input[name=endf-rewriter]:checked")?.value || "stanef",
        mf2: String(mf2.checked), timeout: timeout.value.trim() }));
    box.insertBefore(h("div", {}, h("p.hint", {}, "A program that rewrites runs on a copy; its result is a "
      + "candidate, compared with the file section by section, and replaces nothing until you keep it."), ...radios,
    h("div.inline-form", {}, h("label", {}, mf2, " STANEF: the numbers of MF=2 too"), " ", run)), correctPanel);
  }
  const correctJob = shown("endf_standardize");
  if (correctJob) await showJob(ctx, correctPanel, correctJob.id, drawCandidate);
}

/** The codes pyazr knows: found (path) or not, what each does and what it does to numbers; how they are found. */
function codesView(codes) {
  const found = Object.fromEntries(codes.found.map((c) => [c.name, c]));
  const rows = Object.entries(codes.known).map(([name, k]) => ({ name, ...k, path: found[name]?.path || "" }));
  const table = renderTable(document, {
    columns: [{ key: "name", title: "program" }, { key: "kind", title: "kind" },
      { key: "path", title: "found", render: (r) => h("span", { class: r.path ? "mono" : "hint" },
        r.path || "not found") },
      { key: "does", title: "what it does", wrap: true }, { key: "numbers", title: "what it does to numbers",
        wrap: true }],
    rows, rowKey: (r) => r.name, withFilter: false, label: "checking codes" });
  return h("details.endf-codes", { open: codes.found.length ? null : true },
    h("summary", {}, codes.found.length ? `${codes.found.length} checking codes found: `
      + codes.found.map((c) => c.name).join(", ") : "No checking code found: the module's own format check runs "
      + "alone"), h("p.hint", {}, codes.how), table.element);
}

async function startJob(ctx, el, file, kind, method, settings) {
  try {
    const started = await ctx.api.post("/api/jobs", { kind, method, settings });
    ctx.status(`${kind} job ${started.id} started on ${file}`, "ok");
    go(ctx, routeOf({ file, job: started.id }));
  } catch (e) {
    ctx.showError(e, el);
  }
}

/** A job of this tab: its progress while it runs (one long poll in flight), then `draw(ctx, panel, job)`. */
async function showJob(ctx, panel, id, draw) {
  let job;
  try {
    job = await ctx.api.get(pathOf("/api/jobs", id));
  } catch (e) {
    ctx.showError(e, panel);
    return;
  }
  const paint = (j) => {
    panel.textContent = "";
    if (!isFinal(j.state)) {
      panel.appendChild(h("p.endf-running", {}, `job ${j.id}: ${progressView(j).words}…`));
      return;
    }
    if (j.state !== "finished" || !j.result) {
      panel.appendChild(h("div.warn-box", { role: "alert" }, `job ${j.id} ${j.state}: ${j.result?.message
        || "no result"}`));
      return;
    }
    draw(ctx, panel, j);
  };
  paint(job);
  while (panel.isConnected && !isFinal(job.state)) {
    try {
      job = await ctx.api.get(`${pathOf("/api/jobs", id)}?wait=1`, { background: true });
    } catch (e) {
      ctx.showError(e, panel);
      return;
    }
    if (panel.isConnected) paint(job);
    if (!isFinal(job.state)) await new Promise((done) => setTimeout(done, 200));
  }
}

const findingColumns = (own) => [
  { key: "severity", title: "severity", render: (r) => h(`span.severity.${r.severity}`, {}, r.severity) },
  ...(own ? [{ key: "line", title: "line", numeric: true }] : []),
  { key: "mf", title: "MF", doc: "the file (MF) of the finding; whole file: a finding of no section" },
  { key: "mt", title: "MT", numeric: true },
  ...(own ? [] : [{ key: "record", title: "record", numeric: true }]),
  { key: "message", title: own ? "message" : "the code's text", wrap: true },
  own ? { key: "page", title: "ENDF-102" } : { key: "meaning", title: "what it means (pyazr.endf)", wrap: true },
];

/** A finished check: the counts, the own findings, then per code its findings, the counts its report states and
 * the raw report one click away. */
function drawCheck(ctx, panel, job) {
  const c = job.result.check;
  panel.appendChild(h("p.endf-check-summary", { data: { job: job.id } }, h("strong", {}, `${c.file}: `),
    job.result.message, ` (job ${job.id})`));
  const own = c.own;
  panel.appendChild(h("h3", {}, `pyazr.endf's own format check: ${countText(own.counts)}`));
  panel.appendChild(own.findings.length ? renderTable(document, { columns: findingColumns(true),
    rows: findingRows(own.findings), rowKey: (r) => r.k, withFilter: false, label: "own findings" }).element
    : h("p.hint", {}, "No finding."));
  if (!c.reports.length) {
    panel.appendChild(h("p.hint", {}, "No checking code was run (none found or none ticked)."));
  }
  for (const r of c.reports) {
    const stated = r.totals ? `; the report counts ${r.totals[0]} errors, ${r.totals[1]} warnings` : "";
    panel.appendChild(h("h3", {}, `${r.program.toUpperCase()}${r.version ? ` ${r.version}` : ""} (exit `
      + `${r.exit_status ?? "none: stopped at the time limit"}): ${countText(r.counts)} read${stated}`));
    panel.appendChild(r.findings.length ? renderTable(document, { columns: findingColumns(false),
      rows: findingRows(r.findings), rowKey: (x) => x.k, withFilter: false, label: `${r.program} findings` }).element
      : h("p.hint", {}, "No finding."));
    panel.appendChild(h("details.endf-report", {}, h("summary", {}, `the raw report of ${r.program} (${r.path})`),
      h("pre.endf-text", {}, r.report || r.terminal)));
  }
}

/** A finished rewrite: what the program does to numbers, the candidate's difference per section, and the two ways
 * to keep it (a new file; the original replaced after a second question). */
function drawCandidate(ctx, panel, job) {
  const c = job.result.candidate;
  const rows = c.diff.map((d) => ({ ...d, lines: `${d.lines_before} → ${d.lines_after}` }));
  panel.appendChild(h("p", {}, h("strong", {}, `${c.program} on ${c.file}: `), job.result.message,
    ` (job ${job.id}). What it does to numbers: ${c.numbers}`));
  panel.appendChild(renderTable(document, {
    columns: [{ key: "mf", title: "MF", numeric: true }, { key: "mt", title: "MT", numeric: true },
      { key: "kind", title: "change" }, { key: "numbers_changed", title: "numbers changed", numeric: true },
      { key: "max_relative_change", title: "largest relative change", numeric: true, digits: 3 },
      { key: "digits_lost", title: "digits lost", numeric: true }, { key: "lines", title: "lines" },
      { key: "example", title: "example", wrap: true }],
    rows, rowKey: (r) => `${r.mf}/${r.mt}`, withFilter: false, label: "candidate's difference",
    caption: "The candidate against the file, section by section (Candidate.diff())" }).element);
  panel.appendChild(h("details.endf-report", {}, h("summary", {}, "the program's terminal text and listing"),
    h("pre.endf-text", {}, c.report || "(none)")));
  const kept = job.kept || [];
  for (const k of kept) panel.appendChild(h("p.ok-box", {}, `kept as ${k.path}${k.replaced ? " (the original "
    + "replaced)" : ""}`));
  const name = h("input.mono", { type: "text", size: 36, "aria-label": "new file",
    placeholder: `beside ${c.file}: pyazr's name for it` });
  const keep = button("Keep as a new file", "Write the candidate to a new file (Candidate.accept(to=...)); an "
    + "existing file is not replaced", () => accept(ctx, panel, job, { to: name.value.trim() || null }));
  const confirmBox = h("span.endf-confirm");
  const replace = button(`Replace ${c.file}…`, "Replace the original with the candidate (asks again)", () => {
    confirmBox.textContent = "";
    confirmBox.append(h("span.warn-box", {}, `Replace ${c.file} with the candidate? The original is not kept.`),
      " ", button("Yes, replace", "Candidate.accept(to=the file, overwrite=True)", () => accept(ctx, panel, job,
        { to: c.file, overwrite: true }), "danger"), " ", button("Cancel", "Keep the original",
        () => { confirmBox.textContent = ""; }));
  });
  panel.appendChild(h("div.inline-form.endf-keep", {}, name, " ", keep, " ", replace, " ", confirmBox));
}

async function accept(ctx, panel, job, body) {
  try {
    const answer = await ctx.api.post("/api/endf/accept", { job_id: job.id, ...body });
    ctx.status(`the candidate of job ${job.id} kept as ${answer.path}${answer.replaced ? " (the original replaced)"
      : ""}`, "ok");
    go(ctx, routeOf({ file: ctx.route.params.file, job: job.id }));
  } catch (e) {
    ctx.showError(e, panel);
  }
}
