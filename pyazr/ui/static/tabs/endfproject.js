// The ENDF tab's two ways between the project and an evaluation (#/endf/new, #/endf/export), beside its screens over
// the data and of the check (tabs/endf.js).
//
// New project from this file (#/endf/new?file=a.endf[&job=id]): what MF=2 of the file holds, nothing written
// (POST /api/endf/resonances: Material.resonances()), or the refusal's three parts; Convert, a job (endf_import:
// pyazr.endf.to_project in a directory of the job's own): the report of every conversion, the project's text; Write
// the new project, a new file of the folder (POST /api/endf/keep_project); then Open it here (POST
// /api/project/open: the page serves one project at a time) or `pyazr ui <file>` for a page of its own.
//
// Export this project (#/endf/export[?job=id]): what a file of the project at the working x would lose, before it is
// written (GET /api/endf/export: pyazr.endf.from_session), the MF=1 fields, the range and the grid (the job row's
// settings), Write, a job (endf_write: pyazr.endf.write): the file's sections, grids, the MF=2 check, the losses,
// what is not written, and Check on the file.  The page computes nothing: it shows the server's answers.

import { renderTable } from "../tables.js";
import { button, h, section } from "../ui.js";
import { SETTING_WIDGETS } from "./fitjob.js";

/** The rows of a new project's report (Conversion as the server sent it): pairs and levels. */
export function reportRows(report) {
  const pairs = (report.pairs || []).map((p) => ({ ...p }));
  const levels = (report.levels || []).map((l, i) => ({ k: i, jpi: l.jpi, er_eV: l.er_eV, formal_MeV: l.formal_MeV,
    energy_MeV: l.energy_MeV, status: l.status, background: l.background ? "above EH: background" : "",
    values: l.values.map((v, j) => `${v} ${l.kinds[j] === "width" ? "eV" : "fm^-1/2"}`).join(", ") }));
  return { pairs, levels };
}

/** A refusal's three parts in a box. */
export function refusalView(parts, title = "") {
  return h("div.warn-box.endf-refusal", { role: "note" }, title ? h("strong", {}, `${title}: `) : null,
    parts.what || parts.message, parts.allowed ? h("div", {}, h("span.k", {}, "Allowed: "), parts.allowed) : null,
    parts.fix ? h("div", {}, h("span.k", {}, "Fix: "), parts.fix) : null);
}

/** The job's settings as widgets (the job table's row: GET /api/jobs/kinds), with values given (`values`): the
 * setting's words (its doc) beside its box, its key on hover. */
function settingWidgets(row, values) {
  return row.settings.map((s) => {
    const w = SETTING_WIDGETS[s.type](s, values[s.name]);
    return { s, w, el: [h("label.setting-words", { title: s.name }, s.doc), w.el] };
  });
}

async function startJob(ctx, el, kind, method, settings, route, go) {
  try {
    const started = await ctx.api.post("/api/jobs", { kind, method, settings });
    ctx.status(`${kind} job ${started.id} started`, "ok");
    go(ctx, route(started.id));
  } catch (e) {
    ctx.showError(e, el);
  }
}

const shownJob = (jobs, kind, id, match = () => true) => {
  const mine = jobs.filter((j) => j.kind === kind && match(j)).sort((a, b) => (a.started || 0) - (b.started || 0));
  return mine.find((j) => j.id === id) || mine.slice(-1)[0] || null;
};

// -- a new project from a file --------------------------------------------------------------------------------------

export async function renderNew(ctx, el, listing, params, { showJob, go }) {
  const file = params.file || "";
  const route = (more = {}) => `#/endf/new?${new URLSearchParams({ file, ...more })}`;
  const choose = h("select", { "aria-label": "ENDF-6 file" }, [h("option", { value: "" }, "choose a file"),
    ...listing.files.map((f) => h("option", { value: f.path }, f.path))]);
  if (file && !listing.files.some((f) => f.path === file)) choose.appendChild(h("option", { value: file }, file));
  choose.value = file;
  choose.addEventListener("change", () => go(ctx, `#/endf/new?file=${encodeURIComponent(choose.value)}`));
  el.appendChild(section("New project from an ENDF-6 file", h("p.hint", {}, "The file's resonance parameters "
    + "(MF=2, LRF=7, R-Matrix Limited) as a new project in the Brune basis, every level and width fixed "
    + "(pyazr.endf.to_project). What the file holds and the report of every conversion are shown before anything "
    + "is written; the new project is a new file of the folder, never the open one."),
  h("div.inline-form", {}, h("label", {}, "file ", choose))));
  if (!file) {
    el.appendChild(h("p.hint", {}, "Choose a file."));
    return;
  }
  const [table, jobs] = await Promise.all([ctx.api.get("/api/jobs/kinds"), ctx.api.get("/api/jobs")]);
  const row = table.find((r) => r.kind === "endf_import");
  const job = shownJob(jobs, "endf_import", params.job, (j) => j.settings?.path === file);
  const values = { ...(job?.settings || {}), path: file, isotope: params.isotope ?? job?.settings?.isotope ?? 0,
    energy_range: params.range ?? job?.settings?.energy_range ?? 0 };
  const holds = h("div.endf-holds", {}, "…");
  el.appendChild(section("What MF=2 holds (nothing written)", holds));
  const widgets = settingWidgets(row, values).filter((x) => x.s.name !== "path");
  // another isotope or range of MF=2: what it holds is shown again (the route keeps the choice)
  const read = (name) => widgets.find((x) => x.s.name === name).w.read();
  for (const { w } of widgets.filter((x) => ["isotope", "energy_range"].includes(x.s.name))) {
    w.el.addEventListener("change", () => go(ctx, route({ isotope: read("isotope"), range: read("energy_range") })));
  }
  const convert = button("Convert", "Make the new project in a directory of the job's own (a job: "
    + "pyazr.endf.to_project); nothing is written into the folder", () => startJob(ctx, el, "endf_import",
    "to_project", { path: file, ...Object.fromEntries(widgets.map(({ s, w }) => [s.name, w.read()])) },
    (id) => route({ job: id }), go), "primary");
  el.appendChild(section("Convert", h("div.endf-settings", {}, ...widgets.flatMap((x) => x.el)),
    h("div.form-actions", {}, convert)));
  const result = h("div.endf-import-job", {}, h("p.hint", {}, "Convert shows the report of every conversion here; "
    + "nothing is written into the folder until you write the new project."));
  el.appendChild(section("The new project", result));
  try {
    const r = await ctx.api.post("/api/endf/resonances", { path: file, isotope: Number(values.isotope) || 0,
      energy_range: Number(values.energy_range) || 0 }, { quiet: true });
    holds.textContent = "";
    if (r.refused) {
      holds.appendChild(refusalView(r.refused, "Not convertible"));
      convert.disabled = true;
    } else {
      holds.append(h("p", {}, `EL ${r.el_eV} eV – EH ${r.eh_eV} eV; incident pair ${r.entrance}; `
        + `${r.groups.length} spin groups, ${r.groups.reduce((n, g) => n + g.levels, 0)} levels.`),
      renderTable(document, { rows: r.pairs, rowKey: (p) => String(p.pair), withFilter: false, label: "pairs",
        columns: [{ key: "pair", title: "pair", numeric: true }, { key: "mt", title: "MT", numeric: true },
          { key: "reaction", title: "reaction" }, { key: "q_eV", title: "Q", unit: "eV", numeric: true }] }).element,
      h("details.more", {}, h("summary", {}, "the record as read"), h("pre.endf-text", {}, r.text)));
    }
  } catch (e) {
    ctx.showError(e, holds);
  }
  if (job) await showJob(ctx, result, job.id, drawImport);
}

function drawImport(ctx, panel, job) {
  const p = job.result.project;
  const { pairs, levels } = reportRows(p.report);
  panel.append(h("h3", {}, `The new project ${p.name} (job ${job.id}): ${job.result.message}`),
    h("p", {}, `Separation energy of the incident pair ${p.report.separation_MeV} MeV (${p.report.separation_source})`),
    h("details.more", { open: true }, h("summary", {}, "the report as pyazr prints it (the masses against AME2020 in "
      + "keV)"), h("pre.endf-text", {}, p.report_text)),
    renderTable(document, { rows: pairs, rowKey: (r) => String(r.key), withFilter: false, label: "new pairs",
      caption: "The pairs written", columns: [{ key: "key", title: "key", numeric: true },
        { key: "mt", title: "MT", numeric: true }, { key: "q_eV", title: "Q", unit: "eV", numeric: true },
        { key: "separation_MeV", title: "separation energy", unit: "MeV", numeric: true },
        { key: "excitation_MeV", title: "excitation", unit: "MeV", numeric: true }] }).element,
    renderTable(document, { rows: levels, rowKey: (r) => r.k, withFilter: levels.length > 8, label: "new levels",
      caption: "The levels written (fixed): the file's ER, the formal and the Brune energy, the observed values",
      columns: [{ key: "jpi", title: "Jπ" }, { key: "er_eV", title: "ER", unit: "eV, lab", numeric: true },
        { key: "formal_MeV", title: "formal E_x", unit: "MeV", numeric: true },
        { key: "energy_MeV", title: "E_x", unit: "MeV", numeric: true }, { key: "status", title: "Brune iteration" },
        { key: "values", title: "widths and ANCs", wrap: true }, { key: "background", title: "" }] }).element);
  for (const s of p.report.skipped) panel.appendChild(h("p.warn-box", {}, `Not written: ${s}`));
  for (const n of p.report.notes) panel.appendChild(h("p.hint", {}, `Note: ${n}`));
  panel.appendChild(h("details.more", {}, h("summary", {}, "the new project's text"), h("pre.endf-text", {}, p.text)));
  const kept = job.kept || [];
  const name = h("input.mono", { type: "text", size: 30, value: p.name, "aria-label": "the new project's file" });
  const keep = button("Write the new project", "A new file of the project folder (POST /api/endf/keep_project); "
    + "an existing file is never replaced", async () => {
    try {
      const answer = await ctx.api.post("/api/endf/keep_project", { job_id: job.id, to: name.value.trim() || null });
      ctx.status(`the new project is ${answer.path}`, "ok");
      ctx.navigate(location.hash);
    } catch (e) {
      ctx.showError(e, panel);
    }
  }, "primary");
  panel.appendChild(h("div.inline-form.endf-keep", {}, name, " ", keep));
  for (const k of kept) panel.appendChild(openOffer(ctx, k));
}

/** A kept project: where it is, and how to open it (this page serves one project: Open it here replaces the open
 * project; a second page is `pyazr ui <file>`). */
function openOffer(ctx, kept) {
  const unsaved = ctx.project?.unsaved;
  const open = button(`Open ${kept.path} here`, unsaved ? "The open project has unsaved edits: save them first"
    : "This page then serves the new project (POST /api/project/open); the open one stays as saved", async () => {
    try {
      await ctx.api.post("/api/project/open", { path: kept.open_path });
      location.hash = "#/project";
      location.reload();
    } catch (e) {
      ctx.showError(e);
    }
  });
  open.disabled = Boolean(unsaved);
  return h("div.ok-box.endf-kept", {}, h("strong", {}, `Written: ${kept.path}. `), open,
    unsaved ? h("span.hint", {}, " (save the open project first: it has unsaved edits)") : null,
    h("div.hint", {}, "Or beside this page, from a terminal in the project's folder: ", h("code", {},
      `pyazr ui ${kept.path}`)));
}

// -- export ---------------------------------------------------------------------------------------------------------

export async function renderExport(ctx, el, params, { showJob, routeOf, go }) {
  const [table, jobs] = await Promise.all([ctx.api.get("/api/jobs/kinds"), ctx.api.get("/api/jobs")]);
  const row = table.find((r) => r.kind === "endf_write");
  const job = shownJob(jobs, "endf_write", params.job);
  const head = section("Export this project as an ENDF-6 file", h("p.hint", {}, "The model at the working x as "
    + "an ENDF-6 material (pyazr.endf.write): MF=1, MF=2 (LRF=7, the particle channels in the formal basis), MF=3 of "
    + "every exit, of capture and of elastic scattering on a grid refined to the tolerance, MF=6 (the angular "
    + "distributions on the same energies), a new file of the project folder. Writing takes a minute or more. "
    + "What it loses is said before it is written and in the result; then Check the file."));
  el.appendChild(head);
  let preview;
  try {
    preview = await ctx.api.get("/api/endf/export", { quiet: true });
  } catch (e) {
    ctx.showError(e, head);
    return;
  }
  const lost = h("div.endf-losses");
  if (preview.refused) {
    lost.appendChild(refusalView(preview.refused, "Not written"));
  } else {
    lost.append(h("h3", {}, "What the file will not hold"),
      h("ul", {}, ...preview.losses.map((t) => h("li", {}, t)), ...preview.not_written.map((t) => h("li", {}, t))),
      ...(preview.photon_widths.length ? [renderTable(document, { rows: preview.photon_widths.map((w, i) => ({ k: i,
        ...w })), rowKey: (r) => r.k, withFilter: false, label: "photon widths",
      caption: "The photon widths, listed in the MF=1 text (MF=2 holds the particle channels)",
      columns: [{ key: "jpi", title: "Jπ" }, { key: "energy_MeV", title: "E_x", unit: "MeV", numeric: true },
        { key: "pair", title: "to pair", numeric: true }, { key: "multipole", title: "multipole" },
        { key: "width_eV", title: "Γ_γ", unit: "eV", numeric: true, digits: 6 }] }).element] : []));
  }
  el.appendChild(section("Before writing", lost));
  const range = preview.energy_range_MeV || [];
  const widgets = settingWidgets(row, { ...(job?.settings || {}), name: preview.name });
  for (const { s, w } of widgets) {
    if (s.name === "e_low" && range.length) w.el.placeholder = `${range[0]} (the data's)`;
    if (s.name === "e_high" && range.length) w.el.placeholder = `${range[1]} (the data's)`;
  }
  const write = button("Write", "Write the file (a job: pyazr.endf.write); an existing file is never replaced",
    () => startJob(ctx, el, "endf_write", "writer", Object.fromEntries(widgets.map(({ s, w }) => [s.name,
      w.read()])), (id) => `#/endf/export?job=${id}`, go), "primary");
  write.disabled = Boolean(preview.refused);
  el.appendChild(section("Write", h("div.endf-settings", {}, ...widgets.flatMap((x) => x.el)),
    h("div.form-actions", {}, write)));
  const result = h("div.endf-write-job", {}, h("p.hint", {}, "Write shows the file's sections and what it does "
    + "not hold here."));
  el.appendChild(section("The file", result));
  if (job) await showJob(ctx, result, job.id, (c, panel, j) => drawWritten(c, panel, j, routeOf));
}

function drawWritten(ctx, panel, job, routeOf) {
  const w = job.result.written;
  const rows = w.sections.map(([mf, mt, lines]) => ({ mf, mt, lines, grid: w.grids[mt] ? `${w.grids[mt][0]} points, `
    + `${w.grids[mt][1]} – ${w.grids[mt][2]} eV` : "", check: w.mf2_check[mt] ?? "" }));
  panel.append(h("h3", {}, `${w.path} (job ${job.id})`, " ", h("a", { href: routeOf({ file: w.path }) }, "Check…")),
    renderTable(document, { rows, rowKey: (r) => `${r.mf}/${r.mt}`, withFilter: false, label: "sections written",
      columns: [{ key: "mf", title: "MF", numeric: true }, { key: "mt", title: "MT", numeric: true },
        { key: "lines", title: "lines", numeric: true }, { key: "grid", title: "grid of MF=3" },
        { key: "check", title: "MF=2 reproduces MF=3 to", numeric: true, digits: 2 }] }).element,
    h("h3", {}, "Not in the file"), h("ul", {}, ...w.losses.map((t) => h("li", {}, t)),
      ...w.not_written.map((t) => h("li", {}, t))));
}
