// The prior builder (#/prior), opened by "Prior on the parameters…" in the Fit tab: a multivariate Gaussian prior
// on the fit parameters (the <config> key priorCovariance) made from parts.  The parameter table with check boxes
// (the Fit tab's rows: name, kind, unit, x, free or fixed, limits, a prior of its own); the parts as rows, each with
// the fields of its kind as pyazr.priors.PART_KINDS states them (GET /api/prior sends the table); Build (POST
// /api/prior/build: pyazr builds and checks the prior) shows the centres, sigmas, the correlation matrix as a heat
// map, the engine's refusal and the notes; "Write and use" (POST /api/prior/write: the file and the key in one history
// entry, one Undo takes back both); "Remove the prior" (PATCH /api/project/config, the key set to null).  The parts
// are kept in the file the builder writes (its #pyazr.prior line), and the builder reopens with them (#/prior?build=1
// builds them at once); the page computes nothing of the prior.

import { parameterRows } from "../fitting.js";
import { formatNumber } from "../tables.js";
import { button, h, section, themeColours } from "../ui.js";

/** The tab is reached from the Fit tab; it has no link of its own in the header. */
export const shown = () => false;

export async function render(ctx, el) {
  let opened, params;
  try {
    [opened, params] = await Promise.all([ctx.api.get("/api/prior"), ctx.api.get("/api/parameters")]);
  } catch (e) {
    ctx.showError(e, el);
    return;
  }
  const view = ctx.priorView || (ctx.priorView = { parts: opened.parts.map((p) => ({ ...p })),
    path: opened.proposed, checked: new Set(), filter: "", freeOnly: true, summary: null });
  const rows = parameterRows(params);

  el.appendChild(section("Prior on the parameters",
    h("p", {}, "A multivariate Gaussian prior adds (p − μ)ᵀC⁻¹(p − μ) to the objective over the parameters it "
      + "lists, in the fit vector's units (the x column). μ and C are made from the parts below. ",
      h("a", { href: "#/fit" }, "→ back to the Fit tab")),
    h("p.hint", {}, opened.named
      ? `The project names ${opened.named} (${opened.key})` + (opened.missing ? ": the file does not exist."
        : opened.hand_written ? ": written by hand or by prior_from_fit, it has no parts; writing here replaces it."
          : "; the builder opened with its parts.")
      : `The project names no prior file (${opened.key} is not set).`),
    opened.named ? h("div.inline-form", {}, button("Remove the prior", `Remove the key ${opened.key} from <config> `
      + "(the file stays; Undo puts the key back)", () => ctx.edit("remove the prior", () => ctx.api.patch(
      "/api/project/config", { entries: { [opened.key]: null } })).catch(() => {}))) : null));

  el.appendChild(parameterSection(view, rows));
  const partsBox = h("div.prior-parts");
  const drawParts = () => {
    partsBox.textContent = "";
    view.parts.forEach((part, i) => partsBox.appendChild(partRow(opened, view, rows, part, i, drawParts)));
    const kinds = Object.keys(opened.kinds);
    const add = h("select", { "aria-label": "kind of the new part" }, kinds.map((k) =>
      h("option", { value: k, title: opened.kinds[k].doc }, k)));
    partsBox.appendChild(h("div.inline-form", {}, add, button("+ part", "Add a part of the kind chosen; a "
      + "required names field takes the parameters checked above", () => {
      const part = { kind: add.value };
      const names = opened.kinds[add.value].fields.find((f) => f.type === "names" && f.default === null);
      if (names && view.checked.size) part[names.name] = [...view.checked];
      view.parts.push(part);
      drawParts();
    })));
  };
  drawParts();
  const out = h("div.prior-result");
  const build = button("Build", "Build and check the prior (pyazr.priors.build)", async () => {
    out.textContent = "…";
    try {
      view.summary = await ctx.api.post("/api/prior/build", { parts: view.parts });
      out.textContent = "";
      drawSummary(ctx, out, view.summary);
    } catch (e) {
      out.textContent = "";
      ctx.showError(e, out);
    }
  });
  const path = h("input.mono", { type: "text", size: 32, value: view.path, "aria-label": "prior file",
    title: "the file written, relative to the project's directory" });
  path.addEventListener("change", () => { view.path = path.value.trim(); });
  const write = button("Write and use", `Write the file and set ${opened.key} to it (one Undo takes back both)`,
    () => ctx.edit(`prior ${view.path}`, () => ctx.api.post("/api/prior/write",
      { parts: view.parts, path: view.path, rev: ctx.project.rev }), out).catch(() => {}));
  el.appendChild(section("Parts", partsBox, h("div.inline-form", {}, build, h("label", {}, "file ", path), write)));
  el.appendChild(out);
  if (view.summary) drawSummary(ctx, out, view.summary);
  else if (ctx.route.params.build && view.parts.length) build.click();  // #/prior?build=1: a link to the built prior
}

/** The parameter table with a check box per row (the Fit tab's rows), a text filter and "free only". */
function parameterSection(view, rows) {
  const holder = h("div.table-wrap");
  const count = h("span.hint");
  const filter = h("input", { type: "search", value: view.filter, placeholder: "filter", "aria-label":
    "filter the parameters" });
  const freeOnly = h("input", { type: "checkbox", "aria-label": "free parameters only" });
  freeOnly.checked = view.freeOnly;
  const draw = () => {
    const words = view.filter.toLowerCase();
    const shownRows = rows.filter((r) => (!view.freeOnly || r.free) && (!words || r.name.toLowerCase()
      .includes(words) || r.kind.includes(words)));
    count.textContent = `${view.checked.size} checked; ${shownRows.length} of ${rows.length} shown`;
    holder.textContent = "";
    holder.appendChild(h("table.grid", {}, h("caption", {}, "The session's parameters: check those a new part "
      + "takes"), h("thead", {}, h("tr", {}, ["", "name", "kind", "unit", "x", "free", "limits", "own prior"].map(
      (t) => h("th", {}, t)))), h("tbody", {}, shownRows.map((r) => {
      const box = h("input", { type: "checkbox", "aria-label": `check ${r.name}` });
      box.checked = view.checked.has(r.name);
      box.addEventListener("change", () => {
        if (box.checked) view.checked.add(r.name); else view.checked.delete(r.name);
        count.textContent = `${view.checked.size} checked; ${shownRows.length} of ${rows.length} shown`;
      });
      return h("tr", {}, h("td", {}, box), h("td.mono", {}, r.name), h("td", {}, r.kind), h("td", {}, r.unit),
        h("td.num", {}, r.x === null ? (r.physical === null ? "" : formatNumber(r.physical, 8)) : formatNumber(r.x,
          8)), h("td", {}, r.free ? "free" : "fixed"), h("td", {}, r.bounds), h("td", {}, r.prior));
    }))));
  };
  filter.addEventListener("input", () => { view.filter = filter.value; draw(); });
  freeOnly.addEventListener("change", () => { view.freeOnly = freeOnly.checked; draw(); });
  draw();
  return section("Parameters", h("div.inline-form", {}, filter, h("label", {}, freeOnly, " free only"), count,
    button("Clear", "Uncheck every parameter", () => { view.checked.clear(); draw(); })), holder,
  h("p.hint", {}, "x: the fit vector's value of a free parameter (energies in MeV, reduced-width amplitudes in "
    + "MeV^1/2, nuisances in their own units); a fixed parameter shows its value as read, and the engine "
    + "marginalizes it in a prior."));
}

/** One part as a row: its kind, a control per field of the kind (with unit and doc), Remove. */
function partRow(opened, view, rows, part, i, redraw) {
  const kind = opened.kinds[part.kind];
  const row = h("div.prior-part.panel", {}, h("strong", { title: kind?.doc || "" }, `${i + 1}. ${part.kind}`),
    h("span.hint", {}, ` ${kind?.doc || ""}`));
  const fields = h("div.inline-form");
  for (const f of kind?.fields || []) {
    fields.append(h("label", { title: f.doc }, `${f.name}${f.unit && f.unit !== "unit" ? ` (${f.unit})` : ""} `,
      fieldInput(opened, view, rows, part, f, redraw)));
  }
  row.append(fields, h("div.inline-form", {}, button("Remove", "Remove this part", () => {
    view.parts.splice(i, 1);
    redraw();
  })));
  return row;
}

/** The control of one field, by its type in the kinds' table: a file chosen among the project's, the checked names,
 * a number box, pairs chosen from the parameters with a coefficient; a matrix (made in Python) is shown. */
function fieldInput(opened, view, rows, part, f, redraw) {
  const set = (v) => { if (v === "" || v === undefined || (Array.isArray(v) && !v.length)) delete part[f.name];
  else part[f.name] = v; };
  if (f.type.startsWith("path")) {
    const suffixes = (f.type.split(":")[1] || "").split(",").filter(Boolean);
    const files = opened.files.filter((p) => !suffixes.length || suffixes.some((s) => p.endsWith(s)));
    const select = h("select", { "aria-label": f.name }, h("option", { value: "" }, "—"),
      files.map((p) => h("option", { value: p }, p)));
    if (part[f.name] && !files.includes(part[f.name])) select.appendChild(h("option", { value: part[f.name] },
      part[f.name]));
    select.value = part[f.name] || "";
    select.addEventListener("change", () => set(select.value));
    return select;
  }
  if (f.type === "names") {
    const names = part[f.name] || [];
    return h("span", {}, h("span.mono", {}, names.length ? names.join(", ") : (f.default === null ? "(none: "
      + "required)" : "(empty)")), " ", button("take the checked", "Set this field to the parameters checked in "
      + "the table", () => { set([...view.checked]); redraw(); }), names.length ? button("×", "Empty this field",
      () => { set([]); redraw(); }) : null);
  }
  if (f.type === "pairs") {
    const pairs = part[f.name] || [];
    const a = h("select", { "aria-label": "first parameter" }, rows.map((r) => h("option", { value: r.name },
      r.name)));
    const b = h("select", { "aria-label": "second parameter" }, rows.map((r) => h("option", { value: r.name },
      r.name)));
    const rho = h("input", { type: "text", size: 6, "aria-label": "rho", placeholder: "rho" });
    return h("span", {}, pairs.map((p, k) => h("span.mono", {}, `${p[0]} – ${p[1]}: ${p[2]} `,
      button("×", "Remove this pair", () => { pairs.splice(k, 1); set(pairs); redraw(); }), " ")), a, b, rho,
    button("+ pair", "Add the pair with the coefficient typed", () => {
      set([...pairs, [a.value, b.value, rho.value.trim()]]);
      redraw();
    }));
  }
  if (f.type === "matrix") {
    return h("span.hint", {}, part[f.name] ? `${part[f.name].length} × ${part[f.name].length} (made in Python)`
      : "made in Python only; the page adds pairs");
  }
  const input = h("input", { type: "text", size: 10, value: part[f.name] ?? "", "aria-label": f.name, title: f.doc,
    placeholder: f.default === null ? "required" : "" });
  input.addEventListener("change", () => set(input.value.trim()));
  return input;
}

/** The built prior: the refusal or the checks, the notes, the table of centres and sigmas, the correlation. */
function drawSummary(ctx, el, s) {
  if (s.refusal) {
    el.appendChild(h("div.error-box", { role: "alert" }, h("strong", {}, "The engine would refuse it: "),
      s.refusal.what, h("div", {}, `Allowed: ${s.refusal.allowed}. ${s.refusal.fix}.`),
      s.refusal_parts.length ? h("div", {}, `Parts concerned: ${s.refusal_parts.map((p) => `${p + 1} (${s.parts[p]
        .kind})`).join(", ")}`) : null));
  } else {
    const c = s.checks;
    el.appendChild(h("div.ok-box", {}, `Positive definite on ${s.names.length} parameter(s), parametrization `
      + `${s.basis_name}: eigenvalues ${c.smallest_eigenvalue.toPrecision(4)} to `
      + `${c.largest_eigenvalue.toPrecision(4)}, condition number ${c.condition_number.toPrecision(4)}.`));
  }
  if (s.notes.length) el.appendChild(section("What the engine does with it", h("ul", {}, s.notes.map((n) =>
    h("li", {}, n)))));
  el.appendChild(section("Centres and sigmas", h("table.grid", {}, h("thead", {}, h("tr", {}, ["name", "part",
    "unit", "centre", "sigma", ""].map((t) => h("th", {}, t)))), h("tbody", {}, s.names.map((n, i) => h("tr", {},
    h("td.mono", {}, n), h("td", {}, `${s.part_of[i] + 1} ${s.parts[s.part_of[i]].kind}`), h("td", {}, s.units[i]),
    h("td.num", {}, formatNumber(s.mean[i], 10)), h("td.num", {}, formatNumber(s.sigma[i], 4)),
    h("td", {}, s.fixed[i] ? "fixed: marginalized" : "")))))));
  const heat = h("div.plot-area.prior-heat", { role: "img", "aria-label": "correlation of the prior" });
  el.appendChild(section("Correlation", heat));
  const Plotly = globalThis.Plotly;
  if (!Plotly || ctx.route.params.plot === "table") {
    heat.appendChild(h("table.grid", {}, h("thead", {}, h("tr", {}, h("th", {}, ""), s.names.map((n) =>
      h("th.mono", {}, n)))), h("tbody", {}, s.correlation.map((r, i) => h("tr", {}, h("th.mono", {}, s.names[i]),
      r.map((v) => h("td.num", {}, formatNumber(v, 3))))))));
  } else {
    const t = themeColours();
    Plotly.newPlot(heat, [{ type: "heatmap", z: s.correlation, x: s.names, y: s.names, zmin: -1, zmax: 1,
      colorscale: "RdBu", reversescale: true, colorbar: { title: { text: "ρ" } },
      hovertemplate: "%{y} – %{x}: ρ = %{z:.4f}<extra></extra>" }],
    { paper_bgcolor: t.bg, plot_bgcolor: t.bg, font: { color: t.fg }, margin: { l: 140, r: 20, t: 20, b: 120 },
      height: Math.min(700, 260 + 24 * s.names.length), yaxis: { autorange: "reversed", scaleanchor: "x" } },
    { displaylogo: false, responsive: true });
  }
  if (s.text) el.appendChild(h("details", {}, h("summary", {}, "The file's text"), h("pre.mono", {}, s.text)));
}
