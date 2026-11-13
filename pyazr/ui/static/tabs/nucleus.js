// The Nucleus tab (docs/upgrade/web_ui.md section 6.2): the particle pairs (edited through set_pair), the levels
// with their channels (a width, an energy, a fixed flag edited through the project's routes: each edit is one entry
// of the history that Undo takes back, and the header shows the χ² the engine computes after it), and the records of
// the sections NUCLEUS_SECTIONS (γ decays, potentials) through the generated forms.
//
// Adding and removing pairs waits for pyazr (Project.add_pair, remove_pair): the page has no button for it.

import { pathOf } from "../api.js";
import { formatNumber, renderTable } from "../tables.js";
import { button, h, section } from "../ui.js";
import { NUCLEUS_SECTIONS, newRecordPanel, recordPanel, recordRoute } from "./records.js";

// The editable numbers of a pair, with the units the route documents (PATCH /api/project/pairs/{key}).
const PAIR_FIELDS = [["radius", "fm"], ["light_mass", "u"], ["heavy_mass", "u"], ["separation_energy", "MeV"],
  ["light_excitation", "MeV"], ["heavy_excitation", "MeV"]];

export async function render(ctx, el, args) {
  const [pairs, levels, records, problems] = await Promise.all([ctx.api.get("/api/project/pairs"),
    ctx.api.get("/api/project/levels"), ctx.api.get("/api/project/records"), ctx.api.get("/api/project/validate")]);
  ctx.pairs = pairs;
  const pairLabel = Object.fromEntries(pairs.map((p) => [p.key, p.label]));
  const msg = h("div.messages");
  el.appendChild(msg);

  // -- pairs
  const pairRows = pairs.map((p) => ({ ...p, j_light: p.two_j_light, j_heavy: p.two_j_heavy }));
  const pairTable = renderTable(document, {
    rowKey: (r) => String(r.key), label: "pairs", withFilter: false,
    caption: "Particle pairs (adding and removing pairs waits for pyazr: Project.add_pair, remove_pair)",
    columns: [{ key: "key", title: "key", numeric: true }, { key: "label", title: "pair" },
      { key: "kind", title: "kind", numeric: true, doc: "0 particle, 10 photon, 20 beta (the engine's pair type)" },
      { key: "j_light", title: "2J₁", numeric: true }, { key: "j_heavy", title: "2J₂", numeric: true },
      { key: "light_mass", title: "M₁", unit: "u", numeric: true }, { key: "heavy_mass", title: "M₂", unit: "u",
        numeric: true }, { key: "separation_energy", title: "S", unit: "MeV", numeric: true },
      { key: "heavy_excitation", title: "E_x2", unit: "MeV", numeric: true },
      { key: "radius", title: "a", unit: "fm", numeric: true }],
    rows: pairRows, onSelect: (r) => showPair(r),
  });
  const pairEdit = h("div.inline-form");
  el.appendChild(section("Pairs", pairTable.element, pairEdit));

  function showPair(p) {
    pairEdit.textContent = "";
    const inputs = PAIR_FIELDS.map(([k, unit]) => {
      const input = h("input.num", { type: "text", value: String(p[k]), size: 12, "aria-label": `${k} (${unit})` });
      pairEdit.appendChild(h("label.ctl", {}, h("span.ctl-label", {}, `${k.replace(/_/g, " ")} (${unit})`), input));
      return [k, input];
    });
    pairEdit.appendChild(button(`Apply to pair ${p.key}`, "set_pair with the changed values", () => {
      const changed = Object.fromEntries(inputs.filter(([k, i]) => i.value.trim() !== String(p[k]))
        .map(([k, i]) => [k, Number(i.value)]));
      if (!Object.keys(changed).length) return;
      ctx.edit(`set pair ${p.key}`, () => ctx.api.patch(pathOf("/api/project/pairs", p.key), changed), msg)
        .catch(() => {});
    }, "primary"));
  }

  // -- levels
  const selectedLevel = args[0] === "level" ? levels.find((l) => l.label === args[1]) : null;
  const levelRows = levels.map((l) => ({ ...l, n: l.channels.length,
    on: l.active ? "on" : "off", efix: l.fixed ? "fixed" : "free" }));
  const levelTable = renderTable(document, {
    rowKey: (r) => r.label, label: "levels", selected: selectedLevel?.label, caption: "Levels of <levels>",
    columns: [{ key: "jpi", title: "Jπ" }, { key: "energy", title: "E_λ", unit: "MeV", numeric: true,
      value: (r) => r.energy, format: (v, r) => r.energy_text },
    { key: "efix", title: "energy" }, { key: "on", title: "active" }, { key: "n", title: "channels", numeric: true },
    { key: "nucleus", title: "nucleus", format: (v) => (Array.isArray(v) ? `Z ${v[0]}, A ${v[1]}` : "") },
    { key: "label", title: "label" }],
    rows: levelRows,
    onSelect: (r) => ctx.navigate(`#/nucleus/level/${encodeURIComponent(r.label)}`),
  });
  el.appendChild(section("Levels", levelTable.element));
  if (selectedLevel) el.appendChild(levelPanel(ctx, selectedLevel, pairLabel, msg));

  // -- the records of the sections the Nucleus tab claims (γ decays, potentials)
  const kinds = NUCLEUS_SECTIONS.flatMap((s) => ctx.schema.sections[s] || []);
  const mine = records.filter((r) => kinds.includes(r.kind));
  const recTable = renderTable(document, {
    rowKey: (r) => r.name, label: "records", selected: args[0], withFilter: false,
    caption: `Records of ${NUCLEUS_SECTIONS.join(", ")}: ${kinds.join(", ")}`,
    columns: [{ key: "name", title: "name" }, { key: "kind", title: "kind" }, { key: "line", title: "line",
      numeric: true }],
    rows: mine, onSelect: (r) => ctx.navigate(recordRoute(ctx.schema, r)),
  });
  el.appendChild(section("γ decays and potentials", recTable.element,
    button("+ new", `Add a record of ${kinds.join(", ")}`, () => ctx.navigate("#/nucleus/+new"))));
  if (args[0] === "+new") {
    el.appendChild(await newRecordPanel(ctx, kinds, records,
      { onAdded: (r) => ctx.navigate(recordRoute(ctx.schema, r)) }));
  }
  const record = args[0] && args[0] !== "level" ? mine.find((r) => r.name === args[0]) : null;
  if (record) {
    el.appendChild(await recordPanel(ctx, record, records, problems, { onGone: () => ctx.navigate("#/nucleus") }));
  }
}

/** One level: its energy, fixed and active flags, and its channels with width and fixed flag, each editable. */
function levelPanel(ctx, level, pairLabel, msg) {
  const path = pathOf("/api/project/levels", level.label);
  const box = h("section.panel.level", { "aria-label": `level ${level.label}` });
  const energy = h("input.num", { type: "text", value: level.energy_text, size: 12, "aria-label": "energy (MeV)" });
  const fixed = h("input", { type: "checkbox", "aria-label": "energy fixed" });
  fixed.checked = level.fixed;
  const active = h("input", { type: "checkbox", "aria-label": "active" });
  active.checked = level.active;
  const applyLevel = () => {
    const changed = {};
    if (energy.value.trim() !== level.energy_text) changed.energy = Number(energy.value);
    if (fixed.checked !== level.fixed) changed.fixed = fixed.checked;
    if (active.checked !== level.active) changed.active = active.checked;
    if (Object.keys(changed).length) chi2Edit(ctx, `set level ${level.label}`, () => ctx.api.patch(path, changed), msg);
  };
  energy.addEventListener("keydown", (e) => { if (e.key === "Enter") applyLevel(); });
  box.append(h("h2", {}, `Level ${level.jpi} at ${level.energy_text} MeV`),
    h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "E_λ (MeV)"), energy),
      h("label.ctl", {}, fixed, h("span.ctl-label", {}, "energy fixed")),
      h("label.ctl", {}, active, h("span.ctl-label", {}, "active")),
      button("Apply", "set_level with the changed values (Enter in the energy box)", applyLevel, "primary")));

  const rows = level.channels.map((c) => ({ ...c, pair_label: pairLabel[c.pair] ?? c.pair,
    key: `${c.pair}/${c.l}/${c.s}` }));
  const columns = [
    { key: "channel", title: "ch", numeric: true }, { key: "pair_label", title: "pair" },
    { key: "l", title: "l", numeric: true }, { key: "s", title: "s", numeric: true },
    { key: "width", title: "width as written", numeric: true, sortable: false,
      doc: "the row's convention: a partial width Γ (eV), an ANC (fm^-1/2) or an amplitude (MeV^1/2)",
      render: (r) => widthEditor(r) },
    { key: "fixed", title: "fixed", sortable: false, render: (r) => fixedEditor(r) },
    { key: "amplitude", title: "amplitude", format: (v) => (v ? "yes" : "") },
    { key: "line", title: "line", numeric: true },
  ];
  const table = renderTable(document, { columns, rows, rowKey: (r) => r.key, withFilter: false, label: "channels",
    caption: "Channels: Enter in a width box sends the width; the header shows the new χ²" });
  box.appendChild(table.element);
  box.appendChild(h("p.hint", {}, "A width is written as the file has it: Γ in eV for a particle channel above "
    + "threshold, an ANC in fm^-1/2 below it, an amplitude in MeV^1/2 where the row says so. The engine's "
    + "physical values and Wigner limits come with the Fit tab."));
  return box;

  function channelPath(r) {
    return `${path}/channels/${r.pair}/${r.l}/${r.s}`;
  }
  function widthEditor(r) {
    const input = h("input.num", { type: "text", value: r.width_text, size: 14,
      "aria-label": `width of channel ${r.channel}` });
    let last = r.width_text;  // Enter and the change that follows it send once
    const send = () => {
      if (input.value.trim() === last) return;
      last = input.value.trim();
      chi2Edit(ctx, `width of ${level.label} channel ${r.channel}`,
        () => ctx.api.patch(channelPath(r), { width: Number(input.value) }), msg);
    };
    input.addEventListener("keydown", (e) => {
      if (e.key === "Enter") { e.preventDefault(); send(); } else if (e.key === "Escape") input.value = r.width_text;
    });
    input.addEventListener("change", send);
    return input;
  }
  function fixedEditor(r) {
    const box2 = h("input", { type: "checkbox", "aria-label": `channel ${r.channel} fixed` });
    box2.checked = r.fixed;
    box2.addEventListener("change", () => chi2Edit(ctx, `fix ${level.label} channel ${r.channel}`,
      () => ctx.api.patch(channelPath(r), { fixed: box2.checked }), msg));
    return box2;
  }
}

/** An edit after which the status line names the χ² before and after (both the engine's, as sent). */
function chi2Edit(ctx, what, run, msg) {
  const before = ctx.evaluation?.chi2;
  return ctx.edit(what, run, msg).then(() => {
    const after = ctx.evaluation?.chi2;
    if (before !== undefined && after !== undefined) {
      ctx.status(`${what}: χ² ${formatNumber(before, 10)} → ${formatNumber(after, 10)} (Undo: Ctrl+Z)`, "ok");
    }
  }, () => {});
}
