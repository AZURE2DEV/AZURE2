// The Nucleus tab (docs/upgrade/web_ui.md section 6.2): the particle pairs (added with add_pair, edited with
// set_pair, the THM binding energy among their numbers, removed with remove_pair), the levels with their channels (a
// width, an energy, a fixed flag, whether a width is an amplitude; a new level with the channels allowed_channels
// gives), and the records of the sections NUCLEUS_SECTIONS (γ decays, potentials) through the generated forms.  Every
// edit is a Project call through the project's routes: one entry of the history that Undo takes back, the header
// showing the χ² the engine computes after it, a refusal shown in the engine's words.

import { pathOf } from "../api.js";
import { onlyColumns, pairColumns } from "../lists.js";
import { formatNumber, renderTable } from "../tables.js";
import { button, h, section } from "../ui.js";
import { NUCLEUS_SECTIONS, newRecordPanel, recordPanel, recordRoute } from "./records.js";

// A pair's numbers and a new pair's partners are described by pyazr's table (kinds()["pairs"]: `fields` with title,
// unit and doc, `partner_fields` of add_pair's dicts {J, parity, Z, mass, g}, `kinds`); the fields a form offers are
// the routes' arguments (GET /api/describe), so a field a route gains is shown with its name where the table lacks
// it.
/** A field's words for a form: {title, unit, doc} from pyazr's table, or the name alone. */
export function pairWords(table, name) {
  const f = table?.fields?.[name] || table?.arguments?.[name];
  return f ? { title: f.title, unit: f.unit || "", doc: f.doc || "" } : { title: name.replace(/_/g, " "), unit: "",
    doc: "" };
}
const STRUCTURED = new Set(["kind", "light", "heavy"]);
/** The fields of a new pair's partner (add_pair's dicts), from pyazr's table: [[key, {title, unit, doc, values?}]]. */
export const partnerFields = (table) => Object.entries(table?.partner_fields || {});

/** The arguments of a route of the server's table (GET /api/describe: routes[].arguments), without the path's own
 * and the revision. */
export function routeArguments(routes, method, path) {
  const r = (routes || []).find((x) => x.method === method && x.path === path);
  const own = new Set((path.match(/\{(\w+)/g) || []).map((m) => m.slice(1)));
  return r ? r.arguments.map((a) => a.split(":")[0]).filter((a) => !own.has(a) && a !== "rev") : [];
}

/** A number typed as a number; any other text as typed (pyazr converts it or refuses it with its words). */
const valueOf = (t) => (t !== "" && Number.isFinite(Number(t)) ? Number(t) : t);

/** The body of POST /api/project/pairs from a new pair's texts ({kind, "light.J", ..., field: text}): the partners
 * as add_pair's dicts (none for a photon's light partner), the route's other arguments where typed. */
export function newPairRequest(texts, fields, table) {
  const body = { kind: texts.kind || (table?.kinds || ["particle"])[0] };
  for (const side of ["light", "heavy"]) {
    if (side === "light" && body.kind === "photon") continue;
    const partner = {};
    for (const [k] of partnerFields(table)) {
      const t = String(texts[`${side}.${k}`] ?? "").trim();
      if (t !== "") partner[k] = valueOf(t);
    }
    body[side] = partner;
  }
  for (const f of fields) {
    const t = String(texts[f] ?? "").trim();
    if (!STRUCTURED.has(f) && t !== "") body[f] = valueOf(t);
  }
  return body;
}

/** The body of POST /api/project/levels: J, parity, E (MeV), whether its energy is fixed, and the channels checked
 * among those allowed_channels gave ({pair, l, s}, width 0, fixed) -- or "allowed" when every one is checked
 * (add_level then gives a J^π the project has its group's channels). */
export function newLevelRequest({ J, parity, energy, fixed }, allowed, checked) {
  const chosen = allowed.filter((c, i) => checked[i]);
  return { J: valueOf(String(J).trim()), parity: Number(parity), energy: valueOf(String(energy).trim()),
    fixed: Boolean(fixed), channels: chosen.length === allowed.length && allowed.length ? "allowed"
      : chosen.map((c) => ({ pair: c.pair, l: c.l, s: c.s, width: 0, fixed: true })) };
}

export async function render(ctx, el, args) {
  const [pairs, levels, records, problems] = await Promise.all([ctx.api.get("/api/project/pairs"),
    ctx.api.get("/api/project/levels"), ctx.api.get("/api/project/records"), ctx.api.get("/api/project/validate")]);
  ctx.pairs = pairs;
  const pairLabel = Object.fromEntries(pairs.map((p) => [p.key, p.label]));
  const msg = h("div.messages");
  el.appendChild(msg);

  // -- pairs
  const routes = ctx.describe?.routes || [];
  if (!ctx.kinds) ctx.kinds = await ctx.api.get("/api/kinds");
  const pairTable = renderTable(document, {
    rowKey: (r) => String(r.key), label: "pairs", withFilter: false, selected: args[0] === "pair" ? args[1] : null,
    caption: "Particle pairs of <levels>", columns: ctx.expert ? pairColumns(ctx.kinds.pairs?.fields, pairs)
      : onlyColumns(pairColumns(ctx.kinds.pairs?.fields, pairs), ctx.describe?.first_view?.pairs),
    rows: pairs, onSelect: (r) => ctx.navigate(`#/nucleus/pair/${r.key}`),
  });
  const pairEdit = h("div.pair-edit");
  const pairsBox = section("Pairs", pairTable.element,
    pairs.length ? null : h("p.hint.empty", {}, "No particle pairs yet. A pair is the two particles of a channel "
      + "(the entrance, an exit, or a photon and the final nucleus); the levels' channels are made from them."),
    button("+ pair", "Add a pair (POST /api/project/pairs: add_pair)", () => ctx.navigate("#/nucleus/pair/+new"),
      pairs.length ? "" : "primary"),
    pairEdit);
  const selectedPair = args[0] === "pair" ? pairs.find((p) => String(p.key) === args[1]) : null;
  if (selectedPair) showPair(selectedPair);
  else if (args[0] === "pair" && args[1] === "+new") pairEdit.appendChild(newPairForm());

  function showPair(p) {
    const fields = routeArguments(routes, "PATCH", "/api/project/pairs/{key}");
    const inputs = fields.map((k) => {
      const w = pairWords(ctx.kinds.pairs, k);
      const input = h("input.num", { type: "text", value: p[k] === undefined ? "" : String(p[k]), size: 12,
        "aria-label": w.unit ? `${k} (${w.unit})` : k });
      return [k, input, h("label.ctl", { title: `${w.doc} (${k})` }, h("span.ctl-label", {},
        `${w.title}${w.unit ? ` (${w.unit})` : ""}`), input)];
    });
    const apply = button(`Apply to pair ${p.key}`, "set_pair with the changed values", () => {
      const changed = Object.fromEntries(inputs.filter(([k, i]) => i.value.trim() !== String(p[k] ?? ""))
        .map(([k, i]) => [k, valueOf(i.value.trim())]));
      if (!Object.keys(changed).length) return;
      ctx.edit(`set pair ${p.key}`, () => ctx.api.patch(pathOf("/api/project/pairs", p.key), changed), msg)
        .catch(() => {});
    }, "primary");
    const removal = h("div.form-actions");
    const removeWith = button(`Remove pair ${p.key} and its channels`, "remove_pair(remove_channels=True): the "
      + "channel rows of the pair go with it, and with them their fit parameters", () => ctx.edit(
      `remove pair ${p.key} and its channels`, () => ctx.api.del(pathOf("/api/project/pairs", p.key),
        { remove_channels: true }), msg).then(() => ctx.navigate("#/nucleus"), () => {}));
    removeWith.dataset.action = "remove-pair-channels";
    const remove = button(`Remove pair ${p.key}`, "remove_pair: refused while a level has a channel of it (the "
      + "refusal names them; removing the channels too is the next button)", () => ctx.edit(`remove pair ${p.key}`,
      () => ctx.api.del(pathOf("/api/project/pairs", p.key), {}), msg).then(() => ctx.navigate("#/nucleus"), () => {
      if (!removal.contains(removeWith)) removal.appendChild(removeWith);
    }));
    remove.dataset.action = "remove-pair";
    removal.appendChild(remove);
    pairEdit.append(h("h3", {}, `Pair ${p.key}: ${p.label}`),
      h("div.inline-form", {}, ...inputs.map(([, , label]) => label), apply), removal);
  }

  function newPairForm() {
    const fields = routeArguments(routes, "POST", "/api/project/pairs");
    const texts = {};
    const box = (name, label, title, value = "", words = label) => {
      texts[name] = value;
      const input = h("input.num", { type: "text", value, size: 10, "aria-label": label, title });
      input.addEventListener("input", () => { texts[name] = input.value; });
      return h("label.ctl", { title }, h("span.ctl-label", {}, words), input);
    };
    const table = ctx.kinds.pairs;
    const kind = h("select", { "aria-label": "kind" }, (table?.kinds || []).map((k) => h("option", { value: k }, k)));
    texts.kind = kind.value;
    kind.addEventListener("change", () => { texts.kind = kind.value; });
    // a choice where the table lists the values (the parity: + or -), a box otherwise; labelled in the table's words
    const choice = (name, label, f) => {
      const select = h("select", { "aria-label": label, title: f.doc }, f.values.map((v, i) => h("option",
        { value: String(v) }, f.labels?.[i] ?? String(v))));
      texts[name] = select.value;
      select.addEventListener("change", () => { texts[name] = select.value; });
      return h("label.ctl", { title: f.doc }, h("span.ctl-label", {}, f.title), select);
    };
    const partner = (side) => h("div.inline-form", {}, h("span.k", {}, `${side} partner `),
      ...partnerFields(table).map(([k, f]) => (f.values ? choice(`${side}.${k}`, `${side} ${k}`, f)
        : box(`${side}.${k}`, `${side} ${k}`, `${f.doc}${f.unit ? ` (${f.unit})` : ""}`, "",
          `${f.title}${f.unit ? ` (${f.unit})` : ""}`))));
    const others = fields.filter((f) => !STRUCTURED.has(f)).map((f) => {
      const w = pairWords(table, f);
      return box(f, w.unit ? `${f.replace(/_/g, " ")} (${w.unit})` : f.replace(/_/g, " "), `${w.doc} (${f})`, "",
        `${w.title}${w.unit ? ` (${w.unit})` : ""}`);
    });
    const add = button("Add pair", "add_pair (POST /api/project/pairs): the engine checks every number", () => ctx.edit(
      "add a pair", () => ctx.api.post("/api/project/pairs", newPairRequest(texts, fields, table)), msg)
      .then((answer) => ctx.navigate(`#/nucleus/pair/${answer?.result ?? ""}`), () => {}), "primary");
    add.dataset.action = "add-pair";
    return h("div.new-pair", {}, h("h3", {}, "New pair"), h("p.hint", {}, "A pair no level names yet is kept in "
      + "<levels> by a row of its own (inactive); a level that names it takes its place. A photon pair has no light "
      + "partner: the heavy one is the compound nucleus in its final level."),
    h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "kind"), kind)), partner("light"),
    partner("heavy"), h("div.inline-form", {}, ...others), h("div.form-actions", {}, add));
  }

  // -- levels
  const selectedLevel = args[0] === "level" ? levels.find((l) => l.label === args[1]) : null;
  const newLevel = args[0] === "level" && args[1] === "+new";
  const levelRows = levels.map((l) => ({ ...l, n: l.channels.length,
    on: l.active ? "on" : "off", efix: l.fixed ? "fixed" : "free" }));
  const levelTable = renderTable(document, {
    rowKey: (r) => r.label, label: "levels", selected: selectedLevel?.label, caption: "Levels of <levels>",
    columns: [{ key: "jpi", title: "Jπ" }, { key: "energy", title: "E_λ", unit: "MeV", numeric: true,
      value: (r) => r.energy, format: (v, r) => r.energy_text },
    { key: "efix", title: "energy" }, { key: "on", title: "active" }, { key: "n", title: "channels", numeric: true },
    { key: "nucleus", title: "nucleus", format: (v) => (Array.isArray(v) ? `Z ${v[0]}, A ${v[1]}` : "") }],
    rows: levelRows,
    onSelect: (r) => ctx.navigate(`#/nucleus/level/${encodeURIComponent(r.label)}`),
  });
  const levelsBox = h("div", {}, section("Levels", levelTable.element,
    levels.length ? null : h("p.hint.empty", {}, pairs.length ? "No levels yet. A level is a state of the compound "
      + "nucleus with its J^π and energy; its channels are the pairs' that J^π allows, offered when you add it."
      : "No levels yet: add the particle pairs first (Particle Pairs), then the levels, whose channels come from "
        + "them."),
    button("+ level", "Add a level with the channels its J^π allows (add_level)",
      () => ctx.navigate("#/nucleus/level/+new"), levels.length || !pairs.length ? "" : "primary")));
  if (selectedLevel) levelsBox.appendChild(levelPanel(ctx, selectedLevel, pairLabel, msg));
  if (newLevel) levelsBox.appendChild(newLevelPanel(ctx, pairLabel, msg));

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
  const recordsBox = h("div", {}, section("γ decays and potentials", recTable.element,
    button("+ new", `Add a record of ${kinds.join(", ")}`, () => ctx.navigate("#/nucleus/+new"))));
  if (args[0] === "+new") {
    recordsBox.appendChild(await newRecordPanel(ctx, kinds, records,
      { onAdded: (r) => ctx.navigate(recordRoute(ctx.schema, r)) }));
  }
  const record = args[0] && args[0] !== "level" && args[0] !== "pair" ? mine.find((r) => r.name === args[0]) : null;
  if (record) {
    recordsBox.appendChild(await recordPanel(ctx, record, records, problems,
      { onGone: () => ctx.navigate("#/nucleus") }));
  }
  // what was chosen comes first: a record (or a new one) above the levels, a pair above both; the pairs, rarely
  // edited, last otherwise
  if (!ctx.expert) {
    // the first view: one screen per task (layout.js FIRST_VIEW); the γ decays and potentials one click away
    if (args[0] === "pair") { el.append(pairsBox); return; }
    const more = h("details.more", {}, h("summary", {}, `More: γ decays and potentials (${mine.length})`),
      recordsBox);
    if (record || args[0] === "+new") more.open = true;
    el.append(...(record || args[0] === "+new" ? [more, levelsBox] : [levelsBox, more]));
    return;
  }
  const recordFirst = Boolean(record) || args[0] === "+new";
  const boxes = recordFirst ? [recordsBox, levelsBox] : [levelsBox, recordsBox];
  el.append(...(args[0] === "pair" ? [pairsBox, ...boxes] : [...boxes, pairsBox]));
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
    { key: "amplitude", title: "amplitude", sortable: false, render: (r) => amplitudeEditor(r),
      doc: "whether the width column is a reduced-width amplitude (MeV^1/2) instead of a partial width or ANC (a "
        + "particle channel only; the value is not converted: type the width in the new convention)" },
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
  function amplitudeEditor(r) {
    const flag = h("input", { type: "checkbox", "aria-label": `channel ${r.channel} width is an amplitude` });
    flag.checked = Boolean(r.amplitude);
    flag.addEventListener("change", () => chi2Edit(ctx, `amplitude flag of ${level.label} channel ${r.channel}`,
      () => ctx.api.patch(channelPath(r), { width_is_amplitude: flag.checked }), msg));
    return flag;
  }
  function fixedEditor(r) {
    const box2 = h("input", { type: "checkbox", "aria-label": `channel ${r.channel} fixed` });
    box2.checked = r.fixed;
    box2.addEventListener("change", () => chi2Edit(ctx, `fix ${level.label} channel ${r.channel}`,
      () => ctx.api.patch(channelPath(r), { fixed: box2.checked }), msg));
    return box2;
  }
}

/** A new level: J, parity, E and whether its energy is fixed; the channels its J^π allows with the project's pairs
 * (GET /api/project/allowed_channels: the engine's rule), each to be checked; Add sends add_level. */
function newLevelPanel(ctx, pairLabel, msg) {
  const box = h("section.panel.new-level", { "aria-label": "new level" }, h("h2", {}, "New level"));
  const J = h("input.num", { type: "text", value: "0.5", size: 6, "aria-label": "J",
    title: "the spin as a number: 0.5, 1, 1.5, ..." });
  const parity = h("select", { "aria-label": "parity" }, h("option", { value: "1" }, "+"),
    h("option", { value: "-1" }, "−"));
  const energy = h("input.num", { type: "text", value: "", size: 12, "aria-label": "E_λ (MeV)" });
  const lmax = h("input.num", { type: "text", value: "2", size: 4, "aria-label": "l max" });
  const fixed = h("input", { type: "checkbox", "aria-label": "energy fixed" });
  fixed.checked = true;
  const list = h("div.allowed-channels");
  let allowed = [];
  let checks = [];
  const jValue = () => valueOf(J.value.trim());
  const channels = button("Channels", "The channels this J^π has with the project's pairs (allowed_channels)",
    async () => {
      list.textContent = "";
      try {
        allowed = await ctx.api.get(`/api/project/allowed_channels?J=${encodeURIComponent(jValue())}&parity=`
          + `${parity.value}&l_max=${encodeURIComponent(lmax.value.trim())}`, { quiet: true });
      } catch (e) {
        allowed = [];
        ctx.showError(e, list);
        return;
      }
      checks = allowed.map((c) => {
        const input = h("input", { type: "checkbox", "aria-label": `channel pair ${c.pair} l ${c.l} s ${c.s}` });
        input.checked = true;
        return input;
      });
      list.append(h("p.hint", {}, `${allowed.length} channel(s), each added with width 0, fixed; with every one `
        + "checked, a J^π the project has gets its group's channels."),
      ...allowed.map((c, i) => h("label.choice", {}, checks[i], ` pair ${c.pair} (${pairLabel[c.pair] ?? ""}), `
        + `l ${c.l}, s ${c.s}`)));
    });
  channels.dataset.action = "allowed-channels";
  const add = button("Add level", "add_level (POST /api/project/levels)", () => {
    const body = newLevelRequest({ J: jValue(), parity: parity.value, energy: energy.value, fixed: fixed.checked },
      allowed, checks.map((c) => c.checked));
    ctx.edit("add a level", () => ctx.api.post("/api/project/levels", body), msg)
      .then(() => ctx.navigate("#/nucleus"), () => {});
  }, "primary");
  add.dataset.action = "add-level";
  box.append(h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "J"), J),
    h("label.ctl", {}, h("span.ctl-label", {}, "parity"), parity),
    h("label.ctl", {}, h("span.ctl-label", {}, "E_λ (MeV)"), energy),
    h("label.ctl", {}, fixed, h("span.ctl-label", {}, "energy fixed")),
    h("label.ctl", {}, h("span.ctl-label", {}, "l max"), lmax), channels), list, h("div.form-actions", {}, add));
  return box;
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
