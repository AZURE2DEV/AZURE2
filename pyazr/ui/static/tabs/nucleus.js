// The Nucleus tab (docs/upgrade/web_ui.md section 6.2): the particle pairs (added with add_pair, edited with
// set_pair, the THM binding energy among their numbers, removed with remove_pair), the levels with their channels (a
// width, an energy, a fixed flag, whether a width is an amplitude; a new level with the channels allowed_channels
// gives), and the records of the sections NUCLEUS_SECTIONS (γ decays, potentials) through the generated forms.  Every
// edit is a Project call through the project's routes: one entry of the history that Undo takes back, the header
// showing the χ² the engine computes after it, a refusal shown in the engine's words.

import { pathOf } from "../api.js";
import { SECTION_TITLES } from "../layout.js";
import { onlyColumns, pairColumns } from "../lists.js";
import { renderTable } from "../tables.js";
import { button, h, section } from "../ui.js";
import { NUCLEUS_SECTIONS, newRecordPanel, recordPanel, recordRoute } from "./records.js";
import { levelPanel, newLevelPanel, valueOf } from "./levels.js";

// the rules of a level's requests live with its panels; they are part of this tab's interface
export { newLevelRequest, settingOf, settingRequest, settingText } from "./levels.js";

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
  const selectedPair = args[0] === "pair" ? pairs.find((p) => String(p.key) === args[1]) : null;
  const pairEdit = h("section.panel.detail.pair-edit", { "aria-label": "the selected pair" });
  // Add and Remove at the table (the old program's + and -); the selected pair's numbers beside it
  const addPair = button("Add pair", "Add a pair (POST /api/project/pairs: add_pair)",
    () => ctx.navigate("#/nucleus/pair/+new"), pairs.length ? "" : "primary");
  addPair.dataset.action = "tool-add-pair";
  const removePair = button("Remove pair", selectedPair ? `Remove pair ${selectedPair.key} (refused while a level `
    + "has a channel of it; the pair's panel then offers to remove its channels too)" : "Select a pair first",
  () => ctx.edit(`remove pair ${selectedPair.key}`, () => ctx.api.del(pathOf("/api/project/pairs", selectedPair.key),
    {}), msg).then(() => ctx.navigate("#/nucleus/pair"), () => {}), "danger");
  removePair.disabled = !selectedPair;
  removePair.dataset.action = "tool-remove-pair";
  const pairsBox = section("Particle Pairs", h("div.table-tools", {}, addPair, removePair,
    h("span.hint", {}, "select a row to edit it")), pairTable.element,
  pairs.length ? null : h("p.hint.empty", {}, "No particle pairs yet. A pair is the two particles of a channel "
    + "(the entrance, an exit, or a photon and the final nucleus); the levels' channels are made from them."));
  if (selectedPair) showPair(selectedPair);
  else if (args[0] === "pair" && args[1] === "+new") pairEdit.appendChild(newPairForm());
  else {
    pairEdit.appendChild(h("p.detail-empty", {}, pairs.length ? "Choose a pair in the list to see and edit its "
      + "numbers, or Add pair." : "Add the first pair: the entrance pair, then every exit pair."));
  }

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
    // the partners by name: charges, masses and the separation energy from the mass table (GET /api/masses)
    const lightName = h("input", { type: "text", size: 6, "aria-label": "light nuclide", placeholder: "p" });
    const heavyName = h("input", { type: "text", size: 6, "aria-label": "heavy nuclide", placeholder: "12C" });
    const fill = button("Fill from the mass table", "Z, the atomic masses and the separation energy from AME2020 "
      + "(pyazr.sources.masses); J and parity are not in it", async () => {
      try {
        const q = new URLSearchParams({ heavy: heavyName.value.trim() });
        if (lightName.value.trim()) q.set("light", lightName.value.trim());
        const m = await ctx.api.get(`/api/masses?${q}`);
        const set = (label, value) => {
          const input = boxOf(label);
          if (!input || value === undefined) return;
          input.value = String(value);
          input.dispatchEvent(new Event("input", { bubbles: true }));
        };
        set("light Z", m.Z1); set("light mass", m.M1); set("heavy Z", m.Z2); set("heavy mass", m.M2);
        set("separation energy (MeV)", m.S);
        ctx.status(`from ${m.table}: Z and masses${m.S !== undefined ? ", S = " + m.S + " MeV" : ""}; J and `
          + "parity still to choose", "ok");
      } catch (e) {
        ctx.showError(e, form);
      }
    });
    fill.dataset.action = "fill-masses";
    const form = h("div.new-pair", {}, h("h3", {}, "New pair"), h("p.hint", {}, "A pair no level names yet is kept "
      + "in <levels> by a row of its own (inactive); a level that names it takes its place. A photon pair has no "
      + "light partner: the heavy one is the compound nucleus in its final level."),
    h("div.inline-form", {}, h("label.ctl", {}, h("span.ctl-label", {}, "kind"), kind)),
    h("div.inline-form", {}, h("span.k", {}, "by name "), h("label.ctl", {}, h("span.ctl-label", {}, "light"),
      lightName), h("label.ctl", {}, h("span.ctl-label", {}, "heavy"), heavyName), fill),
    partner("light"), partner("heavy"), h("div.inline-form", {}, ...others), h("div.form-actions", {}, add));
    const boxOf = (label) => form.querySelector(`[aria-label="${label}"]`);
    return form;
  }

  // -- levels
  const selectedLevel = args[0] === "level" ? levels.find((l) => l.label === args[1]) : null;
  const newLevel = args[0] === "level" && args[1] === "+new";
  ctx.levelLabels = levels.map((l) => l.label);
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
  const addLevel = button("Add level", "Add a level with the channels its J^π allows (add_level)",
    () => ctx.navigate("#/nucleus/level/+new"), levels.length || !pairs.length ? "" : "primary");
  addLevel.dataset.action = "tool-add-level";
  const removeLevel = button("Remove level", selectedLevel ? `Remove the level ${selectedLevel.label} and its `
    + "channels (remove_level; Undo takes it back)" : "Select a level first", () => {
    if (!window.confirm(`Remove the level ${selectedLevel.jpi} at ${selectedLevel.energy_text} MeV and its `
      + `${selectedLevel.channels.length} channel(s)? (Undo takes it back.)`)) return;
    ctx.edit(`remove level ${selectedLevel.label}`, () => ctx.api.del(pathOf("/api/project/levels",
      selectedLevel.label), {}), msg).then(() => ctx.navigate("#/nucleus/level"), () => {});
  }, "danger");
  removeLevel.disabled = !selectedLevel;
  removeLevel.dataset.action = "remove-level";
  const levelsBox = section("Levels", h("div.table-tools", {}, addLevel, removeLevel,
    h("span.hint", {}, "select a row: its channels open beside the list")), levelTable.element,
  levels.length ? null : h("p.hint.empty", {}, pairs.length ? "No levels yet. A level is a state of the compound "
    + "nucleus with its J^π and energy; its channels are the pairs' that J^π allows, offered when you add it."
    : "No levels yet: add the particle pairs first (Particle Pairs), then the levels, whose channels come from "
      + "them."));
  const levelDetail = h("div.detail", { "aria-label": "the selected level" });
  if (selectedLevel) {
    const settings = await ctx.api.get("/api/project/parameter_settings");
    levelDetail.appendChild(levelPanel(ctx, selectedLevel, pairLabel, msg, settings));
  } else if (newLevel) levelDetail.appendChild(newLevelPanel(ctx, pairLabel, msg));
  else {
    levelDetail.appendChild(h("p.detail-empty", {}, levels.length ? "Choose a level in the list: its energy and "
      + "its channels, with their widths, open here." : "Add level: J, parity and energy; the channels its J^π "
      + "allows are offered."));
  }

  // -- the records of the sections the Nucleus tab claims (γ decays, potentials, channel radii, isospin)
  const kinds = NUCLEUS_SECTIONS.flatMap((s) => ctx.schema.sections[s] || []);
  const words = NUCLEUS_SECTIONS.filter((s) => ctx.schema.sections[s]).map((s) => SECTION_TITLES[s] || s);
  const mine = records.filter((r) => kinds.includes(r.kind));
  const recTable = renderTable(document, {
    rowKey: (r) => r.name, label: "records", selected: args[0], withFilter: false,
    caption: `Records of ${words.join(", ")}`,
    columns: [{ key: "name", title: "name" }, { key: "kind", title: "kind" }, { key: "line", title: "line",
      numeric: true }],
    rows: mine, onSelect: (r) => ctx.navigate(recordRoute(ctx.schema, r)),
  });
  const recordsBox = h("div", {}, section(words.join(", "), h("div.table-tools", {}, button("Add…",
    `Add a record of ${kinds.join(", ")}`, () => ctx.navigate("#/nucleus/+new"))), recTable.element));
  if (args[0] === "+new") {
    recordsBox.appendChild(await newRecordPanel(ctx, kinds, records,
      { onAdded: (r) => ctx.navigate(recordRoute(ctx.schema, r)) }));
  }
  const record = args[0] && args[0] !== "level" && args[0] !== "pair" ? mine.find((r) => r.name === args[0]) : null;
  if (record) {
    recordsBox.appendChild(await recordPanel(ctx, record, records, problems,
      { onGone: () => ctx.navigate("#/nucleus") }));
  }
  // one screen per tab: the pairs with the selected pair beside them; the levels with the selected level's channels
  // beside them, the further records one click away (open in the expert's depth, and first when one is chosen)
  if (args[0] === "pair") { el.append(h("div.master-detail", {}, pairsBox, pairEdit)); return; }
  const more = h("details.more", { data: { show: "records" } }, h("summary", {}, `More: ${words.join(", ")} `
    + `(${mine.length})`), recordsBox);
  const chosen = Boolean(record) || args[0] === "+new";
  if (chosen || ctx.expert) more.open = true;
  const main = h("div.master-detail", {}, levelsBox, levelDetail);
  el.append(...(chosen ? [more, main] : [main, more]));
}
