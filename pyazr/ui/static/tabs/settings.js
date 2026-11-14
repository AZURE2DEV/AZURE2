// The run's settings (the Project tab): every entry of <config> the engine knows, as a row with its control, built
// from the engine's table (GET /api/kinds: kinds.config, read from io/src/ConfigSection.cpp): a switch for a key that
// is 0 or 1, a choice for a positional entry with values, a number box with its bound, a text box for a path or a
// directory.  Each row says the setting in the old program's words, its key, what it does, what the engine does
// when the project does not set it, and whether the project sets it.  A change is sent at once (PATCH
// /api/project/config: one history entry, Undo takes it back); the engine refuses a value in its words.  The
// rarely set rows (the table's `rare`) wait under "More settings" unless the project sets them.

import { button, h, section } from "../ui.js";

const SWITCH = ["0", "1"];

/** Whether an entry is a switch (a key of value 0 or 1). */
export const isSwitch = (field) => field.type === "enum" && field.values?.length === 2
  && field.values.every((v, i) => v === SWITCH[i]);

/**
 * The rows of the settings, in the table's order: {name, field, value (the project's text or null), set (the
 * project sets it: a key present, a positional entry other than its first value), inUse (the value the engine uses:
 * the project's, else the option's default, else null), group, rare (folded: rare and not set), shown (matches the
 * query)}.  `config`: GET /api/project/config ({name: value}); `query`: words of the search box.
 */
export function settingRows(fields, config, query = "") {
  const words = String(query).toLowerCase().split(/\s+/).filter(Boolean);
  return Object.entries(fields || {}).map(([name, field]) => {
    const value = Object.hasOwn(config || {}, name) ? String(config[name]) : null;
    const set = field.positional ? value !== null && field.values?.length > 0 && value !== field.values[0]
      : value !== null;
    const inUse = value ?? (field.default !== undefined ? String(field.default) : null);
    const text = `${name} ${field.title} ${field.doc} ${field.group}`.toLowerCase();
    return { name, field, value, set, inUse, group: field.group || "", rare: Boolean(field.rare) && !set,
      shown: words.every((w) => text.includes(w)) };
  });
}

/** The groups of the rows in the order they first appear (the table's order: the old Runtime Options first). */
export function groupsOf(rows) {
  return [...new Set(rows.map((r) => r.group))];
}

/** The value a control sends for a switch ticked or not. */
export const switchValue = (checked) => (checked ? "1" : "0");

/** The section: a search box, the groups of rows, the rare rows under "More settings", the free entries behind
 * "Advanced". `send(entries)` sends {name: value | null} through the config route. */
export function settingsSection(ctx, fields, config, send, advanced) {
  const box = section("Run settings");
  const search = h("input.table-filter", { type: "search", placeholder: "search the settings (/)",
    "aria-label": "search the settings" });
  const body = h("div.settings");
  box.append(h("p.hint", {}, "The settings of a run, as the engine's table lists them; a change is applied at once "
    + "(Undo takes it back). The old program's Runtime Options come first."), search, body);
  const draw = () => {
    body.textContent = "";
    const rows = settingRows(fields, config, search.value);
    const often = rows.filter((r) => r.shown && (ctx.expert || !r.rare));
    const rare = rows.filter((r) => r.shown && !ctx.expert && r.rare);
    body.appendChild(groupTables(often, send));
    if (rare.length) {
      const more = h("details.more", {}, h("summary", {}, `More settings (${rare.length}): `
        + `${groupsOf(rare).join(", ")}`), groupTables(rare, send));
      if (search.value.trim()) more.open = true;
      body.appendChild(more);
    }
    if (!often.length && !rare.length) body.appendChild(h("p.hint", {}, "No setting matches the search."));
  };
  search.addEventListener("input", draw);
  draw();
  if (advanced) box.appendChild(h("details.more", {}, h("summary", {}, "Advanced: the entries as text (any key)"),
    advanced));
  return box;
}

function groupTables(rows, send) {
  const out = h("div");
  for (const g of groupsOf(rows)) {
    const tbody = h("tbody");
    for (const r of rows.filter((x) => x.group === g)) tbody.appendChild(settingRow(r, send));
    out.appendChild(h("table.settings-group", { "aria-label": `${g} settings` },
      h("caption", {}, g), tbody));
  }
  return out;
}

/** One row: the title and key, the control, what the engine uses when unset, the meaning, set or not. */
function settingRow(r, send) {
  const f = r.field;
  const label = `${f.title} (${r.name})`;
  let control;
  if (f.readonly) {
    control = h("span.mono", {}, r.value ?? "");
  } else if (isSwitch(f)) {
    control = h("input", { type: "checkbox", "aria-label": label });
    control.checked = r.inUse === "1" || r.inUse === "true";
    control.addEventListener("change", () => send({ [r.name]: switchValue(control.checked) }));
  } else if (f.type === "enum") {
    control = h("select", { "aria-label": label }, (f.values || []).map((v) => h("option", { value: v }, v)));
    if (r.inUse !== null && !(f.values || []).includes(r.inUse)) control.appendChild(h("option", { value: r.inUse },
      r.inUse));
    control.value = r.inUse ?? "";
    control.addEventListener("change", () => send({ [r.name]: control.value }));
  } else {
    const numeric = f.type === "integer" || f.type === "number";
    control = h(numeric ? "input.num" : "input.mono", { type: "text", value: r.value ?? "", size: numeric ? 8 : 24,
      "aria-label": label, placeholder: numeric && f.min !== undefined ? `≥ ${f.min}` : "" });
    const apply = () => {
      const t = control.value.trim();
      if (t === (r.value ?? "")) return;
      send({ [r.name]: t === "" && !f.positional ? null : t });
    };
    control.addEventListener("keydown", (e) => { if (e.key === "Enter") { e.preventDefault(); apply(); } });
    control.addEventListener("change", apply);
  }
  const state = r.field.positional ? "" : (r.value !== null ? "set" : "not set");
  const reset = !f.positional && r.value !== null && !f.readonly
    ? button("×", `Remove ${r.name} from the project: the engine then does what it does without it`,
      () => send({ [r.name]: null }), "quiet") : null;
  const unset = f.absent ? `without it: ${f.absent}` : (f.default !== undefined ? `default ${f.default}` : "");
  return h("tr", { data: { setting: r.name } },
    h("th", {}, h("div", {}, f.title), h("div.mono.hint", {}, r.name)),
    h("td", {}, control, reset),
    h("td.hint", {}, f.doc, unset ? h("div", {}, unset) : null),
    h("td.hint", {}, state));
}
