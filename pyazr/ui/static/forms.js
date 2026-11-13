// Forms generated from the engine's schema (pyazr.schema(); docs/upgrade/web_ui.md sections 1.3 and 10).
//
// WIDGETS is the one table of widgets: one row per value type of the schema.  A row turns a field (or a part of a
// structured field) and its text into controls (inputs and drop-downs, described as plain data) and composes the
// controls' values back into the field's text.  The text is the project's own grammar; the page splits it into
// tokens to show them in separate boxes and joins them again, nothing more: whether a value is valid is the engine's
// decision when the edit is sent (section 5), and a text a widget cannot split is shown whole in one text box.
// A new value type of the engine needs one row here and its test (tests/pyazr/ui/js/forms.test.mjs fails until
// every type of the schema has a row).
//
// describeForm(schema, kind, fields, ctx) gives the whole form as data: one entry per field, in the schema's order
// unless the layout hints (layout.js) put some on one line; formview.js draws it.  No DOM here (node tests it).

import { LAYOUT } from "./layout.js";

/** Split a field's text into tokens; a quoted token keeps its quotes (formulas, paths with blanks). */
export function tokenize(text) {
  const out = [];
  const re = /"[^"]*"|\S+/g;
  let m;
  while ((m = re.exec(String(text ?? ""))) !== null) out.push(m[0]);
  return out;
}

const join = (parts) => parts.filter((p) => p !== undefined && p !== null && String(p).trim() !== "").join(" ");

/** The units a field or part may be written in: its own list, or the quantity's for a dimension "@quantity". */
export function unitsOf(field, ctx) {
  if (Array.isArray(field.units) && field.units.length) return field.units;
  if (field.dimension === "@quantity" && ctx.quantity && ctx.schema?.quantities?.[ctx.quantity]) {
    return Object.keys(ctx.schema.quantities[ctx.quantity].units || {});
  }
  if (field.dimension && ctx.schema?.units) {
    return Object.entries(ctx.schema.units).filter(([, u]) => u.dimension === field.dimension).map(([k]) => k);
  }
  return [];
}

const input = (name, value, extra = {}) => ({ kind: "input", name, value: value ?? "", ...extra });
const select = (name, value, values, extra = {}) => {
  const options = values.map((v) => (typeof v === "object" ? v : { value: v, label: v === "" ? "—" : v }));
  if (value && !options.some((o) => o.value === value)) options.push({ value, label: value });
  return { kind: "select", name, value: value ?? "", options, ...extra };
};
/** A text box holding the whole text: the fallback of every widget, and the widget of free-form types. */
const whole = (field, text, extra = {}) => ({
  controls: [input("text", text, { placeholder: field.grammar || "", wide: true, ...extra })],
  compose: (v) => (v.text ?? "").trim(),
});

/** Names of the project's records of the given kinds (for references). */
function recordNames(ctx, kinds) {
  return (ctx.records || []).filter((r) => (kinds || []).includes(r.kind)).map((r) => r.name);
}

// -- the widget table ------------------------------------------------------------------------------------------
// Each row: (field, text, ctx) -> {controls, compose(values) -> text}.  `values` maps a control's name to its value.

export const WIDGETS = {
  number: (f, t) => ({
    controls: [input("v", t, { numeric: true, placeholder: "number" })],
    compose: (v) => v.v.trim(),
  }),
  integer: (f, t) => ({
    controls: [input("v", t, { numeric: true, placeholder: f.min !== undefined ? `integer ≥ ${f.min}` : "integer" })],
    compose: (v) => v.v.trim(),
  }),
  number_list: (f, t) => whole(f, t, { numeric: true }),
  unit: (f, t, ctx) => {
    const units = unitsOf(f, ctx);
    if (!units.length) return whole(f, t);
    return { controls: [select("u", t, ["", ...units])], compose: (v) => v.u };
  },
  number_with_unit: (f, t, ctx) => {
    const tok = tokenize(t);
    if (tok.length > 2) return whole(f, t);
    const units = unitsOf(f, ctx);
    return {
      controls: [input("v", tok[0], { numeric: true, placeholder: "number" }),
        select("u", tok[1] ?? "", ["", ...units])],
      compose: (v) => join([v.v, v.u]),
    };
  },
  range: (f, t, ctx) => {
    const tok = tokenize(t);
    if (tok.length !== 0 && tok.length !== 3) return whole(f, t);
    return {
      controls: [input("lo", tok[0], { numeric: true, placeholder: "low" }),
        input("hi", tok[1], { numeric: true, placeholder: "high", label: "–" }),
        select("u", tok[2] ?? "", ["", ...unitsOf(f, ctx)])],
      compose: (v) => join([v.lo, v.hi, v.u]),
    };
  },
  nuisance: nuisanceWidget,
  enum: (f, t) => ({
    controls: [select("v", t, [{ value: "", label: f.default !== undefined ? `— (${f.default})` : "—" },
      ...(f.values || [])])],
    compose: (v) => v.v,
  }),
  reference: (f, t, ctx) => {
    const [name, epoch] = String(t ?? "").split("@");
    const controls = [select("v", name, ["", ...recordNames(ctx, f.refers_to)], { reference: f.refers_to })];
    if (f.epoch) controls.push(input("epoch", epoch ?? "", { placeholder: "epoch", label: "@" }));
    return { controls, compose: (v) => (v.v ? (v.epoch ? `${v.v}@${v.epoch}` : v.v) : "") };
  },
  name: (f, t) => whole(f, t),
  path: (f, t) => whole(f, t, { path: true }),
  text: (f, t) => whole(f, t),
  formula: (f, t) => whole(f, t, { mono: true }),
  spin: (f, t, ctx) => whole(f, t, { list: ctx.spins || [] }),
  spin_parity: (f, t, ctx) => whole(f, t, { list: ctx.spinParities || [] }),
  nucleus: (f, t) => whole(f, t),
  pair: (f, t, ctx) => whole(f, t, { list: ctx.pairLabels || [] }),
  composition: (f, t) => whole(f, t),
  reaction_expr: (f, t, ctx) => whole(f, t, { list: ctx.reactions || [] }),
  columns: (f, t, ctx) => whole(f, t, { mono: true,
    hint: `roles: ${Object.keys(ctx.schema?.column_roles || {}).join(" ")}` }),
  grid: (f, t, ctx) => whole(f, t, { hint: `kinds: ${(ctx.schema?.lists?.grid_kinds || []).join(", ")}` }),
  structured: structuredWidget,
};

/** The modes a nuisance may take, from the engine's list (profiled only where the field allows it). */
function nuisanceModes(f, ctx) {
  const modes = ctx.schema?.lists?.nuisance_modes || [];
  return modes.filter((m) => m !== "profiled" || f.profiled);
}

function nuisanceWidget(f, t, ctx) {
  const tok = tokenize(t);
  const shared = recordNames(ctx, f.or_refers_to);
  const modes = nuisanceModes(f, ctx);
  const own = { controls: [], compose: null };
  // a shared record's name, or "profiled" alone, or: value [unit] mode [prior c w] [grid_max g]
  let value = "", unit = "", mode = "", centre = "", width = "", gmax = "", ref = "";
  if (tok.length === 1 && (shared.includes(tok[0]) || (f.or_refers_to && !modes.includes(tok[0]) && isNaN(+tok[0])))) {
    ref = tok[0];
  } else if (tok.length === 1 && tok[0] === "profiled") {
    mode = "profiled";
  } else if (tok.length) {
    let i = 0;
    value = tok[i++];
    if (i < tok.length && !modes.includes(tok[i])) unit = tok[i++];
    if (i < tok.length && modes.includes(tok[i])) mode = tok[i++];
    while (i < tok.length) {
      if (tok[i] === "prior" && i + 2 < tok.length) { centre = tok[i + 1] ?? ""; width = tok[i + 2] ?? ""; i += 3; }
      else if (tok[i] === "grid_max") { gmax = tok[i + 1] ?? ""; i += 2; }
      else return whole(f, t);
    }
  }
  const units = unitsOf(f, ctx);
  if (f.or_refers_to) {
    own.controls.push(select("ref", ref, [{ value: "", label: "own value" },
      ...shared.map((s) => ({ value: s, label: `shared: ${s}` }))], { rerender: true, reference: f.or_refers_to }));
  }
  if (!ref) {
    own.controls.push(input("value", value, { numeric: true, placeholder: "value", disabled: mode === "profiled" }));
    if (units.length || unit) own.controls.push(select("unit", unit, ["", ...units]));
    own.controls.push(select("mode", mode, ["", ...modes], { label: "mode", rerender: true }));
    own.controls.push(input("centre", centre, { numeric: true, placeholder: "centre", label: "prior" }));
    own.controls.push(input("width", width, { numeric: true, placeholder: "width or %", label: "±" }));
    own.controls.push(input("gmax", gmax, { numeric: true, placeholder: "grid max", label: "grid max" }));
  }
  own.compose = (v) => {
    if (v.ref) return v.ref;
    if (v.mode === "profiled") return "profiled";
    const prior = (v.centre || v.width) ? ["prior", v.centre, v.width] : [];
    return join([v.value, v.unit, v.mode, ...prior, ...(v.gmax ? ["grid_max", v.gmax] : [])]);
  };
  return own;
}

/** Tokens a part of a structured field takes when it is not the last positional part. */
const ARITY = { number_with_unit: 2, range: 3 };

function structuredWidget(f, t, ctx) {
  const alts = f.alternatives || [];
  const tok = tokenize(t);
  // the alternative named by the first token, else the one without a keyword (an empty text: the first)
  let alt = alts.find((a) => a.keyword && a.keyword === tok[0]);
  let rest = tok;
  if (alt) rest = tok.slice(1);
  else alt = alts.find((a) => !a.keyword) || (tok.length ? null : alts[0]);
  if (!alt) return whole(f, t);
  // the tokens of each part: keyword parts after their name, positional parts in order before the first keyword
  const keywords = new Set(alt.parts.filter((p) => p.keyword).map((p) => p.name));
  const positional = alt.parts.filter((p) => !p.keyword);
  const texts = {};
  let i = 0;
  for (let k = 0; k < positional.length; k += 1) {
    const p = positional[k];
    const last = k === positional.length - 1;
    const take = [];
    const room = () => last || take.length < (ARITY[p.type] || 1);
    while (i < rest.length && !keywords.has(rest[i]) && room()) take.push(rest[i++]);
    texts[p.name] = take.join(" ");
  }
  while (i < rest.length) {
    if (!keywords.has(rest[i])) return whole(f, t);
    const name = rest[i++];
    const take = [];
    while (i < rest.length && !keywords.has(rest[i])) take.push(rest[i++]);
    texts[name] = take.join(" ");
  }
  const controls = [];
  if (alts.length > 1) {
    controls.push(select("alt", alt.keyword, alts.map((a) => ({ value: a.keyword,
      label: a.keyword || a.parts.map((p) => p.name).join(" ") || "(none)", title: a.doc || "" })),
    { rerender: true, label: "" }));
  }
  const parts = alt.parts.map((p) => {
    const row = WIDGETS[p.type] || ((pf, pt) => whole(pf, pt));
    const w = row(p, texts[p.name] ?? "", ctx);
    w.controls.forEach((c, n) => {
      c.name = `${p.name}.${c.name}`;
      if (n === 0) c.label = p.keyword || positional.length > 1 || alt.parts.length > 1 ? p.name : (c.label || "");
      c.title = p.doc || p.grammar || "";
    });
    controls.push(...w.controls);
    return { p, w };
  });
  return {
    controls,
    compose: (v) => {
      const keyword = alts.length > 1 ? (v.alt ?? alt.keyword) : alt.keyword;
      const chosen = alts.find((a) => a.keyword === keyword) || alt;
      if (chosen !== alt) return keyword;  // another alternative: its parts start empty
      const words = [keyword];
      for (const { p, w } of parts) {
        const sub = {};
        for (const [k, val] of Object.entries(v)) if (k.startsWith(`${p.name}.`)) sub[k.slice(p.name.length + 1)] = val;
        const text = w.compose(sub);
        if (text) words.push(p.keyword ? `${p.name} ${text}` : text);
      }
      return join(words);
    },
  };
}


/** The widget of a field: its row of WIDGETS, or a text box for a type the table lacks (never fails). */
export function widgetOf(field, text, ctx) {
  const row = WIDGETS[field.type];
  return row ? row(field, text ?? "", ctx) : whole(field, text ?? "");
}

/** The initial values of a widget's controls, by name (what compose() takes). */
export function valuesOf(widget) {
  return Object.fromEntries(widget.controls.map((c) => [c.name, c.value ?? ""]));
}

// -- whether a field applies --------------------------------------------------------------------------------------

/**
 * Whether the field is shown, and the condition's text when the page cannot tell (the engine decides then):
 * {show, note}.  A field that holds a value is always shown.  `by_quantity` fields show when the record's quantity
 * row lists them; `applies_when` of the kinds field_present and field_is are evaluated on the record's fields.
 */
export function applies(key, field, fields, ctx) {
  if (fields[key] !== undefined && fields[key] !== "") return { show: true, note: "" };
  if (field.by_quantity) {
    const q = ctx.quantity && ctx.schema?.quantities?.[ctx.quantity];
    if (q && !(q.fields && key in q.fields)) return { show: false, note: "" };
    const when = q?.fields?.[key]?.when;
    if (when) return condition(when, fields);
  }
  if (field.applies_when) return condition(field.applies_when, fields);
  return { show: true, note: "" };
}

function condition(when, fields) {
  const present = (k) => fields[k] !== undefined && fields[k] !== "";
  if (when.if === "field_present") {
    const yes = present(when.field);
    return { show: when.not ? !yes : yes, note: "" };
  }
  if (when.if === "field_is") {
    const first = tokenize(Array.isArray(fields[when.field]) ? fields[when.field][0] : fields[when.field])[0];
    return { show: (when.in || []).includes(first), note: "" };
  }
  return { show: true, note: when.text || when.if };
}

// -- the whole form ---------------------------------------------------------------------------------------------

/** A field's texts as a list: a repeatable field holds one text per occurrence. */
export function textsOf(value) {
  if (value === undefined || value === null) return [];
  return Array.isArray(value) ? value.map(String) : [String(value)];
}

/**
 * The form of a record of `kind` as data: {kind, title, doc, lines: [[entry, ...], ...], hidden: [keys]}.
 * An entry: {key, field, title, type, required, repeatable, doc, grammar, note, texts, widgets, diagnostics}.
 * `fields` are the record's fields as the server sends them; `ctx`: {schema, records, quantity, pairLabels,
 * reactions, spins, spinParities, diagnostics: [{record, field, ...}], name}.
 */
export function describeForm(schema, kind, fields, ctx = {}) {
  const rec = schema.records[kind];
  if (!rec) throw new Error(`the schema has no record kind ${kind}`);
  const c = { ...ctx, schema, quantity: ctx.quantity ?? firstToken(fields.quantity) };
  const entries = {};
  const hidden = [];
  for (const [key, field] of Object.entries(rec.fields)) {
    const a = applies(key, field, fields, c);
    if (!a.show) { hidden.push(key); continue; }
    const texts = textsOf(fields[key]);
    const shown = texts.length ? texts : (field.repeatable ? [] : [""]);  // "+" adds a repeatable one
    entries[key] = {
      key, field, title: field.title || key, type: field.type, required: !!field.required,
      repeatable: !!field.repeatable, doc: field.doc || "", grammar: field.grammar || "", note: a.note,
      texts, widgets: shown.map((t) => widgetOf(field, t, c)),
      diagnostics: (ctx.diagnostics || []).filter((d) => d.field === key || (d.field || "").split(".")[0] === key),
    };
  }
  return { kind, title: rec.title || kind, doc: rec.doc || "", lines: arrange(kind, entries), hidden };
}

function firstToken(v) {
  return tokenize(Array.isArray(v) ? v[0] : v)[0];
}

/** Lines of entries: the layout hints' lines first where their fields are shown, then one field per line. */
export function arrange(kind, entries) {
  const used = new Set();
  const lines = [];
  for (const line of LAYOUT[kind] || []) {
    const row = line.filter((k) => entries[k] && !used.has(k));
    row.forEach((k) => used.add(k));
    if (row.length) lines.push({ at: Math.min(...row.map((k) => Object.keys(entries).indexOf(k))), row });
  }
  Object.keys(entries).forEach((k, at) => { if (!used.has(k)) lines.push({ at, row: [k] }); });
  return lines.sort((a, b) => a.at - b.at).map((l) => l.row.map((k) => entries[k]));
}

/**
 * The fields of an edit: for each entry whose composed texts differ from the record's, the new value (a text, a list
 * of texts for a repeatable field, or null to remove the field).  `composed` maps a key to its list of texts.
 */
export function changedFields(form, composed) {
  const out = {};
  for (const line of form.lines) {
    for (const e of line) {
      const now = (composed[e.key] || []).map((s) => s.trim()).filter((s) => s !== "");
      const before = e.texts.map((s) => s.trim());
      if (now.join("\n") === before.join("\n")) continue;
      out[e.key] = now.length === 0 ? null : (e.repeatable ? now : now[0]);
    }
  }
  return out;
}
