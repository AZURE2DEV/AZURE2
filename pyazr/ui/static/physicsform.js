// The Physics tab's rules without a document (docs/upgrade/web_ui.md section 6.6): a row of a function table (the
// physics table, the stopping table, the session's functions; all of one shape) becomes a form of the existing form
// generator (forms.js), the form's texts become a request, and the answer becomes plots, a status line and a CSV.
//
// Nothing here is physics: the page sends what was typed (numbers parsed from the text, or the text itself where it
// is not a number: the server refuses it with what / allowed / fix), the axis as {from, to, points, spacing} which the
// server turns into an array, and draws the numbers it answers.  A new row of a table needs no line here; a new kind
// of argument is one row of ARGUMENT_KINDS (tests/pyazr/ui/js/physics.test.mjs fails until it exists).

import { describeForm, tokenize } from "./forms.js";
import { csvNumber, seriesStyle } from "./plots.js";

/** The server's bound of elements per request (pyazr/ui/physics.py MAX_ELEMENTS; the node test compares them). */
export const MAX_ELEMENTS = 20000;
/** The most lines (combinations of several values) one plot draws; the points of a first plot. */
export const MAX_LINES = 24;
export const DEFAULT_POINTS = 200;
/** The keys of the axis's points and spacing in the form (not names an argument can have). */
export const POINTS = "@points";
export const SPACING = "@spacing";

const numberOrText = (t) => (t !== "" && Number.isFinite(Number(t)) ? Number(t) : t);
const several = "several values separated by blanks give one line each";

/**
 * The kinds of argument: per kind, the field of the form (a value type of forms.js's WIDGETS) and how its text becomes
 * the values of a request (a list: one line per value).  `ctx.pairs`: the project's pairs ({key, label}).
 */
export const ARGUMENT_KINDS = {
  real: { field: (a) => ({ type: "number_list", doc: `${a.doc} (${a.domain}); ${several}` }),
    values: (t) => tokenize(t).map(numberOrText) },
  integer: { field: (a) => ({ type: "integer", doc: `${a.doc} (${a.domain}); ${several}` }),
    values: (t) => tokenize(t).map(numberOrText) },
  choice: { field: (a) => ({ type: "enum", default: a.default,
    values: [...a.choices, { value: a.choices.join(" "), label: "all (one line each)" }] }),
  values: (t) => tokenize(t) },
  flag: { field: (a) => ({ type: "enum", default: a.default,
    values: [{ value: "0", label: "0 (no)" }, { value: "1", label: "1 (yes)" }, { value: "0 1", label: "both" }] }),
  values: (t) => tokenize(t).map(numberOrText) },
  text: { field: (a) => ({ type: "text", doc: `${a.doc} (one text)` }),
    values: (t) => [String(t ?? "").trim()].filter((v) => v !== "") },
  pair: { field: (a, ctx) => ({ type: "enum", values: (ctx.pairs || []).map((p) => ({ value: String(p.key),
    label: `${p.key}: ${p.label}` })) }),
  values: (t) => tokenize(t).map(numberOrText) },
};

/** The arguments a plot may run over: those with a sweep (the row's axis first). */
export function axisChoices(row) {
  const withSweep = row.arguments.filter((a) => a.sweep);
  return [...withSweep.filter((a) => a.name === row.axis), ...withSweep.filter((a) => a.name !== row.axis)];
}

const argOf = (row, name) => row.arguments.find((a) => a.name === name);
const numberText = (x) => (x === null || x === undefined ? "" : String(x));
const isInteger = (a) => a.kind === "integer";

/** The form's start texts: each argument's default or example, the axis from its sweep with DEFAULT_POINTS points. */
export function startTexts(row, axis, ctx = {}) {
  const texts = {};
  for (const a of row.arguments) {
    if (a.name === axis) continue;
    if (a.kind === "pair") texts[a.name] = String((ctx.pairs || [])[0]?.key ?? "");
    else if (a.default !== null && a.default !== undefined) {
      texts[a.name] = a.kind === "flag" ? (a.default ? "1" : "0") : numberText(a.default);
    } else texts[a.name] = numberText(a.example);
  }
  const w = argOf(row, axis)?.sweep || { lower: "", upper: "", spacing: "linear" };
  texts[axis] = `${numberText(w.lower)} ${numberText(w.upper)}`;
  texts[POINTS] = String(DEFAULT_POINTS);
  texts[SPACING] = isInteger(argOf(row, axis) || {}) ? "integer" : w.spacing;
  return texts;
}

/** The row as a record kind of a schema of its own, for describeForm: one field per argument, the axis a range with
 * its points and spacing.  Returns {schema, kind}. */
export function formSchema(row, axis, ctx = {}) {
  const fields = {};
  for (const a of row.arguments) {
    if (a.name === axis) {
      fields[a.name] = { type: "range", title: `${a.name} (axis)`,
        units: a.unit ? [a.unit] : [],
        doc: `${a.doc} (${a.domain}): the plot's horizontal axis, from and to`, required: true };
      if (!isInteger(a)) {
        fields[POINTS] = { type: "integer", title: "points",
          doc: `points of the axis; at most ${MAX_ELEMENTS} elements per request (points × lines)` };
      }
      fields[SPACING] = { type: "enum", title: "spacing", values: isInteger(a) ? ["integer"] : ["linear", "log"],
        doc: "linear, log (from and to > 0), or every whole number from to to" };
      continue;
    }
    const kind = ARGUMENT_KINDS[a.kind];
    const field = kind ? kind.field(a, ctx) : { type: "text" };
    fields[a.name] = { title: a.unit ? `${a.name} (${a.unit})` : a.name, doc: `${a.doc} (${a.domain})`,
      required: a.default === null || a.default === undefined, ...field };
  }
  const kind = `${row.name}()`;
  return { kind, schema: { records: { [kind]: { title: row.name, doc: row.summary, fields } } } };
}

/** The form of a row as data (forms.js's describeForm). */
export function rowForm(row, axis, texts, ctx = {}) {
  const { schema, kind } = formSchema(row, axis, ctx);
  return describeForm(schema, kind, texts, ctx);
}

/** Wrap v into n levels of arrays. */
const nest = (v, n) => (n <= 0 ? v : nest([v], n - 1));

/**
 * The request of the form's texts: {arguments, axis, lines: [{label, values}], dims: [names], error}.  Each argument
 * with several values is a broadcast dimension of its own, in the order of the row; the axis is the last.  An empty
 * text leaves the argument to its default (not sent).
 */
export function requestOf(row, axis, texts) {
  const args = {};
  const lists = [];
  for (const a of row.arguments) {
    if (a.name === axis) continue;
    const kind = ARGUMENT_KINDS[a.kind];
    const values = kind ? kind.values(texts[a.name] ?? "") : [texts[a.name] ?? ""].filter((t) => t !== "");
    if (!values.length) continue;
    if (values.length === 1) args[a.name] = values[0];
    else lists.push({ name: a.name, values });
  }
  lists.forEach((l, j) => { args[l.name] = l.values.map((v) => nest(v, lists.length - j)); });
  const [from, to] = tokenize(texts[axis] ?? "");
  const spacing = texts[SPACING] || "linear";
  args[axis] = { from: numberOrText(from ?? ""), to: numberOrText(to ?? ""), spacing };
  if (spacing !== "integer") args[axis].points = numberOrText(texts[POINTS] ?? "");
  const lines = combinations(lists);
  const error = lines.length > MAX_LINES ? `${lines.length} lines: at most ${MAX_LINES} in one plot` : "";
  return { arguments: args, axis, lines, dims: lists.map((l) => l.name), error };
}

/** Every combination of the values of the lists, in the order of the broadcast (the last list varies fastest). */
function combinations(lists) {
  let out = [{ label: "", values: {} }];
  for (const l of lists) {
    out = out.flatMap((c) => l.values.map((v) => ({ label: [c.label, `${l.name} = ${v}`].filter(Boolean).join(", "),
      values: { ...c.values, [l.name]: v } })));
  }
  return out;
}

/** A nested list (or a number) in C order as a flat array. */
export const flat = (v) => (Array.isArray(v) ? v.flatMap(flat) : [v]);

/** A returned value as a number for a plot: null where it is not a finite number (the server sends "nan"). */
const plotted = (v) => (typeof v === "number" && Number.isFinite(v) ? v : null);

/** The status line of an answer: "N of M elements computed", then each other status's words with its count. */
export function statusLine(answer) {
  const codes = flat(answer.status);
  const counts = new Map();
  for (const c of codes) counts.set(c, (counts.get(c) || 0) + 1);
  const ok = counts.get(0) || 0;
  const rest = [...counts.entries()].filter(([c]) => c !== 0).sort((a, b) => b[1] - a[1])
    .map(([c, n]) => `${answer.status_words?.[String(c)] ?? `status ${c}`}: ${n}`);
  return { computed: ok, elements: codes.length,
    text: [`${ok} of ${codes.length} elements computed`, ...rest].join("; ") };
}

/** The units of the returned fields in the order of the row, each with its fields: [{unit, fields}]. */
export function unitGroups(row, shown) {
  const out = [];
  for (const r of row.returns) {
    if (shown && !shown.has(r.name)) continue;
    let g = out.find((x) => x.unit === r.unit);
    if (!g) out.push(g = { unit: r.unit, fields: [] });
    g.fields.push(r.name);
  }
  return out;
}

/**
 * The plots of an answer: one per unit of the fields shown, one trace per field and line; an element whose status is
 * not 0 is null (a gap), never a number.  `view`: {shown (Set of fields), logx, logy}; `request`: requestOf's.
 */
export function physicsPlots(row, answer, request, view = {}) {
  const status = flat(answer.status);
  const x = flat(answer.arguments[request.axis] ?? []);
  const n = x.length;
  const axisArg = argOf(row, request.axis) || { name: request.axis, unit: "" };
  const lines = request.lines.length ? request.lines : [{ label: "", values: {} }];
  const gaps = status.filter((c) => c !== 0).length;
  const plots = unitGroups(row, view.shown).map((g) => {
    const traces = [];
    g.fields.forEach((f) => {
      const y = flat(answer.returns[f]);
      lines.forEach((line, li) => {
        const ys = [];
        for (let t = 0; t < n; t += 1) {
          const k = li * n + t;
          ys.push(status[k] === 0 ? plotted(y[k]) : null);
        }
        traces.push({ role: "line", name: line.label ? `${f}, ${line.label}` : f, x, y: ys,
          color: seriesStyle(traces.length).color });
      });
    });
    return { title: `${row.name}: ${g.fields.join(", ")} against ${axisArg.name}`,
      x: { title: axisArg.name, unit: axisArg.unit, log: !!view.logx },
      y: { title: g.fields.join(", "), unit: g.unit, log: !!view.logy }, traces, note: "" };
  });
  const notes = [];
  if (gaps) notes.push(`${gaps} elements not computed are gaps (see the status line)`);
  if (view.logy) {
    const dropped = plots.reduce((s, p) => s + p.traces.reduce((u, t) => u + t.y.filter((v) => v !== null
      && v <= 0).length, 0), 0);
    if (dropped) notes.push(`${dropped} values ≤ 0 are not drawn on the log y axis`);
  }
  if (view.logx && x.some((v) => !(v > 0))) notes.push("points at x ≤ 0 are not drawn on the log x axis");
  plots.forEach((p) => { p.note = notes.join("; "); });
  return plots;
}

/** The value of a (nested) argument at multi-index `idx` of the broadcast shape (numpy's rules: aligned right). */
export function valueAt(value, idx) {
  if (!Array.isArray(value)) return value;
  let depth = 0;
  for (let v = value; Array.isArray(v); v = v[0]) depth += 1;
  let v = value;
  for (const i of idx.slice(idx.length - depth)) v = v.length === 1 ? v[0] : v[i];
  return v;
}

/** The multi-index of flat index k in `shape` (C order). */
function unravel(k, shape) {
  const idx = new Array(shape.length);
  for (let d = shape.length - 1; d >= 0; d -= 1) { idx[d] = k % shape[d]; k = Math.floor(k / shape[d]); }
  return idx;
}

/** The CSV of an answer: one line per element, the arguments at it, every returned field (17 significant digits,
 * "nan" where not computed) and the status; units in the header.  `notes`: lines written first, each after "# ". */
export function answerCsv(row, answer, notes = []) {
  const quote = (s) => (/[",\n]/.test(String(s)) ? `"${String(s).replace(/"/g, '""')}"` : String(s));
  const head = (name, unit) => (unit ? `${name} (${unit})` : name);
  const columns = [...row.arguments.map((a) => head(a.name, a.unit)),
    ...row.returns.map((r) => head(r.name, r.unit)), "status"];
  const shape = answer.shape || [];
  const fields = row.returns.map((r) => flat(answer.returns[r.name]));
  const status = flat(answer.status);
  const out = [...notes.map((t) => `# ${t}`), columns.map(quote).join(",")];
  for (let k = 0; k < status.length; k += 1) {
    const idx = unravel(k, shape);
    const cells = [...row.arguments.map((a) => valueAt(answer.arguments[a.name], idx)),
      ...fields.map((f) => f[k]), status[k]];
    out.push(cells.map((v) => quote(csvNumber(v === null ? "" : v))).join(","));
  }
  return `${out.join("\n")}\n`;
}

/** The view of a row from hash parameters: texts (v.<argument>), axis, the fields shown, logx, logy, table. */
export function applyParams(row, state, params, ctx = {}) {
  const axes = axisChoices(row).map((a) => a.name);
  if (params.axis && axes.includes(params.axis) && params.axis !== state.axis) {
    state.axis = params.axis;
    state.texts = startTexts(row, state.axis, ctx);
  }
  const t = state.texts;
  const [from, to] = tokenize(t[state.axis] ?? "");
  if (params.from !== undefined || params.to !== undefined) {
    t[state.axis] = `${params.from ?? from ?? ""} ${params.to ?? to ?? ""}`;
  }
  if (params.points !== undefined) t[POINTS] = params.points;
  if (params.spacing !== undefined) {
    t[SPACING] = params.spacing;
    state.logx = params.spacing === "log";  // the axis as spaced, unless logx is given
  }
  for (const [k, v] of Object.entries(params)) {
    if (k.startsWith("v.") && argOf(row, k.slice(2))) t[k.slice(2)] = v;
  }
  if (params.fields !== undefined) state.shown = new Set(params.fields.split(",").filter(Boolean));
  for (const f of ["logx", "logy"]) if (params[f] !== undefined) state[f] = params[f] === "1";
  if (params.plot !== undefined) state.table = params.plot === "table";
  return state;
}

/** The page state of a row: its axis, texts, fields shown and log choices (kept while the page lives). */
export function newState(row, ctx = {}) {
  const axis = axisChoices(row)[0]?.name || row.axis;
  const w = argOf(row, axis)?.sweep;
  return { axis, texts: startTexts(row, axis, ctx), shown: new Set(row.returns.map((r) => r.name)),
    logx: w?.spacing === "log", logy: false, table: false };
}

/** The rows grouped as the tables group them, filtered by a search text (name, summary, group): [{group, rows}]. */
export function groupedRows(rows, search = "") {
  const q = search.trim().toLowerCase();
  const out = [];
  for (const r of rows) {
    if (q && ![r.name, r.summary, r.group].some((s) => String(s || "").toLowerCase().includes(q))) continue;
    let g = out.find((x) => x.group === r.group);
    if (!g) out.push(g = { group: r.group, rows: [] });
    g.rows.push(r);
  }
  return out;
}
