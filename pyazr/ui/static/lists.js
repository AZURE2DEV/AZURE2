// The lists of records in words (docs/upgrade/web_ui.md section 11.14): one list per record kind, each under the
// kind's title, its column heads the engine's and pyazr's words.  A kind the engine describes (a data set: the
// records `GET /api/project/describe` lists) shows its description's entries that pyazr's table marks for the list
// (kinds()["description"]: reaction, observable, file, range, points, ...); any other kind shows the fields its
// records set (but those the schema marks rare), headed by the schema's titles.  A record that other records name
// says which (`GET /api/project/references`).  No column is written here per kind: a new record kind gets its list
// with no line of this file.
//
// Every function but none of the elements needs a document: node's test runner loads them
// (tests/pyazr/ui/js/lists.test.mjs).

import { formatNumber } from "./tables.js";

/** The entries of a description as text: a range as "lo – hi" (one number when both are equal), a list joined,
 * a boolean as yes / no. */
export function describedText(v) {
  if (v === null || v === undefined || v === "") return "";
  if (Array.isArray(v)) {
    if (v.length === 2 && v.every((x) => typeof x === "number")) {
      return v[0] === v[1] ? formatNumber(v[0]) : `${formatNumber(v[0])} – ${formatNumber(v[1])}`;
    }
    return v.map((x) => (typeof x === "object" ? JSON.stringify(x) : String(x))).join("; ");
  }
  if (typeof v === "boolean") return v ? "yes" : "no";
  return v;
}

/** The text of a record's field in a list: a repeated field's values joined with "; ". */
export const fieldText = (v) => (Array.isArray(v) ? v.join("; ") : v ?? "");

/** For each record name, the records that name it, as "record (field)" texts: from GET /api/project/references. */
export function usersOf(references) {
  const users = {};
  for (const r of references || []) {
    (users[r.name] = users[r.name] || []).push(r.field ? `${r.record} (${r.field})` : r.record);
  }
  return users;
}

/** The normalization of a described record as the list says it: the value in use (the evaluation's), whether it is
 * fitted (the status of the description's parameter of kind "norm") and the shared record it comes from, if any. */
export function normText(description, normUsed) {
  const p = (description?.parameters || []).find((x) => x.kind === "norm");
  const value = typeof normUsed === "number" ? formatNumber(normUsed, 6) : "";
  if (!p) return value;
  const shared = !p.name.startsWith(`${description.name}.`) ? p.name.split(".")[0] : "";
  return [value, p.status === "fit" ? "fitted" : p.status, shared ? `shared ${shared}` : ""].filter(Boolean)
    .join(" · ");
}

/**
 * The columns of the list of one record kind, {key, title, doc, value(row), clip?, numeric?}:
 *   schema: the engine's schema; kind: the record kind; rows: its records ({name, fields, ...});
 *   described: {name: description} of the records the engine describes; descriptionFields: kinds()["description"]
 *   .fields; users: usersOf(references); evaluated: {name: the evaluation's entry}.
 */
export function listColumns({ schema, kind, rows, described = {}, descriptionFields = {}, users = {},
  evaluated = {} }) {
  const fields = schema.records[kind]?.fields || {};
  const columns = [{ key: "name", title: "Name", doc: "the record's name" }];
  const isDescribed = rows.some((r) => described[r.name]);
  const titled = (key, f) => `${f?.doc || ""}${f?.doc ? " " : ""}(${key})`;
  if (isDescribed) {
    for (const [key, f] of Object.entries(descriptionFields)) {
      if (!f.list || !rows.some((r) => describedText(described[r.name]?.[key]) !== "")) continue;
      const numeric = rows.every((r) => typeof (described[r.name]?.[key] ?? 0) === "number");
      columns.push({ key: `d.${key}`, title: f.title, unit: f.unit || "", doc: titled(key, f), numeric,
        clip: !numeric, value: (r) => describedText(described[r.name]?.[key]) });
    }
    if (fields.norm) {
      columns.push({ key: "norm", title: fields.norm.title, doc: "the normalization in use; fitted or fixed; "
        + "the shared record it comes from (norm)", clip: true,
      value: (r) => normText(described[r.name], evaluated[r.name]?.norm) });
    }
  } else {
    for (const [key, f] of Object.entries(fields)) {  // a rare field stays in the record's form
      if (f.rare || !rows.some((r) => fieldText(r.fields?.[key]) !== "")) continue;
      columns.push({ key: `f.${key}`, title: f.title || key, doc: titled(key, f), clip: true,
        value: (r) => fieldText(r.fields?.[key]) });
    }
  }
  if (rows.some((r) => evaluated[r.name])) {
    columns.push({ key: "chi2", title: "χ²", numeric: true, digits: 6, doc: "the program's χ² of the data set",
      value: (r) => evaluated[r.name]?.chi2 });
    if (!isDescribed) columns.push({ key: "n", title: "N", numeric: true, doc: "points used",
      value: (r) => evaluated[r.name]?.n });
  }
  if (rows.some((r) => users[r.name])) {
    columns.push({ key: "users", title: "Used by", clip: true, doc: "the records that name this one (field)",
      value: (r) => (users[r.name] || []).join(", ") });
  }
  return columns;
}

/** The order of a group's lists: the kinds whose records hold a quantity (data sets, extrapolations, model data)
 * first, then the records they name, each in the schema's order. */
export function listKinds(schema, kinds) {
  const holds = (k) => Boolean(schema.records[k]?.fields?.quantity);
  return [...kinds.filter(holds), ...kinds.filter((k) => !holds(k))];
}

/** The columns of the list of particle pairs from pyazr's table of a pair's fields (kinds()["pairs"].fields): the
 * key, then the fields it marks for the list, in its order; a field marked rare only where a pair sets it (a value
 * other than 0 or empty). */
export function pairColumns(fields, pairs) {
  const set = (v) => v !== 0 && v !== "" && v !== null && v !== undefined;
  const columns = [];
  for (const [key, f] of Object.entries(fields || {})) {
    if (key !== "key" && !f.list) continue;
    if (f.rare && !pairs.some((p) => set(p[key]))) continue;
    const numeric = pairs.every((p) => typeof p[key] === "number");
    columns.push({ key, title: f.title || key, unit: f.unit || "", doc: `${f.doc || ""} (${key})`, numeric });
  }
  return columns;
}
