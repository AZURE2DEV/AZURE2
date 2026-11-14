// Records of the project in the page: the panel of one record (its generated form, the engine's description beside
// it, the engine's messages under the fields, duplicate / rename / delete) and the form of a new record.  Used by the
// Data tab (every section of the schema) and the Nucleus tab (the sections it claims, NUCLEUS_SECTIONS).

import { describeForm } from "../forms.js";
import { renderForm } from "../formview.js";
import { button, diagnosticLine, h } from "../ui.js";
import { pathOf } from "../api.js";

/** The sections of the project file the Nucleus tab shows (by the schema's section names); the Data tab shows all
 * the others, so a new section of the schema becomes a group of the Data tab with no edit here. */
export const NUCLEUS_SECTIONS = ["potentials", "gammaDecays"];

/** The tab that shows a record of `kind`: "nucleus" for the kinds of NUCLEUS_SECTIONS, "data" otherwise. */
export function tabOfKind(schema, kind) {
  const sections = schema.records[kind]?.sections || [];
  return sections.some((s) => NUCLEUS_SECTIONS.includes(s)) ? "nucleus" : "data";
}

/** The hash route of a record. */
export function recordRoute(schema, record) {
  return `#/${tabOfKind(schema, record.kind)}/${encodeURIComponent(record.name)}`;
}

/** The context the forms need: the schema, the records (for references), pair labels, spins, diagnostics. */
export async function formContext(ctx, records, diagnostics = []) {
  if (!ctx.pairs) ctx.pairs = await ctx.api.get("/api/project/pairs");
  const labels = ctx.pairs.map((p) => p.label);
  const spins = ["0", "1/2", "1", "3/2", "2", "5/2", "3", "7/2", "4"];
  return {
    schema: ctx.schema, records, pairLabels: labels, spins,
    spinParities: spins.flatMap((j) => [`${j}+`, `${j}-`]), diagnostics,
    foldRare: !ctx.expert, showMeaning: !ctx.expert,
  };
}

/** The panel of one record: header with actions, form beside the engine's description, the record's messages. */
export async function recordPanel(ctx, record, records, problems, { onGone } = {}) {
  const mine = problems.filter((d) => d.record === record.name);
  const fctx = await formContext(ctx, records, mine);
  const form = describeForm(ctx.schema, record.kind, record.fields, fctx);
  const panel = h("section.panel.record", { "aria-label": `${record.kind} ${record.name}` });
  const msg = h("div.messages");
  const actionsBox = h("span.inline-action");
  const head = h("div.record-head", {}, h("h2", {}, h("span.kind", {}, record.kind), " ", record.name),
    h("span.hint", {}, `line ${record.line} · ${form.doc}`), h("span.spacer"),
    button("Duplicate", "A copy of this record under a new name", () => ask("Duplicate as", `${record.name}_copy`,
      (name) => ctx.api.post("/api/project/records", { kind: record.kind, name, fields: record.fields }),
      (name) => ctx.navigate(recordRoute(ctx.schema, { kind: record.kind, name })))),
    button("Rename", "Rename the record and every reference to it", () => ask("New name", record.name,
      (name) => ctx.api.post(`${pathOf("/api/project/records", record.name)}/rename`, { new: name }),
      (name) => ctx.navigate(recordRoute(ctx.schema, { kind: record.kind, name })))),
    button("Delete", "Remove the record (refused while another refers to it)", () => confirmDelete()),
    actionsBox);
  panel.append(head, msg);
  for (const d of mine.filter((x) => !x.field)) msg.appendChild(diagnosticLine(d));

  const formEl = renderForm(document, form, fctx, {
    onApply: (changed) => ctx.edit(`set ${record.name}`, () => ctx.api.patch(pathOf("/api/project/records",
      record.name), { fields: changed }), msg),
  });
  formEl.onrevert = () => ctx.navigate(location.hash);
  const described = h("aside.describe", { "aria-label": "the engine's description" }, h("h3", {}, "The engine reads"));
  panel.appendChild(h("div.form-and-description", {}, formEl, described));
  ctx.api.get(pathOf("/api/project/describe", record.name), { quiet: true })
    .then((d) => described.appendChild(descriptionList(d)))
    .catch((e) => described.appendChild(h("p.hint", {}, e.view ? e.view.what : String(e))));
  return panel;

  function ask(label, initial, run, then) {
    actionsBox.textContent = "";
    const input = h("input", { type: "text", value: initial, "aria-label": label, size: 18 });
    const go = () => ctx.edit(`${label.toLowerCase()} ${record.name}`, () => run(input.value.trim()), msg)
      .then(() => then(input.value.trim()), () => {});
    input.addEventListener("keydown", (e) => {
      if (e.key === "Enter") { e.preventDefault(); go(); } else if (e.key === "Escape") actionsBox.textContent = "";
    });
    actionsBox.append(h("span.ctl-label", {}, label), input, button("OK", "Enter", go),
      button("Cancel", "Esc", () => { actionsBox.textContent = ""; }));
    input.focus();
    input.select();
  }
  function confirmDelete() {
    actionsBox.textContent = "";
    const remove = () => ctx.edit(`delete ${record.name}`,
      () => ctx.api.del(pathOf("/api/project/records", record.name)), msg).then(() => onGone && onGone(), () => {});
    const yes = button(`Delete ${record.name}`, "Remove it (Undo brings it back)", remove, "danger");
    actionsBox.append(yes, button("Cancel", "Keep it", () => { actionsBox.textContent = ""; }));
    yes.focus();
  }
}

/** The engine's description of a record (a JSON object) as a definition list; nested values compactly. */
export function descriptionList(d) {
  const dl = h("dl.description");
  for (const [k, v] of Object.entries(d || {})) {
    if (v === "" || v === null || (Array.isArray(v) && !v.length)) continue;
    const path = typeof v === "string" && isPath(v);
    dl.append(h("dt", {}, k.replace(/_/g, " ")), h(path ? "dd.path" : "dd", {}, path ? pathElement(v) : compact(v)));
  }
  return dl;
}

/** A text that is a path: one word with a "/" in it. */
export const isPath = (v) => /\//.test(v) && !/\s/.test(v);

/** A path that breaks across lines after a "/" only (never inside a file name or its extension). */
export function pathElement(text) {
  const span = h("span.path");
  text.split("/").forEach((part, i, all) => {
    span.append(part + (i < all.length - 1 ? "/" : ""));
    if (i < all.length - 1) span.appendChild(document.createElement("wbr"));
  });
  return span;
}

function compact(v) {
  if (Array.isArray(v)) {
    if (v.every((x) => typeof x !== "object" || x === null)) return v.join(", ");
    return h("ul", {}, v.map((x) => h("li", {}, compact(x))));
  }
  if (v && typeof v === "object") {
    return Object.entries(v).map(([k, x]) => `${k} ${typeof x === "object" ? JSON.stringify(x) : x}`).join(" · ");
  }
  return String(v);
}

/** The form of a new record of one of `kinds`: a kind drop-down, a name, the kind's empty form; POST on apply. */
export async function newRecordPanel(ctx, kinds, records, { onAdded } = {}) {
  const panel = h("section.panel.record.new-record");
  const kindSelect = h("select", { "aria-label": "kind of the new record" },
    kinds.map((k) => h("option", { value: k }, `${k}: ${ctx.schema.records[k]?.title || k}`)));
  const body = h("div");
  const msg = h("div.messages");
  panel.append(h("div.record-head", {}, h("h2", {}, "New record"), kindSelect), msg, body);
  const fctx = await formContext(ctx, records);
  const draw = () => {
    body.textContent = "";
    const kind = kindSelect.value;
    const form = describeForm(ctx.schema, kind, {}, fctx);
    body.appendChild(h("p.hint", {}, form.doc));
    body.appendChild(renderForm(document, form, fctx, {
      nameInput: true, applyLabel: "Add",
      onApply: (name, fields) => ctx.edit(`add ${kind} ${name}`, () => ctx.api.post("/api/project/records",
        { kind, name, fields }), msg).then(() => onAdded && onAdded({ kind, name })),
    }));
  };
  kindSelect.addEventListener("change", draw);
  draw();
  return panel;
}
