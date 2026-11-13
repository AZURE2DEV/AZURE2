// Drawing a generated form (forms.js describes it as data) in the browser: labels, the widgets' controls, the
// engine's diagnostics under the field they concern, Apply and Revert.  Every control is a native input or select,
// so the form is reachable by keyboard (Tab between controls, Ctrl+Enter applies, Esc reverts).

import { valuesOf, widgetOf, changedFields } from "./forms.js";

let uid = 0;

/** One control described by a widget, as an element; `onChange` after each change. */
function controlElement(doc, c, onChange, ctx) {
  const wrap = doc.createElement("span");
  wrap.className = "ctl";
  if (c.label) {
    const l = doc.createElement("span");
    l.className = "ctl-label";
    l.textContent = c.label;
    wrap.appendChild(l);
  }
  let el;
  if (c.kind === "select") {
    el = doc.createElement("select");
    for (const o of c.options) {
      const opt = doc.createElement("option");
      opt.value = o.value;
      opt.textContent = o.label ?? o.value;
      if (o.title) opt.title = o.title;
      if (o.disabled) opt.disabled = true;
      el.appendChild(opt);
    }
    el.value = c.value ?? "";
  } else {
    el = doc.createElement("input");
    el.type = "text";
    el.value = c.value ?? "";
    if (c.placeholder) el.placeholder = c.placeholder;
    if (c.numeric) el.classList.add("num");
    if (c.wide) el.classList.add("wide");
    if (c.mono) el.classList.add("mono");
    if (c.list && c.list.length) {
      const id = `dl${++uid}`;
      const dl = doc.createElement("datalist");
      dl.id = id;
      for (const v of c.list) { const o = doc.createElement("option"); o.value = v; dl.appendChild(o); }
      wrap.appendChild(dl);
      el.setAttribute("list", id);
    }
    const chars = String(c.value || c.placeholder || "").length;
    el.style.width = `${Math.max(c.numeric ? 7 : 9, Math.min(c.wide ? 64 : 16, chars + 3))}ch`;
  }
  if (c.disabled) el.disabled = true;
  if (c.title) el.title = c.title;
  el.setAttribute("aria-label", c.ariaLabel || c.label || c.placeholder || c.name);
  el.dataset.name = c.name;
  el.addEventListener(c.kind === "select" ? "change" : "input", () => onChange(c, el.value));
  wrap.appendChild(el);
  if (c.hint) {
    const h = doc.createElement("span");
    h.className = "hint";
    h.textContent = c.hint;
    wrap.appendChild(h);
  }
  return wrap;
}

/**
 * The element of a form.  `form` is describeForm's result; `ctx` the context it was made with;
 * onApply(changed) gets the changed fields ({key: text | [texts] | null}); `actions` are extra buttons
 * ({label, title, run}); `nameInput`: a name box first (a new record).
 */
export function renderForm(doc, form, ctx, { onApply, actions = [], nameInput = false, applyLabel = "Apply" } = {}) {
  const root = doc.createElement("form");
  root.className = "record-form";
  root.setAttribute("aria-label", `${form.title} form`);
  root.addEventListener("submit", (e) => e.preventDefault());
  // the state of every widget instance: key -> [{widget, values}]
  const state = {};
  let nameBox = null;
  if (nameInput) {
    const line = doc.createElement("div");
    line.className = "form-line";
    const lab = doc.createElement("label");
    lab.className = "field";
    lab.innerHTML = '<span class="field-title">Name <span class="req">*</span></span>';
    nameBox = doc.createElement("input");
    nameBox.type = "text";
    nameBox.placeholder = "name";
    nameBox.setAttribute("aria-label", "name of the new record");
    lab.appendChild(nameBox);
    line.appendChild(lab);
    root.appendChild(line);
  }
  for (const line of form.lines) {
    const row = doc.createElement("div");
    row.className = "form-line";
    for (const entry of line) row.appendChild(fieldElement(doc, entry, ctx, state));
    root.appendChild(row);
  }
  if (form.hidden.length) {
    const h = doc.createElement("div");
    h.className = "form-hidden";
    h.textContent = `Not shown (they do not apply to this record as it stands): ${form.hidden.join(", ")}`;
    root.appendChild(h);
  }
  const bar = doc.createElement("div");
  bar.className = "form-actions";
  const apply = button(doc, applyLabel, "Send the changed fields to the engine (Ctrl+Enter)", () => submit());
  apply.classList.add("primary");
  bar.appendChild(apply);
  bar.appendChild(button(doc, "Revert", "Show the record as the server has it (Esc)", () => revert()));
  for (const a of actions) bar.appendChild(button(doc, a.label, a.title || "", a.run));
  const status = doc.createElement("span");
  status.className = "form-status";
  status.setAttribute("role", "status");
  bar.appendChild(status);
  root.appendChild(bar);
  root.addEventListener("keydown", (e) => {
    if (e.key === "Enter" && (e.ctrlKey || e.metaKey)) { e.preventDefault(); submit(); }
    else if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); revert(); }
  });

  function composed() {
    const out = {};
    for (const [key, list] of Object.entries(state)) out[key] = list.map((s) => s.widget.compose(s.values));
    return out;
  }
  function submit() {
    const changed = changedFields(form, composed());
    if (!nameInput && !Object.keys(changed).length) { status.textContent = "nothing changed"; return; }
    status.textContent = "sending…";
    const args = nameInput ? [nameBox.value.trim(), changed] : [changed];
    Promise.resolve(onApply && onApply(...args)).then(() => { status.textContent = ""; },
      (e) => { status.textContent = e?.view?.what ? "refused (see the messages)" : String(e); });
  }
  function revert() {
    if (root.onrevert) root.onrevert();
  }
  root.setStatus = (t) => { status.textContent = t; };
  return root;
}

function button(doc, label, title, run) {
  const b = doc.createElement("button");
  b.type = "button";
  b.textContent = label;
  b.title = title;
  b.addEventListener("click", run);
  return b;
}

/** One field: its title, one line of controls per occurrence, "+" for a repeatable one, its messages. */
function fieldElement(doc, entry, ctx, state) {
  const box = doc.createElement("div");
  box.className = `field type-${entry.type}`;
  box.dataset.field = entry.key;
  const title = doc.createElement("div");
  title.className = "field-title";
  title.textContent = entry.title;
  title.title = `${entry.key}: ${entry.doc}\n${entry.grammar}`;
  if (entry.required) {
    const r = doc.createElement("span");
    r.className = "req";
    r.textContent = " *";
    r.title = "required";
    title.appendChild(r);
  }
  box.appendChild(title);
  const occurrences = doc.createElement("div");
  occurrences.className = "occurrences";
  box.appendChild(occurrences);
  state[entry.key] = [];
  const add = (widget) => {
    const slot = { widget, values: valuesOf(widget) };
    state[entry.key].push(slot);
    const line = doc.createElement("div");
    line.className = "occurrence";
    occurrences.appendChild(line);
    const draw = () => {
      line.textContent = "";
      for (const c of slot.widget.controls) {
        line.appendChild(controlElement(doc, c, (ctl, value) => {
          slot.values[ctl.name] = value;
          if (ctl.rerender) {
            slot.widget = widgetOf(entry.field, slot.widget.compose(slot.values), ctx);
            slot.values = valuesOf(slot.widget);
            draw();
          }
        }, ctx));
      }
      if (entry.repeatable) {
        line.appendChild(button(doc, "−", `remove this ${entry.title}`, () => {
          state[entry.key] = state[entry.key].filter((s) => s !== slot);
          line.remove();
        }));
      }
    };
    draw();
  };
  entry.widgets.forEach(add);
  if (entry.repeatable) {
    box.appendChild(button(doc, `+ ${entry.title}`, `add a ${entry.title}`, () => add(widgetOf(entry.field, "", ctx))));
  }
  if (entry.note) {
    const n = doc.createElement("div");
    n.className = "field-note";
    n.textContent = `applies when ${entry.note}`;
    box.appendChild(n);
  }
  for (const d of entry.diagnostics) {
    const m = doc.createElement("div");
    m.className = `diag ${d.severity}`;
    m.textContent = `${d.code}: ${d.text}`;
    box.appendChild(m);
  }
  return box;
}
