// Tables of the page (docs/upgrade/web_ui.md section 7): sortable by any column, filtered by a text, moved through
// with the keyboard.  Numbers are right-aligned with the unit in the column head; a number is formatted, never
// altered (formatNumber keeps every digit the server sent unless the column asks for a fixed number of significant
// digits for display, and the cell's title then holds the number as sent).
//
// sortRows, filterRows, formatNumber and columnHead need no document: node's test runner loads them
// (tests/pyazr/ui/js/tables.test.mjs).  renderTable builds the element in the browser.

/** A cell's value for sorting: numbers as numbers, the strings "inf"/"-inf"/"nan" of the server as such. */
export function sortKey(v) {
  if (v === null || v === undefined || v === "") return { kind: 2, v: 0 };
  if (typeof v === "number") return { kind: 0, v };
  if (v === "inf") return { kind: 0, v: Infinity };
  if (v === "-inf") return { kind: 0, v: -Infinity };
  if (v === "nan") return { kind: 1, v: 0 };
  if (typeof v === "boolean") return { kind: 0, v: v ? 1 : 0 };
  return { kind: 1, v: String(v) };
}

const collator = typeof Intl !== "undefined" ? new Intl.Collator(undefined, { numeric: true, sensitivity: "base" })
  : null;

/** Compare two cell values: numbers before text, text in natural order (s2 before s10), empty cells last. */
export function compareValues(a, b) {
  const ka = sortKey(a), kb = sortKey(b);
  if (ka.kind !== kb.kind) return ka.kind - kb.kind;
  if (ka.kind === 0) return ka.v - kb.v;
  if (ka.kind === 1) return collator ? collator.compare(ka.v, kb.v) : (ka.v < kb.v ? -1 : ka.v > kb.v ? 1 : 0);
  return 0;
}

/** The rows sorted by `key` (ascending unless `descending`); a stable sort, the input left as it is. */
export function sortRows(rows, key, descending = false) {
  if (!key) return rows.slice();
  const sign = descending ? -1 : 1;
  return rows.map((r, i) => [r, i]).sort((x, y) => {
    const a = valueOf(x[0], key), b = valueOf(y[0], key);
    const empty = (sortKey(a).kind === 2) - (sortKey(b).kind === 2);  // empty cells last in both directions
    const c = empty || sign * compareValues(a, b);
    return c !== 0 ? c : x[1] - y[1];
  }).map((p) => p[0]);
}

function valueOf(row, key) {
  return typeof key === "function" ? key(row) : row[key];
}

/**
 * The rows whose cells match every word of `text` (case-insensitive substrings); a word `col:value` matches only
 * the column whose key or title is `col`.  `columns` are the table's column descriptions.
 */
export function filterRows(rows, text, columns) {
  const words = String(text || "").trim().toLowerCase().split(/\s+/).filter(Boolean);
  if (!words.length) return rows.slice();
  return rows.filter((row) => words.every((w) => {
    const [col, val] = w.includes(":") ? w.split(/:(.*)/s) : [null, w];
    const cols = col ? columns.filter((c) => c.key.toLowerCase() === col || (c.title || "").toLowerCase() === col)
      : columns;
    return cols.some((c) => cellText(row, c).toLowerCase().includes(val));
  }));
}

/** The text a cell shows. */
export function cellText(row, column) {
  const v = column.value ? column.value(row) : row[column.key];
  if (column.format) return column.format(v, row);
  if (typeof v === "number") return formatNumber(v, column.digits);
  if (v === null || v === undefined) return "";
  if (typeof v === "boolean") return v ? "yes" : "no";
  return String(v);
}

/**
 * A number as text: as JavaScript writes it (every digit kept), or with `digits` significant digits for display.
 * Very large or small magnitudes use an exponent.  The server's "inf", "-inf", "nan" pass as they are.
 */
export function formatNumber(v, digits) {
  if (typeof v !== "number") return v === null || v === undefined ? "" : String(v);
  if (!Number.isFinite(v)) return Number.isNaN(v) ? "nan" : (v > 0 ? "inf" : "-inf");
  if (!digits) return String(v);
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e6 || a < 1e-4) return v.toExponential(digits - 1).replace(/e\+?/, "e");
  return String(Number(v.toPrecision(digits)));
}

/** A column's head: its title and, where it has one, its unit in brackets. */
export function columnHead(column) {
  return column.unit ? `${column.title} (${column.unit})` : column.title;
}

/**
 * Build a table element.
 *   columns: [{key, title, unit, numeric, digits, format(v,row), value(row), render(row) -> Node, sortable}]
 *   rows: the data; rowKey(row): a row's identity; selected: the key of the selected row
 *   onSelect(row): a row chosen by click, Enter or the arrow keys; filter: the initial filter text
 *   caption: the table's caption (read by screen readers, shown above it)
 * Returns {element, setRows(rows), focus(), filterInput}.
 */
export function renderTable(doc, { columns, rows, rowKey = (r) => r.name, selected = null, onSelect = null,
  filter = "", caption = "", withFilter = true, label = "" }) {
  const wrap = doc.createElement("div");
  wrap.className = "table-wrap";
  let sortCol = null, descending = false, data = rows.slice(), current = selected, filterText = filter;
  let input = null;
  if (withFilter) {
    input = doc.createElement("input");
    input.type = "search";
    input.className = "table-filter";
    input.placeholder = "filter (/)";
    input.setAttribute("aria-label", `filter ${label || caption || "the table"}`);
    input.value = filter;
    input.addEventListener("input", () => { filterText = input.value; draw(); });
    input.addEventListener("keydown", (e) => {
      if (e.key === "ArrowDown") { e.preventDefault(); move(1); body.focus(); }
    });
    wrap.appendChild(input);
  }
  const table = doc.createElement("table");
  table.className = "grid";
  if (caption) {
    const cap = doc.createElement("caption");
    cap.textContent = caption;
    table.appendChild(cap);
  }
  const head = doc.createElement("thead");
  const hr = doc.createElement("tr");
  columns.forEach((c) => {
    const th = doc.createElement("th");
    th.textContent = columnHead(c);
    if (c.numeric) th.className = "num";
    if (c.doc) th.title = c.doc;
    if (c.sortable !== false) {
      th.tabIndex = 0;
      th.setAttribute("aria-sort", "none");
      const sort = () => {
        descending = sortCol === c ? !descending : false;
        sortCol = c;
        head.querySelectorAll("th").forEach((h) => h.setAttribute("aria-sort", "none"));
        th.setAttribute("aria-sort", descending ? "descending" : "ascending");
        draw();
      };
      th.addEventListener("click", sort);
      th.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.key === " ") { e.preventDefault(); sort(); }
      });
    }
    hr.appendChild(th);
  });
  head.appendChild(hr);
  table.appendChild(head);
  const body = doc.createElement("tbody");
  body.tabIndex = 0;
  body.setAttribute("aria-label", `${label || caption || "rows"}: arrow keys move, Enter opens`);
  table.appendChild(body);
  wrap.appendChild(table);
  const count = doc.createElement("div");
  count.className = "table-count";
  wrap.appendChild(count);

  let shown = [];
  function draw() {
    const key = sortCol ? (sortCol.value || sortCol.key) : null;
    shown = sortRows(filterRows(data, filterText, columns), key, descending);
    body.textContent = "";
    for (const row of shown) {
      const tr = doc.createElement("tr");
      const k = rowKey(row);
      tr.dataset.key = k;
      if (k === current) tr.className = "selected";
      if (row._class) tr.classList.add(row._class);
      for (const c of columns) {
        const td = doc.createElement("td");
        if (c.numeric) td.className = "num";
        if (c.wrap) td.className = "wrap";
        if (c.render) {
          const node = c.render(row);
          if (node) td.appendChild(node);
        } else {
          td.textContent = cellText(row, c);
          const raw = c.value ? c.value(row) : row[c.key];
          if (c.digits && typeof raw === "number") td.title = String(raw);
        }
        tr.appendChild(td);
      }
      if (onSelect) tr.addEventListener("click", (e) => {
        if (e.target.closest("button, input, select, a")) return;
        choose(row);
      });
      body.appendChild(tr);
    }
    count.textContent = shown.length === data.length ? `${data.length} rows` : `${shown.length} of ${data.length} rows`;
  }
  function choose(row) {
    current = rowKey(row);
    body.querySelectorAll("tr").forEach((tr) => tr.classList.toggle("selected", tr.dataset.key === current));
    if (onSelect) onSelect(row);
  }
  function move(step) {
    if (!shown.length) return;
    let i = shown.findIndex((r) => rowKey(r) === current);
    i = i < 0 ? 0 : Math.max(0, Math.min(shown.length - 1, i + step));
    current = rowKey(shown[i]);
    body.querySelectorAll("tr").forEach((tr) => {
      const on = tr.dataset.key === current;
      tr.classList.toggle("selected", on);
      if (on && tr.scrollIntoView) tr.scrollIntoView({ block: "nearest" });
    });
  }
  body.addEventListener("keydown", (e) => {
    if (e.key === "ArrowDown" || e.key === "j") { e.preventDefault(); move(1); }
    else if (e.key === "ArrowUp" || e.key === "k") { e.preventDefault(); move(-1); }
    else if (e.key === "Enter") {
      const row = shown.find((r) => rowKey(r) === current);
      if (row) { e.preventDefault(); choose(row); }
    }
  });
  draw();
  return {
    element: wrap,
    filterInput: input,
    setRows(r) { data = r.slice(); draw(); },
    focus() { body.focus(); },
  };
}
