// The page's access to the server's /api (docs/upgrade/web_ui.md sections 3.3, 4.3 and 5).
//
// The page holds no state of the project: it reads what the server answers and sends every edit with the revision
// it was made on (`rev`).  An edit refused because another client changed the project first (412) makes the page
// re-read; any other refusal comes back as an ApiError whose `view` has the three parts of section 4.3: what was
// wrong, what is allowed, how to fix it, and the engine's diagnostics (code, severity, record, field, line, text).
//
// No DOM here: node's test runner loads this module with a fetch of its own (tests/pyazr/ui/js/api.test.mjs).

/** The parts of a refusal the page shows: status, title, what, allowed, fix, diagnostics, names known. */
export function errorView(status, body) {
  const b = body && typeof body === "object" ? body : { message: String(body ?? "") };
  const diagnostics = Array.isArray(b.diagnostics) ? b.diagnostics.map(diagnosticView) : [];
  const what = b.what || b.message || (status ? `HTTP ${status}` : "no answer from the server");
  return {
    status: status ?? 0,
    title: b.error || (status ? `HTTP ${status}` : "Network error"),
    what: String(what),
    allowed: textOf(b.allowed),
    fix: textOf(b.fix),
    diagnostics,
    known: Array.isArray(b.known) ? b.known : [],
    job: b.job ?? null,
    rev: b.rev ?? b.current ?? null,
  };
}

/** One diagnostic of the engine, with the fields the page places it by (record, field, line). */
export function diagnosticView(d) {
  return {
    code: d.code ?? "", severity: d.severity ?? "error", record: d.record ?? "", field: d.field ?? "",
    line: d.line ?? null, text: d.text ?? d.message ?? "", found: d.found ?? "",
  };
}

function textOf(v) {
  if (v === null || v === undefined || v === "") return "";
  return Array.isArray(v) ? v.join(", ") : String(v);
}

/** A refusal or failure of a request; `view` is errorView's result. */
export class ApiError extends Error {
  constructor(status, body, request) {
    const view = errorView(status, body);
    super(view.what);
    this.name = "ApiError";
    this.status = status;
    this.body = body;
    this.view = view;
    this.request = request;
  }
}

/** The methods that change the project and therefore carry the revision. */
const EDITS = new Set(["POST", "PUT", "PATCH", "DELETE"]);

/** Paths whose POST is not a project edit (no revision needed): evaluations and jobs read or run. */
const NOT_EDITS = [/^\/api\/evaluate$/, /^\/api\/transform$/, /^\/api\/jobs/];

export function isEdit(method, path) {
  return EDITS.has(method) && !NOT_EDITS.some((re) => re.test(path.split("?")[0]));
}

/** Join path segments, each encoded (level labels hold blanks and "/", record names are plain). */
export function pathOf(...parts) {
  return parts.map((p, i) => (i === 0 ? String(p) : encodeURIComponent(String(p)))).join("/");
}

/**
 * The client of one server.
 *   fetchFn: the fetch to use (the browser's; a test's in node)
 *   onRevision(rev): called when an answer carries a newer revision than the page last saw
 *   onStale(view): called after a 412 (another client edited first), before the ApiError is thrown
 *   onActivity(n): the number of requests in flight changed (the header's busy mark)
 */
export function createApi({ fetchFn = globalThis.fetch?.bind(globalThis), base = "", onRevision, onStale,
  onActivity } = {}) {
  let rev = null;
  let inFlight = 0;
  const log = [];  // every failed request (the self-test reports them): {method, path, status, what}

  function seen(r) {
    if (typeof r === "number" && r !== rev) {
      rev = r;
      if (onRevision) onRevision(r);
    }
  }

  // quiet: a refusal the page expects and shows itself (not listed among the failures the self-test reports);
  // background: a request that follows something (a job's long poll) and does not light the header's busy mark
  async function request(method, path, body, { text = false, quiet = false, background = false } = {}) {
    const init = { method, credentials: "same-origin", headers: { Accept: text ? "text/plain" : "application/json" } };
    let payload = body;
    if (isEdit(method, path)) {
      payload = Object.assign({}, body || {});
      if (payload.rev === undefined && rev !== null) payload.rev = rev;
    }
    if (payload !== undefined && method !== "GET") {
      init.headers["Content-Type"] = "application/json";
      init.body = JSON.stringify(payload);
    }
    const counted = !background;
    if (counted) { inFlight += 1; if (onActivity) onActivity(inFlight); }
    let response;
    try {
      response = await fetchFn(base + path, init);
    } catch (e) {
      log.push({ method, path, status: 0, what: String(e) });
      throw new ApiError(0, { error: "Network error", message: `${method} ${path}: ${e}` }, { method, path });
    } finally {
      if (counted) { inFlight -= 1; if (onActivity) onActivity(inFlight); }
    }
    const raw = await response.text();
    let data = raw;
    if (!text || !response.ok) {
      try { data = raw ? JSON.parse(raw) : {}; } catch { data = { message: raw }; }
    }
    if (!response.ok) {
      const err = new ApiError(response.status, data, { method, path });
      if (!quiet) log.push({ method, path, status: response.status, what: err.view.what });
      if (response.status === 412) {
        if (typeof err.view.rev === "number") rev = err.view.rev;
        if (onStale) onStale(err.view);
      }
      throw err;
    }
    if (data && typeof data === "object" && !Array.isArray(data)) seen(data.rev);
    return data;
  }

  return {
    get: (path, options) => request("GET", path, undefined, options),
    text: (path) => request("GET", path, undefined, { text: true }),
    post: (path, body) => request("POST", path, body ?? {}),
    put: (path, body) => request("PUT", path, body ?? {}),
    patch: (path, body) => request("PATCH", path, body ?? {}),
    del: (path, body) => request("DELETE", path, body ?? {}),
    request,
    /** The revision the page last saw; edits carry it. */
    get rev() { return rev; },
    set rev(r) { rev = r; },
    get failures() { return log.slice(); },
  };
}

/**
 * Watch the project's revision: the server's event stream is a later step (its route answers 501), so the page
 * polls GET /api/project every `periodMs` while it is visible.  An answer whose revision differs from the one the
 * page last saw (another client's edit) reaches onRevision through the client, which re-reads the page.
 */
export function watchRevision(api, { periodMs = 3000, isVisible = () => true } = {}) {
  async function tick() {
    if (!isVisible()) return;
    try { await api.get("/api/project"); } catch { /* a failed poll shows in the next real request */ }
  }
  const timer = setInterval(tick, periodMs);
  return () => clearInterval(timer);
}
