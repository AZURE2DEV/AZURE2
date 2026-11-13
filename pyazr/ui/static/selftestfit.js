// The self-test's walk through the Fit tab ("#selftest?fit=1"), as a person would do it with the buttons: a recipe
// fit started with Run, its progress seen in the document, Stop, Accept into project, then the header's Undo.  It
// records the numbers the document showed (the engine's, as sent) for tests/pyazr/ui/smoke_test.py.

import { isFinal } from "./fitting.js";

const tick = (ms) => new Promise((r) => setTimeout(r, ms));

/** Wait until test() is true (polling the document), at most `tries` × 50 ms of the page's time. */
async function until(test, tries = 1200) {
  for (let i = 0; i < tries && !test(); i += 1) await tick(50);
  return test();
}

export async function fitCheck(ctx, draw, join) {
  const out = {};
  const q = (sel) => document.querySelector(sel);
  location.hash = join("#/fit");
  await tick(0);
  await draw();
  if (!ctx.evaluation) { out.skipped = "the engine does not evaluate this project"; return out; }
  out.chi2_before = ctx.evaluation.chi2;
  out.rev_before = ctx.project.rev;
  const kind = q('select[aria-label="job kind"]');
  const method = q('select[aria-label="method"]');
  out.method = method?.value;
  // a long iteration limit: the walk stops the fit itself
  const limit = q('input[aria-label="max_iterations"], input[aria-label="max_calls"]');
  if (limit) { limit.value = "100000"; limit.dispatchEvent(new Event("change")); }
  out.kind = kind?.value;
  q('button[data-action="run"]').click();
  out.progress_seen = await until(() => Number(q(".job-progress")?.dataset.calls || 0) >= 3);
  out.calls_seen = Number(q(".job-progress")?.dataset.calls || 0);
  out.plot = q(".progress-plot")?.querySelector(".plot-area, .plot-fallback") ? "drawn" : null;
  out.log_lines = (q("pre.job-log")?.textContent || "").split("\n").length;
  q('button[data-action="stop"]').click();
  out.ended = await until(() => isFinal(q(".job-state")?.dataset.state) && Boolean(q(".job-result")));
  out.state = q(".job-state")?.dataset.state;
  out.chi2_fit_shown = await until(() => Boolean(q("[data-chi2-fit]")));
  out.chi2_fit = Number(q("[data-chi2-fit]")?.dataset.chi2Fit);
  q('button[data-action="accept"]').click();
  out.accepted = await until(() => Boolean(q("[data-chi2-written]")));
  out.chi2_written = Number(q("[data-chi2-written]")?.dataset.chi2Written);
  out.rev_accepted = ctx.project.rev;
  out.chi2_header_after = ctx.evaluation?.chi2;
  ctx.undoButton.click();
  out.undone = await until(() => ctx.statusLine.textContent.startsWith("undo: done"));
  out.chi2_back = ctx.evaluation?.chi2;
  out.rev_back = ctx.project.rev;
  return out;
}
