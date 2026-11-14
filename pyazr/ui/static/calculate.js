// "Calculate" (docs/upgrade/web_ui.md section 6, the first view): from any tab, one action shows every data set with
// the model and the residuals, the χ² per data set and in total, and what changed since the last calculation.
//
// The page computes no χ²: a calculation is the server's evaluation of the current text (POST /api/evaluate), and a
// change is the difference of two of the engine's numbers.  The calculations seen are page state (section 5): kept
// while the page lives, lost on reload.

/** What a calculation leaves to compare with the next: the revision, the totals and the χ² per data set. */
export function snapshotOf(evaluation, rev) {
  const sets = {};
  for (const d of evaluation?.datasets || []) sets[d.name] = { chi2: d.chi2, n: d.n };
  return { rev, chi2: evaluation?.chi2, objective: evaluation?.objective, sets };
}

/** Whether two snapshots are of the same text at the same parameters (the same revision and the same χ²). */
export function sameCalculation(a, b) {
  return Boolean(a && b) && a.rev === b.rev && a.chi2 === b.chi2 && a.objective === b.objective;
}

/** The calculations seen, with `now` added unless it is the newest already. */
export function remember(calculations, now) {
  const last = calculations[calculations.length - 1];
  return sameCalculation(last, now) ? calculations : [...calculations, now];
}

/** The calculation `now` is compared with: the newest one seen that is not `now` itself, or null. */
export function baselineOf(calculations, now) {
  for (let i = calculations.length - 1; i >= 0; i--) if (!sameCalculation(calculations[i], now)) return calculations[i];
  return null;
}

/**
 * The changes from `before` to `now`: {chi2, objective, sets: [{name, before, now, change, n}]}, the data sets
 * sorted by the size of their change (the largest first; a set only in one of the two has a null side and comes
 * first).  Each change is now − before of the engine's two numbers.
 */
export function changesOf(before, now) {
  const minus = (a, b) => (Number.isFinite(a) && Number.isFinite(b) ? a - b : null);
  const names = [...new Set([...Object.keys(now?.sets || {}), ...Object.keys(before?.sets || {})])];
  const sets = names.map((name) => {
    const a = before?.sets?.[name]?.chi2 ?? null;
    const b = now?.sets?.[name]?.chi2 ?? null;
    return { name, before: a, now: b, change: minus(b, a), n: now?.sets?.[name]?.n ?? before?.sets?.[name]?.n };
  });
  const size = (s) => (s.change === null ? Infinity : Math.abs(s.change));
  sets.sort((x, y) => size(y) - size(x));
  return { chi2: minus(now?.chi2, before?.chi2), objective: minus(now?.objective, before?.objective), sets,
    from: before?.rev, to: now?.rev };
}

/** A change as text with its sign ("+1.25", "−0.5", "0"); "—" where one side is missing. */
export function changeText(change, digits = 6) {
  if (change === null || change === undefined || !Number.isFinite(change)) return "—";
  if (change === 0) return "0";
  const p = Math.abs(change).toPrecision(digits);
  // trailing zeros of the fraction only ("1.25000" -> "1.25", "1.50000e+6" -> "1.5e+6"; "120000" stays)
  const t = p.includes("e") ? p.replace(/\.?0+e/, "e") : (p.includes(".") ? p.replace(/\.?0+$/, "") : p);
  return `${change > 0 ? "+" : "−"}${t}`;
}
