#!/usr/bin/env python3
"""classify_quadrature.py CMP REFTAG FINETAG NEWTAG [basis ...] -- classify every changed number of a new engine's
outputs by the quadrature criterion of docs/upgrade/energy_grid.md section 7.1.

CMP is the results root of tools/runall.sh ($WORK/cmp); the bases default to brune, no-brune and use-park.
REFTAG is a run of the reference engine, FINETAG a run of the same reference engine with every sub-point grid refined
k-fold (`AZURE_REFINE_QUADRATURE=k`: points per width x k, smooth step / k, tail step per distance / k), NEWTAG the
engine under test.  FINETAG may be a comma-separated list, finest first (`tab-fine8,tab-fine4`): each file is compared
with the first of them that has it, and the tag used is printed with the file.  Every file of NEWTAG and REFTAG is read as tools/classify.py reads it (tokens(): commas are blanks, "=" a
token, progress lines of run.log dropped).  A file is

  identical    when classify.compare finds no difference;
  quadrature   when it has the reference's lines and tokens, and every number that differs satisfies
                   |m_new - m_fine| <= max(|m_ref - m_fine|, 1e-10 |m_ref|)
               (the new value is at least as close to the converged value as the reference's), and no number that the
               refinement leaves unchanged (m_ref == m_fine: not folded) has moved;
               this is the strict, per-number rule;
  failure      otherwise: a structural difference, a file missing in a run, a number further from the converged value
               than the reference's, or a number of an unfolded output that moved.

A run.log is classified like any file, but its numbers are counters and settings, not observables; its class is
printed with the suffix "(log)" and counted apart.  The external-capture caches intEC.dat and intEC.dat.sig hold one
amplitude per sub-point (their layout follows the grids) and a signature: they are listed as "cache" and not compared.
A number the refinement leaves unchanged may move by RELATIVE_FLOOR of itself (the last printed digit of a folded
value that does not depend on the grid, as a constant cross section).

The per-file rule (the director's decision of 2026-11-12, docs/upgrade/energy_grid.md section 7.1, for a mesh whose
nodes are elsewhere than the reference's): a file is "quadrature" when its structure is the reference's, no number the
refinement leaves unchanged has moved, and its largest relative distance to the refined run over the folded numbers
(m_ref != m_fine), max |m_new - m_fine| / |m_fine|, is not larger than FILE_FACTOR (1.2) times the same maximum of the
reference, and the root mean square of the same relative distances over the file's folded numbers,
sqrt(mean (m_new - m_fine)^2 / m_fine^2), is not larger than RMS_FACTOR (1.03) times the reference's: the k = 8 refined
run is itself 1/64 of the reference's distance from the converged value, so a ratio of two distances to it is resolved
to about 3 % and no better.  A file whose largest new distance is at most CONVERGED (1e-7, the size by which the
Coulomb table moves a charged-channel number) passes: both engines are converged there.  A file of NAMED, a table of
projects whose difference was examined and accepted (one row each, with the reason), passes with its own factor on
both measures and is printed and counted as "named", never as "quadrature".  Both classes are printed; the per-file
class decides the summary and the exit status, the strict one is information.

Output (standard output, tab-separated): one line per file that is not identical -- basis, file, class, the counts of
numbers that moved (quadrature-conform, unfolded-moved, further-than-reference), the largest |m_new - m_fine| /
|m_ref - m_fine| of the conform ones, the first offending number, the per-file class, the largest and the root mean
square relative distances to the refined run of the new and the reference values, and the refined tag used -- then per
basis
"## <basis>: N files, I identical, Q quadrature, M named, F failure (L logs not identical); strict rule: ...".
Exit 2 on a usage error, 1 when any non-log file fails, else 0.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classify  # noqa: E402  (tools/classify.py, beside this file)

BASES = ["brune", "no-brune", "use-park"]
RELATIVE_FLOOR = 1e-10
FILE_FACTOR = 1.2  # largest relative distance: new <= FILE_FACTOR x reference
RMS_FACTOR = 1.03  # root mean square of the relative distances: new <= RMS_FACTOR x reference (resolution of k = 8)
CONVERGED = 1e-7  # a file whose largest new relative distance is no larger passes
# Projects accepted by name: path fragment -> (factor on the largest and the rms distance, reason).
NAMED = {
    ".res-hybrid_potential/": (1.25, "5 % from the converged value with its own settings in both engines; the mesh is"
                               " 1.20 x further at k = 1, 0.85 x at k = 2, equal from k = 4 (second order, same limit):"
                               " docs/upgrade/energy_grid.md section 7.1"),
}


def number(token):
    return float(token) if classify.NUM.match(token) else None


def classify_file(ref, fine, new, named_factor=None):
    """named_factor: the factor of the file's row of NAMED, if it has one (the per-file class is then "named" where
    the rule alone fails and the factor passes).
    -> (class, conform, unfolded_moved, further, worst_ratio, first_offence, per_file_class, distances), distances
    = [largest new, largest ref, rms new, rms ref] relative distance to the refined run over the folded numbers"""
    if ref is None or new is None:
        return "failure", 0, 0, 0, 0.0, "unreadable or missing", "failure", [0.0] * 4
    if ref == new:
        return "identical", 0, 0, 0, 0.0, "", "identical", [0.0] * 4
    if len(ref) != len(new) or any(len(a) != len(b) for a, b in zip(ref, new)):
        return "failure", 0, 0, 0, 0.0, "lines or tokens differ in number", "failure", [0.0] * 4
    if fine is None or len(fine) != len(ref) or any(len(a) != len(b) for a, b in zip(ref, fine)):
        return "failure", 0, 0, 0, 0.0, "the refined run's file differs in structure", "failure", [0.0] * 4
    conform = unfolded = further = 0
    worst = 0.0
    first = ""
    text_moved = False
    dist = [0.0, 0.0]  # largest |m_new - m_fine| / |m_fine| and |m_ref - m_fine| / |m_fine| over the folded numbers
    squares = [0.0, 0.0]  # sums of the squares of the same relative distances
    folded = 0  # folded numbers (m_ref != m_fine), whether or not the new value moved
    for i, (lr, lf, ln) in enumerate(zip(ref, fine, new)):
        if lr == ln and lr == lf:
            continue
        for j, (tr, tf, tn) in enumerate(zip(lr, lf, ln)):
            if tr == tn and tr == tf:
                continue
            r, f, n = number(tr), number(tf), number(tn)
            if tr == tn:
                if r is not None and f is not None and r != f:  # a folded number the new engine left in place
                    scale = abs(f) if f != 0 else abs(r)
                    folded += 1
                    squares[0] += ((r - f) / scale) ** 2
                    squares[1] += ((r - f) / scale) ** 2
                    dist[0] = max(dist[0], abs(r - f) / scale)
                    dist[1] = max(dist[1], abs(r - f) / scale)
                continue
            if r is None or n is None or f is None:
                if r == n and r is not None:
                    continue  # the same value printed otherwise
                further += 1
                text_moved = True
                first = first or f"line {i + 1} token {j + 1}: text {tr!r} -> {tn!r}"
                continue
            if r != f:
                scale = abs(f) if f != 0 else max(abs(r), abs(n))
                if scale > 0:
                    folded += 1
                    squares[0] += ((n - f) / scale) ** 2
                    squares[1] += ((r - f) / scale) ** 2
                    dist[0] = max(dist[0], abs(n - f) / scale)
                    dist[1] = max(dist[1], abs(r - f) / scale)
            if r == n:
                continue
            if r == f:
                if abs(n - r) <= RELATIVE_FLOOR * abs(r):
                    conform += 1  # the last printed digit of a value the refinement leaves unchanged
                    continue
                unfolded += 1
                first = first or f"line {i + 1} token {j + 1}: {tr} -> {tn} (unchanged by the refinement)"
                continue
            if abs(n - f) <= max(abs(r - f), RELATIVE_FLOOR * abs(r)):
                conform += 1
                worst = max(worst, abs(n - f) / abs(r - f))
            else:
                further += 1
                first = first or f"line {i + 1} token {j + 1}: ref {tr}, fine {tf}, new {tn}"
    rms = [(q / folded) ** 0.5 if folded else 0.0 for q in squares]
    dist += rms
    if unfolded or text_moved:
        per_file = "failure"
    elif dist[0] <= CONVERGED or (dist[0] <= FILE_FACTOR * dist[1] and rms[0] <= RMS_FACTOR * rms[1]):
        per_file = "quadrature"
    elif named_factor is not None and dist[0] <= named_factor * dist[1] and rms[0] <= named_factor * rms[1]:
        per_file = "named"
    else:
        per_file = "failure"
    if unfolded or further:
        return "failure", conform, unfolded, further, worst, first, per_file, dist
    return ("quadrature" if conform else "identical"), conform, 0, 0, worst, "", per_file, dist


def main(argv):
    if len(argv) < 5:
        sys.stderr.write(__doc__)
        return 2
    cmp_root, ref_tag, fine_tag, new_tag = argv[1:5]
    bases = argv[5:] or BASES
    failed = False
    for basis in bases:
        ref_files = classify.collect(os.path.join(cmp_root, basis, ref_tag))
        new_files = classify.collect(os.path.join(cmp_root, basis, new_tag))
        fine_roots = [(t, os.path.join(cmp_root, basis, t)) for t in fine_tag.split(",")]
        counts = {"identical": 0, "quadrature": 0, "named": 0, "failure": 0}
        strict = {"identical": 0, "quadrature": 0, "failure": 0}
        logs = 0
        caches = 0
        for rel in sorted(set(ref_files) | set(new_files)):
            if os.path.basename(rel).startswith("intEC.dat"):
                caches += 1
                continue
            is_log = rel.endswith("run.log")
            ref = classify.tokens(ref_files[rel], is_log) if rel in ref_files else None
            new = classify.tokens(new_files[rel], is_log) if rel in new_files else None
            fine, fine_used = None, "-"
            for t, root in fine_roots:
                if os.path.exists(os.path.join(root, rel)):
                    fine, fine_used = classify.tokens(os.path.join(root, rel), is_log), t
                    break
            cls, conform, unfolded, further, worst, first, per_file, dist = classify_file(
                ref, fine, new, next((factor for fragment, (factor, _) in NAMED.items() if fragment in rel), None))
            if is_log:
                if cls != "identical":
                    logs += 1
                    print(f"{basis}\t{rel}\t{cls} (log)\t{conform}\t{unfolded}\t{further}\t{worst:.3g}\t{first}")
                continue
            counts[per_file] += 1
            strict[cls] += 1
            if cls != "identical":
                print(f"{basis}\t{rel}\t{cls}\t{conform}\t{unfolded}\t{further}\t{worst:.3g}\t{first}\t"
                      f"file: {per_file}\tnew {dist[0]:.3g} ref {dist[1]:.3g}\trms new {dist[2]:.3g} ref {dist[3]:.3g}\t"
                      f"{fine_used}")
            failed = failed or per_file == "failure"
        n = sum(counts.values())
        print(f"## {basis}: {n} files, {counts['identical']} identical, {counts['quadrature']} quadrature, "
              f"{counts['named']} named, {counts['failure']} failure ({logs} logs not identical, {caches} caches not compared); strict rule: "
              f"{strict['identical']} identical, {strict['quadrature']} quadrature, {strict['failure']} failure")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
