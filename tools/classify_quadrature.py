#!/usr/bin/env python3
"""classify_quadrature.py CMP REFTAG FINETAG NEWTAG [basis ...] -- classify every changed number of a new engine's
outputs by the quadrature criterion of docs/upgrade/energy_grid.md section 7.1.

CMP is the results root of tools/runall.sh ($WORK/cmp); the bases default to brune, no-brune and use-park.
REFTAG is a run of the reference engine, FINETAG a run of the same reference engine with every sub-point grid refined
twofold (`AZURE_REFINE_QUADRATURE=2`: points per width x2, smooth step /2, tail step per distance /2), NEWTAG the engine
under test.  Every file of NEWTAG and REFTAG is read as tools/classify.py reads it (tokens(): commas are blanks, "=" a
token, progress lines of run.log dropped).  A file is

  identical    when classify.compare finds no difference;
  quadrature   when it has the reference's lines and tokens, and every number that differs satisfies
                   |m_new - m_fine| <= max(|m_ref - m_fine|, 1e-10 |m_ref|)
               (the new value is at least as close to the converged value as the reference's), and no number that the
               refinement leaves unchanged (m_ref == m_fine: not folded) has moved;
  failure      otherwise: a structural difference, a file missing in a run, a number further from the converged value
               than the reference's, or a number of an unfolded output that moved.

A run.log is classified like any file, but its numbers are counters and settings, not observables; its class is
printed with the suffix "(log)" and counted apart.  The external-capture caches intEC.dat and intEC.dat.sig hold one
amplitude per sub-point (their layout follows the grids) and a signature: they are listed as "cache" and not compared.
A number the refinement leaves unchanged may move by RELATIVE_FLOOR of itself (the last printed digit of a folded
value that does not depend on the grid, as a constant cross section).

Output (standard output, tab-separated): one line per file that is not identical -- basis, file, class, the counts of
numbers that moved (quadrature-conform, unfolded-moved, further-than-reference), the largest |m_new - m_fine| /
|m_ref - m_fine| of the conform ones, and the first offending number -- then per basis
"## <basis>: N files, I identical, Q quadrature, F failure (L logs not identical)".
Exit 2 on a usage error, 1 when any non-log file fails, else 0.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classify  # noqa: E402  (tools/classify.py, beside this file)

BASES = ["brune", "no-brune", "use-park"]
RELATIVE_FLOOR = 1e-10


def number(token):
    return float(token) if classify.NUM.match(token) else None


def classify_file(ref, fine, new):
    """-> (class, conform, unfolded_moved, further, worst_ratio, first_offence)"""
    if ref is None or new is None:
        return "failure", 0, 0, 0, 0.0, "unreadable or missing"
    if ref == new:
        return "identical", 0, 0, 0, 0.0, ""
    if len(ref) != len(new) or any(len(a) != len(b) for a, b in zip(ref, new)):
        return "failure", 0, 0, 0, 0.0, "lines or tokens differ in number"
    if fine is None or len(fine) != len(ref) or any(len(a) != len(b) for a, b in zip(ref, fine)):
        return "failure", 0, 0, 0, 0.0, "the refined run's file differs in structure"
    conform = unfolded = further = 0
    worst = 0.0
    first = ""
    for i, (lr, lf, ln) in enumerate(zip(ref, fine, new)):
        if lr == ln:
            continue
        for j, (tr, tf, tn) in enumerate(zip(lr, lf, ln)):
            if tr == tn:
                continue
            r, f, n = number(tr), number(tf), number(tn)
            if r is None or n is None or f is None:
                if r == n and r is not None:
                    continue  # the same value printed otherwise
                further += 1
                first = first or f"line {i + 1} token {j + 1}: text {tr!r} -> {tn!r}"
                continue
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
    if unfolded or further:
        return "failure", conform, unfolded, further, worst, first
    return ("quadrature" if conform else "identical"), conform, 0, 0, worst, ""


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
        fine_root = os.path.join(cmp_root, basis, fine_tag)
        counts = {"identical": 0, "quadrature": 0, "failure": 0}
        logs = 0
        caches = 0
        for rel in sorted(set(ref_files) | set(new_files)):
            if os.path.basename(rel).startswith("intEC.dat"):
                caches += 1
                continue
            is_log = rel.endswith("run.log")
            ref = classify.tokens(ref_files[rel], is_log) if rel in ref_files else None
            new = classify.tokens(new_files[rel], is_log) if rel in new_files else None
            fine_path = os.path.join(fine_root, rel)
            fine = classify.tokens(fine_path, is_log) if os.path.exists(fine_path) else None
            cls, conform, unfolded, further, worst, first = classify_file(ref, fine, new)
            if is_log:
                if cls != "identical":
                    logs += 1
                    print(f"{basis}\t{rel}\t{cls} (log)\t{conform}\t{unfolded}\t{further}\t{worst:.3g}\t{first}")
                continue
            counts[cls] += 1
            if cls != "identical":
                print(f"{basis}\t{rel}\t{cls}\t{conform}\t{unfolded}\t{further}\t{worst:.3g}\t{first}")
            failed = failed or cls == "failure"
        n = sum(counts.values())
        print(f"## {basis}: {n} files, {counts['identical']} identical, {counts['quadrature']} quadrature, "
              f"{counts['failure']} failure ({logs} logs not identical, {caches} caches not compared)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
