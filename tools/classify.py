#!/usr/bin/env python3
"""classify.py CMP REFTAG NEWTAG [basis ...] -- compare the outputs of two runs of tools/runall.sh file by file.

CMP is the results root ($WORK/cmp); the bases default to brune, no-brune and use-park.  Every file written by either
run under CMP/<basis>/<TAG>/.res-*/ (output/, checks/ and run.log of each project) is compared after this reading:
carriage returns are line breaks, commas are blanks, "=" is a token of its own, and each line is split at blanks;
lines of run.log with a progress bar (|###|), a percentage ("12%"), "Calculating" or "Progress" are dropped.
A file is *identical* when its lines of tokens are equal, or when, zipping the lines and the tokens of each line, every
pair of tokens is equal as text or as numbers (float("1.0") == float("1")), no line has a token more or less, and the
two files have the same number of lines.  Spacing and the format of a number with the same value do not count.

Output (standard output, tab-separated): one line per file that differs or exists in one run only ("only-ref",
"only-new"), with the largest relative difference over numeric tokens (|a-b| / max(|a|,|b|), sign flips left out),
the same of the absolute values, the counts of differing numbers, sign flips and non-numeric tokens, and the change in
the line count; then per basis "## <basis>: N files, M identical" and one "##   <kind>: ..." line per kind of file
(run.log, a checks/ file, param.par, chiSquared.out, ..., "observable:<suffix>") that has a difference.
Exit 2 on a usage error, else 0.
"""
import os
import re
import sys
from collections import defaultdict

ABSW = 0.0  # the largest relative difference of the absolute values in the last compare() (a sign flip counts 0)

NUM = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$')
PROGRESS = re.compile(r'(\r|\|[#= ]*\||\d+%|Calculating|Progress)')


def tokens(path, is_log):
    try:
        with open(path, 'rb') as f:
            raw = f.read().decode('utf-8', 'replace')
    except OSError:
        return None
    lines = raw.replace('\r', '\n').split('\n')
    if is_log:
        lines = [l for l in lines if not PROGRESS.search(l)]
    return [l.replace(',', ' ').replace('=', ' = ').split() for l in lines]


def compare(a, b):
    """-> (max_rel, n_num_diff, n_sign_flip, n_text_diff, dlines)"""
    worst, nnum, nsign, ntext = 0.0, 0, 0, 0
    global ABSW
    ABSW = 0.0
    for la, lb in zip(a, b):
        if la == lb:
            continue
        for ta, tb in zip(la, lb):
            if ta == tb:
                continue
            if NUM.match(ta) and NUM.match(tb):
                x, y = float(ta), float(tb)
                if x == y:
                    continue
                nnum += 1
                if max(abs(x), abs(y)) > 0:
                    ABSW = max(ABSW, abs(abs(x) - abs(y)) / max(abs(x), abs(y)))
                if x == -y:
                    nsign += 1
                    continue
                d = abs(x - y) / max(abs(x), abs(y))
                worst = max(worst, d)
            else:
                ntext += 1
        if len(la) != len(lb):
            ntext += abs(len(la) - len(lb))
    return worst, nnum, nsign, ntext, len(b) - len(a)


def kind(rel):
    base = os.path.basename(rel)
    if base == 'run.log':
        return 'run.log'
    if rel.startswith('checks/'):
        return 'checks/' + base
    if base in ('param.par', 'param.sav', 'parameters.out', 'covariance.out', 'chiSquared.out',
                'relativistic.out', 'normalizations.out'):
        return base
    if base.startswith('AZUREOut') or base.endswith('.extrap') or base.endswith('.out'):
        return 'observable:' + os.path.splitext(base)[1]
    return 'other'


def collect(root):
    files = {}
    for dp, dn, fn in os.walk(root):
        for f in fn:
            p = os.path.join(dp, f)
            rel = os.path.relpath(p, root)
            if '.res-' not in rel:
                continue
            files[rel] = p
    return files


def main():
    if len(sys.argv) < 4:
        sys.stderr.write(__doc__)
        return 2
    work, ref, new = sys.argv[1:4]
    modes = sys.argv[4:] or ['brune', 'no-brune', 'use-park']
    for m in modes:
        rroot, nroot = os.path.join(work, m, ref), os.path.join(work, m, new)
        rf, nf = collect(rroot), collect(nroot)
        allk = sorted(set(rf) | set(nf))
        summary = defaultdict(lambda: [0, 0, 0.0, 0, 0, 0.0])  # files, differing, max rel, sign flips, text diffs
        ident = 0
        for k in allk:
            # rel = <dir>/.res-<name>/<azr>/<output|checks|run.log>/...
            parts = k.split('/')
            i = [j for j, p in enumerate(parts) if p.startswith('.res-')][0]
            inner = '/'.join(parts[i + 2:])
            kd = kind(inner)
            s = summary[kd]
            s[0] += 1
            if k not in rf or k not in nf:
                s[1] += 1
                print(f"{m}\t{k}\t{'only-new' if k not in rf else 'only-ref'}")
                continue
            a = tokens(rf[k], kd == 'run.log')
            b = tokens(nf[k], kd == 'run.log')
            if a == b:
                ident += 1
                continue
            w, nn, ns, nt, dl = compare(a, b)
            if nn == 0 and nt == 0 and dl == 0:
                ident += 1
                continue
            s[1] += 1
            s[2] = max(s[2], w)
            s[3] += ns
            s[4] += nt
            s[5] = max(s[5], ABSW)
            print(f"{m}\t{k}\tabsrel={ABSW:.3g}\tmaxrel={w:.3g}\tnumdiff={nn}\tsignflip={ns}\ttextdiff={nt}\tdlines={dl}")
        print(f"## {m}: {len(allk)} files, {ident} identical")
        for kd in sorted(summary):
            s = summary[kd]
            if s[1]:
                print(f"##   {kd}: {s[1]}/{s[0]} differ, max rel {s[2]:.3g}, max |abs| rel {s[5]:.3g}, sign flips {s[3]}, text diffs {s[4]}")


if __name__ == '__main__':
    sys.exit(main())
