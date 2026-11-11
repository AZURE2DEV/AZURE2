#!/usr/bin/env python3
"""classify.py CMP REFTAG NEWTAG [basis ...] -- compare the outputs of two runs of tools/runall.sh file by file.

CMP is the results root ($WORK/cmp); the bases default to brune, no-brune and use-park.  Every file written by either
run under CMP/<basis>/<TAG>/.res-*/ (output/, checks/ and run.log of each project) is compared after this reading:
carriage returns are line breaks, commas are blanks, "=" is a token of its own, and each line is split at blanks;
lines of run.log with a progress bar (|###|), a percentage ("12%"), "Calculating" or "Progress" are dropped.
A file is *identical* when its lines of tokens are equal, or when, zipping the lines and the tokens of each line, every
pair of tokens is equal as text or as numbers (float("1.0") == float("1")), no line has a token more or less, and the
two files have the same number of lines.  Spacing and the format of a number with the same value do not count.

A file that cannot be read (in either run) is never identical.

Output (standard output, tab-separated): one line per file that differs, exists in one run only ("only-ref",
"only-new") or cannot be read ("unreadable-ref", "unreadable-new", "unreadable-both"), with the largest relative
difference over numeric tokens (|a-b| / max(|a|,|b|), sign flips left out), the same of the absolute values, the
counts of differing numbers, sign flips and non-numeric tokens, and the change in the line count; then per basis "## <basis>: N files, M identical" and one "##   <kind>: ..." line per kind of file
(run.log, a checks/ file, param.par, chiSquared.out, ..., "observable:<suffix>") that has a difference.
With --class lab-frame-rounding, a file of a project listed in LAB_FRAME_PROJECTS that differs only by rounding is
counted in that class (listed with "lab-frame-rounding" and added to the summary line) instead of differing: the same
lines and tokens, and every differing number within LAB_FRAME_LIMIT relative or within one unit of the last printed
digit of either token (two values that agree to 1e-14 can print one last digit apart).  The class is the comparison
of the residual in the frame of the datum (docs/upgrade/energy_grid.md section 4.4): only the projects with
lab-differential data, whose outputs move by rounding, are listed; every other project must stay identical.
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


# The lab-frame rounding class: the projects (".res-<folder>/<file>.azr" of runall.sh) whose outputs may move by
# rounding when lab-differential data are compared in the lab frame, and the limit on any output number.
LAB_FRAME_LIMIT = 1e-14
# Listed from a probe of every project of tests/ (outside tests/convert/) that loads without generated files: the
# project files with a data row whose dOmega_lab/dOmega_cm is not 1.
LAB_FRAME_PROJECTS = {
    '.res-13N/13N.azr',
    '.res-angle_mapping/angle_mapping.azr',
    '.res-channel_space/channel_space.azr',
    '.res-hybrid_potential/hybrid_potential.azr',
    '.res-identical_entrance_reaction__projects/c12c12.azr',
    '.res-identical_entrance_reaction__projects/dd_pt.azr',
    '.res-relativistic__models/he5_ang.azr',
    '.res-relativistic__models/he5_ang_rel.azr',
    '.res-relativistic__models/li7pn_ang.azr',
    '.res-relativistic__models/li7pn_ang_rel.azr',
    '.res-relativistic__models/li7pn_edge.azr',
    '.res-relativistic__models/pa_ang.azr',
    '.res-relativistic__models/pa_ang_rel.azr',
    '.res-spin_orbit_potential/spin_orbit_potential.azr',
    '.res-sum_lab_frame__models/plain_data.azr',
    '.res-sum_lab_frame__models/sum_data.azr',
    '.res-sum_lab_frame__models/sum_free.azr',
    '.res-upgrade__T1/T1.azr',
    '.res-upgrade__T1/T1_start.azr',
    '.res-upgrade__T2/T2.azr',
    '.res-upgrade__T2/T2_true.azr',
    '.res-upgrade__T3/T3.azr',
    '.res-upgrade__T3/T3_start.azr',
}


def last_place(t):
    """One unit in the last printed digit of the numeric token t ('2.19e+00' -> 0.01, '692' -> 1)."""
    m = re.match(r'^[-+]?(\d*)(?:\.(\d*))?(?:[eE]([-+]?\d+))?$', t)
    decimals = len(m.group(2) or '')
    exponent = int(m.group(3) or 0)
    return 10.0 ** (exponent - decimals)


def within_rounding(a, b, limit):
    """True when the token lines a and b differ only in numbers, each pair within `limit` relative or one unit of
    the last printed digit."""
    if len(a) != len(b):
        return False
    for la, lb in zip(a, b):
        if len(la) != len(lb):
            return False
        for ta, tb in zip(la, lb):
            if ta == tb:
                continue
            if not (NUM.match(ta) and NUM.match(tb)):
                return False
            x, y = float(ta), float(tb)
            # (1 + 1e-6): the difference of two decimal tokens one unit apart is not exactly one unit in binary
            if abs(x - y) > max(limit * max(abs(x), abs(y)), last_place(ta), last_place(tb)) * (1 + 1e-6):
                return False
    return True


def project_of(k):
    """'.res-<folder>/<file>.azr' of a relative path of collect()."""
    parts = k.split('/')
    i = [j for j, p in enumerate(parts) if p.startswith('.res-')][0]
    return '/'.join(parts[i:i + 2])


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
    args = [a for a in sys.argv[1:] if not a.startswith('--class')]
    classes = [a.split('=', 1)[1] if '=' in a else '' for a in sys.argv[1:] if a.startswith('--class')]
    if any(c != 'lab-frame-rounding' for c in classes):
        sys.stderr.write('classify.py: the only class is --class=lab-frame-rounding\n')
        return 2
    lab_frame = bool(classes)
    if len(args) < 3:
        sys.stderr.write(__doc__)
        return 2
    work, ref, new = args[0:3]
    modes = args[3:] or ['brune', 'no-brune', 'use-park']
    for m in modes:
        rroot, nroot = os.path.join(work, m, ref), os.path.join(work, m, new)
        rf, nf = collect(rroot), collect(nroot)
        allk = sorted(set(rf) | set(nf))
        summary = defaultdict(lambda: [0, 0, 0.0, 0, 0, 0.0])  # files, differing, max rel, sign flips, text diffs
        ident = 0
        rounding = 0
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
            if a is None or b is None:
                # a file that cannot be read is never identical, not even to another one that cannot be read
                s[1] += 1
                which = 'both' if a is None and b is None else ('ref' if a is None else 'new')
                print(f"{m}\t{k}\tunreadable-{which}")
                continue
            if a == b:
                ident += 1
                continue
            w, nn, ns, nt, dl = compare(a, b)
            if nn == 0 and nt == 0 and dl == 0:
                ident += 1
                continue
            if lab_frame and project_of(k) in LAB_FRAME_PROJECTS and within_rounding(a, b, LAB_FRAME_LIMIT):
                rounding += 1
                print(f"{m}\t{k}\tlab-frame-rounding\tmaxrel={w:.3g}\tnumdiff={nn}")
                continue
            s[1] += 1
            s[2] = max(s[2], w)
            s[3] += ns
            s[4] += nt
            s[5] = max(s[5], ABSW)
            print(f"{m}\t{k}\tabsrel={ABSW:.3g}\tmaxrel={w:.3g}\tnumdiff={nn}\tsignflip={ns}\ttextdiff={nt}\tdlines={dl}")
        print(f"## {m}: {len(allk)} files, {ident} identical" + (f", {rounding} lab-frame rounding" if lab_frame else ""))
        for kd in sorted(summary):
            s = summary[kd]
            if s[1]:
                print(f"##   {kd}: {s[1]}/{s[0]} differ, max rel {s[2]:.3g}, max |abs| rel {s[5]:.3g}, sign flips {s[3]}, text diffs {s[4]}")


if __name__ == '__main__':
    sys.exit(main())
