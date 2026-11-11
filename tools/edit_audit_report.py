#!/usr/bin/env python3
"""The table of the audit of the checks' edits (tools/edit_audit.sh writes the log this reads).

    python3 tools/edit_audit_report.py <out dir> [--details]

Per check: its verdict, the edits it made (calls of sed/awk/perl that were given a project, data or parameter file
and wrote a changed copy or changed it in place, and python3 calls that changed such a file), the edits that hit
nothing (output equal to input), the AZURE2 runs, the runs whose input equals an earlier run's (the same run twice),
and the runs whose outputs equal those of a run with another input (a variant that made no difference).  An edit is
counted once per place in the script (same tool and arguments, the check's scratch folder written as WORK).  Whether
a variant without a difference is a defect depends on what the check asserts: read the check for each one listed.
"""
import collections
import os
import re
import sys


def norm(args, check):
    a = re.sub(r"\S*/tmp/" + re.escape(check) + r"/[^ /]*", "WORK", args)
    a = re.sub(r"\S*/shimtmp/\S*", "TMP", a)
    return a


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    out = sys.argv[1]
    details = "--details" in sys.argv
    status = {}
    with open(os.path.join(out, "status.tsv")) as f:
        for line in f:
            name, st, sec = line.rstrip("\n").split("\t")
            status[name] = (st, sec)
    edits = collections.defaultdict(dict)   # check -> site -> [hit count, miss count]
    pymiss = collections.defaultdict(dict)
    runs = collections.defaultdict(list)
    with open(os.path.join(out, "edits.tsv"), errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            check, tool, event, args = parts[0], parts[1], parts[2], "\t".join(parts[3:])
            ev = event.split()
            if tool == "AZURE2":
                _, inkey, outkey, st, proj = ev[:5]
                runs[check].append((inkey, outkey, st, proj, args))
                continue
            site = tool + " " + norm(args, check)
            if ev[0] == "INPLACE":
                hit = ev[1] == "changed"
            elif ev[0] == "FILTER" and ev[2] == "project":
                hit = ev[1] == "differs"
            elif ev[0] == "PYTHON":
                hit = ev[1] != "0"
                if not hit:
                    pymiss[check].setdefault(site, 0)
                    pymiss[check][site] += 1
                    continue
            else:
                continue
            c = edits[check].setdefault(site, [0, 0])
            c[0 if hit else 1] += 1

    rows = []
    totals = collections.Counter()
    for check in sorted(status):
        sites = edits.get(check, {})
        n_sites = len(sites)
        miss_always = [s for s, (h, m) in sites.items() if h == 0]
        miss_some = [s for s, (h, m) in sites.items() if h > 0 and m > 0]
        rr = runs.get(check, [])
        seen_in = {}
        repeats = 0
        by_out = collections.defaultdict(set)
        for inkey, outkey, st, proj, args in rr:
            if inkey in seen_in:
                repeats += 1
            seen_in.setdefault(inkey, args)
            if outkey not in ("-", "none"):
                by_out[outkey].add(inkey)
        nodiff = [ks for ks in by_out.values() if len(ks) > 1]
        rows.append((check, status[check][0], n_sites, len(miss_always), len(miss_some), len(rr), repeats,
                     sum(len(k) for k in nodiff), miss_always, miss_some, nodiff, pymiss.get(check, {}), seen_in))
        totals.update(checks=1, sites=n_sites, miss=len(miss_always), some=len(miss_some), runs=len(rr))

    print("| check | exit | edit places | hit nothing | hit only sometimes | runs | same input again | runs in groups "
          "with equal outputs |")
    print("|---|---|---|---|---|---|---|---|")
    for r in rows:
        print("| %s | %s | %d | %d | %d | %d | %d | %d |" % r[:8])
    print()
    print("checks %(checks)d, edit places %(sites)d, hit nothing %(miss)d, hit only sometimes %(some)d, AZURE2 runs "
          "%(runs)d" % totals)
    if details:
        for r in rows:
            check, miss_always, miss_some, nodiff, pym, seen_in = r[0], r[8], r[9], r[10], r[11], r[12]
            if not (miss_always or miss_some or nodiff or pym):
                continue
            print("\n## " + check)
            for s in miss_always:
                print("  HIT NOTHING      " + s[:400])
            for s in miss_some:
                h, m = edits[check][s]
                print("  HIT %d, MISSED %d  %s" % (h, m, s[:400]))
            for s, n in pym.items():
                print("  PYTHON NO CHANGE x%d  %s" % (n, s[:300]))
            for ks in nodiff:
                print("  EQUAL OUTPUTS from %d inputs:" % len(ks))
                for k in sorted(ks):
                    print("      %s  %s" % (k, seen_in[k][:300]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
