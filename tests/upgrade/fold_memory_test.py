#!/usr/bin/env python3
"""The memory of a project with many folded data sets (docs/upgrade/data_model.md section 9).

The project is made here from examples/c12c12_tumino2018 (in the repository): its twelve extrapolations become eight
active excitation functions of 2 to 3 MeV (lab) -- the four THM lines (entrance pair 1, 0.02 MeV steps), folded as in
the example with its beam_E1 spread and a 60-point quadrature, and the four direct lines (entrance pair 6, 0.004 MeV
steps), folded here with the example's own beam_targetInt7 (p0, a1) and beam_targetInt5 (p1, a0) spreads and a
10-point quadrature -- all on the one shared mesh of 24Mg.  The same project with no extrapolation folded is the
baseline.  `AZURE2 extrapolate` runs both at one thread; the peak resident memory of the folded run minus that of the
baseline (getrusage of the child, the same program and Coulomb table in both) is what the folds hold.

The bound, 125 MiB: measured 97.6 MiB on Linux (repeated runs agree to 0.1 MiB; 97 MiB without a Coulomb table),
with about 25 % left for another allocator or C++ library; 100.4 MiB on macOS, where the runs have libmalloc's
cache of freed large blocks switched off (MallocLargeCache=0), because the peak resident memory counts the blocks
that cache keeps: 31 MiB of them here, the setup's vectors grown by push_back and freed again.
What it catches (each measured on this project): the THM
vertex pieces stored for every channel of the compound instead of the entrance pair's, 139 MiB; every slot computing
and holding its own channel rows instead of one table per compound energy (WorkingCopy::FillChannelRows), 229 MiB;
the build before both, 334 MiB.

    fold_memory_test.py --binary AZURE2

Skipped (exit 77) where the resource module is missing (Windows).  About 1 s a run with a Coulomb table, 22 s without.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

try:
    import resource
except ImportError:  # Windows
    print("skip: no resource module (getrusage) on this platform")
    sys.exit(77)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
EXAMPLE = os.path.join(ROOT, "examples", "c12c12_tumino2018")
BOUND_MIB = 125.0
# The direct extrapolations folded here, with the beam of the example's <targets> each takes.
DIRECT_BEAMS = {9: "beam_targetInt7", 10: "beam_targetInt5", 11: "beam_targetInt5", 12: "beam_targetInt7"}

failures = 0


def check(name, ok, detail=""):
    global failures
    print(("  ok    " if ok else "  FAIL  ") + name + (": " + detail if detail else ""))
    if not ok:
        failures += 1


def project(folded):
    """The example's text with every extrapolation on 2-3 MeV lab; the direct ones folded when `folded`, and no fold
    at all when not."""
    out, line_no, inside = [], 0, False
    with open(os.path.join(EXAMPLE, "c12c12_tumino2018.azr")) as f:
        for line in f.read().split("\n"):
            if line.startswith("<extrapolations>"):
                inside = True
            elif line.startswith("</extrapolations>"):
                inside = False
            m = re.match(r"extrapolation x(\d+)$", line)
            if m:
                line_no = int(m.group(1))
            if inside and line.startswith("    energies"):
                step = "0.02" if line_no <= 4 else "0.004"
                out.append("    energies       range 2 3 MeV lab step " + step)
                if folded and line_no in DIRECT_BEAMS:
                    out.append("    beam           " + DIRECT_BEAMS[line_no])
                    out.append("    quadrature     points 10 width_multiplier 20 per_width 50")
                continue
            if inside and not folded and (line.startswith("    beam ") or line.startswith("    quadrature ")):
                continue
            out.append(line)
    return "\n".join(out)


def peak_kb(binary, directory, name):
    """Run `extrapolate` on name.azr in `directory`; the largest peak resident memory of the children so far (kB)
    and the exit status."""
    env = dict(os.environ, OMP_NUM_THREADS="1")
    if sys.platform == "darwin":
        # libmalloc keeps freed blocks above its small sizes resident for reuse (the large cache), and the peak
        # resident memory counts them: 131.3 MiB there against 100.4 MiB with the cache off (macOS 15.7, Intel).
        env["MallocLargeCache"] = "0"
    with open(os.path.join(directory, name + ".log"), "w") as log:
        rc = subprocess.run([binary, "extrapolate", "--no-gui", name + ".azr"], cwd=directory, stdout=log,
                            stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, env=env, timeout=900).returncode
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if sys.platform == "darwin":  # bytes there, kB on Linux
        peak //= 1024
    return peak, rc


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--binary", default=os.environ.get("AZURE2_BIN"))
    parser.add_argument("--keep", action="store_true", help="keep the scratch directory")
    args = parser.parse_args()
    if not args.binary:
        sys.exit("fold_memory_test.py: give --binary or AZURE2_BIN")
    binary = os.path.abspath(args.binary)
    work = tempfile.mkdtemp(prefix="fold_memory.")
    try:
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        shutil.copytree(os.path.join(EXAMPLE, "data"), os.path.join(work, "data"))
        for name, folded in (("unfolded", False), ("folded", True)):
            with open(os.path.join(work, name + ".azr"), "w") as f:
                f.write(project(folded))
        # the smaller run first: RUSAGE_CHILDREN keeps the largest peak of the children waited for
        base, rc0 = peak_kb(binary, work, "unfolded")
        check("the project without folds runs", rc0 == 0, "exit %d" % rc0)
        peak, rc1 = peak_kb(binary, work, "folded")
        check("the folded project runs", rc1 == 0, "exit %d" % rc1)
        written = sorted(f for f in os.listdir(os.path.join(work, "output")) if f.endswith(".extrap"))
        check("eight extrapolation files", len(written) == 8, ", ".join(written))
        check("the folded run's peak is its own", peak > base, "%d kB against %d kB" % (peak, base))
        held = (peak - base) / 1024.0
        check("the folds hold at most %.0f MiB" % BOUND_MIB, held <= BOUND_MIB,
              "%.1f MiB (peak %.1f MiB, without folds %.1f MiB)" % (held, peak / 1024.0, base / 1024.0))
    finally:
        if args.keep:
            print("scratch kept: " + work)
        else:
            shutil.rmtree(work, ignore_errors=True)
    print("fold_memory: %d failure(s)" % failures)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
