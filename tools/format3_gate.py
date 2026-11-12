#!/usr/bin/env python3
"""The converted-tree gate of project format 3 (docs/upgrade/project_format3.md section 8.5): its two Python steps.

  format3_gate.py tree CONVERTED TREE
      A runnable tree from the output of tools/convert_tree.sh: each original whose conversion validated (or whose
      only validation errors are E_DATA_FILE, a data file a test's check.sh writes) is replaced by its .f3.azr under
      the original's name, so its outputs land where the original's did; a project the converter left alone because
      it is format 3 already stays as it is (with the data files it names, which the copy of the tree holds), so its
      outputs are compared too and must be identical; the others are removed and listed with the converter's reason.

  format3_gate.py compare CMP REF NEW CONVERTED
      Every output file of REF (the original tree's run, tools/runall.sh) against NEW (the converted tree's), per
      basis (brune, no-brune, use-park, under CMP).  A file is identical (param.par after the converter's
      parameter-name map, CONVERTED/tests/<dir>/<name>.report.json), or differs in a way named in CLASSES below, or
      differs.  run.log is not compared.  Prints per basis "## <basis>: N files of REF: a identical, b classified,
      c differ, d missing", then every classified and differing file.  Exit status 1 when a file differs or is
      missing.

A difference is accepted only as a row of CLASSES: the rules below say which classes a difference consists of, and
every one of them must be listed for that file of that project.  A new difference therefore fails the gate until a
row with its reason is added here.
"""
import json
import os
import subprocess
import sys

BASES = ("brune", "no-brune", "use-park")

# ---------------------------------------------------------------- the classes (explicit data)
#
# Each row: the class, the output file it applies to ("*": every file of the project), the projects (<test dir>/
# <file>.azr, as tools/runall.sh names them) and why the converted file's output differs from the original's.

CLASSES = [
    {
        "class": "step of a parameter",
        "file": "param.par",
        "projects": ["13N/13N.azr", "hybrid_potential/hybrid_potential.azr",
                     "shift_sqrt_subpoints/models/narrow_sqrt.azr", "shift_sqrt_subpoints/models/noshift.azr",
                     "shift_sqrt_subpoints/models/shiftfn.azr", "shift_sqrt_subpoints/models/sqrt.azr",
                     "thick_target_resonance/thick_target_resonance.azr"],
        "why": "the start step Minuit is given for a fixed zero energy shift differs between the formats; the values "
               "are equal (director's decision 4 on the G3 report)",
    },
    {
        "class": "follower",
        "file": "param.par",
        "projects": ["angle_mapping/angle_mapping.azr", "upgrade/T2/T2.azr", "upgrade/T2/T2_true.azr"],
        "why": "a later member of a shared norm or shift group is written as <record>.<field>[<data set>], a fixed "
               "follower row of the group's parameter (SetupSpecs.h, ParameterNames)",
    },
    {
        "class": "follower listed",
        "file": "normalizations.out",
        "projects": ["upgrade/T2/T2.azr"],
        "why": "the same followers listed as rows of normalizations.out",
    },
    {
        "class": "angle window of an angle-integrated data set",
        "file": "normalizations.out",
        "projects": ["17O/17O.azr"],
        "why": "an angle-integrated data set has no angle window in format 3: the angle columns of its row differ",
    },
    {
        "class": "angle window of an angle-integrated data set",
        "file": "shifts.out",
        "projects": ["secondary_gamma/models/mixing.azr", "secondary_gamma/models/spin2_minus.azr",
                     "secondary_gamma/models/spin2_multil.azr", "secondary_gamma/models/spin2_plus.azr",
                     "secondary_gamma/models/spin2_zero.azr"],
        "why": "as for normalizations.out: the angle columns of an angle-integrated data set's row",
    },
    {
        "class": "sample named by its target",
        "file": "AZUREOut_sample_seg=*.out",
        "projects": ["angle_mapping/angle_mapping.azr", "natural_neon/natural_neon.azr"],
        "why": "the '# Segment' header names the sample by its target record instead of its old number",
    },
    {
        "class": "scattering lengths as model data (label)",
        "file": "chiSquared.out",
        "projects": ["scattering_length/be7p_fit.azr", "scattering_length/he3he4_fit.azr", "upgrade/T3/T3.azr",
                     "upgrade/T3/T3_start.azr"],
        "why": "the old <scatteringLengths> section is converted to <modelData> rows: its chi-squared line is "
               "labelled Model-Data instead of Scattering-Length, with the same values",
    },
    {
        "class": "scattering lengths as model data (table)",
        "file": "scatteringLengths.out",
        "projects": ["scattering_length/be7p_fit.azr", "scattering_length/he3he4_fit.azr", "upgrade/T3/T3.azr",
                     "upgrade/T3/T3_start.azr"],
        "why": "scatteringLengths.out is written without the old section's table of measured values (they are in "
               "modelData.out)",
    },
    {
        "class": "refused: data files absent (a template)",
        "file": "*",
        "projects": ["shift_subpoint_angle/template.azr"],
        "why": "a template whose data files its check.sh writes: the converted file is refused (E_DATA_FILE) where "
               "the old engine skipped the missing files, so it writes no output (director's decision 3 on G3)",
    },
]


def _match(pattern, name):
    if pattern == "*" or pattern == name:
        return True
    if "*" in pattern:
        head, tail = pattern.split("*", 1)
        return name.startswith(head) and name.endswith(tail) and len(name) >= len(head) + len(tail)
    return False


def listed(cls, filename, project):
    """The row of CLASSES for this class, file and project, or None."""
    for row in CLASSES:
        if row["class"] == cls and _match(row["file"], filename) and project in row["projects"]:
            return row
    return None


# ---------------------------------------------------------------- the rules: which classes a difference consists of


def _lines(b):
    return b.decode(errors="replace").splitlines()


def _param_par(a, b, names):
    xs = [[names.get(t[0], t[0])] + t[1:] if t else t for t in (l.split() for l in a)]
    ys = [l.split() for l in b]
    if xs == ys:
        return set()
    known = {t[0] for t in xs if t}
    ys2 = [t for t in ys if not (t and "[" in t[0] and t[0] not in known)]
    if len(xs) != len(ys2):
        return None
    classes = set()
    for p, q in zip(xs, ys2):
        if p == q:
            continue
        if p and q and p[0] == q[0] and "[" in p[0]:
            classes.add("follower")
        elif p and q and p[0] == q[0] and p[1] == q[1]:
            classes.add("step of a parameter")
        else:
            return None
    if len(ys2) != len(ys):
        classes.add("follower")
    return classes


def _rows_by_label(a, b):
    extra = [l for l in b if l not in a]
    missing = {l.split(",")[0]: l for l in a if l not in b}
    classes = set()
    for l in extra:
        k = l.split(",")[0]
        if k in missing:
            u, v = missing.pop(k).split(","), l.split(",")
            if len(u) == len(v) and all(i in (2, 3) or u[i] == v[i] for i in range(len(u))):
                classes.add("angle window of an angle-integrated data set")
                continue
            return None
        classes.add("follower listed")
    return None if missing else classes


def classes_of(filename, old, new, names):
    """The classes a difference of one file consists of: an empty set when the files are equal (param.par after the
    name map), None when no rule explains it."""
    if old == new:
        return set()
    a, b = _lines(old), _lines(new)
    if filename == "param.par":
        return _param_par(a, b, names)
    if filename in ("normalizations.out", "shifts.out"):
        return _rows_by_label(a, b)
    if len(a) == len(b) and all(p == q or (p.startswith("# Segment") and "(sample " in p
                                           and p.split("(sample")[0] == q.split("(sample")[0]) for p, q in zip(a, b)):
        return {"sample named by its target"}
    if filename == "chiSquared.out" and [l.replace("Scattering-Length", "Model-Data") for l in a] == b:
        return {"scattering lengths as model data (label)"}
    if filename == "scatteringLengths.out":
        return {"scattering lengths as model data (table)"}
    return None


# ---------------------------------------------------------------- tree


ALREADY_FORMAT3 = "already format 3"  # the converter's refusal of a format-3 project (its report's "refusal")


def _read_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def not_converted_reason(base):
    """Why the original BASE.azr was not converted: the report's refusal, else the first validation error of the
    conversion, else the ERROR line of the converter's log."""
    report = _read_json(base + ".report.json")
    if isinstance(report, dict) and report.get("refusal"):
        return report["refusal"]
    diags = _read_json(base + ".validate.json")
    if isinstance(diags, list):
        errors = [x for x in diags if isinstance(x, dict) and x.get("severity") == "error"]
        if errors:
            return f"{errors[0].get('code')}: {errors[0].get('message')} ({len(errors)} error(s))"
    try:
        with open(base + ".convert.log", errors="replace") as f:
            for line in f:
                if line.startswith("ERROR"):
                    return line.strip()
    except OSError:
        pass
    return "no conversion, report or log"


def make_tree(conv, tree):
    """The runnable tree (see the module docstring); prints the files kept without data files, the format-3 projects
    kept as they are and those left out."""
    if os.path.exists(tree):
        subprocess.run(["rm", "-rf", tree], check=True)
    os.makedirs(tree)
    subprocess.run(["cp", "-r", os.path.join(conv, "tests"), os.path.join(tree, "tests")], check=True)
    root = os.path.join(tree, "tests")
    originals = []
    for d, _, files in os.walk(root):
        if "output" in os.path.relpath(d, root).split(os.sep):  # outputs of earlier runs are not projects
            continue
        originals += [os.path.join(d, f) for f in files if f.endswith(".azr") and not f.endswith(".f3.azr")]
    for f in sorted(originals):
        rel = "tests/" + os.path.relpath(f, root)
        base = f[:-4]
        converted, validation = base + ".f3.azr", base + ".validate.json"
        kept, without_data = False, False
        if os.path.exists(converted) and os.path.exists(validation):
            try:
                diags = json.load(open(validation))
                kept = all(x["severity"] != "error" or x["code"] == "E_DATA_FILE" for x in diags)
                without_data = any(x["code"] == "E_DATA_FILE" for x in diags)
            except (ValueError, KeyError, TypeError):
                kept = False
        report = _read_json(base + ".report.json")
        if kept:
            if without_data:
                print("kept without its data files: " + rel)
            os.replace(converted, f)
        elif isinstance(report, dict) and report.get("refusal") == ALREADY_FORMAT3 and not os.path.exists(converted):
            print("format 3 already, kept as it is: " + rel)
        else:
            print(f"not in the tree: {rel} ({not_converted_reason(base)})")
            for g in (f, converted):
                if os.path.exists(g):
                    os.remove(g)
    for d, _, files in os.walk(root):
        for f in files:
            if f.endswith(".f3.azr"):
                os.remove(os.path.join(d, f))


# ---------------------------------------------------------------- compare


def compare(cmp_dir, ref, new, conv):
    failed = False
    for base in BASES:
        rdir, ndir = os.path.join(cmp_dir, base, ref), os.path.join(cmp_dir, base, new)
        same = differ = missing = 0
        diffs, classified = [], []
        for res in sorted(os.listdir(rdir)):
            if not res.startswith(".res-"):
                continue
            reldir = res[5:].replace("__", "/")
            for azr in sorted(os.listdir(os.path.join(rdir, res))):
                project = f"{reldir}/{azr}"
                report = os.path.join(conv, "tests", reldir, azr[:-4] + ".report.json")
                names = json.load(open(report))["parameters"] if os.path.exists(report) else {}
                for sub in ("output", "checks"):
                    rd = os.path.join(rdir, res, azr, sub)
                    if not os.path.isdir(rd):
                        continue
                    for root, _, files in os.walk(rd):
                        for f in sorted(files):
                            a = os.path.join(root, f)
                            rel = os.path.relpath(a, os.path.join(rdir, res, azr))
                            b = os.path.join(ndir, res, azr, rel)
                            key = f"{project}/{rel}"
                            if not os.path.exists(b):
                                log = os.path.join(ndir, res, azr, "run.log")
                                refused = os.path.exists(log) and "error E_DATA_FILE" in open(log, errors="replace").read()
                                cls = "refused: data files absent (a template)"
                                if refused and listed(cls, f, project):
                                    classified.append(f"{cls}: {key}")
                                else:
                                    missing += 1
                                    diffs.append("missing " + key)
                                continue
                            with open(a, "rb") as fa, open(b, "rb") as fb:
                                found = classes_of(f, fa.read(), fb.read(), names)
                            if found is not None and not found:
                                same += 1
                            elif found is not None and all(listed(c, f, project) for c in found):
                                classified.append(f"{f}: {', '.join(sorted(found))}: {key}")
                            else:
                                differ += 1
                                why = "" if found is None else f" (class not listed for this file and project: " \
                                                                f"{', '.join(sorted(found))})"
                                diffs.append("differs " + key + why)
        print(f"## {base}: {same + len(classified) + differ + missing} files of {ref}: {same} identical, "
              f"{len(classified)} classified, {differ} differ, {missing} missing")
        for d in classified:
            print(f"   {base}: {d}")
        for d in diffs:
            print(f"   {base}: {d}")
        failed = failed or differ > 0 or missing > 0
    return 1 if failed else 0


def main(argv):
    if len(argv) == 4 and argv[1] == "tree":
        make_tree(argv[2], argv[3])
        return 0
    if len(argv) == 6 and argv[1] == "compare":
        return compare(*argv[2:6])
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
