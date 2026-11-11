#!/usr/bin/env python3
"""Parameter kinds: what a fit parameter is comes from the record its creation site writes, not from its name.

The engine's table of kinds (include/ParameterKind.h, ``_azure2.parameter_kinds()``) has one row per code 0..11
with its key, unit and flags.  Each parameter's record (``Session.parameter_records()``) carries its kind and what
it belongs to (level and channel, segment key, target effect, sample, THM combination).  Checked on projects that
together hold every kind (tests/upgrade/T2: levels, abundance, norms, the three shift kinds, density;
tests/screening_6li/fit.azr: screening; the grad_fit models of tests/doppler_fgm and tests/transmission_yield:
Doppler temperature and sample thickness; tests/18O_p_a_thm with a cbackground=: THM coherent background):

  1. every parameter's name, rebuilt from its record alone by the naming rules of the creation sites, is the
     engine's name -- the record says the same thing the name does, for every kind;
  2. the record's code is the type code of parameter_info, its owner fields agree with parameter_info's numeric
     ones (segment key, effect), and its ordinal counts the parameters of its kind;
  3. pyazr's Parameter.kind is the record's key, an abundance's sample its record's sample name, and a cbkg
     parameter's thm_* fields its combination.

Needs the compiled engine; without it exits 77 (ctest: Skipped).

Run from anywhere:  python3 tests/pyazr/parameter_kinds_test.py
"""
import importlib.util
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import old_tree  # noqa: E402  (the old originals of converted projects)
os.environ.setdefault("OMP_NUM_THREADS", "2")

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    from pyazr import azure2
    from pyazr import _azure2
    from pyazr.azrfile import AzrModel
    from pyazr.parameters import parameter_kinds
except Exception as err:  # pragma: no cover
    print(f"skip: the engine module is not available ({type(err).__name__}: {err})")
    sys.exit(77)


def name_from_record(r):
    """The name the creation sites give a parameter, from its record alone (CNuc::FillMnParams,
    CNuc::AbundanceParamName, WorkingCopy::SegmentParameterName, WorkingCopy::FillMnParams, ThmCoherentParamStem)."""
    k = r["kind"]
    if k == "energy":
        return f"energy_{r['level_number']}"
    if k == "width":
        return f"width_{r['level_number']}_{r['channel']}"
    if k == "abundance":
        sample = f"{r['sample_name']}_" if r["sample_name"] else ""
        return f"abundance_{sample}{r['nucleus_Z']}_{r['nucleus_A']}"
    field = {"norm": "norm", "shift": "energy_shift", "shift_sqrt": "energy_shift_sqrt", "thickness": "thickness"}
    if k in field:
        return f"segment_{r['segment_key']}_{field[k]}"
    if k == "shift_fn":
        return f"segment_{r['segment_key']}_energy_shift_fn{r['fn_term']}"
    if k in ("screening", "density", "doppler"):
        return f"{k}_{r['effect']}"
    if k == "stopping_scale":
        return f"stopping_scale.{r['material']}.{r['projectile']}"
    if k == "cbkg":
        part = ("re0", "im0", "re1", "im1")[r["thm_term"]]
        return f"cbkg_{r['thm_experiment']}_{r['thm_jpi']}_{r['thm_exit_key']}_{r['thm_channels']}_{part}"
    return None


def check_project(label, path, cwd, want_kinds, old_names=True):
    with azure2(path, cwd=cwd) as s:
        sess = s.sess
        records = sess.parameter_records()
        n = len(records)
        names = [sess.params_names(i) for i in range(n)]
        info = np.asarray(sess.parameter_info(), float).reshape(n, -1)
        # a format-3 project names its data-set parameters by its records (scan1.norm): only the others are rebuilt
        rebuilt = [i for i, r in enumerate(records) if old_names or r["kind"] not in ("norm", "shift", "shift_sqrt", "density")]
        bad = [(names[i], name_from_record(records[i])) for i in rebuilt if name_from_record(records[i]) != names[i]]
        check(f"{label}: all {len(rebuilt)} names rebuilt from their records", not bad, bad[:5])
        check(f"{label}: record code == parameter_info type code",
              all(r["code"] == int(round(info[i][0])) for i, r in enumerate(records)))
        seg = [i for i, r in enumerate(records) if "segment_key" in r]
        check(f"{label}: segment key == parameter_info's ({len(seg)} segment parameters)",
              all(int(round(info[i][13])) == records[i]["segment_key"] for i in seg))
        eff = [i for i, r in enumerate(records) if "effect" in r]
        check(f"{label}: effect == parameter_info's ({len(eff)} target-effect parameters)",
              all(int(round(info[i][6])) == records[i]["effect"] for i in eff))
        counts, ordinals_ok = {}, True
        for r in records:
            counts[r["kind"]] = counts.get(r["kind"], 0) + 1
            ordinals_ok &= r["ordinal"] == counts[r["kind"]]
        check(f"{label}: ordinals count each kind", ordinals_ok)
        check(f"{label}: holds {sorted(want_kinds)}", want_kinds <= set(counts), sorted(counts))
        params = s.parameters
        check(f"{label}: Parameter.kind is the record's key",
              [p.kind for p in params] == [r["kind"] for r in records])
        ab = [(p, r) for p, r in zip(params, records) if r["kind"] == "abundance"]
        check(f"{label}: Parameter.sample is the record's sample name",
              all(p.sample == r["sample_name"] for p, r in ab))
        cb = [(p, r) for p, r in zip(params, records) if r["kind"] == "cbkg"]
        check(f"{label}: cbkg Parameter.thm_* are the record's",
              all((p.thm_experiment, p.thm_jpi, p.thm_exit_key, p.thm_channels, p.thm_term) ==
                  (r["thm_experiment"], r["thm_jpi"], r["thm_exit_key"], r["thm_channels"], r["thm_term"])
                  for p, r in cb))
        return counts


def copy_project(src_dir, dst):
    # the old originals of converted projects: the names are rebuilt here by the old scheme (segment_<k>_norm, ...)
    old_tree.copytree(src_dir, dst, ignore=shutil.ignore_patterns("output", "checks", "run.log"))
    os.makedirs(os.path.join(dst, "output"), exist_ok=True)
    os.makedirs(os.path.join(dst, "checks"), exist_ok=True)
    return dst


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# -- the table ----------------------------------------------------------------------------------------------------
print("the table of kinds")
kinds = parameter_kinds()
check("codes 0..12 in order", [k["code"] for k in kinds] == list(range(13)), [k["code"] for k in kinds])
keys = [k["key"] for k in kinds]
check("keys", keys == ["energy", "width", "norm", "shift", "shift_sqrt", "shift_fn", "abundance", "screening",
                       "density", "doppler", "thickness", "cbkg", "stopping_scale"], keys)
check("every row has a unit, owner and doc", all(k["unit"] and k["owner"] and k["doc"] for k in kinds))
check("own_prior: the nuisance kinds 2..10 and 12",
      [k["key"] for k in kinds if k["own_prior"]] == keys[2:11] + ["stopping_scale"])
check("segment_field: norm, the shifts, thickness",
      [k["key"] for k in kinds if k["segment_field"]] == ["norm", "shift", "shift_sqrt", "shift_fn", "thickness"])
check("vary_from_project: the shifts, screening, density, Doppler, thickness, stopping scale",
      [k["key"] for k in kinds if k["vary_from_project"]] ==
      ["shift", "shift_sqrt", "shift_fn", "screening", "density", "doppler", "thickness", "stopping_scale"])

seen = set()
with tempfile.TemporaryDirectory() as tmp:
    print("tests/upgrade/T2")
    d = copy_project(os.path.join(ROOT, "tests", "upgrade", "T2"), os.path.join(tmp, "T2"))
    seen |= set(check_project("T2", os.path.join(d, "T2.azr"), d,
                              {"energy", "width", "abundance", "shift", "shift_sqrt", "shift_fn", "density"}))

    print("tests/screening_6li/fit.azr")
    d = copy_project(os.path.join(ROOT, "tests", "screening_6li"), os.path.join(tmp, "scr"))
    seen |= set(check_project("screening_6li", os.path.join(d, "fit.azr"), d, {"screening", "norm"}))

    for label, suite, maker, kind in (("doppler_fgm", "doppler_fgm", "main", "doppler"),
                                      ("transmission_yield", "transmission_yield", "n12c_models", "thickness")):
        print(f"tests/{suite} grad_fit")
        gen = os.path.join(tmp, suite)
        os.makedirs(gen)
        # the old generator and the old n12C (tests/convert/fixtures): these models are the old file class's
        mm = load_module(suite + "_models", old_tree.original(os.path.join(ROOT, "tests", suite, "make_models.py")))
        getattr(mm, maker)(old_tree.original(os.path.join(ROOT, "tests", "reaction_rate", "models", "n12C.azr")), gen)
        work = os.path.join(gen, "run")
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        shutil.copy(os.path.join(gen, "grad_fit.azr"), work)
        os.symlink(os.path.join(gen, "data"), os.path.join(work, "data"))
        seen |= set(check_project(label, os.path.join(work, "grad_fit.azr"), work, {kind}))

    print("tests/18O_p_a_thm with a linear cbackground=")
    d = copy_project(os.path.join(ROOT, "tests", "18O_p_a_thm"), os.path.join(tmp, "thm"))
    m = AzrModel.from_file(os.path.join(d, "18O_p_a_thm.azr"))
    m.set_thm_experiment("E", [1, 2], cbackground="1/2+:2:0.5,0,0.5,1:linear=0.25,-0.5f,0,1e-3")
    m.write(os.path.join(d, "cb.azr"))
    seen |= set(check_project("18O_p_a_thm", os.path.join(d, "cb.azr"), d, {"cbkg"}))

    print("tests/stopping_yield/models/fit.azr: two targets of TiN, one stopping scale")
    d = os.path.join(tmp, "stop")
    shutil.copytree(os.path.join(ROOT, "tests", "stopping_yield", "models"), d)
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub))
    kinds_sy = check_project("stopping_yield", os.path.join(d, "fit.azr"), d, {"stopping_scale", "density"}, old_names=False)
    check("stopping_yield: one stopping_scale for the two targets", kinds_sy.get("stopping_scale") == 1, kinds_sy)
    seen |= set(kinds_sy)

check("the projects hold every kind", seen == set(keys), sorted(set(keys) - seen))

if failures:
    print(f"{len(failures)} failure(s)")
    sys.exit(1)
print("parameter kinds: ok")
