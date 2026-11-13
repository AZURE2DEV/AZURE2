#!/usr/bin/env python3
"""Parameter kinds: what a fit parameter is comes from the record its creation site writes, not from its name.

The engine's table of kinds (include/ParameterKind.h, ``_azure2.parameter_kinds()``) has one row per code 0..11
with its key, unit and flags.  Each parameter's record (``Session.parameter_records()``) carries its kind and what
it belongs to (level and channel, segment key, target effect, sample, THM combination).  Checked on projects that
together hold every kind (tests/upgrade/T2: levels, abundance, norms, the three shift kinds, density;
tests/screening_6li/fit.azr: screening; the grad_fit models of tests/doppler_fgm and tests/transmission_yield:
Doppler temperature and sample thickness; tests/18O_p_a_thm with a cbackground=: THM coherent background; all in
project format 3, whose parameter names are made from the records' names, docs/upgrade/project_format3.md):

  1. every parameter's name, rebuilt from its record and the project's records (the data set of its segment key,
     the normalization, energy scale or target that data set names; the sample of an abundance; the THM
     combination) by the naming rules of format 3, is the engine's name -- the record says the same thing the
     name does, for every kind;
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
sys.path.insert(0, HERE)
import thm_records  # noqa: E402  (a <thm> experiment and its members, on the records)
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
    from pyazr import Project
    from pyazr.parameters import parameter_kinds
except Exception as err:  # pragma: no cover
    print(f"skip: the engine module is not available ({type(err).__name__}: {err})")
    sys.exit(77)


def name_from_record(r, project):
    """The name format 3 gives a parameter, from its record and the project's records (experiment::ParameterNamesOf):
    a level's as in CNuc::FillMnParams; a data set's own nuisance <data set>.<field>; a nuisance of a record shared by
    several data sets <record>.<field> for the first of them, <record>.<field>[<data set>] for the others; an
    abundance <sample>.isotope.<isotope>; a target effect <target>.<field>; a THM coherent background as
    ThmCoherentParamStem."""
    k = r["kind"]
    if k == "energy":
        return f"energy_{r['level_number']}"
    if k == "width":
        return f"width_{r['level_number']}_{r['channel']}"
    if k == "radiative":  # <record>.<part>[.<final state>] (GammaDecaySet::AddParameter)
        return f"{r['decay_record']}.{r['decay_part']}" + (f".{r['final_state']}" if r["final_state"] else "")
    if k == "apparatus":  # the owner's values joined by '.' (ApparatusDeclaration::Name)
        return ".".join(v for _, v in r["owner"])
    if k == "cbkg":
        part = ("re0", "im0", "re1", "im1")[r["thm_term"]]
        return f"cbkg_{r['thm_experiment']}_{r['thm_jpi']}_{r['thm_exit_key']}_{r['thm_channels']}_{part}"
    if k == "abundance":
        isotopes = project[r["sample_name"]].get("isotope", [])
        isotopes = isotopes if isinstance(isotopes, list) else [isotopes]
        return [f"{r['sample_name']}.isotope.{i.split()[0]}" for i in isotopes]  # one of the sample's isotopes
    if k in ("screening", "density", "doppler"):
        field = "areal_density" if k == "density" else k
        return [f"{n}.{field}" for n, t in project.records().items() if t.record_kind == "target" and field in t]
    sets = [n for n, x in project.records().items() if x.record_kind == "dataset"]
    ds = sets[r["segment_key"] - 1]

    def shared(record, field, named_by):  # the record's nuisance FIELD; NAMED_BY: the data sets' field naming it
        users = [n for n in sets if record in _values(project[n], named_by)]
        return f"{record}.{field}" + ("" if users[0] == ds else f"[{ds}]")

    if k == "norm":
        n = project[ds]["norm"]
        return shared(n, "value", "norm") if n in project else f"{ds}.norm"
    if k == "thickness":
        return shared(project[ds]["target"], "areal_density", "target")
    term = {"shift": "offset", "shift_sqrt": "sqrt"}.get(k)
    scales = [v for v in _values(project[ds], "shift") if v in project]
    if k == "shift_fn":
        params = _values(project[scales[0]], "param")
        return f"{scales[0]}.{params[r['fn_term']].split()[0]}"
    holder = [g for g in scales if term in project[g]]
    return shared(holder[0], term, "shift") if holder else f"{ds}.shift.{term}"


def _values(record, field):
    v = record.get(field, [])
    return v if isinstance(v, list) else [v]


def check_project(label, path, cwd, want_kinds):
    project = Project.open(path)
    with azure2(path, cwd=cwd) as s:
        sess = s.sess
        records = sess.parameter_records()
        n = len(records)
        names = [sess.params_names(i) for i in range(n)]
        info = np.asarray(sess.parameter_info(), float).reshape(n, -1)
        rebuilt = [name_from_record(r, project) for r in records]
        bad = [(names[i], want) for i, want in enumerate(rebuilt)
               if (names[i] not in want if isinstance(want, list) else names[i] != want)]
        effects = {}  # a target-effect number is one target: the same name for the same effect
        for i, r in enumerate(records):
            if "effect" in r and isinstance(rebuilt[i], list):
                effects.setdefault((r["kind"], r["effect"]), set()).add(names[i])
        bad += [e for e, n in effects.items() if len(n) != 1]
        check(f"{label}: all {n} names rebuilt from their records", not bad, bad[:5])
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
    shutil.copytree(src_dir, dst, ignore=shutil.ignore_patterns("output", "checks", "run.log"))
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
check("codes 0..13 in order", [k["code"] for k in kinds] == list(range(14)), [k["code"] for k in kinds])
keys = [k["key"] for k in kinds]
check("keys", keys == ["energy", "width", "norm", "shift", "shift_sqrt", "shift_fn", "abundance", "screening",
                       "density", "doppler", "thickness", "cbkg", "apparatus", "radiative"], keys)
check("every row has a unit, owner and doc", all(k["unit"] and k["owner"] and k["doc"] for k in kinds))
check("own_prior: the nuisance kinds 2..10, 12 and 13",
      [k["key"] for k in kinds if k["own_prior"]] == keys[2:11] + ["apparatus", "radiative"])
check("segment_field: norm, the shifts, thickness",
      [k["key"] for k in kinds if k["segment_field"]] == ["norm", "shift", "shift_sqrt", "shift_fn", "thickness"])
check("vary_from_project: the shifts, screening, density, Doppler, thickness, apparatus parameters",
      [k["key"] for k in kinds if k["vary_from_project"]] ==
      ["shift", "shift_sqrt", "shift_fn", "screening", "density", "doppler", "thickness", "apparatus", "radiative"])

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
        # the suite's generator (format 3) on the tree's n12C
        mm = load_module(suite + "_models", os.path.join(ROOT, "tests", suite, "make_models.py"))
        getattr(mm, maker)(os.path.join(ROOT, "tests", "reaction_rate", "models", "n12C.azr"), gen)
        work = os.path.join(gen, "run")
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        shutil.copy(os.path.join(gen, "grad_fit.azr"), work)
        os.symlink(os.path.join(gen, "data"), os.path.join(work, "data"))
        seen |= set(check_project(label, os.path.join(work, "grad_fit.azr"), work, {kind}))

    print("tests/18O_p_a_thm with a linear cbackground=")
    d = copy_project(os.path.join(ROOT, "tests", "18O_p_a_thm"), os.path.join(tmp, "thm"))
    # experiment E of both data sets (they name it in their thm field) with the background
    text = thm_records.thm_set(open(os.path.join(d, "18O_p_a_thm.azr")).read(),
                               "experiment[E] cbackground=1/2+:2:0.5,0,0.5,1:linear=0.25,-0.5f,0,1e-3",
                               thm_records.datasets(open(os.path.join(d, "18O_p_a_thm.azr")).read()))
    open(os.path.join(d, "cb.azr"), "w").write(text)
    seen |= set(check_project("18O_p_a_thm", os.path.join(d, "cb.azr"), d, {"cbkg"}))

    print("tests/stopping_yield/models/fit.azr: two targets of TiN, one stopping scale")
    d = os.path.join(tmp, "stop")
    shutil.copytree(os.path.join(ROOT, "tests", "stopping_yield", "models"), d)
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub))
    kinds_sy = check_project("stopping_yield", os.path.join(d, "fit.azr"), d, {"apparatus", "density"})
    check("stopping_yield: one stopping scale (kind apparatus) for the two targets", kinds_sy.get("apparatus") == 1, kinds_sy)
    seen |= set(kinds_sy)

    print("tests/rectangular_resolution/models/fit.azr: one fitted width of a detector named by two data sets")
    d = os.path.join(tmp, "rect")
    shutil.copytree(os.path.join(ROOT, "tests", "rectangular_resolution", "models"), d)
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub))
    kinds_rr = check_project("rectangular_resolution", os.path.join(d, "fit.azr"), d, {"apparatus"})
    check("rectangular_resolution: one apparatus parameter for the two data sets", kinds_rr.get("apparatus") == 1, kinds_rr)

    print("tests/gamma_decay/models/record.azr: a level written as Gamma_gamma, branchings and mixing ratios")
    d = copy_project(os.path.join(ROOT, "tests", "gamma_decay", "models"), os.path.join(tmp, "gamma"))
    kinds_gd = check_project("gamma_decay", os.path.join(d, "record.azr"), d, {"radiative", "width"})
    check("gamma_decay: one radiative parameter per photon channel of level 5 (8)", kinds_gd.get("radiative") == 8, kinds_gd)
    seen |= set(kinds_gd)

check("the projects hold every kind", seen == set(keys), sorted(set(keys) - seen))

if failures:
    print(f"{len(failures)} failure(s)")
    sys.exit(1)
print("parameter kinds: ok")
