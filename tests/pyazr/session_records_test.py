#!/usr/bin/env python3
"""pyazr on a format-3 project without the old file class (brief: validation follow-ups):

  1. a calculation session on a format-3 project (T2's conversion): datasets, active_datasets and evaluate().data come
     from the engine's description of the built project: the data sets' record names in file order, one per segment
     the engine evaluates, rows kept as --describe counts them; the chi^2 per data set adds up to the session's chi^2
     and equals that of the tree's T2 (tests/upgrade/T2, converted with exact grids; the same data); extrapolations of a format-3 project likewise;
  2. Project.rename takes the parameter names from the engine: a target's isotope abundance (natN.isotope.15N) and a
     shared normalization's follower name (D.value[s03_diff150]) in <parameterSettings> follow, a word that only
     contains the name does not (and a warning lists its row); an inactive data set's <parameterSettings> and <mcmc>
     rows follow its rename and name the engine's parameters once it is switched on;
  3. Project.set keeps a whole-line comment between two repeats of the field it replaces;
  4. records_from_csv: data sets declared in a sheet, each row checked by the engine against the project, every
     refusal at once with its row; a valid sheet's records are added and the project validates.

Needs the compiled module.  Run from anywhere:  python3 tests/pyazr/session_records_test.py
"""
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
T2 = os.path.join(ROOT, "tests", "upgrade", "T2")
FIXTURES = os.path.join(ROOT, "convert", "tests", "fixtures")
RTH = os.path.join(ROOT, "tests", "ratio_to_rutherford", "models")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
try:
    import pyazr  # noqa: E402
    from pyazr import Session, Project, FormatError, records_from_csv  # noqa: E402
    pyazr.Session  # the engine module is loaded on first use
    from pyazr import _azure2  # noqa: F401,E402
except ImportError as err:
    print(f"skip: pyazr or its engine module cannot be imported ({err})")
    sys.exit(77)


def t2_dir(work, name):
    d = os.path.join(work, name)
    os.makedirs(os.path.join(d, "output"))
    os.symlink(os.path.join(T2, "data"), os.path.join(d, "data"))
    shutil.copy(os.path.join(T2, "cov_diff.txt"), d)
    shutil.copy(os.path.join(FIXTURES, "T2.prior_cov.f3.txt"), os.path.join(d, "prior_cov.f3.txt"))
    shutil.copy(os.path.join(FIXTURES, "T2.f3.azr"), os.path.join(d, "T2.azr"))
    return os.path.join(d, "T2.azr")


def session(work):
    print("1. the session's data sets on format 3")
    f3 = t2_dir(work, "T2f3")
    names = [n for n, r in Project.open(f3).records().items() if r.record_kind == "dataset"]
    with Session(f3) as m:
        ds, active = m.datasets, m.active_datasets
        chi = {d.name: (d.chi2, d.n) for d in m.evaluate().data}
        total = m.evaluate(m.parameters.x0).chi2
        described = pyazr.describe(f3, as_json=True)["datasets"] if os.getcwd() == os.path.dirname(f3) else None
        check("datasets: the engine's description, the record names in file order",
              type(ds).__name__ == "DescribedDatasets" and [d.name for d in ds] == names, [d.name for d in ds])
        check("active_datasets: one per segment the engine evaluates", len(active) == m.nsegments,
              (len(active), m.nsegments))
        check("rows kept: the energies the engine evaluates",
              all(d.rows_kept == len(m.energies[i]) for i, d in enumerate(active)))
        check("evaluate().data: named by record name, adding up to the session's chi^2",
              list(chi) == [d.name for d in active] and abs(sum(v[0] for v in chi.values()) - total) <= 1e-9 * total,
              (sum(v[0] for v in chi.values()), total))
        del described
    tree = shutil.copytree(T2, os.path.join(work, "T2tree"))
    os.makedirs(os.path.join(tree, "output"), exist_ok=True)
    with Session(os.path.join(tree, "T2.azr"), cwd=tree) as m:
        old = {d.name: (d.chi2, d.n) for d in m.evaluate().data}
    check("data[i].chi2 equals the tree's T2's per data set",
          set(old) == set(chi) and all(abs(old[k][0] - chi[k][0]) <= 1e-9 * max(1.0, old[k][0]) and old[k][1] == chi[k][1]
                                       for k in old), [(k, old[k], chi.get(k)) for k in old][:3])
    d = os.path.join(work, "rth")
    shutil.copytree(RTH, d)
    os.makedirs(os.path.join(d, "output"), exist_ok=True)
    with Session(os.path.join(d, "aC_hs.azr"), extrapolations=True) as m:
        x = m.extrapolations
        names = [n for n, r in Project.open(os.path.join(d, "aC_hs.azr")).records().items()
                 if r.record_kind == "extrapolation"]
        check("extrapolations: the engine's description, by record name, one per segment",
              type(x).__name__ == "DescribedExtrapolations" and [e.name for e in x] == names
              and len(x.active) == m.nsegments, ([e.name for e in x], names, m.nsegments))


def rename(work):
    print("2. rename with the engine's parameter names")
    f3 = t2_dir(work, "T2rename")
    with open(f3) as f:
        text = f.read().replace("<parameterSettings>\n", "<parameterSettings>\n"
                                "natN.isotope.15N 0.0037 0 0 0.0002 0 1 abundance 0\n"
                                "D.value[s03_diff150] 1 0 0 0.05 0 1 norm 0\n"
                                "natN_note 1 0 0 0 0 1 note 0\n")
    with open(f3, "w") as f:
        f.write(text)
    p = Project.open(f3)
    p.rename("natN", "natX")
    p.rename("s03_diff150", "diff150")
    body = p.text.split("<parameterSettings>")[1].split("</parameterSettings>")[0]
    check("the isotope's abundance follows the target (natX.isotope.15N)", "\nnatX.isotope.15N 0.0037" in body, body)
    check("the follower's name follows the data set (D.value[diff150])", "\nD.value[diff150] 1" in body, body)
    check("a word that only contains the name is left alone", "\nnatN_note 1" in body, body)
    check("the renamed project validates", p.validate() == [], p.validate())


def rename_inactive(work):
    print("2b. rename of an inactive data set: its rows follow, and apply once it is switched on")
    f3 = t2_dir(work, "T2inactive")
    with open(f3) as f:
        text = f.read().replace("<parameterSettings>\n", "<parameterSettings>\n"
                                "s04_thick.norm 1 0 0 0.03 0 1 norm 0\n"
                                "natN_note 1 0 0 0 0 1 note 0\n")
    text = text.rstrip("\n") + "\n<mcmc>\nnwalkers 10\n<parameters>\ns04_thick.norm 1 0.03\n</parameters>\n</mcmc>\n"
    with open(f3, "w") as f:
        f.write(text)
    p = Project.open(f3)
    p["s04_thick"].set(active="no")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        p.rename("s04_thick", "thick4")
        p.rename("natN", "natY")  # a row that contains the name but names no parameter: a warning lists it
    p["thick4"].set(active="yes")
    rows = [l.split()[0] for l in p.text.splitlines() if l.startswith(("thick4.", "s04_thick."))]
    names = pyazr._azure2.parameter_names(p.to_bytes())
    check("the inactive data set's <parameterSettings> and <mcmc> rows follow the rename", rows == ["thick4.norm"] * 2, rows)
    check("activated: the rows name parameters the engine gives", all(r in names for r in rows), (rows, names[:8]))
    check("a row the engine cannot name is listed in a warning",
          any("natN_note" in str(w.message) for w in caught), [str(w.message) for w in caught])


def comments(work):
    print("3. set keeps comments between repeats")
    f3 = t2_dir(work, "T2set")
    with open(f3) as f:
        text = f.read().replace("    shift          offset 0 MeV fit prior 0 0.003\n",
                                "    shift          offset 0 MeV fit prior 0 0.003\n    # the sqrt term below\n"
                                "    shift          sqrt 0 MeV^1/2 fixed\n", 1)
    with open(f3, "w") as f:
        f.write(text)
    p = Project.open(f3)
    p["s01_int_thin"].set(shift=["offset 0.001 MeV fit prior 0 0.003", "sqrt 0 MeV^1/2 fixed"])
    lines = p.text.splitlines()
    i = next(k for k, l in enumerate(lines) if l.strip().startswith("shift") and "0.001" in l)
    check("the comment stays between the two repeats", lines[i + 1].strip() == "# the sqrt term below" and
          lines[i + 2].strip().startswith("shift") and "sqrt" in lines[i + 2], lines[i:i + 3])


def csv_records(work):
    print("4. records_from_csv")
    f3 = t2_dir(work, "T2csv")
    p = Project.open(f3)
    sheet = os.path.join(work, "datasets.csv")
    with open(sheet, "w") as f:
        f.write("name,reaction,quantity,unit,file,columns,norm\n"
                "csv1,\"15N(p,a)12C\",sigma,b,data/s01_int_thin.dat,E_lab[MeV] _ value error,1 fixed\n"
                "csv2,\"15N(p,a)12C\",sigma,parsec,data/s01_int_thin.dat,E_lab[MeV] _ value error,1 fixed\n"
                "# a comment row\n"
                "csv3,\"15N(p,a)12C\",sigma,b,data/s01_int_thin.dat,E_lab[MeV] _ value error,D\n"
                "csv4,\"15N(p,a)12C\",sigma,b,data/s01_int_thin.dat,E_lab[MeV] _ value error,\n")
    try:
        records_from_csv(sheet, project=p)
        e = None
    except FormatError as err:
        e = err
    rows = sorted({d["row"] for d in e.diagnostics}) if e else []
    check("every refused row at once, with its row and the engine's codes", e is not None and rows == [3, 6] and
          all(d["message"].startswith(f"row {d['row']}:") for d in e.diagnostics), (rows, e and e.codes))
    with open(sheet, "w") as f:
        f.write("name,reaction,quantity,unit,file,columns,norm\n"
                "csv1,\"15N(p,a)12C\",sigma,b,data/s01_int_thin.dat,E_lab[MeV] _ value error,1 fixed\n"
                "csv3,\"15N(p,a)12C\",sigma,b,data/s01_int_thin.dat,E_lab[MeV] _ value error,1.1 fit prior 1 5%\n")
    records = records_from_csv(sheet, project=p)
    for kind, name, fields in records:
        p.add(kind, name, **fields)
    check("a valid sheet: its records added, the project validates",
          [r[1] for r in records] == ["csv1", "csv3"] and "csv3" in p and p.validate() == [], p.validate())


def main():
    work = tempfile.mkdtemp(prefix="session_records.")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            session(work)
            rename(work)
            rename_inactive(work)
            comments(work)
            csv_records(work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print(f"{len(failures)} failure(s)" if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
