#!/usr/bin/env python3
"""pyazr.Project: a format-3 project as records (docs/upgrade/project_format3.md sections 4 and 6).

  1. round trip: every project of tests/ that pyazr.convert_azr(..., to=3, exact_grids=True) converts, and every
     format-3 file of tests/, gives its own bytes through Project.open(f).save(g);
  2. one edit changes one line: set(norm=...) on a data set of T2's conversion; the diff is that line;
  3. rename updates the header, every reference (the validator's), the <parameterSettings> rows and the names of the
     prior-covariance file (written by save); remove of a referenced record is refused naming the referrers (the
     engine's E_UNKNOWN_REFERENCE), of an unreferenced one leaves a valid project;
  4. every refusal is the engine's: a missing required field (E_FIELD_REQUIRED), a unit of another dimension
     (E_UNIT_NOT_ALLOWED), a reference to no record (E_UNKNOWN_REFERENCE), a name in use (E_DUPLICATE_NAME); the
     project is unchanged after a refusal;
  5. user stories built from scratch in Python (Project.new on the committed file's nuclear part, then add_*) give the
     bytes of the committed project's outputs (menu 1 of the CLI, every file of output/): story (a) yield per charge
     (tests/units_scales/models/yield_story.azr), story (b) ratio to Rutherford (tests/ratio_to_rutherford/models/
     aC_story.azr), a potential with a spin-orbit term (tests/spin_orbit_potential/spin_orbit_potential.azr); and an
     S-factor data set (no committed project has one): built in Python against the same records written as file
     text by this script, on the nuclear part of tests/spin_orbit_potential/models/capture.azr;
  6. help(Project.add_dataset) lists every field key of the schema's dataset record (generated, not written);
  7. describe(name) is the engine's description of that record; a target is refused naming its users;
  8. session(): a project with unsaved edits is saved and run (the chi^2 of its data set as the CLI's).

Needs the compiled module and AZURE2_BIN (ctest sets it).  Run:  AZURE2_BIN=<build>/src/AZURE2 python3 project_test.py
"""
import difflib
import filecmp
import glob
import inspect
import os
import re
import shutil
import subprocess
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
TESTS = os.path.join(ROOT, "tests")
FIXTURES = os.path.join(ROOT, "convert", "tests", "fixtures")
BIN = os.environ.get("AZURE2_BIN")
sys.path.insert(0, ROOT)

import pyazr  # noqa: E402
from pyazr import Project, Nuisance, FormatError  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def refused(call):
    """The FormatError a call raises, or None."""
    try:
        call()
    except FormatError as e:
        return e
    return None


def run_cli(directory, azr):
    """Menu 1 (calculate with data) of the CLI on directory/azr, one thread; the output directory's files."""
    for d in ("output", "checks"):
        shutil.rmtree(os.path.join(directory, d), ignore_errors=True)
        os.makedirs(os.path.join(directory, d))
    subprocess.run([BIN, "--no-gui", "--no-readline", azr], input="1\n\n\n7\n", text=True, cwd=directory,
                   capture_output=True, env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
    out = os.path.join(directory, "output")
    return out, sorted(os.listdir(out))


def same_outputs(a, b):
    """(identical, files, differing): every file of output/ of two runs, byte for byte."""
    oa, fa = a
    ob, fb = b
    if fa != fb or not fa:
        return False, len(fa), [f"files {fa} != {fb}"]
    bad = [f for f in fa if not filecmp.cmp(os.path.join(oa, f), os.path.join(ob, f), shallow=False)]
    return not bad, len(fa), bad


def workdir(parent, name, data):
    d = os.path.join(parent, name)
    os.makedirs(d)
    if data:
        os.symlink(data, os.path.join(d, "data"))
    return d


def story(tmp, name, committed, build):
    """Build the story in Python on the committed file's nuclear part; run both; compare the outputs."""
    src_dir = os.path.dirname(committed)
    a = workdir(tmp, name + "_committed", os.path.join(src_dir, "data"))
    shutil.copy(committed, os.path.join(a, "m.azr"))
    b = workdir(tmp, name + "_python", os.path.join(src_dir, "data"))
    p = Project.new(os.path.join(b, "m.azr"), nuclear_part=committed)
    build(p)
    p.save()
    same, n, bad = same_outputs(run_cli(a, "m.azr"), run_cli(b, "m.azr"))
    check(f"story {name}: built in Python, the {n} output files of the committed project byte for byte", same, bad)
    return p


def main():
    if not BIN:
        print("AZURE2_BIN is not set; skipped")
        return 77
    tmp = tempfile.mkdtemp(prefix="project_test.")
    try:
        sections(tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"{len(failures)} failure(s)" if failures else "PASS")
    return 1 if failures else 0


def sections(tmp):
    print("1. round trip")
    conv = os.path.join(tmp, "converted")
    sources = sorted(glob.glob(os.path.join(TESTS, "**", "*.azr"), recursive=True))
    files, not_converted = [], 0
    for i, f in enumerate(sources):
        with open(f, "rb") as h:
            fmt = pyazr._azure2.project_format(h.read())
        if fmt == 3:
            files.append(f)
            continue
        out = os.path.join(conv, f"{i:03d}_" + os.path.basename(f))
        os.makedirs(os.path.dirname(out), exist_ok=True)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pyazr.convert_azr(f, out, to=3, exact_grids=True)
            files.append(out)
        except (ValueError, FileExistsError):
            not_converted += 1
    bad = []
    for f in files:
        g = os.path.join(tmp, "roundtrip.azr")
        Project.open(f).save(g)
        with open(f, "rb") as x, open(g, "rb") as y:
            if x.read() != y.read():
                bad.append(f)
    check(f"Project.open(f).save(g) == f for {len(files)} format-3 files ({len(sources)} projects of tests/, "
          f"{not_converted} the converter refuses without answers)", not bad and len(files) > 150, bad[:5])
    # The converter ends every line with LF: a variant with CR LF and no final line ending, unedited and edited.
    with open(os.path.join(FIXTURES, "T2.f3.azr"), "rb") as f:
        crlf = f.read().replace(b"\n", b"\r\n").rstrip(b"\r\n")
    v = os.path.join(tmp, "crlf.azr")
    with open(v, "wb") as f:
        f.write(crlf)
    q = Project.open(v)
    q.save(os.path.join(tmp, "crlf2.azr"))
    check("CR LF and no final line ending: written back byte for byte", open(os.path.join(tmp, "crlf2.azr"), "rb").read() == crlf)
    q["s01_int_thin"].set(norm="1.05 fit prior 1 5%")
    edited = q.to_bytes()
    check("  an edit keeps the line's CR LF and changes nothing else",
          edited == crlf.replace(b"    norm           1 fit prior 1 5%\r\n", b"    norm           1.05 fit prior 1 5%\r\n", 1)
          and edited != crlf)

    print("2. one edit, one line")
    t2 = workdir(tmp, "T2", os.path.join(TESTS, "upgrade", "T2", "data"))
    shutil.copy(os.path.join(FIXTURES, "T2.f3.azr"), os.path.join(t2, "T2.azr"))
    shutil.copy(os.path.join(FIXTURES, "T2.prior_cov.f3.txt"), os.path.join(t2, "prior_cov.f3.txt"))
    shutil.copy(os.path.join(TESTS, "upgrade", "T2", "cov_diff.txt"), t2)
    p = Project.open(os.path.join(t2, "T2.azr"))
    before = p.text.splitlines(True)
    p["s01_int_thin"].set(norm="1.05 fit prior 1 5%")
    diff = [d for d in difflib.unified_diff(before, p.text.splitlines(True), n=0) if d[:1] in "+-" and d[:3] not in ("+++", "---")]
    check("set(norm=...) on s01_int_thin: the diff is one line", diff == ["-    norm           1 fit prior 1 5%\n",
                                                                         "+    norm           1.05 fit prior 1 5%\n"], diff)
    check("the record reads the new text", p["s01_int_thin"].norm == "1.05 fit prior 1 5%", p["s01_int_thin"].norm)
    check("a field not given reads None; a target's 'kind' is its field", p["s01_int_thin"].angle_window is None
          and p["tgt_s01_int_thin"].kind == "layer" and p["tgt_s01_int_thin"].record_kind == "target")

    print("3. rename and remove")
    t2r = workdir(tmp, "T2rename", os.path.join(TESTS, "upgrade", "T2", "data"))  # rename rewrites the covariance file
    for f in ("T2.azr", "prior_cov.f3.txt", "cov_diff.txt"):
        shutil.copy(os.path.join(t2, f), t2r)
    with open(os.path.join(t2r, "T2.azr")) as f:  # a <parameterSettings> row naming a parameter of s07_profile
        text = f.read().replace("<parameterSettings>\n", "<parameterSettings>\ns07_profile.norm 0.94 0 0 0.01 0 1 norm 0\n")
    with open(os.path.join(t2r, "T2p.azr"), "w") as f:
        f.write(text)
    p = Project.open(os.path.join(t2r, "T2p.azr"))
    p.rename("s07_profile", "profile7")
    text = p.text
    left = re.findall(r"(?<![\w/-])s07_profile\b(?!\.dat)", text)
    check("rename: the header renamed and no name or parameter of 's07_profile' left", "dataset profile7" in text and
          not left and "\nprofile7.norm 0.94" in text, [l for l in text.splitlines() if "s07_profile" in l])
    check("rename: the data file named after the record is not touched", "data/s07_profile.dat" in text)
    p.save(os.path.join(t2r, "T2r.azr"))
    with open(os.path.join(t2r, "prior_cov.f3.txt")) as f:
        names = [l for l in f if l.strip() and not l.startswith("#")][0].split()
    check("rename: the prior-covariance file's names follow (written by save)", names == ["profile7.norm", "s08_sum.norm"],
          names)
    check("rename: the project validates", p.validate() == [], p.validate())
    p.rename("D", "Dnorm")
    users = [r for r in pyazr._azure2.record_references(p.to_bytes()) if r["name"] == b"Dnorm"]
    check("rename of a normalization: every reference renamed (2)", len(users) == 2 and "Dnorm" in p["s02_diff090"].norm)
    e = refused(lambda: p.remove("Dnorm"))
    check("remove of a referenced record: refused, naming the referrers", e is not None and
          "s02_diff090.norm" in str(e) and "s03_diff150.norm" in str(e) and set(e.codes) == {"E_UNKNOWN_REFERENCE"},
          str(e))
    check("  and the project is unchanged", "Dnorm" in p)
    p.remove("s13_thm_oes")
    check("remove of an unreferenced data set: gone, the project valid", "s13_thm_oes" not in p and p.validate() == [],
          p.validate())
    e = refused(lambda: p.rename("s01_int_thin", "s02_diff090"))
    check("rename to a name in use: E_DUPLICATE_NAME", e is not None and "E_DUPLICATE_NAME" in e.codes, str(e))

    print("4. the engine's refusals")
    p = Project.open(os.path.join(t2, "T2.azr"))
    snapshot = p.to_bytes()
    common = dict(reaction="15N(p,a)12C", quantity="sigma", file="data/s01_int_thin.dat",
                  columns="E_lab[MeV] _ value error", norm=Nuisance(1, "fixed"))
    e = refused(lambda: p.add_dataset("x1", **common))
    check("a missing required field (unit): E_FIELD_REQUIRED", e is not None and e.codes == ["E_FIELD_REQUIRED"] and
          e.diagnostics[0]["field"] == "unit", str(e))
    e = refused(lambda: p.add_dataset("x2", unit="MeV", **common))
    check("a unit of another dimension (MeV for sigma): E_UNIT_NOT_ALLOWED", e is not None and
          e.codes == ["E_UNIT_NOT_ALLOWED"], str(e))
    e = refused(lambda: p.add_dataset("x3", unit="mb", target="nowhere", **common))
    check("a reference to no record: E_UNKNOWN_REFERENCE", e is not None and e.codes == ["E_UNKNOWN_REFERENCE"], str(e))
    e = refused(lambda: p["s01_int_thin"].set(norm="1 fitted"))
    check("set: a nuisance of the wrong grammar is refused with the engine's code", e is not None and e.codes and
          all(c.startswith("E_") for c in e.codes), str(e))
    e = refused(lambda: p.add_dataset("D", unit="mb", **common))
    check("a name in use: E_DUPLICATE_NAME", e is not None and "E_DUPLICATE_NAME" in e.codes, str(e))
    check("the project is unchanged after every refusal", p.to_bytes() == snapshot)
    p.add_dataset("x4", unit="mb", **common)
    check("the same record with its unit is added", p["x4"].unit == "mb" and p.validate() == [], p.validate())

    print("5. user stories from scratch")
    os.environ.setdefault("OMP_NUM_THREADS", "1")

    def story_a(p):
        p.add_target("N15T", kind="infinite", active="15N", composition="N 1",
                     stopping='formula "a0*x^a1" a0=4.5e-20 a1=-0.6 MeV*cm2/atom', straggling="none")
        p.add_beam("pbeam", spread=(0.3, "keV", "lab"))
        y = dict(reaction="15N(p,a)12C", quantity="yield", unit="counts/uC", efficiency=Nuisance(0.0182, "fixed"),
                 file="data/yield_curve.dat", columns="E_lab[keV] _ value error",
                 norm=Nuisance(1, "fit", prior=(1, "3%")), target="N15T", beam="pbeam", quadrature="points 60")
        e = refused(lambda: p.add_dataset("Y_curve", **y))
        check("story (a): without charge_state the engine refuses (E_FIELD_REQUIRED, charge_state)",
              e is not None and e.codes == ["E_FIELD_REQUIRED"] and e.diagnostics[0]["field"] == "charge_state", str(e))
        p.add_dataset("Y_curve", charge_state=1, **y)
        p.add_extrapolation("Y_thick", reaction="15N(p,a)12C", quantity="yield", unit="1/particle",
                            energies="range 0.6 1.4 MeV lab step 0.01", target="N15T", beam="pbeam",
                            quadrature="points 60")

    story(tmp, "a_yield", os.path.join(TESTS, "units_scales", "models", "yield_story.azr"), story_a)

    def story_b(p):
        p.add_dataset("aC_RTH", reaction="12C(a,a)12C", quantity="ratio_to_rutherford", unit=1, file="data/aC_rth.dat",
                      columns="E_lab[MeV] theta_lab[deg] value error", norm="1 fit prior 1 3%")
        p.add_extrapolation("aC_RTH_curve", reaction="12C(a,a)12C", quantity="ratio_to_rutherford", unit=1,
                            energies="range 3.0 6.0 MeV lab step 0.05", angles="list 30 60 90 120 150 deg lab")
        e = refused(lambda: p.add_dataset("aa", reaction="4He(a,a)4He", quantity="ratio_to_rutherford", unit=1,
                                          file="data/aC_rth.dat", columns="E_lab[MeV] theta_lab[deg] value error",
                                          norm="1 fixed"))
        check("story (b): 4He(a,a)4He as ratio_to_rutherford is refused by the engine", e is not None and e.codes, str(e))

    story(tmp, "b_rutherford", os.path.join(TESTS, "ratio_to_rutherford", "models", "aC_story.azr"), story_b)

    def story_so(p):
        p.add_potential("p12C", pair="p+12C", coulomb_radius=(3.5, "fm"),
                        central="woods_saxon depth 20 MeV radius 3.6 fm diffuseness 0.6 fm",
                        spin_orbit="woods_saxon depth 6 MeV radius 3.6 fm diffuseness 0.6 fm on light operator 2l.s",
                        matching_radius=(15, "fm"), step=(0.0275, "fm"))
        p.add_dataset("meyer_144", reaction="12C(p,p)12C", quantity="dsigma_domega", unit="b/sr",
                      file="data/meyer_144.dat", columns="E_lab[MeV] theta_lab[deg] value error",
                      energy_window=(0, 5, "MeV"), angle_window=(144.1, 144.1, "deg"), norm="1 fixed")

    story(tmp, "spin_orbit", os.path.join(TESTS, "spin_orbit_potential", "spin_orbit_potential.azr"), story_so)

    # An S-factor data set: no committed project has one; the Python-built project against the same records as text.
    capture = os.path.join(TESTS, "spin_orbit_potential", "models", "capture.azr")
    d = workdir(tmp, "sfactor_python", None)
    os.makedirs(os.path.join(d, "data"))
    rows = [l.split() for l in open(os.path.join(os.path.dirname(capture), "data", "skowronski_luna_hpge.dat"))
            if l.strip() and not l.startswith("#")]
    conv = pyazr.physics.sfactor_conversion(1, 6, 1.00728, 12.0, [float(r[0]) * 12 / 13.00728 for r in rows])
    with open(os.path.join(d, "data", "s.dat"), "w") as f:
        for r, c in zip(rows, conv):
            f.write(f"{r[0]} 0 {float(r[2]) * float(c) * 1e3!r} {float(r[3]) * float(c) * 1e3!r}\n")
    p = Project.new(os.path.join(d, "m.azr"), nuclear_part=capture)
    nuclear = p.text
    p.add_dataset("S_luna", reaction="12C(p,g)13N", quantity="sfactor", unit="keV*b", file="data/s.dat",
                  columns="E_lab[MeV] _ value error", norm=Nuisance(1, "fit", prior=(1, "10%")),
                  shift=("offset", Nuisance(0, "fit", unit="keV", prior=(0, 2))))
    p.save()
    h = workdir(tmp, "sfactor_text", os.path.join(d, "data"))
    records = ("dataset S_luna\n    reaction       12C(p,g)13N\n    quantity       sfactor\n    unit           keV*b\n"
               "    file           data/s.dat\n    columns        E_lab[MeV] _ value error\n"
               "    norm           1 fit prior 1 10%\n    shift          offset 0 keV fit prior 0 2\n")
    with open(os.path.join(h, "m.azr"), "w") as f:
        f.write(re.sub(r"<data>\n", "<data>\n" + records, nuclear, count=1))
    check("S factor: the Python-built record is the hand-written one", p["S_luna"].quantity == "sfactor")
    same, n, bad = same_outputs(run_cli(h, "m.azr"), run_cli(d, "m.azr"))
    check(f"story sfactor: built in Python, the {n} output files of the same records written as text byte for byte",
          same, bad)

    print("6. help() from the table")
    import pydoc
    text = pydoc.render_doc(Project.add_dataset, renderer=pydoc.plaintext)
    keys = list(pyazr.schema()["records"]["dataset"]["fields"])
    sig = str(inspect.signature(Project.add_dataset))
    missing = [k for k in keys if f"\n    {k}: " not in Project.add_dataset.__doc__ or not re.search(rf"\b{k}\b", sig)
               or not re.search(rf"^ *{re.escape(k)}: ", text, re.M)]
    check(f"help(add_dataset): each of the {len(keys)} field keys of the schema's dataset in the signature and with its "
          f"own line of help", not missing, missing)
    f = pyazr.schema()["records"]["dataset"]["fields"]["charge_state"]
    check("  the help line is the table's (charge_state: title, grammar)", f"    charge_state: {f['title']}; {f['grammar']}"
          in Project.add_dataset.__doc__)
    check("every record kind of the schema has its add_<kind>",
          all(hasattr(Project, "add_" + k) for k in pyazr.schema()["records"]))

    print("7. describe")
    p = Project.open(os.path.join(t2, "T2.azr"))
    one = p.describe("s09_ratio")
    check("describe(name): the engine's block of that record only", one.startswith("dataset s09_ratio") and
          "dataset s01_int_thin" not in one, one[:200])
    j = p.describe("s09_ratio", as_json=True)
    check("describe(name, as_json=True): its entry", j.get("name") == "s09_ratio" and j.get("rows_read") == 20, j)
    try:
        p.describe("D")
        check("describe of a normalization is refused", False)
    except ValueError as e:
        check("describe of a normalization is refused naming its users", "s02_diff090" in str(e), str(e))

    print("8. session")
    p = Project.open(os.path.join(t2, "T2.azr"))
    p["s01_int_thin"].set(norm="1.05 fixed")
    os.makedirs(os.path.join(t2, "output"), exist_ok=True)
    with p.session() as m:
        chi2 = m.evaluate().chi2
    check("session(): the unsaved edit is saved first", not p.unsaved)
    out, _ = run_cli(t2, "T2.azr")
    cli = float(re.findall(r"Total-Chi-Squared:\s*(\S+)", open(os.path.join(out, "chiSquared.out")).read())[-1])
    check("session(): the chi^2 of the edited project is the CLI's (the 6 digits chiSquared.out prints)",
          abs(chi2 / cli - 1) < 1e-5 and "1.05 fixed" in open(p.path).read(), (chi2, cli))


if __name__ == "__main__":
    sys.exit(main())
