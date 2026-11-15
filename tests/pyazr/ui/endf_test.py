#!/usr/bin/env python3
"""The ENDF tab's routes and jobs against the pyazr.endf calls made directly (T1 with the two ENDF-6 files of
tests/pyazr/data/endf in its folder, the stand-in of the checking codes of tests/pyazr/endf_standin):

  1. GET /api/endf/files lists exactly the two ENDF-6 files (one in a subdirectory), not the project's data files;
  2. POST /api/endf/summary answers Material.summary().to_dict(); a path outside the root is refused (403), a file
     that is not ENDF-6 with pyazr's FormatError (422);
  3. POST /api/endf/curves: per data set the values of Material.curves(data, session=s) at the text's x0 bit for
     bit, the ratios to the datum times its normalization and to the model; where the file has nothing, the
     refusal's three parts (the quantity, the pair of another file); no file named is refused;
  4. GET /api/endf/checkers: the stand-in codes found through AZURE2_ENDF_UTILITIES, each with what it does to
     numbers (pyazr.endf.PROGRAMS); none found without the variable, with the sentence saying how codes are found;
  5. the job endf_check equals pyazr.endf.check on the same file (findings, counts, meanings of the known false
     alarms); a path outside the root is refused before the job starts;
  6. the job endf_standardize gives Candidate.diff(); POST /api/endf/accept writes a new file beside the original,
     refuses to write over an existing file without overwrite, replaces the original only with overwrite, and
     refuses a target outside the root and a job without a candidate.

Run: python3 tests/pyazr/ui/endf_test.py
"""
import os
import shutil
import sys

try:
    import pytest
except ImportError:  # ctest counts exit 77 as skipped
    print("skip: the package pytest is not installed (pip install pyazr[web])")
    sys.exit(77)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client, wait_job  # noqa: E402

DATA = os.path.join(uitest.ROOT, "tests", "pyazr", "data", "endf")
STANDIN = os.path.join(uitest.ROOT, "tests", "pyazr", "endf_standin", "standin.py")
HAND, EXTRACT = "p-C012_hand.endf", os.path.join("endf", "h-002_He_004_extract.endf")


@pytest.fixture()
def t1(tmp_path, monkeypatch):
    path = copy_project("T1", tmp_path)
    folder = os.path.dirname(path)
    shutil.copyfile(os.path.join(DATA, "p-C012_hand.endf"), os.path.join(folder, HAND))
    os.makedirs(os.path.join(folder, "endf"))
    shutil.copyfile(os.path.join(DATA, "h-002_He_004_extract.endf"), os.path.join(folder, EXTRACT))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name in ("checkr", "fizcon", "psyche", "stanef"):
        os.symlink(STANDIN, bin_dir / name)
    monkeypatch.delenv("AZURE2_PREPRO", raising=False)
    monkeypatch.setenv("AZURE2_ENDF_UTILITIES", str(bin_dir))
    monkeypatch.setenv("PATH", "/usr/bin:/bin")  # no code of the machine's PATH is found
    monkeypatch.setenv("ENDF_STANDIN_ROLE", "findings")
    with make_client(path) as client:
        client.path, client.folder, client.bin = path, folder, str(bin_dir)
        yield client


def test_files_and_summary(t1):
    import pyazr.endf as E
    listed = t1.get("/api/endf/files").json()
    assert [f["path"] for f in listed["files"]] == [EXTRACT, HAND], listed
    r = t1.post("/api/endf/summary", json={"path": HAND})
    assert r.status_code == 200, r.text
    assert r.json() == dict(E.read(os.path.join(t1.folder, HAND)).summary().to_dict(), name=HAND)
    assert t1.post("/api/endf/summary", json={"path": "../../outside.endf"}).status_code == 403
    bad = t1.post("/api/endf/summary", json={"path": "T1.azr"})
    assert bad.status_code == 422 and bad.json()["error"] == "FormatError", bad.text


def test_curves_against_material_curves(t1):
    import numpy as np
    import pyazr
    import pyazr.endf as E
    r = t1.post("/api/endf/curves", json={"files": [HAND, EXTRACT]})
    assert r.status_code == 200, r.text
    answer = r.json()
    assert answer["files"] == [{"name": HAND, "refused": None}, {"name": EXTRACT, "refused": None}]
    sets = {d["name"]: d for d in answer["datasets"]}
    hand = E.read(os.path.join(t1.folder, HAND))
    drawn = refused = 0
    with pyazr.Session(t1.path) as s:
        ev = s.evaluate(np.asarray(s.parameters.x0, float))
        for d in ev.data:
            entry = sets[d.name]
            try:
                direct = hand.curves(d, session=s)
            except pyazr.FormatRefused as e:
                got = entry["curves"][HAND]["refused"]
                assert (got["what"], got["allowed"], got["fix"]) == (e.what, e.allowed, e.fix), got
                refused += 1
                continue
            got = entry["curves"][HAND]
            want = [float(v) if np.isfinite(v) else None for v in direct.values]
            assert got["values"] == want and got["source"] == direct.source, d.name
            assert got["coverage"] == direct.coverage.tolist() and got["unit"] == entry["engine_unit"]
            for k, v in enumerate(got["values"]):
                if v is None:
                    assert got["ratio_to_data"][k] is None and got["ratio_to_model"][k] is None
                else:
                    assert got["ratio_to_data"][k] == v / entry["value_norm"][k]
                    assert got["ratio_to_model"][k] == v / entry["model"][k]
                    assert entry["model"][k] == float(d.model[k])
            assert "refused" in entry["curves"][EXTRACT]
            assert "projectile + target" in entry["curves"][EXTRACT]["refused"]["what"]
            drawn += 1
    assert drawn == 3 and refused == 7, (drawn, refused)
    none = t1.post("/api/endf/curves", json={"files": []})
    assert none.status_code == 422 and none.json()["allowed"], none.text
    assert t1.post("/api/endf/curves", json={"files": ["../../x.endf"]}).status_code == 403


def test_checkers(t1, monkeypatch):
    import pyazr.endf as E
    found = t1.get("/api/endf/checkers").json()
    assert [c["name"] for c in found["found"]] == ["checkr", "fizcon", "psyche", "stanef"], found
    stanef = next(c for c in found["found"] if c["name"] == "stanef")
    assert stanef["numbers"] == E.PROGRAMS["stanef"][3] and "7 significant digits" in stanef["numbers"]
    assert found["known"].keys() == E.PROGRAMS.keys() and found["search"]["AZURE2_ENDF_UTILITIES"] == t1.bin
    monkeypatch.delenv("AZURE2_ENDF_UTILITIES")
    none = t1.get("/api/endf/checkers").json()
    assert none["found"] == [] and "AZURE2_ENDF_UTILITIES" in none["how"] and "PATH" in none["how"]


def test_check_job_against_check(t1):
    import pyazr.endf as E
    r = t1.post("/api/jobs", json={"kind": "endf_check", "method": "codes", "settings": {"path": EXTRACT}})
    assert r.status_code == 200, r.text
    job = wait_job(t1, r.json()["id"])
    assert job["state"] == "finished", job
    got = job["result"]["check"]
    # the page's check includes the check of the resonance parameters (MF=2 against MF=3, masses and Q, the text)
    direct = E.check(os.path.join(t1.folder, EXTRACT), resonances=True).to_dict()
    assert got["own"]["findings"] == direct["own"]["findings"]
    assert got["own"]["counts"] == {"error": 3, "warning": 6, "note": 2}
    assert {"RECONSTRUCTION", "MASS", "Q_MASSES", "TEXT"} <= {f["code"] for f in got["own"]["findings"]}
    assert [x["program"] for x in got["reports"]] == ["checkr", "fizcon", "psyche"]
    for g, d in zip(got["reports"], direct["reports"]):
        assert g["findings"] == d["findings"] and g["totals"] == (list(d["totals"]) if d["totals"] else None)
    checkr = got["reports"][0]
    assert checkr["counts"] == {"error": 3, "warning": 0, "note": 1} and checkr["totals"] == [3, 0]
    lrp = next(f for f in checkr["findings"] if "LRP=2" in f["message"])
    assert lrp["meaning"].startswith("ENDF-102 p. 42 reserves LRP=2"), lrp
    ap = got["reports"][2]["findings"][0]
    assert "AP NOT IN RANGE" in ap["message"] and "not a radius" in ap["meaning"], ap
    one = t1.post("/api/jobs", json={"kind": "endf_check", "method": "codes",
                                     "settings": {"path": HAND, "codes": "checkr"}})
    assert [x["program"] for x in wait_job(t1, one.json()["id"])["result"]["check"]["reports"]] == ["checkr"]
    outside = t1.post("/api/jobs", json={"kind": "endf_check", "method": "codes",
                                         "settings": {"path": "../../x.endf"}})
    assert outside.status_code == 403, outside.text


def test_standardize_and_accept(t1, monkeypatch):
    import pyazr.endf as E
    monkeypatch.setenv("ENDF_STANDIN_ROLE", "lossy")
    original = os.path.join(t1.folder, EXTRACT)
    with open(original) as f:
        before = f.read()
    r = t1.post("/api/jobs", json={"kind": "endf_standardize", "method": "program",
                                   "settings": {"path": EXTRACT, "program": "stanef"}})
    job = wait_job(t1, r.json()["id"])
    assert job["state"] == "finished", job
    cand = job["result"]["candidate"]
    direct = E.standardize(original, program="stanef")
    import dataclasses
    assert cand["diff"] == [dataclasses.asdict(c) for c in direct.diff()]
    assert any(c["kind"] == "numbers" for c in cand["diff"]) and cand["numbers"] == E.PROGRAMS["stanef"][3]
    with open(original) as f:
        assert f.read() == before  # nothing replaced by the run
    kept = t1.post("/api/endf/accept", json={"job_id": job["id"]})
    assert kept.status_code == 200, kept.text
    assert kept.json() == {"path": os.path.join("endf", "h-002_He_004_extract.stanef.endf"), "replaced": False,
                           "overwrite": False}
    with open(os.path.join(t1.folder, "endf", "h-002_He_004_extract.stanef.endf")) as f, open(cand["path"]) as g:
        assert f.read() == g.read()
    again = t1.post("/api/endf/accept", json={"job_id": job["id"]})
    assert again.status_code == 422 and again.json()["allowed"], again.text
    named = t1.post("/api/endf/accept", json={"job_id": job["id"], "to": "kept/new.endf"})
    assert named.status_code == 422 and "no directory" in named.json()["what"], named.text
    assert t1.post("/api/endf/accept", json={"job_id": job["id"], "to": "../../x.endf"}).status_code == 403
    project = os.path.basename(t1.path)
    not_endf = t1.post("/api/endf/accept", json={"job_id": job["id"], "to": project, "overwrite": True})
    assert not_endf.status_code == 422 and "not an ENDF-6 file" in not_endf.json()["what"], not_endf.text
    no_overwrite = t1.post("/api/endf/accept", json={"job_id": job["id"], "to": EXTRACT})
    assert no_overwrite.status_code == 422
    with open(original) as f:
        assert f.read() == before
    replaced = t1.post("/api/endf/accept", json={"job_id": job["id"], "to": EXTRACT, "overwrite": True})
    assert replaced.status_code == 200 and replaced.json()["replaced"] is True, replaced.text
    with open(original) as f:
        assert f.read() != before
    assert len(t1.get(f"/api/jobs/{job['id']}").json()["kept"]) == 2
    check = t1.post("/api/jobs", json={"kind": "endf_check", "method": "codes", "settings": {"path": HAND}})
    check_job = wait_job(t1, check.json()["id"])
    refused = t1.post("/api/endf/accept", json={"job_id": check_job["id"]})
    assert refused.status_code == 422 and "no candidate" in refused.json()["what"], refused.text


if __name__ == "__main__":
    uitest.run(__file__)
