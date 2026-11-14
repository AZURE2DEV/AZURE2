#!/usr/bin/env python3
"""The ENDF tab's two ways between a project and an evaluation, against the pyazr.endf calls made directly (T1 with
the ENDF/B-VIII.1 3He + 4He extract and the hand-written p + 12C file of tests/pyazr/data/endf in its folder):

  1. POST /api/endf/resonances: what MF=2 holds (the record's text, the pairs, the spin groups) equals
     ``read(path).resonances()``; a file without LRF=7 parameters and an isotope the file lacks are refused with the
     three parts; nothing is written;
  2. the job endf_import gives the report of ``pyazr.endf.to_project`` made directly (pairs, levels, masses against
     AME2020, notes) and the same project text; the project folder is unchanged until the project is kept;
  3. POST /api/endf/keep_project writes the job's project as a new file of the folder; it refuses an existing file
     (the open project among them), a name without .azr, a path outside the root (403), a job without a project;
     the kept project then opens as the server's project (POST /api/project/open);
  4. GET /api/endf/export: the losses and photon widths of ``pyazr.endf.from_session(s, x0)``, the data's energy
     range, a file name not in use;
  5. the job endf_write writes the file ``pyazr.endf.write`` writes (the same bytes, the same sections, grids and
     MF=2 check) as a new file of the folder, listed by GET /api/endf/files; a second write to the same name fails
     with the writer's refusal, a name outside the root is refused before the job starts.

Run: python3 tests/pyazr/ui/endf_project_test.py
"""
import dataclasses
import os
import shutil
import sys
import tempfile
import warnings

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client, wait_job  # noqa: E402

DATA = os.path.join(uitest.ROOT, "tests", "pyazr", "data", "endf")
HAND, EXTRACT = "p-C012_hand.endf", "h-002_He_004_extract.endf"


@pytest.fixture()
def t1(tmp_path):
    path = copy_project("T1", tmp_path)
    folder = os.path.dirname(path)
    for name in (HAND, EXTRACT):
        shutil.copyfile(os.path.join(DATA, name), os.path.join(folder, name))
    with make_client(path) as client:
        client.path, client.folder = path, folder
        yield client


def listing(folder):
    """The user's files of the folder (not the program's output, not the server's hidden files)."""
    out = []
    for d, _, files in os.walk(folder):
        for f in files:
            rel = os.path.relpath(os.path.join(d, f), folder)
            if not any(part.startswith(".") or part == "output" for part in rel.split(os.sep)):
                out.append(rel)
    return sorted(out)


def test_resonances_and_new_project(t1):
    import pyazr.endf as E
    before = listing(t1.folder)
    r = t1.post("/api/endf/resonances", json={"path": EXTRACT})
    assert r.status_code == 200, r.text
    answer = r.json()
    record = E.read(os.path.join(t1.folder, EXTRACT)).resonances()
    assert answer["refused"] is None and answer["text"] == str(record)
    assert [p["mt"] for p in answer["pairs"]] == [p.mt for p in record.pairs]
    assert [(g["jpi"], g["levels"]) for g in answer["groups"]] == [(g.jpi, len(g.levels)) for g in record.groups]
    for body, word in (({"path": HAND}, ""), ({"path": EXTRACT, "isotope": 3}, "isotope")):
        refused = t1.post("/api/endf/resonances", json=body).json()["refused"]
        assert refused and refused["what"] and refused["allowed"] and refused["fix"], refused
        assert word in refused["what"]
    assert t1.post("/api/endf/resonances", json={"path": "../../outside.endf"}).status_code == 403

    started = t1.post("/api/jobs", json={"kind": "endf_import", "method": "to_project",
                                         "settings": {"path": EXTRACT, "name": "he.azr"}})
    assert started.status_code == 200, started.text
    job = wait_job(t1, started.json()["id"])
    assert job["state"] == "finished", job
    result = job["result"]["project"]
    assert listing(t1.folder) == before, "nothing is written into the folder before the project is kept"
    work = tempfile.mkdtemp()
    m = E.read(os.path.join(t1.folder, EXTRACT))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        project, conversion = E.to_project(m.resonances(), os.path.join(work, "he.azr"), material=m)
    expected = dataclasses.asdict(conversion)
    assert {k: v for k, v in result["report"].items() if k != "path"} == \
        {k: (list(v) if isinstance(v, tuple) else v) for k, v in uitest_plain(expected).items() if k != "path"}
    assert result["text"] == project.text and result["report_text"].splitlines()[1:] == \
        str(conversion).splitlines()[1:]

    for body, status in (({"job_id": job["id"], "to": "T1.azr"}, 422), ({"job_id": job["id"], "to": "he.txt"}, 422),
                         ({"job_id": job["id"], "to": "../../he.azr"}, 403)):
        r = t1.post("/api/endf/keep_project", json=body)
        assert r.status_code == status, (body, r.text)
        assert status == 403 or r.json()["what"]
    kept = t1.post("/api/endf/keep_project", json={"job_id": job["id"]})
    assert kept.status_code == 200 and kept.json()["path"] == "he.azr", kept.text
    with open(os.path.join(t1.folder, "he.azr")) as f:
        assert f.read() == project.text
    again = t1.post("/api/endf/keep_project", json={"job_id": job["id"]})
    assert again.status_code == 422 and "exists" in again.json()["what"]
    other = t1.post("/api/jobs", json={"kind": "endf_check", "method": "codes", "settings": {"path": HAND}}).json()
    wait_job(t1, other["id"])
    assert t1.post("/api/endf/keep_project", json={"job_id": other["id"]}).status_code == 422
    opened = t1.post("/api/project/open", json={"path": kept.json()["open_path"]})
    assert opened.status_code == 200 and opened.json()["path"].endswith("he.azr"), opened.text
    shutil.rmtree(work, ignore_errors=True)


def uitest_plain(value):
    """A value as JSON gives it back: tuples as lists, recursively."""
    if isinstance(value, dict):
        return {k: uitest_plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [uitest_plain(v) for v in value]
    return value


def test_export(t1):
    import numpy as np
    import pyazr
    import pyazr.endf as E
    warnings.simplefilter("ignore", pyazr.EngineWarning)
    preview = t1.get("/api/endf/export").json()
    s = pyazr.Session(t1.path)
    x0 = np.asarray(s.parameters.x0, float)
    _, losses = E.from_session(s, x0)
    assert preview["refused"] is None and preview["losses"] == list(losses.items)
    assert [w["width_eV"] for w in preview["photon_widths"]] == [w[-1] for w in losses.photon_widths]
    lab = [e for d in s.evaluate(x0).data for e in d.extra["energy_lab"]]
    assert preview["energy_range_MeV"] == [min(lab), max(lab)] and preview["name"] == "T1.endf"
    assert preview["not_written"] == list(E.writer.NOT_WRITTEN)

    settings = {"name": "T1.endf", "authors": "pyazr tests", "lab": "TEST", "mat": 625}
    r = t1.post("/api/jobs", json={"kind": "endf_write", "method": "writer", "settings": settings})
    assert r.status_code == 200, r.text
    job = wait_job(t1, r.json()["id"])
    assert job["state"] == "finished", job
    written = job["result"]["written"]
    direct = os.path.join(tempfile.mkdtemp(), "T1.endf")
    w = E.write(s, direct, x0, header={"AUTH": "pyazr tests", "ALAB": "TEST", "MAT": 625})
    with open(os.path.join(t1.folder, "T1.endf")) as a, open(direct) as b:
        assert a.read() == b.read(), "the job's file is the direct call's"
    assert written["path"] == "T1.endf" and written["sections"] == [list(x) for x in w.sections]
    assert written["mf2_check"] == {str(k): v for k, v in w.mf2_check.items()}
    assert written["losses"] == list(w.losses) and written["not_written"] == list(w.not_written)
    assert "T1.endf" in [f["path"] for f in t1.get("/api/endf/files").json()["files"]]
    assert t1.get("/api/endf/export").json()["name"] == "T1_2.endf"
    again = wait_job(t1, t1.post("/api/jobs", json={"kind": "endf_write", "method": "writer",
                                                    "settings": settings}).json()["id"])
    assert again["state"] == "failed" and "exists" in again["result"]["message"], again["result"]
    outside = t1.post("/api/jobs", json={"kind": "endf_write", "method": "writer",
                                         "settings": {"name": "../../T1.endf"}})
    assert outside.status_code == 403, outside.text
    s.close()


if __name__ == "__main__":
    uitest.run(__file__)
