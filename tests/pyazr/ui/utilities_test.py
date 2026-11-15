#!/usr/bin/env python3
"""The page's utilities through their routes, against the pyazr calls made directly (T1):

  1. the covariance builder: GET /api/covariance/{name} proposes a file and the parts to start from; POST
     /api/covariance/build answers pyazr.covariance.build's numbers; a scale part alone comes back refused, naming
     it; POST /api/covariance/write writes the file pyazr writes and sets the field (the text equals Project.set's),
     the objective then equals a session's on that text; the builder reopens with the parts; one Undo removes the
     file and restores the text byte for byte, Redo writes the same bytes again; a refused matrix is not written;
  2. the rate job: the rate of 12C(p,g) at the text's x0 equals pyazr.rate's at the same vector bit for bit, with the
     Gamow window of pyazr.rate.gamow_window, the CSV in the output directory and the ratios; the libraries' files
     refused off their grids (STARLIB without a band, REACLIB off 0.01-10 GK) and, on the table grid, the REACLIB
     file pyazr.rate_formats writes from the same rate; an unknown grid and a
     covariance band without a fit are refused with their three parts.

Run: python3 tests/pyazr/ui/utilities_test.py
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client, wait_job  # noqa: E402

PARTS = [{"kind": "own"}, {"kind": "scale", "percent": 4}, {"kind": "length", "percent": 2, "length": 0.3}]


@pytest.fixture()
def t1(tmp_path):
    path = copy_project("T1", tmp_path)
    with make_client(path) as client:
        client.path = path
        yield client


def rev(client):
    return client.get("/api/project").json()["rev"]


def test_covariance_builder(t1):
    import numpy as np
    import pyazr
    from pyazr import covariance
    directory = os.path.dirname(t1.path)
    opened = t1.get("/api/covariance/el_lab090").json()
    assert opened["proposed"] == "cov_el_lab090.txt" and opened["parts"] == [{"kind": "own"}], opened
    assert opened["covariance"] is None and set(opened["kinds"]) == set(covariance.PART_KINDS), opened
    r = t1.post("/api/covariance/build", json={"datasets": ["el_lab090"], "parts": PARTS})
    assert r.status_code == 200, r.text
    built = r.json()
    direct = covariance.build(pyazr.Project.open(t1.path), ["el_lab090"], PARTS)
    assert built["refusal"] is None and built["rows"]["total_error"] == direct.total_error.tolist()
    assert built["correlation"]["matrix"] == direct.correlation.tolist()
    alone = t1.post("/api/covariance/build", json={"datasets": ["el_lab090"],
                                                    "parts": [{"kind": "scale", "percent": 4}]}).json()
    assert alone["refusal"] and alone["refusal_parts"] == [0], alone
    bad = t1.post("/api/covariance/build", json={"datasets": ["no_such"], "parts": PARTS})
    assert bad.status_code == 422 and "el_lab090" in bad.json()["allowed"], bad.text

    before = t1.get("/api/project/text").text
    refused = t1.post("/api/covariance/write", json={"datasets": ["el_lab090"], "path": "cov_x.txt",
                                                      "parts": [{"kind": "scale", "percent": 4}], "rev": rev(t1)})
    assert refused.status_code == 422 and not os.path.exists(os.path.join(directory, "cov_x.txt")), refused.text
    r = t1.post("/api/covariance/write", json={"datasets": ["el_lab090"], "parts": PARTS,
                                               "path": "cov_el_lab090.txt", "rev": rev(t1)})
    assert r.status_code == 200, r.text
    file = os.path.join(directory, "cov_el_lab090.txt")
    with open(file) as f:
        written = f.read()
    assert written == direct.text()
    expect = pyazr.Project.open(t1.path)
    expect.set("el_lab090", covariance="cov_el_lab090.txt")
    after = t1.get("/api/project/text").text
    assert after == expect.text
    again = t1.get("/api/covariance/el_lab090").json()
    assert again["covariance"] == "cov_el_lab090.txt" and not again["hand_written"]
    assert again["parts"] == direct.parts, again["parts"]
    objective = t1.post("/api/evaluate", json={}).json()["objective"]
    saved = os.path.join(directory, "with_cov.azr")
    expect.save(saved)
    with pyazr.Session(saved) as s:
        assert objective == s.objective()(np.asarray(s.parameters.x0, float))
    assert t1.post("/api/project/undo", json={"rev": rev(t1)}).status_code == 200
    assert t1.get("/api/project/text").text == before and not os.path.exists(file)
    assert t1.post("/api/project/redo", json={"rev": rev(t1)}).status_code == 200
    with open(file) as f:
        assert f.read() == written
    assert t1.get("/api/project/text").text == after


def test_rate_job(t1):
    import numpy as np
    from pyazr import rate
    settings = {"entrance": 1, "exit": 2, "grid": "custom", "t_min": 0.05, "t_max": 0.5, "t_step": 0.15}
    r = t1.post("/api/jobs", json={"kind": "rate", "method": "engine", "settings": settings, "x0": "text"})
    assert r.status_code == 200, r.text
    job = wait_job(t1, r.json()["id"])
    assert job["state"] == "finished", job
    # the engine's lines (pyazr's logger pyazr.engine) reach the job's log, as the program prints them
    assert any("Integrating reaction rate over 4 temperatures" in line for line in job["log"]), job["log"]
    res = job["result"]
    # the job's parameters are a data session's x0 (the start of every job), given to the rate by name; a rate
    # session's own transformation of the .azr's capture widths differs from it at 1e-10 (the program's rate takes
    # that one: tests/pyazr/reaction_rate_test.py)
    import pyazr
    with pyazr.Session(t1.path) as s:
        x0 = dict(zip([p.name for p in s.parameters.free], np.asarray(s.parameters.x0, float)))
    with rate.RateSession(t1.path, 1, 2) as rs:
        direct = rs.rate(x0, t_min=0.05, t_max=0.5, t_step=0.15)
    assert res["rate"]["rate"] == direct.rate.tolist() and res["rate"]["T9"] == direct.T9.tolist(), res["rate"]
    assert res["rate"]["units"]["rate"] == "cm^3 mol^-1 s^-1" and res["rate"]["reaction"] == "12C(p,g)13N"
    g = rate.gamow_window(t1.path, 1, direct.T9)
    assert res["gamow"]["E0"] == g["E0"].tolist()
    assert res["ratios"]["rate"] == [1.0] * direct.T9.size and "band" not in res
    assert os.path.isfile(res["csv"]) and "units: T9 in GK" in res["table"]
    formats = res["formats"]
    assert formats["species"]["reactants"] == ["p", "c12"] and formats["species"]["products"] == ["n13"]
    assert "no band" in formats["starlib"]["what"] and formats["starlib"]["fix"], formats["starlib"]
    assert "a REACLIB set stands for 0.01 to 10 GK" in formats["reaclib"]["what"], formats["reaclib"]
    # on the table grid the REACLIB file is the one pyazr.rate_formats writes from the same rate
    from pyazr import rate_formats
    r = t1.post("/api/jobs", json={"kind": "rate", "method": "engine", "x0": "text",
                                   "settings": {"entrance": 1, "exit": 2, "grid": "table", "reaclib_flag": "r"}})
    table_job = wait_job(t1, r.json()["id"])["result"]
    with rate.RateSession(t1.path, 1, 2) as rs:
        table = rs.rate(x0, grid="table")
    path = os.path.join(os.path.dirname(t1.path), "direct.reaclib")
    rate_formats.write_reaclib(path, rate_formats.fit_reaclib(table), rate_formats.species_of(t1.path, 1, 2),
                               flag="r")
    with open(path) as f:
        assert table_job["formats"]["reaclib"]["text"] == f.read(), table_job["formats"]["reaclib"]
    assert table_job["formats"]["reaclib"]["sets"] >= 1 and table_job["formats"]["reaclib"]["max_deviation"] <= 0.05
    for bad, words in (({"grid": "nope"}, "no temperature grid"), ({"band": "covariance"}, "no covariance")):
        r = t1.post("/api/jobs", json={"kind": "rate", "method": "engine", "settings": dict(settings, **bad)})
        failed = wait_job(t1, r.json()["id"])
        assert failed["state"] == "failed" and words in failed["result"]["message"], failed["result"]
    assert np.all(direct.rate > 0)


if __name__ == "__main__":
    uitest.run(__file__)
