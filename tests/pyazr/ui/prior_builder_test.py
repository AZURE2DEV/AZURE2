#!/usr/bin/env python3
"""The prior builder through its routes, against the pyazr calls made directly (T1):

  GET /api/prior opens with the key unset, no parts, the proposed file, pyazr.priors.kinds_table() and the project's
  files; POST /api/prior/build answers pyazr.priors.build's summary on the session's parameter table (names, centres,
  covariance, checks, notes, the file's text); a name the table lacks is refused with what / allowed / fix; a matrix
  that is not positive definite comes back refused, naming its parts, and is not written; POST /api/prior/write
  writes the file pyazr writes and sets priorCovariance (the text equals Project.set_config's), the prior's term in
  the server's evaluation equals a session's on that text; the builder reopens with the parts; one Undo removes the
  file and restores the text byte for byte, Redo writes the same bytes; a second write under the same name with
  other parts is one history entry whose Undo puts the first file back; "Remove the prior" (the config route) takes
  the key out and Undo puts it back.

Run: python3 tests/pyazr/ui/prior_builder_test.py
"""
import os
import sys

try:
    import pytest
except ImportError:  # ctest counts exit 77 as skipped
    print("skip: the package pytest is not installed (pip install pyazr[web])")
    sys.exit(77)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client  # noqa: E402


@pytest.fixture()
def t1(tmp_path):
    path = copy_project("T1", tmp_path)
    with make_client(path) as client:
        client.path = path
        yield client


def rev(client):
    return client.get("/api/project").json()["rev"]


def test_prior_builder(t1):
    import numpy as np
    import pyazr
    from pyazr import priors
    directory = os.path.dirname(t1.path)
    opened = t1.get("/api/prior").json()
    assert opened["named"] is None and opened["parts"] == [] and opened["proposed"] == "prior_parameters.txt"
    assert opened["kinds"] == priors.kinds_table() and "T1.azr" in opened["files"], opened

    with pyazr.Session(t1.path) as s:
        table = s.parameters
        free = table.free
        norms = [p.name for p in free if p.kind == "energy"][:2]   # two free level energies (T1 fits no norm)
        widths = [p.name for p in free if p.kind == "width" and p.start != 0][:2]
        assert len(norms) == 2 and len(widths) == 2, (norms, widths)
        fixed = next(p.name for p in table if p.fixed and p.kind == "energy")
        parts = [{"kind": "width", "names": norms, "percent": 5},
                 {"kind": "width", "names": [fixed], "sigma": 0.02},
                 {"kind": "width", "names": widths, "percent": 10},
                 {"kind": "correlation", "pairs": [[widths[0], norms[0], "0.5"], [widths[1], widths[0], -0.3]]}]
        direct = priors.build(table, parts, directory=directory)
    t1.get("/api/parameters")   # the page reads the table first: it must not change the rows the build reads
    r = t1.post("/api/prior/build", json={"parts": parts})
    assert r.status_code == 200, r.text
    built = r.json()
    want = direct.summary()
    for key in ("names", "mean", "sigma", "covariance", "correlation", "checks", "notes", "text", "basis", "fixed"):
        assert built[key] == want[key], (key, built[key], want[key])
    assert built["refusal"] is None

    bad = t1.post("/api/prior/build", json={"parts": [{"kind": "width", "names": ["no_such"], "sigma": 1}]})
    body = bad.json()
    assert bad.status_code == 404 and body["what"] and body["allowed"] and body["fix"], bad.text
    singular = [{"kind": "width", "names": norms + widths[:1], "sigma": 0.1},
                {"kind": "correlation", "pairs": [[norms[0], norms[1], 0.9], [norms[0], widths[0], 0.9],
                                                  [norms[1], widths[0], -0.9]]}]
    refused = t1.post("/api/prior/build", json={"parts": singular}).json()
    assert refused["refusal"] and refused["refusal_parts"] == [0, 1] and refused["text"] is None, refused
    before = t1.get("/api/project/text").text
    w = t1.post("/api/prior/write", json={"parts": singular, "path": "p_x.txt", "rev": rev(t1)})
    assert w.status_code == 422 and not os.path.exists(os.path.join(directory, "p_x.txt")), w.text
    assert t1.get("/api/project/text").text == before

    w = t1.post("/api/prior/write", json={"parts": parts, "path": "prior_parameters.txt", "rev": rev(t1)})
    assert w.status_code == 200, w.text
    file = os.path.join(directory, "prior_parameters.txt")
    with open(file) as f:
        written = f.read()
    assert written == direct.text()
    expect = pyazr.Project.open(t1.path)
    expect.set_config(priorCovariance="prior_parameters.txt")
    after = t1.get("/api/project/text").text
    assert after == expect.text
    again = t1.get("/api/prior").json()
    assert again["named"] == "prior_parameters.txt" and not again["hand_written"]
    assert again["parts"] == direct.parts, again["parts"]
    moved = {norms[0]: direct.mean[0] + 0.01, widths[0]: direct.mean[3] * 1.05}
    evaluation = t1.post("/api/evaluate", json={"x": moved}).json()
    saved = os.path.join(directory, "with_prior.azr")
    expect.save(saved)
    with pyazr.Session(saved) as s:
        x = s.parameters.put(s.parameters.x0, list(moved), list(moved.values()))
        assert evaluation["objective"] == s.objective()(x)
        term = s.evaluate(x).terms["prior covariance"]
    kept = [i for i, f in enumerate(direct.fixed) if not f]   # the engine marginalizes the fixed parameter
    assert len(kept) == len(direct.names) - 1 and np.all(np.isfinite(direct.mean)), direct.mean
    d = np.zeros(len(direct.names))
    d[0], d[3] = 0.01, direct.mean[3] * 0.05
    d, C = d[kept], direct.cov[np.ix_(kept, kept)]
    by_hand = float(d @ np.linalg.solve(C, d))
    assert evaluation["terms"]["prior covariance"] == term and abs(term - by_hand) <= 1e-10 * by_hand, (term, by_hand)

    # a second write under the same name: one entry, whose undo puts the first file back
    other = [{"kind": "width", "names": norms, "sigma": 0.2}]
    entries = len(t1.get("/api/project/history").json())
    w = t1.post("/api/prior/write", json={"parts": other, "path": "prior_parameters.txt", "rev": rev(t1)})
    assert w.status_code == 200, w.text
    assert len(t1.get("/api/project/history").json()) == entries + 1
    assert priors.read_parts(file)["parts"][0]["sigma"] == 0.2
    assert t1.post("/api/project/undo", json={"rev": rev(t1)}).status_code == 200
    with open(file) as f:
        assert f.read() == written

    # undo the first write: the file removed, the text as before; redo: the same bytes
    assert t1.post("/api/project/undo", json={"rev": rev(t1)}).status_code == 200
    assert t1.get("/api/project/text").text == before and not os.path.exists(file)
    assert t1.post("/api/project/redo", json={"rev": rev(t1)}).status_code == 200
    with open(file) as f:
        assert f.read() == written
    assert t1.get("/api/project/text").text == after

    # "Remove the prior": the config route; undo puts the key back
    r = t1.patch("/api/project/config", json={"entries": {"priorCovariance": None}, "rev": rev(t1)})
    assert r.status_code == 200, r.text
    assert t1.get("/api/project/text").text == before and os.path.exists(file)
    assert t1.post("/api/project/undo", json={"rev": rev(t1)}).status_code == 200
    assert t1.get("/api/project/text").text == after


if __name__ == "__main__":
    uitest.run(__file__)
