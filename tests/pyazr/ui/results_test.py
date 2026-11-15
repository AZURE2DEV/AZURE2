"""The Results payload of the web UI (GET /api/evaluation, POST and GET /api/outputs) against the program's own output
files (docs/upgrade/web_ui.md section 6.5).

- T2 (thin and thick targets, a data covariance joining two data sets, a sum, a ratio, a sample, THM) and T1 (angular
  distributions of elastic scattering and capture, analyzing powers, phase shifts, five extrapolations): every data
  set's points (c.m. energy, excitation energy, c.m. angle), its model, its datum and error times the norm, and the S
  factor where it is offered equal the columns of the program's AZUREOut_* (menu 1) to their eleven digits; every
  extrapolation equals the program's AZUREOut_*.extrap (menu 3); the chi-squared table is the program's
  chiSquared.out (rows and Total- line); the squares of the residuals add up to each data set's chi-squared, and
  for the two data sets of T2's data covariance to their sum.
- T3: the model data equal pyazr's model_data; no S factor for a neutral entrance.
- The program's output files written at the working x into the project's output directory and listed.

Run: python3 tests/pyazr/ui/results_test.py   (AZURE2_BIN: the program; ctest gives it)
"""
import os
import subprocess
import sys

try:
    import pytest
except ImportError:  # ctest counts exit 77 as skipped
    print("skip: the package pytest is not installed (pip install pyazr[web])")
    sys.exit(77)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client  # noqa: E402

np = pytest.importorskip("numpy")

REL = 1e-9  # AZUREOut_* columns have eleven significant digits (precision 10, scientific)


def run_program(azr, command):
    """The program on ``azr`` by ``command`` (calculate: with data, extrapolate: extrapolations); its output directory."""
    proc = subprocess.run([os.environ["AZURE2_BIN"], command, os.path.basename(azr)],
                          stdin=subprocess.DEVNULL, text=True, cwd=os.path.dirname(azr), capture_output=True,
                          timeout=600)
    out = os.path.join(os.path.dirname(azr), "output")
    assert any(n.startswith("AZUREOut_") for n in os.listdir(out)), proc.stdout[-2000:]
    return out


def blocks(directory, suffix):
    """Every block of numbers (one per segment, blank lines between) of the AZUREOut_* files ending in ``suffix``."""
    out = []
    for name in sorted(os.listdir(directory)):
        if not (name.startswith("AZUREOut_") and name.endswith(suffix)) or "sample_seg" in name:
            continue
        rows = []
        with open(os.path.join(directory, name)) as f:
            for line in list(f) + [""]:
                cells = line.split()
                if not cells:
                    if rows:
                        out.append((name, np.array(rows)))
                    rows = []
                    continue
                rows.append([float(c) for c in cells])
    return out


def close(a, b, rel=REL):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return a.shape == b.shape and bool(np.all(np.abs(a - b) <= rel * np.abs(b) + 1e-300))


def match(entry, candidates, with_data):
    """The block of the program's files that holds ``entry``'s points and agrees in every column it shares."""
    for name, b in candidates:
        if b.shape[0] != len(entry["e_cm"]) or b.shape[1] < 5:
            continue
        cols = [(entry["e_cm"], b[:, 0]), (entry["e_x"], b[:, 1]), (entry["angle_cm"], b[:, 2]),
                (entry["model"], b[:, 3])]
        if "s_factor" in entry:
            cols.append((entry["s_factor"]["model"], b[:, 4]))
        if with_data:
            cols += [(entry["value_norm"], b[:, 5]), (entry["error_norm"], b[:, 6])]
            if "s_factor" in entry:
                cols += [(entry["s_factor"]["value"], b[:, 7]), (entry["s_factor"]["error"], b[:, 8])]
        if all(close(a, c) for a, c in cols):
            return name
    return None


def program_chi_squared(directory):
    from pyazr.ui.engine import _read_chi_squared
    return _read_chi_squared(os.path.join(directory, "chiSquared.out"))


@pytest.mark.parametrize("name", ["T1", "T2"])
def test_payload_is_the_programs_output(name, tmp_path):
    program = run_program(copy_project(name, tmp_path / "cli"), "calculate")
    rows, totals = program_chi_squared(program)
    with make_client(copy_project(name, tmp_path / "api")) as client:
        r = client.get("/api/evaluation")
        assert r.status_code == 200, r.text
        e = r.json()
    # the chi-squared table: the program's rows and Total- line
    assert [(d["chi2"], d["n"], d["norm"], d["norm_chi2"]) for d in e["datasets"]] == \
        [(r["chi2"], r["n"], r["norm"], r["norm_chi2"]) for r in rows]
    assert e["totals"] == totals
    # every data set's arrays: the program's AZUREOut_* columns
    found = blocks(program, ".out")
    for d in e["datasets"]:
        assert match(d, found, with_data=True), f"{d['name']}: no block of the program's files agrees"
    # the residuals: their squares are the data set's chi-squared, except where a data covariance joins sets
    joined = {"s02_diff090", "s03_diff150"} if name == "T2" else set()
    for d in e["datasets"]:
        r2 = float(np.sum(np.asarray(d["residual"], float) ** 2))
        assert abs(r2 - d["residual_chi2"]) <= 1e-12 * max(r2, 1.0), d["name"]
        if d["name"] not in joined:
            assert abs(d["residual_chi2"] - d["chi2"]) <= 5e-6 * d["chi2"], (d["name"], d["residual_chi2"], d["chi2"])
    if joined:
        group = [d for d in e["datasets"] if d["name"] in joined]
        assert abs(sum(d["residual_chi2"] for d in group) - sum(d["chi2"] for d in group)) <= 1e-5 * 513.9
        assert not all(abs(d["residual_chi2"] - d["chi2"]) <= 1e-4 * d["chi2"] for d in group)  # two splits
    assert "residual" in e and "n*datum" in e["residual"]["definition"]
    if name == "T1":  # the extrapolations: the program's extrapolate
        extrap = blocks(run_program(copy_project(name, tmp_path / "cli3"), "extrapolate"), ".extrap")
        assert len(e["extrapolations"]) == 5
        for x in e["extrapolations"]:
            assert match(x, extrap, with_data=False), f"{x['name']}: no block of the program's .extrap files agrees"
        assert any("s_factor" in d for d in e["datasets"]) and any("s_factor" in x for x in e["extrapolations"])


def test_model_data_and_neutral_entrance(tmp_path):
    import dataclasses

    import pyazr
    azr = copy_project("T3", tmp_path)
    with make_client(azr) as client:
        e = client.get("/api/evaluation").json()
        reactions = {r["name"]: r["fields"].get("reaction", "") for r in client.get("/api/project/records").json()}
    session = pyazr.Session(azr)  # at the text's x, as the server's working x before any edit
    want = [dataclasses.asdict(m) for m in pyazr.model_data(session)]
    session.close()
    assert e["model_data"] and e["model_data"] == want
    # the S factor where the engine unit is a cross section and the entrance charged: 7Li + p, not 7Be + n
    offered = sorted(reactions[d["name"]] for d in e["datasets"] + e["extrapolations"] if "s_factor" in d)
    assert offered and all(r.startswith("7Li(p,") for r in offered), offered
    assert any(r.startswith("7Be(n,") for r in reactions.values())


def test_outputs_written_into_the_projects_output_directory(tmp_path):
    azr = copy_project("T2", tmp_path / "api")
    program = run_program(copy_project("T2", tmp_path / "cli"), "calculate")
    with make_client(azr) as client:
        assert client.get("/api/outputs").json()["files"] == []
        r = client.post("/api/outputs", json={})
        assert r.status_code == 200, r.text
        written = r.json()
        assert written["directory"] == "output" and "chiSquared.out" in written["files"]
        assert any(n.startswith("AZUREOut_") for n in written["files"])
        listed = [f["name"] for f in client.get("/api/outputs").json()["files"]]
        assert set(written["files"]) <= set(listed)
        head = client.get("/api/files/output/chiSquared.out", params={"rows": 3}).json()
        assert head["lines"][0].startswith("Segment#")
    with open(os.path.join(os.path.dirname(azr), "output", "chiSquared.out")) as a, \
            open(os.path.join(program, "chiSquared.out")) as b:
        assert a.read() == b.read()


if __name__ == "__main__":
    uitest.run(__file__)
