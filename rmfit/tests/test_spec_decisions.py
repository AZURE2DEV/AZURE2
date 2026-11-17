"""spec.py + decisions.py (2026-11-07): engine-free, seconds on the login node.  Uses the 13N test model as the .azr."""
import json, os, shutil, sys, tempfile, types
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_spec_decisions.py)")
    sys.exit(77)
from rmfit import AZURE2_ROOT
from rmfit.spec import FitSpec, Finding
from rmfit.decisions import DecisionQueue
from rmfit.campaign import Campaign

AZR = os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr")


def codes(fs):
    return sorted({f.code for f in fs})


def test_spec(tmp):
    from pyazr.datasets import SegmentSet
    segs = SegmentSet.from_file(AZR)
    s0 = segs[0]
    name = os.path.basename(s0.data_file)
    spec = {
        "reaction": "12C+p", "parametrization": "park", "nshards": 3,
        "engine": {"flags": {"use_long_wavelength": False, "ignore_externals": True}},
        "model": {"pairs_required": [{"key": 1, "label": "p0"}, {"key": 99, "label": "nonexistent"}]},
        "data_sets": [
            {"file": name, "role": "anchor", "norm": {"nominal": 1.0, "fixed": True}},          # 13N: norms fixed -> ok
            {"file": name, "role": "excluded"},                                                 # but it is active -> HIGH
            {"file": "no_such_file", "role": "free"},                                           # missing -> MEDIUM
            {"file": name, "role": "free", "norm": {"fixed": False}, "window": [0.0, 0.5]},     # fixed in .azr; window
        ],
        "bounds": {"bg_energy": 15.0, "e_window": 0.5, "not_a_field": 1},
        "search": {"top_k": 3},
        "pre_authorized": {"adopt_level": {"max_gamma_kev": 50, "forbid_near_same_jpi": True}},
        "standing_decisions": [{"date": "2026-11-07", "who": "test", "text": "x"}],
        "bogus_section": 1,
    }
    sp = os.path.join(tmp, "spec.json")
    json.dump(spec, open(sp, "w"))
    fs = FitSpec.load(sp)
    assert fs.sha and len(fs.sha) == 64
    v = codes(fs.validate())
    assert "spec-unknown-key" in v and "spec-bounds-keys" in v, v
    d = fs.derive_campaign_cfg()
    assert d["flags"]["use_park"] is True and d["flags"]["use_brune"] is False and d["flags"]["ignore_externals"] is True, d
    assert d["policy"] == {"bg_energy": 15.0, "e_window": 0.5} and d["search"] == {"top_k": 3} and d["nshards"] == 3, d
    f = fs.check(AZR)
    c = codes(f)
    for want in ("spec-excluded-active", "spec-dataset-missing", "spec-pair-missing", "spec-norm-free"):
        assert want in c, (want, c)
    assert not any(x.code == "spec-norm-fixed" for x in f), "anchor with fixed norm must pass on the 13N model"
    assert any(x.code == "spec-window" and x.severity == "MEDIUM" for x in f), c
    # against a campaign.json with a different flag and pin
    cfg = {"flags": {"use_park": False}, "policy": {"bg_energy": 11.5}, "engine": {"so_md5": "0" * 32}, "spec_sha": "x"}
    fs.data["engine"]["so_md5"] = "f" * 32
    c2 = codes(fs.check(AZR, campaign_cfg=cfg))
    for want in ("spec-flags", "spec-bounds", "spec-engine", "spec-changed"):
        assert want in c2, (want, c2)
    md = fs.render_md()
    assert "Standing decisions" in md and name in md and "pre_authorized" in md
    print("spec: ok")
    return fs


def test_adopt_gate(fs):
    stub = types.SimpleNamespace(spec=fs)
    allows = Campaign._spec_allows_adopt
    assert allows(stub, "1/2+", None, 10e3, 100.0, 1e5, [], None) is True          # 10 keV < 50, not near
    assert allows(stub, "1/2+", None, 80e3, 100.0, 1e5, [], None) is False         # too broad
    assert allows(stub, "1/2+", None, 10e3, 100.0, 1e5, ["1/2+@3.5"], None) is False  # near same J^pi
    fs.data["pre_authorized"]["adopt_level"] = False
    assert allows(stub, "1/2+", None, 1e3, 1e4, 1e5, [], None) is False
    fs.data["pre_authorized"]["adopt_level"] = True
    assert allows(stub, "1/2+", None, 1e3, 1e4, 1e5, [], None) is True
    assert allows(types.SimpleNamespace(spec=None), "1/2+", None, 1e9, 0, 1, ["x"], 1.0) is True   # no spec: old behaviour
    print("adopt gate: ok")


def test_decisions(tmp):
    q = DecisionQueue(os.path.join(tmp, "DECISIONS.md"))
    assert q.pending() == []
    d1 = q.add("adopt level 6+ 15.23", "adopt the 6+?", evidence="fig.png", options=("adopt", "reject"), default="wait")
    d2 = q.add("drop-row", "drop the spiked row?", options=("yes", "no"), default="no", prose="Row 17 of Kunz is 5 sigma off.")
    assert d1 == "D-001" and d2 == "D-002"
    p = q.pending()
    assert [d.id for d in p] == ["D-001", "D-002"] and p[0]["options"] == "adopt | reject" and p[1]["prose"].startswith("Row 17")
    q.answer("D-001", "adopt", by="the evaluator", note="J undecided, test 2+ too")
    assert [d.id for d in q.pending()] == ["D-002"]
    got = q.get("D-001")
    assert got["status"] == "answered" and got["answer"] == "adopt" and got["answered_by"] == "the evaluator"
    done = q.apply_defaults()
    assert [d.id for d in done] == ["D-002"] and q.get("D-002")["answer"] == "no" and q.pending() == []
    text = open(q.path).read()
    assert "## D-001 [answered]" in text and "## D-002 [answered]" in text
    # re-read round trip keeps every field
    q2 = DecisionQueue(q.path)
    assert q2.get("D-002")["note"] == "default_if_unanswered applied"
    print("decisions: ok")


if __name__ == "__main__":
    tmp = tempfile.mkdtemp(prefix="rmfit_spec_")
    try:
        fs = test_spec(tmp)
        test_adopt_gate(fs)
        test_decisions(tmp)
        print("test_spec_decisions: PASS")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
