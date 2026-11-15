"""audit_subagent.py (2026-11-07): the three idle-time review prompts and the FINDINGS.md parser.  Engine-free, seconds."""
import os, shutil, sys, tempfile, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit import AZURE2_ROOT
from rmfit import audit_subagent as asub
from rmfit.audit import audit_model
from rmfit.spec import FitSpec
from pyazr import AzrModel

SRC = os.path.join(AZURE2_ROOT, "tests", "13N")
FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""))
    if not ok:
        FAILS.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_asub_")
    try:
        shutil.copytree(os.path.join(SRC, "data"), os.path.join(tmp, "data"))
        os.makedirs(os.path.join(tmp, "output"))
        m = AzrModel.from_file(os.path.join(SRC, "13N.azr"))
        m.set_output_dir(os.path.join(tmp, "output"))
        azr = m.write(os.path.join(tmp, "13N.azr"))
        open(os.path.join(tmp, "readme"), "w").write("test campaign -- the evaluator's framing line\nsecond line\n")
        spec = {"name": "13N test", "data_sets": [{"name": "capture", "file": "13N", "role": "anchor", "norm": {"fixed": True, "nominal": 1.0}}],
                "standing_decisions": [{"date": "2026-11-07", "who": "test", "text": "norms fixed at 1 for the test"}]}
        json.dump(spec, open(os.path.join(tmp, "spec.json"), "w"))
        fs = FitSpec.load(os.path.join(tmp, "spec.json"))
        machine = audit_model(azr, campaign_dir=tmp, run_target_energies=False)
        for kind in asub.KINDS:
            r = asub.write_subagent_audit(kind, azr, tmp, fs, machine)
            txt = open(r["prompt"]).read()
            check(f"{kind}: prompt written", os.path.exists(r["prompt"]) and len(txt) > 2000, str(len(txt)))
            check(f"{kind}: no unfilled $placeholders", "$" not in txt.replace("$(", ""), [ln for ln in txt.splitlines() if "$" in ln][:2])
            check(f"{kind}: segment table present", "| key | active |" in txt and "| 1 | 1 |" in txt)
            check(f"{kind}: standing decision carried", "norms fixed at 1 for the test" in txt)
            check(f"{kind}: readme head carried", "the evaluator's framing line" in txt)
            check(f"{kind}: report path and rules", r["report"] in txt and "READ-ONLY" in txt)
            check(f"{kind}: skeleton FINDINGS.md", os.path.exists(r["report"]) and asub.parse_findings(r["report"]) == [])
        txt = open(os.path.join(tmp, "audits", "literature", "PROMPT.md")).read()
        check("literature: pairs with thresholds", "Compound nucleus 13N" in txt and "threshold at Ex" in txt)
        txt = open(os.path.join(tmp, "audits", "data", "PROMPT.md")).read()
        check("data: active files listed with row counts", "rows; segments" in txt)
        # a subagent writes its report: parse, summarize, mark
        rep = os.path.join(tmp, "audits", "data", "FINDINGS.md")
        open(rep, "w").write("""# data audit: test, 2026-11-07

one paragraph.

## Findings

| # | Sev. | Status | Where | What | Evidence | Why it matters | Cheapest check |
|---|------|--------|-------|------|----------|----------------|----------------|
| 1 | **HIGH** | candidate | `x.dat:4007` seg 3 | a 0.02 % error | row 4007: 2.92E-04 vs 3.6E-02 | pins the fit | re-evaluate seg 3 |
| 2 | MEDIUM | candidate | `y.dat` | last row repeated | lines 156-157 equal | double weight | drop one |
| 3 | LOW | `z.dat` | no status column here | n/a | n/a | n/a |

## Checked and clean
- integrity
""")
        rows = asub.parse_findings(rep)
        check("parse: three rows", len(rows) == 3, str(len(rows)))
        check("parse: bold severity stripped", rows[0]["severity"] == "HIGH" and rows[0]["status"] == "candidate")
        check("parse: missing status -> candidate, columns shift", rows[2]["status"] == "candidate" and rows[2]["where"] == "`z.dat`", str(rows[2]))
        check("set_status: existing cell", asub.set_status(rep, 1, "confirmed", "line checked") and asub.parse_findings(rep)[0]["status"] == "confirmed")
        check("set_status: inserted cell", asub.set_status(rep, 3, "rejected") and asub.parse_findings(rep)[2]["status"] == "rejected")
        check("set_status: unknown row", asub.set_status(rep, 9, "confirmed") is False)
        sm = asub.findings_summary(tmp)
        d = next(x for x in sm if x["kind"] == "data")
        check("summary: counts", d["n"] == 3 and d["by_status"]["candidate"] == 1 and d["by_severity"]["HIGH"] == 1, str(d["by_status"]))
        check("summary: unverified is the MEDIUM", [r["n"] for r in d["unverified"]] == [2])
        # the second write keeps a filled FINDINGS.md
        asub.write_subagent_audit("data", azr, tmp, fs, machine)
        check("rewrite keeps a filled report", len(asub.parse_findings(rep)) == 3)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("ALL OK" if not FAILS else f"FAILED: {len(FAILS)}: " + "; ".join(FAILS))
    return 0 if not FAILS else 1


if __name__ == "__main__":
    raise SystemExit(main())
