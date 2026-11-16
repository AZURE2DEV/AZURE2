#!/usr/bin/env python3
"""pyazr.schema(): the table of project format 3 (docs/upgrade/project_format3.md section 6.1).

  1. pyazr.schema() is a dict equal to the committed tests/schema/schema_v3.json, schema_version "3.0";
  2. _azure2.schema_json() is the text of that file (without its final newline), so the binding and the file agree
     byte for byte; with AZURE2_BIN set, `AZURE2 --schema` prints exactly the file as well;
  3. the content an agent or a form needs is there: the record kinds of each section, every field with a type and a
     grammar, the conditional fields of the user stories (charge_state with a per-charge unit), units with factors,
     and every diagnostic code with its meaning.
Needs the compiled engine; skips (exit 77) without it.  Run from anywhere:  python3 tests/pyazr/schema_test.py
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
GOLDEN = os.path.join(ROOT, "tests", "schema", "schema_v3.json")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    import pyazr
    from pyazr import _azure2
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

with open(GOLDEN, encoding="utf-8") as f:
    golden_text = f.read()
s = pyazr.schema()
check("pyazr.schema() is a dict", isinstance(s, dict))
check("schema_version 3.0", s.get("schema_version") == "3.0")
check("equals tests/schema/schema_v3.json", s == json.loads(golden_text))
check("_azure2.schema_json() + newline is the committed file", _azure2.schema_json() + "\n" == golden_text)

binary = os.environ.get("AZURE2_BIN")
if binary:
    out = subprocess.run([binary, "--schema"], capture_output=True, check=False)
    check("AZURE2 --schema prints the committed file", out.returncode == 0 and out.stdout.decode("utf-8") == golden_text,
          f"exit {out.returncode}")
else:
    print("  (AZURE2_BIN not set: the CLI comparison is skipped)")

check("sections", s["sections"] == {"data": ["dataset", "normalization", "energy_scale"],
                                    "targets": ["target", "beam", "detector"],
                                    "extrapolations": ["extrapolation"], "modelData": ["datum"],
                                    "potentials": ["potential"], "gammaDecays": ["gamma_decay"], "channelRadii": ["channel_radius"],
                                    "isospin": ["isospin_reference", "isospin_shift", "isospin_level", "isospin_mix",
                                                "isospin_factor", "isospin_pair"]})
pot = s["records"]["potential"]["fields"]
check("potential: pair, coulomb_radius and central required; matching_radius and step with their defaults",
      pot["pair"].get("required") is True and pot["coulomb_radius"].get("required") is True and
      pot["central"].get("required") is True and pot["matching_radius"].get("default") == "15 fm" and
      pot["step"].get("default") == "0.0275 fm", str({k: (v.get("required"), v.get("default")) for k, v in pot.items()}))
for kind, rec in s["records"].items():
    for key, field in rec["fields"].items():
        if not ({"type", "grammar", "doc"} <= set(field)):
            check(f"{kind}.{key} has type, grammar and doc", False)
yield_fields = s["quantities"]["yield"]["fields"]
check("yield: charge_state required with a per-charge unit",
      yield_fields["charge_state"]["required"] and set(yield_fields["charge_state"]["when"]["in"]) == {"per_charge", "counts_per_charge"})
check("unit factors", s["units"]["mb"]["factor"] == 1e-3 and s["units"]["mb"]["engine_unit"] == "b"
      and s["quantities"]["resonance_strength"]["units"]["keV"] == 1000.0)
check("ratio_to_rutherford refuses identical particles",
      any(r["rule"] == "identical" for r in s["quantities"]["ratio_to_rutherford"]["refusals"]))
check("diagnostic codes with meanings", len(s["diagnostic_codes"]) >= 40
      and all(c["meaning"] for c in s["diagnostic_codes"].values()) and "E_FIELD_REQUIRED" in s["diagnostic_codes"])

if failures:
    print(f"FAILED: {len(failures)} check(s)")
    sys.exit(1)
print("schema: all checks passed")
