"""Writes the p+12C elastic models of this check (run once; the .azr files are committed).

All share the levels of tests/hard_sphere/models/no_dummies.azr (1/2+ 2.365, 3/2- 3.502, 5/2+ 3.547 MeV in 13N,
includeAllJPi) and one test segment: elastic differential cross section at 120 deg, 1.0-3.0 MeV lab in 1 keV steps.
They differ only in the <targetInt> line:

  bare             no experimental effect
  fixed            Gaussian beam spread, the constant-sigma field: 30 keV (lab)
  legacy_cm_const  the expression form without the sigmalab keyword, i.e. the pre-2026-10 convention: expression of
                   the c.m. energy, value in c.m. MeV; a constant 30 keV * F (F = 12/(12+1.00782)) -- the same beam
                   spread as `fixed`
  legacy_cm_E      the same convention with sigma_cm(E_cm) = a0 * E_cm^1.5, a0 = 0.01
  lab_const        the expression form with the sigmalab keyword: sigma_lab(E_lab) = 0.030, must equal `fixed`
  lab_E            sigmalab with sigma_lab(E_lab) = a0 * sqrt(F) * E_lab^1.5, the lab form of legacy_cm_E
  tgt              target integration alone: 5e17 atoms/cm2, stopping cross section 4e-20 MeV cm2 (~20 keV thick)
  tgt_fixed / tgt_legacy_const / tgt_lab_const   target integration plus the beam spread, three ways
  tgt_straggle     target integration with straggling (0.04 keV^1/2 per keV^1/2), no beam spread
"""
import os, re
F = 12.0 / (12.0 + 1.00782)
TAIL = "0 0.04 20 50"      # straggling off, its coefficient, resonance-width multiplier, points per width
SEG = '1 "1" 50'
LINES = {
    "bare": None,
    "fixed":           f'{SEG} 1 0.030 0 0 "" 0 0 0 0 "" 0 {TAIL}',
    "legacy_cm_const": f'{SEG} 0 0 0 0 "" 0 0 0 1 "a0" 1 {0.030 * F:.10g} {TAIL}',
    "legacy_cm_E":     f'{SEG} 0 0 0 0 "" 0 0 0 1 "a0*x^1.5" 1 0.01 {TAIL}',
    "lab_const":       f'{SEG} 0 0 0 0 "" 0 0 0 1 "a0" 1 0.030 {TAIL} sigmalab',
    "lab_E":           f'{SEG} 0 0 0 0 "" 0 0 0 1 "a0*x^1.5" 1 {0.01 * F ** 0.5:.12g} {TAIL} sigmalab',
    "tgt":             f'{SEG} 0 0 1 5e17 "a0" 1 4e-20 0 0 0 "" 0 {TAIL}',
    "tgt_fixed":       f'{SEG} 1 0.030 1 5e17 "a0" 1 4e-20 0 0 0 "" 0 {TAIL}',
    "tgt_legacy_const": f'{SEG} 0 0 1 5e17 "a0" 1 4e-20 0 0 1 "a0" 1 {0.030 * F:.10g} {TAIL}',
    "tgt_lab_const":   f'{SEG} 0 0 1 5e17 "a0" 1 4e-20 0 0 1 "a0" 1 0.030 {TAIL} sigmalab',
    "tgt_straggle":    f'{SEG} 0 0 1 5e17 "a0" 1 4e-20 0 0 0 "" 0 1 0.04 20 50',
}
base = open(os.path.join(os.path.dirname(__file__), "..", "hard_sphere", "models", "no_dummies.azr")).read()
base = re.sub(r"<segmentsTest>.*?</segmentsTest>", "<segmentsTest>\n1 1 1 1.0 3.0 0.002 120 120 1 1\n</segmentsTest>", base, flags=re.S)
for name, line in LINES.items():
    text = re.sub(r"<targetInt>\n</targetInt>", "<targetInt>\n" + (line + "\n" if line else "") + "</targetInt>", base)
    assert ("<targetInt>\n" + (line or "")) in text, name
    with open(os.path.join("models", name + ".azr"), "w") as f: f.write(text)
print("written", len(LINES), "models")
