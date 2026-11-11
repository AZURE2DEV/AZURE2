# tests/convolution_unify

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show dad58dcd:tests/convolution_unify/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Writes the p+12C elastic models of this check (run once; the .azr files are committed).

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
```
