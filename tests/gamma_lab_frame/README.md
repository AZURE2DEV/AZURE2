# gamma_lab_frame — lab-to-c.m. conversion of capture angular distributions

Capture data given in the lab frame (`isDiff 1` with a photon exit pair) are converted to the c.m. frame by
`Unit::ConvertLabAngleGammas` and `Unit::ConvertCrossSectionGammas`, using the photon aberration with the c.m.
velocity β:

- cos θ* = (cos θ − β)/(1 − β cos θ)
- dσ/dΩ* = dσ/dΩ_lab · (1 − β cos θ)²/(1 − β²)

The second factor is the solid-angle Jacobian d cos θ/d cos θ*.

`make_models.py` (deleted; see "How the models were made") writes a ¹²C(α,γ)¹⁶O(g.s.) external-capture model (E1 + E2, strongly asymmetric). `check.sh` does
three things:

1. It computes the c.m. differential cross section at E_lab = 2, 3, 4.5, 6 MeV and θ* = 10°–170° (`cm.azr`).
2. It boosts those values to the lab with `awk`, using the photon's four-vector transformation. No AZURE2 code is
   used for this step.
3. It reads the boosted values back as lab-frame data (`lab.azr`).

It checks that:

- (a) χ² = 0: every point returns to its θ* and σ*;
- (b) the lab-frame output file (`kinflag.dat` = 1) shows the model equal to the data, i.e. the fit is taken back
  to the lab with the same factor.

**What it caught (fixed 2026-10-23).** The cross-section factor was (1 − β)²/(1 − β cos θ). Relative to the correct
Jacobian this is ≈ 1 − 2β + 3β cos θ: an overall −2β plus a P1-shaped distortion. For ¹²C(α,γ) at E_lab = 3 MeV
(β = 0.010) it is +0.9 % at 13°, −1.9 % at 89° and −4.8 % at 167°. Here (a) gave χ² = 31,722 over 28 points.
The kinematic factors for the lab-frame output were also not set, so those files printed the c.m. model next to the
lab data; that is (b). The angle conversion itself was right.

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show f389d46d:tests/gamma_lab_frame/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Models for tests/gamma_lab_frame: 12C(a,g)16O(g.s.) external capture only (E1 from l = 1 and E2 from l = 2, ANC
100 fm^-1/2, a = 5.43 fm, one 30 MeV level with a 1e-6 eV alpha width in 1- and in 2+), whose E1-E2 interference
makes the angular distribution strongly asymmetric.

  cm.azr  : one c.m. differential segment (isDiff 4) at E_LAB x ANGLES (data/cm.dat, dummy cross sections)
  lab.azr : the same reaction as a lab differential segment (isDiff 1) reading data/lab.dat, which check.sh writes
            from cm's output with an explicit Lorentz boost of the photon

Run from this directory: python3 make_models.py   (writes models/)
```
