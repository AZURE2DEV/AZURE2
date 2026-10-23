# gamma_lab_frame — lab-to-c.m. conversion of capture angular distributions

Capture data given in the lab frame (`isDiff 1` with a photon exit pair) are converted to the c.m. frame by
`EPoint::ConvertLabAngleGammas` and `EPoint::ConvertCrossSectionGammas`, using the photon aberration with the c.m.
velocity β:

- cos θ* = (cos θ − β)/(1 − β cos θ)
- dσ/dΩ* = dσ/dΩ_lab · (1 − β cos θ)²/(1 − β²)

The second factor is the solid-angle Jacobian d cos θ/d cos θ*.

`make_models.py` writes a ¹²C(α,γ)¹⁶O(g.s.) external-capture model (E1 + E2, strongly asymmetric). `check.sh` does
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
