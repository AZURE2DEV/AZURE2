# 7Li_p_ay — analyzing power on a target with spin

The vector analyzing power of 7Li(p,p)7Li: a **spin-3/2 target**, so the
channel spins are 1 and 2 and the projectile's spin projection only exists
after the channel-spin index is decomposed,

    M_{out; m1 m2} = sum_s <j1 m1 j2 m2 | s, m1+m2> A_{out; s, m1+m2}.

That decomposition is the general-target-spin path of
`AmplitudeMatrix::AnalyzingPowerAy`. Every other analyzing-power test is
12C+p, where the target spin is zero, the decomposition collapses, and the
path this project exists for never runs.

## The model and the data

The model is the 8Be evaluation of Paneru et al., PRC 111 (2025) 064609
(Table IV parameters; the same model as the validated `7Li_p_paper.azr`),
with two analyzing-power segments added:

| # | data | theta_lab | theta_cm | points |
|---|---|---|---|---|
| 1 | `deineko_130.dat` | 130 deg | 136.32 | 35 |
| 2 | `deineko_160.dat` | 160 deg | 162.82 | 37 |

Data: A. S. Deineko et al. (1974), EXFOR entry **F0285004**, E_lab =
1.83–3.02 MeV. EXFOR gives lab angles and no experimental uncertainties
(the DATA-ERR field is empty; only a digitizing error is recorded), so the
files carry the lab->cm converted angle and an assigned **absolute**
uncertainty of 0.03 — a typical polarimeter accuracy of the period. Ay
uncertainties must be absolute, never relative: a point at a zero crossing
is not infinitely precise.

## What the recorded chi-squared means

The Paneru parameters were fit **without any polarization data**, and Ay is
pure relative phase, so the recorded chi2 = 3955 over 72 points is not a
fit — it pins the calculation against drift, exactly as `13N_capture_ay`
does. The calculation does place structure at the right energies (Ex = 18.9
and 19.25 MeV, the 2- and 3+ states) and stays inside |Ay| <= 1.

The interesting physics: Table IV's width *signs* barely move the cross
section (per-level overall signs are gauge) but move Ay directly, since it
is built from channel-spin off-diagonal interference. This dataset is the
kind of measurement that would discriminate them.
