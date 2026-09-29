# thm_distortion — the distortion factor R(E) of a THM experiment

`check.sh` runs copies of `tests/18O_p_a_thm` (two THM segments with free
norms; the real 2H(18O,α15N)n at 54 MeV has a neutron spectator, a charged one
is made up as 18O(3He,α15N)d at 115 MeV) and a one-segment, unfolded copy of
`examples/c12c12_tumino2018` (12C(14N,α/p)d at 30 MeV) with
`experiment[<name>] ... distortion=...` in the `<thm>` block
(docs/source/theory/thm_implementation.rst, "Distortion factor R(E)").

| case | expectation |
|---|---|
| (a) `distortion=none` | output files byte-identical to the same line without the key |
| (b) neutral spectator, `distortion=optical opticalAA=plane` | nothing is distorted: model ratio 1 at all 118 points (measured 5e-11); with point Coulomb in d + 18O (η = 0.73) R = 0.93–1.07 |
| (c) lab 30°, no folding | model(on)/model(off) at the lowest and highest point equals the R that `thm_experiments.out` evaluates there directly (1e-5); `distortionRef=` a point energy gives R = 1 there (1e-6) |
| (d) 12C(14N,d), points at 0.8–2.55 MeV, E_ref = 2.664 MeV | 1/R against Mukhamedzhanov & Pang PRC 99 (2019) Fig. 10 (defaults: rms 0.031, max 0.049 dex) and Mukhamedzhanov arXiv:2609.04498 Fig. 9 (`distortionRatio=dw`, `Ebeam=30.11`: rms 0.010, max 0.023 dex); limits rms 0.035, max 0.06 dex (digitisation 0.03 dex) |
| (e) `distortion=optical` with ten zeros | chiSquared.out, normalizations.out and the model byte-identical to `distortion=coulomb`; a Woods-Saxon potential changes the model |
| (f) refusals and warnings | no kinematics, bad values of every key, optical keys without `optical`, detail keys without a computed distortion, a malformed `theta`, a key twice, a missing or too short table, E_sF ≤ 0, an unreachable lab angle, E_ref beyond the spectator threshold; the B(x+s) warning when the masses (5.493 MeV for 3He) and field 32 (2.2246) disagree, none when they agree |

The amplitude itself is checked against an independent mpmath/scipy
evaluation in `tests/reference` (ctest `thm_distortion`), the Python side in
`tests/pyazr/thm_distortion_test.py`.
