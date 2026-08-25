# 7Li_p_a — the THM validity test of 7Li(p,alpha), pinned

The half-off-shell excitation function of 7Li(p,alpha)4He extracted by the
Trojan Horse method from d(7Li,alpha alpha)n breakup:

A. Tumino et al., Eur. Phys. J. A 27 (S1) 243 (2006) — the THM *validity
test* of this reaction — EXFOR **O1653002**, 66 points,
E_cm = 0.076-6.87 MeV, arbitrary units, theta_cm = 50-70 deg.

Model: the validated Paneru et al. 8Be evaluation (PRC 111, 064609,
Table IV; the same base as `6Li_d` and `7Li_p_ay`), entrance pair 5
(7Li+p), exit the identical-boson alpha+alpha pair. The Trojan horse is
the deuteron, so the entrance-pair channel lines carry B = 2.2246 MeV in
the optional 33-field column — the complementary case to `6Li_d`, whose
carrier is 6Li at 1.4735 MeV. A 30 keV (sigma) Gaussian resolution is
folded through the engine's sub-point convolution on the THM segment.

EXFOR carries no experimental uncertainties for this entry (only a
digitizing error), so the data files hold an assigned 10% (0.01 absolute
floor); the c.m. energies were converted to the lab convention the data
files use. The segment normalization is the THM arbitrary scale, solved
analytically against this model and written into the file.

The recorded chi2 = 4952.9 over 66 points is a **pin, not a fit**: the
Paneru parameters come from on-shell data alone, and the half-off-shell
observable weights the same poles differently. The two-peak structure
(E_cm = 2.6 and 5.1 MeV) sits at the right energies; the misfit
concentrates below 0.3 MeV, where the HOES observable rises toward
threshold against falling data, and in the relative peak heights — the
behaviour a joint THM+direct refit corrects (see the analysis notes).
