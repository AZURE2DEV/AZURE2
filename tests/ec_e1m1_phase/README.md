# ec_e1m1_phase

E1-M1 interference of external capture in a differential cross section, against an independent
m-scheme calculation (audit 2026-10-25, 04_external_capture_photons S5, checklist item 26).

The E external-capture amplitude carries `-i i^(li+L-lf)` and the M1 one `+i` (`ECIntegral.cpp`).
Angle-integrated cross sections never see the relative phase, so the E-only and M-only checks of
`tests/ec_m1_exact` cannot tell whether E1-M1 interference in an angular distribution has the right
sign. `truth_dsdo.py` computes dsigma/dOmega of p + 16O capture from the dipole interaction
`-d.E - mu.B` (the plane wave's E and B carry the same factor, so the amplitude is proportional to
`d.eps* + mu.(khat x eps*)`) between the Coulomb + hard-sphere scattering state, beam along z, each
partial wave with its Coulomb and hard-sphere phase, and the ANC-normalized bound level
`|(lf 1/2) Jf Mf>`, summed over magnetic substates and photon helicities. No multipole phase
convention enters. Its golden-rule constant is fixed by the E1-only angle integral against the
`(l s)J` reduced-matrix-element cross section of `tests/ec_m1_exact/truth.py`, and the M1-only
integral is checked against the same formula (orbital + spin) to 1e-10, so the interference term
carries no free factor.

Two final states (`make_models.py` (deleted; see "How the models were made"), models committed):

- **s**: 17F(0.4953, 1/2+), l_f = 0, ANC 84: E1 from l = 1, M1 (spin term) from l = 0. With a
  spin-0 target and an unpolarized beam these cannot interfere: the E1 amplitude is spin-diagonal
  and the spin operator is traceless, so dsigma/dOmega(E1 + M1) = E1 + M1 exactly. AZURE2 reproduces
  that, a check of its unpolarized-beam algebra.
- **p**: a fictitious 1/2- state at E_x = 0.4996 MeV, l_f = 1, ANC 20: E1 from l = 0 and 2, M1
  (orbital + spin) from l = 1. The orbital M1 amplitude is spin-diagonal and interferes with E1 in
  the odd Legendre orders: a forward/backward asymmetry of 16 % at 1 MeV. (An ANC of 84 here
  would put 130 % of the state outside 5 fm; AZURE2 then warns "Denominator less than zero while
  transforming" and its cross section stops scaling as ANC squared.)

## Checks (`check.sh path/to/AZURE2`)

For each case, with the photon pair's multipolarity mask 1 (E1), 2 (M1), 3 (both) and one c.m.
differential capture segment at E_cm = 0.3, 1.0, 2.0 MeV and seven angles:

- (a) E1 only and (b) M1 only equal the truth at every point (2e-5 of the row maximum): the
  differential capture machinery and units;
- (c) E1 + M1 equals the truth at every point: the relative phase of the E and M amplitudes.

`expected/dsdo_{s,p}.dat` come from `truth_dsdo.py` (mpmath, numpy, scipy, sympy; a few minutes).

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show 3983cd76:tests/ec_e1m1_phase/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Models for ec_e1m1_phase (committed): p + 16O external capture with E1 and M1 into

  s : 17F(0.4953, 1/2+), the 2s_1/2 state of tests/ec_m1_exact (ANC 84 fm^-1/2, bound by 104 keV): E1 from l = 1
      (J = 1/2-, 3/2-), M1 from l = 0 (J = 1/2+, the spin term alone).  For a spin-0 target and an unpolarized beam the
      two cannot interfere in dsigma/dOmega: the E1 amplitude is spin-diagonal and the M1 spin operator is traceless.
  p : a fictitious 1/2- state at E_x = 0.4996 MeV (l_f = 1, bound by 100 keV, ANC 20): E1 from l = 0 (J = 1/2+) and
      l = 2 (J = 3/2+), M1 from l = 1 (J = 1/2-, 3/2-) with the orbital term g_l L and the spin term g_p S_p; the
      orbital M1 amplitude is spin-diagonal and interferes with E1 (odd Legendre orders).

Each with the photon pair's multipolarity mask 1 (E1), 2 (M1) and 3 (both), and one c.m. differential capture
segment (isDiff 4) at the energies and angles of truth_dsdo.py.  Entrance J-groups hold one 30 MeV level with a
1e-6 eV proton width (negligible internal amplitude), as in tests/ec_m1_exact.
```
