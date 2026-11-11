# scattering_length — s-wave scattering lengths at each pair's threshold

The project files are format 3 (`docs/upgrade/format3_tree.md`); the old lines quoted on this page are those of the
originals in `tests/convert/fixtures/tests/scattering_length/`.

`check.sh` (run by `run_tests.sh`) tests `ScatteringLength` (src/ScatteringLength.cpp): the s-wave
scattering length a0 = -a [M/(x^2 K1^2) - 2 I1/(x^2 K1)], M = gamma^T A gamma at E = 0
(Paneru et al., PRC 99, 045807 (2019), Eq. 7), written to `output/scatteringLengths.out`, and the
`<scatteringLengths>` chi2 term.  Theory: docs/source/theory/scattering_length.rst.

Models (built by `make_models.py` (deleted; see "How the models were made"), which needs mpmath; widths are observed Brune widths
Gamma = 2 P gamma~^2 / (1 + sum gamma~^2 S') so that AZURE2 recovers the published gamma~):

| model | content | reference a0 (fm) |
|---|---|---|
| `he3he4.azr` | 3He+4He, Paneru et al., PRC 109, 015802 (2024) Table V: 1/2+ background, Ex 14.000, gamma~ 1.683, a = 4.2 | 33.099574 (paper: 33.10) |
| `be7p.azr` | 7Be+p, Paneru 2019 Table II, a = 4.3: 2- (Ex 3.520, gamma~ 1.700), 1- background (Ex 9.000, gamma~ 1.433 and -1.822 in 7Be*(0.429)+p, closed) | a01 17.06776438, a02 -5.018395912 |
| `n16o.azr` | n+16O, one 1/2+ level (Ex 6.0, gamma~ 0.8, a = 4.0) | 4 (1 - 0.64/1.8569) = 2.6213582 |
| `he3he4_fit.azr` | he3he4 with E and width free, phase-shift data from he3he4 (`data/he3he4_phase.dat`), constraint 31.0 +- 0.3 | |
| `be7p_fit.azr` | be7p with every parameter free, a second 1- level, phase-shift data (`data/be7p_phase_*.dat`), every kind of constraint line (a01, a02, coherent, complex, inactive) | constraint chi2 48.3471 |

The charged references are independent of AZURE2: mpmath Coulomb functions and the Brune level matrix at
E = 0 (`7Be+p/10-20-26_scattering_length/proto_*.py` in the fit archive), confirmed there by the
effective-range extrapolation of the same levels' phase shifts.  Paneru 2019 quotes a02 = -3.18 fm,
which does not follow from its printed Table II (see the theory page).

Checks:
1. a0 against the references (1e-5; neutral 1e-7); the 7Be*+p 1- length is complex.
2. a0 from AZURE2's own phase shifts by the Coulomb-modified effective-range expansion (`ere.py`,
   standard library python3; skipped without python3), within 2e-4 (observed <= 4e-5).
3. `--no-brune` gives the same lengths (1e-8).
4. The constraint lines: coherent = 3/8 a01 + 5/8 a02, chi2 total = sum of lines (inactive 0) = the
   `Total-Scattering-Length-Chi-Squared` of chiSquared.out = 48.3471.
5. `AZURE_GRAD_CHECK` on be7p_fit (8 parameters) in the Brune, standard and Park parametrizations:
   analytic vs finite differences < 1e-6 (observed 4e-10..7e-8); Park < 1e-3 (4e-8 with the Coulomb
   table, 2.5e-4 with the exact Coulomb routine as in CI: Park's dS/dE at the level energy is a numerical
   derivative and its noise enters the finite differences, as in tests/pyazr/park_gradient_test.py).
6. A `--use-lm` fit with the constraints reaches data + constraint chi2 < 6 (from 227 + 48).
7. Malformed block lines are refused with the line named.

Runs in about 6 s on two threads.

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show 8506786e:tests/scattering_length/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Builds the scattering-length test models (committed; rerun only to change them).

he3he4.azr  3He+4He, Paneru et al., PRC 109, 015802 (2024) Table V: the 1/2+ background level
            (Ex 14.000 MeV, gamma~ 1.683 MeV^1/2, a = 4.2 fm), the only s-wave term.  a0 = 33.10 fm.
be7p.azr    7Be+p, Paneru et al., PRC 99, 045807 (2019) Table II, a = 4.3 fm: the 2- level (Ex 3.520,
            gamma~(S=2) 1.700) and the 1- background (Ex 9.000, gamma~_el(S=1) 1.433, gamma~ in
            7Be*(0.429)+p (S=1) -1.822, closed at the 7Be+p threshold).
n16o.azr    n+16O with one 1/2+ level (Ex 6.0 MeV in 17O, gamma~ 0.8 MeV^1/2, a = 4.0 fm): neutral, so
            a0 = a (1 - gamma~^2 / E~) in closed form (S = 0 for l = 0).
Widths are written as observed (Brune) partial widths in eV, Gamma~_c = 2 P_c gamma~_c^2 / (1 + sum gamma~^2 S'),
so that AZURE2 recovers the published gamma~ (parameters.out g_int); each model has
a phase-shift test segment at low energy for the effective-range check of check.sh.
Needs mpmath.
```
