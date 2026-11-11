# screening_6li — electron screening benchmark on 6Li(p,α)3He

A benchmark of the `screening Ue vary sigma` block of a `<targetInt>` line
(see `tests/electron_screening` for the unit checks) on measured data with
published screening potentials.

**Data.** Cruz et al., Phys. Lett. B 624, 181 (2005), EXFOR O1282005–007:
6Li(p,α)3He S-factors at E_c.m. = 24–82 keV on three targets, with the
screening potentials the authors derived:

| target | EXFOR | N | published U_e (eV) |
|---|---|---|---|
| Li metal | O1282005 | 9 | 1320 ± 70 |
| Li₂WO₄ (insulator) | O1282006 | 7 | 320 ± 110 |
| PdLi alloy | O1282007 | 7 | 3760 ± 260 |

EXFOR tags the energy column `EN` (laboratory), but it is the c.m. energy:
the paper's 30–100 keV proton beam gives these c.m. energies, and read as
laboratory energies the metals come out 26–27 % below the published U_e.
`make_models.py` (deleted; see "How the models were made") converts S to σ = S/E · exp(−2πη) with AZURE2's constants
and the model's masses (1.00782, 6.0151), and writes `data/cruz_*.dat`
(E_lab in MeV, 0, σ and its statistical error in b).

**Model.** J. Skowronski's 7Be R-matrix evaluation (3He+α, p+6Li, p+6Li\*,
capture; Brune parametrization), best fit written into `<levels>` with pyazr
`save_fit`, every level energy and width fixed. Its own data are removed; it
has no 6Li(p,α) points below E_p = 100 keV, so the Cruz sets are out of
sample: the test asks whether the bare low-energy extrapolation of a fit to
higher-energy data, times exp(πηU_e/E), describes them, and with which U_e.
Default CLI flags (Brune, externals, long wavelength).

| model | what it is |
|---|---|
| `bare.azr` | U_e = 0 on all three segments, norms 1 |
| `scr.azr` | U_e and norms at the benchmark minimum, fixed |
| `fit.azr` | U_e (start 500 eV) and the three norms free, no priors |

**Benchmark** (AZURE2 dev 334f104; χ² with the norm free and Δχ² = 1
intervals):

| target | U_e (eV) | norm | χ² | χ² at published U_e | χ² at U_e = 0 |
|---|---|---|---|---|---|
| Li metal | 1211 ± 65 | 1.015 | 12.49 | 15.33 | 311.5 |
| Li₂WO₄ | 323 ± 98 | 1.010 | 2.95 | 2.95 | 13.4 |
| PdLi | 3433 ± 231 | 1.035 | 8.26 | 10.28 | 201.8 |

The insulator reproduces the published value; the metals come out 8–9 %
lower, about 1σ, which is the size of difference expected from a different
bare S(E) (the published analysis did not use this fit). The adiabatic limit
for p + Li is 182 eV (NIST total electron binding energies).

`check.sh path/to/AZURE2` checks:

1. screened / bare = exp(πηU_e/E) at all 23 points to 1e-8;
2. the χ² of each segment at the benchmark minimum to 1e-4;
3. MIGRAD and `--use-gradient`, with only U_e and the norms free, find the
   minimum from U_e = 500 eV, each U_e within 5 % of its uncertainty (MIGRAD
   stops on its EDM tolerance, which leaves the weakly constrained insulator
   value about 1 eV short with the analytic gradient).

The full study (both readings of the energy column, engine-vs-profile
checks, the pyazr fit and its gradient) is in the deBoer R-matrix archive,
`6Li+p/3He+a/AZURE2_evaluation/10-21-26_screening_benchmark/`.

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show 8506786e:tests/screening_6li/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Models and data of tests/screening_6li: electron screening of 6Li(p,a)3He, Cruz et al. 2005 (Phys. Lett. B 624, 181;
EXFOR O1282005-007) on a Li2WO4 insulator, Li metal and a PdLi alloy, against J. Skowronski's 7Be R-matrix evaluation held
fixed.

The levels are that evaluation's best fit (Supplemental 7Be.azr + param.sav, written into <levels> with pyazr save_fit and
checked on re-read), every energy and width fixed; the evaluation's own segments are removed (it has no 6Li(p,a) data below
E_p = 100 keV, so these points are out of sample).  Only the three Cruz sets remain, each with its own <targetInt> line
carrying a screening block.  Written next to this script:

  data/cruz_{metal,insulator,pdli}.dat   E_lab (MeV), 0, sigma (b), dsigma (b) from the EXFOR S-factors
  bare.azr     U_e = 0 on all three, norms 1
  scr.azr      U_e and norms at the benchmark minimum, fixed
  fit.azr      U_e (start 500 eV) and the norms free, no priors

Run:  python3 tests/screening_6li/make_models.py path/to/7Be_fit.azr
(the snapshot is R-matrix/6Li+p/3He+a/AZURE2_evaluation/10-21-26_screening_benchmark/7Be_fit.azr in the deBoer archive.)
```
