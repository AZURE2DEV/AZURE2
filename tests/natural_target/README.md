# natural_target: several compound nuclei in one project

`check.sh` (run by `tests/run_tests.sh`) tests natural-abundance targets: two p + C nuclei (13N from 12C, 14N from 13C)
in one file, a target with isotopes 12C and 13C, abundance-weighted sums, abundances as fit parameters
(`<target>.isotope.<nucleus>`), and the refusals. What each part checks is in the script's header.

## The models (`models/`, format 3)

They were written once by a generator (removed: it edited old keyword lines) and converted with
`AZURE2 --convert ... --to 3 --exact-grids`; the originals are in `tests/convert/fixtures/tests/natural_target/models/`.
Both nuclei are built on the three-level p + 12C skeleton of `tests/convolution_unify/models/fixed.azr`:

| file | what |
|---|---|
| `A.azr` | p + 12C -> 13N (pair key 1), the skeleton as it is |
| `B.azr` | p + 13C -> 14N (pair key 2): the 13C masses and spin, S_p(14N) = 7.5506 MeV, levels at the same c.m. energies so the resonances sit in the same window, other widths |
| `M.azr` | A's and B's levels in one file, target `natural` with 12C 0.6 and 13C 0.4 fixed; extrapolations x1 (pair 1 alone), x2 (pair 2 alone), x3 = the sum 12C(p,p)12C + 13C(p,p)13C on the natural target (= 0.6 x1 + 0.4 x2), x4 = pair 2 on the natural target (= 0.4 x2) |
| `A_thick.azr`, `B_thick.azr` | A and B under a target-integration effect (no beam spread) |
| `M_thick.azr` | the same for M: the old file's one `<targetInt>` line served x1, x2 (no sample) and x3, x4 (the natural sample); in format 3 the layer is written twice, `tgt_targetInt1` for x1, x2 and `tgt_targetInt1_2` with the isotopes for x3, x4 (the converter's copy of a target with no fitted parameter; outputs equal to the old file's) |
| `M_samples.azr` | M plus a second target `enr` (0.1 / 0.9, its own sum-to-one width): x5 = the sum on `enr` = 0.1 x1 + 0.9 x2, x6 = pair 2 on `enr` = 0.9 x2 |
| `M_samples_data.azr` | two data sets on two targets, four abundances free: the fit must find 0.6 / 0.4 and 0.1 / 0.9 |
| `M_data.azr` | x3 as a data set reading `data.dat` (written by `check.sh` from A and B), abundances 0.5 / 0.5 free with an `isotope_sum` prior of width 0.01: the fit must find 0.6 / 0.4 |
| `bad_level.azr`, `bad_ratio.azr` | a level with channels of both nuclei; a ratio on a natural target: refused by the run |

`bad_key.azr` (a component naming pair key 9) and `bad_sample.azr` (a sample `<abundances>` lacks) have no format-3
form; their old files are kept in the fixtures, and the converter refuses them with the old engine's reason
(`check.sh` part (d)).

`data.dat` and `data_enr.dat` are not committed: `check.sh` writes them (0.6 sigma_A + 0.4 sigma_B and 0.1 sigma_A +
0.9 sigma_B on the 2 keV lab grid of the extrapolations, 0.1 % errors), so the data sets of `M_data.azr` and
`M_samples_data.azr` cannot be set up in the tree itself (a listed class of `tools/format3_gate.py`).
