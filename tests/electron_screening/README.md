# electron_screening — the screening block of a `<targetInt>` line

A `<targetInt>` line may carry `screening Ue vary sigma` (U_e and sigma in eV):
the bare cross section of its segments, at every point and every sub-point of
the line's other effects, is multiplied by f(E) = exp(pi eta U_e / E), with eta
the Sommerfeld parameter of the entrance pair and E the c.m. energy. With
`vary` 1, U_e is the fit parameter `screening_<n>`; `sigma` > 0 adds the prior
((U_e − U_e⁰)/sigma)² to chi-squared.

The models are `tests/13N` (12C(p,γ)13N) with only the low-energy capture
segments 6 and 7 (LUNA HPGe and BGO) and 8 (FELS) active and 13N's own effects
switched off; `make_models.py` writes them, and `check.sh` links in
`tests/13N/data`.

| model | what it is |
|---|---|
| `bare.azr` | no screening |
| `scr.azr` | U_e = 1000 eV, fixed, on segments 6 and 7 |
| `conv.azr` | the same inside a Gaussian beam spread of 1e-6 MeV |
| `grad.azr` | U_e free with a 200 eV prior, inside a 2 keV beam spread |
| `bad.azr` | a screening block without its sigma |

`check.sh path/to/AZURE2` checks:

1. screened / bare = exp(pi eta U_e / E) on segments 6 and 7 to 1e-9 (AZURE2's
   constants and the pair's masses), segment 8 unchanged;
2. inside a 1e-6 MeV beam spread the result equals plain screening (1e-5), so
   the factor reaches the sub-points;
3. `--use-gradient` with `AZURE_GRAD_CHECK=exit`: d chi2 / d U_e agrees with
   finite differences in both the adjoint and the Jacobian path;
4. data made from `scr.azr` (1 % errors, every level fixed): MIGRAD with the
   analytic gradient and `--use-lm` recover U_e = 1000 eV from 300 eV to 1e-3;
   a prior of 300 ± 50 eV pulls the fit between the two, and
   `Total-Screening-Chi-Squared` equals ((U_e − 300)/50)²;
5. the malformed block is refused with a message.

The pyazr side (helpers, parameter record, penalties, gradient and Jacobian
column) is `tests/pyazr/screening_test.py`; the GUI round trip and the
adiabatic limit are in `tests/gui/target_int_tab_test.cpp`.
