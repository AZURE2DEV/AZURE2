# The owner's answers to `FOR_THE_OWNER.md`

R. J. deBoer's answers of 2026-11-18 to the questions `FOR_THE_OWNER.md` puts to the owner (sections 3, 7, 8 and 10
as they stand at `a405a9e9`). Each row names the place of the question, the answer, and whether the answer is this
branch's proposal. Five answers repeat decisions given earlier (2026-11-15 and 2026-11-16) and are marked so.
An answer says what is wanted; none of the code changes that follow from it is made by this page.

## Section 3: where the numbers differ from `dev`

| Row | Question | Answer |
|---|---|---|
| 1 | Factor 2 for identical entrance particles | As proposed: the factor 2 stays on every reaction cross section (decided 2026-11-16). |
| 2 | Meaning of column 32 of the levels block | As proposed: the THM binding energy (decided 2026-11-16). |
| 3 | Boundary condition of the standard basis | As proposed: B_c at the lowest level of the Jπ group (decided 2026-11-16). |
| 4 | Two levels of one Jπ very close in energy | As proposed: 1e-6 MeV with the mean S′, on `dev` too. |
| 5 | Where the heavy-ion Coulomb repair lives | Keep it in `physics/coulomb/coulsm` as a listed local change. |
| 8 | Acceptance rule for a change of the folding quadrature | As proposed: 1.2 times the largest and 1.03 times the r.m.s. distance of the old engine. |
| 13 | Capture as a later term of a sum or ratio | As proposed: the external part is computed for every term. |
| 14 | One final state on two photon pairs | As proposed: such a project is refused with a message. |
| 15 | Analytic gradient of an E1-only or E2-only sum | As proposed: take the fix. |
| 16 | A MIGRAD fit stopped by its call limit | As proposed: end at the lowest objective seen. |

## Section 7: the one request

| Question | Answer |
|---|---|
| Write the next features against this branch | Agreed. `dev` gets no further features (in effect since 2026-11-15); new work is cut from `feature/upgrade`. |

## Section 8: what is asked of the owner

| Question | Answer |
|---|---|
| `--gsl-coul` on `dev` | As proposed: retire the flag on `dev` too. |
| The covariance written to ENDF MF=32 | The fit's covariance carried through the transformation (not HESSE in a formal-basis session, which is what is built). |
| The two print-only debug switches (`AZURE_GRAD_DEBUG`, `AZURE_CHANNEL_DEBUG`) | Drop both. |
| `rmfit-campaign` | As proposed: the skill lives with rmfit; the ported copy here is dropped. |
| The logo | `dev`'s new logo. |
| The format of dev's new blocks | As proposed: `channel_radius`, `isospin_*`, `spin_observable` with an operator field (decided 2026-11-16). |
| What an energy shift means for three kernels | All three should follow the shifted energy. |
| Fitted target parameters in the band | Include them in the band. |
| Minuit2's call counter | Apply the two-line change to the bundled Minuit2. |
| Two constants kept so that no result moves | Update both to the current values and accept that results move. |
| Publishing the Coulomb table | Do not publish; users build it locally (about 5 minutes). |
| A Windows wheel | MinGW, as the Windows program is built now. |
| SRIM's two data files | Do not ship them; users supply their own SRIM files. |
| `USE_READLINE` in `.github/workflows/build.yml` | Remove it from the workflow. |

## Section 10: what is not finished

| Item | Question | Answer |
|---|---|---|
| 1 | Which tests CI runs on each platform | Everything on Linux; a shorter set on Windows and macOS. |
| 1 | Line endings in the Windows CI checkout | As proposed: keep the checkout as committed. |
| 8 | `numcmc/` and the GSL fitter | Remove both. |
| 9 | Weisskopf-unit coefficients | Switch to the closed form. |
