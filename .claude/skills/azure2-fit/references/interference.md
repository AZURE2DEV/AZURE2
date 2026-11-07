# Interference, branches and parametrization

Amplitudes from several levels feeding one channel add before squaring, so the relative signs and magnitudes of the
reduced-width amplitudes γ_λc (MeV^½) decide the cross section between and far from resonances. The objective is not
convex in those signs, so a fit can converge to a minimum on the wrong interference branch. This page states the
physics that matters for fitting, how to recognise a wrong branch, and the parametrization choice. Equations:
`references/formalism.md`. The level-scheme picture (which levels to include): `azure2-setup`.

## Same-J^π interference (energy dependence)

Two levels of the same J^π contribute amplitudes to the same channel. The cross section carries a cross term of size
about 2√(σ₁σ₂), so a weak level can change the cross section far from its own peak through interference with a strong
one. The relative sign of the reduced-width amplitudes decides whether they add (constructive) or cancel (destructive)
in the valley between them and in the tails. In ¹²C(α,γ) the constructive and destructive E1 choices change the
extrapolated S factor by a large factor and differ by hundreds in χ². Theory rarely predicts the sign; the data decide
it, and capture data between resonances are the decisive constraint. Test both relative signs of interfering
same-J^π levels; the off-resonance regions and the low-energy tail are where the sign shows.

## Different-J^π interference (angle dependence)

Amplitudes of different J^π do not interfere in the angle-integrated cross section, but they do in the angular
distribution, producing odd Legendre terms (E1–E2 interference gives a forward–backward asymmetry). This is how
multipoles are separated.

## The ghost anomaly

A broad resonance just above a threshold has a width Γ_c(E) = 2P_c(E)γ²_c that the penetrability P_c throttles to
nearly zero at threshold and then raises steeply. The line shape is skewed, the peak is pushed up, and a secondary
"ghost" bump can appear away from the nominal energy (the ⁸Be ground state in α+α; Barker and Treacy, Nucl. Phys. 38
(1962) 33). It follows from the energy dependence of P_c and S_c; do not "fix" it by adding a level.

## Wrong branches

Interference makes the objective non-convex in the interfering amplitudes' signs and relative magnitudes, so several
well-separated local minima can exist, each a different combination of constructive and destructive interference. This
holds for any amplitudes feeding one final channel: background poles against a resonance, two resonances of one J^π,
internal against external capture, several background poles, in particle and γ channels alike. A least-squares fit
(LM, TRF) started with every free amplitude at a small same-sign seed walks to the nearest branch. Reaching a better
branch that needs a larger magnitude or a sign flip means crossing worse χ² first, which a gradient method does not
do; a revert-on-bound guard built against numerical runaway cannot tell this from progress and stops it.

- **The tell.** Independent attempts (different seeds, subsets, guard bounds) keep reverting with the same parameters
  pushing toward the same large value. The recurrence is diagnostic: the optimizer is pointing at the other branch.
- **The remedy is a starting point, not a guard.** Seed the interfering amplitudes at a large-magnitude, sign-varied
  point (from a previous fit of the same or a similar reaction, or by trying several sign and magnitude combinations)
  and let the fit descend into that branch. Worked case: a ¹²C(p,γ) fit reverted three times pushing two background γ
  widths to 4–7 × 10⁵ eV; recalculating, without fitting, from a hand-picked point with those channels at
  200–1300 eV and one sign flipped dropped one segment's χ²/N from 79.8 to 17.2.
- When standing up interfering amplitudes for a new reaction, do not assume the small-positive-seed branch is the best
  one. If a fit of the same or a similar reaction exists, seed from its converged values, signs included. Systematic
  sign rounds over a whole model are rmfit's job (`rmfit-campaign`).

### Guard or runaway

When a fit reverts at a guard boundary, check how it stopped. In the same ¹²C(p,γ) fit, one background E2 channel then
reverted against a 50 W.u. cap: its optimizer had terminated on `gtol`, converged in that direction at 66.7–83.6 W.u.,
unlike the earlier true runaways (hundreds to hundreds of thousands of W.u.). With the cap at 100 W.u. the fit had no
reverts and the channel settled at 77.2 W.u., inside the range the reverted run had indicated.

- `gtol`/`xtol`/`ftol` termination with the guard flagging afterwards: a real local minimum just past the bound.
  Raising the bound is then a modelling decision (how large a W.u. value is reasonable for that pole).
- `max_nfev` exhaustion, or a value at many multiples of the guard: a runaway. Raising the bound would only hide a fit
  that has not converged.

## Formal, observed, Brune and Park parameters

The standard R-matrix parameters E_λ (MeV) and γ_λc depend on the boundary constants B_c, so they are not the observed
resonance energies and widths, and B_c = S_c(E_λ) can be chosen for only one level at a time. Brune's alternative
parametrization (Phys. Rev. C 66 (2002) 044611) writes the collision matrix in terms of the observed energies Ẽ_λ and
on-resonance amplitudes γ̃_λc of all levels: B_c does not appear, every parameter has a direct physical meaning, and
literature values go in and come out without conversion. It gives the same cross sections; the channel radii must
still be chosen. Prefer it for fits that draw on tabulated level data (`--use-brune`). The Park basis (`--use-park`)
fits observed amplitudes per channel. Basis definitions, conversions between them and the cases where a standard
equivalent does not exist: `azure2-exchange`.

## Pitfalls

- Using observed widths as formal widths, or the reverse, without the conversion (the factor
  1 + Σ_c γ²_λc dS_c/dE is often about 30 % for broad levels). The Brune basis sidesteps it.
- A wrong J^π: coupling a level to channels or multipoles it cannot reach.
- Forgetting a sub-threshold state that dominates the low-energy extrapolation.
- Not testing both interference signs of same-J^π levels: the fit can look right on resonance and be far off in the
  extrapolation.
- Over-parametrizing with background poles and reading physics into them; ignoring the radius–background correlation
  (`references/background-poles.md`).
- Chasing the lowest χ² rather than the most physical fit. Data uncertainties are rarely fully right, and the
  treatment of normalizations and systematics (Peelle's pertinent puzzle) matters.
- Treating the channel radius as a physical radius, or changing it without re-examining the background poles.
- Forgetting that γ widths set the peak height, not the total width.
- Judging a fit by a total χ²/N per data set instead of point by point (`azure2-check`).
