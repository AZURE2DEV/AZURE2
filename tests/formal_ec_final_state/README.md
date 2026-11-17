# formal_ec_final_state — an external-capture final state in the formal parametrization

dev's `tests/formal_ec_final_state` (dev `081dddf9`, `630dff75`), brought here on 2026-11-15: the four models
converted with `AZURE2 --convert` (no note, no question), the script rewritten to commands (`extrapolate`,
`calculate`, `fit`), the pseudo-data set as a format-3 `dataset` record. The level lines are edited with dev's awk
programs unchanged (the `<levels>` block keeps its columns).

On dev the formal parametrization took the working energy and amplitude of every final state from the first level of
its J^π group, so 12C(p,γ0) depended on the order of the level lines (up to 9.7 times). This branch has its own
version (the boundary condition at the lowest level of a group; the final states' amplitudes carried by
`CompoundSystem`), and `docs/upgrade/dev_features.md` section 1 asked for dev's checks to be run against it.

## Result

50 checks pass. 45 of dev's 51 passed unchanged; the six others met two differences of design, not of numbers:

| dev's checks | here | what the script now checks |
|---|---|---|
| `--no-transform`: a NOTE names the level taken as the final state when none is within 1 keV of the pair's energy; a WARNING when the pole is off that energy; the distance is reported (4 checks) | the boundary condition is at the lowest level whatever the order, with a Note when that is not the first listed one; a final state moved 50 keV off its pair's energy is refused ("needs its final state ... none is listed within 1 keV") | the Note with the ground state second and not with it first; the refusal, in both orders, with no cross section written (4 checks) |
| the analytic gradient is used with every parameter of the final state's group free, and the background pole's width has a gradient there (2 checks) | the final state's amplitude is not differentiated analytically (`AZUREGrad.cpp`, `HasFormalECFinalBlocks`): the gradient is by finite differences, the run says so, and MIGRAD seeds numerically ("Analytic G2 unavailable for this model") | the announcement (either route is accepted, a silent third is not), and a formal fit of that model returns to the data (220.515 to below 1) (2 checks) |

The background pole's width has no gradient here in either order: with the boundary condition at the ground state
the pole does not enter its amplitude (on dev that is the case "ground state first").

## Open

- **The analytic derivative of the final state's amplitude** (dev: `CNuc::FormalECFinalPullback`) is not here. A formal
  fit with the final state's group free costs one finite-difference gradient per step and a numerical seed (about 210
  calls for 16 parameters in this model before the first step).
- **Refusal or warning** for a formal final state off its pair's energy under `--no-transform`: this branch refuses,
  dev warns and goes on. With the question of the boundary condition's level (`FOR_THE_OWNER.md`), it is the owner's.

## Checks (`check.sh path/to/AZURE2`)

Two- (A) and three-level (B) 1/2⁻ groups of 13N, the ground state first or not: order independence and formal = Brune
for capture and elastic scattering; the sign of a negative ANC; a static moment on the ground state; the γ widths
reported by a formal run; `--no-transform` round trips; a zero ANC; pseudo-data χ² in both orders and bases, the
gradient check, and two formal fits. About a minute.
