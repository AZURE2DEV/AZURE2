# thm_lineshape — Coulomb line shape of the THM spectator

`check.sh` runs `tests/18O_p_a_thm` (two 1/2+ levels of 19F, two THM segments
with free norms, ~0.1 s a run) with `experiment[<name>] ... lineshape=on|off`
in a `<thm>` block appended to a temporary copy of the .azr
(docs/source/theory/thm_implementation.rst, "Coulomb line shape").  The real
reaction, 2H(18O,α15N)n, has a neutron spectator; a charged one is made up as
18O(3He,α15N)d at 115 MeV (E_sF ≈ 10.4 MeV, ζ ≈ −0.14).

| case | expectation |
|---|---|
| (a) `lineshape=off` | output files byte-identical to the same line without the key; one-segment experiments with `lineshape=off` byte-identical to no `<thm>` block |
| (b) spectator `n` (2H(18O,α15N)n at 54 MeV), `lineshape=on` | byte-identical to off; ζ = 0 in `thm_experiments.out` |
| (c) one level, no folding | model(on)/model(off) = exp[2ζ arctan(2(E_λ−E)/Γ_λ)] at all 118 points to 1e-6 (measured 3e-9), ζ computed here from the masses, Γ_λ from `parameters.out`; the same with the second level kept (it interferes through the level matrix) but its α width zero: the factor belongs to the decaying level |
| (d) one level, 0.5 keV grid | ζ < 0: the peak moves up (0.592 → 0.605 MeV), the model is lowered below E_λ and raised above |
| (e) refusals | `lineshape=on` without kinematics, `lineshape=1`, E_sF ≤ 0 at a data point (Ebeam = 40 MeV), `--use-rmc` (no Brune) |

`tests/pyazr/thm_lineshape_test.py` checks the same through pyazr
(`AzrModel.set_thm_experiment(..., lineshape=True)`, CLI == session,
`session.thm_lineshape` against the formulas).
