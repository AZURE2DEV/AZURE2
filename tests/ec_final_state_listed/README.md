# ec_final_state_listed — an external-capture final state must be a level of the project

The user guide ("Gamma-Ray Calculations") asks for one level per gamma-ray pair: that level carries the bound
state's ANC and gamma widths. Before 2026-10-26, `CNuc::ParseExternalCapture` silently created a missing final
state as a hidden level: energy fixed at the pair's E_x, not in the R-matrix, and a FREE reduced-width amplitude
of 0.1 MeV^1/2 in every particle channel, which appeared in no `.azr` line and not in the GUI. Because the level
was appended only when external capture was on (CLI: modes with a capture segment; API: always), every later
level's `energy_N` / `width_N_c` name moved by one between runs of the same project, so a `param.sav` matched
by name landed on the wrong levels with no message (2026-10-25 audit, item 4; reproduced on 13N, where the
2.37 MeV 1/2+ level is `energy_3` with external capture and `energy_2` without).

The engine now refuses, when an active data or test segment (or the reaction-rate pair) exits into that pair, with

    ERROR: external capture into pair 3 needs its final state, J^pi = 0.5+ at E_x = 0.4953 MeV, as a level of the
    project (with its ANC and gamma widths); none is listed within 1 keV in that J^pi group.

(or "... and that J^pi group does not exist."), and `Fill` fails, so nothing is calculated. A capture pair that no
segment uses (a level-isolation or scratch model whose final-state lines are switched off, as in the 11B+a figure
decompositions) gets the same text as a WARNING, no external capture into that pair and no hidden level, and the
run goes on. The CLI only switches external capture on in a mode with a capture segment, so there the WARNING
branch is never reached; `tests/pyazr/ec_final_state_test.py` covers it through the API, which always has it on.

Models are p + 16O -> 17F(0.495 MeV, 1/2+) + gamma with E1 external capture (`tests/ec_m1_exact/models/e1.azr`):

| model | levels | result |
|---|---|---|
| listed | as e1.azr | runs; chi2 against `expected/listed_chi2` (1e-3) |
| missing_level | the 0.495 MeV lines removed, the 1/2+ group kept through its 30 MeV level | refused, no output |
| missing_group | every 1/2+ line removed | refused, no output |
| unused | missing_level, but the only segment is elastic 16O(p,p) at 90 deg (`data/pp90.dat`) | runs, three levels in param.par |

Before the fix both "missing" models ran (chi2 written, a hidden `energy_2` / `width_2_1` in the parameter list),
and the API gave the unused model four levels.
To make such a project run again, add the final state to `<levels>` at the pair's E_x with its ANC in the
particle channel(s) and its gamma-width line, as the GUI does.
