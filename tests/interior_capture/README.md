# interior_capture

Capture from the interior of a hybrid potential through the whole engine (`docs/upgrade/interior_potential.md`
section 6). `check.sh path/to/AZURE2` (ctest `interior_capture`) runs `models/direct.azr` — alpha + 12C with the
Woods-Saxon of Koonin, Tombrello and Fox in every l, `interior_boundary 0`, E1 capture into the potential's own 0+
bound state, no 1- level — at the channel radii 4.0, 5.35 and 6.5 fm: the cross section is the potential model's
direct capture and must not depend on the radius. The header of `check.sh` lists the checks.

The radial integrals are checked against an independent calculation by ctest `upgrade_interior_capture`
(`tests/upgrade/interior_capture_test.cpp`, reference `tests/upgrade/interior_capture/reference.csv` from
`make_reference.py`).

The final state's energy, E_x = 7.16191 - 2.427561455 MeV, is the bound state of this potential for the masses of
the project (4.0026 and 12 u) as `make_reference.py` finds it (row `bound`).

`models/unbound.azr` is the same potential with E1 capture from l = 2 and 4 into the potential's own l = 3
resonance, given as a 3- level at E_x = 10.4229 MeV with a width of 173.042 keV: the state `AZURE2` reports for this
potential at 5.35 fm (the set-up warning that the level lies on a resonance of the potential, which is the intent
here). `temperatures.dat` holds the temperatures of the rate check, whose reference is made in the run from a dense
extrapolation by `tests/reaction_rate/truth.py`.
