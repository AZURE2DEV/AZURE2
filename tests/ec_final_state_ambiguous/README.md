# ec_final_state_ambiguous

An external-capture final state that matches more than one level is refused (audit 2026-10-25,
03_rmatrix_core F9, checklist item 28).

`CNuc::ParseExternalCapture` attaches a capture pair's final state to the level of the pair's J^pi
group whose energy lies within 1 keV of the pair's excitation energy (`JGroup::IsLevel`, an
energy-only match that returns the first fit). Since e239a12 two levels of one J^pi may sit within
1 keV of each other (they carry different level IDs), so a pair's energy can match both: the first
one in the file took the ANC silently and the other was ignored. Now `ParseExternalCapture` counts
the candidates and stops with an ERROR naming the pair, the J^pi, the energy and the count.

## Checks (`check.sh path/to/AZURE2`)

`tests/13N` with a second 1/2- level 0.8 keV above the ground state (the external-capture final
state of its gamma pair), with its own ANC line:

- (a) the plain 13N reproduces its expected total;
- (b) the duplicate-level model stops with the ERROR and writes no chiSquared.out;
- (c) the duplicate moved to 2 keV, outside the window, runs.

dev b5faec4: (b) runs, attaching the capture to the first level.
