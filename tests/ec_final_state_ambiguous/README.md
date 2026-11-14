# ec_final_state_ambiguous

An external-capture final state that matches more than one level is refused (audit 2026-10-25,
03_rmatrix_core F9, checklist item 28).

`CompoundSystem::ParseExternalCapture` attaches a capture pair's final state to the level of the pair's J^pi
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

## Two photon pairs with one final state

The other way round: `tests/13N` with a copy of its gamma pair under a new key (pair 3) and a copy of its first
capture data set into that pair.  Both pairs name the ground state.  `ParseExternalCapture` gave the level the
second pair's capture parameters, and capture into pair 2 ran without its external part and without a message
(1.95e-9 against 4.32e-9 at the first point).  Now it stops with an ERROR naming both pair keys and the level
(`FinalStateOfAnotherPair`, `src/ExternalCaptureAmplitudes.cpp`):

- (d) the two-pair model writes the ERROR and no chiSquared.out;
- (e) its format-3 conversion (`AZURE2 --convert ... --to 3`) converts and is refused the same way when run.

