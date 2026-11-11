# ec_exit_position — the external-capture exit is found by the photon pair's number

The pathways from an entrance pair are listed exit by exit, in pair order, but only for the pairs the entrance has a
pathway to: a beta pair and a pair of another compound nucleus are skipped. The external-capture code (the amplitudes
at set-up, the `intEC.dat` file and its signature, the component segments, the amplitude at a shifted energy) used to
take the photon pair's *number* as the exit's *position* (`docs/upgrade/cnuc_split.md` §11.9, defect 1). With a
skipped pair numbered below the photon pair it read another exit or past the end of the list. Now the exit is looked
up by the pair number (`include/ExternalCaptureExit.h`), and an entrance without a pathway to the photon pair is
refused with a message naming the data segment and the two pair keys.

## Models (`models/`, written by hand)

Both hold the levels of `tests/13N` verbatim and a 1/2+ level of 13C at 2 MeV with one n + 12C channel (pair key 3;
energy and width fixed, no data). The 13C level is another compound nucleus, so p + 12C has no pathway to n + 12C.
The data sets (files copied from `tests/13N/data`): the LUNA HPGe capture cross section three times (as the
reaction, as `[total]`, and as the sum `12C(p,g)13N + 0.5*12C(p,g)13N`, which gives a component segment), the
Kettner 55° angular distribution, and the Meyer 84° elastic scattering.

| model | pairs (number: key) | photon exit of p + 12C |
|---|---|---|
| `neutron_pair_first.azr` (13C level first) | 1: n + 12C (3), 2: p + 12C (1), 3: photon (2) | position 2, pair 3 |
| `photon_pair_first.azr` (13C level last) | 1: p + 12C (1), 2: photon (2), 3: n + 12C (3) | position 2, pair 2 |

The two describe the same physics; only the order of the pairs differs. The pathways check (`file`) is on, so each run
writes `checks/pathways.chk`, whose external rows show the exit's position ("Decay #") beside its pair ("Decay Pair #").

## Checks (`check.sh path/to/AZURE2`)

1. The check files: the photon exit is position 2 of pair 3 (`neutron_pair_first`) and position 2 of pair 2
   (`photon_pair_first`).
2. Both models run.
3. `neutron_pair_first`'s `intEC.dat` holds 1560 amplitudes ((32 + 32 + 32 + 32 + 67) points × 8 pathways), and it,
   its signature, every `AZUREOut_*` file and `chiSquared.out` equal those of `photon_pair_first` byte for byte. The
   expected values come from the model whose positions equal the pair numbers, which the defective look-up read
   correctly. `param.par`, `parameters.out` and `scatteringLengths.out` list the levels in file order and are not
   compared.

The refusal is checked by `tests/upgrade/external_capture_exit_test.cpp` (ctest `upgrade_external_capture_exit`) on
the same model: asking for the exit from p + 12C to n + 12C throws with "data segment #7" and "from pair 1 to pair 3".

Before the fix (base `feature/upgrade-4-data` at 2363abac): `neutron_pair_first` wrote an empty `intEC.dat` (the exit
read past the end of the list had no pathways) and died with SIGSEGV (exit 139) in the A-matrix calculation;
`photon_pair_first` ran and gives the outputs above.
