# An extract of an ENDF/B-VIII.1 evaluation (test data of `pyazr.endf`)

`h-002_He_004_extract.endf`: records of the ENDF/B-VIII.1 evaluation of 3He + 4He (incident-helion sublibrary,
MAT 228; I.J. Thompson (LLNL) and the INDEN light-element collaboration, EVAL-OCT22), kept for the tests of
`pyazr.endf` with the PI's agreement (2026-11-14).

- **Source**: NNDC, ENDF/B-VIII.1 release, helion sublibrary archive (`helium3s`, MD5
  `402f3675ca8fb2b4dc02624a1de7ca3a` as published); the file `h-002_He_004.endf` of that archive (257 564 bytes, MD5
  `badff4280df4053830547ee2e6d0d581`).
- **Records kept, unchanged, in their order** (182 lines, 13 680 bytes, MD5 `210938d3e02baa26f114bc19a1ab83e9`):
  the tape's first line (TPID); MF=1/MT=451 whole (header, text and directory, 51 lines) with its SEND and FEND;
  MF=2/MT=151 whole (84 lines) with its SEND and FEND; MF=3/MT=600 (36 lines) with its SEND and the FEND of MF=3;
  MEND and TEND.  **Left out**: MF=3/MT=2, MF=6/MT=2 and MF=6/MT=600 (and the SEND of MF=3/MT=2, the FEND of MF=6).
  The directory of MF=1 still lists the sections left out: the format check reports exactly that, and the tests
  expect it.
- **Citation** asked by the archive's LICENSE: D.A. Brown, M.B. Chadwick, R. Capote, et al., "ENDF/B-VIII.0: The 8th
  Major Release of the Nuclear Reaction Data Library with CIELO-project Cross Sections, New Standards and Thermal
  Scattering Data", Nuclear Data Sheets 148, 1-142 (2018).
