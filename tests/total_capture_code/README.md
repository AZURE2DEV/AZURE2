# total_capture_code

Total capture is the exit key `-1` together with the total-capture data type
(audit 2026-10-25, 01_kinematics_frames B2, checklist item 32).

The reader expands a total-capture segment (exit key `-1` and data type 3 on a
data line, 4 on a test line) into one term per capture pair. The key `-1` with
another data type left a segment whose exit pair does not exist:
`GetPair(GetPairNumFromKey(-1))` read `pairs_[-1]`, and in mode 3 the
external-capture lookup for that point crashed (SIGSEGV on `tests/13N`; a model
without external capture gave meaningless numbers instead). The two halves are
now checked against each other and the mismatch is refused; the segment-level
lookups use `CNuc::GetPairFromKey`, which returns null for a key no pair has.

## Checks (`check.sh path/to/AZURE2`)

- (a) `tests/13N` plus a test line `1 1 -1 0.5 2.0 0.5 0 0 0 0` (exit `-1`,
  type 0), mode 3: refused, no output;
- (b) a data line with exit `-1` and type 0 (angle-integrated), mode 1:
  refused, no `chiSquared.out`;
- (c) the proper test line `1 1 -1 0.5 2.0 0.5 0 0 0 4` beside
  `1 1 2 ... 0`: `AZUREOut_aa=1_TOTAL_CAPTURE.extrap` equals
  `AZUREOut_aa=1_R=2.extrap` point by point (13N has one capture pair).

dev 510e6d1: (a) and (b) segmentation fault (exit 139) with no output.
