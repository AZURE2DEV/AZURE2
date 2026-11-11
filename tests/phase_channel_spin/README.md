# phase_channel_spin

A phase-shift segment is one (J, l, s) partial wave: the `spin <s>` keyword, and an ambiguous
(J, l) is refused.

A phase-shift segment is specified by J and l. Up to dev 62c7ed2 `GenMatrixFunc` accumulated the
diagonal U-matrix elements of every channel spin with that (J, l) in one `TempTMatrix` and reported
half the argument of the sum (audit 2026-10-25, checklist t3-23). For a spin-0 target, or any
(J, l) with a single channel spin, that is the phase shift; for a target with spin (p + 3He 1P1 and
3P1, p + 7Be J = 1 l = 1 with s = 1 and 2) it is the phase of neither wave, and the
`assert(NumTempTMatrices() <= 1)` could not see it because both spins landed in one entry.

Now a segment line (data or test) may end with `spin <s>` (SegLine / ExtrapLine keyword blocks,
`ESegment::GetPhaseS`, `EPoint::GetPhaseS`); the phase branch keeps only that channel spin, and
`EData::Initialize` refuses a phase-shift segment whose (J, l) has several channel spins in the
level scheme without the keyword, or whose spin has no channel, naming the spins. The GUI segment
dialogs have an **s** field next to J and l (`SegmentKeywords::phaseSpin / withPhaseSpin`), pyazr
reads the keyword into `Segment.phase_S` / `TestSegment.phase_S`, and older engines ignore it.

## Checks (`check.sh path/to/AZURE2`)

Model from `tests/scattering_length/be7p_fit.azr` (p + 7Be, s = 1 or 2) plus two 1+ levels with
l = 1 channels: level A (Ex 1.2 MeV) couples only to s = 1 and level B (Ex 1.7 MeV) only to s = 2,
so each wave is independent of the other level. The truth for each spin is the same model with the
other spin's channel lines removed, which is unambiguous and the same on every build. One
phase-shift data segment J = 1, l = 1 at lab 0.3 to 2.0 MeV.

- (a) `spin 1` on the full model equals the s = 1-only model at every energy, (b) likewise `spin 2`;
- (c) the two waves differ by up to 69 degrees, so (a) and (b) have teeth;
- (d) the full model without `spin` stops with an ERROR naming s = 1 and 2 and writes no chiSquared.out;
- (e) `spin 1.5` stops with an ERROR.

dev 62c7ed2: (a) and (b) fail by 34 degrees (the keyword is ignored and both spins summed), (d) and
(e) run. `tests/scattering_length` (single-spin (J, l) segments) is unchanged.
