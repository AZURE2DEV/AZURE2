# progress_bars

Progress bars of the external-capture amplitude loops (2026-11-06).

`EData::CalculateECAmplitudes` (data segments) and `EData::InitializeComponentSegments` (component segments)
compute the EC amplitudes of a segment's points in an OpenMP loop and redraw a 25-cell bar
(`\r\tSegment #<key> [*********                ] 36%`). Up to dev 0119771 every thread incremented the shared
point counter without synchronization and wrote its bar to the output stream piece by piece, without a lock. The
bars of several threads interleaved character by character (`Segment #\r\tSegment #\r\tSegment #1111   [****`),
and a fragment could land on a line of its own after another thread's newline. This garbled the log and made the
`tests/light_excitation` log compare flaky on the macOS and Windows runners (dev 016bcfe). It was also a data race
on the GUI's output stream, which is not thread-safe. Now the counter is atomic and one thread at a time draws a
whole bar, under the critical section `ec_progress`, which the loops' error messages share (`ECProgressRedraw`).
The bars' format and the results are unchanged.

Models: `tests/13N` (external capture on seven capture segments), and the same with segment 1 a SUM of its own
term and one identical capture component (key -1000; exercises the component loop), mode 1.

| check | what |
|---|---|
| (a) | at `THREADS` (8) threads, `RUNS` (3) runs per model: every progress line, split at `\r`, is a sequence of whole bars of one segment ending at 100 % |
| (b) | no other line is made only of bar characters |
| (c) | one finished bar per EC segment, in segment order, and the single-thread chi2 |

The old build fails (a) and (c) on every 8-thread run of the plain model: segment 1's line holds torn fragments.
The race is timing-dependent, so on a slow machine one run might happen not to interleave; three runs make a
false pass unlikely, and a pass on the old build is not a pass of this check.
The component loop computes only the component's own points (here a few), so it seldom races even on the old
build; its bar is checked for form.

    ./tests/progress_bars/check.sh path/to/AZURE2      # ~100 s; THREADS=8 RUNS=3
