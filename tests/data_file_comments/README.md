# data_file_comments

`tests/7Li_p_ay` with its data files rewritten to carry a `#` header, blank
lines, an indented `#` comment mid-file, and CRLF endings on some rows.

The expected chi-squared is the unmodified 7Li_p_ay reference: comments and
blank lines must be ignored and the numbers must not move.

Before the line-based reader, any of these hung AZURE2 forever in
`ESegment::Fill` -- a line that failed `stream >>` set failbit, eof was never
reached, and the loop spun at 100% CPU with no message.  If this test ever
hangs, that regression is back; the harness timeout turns it into a FAIL.

Re-pinned 2026-09-26: the shift-function fix b6cc41b removed ~1e-9 noise from S(E) that ShftFunc::EnergyDerivative amplified into dS/dE, which enters the Brune transformation of sub-threshold levels. As in 7Li_p_ay: 3955.30 -> 3955.32, segments by up to 1.7e-5.

Re-pinned 2026-09-28: the dS/dE fix (ChannelFunc::DerivativeStep) starts the central difference of the shift function from a quarter of the distance to threshold, at most 1 keV, instead of 1 eV, where it was round-off dominated (~1e-7..1e-5 relative, and different with and without FMA -- Linux and Windows CI differed from an x86-64-v3 build by up to 1.6e-3 in tests/thm_options). As in 7Li_p_ay: 3955.32 -> 3955.41.
