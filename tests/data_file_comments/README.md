# data_file_comments

`tests/7Li_p_ay` with its data files rewritten to carry a `#` header, blank
lines, an indented `#` comment mid-file, and CRLF endings on some rows.

The expected chi-squared is the unmodified 7Li_p_ay reference: comments and
blank lines must be ignored and the numbers must not move.

Before the line-based reader, any of these hung AZURE2 forever in
`ESegment::Fill` -- a line that failed `stream >>` set failbit, eof was never
reached, and the loop spun at 100% CPU with no message.  If this test ever
hangs, that regression is back; the harness timeout turns it into a FAIL.
