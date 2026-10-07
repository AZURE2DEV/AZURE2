# Sourced by the tests/*/check.sh scripts that use them: the reporting and
# comparison helpers they all had a copy of.  A script sets $fail (0) and,
# for ran/model, $WORK and $OUT (the output file of the segment) first.

# ok MESSAGE / bad MESSAGE: one result line; bad marks the check failed.
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }
# same A B: both files exist and are equal up to \r (Windows line ends).
same() { [ -f "$1" ] && [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
# ran NAME: the run in $WORK/NAME ended with status 0 and wrote $OUT; else a
# failure with the end of its log.
ran() { [ "$(cat "$WORK/$1/status")" = 0 ] && [ -f "$WORK/$1/output/$OUT" ] || { bad "run $1 failed"; tail -5 "$WORK/$1/log" | sed 's/^/        /'; return 1; }; }
# model NAME: "E model" per point of $OUT (columns 1 and 4).
model() { tr -d '\r' < "$WORK/$1/output/$OUT" | awk 'NF > 4 { print $1, $4 }'; }
# counted N: N is a positive number (a comparison that compared something).
counted() { awk -v n="$1" 'BEGIN { exit !(n + 0 > 0) }'; }
