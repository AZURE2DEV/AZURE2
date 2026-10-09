---
name: azure2-gui-review
description: Open the AZURE2 GUI so a human can visually check a fit or a model edit. Use whenever the task is "show me in the GUI", "let me look at it in the GUI", "open the GUI", or any request to eyeball segments, curves or the level scheme before approving a fit or a refit. Covers the recalculate-first rule, what the GUI actually reads, and the edit/close ordering.
---

# Showing a result in the AZURE2 GUI

## The rule that keeps getting broken

**The GUI plots `output/`, not the `.azr`. Recalculate before you open it, every
time, or the person is looking at stale curves.**

Opening the GUI on a `.azr` you just edited shows the *old* results, because the
plot tab reads the `AZUREOut_*.out` files sitting in the output directory. Those
files were written by whatever run last happened — possibly a different segment
set, a different parameter file, even a different binary. The segment list and
level tabs update from the `.azr` immediately, so the window looks correct while
the curves silently disagree with it. That mismatch is worse than no plot at all:
it invites a review decision based on a picture of something else.

There is no warning. Nothing in the GUI marks the plots as out of date.

So the sequence is always:

1. finish every `.azr` edit
2. **close any running GUI**
3. run a calculate (CLI mode 1)
4. confirm it wrote fresh output
5. open the GUI

Never reorder these. In particular do not leave a GUI open while editing the
`.azr` on disk — the GUI holds its own in-memory copy and can write it back over
your edit on exit, and the binary stays mapped (an attempt to replace
`~/bin/AZURE2` while a GUI is running fails with "Text file busy").

## Doing it

```bash
cd <fit directory>

# 2. close any GUI
pkill -f "bin/AZURE2 <name>.azr"; sleep 1
ps -eo pid,cmd | grep "[b]in/AZURE2 <name>.azr" || echo "GUI closed"

# 3. calculate -- menu 1, external parameter file, blank EC file, exit
~/bin/AZURE2 <project flags> --no-gui --no-readline <name>.azr <<EOF
1
output/param.sav

7
EOF

# 4. confirm
grep Total-Chi output/chiSquared.out
ls -la output/AZUREOut_*.out | head       # timestamps must be from just now

# 5. open
nohup ~/bin/AZURE2 <name>.azr > azure2_gui.log 2>&1 &
sleep 2; ps -eo pid,cmd | grep "[b]in/AZURE2 <name>.azr"
```

`pgrep -f` will match your own shell wrapper as well as the GUI; `ps -eo
pid,cmd | grep "[b]in/AZURE2"` is the check that actually tells you whether the
process is alive.

## Open it on the screen the person is looking at

A Claude session inside tmux keeps the `DISPLAY` it had when it started. After the
person logs in again, that X forwarding is gone, and the GUI dies at once with
"The X11 connection broke (error 1). Did the X11 server die?" in its log. Nothing
appears, and `nohup` hides it. Seen 2026-10-08: the session had `localhost:12.0`,
whose port 6012 was no longer listening.

Find the display of the login that is showing this session, and pass it explicitly:

```bash
tmux list-clients -F '#{client_tty} #{session_name}'      # which login tty shows which session
p=$(ps -t pts/0 -o pid= | head -1)                          # the shell on that tty
tr '\0' '\n' < /proc/$p/environ | grep '^DISPLAY='          # its DISPLAY
ss -ltn | grep 127.0.0.1:60                                 # live forwards: 60NN = localhost:NN
DISPLAY=localhost:NN.0 nohup <binary> <name>.azr > azure2_gui.log 2>&1 &
sleep 5; cat azure2_gui.log                                 # must be empty of X11 errors
```

## Things that make the picture wrong in ways the GUI will not tell you

**Pass the project's own CLI flags.** CLI mode does not read Runtime Options
from the `.azr` — `--use-brune`, `--ignore-externals`, `--no-long-wavelength`
and friends only take effect if passed on the command line. Take them from the
project's job script (`run_crc_*`), not from memory and not from a default. A
calculate run under different flags than the fit produces plot curves computed
with different physics from what the `.azr` describes, and the GUI shows them
without complaint.

**Use the same binary the review is about.** `~/bin/AZURE2` gets rebuilt; the
GUI and the calculate must be the same build or the curves and the person's
mental model diverge. Check `ls -la ~/bin/AZURE2` and say which build is in
play when handing over.

**Give it the right parameter file.** A calculate with a *blank* external
parameter file uses the `<levels>` block; with `output/param.sav` it uses the
fitted values including normalizations and energy shifts. These can differ by a
large factor (on one 11B+a model, 127,588 vs 96,811). Blank is the right choice
only when the point is to verify what `<levels>` itself contains -- and note that
"blank" does not mean norm = 1: it means the `dataNorm`/`energyShift` fields stored in
`<segmentsData>`. Those are the NOMINAL values (the systematic-error penalty is
measured from them); no fit and no `save_fit` ever updates them, and they must not
be overwritten with fitted values (azure2-eval skill, "Normalizations live in two
places"). The fitted values are only in the `.sav`.

**Back up `output/` before recalculating** when the existing contents are a
record someone still cites — a published baseline, the numbers in a manuscript,
anything you have been quoting in conversation. The calculate overwrites
`chiSquared.out`, `normalizations.out`, `parameters.out` and every
`AZUREOut_*.out` in place. Copy them to a sibling directory first; they are
small (~1 MB) and may not be reproducible if the binary has since changed.
`.extrap` and `.acoeff` files are not touched by a data-mode calculate.

## Budget the time

A plain mode-1 calculate is not instant on a large multichannel model with
target integrations. Measured on the 11B+a fit (107 segments, ~4,800 points,
`<targetInt>` entries up to 1,000 integration points):

| | |
|---|---|
| fresh directory, empty `output/` | **~33 min** |
| live directory with populated `output/` | ~5 min |

Tell the person the calculate is running rather than letting them wait at a
closed GUI, and queue it (`qsub`) rather than running it on a login node — these
runs pull several GB and effectively one core regardless of `OMP_NUM_THREADS`.

## Before handing the window over

Say what they are looking at, because the GUI itself does not:

- which `.azr`, and what changed in it since last time
- which binary (date), and which parameter file the curves came from
- the total χ² from the calculate you just ran
- which tabs are worth checking for this particular change

If the change was to `<segmentsData>`, name the segment **keys** you touched and
what file each one holds. A segment's key is its 1-based position among *all*
lines in `<segmentsData>`, counting inactive ones — the same key `<targetInt>`
binds to, and the same key `chiSquared.out` and `normalizations.out` print.
Those two output files are the authoritative cross-check on any claim about
which segment is which; quote them rather than a count you did by eye.
