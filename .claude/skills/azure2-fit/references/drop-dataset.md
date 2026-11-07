# Dropping a data set from a fit, and bringing it back

Removing one or more `<segmentsData>` segments from a fit (leave-one-out, uncertainty runs, a speed-up variant, a
set judged unusable) must change exactly that segment's contribution to the objective and nothing else. This page
gives the rule (deactivate, never delete), where the change must be recorded, how to make it, and the numerical check
that proves it did what was intended.

## The rule: deactivate, never delete

Set the segment's include flag (first token of its `<segmentsData>` line) from `1` to `0` and leave the line where it
is. Everything else in the project is keyed to a segment's position:

| keyed by segment position | what happens if a line is deleted |
|---|---|
| `segment_<key>_norm`, `segment_<key>_energy_shift` in `param.sav` | every later segment silently takes its neighbour's norm and shift |
| `<targetInt>` key lists (`"60-63"`, `"9-12,107"`) | target and resolution functions land on the wrong data |
| `chiSquared.out`, `normalizations.out`, `shifts.out` row labels | every table quoted so far changes numbering |
| readme notes, plot scripts, watcher scripts | refer to the wrong data with no error |

A key is the 1-based position among all lines, inactive ones included. A deactivated line keeps its key and everyone
else's. If a line really must be deleted, remap every per-segment parameter family (`references/parameters-and-sav.md`,
"Segment keys").

With the line only deactivated, an existing `param.sav` is safe to reuse as is: `AZUREParams::ReadUserParameters`
matches entries by name, so the dropped segment's `segment_<key>_*` lines match nothing and are ignored. Do not clean
them out of the file: a filtered file is one more thing that can be wrong, and leaving them makes re-activation a
one-character change.

## The drop must reach the model being run

A fit directory accumulates variant `.azr` files (uncertainty model, extrapolation model, polarization variant,
campaign seed, GUI copy). Deactivating a segment in one of them does nothing for the others, and a variant built later
from the master inherits the master's flag. This is how a segment believed dropped returns in a day-long refit. Before
submitting anything, list which files have the segment active:

```bash
# which .azr files in this fit have segment KEY active?
KEY=103
for f in *.azr */*.azr; do python3 - "$f" $KEY <<'PY'
import re,io,sys
c=io.open(sys.argv[1],encoding='latin-1').read(); k=int(sys.argv[2])
m=re.search(r'<segmentsData>(.*?)</segmentsData>',c,re.DOTALL)
rows=[l.split() for l in m.group(1).split('\n') if l.split()] if m else []
if len(rows)>=k:
    t=rows[k-1]; f=[x for x in t if x.endswith('.dat')]
    print("%-40s key %d active=%s  %s"%(sys.argv[1],k,t[0],f[0] if f else '?'))
PY
done
```

Decide where the drop is recorded as the fit's state. If it is a property of the analysis ("this data set is not
used"), change the master `.azr` (in a new dated directory, or after a `.bak_pre_*` copy) and say so in the readme,
so every later variant inherits it. If it is a property of one run ("off for the uncertainty runs only"), it belongs in
that run's `.azr`, and the readme says that the master still has it on.

## Doing it

Edit by position and assert what you are touching. `.azr` files are Latin-1; read and write them as such.

```python
import re, io
KEY, EXPECT = 103, "nToF_2023_Ecorr"          # key and a substring of its data file
c = io.open("model.azr", encoding="latin-1").read()
m = re.search(r'(<segmentsData>)(.*?)(</segmentsData>)', c, re.DOTALL)
rows = m.group(2).split('\n'); key = 0
for i, l in enumerate(rows):
    if l.split():
        key += 1
        if key == KEY:
            assert EXPECT in l and l.split()[0] == '1', l      # right line, currently active
            j = l.index('1'); rows[i] = l[:j] + '0' + l[j+1:]  # flip the flag only, keep spacing
io.open("model_noKEY.azr", "w", encoding="latin-1").write(
    c[:m.start(2)] + '\n'.join(rows) + c[m.end(2):])
```

In pyazr the same edit is `mdl.set_segment_active("<name>", False)`; when a data file appears in several segments
(angles, energy windows), check afterwards that exactly the intended segments changed.

What needs no change, and what might:

- **`<targetInt>`**: leave it. An entry bound only to inactive segments costs nothing. An entry shared with segments
  that stay active (`"6,102"`) must stay.
- **`<segmentsTest>`**: independent of data segments.
- **External-capture cache**: `output/intEC.dat` depends on the active segment set only if the dropped segment was a
  capture segment. If it was, delete `intEC.*` in the output directory before running.
- **Output directory**: give the variant its own (first path line of the `.azr`), so it cannot overwrite the reference
  run it is compared against.

## Verify: the difference is exactly the dropped segment

Run mode 1 of the new `.azr` with the same parameter file as the reference and compare `chiSquared.out`. Use the same
file, not a blank one: with a blank file both runs use the nominal norms, not the fitted ones
(`references/normalizations-and-shifts.md`).

```python
def chi(f):
    o = {}
    for l in open(f):
        p = l.split(',')
        if p[0].strip().isdigit(): o[int(p[0])] = (float(p[1]), int(p[2]), float(p[3]))
    return o
ref, new = chi("output/chiSquared.out"), chi("output_noKEY/chiSquared.out")
assert KEY in ref and KEY not in new
tot = lambda x: sum(v[0] for v in x.values())
print("total %.3f -> %.3f, expected %.3f" % (tot(ref), tot(new), tot(ref) - ref[KEY][0]))
print("max change in a remaining segment:", max(abs(new[k][0]-ref[k][0]) for k in new))
print("max norm change:", max(abs(new[k][2]-ref[k][2]) for k in new))
```

Pass criteria: the dropped key is absent; the segment count and `Total-N` fall by exactly that segment's; the total
falls by exactly that segment's χ²; no remaining segment and no normalization moves beyond the file's printing
precision. Worked case (key 103, 45 points, χ² = 44.770): 75 → 74 segments, 4654 → 4609 points, total
10,790.481 → 10,745.712 against an expected 10,745.710; largest change elsewhere 0.011, largest norm change 10⁻⁷.

Anything else means the edit touched more than the flag: a shifted key, a different parameter file, different CLI
flags, or a different binary. If the dropped segment had a penalized norm (`normErr != 0`), its `Norm-Chi-Squared`
leaves the total too; compare `Total-Norm-Chi-Squared` as well.

## The refit afterwards

- Modes 2 and 3 ask `Calculate cross-section uncertainty band? (y/n)` before the parameter-file prompt. Scripted
  input is `2 / n / param.sav / <blank EC> / 7`. Confirm `Reading User Parameter File...` in the log; `Creating New
  param.par File...` means the file was ignored and the fit started from `<levels>`.
- The first χ² the fit prints must equal the verified mode-1 total. If it does not, stop the job: it did not start
  where you think.
- Flags (`--use-brune`, `--ignore-externals`, ...) are CLI-only: copy them from the reference run's job script.
- Run it on a batch node, not on a login node.

## Record it

In the fit's readme: which key, which data file, why it was dropped, in which `.azr` files it is off (and in which it
is deliberately still on), and the verification numbers. "Segment 103 disabled" without saying where is how the
confusion starts.

## Bringing a data set back

Flip the flag to `1` in the same place. Its `segment_<key>_*` entries are still in `param.sav` if they were left
alone, so it returns with its last fitted norm and shift. Verify in reverse: the total rises by that segment's own χ²
at those parameters.
