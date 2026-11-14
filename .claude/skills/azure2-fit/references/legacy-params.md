# Old-format parameter files

A parameter file (`param.sav`, `param.par`, `param.fit`) written by an old build names its level parameters
differently from the current engine. The engine and pyazr read parameter files by exact name and ignore unknown
names without a message, so an old file applied to a current build loads only its `segment_*` lines: every level
stays at its `.azr` value. This page says how to detect such a file, convert it, and verify the conversion.

## The two naming schemes

| old builds | current builds (`CNuc::FillMnParams`) |
|---|---|
| `j=<J>_la=<L>_energy` | `energy_<N>` |
| `j=<J>_la=<L>_ch=<C>_rwa` | `width_<N>_<C>` |

`N` counts levels in J-group order (j ascending, then la within the group): the same double loop that wrote the old
names, so the mapping is deterministic. `C` is unchanged, and so are the values (level energy in MeV, reduced-width
amplitude). `segment_<k>_norm` lines have the same name in both.

What the failure looks like: the log prints `Reading User Parameter File...`, nothing reports a mismatch, and the
result is that of the `.azr` levels. pyazr has nearly the same blind spot: `pyazr.Session(file, parameter_file=path)`
only warns (`W_PARAMETER_FILE_NAMES`, the names it skipped), `s.parameters.read(path)` says nothing, and both raise
only when *no* name matches; the `segment_*` lines always match.

Scale of the problem: in one group's archive surveyed in 2026, 6,997 of 11,073 parameter files were old format. One
fit with archived `chiSquared.out` total 10,422 gave 133,516 with its own unconverted `param.sav`; its α₀ reaction
rate at T₉ = 1 came out at the `.azr` levels' value (9.41 × 10⁵ cm³ mol⁻¹ s⁻¹) instead of the fit's (9.70 × 10⁵).

Use this page before any `calculate`, `extrapolate`, `minos` or `rate` run, pyazr `parameter_file=` or `parameters.read` call, or refit that starts
from an archived fit's parameter file, and whenever an archived `chiSquared.out`, rate or extrapolation will not
reproduce.

## Detect

```bash
grep -v '^#' output/param.sav | head -1          # an old file starts with  j=1_la=1_energy
```

A whole directory tree (prints the old-format files):

```bash
find . -name '*.sav' -o -name 'param.par' -o -name 'param.fit' | while read f; do
  case "$(grep -v '^#' "$f" | head -1 | awk '{print $1}')" in j=*) echo "$f";; esac; done
```

## Convert

```bash
python3 "$AZURE2_SRC/.claude/skills/azure2-fit/scripts/convert_sav_names.py" \
    output/param.sav output/param_newnames.sav --check <fresh>/output/param.par
```

- Write to a new file. Archive directories are never overwritten in place; put the converted file in the new dated
  fit directory.
- Value columns are copied byte for byte. Only the names change.
- Exit status: 0 converted (and check passed), 1 no old names (nothing written), 2 check failed.
- The converted file carries no `#parametrization` tag, like the original. An untagged file read by a `--use-park`
  run is taken as Brune amplitudes (`azure2-exchange`). Run in the formalism the original fit used, taken from the
  directory's job script or readme.

## Verify (all three, in order)

1. **Names and structure: `--check`.** `<fresh>/output/param.par` must be written by the current build from the same
   `.azr`. Make it in a scratch copy so the archive's `output/` is not touched:
   ```bash
   mkdir -p $SCR/output $SCR/checks && cp -rp model.azr data $SCR/ && cd $SCR
   "$AZURE2_BINARY" calculate <the fit's flags> model.azr
   ```
   The log must say `Creating New param.par File...`. The check requires every level name to match both ways and
   every width that is exactly 0 in one file to be 0 in the other. Widths that are 0 in the `.azr` are fixed and
   cannot become nonzero in a fit, so a pattern mismatch means the levels or channels pair up differently: the
   `.azr` structure changed after the fit, or the file belongs to another model. A deliberately wrong pairing (a
   parameter file of an older model of the same reaction) failed with 257 of 331 names in common and dozens of
   pattern mismatches.

2. **χ² against the archived `chiSquared.out`.** Run `calculate` with the converted file (scratch copy again) and compare
   per segment. First check that the archived `param.sav`, `chiSquared.out` and `normalizations.out` share a
   timestamp (`ls --time-style=full-iso`): written together means they came from one run. Check also that the `.azr`
   is not newer than them; an `.azr` edited after the fit cannot be expected to reproduce it.

3. **What to expect.** A perfect match is unlikely across a year or more of engine changes. One worked case (10B+p,
   current binary against a fit three years older):
   - unconverted 133,516 → converted 23,776, against 10,422 archived;
   - segments without target effects (p,p), (p,α₀), (p,α₁) agreed to 1–7 %, and the α₀ S factor to 1–3 % over
     E_c.m. = 5 keV–2.5 MeV;
   - the large differences were all `<targetInt>` segments (the target-effect rework) plus one (p,α₁) set about
     4 % low;
   - `--gsl-coul` did not close the gap (χ² 142,804);
   - the `.azr` had been edited four days after the fit, which confounds the comparison.

   Report the residual to the user; do not tune it away. A segment off by orders of magnitude (not percent) after
   conversion means the pairing is wrong: return to step 1.

## Related

- Reading `.sav` files by name in your own scripts, and remapping after a level or segment edit (a same-format remap):
  `references/parameters-and-sav.md`.
