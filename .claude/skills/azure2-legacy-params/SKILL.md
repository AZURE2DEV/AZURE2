---
name: azure2-legacy-params
description: Reuse an AZURE2 parameter file (param.sav, param.par, param.fit) written by a build older than August 2025 — detect the old j=<J>_la=<L>_energy / j=..._ch=..._rwa names, convert them to the current energy_<N> / width_<N>_<C> names, and verify the converted file against the model. Use before ANY mode 1/3/4/5 run, pyazr load_sav/full_rwa_from_sav call, or refit that starts from an archived fit's parameter file; whenever an archived chiSquared.out, rate or extrapolation will not reproduce; or when a run "used" a .sav but the result looks like the .azr's own levels.
---

# Old-format AZURE2 parameter files

## The trap

Commit `dcfdb34` (2025-08-01, "Added energy shifts…") renamed every level
parameter in `CNuc::FillMnParams`:

| before 2025-08-01 | now |
|---|---|
| `j=<J>_la=<L>_energy` | `energy_<N>` |
| `j=<J>_la=<L>_ch=<C>_rwa` | `width_<N>_<C>` |

`N` counts levels in J-group order (j ascending, then la within the group).
That is the same double loop that wrote the old names, so the mapping is
deterministic. `C` is unchanged, and so are the values (level energy, reduced-width
amplitude). `segment_<k>_norm` lines are named the same in both.

`AZUREParams::ReadUserParameters` matches a parameter file **by exact name**
and says nothing about names it does not recognise. An old file therefore
loads only its `segment_*` lines. **Every level stays at its `.azr` value**,
and the log still prints `Reading User Parameter File...`. Nothing in the
output says the fit was not applied.

pyazr has the same blind spot. `full_rwa_from_sav` (and `load_sav`, which uses
it) raises only when *no* name matches, and the `segment_*` lines always match.

**Scope (2026-10-14):** of 11,073 `*.sav`/`param.par`/`param.fit` files under
`~/R-matrix` (find depth 5), **6,997 are old format** and 2,651 are current.
Most archived fits before 2025-08 are affected.

**Observed (10B+p, 2026-10-14):** the 7-11-22 fit (`chiSquared.out` total
10,422) gave 133,516 with its own `param.sav` and today's binary. Its α₀
rate, rerun the same way, came out at the `.azr` levels' value (T9=1:
9.41e5) instead of the fit's (9.70e5). Details are in
`~/R-matrix/10B+p/10-14-26_reaction_rate/readme`.

## Detect

```bash
grep -v '^#' output/param.sav | head -1          # old format starts with  j=1_la=1_energy
```

Whole directory tree (prints the old-format files):

```bash
find . -name '*.sav' -o -name 'param.par' -o -name 'param.fit' | while read f; do
  case "$(grep -v '^#' "$f" | head -1 | awk '{print $1}')" in j=*) echo "$f";; esac; done
```

## Convert

```bash
python3 ~/AZURE2/.claude/skills/azure2-legacy-params/convert_sav_names.py \
    output/param.sav output/param_newnames.sav --check <fresh>/output/param.par
```

- Write to a **new** file. Archive directories are never overwritten in place.
  Put the converted file in the new dated fit directory.
- Value columns are copied byte for byte. Only the names change.
- Exit status: 0 = converted (and check passed), 1 = no old names (nothing
  written), 2 = check failed.
- The converted file carries no `#parametrization` tag, like the original. An
  untagged file read by a `--use-park` run is taken as Brune amplitudes (see
  azure2-eval). Run in the same formalism the original fit used, taken from the
  directory's job script (`run_crc_*`) or `readme`.

## Verify (all three, in order)

1. **Names and structure: `--check`.** `<fresh>/output/param.par` must be
   written by the *current* build from the *same* `.azr`. Make it in a scratch
   copy so the archive's `output/` is not touched:
   ```bash
   mkdir -p $SCR/output $SCR/checks && cp -rp model.azr data $SCR/ && cd $SCR
   printf '1\n\n\n\n7\n' | ~/bin/AZURE2 <the fit's flags> --no-gui --no-readline model.azr
   ```
   The log must say `Creating New param.par File...`. The check requires every
   level name to match both ways and every width that is exactly 0 in one file
   to be 0 in the other. Widths that are 0 in the `.azr` are fixed and cannot
   become nonzero in a fit. A pattern mismatch therefore means the levels or
   channels pair up differently: the `.azr` structure changed after the fit, or
   the file belongs to another model. A deliberate wrong pairing (10B+p
   9-8-16 `.sav` against the 7-11-22 model) fails with 257/331 names in common
   and dozens of pattern mismatches.

2. **χ² against the archived `chiSquared.out`.** Run mode 1 with the converted
   file (scratch copy again) and compare per segment. First check that the
   archived `param.sav`, `chiSquared.out` and `normalizations.out` share a
   timestamp (`ls --time-style=full-iso`): written together means they came from
   the same run. Also check that the `.azr` is not newer than them; an `.azr`
   edited after the fit cannot be expected to reproduce it.

3. **What to expect.** A perfect match is unlikely across a year or more of
   code changes. On 10B+p (2026-10-14, binary of 2026-10-13):
   - unconverted 133,516 → converted 23,776, against 10,422 in 2022;
   - segments without target effects (pp, pα₀, pα₁) agreed to 1–7%, and the
     α₀ S-factor to 1–3% over 5 keV–2.5 MeV;
   - the large differences were all `<targetInt>` segments (the 2026
     target-effect rework) plus one pα₁ set about 4% low;
   - `--gsl-coul` did not close the gap (χ² 142,804);
   - the `.azr` there had also been edited four days after the fit, which
     muddies the comparison.

   Report the residual to the user. Do not tune it away. A segment that is
   off by orders of magnitude (not percent) after conversion means the pairing
   is wrong: go back to step 1.

## Related

- Reading `.sav` files by name in your own scripts: azure2-eval, "Read `.sav`
  files by parameter NAME".
- Renumbering after adding or removing a level or a `<segmentsData>` line is a
  different, same-format remap: azure2-eval, "How the names are assigned".
