---
name: nds-explorer
description: Navigate the IAEA EXFOR database and NDS nuclear-data services to find experimental data — level structures, integrated cross sections, differential cross sections, analyzing powers, yields, gamma transitions — and convert them into the formats AZURE2 expects, including resolving the bibliographic reference behind every dataset. Use whenever the task involves finding published cross-section/scattering data, level schemes, or reaction data for an R-matrix evaluation, or fetching the paper behind an EXFOR entry. Works entirely through the pyazr.nds module.
---

# Navigating EXFOR and NDS for nuclear data

This skill covers the two IAEA services that supply the experimental data for
R-matrix evaluations and how to get that data into AZURE2. Everything is
wrapped in **`pyazr.nds`** (imported as `from pyazr import nds`), which talks
straight to the web APIs — no local database, needs a network connection, and
no dependency beyond NumPy and the standard library.

`pyazr.nds` is the *only* EXFOR client in the package. `gui/src/ExforData.cpp`
is its counterpart inside the Qt setup utility; the two are independent
implementations of the same Web-API, so a parsing rule learned by either
belongs in both.

- **EXFOR** — the experimental reaction database: cross sections, differential
  cross sections, analyzing powers, yields, etc., one entry per measurement,
  with the full bibliography.
- **LiveChart / ENSDF** — evaluated level schemes, gamma transitions and
  ground-state properties.

For the AZURE2 side of converting fetched data into segments, levels and `.azr`
edits, load the **`azure2-eval`** skill as well — this skill hands data to it.

## The two services at a glance

| Service | Base URL | What it gives | pyazr.nds wrapper |
|---|---|---|---|
| EXFOR search | `https://nds.iaea.org/exfor/x4list` | matching datasets (ID, reaction, points, E/angle range, author, reference) | `search_exfor()` |
| EXFOR data | `https://nds.iaea.org/exfor/x4get?DatasetID=…&op=csv` | one dataset's numbers with units | `fetch_exfor()` |
| EXFOR bibliography | `https://nds.iaea.org/exfor/x4get?sub=ENTRY` | title, authors, journal reference | `reference()` |
| EXFOR master files (local fallback) | a clone of `https://github.com/IAEA-NDS/exfor_master` at `$AZURE2_EXFOR_MASTER` or `~/.azure2/exfor_master` | every entry as raw EXFOR (`exforall/<3 chars>/<entry>.x4`) | `fetch_exfor(source="local")`, `reference(source="local")` |
| LiveChart levels | `https://www-nds.iaea.org/relnsd/v1/data?fields=levels&nuclides=…` | level scheme CSV | `fetch_levels()` |
| NNDC ENSDF API (fallback) | `https://www.nndc.bnl.gov/ensdf-api/levels?nuclide=…&document-type=adopted` | the same adopted levels, JSON | `fetch_levels(source="nndc")` |
| LiveChart gammas | `…?fields=gammas&nuclides=…` | gamma transitions | `fetch_gammas()` |
| LiveChart ground state | `…?fields=ground_states&nuclides=…` | masses, J^pi, separations | `fetch_ground_state()` |
| CrossRef | `https://api.crossref.org/works` | DOI for a paper | `resolve_doi()` |

## Searching EXFOR

`search_exfor(target, reaction, quantity, author, accnum, limit)` — all
filters optional, give at least one. Naming rules that will bite you:

- **`target` uses EXFOR names**: `"C-13"`, `"12C"`, `"U-235"` (hyphen, no
  isotope before the symbol for Z>~2). LiveChart's `nuclides=` uses the
  *opposite* convention, lowercase without hyphen: `"14n"`.
- **`reaction`** is `projectile,exit`: `"p,g"`, `"p,el"`, `"n,f"`, `"d,p"`.
  (this is the EXFOR *REAC* shorthand, not the full `6-C-13(P,G)7-N-14` code).
- **`quantity`** codes that actually work (all verified live):
  | want | quantity= | reaction code you'll see |
  |---|---|---|
  | total cross section | `SIG` | `…,,SIG,,SFC` for an S-factor capture set |
  | differential (angle-resolved) | `DA` | `6-C-13(P,EL)6-C-13,,DA` |
  | analyzing power | `pol` | `…,,POL/DA` (note: `AP` returns **nothing**) |
  | yields (fission etc.) | `FY` | `(N,F),,FY` |
  | elastic | — | reaction `p,el` |
- Result objects carry `dataset_id`, `reaction`, `npoints`, `en_min`/`en_max`
  (eV), `an_min`/`an_max` (deg, 0 if not angular), `author`, `reference`,
  `year`. The `reference` field is a human citation; the DOI needs `reference()`
  + `resolve_doi()`.

```python
from pyazr import nds as N
hits = N.search_exfor(target="C-13", reaction="p,g", quantity="SIG", limit=10)
for h in hits: print(h)   # <ExforDataset O2599004 6-C-13(P,G)7-N-14,,SIG,,SFC n=31 …>
```

A target the EXFOR search can't name (exotic, or the search is too broad)
can often still be reached by `accnum=` if you know an entry number, or by
searching the *product* via reaction `…,g` plus a broad target.

## Fetching a dataset

`fetch_exfor(dataset_id, plus=0)` → `ExforData`. `plus=0` keeps the EXFOR
native units, which is what you want for AZURE2 conversion (it reads units off
the CSV header):

- **S-factor capture** (`SIG,,SFC`): `DATA` in `B*KEV` (or `B*EV`), energy as
  `EN-CM (KEV)`, systematic error `ERR-SYS (PER-CENT)`. AZURE2 needs barns →
  `to_azr` converts via the Sommerfeld factor automatically.
- **Differential** (`DA`): `DATA-CM (NB/SR)`, `EN (MEV)` (lab), `ANG-CM (ADEG)`.
  `to_azr` scales `NB/SR → b/sr` (÷1e9) and passes angles through.
- **Analyzing power** (`POL/DA`): `DATA (NO-DIM)`, error column is `DATA-ERR`
  (not `ERR-S`), `ANG-CM (ADEG)`. Dimensionless; `to_azr` leaves values alone.
- **Ratio to Rutherford** (`DA,,RTH`): also `DATA (NO-DIM)` with an angle
  column — *indistinguishable from an analyzing power once fetched*, because
  `x4get` drops the `,,RTH` suffix that `x4list` reports. `to_azr` warns and
  passes the values through; pass `rutherford=True` to multiply by the Coulomb
  cross section (→ b/sr), or `rutherford=False` to accept the ratio and
  silence the warning. **Check the search result's reaction code** — that is
  the only place the marker survives. The cross section written is in the
  frame of the angle column: for a lab `ANG` column `to_azr` converts the
  angle to the c.m. frame (`elastic_lab_to_cm`, the entrance masses — pass
  `target_mass`/`projectile_mass` for exact ones) and divides the c.m.
  Rutherford value by the lab→c.m. Jacobian, so AZURE2 reads back
  ratio × σ_R,cm(θ_cm). **Before 2026-10-26 (audit item 8) it wrote the c.m.
  Rutherford at the lab angle** (1.13× off for p+16O, 1.8–2.2× for α+12C) and
  `_guess_observable` called every angle column `differential-cm`; a lab
  column is `differential` now, and a ratio is no longer guessed as an
  analyzing power. Any data file made from a `,,RTH` set with a lab `ANG`
  column before that date must be regenerated (`tests/pyazr/nds_rutherford_test.py`).
  Identical particles (ratio to Mott) get a warning and ratio × Rutherford.
  The GUI's own EXFOR dialog (gui/src/ExforData.cpp, "Retrieve Experimental Data
  from EXFOR") had the same lab-angle mistake and is fixed the same way (2026-10-27,
  tests/gui/exfor_rutherford_test); data files it made from a lab-angle `,,RTH`
  set before that date (12C+p/8-9-26_claude_learns_12C+p: Milne 1954, Jiang 2004,
  if their EXFOR angle column is `ANG`) need regenerating too.
- `plus=1` gives computational units (eV, `B*EV`, `B/SR`); `plus=2` the
  universal grid. Not needed for AZURE2; `plus=0` is the default.

Units are read by prefix, not from a fixed table, so the many spellings EXFOR
uses for the same quantity (`B`, `MB`, `MICRO-B`, `MU-B`, `UB/SR`, `NB/SR`,
`KB`, `B*KEV`, …) all convert. An unrecognised unit no longer falls through to
"assume it is already right", which was silently wrong by powers of ten.

The parsed object exposes `reaction`, `year`, `author`, `projectile`,
`target`, `exit`, `arrays` (the numeric columns keyed by plain name: `DATA`,
`ERR-S`, `ERR-SYS`, `EN-CM`/`EN`, `ANG-CM`/`ANG`) and `columns` (name + unit).

## Converting to AZURE2 (the important part)

`ExforData.to_azr(data_dir, entrance, exit, observable=None, …)` writes an
AZURE2 data file (`energy angle crossSection error`, **lab frame**) and returns
the exact keyword arguments for `AzrModel.add_data_segment(**kw)`. Frames are
handled for you:

- **`EN-CM` → lab** via `E_lab = E_cm·(m_t + m_p)/m_t` (masses come from the
  reaction code; override with `target_mass=`/`projectile_mass=`).
- **`EN` is assumed lab** already — passed through.
- **S-factor (`B*KEV`/`B*EV`) → barns** with `sfactor_to_cross_section`, using
  AZURE2's own Sommerfeld constant (`0.157488`), so feeding the resulting barns
  reproduces the reported S-factor exactly (verified against EPoint.cpp).
- Angle is **written through unchanged**. Pick the `observable` to match its
  frame (see below). `cross_section_scale`/`energy_scale`/`angle_scale` are the
  escape hatches for manual fix-ups.

Choose `observable` from the data's frame:

| EXFOR gives | use observable= | AZURE2 code |
|---|---|---|
| `EN-CM` + `ANG-CM`, capture | `"total-capture"` (exit is a gamma) | 3 |
| `EN-CM`, no angle | `"angle-integrated"` | 0 |
| `ANG-CM` angle | `"differential-cm"` (CM frame in the file) | 4 |
| lab `ANG` angle | `"differential"` | 1 |
| `POL`/`NO-DIM` | `"analyzing-power"` (angle is CM) | 7 |

`to_azr` also derives `energy_min/max` (lab MeV) and `angle_min/max` from the
converted points and folds the EXFOR `ERR-SYS` percent into `norm_error` when
present. If `observable` is omitted it guesses from the reaction code
(capture → `total-capture`, `POL` → `analyzing-power`, angular → `differential-cm`).

**Uncertainties.** The per-point error is taken from `ERR-S`, `DATA-ERR` or
`ERR-T`, whichever is present, and a column marked `PER-CENT` is applied as a
percentage of the value rather than as an absolute one — getting that backwards
turns "6 per cent" into "6 millibarn/sr". A point EXFOR gives no error for is
assigned 5% of its value, because a zero error reaches AZURE2 as an infinitely
precise measurement and would dominate the χ².

```python
d = N.fetch_exfor("O2599004")
kw = d.to_azr("run/data", entrance=1, exit=2, observable="total-capture")
AzrModel.from_file("13N.azr").add_data_segment(**kw).write("13N_new.azr")
```

**The usual segment-editing rule applies**: adding/removing data segments
changes the integration grid, so delete `output/intEC.dat` and
`output/intEC.extrap` before running the edited model (see azure2-eval).

## Level structure, gammas, ground state

`fetch_levels("14n")` → list of `Level` (`energy` keV, `.energy_mev`, `jp`,
`.spin` float, `.parity` ±1, `half_life`; from the NNDC fallback also `width_kev`, and `source` says which
service answered). This is the ENSDF-evaluated scheme
and is what you seed an R-matrix J-group structure from. `fetch_gammas("14n")`
gives transitions (energy keV, multipolarity, intensity) for decay-scheme
cross-checks. `fetch_ground_state("14n")` gives J^pi, spin/parity, mass excess,
neutron/proton separation energies (Sn/Sp, keV) and radius — useful for
checking thresholds and pair energies when building the `.azr`.

```python
for lv in N.fetch_levels("14n"):
    print(lv.energy_mev, lv.jp, lv.half_life)
gs = N.fetch_ground_state("14n")   # -> Sn=10553 keV, Sp=7550 keV, J^pi=1+
```

## Getting the paper behind a dataset

EXFOR's bibliography is per **entry**, which is the first 5 characters of the
dataset ID (`O2599004` → entry `O2599`).

```python
ref = N.reference("O2599004")
# ref.title, ref.authors, ref.reference ("Jour: … (2023)"), ref.x4ref ("J,PRL,131,162701,2023"), ref.year
doi = N.resolve_doi(ref)           # -> "10.1103/physrevlett.131.162701"
# open the paper:
#   https://doi.org/<doi>
```

`resolve_doi` first tries the structured `x4ref` journal code against known
DOI patterns (APS journals — `aps_doi("PRL","131","162701")` is
`10.1103/physrevlett.131.162701`, verified) then falls back to a CrossRef
title search. For non-APS journals give it the title and year and it does its
best; the `reference()` string is always there to cite manually. Unpaywalled
copies often live on arXiv — searching the title there is a good follow-up.

## Frames: EXFOR's target is not always AZURE2's target

AZURE2 converts lab to c.m. with `E_cm = E_lab · m2/(m1+m2)`, where `m1` is the
**light** member of the pair (`EPoint.cpp:622`), and its angle transform
`asin(m1/m2 · sin θ)` is only defined for `m1 < m2`. EXFOR names target and
projectile by how the experiment was run, which is often the other way round —
`3He(α,γ)7Be` is compiled as `2-HE-3(A,G)4-BE-7`, target ³He, projectile α.
Letting `to_azr` infer the masses there writes lab energies too large by
`(3+4)/3 ÷ (3+4)/4 = 4/3`.

**State the pair's light/heavy assignment on every fetch**, not only when it
looks reversed:

```python
d.to_azr(data_dir, entrance=1, exit=-1, observable="total-capture",
         projectile_mass=nuclear_mass("3He"),   # AZURE2's m1
         target_mass=nuclear_mass("4He"))       # AZURE2's m2
```

The masses do not change the S-factor→barn conversion (the reduced mass is
symmetric), only the energy axis. Getting it wrong displaces a dataset by tens
of per cent in energy and shows up as a resonance that will not sit still.

## Two conversion bugs fixed here (August 2026)

Both were silent and both changed data by orders of magnitude. They are fixed
in `pyazr/nds.py`; the reasoning is recorded because the same traps exist in
any hand-rolled EXFOR reader, including `gui/src/ExforData.cpp`.

**1. The S-factor unit does not name the abscissa's unit.** `_energy_column`
used to override the declared `EN`/`EN-CM` unit whenever the `DATA` unit was an
S-factor unit with a different energy part — reading `B*KEV` as "the energies
are in keV". But `B*KEV` is a unit of the *product* barn × energy, chosen so
the tabulated S values read well, and it says nothing about the energy column.
Checked against the independent lab ranges `x4list` reports, the declared
column unit was right and the S-factor unit different in every case that fired:
`C1610002` (Brown 2007, `EN-CM/MEV` with `DATA/B*KEV`), `O2521002` (Piatti
2020, `EN-CM/KEV` with `DATA/B*EV`), `O1590003` (Cruz 2008, `EN/KEV` with
`DATA/B*MEV`). The declared unit now wins and the disagreement is a warning.

**2. `za()` could not read a light-particle projectile, and the charge fell
back silently.** `za()` matched only `Z-SYM-A` codes, so the `A`, `HE3`, `D`,
`T` shorthands EXFOR writes inside the reaction parentheses returned `None`,
and `_data_column` then assumed `z1 = 1`, `z2 = 6` — p + ¹²C. Since
`σ = S·exp(−2πη)/E` with `η ∝ Z₁Z₂`, halving `Z₁Z₂` for ³He+α multiplies the
cross section by **≈5000** at 93 keV. The LUNA `3He(α,γ)7Be` S-factor set
`O15160021` came out 5000× above the `O17090021` cross sections from the same
measurement; after the fix the two agree to 0.3%. `za()` now resolves the
shorthands and the S-factor path raises rather than assuming a charge.

**The lesson generalises: cross-check every converted dataset against an
independent one.** A capture measurement compiled twice, once as `SIG` and once
as `SIG,,SFC`, is the cheapest possible check and catches both bugs at once.

## Common failure modes

- **IAEA answers with a Cloudflare challenge** — HTTP 403 with a "Just a moment..." HTML page (or that page
  with HTTP 200), intermittently, for EXFOR (`x4list`/`x4get`), the dataexplorer API and LiveChart alike.
  Seen from the Notre Dame CRC cluster on 2026-10-15: LiveChart answered one request and refused the next
  ones for minutes.  `fetch_levels` falls back to the NNDC ENSDF API by itself (a `UserWarning` says so):
  same ENSDF adopted levels, checked against LiveChart on 14N and 16O (same levels, energies and
  uncertainties), written in LiveChart's spellings (`jp` "2+ 1+", `half_life` the bare number), plus
  `Level.width_kev` and `Level.source`.  Only the trailing zeros of a half-life text (LiveChart "4.70",
  NNDC "4.7") and `half_life_sec` at ~1e-7 for width-only levels can differ.  Force a service with
  `source="livechart"` / `"nndc"`.  `fetch_exfor` and `reference` fall back to a local clone of the
  EXFOR master files (`$AZURE2_EXFOR_MASTER` or `~/.azure2/exfor_master`; on the CRC cluster
  `/groups/rdeboer1/user/rdeboer1/exfor_master`, linked from `~/.azure2/exfor_master`, cloned 2026-10-22,
  refresh with `git pull origin main`).  The fallback warns; `source="local"` skips the web.  It gives
  plus=0 only, needs the pointer for a multi-dataset subentry (`"C1010003.2"`, the error lists them), and
  carries the COMMON values (an angle, a per-cent error) as columns as x4get does.  Checked against the
  dataexplorer API's own parse of 439 entries: 1,596/1,596 tables and 1,663/1,663 reaction codes identical.
  On 2026-10-22 x4get AND the dataexplorer API were both behind the challenge, so the dataexplorer
  work-around no longer helps there.  `search_exfor`, `fetch_gammas` and `fetch_ground_state` have no
  fallback (search the mirror by hand: `grep -l "6-C-13(P,G)" exforall/*/*.x4`) -- retry later.  NNDC also
  serves ENDF evaluations at `https://www.nndc.bnl.gov/endf-api/` (pointwise cross sections; its listed
  SHA-1 for the main GNDS .xml files does not match the served files, the ENDF-6 ones do; 2026-10-15).

- **`x4list` returns the HTML "Select" page, not JSON** — the request was
  malformed (e.g. `quantity="AP"` instead of `pol`). `search_exfor` raises with
  the HTML text; re-check the quantity/reaction codes.
- **`x4list` energies are LAB, `x4get` energies are whatever the column says**
  — usually `EN-CM`. The search result's `en_min`/`en_max` and the fetched
  `EN-CM` array therefore differ by the kinematic factor, which is not a bug
  but is a good consistency check on a fetch.
- **Missing uncertainties are invented.** A point with no error column gets 5%
  of its value (`_DEFAULT_REL_ERROR`), because zero would reach AZURE2 as an
  infinitely precise measurement. Datasets that trigger this must be reported
  as such — see the `diagnose()` helper in
  `evaluations/_tools/exfor_to_azr.py`, which also flags an `ERR-SYS` given in
  absolute units instead of per cent (`C2032002`, Kontos 2013) and points with
  non-positive values.
- **Wrong frame silently** — mixing a CM-angle dataset with
  `observable="differential"` makes AZURE2 convert the angle a second time.
  Use `differential-cm` for `ANG-CM` EXFOR data, `differential` for lab `ANG`.
- **S-factor vs barns** — a `SIG,,SFC` dataset fed as if it were already barns
  is wrong by many orders of magnitude. `to_azr` handles it; hand-editing does
  not (see `sfactor_to_cross_section`).
- **Data edited ⇒ stale EC cache** — delete `output/intEC.dat`/`.extrap` or
  the run silently reuses integrals for the old grids.
- **`search_exfor` with no filters** — returns an empty list; give at least one
  filter to narrow the hit set.

## See also

- `pyazr/examples/exfor_fetch.py` — the whole flow end to end.
- `pyazr/nds.py` — source of the wrappers (unit rules, DOI logic). Also a
  command line for a quick look: `python -m pyazr.nds search --target C-13
  --reaction p,g --quantity SIG`, `… download <DatasetID> -o file.dat`,
  `… reference <DatasetID>`.
- `gui/src/ExforData.cpp` — the GUI's own client for the same service.
- EXFOR Web-API manual: `https://nds.iaea.org/exfor/x4guide/API/`.
- For turning fetched data into fitted resonances, levels, χ²: load
  `azure2-eval`.
