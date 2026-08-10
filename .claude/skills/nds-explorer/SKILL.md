---
name: nds-explorer
description: Navigate the IAEA EXFOR database and NDS nuclear-data services to find experimental data — level structures, integrated cross sections, differential cross sections, analyzing powers, yields, gamma transitions — and convert them into the formats AZURE2 expects, including resolving the bibliographic reference behind every dataset. Use whenever the task involves finding published cross-section/scattering data, level schemes, or reaction data for an R-matrix evaluation, or fetching the paper behind an EXFOR entry. Works entirely through the pyazr.nds module.
---

# Navigating EXFOR and NDS for nuclear data

This skill covers the two IAEA services that supply the experimental data for
R-matrix evaluations and how to get that data into AZURE2. Everything is
wrapped in **`pyazr.nds`** (imported as `from pyazr import nds`), which talks
straight to the web APIs — no local database, needs a network connection.

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
| LiveChart levels | `https://www-nds.iaea.org/relnsd/v1/data?fields=levels&nuclides=…` | level scheme CSV | `fetch_levels()` |
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
- `plus=1` gives computational units (eV, `B*EV`, `B/SR`); `plus=2` the
  universal grid. Not needed for AZURE2; `plus=0` is the default.

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
`.spin` float, `.parity` ±1, `half_life`). This is the ENSDF-evaluated scheme
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

## Common failure modes

- **`x4list` returns the HTML "Select" page, not JSON** — the request was
  malformed (e.g. `quantity="AP"` instead of `pol`). `search_exfor` raises with
  the HTML text; re-check the quantity/reaction codes.
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
- `pyazr/nds.py` — source of the wrappers (units tables, DOI logic).
- EXFOR Web-API manual: `https://nds.iaea.org/exfor/x4guide/API/`.
- For turning fetched data into fitted resonances, levels, χ²: load
  `azure2-eval`.
