# EXFOR and NDS through `pyazr.nds`

This reference covers the two IAEA services that supply experimental data and level schemes for an evaluation, and
the conversion of their records into AZURE2 data files and segments. Everything is wrapped in `pyazr.nds`
(`from pyazr import nds as N`), which calls the web APIs directly: no local database, a network connection, no
dependency beyond NumPy. `gui/src/ExforData.cpp` is an independent client of the same Web API inside the Qt GUI; a
parsing rule learned in one belongs in both. Worked example: `pyazr/examples/exfor_fetch.py`. Command line:
`python -m pyazr.nds search --target C-13 --reaction p,g --quantity SIG`, `... download <DatasetID> -o file.dat`,
`... reference <DatasetID>`. API manual: `https://nds.iaea.org/exfor/x4guide/API/`.

## 1. Services

| Service | URL | Gives | wrapper |
|---|---|---|---|
| EXFOR search | `https://nds.iaea.org/exfor/x4list` | datasets (ID, reaction, points, E/θ range, author, reference) | `search_exfor()` |
| EXFOR data | `https://nds.iaea.org/exfor/x4get?DatasetID=…&op=csv` | one dataset with units | `fetch_exfor()` |
| EXFOR bibliography | `https://nds.iaea.org/exfor/x4get?sub=ENTRY` | title, authors, journal | `reference()` |
| EXFOR master files | clone of `https://github.com/IAEA-NDS/exfor_master` at `$AZURE2_EXFOR_MASTER` or `~/.azure2/exfor_master` | every entry as raw EXFOR (`exforall/<3 chars>/<entry>.x4`) | `fetch_exfor(source="local")`, `reference(source="local")` |
| LiveChart levels | `https://www-nds.iaea.org/relnsd/v1/data?fields=levels&nuclides=…` | level scheme | `fetch_levels()` |
| NNDC ENSDF API | `https://www.nndc.bnl.gov/ensdf-api/levels?nuclide=…&document-type=adopted` | the same adopted levels | `fetch_levels(source="nndc")` |
| LiveChart γ rays | `…?fields=gammas&nuclides=…` | transitions | `fetch_gammas()` |
| LiveChart ground state | `…?fields=ground_states&nuclides=…` | masses, J^π, separation energies | `fetch_ground_state()` |
| CrossRef | `https://api.crossref.org/works` | DOI | `resolve_doi()` |

## 2. Searching

`search_exfor(target, reaction, quantity, author, accnum, limit)`; give at least one filter (none returns an empty
list).

- `target` uses EXFOR names: `"C-13"`, `"12C"`, `"U-235"`. LiveChart's `nuclides=` uses lowercase without hyphen:
  `"14n"`.
- `reaction` is the REAC shorthand `projectile,exit`: `"p,g"`, `"p,el"`, `"n,f"`, `"d,p"`.
- `quantity`: `SIG` (integrated; `…,,SIG,,SFC` for an S-factor capture set), `DA` (differential), `pol` (analysing
  power, `…,,POL/DA`; `AP` returns nothing), `FY` (yields).
- Results carry `dataset_id`, `reaction`, `npoints`, `en_min`/`en_max` (eV, **lab**), `an_min`/`an_max` (deg),
  `author`, `reference`, `year`.
- A target the search cannot name is often reached by `accnum=` or by searching the product (reaction `…,g` plus a
  broad target).

## 3. Fetching

`fetch_exfor(dataset_id, plus=0)` → `ExforData` in EXFOR native units (what the conversion needs; `plus=1`
computational units, `plus=2` the universal grid). It exposes `reaction`, `year`, `author`, `projectile`, `target`,
`exit`, `arrays` (`DATA`, `ERR-S`, `ERR-SYS`, `EN-CM`/`EN`, `ANG-CM`/`ANG`) and `columns`.

- S-factor capture (`SIG,,SFC`): `DATA` in `B*KEV` or `B*EV`, energy `EN-CM`, `ERR-SYS` in per cent.
- Differential (`DA`): `DATA-CM (NB/SR)`, `EN (MEV)` lab, `ANG-CM (ADEG)`.
- Analysing power (`POL/DA`): `DATA (NO-DIM)`, error column `DATA-ERR`, `ANG-CM`.
- Ratio to Rutherford (`DA,,RTH`): `DATA (NO-DIM)` with an angle column, indistinguishable from an analysing power
  once fetched, because `x4get` drops the `,,RTH` suffix that `x4list` reports. The search result's reaction code is
  the only place the marker survives.
- Units are read by prefix (`B`, `MB`, `MICRO-B`, `MU-B`, `UB/SR`, `NB/SR`, `KB`, `B*KEV`, …); an unrecognised unit
  raises instead of passing through.
- The S-factor unit does not name the abscissa's unit: `B*KEV` is a unit of the product barn × energy. The declared
  `EN`/`EN-CM` unit wins and a disagreement is a warning.
- `x4list` energies are lab; `x4get` energies are what the column says (usually `EN-CM`). The two differ by the
  kinematic factor, a useful consistency check.

## 4. Converting to AZURE2

`ExforData.to_azr(data_dir, entrance, exit, observable=None, …)` writes a data file (E_lab, θ, σ, δσ) and returns
keyword arguments for `AzrModel.add_data_segment(**kw)`:

```python
d = N.fetch_exfor("O2599004")
kw = d.to_azr("run/data", entrance=1, exit=2, observable="total-capture")
AzrModel.from_file("13N.azr").add_data_segment(**kw).write("13N_new.azr")
```

- `EN-CM` → lab by E_lab = E_c.m.(m_t + m_p)/m_t; `EN` is passed through as lab.
- S factor → σ with `sfactor_to_cross_section`, using AZURE2's own Sommerfeld constant (0.157488), so the barns
  reproduce the reported S exactly.
- Angles are written unchanged; choose `observable` to match their frame:

| EXFOR gives | observable= | AZURE2 code |
|---|---|---|
| capture, `EN-CM` | `"total-capture"` (exit is a γ pair) | 3 |
| `EN-CM`, no angle | `"angle-integrated"` | 0 |
| `ANG-CM` | `"differential-cm"` | 4 |
| lab `ANG` | `"differential"` | 1 |
| `POL`/`NO-DIM` | `"analyzing-power"` (angle c.m.) | 7 |

  Without `observable` it guesses from the reaction code (capture → total-capture, `POL` → analyzing-power, `ANG-CM`
  → differential-cm, lab `ANG` → differential). A c.m.-angle set given as `differential` is converted a second time.
- `energy_min/max` (lab MeV) and `angle_min/max` are derived from the points; `ERR-SYS` in per cent goes into
  `norm_error`. `cross_section_scale`, `energy_scale`, `angle_scale` are the manual fix-ups.
- Uncertainties come from `ERR-S`, `DATA-ERR` or `ERR-T`; a `PER-CENT` column is applied as a percentage. A point
  without an error gets 5 % of its value (`_DEFAULT_REL_ERROR`), since zero would be an infinitely precise point.
  Report every data set that triggered this. An `ERR-SYS` given in absolute units instead of per cent, and
  non-positive values, need a manual check.
- **Ratio to Rutherford**: `rutherford=True` multiplies by the Coulomb cross section (→ b/sr), `rutherford=False`
  accepts the ratio. For a lab `ANG` column `to_azr` converts the angle to c.m. (`elastic_lab_to_cm`, entrance
  masses; pass `target_mass`/`projectile_mass` for exact ones) and divides the c.m. Rutherford value by the lab→c.m.
  Jacobian, so AZURE2 reads back ratio × σ_R,c.m.(θ_c.m.). Identical particles (ratio to Mott) get a warning and
  ratio × Rutherford. Test: `tests/pyazr/nds_rutherford_test.py`; the GUI dialog: `tests/gui/exfor_rutherford_test.cpp`.
- builds before 2026-10-27 (expires: 2027-02-01): `to_azr` and the GUI EXFOR dialog multiplied a `,,RTH` ratio by
  the c.m. Rutherford cross section at the **lab** angle (1.13× off for p+¹⁶O, 1.8–2.2× for α+¹²C) and called every
  angle column c.m. Regenerate data files made from a lab-angle `,,RTH` set with such a build.
- **The light/heavy assignment.** AZURE2 converts lab to c.m. with E_c.m. = E_lab·m₂/(m₁ + m₂), m₁ the light member
  (`src/EPoint.cpp`), and its angle transform asin(m₁/m₂ sin θ) needs m₁ < m₂. EXFOR names target and projectile as
  the experiment ran: ³He(α,γ)⁷Be is `2-HE-3(A,G)4-BE-7`, target ³He. Letting `to_azr` infer the masses there writes
  E_lab too large by 4/3. State the assignment on every fetch, with atomic masses from `pyazr.masses`:

  ```python
  from pyazr.masses import atomic_mass
  d.to_azr(data_dir, entrance=1, exit=-1, observable="total-capture",
           projectile_mass=atomic_mass("3He"),   # AZURE2's m1
           target_mass=atomic_mass("4He"))       # AZURE2's m2
  ```

  The masses do not change the S → σ conversion (μ is symmetric), only the energy axis.
- **Charges of light projectiles.** `za()` resolves the `A`, `HE3`, `D`, `T` shorthands, and the S-factor path
  raises rather than assume a charge. σ = S exp(−2πη)/E with η ∝ Z₁Z₂, so a wrong charge is a factor ≈ 5000 at
  E_c.m. = 93 keV for ³He+α.
- **Cross-check every converted data set against an independent one.** A capture measurement compiled twice, as
  `SIG` and as `SIG,,SFC`, should agree after conversion (0.3 % for the LUNA ³He(α,γ)⁷Be pair).
- Editing data segments makes the cached EC integrals stale: delete `output/intEC.dat` and `output/intEC.extrap`
  before running the edited model (azure2-calculate).

## 5. Level schemes, γ rays, ground states

`fetch_levels("14n")` → `Level` objects: `energy` (keV), `energy_mev`, `jp`, `spin`, `parity` (±1), `half_life`;
from the NNDC fallback also `width_kev`; `source` names the service. This adopted ENSDF scheme seeds the J-groups.
`fetch_gammas` gives transitions (E_γ keV, multipolarity, intensity) for decay-scheme checks. `fetch_ground_state`
gives J^π, mass excess, S_n and S_p (keV), radius.

```python
for lv in N.fetch_levels("14n"):
    print(lv.energy_mev, lv.jp, lv.half_life)
gs = N.fetch_ground_state("14n")   # S_n = 10553 keV, S_p = 7550 keV, J^pi = 1+
```

Masses and separation energies for the `.azr` come from the shipped AME2020 table (`pyazr.masses`), not from
LiveChart (SKILL.md §2).

## 6. The paper behind a data set

The bibliography is per entry, the first five characters of the dataset ID (`O2599004` → `O2599`).

```python
ref = N.reference("O2599004")      # .title .authors .reference .x4ref ("J,PRL,131,162701,2023") .year
doi = N.resolve_doi(ref)           # "10.1103/physrevlett.131.162701"; open https://doi.org/<doi>
```

`resolve_doi` tries the `x4ref` journal code against known DOI patterns (APS: `aps_doi("PRL", "131", "162701")`),
then a CrossRef title search. The `reference()` string can always be cited by hand; arXiv often has a free copy.

## 7. Service failures

- **Cloudflare challenge**: HTTP 403 (or 200) with a "Just a moment..." page, intermittently, for EXFOR, the data
  explorer API and LiveChart. `fetch_levels` falls back to the NNDC ENSDF API with a `UserWarning` (same adopted
  levels in LiveChart's spellings; only trailing zeros of half-life text and `half_life_sec` of width-only levels can
  differ). Force a service with `source="livechart"` or `"nndc"`.
- `fetch_exfor` and `reference` fall back to the local master files (`$AZURE2_EXFOR_MASTER` or
  `~/.azure2/exfor_master`; refresh with `git pull origin main` in the clone). The fallback warns; `source="local"`
  skips the web. It gives `plus=0` only, needs the pointer for a multi-dataset subentry (`"C1010003.2"`; the error
  lists them), and carries COMMON values as columns, as `x4get` does. On 439 entries it matched the data explorer's
  parse: 1,596/1,596 tables and 1,663/1,663 reaction codes.
- `search_exfor`, `fetch_gammas` and `fetch_ground_state` have no fallback: search the mirror by hand
  (`grep -l "6-C-13(P,G)" exforall/*/*.x4`) or retry later.
- `x4list` answering with the HTML "Select" page means a malformed request (e.g. `quantity="AP"`): `search_exfor`
  raises with the HTML text.
- NNDC also serves ENDF evaluations at `https://www.nndc.bnl.gov/endf-api/` (azure2-exchange); its listed SHA-1 of the
  GNDS `.xml` files did not match the served files, the ENDF-6 ones did.
