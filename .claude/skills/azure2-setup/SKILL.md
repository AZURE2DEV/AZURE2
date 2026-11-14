---
name: azure2-setup
description: Build or audit an AZURE2 project (format 3) before any fit - compound nucleus and thresholds, particle and photon pairs with AME masses and separation energies, levels and channels from ENSDF, channel radii, data sets as records (which quantity, frame, units, totals, ratios, natural targets, model data), EXFOR data (pyazr.sources.nds), conventions and the first checks (validate, describe, audit script). Use for a new reaction, a new data set, or an inherited project. Calculating is azure2-calculate; fitting azure2-fit; target effects azure2-targets; THM azure2-thm.
---

# Setting up an AZURE2 project

An evaluation starts from one `.azr` file (project format 3): `<config>` (options as `key=value`), `<levels>` (the
compound nucleus: particle pairs with masses m₁, m₂ in u and separation energies S in MeV; levels with J^π and E_x in
MeV; channels with pair, l, s, the channel radius a_c in fm and a width: Γ in eV for an open channel, an ANC in
fm^-1/2 for a closed one), and records: data sets, normalizations, energy scales, extrapolations, targets, beams,
detectors, model data, potentials, γ-decay records. Every later number depends on these inputs, and several setup
errors load without a message. This skill gives the order of building a project, the conventions, and the checks.
The grammar of every record and field: `docs/source/reference/format3.rst` (`AZURE2 --schema` prints the table).

## When to use, and when not

- A new compound nucleus or reaction, a new data set, EXFOR and ENSDF lookups, the `<levels>` block and the records,
  auditing an inherited project. An old-format file is refused (the message names `AZURE2 --convert OLD NEW`, or
  `pyazr.convert_azr(old, new)`): convert, then audit the conversion's notes
  (`docs/source/user_guide/older_projects.rst`).
- Calculations, extrapolations, rates: azure2-calculate. Fits, editing a fitted model, norms and shifts: azure2-fit.
  Trust: azure2-check. Data covariance files: azure2-uncertainty. Target, beam, detector: azure2-targets. THM:
  azure2-thm. Other codes and bases: azure2-exchange.

## Scripts

| Script | Does |
|---|---|
| `scripts/audit_project.py` | validation, pairs against AME2020, radii, data sets and rows kept, options; flags |
| `scripts/selftest.py` | the audit on `tests/upgrade/T1` and on copies with an S error and a mass number |

`python3 .claude/skills/azure2-setup/scripts/audit_project.py model.azr` — no flag is the start of a fit, not a proof
of a correct model.

## Procedure

1. **The compound nucleus and its thresholds.** Every partition open in the energy range of the data, and one photon
   pair per final state reached by capture. Level schemes and ground states: `pyazr.sources.nds`
   (references/nds-exfor.md).
2. **Pairs: masses and S from one table.** Atomic masses (u, electrons included) from AME2020 and S = m₁ + m₂ −
   M(compound) from the same table: `pyazr.sources.masses.fill_masses(project)`, then the audit (§ "Masses").
3. **Channels and radii.** The allowed (l, s) for each J^π (references/azr-levels.md §5;
   `Project.allowed_channels(J, parity, l_max=)` lists them by the Qt rule, `Project.add_pair(...)` adds a pair). Start a_c at
   1.4 (A₁^⅓ + A₂^⅓) fm; it is a model choice to scan with refits (azure2-check, `radius_scan.py`).
4. **Levels.** Known resonances in and near the data range, sub-threshold states whose tails reach it, background
   poles above (references/rmatrix-picture.md §3–4). `Project.add_level(J, parity, E_x, [dict(pair=, l=, s=,
   width=)])`; the width is observed Γ (eV) or an ANC, never a reduced-width amplitude. Every external-capture final
   state is a level. `includeAllJPi=1` in `<config>` for elastic data.
5. **Data.** One `dataset` record per measured quantity (references/data-kinds.md: the quantity, the columns with
   their frame and unit, totals and sums, natural targets, provenance). From EXFOR: references/nds-exfor.md. Each data
   set's error model (statistical, common, shape) and covariance file: azure2-uncertainty.
6. **Conventions.** Non-relativistic or relativistic kinematics (§ "Kinematics"), once, before the first fit; the
   options (`use_brune`, `use_long_wavelength`, ...) as `<config>` keys so that every front end runs the same model
   (`Project.set_config(use_long_wavelength=False)`). Keep the project in its own directory with `output/` and
   `checks/` beside it.
7. **First checks** (minutes, no fit):
   - `scripts/audit_project.py model.azr` (or `AZURE2 --validate model.azr`, exit 2 with an error);
   - `AZURE2 --describe model.azr`: each data set's observable, reaction tree, rows kept; `--print-predictions`;
   - a calculation at the file's own parameters (azure2-calculate), each data set plotted against it in the same
     frame (azure2-plot).

## Rules and pitfalls

### Masses and separation energies

- **Particle masses are atomic masses in u, not mass numbers.** The mass columns of `<levels>` are doubles in u;
  the reduced mass, lab↔c.m. conversion, penetrabilities and phase space all use it. An integer mass (`4` for an α
  instead of 4.002602) degrades every kinematic quantity of the fit.
- **Masses and S are checked against each other at load.** For two particle pairs of one compound nucleus,
  S_b − S_a must equal (m₁b + m₂b − m₁a − m₂a)·u. Thresholds come from S, kinematics from the masses; a mismatch
  costs about the mass error in the kinematics but exposes typos and mixed conventions. The engine warns above
  1 keV (three listed, then a count), adds an electron-mass hint near −0.511 MeV multiples (a bare proton 1.00728 u
  among atomic masses), and prints one note for projects with mass numbers as masses. γ and β pairs are skipped.
  In an archive of 5,705 projects with two or more particle pairs, 1,138 agreed to 0.1 keV; 4-decimal masses give
  tens of keV. Test: `tests/mass_separation_check`.
- **AME2020 is shipped** (`pyazr/data/mass_1.mas20`): `pyazr.sources.masses.atomic_mass("12C")`,
  `separation_energy("p", "12C")`, `pair_masses("a", "12C")`. `masses.fill_masses(project, pairs=None,
  separation=True)` sets the masses and S of each pair through `Project.set_pair`; particles are identified by Z and
  round(M), so mass numbers and bare protons are recognised; γ pairs get M(compound) and keep S = 0.
- **Filling S moves thresholds.** Never fill an archived fit in place: copy it, fill, compare χ², and leave any refit
  to the evaluator.
- Behaviour: `docs/source/user_guide/particle_pairs.rst`.

### Frames, units and directories

- **Every output and API result is c.m.**; a data file states the frame of each column (`E_lab[MeV]`, `E_cm[keV]`,
  `theta_lab[deg]`, `theta_cm[deg]`); an extrapolation grid states its frame (`MeV lab`, `MeV cm`).
- **Relative paths are relative to the project file** (`output/`, `checks/`, data, covariance and table files);
  pyazr sets `cwd` to the file's directory. `output/` and `checks/` must exist.
- **A lab angle the forward-branch kinematics cannot reach is refused** (`ERROR: ... beyond the kinematic maximum`):
  past the maximum lab angle (an endothermic exit near threshold, elastic scattering with m₁ > m₂: arcsin(m₂/m₁)) or
  below the exit threshold. Give such points in c.m. or drop them. Test: `tests/kinematic_limits`.

### Kinematics: non-relativistic unless `relativistic=1`

- Without the `<config>` key the kinematics is non-relativistic. With `relativistic=1`
  (`Project.set_config(relativistic=True)`) k, η, lab↔c.m. energies, lab angles and cross sections come from the
  invariant mass (the LANL/EDA prescription, ENDF KRL = 1; `include/Kinematics.h`). Level energies stay excitation
  energies; `output/relativistic.out` lists them in EDA's E(s) convention (E_s − E_x = T₀²/(2m₀), 43 keV at T₀ = 20 MeV
  for A = 5).
- It matters above E_lab ≈ 10 MeV and above all between partitions with a large Q: in ⁵He an n+⁴He point at
  E_n = 20 MeV moves 27 keV relative to d+t.
- It is a model choice: a project fitted with it is evaluated with it (the key travels with the file), and switching
  it moves χ². It is refused with a `potential` record.
- The relativistic lab→c.m. map is an exact boost; it agrees with an independent reaction-kinematics calculation to
  1e-7 deg and 4e-9 in the Jacobian up to 200 MeV when the `.azr` carries full-precision masses and S (6-decimal S
  values gave 4e-4 deg near thresholds). The relativistic maximum lab angle decides which angles are reachable; it
  is 0.3–0.6° smaller than the non-relativistic one just above the ⁷Li(p,n) and ¹²C(α,p) thresholds.
- Remaining linearisations in relativistic mode: target-effect kernels, sub-point offsets and energy shifts convert
  lab→c.m. with m₂/(m₁ + m₂) instead of m₂/√s, a relative error of about 4e-4 on offsets (tens of eV over a 50 keV
  target).
- rmfit's convolution pieces use the non-relativistic map: set `conv_grid=False` for a relativistic model.
- Tests: `tests/relativistic`. Behaviour: `docs/source/user_guide/configure_menu.rst`,
  `docs/source/reference/data_formats.rst`.

### Levels and channels

- **A J^π is in the calculation only if an active level carries it, unless `includeAllJPi=1`** in `<config>`
  (references/azr-levels.md §6): without it an elastic curve misses the hard-sphere scattering of every J^π without a
  level (p + ¹²C elastic off by up to 14 %). `audit_project.py` notes it.
- **The `gamma` field is a physical width** (references/azr-levels.md §4).
- **Every external-capture final state is a level, and only one level of its J^π lies within 1 keV of the pair's
  E_x** (references/azr-levels.md §7).
- **An excited light particle** is column 32 (`e1`, MeV, written only when nonzero); the pair's total excitation is
  e1 + e2. Identical particles need equal Z, mass, J, parity, g-factor and excitation (α + α* is not identical;
  e1 = e2 is). A UPOS segment needs the excitation on the heavy partner. Engines that read 31 columns ignore column
  32 silently and use the wrong threshold. `Project.set_pair(key, light_excitation=E)`. Test: `tests/light_excitation`.
- **J is stored as J; channel spin and l are stored doubled** (2S, 2L) (references/azr-levels.md §2).
- **All levels of one J^π share one channel set**; `Project.add_level` clones the group's channels.

### Data sets

- Read `AZURE2 --describe model.azr` for every data set before modelling it: what it computes, rows kept.
- Totals, sums, ratios, natural targets, secondary γ rays: references/data-kinds.md §2–4.
- Parameter files are matched by name (`<data set>.norm`); one written before a conversion is renamed with
  `AZURE2 --convert-params`.

## References

- references/azr-levels.md: the `<levels>` columns, rules the engine applies, channel enumeration, `includeAllJPi`,
  capture final states, verifying a new file.
- references/data-kinds.md: the data set record, which quantity, totals and sums, natural targets, provenance.
- references/nds-exfor.md: `pyazr.sources.nds` search, fetch, conversion, levels, bibliography, service failures.
- references/rmatrix-picture.md: what the parameters mean and which levels to include.
- references/project-from-python.md: a project written from Python (`pyazr.Project`).
- references/key-papers.md: annotated primary literature.

Sphinx pages: `docs/source/user_guide/particle_pairs.rst`, `levels_channels.rst`, `segments.rst`,
`project_from_python.rst`, `docs/source/reference/format3.rst`, `docs/source/examples/natural_neon.rst`. Examples:
`pyazr/examples/exfor_fetch.py`, `edit_scheme.py`, `print_scheme.py`, `project_from_python.py`. Test projects to copy
from: `tests/upgrade/T1` (capture and scattering), `tests/13N`, `tests/natural_neon`.
