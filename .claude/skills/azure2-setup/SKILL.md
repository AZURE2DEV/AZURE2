---
name: azure2-setup
description: Build or audit an AZURE2 .azr project before any fit - compound nucleus and thresholds, particle and photon pairs with AME masses and separation energies, levels and channels from ENSDF, channel radii, EXFOR data converted into segments (pyazr.nds), what each segment measures (totals, natural targets), frames and units, and the first checks. Use for a new reaction, a new data set, reading or writing level lines or data-set records, or EXFOR/ENSDF lookups. Running a model is azure2-calculate; fitting azure2-fit; target effects azure2-targets; THM projects azure2-thm.
---

# Setting up an AZURE2 project

An AZURE2 evaluation starts from one `.azr` file that states the compound nucleus, its partitions (particle pairs
with masses m₁, m₂ in u and separation energies S in MeV), the levels (J^π, excitation energy E_x in MeV, partial
widths Γ in eV or ANCs in fm^-1/2) and channels (pair, 2S, 2L, channel radius a_c in fm), and the data segments
(lab-frame data files with their observable codes and tails). Every later number depends on these inputs, and most
setup errors load without a message. This skill gives the order in which to build a project, the conventions it
must follow, and the checks that show it is what was intended.

## When to use, and when not

- Use it for a new compound nucleus or reaction, for adding a data set, for fetching data or level schemes from
  EXFOR and ENSDF, for reading or hand-writing `<levels>` lines and `dataset` records, and for auditing an inherited
  project before trusting it (an old `.azr` is refused: convert it with `AZURE2 --convert OLD NEW`, then audit the
  conversion and its report; `docs/source/user_guide/older_projects.rst`).
- Running calculations, extrapolations, rates and decompositions: azure2-calculate. Fitting, editing levels of a
  fitted model, norms and shifts: azure2-fit. Whether a result can be trusted: azure2-check. Covariance files for
  data: azure2-uncertainty. Target thickness, resolution and Doppler: azure2-targets. THM data: azure2-thm.
  Format conversion (GNDS, ENDF, other codes): azure2-exchange.

## Procedure

1. **The compound nucleus and its thresholds.** List every partition open in the energy range of the data, and the
   photon pairs (one per final state reached by capture). Fetch the adopted level scheme and ground-state data with
   `pyazr.nds` (references/nds-exfor.md §5).
2. **Pairs: masses and separation energies from one table.** Enter atomic masses (u, electrons included) from
   AME2020 and S = m₁ + m₂ − M(compound) from the same table: `pyazr.masses.fill_masses(project)`. Then `masses.mass_check(project)` (§ Rules, "Masses").
3. **Channels and radii.** Enumerate the allowed (L, S) for each J^π (references/azr-levels.md §5). Start the
   channel radius at a_c = 1.4 (A₁^⅓ + A₂^⅓) fm. a_c is a model parameter to be scanned (azure2-check), not a
   physical radius.
4. **Levels.** Choose them by references/rmatrix-picture.md §3–4: known resonances in and near the data range,
   sub-threshold states whose tails reach it, background poles above. The `gamma` field carries observed Γ (eV) for
   an open channel and an ANC (fm^-1/2) for a closed one, never a reduced-width amplitude. Add every external-capture
   final state as a level. Add `includeAllJPi=1` to `<config>` for elastic data (references/azr-levels.md §6). A new
   compound nucleus that has no project yet is written line by line (references/azr-levels.md §1–3).
5. **Data.** Fetch and convert with `ExforData.to_azr` (references/nds-exfor.md §4), stating the light/heavy masses
   on every call. Data files are lab frame, forward kinematics (references/azr-segments.md §5). Decide what each
   segment measures: a total over channels, a ratio, a natural target, a secondary γ ray
   (references/data-definition.md, references/azr-segments.md §2). Give every data set a covariance file
   (azure2-uncertainty).
6. **Conventions.** Choose non-relativistic or relativistic kinematics (§ Rules, "Kinematics") once, before the first
   fit. Keep the project in its own directory with `output/` and `checks/` beside it.
7. **First checks** (minutes, no fit):
   - CLI mode 3 with no data, then a mode-1 calculation from the file's own parameters (references/azr-levels.md §8);
   - `pyazr.masses.mass_check(project)`: no mismatch above 1 keV;
   - the segment summary in the run log: every segment filled, every keyword recognised;
   - the tail audit and the shared-file audit (references/azr-segments.md §2, references/data-definition.md §2);
   - each data file plotted against the first calculation in the same frame (azure2-plot).

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
- **AME2020 is shipped.** `pyazr/data/mass_1.mas20`: `pyazr.masses.atomic_mass("12C")`,
  `separation_energy("p", "12C")`, `pair_masses("a", "12C")`. `masses.fill_masses(project, pairs=None, separation=True)`
  sets only the masses and S of each pair through `Project.set_pair` (and the entrance-S copy); particles are identified by Z and round(M), so
  mass numbers and bare protons are recognised; γ pairs get M(compound) and keep S = 0.
- **Filling S moves thresholds.** Never fill an archived fit in place: copy it, fill, compare χ², and leave any refit
  to the evaluator.
- Behaviour: `docs/source/user_guide/particle_pairs.rst`.

### Frames, units and directories

- **Input is lab frame, forward kinematics (light particle = projectile). Every output file and API result is c.m.**
  Never compare across the two without converting.
- **`add_extrapolation(e_min, e_max, e_step)` takes lab energies by default** (`frame="lab"`; `frame="cm"`
  multiplies by (m_beam + m_target)/m_target), while `evaluate(x).data[i].energies` are c.m. energies. Comparing a returned
  c.m. column with the lab e_max looks like an extrapolation that stops early at e_max·m_target/(m_beam + m_target);
  there is no ceiling, only a frame mismatch.
- **Run from the directory that holds the `.azr`.** The file stores `output/`, `checks/` and data paths relative to
  itself; pyazr sets `cwd=` to the file's directory by default. `output/` must exist.
- **A lab angle the forward-branch kinematics cannot reach is refused.** The lab→c.m. map of a lab-differential
  segment (code 1) covers the forward branch only. Past the maximum lab angle (an endothermic exit near threshold:
  ⁷Li(p,n) at E_p = 1.89–1.90 MeV has θ_max = 22–39°; elastic scattering with m₁ > m₂: arcsin(m₂/m₁)) or below the
  exit threshold, the run stops with `ERROR: segment #k: lab angle X deg at E_lab = Y MeV: beyond the kinematic
  maximum ...`. Give such points in c.m. (codes 4 or 5) or drop them. An energy shift that moves a point past the
  limit during a fit warns and keeps the last valid c.m. angle. Recoil detection is not an option (convert by hand);
  A_y angles (code 7) are c.m.; a UPOS lab-differential segment is converted with particle b's kinematics (warned).
  Test: `tests/kinematic_limits`. builds before 2026-10-27 (expires: 2027-02-01): such points became `nan` in
  `chiSquared.out` with exit status 0.

### Kinematics: non-relativistic unless `relativistic=1`

- Without the `<config>` key the kinematics is non-relativistic. With `relativistic=1` (pyazr `Project.set_config(relativistic=True)`) k, η, lab↔c.m. energies, lab angles and cross sections come
  from the invariant mass (the LANL/EDA prescription, ENDF KRL = 1; `include/Kinematics.h`). Level energies stay
  excitation energies; `output/relativistic.out` lists them in EDA's E(s) convention (E_s − E_x = T₀²/(2m₀), 43 keV
  at T₀ = 20 MeV for A = 5).
- It matters above E_lab ≈ 10 MeV and above all between partitions with a large Q: in ⁵He an n+⁴He point at
  E_n = 20 MeV moves 27 keV relative to d+t.
- It is a model choice: a project fitted with it is evaluated with it (the key travels with the file), and switching
  it moves χ². It is refused with a `<potential>` block.
- The relativistic lab→c.m. map is an exact boost; it agrees with an independent reaction-kinematics calculation to
  1e-7 deg and 4e-9 in the Jacobian up to 200 MeV when the `.azr` carries full-precision masses and S (6-decimal S
  values gave 4e-4 deg near thresholds). The relativistic maximum lab angle decides which angles are reachable; it
  is 0.3–0.6° smaller than the non-relativistic one just above the ⁷Li(p,n) and ¹²C(α,p) thresholds.
  builds before 2026-10-29 (expires: 2027-02-01): the non-relativistic limit decided, and a point between the two
  limits got a non-relativistic angle and Jacobian (9.5 % off in the Jacobian for ⁷Li(p,n) at 1.90 MeV, 44.9°).
- Remaining linearisations in relativistic mode: target-effect kernels, sub-point offsets and energy shifts convert
  lab→c.m. with m₂/(m₁ + m₂) instead of m₂/√s, a relative error of about 4e-4 on offsets (tens of eV over a 50 keV
  target).
- rmfit's convolution pieces use the non-relativistic map: set `conv_grid=False` for a relativistic model.
- Tests: `tests/relativistic`. Behaviour: `docs/source/user_guide/configure_menu.rst`,
  `docs/source/reference/data_formats.rst`.

### Levels and channels

- **A J^π is in the calculation only if an active level carries it, unless `includeAllJPi=1`**
  (references/azr-levels.md §6). Check the Levels tab status line or `grep includeAllJPi` before trusting an elastic
  curve from an older project.
- **The `gamma` field is a physical width** (references/azr-levels.md §4).
- **Every external-capture final state is a level, and only one level of its J^π lies within 1 keV of the pair's
  E_x** (references/azr-levels.md §7).
- **An excited light particle** is column 32 (`e1`, MeV, written only when nonzero); the pair's total excitation is
  e1 + e2. Identical particles need equal Z, mass, J, parity, g-factor and excitation (α + α* is not identical;
  e1 = e2 is). A UPOS segment needs the excitation on the heavy partner. Engines that read 31 columns ignore column
  32 silently and use the wrong threshold. pyazr: `Pair.light_excitation`. Test: `tests/light_excitation`.
- **J is stored as J; channel spin and l are stored doubled** (2S, 2L) (references/azr-levels.md §2).
- **All levels of one J^π share one channel set**; the line format otherwise refuses the file.

### Segments and data

- Read every segment's tail before modelling it (references/azr-segments.md §2).
- A data set names its target, beam and detector records; parameter-file entries are matched by name
  (`<data set>.norm`, `<target>.areal_density`); a parameter file of the old project is renamed with
  `AZURE2 --convert-params` (references/azr-segments.md §4 for old files).
- A total is summed over every channel it measures; same-file lines are not summed; prefer `total 1` for neutron
  totals (references/data-definition.md §1–2).
- Natural targets are one SUM segment with `sample`, built with `merge_projects` (references/data-definition.md §3).
- Data converted from EXFOR need the frame of their angle column and an explicit light/heavy assignment
  (references/nds-exfor.md §4).

## References

- references/azr-levels.md: the 31 (+2) level-line fields, rules the engine applies, channel enumeration,
  `includeAllJPi`, capture final states, verifying a new file.
- references/azr-segments.md: blocks of the `.azr`, segment line and tail keywords, test segments, key matching,
  the data file.
- references/data-definition.md: totals, shared data files, natural targets and samples.
- references/nds-exfor.md: `pyazr.nds` search, fetch, conversion, levels, bibliography, service failures.
- references/rmatrix-picture.md: what the parameters mean and which levels to include.
- references/key-papers.md: annotated primary literature.
- references/project-from-python.md: a format-3 project written from Python (`pyazr.Project`): schema, records, the
  engine's diagnostics, describe, calculate.

Sphinx pages: `docs/source/user_guide/particle_pairs.rst`, `docs/source/user_guide/levels_channels.rst`,
`docs/source/user_guide/segments.rst`, `docs/source/reference/data_formats.rst`, `docs/source/user_guide/pyazr.rst`,
`docs/source/examples/natural_neon.rst`. Examples: `pyazr/examples/exfor_fetch.py`, `pyazr/examples/edit_scheme.py`,
`pyazr/examples/print_scheme.py`, `pyazr/examples/project_from_python.py`. Test projects to copy from: `tests/13N`, `tests/natural_neon`.
