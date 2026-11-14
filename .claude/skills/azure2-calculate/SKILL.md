---
name: azure2-calculate
description: Compute AZURE2 observables from a format-3 project and a parameter file without changing the parameters - cross sections, S factors, differential cross sections, A_y, Legendre coefficients, phase shifts and scattering lengths on the data or on any grid (extrapolation records), S(0), reaction rates, level/interference/external-capture decomposition, adaptive tables; the CLI menu and pyazr sessions. Use for "calculate", "extrapolate", "rate", "S factor", "decompose". Fitting: azure2-fit; trust checks: azure2-check; bands: azure2-uncertainty; figures: azure2-plot.
---

# Calculating observables with AZURE2

A calculation takes a project (`.azr`, format 3: `<levels>`, data sets and extrapolations as records) and a parameter
vector (the project's own `<levels>` values, or a fit's `param.sav` matched by name) and returns, per point, the c.m.
observable: σ (b), dσ/dΩ_c.m. (b/sr), S (MeV b), A_y, Legendre coefficients, a phase shift (deg). This skill runs that
calculation on the data points, on grids of your choice, as a rate N_A⟨σv⟩ (cm³ mol⁻¹ s⁻¹), or as a decomposition.

| Task | Skill |
|---|---|
| any forward evaluation, extrapolation, rate, decomposition, output files | this one |
| build a project, masses, levels, data records | azure2-setup |
| move parameters (fits, priors, editing levels) | azure2-fit |
| decide whether a result can be trusted | azure2-check |
| bands, covariance, sampling | azure2-uncertainty |
| figures | azure2-plot |
| target, beam, detector effects | azure2-targets |
| THM data | azure2-thm |

Behaviour: `docs/source/user_guide/calculate.rst`, `docs/source/user_guide/pyazr.rst` ("The objective and the
evaluation", "Frames and units"), `docs/source/reference/command_line.rst`, `docs/source/reference/output_files.rst`.
Worked examples: `pyazr/examples/` (`references/examples.md`).

## Scripts

| Script | Does |
|---|---|
| `scripts/extrapolate.py` | one bare extrapolation record on any grid: E_c.m., E_x, θ, value, S |
| `scripts/decompose.py` | full, level off, level only, bare; interference; the pairwise term of two levels |
| `scripts/selftest.py` | runs both on `tests/upgrade/T1` |

```bash
S=.claude/skills/azure2-calculate/scripts
python3 $S/extrapolate.py model.azr --params output/param.sav --reaction "12C(p,g)13N" \
    --energies "range 0.001 0.6 MeV cm points 200 log" --out sfactor.txt
python3 $S/decompose.py model.azr --params output/param.sav --level "1/2+ 2.365" --with "3/2- 3.502" \
    --dataset cap_g0 --out decomposition.txt
```

## Two ways to run the engine

- **The program** (`AZURE2 --no-gui --no-readline project.azr`, the menu piped): one-shot calculate (1), fit (2),
  calculate the extrapolations (3), MINOS (4), rate (5), MCMC (6). Use it when its output files are the product.
- **pyazr** (`pyazr.Session`): the engine in the interpreter, for everything else (scans, decompositions, samplers,
  figures). Installed by `pip install .` from the repository; the Coulomb table once with
  `pyazr coulomb-table fetch` (or `build`).

## Procedure A: the program

1. **The options are the project's.** A project states its options as `<config>` keys (`use_brune=0`,
   `use_long_wavelength=0`, ...); a flag given wins over the key. `AZURE2 --describe project.azr` prints every option
   with its value and source: read it before trusting a run of an inherited project. Without keys, take the flags
   from the project's job scripts, never from defaults. Bases and their conversion: azure2-exchange.
2. **Pass flags as words.** zsh does not split an unquoted variable (`F="--use-brune --no-readline"; AZURE2 $F` passes
   one unknown argument, which is ignored): use an array or a shell function.
3. **Create `output/` and `checks/`** beside the project (paths in `<config>` are relative to the project file).
4. **Answer the prompts in order** (`command_line.rst`, "Scripted input"): menu → (2, 3, 4) the band question →
   parameter file (blank: the project's `<levels>`) → capture-amplitude file (blank: compute) → `7`.

   ```bash
   az() { AZURE2 --no-gui --no-readline "$@"; }
   printf '1\n\n\n7\n'                     | az project.azr   # calculate at the project's values
   printf '1\noutput/param.sav\n\n7\n'     | az project.azr   # calculate at a fit (names matched)
   printf '3\nn\noutput/param.sav\n\n7\n'  | az project.azr   # calculate the extrapolation records
   printf '5\noutput/param.sav\n1\n-1\nyes\ntemperatures.dat\n' | az project.azr   # rate, pair 1, total capture
   ```

   A y/n prompt that receives anything else stops a piped run with exit 1; the log says
   `Reading User Parameter File...` when a parameter file was applied.
5. **Check the exit status and `output/chiSquared.out`** (its `Total-Chi-Squared` line is the quickest check).

## Procedure B: pyazr

```python
import pyazr
with pyazr.Session("project.azr") as s:          # cwd defaults to the project's directory
    x = s.parameters.read("output/param.sav")    # the fit vector, matched by name (or s.parameters.x0)
    ev = s.evaluate(x)                           # every number of one evaluation
    d = ev.data[0]                               # per data set: model, energies (c.m.), angles, values, errors
    print(ev.chi2, d.name, d.model[:3])
```

- **Options**: `Session(file, use_brune=..., use_long_wavelength=..., ...)`; a keyword wins over the `<config>` key.
  `s.options` holds the values in effect and `s.option_sources` where each came from. Pass the project's own options
  when comparing with an archived CLI run (pyazr's default is `use_long_wavelength=True`).
- **One vector**: the fit vector `x` (reduced-width amplitudes in the session's basis, energies in MeV, norms, shifts,
  ...); `s.parameters` names its entries (`print(s.parameters.free.table())`). A vector of one basis means nothing in
  another (`s.basis`: 0 standard, 1 Brune, 2 Park). Physical values (Γ in eV, ANCs): `s.transform(x)`, one per row
  of `s.parameters`.
- **Data or extrapolations**: a session evaluates either its data sets (default) or its extrapolation records
  (`Session(file, extrapolations=True)`); open one session of each kind rather than switching. In extrapolation mode only
  active records are evaluated and every χ² field is `None`; find a record by `d.name`, not by position.
- **Threads**: every evaluation uses `OMP_NUM_THREADS` threads; two sessions may evaluate on two threads at once; for
  more cores give each worker process its own session (`pyazr.rst`, "Running in parallel").
- **Never point an exploratory session at a live fit's output directory**: a session whose options or points differ
  from the cached ones rewrites `intEC.dat`/`intEC.extrap` there. Give the copy its own directory
  (`Project.set_config(output_dir="out2/")`, and create it).

### Observables

```python
d = s.evaluate(x).data[i]
d.model                                  # σ (b), dσ/dΩ_c.m. (b/sr), A_y, phase (deg): as χ² uses it
d.model * d.extra["sfactor_conversion"]  # S factor (MeV b)
d.energies, d.excitation_energies        # E_c.m. and E_x (MeV)
d.angles                                 # θ_c.m. (deg)
d.extra["legendre"]                      # Legendre coefficients of an angular distribution, per point
d.extra                                  # every per-point number: pyazr.kinds()["data_outputs"]
s.write_outputs(x)                  # AZUREOut_*, chiSquared.out into the output directory
```

E_x is the only axis common to all entrance pairs. A capture data set of `multipole E1` gives the E1 part as `model`.
Special observables (A_y, Legendre data, phase shifts, sums and ratios, total capture, scattering lengths):
`references/observables.md`.

## Extrapolation

An `extrapolation` record (`<extrapolations>`) is a grid with a reaction and a quantity; `scripts/extrapolate.py` adds
one to a copy of the project and evaluates it. By hand:

```python
p = pyazr.Project.open("model.azr")
p.add("extrapolation", "sfac", reaction="12C(p,g)13N", quantity="sigma", unit="b",
      energies="range 0.001 0.6 MeV cm points 200 log")
p.save("model_grid.azr")                  # beside the original: relative paths hold
with pyazr.Session("model_grid.azr", extrapolations=True) as s:
    d = next(d for d in s.evaluate(s.parameters.read("output/param.sav")).data if d.name == "sfac")
```

1. **A grid states its frame** (`MeV lab` or `MeV cm`); returned energies are always c.m. Grid grammar:
   `docs/source/reference/format3.rst` ("extrapolation"); `pyazr.extrapolation_grid(text, masses=...)` shows it.
2. **Bare or folded**: an extrapolation is folded only by the `target`, `beam` or `detector` it names. For a bare σ
   or S, name none (`Project.set(name, target=None)` removes one).
3. **Fitted norms and shifts do not reach extrapolations**: the parameters of data sets are not in a run without
   data, and a parameter file's data-set entries are skipped there.
4. **Resolve narrow resonances**: 0.25–0.5 keV steps across a keV-wide level, or an adaptive table (below). A grid
   too coarse misdraws a peak by up to a factor 10.
5. **S(0)**: a log grid down to about 1 keV; read S at the lowest points and check that a lower start does not move
   it. Its band: azure2-uncertainty.
6. **A secondary γ ray (`secondary_gamma` of a data set) is not an extrapolation of the primary reaction**: to draw
   such a curve, evaluate a copy of the data set on a pseudo-data file of the grid energies.

## Decomposing a cross section

At fixed parameters, zero the free widths of chosen levels in `x` (`P.of_level(label)`, `P.put`) and evaluate
(`scripts/decompose.py` does all four):

```python
P = s.parameters
w = [p for p in P.of_level("1/2+ 2.365").of_kind("width") if not p.fixed]
off = s.evaluate(P.put(x, w, 0.0)).data[i].model        # the model without the level
```

- `full − off` is everything the level does (resonance plus interference); `only − bare` the bare resonance; their
  difference the interference.
- Levels of different J^π do not interfere in an angle-integrated σ: the pairwise term
  `only{T,L} − only{T} − only{L} + bare` is zero there (the self-test checks it). A decomposition that breaks this
  zeroed the wrong parameters.
- **Keep the final states' ANCs**: external capture scales with the ANC of the bound final state, which is a width
  of that level. "bare" with every width zeroed has no direct capture left; `decompose.py` keeps the bound levels
  unless `--zero-bound`.
- Fixed widths are not entries of `x`: zero them in a copy of the project (`Project.set_channel(level, pair, l, s,
  width=0.0)`) and open a new session.

## Reaction rates

- **Menu 5** integrates σ(E_c.m.) on a grid that follows the parameters (48 points per decade from 0.1 eV to 50 MeV,
  plus a cluster around every level), writes `reactionrates.out` (N_A⟨σv⟩ against T9) and `cross_section_grid.dat`
  (plot it when a rate looks odd). Exit pair `-1` is the total radiative capture (every capture pair of the
  entrance's compound nucleus). Narrow resonances need no hand-added terms (a 22 eV level to 1e-4,
  `tests/reaction_rate`).
- **In Python**: `pyazr.rate` runs the engine's rate in-process (`reaction_rate(file, entrance, exit, grid=...)`,
  `RateSession` for many vectors, `rate_band` from a covariance or from draws, `gamow_window`); the program's
  digits (`tests/pyazr/reaction_rate_test.py`); `pyazr/examples/reaction_rate.py` runs it.

## Adaptive cross-section tables

`pyazr.tabulate(azr_file, pairs=[(entrance, exit)], e_min, e_max, rel_tol=5e-3, frame="cm")` refines a grid until
log-linear interpolation reproduces the engine to `rel_tol`; `.save(path)` writes `E_cm E_lab sigma`. Central values
only. Example: `pyazr/examples/tabulate_cross_sections.py`.

## Caches

- `output/intEC.dat` and `intEC.extrap` hold external-capture integrals; each carries a signature of its grid and
  options (`.sig`) and is recomputed at session start when it no longer matches. A running session does not re-check:
  after editing data or grids, open a new session (or `s.recalculate_external_capture()`). Safe to delete.
- The universal Coulomb table (`pyazr coulomb-table fetch`, or `AZURE2 --build-coulomb-table`;
  `$AZURE2_COULOMB_TABLE` names the file) speeds charged-particle runs about 4× without changing printed digits.

## Pitfalls by topic

- Runtime options that change the physics (route, capture operator, hybrid potential, Coulomb functions, tables):
  `references/physics-options.md`.
- Output columns, the norm in them, the χ²-optimal norm: `references/outputs.md`.
- Results made by older builds: azure2-check ("Constants and logs").
- Every figure is checked against the engine's model at the data points (azure2-plot).
