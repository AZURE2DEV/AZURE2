---
name: azure2-calculate
description: Compute AZURE2 observables from an existing .azr and parameter file -- cross sections, S factors, differential cross sections, analyzing powers, Legendre coefficients, phase shifts and scattering lengths at any E and angle; extrapolation (mode 3, extrap_mode), S(0), reaction rates (mode 5), level/interference/external-capture decomposition, adaptive cross-section tables, the CLI menu and pyazr sessions. Use for "calculate", "extrapolate", "rate", "S factor", "decompose". For fitting use azure2-fit; trust checks azure2-check; bands azure2-uncertainty; figures azure2-plot.
---

# Calculating observables with AZURE2

An AZURE2 calculation takes a project (`.azr`: particle pairs, levels and channels, data and test segments) and a
parameter vector (the `<levels>` block, or a fit's `param.sav`) and returns, per segment and point, the c.m.
observable: σ (b), dσ/dΩ_c.m. (b/sr), the S factor (MeV b), A_y, Legendre coefficients a_k, or a phase shift
(deg). This skill says how to run that calculation without changing the parameters: on the data points (data mode),
on user grids (extrapolation), as a thermonuclear rate N_A⟨σv⟩ (cm³ mol⁻¹ s⁻¹), or as a decomposition into level,
interference and external-capture parts.

## When to use, and when a sibling instead

| Task | Skill |
|---|---|
| any forward evaluation, extrapolation, rate, decomposition, output files | this one |
| build a new project, masses, levels, EXFOR data, segment syntax | azure2-setup |
| move parameters (fits, free/fixed, priors, editing levels, `.sav` names) | azure2-fit |
| decide whether a result can be trusted (pulls, χ²/ν, Wigner limits, radius) | azure2-check |
| uncertainty bands, covariance, sampling | azure2-uncertainty |
| figures, and showing a result to a person | azure2-plot |
| target integration, resolution, Doppler, transmission | azure2-targets |
| THM (HOES) segments | azure2-thm |

## Two ways to run the engine

- **CLI** (`$AZURE2_BINARY --no-gui --no-readline project.azr`): one-shot calculate, fit, extrapolate, MINOS, rate or
  MCMC run that writes the standard output files. Use it when the files are the product.
- **pyazr** (`from pyazr import azure2, Project`; format-3 projects only): the engine compiled into the extension module `pyazr/_azure2`,
  running in the interpreter. Use it for all evaluation work: χ², scans, decompositions, custom fitters and
  samplers. pyazr does not call the CLI binary; rebuilding `_azure2` is the update. Its version is
  `pyazr.__version__`.

Install pyazr from the repository root: `pip install -e .` (core: numpy, mpmath, scipy), `".[examples]"` adds
matplotlib, emcee, multiprocess; `".[all]"` adds zeus-mcmc. The C++ module is not built by pip: build `_azure2` with
CMake (`USE_API=ON`, the default), which places it in `pyazr/`.

Behaviour reference: `docs/source/reference/command_line.rst` (flags, prompts, exit status),
`docs/source/reference/output_files.rst`, `docs/source/user_guide/calculate.rst`, `docs/source/user_guide/pyazr.rst`,
`docs/source/user_guide/segments.rst`. Worked examples: `pyazr/examples/` (catalogue in `references/examples.md`).
Runnable projects: `tests/13N`, `tests/13N_capture_ay`, `tests/hybrid_potential`.

## Procedure A: the CLI

1. **Pass the runtime options as flags, every time.** The CLI does not read the Runtime Options saved in the `.azr`.
   Take the flags from the project's own job scripts, never from memory or defaults. The parametrization is a flag:
   Brune is the default; `--use-park` selects Park's observed amplitudes; `--no-brune` the Lane–Thomas basis with
   constant boundary conditions; `--use-rmc` the Reich–Moore formalism (switches Brune and Park off; for (n,γ) only).
   Conventions and conversions between bases: azure2-exchange.
2. **Write the flags out.** zsh does not word-split an unquoted variable: `F="--no-readline --use-brune"; AZURE2 $F x.azr`
   passes one argument, which AZURE2 reports as an unknown option and ignores, so the run uses other physics. Use a
   zsh array (`F=(--no-readline --use-brune)`) or `${=F}`; csh job scripts with literal flags are fine.
3. **Create `checks/` and `output/`** in a new run directory. A missing `checks/` stops the run before it starts
   ("Could not find checks directory"), and the job still ends normally with no output.
4. **Answer the prompts in order** (`docs/source/reference/command_line.rst`, "Scripted input"):
   menu choice → (modes 2, 3, 4) "Calculate cross-section uncertainty band? (y/n)", and after a yes in a fit mode the
   reduced-χ² scaling question → external parameter file (blank: `<levels>`) → external-capture amplitude file
   (blank: compute) → `7` to exit. Both file prompts appear whether or not the model has capture.

   ```bash
   az() { "$AZURE2_BINARY" <flags> --no-gui --no-readline "$@"; }   # a function, not a string variable (step 2)
   printf '1\n\n\n7\n'                  | az project.azr  # calculate from the .azr <levels>, both files blank
   printf '1\noutput/param.sav\n\n7\n'    | az project.azr  # calculate at a fit's parameters
   printf '2\nn\n\n\n7\n'               | az project.azr  # fit from <levels> (band: no)
   printf '2\nn\noutput/param.sav\n\n7\n' | az project.azr  # fit from saved parameters
   printf '3\nn\noutput/param.sav\n\n7\n' | az extrap.azr   # extrapolate: a copy with <targetInt> inactive
   ```

   Give a saved `output/param.sav` at the parameter-file prompt to start from the best-fit formal parameters, or a
   blank line to start from the `.azr` levels.

   `--covariance-band` (with `--scale-covariance`) answers the band question yes without a prompt.
5. **Confirm the parameter file was read.** The log prints `Reading User Parameter File...` when it was, and
   `Creating New param.par File...` when the run started from `<levels>`. Current builds stop a piped run with exit
   status 1 when a y/n prompt receives a file name; builds before 2026-11-01 (expires: 2027-02-01) took any other
   answer as "no" and silently ran from `<levels>`.
6. **Check the exit status and `output/chiSquared.out`.** Exit status 1 follows "Calculation was aborted." or a
   refused project. The last line of `chiSquared.out` (total χ²) is the quickest scalar check.

| Mode | What it does | Main output |
|---|---|---|
| 1 | calculate the data segments | `AZUREOut_*.out`, `chiSquared.out` |
| 2 | fit (Minuit2), see azure2-fit | `param.sav`, `parameters.out` |
| 3 | calculate the test segments (extrapolate) | `AZUREOut_*.extrap` |
| 4 | MINOS errors | `param.errors`, `covariance_matrix.out` |
| 5 | reaction rate | `reactionrates.out`, `cross_section_grid.dat` |
| 6 | MCMC sampling, see azure2-uncertainty | `samples.mcmc` |
| 7 | exit | |

## Procedure B: pyazr

```python
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")      # before importing numpy
import numpy as np
from pyazr import azure2

with azure2("project.azr", cwd=HERE) as m:          # a session is a context manager: close it
    x = m.load_sav("output/param.sav")              # free rwa vector of a fit (load by NAME)
    chi2 = np.sum(m.calculate_chi2_rwa(x))
```

- **Pass the same physics options as the project's CLI fits.** pyazr's defaults can differ from an archive's job
  scripts: `azure2()` defaults to `use_long_wavelength=True`, while many archives run `--no-long-wavelength`. Pass
  `use_long_wavelength=False`, `use_brune`/`use_park`, `ignore_externals`, … explicitly to match.
- **One session per process for parallel work.** Several sessions may be open in one process (each is an independent
  engine that enters its own directory per call), but the engine is not thread-safe. Create the session at module
  level of the worker module; every pool worker then owns one, under `spawn` or `fork`. Builds before 2026-10-07 (expires: 2027-02-01)
   shared the external-capture amplitude cache between sessions of one process: the first
  session computed wrong capture cross sections (up to 28 %) once a second session on the same model existed;
  multi-session capture results from such a build are suspect (test `tests/pyazr/shared_state_test.py`).
- **A session parses its `.azr` once** and re-reads it when the file changes on disk or `configure()` runs. In builds
  before 2026-11-03 (expires: 2027-02-01) every `calculate_rwa` re-parsed the file per segment (3–75 ms per call on
  NFS); look for this before blaming the worker pool.
- `m.datasets`, `m.active_datasets` and `m.extrapolations` return copies; editing them changes nothing in the
  session.
- **Never point an exploratory session at a live fit's `output/`.** Opening a session whose runtime options or
  segment count differ from the cached ones recomputes and overwrites `output/intEC.dat`/`intEC.extrap` there, even
  for a calculation. Work on a copy with its own output directory (`Project.set_config(output_dir="out2/")`). To repair, rerun a
  CLI mode-1 calculate under the project's flags and confirm `chiSquared.out` reads as before.

### Two parameter vectors

- `*_rwa` methods (`calculate_rwa`, `calculate_chi2_rwa`, `calculate_sfactor_rwa`, `chi2_and_grad`,
  `residual_jacobian`) take the **reduced-width-amplitude** free vector `m.params_rwa`: the fit space and the only one
  with analytic derivatives. Default to it. Its content depends on the basis (`m.basis`: 1 Brune, 2 Park,
  0 standard); a vector from one basis is meaningless in another.
- Plain methods (`calculate`, `calculate_chi2`) take the **physical** free vector `m.params` (E_λ in MeV, partial
  widths Γ in eV).
- `m.transform_rwa(x)` maps rwa → physical for the leading R-matrix block only (energies and widths) and returns that
  block, indexed by `free_index`, whatever the input length. Never feed its result to `calculate`,
  `calculate_sfactor` or `calculate_chi2`: they take the full physical free vector. Current builds raise
  `AZURE2Error: ... expected the free parameter vector (N values ...)`; builds before 2026-10-08 (expires: 2027-02-01)
   accepted a short vector silently, leaving norms and shifts at whatever the engine last held (σ and S
  right, χ² meaningless, non-monotonic in the scanned parameter).
- Both vectors hold only free parameters, in `.azr` order: `p.free_index` is the position in the vector, `p.index`
  the position among all parameters (`m.fixed_params`, `param.sav` lines). Parameter names and `.sav` handling:
  azure2-fit.
- `calculate_rwa` returns, per segment, the quantity χ² is built from: the E1 or E2 component for an
  `angle-integrated-E1`/`-E2` segment, the total otherwise.

### Data mode and extrapolation mode

`m.data_mode()` (default) evaluates the `<segmentsData>` segments; `m.extrap_mode()` the `<segmentsTest>` grids. Both
re-initialise the session; switch rarely and batch the work.

- Data mode: segment `i` ↔ `m.datasets[i]` (the i-th `<segmentsData>` line); `len(m.energies[i])` is its point count.
- Extrapolation mode: only active test segments are returned, in file order: segment `i` ↔
  `m.extrapolations.active[i]`, not `m.extrapolations[i]`.

### Observables

```python
m.calculate_rwa(x)                  # σ (b) or dσ/dΩ_c.m. (b/sr) per segment
m.calculate_sfactor_rwa(x)          # S factor (MeV b)
m.calculate_energies(x)             # E_c.m. (MeV) of the calculated points
m.calculate_excitation_energy(x)    # E_x of the compound nucleus (MeV)
m.calculate_angles(x)               # θ_c.m. (deg); differ from the lab angles in the .azr
m.calculate_analyzing_power_rwa(x)  # A_y, analyzing-power segments
m.calculate_angular_dists_rwa(x)    # Legendre coefficients, angular-distribution segments
m.write_output_files(x)             # AZUREOut_*, chiSquared.out into the existing output directory
```

E_x is the only axis common to all entrance pairs; E_c.m. = E_x − (threshold of the pair). `write_output_files` runs
the χ² evaluation first in data mode, so `chiSquared.out` is filled; `param.par` and `parameters.out` are not written
by it (use `m.transform_rwa(x)` for the physical values `parameters.out` would hold). Special observables (A_y,
Legendre coefficients as data, phase shifts, composite SUM/RATIO segments, total capture, scattering lengths):
`references/observables.md`.

## Extrapolation (mode 3, `extrap_mode()`)

1. **Deactivate `<targetInt>` for bare cross sections.** `<targetInt>` entries bind by segment number and apply to
   `<segmentsTest>` lines too: test segment k inherits the convolution of any active entry that lists key k. The
   folded curve looks plausible. Run mode 3 on a copy of the `.azr` with every `<targetInt>` line inactive (first
   field `0`); keep mode-1 χ² checks on the original. A mode-3 run that is slow or memory-hungry is convolving (one
   model: 13k points in 224 s and 65k points at 30 GB with `<targetInt>` active, 1 s with it off). To get the folded
   shape on a fine grid instead, give a data segment a dummy fine-grid data file and run mode 1.
2. **Delete `output/intEC.extrap` (and `intEC.dat`) before every extrapolation run**, not only when the grid changes.
   A session that initialises in data mode and then calls `extrap_mode()` reuses the file on disk, whose amplitudes
   do not match the new data+test integration set. Only capture reads these integrals, so particle exits stay right
   while capture S moves (one model: S(p,γ)(10 keV) 0.0912 keV b with fresh caches, 2.4 × 10⁸ keV b reused). Put
   the deletion in the script:

   ```python
   for f in ("intEC.dat", "intEC.extrap"):
       p = os.path.join(eval_dir, "output", f)
       if os.path.exists(p):
           os.remove(p)
   ```

3. **Fitted shifts and norms do not reach test grids.** A `param.sav` carries `segment_<k>_energy_shift`/`_norm` of
   data segment k. Current CLI and `load_sav` skip `segment_*` entries in a run without data.
   Builds before 2026-10-26 (expires: 2027-02-01) applied data segment k's shift to test segment k, so the printed energies were
   E + shift; strip the `segment_*` lines when using such a build.
4. **State the frame of every grid.** An `extrapolation` record's `energies` and `angles` name their frame (`lab` or
   `cm`); the engine converts a `cm` grid with the masses of `<levels>`:

   ```python
   p = Project.open("3He_a.azr")
   p.add_extrapolation("capture", reaction=p.reaction(1, -1), quantity="sigma", unit="b",
                       energies="range 0.10 8.5 MeV cm step 0.05")
   p.add_extrapolation("el90", reaction=p.reaction(1, 1), quantity="dsigma_domega", unit="b/sr",
                       energies="range 0.05 6.0 MeV lab step 0.05", angles="list 90 deg cm")
   p.save("3He_a_grids.azr")
   ```

   Quantities: `sigma`, `dsigma_domega`, `analyzing_power`, `legendre` (with `expansion`, `orders`), `phase_shift`,
   ... (`pyazr.schema()` lists them with their fields); `p.reaction(entrance, -1)` is every capture pair summed.
5. **Resolve narrow resonances.** A 5 keV grid misdraws keV-wide peaks by up to a factor 11; use 0.25–0.5 keV steps,
   or an adaptive table (below).

## S factor

AZURE2 outputs S directly; never re-derive the Sommerfeld conversion by hand (it is the usual way to evaluate the
factor at the wrong energy).

| Want | Use |
|---|---|
| model S per segment | `m.calculate_sfactor_rwa(x)` (`calculate_sfactor` for a physical vector) |
| data S per point | `m.sfactor[i]`, `m.sfactor_err[i]` |
| from a finished run | S columns of `AZUREOut_*.out` (fit S, data S, data S error) and of `.extrap` |

`calculate_sfactor_rwa` is σ times a factor that depends only on E, so linear combinations of cross sections you
built yourself (a decomposition, an E1/E2 component) convert with `conv = S_full / σ_full` (guard zeros). Units:
MeV b; ×10³ for keV b. S(0) is read from an extrapolation grid that starts near E_c.m. = 0, with its band from
azure2-uncertainty (`extrapolation_bands(..., quantity="sfactor")`).

## Decomposing a cross section

At the fitted parameters, without refitting, so every term is that component's contribution:

```python
m.extrap_mode()
ex   = m.calculate_excitation_energy(best)
full = m.calculate_rwa(best)
off  = m.calculate_rwa(m.without_level(best, jpi="5/2-", energy=6.588))
only = m.calculate_rwa(m.only_level(best,    jpi="5/2-", energy=6.588))
bare_x = best.copy()                                  # every free width zeroed: hard sphere + external capture
for p in m.parameters.widths:
    if p.free_index is not None:
        bare_x[p.free_index] = 0.0
bare = m.calculate_rwa(bare_x)
interference = (full - off) - (only - bare)
```

- `full − off` is everything the level does (resonance plus its interference); `only − bare` is the bare resonance on
  the non-resonant background; their difference is the interference.
- Interference is block-diagonal in J^π: only levels of the same J^π interfere. A decomposition that violates this
  is wrong.
- Pairwise, levels T and L: `only{T,L} − only{T} − only{L} + bare = 2 Re(A_T A_L*)`. Integrate its modulus over E_x
  to rank a level's interference partners.
- Capture: zeroing only the γ-channel widths (`p.radiation_type in ("E", "M")`) leaves every particle channel, and
  thus the scattering solution, unchanged; zeroing all γ widths leaves the external (direct) capture.
- Script: `pyazr/examples/decompose_cross_section.py`.

## Reaction rates

- **Mode 5** integrates σ(E_c.m.) on a grid that follows the evaluated parameters (resonance clusters at the levels
  of the parameter file) from 0.1 eV, in ln E, and writes `reactionrates.out` (N_A⟨σv⟩ in cm³ mol⁻¹ s⁻¹ versus T9)
  and `cross_section_grid.dat` (E_c.m., σ, S: plot it when a rate looks odd). The grid has 48 points per decade from
  0.1 eV to 50 MeV plus a cluster from 1/32 to 64 widths around every level; ln σ is splined against ln E. A 22 eV
  level is integrated to 10⁻⁴ (`tests/reaction_rate`), so narrow resonances need no hand-added terms.
  There is no other path: `useAdaptiveGrid` is not a key of format 3 (a file with it is refused, `E_CONFIG_KEY_REMOVED`).
  Behaviour: `docs/source/user_guide/calculate.rst`, troubleshooting in
  `docs/source/reference/troubleshooting.rst`. Apply the extrapolation rules above (no `<targetInt>`, fresh caches).
- **Total radiative capture rate**: give the exit pair `-1` at the "Reaction Rate Exit Pair" prompt. The rate integrates the sum over every capture pair (pType 10) of the entrance pair's
  compound nucleus; the log lists the pairs. Builds before 2026-11-04 (expires: 2027-02-01) refuse `-1`: run each
  capture pair and add the rates (the rate is linear in σ).
- **In Python**, for a rate per posterior sample or per level: `pyazr/examples/reaction_rate.py`. Its grid must
  resolve the Gamow window at the lowest temperature and every resonance narrower than the spacing.
- Rates from builds before 2026-10-26 (expires: 2027-02-01) are low by 0.8–1.2 % for non-resonant charged-particle
  rates at T9 ≤ 0.1, up to 49 % for neutron (1/v) entrances at T9 = 0.01, and can miss a narrow level moved by the
  parameter file: recompute them.

## Adaptive cross-section tables

`pyazr.tabulate(azr_file, pairs=[(entrance, exit), ...], e_min, e_max, rel_tol=5e-3, frame="cm")` builds σ(E) on a
non-uniform grid by bisection until log-linear interpolation reproduces the engine to `rel_tol`: knots gather at
resonances (a few hundred to a few thousand points where a uniform grid of equal fidelity needs ~10⁵). Each returned
`TabulatedCrossSection` interpolates when called; `.save(path)` writes `E_cm E_lab sigma`. Central values only; bands
come from azure2-uncertainty. Example: `pyazr/examples/tabulate_cross_sections.py`.

## Caches

- `output/intEC.dat` (data segments) and `output/intEC.extrap` (test segments) hold external-capture integrals. Both
  are safe to delete (they cost time to rebuild); delete both after any change of segments, grids or of a flag that
  enters the capture operator. Details and signatures: `references/physics-options.md`.
- The universal Coulomb table (`~/.azure2/coulomb_table.bin`, built once with `AZURE2 --build-coulomb-table`, 88 MB,
  ~5 min on 4 threads) speeds charged-particle runs about 4× and does not change results at the printed digits.

## Pitfalls by topic

- **Runtime options**: R-matrix route, ignore-externals, long-wavelength and convection-current capture operator,
  hybrid potential, Coulomb and EC tables, `--gsl-coul` (do not use): `references/physics-options.md`.
- **Output columns and normalization**: the norm multiplies the data and error columns, not the fit column; the
  χ²-optimal norm n* from one output file; `param.fit` semantics: `references/outputs.md`.
- **Results from older builds**: what each build range computed differently, and what to redo:
  `references/older-builds.md`.
- **Figures**: compare every plotted curve with the fit column of `AZUREOut_*.out` (azure2-plot).
