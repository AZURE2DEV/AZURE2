---
name: azure2-eval
description: Run AZURE2 R-matrix evaluations of nuclear reaction/scattering data — calculate cross sections & chi-squared, fit levels/widths, add/remove levels, decompose a cross section into level/interference/external contributions, extrapolate, compute reaction rates, or drive AZURE2 from Python via pyazr. Use whenever the task involves an .azr project file, AZURE2 levels/channels/segments, R-matrix fitting, S-factors, resonance significance tests, or the worked examples in pyazr/examples/.
---

# Running AZURE2 evaluations

AZURE2 is a multi-channel, multi-level R-matrix code. An evaluation is driven by
a single `.azr` project file that describes the compound nucleus, particle
pairs, levels/channels, and data segments. Two headless ways to run it: the
**interactive CLI** (`--no-gui`) and the **pyazr** Python API. Use the CLI for a
one-shot calculate/fit/extrapolate/rate that writes the standard output files;
use **pyazr for all evaluation work** — χ², parameter scans, level add/remove
tests, cross-section decomposition, custom fitters and samplers.

**pyazr does not run the binary.** The engine is compiled into the extension
module `pyazr/_azure2` and runs in your interpreter, so there is no path to
resolve and nothing to keep in step — a rebuild of `_azure2` *is* the update.
The CLI binary (`build/src/AZURE2`) is a separate artefact, used only for
Workflow A.

`pyazr/` lives at the repo root and is mirrored into evaluation directories, so
`from pyazr import ...` works when cwd is either. Its version is
`pyazr.__version__`; do not assume one.

Reference docs: `docs/source/` — `reference/` covers command_line,
data_formats, output_files; `user_guide/` covers levels_channels, segments,
particle_pairs, fitting, mcmc, pyazr. Worked examples: `pyazr/examples/`;
runnable projects: `tests/13N`, `tests/13N_capture_ay`, `tests/hybrid_potential`.

## Golden rules

- **Input is LAB frame, forward kinematics** (light particle = projectile).
  **All output files and API results are CENTER-OF-MASS.** Never mix them.
  This includes `add_extrapolation(e_min, e_max, e_step)` — those are **lab**
  energies (multiply c.m. by `(m_beam + m_target)/m_target`), while the energies
  that come *back* from `calculate_energies` are c.m.
- Always run **from the directory that contains the `.azr` file** (pyazr does
  this for you via `cwd=`, defaulting to the `.azr`'s directory), because the
  `.azr` stores `output/`, `checks/` and data paths *relative to itself*.
- The `output/` directory must exist; AZURE2 writes results there.
- **Delete `output/intEC.extrap` whenever the `<segmentsTest>` grid changes.**
  AZURE2 caches external-capture integrals there and silently reuses them on a
  different grid — this corrupts capture cross sections *and* data-mode χ²
  (integrals for data and test segments are built together at INITIALIZE).
  `intEC.dat` is the data-segment equivalent and is safe while `<segmentsData>`
  is untouched. Both are safe to delete; they just cost time to rebuild.
- **Several sessions can be open at once**, each an independent engine that
  enters its own directory per call — that is how `save_fit` verifies what it
  wrote. The engine is not *thread*-safe, so parallelism is one session per
  process, not per thread.
- CLI mode does **not** read Runtime Options from the `.azr` — pass them as
  flags every time (`--gsl-coul`, `--ignore-externals`, …). Note the Brune
  parameterization is **on by default** and there is no flag to turn it off;
  `--use-rmc` selects the mutually exclusive RMC formalism, and pyazr takes
  `use_brune=False` directly.

## Workflow A — interactive CLI (one-shot runs)

Prompt order: **(1) menu choice → (2) external parameter file → (3) external
capture amplitude file → (4) mode-specific prompts**, then `7` to exit. Both
file prompts appear whether or not the model has capture, and a recipe that
feeds only one leaves the rest of the answers off by a line. Pipe stdin to run
non-interactively, and pass `--no-readline` so readline does not fight the
pipe. `tests/run_tests.sh` is the working reference.

| # | Mode |
|---|------|
| 1 | Calculate Segments From Data (`AZUREOut_*.out`, `chiSquared.out`) |
| 2 | Fit Segments From Data (Minuit2 → `param.sav`, `parameters.out`) |
| 3 | Calculate Segments Without Data / extrapolate (`AZUREOut_*.extrap`) |
| 4 | MINOS Error Analysis (`param.errors`, `covariance_matrix.out`) |
| 5 | Calculate Reaction Rate (needs temps → `reactionrates.dat`) |
| 6 | MCMC Bayesian Sampling (`samples.mcmc`) |
| 7 | Exit |

The **external parameter file** prompt: give a saved `output/param.sav` to start
from those best-fit formal parameters, or leave **blank** to build fresh from
the `.azr` levels.

```bash
# Calculate from the .azr's own parameters: menu, both file prompts blank, exit.
printf '1\n\n\n7\n' | AZURE2 --no-gui --no-readline 7Be.azr

# Calculate with saved best-fit params (parameter file given, EC file blank):
printf '1\noutput/param.sav\n\n7\n' | AZURE2 --no-gui --no-readline 7Be.azr

# Fit fresh from the .azr levels. Mode 2 then asks about the cross-section
# uncertainty band (y/n) and, if yes, reduced-chi2 scaling (y/n):
printf '2\n\n\nn\n7\n' | AZURE2 --no-gui --no-readline 7Be.azr

# Extrapolate (no data) using saved params:
printf '3\noutput/param.sav\n\nn\n7\n' | AZURE2 --no-gui --no-readline 7Be.azr
```

Non-interactive band control: `--covariance-band` (+ `--scale-covariance`) skips
the y/n prompts. Other flags: `--use-brune`, `--gsl-coul`, `--ignore-externals`,
`--use-rmc`, `--no-transform`; `--help` lists all.

## Workflow B — pyazr

pyazr embeds the AZURE2 engine in-process: the R-matrix code is compiled into
the pybind11 extension module `pyazr/_azure2`, and an `azure2()` object is a
real `AZUREAPI` session living in the interpreter. No subprocesses, no sockets,
no instance pool. A session is a context manager; always close it — that frees
the compound nucleus and data, which are several MB per model.

Install it from the repository root (it is not on PyPI):

```bash
pip install -e .                # core: numpy, mpmath, scipy
pip install -e ".[examples]"    # + matplotlib, emcee, multiprocess
pip install -e ".[all]"         # + zeus-mcmc
```

`-e` is the useful form, since the package lives in the checkout. The engine is
C++ and is **not** installed by pip — build the `_azure2` module with CMake
(`USE_API=ON`, the default), which lands it in `pyazr/`; a `pip install` ships
it as package data. `import pyazr` from the repository root picks up the built
module.

```python
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")   # set BEFORE importing numpy
import numpy as np
from pyazr import azure2, AzrModel

with azure2("7Be.azr", cwd=HERE) as m:
    best = np.asarray(m.params_rwa, float)      # free parameter vector
    chi2 = np.sum(m.calculate_chi2_rwa(best))
```

Every `azure2()` object is an independent engine, and several can be open at
once — each enters its own directory per call, so they never disturb each other
or your cwd. The engine is *not* thread-safe, so parallelism is per process:
construct the session at module level of the worker module and every pool
worker gets its own, under either `spawn` or `fork`.

### Two parameter conventions

- **`*_rwa` methods** (`calculate_rwa`, `calculate_chi2_rwa`,
  `calculate_sfactor_rwa`, `chi2_and_grad`, `residual_jacobian`) take the
  **reduced-width-amplitude** vector `m.params_rwa`. This is the natural fit
  space and the only one with analytic derivatives. **Default to it.**
- **plain methods** (`calculate`, `calculate_chi2`) take the transformed
  **physical** vector `m.params` (level energies in MeV, partial widths in eV).
- `m.transform_rwa(x)` maps rwa → physical; it takes either the full free
  vector or just the leading R-matrix block (energies + widths), so the
  norm/shift tail can be sliced off:
  `nR = 1 + max(p.free_index for p in m.parameters if p.kind in ("energy","width") and not p.fixed and p.free_index is not None)`.

Both vectors hold only the **free** parameters, in `.azr` order:
`p.free_index` is the position in that vector, `p.index` the position among all
parameters (`m.fixed_params`, `param.sav` lines).

### Data mode vs extrapolation mode

`m.data_mode()` (default) evaluates the `<segmentsData>` segments — this is what
χ² uses. `m.extrap_mode()` switches to the `<segmentsTest>` grids. Both
re-INITIALIZE every instance, so switch sparingly and batch the work.

**Segment indexing** (a frequent source of silent misalignment):

- data mode: segment `i` ↔ `m.datasets[i]` (the i-th `<segmentsData>` line);
  `m.nsegments` and `len(m.energies[i])` are its point count.
- extrap mode: **only ACTIVE test segments are returned**, in file order, so
  segment `i` ↔ `m.extrapolations.active[i]` — not `m.extrapolations[i]`.
  (7Be: 113 test segments declared, 108 active and returned.)

### Inspecting a model

```python
print(m.level_scheme)          # pairs → J-groups → levels → channels (read-only)
print(m.parameters.table())    # idx, name, kind, free, J^pi, E, L, S, pair, rad, value, wigner
print(m.datasets.table())      # per data segment: file, reaction, observable, E range, norm err
print(m.extrapolations.table())
m.pairs                        # PairSet: masses, charges, spins, sep/excitation E, channel radius
m.parameters.free / .widths / .energies / .norms / .shifts
m.parameters.by_physical_level()   # {LevelKey: ParameterSet} — one entry per resonance/pole
m.physical_levels()                # the LevelKeys, ordered (jgroup, level)
m.datasets.sys_errors(vary_only=True)   # per-segment normalization systematics (fractional)
```

A `LevelKey` prints as `5/2-#2@6.588MeV`; `(jgroup, level)` is its identity
(AZURE2 restarts level numbering inside every J-group).

### χ², residuals, and AZURE2's objective

```python
chi2   = np.sum(m.calculate_chi2_rwa(x))    # total DATA chi-squared
val, g = m.chi2_and_grad(x)                 # analytic gradient (energies/widths/norms)
r, J   = m.residual_jacobian(x)             # r_i standardized: sum(r**2) == chi2

m.segment_chi2(x)      # one entry per <segmentsData> segment; sums to chi2
m.dataset_chi2(x)      # {name: (chi2, points, segments)} -- per experiment
m.penalties(x)         # {"norm": array, "shift": array}, per segment
m.objective(x)         # chi2 + penalties: what AZURE2's own fit minimizes
```

`residual_jacobian` costs ~2 forward evaluations for the whole Jacobian
(reverse-mode adjoint) — that is what makes a Gauss-Newton / LM fit in Python
cheaper than Minuit's numerical gradients. Energy-shift columns come back zero;
an analytically unsupported segment raises.

**Fit against `objective`, not `calculate_chi2_rwa`.** The latter is the *data*
χ² only. `AZURECalc::operator()` also adds, per segment,
`((norm − nominal)/(nominal/100 · norm_error))²` and, for a free energy shift,
`((shift − nominal)/shift_error)²`. Minimize the bare χ² and the normalizations
drift to absorb every discrepancy — a "better" number AZURE2 would never have
found (−480 on the 7Be model). Note the denominator uses the **nominal**
normalization, and that `norm_error` is a **percentage**.

`m.objective(x)` is exactly what `chiSquared.out` reports as
`Total-Chi-Squared` + `Total-Norm-Chi-Squared`; `tests/pyazr/objective_test.py`
cross-checks it against the engine rather than against the formula.

For a least-squares fit that needs the penalties as residual *rows* (so LM sees
their Jacobian), append them with their constant derivatives:

```python
pen = [(p.free_index, d.norm, d.norm / 100.0 * d.norm_error)
       for p in m.parameters.norms
       for d in [m.datasets.by_key(p.segment_key)] if d.norm_error > 0]
Jpen = np.zeros((len(pen), nfree))
for k, (fi, _, s) in enumerate(pen): Jpen[k, fi] = 1.0 / s
resid = np.concatenate([r, [(z[fi] - n0) / s for fi, n0, s in pen]])
```

### Writing the run's output files

```python
m.write_output_files()        # AZUREOut_*, chiSquared.out, into output/
m.write_output_files(x_best)  # at a particular parameter vector
```

The files a colleague or the GUI's plot tab expects, without re-running the
binary. `output/` must already exist. In data mode this runs the χ²
evaluation first — `chiSquared.out` reports the per-segment χ² *stored on each
segment*, and only that evaluation sets it, so writing after a bare forward
pass gives a well-formed file full of zeros.

`param.par` and `parameters.out` do not come from this call: the first is
written at initialization, the second by `CNuc::PrintTransformParams` at the
end of a CLI run. Use `m.transform_rwa(x)` for the numbers `parameters.out`
would carry.

### Free-vector helpers

`Parameter.free_index` is a parameter's position in the free vector; threading
that by hand is the usual source of silent misalignment.

```python
m.n_rmatrix                       # where the R-matrix block ends: x[:m.n_rmatrix]
w = m.parameters.widths
w.indices()                       # their free-vector positions
w.take(x)                         # their entries of x
w.put(x, 0.0)                     # a copy of x with every width zeroed
m.datasets.by_key(p.segment_key)  # the segment a norm/shift belongs to
```

`by_key` rather than `datasets[key - 1]`: a segment key counts the *inactive*
segments too, so it is not an index.

### Calculating observables

```python
m.calculate_rwa(x)                  # cross section per segment (b or b/sr)
m.calculate_sfactor_rwa(x)          # S-factor per segment (MeV b)
m.calculate_energies(x)             # c.m. energies of the calculated points
m.calculate_excitation_energy(x)    # compound-nucleus Ex — the common axis across channels
m.calculate_angles(x)               # c.m. angles (differ from the lab angles in the .azr)
m.calculate_analyzing_power_rwa(x)  # vector A_y, for analyzing-power segments
m.calculate_angular_dists_rwa(x)    # Legendre coefficients, for ang-dist segments
```

Plot different reaction channels against **excitation energy**; it is the only
axis shared by all entrance pairs. `E_cm = Ex − threshold`.

### Analyzing power (observable 7)

The vector analyzing power `A_y` for a spin-1/2 projectile. Declare a segment
with `observable="analyzing-power"` (code 7 in the raw `.azr`); it is reported
**in place of** the cross section, so χ², fitting, plotting and output files
need no special handling.

```python
ay = m.calculate_analyzing_power_rwa(m.params_rwa)
```

Works for a spin-1/2 **projectile** on a target of any spin, in **both particle
and capture** exit channels — by two different routes.

*Particle exits* use Seyler's channel-spin amplitude matrix. The amplitudes come
out of the R-matrix in the channel-spin basis, so for a target with spin the
entrance index is decomposed into projectile and target projections before
`sigma_y` is applied to the projectile alone. For a spin-0 target (¹²C) the
channel spin *is* the projectile's and no decomposition is needed; for a
spin-1/2 target (¹⁵N) the channel spins are 0 and 1 and never 1/2.

*Capture exits* have no such amplitude matrix and use the Legendre coefficients
of Seyler & Weller, PRC 20 (1979) 453, Eqs. (20) and (21):
`A_y = sum_k b_k P_k^1(cos t) / sum_k a_k P_k(cos t)`, built by
`CNuc::CalcCaptureAnalyzingPower` (lazily, on first use — the 9-j symbols cost
time) and evaluated by `GenMatrixFunc::CalculateCaptureAnalyzingPower`. External
capture is included; a term may pair two internal pathways, two external ones or
one of each. The analytic adjoint is in `AMatrixFunc::PointAdjoint`.

The reason capture needs its own pathway table: `a_k` forces `s = s'`, so
`CalcAngularDists` only ever pairs pathways inside one KGroup. `b_k` does not,
and those channel-spin off-diagonal terms are exactly what the polarization
observable adds. `AZURE2_CAPPOL_DEBUG=1` prints, per point, the check that
`sum_k a_k P_k` reproduces AZURE2's own differential capture cross section, up
to `400*pi/(geom*I1I2)` — it agrees to machine precision. Reference case and the
validation against the paper's worked example: `tests/13N_capture_ay`.

Four things that will bite:

- **Angles are centre-of-mass** in the data file, unlike an ordinary
  differential segment. Columns are `E_lab · theta_cm · A_y · dA_y`.
- **A_y is a dimensionless ratio in [-1,1] and is often negative.** Do not put a
  relative uncertainty on it — a point near a zero crossing does not have a
  small uncertainty, and doing so gives those points enormous weight and wrecks
  the fit. Use a roughly constant absolute uncertainty.
- **Leave `vary_norm` off.** A normalization factor is meaningless for a ratio.
- **Use segments with no target integration** when comparing against
  thin-target data. `A_y` averaged over a target is cross-section weighted,
  `<A_y> = ∫A_y σ dE / ∫σ dE`, and since Rutherford σ diverges at low energy
  where A_y ≈ 0, a thick target drives the average to ~1e-6. That is physics,
  not a bug.

Analytic derivatives cover A_y in both channels, *except* for a point that also
carries target integration — that one reports itself unsupported, which drops
the analytic Jacobian for the whole fit back to numerical. Never wrong, just
slower. For capture there is a second reason to avoid target integration: the
geometric attenuation coefficients AZURE2 folds into `P_k` have no counterpart
for the `P_k^1` of the numerator.

Worked comparison against measured data: `tests/13N`, segments 11–16 (Baumann
1992). Capture: `tests/13N_capture_ay`. Formalism and implementation:
`docs/source/theory/polarization_*.rst`.

## Editing the model: adding and removing levels

AZURE2 reads its model from the file, so an edited scheme is applied by
**writing a new `.azr` and launching a fresh instance from it**. `AzrModel`
parses only the `<levels>` block and re-emits everything else verbatim; the
original file is never modified.

```python
from pyazr import AzrModel, azure2

mdl = AzrModel.from_file("7Be.azr")
print(mdl)                                        # J^pi / energies / channels / fixed flags
mdl.find(jpi="5/2+", energy=10.253, tol=2e-2)     # -> [AzrLevel]

mdl.remove_level(jpi="1/2+", energy=20)           # drop a background pole entirely
mdl.add_level(J=1.5, parity=+1, energy=8.6,       # add a 3/2+ resonance
              channels=[dict(pair=1, L=2, S=0.5, gamma=1000.0, fixed=False),
                        dict(pair=2, L=1, S=0.5, gamma=0.1)],
              level_fixed=False)                  # level_fixed=False -> energy is a fit parameter
mdl.deactivate_level(jpi="7/2-", energy=4.572)    # keep in file, zero+fix every gamma
path = mdl.write("_test.azr")

with azure2(path, cwd=HERE) as m:
    ...
```

Rules that will bite you:

- **All levels of one J^π share one channel set** (an R-matrix J-group
  requirement). If a level of that `(J, parity)` already exists, `add_level`
  clones that group's exact channel structure and applies your `gamma`/`fixed`
  only to matching `(pair, L, S)`; unmentioned channels are added inert
  (`gamma=0`, fixed). A spec matching none of them raises and lists the allowed
  channels. Only for a brand-new J^π do your channels define the group.
- For the same reason `remove_channel(jpi, pair, L, S)` removes that channel
  from **every level of the group**; removing it from one level makes AZURE2
  reject the file.
- `add_level`/`add_channel` can only reference a **pair that already exists**
  in the file — new particle pairs must be added in the GUI.
- `gamma` in the file is the **reduced width / ANC input**, not the eV partial
  width the GUI shows. Seed a new channel with a small nonzero value (fits from
  exactly 0 have zero gradient) and set `fixed=False` to free it.
- Levels are renumbered (`levelID`) automatically on every edit.

**Removing a level: file-level vs runtime.** Two different tools:

| | effect | use for |
|---|---|---|
| `AzrModel.remove_level` | level gone from the file; parameter vector shrinks | a persisted reduced model to refit |
| `AzrModel.deactivate_level` | stays in the file, all γ = 0 and fixed | persisted "level off" that keeps the J-group's channel set |
| `azure2.without_level(x, ...)` | returns a copy of `x` with that level's free reduced widths zeroed | *frozen* Δχ² probes, no refit, same parameter vector |
| `azure2.only_level(x, ...)` | zeroes every *other* level's widths | the bare contribution of one resonance |

Zeroing reduced widths decouples a level from every channel, so it contributes
exactly nothing — provided its *fixed* widths are already zero (check with
`m.parameters.widths`; in the 7Be model all fixed widths are 0).

Other `.azr` edits without touching the GUI:

```python
mdl.channel_radii()                                      # {pair: radius fm}
mdl.set_channel_radius(1, 5.0)                           # 3He+alpha, all its lines
mdl.set_segment_norm("Toth", vary=True, sys_error=8.0)   # percent, as stored
mdl.set_segment_active("Spiger", False)                  # drop a dataset from the fit
mdl.set_segment_datafile("Elwyn-F0012002", "data/new.dat")
mdl.add_data_segment("data/roughton.dat", entrance=1, exit=2,
                     observable="total-capture", energy_min=0.3,
                     energy_max=2.3, norm_error=5.0)     # brand-new dataset
mdl.remove_data_segments("artemov.dat")                  # or clear_data_segments()
mdl.set_extrapolations([...]) / add_extrapolation(...) / clear_extrapolations()
m.save_fit("fitted.azr", x_best)   # fit -> a .azr + its param.sav, verified
```

`add_data_segment` observables include `analyzing-power` (code 7), and
`pyazr/examples/exfor_fetch.py` + the `nds-explorer` skill show how to fetch
real datasets from EXFOR/NDS and feed them into `add_data_segment`.

**Adding/removing data segments invalidates the EC integral cache.** After any
`add_data_segment` / `remove_data_segments` / `clear_data_segments` /
`set_extrapolations` edit, delete `output/intEC.dat` and `output/intEC.extrap`
before the next run, or give the edited model its own output dir
(`mdl.set_output_dir(...)`); otherwise AZURE2 silently reuses integrals
computed for the *old* grids. In a live session,
`azr.recalculate_external_capture()` forces a recompute. See
`pyazr/examples/edit_model.py`.

**`m.save_fit(path, x=None)` is how you snapshot a fit.** It writes the `.azr`,
writes the companion `param.sav`, and verifies the result — reopening what it
wrote and comparing every R-matrix value against the fit. If they disagree it
removes both files and raises, so a snapshot you still have is one that reads
back as the fit it came from. `path` is always explicit; nothing is written in
place. Returns `(azr_path, sav_path)`.

```python
azr, sav = m.save_fit("7Be_fit.azr")          # current parameters
azr, sav = m.save_fit("7Be_fit.azr", x_best)  # or an explicit free vector
```

Three things it handles that used to be the caller's problem:

**The `<levels>` `gamma` field is NOT a reduced-width amplitude.** It holds the
physical value `parameters.out` prints: Γ in **eV** for an open particle
channel, an **ANC in fm^-1/2** for a closed (sub-threshold) one, Γ_γ in eV for a
photon channel — i.e. exactly `transform_rwa(x)`. In the 7Be model the two
differ by factors of 10² to 10⁷, and a file written with the rwa loads without
complaint and is wrong.

**A `Parameter`'s `pair` is not the file's pair key.** It is the *engine's*
number, which counts particle pairs in the order `<levels>` first mentions them.
On the 8Be model engine pair 1 is file key 2 and file key 1 is engine pair 6 —
match one against the other and every width lands on the wrong channel. Calling
`AzrModel.apply_fit` directly? Pass `pairs=m.pairs`, or it cannot translate.

**Normalizations and energy shifts are not in `<levels>`.** A calculate run on a
bare `.azr` uses 1.0 for every dataset, so its χ² sits *above* the fit's by
whatever they were absorbing (3He: 166 fitted, 769 from the snapshot). That is
what the companion `param.sav` carries — hand it to AZURE2 as the external
parameter file and the model is whole.

`AzrModel.apply_fit` is the lower-level half if you need it: it takes
`pairs=`, matches levels on the engine's own `(jgroup, level)` via
`engine_level_keys()`, and refuses (rather than silently skipping) any parameter
it cannot place. It does not verify — `save_fit` does that.

`pyazr/examples/save_fit_to_azr.py` does the whole thing — loads a fit from
`param.sav` or an `.npz`, optionally re-radiuses it (`--radius 1=4.70`, dropping
the stale `intEC` caches), writes the `.azr` plus a companion `param.sav` with
the norms, and fails loudly if the result does not round-trip.

### What a snapshot still cannot carry

`save_fit` verifies the `.azr` it writes, so a mismatched snapshot no longer
reaches you silently — it raises and removes the file. Two things remain true
of the `.azr` itself:

- **Normalizations do not live in `<levels>`.** A calculate run on a snapshot
  uses 1.0 for every dataset, so its `chiSquared.out` sits *above* the fit's by
  whatever the normalizations were absorbing — 3He: 166 fitted, 769 from the
  snapshot. Write them alongside (`output/normalizations.out`) or supply a
  `param.sav`.
- **The check dumps are keywords, not filenames.** `<config>` accepts only
  `none`, `screen` or `file` (`Config::ReadConfigFile`); anything else silently
  leaves the check off and `checks/` stays empty. Write `file`.

### Defining extrapolation grids

```python
mdl.set_extrapolations([
    dict(entrance=1, exit=-1, e_min=0.10/cm2lab, e_max=8.5/cm2lab,
         e_step=0.05/cm2lab, observable="total-capture"),
    dict(entrance=1, exit=1, e_min=0.05, e_max=6.0, e_step=0.05,
         observable="differential-cm", angle=90.0),
])
```

`exit=-1` means summed/total (capture). `observable` ∈ `angle-integrated`,
`differential`, `differential-cm`, `total-capture`, `angular-distribution`
(needs `order=`), `phase-shift` (needs `phase_J=`, `phase_L=`). Energies are
**lab** (`cm2lab = m_target/(m_beam+m_target)`; ³He+α 0.5703, p+⁶Li 0.8565 —
divide c.m. by it). Remember to delete `output/intEC.extrap` afterwards.

Note the `isDiff` codes differ between the two blocks: in `<segmentsData>` 3 is
total-capture and 4 is differential-cm; in `<segmentsTest>` 3 is angular
distribution, 4 total-capture, 5 differential-cm. pyazr handles this — hand-edits
must not.

## Fitting from Python

**pyazr does not fit.** It gives you the model, the χ², the objective and an
analytic Jacobian; the minimizer is yours. That is deliberate rather than a
gap: `residual_jacobian` returns the whole Jacobian for ~2 forward evaluations
(reverse-mode adjoint), which beats Minuit's numerical gradients on this
problem, and binding `MnMigrad`/`MnMinos` would mean mirroring AZURE2's
objective, parameter limits and Wigner bounds inside the API — two "AZURE2
fits" that can silently disagree. The CLI (modes 2 and 4) remains the reference
implementation.

Two things a hand-rolled fit must do to match AZURE2:

- **Minimize `m.objective(x)`**, not `calculate_chi2_rwa` — see the χ² section
  above. Bare χ² lets the normalizations absorb everything.
- **Respect the Wigner limits** if you care about them; `m.wigner_widths()`
  gives the bound per channel. The engine's `ParameterLimitsManager` enforces
  them, a Python fit will not unless you add them as bounds.

`scipy.optimize.least_squares` with the analytic Jacobian, fitting a **subset**
at frozen everything-else — the usual evaluation move:

```python
from scipy.optimize import least_squares

fit = np.array([p.free_index for p in m.parameters.widths
                if not p.fixed and p.radiation_type in ("E", "M")], int)

def resid(z):
    v = np.array(best, float); v[fit] = z
    return m.residual_jacobian(v)[0]

def jac(z):
    v = np.array(best, float); v[fit] = z
    return np.asarray(m.residual_jacobian(v)[1])[:, fit]

sol = least_squares(resid, best[fit], jac=jac, method="lm",
                    xtol=1e-14, ftol=1e-14, max_nfev=2000)
x = np.array(best, float); x[fit] = sol.x
```

For asymmetric (MINOS-style) errors, point `iminuit` at `m.objective` and
`m.chi2_and_grad` — it gives you `minos()` against the same objective the CLI
uses, without a second implementation of it inside AZURE2.

Fit in **rwa space** (that is where the Jacobian is exact) and convert for
reporting with `transform_rwa`. Free the relevant normalizations too
(`m.parameters.norms`; `m.datasets.by_key(p.segment_key)` is the segment) when
a subset fit would otherwise be absorbed by them.

Sanity-check the fit space: `transform_rwa(best[:nR])[p.free_index]` must equal
`p.value` for every width.

Save the result with `m.save_fit("fitted.azr", x)` — see the editing section.

For MCMC, `pyazr/examples/fit_emcee.py` and `fit_zeus.py` show the pattern —
one in-process engine per pool worker (constructed at module level), priors
chosen per `p.kind`, and
`log_prob = -0.5*(chi2 + offset) + log_prior` with
`offset = Σ log(2π σ²)`.

### Three traps that cost a whole evaluation

Found while building a BBN/pp-chain campaign; each produced a wrong answer
*silently*. The `evaluations/_tools/...` files cited below live in that analysis
repository, not in this one — the descriptions stand without them.

#### 1. `intEC.extrap` must be deleted before every extrapolation run

The golden rule above says to delete it "whenever the `<segmentsTest>` grid
changes". **That is not sufficient.** A session that initialises in data mode
and then calls `extrap_mode()` reuses whatever `intEC.extrap` is on disk, and
those amplitudes do not match the combined data+test integration set the second
INITIALIZE builds. Three identical runs of the same 7Be analysis script:

| run | caches | S₃₄(10 keV) | S(p,γ)(10 keV) | S(p,α)(10 keV) |
|---|---|---|---|---|
| 1 | deleted | 0.5402 keV b | 0.0912 keV b | 3254 keV b |
| 2 | deleted | 0.5402 keV b | 0.0912 keV b | 3254 keV b |
| 3 | reused | **1.055** keV b | **2.4e8** keV b | 3254 keV b |

The particle-exit grid is identical in all three — that is the signature, since
only capture reads these integrals. Rebuilding costs seconds. Put the deletion
in the script, not in your memory:

```python
for f in ("intEC.dat", "intEC.extrap"):
    p = os.path.join(eval_dir, "output", f)
    if os.path.exists(p): os.remove(p)
```

(`evaluations/_tools/uq.py:clear_ec_cache`.)

#### 2. A zero Jacobian column freezes `scipy.least_squares`

`residual_jacobian` returns an identically zero column for any parameter no
active dataset can see — a γ width in a model with no capture data, a channel
of a level that nothing populates. With `x_scale="jac"` that column's scale is
zero, and TRF reports convergence after **one** function evaluation while the
gradient at other parameters is still in the hundreds. On the 8Be model this
froze every normalization at exactly 1.0 for the entire run and left
χ²/N = 5.6; dropping the dead columns first gave χ²/N = 1.01 from the same
seeds.

```python
J = fitter.jacobian(x0)[:, cols]
cols = cols[np.max(np.abs(J), axis=0) > 0]      # before least_squares
```

A norm penalty of *exactly* 0.0 alongside free normalizations is the tell:
`m.penalties(x)["norm"].sum()` says so directly.

#### 3. `residual_jacobian` raises, and the exception aborts the fit

A trust-region step can put a reduced width where the Coulomb functions
overflow (`RuntimeError: z is not finite in log_Gamma`) or where the level
matrix is singular. Catch it and return a large residual with a zero Jacobian:
the optimiser then shrinks the region, which is what it would have concluded
from a finite but terrible χ². Letting it propagate loses the whole fit.

### Fitting a model built from hand-chosen seeds

A fresh `.azr` is far from any minimum and a single global fit over sixty
parameters walks into a bad one. Free the parameters in the order the data
constrain them, warm-starting each stage, and run the sequence two or three
times:

1. the particle widths of the pair that carries the dominant reaction;
2. the other particle pairs;
3. the γ widths, at fixed particle widths;
4. the background poles;
5. everything, including level energies and normalizations.

`evaluations/_tools/fitting.py` (`Fitter.stages`) implements this with the
penalty rows, the dead-column filter and the exception guard already in place.

## Decomposing a cross section

Everything below is done at the **fitted point, without refitting**, so each
number is that component's raw contribution.

```python
m.extrap_mode()
ex   = m.calculate_excitation_energy(best)     # common Ex axis
full = m.calculate_rwa(best)
off  = m.calculate_rwa(m.without_level(best, jpi="5/2-", energy=6.588))
only = m.calculate_rwa(m.only_level(best,    jpi="5/2-", energy=6.588))
bg   = m.calculate_rwa(zero_all_widths(best))  # non-resonant / hard-sphere / external
interference = (full - off) - (only - bg)
```

- `full − off` = everything the level does (resonant + its interference).
- `only − bg` = the bare resonance on the non-resonant background.
- The difference is pure interference. It is **block-diagonal in J^π** — only
  same-J^π levels interfere; that is a correctness check on any decomposition.
- Pairwise, level T against level L:
  `cross = only{T,L} − only{T} − only{L} + bg = 2 Re(A_T A_L*)`. Integrate
  `|cross|` over Ex to rank which partners a level interferes with.

**Capture-specific.** Zeroing only the **γ-channel** widths
(`p.radiation_type in ("E","M")`) leaves every particle channel — and therefore
the whole scattering solution — untouched, so it decomposes the capture
amplitude cleanly. Zero *all* γ widths at once and what remains is the
**external (direct) capture**.

**S-factor.** `calculate_sfactor_rwa` = cross section × an energy-only
conversion factor, so linear combinations of cross sections may be converted
after the fact: `conv = sfactor_full / xs_full` (guard the zeros) and multiply
each curve by it. Units: MeV b — ×10³ for keV b, ×10⁶ for eV b.

## Evaluation recipes

**Is a level needed?** Frozen Δχ² ranking — recompute χ² with each level in turn
switched off (`without_level`) at fixed parameters. Then, for the candidates,
refit with the level removed (`AzrModel.remove_level` / `deactivate_level` +
`least_squares`) — the frozen Δχ² is an upper bound on the level's importance,
the refit Δχ² is the honest one. Report both, plus which experiments' χ² move
(per-dataset partition above): a level justified by one dataset is a different
claim from one justified by all of them.

**Changing the channel radius.** `m.set_channel_radius(pair, r)` on a live
instance, or `AzrModel.set_channel_radius(pair, r)` + `write()` to persist one.
The radius is the matching surface, so AZURE2 rebuilds the compound nucleus and
data, redoes penetrabilities / shift functions / boundary conditions / Wigner
limits, and recomputes **every external-capture integral** (the API clears
`USE_PREVIOUS_INTEGRALS` for the rebuild, so `output/intEC*` is not reused).
A reduced width means something different at a different radius — **the model is
no longer fitted, so refit before reading anything off it.** For a scan, give
each radius its own directory (own `output/`, `data` symlinked) and walk the
ladder outward, warm-starting each radius from its converged neighbour rather
than from the original fit; starting them all from one point biases the ends.

**Is a new level/channel warranted?** Add it (`add_level`, or free one more
channel of an existing level), refit the subset it can affect, and weigh Δχ²
against the parameters spent and against physical bounds: θ² ≤ 1 for particle
channels, a few W.u. for γ channels. Sweep candidates in a loop; several sessions may be open at once.

**Dimensionless widths.** `m.dimensionless_widths(x)` (pyazr ≥ 2.4) does the
whole job — θ² for every particle channel, W.u. for every γ channel, at any
parameter vector:

```python
t = m.dimensionless_widths(best)          # -> WidthTable of ChannelWidth rows
print(t.photons.nonzero.table())          # .particles / .photons / .free / .nonzero
[c for c in t.particles if c.theta2 and c.theta2 > 1]
```

Each row carries `theta2` (= Γ/Γ_W), `theta2_formal` (= γ²/γ²_W),
`wigner_width` (Γ_W, eV), `wigner_gamma2` (γ²_W, MeV), and for photons
`e_gamma`, `weisskopf`, `wu`. Closed (sub-threshold) channels have `is_open =
False`, no Γ_W — AZURE2 reports an ANC there, not a width. `m.wigner_widths(x)`
returns just `{free_index: Γ_W}`; `weisskopf_width(rad, L, E_γ, A)` and
`m.mass_number` are exposed separately. Example:
`pyazr/examples/dimensionless_widths.py <azr> [--params output/param.sav]`.

The reason there are two answers: AZURE2 exposes the Wigner limit in two
non-interchangeable forms.

1. **GUI / width form**, in eV next to each channel: `Γ_W = 2 P_l(E_r) γ²_W`,
   energy-dependent. **θ² = Γ/Γ_W** — prefer this; it is checkable one channel
   at a time in the GUI.
2. **API form**, `Parameter.wigner_limit`: `γ²_W = ħ²/(μ a²)` in **MeV**,
   already squared and energy-independent. Then θ² = γ²_formal/γ²_W with the
   formal reduced width (`g_int` in `parameters.out`, squared).

Never divide by `wigner_limit` twice. The two differ by the level-shift factor
`1/(1 + Σ_c γ_c² dS_c/dE)`, negligible for narrow levels but a factor 70+ for
the 30 MeV background poles. AZURE2 reports **observed** partial widths,
`Γ_c = 2 P_c γ_c²/(1 + Σ γ² dS/dE)`. Some references use the Teichmann–Wigner
sum-rule unit `3ħ²/(2μa²)`, 1.5× larger (`teichmann_wigner()`).

`Γ_W` needs `P_l(E_r)`, which the API does not expose; `wigner_widths` gets it
out of AZURE2 itself by transforming a probe vector with every rwa set to a tiny
`eps`, so the shift denominator is 1 and `Γ_c(eps) = 2 P_c eps²`. One extra
`TRANSFORM_RWA`, at AZURE2's own radii and Coulomb functions — it reproduces the
GUI's Γ_W (7/2⁻ 4.57 ³He+α: 3.487e5 eV, θ² = 0.494).

**γ-ray strengths.** W.u. against the Weisskopf estimates E1 `6.8e-2 A^{2/3}
E³`, E2 `4.9e-8 A^{4/3} E⁵`, M1 `2.1e-2 E³` eV with E in MeV (E1–E4/M1–M4 in
`pyazr/widths.py`). Note Γ is **not** ∝ rwa² for a capture channel — the
external-capture term is linear in the amplitude — so convert with
`transform_rwa`, never by squaring.

## `.azr` file anatomy

Plain-text, section-tagged; prefer the GUI or `AzrModel` over hand edits.

- `<config>` — A-matrix flag, `output/` and `checks/` dirs, check toggles.
- `<levels>` — one line **per channel of each level**, 31 fields matching
  `NucLine` (`include/NucLine.h`); the file stores `2J`, `2S`, `2L` as integers.
  Blank line between levels; `levelID` groups them.
- `<segmentsData>` — one line per data segment: `isActive entranceKey exitKey
  minE maxE minA maxA isDiff [phaseJ phaseL] dataNorm varyNorm dataNormError
  [energyShift …] dataFile`. A `+10` on `isDiff` marks a THM/HOES segment.
- `<segmentsTest>` — extrapolation grids (see above).
- `<targetInt>` — target/experimental effects (integration, convolution).
- `<parameterSettings>` — free/fixed, limits, nuisance, category, Minuit index.
- `<mcmc>` — walkers, steps, threads.

## Data file format (`.dat`)

Four whitespace-delimited columns, **lab frame, forward kinematics**:
energy (MeV) · angle (deg) · cross section (b, or b/sr if differential) ·
uncertainty. The angle column is required even for angle-integrated data (dummy
value). Phase-shift data uses the same layout with the phase (deg) in column 3.
Analyzing-power data uses it with `A_y` (dimensionless, signed) in column 3 —
and its **angle column is centre-of-mass**, not lab.

## Output files (all center-of-mass)

- `AZUREOut_aa=<in>_R=<out>.out` — 9 cols: cm E, excitation E, cm angle, **fit**
  σ, fit S, **data** σ, data σ err, data S, data S err. For an analyzing-power
  segment, cols 4 and 6 hold `A_y` instead of σ (dimensionless, may be negative). `TOTAL_CAPTURE` in place
  of `R=<out>` for summed capture. `.band` files carry the covariance band.
- `AZUREOut_*.extrap` — 5 cols: cm E, excitation E, cm angle, σ, S (mode 3).
- `chiSquared.out` — per-segment χ²/N and norms; last line total χ². **The
  quickest scalar check that a run succeeded.**
- `param.par` initial / `param.sav` best-fit formal params (reload as the
  external parameter file); `parameters.out` physical/observable params.
- `normalizations.out` — fitted segment norms (auto-loaded with `param.sav`).
- `param.errors`, `covariance_matrix.out` — MINOS (mode 4).
- `reactionrates.dat` (mode 5); `samples.mcmc` (mode 6).
- `intEC.dat` / `intEC.extrap` — external-capture integral caches; see the
  golden rule above.

## Verifying a run

`cat output/chiSquared.out` — a finite total χ² and the expected N per segment.
From pyazr: `np.sum(m.calculate_chi2_rwa(m.params_rwa))` must reproduce the
model's known total (record it whenever the `.azr` changes; a jump concentrated
in the capture datasets almost always means a stale `intEC.extrap`). After a
fit, check `param.sav` updated and compare `parameters.out` widths — and their
θ² — against expectations.

## Examples shipped with pyazr

In `pyazr/examples/` — run `ls` there rather than trusting this list to stay
complete:

| | |
|---|---|
| model & scheme | `print_scheme.py`, `edit_scheme.py`, `edit_model.py`, `deactivate_level.py` |
| widths | `transform_widths.py`, `dimensionless_widths.py` |
| observables | `angular_distribution.py`, `sfactor_extrapolation.py`, `decompose_cross_section.py`, `reaction_rate.py` |
| fitting & UQ | `fit_emcee.py`, `fit_zeus.py`, `per_dataset_chi2.py`, `sensitivities.py`, `uncertainty_band.py` |
| snapshots | `save_fit_to_azr.py` |
| model internals | `coulomb_functions.py`, `ec_integrals.py`, `channel_radius_scan.py`, `nuclear_potential.py` |
| data | `exfor_fetch.py` |
