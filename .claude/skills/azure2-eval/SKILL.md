---
name: azure2-eval
description: Run AZURE2 R-matrix evaluations of nuclear reaction/scattering data — calculate cross sections & chi-squared, fit levels/widths, add/remove levels, decompose a cross section into level/interference/external contributions, extrapolate, compute reaction rates, or drive AZURE2 from Python via pyazr. Use whenever the task involves an .azr project file, AZURE2 levels/channels/segments, R-matrix fitting, S-factors, resonance significance tests, or the o17_guardo (THM) example in this repo.
---

# Running AZURE2 evaluations

AZURE2 is a multi-channel, multi-level R-matrix code. An evaluation is driven by
a single `.azr` project file that describes the compound nucleus, particle
pairs, levels/channels, and data segments. Two headless ways to run it: the
**interactive CLI** (`--no-gui`) and the **pyazr** Python API. Use the CLI for a
one-shot calculate/fit/extrapolate/rate that writes the standard output files;
use **pyazr for all evaluation work** — χ², parameter scans, level add/remove
tests, cross-section decomposition, custom fitters and samplers.

The binary is on `$PATH` as `AZURE2` (also `build/src/AZURE2` in the source
repo); pyazr resolves it via `$AZURE2_BINARY` → `<repo>/build/src/AZURE2` →
`$PATH`, so `binary=` is usually unnecessary. **The two can drift**:
`/usr/local/bin/AZURE2` is a copy, not a symlink, so a `make AZURE2` in the
source repo does not reach it — `sudo cp build/src/AZURE2 /usr/local/bin/` after
every rebuild, or export `$AZURE2_BINARY`. `pyazr/` lives at the repo root
(`/Users/kuba/Desktop/R-Matrix/Codes/AZURE2/pyazr`, v2.3.0) and is mirrored into
evaluation directories, so `from pyazr import ...` works when cwd is either.

Full reference docs: `docs/_build/html/_sources/` — `reference/` covers
command_line, data_formats, output_files; `user_guide/` covers levels_channels,
segments, particle_pairs, fitting, mcmc. Worked example: `examples/o17_guardo/`.

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
- **One AZURE2 session per process.** Opening a second `azure2()` in the same
  interpreter desyncs `GET_PARAMS_INFO`. Sweep variants from a shell loop, or
  build one temp `.azr` that carries everything the run needs.
- CLI mode does **not** read Runtime Options from the `.azr` — pass them as
  flags every time (`--use-brune`, `--gsl-coul`, …).

## Workflow A — interactive CLI (one-shot runs)

Prompt order: **(1) menu choice → (2) external parameter file → (3) mode-specific
prompts**. Pipe stdin to run non-interactively.

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
# Calculate with saved best-fit params; writes chiSquared.out + AZUREOut_*.out
printf '1\noutput/param.sav\n' | AZURE2 --no-gui 7Be.azr

# Fit fresh from the .azr levels. Mode 2 then asks about the cross-section
# uncertainty band (y/n) and, if yes, reduced-chi2 scaling (y/n):
printf '2\n\nn\n' | AZURE2 --no-gui 7Be.azr

# Extrapolate (no data) using saved params:
printf '3\noutput/param.sav\nn\n' | AZURE2 --no-gui 7Be.azr
```

Non-interactive band control: `--covariance-band` (+ `--scale-covariance`) skips
the y/n prompts. Other flags: `--use-brune`, `--gsl-coul`, `--ignore-externals`,
`--use-rmc`, `--no-transform`; `--help` lists all.

## Workflow B — pyazr

pyazr spawns headless `AZURE2 --no-gui --use-api` processes and talks to them
over a socket. A session is a context manager; always close it.

Install it from the repository root (it is not on PyPI):

```bash
pip install -e .                # core: numpy, mpmath, scipy
pip install -e ".[examples]"    # + matplotlib, emcee, multiprocess
pip install -e ".[all]"         # + zeus-mcmc
```

`-e` is the useful form, since the package lives in the checkout. AZURE2 itself
is C++ and is **not** installed by pip — build it separately. pyazr locates the
binary via `$AZURE2_BINARY`, then `build/src/AZURE2`, then `$PATH`, and it must
be built with `USE_API=ON` (the default). Installing also puts `pyazr-cleanup`
on the path, which reaps API instances orphaned by an interpreter that died
before it could close them.

```python
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")   # set BEFORE importing numpy
import numpy as np
from pyazr import azure2, AzrModel

with azure2("7Be.azr", nprocs=1, cwd=HERE) as m:
    best = np.asarray(m.params_rwa, float)      # free parameter vector
    chi2 = np.sum(m.calculate_chi2_rwa(best))
```

`nprocs=N` spawns N independent instances (`proc=i` selects one) — for one χ²
per walker in emcee/zeus. Ports are OS-assigned, so parallel sessions never
collide.

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

### χ², residuals, per-dataset χ²

```python
chi2  = np.sum(m.calculate_chi2_rwa(x))     # total only
val, g = m.chi2_and_grad(x)                 # analytic gradient (energies/widths/norms)
r, J   = m.residual_jacobian(x)             # r_i standardized: sum(r**2) == chi2
```

`residual_jacobian` costs ~2 forward evaluations for the whole Jacobian
(reverse-mode adjoint) — use it for Gauss-Newton/LM fits *and* to get the
**per-segment χ² the API does not expose**:

```python
seglens = [len(m.energies[i]) for i in range(m.nsegments)]   # data mode
idx = np.cumsum([0] + seglens)
seg_chi2 = np.array([np.sum(r[idx[i]:idx[i+1]]**2) for i in range(m.nsegments)])
# aggregate by experiment: m.datasets[i].name
```

Energy-shift columns come back zero; an analytically unsupported segment raises.

**`chi2`/`residual_jacobian` are the DATA χ² only.** AZURE2's own fit objective
(`AZURECalc::operator()`) also carries a penalty per free normalization,
`((norm − nominal)/(nominal·sys_err))²`, and one per free energy shift. Minimize
the API's residuals alone and the norms drift to absorb every discrepancy — a
"better" χ² that AZURE2 would never have found (−480 on the 7Be model). Append
the penalty rows, with their constant Jacobian rows:

```python
pen = [(p.free_index, d.norm, d.norm * d.norm_error / 100.0)
       for p in m.parameters.norms
       for d in [m.datasets[p.segment_key - 1]] if d.norm_error > 0]
Jpen = np.zeros((len(pen), nfree))
for k, (fi, _, s) in enumerate(pen): Jpen[k, fi] = 1.0 / s
resid = np.concatenate([r, [(z[fi] - n0)/s for fi, n0, s in pen]])
```

`chiSquared.out`'s `Total-Chi-Squared` is likewise the data part; the penalty is
the separate `Total-Norm-Chi-Squared`.

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

Analytic derivatives cover A_y, *except* for a point that also carries target
integration — that one reports itself unsupported, which drops the analytic
Jacobian for the whole fit back to numerical. Never wrong, just slower.

Worked comparison against measured data: `tests/13N`, segments 11–16 (Baumann
1992). Formalism and implementation: `docs/source/theory/polarization_*.rst`.

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
mdl.set_extrapolations([...]) / add_extrapolation(...) / clear_extrapolations()
mdl.apply_fit(m.parameters, x_best, transform=m.transform_rwa)   # fit -> <levels>
```

`apply_fit` turns a fit result into a real `.azr` snapshot — the reference to
reuse and hand to the GUI. Matching is by (2J, parity, input energy, pair/L/S).

**The `<levels>` `gamma` field is NOT a reduced-width amplitude.** It holds the
physical value `parameters.out` prints: Γ in **eV** for an open particle
channel, an **ANC in fm^-1/2** for a closed (sub-threshold) one, Γ_γ in eV for a
photon channel — i.e. exactly `transform_rwa(x)`. In the 7Be model the two
differ by factors of 10² to 10⁷, and a file written with the rwa loads without
complaint and is wrong. So `apply_fit` requires the conversion to be explicit —
`transform=m.transform_rwa`, or `physical=True` if you converted already — and
raises if given neither. Level energies need no conversion.

`apply_fit` covers `<levels>` only: normalizations do not live there, so a fit
that moved them needs its `param.sav` alongside. Verify a written file by
reloading it and comparing `transform_rwa(m.params_rwa)` against the fit's
physical vector — every R-matrix entry must match.

`pyazr/examples/save_fit_to_azr.py` does the whole thing — loads a fit from
`param.sav` or an `.npz`, optionally re-radiuses it (`--radius 1=4.70`, dropping
the stale `intEC` caches), writes the `.azr` plus a companion `param.sav` with
the norms, and fails loudly if the result does not round-trip.

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

Minuit (CLI mode 2) fits everything free in the `.azr`. For evaluation work you
usually want to fit a **subset** at frozen everything-else, which pyazr does with
`scipy.optimize.least_squares` and the analytic Jacobian:

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

Fit in **rwa space** (that is where the Jacobian is exact) and convert for
reporting with `transform_rwa`. Free the relevant normalizations too
(`m.parameters.norms`, `p.segment_key - 1` indexes `m.datasets`) when a subset
fit would otherwise be absorbed by them.

Sanity-check the fit space: `transform_rwa(best[:nR])[p.free_index]` must equal
`p.value` for every width.

For MCMC, `pyazr/examples/fit_emcee.py` and `fit_zeus.py` show the pattern —
`nprocs` instances, one per pool worker, priors chosen per `p.kind`, and
`log_prob = -0.5*(chi2 + offset) + log_prior` with
`offset = Σ log(2π σ²)`.

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
channels, a few W.u. for γ channels. Sweep candidates from a shell loop — one
session per process.

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

In `pyazr/examples/`: `angular_distribution.py`, `print_scheme.py`,
`edit_scheme.py`, `deactivate_level.py`, `transform_widths.py`,
`dimensionless_widths.py`, `save_fit_to_azr.py`, `uncertainty_band.py`,
`fit_emcee.py`, `fit_zeus.py`.
