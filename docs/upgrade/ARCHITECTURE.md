# AZURE2 upgrade: target architecture

Status: agreed direction (Jakub Skowroński, 2026-11-07); not yet discussed with James deBoer's group in detail.
Companion files: `ROADMAP.md` (phases and tasks), `THM.md` (Trojan horse method), `RULES.md` (rules for `CLAUDE.md`),
`HANDOFF.md` (current state). Figures: `figures/structure.html` (open in a browser).

AZURE2 stays a C++ code built for speed. Python is used only to set projects up (library readers), for the web UI and
for users' own scripts.

## 1. The calculation in three lines

```
O_n(θ)     = observables at grid energy E_n from the collision matrix U(E_n; E_λ(θ), γ_λc(θ))
m_i(θ, ν)  = s_k(ν) · Σ_n K_in(ν) · O_n(θ)
χ²(θ, ν)   = (m − d)ᵀ C⁻¹ (m − d) + Σ priors(θ, ν)
```

| Symbol | What it is | Module | Computed |
|---|---|---|---|
| θ | fitted level parameters (observed energies, widths or ANCs) | inference | by the fitter |
| E_λ, γ_λc | formal level parameters in the one internal convention | physics | each evaluation |
| E_n | grid of compound-nucleus excitation energies covering every data point and its target/resolution window | experiment | once per fit stage |
| U(E_n) | collision matrix: one level-matrix solve per Jπ block per grid energy, shared by every reaction, angle and window | physics | each evaluation |
| O_n | observables at E_n: integrated σ, Legendre coefficients, nuclear amplitudes (for Coulomb interference), THM half-off-shell amplitudes | physics | each evaluation |
| K(ν) | apparatus kernel, sparse: target, stopping and straggling, resolution, Doppler, energy shift, angle weights, THM fold | experiment | when ν changes |
| ν, s_k | nuisances: normalizations, shifts, target density and stoichiometry, stopping scale, Doppler T | inference | by the fitter |
| d, C | data and covariance | experiment | once |

The gradient is the same three lines differentiated: ∂m/∂θ = s·K·∂O/∂θ, and ∂m/∂ν from K(ν). Nothing else.

**Exact energies (decision, Jakub, 2026-11-07).** As in AZURE2 today, every observable is calculated at the energy
where it is compared with data; nothing is interpolated by default. The grid E_n is therefore the set of evaluation
energies, not a mesh: (i) the data energies themselves, so a point without shift or fold has K_in = δ at its own
energy; (ii) for a free energy shift, the shifted energies, which move with ν — when the fitter changes a shift, the
channel functions of that measurement's energies are recomputed (they depend on the channel energy only, not on θ),
and ∂m/∂ν comes from S′, P′ and dO/dE at those energies; (iii) for a fold (target, resolution, Doppler, THM),
quadrature nodes placed around each level as `AdaptiveIntegrationGrid` does now, with the tolerance stated and
tested. An **interpolated mode** (observables on a fixed mesh, shifts and folds as interpolation weights) is an
opt-in speed-up for large data sets or sampling; it runs only with a check that reports its largest difference from
the exact mode on the project's own data, and T2 pins the two modes against each other.

## 2. Modules and the one-way rule

```
front ends (azure2 CLI, pyazr, web UI, user scripts)
   │                 │
   ▼                 ▼
inference/         io/            ← io builds physics and experiment objects from files and libraries
   │   ╲             │   ╲
   ▼    ╲            ▼    ╲
experiment/ ◄────────┘     ╲
   │                        ▼
   └──────────────────► physics/       ← uses nothing
```

- `physics/` uses nothing. `experiment/` uses physics (masses, kinematics). `inference/` uses physics and
  experiment. `io/` uses physics and experiment. Front ends use inference and io. `support/` (thread budget, per-task
  scratch, profiling) is used by all and uses nothing.
- Each module is its own CMake library target and links only the ones below it: a wrong dependency does not compile.
- Only the fitter writes θ and ν. During an evaluation everything else is read-only, so any number of threads is safe.

Why this is not spaghetti: spaghetti is (1) circular dependencies, (2) state written by several parts, (3) effects
not visible from the call. The one-way rule removes (1), the single feedback loop removes (2), and "each module
computes one symbol" removes (3). Today's `EPoint` holds the datum, its channel tables, its results and pointers to
other points; hosts write into mapped points; every thread deep-copies `EData` and shares parts back (copy-on-write,
node cache). That cycle through shared state produced the races and memory growth of 2026-08..11 (~55 % of ~126 dev
fixes touched `EData`/`EPoint`/`ESegment`; at least 21 fixes repaired bugs introduced by earlier perf commits).

## 3. Module contents (physical objects; physics stages are functions, not classes)

### physics/
- `System`: partitions α (masses, charges, spins, separation energies), channels c = (α, kind, l, s) with radius
  a_c, Jπ blocks. Immutable after loading.
- `JpiBlock`: arrays E[λ], γ[λ, c] (no Level/Channel objects with pointers), and the boundary value B_c of each of
  its channels: the engine sets B_c = S_c(E_1) per Jπ group, so one (α, l, s) carries a different B_c in each block
  (`physics_system.md`). Contiguous, cheap to copy.
- `channel_functions(α, l, E)` → P, S, φ, and Coulomb phase σ_l; from **pycoulomb** (Jakub's library, see §6).
  Indexed by (α, l), not by Jπ × channel (today's `EDependentTables::Slot(jGroup, channel)` stores the same numbers
  several times). Spin–orbit hybrid-potential channels are indexed (α, l, j).
- `level_parameters(θ, convention, basis)` → E_λ, γ_λc: observed ↔ formal (Brune, Park, none), angular-momentum and
  phase conventions (§5). Brune needs S_c(E_λ), which moves with θ: evaluated here, a few per evaluation.
- `collision_matrix(block, tables at E_n)` → A-matrix and U. One solve per Jπ block per grid energy.
- `external_capture`: radial integrals per grid energy, tabulated once per fit stage as unit amplitudes and scaled by
  the ANC or γ (they depend on the fit parameters only through that scale).
- `observables(U, geometry)` → O_n: integrated σ; Legendre coefficients B_L for angular distributions (coupling
  coefficients precomputed per reaction); nuclear amplitudes per (l, s, J) for elastic scattering, because Rutherford
  and Coulomb–nuclear interference are not a finite Legendre series (f_C(θ) enters through the kernel row); A_y as a
  numerator series folded separately from σ; phase shifts; total cross section (optical theorem); β-delayed spectra;
  THM half-off-shell amplitudes (`THM.md`); identical particles through symmetrised coefficients.
- `transfer/` (DWBA, for THM vertices): bound states, distorted waves, form factors (`THM.md` §3).
- `potentials/`: Woods–Saxon, Coulomb, global optical potentials (for DWBA and hybrid channels).
- `extrapolation`, `rates`: the same chain on a test grid; rates integrate σ on the grid and add the analytic
  narrow-resonance term for resonances narrower than the grid can follow.

### experiment/
- `Measurement`: data arrays (E, θ, y, δy) already converted to c.m. at load (frames are a load-time function),
  observable kind, reaction, apparatus settings, provenance.
- `Material`: elements and stoichiometry, density, phase (solid/gas), depth profile (implanted targets).
- `stopping/`: AZURE2's own electronic-stopping evaluation and straggling (§7).
- `EnergyGrid`: the evaluation energies E_n (§1, exact energies): data energies, shifted energies (updated when a
  shift moves; only that measurement's channel functions are recomputed), and fold quadrature nodes around the
  current levels (each width resolved), re-planned between stages. The interpolated mode adds a fixed mesh.
- Kernel builders, each returning sparse weights on E_n: target (thick/thin/implanted/gas, with straggling propagation),
  resolution (Gaussian, beam profile, UDR), Doppler (thermal, FGM), energy shift (constant, √E, f(E), groups: by
  default the shifted energies are evaluation points; in the interpolated mode, interpolation weights on the mesh), composites (SUM, RATIO, natural targets as
  rows), S-factor conversion (diagonal), THM fold (`THM.md`).

### inference/
- `ParameterVector`: names, physical ↔ fit transform, limits, free/fixed.
- `Objective`: one list of terms (data χ² per measurement, normalization and shift priors, stopping-scale priors,
  parameter priors, prior covariance); value, gradient, residuals, Jacobian. Used by the CLI fit, pyazr, outputs.
- `Fitter`: Minuit2 (MIGRAD/HESSE) and Levenberg–Marquardt on the residual Jacobian. Sampling (MCMC) is done by users'
  scripts on the objective's value and gradient (pyazr does not own fitting).

### io/
- `Project`: the single `.azr` reader/writer (today two exist: engine, pyazr; the Qt GUI's is gone, decision 3). Old files keep loading.
- Library readers live in pyazr `sources/` (Python, run at project setup; the C++ core reads only the project):
  AME/NUBASE (masses), ENSDF (level suggestions, tentative Jπ flagged), EXFOR X4 master files (measurements with
  units, frames, error split into covariance terms, DOI; unknown quantities refused), IAEA measured stopping (for the
  stopping evaluation). No network at run time; `pyazr libraries update` is the only networked command. Adopted
  records are written into the project with provenance (library, record, version, checksum).

## 4. What runs when

- Once per run: `System`, `Measurement`s, d, C.
- Once per fit stage (read-only afterwards): E_n, channel functions on E_n, angular coupling, capture unit
  amplitudes, THM vertex tables, kernel builders' fixed parts.
- When a free shift ν changes (before the evaluation, single-threaded hand-over): the shifted energies of that
  measurement and their channel functions are recomputed into a new table, which then replaces the old one; no table is
  written while an evaluation reads it.
- Each evaluation (depends only on θ, ν): level parameters → U(E_n) → O_n → m = s·K(ν)·O → χ².
- Threads split the grid energies into blocks; each block writes only its slice of O. Reductions in fixed order
  (results independent of thread count).

## 5. Angular-momentum and phase conventions

The core computes in one internal convention, written down once in the Sphinx docs: channel-spin coupling (l, s),
s = j1 + j2, Lane–Thomas phases. A `Convention` record in the System says what a project's parameters mean:
coupling scheme ((l, s) or jj), coupling order (j1 + j2 or j2 + j1), iˡ factors, sign of reduced widths, Coulomb-phase
reference, E/M multipole phases of photon channels. `level_parameters` converts at that boundary only: an orthogonal
block matrix of 6j symbols and signs on the γ columns, with its derivative. Test: every test project in every supported
convention gives identical observables to rounding (only parameter values change); a missing channel in a jj ↔ (l, s)
recoupling is reported. Care goes into interference signs (level–level, level–external capture), pinned by tests.

## 6. Coulomb functions: pycoulomb

The Coulomb functions are physics, not an external dependency: they live in `physics/coulomb/`, part of the
`physics/` target, taken from Jakub's **pycoulomb** C++ sources (github.com/skowrons94/pycoulomb, source commit
recorded in the directory's README; decision 6, 2026-11-07). One interface, `coulomb(l, η, ρ, method)` → F, F′, G,
G′, σ_l, and the Whittaker functions for bound channels; each algorithm one file behind it (`coulsm` = AZURE2 since
2025-09, the default; `coul` = AZURE2 before; `coulfg`; COULCC when its C++ version is done). `channel_functions`
builds P, S, φ on it. It replaces `external/coul`, `CoulFunc`, `CoulFuncCache`, `WhitFunc` and the hybrid-potential
code; `external/coul` is deleted. The GSL algorithm is not in `physics/` (no GSL there); it stays a cross-check in
the tests. The algorithms cross-check each other as an identity test.
Prerequisites in pycoulomb (owner: Jakub; see `ROADMAP.md` Phase 1): licence file, library target without Python,
stateless thread-safe call, Whittaker functions, no committed build products, tests against mpmath.

## 7. Stopping and straggling (no SRIM, no ERYA, no LunaSRIM)

- Today: a user formula, or ERYA, which ships SRIM's coefficient file `SRIM2013.xml` (SRIM's non-commercial terms,
  inside a GPL code). Straggling is c·√ΔE or plain Bohr. Both go.
- Electronic stopping, protons and α: AZURE2's own evaluation per element in the ICRU 49 functional form
  (velocity-proportional at low energy, Bethe with shell, Barkas and Bloch corrections at high energy), coefficients
  and covariance fitted to the IAEA measured-stopping database; gas and solid separately where data exist.
  d, t from protons and ³He from α at the same velocity; heavier ions by effective-charge scaling, with a larger
  stated uncertainty. Nuclear stopping: ZBL universal. Compounds: Bragg additivity plus a core-and-bond correction
  where measured compound data show one.
- Each evaluation carries its uncertainty, used as the prior of a stopping-scale nuisance per material.
- Straggling: Bohr with Chu's correction and the charge-state term (Yang, O'Connor, Wang 1991), propagated through
  the target, σ²(x) = S(E(x))² ∫ (dΩ²/dx′)/S(E(x′))² dx′ (Tschalär); beam spread in quadrature.
- User-measured tables or formulas override, with provenance.
- Before any third-party table is shipped (NIST PSTAR/ASTAR, IAEA data), check its redistribution terms; until then
  they are used for fitting and comparison only.

## 8. Physics explorer (web UI) — every physics function is plottable

Every physics function is a free function exposed in the binding with array inputs, so the web UI's Physics tab plots
it with no physics of its own: Coulomb F, G, F′, G′ versus ρ (all pycoulomb methods side by side), P_l, S_l, φ_l and
σ_l versus E, Whittaker functions, Wigner limits and single-particle widths, level-matrix poles and eigenphases, R and
U matrix elements, phase shifts, THM vertex functions W_l(p) and momentum distributions, DWBA bound states and
distorted waves, stopping and straggling versus E and depth, kernel rows K_i(E_n), and the grid E_n itself.

## 9. Extension points (one file, one symbol, one test project)

| You add | It changes | You write | Tested in |
|---|---|---|---|
| an observable | O | how it follows from U | the test project with that observable |
| an apparatus effect | K | its kernel weights on the grid | T2 |
| a parameterization or convention | E_λ, γ_λc | the map from θ, with derivative | all (observables unchanged) |
| a channel type | channel functions | P, S, φ for it (photon, β, potential) | T1 or T3 |
| a THM vertex model | O (THM) | the vertex function W(p) | T2 THM part |
| a prior | χ² | value and gradient | T2 |
| a data source | d, C or the System | a reader with provenance | a fixture from that library |
| a file keyword | the project | one field | round trip of the test projects |

## 10. Tests

Three synthetic test projects, generated from known parameters (a fit must recover them), small, seconds to run:

- **T1 capture and scattering, p + ¹²C → ¹³N**: elastic dσ/dΩ at lab angles, (p,γ) to ground and excited states with
  external capture E1/E2/M1, A_y, phase shifts, S factor, rate, Brune/Park/formal and convention round trips,
  extrapolation.
- **T2 reactions through an apparatus, p + ¹⁵N → ¹⁶O**: (p,α₀) integrated and differential; thick, implanted and gas
  targets with stopping and straggling; resolution, profile, Doppler; shifts (constant, √E, f(E), groups);
  normalization groups and priors; covariance; SUM and RATIO; natural target; THM data (²H(¹⁵N,α¹²C)n, HOES and OES
  points, the La Cognata 2009 case); all at 8 threads.
- **T3 neutrons, identical particles, special channels, ⁸Be**: α + α (identical bosons), p + ⁷Li, n + ⁷Be, β-delayed α,
  optical theorem, scattering length, hybrid potential channel.

Identity checks (no data): Coulomb functions vs mpmath and Wronskian; U unitary with elastic only; optical theorem;
isolated narrow level → Breit–Wigner area; constant σ survives every kernel; Gaussian∘Gaussian; Brune ↔ formal ↔ Park
round trip; convention invariance; gradient = finite differences; thread-count independence; stopping thin-layer limit,
Bohr limit of straggling, probability conservation of kernel weights; THM Green identity and DWBA vs FRESCO
(`THM.md`).

During migration today's engine provides golden outputs for T1–T3; the ~96 current test folders stay as oracles until
the phase that replaces their code is merged.

## 11. Front ends

- **pyazr 3.0**: modules named like the core (`project`, `system`, `data`, `objective`, `evaluate`, `physics` (the
  explorer functions), `sources`, `bands`, `plot`, `ui`, `cli`); read-only numpy views of setup data; no fitting
  framework (recipes for iminuit, scipy, emcee in the docs); no physics in Python except library parsing and
  presentation (any Python formula tested against the core); commands `pyazr ui | import exfor|ensdf | stopping |
  libraries update`; wheels (scikit-build-core + cibuildwheel) for Linux/macOS/Windows carrying the core and the
  `azure2` CLI, replacing the Qt installers; today's `azure2()` session kept one release as a warning shim.
- **Web UI** (replaces Qt): `pyazr ui project.azr`, FastAPI + uvicorn (`pyazr[web]`), one static vanilla-JS page with
  plotly.js vendored for offline use (the stack of Jakub's LunaSRIM web app, without depending on it); thin `/api`
  mapping one-to-one onto pyazr objects; long jobs in worker processes with progress and stop; the project is the state
  (saved as `.azr`, undo = history). Tabs: Project · Nucleus · Data · Fit · Results · Physics (§8) · THM (only when the
  project enables it). Tested with pytest on `/api` on T1–T3, a few browser smoke tests, and the Qt parity checklist
  (`docs/upgrade/qt_parity.md`).
- **CLI**: `azure2` takes flags and a run file; no readline menus.

## 12. Out of the tree at the end

Qt/Qwt GUI (removed already, decision 3), ERYA (and `SRIM2013.xml`), numcmc, readline menus, the GSL fitter, per-thread `EData`/`CNuc` clones,
point mapping, copy-on-write tables, `NodeTableCache`, process-wide memos, `external/coul` (→ pycoulomb). Kept
externals: Minuit2, GSL (quadrature and special functions only), pybind11, pycoulomb.
