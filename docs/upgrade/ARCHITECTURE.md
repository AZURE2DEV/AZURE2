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

As built (Phase 5, second branch, `src/GradEnergyTasks.cpp`): the gradient runs over the same energy tasks as the
point pass. The point pass gives every unit's record and every point's model; each point's cotangent goes back
through its expression tree and its fold to cotangents of its units (the transpose of K); per energy task the
level matrix is filled and factored once, the T-matrix cotangents of the units of each T-matrix group are summed,
and the bilinear terms g_c^T A g_c' are back-propagated once per J-group (∂O/∂θ, transposed). The sums are formed
in an order fixed by the tasks, not by the threads (64 interleaved chunks of tasks, each summed in task order, then
the chunks in order; rows task by task), so a gradient and a fit give the same bits at any thread count. Why:
floating-point addition is not associative, and a fit amplifies a difference of 1e-16 in one gradient into a different
path (seen: two 3-thread runs of the same T2 fit ending apart before this; a 30-call natural_neon fit moving 6e-5 in a
parameter from a 1e-15 change of the summation order).

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
computes one symbol" removes (3). Until 2026-11 one data point (`EPoint`) held the datum, its channel tables, its
results and pointers to other points; hosts wrote into mapped points; every thread deep-copied the data (`EData`) and
shared parts back.  That cycle through shared state produced the races and memory growth of 2026-08..11 (~55 % of ~126
dev fixes touched `EData`/`EPoint`/`ESegment`; at least 21 fixes repaired bugs introduced by earlier perf commits).
The data path now has three types (`data_structure.md`): what was measured (`experiment::Measurement`), the data set
as the set-up leaves it (`DataSet`, `include/DataSet.h`: segments with their settings and terms, target effects,
covariance groups, parameter layout; immutable, one object shared by every working copy through a
`shared_ptr<const DataSet>`), and the working copy (`WorkingCopy`, `include/WorkingCopy.h`: per segment a
`SegmentState` of norms, shifts, frames and shift tables, the present values of the target-effect parameters, the
calculation plan and one flat vector of calculation `Unit`s; copied by a plain copy).  An evaluation writes only its
working copy and the caller's `Results`.  Its point pass (`WorkingCopy::CalculatePointsInParallel`) runs the plan's
energy tasks (`CalculationPlan::EnergyTasks`, since Phase 5): one task per compound energy E_x (nucleus and the bits of
`AMatrixFunc::InEnergy`) with every unit at it -- points without a fold and sub-points alike -- one level-matrix fill
and factorisation, one T-matrix per (entrance, exit) pair and external-capture amplitude source, one observable per unit (`azure2::SolveAtOneEnergy`); then
the folds, on values already computed (`Unit::FoldSubPoints`).  Each task writes only its units' records.  The plan
has one task per point (the point mapping, which made one point calculate others at its energy, is deleted).  Beside
the units stands the stage's node table (`include/NodeTable.h`: one node per nucleus and bits of E_x with its channel
rows, built at the end of the set-up, shared by every copy) with, per plan task, one slot per unit naming its node;
nodes the hand-over moves off the stage table (re-placed grids, shifted points) live in the working copy's own overflow
table (`NodeOverflow`, placed by `WorkingCopy::PlaceOverflowNodes`).  The energy tasks are grouped by the slots' nodes;
the slots of the folded points (`Site`s) carry the rest of the calculation, and `AZURE_NODE_CHECK=1` shows at every
hand-over and point pass that every unit is at the node its slot names.  Folded points on the shared mesh share their
slots (`include/SlotGroup.h`): one site per node and sharing group, each point a range of its group's block, the slot
written by its first reader (`WorkingCopy::WritesSlot`: a plan task's units are its point's written slots) and read by
all; `AZURE_NODE_CHECK` checks the groups too (`WorkingCopy::CheckSlotGroups`).

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
- `transfer/` (DWBA, for THM vertices; `physics/transfer`, `transfer.md`): bound states, distorted waves, form
  factors (`THM.md` §3), the zero-range THM weight and the global optical sets (one table, `GlobalOptical.h`). The
  engine's `src/ThmDistortion.cpp` (R(E)) and `src/ThmDwVertex.cpp` (its waves) call it; `src/ThmOptical.cpp` is
  gone.
- `potentials/`: Woods–Saxon, Coulomb, global optical potentials (for DWBA and hybrid channels).
- `extrapolation`, `rates`: the same chain on a test grid; rates integrate σ on the grid and add the analytic
  narrow-resonance term for resonances narrower than the grid can follow.

### experiment/
- `Measurement`: data arrays (E, θ, y, δy) as read, in the frame they were given in; the model is compared with them in
  that frame, r = (f/J − n y)/(n δy), the frame of a row (`MeasurementFrame`) holding only its kinematics and J, the
  **model expression** that predicts them (below), apparatus settings, provenance.
- **Model expressions (decision 12, Jakub, 2026-11-07): what a data set measures is written, not coded.** A
  measurement's prediction is a small expression tree, so a new kind of data is a new combination, not new engine code.
  - Leaves are `Observable` selectors: a reaction (entrance → exit partition, exit excitation), an observable kind
    (σ, dσ/dΩ, Legendre coefficient a_L, A_y, phase shift, S factor, β-delayed spectrum, THM yield), and an optional
    **amplitude filter** — multipolarity (E1, E2, M1, …), Jπ, l, channel spin, internal or external part. A filter
    selects amplitudes *before* they are squared: "E1 only" is |Σ_{E1} T|², with the E1 interferences kept and the
    E1–E2 cross term excluded; a sum of filtered leaves is therefore not the unfiltered observable, and the docs and
    tests say so.
  - Nodes combine leaves: `Sum` (incoherent, e.g. several final states or γ lines), `Ratio`, `Scale` (fixed factors,
    branching ratios), `Weighted` (isotopic abundances of a natural target), and the apparatus folds `Yield` (target
    thickness and stopping), `Resolution`, `Doppler`, `Shift`, each a kernel on its sub-expression.
  - The same tree gives the value and, by the chain rule through each node, the gradient; each node is one class with
    `value` and `derivative`, tested alone.
  - A project is format 3 only: its records (data sets, extrapolations, targets, beams, detectors, model data) are
    read by `io/` and built by `experiment/`'s table into these trees; a new data type needs at most one new leaf
    filter or node. An older `.azr` (segment lines, isDiff codes, SUM/RATIO composites) is refused with one message
    (`io::NotFormat3Message`, exit 3) and converted by `convert/` (`AZURE2 --convert`), the only code that knows the
    old grammar; the engine library links nothing of it (repo_consistency check 34).
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

## 7. Stopping and straggling: SRIM-2013 through LunaSRIM (no ERYA)

- Before: a user formula, or ERYA, which ships SRIM's coefficient file `SRIM2013.xml` inside a GPL code and evaluates
  it with one formula at all energies. Straggling was c·√ΔE or plain Bohr. ERYA is removed (branch
  `stopping-wire`, 2026-11-11), where a target's `stopping evaluation` calls `stopping/`.
- **Decision of the PI, 2026-11-11:** "I would use the SRIM stopping powers from my LunaSRIM repository. They are free
  to use as long as SRIM is cited." This replaces the earlier plan of an own evaluation fitted to the IAEA database.
  So the stopping library `stopping/` (`docs/upgrade/stopping.md`) computes SRIM-2013's stopping powers, ported to
  C++ from LunaSRIM (github.com/skowrons94/LunaSRIM, MIT; a re-implementation of SRIM-2013 that reproduces its
  tables), with SRIM's data files `SCOEF03.DAT` and `SNUC03.dat` (`stopping/data/srim/`, not under AZURE2's licence;
  whether the public repository ships them is the owner's decision, the library can read them from a directory):
  electronic stopping of any ion Z₁ = 1–92 in any element, solid or gas with SRIM's phase factor; SRIM's
  pair-specific nuclear stopping; compounds by Bragg's rule with SRIM's compound correction (a number the user
  gives). Agreement with LunaSRIM: to rounding (10⁻¹⁰ in the tests). SRIM is cited in the data directory, in the
  provenance text and in every output file of `azure2-stopping`.
- Each evaluation carries its uncertainty — SRIM's published accuracy for the ion (4.0 % H, 3.9 % He, 4.8 % Li,
  5.8 % heavier) — used as the prior of a stopping-scale nuisance per material.
- Straggling: TRIM's (Bohr with TRIM's relativistic factor and low-energy limit), plain Bohr selectable, propagated
  through the target, σ²(x) = S(E(x))² ∫ (dΩ²/dx′)/S(E(x′))² dx′ (Tschalär); beam spread in quadrature.
- User-measured tables or formulas override, with provenance.
- Before any other third-party table is shipped (NIST PSTAR/ASTAR, IAEA data), check its redistribution terms; until
  then they are used for comparison only.

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
| an apparatus effect | K | its kernel: a fold template (value and weights) and a support, one row of `Kernels()` (`experiment/src/Kernels.cpp`) | `experiment/tests/kernels_test.cpp`, T2 |
| a parameterization or convention | E_λ, γ_λc | the map from θ, with derivative | all (observables unchanged) |
| a channel type | channel functions | P, S, φ for it (photon, β, potential) | T1 or T3 |
| a THM vertex model | O (THM) | the vertex function W(p) | T2 THM part |
| a prior | χ² | value and gradient | T2 |
| a data source | d, C or the System | a reader with provenance | a fixture from that library |
| a file keyword | the project | one field | round trip of the test projects |
| a kind of data (sum, ratio, yield, E1 part, …) | a model expression | one leaf filter or node, if not a combination of existing ones | the test project with that data |

**State.** Two kinds only.

- *Constant tables*, built once and never written afterwards: masses, constants, element symbols, the format-3 schema
  and its quantity rows, the list of χ² terms. Each is reached through a function returning a `const` reference. This
  is the only form of singleton in the code.
- *Everything an evaluation writes* belongs to the caller and is passed in: the `Results` buffer, the per-stage
  tables, per-task scratch. No function-level `static` or `thread_local` holds a value that a calculation reads back,
  and no cache is keyed by "the last call". Such state made calculations depend on the thread count (the Coulomb memo,
  removed in Phase 1).

**State: the rule, the allow-list, the check** (Phase 6, `evaluation.md`).  The only process-wide objects are
immutable tables built before the first evaluation; every other state is owned by the caller and passed down.  Where
the code keeps something between calls, it is one of three owned kinds:
- the energy task's scratch, `azure2::TaskScratch` (`include/TaskScratch.h`), one per thread in the caller's
  `Results::scratch`: the level-matrix solvers with their buffers and the channel-space route's per-J-group values;
- the run's scratch, `azure2::RunScratch` (`include/RunScratch.h`), held by the run's `Config` (`Config::scratch`):
  per thread the external-capture integrals' radial values and workspace and the THM vertex boundary's S_c; shared by
  the run's threads behind their locks the external-capture term table (`ECTermStore`), Park's shift-derivative memo
  and the `--no-ec-table` amplitude store;
- the run's once-only warnings, `Config::printedOnce` (`Config::FirstTime(key)`), so a second run in the same process
  prints them again.
Each memo in them is keyed by everything its value depends on, so a hit returns what the computation gives.  The
settings a channel function reads deep inside a header (relativistic kinematics, the Coulomb table) are given to each
particle pair when the run builds it (`PPair::FunctionSettings`); no process-wide pointer to a `Config` exists.  Two
sessions may evaluate at the same time on two threads (`tests/pyazr/concurrent_sessions_test.py`).  The
check (`tests/consistency/repo_consistency_test.py` check 38, scanner `tests/consistency/state_scan.py`) refuses any
`thread_local`, non-`const` static or namespace-scope variable and `mutable` member of the engine's sources that is
not a row of `tests/consistency/allowed_state.tsv` (file, symbol, kind, reason); a row whose symbol is gone fails too.
Rows marked `table:` stay (the Coulomb tables read from files, GSL's error handler, locks of process-wide streams);
rows marked `pending (milestone N)` name the milestone of `evaluation.md` §5 that removes them.  The test that two
calculations in one process do not see each other: `tests/pyazr/session_isolation_test.py`.

**Extension.** Each row of the table above is an explicit table in one source file, whose rows hold the name, what
the front ends need (units, fields, documentation) and the function that does the work: the quantity table of
`experiment/` (`project_format3.md` §7), the χ² terms of `src/Objective.cpp`, the apparatus kernels (`Kernels()`,
`experiment/src/Kernels.cpp`, a row per kernel file; a detector-resolution row also generates its format-3 form) and
the stopping models of Phase 5. A new feature is one row and one test; `--schema`, the documentation tables and the GUI forms are
generated from the same rows, and a consistency check fails when a row is bypassed. Not used, on purpose:

- classes that register themselves from a static object in their own `.cpp` (initialisation order is undefined, a
  static library can drop them at link time, and no file lists what exists);
- abstract-factory hierarchies and injection frameworks (indirection that removes no duplicated formula).

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

**Comparing the outputs of two engines.** A change that must not move a number (a move of code, a split of a file, a
new reader of the same input) is checked by running every project of the test tree with both engines and comparing
every file they write. Three commands, from the repository root:

```bash
export SRC_TESTS=$PWD/tests WORK=<scratch> JOBS=<runs at once>
tools/runall.sh <reference engine, absolute path> tab-ref && tools/runall.sh <new engine, absolute path> tab-new
tools/cls.sh tab-ref tab-new; grep '^## ' $WORK/cls-tab-new.tsv
```

- `tools/runall.sh ENGINE TAG [MODE ...]` runs each folder of `$SRC_TESTS` that holds an `.azr` (`output/` folders
  left out) as one job of `tools/run_one.sh`: menu 1 (calculate with data), one OpenMP thread, in a private copy of
  the folder's top-level test folder, in the three bases ("" Brune, `--no-brune`, `--use-park`), `JOBS` jobs at once
  (default 1; the caller states its thread budget). A run is cut off after 900 s, its log after 2 MB. The files land in
  `$WORK/cmp/<basis>/<TAG>/.res-<folder>/<file>.azr/`: `output/`, `checks/` and `run.log` (`<basis>` is `brune`,
  `no-brune` or `use-park`; `/` in the folder is written `__`). Every run gets a fresh tag.
- `tools/cls.sh REFTAG NEWTAG` writes `$WORK/cls-NEWTAG.tsv` with `tools/classify.py`: one line per file that differs
  or exists in one run only (largest relative difference of the numbers, counts of numbers, sign flips and words that
  differ, change of the line count), then per basis `## <basis>: N files, M identical` and a line per kind of file
  with a difference. *Identical* means: each file is read as lines of tokens (carriage returns are line breaks, commas
  blanks, `=` a token of its own, split at blanks; the lines of `run.log` with a progress bar, a percentage,
  "Calculating" or "Progress" dropped), and the two files have the same number of lines, the same number of tokens on
  each line, and each pair of tokens is equal as text or, when both are numbers, equal in value (`1.0` and `1`).
  Spacing and the spelling of a number do not count; any change of a value, however small, does.
- The converted-tree gate of format 3 (`tools/format3_gate.sh`, `project_format3.md` section 8.5) runs the same
  `runall.sh` on a converted copy of the tree and compares with its own rules (`tools/format3_gate.py`).

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

Qt/Qwt GUI (removed already, decision 3), ERYA (and `SRIM2013.xml`), numcmc, readline menus, the GSL fitter, per-thread `WorkingCopy`/`CNuc` clones,
point mapping (deleted 2026-11-11), copy-on-write tables, `NodeTableCache`, process-wide memos, `external/coul` (→ pycoulomb). Kept
externals: Minuit2, GSL (quadrature and special functions only), pybind11, pycoulomb.
