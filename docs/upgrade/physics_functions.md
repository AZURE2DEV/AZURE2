# The collision matrix and the observables as functions of `physics/`

Status: stage 1 (inventory) written 2026-11-15 on `feature/upgrade-1-physics-functions`, cut from `feature/upgrade`
at `486d746e`. Section 4, "Done", says what of the list below exists; until a line stands there, it is a
plan. Goal (`ARCHITECTURE.md` section 3, `FOR_THE_OWNER.md` section 10 item 5): the arithmetic of the collision
matrix and of the observables lives in `physics/` as functions of arrays and plain structs, with unit tests that
need no project file; the engine's classes collect the arrays and call. Every output file stays identical: the
arithmetic is moved, with its order of operations, not rewritten.

## 1. What the engine computes today, and where

Sizes at the start (lines): `src/AMatrixFunc.cpp` 2766, `include/AMatrixFunc.h` 415, `src/GenMatrixFunc.cpp` 756,
`src/ChannelSpace.cpp` 867, `src/LevelMatrixSolver.cpp` 112, `src/THMMatrixFunc.cpp` 863, `src/RMatrixFunc.cpp` 462,
`src/SiteObservables.cpp` 375; `physics/src/*.cpp` 2658.

At each compound energy every J group of the active nucleus takes one of two routes (`AMatrixFunc::FillMatrices`):

- the **channel route** (the default when the group has fewer coupled channels than levels): the level matrix is
  never formed; with the formal levels of the group (`ChannelSpace::Update`: a generalized eigenproblem per
  parameter set) the width bilinear forms follow from a solve in channel space (`SolveChannelSpace`);
- the **level route** (few levels, `AZURE_LEVEL_SPACE`, or a channel-route refusal near a formal pole): the
  inverse level matrix is filled (`FillLevelMatrix`) and factorized (`LevelMatrixSolver`).

Both serve the same readers: `GetUBilinear` (g_c^T A g_c'), `GetChannelCaptureBilinear`, the adjoint solves.

- `AMatrixFunc::FillLevelMatrix`: inverse level matrix of one J group at E (standard, Brune, Park; RMC term). Reads:
  `LevelState` (E, gamma, S(E_l), S', S'', Park norm), `Site::GetLoElement`, channel type and B_c,
  `Config::paramMask`. Writes: `a_inv_matrices_`.
- `AMatrixFunc::InvertMatrices`: LU of each level-route group. Reads: `a_inv_matrices_`. Writes: `solvers_`.
- `LevelMatrixSolver`: complex LU with partial pivoting; solve; bilinear form; inverse on request. Reads: a matrix.
  Writes: its factors.
- `ChannelSpace::Update`: formal levels of a group: K, N, generalized eigenproblem (extended precision), formal
  widths; or the real pencil. Reads: `LevelState`, `JGroup`, `Config::paramMask`. Writes: `FormalLevels` (kept per
  thread, keyed by level stamps).
- `AMatrixFunc::SolveChannelSpace`, `SolveFormal`, `SolvePencil`: R = Gamma^T D^-1 Gamma, deflation near a pole, LU of
  the bordered system 1 - Lambda R. Reads: `FormalLevels`, `Site::GetLoElement`. Writes: `ChannelState`.
- `SolveBordered`, `FColumn`, `FormalApply`, `ApplyA`, `ChannelBilinear`, `FormalSolve`: columns of G^T A G, A u, u^T
  A v on the channel route. Reads: `FormalLevels`, `ChannelState`. Writes: `ChannelState`.
- `GetGammaVector`, `GetChannelSolve`, `GetUBilinear`, `GetChannelCaptureBilinear`: g_c, A g_c, g_c^T A g_c' (per
  energy memo). Reads: `LevelState`, the factors. Writes: memo vectors.
- `AMatrixFunc::CalculateTMatrix`: T per pathway: phases, sqrt(P P'), 2i, the bilinear form; external capture
  amplitude times its normalisation. Reads: pathway table, `Site` tables (phases, penetrabilities, EC amplitudes),
  `LevelState`. Writes: `tmatrix_`, `ec_tmatrix_`.
- `GenMatrixFunc::CalculateBareCrossSection`: integrated cross section (sum over (J, l, l') of the squared coherent
  sums), total by the optical theorem, beta spectra, phase shift. Reads: T, pathway table, `Site` (geometry factor,
  kind), `PPair`. Writes: the site's record through `Site::Set...`.
- `BuildOrderProducts`, `BuildCoulombInterference`: sum of Z T1 T2* per Legendre order; sum of statSpin T* per (s, l).
  Reads: T, pathway table. Writes: `rt_wl_`, `rt_it_`.
- `LegendreCoefficientsOf`, `ElasticCoefficientsOf`, `AnalyzingPowerCoefficientsOf`, `GetRk`: B_L; elastic nuclear and
  interference coefficients; A_y coefficients; R_k of an unobserved primary. Reads: the two above, `Site` (geometry
  factor, secondary decay), `PPair`. Writes: returns the coefficients.
- `THMMatrixFunc::CalculateTHMCrossSection`: half-off-shell amplitude and THM cross section (channel form). Reads: T,
  the THM tables of the site. Writes: the site's record.
- `RMatrixFunc`: the channel-matrix route of the old code (R-matrix, kept as a cross-check). Reads: as `AMatrixFunc`.
  Writes: T.
- `SiteObservables.cpp` (`Site::ValueFor`, `AnalyzingPowerFor`, `ValuesWith`, `Write...`): contracts the stored
  coefficients with a row's angle (Legendre sums, Coulomb amplitude, A_y ratio); record layout. Reads: the site's
  record, the row. Writes: values.
- `AMatrixFunc::PointTBar`, `TBarAdjoint`, `LevelMatrixAdjoint`, `BilinearAdjoint`, `FormalBilinearAdjoint`,
  `LevelAdjoint`, `LevelBasis` (about 1,700 lines): reverse mode of all of the above. Reads: as above, the cotangents.
  Writes: `GradAccum`.

Already in `physics/`: `JpiBlock`, `System`, `Pathways` (the pathway table from a `System`), `ChannelFunctions`
(P, S, phases), `LevelPairShift` (the divided differences of S between two levels), `ObservableCoupling` (Z
coefficients, interference lists), `CoulombScattering` (the Coulomb amplitude), `LevelParameters`, `Widths`,
`ExternalCapture`, `SymmetricEigen`, `AngularCoupling`.

## 2. The functions to make

Indices count from 0 in `physics/`; matrices are row-major arrays unless a type says otherwise.

Stage 2, `azure2/physics/CollisionMatrix.h` (+ `LevelMatrixSolver.h`, `ChannelSolve.h`):

```
struct LevelMatrixTables {            // one J-pi block at one energy; n levels, nc channels
  int n, nc;
  const double *E, *gamma;            // E[l], gamma[l*nc+c]
  const double *shift, *shift_d1, *shift_d2;   // S_c(E_l) and its derivatives (Brune, Park)
  const double *park_norm;            // J_ll (Park)
  const std::complex<double> *lo;     // L_c(E) - B_c = S_c(E) - B_c + i P_c
  const double *boundary;             // B_c
  const char *rad_type;               // 'P' particle, 'E'/'M' photon, ...
  bool brune, park, rmc;
};
void fill_level_matrix(const LevelMatrixTables &, double E, std::vector<std::vector<std::complex<double>>> &M);
class LevelMatrixSolver;              // moved whole: Decompose, Solve, SolveInto, Bilinear, Inverse
std::complex<double> width_bilinear(const double *g, const std::complex<double> *Ag, std::size_t n);
std::complex<double> pathway_amplitude(sqrt_pen_in, sqrt_pen_out, bilinear);   // 2i sqrt(P P') g^T A g'
struct ChannelPhases { std::complex<double> coulomb, hard_sphere; };
std::complex<double> collision_element(diagonal, in, out, amplitude);          // U
std::complex<double> t_matrix_element(diagonal, in, out, amplitude);           // T of an internal pathway
rmc_t_matrix_element, external_capture_element, channel_capture_factor         // the capture elements
struct CollisionMatrix { int n, nc; bool singular; std::vector<std::complex<double>> A, U, T; };
CollisionMatrix collision_matrix(const LevelMatrixTables &, double E, const ChannelTablesAtE &);  // whole block
// the channel route at one energy, on the formal levels of a block (namespace channel_space):
struct FormalBlock;  struct ChannelSolve;  enum class Route { kTaken, kPole, kPivot };
Route solve_channel_space(const FormalBlock &, const std::complex<double> *lo, double E, ChannelSolve &);
solve_formal, solve_pencil, solve_bordered, compute_f_column, f_column, formal_apply, apply_a, channel_bilinear
lu_factor, lu_solve, hessenberg_lu, hessenberg_solve, lu_factor_real, lu_solve_real (from ChannelSpace)
```

Stage 3, `azure2/physics/Observables.h` (header only; each function adds into the caller's sum):

```
struct PartialWaveSum { double J; int l, lp; std::complex<double> T; };   // the coherent sum of one (J, l, l')
struct PartialWaveSums { rows; clear(); add(J, l, lp, T); };              // of one channel-spin group
void add_partial_wave_cross_section(sums, geom, i1i2, std::complex<double> &sum);        // fm^2
void subtract_partial_wave_cross_section(...);  double rmc_elastic_term(geom, J, i1i2, T);   // R-matrix capture
void add_beta_delayed_term(T, sum);             double cross_section_barns(sum);
void add_optical_theorem_term(J, T, double &total);   double total_cross_section(geom, i1i2, total);
std::complex<double> collision_from_t_matrix(exp_coulomb_phase, T);   double phase_shift_degrees(U, identical);
void add_interference_product(z, T1, T2, std::complex<double> &b);    // into the Legendre order's sum
void scale_legendre_coefficients(std::vector<double> &c, geom, i1i2, symmetry);          // B_L in b/sr
double upos_mixed_rk(R_L, R_LLp, R_LpLp, delta);   double upos_term(z_upos, T1, T2, l, R_l);
struct CoulombInterference { terms; clear(); add(s, l, stat_spin, T); };  // per (channel spin, l)
struct ElasticInterference { double coulomb, cross; std::vector<std::complex<double>> direct, exchange; };
int elastic_interference(terms, geom, identical, spin_dependent, spin, max_order, ElasticInterference &);
void add_capture_ay_term(ak, bk, T1, T2, double &denominator, double &numerator);
```

The engine keeps: the traversal of its pathway table, the readers of `Site`, the record writes, the memo of an
energy, the choice of the route, the profile counters.

## 3. What cannot move as it stands, and why

- **The formal levels' set-up** (`ChannelSpace::Update`, 280 lines, with `GeneralizedEigen`, `Jacobi`,
  `HessenbergTriangular`, 200 lines): arithmetic in `long double`, interleaved with the per-thread memo keyed by
  level stamps and with a debug print read from the environment. It can move as `formal_levels(arrays, flags)`
  once the gather, the key and the print are separated from the numbers; it is the last item of stage 2 and is
  reported open if the stage's time is used up before it.
- **The per-energy memo** (which channel solves exist, the inverse on request) and **the route choice** are state
  of an evaluation, not physics: they stay in `AMatrixFunc`.
- **`Site::GetLoElement` inside the level-matrix loop** becomes one gather per group and energy (fewer calls, same
  numbers).
- **The reverse mode** (stage 4): its routines read the same memo and the engine's free-parameter mask; they are
  derivatives of the functions above and move only after those are verified.
- **`SiteObservables.cpp`** is the record layout of `Site`, which the split of `Site` (another branch) is changing
  now: its contractions (Legendre sum, Coulomb amplitude times coefficients, A_y ratio) are candidates for
  `physics/`, but are not touched here.
- **`THMMatrixFunc`** and **`RMatrixFunc`**: not in this task's stages; listed for completeness.
- A function moved to another translation unit keeps its bits: the build has no `-march`, no fast-math and no
  contraction, so the result depends on the order and the types of the operations only. Speed can change where a
  call is no longer inlined; the kernels that run per level pair or per pathway are therefore header-inline.

## 4. Done

### Stage 2: the collision matrix (2026-11-15)

In `physics/` now (namespace `azure2::physics`; nothing of the engine is included):

- `LevelMatrixSolver.h`, `src/LevelMatrixSolver.cpp`: the class, moved whole from `include/`, `src/`.
- `CollisionMatrix.h`, `src/CollisionMatrix.cpp`: `LevelMatrixTables`, `fill_level_matrix` (the loop of the old
  `AMatrixFunc::FillLevelMatrix`), `width_bilinear`, `pathway_amplitude`, `ChannelPhases`, `collision_element`,
  `t_matrix_element`, `rmc_t_matrix_element`, `external_capture_element`, `channel_capture_factor`, and
  `collision_matrix` (a whole block on the level route: A, U, T for every channel pair; the call for a Python user
  or an agent, not used by the engine, which evaluates per pathway).
- `ChannelSolve.h`, `src/ChannelSolve.cpp` (namespace `channel_space`): `FormalBlock`, `ChannelSolve`, `Route`,
  `solve_channel_space`, `solve_formal`, `solve_pencil`, `solve_bordered`, `compute_f_column`, `f_column`,
  `formal_apply`, `apply_a`, `channel_bilinear` (bodies of the `AMatrixFunc` routines of those names), the six LU
  routines of `ChannelSpace`, and the six limits (`kMaxPoleRatio` ...).
- Tests, no project file: `physics/tests/collision_matrix_test.cpp` (29 checks: one level in closed form, |U| = 1,
  symmetry and unitarity for three levels and three channels, the Brune and Park elements, the Reich-Moore term,
  the decoupled level, the cutoff, the element functions) and `physics/tests/channel_solve_test.cpp` (11 checks:
  the channel route against the level route between levels, at a level energy and 1e-9 MeV from it, the refusal,
  the LU routines).

The engine after it:

- `AMatrixFunc::FillLevelMatrix` gathers the arrays of a group (`LevelMatrixGather level_gather_`) and calls
  `fill_level_matrix`; `AddAInvMatrixElement` is gone.
- `AMatrixFunc::SolveChannelSpace` gathers the Lo-matrix elements, calls `solve_channel_space` and counts the
  outcome for the profile; `SolvePencil` and `SolveFormal` are gone; `SolveBordered`, `FColumn`, `FormalApply`,
  `ApplyA`, `ChannelBilinear` are one-line calls (the reverse mode still uses those names).
- `AMatrixFunc::CalculateTMatrix` keeps the traversal of the pathways and the reads of the site, and calls the
  element functions.
- `ChannelState` is `channel_space::ChannelSolve`; `FormalLevels` is a `channel_space::FormalBlock` with the
  engine's key; `ChannelSpace::LUFactor` ... are inline calls.

Sizes (lines): `src/AMatrixFunc.cpp` 2766 -> 2380, `include/AMatrixFunc.h` 415 -> 387, `src/ChannelSpace.cpp`
867 -> 697, `src/LevelMatrixSolver.cpp` 112 -> gone; `physics/src/*.cpp` 2658 -> 3416.

Not done in stage 2: the formal levels' set-up (`ChannelSpace::Update` and its eigen-solvers) is still in the
engine (section 3, first item).

### Stage 3: the observables (2026-11-16)

In `physics/` now: `Observables.h` (header only, 195 lines; the list of section 2). Each function is the
expression of `GenMatrixFunc` it replaces, with the operations in the same order, and adds into a sum the caller
holds, so the order of the additions over pathways and groups is the engine's as before. Test, no project file:
`physics/tests/observables_test.cpp` (36 checks: a one-level, two-channel block through `collision_matrix` gives
the Breit-Wigner cross section, and its total by the optical theorem is the elastic plus the reaction one; rows
by (J, l, l'); the phase shift and its wrap; 4 pi B_0 = sigma for an s wave with `particle_interference` as the
coefficient; the direct and exchange coefficients of the Coulomb-nuclear interference, the cross term's sign, the
order limit; the capture A_y sums; the unobserved-primary, R-matrix-capture and beta-delayed terms).

The engine after it (`src/GenMatrixFunc.cpp`, `include/GenMatrixFunc.h` only):

- `CalculateBareCrossSection` keeps the walk over the exit's spin groups and pathways, the reads of the site and
  of the pair, and the writes of the record; the R-matrix-capture, total, angle-integrated and phase-shift
  branches call the functions. The `TempTMatrix` structure and its eighteen member functions (`NewTempTMatrix`,
  `IsTempTMatrixE1` ...) are gone: three `PartialWaveSums` members take their place.
- `BuildOrderProducts`, `LegendreCoefficientsOf`, `ElasticCoefficientsOf`, `BuildCoulombInterference` and the
  capture branch of `AnalyzingPowerCoefficientsOf` call the functions; `rt_it_` is a `CoulombInterference`.

Sizes (lines): `src/GenMatrixFunc.cpp` 756 -> 549, `include/GenMatrixFunc.h` 207 -> 152.

Not done in stage 3:

- `GenMatrixFunc::GetRk` (the R_k of a secondary decay) stays: it calls the engine's `AngCoeff`, which takes its
  3j and 6j symbols from GSL, while `physics/AngularCoupling.h` computes its own; the two agree to round-off, not
  to the last bit, so the move would change output digits.  (Since 2026-11-16 `AngCoeff` calls
  `physics/AngularCoupling.h` itself, `cnuc_split.md` section 9, so this reason is gone; `GetRk` is still engine
  code.)
- The analysing power of a particle exit stays in `Polarization::AmplitudeMatrix` (`PolarizationFunc.h`, an engine
  class over the compound system).
- The contractions with a reader's angle (`SiteObservables.cpp`: the Legendre sum, the Coulomb amplitude times
  the coefficients, the A_y ratio) are not touched (section 3).
- There is no single `observables(U, geometry)` call that returns every observable of a reaction: a caller
  outside the engine walks its own pathways (`Pathways.h` lists them for a `System`) and calls the functions
  above, as `observables_test.cpp` does for one block.

Checked for stages 2 and 3 alike (each on its own commit): every output file of every project in the three bases
identical to the reference tables (887, 851 and 887 files; none differs), eight fits identical to the base build
in chi-squared and output, the quick test set and the 172 regression scripts. Instructions at one thread against
the base build, after stage 3: a calculation -0.8 %, two gradient fits -1.9 % and -0.6 %, a sampler run -1.3 %.

### Stage 4: the collision matrix in pyazr (2026-11-16); the gradient not moved

Done: `pyazr.physics.collision_matrix(level_energies, reduced_widths, energy, penetrability=, shift=, boundary=,
hard_sphere_phase=, coulomb_phase=)` returns A, U and T of one block of formal levels and particle channels at one
energy as complex arrays, without a session or a project. The door is `api/src/pybind_collision_matrix.cpp` (it
checks the shapes, forms L_c - B_c and the two phase factors and calls `collision_matrix`; no physics); the Python
function is in `pyazr/physics.py`, outside the table of rows because its arguments are matrices. Test:
`tests/pyazr/collision_matrix_test.py` (ctest `pyazr_collision_matrix`), against the level matrix written in numpy,
unitarity and symmetry of U, the one-level phase shift, the Breit-Wigner form, hard-sphere scattering without a
level, a closed channel, and the six errors. The channel functions are the caller's (for a pair: the rows
`penetrability`, `shift`, `hard_sphere_phase`, `coulomb_phase`).

Checked on the commit of this stage and again on the tree it was merged into (with the grouping of a site's
values, `evaluation.md` section 7.8): the build without a warning, every output file of every project in the three
bases identical (887, 851, 887), the eight fits identical to the base build, 314 quick tests, the 172 regression
scripts. Besides the test, by hand on the merged build: a one-level, one-channel case against its closed form and
a block of five levels and three channels against a direct numpy inverse (off a level, 1e-7 MeV beside one, and
with a closed channel) agree to rounding (below 1e-15), U unitary and symmetric.

Not offered by the door at stage 4, though `fill_level_matrix` takes them: Brune's and Park's parameters, photon
channels and R-matrix capture (offered since stage 5, below). The observables' sums of `Observables.h` are not
bound: a Python caller forms a cross section from T with numpy.

Not done: **the analytic gradient.** The adjoints of the functions of stages 2 and 3 are `AMatrixFunc::PointTBar`
(cross section to T), `TBarAdjoint`, `LevelMatrixAdjoint`, `BilinearAdjoint`, `FormalBilinearAdjoint`,
`ResolveFormalCotangents` and `LevelAdjoint` (T to the level parameters), about 1,500 of the 2,380 lines of
`src/AMatrixFunc.cpp`. They read the engine's free-parameter mask, the memo of the channel route and the pathway
table in one pass; to move them, each needs the same split as its forward function (the arithmetic of one pathway
as a function, the walk in the engine), with every gradient output identical. That is a stage of its own size.

### Stage 5: the door takes the three parametrizations and photon channels (2026-11-16)

`pyazr.physics.collision_matrix` (pyazr 3.0.0.dev4) has six more keywords; the door names the members of
`LevelMatrixTables`, no formula is added and no engine file is touched:

| keyword | shape, unit | read by | member of `LevelMatrixTables` |
|---|---|---|---|
| `parametrization` | `"formal"` (default), `"brune"`, `"park"` | | `brune`, `park` (Park sets both, as the engine) |
| `boundary` | channels, dimensionless | formal only; refused otherwise | `boundary`, and B_c of `lo` |
| `shift_at_levels` | levels x channels, dimensionless | Brune, Park | `shift`: S_c(E_l) |
| `shift_derivative_at_levels_per_MeV` | levels x channels, MeV^-1 | Brune, Park; coincident levels only | `shift_d1` |
| `park_norms` | levels, dimensionless | Park | `park_norm`: J_l = 1 - sum_c gamma_lc^2 dS_c/dE at E_l |
| `photon_channels` | channels, bool | all | `rad_type` ('E' for True), `lo` = 0 |
| `use_rmc` | bool | formal, with a photon channel | `rmc` |

Choices made in the door, each stated in the docstring:

- B_c is not an argument of Brune's and Park's matrices (S_c(E_l) stands where B_c stood), so `boundary` is refused
  there instead of being read and cancelled; the door forms L_c with B_c = 0 for them, which removes the round-off
  of adding and subtracting B_c (the engine adds its B_c and subtracts it; the two agree to rounding).
- The derivative dS_c/dE enters an element only for two levels closer than `kCoincidentLevelGap_MeV` (1e-6 MeV,
  `LevelPairShift.h`): it is asked for only then, and its absence is an error naming the two levels. The second
  derivative enters no element (the gradients read it) and is not an argument.
- Park's J_l is the caller's: forming it in the door would be a formula there. `Session.park_norms()` gives it for
  a session's levels; the docstring gives the sum.
- A photon channel has L_c - B_c = 0 and unit phase factors in the engine's channel tables; the door refuses a
  shift, boundary value or phase other than 0 for it, so that a number that is not used is never accepted. Its
  `penetrability` is the caller's factor: the engine's is (E_gamma / hbar c)^(2L+1), the row
  `photon_penetrability`, or 1 with `use_rmc`. A beta-decay channel of a block is passed the same way (L_c = 0).
- With `use_rmc` the door returns T by the same formula as otherwise; the engine's R-matrix-capture element of a
  photon pathway is -i times the resonant amplitude (`rmc_t_matrix_element`), which a caller forms from U.

Tests (`tests/pyazr/`):

- `collision_matrix_test.py` (numpy, no session; 46 checks): A against the inverse of Brune's Eq. (33) and of
  Park's Eqs. (21), (22) written out in numpy (1e-15); three formal levels in two channels taken to Brune's and
  Park's parameters by numpy (the eigenvalue problem of Brune's paper, with shift functions written in the test)
  give the U of the formal levels at eight energies (largest difference 1.4e-14); two levels 0.5 eV apart against
  the divided difference; one Park level in closed form; one level with a proton and a photon channel in the
  capture form, with and without the Reich-Moore term; the 18 errors.
- `collision_matrix_session_test.py` (ctest `pyazr_collision_matrix_session`, 12 s): on T1, T2 and T3, in each
  parametrization (a session with its defaults, with `use_brune=False`, with `use_park=True`), the arrays of every
  block are read off the session (the parameter records' `start`, `channel_functions`, `park_norms`, B_c = S_c at
  the lowest level) and an observable formed in numpy from U or T is compared with `Session.evaluate`:

  | project | observable | points | largest difference (Brune / formal / Park) |
  |---|---|---|---|
  | T1 | phase shift p + 12C, 1/2+ (photon channels in the block) | 35 | 1.5e-10 deg in each |
  | T1, no external capture | sigma(p,gamma_0), sigma(p,gamma_1) | 48, 38 | 1.5e-14 / 1.3e-15 / 9.1e-15 relative |
  | T2 | sigma(p,alpha_1) / sigma(p,alpha_0) | 20 | 2.3e-15 / 1.1e-15 / 1.2e-15 |
  | T3 | phase shift alpha + alpha, 2+ (hybrid potential, a beta channel) | 22 | 1.4e-14 / 1.4e-14 / 1.8e-14 deg |
  | T3 | sigma_total n + 7Be; 7Be(n,p) | 40, 40 | 1.6e-13 in each |
  | T3 | 7Li(p,n) | 51 | 1.4e-15 / 1.4e-15 / 1.6e-15 |

  The capture cross sections are compared on a copy of T1 whose photon rows take no external capture (the multipole
  mask of the row set to 0; the test asserts that the copy's cross sections differ from T1's), because the function
  computes the R-matrix's photon channels and not the external-capture amplitudes. The 1.5e-10 deg of T1 is the
  same in the three parametrizations and belongs to the Coulomb table: in a session with
  `use_coulomb_table=False` the same comparison gives 7e-13 deg, so with the table `Session.channel_functions`
  and the evaluation do not hand out P, S and phi equal to the last digits (not looked into further; both are
  within the table's stated accuracy of 1e-7). Shown to fail: B_c taken at the
  highest level fails the eight formal comparisons; S_c(E_l) off by 1e-6 the eight of Brune and the eight of
  Park; J = 1 the eight of Park; the photon factor of another multipolarity the capture comparisons of each.

`pyazr.physics.optical_scattering` (new, `api/src/pybind_optical.cpp`) is the door to `physics/transfer`'s
distorted waves: `SolveDistortedWaves`, `ElasticOverRutherford` / `ElasticCrossSection`, `ReactionCrossSection`.
Test `optical_scattering_test.py` against FRESCO 3.5's reference files of `physics/tests/data/fresco` (p + 15N,
d + 15N, 6Li + 12C): S to 1.0e-5, sigma / sigma_R to 4.8e-5, the reaction cross section to 2.8e-7 relative.

Not done: the Weisskopf width (no function in `physics/`: it is the one formula of `pyazr/widths.py`); bound
states, DWBA cross sections and continuum bins of `physics/transfer` (each takes several potentials and a
reaction; no door); the Physics tab does not show the two functions (they are not rows of the table).
