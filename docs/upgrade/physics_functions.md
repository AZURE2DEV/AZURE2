# The collision matrix and the observables as functions of `physics/`

Status: stage 1 (inventory) written 2026-11-15 on `feature/upgrade-1-physics-functions`, cut from `feature/upgrade`
at `486d746e`. The sections "Done" at the end say what of the list below exists; until a line stands there, it is a
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
std::complex<double> width_bilinear(const double *g, const std::complex<double> *Ag, int n);
std::complex<double> collision_element(...phases, sqrt penetrabilities, bilinear form...);   // U-matrix part of T
std::complex<double> t_matrix_element(...);                   // T of an internal pathway
struct CollisionMatrix { std::vector<std::complex<double>> A, U; int n, nc; };
CollisionMatrix collision_matrix(const LevelMatrixTables &, double E, const ChannelTablesAtE &);  // whole block
// the channel route at one energy, on the formal levels of a block:
struct FormalBlock;  struct ChannelSolve;
ChannelRoute solve_channel_space(const FormalBlock &, const std::complex<double> *lo, double E, ChannelSolve &);
void solve_bordered(...);  const std::complex<double> *f_column(...);  void formal_apply(...);  void apply_a(...);
double lu_factor(...); void lu_solve(...); hessenberg and real variants (from ChannelSpace)
```

Stage 3, `azure2/physics/Observables.h`:

```
struct PartialSum { double J; int l, lp; std::complex<double> T; };     // the coherent sum of one (J, l, l')
void add_partial(std::vector<PartialSum> &, double J, int l, int lp, std::complex<double> T);
double integrated_cross_section(const std::vector<PartialSum> &, double geom, double i1i2);   // b
double total_cross_section_term(...);                                   // optical theorem, per diagonal pathway
double phase_shift_degrees(std::complex<double> U, bool identical);
void legendre_from_products(std::vector<double> &c, double geom, double i1i2, double symmetry);
ElasticInterference elastic_interference(...);                          // direct and exchange coefficients
double upos_rk(double j2f, double L, double Lp, double Ic, int order);
capture A_y numerator and denominator per order
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

Nothing yet (stage 1 is this document).
