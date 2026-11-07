# Coulomb functions and channel functions (`physics/coulomb/`, `channel_functions`)

Every particle channel of the R-matrix enters the collision matrix through four numbers at the channel energy: the
penetrability P_c, the shift function S_c, its energy derivative S′_c (Brune and Park bases) and the hard-sphere phase
φ_c, plus the Coulomb phase σ_l that relates the amplitudes of different l. All of them follow from the regular and
irregular Coulomb functions F_l, G_l and their derivatives at the channel radius, or, below threshold, from the
Whittaker function. This note fixes the interface of the part that computes them (`physics/coulomb/`, from Jakub
Skowroński's pycoulomb), the interface of `channel_functions` built on it, what in the engine they replace, the order
in which they replace it, and the tests. ARCHITECTURE.md §1 (exact energies, decision 11), §3, §4 and §6 set the
frame; this is step 3 of `cnuc_split.md`.

## 1. Coulomb functions

### 1.1 Definition

For a partition α of charge numbers Z₁, Z₂ (particles 1 and 2) and reduced mass μ (u), at channel energy E_c (MeV,
c.m., above the threshold of the pair in its states), the radial functions solve

    u″(ρ) + [1 − 2η/ρ − l(l+1)/ρ²] u(ρ) = 0,                                                   (1)

with ρ = k r (k in fm⁻¹, r in fm) and the Sommerfeld parameter η = Z₁Z₂ e² μ / (ħ² k). η > 0 is a repulsive pair,
η = 0 a neutral one (F, G are then ρ j_l(ρ) and −ρ y_l(ρ)). The normalisation is that of Abramowitz and Stegun ch. 14
and DLMF §33.2: F_l ~ sin θ_l, G_l ~ cos θ_l for ρ → ∞, θ_l = ρ − η ln 2ρ − lπ/2 + σ_l, with the Wronskian

    F′_l G_l − F_l G′_l = 1                                                                  (2)

(primes are d/dρ). The Coulomb phase is σ_l(η) = arg Γ(l + 1 + iη), the continuous branch, σ_0(0) = 0 and
σ_l = σ_{l−1} + arctan(η/l). The engine uses only e^{i(σ_l − σ_0)} (`conventions.rst`, Eq. conv-phases).

### 1.2 Interface

Only the dimensionless arguments enter; the kinematics (η and ρ from E_c, a_c and the masses) belong to the caller
(section 3.2). Namespace `azure2::physics::coulomb`, header `azure2/physics/coulomb/Coulomb.h`:

    enum class Method { coulsm, coul, coulfg, coulcc };
    enum class Status { ok, invalid_argument, failed };
    struct CoulombWaves { double F, Fp, G, Gp, sigma; Status status; };
    CoulombWaves coulomb(int l, double eta, double rho, Method method = Method::coulsm) noexcept;
    double coulomb_phase(int l, double eta) noexcept;
    const char *method_name(Method); bool method_from_name(const std::string &, Method &);

The rules are pycoulomb's (its `DESIGN.md` §2–3): a pure function of its arguments, no object, no static or global
mutable state, no cache, no printing, no exception leaving it; the same arguments give the same bits in any thread
and in any order of calls. Arguments are checked first (0 ≤ l ≤ 10⁴, η and ρ finite, ρ > 0); an algorithm that
throws, reports failure or returns a non-finite value gives `failed`. On any status but `ok` every number is NaN, so
a caller that ignores the status cannot use a wrong value without seeing NaN.

### 1.3 Algorithms, one file each behind the same signature

| `Method` | source (pycoulomb) | what it is |
|---|---|---|
| `coulsm` (default) | `libcoulomb/coulsm/` | N. Michel's `cwfcomp` (Comput. Phys. Commun. 176 (2007) 232) with the patch from the SAMMY side; AZURE2's `external/coul` since 2025-09-05. The default because it reproduces the engine's numbers bit for bit. |
| `coul` | `libcoulomb/coul/` | `cwfcomp` as in AZURE2 before 2025-09. |
| `coulfg` | `libcoulomb/coulfg/` | A. R. Barnett's COULFG (Comput. Phys. Commun. 27 (1982) 147), C translation; F′, G′ by the recurrence A&S 14.2.1 from l + 1. |
| `coulcc` | branch `coulcc-cpp`, `coulcc/` | COULCC (Thompson and Barnett, Comput. Phys. Commun. 36 (1985) 363), the C++17 translation bit-identical to the Fortran; called with real ρ, η and λ = l. |

Each algorithm is one adaptor, `Status <method>(int l, double eta, double rho, CoulombWaves &out)`, in
`physics/coulomb/src/method_<name>.cpp`; `Coulomb.cpp` checks the arguments, dispatches, sets σ_l and applies the
NaN rule. The third-party sources sit unchanged under `physics/coulomb/<method>/` apart from a namespace and the
`M_*` guards for MinGW/MSVC; `README.md` there records the pycoulomb commits and each local change. COULCC's complex
interface (complex l, η, ρ; Bessel functions) is not exposed: AZURE2 has no use for it yet.

GSL's `gsl_sf_coulomb_wave_FG_e` is not a method of `physics/coulomb/` (physics/ links nothing). It stays as an
independent cross-check in a test outside physics/ (section 5), and the engine's `--gsl-coul` path stays in the
engine (section 4, open question).

### 1.4 Whittaker functions (bound channels)

For E_c < 0 the decaying solution of (1) with κ = |k| and η_b = Z₁Z₂ e² μ / (ħ² κ) is the Whittaker function

    W_l(ρ) = W_{−η_b, l+½}(2ρ) = e^{−ρ} (2ρ)^{l+1} U(l + 1 + η_b, 2l + 2, 2ρ),   ρ = κ r,           (3)

and the bound-channel shift function is S_c = ρ W′/W at ρ = κ a_c. The interface is pycoulomb's: `whittaker(l, η_b,
ρ)` → W, W′, status, and `whittaker_scaled(l, η_b, ρ)` → mantissa, log₁₀ scale, ρW′/W, which stays finite just below a
charged threshold, where η_b grows like E_b^−½ and W underflows. pycoulomb computes U with GSL's `gsl_sf_hyperg_U`,
which physics/ cannot use. **This branch does not move the Whittaker functions**: `WhitFunc` and `ShftFunc` stay in
the engine until a GSL-free U is in `physics/coulomb/` (proposed: COULCC at imaginary argument, which gives H⁺ and
hence W up to a known factor, checked against mpmath and GSL), with a comparison of the bound-state outputs, which will
move in their last digits (section 6).

## 2. Channel functions

### 2.1 Definition

At channel energy E_c > 0 and channel radius a_c (fm), ρ = k a_c,

    P_c = ρ / (F_l² + G_l²),   S_c = ρ (F_l F′_l + G_l G′_l) / (F_l² + G_l²),
    e^{−iφ_c} = (G_l − i F_l) / (F_l² + G_l²)^½,   ω_c = σ_l − σ_0,                                (4)

all dimensionless, and S′_c = dS_c/dE_c (MeV⁻¹). The formulas and their order of operations are the engine's
(`CoulFunc.cpp`, `conventions.rst` Eqs. conv-ps, conv-phases), including the two guards for η ≳ 110 just above a
charged threshold, where F² + G² overflows: S and e^{−iφ} are then formed with F, G scaled by hypot(F, G), and only
then, so that every value the plain formula gives is unchanged. S′_c is the Richardson-extrapolated central
difference of S with step h = 10⁻³ E_c, (4 D(h/2) − D(h))/3, D(h) = [S(E + h) − S(E − h)]/2h; its error is about
10⁻⁸ relative, set by the 10⁻¹⁰ accuracy of S from `coulsm` (`cnuc_split.md` §5).

### 2.2 Interface

    struct ChannelKinematics { int z1z2; double mu_u, m1_MeV, m2_MeV; bool relativistic; };
    double eta(const ChannelKinematics &, double E_c);       double rho(…, double a_fm, double E_c);
    struct OpenChannelValues { double P, S, dS_dE; std::complex<double> exp_minus_i_phi; double sigma; Status status; };
    OpenChannelValues channel_functions(const ChannelKinematics &, int l, double a_fm, double E_c,
                                        coulomb::Method = coulsm, ChannelQuantities = all);

`eta` and `rho` are the engine's `Kinematics::Channel::Eta/Rho` (non-relativistic, or the relativistic forms of
Haeck et al., LA-UR-19-31084), moved into physics/ with the constants they use, so that the engine and physics/ share
one formula; the engine's `Kinematics` calls them. The same function is exposed in the binding with array inputs
(ARCHITECTURE.md §8) through the existing `coulomb_functions` call.

`channel_functions` covers the open channel, E_c > 0 with η ≤ 100. Below threshold, at threshold and in the
"Coulomb threshold" band η > 100 (P = 0 in double precision, S continued through E_c = 0), the engine's
`ChannelFunc` keeps its rules on the Whittaker branch until section 1.4 is done; it then becomes one function of
physics/ for all E_c.

### 2.3 The table

ARCHITECTURE.md §3–4 call for a read-only map of P, S, S′, φ, σ_l per (partition α, l, channel energy), built once per
evaluation-energy set and shared by every reaction, angle and J^π block (today's `EDependentTables::Slot(jGroup,
channel)` stores one (α, l) once per J^π group and channel spin). Its contract:

- keyed by (α, l, a_c) and the index of the evaluation energy E_n; values are those of `channel_functions` at
  E_c = (E_n − S_α) − e_α (MeV, c.m.; the order of `cnuc_split.md` §5);
- built in parallel, each thread writing only its own rows, so the values do not depend on the thread count;
- read-only once built; when a free energy shift ν moves the energies of a measurement, a new table for those
  energies is built single-threaded between evaluations and replaces the old one (decision 11); nothing writes into a
  table that an evaluation reads;
- no other cache: neither a process-wide memo nor a tolerance on the energy.

The table replaces the Coulomb rows of `EPoint::EDependentTables` (penetrability, shift, hard-sphere and Coulomb
phase), their copy-on-write sharing and `NodeTableCache` for those rows. That replacement changes who owns the
per-point energies (Phase 4, `EnergyGrid`), and is not part of this branch (section 6).

## 3. What this replaces

| today | role | becomes |
|---|---|---|
| `external/coul` (`cwfcomp`, `cwfcomp_accurate`) | the Coulomb routine | `physics/coulomb/coulsm/` (and `coul/`); deleted |
| direct `Coulomb_wave_functions_accurate` calls in `CoulFunc` (three: plain, hybrid boundary values, hybrid fallbacks), `CoulombTable::ExactOpen`, `ThmDistortion` | five copies of one call | `coulomb::coulomb(l, η, ρ, coulsm)` |
| `CoulFuncCache` (process-wide, per-thread memo; energies matched to 10⁻¹² MeV) | speed | deleted, no replacement in this branch; the table (§2.3) later |
| `CoulFunc::Penetrability`, `PEShift`, `HardSpherePhase`, `PEShift_dE` | Eq. (4) | thin calls of `channel_functions` (the Coulomb-table and hybrid paths stay in `CoulFunc`) |
| `Kinematics::Channel::Eta`, `Rho` | η, ρ | call physics/ |
| `WhitFunc`, `ShftFunc` | bound channels | stay (section 1.4) |
| `EDependentTables` Coulomb rows | per-point store | stay (section 2.3) |
| `g_shiftDerivMemo` (exact-argument memo of S′, S″ at level energies) | speed | stays; becomes the per-stage table of level energies with step 3b |

The memo's tolerance makes the engine's last digits depend on the order of the calls: an energy within 10⁻¹² MeV of
one computed earlier in the same thread gets the earlier one's functions. Without it every point gets its own values,
so outputs may move at the level allowed by S (10⁻¹⁰) and S′ (10⁻⁸) for an energy change of 10⁻¹² MeV, and the
pyazr/CLI identity checks of THM (HANDOFF, director's entry of 2026-11-07), which suspect this call-order dependence,
are expected to pass.

## 4. Engine and front ends

The engine always calls `Method::coulsm` in this branch: no new option. `--gsl-coul` keeps its present meaning (GSL's
routine, called from `CoulFunc`, which may use GSL), because removing or remapping a CLI flag is a user-visible change
to be decided (section 6). pyazr's `cache_stats()` keeps its keys; the six Coulomb-memo counters are 0, since there is
no memo, and its docstring says so.

## 5. Tests

- `physics_coulomb_unit` (physics/ only): every method against mpmath at 40 digits on pycoulomb's reference set
  (`tests/data/coulomb.csv`, copied with its generator), with the tolerance each method meets there stated in the
  test; the Wronskian (2) to 10⁻¹⁰; σ_l against the arctan recurrence; the four methods against each other where
  each is accurate (identity test); argument checks (NaN and status); bit-identical results from 8 threads calling
  in different orders. It fails if a method is wired to the wrong source, if a sign or normalisation is changed, or
  if any state is shared between calls.
- `physics_channel_functions_unit`: P, S, φ against Eq. (4) evaluated with mpmath's F, G; the neutral closed forms
  (l = 0: P = ρ, S = 0, φ = ρ; l = 1: P = ρ³/(1 + ρ²), S = −1/(1 + ρ²), φ = ρ − arctan ρ); S′ against the derivative
  of the closed form; η and ρ against the non-relativistic formulas and the relativistic closed forms of
  `Kinematics.h`.
- `coulomb_gsl_crosscheck` (outside physics/, links GSL): `coulsm` against `gsl_sf_coulomb_wave_FG_e` where GSL is
  accurate.
- Engine: every `.azr` under `tests/` (179) in the Brune, formal (`--no-brune`) and Park (`--use-park`) bases, before
  and after, every output file compared number by number. The routing commit must be bit-identical; the memo's
  removal may move numbers at ≤ 10⁻⁹ relative on observables, each changed file classified; anything else stops the
  step. `tests/upgrade/check.sh`, `ctest -LE slow`, `tests/run_tests.sh`, 1 vs 8 threads identical, T2 benchmark.

Code: `physics/coulomb/` (sources, adaptors, `Coulomb.h`), `physics/src/ChannelFunctions.cpp`,
`physics/include/azure2/physics/ChannelFunctions.h`, `physics/tests/coulomb_*`, `src/CoulFunc.cpp`, `include/Kinematics.h`.

## 6. Order and what is left

1. This note. 2. `physics/coulomb/` with its tests. 3. The engine calls it; `external/coul` deleted (bit-identical).
4. `channel_functions` and the kinematics in physics/; `CoulFunc` calls them (bit-identical). 5. `CoulFuncCache`
deleted (outputs classified). Left for step 3b (`cnuc_split.md` §3): the Whittaker functions without GSL and the
deletion of `WhitFunc`, `ShftFunc` and `ChannelFunc` (bound branch into `channel_functions`); the table of §2.3 in
place of the `EDependentTables` rows and `g_shiftDerivMemo`; the hybrid-potential channels ((α, l, j) entries of the
table); `CoulombTable` (an interpolation of these same functions, opt-in) and `--gsl-coul`, which are user-visible
options to keep, remap to a `Method`, or retire.
