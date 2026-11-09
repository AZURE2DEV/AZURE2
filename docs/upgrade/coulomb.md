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
which physics/ cannot use.

**Design (2026-11-08, step 3b (i)).** Three ways were compared on 990 points (l = 0–8; η_b = 0, 0.01, 0.3, 1, 2, 3,
5, 10, 30, 100, 300; ρ = 10⁻³–50) against mpmath's `whitw` at 40 digits:

| method | worst relative error of W | points without a value |
|---|---|---|
| GSL `gsl_sf_hyperg_U` (the engine's `WhitFunc`) | 1.2 × 10⁻⁹ (l = 0, η_b = 100, ρ = 10⁻³) | 90 (GSL error, overflow or underflow) |
| COULCC at ρ = iκr, η = −iη_b (the earlier proposal) | — | 0/0 at η_b = l + 1, l + 2, … (below) |
| the integral representation below (Python prototype) | 4.6 × 10⁻¹³ (l = 5, η_b = 300, ρ = 10⁻³) | 0 |

COULCC was dropped: H⁺_l(−iη_b, iκr) = e^{−i(l+η_b)π/2} [Γ(l+1+η_b)/Γ(l+1−η_b)]^½ W_l(κr) (DLMF 33.2.7; checked
with mpmath to 30 digits), and the factor vanishes at η_b = l + 1, l + 2, …, so W = H⁺/factor loses digits there,
however accurate H⁺ is. Its logarithmic derivative is normalization-free, but the ANC conversions and N_f need W
itself.

**Method.** For a = l + 1 + η_b ≥ 1 (η_b ≥ 0: a neutral or repulsive pair; nuclear pairs are never attractive),
b = 2l + 2, z = 2ρ:

    U(a, b, z) = (1/Γ(a)) ∫₀^∞ e^{−zt} t^{a−1} (1+t)^{b−a−1} dt,    dU/dz = −a U(a+1, b+1, z)   (DLMF 13.4.4, 13.3.22)

so that ln W = −ρ + (l+1) ln 2ρ − ln Γ(a) + ln I(a, b, 2ρ) and ρW′/W = −ρ + (l+1) − 2aρ I(a+1, b+1, 2ρ)/I(a, b, 2ρ)
(primes d/dρ). The integrand is positive and has one maximum t₀ (a quadratic in t); with t = t₀ exp[(π/2) sinh s]
(double-exponential in both tails) the trapezoidal sum in s converges geometrically. The sum is formed in
logarithms (largest term factored out), so W is never formed before the end: `whittaker_scaled` returns mantissa
and power of ten, and `whittaker` returns W, W′ when they are inside the double range and `out_of_range` (NaN)
otherwise. The step is halved from 1/2 until two sums agree to 4ε; terms are summed while they exceed e⁻⁴⁰ of the
largest. No GSL, no table, no state: the same arguments give the same bits in any thread.

**Interface** (`Whittaker.h` in `physics/coulomb`, beside `Coulomb.h`):

    struct WhittakerWave { double W, Wp; Status status; };                 // W_l(ρ) = W_{−η_b, l+½}(2ρ), dW/dρ
    struct WhittakerScaled { double mantissa; long log10_scale; double rho_dlog; Status status; };
    WhittakerWave whittaker(int l, double eta_b, double rho) noexcept;    // out_of_range outside the double range
    WhittakerScaled whittaker_scaled(int l, double eta_b, double rho) noexcept;   // W = mantissa·10^scale, ρW′/W

`invalid_argument` for l < 0 or > 10000, η_b < 0 or not finite, ρ ≤ 0 or not finite. The binding exposes both with
array inputs.

**Tests.** `physics_whittaker` against a committed mpmath reference (`physics/tests/data/whittaker.csv`, made by
`make_whittaker_reference.py`, 40 digits: l = 0–10, η_b = 0–1000, ρ = 10⁻⁴–200, the scaled form where W leaves the
double range), with the tolerance set from the measured error; ρW′/W against mpmath's derivative; the same bits from
eight threads; `invalid_argument` cases. Outside physics/, a reference test compares with GSL (where GSL gives a
value) and times both. Speed decides the engine part: `ECIntegral` spends about 85 % of its time in GSL's U.

**Part 1 verified (2026-11-09; build of bf8622e2, GCC 11.5 Release, pyazr on).** No warnings. `ctest -LE slow`
114/114, with `physics_whittaker`, `whittaker_gsl` and `pyazr_whittaker_wave` (the binding) passing in the full
build; 1 vs 8 threads byte-identical (T1–T3, 13N, three bases); the 179 projects × 3 bases against a18fe89e:
2040/2040 files identical (the engine does not call the new functions yet). `run_tests.sh` was not run: no engine
path changed.

**Engine (step 3b (i), second part, after the timing).** `WhitFunc` becomes a thin adapter over `whittaker` /
`whittaker_scaled` (it keeps η_b and ρ from `Kinematics`); `ShftFunc`'s bound branch takes ρW′/W in place of
`CoulombTable::ExactBound` (GSL's U ratio, trusted for η_b ≤ 1000 and η_b ≤ 50 or z ≥ 0.05, linearly interpolated
to S(0) closer to threshold; the numerical derivative of W is only the opt-in `AZURE_BOUND_SHIFT_NUMERICAL` route).
The new function holds 10⁻¹² over that band and below it, so whether the interpolation band narrows (outputs near
threshold move) is decided with its own comparison. `CoulombTable`, `ThmNumerics` and `ThmFunc` call the same
function; GSL's U leaves the engine. In the data path `ECIntegral` evaluates W at every radius its adaptive rule
samples (13× GSL's time per call); the rate path already tabulates the final-state W once per pathway on a fixed
grid, and the per-integral table of step 3b (ii) does the same for the data path. Bound-state outputs (ANC ↔ γ conversions, N_f, external capture to bound states, sub-threshold shifts) move
in their last digits; the 179-project comparison classifies every change, and the T1–T3 truth checks must hold.

**Engine part as built (2026-11-09, branch `feature/upgrade-1-whittaker-engine`).** Every call of GSL's U in the
engine goes through the two physics functions, with the engine's η_b and ρ = κr unchanged (non-relativistic or
`Kinematics` relativistic forms):

| engine call | was | becomes |
|---|---|---|
| `WhitFunc::operator()` | e^{−z/2} z^{l+1} `gsl_sf_hyperg_U` (z = 2κr) | `whittaker(l, η_b, κr).W`; outside the double range mantissa × 10^scale of `whittaker_scaled`, i.e. 0 on underflow and inf on overflow (IEEE), never a GSL error |
| `WhitFunc::Scaled` | `gsl_sf_hyperg_U_e10_e` | `whittaker_scaled` (mantissa, power of ten) |
| `CoulombTable::ExactBound` (ρW′/W; `ShftFunc`'s bound branch and the table's bound sheet) | ratio of two `gsl_sf_hyperg_U_e10_e` | `whittaker_scaled(...).rho_dlog` |
| `ThmBoundStateTail` | 2κ e^{−κr} `gsl_sf_hyperg_U_e`(1 + η_b, 2, 2κr) | `whittaker(0, η_b, κr).W / r` (same function: W_0/r) |

Invalid arguments (NaN η_b or ρ, ρ ≤ 0) keep the engine's present answer, 0 (the `WhitFunc` FIXME and `Scaled`'s
guard); they are not new failure modes. `ShftFunc`'s interpolation band near threshold (η_b ≤ 1000, and η_b ≤ 50 or
z ≥ 0.05) is kept unchanged here, so that this step changes only the function; narrowing it is a separate,
measured step. `ZeroEnergyLimit` (Bessel K) and the opt-in numerical route stay on GSL: they are not U.

The Coulomb table's bound sheet is filled from `ExactBound`, so its contents change: the table format becomes 2
(`kFormat`), and a format-1 file is refused with the existing warning ("table format 1, this build reads 2") and the
exact functions are used until a table is built. No file is overwritten.

**Comparison for this step.** Reference: the head of `feature/upgrade-1-whittaker` (engine identical to step 4).
Both builds run the 179 projects × 3 bases without a Coulomb table (`AZURE2_COULOMB_TABLE` pointing to a missing
file), so the difference is the function alone; the new build then also with a table of format 2 against itself
without it. Every changed file is classified; bound-state quantities (ANC ↔ γ, N_f, external capture to bound
states, sub-threshold S) are expected to move at ≲ 10⁻⁸ relative (GSL's error, measured above), anything larger is
explained. Timing: `ECIntegral`'s data path on a capture project (T1 or 13N), calc, one thread, before and after.

**Engine part verified (2026-11-09; build of 5ed85322, GCC 11.5 Release, pyazr on).** No warnings; no engine file
includes `gsl_sf_hyperg.h`. Without a table, against bf8622e2: Brune 595/689, formal 559/662, Park 605/689 files
identical; every changed observable within 2.2 × 10⁻⁹ of its column's largest value (worst: `relativistic/he5`);
`param.par`/`parameters.out` at most 9 × 10⁻⁸ relative in printed digits, larger relative changes only for values at
rounding (formal zero-input widths 10⁻¹⁸–10⁻³⁷; cross sections of 10⁻²² at 0° and 180° where they vanish), one
photon width of 0.99999999 eV printed as 1000 meV (`ec_threshold/ng`), and signs of −0.000000; all 20 `intEC.dat`
files change (every token, ~10⁻¹⁰), hence EC signature v6 with the stale-file case (s) of `tests/ec_m1_exact`
(fails without the bump). With a table built by each build (format 1 by bf8622e2, format 2 by this one), the
table-vs-exact differences are of the same size for both (worst `c12c12`, 5.7 and 7.6 × 10⁻⁵, interpolation), and
`c12c12` is identical between the two builds. `ctest -LE slow` 114/114 with a format-2 table (`pyazr_park_dsde`
compares a table session with an exact one and fails when the only table on the machine is refused);
`ec_m1_exact` PASS; `run_tests.sh` 133/133 (1421 s, 11 threads, format-2 table); 1 vs 8 threads byte-identical
(T1–T3, 13N, three bases). **Cost** (one thread, calc, mean of three, own table): 13N 7.10 → 9.14 s, T1 2.62 →
4.73 s; `AZURE_PROFILE` puts it all in set-up (the external-capture table: 1176 integrals, 1320 node values in T1,
recomputed every run, even with `intEC.dat` present), evaluation ~1 ms in both. This is the 13× cost per W call in
`ECIntegral`'s data path; the per-integral W table (step 3b (ii)) removes it and should come before merging.

**Step 3b (ii): the final-state W of the data path, memoized (design, 2026-11-09).** The cost above is the
final-state factor W_{l_f}(κ_f r) of `ECIntegral`: `fgEvaluate` (FW, GW) and the outgoing factor of `WWIntegrand`
evaluate it at every radius GSL's `qagiu` samples, for every node energy of every `ECTable` term, although it does
not depend on the entrance energy. The rate path's fixed Laguerre grid would remove it, but it changes the quadrature
of the data path (outputs move). Instead: `qagiu` maps [a, ∞) to t ∈ (0, 1] and bisects; a subinterval it creates is
the same bits whatever the integrand, and so are its 21 Gauss–Kronrod radii. The integrals of one final state at
different entrance energies therefore sample largely the same radii, and a memo keyed on the exact arguments returns
the bits a new evaluation would.

*Rule.* `WhitFunc::ScaledMemo(l, r, E_b, log10Scale)`: η_b and ρ formed as in `Scaled`, then a thread-local table
keyed on (l, η_b, ρ) — the arguments of `whittaker_scaled` themselves, so no key can be stale and no model or pair
identity enters — returns the stored (mantissa, scale) or computes and stores them. `whittaker_scaled` is a pure
function (same arguments, same bits, any thread), so results are bit-identical to `Scaled`. The table holds at most
2¹⁷ entries per thread (about 10 MB) and is cleared when full. It is used for the final-state factor in
`fgEvaluate` and in `WWIntegrand` (outgoing factor) only; the entrance factor of `WWIntegrand` depends on the
entrance energy and stays `Scaled`. `AZURE_PROFILE` reports hits and misses.

*Expected.* 179 × 3 byte-identical (with and without a table, against the head of this branch); threads identical;
the cost measured again (13N, T1, one thread, calc). If the hit rate is low (the trees of different energies share
few subintervals), the fallback is the fixed grid, as its own design with a comparison.

**Step 3b (ii) verified (2026-11-09; build of 82ca0985, GCC 11.5 Release, pyazr on).** No warnings. 179 × 3 against
5ed85322 (this branch before the memo): with a format-2 table 2040/2040 identical; without a table 1977/2040, the
63 others being `intEC.dat.sig` alone (the reference run without a table was made by 48fb454f, before the v6 bump;
the memo build's signatures equal the v6 build's, with and without a table). 1 vs 8 threads byte-identical (T1–T3,
13N, three bases); `ctest -LE slow` 114/114. Memo on T1 (calc): 394980 hits, 870 computed (99.8 %); 13N 394965 and
495. **Cost** (one thread, calc, mean of three, each build with its own table): T1 2.89 s (`-whittaker` head, GSL)
→ 5.14 s (engine switch) → 2.77 s (memo); 13N 7.60 → 9.71 → 7.08 s. The engine switch no longer costs time.

**The interpolation band near threshold (design, 2026-11-09).** `ShftFunc::operator()` evaluates ρW′/W exactly for
binding B ≥ B_min. Closer to threshold it interpolates linearly between S(−B_min) and S(0) (`ZeroEnergyLimit`, Bessel K).
B_min is the binding where η_b = 1000, or, if larger, the smaller of the bindings where η_b = 50 and where
z = 2κa = 0.05. The second condition described where GSL's U went wrong (off by 10⁻⁶ at η = 111, z = 1.2 × 10⁻³ for
p+p); `whittaker_scaled` has no such weakness. Measured on p+¹²C, p+p, α+¹²C and ¹²C+¹²C, l = 0 and 3, B = 10⁻¹⁴–1 MeV:
- status ok everywhere;
- S tends to S(0) linearly in B;
- below about 10⁻⁹ MeV (η_b ≳ 10⁴) it jitters at 10⁻⁹ relative, beyond the η_b ≤ 1000 its test covers to 10⁻¹²;
- each call takes 2–9 µs, but 16–230 µs at η_b ≳ 10⁷.

So the η_b = 1000 limit stays, and the GSL condition is dropped: B_min = B(η_b = 1000). For p+¹²C, B_min goes from
3.3 × 10⁻⁴ to 8.2 × 10⁻⁷ MeV. In what is left of the band, S − S(0) is linear to about 10⁻¹¹ (its curvature,
measured from the values above, is about 3 × 10⁻⁵ of the linear term at 10⁻³ MeV and shrinks with B).

*Expected.* S changes, without a table, only for bindings between the old and the new B_min. There the linear
interpolation was off by up to about 10⁻⁸ absolute, so outputs near threshold move at that level. dS/dE moves by up to
about 10⁻⁵ relative on a level within 3 × 10⁻⁴ MeV below threshold. The Coulomb table is filled from `ExactBound` (not
`operator()`), so its contents and format do not change. With a table, only the bindings the table does not cover
change. 179 × 3, with and without a table, against the head of this branch: every changed file classified. A change
beyond about 10⁻⁷ relative in an observable is examined before the step goes on. Then ctest, `run_tests.sh` and the
thread check.

**Band narrowed, verified (2026-11-09; build of 37e1c5ba).** 179 × 3 against the merged head db972c84:
- without a table, 2040/2040 identical: no project has a bound channel between the old and the new B_min;
- with a format-2 table, 2039/2040: `ec_threshold` `ng`, R = 3, one number, 1.2 × 10⁻¹¹ of its column (a point the
  table does not cover), class (b).

ctest 117/117, `run_tests.sh` 133/133, 1 vs 8 threads identical.

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
| `CoulFuncCache` (process-wide, per-thread memo; energies matched to 10⁻¹² MeV) | speed | deleted in step 3b, replaced by the channel table and the radial scratch (`channel_table.md`) |
| `CoulFunc::Penetrability`, `PEShift`, `HardSpherePhase`, `PEShift_dE` | Eq. (4) | thin calls of `channel_functions` (the Coulomb-table and hybrid paths stay in `CoulFunc`) |
| `Kinematics::Channel::Eta`, `Rho` | η, ρ | call physics/ |
| `WhitFunc`, `ShftFunc` | bound channels | stay (section 1.4) |
| `EDependentTables` Coulomb rows | per-point store | filled from the channel table (step 3b); the rows themselves stay until `EnergyGrid` (section 2.3) |
| `g_shiftDerivMemo` (exact-argument memo of S′, S″ at level energies) | speed | stays: level energies are fit parameters, not evaluation energies (`channel_table.md` §2) |

The memo's tolerance makes the engine's last digits depend on the order of the calls: an energy within 10⁻¹² MeV of
one computed earlier in the same thread gets the earlier one's functions. Deleting it lets every point get its own
values; outputs then move at the level allowed by S (10⁻¹⁰) and S′ (10⁻⁸) for an energy change of 10⁻¹² MeV, and the
pyazr/CLI identity checks of THM (HANDOFF, director's entry of 2026-11-07) pass. It is also most of the engine's speed
where energies recur, so it is deleted only with the table that replaces it (§7).

## 4. Engine and front ends

The engine always calls `Method::coulsm` in this branch: no new option. `--gsl-coul` keeps its present meaning (GSL's
routine, called from `CoulFunc`, which may use GSL), because removing or remapping a CLI flag is a user-visible change
to be decided (section 6). The memo stays in front of `CoulombWavesAt` as before, so `cache_stats()` is unchanged.

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
  and after, every output file compared number by number. The routing commits must be bit-identical (each changed file
  classified otherwise; anything above 10⁻⁹ relative on observables, unexplained, stops the step). `tests/upgrade/check.sh`, `ctest -LE slow`, `tests/run_tests.sh`, 1 vs 8 threads identical, T2 benchmark.

Code: `physics/coulomb/` (sources, adaptors, `Coulomb.h`), `physics/src/ChannelFunctions.cpp`,
`physics/include/azure2/physics/ChannelFunctions.h`, `physics/tests/coulomb_*`, `src/CoulFunc.cpp`, `include/Kinematics.h`.

## 6. Order and what is left

1. This note. 2. `physics/coulomb/` with its tests. 3. The engine calls it; `external/coul` deleted (bit-identical).
4. `channel_functions` and the kinematics in physics/; `CoulFunc` calls them (bit-identical). 5. `CoulFuncCache`
deleted (outputs classified) — moved to step 3b, after the table (§7). Left for step 3b (`cnuc_split.md` §3): the Whittaker functions without GSL and the
deletion of `WhitFunc`, `ShftFunc` and `ChannelFunc` (bound branch into `channel_functions`); the table of §2.3 in
place of the `EDependentTables` rows and `g_shiftDerivMemo`; the hybrid-potential channels ((α, l, j) entries of the
table); `CoulombTable` (an interpolation of these same functions, opt-in) and `--gsl-coul`, which are user-visible
options to keep, remap to a `Method`, or retire.

## 7. As built (branch `feature/upgrade-1-coulomb`, 2026-11-07)

Steps 1–4 of section 6 are done, one commit each. Step 5 (the memo's deletion) was built and measured, and moved to
the branch `feature/upgrade-1-coulomb-nomemo`, to be re-cut on top of the per-stage table (step 3b): without the
memo the engine is 3–17× slower and `convolution_unify` exceeds its time limit (below).

**Outputs.** Engine before (f549fa23) and after each step, GCC 11.5, Release, one thread, every `.azr` under
`tests/` (179) in three bases, every output file compared number by number (`cnuc_split.md` §6, step 3):

- Routing to `physics/coulomb` (step 3): identical except the projects with a hybrid-potential channel (`T3`,
  `T3_start`, `hybrid_potential`): observables ≤ 9.9 × 10⁻¹¹, one width in `parameters.out` 5.0 × 10⁻⁸ (formal
  basis, the formal → Brune iteration). Cause: `cwfcomp` returns F and G as complex numbers whose imaginary parts
  are not exactly 0 for real arguments (7 of 125 probe points, up to 1.0 × 10⁻¹⁵ relative); the hybrid model formed
  H⁺ = G + iF from the complex values and so carried them into the Numerov boundary values, `coulomb()` returns the
  real parts. Plain Coulomb channels never used the imaginary parts.
- `channel_functions` and the kinematics in physics/ (step 4): bit-identical to step 3 (2040 files).

**Speed of this branch** against f549fa23, one Condor host, both engines alternately, `--no-coulomb-table`, mean of
two: host d32cepyc086, mean of two, before (f549fa23) → this branch: T2 calc 1.56 → 1.36 s (1 thread; one
f549fa23 run at 1.89 s), 0.44 → 0.42 s (8); T2 fit200 80.5 → 81.9 s (1 thread), 36.5 → 37.0 s (8); menu-1
calculation, one thread, `13N` 10.05 → 10.16 s, `natural_neon` 12.35 → 12.45 s; `convolution_unify` `tgt_fixed`
(2 threads) 29.9 → 30.6 s. The routing costs at most 2 % (fit and `tgt_fixed`), within the run-to-run spread of the
calculations; peak RSS within 2 %.

**The memo's deletion, measured (branch `-nomemo`, c3dc22e5).** Observables ≤ 7.0 × 10⁻¹¹ relative, except the
rounding residue of an exact zero (E1 capture at θ = 0°, 1.5 × 10⁻²⁵ → 3.0 × 10⁻²⁵ beside 1.2 × 10⁻⁹ at 60°);
`intEC.dat` of 20 projects ≤ 4.0 × 10⁻¹¹; the THM checks `pyazr_thm_profiled` and `pyazr_thm_fixed_angle` pass (the
memo's call-order dependence was their failure). Speed (host qa-a10-001, same method): T2 calc 1.37 → 6.16 s
(1 thread), 0.46 → 1.40 s (8); fit200 93.1 → 114.4 s, 42.7 → 48.1 s; menu-1 calculation, one thread, `13N` 10.0 →
32.1 s, `beam_profile_kernel` 52.6 → 114.9 s, `natural_neon` 12.3 → 105.6 s, `18O_p_a_thm` 0.14 → 1.39 s;
`convolution_unify` `tgt_fixed` (2 threads) 84 → 1414 s, past its 600 s limit. The memo returned stored functions
wherever an energy recurred: shared point energies, target and resolution sub-points, external-capture radii, a
second pyazr session.

**Next (step 3b).** The table of section 2.3, one evaluation per (α, l, a_c, E_c) of the run's evaluation energies,
shared read-only, and a per-integral table of F, G at the quadrature radii for the external-capture integrals; then
the `-nomemo` commits re-cut on it, with the timings above repeated (they must come back near those of this branch).
