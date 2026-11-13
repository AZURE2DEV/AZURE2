# The radiative decay of a level: total radiative width, branchings and mixing ratios (`gamma_decay`)

A level λ that decays by several γ-ray transitions is known from experiment as a total radiative width Γ_γ (eV) (or a
resonance strength), the branchings b_f to its final states f, and, for a transition with two multipolarities, the
mixing ratio δ_f. The engine's parameters are partial widths Γ_λ,fL (eV) per final state f and multipolarity L,
turned into reduced width amplitudes before the calculation. This note defines a second parametrisation of the same
partial widths, written as one record per level in a project of format 3 and in pyazr, which a fit varies and reads
directly, and the map from it to the engine's parameters. The physics of the engine (level matrix, amplitudes,
observables) is not changed: it receives the same partial widths as before.

## 1. What exists

**The file.** In `<levels>` (31 columns, `data_formats.rst`, read by `io::LevelRow`) each level has one line per
channel; a photon channel's line names the photon pair (column 6: the pair key; the pair's columns 13-17 give the final
level's J^π and E_x), 2s = 2J_f (column 7), 2l = 2L (column 8), and column 12 is the **observed** partial width
Γ_λ,fL (eV) of that multipolarity; E or M follows from the parities. Column 11 fixes it. The sign of the entry is the
sign of the reduced width amplitude: it sets the interference of this level's amplitude to f with the amplitudes of
the other levels (and of external capture) to the same f; it has no meaning for a level alone, except between the two
multipolarities of one f, where it sets the sign of the mixing ratio (section 2.3).

**The engine's parameters.** `CNuc::TransformIn` turns each observed width into a reduced width amplitude
(level_parameters.md, Eqs. 2.1-2.3): in Brune's basis γ̃²_λc = (|Γ_λc|/P_c) N_λ/D_λ with the photon penetrability
P^γ_fL = (E_γ/ħc)^{2L+1} (fm^−(2L+1); 1 with `--use-rmc`), E_γ = E_λ − E_f; the photon channels add nothing to N_λ and
D_λ, so that

    Γ_λ,fL = 2 P^γ_fL(Ẽ_λ) γ̃²_λ,fL / n_λ,     n_λ = 1 + Σ_{c particle} γ̃²_λc S'_c(Ẽ_λ)                          (1.1)

(Γ in MeV; n_λ = 2N_λ/D_λ). In Park's basis the amplitude is the observed one, Γ_λ,fL = 2 P^γ_fL (γ^P_λ,fL)². In the
formal basis (`--no-brune`) the formal photon amplitudes of a J^π block are linear combinations of the Brune amplitudes
of all its levels (`formal_from_brune`). The photon widths stay out of the level matrix and out of the Brune factor
n_λ (`docs/source/theory/conventions.rst`). The **fit parameters** are these amplitudes (MeV^½), one per channel of
each level in the vector's fixed order (`CNuc::FillMnParams`): `energy_<n>`, `width_<n>_<k>` (n the level counted over
all J^π groups, k the channel of its group). `param.par` / `param.sav` list them by name, `parameters.out` gives the
observed widths and amplitudes per level (`CNuc::TransformOut`), `--describe` the records of the file, pyazr the same
through `Session.parameters` (kind `width`).

**External capture.** With external capture on, the width written for a photon channel to a final state with an
external-capture mask is the observed total; `TransformIn` subtracts Re γ^ext (external_capture.md, Eq. 4) and the fit
parameter is the internal amplitude; the observed width is 2P|γ_int + γ^ext|²/n_λ. γ^ext depends on the level's
particle amplitudes, on E_λ and on the final state's amplitudes, so the observed width of such a channel is a function
of parameters outside the photon channel.

**Priors on a parameter.** A Gaussian prior on one parameter is written in `<parameterSettings>` (nuisance flag, the
centre and width of a level width in eV, converted to the amplitude) or as a multivariate prior in a file
(`priorCovariance=`). A literature value of a level property is a `<modelData>` datum (term T9).

**The datum kinds** (`src/ModelData.cpp`, `LevelValue`) compute from the observed widths of the current parameters:
`partial_width` Γ_λ,pair (sum over the pair's channels), `total_width` Γ = Σ over **all** channels that are widths,
photon channels included, `branching_ratio` Γ_exit/Γ with that Γ, `resonance_strength` ω Γ_in Γ_out/Γ with that Γ.
Section 4 compares this with the engine's calculation.

## 2. The parametrisation

### 2.1 Definitions

For a level λ with photon channels to final states f = 1..n (each with one multipolarity L_f, or two, L_f and L_f+1):

    Γ_γ = Σ_f Γ_f,     Γ_f = Γ_f,L + Γ_f,L+1,     b_f = Γ_f / Γ_γ,     Σ_f b_f = 1                            (2.1)

    δ_f² = Γ_f,L+1 / Γ_f,L                                                                                    (2.2)

Γ_γ (eV) is the sum of the level's observed radiative partial widths, exactly the quantity the engine sums when it
computes the observed widths; it is the radiative width of the level alone. It does **not** enter the total width of
the resonance denominator: in the engine (level matrix, Brune factor n_λ) the total width of a resonance is the sum of
its particle widths, and the radiative widths of this record keep that meaning (the PI, 2026-11-12). b_f are the
branchings of the γ decay (not of the total width).

### 2.2 The map to the engine's partial widths

    Γ_f,L   = σ_f       Γ_γ b_f / (1 + δ_f²)
    Γ_f,L+1 = σ_f s(δ_f) Γ_γ b_f δ_f² / (1 + δ_f²)                                                          (2.3)

σ_f = ±1 is the interference sign of the final state (the sign of its L amplitude, kept as it is today: the record's
`sign`, default +), and s(δ) = −sign δ is the sign that the engine's amplitudes give the Krane-Steffen mixing ratio
(section 2.3). For a final state with one multipolarity Γ_f,L = σ_f Γ_γ b_f. With the remainder b_r = 1 − Σ_{f≠r} b_f
(section 2.6) every partial width is a product of the record's numbers: the map has no free convention left.

### 2.3 Sign of the mixing ratio: the numerical check

The sign is fixed against the engine's own γ angular distributions, not taken from a book. A level of ¹⁵O formed by
p + ¹⁴N in one channel (l, s) with unpolarised beam and target decays to one final state by L and L+1. The engine's
dσ/dΩ (c.m., 0-180° in 10° steps, base engine `b-head3`) is expanded in Legendre polynomials and a₂/a₀, a₄/a₀ are
compared
with Krane, Steffen and Wheeler (At. Data Nucl. Data Tables 11 (1973) 351), W(θ) = Σ_k B_k(J) A_k(δ) P_k(cos θ),
A_k(δ) = [F_k(LLJ_fJ) + 2δ F_k(L L+1 J_fJ) + δ² F_k(L+1 L+1 J_fJ)]/(1+δ²), B_k from the m-populations of formation in
(l, s) with m_l = 0 (the normalisation of B_k is fixed by the pure-multipole runs, which carry no sign). Each case is
run
four times: pure L, pure L+1, and Γ_L+1/Γ_L = 0.25 with the L+1 width written positive and negative (L width positive).

| case | J^π → J_f^π (E_f, MeV) | L / L+1 | engine, Γ_L+1 > 0: a₂/a₀, a₄/a₀ | KSW δ = −0.5 | KSW δ = +0.5 |
|---|---|---|---|---|---|
| a | 3/2⁺ → 1/2⁺ (5.183) | M1/E2 | +0.39282, 0 | +0.39282, 0 | −0.99282, 0 |
| b | 3/2⁺ → 5/2⁺ (5.2409) | M1/E2 | −0.62471, 0 | −0.62471, 0 | +0.32186, 0 |
| c | 3/2⁺ → 1/2⁻ (0) | E1/M2 | +0.39282, 0 | +0.39282, 0 | −0.99282, 0 |
| d | 3/2⁺ → 3/2⁺ (6.7931) | M1/E2 | +0.93968, 0 | +0.93968, 0 | −0.29968, 0 |
| e | 5/2⁺ → 3/2⁻ (6.1763) | E1/M2 | +0.53216, +0.13061 | +0.53216, +0.13061 | −1.09053, +0.13061 |
| f | 5/2⁺ → 7/2⁺ (7.2759) | M1/E2 | −0.77752, +0.02177 | −0.77752, +0.02177 | +0.41017, +0.02177 |
| g | 3/2⁻ → 1/2⁺ (5.183), p wave | E1/M2 | +0.39282, 0 | +0.39282, 0 | −0.99282, 0 |

In all seven cases the run with both widths of the same sign is KSW's δ < 0, and with opposite signs δ > 0; the pure
runs and the matched mixed runs agree with KSW to 2×10⁻⁷ in a₂/a₀ and a₄/a₀ (the Legendre fit of 19 angles). Hence

    sign(δ_KSW) = − sign(Γ_f,L) · sign(Γ_f,L+1)                                                               (2.4)

for M1/E2 and E1/M2, both parities of the level, J up and down, J = 3/2 and 5/2. The record's δ **is Krane and
Steffen's** (the convention of ENSDF and of most γ-ray work). Rose and Brink's δ is related to it by a phase that
Krane and Steffen give in their paper; it was not checked numerically here, so the record does not offer it.
E2/M3 (L = 2) was not run. Scripts: `make_cases.py` and `ks_check.py` (sympy's 3j/6j) in the scratch directory of this
branch (`work/gammab/sign/`); the test `tests/gamma_decay/check.sh` carries one case (a) as a regression.

### 2.4 The map to the engine's parameters

The fit parameter of a photon channel is an amplitude, so the record's numbers are turned into amplitudes at every
evaluation, from the partial widths (2.3) and the current values of the level's other parameters: in Brune's basis by
the inverse of (1.1), γ̃_λ,fL = sign(Γ_f,L) [|Γ_f,L| n_λ / (2 P^γ_fL(Ẽ_λ))]^½, with n_λ from the level's current
particle
amplitudes and Ẽ_λ its current energy; in Park's basis γ^P = sign(Γ) [|Γ|/(2P^γ)]^½. Both are the engine's own
conversion, `physics::amplitudes_from_observed` with the particle columns given as amplitudes: one implementation.
At set-up the record's partial widths are written into the level's input widths and `TransformIn` converts them as it
converts any width, so a project with the record and the same project with the partial widths written by hand give
the same calculation (the map is a product: the last digit).

Refused: the formal basis (`--no-brune`): there a level's formal photon amplitude mixes the Brune amplitudes of all
levels of its block, so the record of one level cannot be the parameters of its own slots. Also refused: a photon
channel to the level's own state (a static moment). A photon channel with external capture in the run is mapped as
section 2.7 says.

### 2.5 Jacobian

With q the record's parameters of level λ and p_k = γ_λ,k the amplitudes of its photon channels:

    ∂p_k/∂q = (p_k / 2Γ_k) ∂Γ_k/∂q   (Γ_γ, b; Γ_k ≠ 0; 0 for Γ_k = 0),      ∂Γ_k/∂q from (2.3), analytic
    ∂p_L/∂δ = −p_L δ/(1+δ²),   ∂p_L+1/∂δ = K/(1+δ²)^{3/2},   K = 2^½ p_L+1 at δ = 1    (regular through δ = 0)
    ∂p_k/∂γ_c = p_k γ_c S'_c / n_λ   (particle channel c of λ, Brune; 0 in Park's basis)
    ∂p_k/∂E_λ   by a central difference of the map (step 10⁻⁶ MeV)                                          (2.5)

The gradient of χ² in the record's parameters is the chain rule on the engine's existing gradient in the amplitudes,
g_q = Σ_k g_{p_k} ∂p_k/∂q, and the level's energy and particle amplitudes get g_{p_k} ∂p_k/∂E, ∂p_k/∂γ_c added: the
amplitudes of a photon channel move when Γ_γ is held and the energy or n_λ changes. The covariance of the record's
parameters is the one the fit gives (MINUIT's, or the MCMC's); the error of a derived number (the remainder, a partial
width) is propagated with ∂/∂q of (2.3) to the linearisation.

### 2.6 The constraint Σ b_f = 1, and the boundary

**Choice: n − 1 free branchings and a remainder.** One final state, named by the user, carries the remainder
b_r = 1 − Σ_{f≠r} b_f; it is not a parameter. Reasons: (i) the record then has exactly one parameter per photon channel
of the level (1 + (n−1) + n_mixed = number of photon channels), so the record's parameters take the positions of the
amplitudes in the vector and every path that addresses a level parameter by position (gradient, Jacobian, MCMC,
`param.par`) keeps its layout; (ii) no redundant scale: with all n branchings and a normalisation the χ² is flat along
the scale and the Hessian is singular; (iii) a reader sees the numbers that were measured, and the remainder is the one
the experiment did not measure (often the weak or unobserved branch).

**Boundary.** Each free b_f has MINUIT limits [0, 1], Γ_γ ≥ 0, δ is free. MINUIT reaches a limit through its sine
transformation without leaving it; MINOS reports the limit as the bound of an interval that hits it; the ensemble
MCMC samples in the external coordinates and rejects a proposal outside the limits. The remainder has no box: when
Σ_{f≠r} b_f > 1 the map uses b_r = 0 and the others divided by their sum (no negative width, no NaN), and the objective
adds the wall ((Σ_{f≠r} b_f − 1)/10⁻³)² so that the minimum stays on Σ ≤ 1 (to 10⁻³ of the wall's width); a fit that
ends
with a branching or the remainder at 0 says so in the run log and in `parameters.out`. At b_f = 0 the amplitude is 0
and ∂p/∂b is infinite (p ∝ b^½); the Jacobian is set to 0 there, as for a width written 0. δ = 0 is regular: the L+1
amplitude is linear in δ, and its derivative is taken in that form, not through p/(2Γ) (0/0 at δ = 0), so δ = 0 is an
allowed start; with angle-integrated data alone, which depend on δ², the sign of δ is not determined and δ = 0 is a
stationary point by symmetry: start away from 0 or give an angular distribution.

### 2.7 A channel with external capture

**What the record's width is.** With external capture on, a photon channel k of λ (multipolarity XL, final state f)
has an external part when f is the external-capture level of k's photon pair and f's mask (column 31) covers XL
(`ExternalCaptureFinalState`, `src/CNucExternalCapture.cpp`). Its width in `<levels>` is the observed total, and
`TransformIn` turns it into the internal amplitude γ_int (section 1). The record's partial width Γ_k of (2.3) is the
same observed total: a b_f measured as the fraction of Γ_γ that goes to f is what the record writes, the external part
included. The observed width of the engine is (cnuc_split.md Eq. 8.1)

    Γ_k = 2 P_k |γ_int,k + γ^ext_k|² / n_λ   (Brune; n_λ = 1 in Park's basis),
    γ^ext_k = − Σ_c Σ_{c_f} s_{f c_f} γ_λc I_{c c_f XL}(E_λ, E_f),   s_{f c} = γ̃_fc N_f^−½ K_c          (2.7)

(external_capture.md Eqs. 2 and 4; I the complex radial integral of channel capture, K_c the ANC conversion,
N_f = n0 + Σ_c A_c γ̃_fc² I_fc with A_c = 2 μ_c a_c/ħ², n0 = 1 in Brune's basis and Park's J_f = 1 − Σ_c γ_fc² S'_c(E_f)
in Park's, as `CNuc::CalcShiftFunctions` forms it). γ^ext depends on E_λ and the particle amplitudes of λ, and on E_f
and the particle amplitudes (ANCs) of f; all of them may be parameters of the vector.

**The map**, at every fill (`ToEngine`), for each record channel with an external part:

1. the total amplitude g_k = sign(Γ_k) [|Γ_k| n_λ / (2 P_k)]^½ of section 2.4 (unchanged, `amplitudes_from_observed`);
2. γ^ext_k at the **current** parameters: E_λ and γ_λc from the vector, E_f and γ̃_fc from the vector, N_f from them
   (`ExternalReducedWidthAt`, of which the set-up's `ExternalReducedWidth` is now the call with the level's stored
   values: one implementation);
3. γ_int,k = sign(g_k) (g_k² − (Im γ^ext_k)²)^½ − Re γ^ext_k, by `physics::subtract_external_capture`, the function
   `TransformIn` uses.

So the record's observed widths stay what the record says when a particle width, E_λ or an ANC of f moves; the twin
(the same level with partial widths in `<levels>`) has the same widths at its set-up point and keeps its γ_int fixed
afterwards. Both parametrisations reach every real γ_int (the two signs of g_k give the two branches), so a fit has the
same minimum in both.

**Im γ^ext and the floor.** Im γ^ext enters through the modulus: |γ_int + γ^ext|² ≥ (Im γ^ext)², so no internal
amplitude gives a width below the floor 2P(Im γ^ext)²/n_λ. A record width below it (g_k² < (Im γ^ext_k)²) is handled
as `TransformIn` handles it: γ_int = −Re γ^ext, the observed width is the floor, not the record's. The fit is not
stopped; a warning is printed once per run, and `parameters.out` notes the channel. No wall term is added (choice):
the floor is a property of the physics, not of the parametrisation (the twin cannot go below it either), and the
map is continuous there; on the floor side the record's numbers have no gradient in that channel.

**The gradient's pull-back** (Eq. 2.5 extended). With a = Re γ^ext_k, b = Im γ^ext_k and r = (g_k² − b²)^½,

    R_k = |g_k| / r,   B_k = − sign(g_k) b / r     (R_k = 1, B_k = 0 without external part; 0, 0 on the floor)
    ∂γ_int,k/∂q     = R_k ∂g_k/∂q                                    (the record's numbers; ∂g/∂q of Eq. 2.5)
    ∂γ_int,k/∂γ_λc  = R_k ∂g_k/∂γ_λc − Re ∂γ^ext_k/∂γ_λc + B_k Im ∂γ^ext_k/∂γ_λc
    ∂γ_int,k/∂γ̃_fd = − Re ∂γ^ext_k/∂γ̃_fd + B_k Im ∂γ^ext_k/∂γ̃_fd                                          (2.8)
    ∂γ^ext/∂γ̃_fd   = Σ_c ∂γ^ext/∂s_fc (δ_cd N_f^−½ K_c − s_fc ∂N_f/∂γ̃_fd / (2 N_f)),
    ∂N_f/∂γ̃_fd     = 2 γ̃_fd (A_d I_fd − [Park] S'_d(E_f))

γ^ext is bilinear in γ_λc and s_fc, so ∂γ^ext/∂γ_λc and ∂γ^ext/∂s_fc are exact: Eq. (4) evaluated with unit vectors and
the integrals of the same fill (no further integral). ∂N_f/∂γ̃_fd uses the coefficient A_d of
`physics::bound_state_coefficient`, as the capture term's gradient does (`AMatrixFunc`). The energy E_λ keeps its
central difference of the whole map (step 10⁻⁴ MeV), γ^ext included, as in section 2.5: the radial integrals have no
derivative in the code. E_f has no term: a bound state's energy is never a free parameter (`CNuc::FillMnParams`).
The fill stores, per record channel, g_k, γ^ext_k and the two slopes; the pull-back reads them (no integral), except
the energy's difference (two further maps), taken only when E_λ is free. The gradient's mask asks for the photon
amplitudes' cotangents when any input of the record is free: its numbers, E_λ, γ_λc and the particle amplitudes of
the final states with an external part (`GammaDecaySet::InputIndices`).

**`parameters.out`.** The record's block gives Γ_γ, b_f, δ_f and the partial widths of (2.3); these are the observed
widths `TransformOut` prints for the level's channels (`observed_levels` with γ^ext at the same parameters), to the
rounding of the two evaluations. A channel on the floor is noted with the engine's width.

**Where and at what cost.** γ^ext is evaluated by the fill of the compound (`GammaDecaySet::ToEngine`, called from
`CNuc::FillCompoundFromParams`) with the run's `Config`, which the set refers to (the run's object, not a copy; the
set is set-up state of the run's compound and its clones). The radial integrals are `ECIntegral`'s direct ones with
the run's per-thread scratch (`RunScratch`: the entrance Coulomb values at the quadrature's radii, the Whittaker
values, the GSL workspace; keyed by content), the same evaluation as in `TransformIn` and `TransformOut`. The energy
tables of `ECTable` are not used (choice): their interpolation would make the record's observed widths differ from
`TransformOut`'s and from the twin's by the table's error. No static or process-wide memo is added (check 38). Cost
per fill: one integral per coupled pair (c, c_f) of each record channel with an external part; two more maps per
gradient when E_λ is free. Measured (user-space instructions, one thread, the ¹⁵O project of the test with capture
into three final states, 90 data points): a calculation +0.25 % against the twin, a MIGRAD fit with the gradient
from displaced starts +14 % (89 objective calls against 85). Threads: each fill runs in its thread with its own
scratch, whose values are keyed by
everything they depend on, so the numbers do not depend on the thread count.

**`--ignore-externals`.** A record channel whose set-up width is 0 has no external part, as `TransformIn` decides on
the input (external_capture.md section 3).

## 3. The file and the API

### 3.1 The record

One `gamma_decay` record per level, in a section `<gammaDecays>`:

    <gammaDecays>
    gamma_decay level5
        level        5
        J            1/2+
        total_width  0.0385 eV fit prior 0.040 0.004
        branching    g+15O[6.7931] 0.230 fit prior 0.23 0.02
        branching    g+15O[6.1763] 0.584 fit
        branching    g+15O[6.8594] 0.0014 fixed
        remainder    g+15O[5.183] prior 0.17 0.02
        sign         g+15O[6.1763] -
        mixing       g+15O[6.7931] -0.12 fit
    </gammaDecays>

`level` and `J` name the level as a datum does (its ID, column 9, and its J^π, F13). `total_width` is Γ_γ with an
energy unit; `branching` one final state (the photon pair by its label) and b_f; `remainder` the final state that takes
1 − Σ b_f, with an optional Gaussian `prior <centre> <width>`; `mixing` δ_f (Krane-Steffen) of a final state with two
multipolarities; `sign` the interference sign σ_f of a final state (default +). Every number is
`<value> fixed|fit [prior <centre> <width>[%]]` (the nuisance grammar of format 3). A final state of the level's photon
channels that the record does not list has b_f = 0, fixed; a final state with two multipolarities and no `mixing` has
δ = 0, fixed (pure L). A level with one final state needs no `branching`: its b = 1.

The level's photon lines of `<levels>` stay (they define the channels) with width 0; their fixed flags are not read.
Refusals, with a code each (`reference/format3.rst`): a record level with a non-zero photon width in `<levels>`
(`E_GAMMA_MIXTURE`: the message names the line and says to write 0 there or to delete the record); branchings that are
all fixed and do not sum to 1 (`E_BRANCHING_SUM`, tolerance 10⁻⁶); fixed and given branchings that leave a negative
remainder (`E_BRANCHING_REMAINDER`); a free branching without a remainder (`E_BRANCHING_SUM`); a `mixing` on a final
state with one multipolarity (`E_MIXING_ONE_MULTIPOLE`); a final state that is no photon channel of the level
(`E_NOT_A_CHANNEL`); two records for one level (`E_DUPLICATE_RECORD` of the validator, by name, or
`E_GAMMA_DUPLICATE`). The engine refuses at set-up, with a message, the formal basis (section 2.4).

### 3.2 Parameters and outputs

The record's numbers become parameters of one new kind, `radiative` (code 13), named `<record>.total_width`,
`<record>.branching.<final state>`, `<record>.mixing.<final state>` (`level5.branching.g+15O[6.7931]`), in the vector at
the positions of the level's photon amplitudes, in this order: Γ_γ, the explicit branchings and the mixing ratios in the
order of the channels of the level's J^π group. `param.par`, `param.sav` and pyazr list them by these names; a fixed
number is a fixed parameter. A prior is a term of the objective (`gamma decay priors`, beside the remainder's wall),
including the remainder's prior, which is a prior on 1 − Σ b_f. `parameters.out` adds a block per record with Γ_γ, every
b_f (the remainder included), δ_f and the partial widths they imply; the observed widths of the level's channels are
written as before.

### 3.3 pyazr and a form

pyazr reads and edits the record as any record of format 3 (`Project.add("gamma_decay", ...)`, `project["level5"]`,
`project.set("level5", branching=[...])`), and the session lists its parameters (`Parameter.kind == "radiative"`).
To fit one branching: `project.set("level5", total_width=Nuisance(0.0385, "fit", unit="eV"))`. A form for the web UI:
one panel per level that has radiative channels, a switch "partial widths / total width and branchings"; in the second
mode a row for Γ_γ and one row per final state (label, b_f or "remainder", δ_f where two multipolarities exist, sign),
each number with fixed/fit and an optional prior; the sum of the branchings shown live.

## 4. The total width in the datum kinds and in the calculation

The engine computes a resonance with the total width of its denominator equal to the sum of the particle widths
(`n_λ` and the level matrix contain no photon width): the strength of a narrow resonance in a calculated yield is
ω Γ_γ Γ_p / Σ Γ_particle. The datum kinds of section 1 use Γ = Σ of all observed widths, photon widths included:
a `resonance_strength` datum compares ω Γ_in Γ_out / (Σ Γ_particle + Γ_γ) with its value, so the datum and the folded
yield of the same level are two different numbers, which differ by the factor Γ/Σ Γ_particle = 1 + Γ_γ/Σ Γ_particle.
For the 259.7 keV level of ¹⁴N(p,γ)¹⁵O as written in `tests/thick_target_resonance` (Γ_p = 974.9 eV, Γ_γ = 0.0385 eV)
the factor is 1 + 4.0×10⁻⁵; in the director's check with Γ_p = 9.7 eV it is 0.39 % (`work/director/omegagamma-223120`);
for Γ_p = Γ_γ it is 2. `total_width` and `branching_ratio` follow the same Γ. This is stated in the user section; what
the datum computes is not changed here.

## 5. Decisions for the director

1. **The section and record name** (`<gammaDecays>`, `gamma_decay`) and the field names. Recommendation: keep; they
   read as the physics.
2. **The formal basis**: refused for a project with the record. Recommendation: keep refused (the formal basis is a
   check mode); a level-by-level map there would need the block's whole transformation per evaluation.
3. **External capture on a record's channel**: decided (director, 2026-11-13) to lift it with γ^ext evaluated in the
   map; done in section 2.7 (Brune and Park bases).
4. **The remainder's wall** (width 10⁻³ in Σ b) against a hard refusal during a fit. Recommendation: the wall.
5. **Rose-Brink δ**: not offered (section 2.3). Recommendation: add only if asked, with its own numerical check.

## 6. Where the code goes, and the tests

- `experiment/`: the record kind (`Schema.cpp`: `GammaDecayKind`, section `gammaDecays`), its build against `<levels>`
  (`GammaDecay.h/.cpp`: `BuildGammaDecays`, the refusals), the codes (`Diagnostics.h`).
- engine: `include/GammaDecayMap.h`, `src/GammaDecayMap.cpp`: reading the section (as `ModelData::Read`), the input
  widths at set-up, the parameters (`FillMnParams`), the map at every fill (`FillCompoundFromParams`), its Jacobian and
  the pull-back of a gradient (`GradAccum` before `Scatter`), the priors and the wall (objective term), the report.
  One call each in `CNuc` (fill, parameters, fill from parameters), in the gradient's three scatter sites, and the
  helper `IsLevelSlot(kind)` where the code asked "energy or width" of a parameter. External capture on a record
  channel (section 2.7): `ExternalCaptureFinalState` and `ExternalReducedWidthAt` in `src/CNucExternalCapture.cpp`,
  the record's inputs for the gradient's mask (`InputIndices`).
- Tests (`tests/gamma_decay/check.sh`): a project with partial widths and the same with the record, built by a script
  of the test from (2.3), calculations identical; pseudo-data from known Γ_γ, b_f fitted in both parametrisations
  (equal χ² minima, the errors of b_f equal the propagated ones), a prior's pull; the gradient against central
  differences with `--use-gradient` and `--use-lm`; a branching driven to 0; the refusals; pyazr (built in Python, an
  agent-style read-change-refit).
