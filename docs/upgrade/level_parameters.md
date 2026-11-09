# Level parameters: bases and conventions (`level_parameters`)

The collision matrix of a J^π block is computed from level energies E_λ (MeV, excitation energy E_x of the compound
nucleus) and reduced width amplitudes γ_λc (MeV^½). A project states these numbers in its own form: observed partial
widths Γ (eV) and asymptotic normalization coefficients C (fm^−½) in Brune's or Park's parametrization, or formal
amplitudes, coupled in some angular-momentum scheme with some phase convention. `level_parameters`
(`physics/include/azure2/physics/LevelParameters.h`) converts them to the one internal set: formal R-matrix parameters
with constant boundary values B_c, Lane–Thomas phases, channel-spin coupling s = j₁ + j₂ (light particle first),
J = s + l (`docs/source/theory/conventions.rst`). It is the only place in the code where a conversion between bases or
conventions happens (`RULES.md` §1.6). This page gives the equations of each conversion, their validity, the
tolerances, and what each test would detect. The engine calls it: `CNuc::CalcBoundaryConditions`, `CNuc::TransformIn`,
the formal → Brune iteration and the observed widths of `CNuc::TransformOut`, `CNuc::ConvertAmplitudeBasis` and the
external-capture final states of the formal basis delegate to these functions (`cnuc_split.md`), with this
page's rules for the formal basis (B_c at the lowest level, orientation of the Brune level) and its own eigensolvers; the
convention conversions (section 5) are not used by the engine yet, which reads no convention record. The `CNuc.cpp`
line numbers below are those of e54e6fcc, the code these functions were checked against before the move.

## 1. Symbols

| symbol | unit | meaning |
|---|---|---|
| E_λ, Ẽ_λ | MeV | formal and observed (Brune) level energy, E_x of the compound nucleus |
| E_c = (E − S_α) − (e₀ + e₁) | MeV | channel energy, c.m.; S_α separation energy, e₀, e₁ excitation energies of the pair |
| γ_λc, γ̃_λc, γ^P_λc | MeV^½ | formal, Brune and Park reduced width amplitudes |
| Γ_λc | eV | observed partial width (signed: the sign of the amplitude) |
| C_λc | fm^−½ | asymptotic normalization coefficient of a bound channel (signed) |
| P_c(E) | — | penetrability (particle channel, E_c > 0) |
| P^b_c(E) = μ a_c / (ħ² W_l(2κa_c)²) | MeV^−1 fm^−1 | the ANC conversion factor (E_c < 0); W Whittaker function, κ the bound-state wave number |
| P^γ_c(E) = (E_γ/ħc)^{2L+1} | fm^−(2L+1) | photon channel of multipolarity L (1 with `--use-rmc`) |
| S_c(E), S'_c = dS_c/dE | —, MeV^−1 | shift function and its energy derivative |
| B_c | — | boundary value of channel c in the block |
| J_λ | — | Park's norm of level λ |

The channel functions P, S, S′ are an input of every function (`ChannelFunctions`): `physics/` has no Coulomb
functions until `channel_functions` exists (pycoulomb, `ROADMAP.md` Phase 1). The engine's adapter
(`src/CNucLevelParameters.cpp`) supplies them from `ChannelFunc` and `WhitFunc`: every call returns S, P and S′ of the
channel, with E_c formed as in the table, in the order the collision matrix forms it at a data point (`EPoint.cpp`).
The engine's S is accurate to about 10⁻¹⁰ (`coulsm`) and its S′, a Richardson difference of S with step 10⁻³ E_c, to
about 10⁻⁸; neither is smooth below that level, so a change of E_c by one ulp can move them there (`cnuc_split.md`
section 6, step 1b).

## 2. Observed widths → Brune and Park amplitudes

For each level λ, over the particle channels of the block (the photon widths stay out of the sums, as in the engine):

    N_λ = 1 + Σ_{c ∈ rwa} γ̃²_λc S'_c(Ẽ_λ),     D_λ = 2 − Σ_{c ∉ rwa} |Γ_λc| S'_c(Ẽ_λ) / P_c(Ẽ_λ)            (2.1)

    γ̃²_λc = (|Γ_λc| / P_c) · N_λ / D_λ,       sign γ̃_λc = sign Γ_λc                                        (2.2)

with Γ in MeV (the file's eV divided by 10⁶). A bound channel uses C² in place of |Γ| and P^b in place of P. A channel
flagged `gammaIsRWA` (rwa) carries γ̃ itself. A photon channel uses P^γ (floored at 10⁻¹⁶) and adds nothing to the sums;
a static moment (photon channel to the level's own state, |E − E_f| < 1 keV, same J^π) uses the engine's factor
3J/[4(J+1)μ_N²] (M1) or 60J(2J−1)/[(J+1)(2J+3)] (E2, input × 100 (αħc)^½) and reverses the sign for half-integer J.
A β channel carries its amplitude. Eq. (2.2) is the inverse of Γ_λc = 2P_c γ̃²_λc / (1 + Σ_c γ̃²_λc S'_c), which
follows from Brune's level matrix (PRC 66, 044611 (2002)) for an isolated level; because 1 + Σ γ̃² S' = 2N/D, the
inverse (`observed_from_amplitudes`) is exact.

Park's amplitudes (PRC 104, 064612 (2021)) are the observed ones, Γ = 2P(γ^P)², with norm

    J_λ = 1 − Σ_c (γ^P_λc)² S'_c = D_λ / (2 N_λ) = 1 / (1 + Σ_c γ̃²_λc S'_c),      γ^P_λc = J_λ^½ γ̃_λc         (2.3)

for every channel, rwa and β included (`park_from_brune`, `brune_from_park`). Validity: D_λ > 0 (Brune) and J_λ > 0
(Park); otherwise no real reduced width reproduces the widths, and a warning names the level (the engine's message).
A level exactly at a channel threshold, or just above a Coulomb threshold where P = 0 in double precision, gets
γ = 0 in that channel, with a warning when the input was not 0. Code: `CNuc.cpp:1200-1460` (e54e6fcc).

## 3. Boundary values and external capture

In the formal basis every particle channel of a block has the constant

    B_c = S_c(Ẽ_low),     Ẽ_low = the lowest level energy of the block (MeV, E_x)                                (3.1)

(HANDOFF decision 10; the first of equal energies). Observables do not depend on B_c; the formal parameters do, and with
Ẽ_low they do not depend on the order in which the file lists the levels: a permutation of the listing permutes the
formal levels. Because Ẽ_λ is an eigenvalue of the formal level matrix with B_c = S_c(Ẽ_λ), the lowest level's formal
energy equals its Brune energy. Builds before 2026-11-07 took the first level listed; formal parameters they wrote for
a block listed out of energy order belong to that other B_c, and the engine says so at load in the formal basis
(`cnuc_split.md` section 7). A block without levels has B_c = 0. Non-particle channels copy channel 0's value, as in
`CNuc::CalcBoundaryConditions` (`CNuc.cpp:1975-2003`); they never enter the level matrix.

When external capture is on, the photon amplitude of the file is the total one; the internal amplitude is

    γ^int = sign(γ̃) (γ̃² − (Im γ^ext)²)^½ − Re γ^ext,    or −Re γ^ext if γ̃² < (Im γ^ext)² (with a warning)     (3.2)

(`subtract_external_capture`, `CNuc.cpp:1488-1515`). γ^ext is an input (`CNuc::CalcExternalWidth`); it belongs to
the external-capture module, not to `level_parameters`.

## 4. Brune ↔ formal

With the Brune parameters Ẽ_λ, γ̃_λc and the constant B_c, the matrices (λ ≠ μ in the second lines)

    M_λλ = 1,              M_μλ = −Σ_c [S_c(Ẽ_μ) − S_c(Ẽ_λ)] / (Ẽ_μ − Ẽ_λ) γ̃_λc γ̃_μc
    N_λλ = Ẽ_λ + Σ_c [S_c(Ẽ_λ) − B_c] γ̃²_λc,
    N_μλ = Σ_c {[Ẽ_μ S_c(Ẽ_λ) − Ẽ_λ S_c(Ẽ_μ)] / (Ẽ_μ − Ẽ_λ) − B_c} γ̃_λc γ̃_μc                                   (4.1)

define the generalized eigenproblem N a = E M a, a^T M a = 1 (Brune 2002). Its eigenvalues are the formal energies;
the formal amplitudes are γ_λc = Σ_μ a_μ,k(λ) γ̃_μc. Level λ keeps its position and receives the eigenvalue of rank
k(λ), the rank of Ẽ_λ among the block's energies (`formal_from_brune`, `CNuc.cpp:1532-1601`). Validity: M positive
definite; otherwise no formal parameters with constant B_c exist and the function throws, as the engine stops.
Distinct Ẽ_λ are required.

The eigenproblem leaves the sign of each eigenvector free. Of the two, `formal_from_brune` takes the one whose
coefficient on the level's own Brune state is positive, a_λ,k(λ) > 0: each formal level keeps the orientation of the
Brune level it comes from. The rule does not depend on the eigensolver (`SymmetricEigen.h`'s is the only one), and
where the levels of a block mix weakly (a ≈ 1) the formal amplitudes carry the signs of the input. Within the R-matrix
a level's sign is free (γ^T A γ is invariant), but external capture is not: the photon amplitudes of other levels
that end on a final state f, and the external-capture amplitude ∝ γ_f, carry f's sign together. The rule "largest
|γ_λc| of each formal level positive" (the proposal of decision 10) breaks that link: on `tests/beam_profile_kernel`
it reverses the 6.917 MeV 2⁺ final state, whose largest formal amplitude is a photon one, and the formal-basis capture
to that state changes by a factor of up to 4.5 (`cnuc_split.md` section 6, step 1b). The orientation rule leaves it
unchanged.

The inverse (`brune_from_formal`) is the engine's fixed-point iteration (`CNuc::TransformOut`, `CNuc.cpp:2563-2745`):
for each level, B_c is moved to S_c at the current energy of that level, the matrix
E_μ δ_μμ' − Σ_c ΔB_c γ_μc γ_μ'c is diagonalized, and the iteration stops when the level's energy moves by less than
the tolerance, 10⁻⁶ MeV in the engine (an option here). A level that has not converged after 1000 iterations keeps its
formal values, with a warning. Each eigenvector is oriented as in `formal_from_brune` (its component on its own level
positive), so Brune → formal → Brune returns the input signs.

## 5. Conventions: the matrix on the γ columns

A project's amplitudes γ^proj_λm in its columns m and the internal ones are related by

    γ_λk = Σ_m T_km γ^proj_λm,                                                                                (5.1)

with T block-orthogonal: one block per (partition α, l) of particle channels, ±1 for photon and β channels. T does not
depend on the parameters, so ∂γ/∂γ^proj = T (`ConventionTransform::jacobian`) and the inverse is T^T. Because P_c, S_c
and S'_c depend on (α, l) only and T is orthogonal within (α, l), the Brune normalization of Eq. (2.1) and Park's J are
invariant under T, and Eqs. (2.2)–(4.1) may be applied before or after the recoupling. The entries
(`recoupling_coefficient`, doubled spins in the code):

- **Spin order** (LS): |(j₂ j₁)s⟩ = (−1)^{j₁+j₂−s} |(j₁ j₂)s⟩.
- **Orbital order** (LS): |(l s)J⟩ = (−1)^{l+s−J} |(s l)J⟩.
- **jj coupling**, particle x's spin coupled with l (x = light: j = j₁ + l, J = j + j₂):

      ⟨((j₁ j₂)s, l)J | ((j_x l)j, j_y)J⟩ = (−1)^{4j_x + 3j_y + 3l − s − j} ŝ ĵ {j_x j_y s; J l j},   ŝ = (2s+1)^½   (5.2)

  times (−1)^{j₁+j₂−s} when x is the heavy particle (the internal s is (j₁ j₂)), and times (−1)^{l+j_x−j} for the order
  (l j_x)j. Eq. (5.2) was derived from the standard (ab)e,d ↔ a,(bd)f recoupling (Edmonds, Ch. 6) with two order swaps
  and is checked against m-scheme overlaps (section 7).
- **iˡ factor.** With iˡY_l^m surface functions a code's amplitudes differ from the internal ones, for real amplitudes,
  by T_cc' → i^{l−l'}T_cc'. For channels whose l have equal parity, i^{l−l'} = ε_c ε_c' with ε_c = (−1)^⌊l/2⌋; for
  opposite parity the remaining factor ±i is common to a pair of partitions and cancels from every observable of one
  reaction. The record's iˡ (and i^−l, which gives the same real signs) conversion is therefore the diagonal ε_c.
  This is the definition adopted here; no code with the iˡ convention has been compared with yet.
- **Photon multipole phases**: a source whose E or M amplitudes carry the opposite real sign sets `electric_sign` or
  `magnetic_sign` = −1. The engine's own E/M phases of external capture (`conventions.rst`, Eqs. conv-ec-e, conv-ec-m)
  are not amplitudes and are not converted here.
- **Coulomb-phase reference** (σ_l − σ₀ or σ_l): a phase common to a partition's U; T is the identity.

**Missing channels.** A recoupling is orthogonal only between complete sets: for (α, l) and J the internal set is every
s with |j₁ − j₂| ≤ s ≤ j₁ + j₂ and (s, l, J) a triangle, the jj set every j with (l, j_x, j) and (j, j_y, J) triangles.
`recoupling_gaps` lists every channel with a non-zero coefficient (> 10⁻¹⁴) that the other side lacks, every project
column that is not allowed or has no counterpart, and duplicates; `convention_transform` refuses (throws, naming them).
A block that lists one channel spin of two allowed ones (a model in which the other carries no width) has no jj form
without that channel: the missing one would also add its hard-sphere scattering.

## 6. External capture in the formal basis

A final state f enters external capture through its asymptotic normalization C_f ∝ γ̃_fc N_f^{−½}, the residue of its
own pole: its Brune amplitude γ̃_fc and energy Ẽ_f (`cnuc_split.md` §8). In the formal basis a final state that shares
its J^π group with other levels has a formal amplitude γ_fc that mixes the group's poles, except the lowest level of
the group: with B_c = S_c(Ẽ_low) its formal parameters are its Brune ones (the iteration of §4 moves B_c by 0). Since
2026-11-08 the engine uses γ̃_f and Ẽ_f in every basis: `CNuc::TransformIn` subtracts γ^ext (Eq. 3.2) of every group
while all levels hold their Brune amplitudes, before any group is made formal (until then a later group saw the formal
amplitude of a final state in an earlier one); in the formal basis a final state above the lowest level of its group
enters the capture term with `brune_from_formal` of its group at each evaluation (to 10⁻¹² MeV) and the capture
integrals with its Brune energy at setup. On `tests/beam_profile_kernel` the capture to the 6.049 MeV 0⁺ state (second
of four 0⁺ levels) differed between the bases by up to 0.45 relative and now agrees to the printed 11 digits
(`tests/beam_profile_kernel/check.sh`, with `--no-ignore-externals`). Not covered: a formal-basis fit that moves the
energy of a lowest-level final state away from its setup value (its Brune amplitude then differs from the formal one
at second order); the default `--ignore-externals`, which drops the external capture of a level whose internal photon
amplitude is below 10⁻⁸ MeV^½ in the basis in use and so selects different terms in the two bases.

## 6a. Observed quantities

`observed_levels` (Eq. 8.1 of `cnuc_split.md`) gives, from the Brune (Park) energies and total amplitudes
γ^int + γ^ext of a block's levels, the signed W = 2P|γ^tot|²/N per channel and the value a `<levels>` line holds
(Γ in eV, C in fm^−½, μ, Q, B, or the Brune amplitude of a `gammaIsRWA` column). It is Eq. (2.2) read backwards and
the one definition of these quantities: `CNuc::TransformOut` (`parameters.out`), `GetTransformParams` and
`CheckRadiativeWidths` read it, and `observed_from_amplitudes` is its value without external capture.

## 7. Tests and what they detect

`physics_angular_coupling` (`physics/tests/angular_coupling_test.cpp`): 3j, 6j, Clebsch–Gordan against closed forms
(stretched 3j, {a b c; 0 c b}, tabulated values), orthogonality, the 6j against its four-3j sum, and every recoupling
coefficient of the record against the overlap of the two states built in the m-scheme from Clebsch–Gordan products
(2j₁, 2j₂ ≤ 4, l ≤ 3, eight conventions). A wrong phase in Eq. (5.2), a transposed 6j argument or a wrong factorial
fails it at order 1. Results: On the GCC 11.5 build: 2326 exact values to 1.1 × 10⁻¹⁶; 7151
orthogonality sums to 2.2 × 10⁻¹⁶; 570 6j symbols against the four-3j sum to 1.1 × 10⁻¹⁶; 10 512 recoupling
coefficients in eight conventions against the m-scheme to 8.9 × 10⁻¹⁶.

`physics_symmetric_eigen`: the eigensolvers against their defining identities (A v = e v, V^T V = 1, N v = e M v,
v^T M v = 1) and a 2 × 2 closed form; refusal of an indefinite M. 2 × 2 closed form to 4.4 × 10⁻¹⁶; on 100 matrices with n ≤ 20 and
eigenvalues over six decades, |Av − ev| ≤ 3.7 × 10⁻¹⁵ max|A|, |V^T V − 1| ≤ 5.8 × 10⁻¹⁵, generalized residuals and
normalization ≤ 1.7 × 10⁻¹⁴.

`physics_level_parameters_unit` (synthetic 8Be-like block, non-linear S_c): (1) the identity behind section 4, that
the formal level matrix with constant B_c and Brune's level matrix give the same γ^T A γ at every energy; an error in
N, M, the rank assignment or B_c fails it at > 10⁻³. (2) Round trips observed → Brune → observed, Brune → Park →
Brune, formal → Brune → formal, signs included. (3) Convention invariance on six conventions; (4) reported gaps; (5)
the formal-basis rules: B_c at the lowest level (listed second), the lowest level's formal energy equal to its Brune
energy, the block listed in energy order giving the same formal levels permuted (with B_c at the first listed level
this fails), each formal level oriented like its Brune level, and widths too large for the radii refused. Results
(2026-11-07, GCC 11.5): identity (1) to 3.9 × 10⁻¹⁵ relative; lowest level's formal E = Brune E exactly; listing
order 1.8 × 10⁻¹⁶ (E), 3.6 × 10⁻¹⁶ (γ, signs included); observed → Brune → observed 1.5 × 10⁻¹⁶ and observed → Park →
observed 1.9 × 10⁻¹⁶ element-wise; Brune → Park → Brune 1.5 × 10⁻¹⁷; Brune → formal → Brune 2.2 × 10⁻¹⁵ (E), 1.8 ×
10⁻¹⁵ (γ, signs included) and formal → Brune → formal 2.7 × 10⁻¹⁵ (E, γ) with the iteration tolerance at 10⁻¹³ MeV,
while the engine's 10⁻⁶ MeV leaves 3.8 × 10⁻⁹ relative in E; in six conventions T T^T = 1 to 3.3 × 10⁻¹⁶ and the
formal parameters from the project's widths equal the internal ones to 3.6 × 10⁻¹⁶.

`physics_level_parameters_from_cnuc` (`tests/upgrade/level_parameters_from_cnuc_test.cpp`) on T1, T2, T3: (a) the
file's widths through `level_parameters` against `CNuc` after `Initialize` in the Brune, Park and formal bases, Park →
formal against the formal `CNuc`, formal → Brune against `CNuc::TransformOut`; (b) round trips; (c) every block in six
conventions, written as widths in the project's columns and converted back; (d) T3's 2⁺ block in jj with one p + ⁷Li
channel removed. An error in any conversion step, a channel attached to the wrong column, a wrong B_c or a wrong
level sign fails (a): since step 1b the engine and the test use one eigensolver and one sign rule, so the formal
comparisons allow no sign. Results (2026-11-07, GCC 11.5, relative differences; "identical" counts bit-equal numbers):

| comparison | T1 | T2 | T3 |
|---|---|---|---|
| (a) Brune γ̃ vs `CNuc` (Brune) | 21/21 identical | 20/20 identical | 16/16 identical |
| (a) B_c vs `CalcBoundaryConditions` | 10/10 identical | 8/8 identical | 12/12 identical |
| (a) Park γ^P and J vs `CNuc` (`--use-park`) | identical | identical | identical |
| (a) formal E vs `CNuc` (`--no-brune`) | 7/7 identical | 5/5 identical | 7/7 identical |
| (a) formal γ vs `CNuc`, signs included | 21/21 identical | 20/20 identical | 16/16 identical |
| (a) Park → formal, particle γ | 1.9 × 10⁻¹⁶ | 3.7 × 10⁻¹⁶ | 1.6 × 10⁻¹⁶ |
| (a) formal → Brune vs `TransformOut`: E, γ | identical | identical | identical |
| (b) observed → Brune → observed | 3.1 × 10⁻¹⁶ | 3.6 × 10⁻¹⁶ | 2.4 × 10⁻¹⁶ |
| (b) Brune → Park → Brune | 2.1 × 10⁻¹⁶ | 1.9 × 10⁻¹⁶ | 2.1 × 10⁻¹⁶ |
| (b) formal → Brune → formal (10⁻¹³ MeV): E, γ | 2.4 × 10⁻¹⁵, 5.3 × 10⁻¹⁵ | 2.3 × 10⁻¹⁴, 1.9 × 10⁻¹⁴ | 3.1 × 10⁻¹⁵, 5.0 × 10⁻¹⁶ |
| (c) six conventions: formal γ vs internal columns | ≤ 1.1 × 10⁻¹⁵ | ≤ 1.5 × 10⁻¹⁵ | ≤ 2.7 × 10⁻¹⁶ (jj: 2⁺ block reported, (d)) |

The formal → Brune iteration takes the same number of steps as the engine's (at most 6, 13, 7).

Mutations (each applied alone, all three test programs rebuilt): B_c dropped from N — unit and T1–T3 tests fail; the
6j arguments j_x ↔ j_y — all three fail; the Brune factor N/D inverted — unit and T1–T3 fail; B_c from the last level
— unit and T1–T3 fail; eigenvalues assigned by position instead of rank — the unit test fails (T1–T3 list their levels
in energy order, so they do not see it); the heavy-first jj phase dropped — the m-scheme test fails (a consistent sign
error leaves the round trips and the invariance test unchanged, which is why the m-scheme test exists).

Convention invariance is tested at the parameter level: through `level_parameters` the formal parameters, and so every
observable computed from them, are the same to rounding. It is not tested through the engine, which reads no
convention record and no jj project; that needs the project reader (`io/`) to carry the record.

## 8. Not done

- The bound-state normalization N_f^½ of external-capture final states (`SqrtNFFactor`) is not computed here; it
  belongs to the external-capture module with γ^ext.
- No `.azr` keyword states a convention; the record exists only in `physics/`.
- The binding exposes none of these functions yet (`RULES.md` §1.7): with `channel_functions` available, the
  arrays-in, arrays-out signatures can be bound directly.

Code: `physics/include/azure2/physics/{Convention,LevelParameters,AngularCoupling,SymmetricEigen}.h`,
`physics/src/{LevelParameters,Recoupling,AngularCoupling,SymmetricEigen}.cpp`; adapter `include/CNucLevelParameters.h`,
`src/CNucLevelParameters.cpp`.

## 9. The formal → Brune iteration inside an evaluation (design, 2026-11-08)

**What.** `brune_from_formal` iterates, per level, B_c ← S_c(E) and the rediagonalization of
E − Σ_c (S_c − B_c) γγᵀ until the level's energy moves by at most a tolerance; a level that does not converge in
1000 iterations, or whose eigenproblem fails, keeps its formal values (unit row of T). That is the engine's report
iteration (`TransformOut`, 10⁻⁶ MeV) and stays so. Since step 2 the same function runs inside every formal-basis
evaluation (`ECFinalStatesFromFormal`, `FormalCaptureMasksFromParams`, and `SetUpECFinalBlocks` at setup), at
10⁻¹² MeV. There the iteration can stall at its rounding floor: on `beam_profile_kernel` the 1⁻ level with Brune
energy 20.0 MeV ends at steps of 7.7 × 10⁻¹² MeV (3.9 × 10⁻¹³ relative, about 1700 ε: the rounding of E is amplified
by the slope of the shift function) after a 10⁻⁶ change of one width, its formal values are kept, the capture mask
changes and χ² jumps 40.44 → 42.93 (`external_capture.md` §8). Keeping the formal values is right for a report that
warns; inside a fit it is a silent, discontinuous substitute.

**Rule** (options of `FormalToBruneOptions`; the defaults keep the report iteration unchanged):
1. `relative_tolerance`: converged when the step is at most max(`energy_tolerance_MeV`, `relative_tolerance`·|E|).
   The evaluation callers use 10⁻¹² MeV and 10⁻¹² relative (about 4500 ε): the stall above converges, and the
   formal ↔ Brune capture agreement stays at 10⁻¹¹.
2. `evaluation = true`: a level that has not met the tolerance after the iterations, but whose last step is at most
   `rounding_bound_MeV` (10⁻⁹ MeV), takes its last iterate (status `RoundingLimited`; continuous in the parameters,
   accurate to that step). A level whose eigenproblem fails, or whose last step is larger, gets NaN energy and
   amplitudes (status as now, row of T NaN): the model is undefined there, so the observables and χ² are not finite —
   MCMC rejects the point (its log-probability is tested with `isfinite`), Minuit and pyazr see a failed evaluation.
   No formal values are substituted. The engine warns once per process on stderr with the level and the step.
   At setup (`SetUpECFinalBlocks`) a failure stops the run (`GSLException`), as `TransformIn` does.

**Where.** `physics/` (`brune_from_formal`, `FormalToBruneOptions`, `LevelSet::Status`); the three callers in
`src/CNucExternalCapture.cpp` set the options through one function. Nothing else changes.

**Tests.** `level_parameters_unit_test`: with `max_iterations` = 1, `evaluation` and a rounding bound above the first
step, the last iterate is returned (not the formal values) with status `RoundingLimited`; with a bound below it, NaN;
without `evaluation`, the formal values as before; a relative tolerance of 10⁻³ converges in fewer iterations to an
energy within 10⁻³|E| of the tight one. A pyazr test on `beam_profile_kernel` (formal basis): χ² at the start and
±10⁻⁶ relative on `width_9_2` is smooth (the second difference is below 10⁻³ of the first). Expected outputs: menu 1
of the 179 projects unchanged in Brune and Park; formal at most in the last printed digits, every change classified.

**Verified (2026-11-09; build of 92d1525e, GCC 11.5 Release, pyazr on; local 24-core node).** No warnings.
`ctest -LE slow` 110/111, `physics_level_parameters_unit` and the slow `pyazr_formal_brune_iteration` pass (41.7 s);
`run_tests.sh` 133/133 (1185 s, 11 threads); 1 vs 8 threads byte-identical on T1–T3 and 13N in the three bases. The
179 projects × 3 bases against a18fe89e: Brune 689/689 and Park 689/689 identical; formal 657/662, the five
changes being `beam_profile_kernel` `AZUREOut_aa=1_R=4.out` (two tokens, 1.7 × 10⁻¹¹ relative) and four `run.log`
lines that count Coulomb-table lookups (fewer iterations under the relative tolerance, e.g. 318934 → 300814). No
evaluation of the 179 projects warned or gave NaN.

**Blocker: `pyazr_formal_capture_mask` fails (check (c), formal-basis χ² gradient = NaN).** Its model (13N plus a
1/2⁺ level at 4.5 MeV with Γ_p = 20 keV, `formal_capture_mask_test.py`) sits where the fixed-point iteration of
`brune_from_formal` does not converge for block 1, level 0: at the test's start point it converges, but at ±10⁻⁶ MeV
on `energy_3` (2.3682 MeV) it alternates between E = −0.973 and 4.462 MeV (last step 5.435 MeV after 1000
iterations, a 2-cycle). The rule above then gives NaN, as designed, and the central differences of check (c) are
NaN. On a18fe89e the same four evaluations did not converge either (instrumented build: four "not converged" lines,
same 2-cycle) and silently kept the formal values, so check (c) passed there by comparing two silently substituted
models. This is the case the rule was written for, but it shows a gap in the iteration, not in the rule: a 2-cycle
is not a rounding stall. Options for the PI: (i) a safeguarded iteration in `brune_from_formal` (damping, or a
bracketed root of E − e(E) per level) so that such levels converge, with its own design note and the unit test as
the case; (ii) move the test model away from the 2-cycle and record the 2-cycle as a known NaN region. Not merged
until decided.

**Decision (PI, 2026-11-09): option (i).** The iteration must converge where the root exists; see §10.

## 10. A bracketed root where the fixed-point iteration does not converge (design, 2026-11-09)

**Why the root exists and is unique.** For level λ of rank r (ascending formal energy), the Brune energy is the
solution of g(x) = x − e_r(x) = 0, with e_r(x) the r-th eigenvalue (ascending) of
A(x) = diag(E_formal) − Σ_c (S_c(x) − B_c) γ_c γ_cᵀ (particle channels; γ_c the column of formal amplitudes). The
iteration of §4 is the fixed-point map x ← e_r(x): it re-bases B_c at S_c(x) and rediagonalizes, which leaves the
spectrum of A(x) unchanged. Since dA/dx = −Σ_c S_c′(x) γ_c γ_cᵀ is negative semidefinite whenever every S_c′ ≥ 0 (the
shift function rises with energy, open and closed channels alike), each sorted eigenvalue e_r(x) is non-increasing
(Weyl's monotonicity), so g′(x) ≥ 1: one root per level, always. The fixed-point map converges only where
|e_r′| < 1 at the root. Where e_r′ < −1 the root repels: the iterates jump to either side and end in a 2-cycle
(x_a, x_b) with e_r(x_a) = x_b, e_r(x_b) = x_a — and then g(x_a) = x_a − x_b and g(x_b) = x_b − x_a have opposite
signs, so the 2-cycle brackets the root.

**Measured on the blocker** (`pyazr_formal_capture_mask`, block 1, level 0; debug build of 46efe93e, scan of g on
−3 … 7 MeV in steps of 0.05 MeV with the formal parameters of the four failing evaluations): g increases
everywhere, slope 1.003 … 3.96; one root, at 2.36820 MeV (the 2.368 MeV resonance), e_r′ = −2.88 there; bisection
between the two iterates of the 2-cycle (−0.973, 4.462 MeV) finds it in each of the four evaluations. The centre point
of check (c) "converged" only because the iteration started within the tolerance of the root. The report iteration
of the 179 projects × 3 bases (step 4 build) leaves three levels of `screening_6li` (formal basis; background
levels at 43.9, 122.8, 141.0 MeV, residuals 41–118 MeV) not transformed, in its three files.

**Rule.** In `brune_from_formal`, a level for which the fixed-point iteration ends without meeting the tolerance and
with a last step above `rounding_bound_MeV` (where the status would be `NotConverged`) is solved for the root of g
instead, in the report and the evaluation iteration alike:
1. Bracket: the last two iterates, if g is finite at both and changes sign between them (a 2-cycle); otherwise —
   an iteration that diverges overflows instead of cycling — expand from the level's formal energy towards the sign
   change (g < 0: upwards), starting with a width of 10⁻³ MeV, doubled up to 64 times.
2. Bisection on the bracket until the midpoint equals an end (full double precision) or g = 0; at most about 60
   eigenproblems of the block's size.
3. The level takes the r-th eigenvalue and eigenvector of A at the root: E, the Brune amplitudes Vᵀγ and the row of
   T, the eigenvector oriented with its component on the level's own formal level positive (the orientation rule
   of §4 and `formal_from_brune`).
4. Status `Bracketed` (new), iterations = the fixed-point iterations plus the bisection steps, `residual_MeV` = |g|
   at the root, `boundary_shift` = S_c(root) − B_c. A warning names the level and says it was solved by bracketing
   (the evaluation iteration prints the first one per process, as now; the report prints a NOTE in place of the
   WARNING, without "the input energy and widths are kept").
5. An iteration that diverges overflows: its matrix is no longer finite and the eigenproblem fails. That is the
   same case, solved as in 1 from the formal energy. No bracket (g not finite, or no sign change after the
   expansion), or an eigenproblem that fails on a finite matrix: as now (NaN in an evaluation, formal values with the
   WARNING in the report, a stop at setup).
Levels that converge are untouched: the fallback only runs where the old result was the formal values or NaN, so
all other outputs stay bit-identical. Cost: only for those levels, the 1000 fixed-point iterations plus the bisection.

`FormalToBruneOptions::bracket` (default true) switches the rule off, for the unit tests of the fixed-point branches
(they stop the iteration with `max_iterations` = 1, where a bracket would now be found).

**Where.** `physics/src/LevelParameters.cpp` (`brune_from_formal`), `LevelSet::Status::Bracketed`; the report message
in `src/CNuc.cpp`; `ThrowIfFailed` and `WarnOnce` (`src/CNucExternalCapture.cpp`) unchanged (Bracketed is not a
failure).

**Tests.** `level_parameters_unit_test`: a one-level block made by `formal_from_brune` from a Brune level with
γ²S′(E) > 1 (e′ < −1: the fixed point repels): `brune_from_formal` returns the Brune input to 10⁻¹⁰ relative with
status `Bracketed`, in report and evaluation mode (a linear S: the iterates diverge and overflow, the bracket comes
from the formal energy); the three-level case with `max_iterations` = 1 gives every level Converged or Bracketed at
the tight result (10⁻¹⁰); the existing fixed-point tests with `bracket` = false keep their expectations; channel
functions that are NaN away from the formal energy: no bracket, NaN in an evaluation and the formal values in the
report. `pyazr_formal_capture_mask` (the blocker) passes unchanged. Expected outputs:
179 × 3 identical except the formal basis of `screening_6li` (its three levels now transformed: `parameters.out`
and `run.log`), classified.

**Verified (2026-11-09; build of 5d8dc1e1, GCC 11.5 Release, pyazr on).** No warnings. `ctest -LE slow` 111/111
(the blocker passes: check (c) formal 1.4 × 10⁻⁷, Brune 1.1 × 10⁻⁷), the slow `pyazr_formal_brune_iteration` passes;
`physics_level_parameters_unit` with §6b; `run_tests.sh` 133/133 (1182 s, 11 threads); 1 vs 8 threads
byte-identical (T1–T3, 13N, three bases). 179 × 3 against 46efe93e with the same format-1 table: Brune 689/689 and
Park 689/689 identical, formal 656/662: the six files of `screening_6li` (`parameters.out`, `run.log` of bare, fit,
scr), where the three levels the old iteration left untransformed (formal values 122.75, 43.89, 141.04 MeV printed
under a WARNING) are now bracketed: 30.0000, 11.5670, 30.0000 MeV, every width equal digit for digit to the
Brune-basis run of the same file (the input), and a NOTE in place of the WARNING. Classified as a fix (the round
trip formal → Brune now returns the input).
