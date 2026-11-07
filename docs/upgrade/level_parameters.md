# Level parameters: bases and conventions (`level_parameters`)

The collision matrix of a J^π block is computed from level energies E_λ (MeV, excitation energy E_x of the compound
nucleus) and reduced width amplitudes γ_λc (MeV^½). A project states these numbers in its own form: observed partial
widths Γ (eV) and asymptotic normalization coefficients C (fm^−½) in Brune's or Park's parametrization, or formal
amplitudes, coupled in some angular-momentum scheme with some phase convention. `level_parameters`
(`physics/include/azure2/physics/LevelParameters.h`) converts them to the one internal set: formal R-matrix parameters
with constant boundary values B_c, Lane–Thomas phases, channel-spin coupling s = j₁ + j₂ (light particle first),
J = s + l (`docs/source/theory/conventions.rst`). It is the only place in the code where a conversion between bases or
conventions happens (`RULES.md` §1.6). This page gives the equations of each conversion, their validity, the
tolerances, and what each test would detect. The engine calls it: `CNuc::CalcBoundaryConditions`, `CNuc::TransformIn`
and the formal → Brune iteration of `CNuc::TransformOut` delegate to these functions (`cnuc_split.md`); the
convention conversions (section 5) are not used by the engine yet, which reads no convention record. The `CNuc.cpp`
line numbers below are those of e54e6fcc, the code these functions were checked against before the move.

## 1. Symbols

| symbol | unit | meaning |
|---|---|---|
| E_λ, Ẽ_λ | MeV | formal and observed (Brune) level energy, E_x of the compound nucleus |
| E_c = E − (e₀ + e₁) − S_α | MeV | channel energy, c.m.; S_α separation energy, e₀, e₁ excitation energies of the pair |
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
(`src/CNucLevelParameters.cpp`) supplies them from `ChannelFunc` and `WhitFunc`, with E_c in the order of
`CNuc::TransformIn`.

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

    B_c = S_c(Ẽ_ref),     Ẽ_ref = energy of the first level listed (engine; `BoundaryRule::FirstListedLevel`)   (3.1)

or of the lowest level (`BoundaryRule::LowestLevel`, the proposal of `conventions.rst`, off by default). Observables do
not depend on B_c; the formal parameters do, so with the engine's rule they depend on the order of the levels in the
file. A block without levels has B_c = 0. Non-particle channels copy channel 0's value, as in
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

The overall sign of each formal level's amplitude vector is the eigenvector's, which the eigenproblem leaves free.
`SymmetricEigen.h` fixes each eigenvector's largest component positive; GSL makes another choice. The eigensolver is an
option of `formal_from_brune` and `brune_from_formal` (`GeneralizedEigenSolver`, `FormalToBruneOptions::solver`); the
engine passes GSL's, so its formal parameters keep GSL's signs (`cnuc_split.md` section 4). Observables
are even in that sign (γ^T A γ is invariant). The rule "largest |γ_λc| of a level positive" (the proposal of
`conventions.rst`) is `FormalSignRule::LargestPositive`, off by default.

The inverse (`brune_from_formal`) is the engine's fixed-point iteration (`CNuc::TransformOut`, `CNuc.cpp:2563-2745`):
for each level, B_c is moved to S_c at the current energy of that level, the matrix
E_μ δ_μμ' − Σ_c ΔB_c γ_μc γ_μ'c is diagonalized, and the iteration stops when the level's energy moves by less than
the tolerance, 10⁻⁶ MeV in the engine (an option here). A level that has not converged after 1000 iterations keeps its
formal values, with a warning.

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

## 6. The engine's formal basis and external capture: an ordering effect

`CNuc::TransformIn` converts the J^π groups one after the other and replaces each group's levels by their formal
parameters before the next group is treated (`CNuc.cpp:1464-1607`). The external-capture amplitude γ^ext of a photon
channel (Eq. 3.2) is computed in that loop from the final state's E and γ as they stand: a final state in an earlier
group enters with its formal parameters, one in the same or a later group with its Brune ones. In the Brune and Park
bases the loop reads only the final state's particle amplitudes, which it does not change, so the effect is confined to
`--no-brune` with external capture and a final state in a group of several levels listed before the capturing group.
None of T1–T3 has one: T1's 14 non-zero external amplitudes all go to the ground state, alone in its group, whose
formal and Brune parameters are identical (B_c = S_c(E_1) for a single level), and they are the same in either order
(`physics_level_parameters_from_cnuc` prints the count). The test computes γ^ext in the engine's order, so it reproduces
the engine either way. Whether the engine's order is intended is a decision (HANDOFF); it is not changed here.

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
Brune, formal → Brune → formal. (3) Convention invariance on six conventions; (4) reported gaps; (5) the proposals.
Identity (1) to 4.4 × 10⁻¹⁵ relative (3.9 × 10⁻¹⁵ with B_c from the lowest
level); observed → Brune → observed 1.5 × 10⁻¹⁶ and observed → Park → observed 1.9 × 10⁻¹⁶ element-wise; Brune → Park
→ Brune 1.5 × 10⁻¹⁷; formal → Brune → formal 3.1 × 10⁻¹⁴ (E) and 3.4 × 10⁻¹⁴ (γ) with the iteration tolerance at
10⁻¹³ MeV, while the engine's 10⁻⁶ MeV leaves 2.4 × 10⁻⁸ relative in E; in six conventions T T^T = 1 to 3.3 × 10⁻¹⁶
and the formal parameters from the project's widths equal the internal ones to 4.5 × 10⁻¹⁶.

`physics_level_parameters_from_cnuc` (`tests/upgrade/level_parameters_from_cnuc_test.cpp`) on T1, T2, T3: (a) the
file's widths through `level_parameters` against `CNuc` after `Initialize` in the Brune, Park and formal bases, Park →
formal against the formal `CNuc`, formal → Brune against `CNuc::TransformOut`; (b) round trips; (c) every block in six
conventions, written as widths in the project's columns and converted back; (d) T3's 2⁺ block in jj with one p + ⁷Li
channel removed. An error in any conversion step, a channel attached to the wrong column or a wrong B_c fails (a).
Results (GCC 11.5, relative differences; "identical" counts bit-equal numbers):

| comparison | T1 | T2 | T3 |
|---|---|---|---|
| (a) Brune γ̃ vs `CNuc` (Brune) | 21/21 identical | 20/20 identical | 16/16 identical |
| (a) B_c vs `CalcBoundaryConditions` | 10/10 identical | 8/8 identical | 12/12 identical |
| (a) Park γ^P and J vs `CNuc` (`--use-park`) | identical | identical | identical |
| (a) formal E vs `CNuc` (`--no-brune`) | 3.6 × 10⁻¹⁶ | 6.6 × 10⁻¹⁶ | 2.8 × 10⁻¹⁶ |
| (a) formal γ vs `CNuc`, up to level sign (of max \|γ\|) | 9.3 × 10⁻¹⁷ | 1.0 × 10⁻¹⁴ | 1.1 × 10⁻¹⁶ |
| (a) Park → formal, particle γ | 1.9 × 10⁻¹⁶ | 1.0 × 10⁻¹⁴ | 1.6 × 10⁻¹⁶ |
| (a) formal → Brune vs `TransformOut`: E, γ | 1.6 × 10⁻¹⁴, 3.8 × 10⁻¹⁴ | 5.5 × 10⁻¹⁴, 6.9 × 10⁻¹⁴ | identical |
| (b) observed → Brune → observed | 3.1 × 10⁻¹⁶ | 3.6 × 10⁻¹⁶ | 2.4 × 10⁻¹⁶ |
| (b) Brune → Park → Brune | 2.1 × 10⁻¹⁶ | 1.9 × 10⁻¹⁶ | 2.1 × 10⁻¹⁶ |
| (b) formal → Brune → formal (10⁻¹³ MeV): E, γ | 3.0 × 10⁻¹⁵, 5.6 × 10⁻¹⁵ | 6.4 × 10⁻¹⁴, 4.3 × 10⁻¹⁴ | 3.0 × 10⁻¹⁵, 3.9 × 10⁻¹⁶ |
| (c) six conventions: formal γ vs internal columns | ≤ 1.1 × 10⁻¹⁵ | ≤ 1.5 × 10⁻¹⁵ | ≤ 2.7 × 10⁻¹⁶ (jj: 2⁺ block reported, (d)) |

The formal level vectors whose overall sign differs from GSL's: 1 of 7 (T1), 2 of 5 (T2), 0 of 7 (T3). The
formal → Brune iteration takes the same number of steps as the engine's (at most 6, 13, 7). B_c is bit-identical
although `CalcBoundaryConditions` forms E_c as E − (S + e) and the adapter as E − e − S: on T1–T3 the two orders round
alike. The formal basis, which goes through a different eigensolver, agrees to ≤ 1.0 × 10⁻¹⁴ (the bit-for-bit default
holds for Brune, Park and B_c; ≤ 10⁻¹² relative for formal parameters, as promised).

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
