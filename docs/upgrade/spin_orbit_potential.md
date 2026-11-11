# A spin–orbit term in the hybrid potentials

The hybrid model replaces, for one particle pair, the Coulomb functions F_l, G_l of the external region by the
solutions of a radial equation that contains a nuclear potential besides the Coulomb one; the penetrability P, the
shift function S, the phase φ and the external-capture waves of that pair are then computed from those solutions. The
potential is central today, so the channel functions depend on (α, l, E_c) only. This note states what the model is
(§1), works out what a spin–orbit term does to channels that AZURE2 labels by channel spin (§2), and proposes the
design (§3), the tests (§4), the order of the work (§5) and the decisions the PI has to take (§6). It is the design
note that `cnuc_split.md` step 3b (iii) asks for. No engine code is changed by it. Code references are to commit
070756a1 (`feature/upgrade-4-data`). The numerical checks named here were run with the scripts in
`work/so-note/` (§4.1); everything marked *unchecked* is believed, not verified.

## 1. Today

### 1.1 The radial equation and what is matched

For a pair α with reduced mass μ (MeV, `redmass() · uconv`) and charges Z₁, Z₂, `CoulFunc::computeHybrid`
(`src/CoulFunc.cpp:330-481`) solves, at the channel energy E (MeV, c.m.) and orbital angular momentum l,

  u″(r) = (2μ/ħ²) [V(r) − E] u(r),  V(r) = V_C(r) + V_N(r) + ħ² l(l+1) / (2μ r²),        (1)

with

- V_C the potential of a uniformly charged sphere of radius R_C = 3.5 fm inside R_C and Z₁Z₂e²/r outside, with
  e² = 1.4399764 MeV fm (`include/NuclearPotential.h:189-201`; R_C is the constructor's default,
  `NuclearPotential.h:165-170`, which `CoulFunc.cpp:423-427` does not override);
- V_N either a Woods–Saxon well −V₀ / (1 + exp((r − R)/a)) (`NuclearPotential.h:33-66`) or a Gaussian
  −V₀ exp(−(r/r₀)²) (`NuclearPotential.h:68-99`), V₀ > 0 (MeV) attractive, R, a, r₀ in fm;
- the centrifugal term with ħc = 197.3269804 MeV fm (`NuclearPotential.h:207-211`).

Equation (1) is integrated **inward** by Numerov's method from the matching radius r_m = 15 fm
(`CoulFunc.cpp:25`), where the two outermost nodes are set to the Coulomb H⁺_l = G_l + iF_l at η = Z₁Z₂ α √(μ/2E)
and ρ = kr (`CoulFunc.cpp:363-415`). The step is the one closest to 0.0275 fm that puts the channel radius a_c on a
node (`CoulFunc.cpp:402-405`). The solution H(r) of (1) with H(r) = H⁺_l(kr) for r ≥ r_m is read at a_c (node value),
its ρ-derivative by a five-point difference (`CoulFunc.cpp:451-452`), and stored as G = Re H, F = Im H
(`CoulFunc.cpp:455-458`). Because (1) is real, F and G are each solutions, their Wronskian is that of the Coulomb
functions, and every formula built on F, G (P = ρ/(F² + G²), S = ρ (FF′ + GG′)/(F² + G²), e^{−iφ} = (G − iF)/|H|,
`physics/include/azure2/physics/ChannelFunctions.h:54-58`) keeps its form:

  L(a_c) ≡ ρ H′(a_c)/H(a_c) = S + iP,  φ = arg H(a_c).        (2)

What is matched at a_c is therefore the logarithmic derivative of the outgoing solution *of the potential* (2); the
interior is the R-matrix's as before. The potential acts only between a_c and r_m: the values read at a_c depend on V
for r ≥ a_c − 2 dr only (the Numerov loop runs on to r = 0, `CoulFunc.cpp:437`, but nothing below node n_a − 2 is
read). The Coulomb phase factor e^{i(σ_l − σ₀)} stays the point-Coulomb one (`src/ChannelTable.cpp:31-37`), which is
right because H → H⁺_l beyond r_m. A requested radius beyond r_m returns the plain Coulomb functions
(`CoulFunc.cpp:343-360`), so the external-capture integrand, which runs from a_c to infinity, uses the distorted
waves up to 15 fm and the Coulomb ones beyond.

No paper is cited for the model in the code, its commit (03d9152a, 2025-10-28) or the user documentation. The model is
the "potential in the external region" form of R-matrix matching; I have not found a reference in the repository and
cite none.

### 1.2 Inputs and where they live

- The `<potential>` block of format 1/2 (`src/Config.cpp:302-423`): `useHybridPotential`, `potentialType` (0
  Woods–Saxon, 1 Gaussian), `V0`, `R`, `a`, `r0`; keys before the first `pair=<key>` are the default, each `pair=`
  section starts from the default (`Config.cpp:365-378`). `useHybridMethod` is the OR over pairs (`Config.cpp:421`).
- `NuclearPotentialManager` (`include/NuclearPotentialManager.h:52-223`, `src/NuclearPotentialManager.cpp`): a
  process singleton (`NuclearPotentialManager.cpp:5-8`) holding a default and per-pair settings
  (`NuclearPotentialManager.h:17-30`) and a "tag" counter that changes with each setting (`tagFor`,
  `NuclearPotentialManager.cpp:190-195`). pyazr sessions swap the whole state in and out around each call
  (`exportState`/`importState`, `NuclearPotentialManager.h:193-205`; `api/src/pybind_api.cpp:63-92` (`ConfigScope`), `286-320`).
- `CoulFunc`'s constructor reads the switch, the tag and the potential of its pair (`CoulFunc.cpp:33-36`, `54`) and
  turns the Coulomb table off for that pair (`CoulFunc.cpp:58`). r_m and the step have setters
  (`CoulFunc.cpp:502-528`) that nothing calls: they are constants in practice.
- The tag enters every memo and table keyed on the Coulomb identity of a pair: the radial scratch
  (`CoulFunc.cpp:162-170`, `215`), the channel table's pair signature (`ChannelTable.cpp:51`), the shift-derivative
  memo (`src/CNuc.cpp:1592-1601`, `1711-1713`), the external-capture node table (`src/ECTable.cpp:82`) and the THM
  vertex memo (`src/ThmVertexBoundary.cpp:63-85`). **None of these keys has a j**: they are (pair, l) keys.

### 1.3 Where the hybrid functions enter

| quantity | formula | code path |
|---|---|---|
| P, S, e^{−iφ} of an open channel at a data point | (2) at a_c | `ChannelTable::Compute` (`ChannelTable.cpp:18-40`) → `EPoint::CalcEDependentValues` (`src/EPoint.cpp:1281-1316`): L⁰_c = S − B_c + iP per (J-group, channel) slot |
| level matrix | A⁻¹_{λμ} = (E_λ − E)δ_{λμ} − Σ_c γ_λc L⁰_c γ_μc | diagonal in c: `src/AMatrixFunc.cpp:277-285` |
| collision matrix | U_{cc′} = Ω_cΩ_{c′}[δ_{cc′} + 2i P_c^½ P_{c′}^½ (γᵀAγ)_{cc′}] | `AMatrixFunc.cpp:554-570` (conventions.rst Eq. conv-u) |
| S, P, S′ at a level energy (Brune, Park, ANC) | (2), and its energy derivative | `ChannelFunc::Shift/Penetrability/ShiftDerivative` (`src/ChannelFunc.cpp:39-68`) via `ChannelFunctionsFromCNuc` (`src/CNucLevelParameters.cpp:19-48`) |
| external capture, entrance wave | F, G of (1) at the quadrature radii, √P and e^{−iφ} at a_c | `src/ECIntegral.cpp:45-56`, `357-370`; `src/ECTable.cpp:96-102`; `EPoint.cpp:1762-1780` |
| THM vertex | point-Coulomb O_l only, the potential is not applied | `include/ThmFunc.h:88-89`; S at the level energy through `ChannelFunc` (`ThmVertexBoundary.cpp:84`) |

### 1.4 Limits of today's model

1. **Below threshold the potential is absent.** For E_c < 0, at E_c = 0 and in the closed band η > 100, S comes
   from the Whittaker function or its continuation (`ChannelFunc.cpp:39-43`; `ShftFunc`), which knows no potential;
   so do the ANC ↔ γ conversion (`CNucLevelParameters.cpp:35-40`) and the bound-state factor of external capture
   (`ECIntegral.cpp:50-51`). S of a hybrid pair is therefore discontinuous at threshold. Size, estimated with an
   independent integration (`work/so-note/threshold.py`, not the engine) for the settings of `tests/hybrid_potential`
   (p + ¹²C, V₀ = 20 MeV, R = 3.6 fm, a = 0.6 fm, a_c = 3.6 fm): at E_c = 0.05 MeV the distorted S exceeds the
   Coulomb S by 1.07 (l = 0), 0.92 (l = 1), 0.75 (l = 2), and P by a factor 1.5. A bound level of that pair (the ¹³N
   ground state of the project) is transformed with the Coulomb S while the scattering uses the distorted one. The user
   documentation states the bound-state half of this (`docs/source/user_guide/pyazr.rst:841-855`), not the threshold
   step of S.
2. **The Coulomb part is not the engine's.** e² = 1.4399764 MeV fm (`NuclearPotential.h:191`) differs from
   α ħc = 1.4399645 MeV fm (`physics/include/azure2/physics/Constants.h:8,10`) by 8.3 × 10⁻⁶ relative, while the
   boundary values at r_m use α. With V_N = 0 the hybrid P then differs from the exact Coulomb P by 1.7 × 10⁻⁵ (p +
   ¹²C, l = 1, 0.5 MeV) and S by 3.9 × 10⁻⁶; with e² = α ħc the same Numerov integration agrees to 6 × 10⁻⁹ and
   2 × 10⁻⁸ (`work/so-note/radial.py`, part B). The "about 10⁻⁵, the accuracy of the integration itself"
   (`pyazr.rst:830-833`, `tests/pyazr/hybrid_derivative_test.py:9-11`) is this constant, not the integration.
   R_C = 3.5 fm is fixed for every pair.
3. **Parameters are not fit parameters.** A change rebuilds the model and re-derives the amplitudes
   (`tests/hybrid_potential/README.md:77-85`).
4. **Not relativistic** (refused, `CNuc.cpp:202-207`), **not in the THM vertex** (above).
5. **Imprecise statements in the code and docs**: "inside the channel radius" (`CNuc.cpp:203`,
   `docs/source/user_guide/particle_pairs.rst:439`) — the potential acts outside it; "outward"
   (`include/CoulFunc.h:100`, `tests/hybrid_potential/tolerance:4`) — the integration is inward.

### 1.5 Tests that use the hybrid model

| test | system | what it covers |
|---|---|---|
| `tests/hybrid_potential` | p + ¹²C, j₁ = ½, j₂ = 0, a_c = 3.6 fm, three Meyer elastic segments (dσ/dΩ) with target integration | χ² pinned at 2 × 10⁻² (`tolerance`); per-pair precedence (`check_per_pair.py`) |
| `tests/upgrade/T3` | α + α (both spin 0), a_c = 5.5 fm, WS 10 / 3.0 / 0.6 | α + α elastic in the T3 golden outputs (`tests/upgrade/T3/README.md:13,22-25`) |
| `tests/pyazr/hybrid_derivative_test.py` | p + ¹²C, V₀ = 10⁻⁹ MeV | P, S against the plain Coulomb functions, 10⁻⁴ |
| `tests/pyazr/potential_sessions_test.py` | p + ¹²C | one session's potential does not leak into another |
| `tests/relativistic/check.sh` (g) | — | the refusal with `relativistic=1` |

No test has a capture reaction with a hybrid pair, and none has a hybrid pair whose two spins are both non-zero. For a
spin–orbit term, `hybrid_potential` (spin ½ on spin 0) is the natural case; α + α (T3) has s = 0 and no spin–orbit
coupling at all.

## 2. The physics of the extension

### 2.1 The term and its convention

The potential becomes

  V(r) = V_C(r) + V_N(r) + ħ² l(l+1)/(2μr²) + v_so(r) ⟨2 l·i⟩,  v_so(r) = V_so (ħ/m_π c)² (1/r) df_so/dr,   (3)

with f_so(r) = 1/(1 + exp((r − R_so)/a_so)) (a Woods–Saxon form factor; f_so → 1 inside), V_so (MeV) > 0, and
(ħ/m_π c)² = 2.000 fm² exactly. **i** is the spin of the particle the term acts on (§2.2–2.3). Since df_so/dr < 0,
V_so > 0 lowers j = l + i (the nuclear sign). For spin ½, 2 l·i = l·σ, with eigenvalue l for j = l + ½ and −(l + 1)
for j = l − ½.

Sources consulted: the FRESCO namelist manual (www.fresco.org.uk, input3.4, "&POT", node 8) gives the projectile
(TYPE 3) and target (TYPE 4) spin–orbit form factors, e.g. the Woods–Saxon one −CONLS · V · E/[(1 + E)² a r] with
E = exp(−(r − R)/a) and CONLS = (ħ/m_π c)² = 2.000, "multiplied by j(j+1) − l(l+1) − s(s+1) = 2 l.s, where s =
spin of projectile or target". With f = E/(1 + E), df/dr = −E/[a(1 + E)²], so FRESCO's term is exactly (3). J.
Raynal, arXiv:nucl-th/0312038 (2003), Eq. (2), writes the optical-model spin–orbit term as (2/r)(dV/dr) l·s, the
same factor 2. The factor 2 is the trap: a code or parameter table written with l·s and no 2 has half our V_so for
the same physics. *Unchecked*: the conventions of ECIS, of the global sets (Becchetti–Greenlees, Koning–Delaroche,
RIPL) and whether they use 2.000 fm² or (ħc/m_π±c²)² = 1.9989 fm² — I did not open these sources.

### 2.2 One spinning particle: diagonal

AZURE2's channels are c = (α, l, s) with **s** = **j₁** + **j₂** coupled light particle first and **J** = **s** +
**l** (`docs/source/theory/conventions.rst`, Eq. conv-coupling). If only one particle has spin (i₂ = 0, or i₁ = 0
with the term on the heavy particle), s is that spin, j ≡ J, and the term is diagonal:

  ⟨2 l·s⟩ = J(J+1) − l(l+1) − s(s+1).        (4)

Each (α, l, J) then has its own scalar equation (1) + v_so ⟨2 l·s⟩, its own H_{lJ}, and P, S, φ of (2) depend on
(α, l, J). Nothing else changes in the R-matrix: L⁰ stays diagonal, U keeps its form.

### 2.3 Two spinning particles: the channel spins of one (α, l, J) are coupled

A projectile spin–orbit term l·**i₁** with i₂ ≠ 0 conserves l and J but not s. It is diagonal in the jj basis
|((i₁ l) j, i₂) J⟩, with eigenvalue λ_j = [j(j+1) − l(l+1) − i₁(i₁+1)]/2. The recoupling to the channel-spin basis
is the orthogonal matrix

  T_{sj} = ⟨((i₁ i₂) s, l) J | ((i₁ l) j, i₂) J⟩ = (−1)^{4i₁ + 3i₂ + 3l − s − j} ŝ ĵ {i₁ i₂ s; J l j},        (5)

which the engine already has (`physics/src/Recoupling.cpp:27-36`, `SpinFirstJJ`, checked there against an m-scheme
sum). In the channel-spin basis

  M_{ss′} ≡ ⟨((i₁ i₂)s, l)J | l·i₁ | ((i₁ i₂)s′, l)J⟩ = Σ_j T_{sj} λ_j T_{s′j}.        (6)

`work/so-note/angular.py` builds l·i₁ as a matrix on |l m_l⟩|i₁ m₁⟩|i₂ m₂⟩ and the channel-spin states with
AZURE2's coupling order from Clebsch–Gordan coefficients, and checks (6) against the direct overlap, (5) against the
m-scheme overlap, TᵀT = 1, the eigenvalues of M, the diagonal case (4), the channel-spin operator (below) and
l·i₁ + l·i₂ = l·s: 308 cases (2i₁ = 1–3, 2i₂ = 0–4, l = 0–4, all J), worst deviation 5 × 10⁻¹⁵. Example, p + ⁷Li
(i₁ = ½, i₂ = 3/2), l = 1, J = 1 (s = 1, 2; j = ½, 3/2): T = [[−0.4082, −0.9129], [−0.9129, 0.4082]],
M = [[0.25, −0.5590], [−0.5590, −0.75]].

**The coupled radial problem.** With the channel-spin components u_s(r) of one (α, l, J),

  u_s″ = (2μ/ħ²) Σ_{s′} { [V_C + V_N + ħ²l(l+1)/(2μr²) − E] δ_{ss′} + v_so(r) 2M_{ss′} } u_{s′}.        (7)

All s share k and η (same partition, same l), so the matrix outgoing solution with 𝐇(r) = H⁺_l(kr) 𝟙 for r ≥ r_m is

  𝐇(r) = T diag(H_j(r)) Tᵀ,        (8)

where H_j solves the scalar problem with λ_j: T is constant and orthogonal, so (8) solves (7), and at r ≥ r_m
𝐇 = H⁺_l T Tᵀ = H⁺_l 𝟙. The channel functions become symmetric matrices over the s of one (α, l, J):

  𝐋 = ρ 𝐇′ 𝐇⁻¹ = T diag(S_j + iP_j) Tᵀ ≡ 𝐒 + i𝐏,  Ω_hs = T diag(e^{−iφ_j}) Tᵀ,  𝐇 = √ρ 𝐏^{−½} Ω_hs^{−1}.   (9)

**The collision matrix.** In the external region the matrix solution for unit incoming flux in each channel is
𝐔(r) = 𝐈(r) − 𝐎(r) U with 𝐎 = 𝐇 and 𝐈 = 𝐇* (the potential is real). The R-matrix condition at a_c,
𝐔 = 𝐑 (ρ𝐔′ − B𝐔) (𝐑 = Σ_λ γ_λγ_λᵀ/(E_λ − E), dimensionless with the usual √(ħ²/2μa) factors, which are common
to the coupled channels), gives [1 − 𝐑(𝐋 − B)] 𝐎 U = [1 − 𝐑(𝐋* − B)] 𝐈, so

  U = 𝐎⁻¹ 𝐈 + 2i 𝐎⁻¹ [1 − 𝐑𝐋⁰]⁻¹ 𝐑 𝐏 𝐈,  𝐋⁰ = 𝐒 − B + i𝐏,

and with (9), since 𝐏, Ω_hs and 𝐒 commute (all are functions of the same T),

  U = Ω [𝟙 + 2i 𝐏^½ γᵀ A γ 𝐏^½] Ω,  (A⁻¹)_{λμ} = (E_λ − E) δ_{λμ} − γ_λᵀ 𝐋⁰ γ_μ,  Ω = e^{i(σ_l − σ₀)} Ω_hs.   (10)

This is conventions.rst Eq. conv-u with the diagonal P_c^½, Ω_c and L⁰_c replaced by blocks. `work/so-note/coupled.py`
checks it on p + ⁷Li, l = 1, J = 1 (WS 30 MeV / 2.5 fm / 0.65 fm, V_so = 5 MeV, a_c = 4.0 fm, E = 1.0 MeV): the 2 × 2
system (7) integrated directly agrees with (8) built from two scalar runs to 1.4 × 10⁻¹³; 𝐋 is symmetric to
1.5 × 10⁻¹⁷ and 𝐏 positive definite; the factorisation (9) holds to 1.1 × 10⁻¹³; U of (10) with one level is unitary
and symmetric to 4 × 10⁻¹⁶.

**What "the channel function of channel c" then means.** There is none: a channel spin s of an (α, l, J) with two
allowed s has a row of a 2 × 2 (in general n_s × n_s) block. Keeping only the diagonal of the block (a "restriction")
is not the solution of any equation: it drops the channel-spin flip that the potential alone produces (|U₁₂|² =
3.2 × 10⁻⁷ above, zero in the restriction) and changes resonant s → s′ transitions — in the example |U₁₂|² by −7 % and
+12 % for two amplitude vectors (coupled.py, part 4). It is refused (§3.4).

Further consequences of the coupled treatment (option (b) below), each derived here only as far as stated:

- **Brune and Park.** S enters the level-parameter transformations through bilinear forms Σ_c γ_λc (S_c − B_c) γ_μc;
  with a block 𝐒 they become γ_λᵀ(𝐒 − B)γ_μ. That Brune's (2002) alternative parameters follow with this substitution
  is my reading of the structure of the derivation (only symmetry of 𝐒 is used), *unchecked* against Brune's paper,
  which I did not consult. `level_parameters` today assumes a diagonal S (`ChannelValues`,
  `physics/include/azure2/physics/LevelParameters.h:43-61`).
- **Observed widths.** From the residue of (10) for one level: Γ_λc = 2 [(𝐏^½ γ̃_λ)_c]² / (1 + γ̃_λᵀ 𝐒′ γ̃_λ); a level
  "of channel spin s" decays into s′ through the off-diagonal 𝐏^½.
- **External capture.** The entrance wave in channel s is the column s of 𝐔(r); an E/M multipole that conserves
  channel spin connects the final channel s_f to entrance s_i through 𝐔_{s_f s_i}(r), so the radial integrals become
  matrices over (s_f, s_i).
- **Incomplete channel lists.** (7) needs every s of the (α, l, J). Projects list channels selectively: T3's 2⁺ block
  lists only s = 2 of the two p + ⁷Li l = 1 channel spins (conventions.rst, end of "Conventions"). The coupled
  treatment must then add the missing s as channels (with zero widths) or refuse.

**A spin–orbit term on the channel spin** (l·**s**, s = i₁ + i₂) is diagonal in (α, l, s, J) with
⟨2 l·s⟩ = J(J+1) − l(l+1) − s(s+1) (checked in angular.py). It is a different, phenomenological operator: it equals a
projectile term plus a target term of the same strength and shape (l·i₁ + l·i₂ = l·s, checked in angular.py), which
no optical-model parameter set describes; a projectile spin–orbit strength from the literature does not carry over to
it. For i₂ = 0 it coincides with §2.2.

**The options.** (a) Support the diagonal case only: exactly one particle of the pair has spin, and the term acts on
that spin; the key is (α, l, J). Optionally, as a separate, explicitly named operator, the channel-spin term. (b) The
full coupled treatment (7)–(10), with matrix blocks in the channel functions, the level matrix, the collision matrix,
the level-parameter transformations, the observed widths and external capture. **Recommendation: (a) without the
channel-spin operator as the first step.** It covers the cases a spin–orbit potential is usually wanted for (n, p, ³He,
t on even–even targets; α, ¹²C on any target with the term on the target spin; d on spin-0 targets), changes no
data structure of the R-matrix (L⁰ stays diagonal), and has independent tests (§4). (b) touches the level matrix,
`level_parameters` and external capture, each of which has its own design note and bit-identity gate; it should
follow only if a project needs it (p + ⁷Li, n + ⁷Be, p + ¹⁵N, p + ¹⁹F …). A pair with two non-zero spins and a
spin–orbit term is refused in step (a) with a diagnostic that names the pair and this note (§3.4).

### 2.4 What becomes j-dependent

In option (a), j ≡ J for every channel of the pair. Per quantity:

| quantity | equation | j-dependent | can today's path carry (α, l, J)? |
|---|---|---|---|
| P, S, e^{−iφ} (open) | (2) with (3) | yes | `ChannelTable::Key` is (pair, l, E) (`include/ChannelTable.h:23-33`): needs 2J; the slot loop knows J (`EPoint.cpp:1289-1316`) |
| e^{i(σ_l − σ₀)} | point Coulomb | no | — |
| S′, S″ at E_λ (Brune, Park) | d/dE of (2) | yes | `ChannelFunc` takes l only (`ChannelFunc.h:52-57`); `ShiftDerivKey` has no j (`CNuc.cpp:1592-1601`); the physics callback `ChannelFunctions(b, k, E)` already has the block's J (`LevelParameters.h:61`) |
| B_c | S at the lowest level | yes | held per block channel already (`physics_system.md`) |
| boundary condition at a_c | (2) | yes | as P, S |
| external capture, entrance | F, G of (1) + (3) on [a_c, r_m] | yes | `RadialBlockKey` (`CoulFunc.cpp:162-170`) needs 2J; `ECTable` `TermKey` already carries j_i, j_f (`ECTable.cpp:76`) |
| bound state, closed channels, ANC | Whittaker today (§1.4) | not today | below |
| THM vertex | point Coulomb today | not today | refused with a spin–orbit pair in step (a) |

**Bound states and closed channels.** The natural definition with a potential in the external region is the
decaying solution W̃_{lJ}(r) of (1) + (3) at E = −B with W̃ = W_{−η,l+½}(2κr) for r ≥ r_m, obtained by the same
inward integration (stable: the decaying solution grows inward). The asymptotic normalisation coefficient keeps its
meaning — the bound state is C W(2κr) at large r — and what changes is the factor at a_c: S(E < 0) = ρ W̃′/W̃, the ANC
↔ γ conversion μa/(ħ² W̃(a)²) (`CNucLevelParameters.cpp:35-40`) and the bound-state factor of external capture. Today
none of these is distorted, central or not (§1.4, item 1); with a spin–orbit term the omission would also hide its j
dependence below threshold. This is independent of the spin–orbit term and is proposed as its own step (§5, step 2).

## 3. Design (option (a))

### 3.1 The table key

`ChannelTable::Key` becomes (pair, l, two_j, E_c) with `two_j = 2J` for a channel of a pair whose potential has a
spin–orbit term and `two_j = kNoSpinOrbit = −1` for every other channel. For all present projects every key gets −1,
the order within a pair is unchanged and each row is the same function of the same arguments, so the tables and
outputs keep their bits. The same field enters, with the same rule, `RadialBlockKey` (`CoulFunc.cpp:162-170`),
`ShiftDerivKey` (`CNuc.cpp:1592-1601`) and `ECTable`'s `TermKey` (through the radial scratch; j_i is there). The
radial routine receives the number ⟨2 l·i⟩ of (4) rather than j, computed once from (l, 2J, 2i) where the key is
formed; the term is added last in the sum of (1) and not at all when ⟨2 l·i⟩ V_so = 0, so a potential without
spin–orbit evaluates V(r) with today's operations in today's order (radial.py, part B, checks that V(r) is
bit-identical on the grid).

### 3.2 Where the potential lives

`NuclearPotentialManager` is a process singleton with mutable state, tags from a counter and a state swap per pyazr
call: under ARCHITECTURE.md §9 it is not a constant table and not caller-owned state, and §9 lists it for removal
(`ARCHITECTURE.md:221`). It should become a value:

- `physics/potentials` (ARCHITECTURE.md §3): `struct FormFactor { Shape shape; double radius_fm, diffuseness_fm,
  width_fm; }`, `struct SpinOrbitTerm { double depth_MeV; FormFactor form; Particle acts_on; }`,
  `struct HybridPotential { double coulomb_radius_fm; double central_depth_MeV; FormFactor central;
  std::optional<SpinOrbitTerm> spin_orbit; double matching_radius_fm, step_fm; }` — immutable, compared by value,
  with the radial solver a free function of (HybridPotential, channel kinematics, l, ⟨2 l·i⟩, a_c, E).
- The System holds an `std::optional<HybridPotential>` per partition (physics_system.md §5 lists "the nuclear
  potential itself" as not yet included). The tag disappears: the pair signature of the channel table and the memo
  keys compare the potential's values.
- A pyazr session owns its System, so `ConfigScope`'s potential swap (`pybind_api.cpp:75-91`) and
  `exportState`/`importState` go with the manager.

This move is bit-identical and should precede the spin–orbit term (§5, step 1).

### 3.3 The grammar (format 3)

Format 3 has no potential record yet (`project_format3.md` does not mention one). Proposed, following its principles
(records with names, units stated, pairs by label, no number-changing defaults, one table of words):

```
<potentials>
potential p12C
    pair            p+12C                      # the pair label of <levels>; pair_key <key> only where ambiguous (D2/D3)
    coulomb_radius  3.5 fm                     # required (today's value is implicit)
    central         woods_saxon depth 20 MeV radius 3.6 fm diffuseness 0.6 fm
    spin_orbit      woods_saxon depth 6 MeV radius 3.6 fm diffuseness 0.6 fm on light operator 2l.s
</potentials>
```

- `central`: `woods_saxon depth radius diffuseness` or `gaussian depth width`, depth > 0 attractive, as today.
- `spin_orbit` (optional; absent = no term): Eq. (3), form factor `woods_saxon` (and `gaussian`, FRESCO's shape 2, if
  wanted), `on light|heavy` required, `operator 2l.s` required — the only value accepted; it is written so that a user
  copying an l·s-convention strength sees the factor (§6, Q2).
- `matching_radius 15 fm`, `step 0.0275 fm`: numerical settings with logged defaults (principle 6 allows those).
- Every word is a row of the schema table (`project_format3.md` §7), so `--schema`, pyazr and the web UI get it from
  there. The converter writes the record from a `<potential>` block (with `coulomb_radius 3.5 fm`), so a converted
  project gives today's bits.

### 3.4 Diagnostics

| code (proposed) | when | message names |
|---|---|---|
| `E_SO_TWO_SPINS` | spin–orbit on a pair with both spins non-zero (option (b) not implemented) | pair, both spins, this note |
| `E_SO_ZERO_SPIN` | `on` names a particle of spin 0 (the term would vanish) | pair, particle |
| `E_SO_RELATIVISTIC` | with `relativistic=1` (as today, `CNuc.cpp:202-207`) | — |
| `E_SO_THM` | the pair is the entrance pair of THM data (the vertex uses point-Coulomb waves) | pair, data set |
| `E_SO_NEEDS_J` | an API call for channel functions of a spin–orbit pair without J | pair |
| `W_SO_BELOW_THRESHOLD` | only if step 2 is not merged first: a sub-threshold level or closed channel of the pair (the term is ignored there) | pair, level |

### 3.5 pyazr and the web UI

- `set_nuclear_potential(pair, ..., spin_orbit=None | {"depth": MeV, "radius": fm, "diffuseness": fm, "on": "light"})`
  and the `NuclearPotentialSetting` view gain the term; `coulomb_functions(pair, energies, L, J=None)` takes J, required
  (`E_SO_NEEDS_J`) for a spin–orbit pair and refused otherwise.
- The web UI's form is generated from the schema rows (ARCHITECTURE.md §9, §11); `qt_parity.md` §8 gains the fields.
- pyazr's `params_rwa` caveat (the vector is re-derived when the potential changes) applies unchanged.

## 4. Tests with independent answers

### 4.1 Run for this note (`work/so-note/`, Python 3.12, numpy 2.5, scipy 1.18, mpmath 1.4; one thread)

| script | method | result |
|---|---|---|
| `angular.py` | m-scheme: l·i₁, l·i₂, l·s as matrices on the product space; channel-spin and jj states from Clebsch–Gordan (Racah's formula, exact rationals); (5) with a Racah 6-j | 308 cases, worst 5.3 × 10⁻¹⁵ (tolerance 10⁻¹²). Injected error: dropping (−1)^{3l} from (5) fails 2.0 |
| `radial.py` A | engine's Numerov (transcription of `computeHybrid`) vs scipy DOP853 (rtol 10⁻¹²) with mpmath boundary values; p + ¹²C, WS 20/3.6/0.6, V_so = 6 MeV, l = 1, 2, both j, E = 0.5–3 MeV | |ΔS| ≤ 7.8 × 10⁻⁸, ΔP/P ≤ 9.1 × 10⁻⁸, |Δφ| ≤ 5 × 10⁻⁸ rad |
| `radial.py` B | V₀ = V_so = 0 against exact Coulomb functions | 1.7 × 10⁻⁵ (engine's e²), 6 × 10⁻⁹ (α ħc); spin–orbit term with ⟨2l·s⟩ = 0 leaves V(r) bit-identical |
| `radial.py` C | closed form: neutral square well 30 MeV, R = 4.5 fm, with the Thomas term of its sharp edge, which is a δ shell −V_so (ħ/m_πc)² ⟨2l·s⟩ δ(r − R)/R (a jump of u′), solved with Riccati–Bessel functions; against DOP853 on Woods–Saxon edges of diffuseness d → 0, Richardson in d (the kink of u makes the error O(d)) | l = 1, 2, both j: 4–20 × 10⁻⁵ relative after Richardson (raw: 2–110 × 10⁻³ at d = 2.5 × 10⁻³ fm, halving with d). Injected error (l·s instead of 2l·s in the closed form): 0.30 |
| `coupled.py` | §2.3: coupled 2 × 2 integration vs (8); (9); (10) | 1.4 × 10⁻¹³; 1.5 × 10⁻¹⁷; 1.1 × 10⁻¹³; unitarity 4 × 10⁻¹⁶ |
| `threshold.py` | §1.4 item 1 | ΔS(0⁺) ≈ 0.75–1.07 |
| `engine_check.py` | the reference engine (`work/b-trial2`, pyazr `coulomb_functions`, one thread) on a copy of `tests/hybrid_potential` vs the transcription in `radial.py` | P to ≤ 1.2 × 10⁻⁸, S to ≤ 3.8 × 10⁻¹⁰, l = 0–2, 0.5–3 MeV: the transcription is the engine's algorithm (the residue is the masses I typed) |

Not run: anything with the engine with a spin–orbit term (it does not exist); a p₃/₂–p₁/₂ splitting against a
literature potential (I have no consulted source of parameters for n/p + ⁴He, ¹²C or ¹⁶O and cite none); the
gradient check.

### 4.2 Tests to build with the code

1. **Angular** (physics unit test): (4) for one spinning particle; (6) against an m-scheme overlap; refusal of two
   spins. Detects a wrong sign or factor 2 in ⟨2 l·i⟩ (fails at 100 %).
2. **δ-shell closed form** (physics unit test): the square well + sharp Thomas edge of radial.py C, evaluated in closed
   form in the test, against the engine's solver on a Woods–Saxon edge with Richardson in d, 10⁻³; and a constant
   `spin_orbit` with ⟨2l·s⟩ = 0 (s-wave) equal to the central result bit for bit. Detects the factor 2, the sign of
   df/dr and an l·s–l·σ mix-up (each fails at > 10 %).
3. **Woods–Saxon + Thomas against an independent integration**: P, S, φ of the engine (pyazr `coulomb_functions`
   with J) against DOP853/mpmath values tabulated by a committed script (radial.py A's method), 10⁻⁶ relative.
4. **V_so absent → today's bits**: the 179 × 3 comparison identical; `hybrid_potential` and T3 byte-identical.
5. **A new project `tests/spin_orbit_potential`**: `hybrid_potential` with a `spin_orbit` record; its χ² pinned, its
   p₁/₂ and p₃/₂ phase shifts (segment phase-shift output) compared with item 3's script. A second variant with the
   ¹²C(p,γ) capture of 13N to test the external-capture path (entrance wave per J).
6. **Unitarity**: |U| of every elastic-only block = 1 to 10⁻¹² with the term on (diagonal case: automatic, but it
   catches a row read with the wrong J).
7. **Gradient**: the analytic gradient against finite differences on project 5 (the S′ path with J).
8. **Thread count and order**: tables built in two orders and at 1 and 8 threads identical (`channel_table_unit`
   extended by two_j).
9. **Diagnostics**: one fixture per code of §3.4.

## 5. Steps

Each step is one branch, mergeable on its own, with the 179 × 3 gate.

1. **The potential as a value** (§3.2): `physics/potentials` with `HybridPotential` and the radial solver;
   `NuclearPotentialManager`, the tag and the session swap removed; the System carries the potential. Outputs:
   identical (the solver is the same arithmetic). Needs: nothing beyond the present head; collides with Phase 6's
   statics work only in the manager, which it deletes.
2. **Distorted closed channels and bound states** (§2.4): W̃ by inward integration, S(E < 0), the ANC conversion
   and the bound-state factor of external capture of hybrid pairs; continuity of S through threshold checked.
   Outputs move for hybrid projects with a sub-threshold level or capture (`hybrid_potential`: its bound ¹³N ground
   state), classified with cause and size; nothing else. Needs the PI's decision Q4. Depends on step 3b (i) of
   `cnuc_split.md` (Whittaker functions in physics/), if that lands first, to write it once.
3. **The format-3 `potential` record** (§3.3) and its converter row, without spin–orbit. [Done 2026-11-10,
   `feature/upgrade-4-potential-record`: `project_format3.md` §3.10a.] Outputs identical. Needs the
   format-3 reader of Phase 4 (`project_format3.md` §7).
4. **The key** (§3.1): `two_j` in `ChannelTable::Key`, `RadialBlockKey`, `ShiftDerivKey`, the scratch, always −1.
   Outputs identical (byte-identical tables in `channel_table_unit`).  [Done 2026-11-10,
   `feature/upgrade-4-spin-orbit`: one function forms the key, `PPair::SpinOrbitTwoJ(J)`; also the THM shift memo
   and the rate grid of `ECIntegral`.  The THM Coulomb term passes `kNoSpinOrbit` (the vertex refuses any potential).]
5. **The spin–orbit term, option (a)**: the term in the solver, ⟨2 l·i⟩ from (l, 2J, 2i), the `spin_orbit` field,
   the diagnostics, pyazr's `J` argument, tests 1–9 of §4.2. Outputs identical unless a spin–orbit term is given.
   [Done 2026-11-10, `feature/upgrade-4-spin-orbit`: `physics/HybridPotential.h` (`SpinOrbitTerm`, `TwoLDotI`,
   `SpinOrbitRadialTerm`; the solvers take ⟨2 l·i⟩); `PPair::SpinOrbitFactor`; the record word (woods_saxon only);
   the codes `E_SO_TWO_SPINS`, `E_SO_ZERO_SPIN`, `E_SO_RELATIVISTIC`, `E_SO_THM`, `E_SO_NEEDS_J` and `E_SO_NO_TERM`
   (a J for a pair without the term) in `experiment/SpinOrbitChecks.h`; pyazr 2.29.0 (2.24.0 on the branch); tests: `physics_spin_orbit`,
   `tests/spin_orbit_potential`, `pyazr_spin_orbit_potential`, `channel_table_unit` with both j.  Not done: the
   web UI / Qt forms; `gaussian` as a spin–orbit form.]
6. **Option (b)**, only on demand: matrix blocks in `EDependentTables`/the per-stage table, `AMatrixFunc`,
   `level_parameters`, observed widths and external capture; its own design note (the A-matrix and the level
   parameters are owned by `level_parameters.md` and the plan's S5 / CNuc step 6b, which should be merged first).

Steps 1, 3 and 4 do not depend on each other; step 5 needs 1, 3 and 4, and should follow 2 (otherwise
`W_SO_BELOW_THRESHOLD` is needed).

## 6. Questions for the PI

**Decided 2026-11-10.** The PI: yes to questions 1, 2, 4 and 5. The director: 3 no; 6 yes; 7 no; 8 yes. The
questions stay below as the record of what was decided.

1. **Target spin.** Support only pairs with one spinning particle first (option (a)), refusing the others, and do the
   coupled treatment only when a project needs it? *Recommendation: yes.*
2. **Convention.** V_so as in FRESCO (Eq. (3): (ħ/m_πc)² = 2.000 fm², operator 2 l·s = l·σ), written as
   `operator 2l.s` in the file so that the factor is visible? *Recommendation: yes; no other operator accepted.*
3. **Channel-spin spin–orbit.** Offer l·**s** (s = channel spin) as a separate, named phenomenological operator for
   pairs with two spins? *Recommendation: no, not in the first step; it does not correspond to any optical-model
   parameter set.*
4. **Bound states and closed channels.** Distort them with the same potential (step 2), which removes today's threshold
   step of S (~1 for `hybrid_potential`) and moves that project's outputs? *Recommendation: yes, before the
   spin–orbit term.*
5. **e² of the hybrid Coulomb potential.** Use α ħc (one set of constants, `repo_consistency` check 10's intent),
   which moves hybrid outputs by ~10⁻⁵, or keep 1.4399764 as a named exception? *Recommendation: α ħc, in step 1 or
   its own commit, outputs classified.*
6. **Coulomb radius.** Keep R_C = 3.5 fm for every pair as the converter's explicit value, and make it a required
   field of the record? *Recommendation: yes.*
7. **Fit parameter.** May V_so (and V₀, R, a) be fit parameters? The per-stage channel table assumes rows that do not
   depend on θ (`channel_table.md` §2); a free potential parameter makes every row of the pair depend on it.
   *Recommendation: not in this work; a potential parameter is scanned (profile) by rebuilding, as now.*
8. **Matching radius and step.** Expose `matching_radius` and `step` as fields with logged defaults (15 fm, 0.0275 fm)?
   *Recommendation: yes; today they are constants with dead setters.*
