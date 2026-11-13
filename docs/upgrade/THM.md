# AZURE2 upgrade: Trojan horse method (THM)

Goal: a THM implementation that is simple, valid and natural in the physics, as Marco La Cognata asked in his review
(2026-11; `~/THM/review/marco_2026-11/` on the analysis machine), with precise DWBA vertices. THM is one more
observable and one more kernel in the three lines of `ARCHITECTURE.md`; it shares the level matrix with the direct
data, so THM and direct measurements are fitted together without special code paths.

## 1. Physics, in the notation of the three lines

Reaction A(a, c C)s with Trojan horse a = (x s), quasi-free (QF) break-up, spectator s. Three-body cross section in
the plane-wave impulse approximation, with KF the kinematic factor and φ(p_s) the x–s momentum distribution:

```
d³σ/(dE_c dΩ_c dΩ_C) = KF · |φ(p_s)|² · dσ^HOES/dΩ        (Ref. [21] = La Cognata et al. ApJ 708 (2010) 796, Eq. 15)
```

The half-off-shell (HOES) cross section of x + A → c + C at E = E_xA (c.m.), modified R-matrix (Ref. [21] Eqs. 16–21):

```
σ_HOES(E) = Σ_J (2J+1) Σ_{c'∈f} K(E) · 2P_{c'}(E_f) · Σ_{(s,l)} | Σ_{λμ} γ_{λc'} A_{λμ}(E) · W_{μ,(s,l)}(p) |²
```

- A_{λμ}(E): **the same level matrix** as the direct data (one solve per grid energy E_n, `ARCHITECTURE.md` §1).
- W_{μ,(s,l)}(p): the **entrance vertex** for level μ, with p = √(2μ(E + B_xs + T_s))/ħ the off-shell x–A momentum.
  It replaces √(2P_c) e^{iφ_c} of the on-shell entrance channel. Vertex models (§3): plane-wave surface term
  M_l(p) = (B − 1) j_l(pa) − pa j_l′(pa); full form factor W = (ħ²/2μ) u(a) M_l + (E_λ − E − B_xs − T_s) ∫₀ᵃ u r j_l dr
  (Green identity, verified to 2.6e-10 on the `thm-review` branch, commit 4882ad0); DWBA vertex (§3).
- K(E): the kinematic factor left after dividing by KF·|φ|² and integrating. Default `lacognata` K = k_f/μ_f
  (Tumino et al. 2021 Eq. 50; ApJ 723 k_f dependence). Open question for Catania (Q1): μ in the numerator (ApJ 723
  Eq. 14) or the denominator (Tumino Eq. 50). Other forms (`triple`, `kf3body`, `lambda32`) are kept only to reproduce
  old projects and are marked legacy.

In the three lines: O_n^THM(θ) is σ_HOES at E_n (an observable), the THM kernel K_THM(ν) folds it to the data, and
s_k is the THM normalization.

## 2. Data, kernel and normalization (keep it simple)

- **Two data kinds, never mixed silently.** HOES excitation functions (arbitrary units; fitted with σ_HOES) and
  on-shell-equivalent (OES) S-factor points (penetrability already divided out by the authors; compared with the
  on-shell S(E)). The project states which one each THM measurement is.
- **Kernel K_THM**, one sparse matrix per measurement: the experimental resolution as a Gaussian with σ_E stated in
  the **c.m.** frame of x + A (a lab value must be converted and the file says which); the QF window (spectator-momentum
  cut) through |φ(p_s)|² weights when the data are not already divided by it; for OES points the fold acts on S(E),
  not σ(E) (La Cognata 2009 folds S with σ = 20 keV; folding σ instead inflates the result ×182 at 50 keV).
- **Normalization.** Either a free scale s_k (profiled analytically) or fixed by a reference resonance (double ratio,
  Ref. [21] Eq. 22: (ωγ)_i = (ω_i/ω_3)·[Γ_{p,i}/dσ_i^tr]/[Γ_{p,3}/dσ_3^tr]·(N_i/N_3)·(ωγ)_3) entered as a prior term
  on s_k. No energy shifts to move THM data onto direct data; no negative backgrounds (a background, if any, is
  non-negative by construction).
- **Identical particles** (¹²C + ¹²C): only even J⁺ levels; (1 + δ₁₂) in the rate; resolution frame stated.
- **Reproduce first.** Every published THM analysis we re-analyse is first reproduced with its own settings (C7 in the
  review plan: ¹⁷O Oliva, ¹⁸O La Cognata 2010, ¹⁵N La Cognata 2009, ¹⁹F Indelicato 2017).

## 3. DWBA vertex: precise and validated

`physics/transfer/` computes the vertex functions W(p) on the THM grid once per fit stage (they do not depend on the
fit parameters, except through B_λ in the per-level basis, handled in the level-parameter step):

1. **Bound state** of x in a: Woods–Saxon (+ Coulomb) well, depth fixed by the binding energy B_xs, node count from the
   Wildermuth rule for a nucleon on a core, given explicitly for clusters; spectroscopic factor or ANC as input.
2. **Level function** of x + A inside the channel radius a: depth fixed by a·u′/u = B_λ at E_λ (as `vertexModel=ff`).
3. **Distorted waves** for a + A (entrance) and s + F* (exit) from optical potentials: global sets as defaults
   (Koning–Delaroche for nucleons, An–Cai for deuterons, Pang (GDP08) for ³He/t, Avrigeanu for α), each overridable.
   Integrated with Numerov through pycoulomb's hybrid solver, matched to Coulomb functions.
4. **Form factor**: post-form DWBA amplitude ⟨χ_s F| V_xs |χ_a A⟩ for transfer to the x + A continuum at energy E,
   first in zero range as a check, then **finite range**; the HOES vertex is its projection on each channel (s, l).

Validation (the independent check required by `CLAUDE.md`):
- plane-wave limit (optical potentials off) reproduces M_l and the full form factor to 1e-6;
- FRESCO benchmarks for each THM system we analyse (d → p + n, ³He → d + p, ⁶Li → α + d, ¹⁴N → ¹²C + d): angular
  distributions and form factors agree to 1 %;
- Green identity and probability conservation of the bound state;
- the ratio of vertices between two levels (e.g. ¹⁸O 20 keV / 144 keV) is stable against the potential choice within
  a stated systematic, reported as a model uncertainty.

Status (2026-11-11, branch `feature/upgrade-7-transfer-1`, `transfer.md`): `physics/transfer` exists as a leaf
library, not yet called by the engine. Done: the bound state of item 1 (depth search, ANC, D0; FRESCO to its printed
digits), the distorted waves of item 3 with explicit potentials (spin 0, ½, 1 with spin-orbit; FRESCO to 1e-4), the
post-form transfer of item 4 to a *bound* state in zero and finite range (FRESCO, ¹⁵N(d,n)¹⁶O and ¹²C(⁶Li,d)¹⁶O, to
0.2 %), and the plane-wave limit of the vertex (M_l and the full form factor to 1e-6). Open: the level function of
item 2 inside a, the transfer to the x + A continuum (FRESCO bins), the global optical sets, the ³He and ¹⁴N
benchmarks, the wiring (`transfer.md` §11–12).

Status (2026-11-12, branch `feature/upgrade-7-transfer-2`, `transfer.md` §14–§18): the transfer to the x + A
continuum exists (scattering states, FRESCO's bins, d²σ/dΩ dE as the narrow-bin limit, finite and zero range), within
0.5 % of FRESCO for ¹⁵N(d,n) to three p + ¹⁵N bins; the interior continuum amplitude is the level function's vertex
(item 2's level function: u_E at E_λ inside a) times A_λ(E), exactly at E_λ and to 0.1–1.9 % at E_λ ± Γ, with an
exterior part of 9 % (a = 5 fm) to 0.5 % (7 fm) (`transfer.md` §16); the prior form and the remnants (post = prior
to 10⁻³). Open: the ⁶Li bins against FRESCO (FRESCO loses accuracy there), ³He and ¹⁴N, global sets, wiring.

Status (2026-11-12, branch `feature/upgrade-7-transfer-3`, `transfer.md` §19–§23): spin-orbit in the distorted waves
of a transfer (both channels; the projectile's iT11 from the same amplitudes), FRESCO to 2e-4 for ¹⁵N(d,n); the ³He
Trojan horse ¹⁸O(³He,d)¹⁹F (g.s. and 5/2⁺, finite and zero range, with and without spin-orbit) to ≤ 8.9e-4 and the
¹⁴N Trojan horse ¹²C(¹⁴N,d)²⁴Mg to ¹²C + ¹²C states (0⁺, 2⁺, finite range) to ≤ 8.3e-4 of FRESCO. The ⁶Li bins remain
open: Frescox 7.4 shows FRESCO 3.5's numerics (§23). Open: ³He and ¹⁴N to the continuum, global sets, wiring.

Status (2026-11-13, branch `feature/upgrade-7-thm-wiring`, `transfer.md` §24): the global optical sets are one table of
the library (`GlobalOptical.h`, `src/ThmOptical.cpp` deleted); the distortion factor R(E) calls the library for the
waves, the tail, the sum and M_PW (`ThmDistortion` keeps kinematics, windows, E_ref, the ln R grid); the DW vertex
takes its waves from the library. Numbers moved ≤ 9e-9. Open: the DW vertex's source S_lm(a) as a library function
(§24.2), spin-orbit terms of the global sets, a THM weight with spin-orbit distortions.

Open questions for Catania (formula sheet Q2): their potential geometry and radius convention (r0·A^{1/3} vs
r0·(A1^{1/3} + A2^{1/3})), spin–orbit, Coulomb radius, node rule for clusters, whether the interaction radius equals
the channel radius, which boundary condition fixes the level function, whether u is normalised inside a, and one
benchmark number (e.g. Γ/(dσ/dΩ) of the ¹⁸O 20 keV level relative to 144 keV).

## 4. What is plottable (Physics tab of the web UI)

M_l(p) and W(p) per level and channel, the momentum distribution |φ(p_s)|², bound-state and level functions u(r),
distorted waves, K(E), the THM kernel rows, and σ_HOES next to the on-shell σ for the same levels.

## 5. Where THM stands today and how it moves

- Today's THM code lives on `origin/thm` (30d3bd9) and the local branch `thm-review` on the analysis machine:
  `backgroundPositive=on` (073fe48 + GUI 08ccea3), `vertexModel=ff` full form factor (4882ad0, engine/pyazr/docs/tests;
  GUI part written, not yet compiled). The THM re-analyses (D1–D6 of the Marco review plan) continue there until the
  upgrade reaches Phase 5.
- In the upgrade, THM is rebuilt in the new structure: σ_HOES as an observable (Phase 5), the THM kernel and the
  double-ratio prior (Phase 5), `physics/transfer` DWBA (Phase 7, can start in parallel because it is a leaf module).
  The `.azr` `<thm>` block stays readable; its options map onto the new objects; legacy options are marked.
- On the shared mesh (Phase 5 branch 10, 2026-11-12): σ_HOES at E_n is a kind of site observable (`kHoes`,
  `include/SiteObservable.h`), THM points are folded on the mesh of their compound nucleus with the kernels' table, the
  level matrix of a node is solved once for THM and direct readers, and the scalar gradient takes THM data through the
  energy tasks (`include/ThmCotangents.h`).  A new vertex model: its vertex at E_n (`ChannelTables.cpp`, the site's
  vertex pieces) and its branch in `ThmEntranceVertices` (`THMMatrixFunc.cpp`) with the adjoint's; nothing of the
  mesh, the slots, the passes or the gradient changes.
- THM tests: the THM part of T2 (synthetic ²H(¹⁵N,α¹²C)n HOES and OES data from known parameters), the Green
  identity, the plane-wave limit and the FRESCO benchmarks.
