# experiment/: what was measured and what the model predicts

This module holds the part of the comparison between the R-matrix model and the data that depends on the
measurement and not on the level parameters: the model expression of a data set (which calculated observables a
measured value corresponds to, and how they combine) and, from step S6, the conversion of a measured point from the
laboratory to the centre-of-mass frame. It is the `azure2_experiment` CMake target, namespace `azure2::experiment`,
and uses `physics/` (constants, kinematics) and the C++ standard library only (`docs/upgrade/data_model.md` §2.5,
decision M1). Three checks keep it so: the target links only `azure2_physics` and sees no engine include directory;
the ctest `experiment_standalone` configures and builds this directory on its own; `repo_consistency` (check 9b)
refuses any other include.

## 1. Model expressions (`Expression.h`, step S8)

A measured point i (beam energy E_lab,i in MeV, lab or c.m. angle θ_i in deg) is compared with a model value m_i.
m_i is the value at i of a tree whose leaves are calculated observables (σ in b, dσ/dΩ_c.m. in b/sr, A_y, a phase
shift in deg, ...), each already folded with its apparatus by the calculation unit that computed it (decision M4),
and whose nodes combine them point by point:

| node | m at point i | cotangents passed to the children (`backward`) |
|---|---|---|
| `Leaf` | v_i, the leaf's value | — |
| `Sum` | v₀ + v₁ + … + v_{n−1}, from v₀ (0 for n = 0) | g to every child |
| `Scale` | (c J_i) v_i; c a constant, J_i a per-point factor (the frame factor J_base/J_k), 1 if absent | g (c J_i) |
| `Weighted` | v_i a; a a parameter (abundance a_{s,n}) | g a to the child; g v_i to a |
| `Ratio` | v_num / v_den (J_den / J_num), NaN where v_den = 0 | g (m/v_num), g (−m/v_den) |
| `Transform` | F(s, s_t) at every sub-point of the base unit, then its kernel | see below |

with g the cotangent of the node (∂χ²/∂m for the root). The transmission is F = T = exp(−n s), the self-shielded
yield F = Y₀ = (s/s_t)(1 − exp(−n s_t)); n is the sample thickness (atoms/b, a parameter), s = Σ_q w_q v_q the
weighted sum of the reaction terms and s_t the total term of a yield, with weights w_q = (c_q J_q,i) a_q.

**Order of operations.** Every value is formed as `src/ESegment.cpp` forms it today (CalculateTheoreticalCrossSection
1698–1822, TransformedValue 1605–1696, Transmission, YieldY0, TermFrameFactor 1487–1493): a SUM starts from the base
term and adds the terms in component order; a term is (c J) v and then times a; a ratio is num/den times (J_den/J_num);
inside a transform the weight w = (c J) a is formed first and then multiplied by v, and s starts from 0. The test
`experiment_expression` compares each node and each composition with these expressions written out, bit for bit, on
4000 points with random 53-bit mantissas.

**backward** is the adjoint of each node in the order `src/AZUREGrad.cpp` uses (GradOnePoint's shared-component route
for SUM and RATIO, GradTransformedPoint for transforms). Bit for bit with those expressions: Sum, Scale, Weighted on a
leaf, Ratio, and every Transform case (point values, sub-points, free-gas Doppler with and without a kernel, the
thickness, abundance and Doppler-temperature cotangents). Not bit for bit: a weighted *and* scaled SUM term (a
sample's component). The tree's reverse sweep forms (g a)(c J) for the leaf and g ((c J) v) for the abundance;
AZUREGrad forms g ((c J) a) (shared route) or (g (c J)) a (per-term route) and (g (c J)) v. The two differ by at most
2 ulp (2.0 × 10⁻¹⁶ relative in the test). No per-node adjoint can give AZUREGrad's products, since they combine the
factors of two nodes before g; S11 classifies the gradient bits this moves.

**Departures from AZUREGrad, deliberate.** At v_num = 0 the ratio's numerator gets g (J_den/J_num)/v_den (AZUREGrad
gives 0, because it forms m/v_num); where v_den = 0 both cotangents are NaN (AZUREGrad gives 0, with a NaN model
value).

**Transform and the kernel.** The transform acts between the sub-point values of its terms and the base unit's
kernel. The library does not know the kernel: the caller passes a `SubPointFold`, one virtual call per point, which
returns the layout at that point (`TransformStage`: no sub-points, or nsp sub-points; an optional linear map of them
before F, the free-gas Doppler broadening that SAMMY applies first, in compressed rows; whether a kernel follows), the
kernel applied to a local vector of values, and the kernel's exact weights for the backward pass. In the engine this
is `EPoint::IntegrateKernelsImpl` with a reader over the local vector, under `DopplerAppliedScope`, and
`EPoint::TargetEffectWeights`; nothing is written into the sub-points (W4 of the note goes). A kernel that is not linear
makes `Backward` return false, as AZUREGrad falls back to finite differences there.

**Storage and evaluation.** A `Tree` is a vector of `Node` values in creation order. A node is created after its
children, so the creation order is a post-order: `Evaluate` is one forward sweep over the vector, one switch per node
and one loop over all points inside it (no recursion, no pointer graph, no virtual call per point), into a node ×
point buffer of the caller's `Workspace`. `Backward` recomputes the values and sweeps the vector backwards. A tree holds
structure and constants only; per-point values arrive in `Inputs`:

- a leaf reads `Inputs::leaves[k]`, a strided view (value at i is `data[i * stride]`), so the engine can point it at
  the slots of the calculation units (`Results::slot`, a unit's quantities interleaved) without copying;
- a transform term reads the same leaf's sub-point values in `Inputs::subLeaves[k]` (a run per point, given by
  offsets, as grids differ between points);
- frame factors J_i are per-point arrays (`Inputs::factors`); abundances, sample thicknesses and Doppler temperatures
  are parameters (`Inputs::params`), each node refers to them by index.

One tree therefore serves every evaluation and every working copy, and it is const during an evaluation. `Natural`
gives the per-nucleus terms of a sample (every `Weighted` node's value added to its nucleus in creation order, and a
transform's w_q v_q), as `SetNaturalTerms` records them now.

**Not here yet.** The tree is not wired into the engine (S10), and no keyword builds one (S9). `ESegment::Transmission`,
`YieldY0` and `TermFrameFactor` call the library's (one formula each).

## 2. Frames (`Frames.h`, step S6)

A measured point is given in the laboratory: beam particle 1 of the entrance partition with kinetic energy E_lab (MeV)
on particle 2 at rest, and for lab-frame differential data the detected particle 1 of the exit partition (a photon
for a capture exit) at θ_lab (deg). The comparison is made in the c.m. frame of the entrance partition, so each point
needs its c.m. energy E_c.m. = E_lab m₂/(m₁ + m₂) (relativistically √s − m₁ − m₂, s = (m₁ + m₂)² + 2 m₂ E_lab), the
excitation energy E_x = E_c.m. + S + E*, its c.m. angle θ_c.m. (elastic, Marion–Young for a reaction, the boost of
the detected particle with relativistic=1, photon aberration) and the Jacobian J = dΩ_lab/dΩ_c.m. that multiplies the
measured value and its uncertainty, dσ/dΩ_c.m. = J dσ/dΩ_lab.

The library holds these formulas as functions of plain numbers (masses in u, energies in MeV, angles in deg; a
`Partition` is m₁, m₂, S, E*), in the engine's order of operations: `BeamToCM`, `DecayToCM`, `ElasticLabToCM`,
`ReactionLabToCM`, `RelativisticLabToCM`, `CMVelocity`, `PhotonLabToCM`, `ParticleJacobian`, `PhotonJacobian`, and
the relativistic energy maps `CMFromLab`, `LabFromCM`, `CMFromDecay` (which `include/Kinematics.h` now calls). A map
that does not reach the lab angle returns false with an `AngleLimit` (the maximum lab angle and the reason).

The engine keeps the state rules: `EPoint`'s conversion functions read the point's lab values and pairs, call the
library and store the results (and the maximum and reason, as before); `ESegment::ConvertPointToCM` is the one
sequence (energy, angle, value) that `ESegment::Fill`, `UpdatePointEnergyWithShift`, the test segments and the
component segments call. Where the lab angle is beyond the kinematic maximum a caller either stops (load) or converts
at the last valid angle (an energy shift mid-fit). `ESegment::CMEnergyFactor` is the one E_c.m./E_lab factor of the
sub-point shift and the target-effect grids.

Tests (`experiment_frames`): the engine's expressions before S6 written out, against the library on 20 000 random
settings each, bit for bit; the elastic closed form tan θ_lab = sin θ_c.m./(cos θ_c.m. + m₁/m₂) (3.8 × 10⁻¹³);
the inverse aberration; every Jacobian against a numerical d cos θ_lab / d cos θ_c.m. of the library's own angle map
(≤ 1.3 × 10⁻¹⁰); the relativistic branch against Lorentz boosts built in the test from four-vectors (E_c.m. from the
invariant mass, θ* boosted to the lab and converted back, 8.5 × 10⁻¹³; J = d cos θ_lab / d cos θ*, 2.9 × 10⁻¹⁰).
