# Phase 5: one energy grid, the apparatus as kernels, lab-frame comparison, finite layers — design

This note states what Phase 5 changes in *where* the engine computes what, measured on the code at commit 2363abac
(the head of `feature/upgrade-4-data`, data classes `DataSet` / `WorkingCopy` / `Unit`). It covers the grid of
compound-nucleus energies E_n on which the collision matrix is solved once for all reactions, the observables O_n on
that grid, the apparatus kernels K that take them to the measured quantities, the comparison of lab data with a
fitted energy shift in the frame of the data (decision M3 (b) of `data_model.md`), the yield of a target layer of
finite thickness, and the wiring of the stopping module. It ends with the gate, the branches and the decisions the
PI must take. "Measured" means run on this host with the commands of §1.6; "inferred" means read from the code and
not run.

The main result of the measurements, in one paragraph: today an evaluation solves the level matrix once per
calculation unit (a data point or a fold sub-point), 207 to 291,847 times per evaluation in the projects measured,
but the units of a project lie on far fewer distinct compound energies — **bit for bit**: 13N has 45,872 solves on
2,728 distinct energies E_x (16.8×), `thick_target_resonance` 48,239 on 2,149 (22×), `18O_p_a_thm` 14,689 on 704
(21×), `natural_neon` 291,847 on 56,703 (5.1×), T2 24,300 on 5,442 (4.5×). The sub-point grids of neighbouring points
already coincide, because `AdaptiveIntegrationGrid` anchors them to a quantized lattice. One solve per distinct
energy is therefore available **without changing any output**, before any new quadrature. A shared mesh at the node
density of today's grids would add at most another factor 2.2 (T2, `natural_neon`) and nothing for 13N, the thick
target and THM. The design below therefore does the bit-identical part first, and the new quadrature later, for
memory and deletion rather than for speed.

## 1. What is computed today, measured

### 1.1 The path of one evaluation

An evaluation (`AZURECalc::operator()`) runs the hand-over (`WorkingCopy::FillNormsFromParams`,
`FillEnergyShiftsFromParams`: new frames for shifted segments, `ReplaceSubPointGrids`, re-energised units, a new
plan), then the point pass (`WorkingCopy::CalculatePointsInParallel`, one OpenMP task per plan task), then the model
pass (`azure2::EvaluateModel`, the expression trees) and the objective. A plan task is a host unit and the units
mapped onto it (`WorkingCopy::MapData`: same entrance pair, same exit pair, same energy, `SharesCalculation`; the
components of a SUM/RATIO/total-capture segment onto their base point, `MapLeaves`). A host with a fold calculates
its sub-points inside `Unit::Calculate` (`src/PointPass.cpp`), each with its own `AMatrixFunc::FillMatrices` →
`InvertMatrices` → `CalculateTMatrix` → `CalculateCrossSection`, then folds them (`Unit::IntegrateKernelsImpl`,
`src/TargetFold.cpp`). Sub-points are never mapped across points: two points whose windows overlap solve the level
matrix twice at the same node.

The level matrix of a Jπ block depends on the compound energy E_x alone (and θ): `AMatrixFunc::FillMatrices`
takes `InEnergy(point)` = E_c.m. + S_α + E_x,α of the entrance pair (of the exit pair for a β-delayed entrance), and
every channel energy in the unit's tables is formed from that number (`Unit::ChannelKeys`,
`src/ChannelTables.cpp`). The channel functions themselves are already de-duplicated by the per-stage channel table
(Phase 1): T2 asks 190,748 rows of which 33,207 are distinct, 13N 367,264 → 13,635. The level-matrix solve is not.

### 1.2 Units, solves, distinct energies (measured)

| project | points | sub-points | units | plan tasks | level-matrix solves per evaluation | distinct E_x of the solves (bit for bit) | solves / distinct | one shared mesh, estimate | distinct (E_x, entrance, exit) |
|---|---|---|---|---|---|---|---|---|---|
| T1 | 442 | 0 | 442 | 206 | 207 | 138 | 1.5 | 137 | 198 |
| T2 | 380 | 24,248 | 24,628 | 360 | 24,300 | 5,442 | 4.5 | ≈ 2,780 | 6,195 |
| T3 | 494 | 0 | 494 | 412 | 421 | 316 | 1.3 | 313 | 412 |
| 13N | 593 | 45,565 | 46,158 | 343 | 45,872 | 2,728 | 16.8 | ≈ 2,790 | 2,731 |
| `natural_neon` (3 nuclei) | 6,271 | 289,578 | 295,849 | 5,367 | 291,847 | 56,703 | 5.1 | ≈ 25,300 | 131,199 |
| `thick_target_resonance` | 56 | 48,237 | 48,293 | 56 | 48,239 | 2,149 | 22.4 | ≈ 2,200 | 15,036 |
| `energy_shift_sqrt` | 120 | 0 | 120 | 120 | 121 | 115 | 1.05 | 114 | 114 |
| `shift_fn` | 168 | 0 | 168 | 144 | 145 | 119 | 1.2 | 118 | 118 |
| `18O_p_a_thm` (THM) | 118 | 14,688 | 14,806 | 118 | 14,689 | 704 | 20.9 | ≈ 710 | 703 |

- *Solves* count `AMatrixFunc::FillMatrices` calls (one fill and one factorisation of every Jπ block of the compound
  nucleus) in the calculation of menu 1. The fit's evaluations do the same number (AZURE_PROFILE "matrix fill": T2
  24,299, 13N 45,871, T1 206, T3 414 per evaluation).
- *One shared mesh* is the node count of one grid per compound nucleus that has, at every energy, the largest node
  density of the sub-point grids covering it (∫ max_g 1/h_g(E) dE over the union of the windows), plus 2 nodes per
  folded point for the kernel's edges, plus the exact energies of the points without a fold. It is an estimate of a
  grid at today's accuracy, not a planned grid (inferred accuracy; measured node densities).
- *Distinct (E_x, entrance, exit)* is the number of T-matrices a shared grid needs (one per energy and reaction).
  `thick_target_resonance` has 7 capture exits, `natural_neon` three nuclei and several exits.

Where the projects without folds (T1, T3, the shift projects) gain little (1.05–1.5×: elastic and capture data at
the same energies in different data sets), the folded projects gain 4.5–22× from de-duplication alone, and the
shared mesh adds a factor 2.0 for T2 and 2.2 for `natural_neon`, ≈ 1 elsewhere (13N, the thick target and THM grids
already coincide on the anchored lattice).

### 1.3 Where an evaluation's time goes (measured, 1 thread)

Per evaluation, median over the evaluations after the first of a MIGRAD run (`AZURE_PROFILE=1`, menu 2,
`--max-calls=3`; stages are nested and thread-summed, ms):

| project | evaluation | point pass | level matrix fill + solve | T-matrix | cross sections | re-energised units (hand-over) |
|---|---|---|---|---|---|---|
| T1 | 1.0 | 1.0 | 0.3 | 0.2 | 0.5 | — |
| T2 | 92.6 | 69.5 | 26.9 | 11.9 | 6.8 | 15.4 ms, 5,822 units |
| T3 | 1.1 | 0.9 | 0.4 | 0.2 | 0.1 | — |
| 13N | 259.5 | 258.0 | 77.4 | 88.5 | 66.3 | — |
| `thick_target_resonance` (no free parameter: the second evaluation) | 800.1 | 794.8 | 473.4 | 169.5 | 121.8 | — |
| `energy_shift_sqrt` | 1.2 | 1.2 | 0.1 | 0.1 | 1.0 | 0.1 ms, 48 units |
| `18O_p_a_thm` (second evaluation) | 26.0 | 22.4 | 8.7 | (THM) | (THM) | — |

By function group (perf, 1 thread; the runs include the set-up, which dominates a single calculation):

| group | 13N fit, 40 calls | T2 fit, 40 calls | thick target, menu 1 | `natural_neon`, menu 1 | `18O_p_a_thm`, menu 1 |
|---|---|---|---|---|---|
| Coulomb and Whittaker functions (set-up: channel table, EC) | 20.6 % | 3.0 % | 0.5 % | 2.0 % | 1.4 % |
| level matrix: fill and solve | 27.4 % | 27.7 % | 31.5 % | 10.2 % | 6.6 % |
| T-matrix | 16.5 % | 5.8 % | 5.8 % | 12.0 % | — |
| observables (cross sections, angular sums) | 6.4 % | 12.4 % | 3.2 % | 4.0 % | 4.0 % |
| THM | — | 4.2 % | — | — | 8.0 % |
| folds and grids | 0.5 % | 5.8 % | 0.5 % | 0.4 % | 0.3 % |
| per-unit channel tables and E-dependent values | 0.9 % | 8.1 % | 22.0 % | 38.5 % | 20.7 % |
| memory: allocation, fills, unit copies | 11.1 % | 10.0 % | 15.4 % | 16.8 % | 40.1 % |
| accessors of CNuc / JGroup / ALevel / Unit | 9.3 % | 15.0 % | 10.3 % | 11.6 % | 9.1 % |
| other (incl. external capture 1.9 % in 13N) | 7.1 % | 7.6 % | 10.5 % | 4.4 % | 9.8 % |

Three things follow. (i) In a fit evaluation the level matrix and the T-matrix take 40–80 % of the point pass, and
both are repeated per unit; (ii) the folds themselves are cheap (≤ 6 %); (iii) in a single calculation of a large
folded project the per-unit set-up (filling 480-byte units and their tables, 296,000 of them in `natural_neon`) is the
largest cost, which a grid with nodes outside the units removes.

### 1.4 What is copied per thread (from the code and `HANDOFF.md`, 2026-11-11)

A fit evaluation works on a pooled copy of the `WorkingCopy` and of `CNuc` (`AZURECalc::GetPooledData`, a pool of
max(4, cores) copies, `fit_threads.md` §1). The copy holds the unit vector (480 bytes per unit) and the per-segment
state; the `DataSet` is shared. Measured on milestone 3: one copy after `Initialize` is 11.87 MB in 95 allocations for
T2 and 142.7 MB in 508 for `natural_neon` — which is the unit count times 480 B (24,628 × 480 B = 11.8 MB; 295,849 ×
480 B = 142 MB). The sub-points are 98 % of the units of both. Each thread also keeps a `thread_local AMatrixFunc` and
`THMMatrixFunc` (buffers sized by the largest Jπ block, small). The re-placed node tables are shared between copies
through `NodeTableCache` (its mutexes are one of the two pieces of shared mutable state left in the data path).

### 1.5 What the numbers mean for the design

1. The grid E_n should first be **the set of distinct compound energies of today's units**, with one solve each.
   That is bit-identical (the same solve on the same bits), needs no tolerance and no gate beyond today's, and gives
   most of the gain (§1.2). Phase 6's "parallel over grid-energy blocks" is the same change: a task is one energy and
   all units at it.
2. A shared mesh replaces the per-point sub-point grids later; it is justified by memory (the sub-point units go:
   142.7 → about 3 MB per copy in `natural_neon`, inferred) and by deleting the per-point grid machinery
   (`ReplaceSubPointGrids`, `NodeTableCache`, sub-point units, automatic tolerance, windows), not by speed.
3. The channel functions are already on an energy table; the node table of the grid is that table.

### 1.6 Commands (scratch scripts kept in `work/p5design/`)

```bash
# environment
source ~/miniforge3/etc/profile.d/conda.sh && conda activate <the upgrade environment>   # UP: the directory that holds the worktrees and work/
export AZURE2_COULOMB_TABLE=$UP/work/coulomb_table_v2.bin
S=$UP/work/p5design; T=<worktree>/tests

# instrumented build (never committed): work/p5design/p5_dump.patch on a copy of 2363abac, build work/b-p5design
#   P5_DUMP=1 prints every unit (P5U), its sub-points (P5SUB) and every level-matrix fill (P5S) to stderr
git archive 2363abac | tar -x -C $S/tree && (cd $S/tree && patch -p1 < ../p5_dump.patch)
cmake -S $S/tree -B work/b-p5design -DCMAKE_BUILD_TYPE=Release -DBUILD_GUI=OFF -DUSE_API=OFF
cmake --build work/b-p5design --parallel 2 --target AZURE2
bash $S/runproj.sh work/b-p5design/src/AZURE2 $T/13N 13N.azr dump-13N P5_DUMP=1      # menu 1, 1 thread
python3 $S/grid_count.py $S/runs/dump-13N                                            # §1.2 columns 2-9
python3 $S/reaction_count.py $S/runs/dump-13N                                        # §1.2 last column

# stage times of fit evaluations (Release build of the same tree, work/b-trial18)
MENU=2 ENGINE_FLAGS="--max-calls=3" bash $S/runproj.sh work/b-trial18/src/AZURE2 $T/13N 13N.azr prof-13N AZURE_PROFILE=1
python3 $S/prof_summary.py $S/runs/prof-13N

# perf by function group
CALLS=40 bash $S/perfproj.sh work/b-trial18/src/AZURE2 $T/13N 13N.azr perf-13N      # MENU=1 for a calculation
python3 $S/perf_groups.py $S/runs/perf-13N
```

The projects: `tests/upgrade/T1..T3`, `tests/13N`, `tests/natural_neon`, `tests/thick_target_resonance`,
`tests/energy_shift_sqrt`, `tests/shift_fn`, `tests/18O_p_a_thm`, each copied into `work/p5design/runs/<tag>/p`.

## 2. The grid

### 2.1 What it is

E_n, n = 1 … N, are compound-nucleus excitation energies (MeV), one list per compound nucleus of the project. A
grid energy carries: the channel energies E_c = E_n − S_α − E_x,α of every particle channel of that nucleus, their
channel functions (P, S, φ, σ_l: the rows of the per-stage channel table), and nothing that depends on a data point.
The collision matrix U(E_n) is solved once per Jπ block per E_n and shared by every reaction, angle and data set.

Exact energies (decision of 2026-11-07, `ARCHITECTURE.md` §1) are kept: the grid contains, bit for bit,

1. the energy of every data point without a fold (its E_x at the present shift);
2. the nodes of every fold (§2.2);
3. for a shifted data set without a fold, its shifted energies, which move with ν (§2.5).

### 2.2 How the fold nodes are planned

First (bit-identical step): the nodes are today's — each point's sub-point grid as `WorkingCopy::SetUpTargetEffects`
and `ReplaceSubPointGrids` place it, the grid being the union of all of them with exact duplicates merged. Nothing
about the quadrature changes.

Then (the shared mesh): one mesh per compound nucleus, planned from

- the **supports** of the kernels: for point i the interval [a_i, b_i] in E_x where its kernel is non-zero (target:
  from the back of the layer minus k σ_eff to the beam energy plus k σ_b; Gaussian: E_0 ± r σ; profile, UDR, Doppler:
  their stated windows; `TargetEffect::GridExtent` today);
- the **levels**: their quantized anchors (`AdaptiveIntegrationGrid::Quantize`: Γ to the nearest power of 1.25, E_λ
  to a multiple of Γ/4) and reach (core of m widths at Γ/p per node, then a geometric tail with step = t·distance up
  to the base step h_0);
- the **numerical settings** of each data set (`quadrature` in format 3: points per width p, width multiplier m, base
  step h_0).

The local step is h(E) = min over the data sets whose supports cover E of their own step h_s(E), and the nodes are
placed on the lattice of the anchors (core nodes at E_anchor + k Γ_q/p), so two data sets with the same settings get
the same nodes, and the mesh is a function of the quantized anchors, the supports and the settings only. A point's
kernel is evaluated on the mesh nodes inside its support; the support's edges are handled by clipping the interval
that contains them (as `GridQuadrature` clips [lo, hi] today), so they add no node.

**Status, branch 7 part 2 milestone 2 (`feature/upgrade-5-shared-mesh`, 2026-11-12): wired.** `WorkingCopy::
SetUpTargetEffects` queues every folded point of an unshifted, non-THM segment with an adaptive grid as a support of
its nucleus's mesh (`src/MeshSlots.cpp`: `QueueMeshPoint`, support [a, b] in E_x = E_cm + S + E_x,pair, settings from
the point's `Unit::SubGridSpec`, a hard lower edge at the integration floor or a thick target's 1 keV), then plans each
nucleus once (`PlanMesh`, anchors `MeshAnchorsOf`: every level as each entrance pair of the nucleus's folded points
sees it, quantized in E_x, from the input widths) and gives each point the nodes of its range (`MeshSlotEnergies`, c.m.
energies descending).  The kernels integrate over the support (`KernelNodes::Low()/High()`), not over the span of the
nodes, whose two outer nodes of a soft edge lie up to one interval outside it.  Off the mesh, with grids of their own
as before: shifted segments (branch 8), THM (branch 10), `useAdaptiveGrid` off, and the σ(Ē) fallback nodes and windows
(milestone 3).

**Status, milestone 3 (`feature/upgrade-5-d3-deletions`, 2026-11-12): the four mechanisms of D3 are deleted.**  The
automatic tolerance, the fold windows (explicit and `auto`), the edge blend and the σ(Ē) fallback are gone from the
engine (`TargetSetup.cpp`, `PointPass.cpp`, `Unit.h`, `TargetEffect.*`, the hand-over, the stopping scale, the Coulomb
reference, the gradient), the target kernel (`KernelTarget.cpp`), format 3 (`fold_windows` is refused with
`E_FIELD_DELETED`), the converter (an old line's ranges, blend and tolerance are dropped with a note) and the docs.
Every point of a folded data set is folded with the whole kernel.  What is left off the mesh, each with a grid of its
own and a sharing group of its own (`SlotGroup.h`): a segment with an energy shift (its sub-points move rigidly with the
shift until branch 8 folds it on fixed nodes), THM (its kinematics on the mesh are branch 10), and every folded point
under `useAdaptiveGrid` off (a uniform grid of `points` intervals, the old option).  A thin target narrower than one
mesh interval is an ordinary folded point: the rule of §2.4 clipped at its support integrates the linear interpolant of
the two enclosing nodes over [E0 − Δ, E0], i.e. the interpolated σ at the depth-averaged energy (exactly so for 1/ε
constant across the target) — the linear-interpolation case of the deleted fallback, with the nodes of the mesh, whose
interval near a level of width Γ is at most Γ/p.  Moved: `tests/target_effect_ranges` (the project of the tree that
used them), χ² 3838.66 → 3760.95, the beam spread now folded outside its old window (largest A_y change 0.046, 1.5 data
σ); its numbers equal the engine before the deletion run without the windows, bit for bit.

### 2.3 Fixed per stage, and what happens when a narrow level moves

Today every evaluation re-places the grids (`ReplaceSubPointGrids`) from the quantized anchors; the grid a copy
holds is a function of θ alone since the exact-comparison fix (`fit_threads.md` §2), and it is piecewise constant in
θ: it changes when a level moves by about Γ/4 or its width by 25 %. In T2's fit this re-energises 5,822 units per
evaluation (15.4 of 92.6 ms; shifts included).

The grid keeps that rule (decision D1, §9): **planned at the start of a fit stage; within the stage the nodes in the
reach of a level whose quantized anchor changed are re-placed in the hand-over**, and only those. The re-placed nodes
are a function of the anchors, so every copy at the same θ has the same mesh. Their channel functions are computed in
the hand-over into the working copy's own node table (from the Coulomb table; single-threaded or over the new nodes,
never during the pass). Kernel rows of the points whose support meets a re-placed region are rebuilt there too.
Between stages the whole mesh is re-planned.

The alternative — a mesh fixed through the stage, wide enough to cover each narrow level's motion — removes the
re-placement but makes a level that moves further than planned lose resolution silently, or forces a re-plan that is
a discontinuity of χ² inside MIGRAD. It is the PI's choice (D1); I recommend keeping today's rule.

**Status (milestone 2).** The hand-over (`WorkingCopy::ReplaceSubPointGrids`) takes the anchors of the observed
parameters of the current values; where they differ bit for bit from the copy's (`WorkingCopy::meshAnchors_`), it plans
the nucleus's mesh again and re-places the slots of the points whose support or present slots meet
`MeshChangedRegion(before, after)`, through the same overflow path as before; the others keep theirs (with
`AZURE_PLAN_CHECK=1` it checks that they would get the same bits, else it throws).  The plan itself is a whole
`PlanMesh` of the nucleus (0.6 ms on T2).

### 2.4 Quadrature rule and error bound

The rule is today's (`src/TargetFold.cpp`, `GridQuadrature`): between adjacent nodes E_j > E_{j+1} the observable is
interpolated linearly, σ_lin(E) = (1 − t) σ_j + t σ_{j+1}; the product K(E) σ_lin(E) (and 1/ε, interpolated
linearly, for a target) is integrated by two-point Gauss–Legendre on each interval, clipped to the support. The
error of one interval of width h is bounded by

```
|∫ K (σ − σ_lin) dE| ≤ (h³/12) · max|K| · max|σ″|          (linear interpolation)
|GL₂ remainder|      ≤ (h⁵/4320) · max|(K σ_lin)⁗|          (exact for cubics)
```

so for a Breit–Wigner peak of width Γ resolved with h = Γ/p the relative error of the peak region is of order 1/p²
(p = 50: 4 × 10⁻⁴), and the tails are controlled by the geometric step (`AdaptiveIntegrationGrid::TailFactor`: area
bias ≈ r² Γ/(2π core)). On the shared mesh, h(E) ≤ h_s(E) for every data set s at every E, so **the bound is never
worse than today's for any point**; the error itself differs (the nodes differ), which is what the gate classifies
(§7). Measured size of today's own quadrature error, for scale: `thick_target_resonance` has the same data at width
multiplier 5 and 20; the yields differ by 2.0 × 10⁻³ (relative) at 253.8 keV below the resonance and by 0.6–0.9 ×
10⁻⁴ on the plateau (`AZUREOut_aa=1_TOTAL_CAPTURE.out`, segments 1 and 2).

### 2.5 Shifted data on the grid

- **Without a fold** the point is an exact evaluation energy: its node moves with ν (as today,
  `Segment::UpdatePointEnergyWithShift`), its channel functions are recomputed in the hand-over, and only its own
  node. ∂m/∂ν = (dO/dE) ∂E_0/∂ν, with dO/dE from the channel-function derivatives (S′, P′) or a difference on that
  node alone.
- **With a fold** the nodes stay; the kernel moves: K_in(ν) = K(E_n; E_0(ν)). The quadrature is the same rule on
  fixed nodes, the χ² is smooth in ν (no re-placed node, no table recomputed), and ∂K/∂ν is analytic (§3). Today the
  sub-points move rigidly with the point (`ShiftFrames.cpp`, since 2026-10-26), which in exact arithmetic is the same
  integral; the numbers move by the quadrature difference only.
- The **interpolated mode** of `ARCHITECTURE.md` §1 (observables on a fixed mesh, an unfolded shifted point as
  interpolation weights) stays opt-in, with its check against the exact mode.

**Design of the shifted fold on the fixed nodes (branch 8 `shift-kernel`, 2026-11-12; the full text is the comment of
`include/ShiftedFold.h`, which owns it).**  Δ_i(ν) = E_cm,i(ν) − E_cm,i(0) is the c.m. shift of point i; the kernel is
translated rigidly, K_i(E; ν) = K_i(E − Δ_i; 0), which is the engine's fold before in exact arithmetic.

1. *The planner's support.*  The unshifted support [a_i, b_i] moved by a centre c_i and widened by a reach R_i on each
   side (never below the integration floor, a hard edge there).  R_i = 3 Σ_k |Δ_i(ν_0 + σ_k e_k) − Δ_i(ν_0)| over the
   varied shift parameters k (a follower: its leader's; a shift-function group: the owner's), σ_k the prior width; a
   varied parameter without a prior adds b_i − a_i; R_i = 0 when nothing of the point's shift varies (a fixed shift
   moves the support once, no node added).  From what the project already says (the priors), no new field.
2. *Leaving the reach.*  c_i = Δ_i(ν_0) + 2R_i round((Δ_i(ν) − Δ_i(ν_0))/(2R_i)) (R_i = 0: c_i = Δ_i(ν)), a function of ν
   alone, so the moved support stays inside the planned one.  A change of some c_i re-plans the nucleus in the
   hand-over with the copy's supports and re-places the groups whose nodes changed, the path of a moved narrow level
   (§2.3); deterministic at any thread count, counted in `run.log`; χ² steps by the quadrature difference there.
3. *The sites a point reads.*  Its slots are its planned range (fixed for the stage); the fold reads the tightest
   sub-range enclosing the moved support [a_i + Δ_i, b_i + Δ_i] (at least two nodes), the kernel clips there; the block
   is not rewritten.  What depends on ν besides the kernel's centre: the stopping the target kernel reads, ε(max(E_n − Δ_i, a_i))
   (the fold before kept each sub-point's stopping at its unshifted energy: the depth profile and T are the point's as
   measured), computed in the gather from the target's stopping — the point's, not the site's; and the c.m. angle of a
   lab-frame differential point (converted at its shifted energy), which its group's sites carry: their Legendre rows
   and Coulomb amplitude are recomputed by the slot's writer when it changes (the channel rows stay), and the group key
   holds the point's nominal lab energy so that every reader of a block has the same angle at every ν.
   (The argument is held at the nominal support's lower edge a_i: an outer node below the moved support only shapes
   the edge interval's interpolant, and at a thick target's floor E_n − Δ_i < 0 is outside the stopping's domain.)
   A segment with external capture whose shift leaves the reach is re-placed like any other: its new sites take their
   capture amplitudes from the tables at their energies (`Site::GetECAmplitudeWithShift`).
4. *Writers.*  Ranges, blocks, writers, plan tasks and the energy pass's task list are those of the planned ranges, fixed
   for the stage; every site of a block is computed (the margin is the cost).
5. *Gradient.*  A shift column of a folded segment is the central difference of the fold on the fixed nodes: at ν ± h
   the points move and are folded again from the node values of ν (the same `FoldValues`), no node recomputed, no level
   matrix solved; shift kinds enter through Δ_i only.  The analytic ∂K/∂E_0 is open for every kernel (it needs the
   kernels' energy arithmetic on a dual type).  A lab-frame differential segment's column still recomputes its sites'
   observables at the moved angle (open for the observables per node, branch 9).

**Measured (branch 8, build `work/b-shift`, 1 thread unless stated).**
- *Moves* (classifier, three bases, against the engine before refined eightfold; T2's `covariance.out` is a named row
  since the director's review): every output file of T2, T2_true and
  `shift_sqrt_subpoints` (const, sqrt, shiftfn, narrow_const, narrow_sqrt) is per-file quadrature and closer to the
  refined value than before (T2 `AZUREOut`: largest 1.1e-3 against 5.4e-3, rms 1.7e-4 against 5.5e-4; T2 χ² of the
  thin-target segment 692.29 → 697.33, refined 697.67; thick target 283.21 → 286.44, refined 286.53), **except T2's
  `covariance.out`**: the correlated χ² of the two lab-frame segments, 513.884 → 513.904 against refined 513.771
  (2.6e-4 against 2.2e-4, 1.15 × the old distance; each of the two segments' own χ² is closer than before).  It is the
  same with every shift fixed (reach 0): a matter of the mesh's nodes for those segments, not of the reach.
- *Cost of the reach in nodes*: T2 3919 energy tasks (distinct nodes solved; before 5157; with every reach 0: 3767) and
  23602 slots (18172 with reach 0); `shift_sqrt_subpoints/sqrt` (fixed shift) 2057 tasks against 6738 before.
- *Step at a re-plan* (the director, `work/director/shift_review/steps.py`, 500-600 values of each of T2's three shift
  parameters): one step, at the re-plan boundary of the offset with prior 0.003 (at +9.0e-3 MeV), 1.3e-2 in a χ² of
  10,500; none inside the reach.  The engine before had steps of the same size elsewhere (1.7e-2 at 1.16e-3 MeV on the
  same parameter, 1e-3 on the sqrt shift at three places).
- *Smoothness* (`pyazr_shift_fold_smooth`): 21 shifts across 0.2 keV lab of a thin layer whose back edge lies in a
  20 eV core lattice: second differences −0.79 … −0.69, third differences 0.002 … 0.009; without the clip at the moved
  support, third differences ±700.
- *Gradient* (one gradient, `AZURE_PROFILE`): T2 158,929 → 48,526 level-matrix solves, 49,828 → 640 channel-function
  sets; `shift_sqrt_subpoints/sqrt` with a varied sqrt shift 26,952 → 8,276 solves (the remainder is the adjoint's
  forward pass), 141,582 → 84 channel-function sets, 2.5 s → 0.6 s.

### 2.6 What replaces the node-table cache

The node table of a stage — per node its channel-table rows, per (point, node) its angle rows (Legendre polynomials,
Coulomb amplitudes, attenuation) — is built once at the stage start and held as `shared_ptr<const>` by every copy,
like the `DataSet`. Nodes re-placed inside a stage (§2.3) live in the copy's own small overflow table, computed in
its hand-over. No copy publishes anything to another; `NodeTableCache` and its mutexes go, and so do the per-unit
copy-on-write tables of sub-points.

**Status, branch 7 part 1 (2026-11-11).** The node table exists beside the sub-point units (`include/NodeTable.h`,
`src/NodeTable.cpp`): one node per (nucleus, route, bits of E_x) of the units of the point pass, with its channel rows
(a holder of the rows the first unit at the node computed: no second copy), built at the end of
`WorkingCopy::Initialize`, held `shared_ptr<const>` by every copy.

**Status, branch 7 part 2 (`feature/upgrade-5-grid-readers`, 2026-11-12).** Every slot is on a table: each plan task
has one slot per unit it calculates (its point, or its sub-points), and each slot names a node of the stage table or of
the working copy's **overflow table** (`NodeOverflow`), which the hand-over places from nothing for the tasks whose
units moved (re-placed grids, `ReplaceSubPointGrids`; shifted points, `Segment::UpdatePointEnergyWithShift`) and the
tasks already on it (`WorkingCopy::PlaceOverflowNodes`, end of `FillEnergyShiftsFromParamsImpl`).  The overflow table
is a function of the vector handed over alone; its rows are let go before units are re-energised, so the
copy-on-write of a unit's rows is not triggered by it.  `AZURE_NODE_CHECK=1` checks every unit of every task at every
hand-over and every point pass: 0 off the tables, 0 on a wrong node (test `upgrade_node_check`).  The seam for the
shift-kernel branch is the comment of `NodeOverflow`: what branch 8 deletes there (the rigid move of sub-points and
their overflow entries).  First reader moved: the energy tasks are grouped by the slots' nodes (nucleus, route, E_x
from the node tables, `CalculationPlan::EnergyTasks`; `EnergyTask::node`), no longer by recomputing E_x per unit.
Not moved (the sub-point units still carry them): the observables the energy pass writes (unit records), the kernels'
gather (`PointFold.cpp`: c.m. energy and stopping per sub-point), the result records, the gradient's unit cotangents,
the external-capture amplitudes and `intEC.dat`, THM `thm_ps_`; not deleted: sub-point units, `SubPointAccess.h`,
`NodeTableCache.h`, the copy-on-write tables.  What a slot must carry for them (the c.m. energy of the entrance pair,
which is not E_x minus a constant to the last bit; the stopping, which a fitted stopping scale moves per copy; angle
rows, Coulomb amplitude, screening, EC amplitudes) is listed in `HANDOFF.md`.

**Status, branch 7 part 3 (`feature/upgrade-5-grid-readers2`, 2026-11-12).** The sub-point units are gone: a slot of a
folded point is a `Site` (`include/Site.h`: what depends on the point at one energy), held in one vector of the working
copy (`WorkingCopy::SlotAt`, each point's a range); a data point (`Unit`) is a site with its row, fold, plan task and
model record. The observables, the energy pass, the folds, the gradient's cotangents and the external capture read
sites by one index (`WorkingCopy::SiteAt`: points, then slots). `SubPointAccess.h` is deleted. Not deleted:
`NodeTableCache.h` (re-placed slots share their rows across copies through it until the shared mesh stops nodes from
moving with the levels) and the slots' copy-on-write rows (the node rows are holders of them).

**Status (milestone 2).** `NodeTableCache` is kept: the shifted folded points still have grids of their own that the
hand-over re-places from the shifted anchors at every move of a shift, and mesh slots re-placed after an anchor change
are the same nodes in every copy at the same parameters; the cache lets a copy adopt the rows another copy computed for
them (on T2 most re-placed slots take another copy's rows).  It goes with branch 8, when no node moves with a shift.

**Status, branch 8 (`feature/upgrade-5-shift-kernel`, 2026-11-12): deleted.**  `NodeTableCache.h`, its mutexes and
`SharedTable::Weak` are gone; a copy that re-places a block after an anchor change computes the rows of its new nodes
itself (only the nodes near a moved narrow level move on the shared mesh).  `NodeOverflow` and `PlaceOverflowNodes`
stay: their users are the re-placed blocks of a moved narrow level and THM's sub-points, which still move with a
shift (`Segment::MoveOwnGridSlots`, branch 10).

**Status, shared sites (`feature/upgrade-5-shared-sites`, 2026-11-12).** The folded points of a data set that need the
cross section at one mesh node with everything else equal share one `Site` there (`include/SlotGroup.h`: the key --
segment and term, c.m. angle bits, Legendre order of a coefficient row, target effect, stopping factor, frame shift,
screening -- and the audit behind it).  A group's sites are one block of `WorkingCopy::SlotAt` in node order; a point's
slots are a range of it, ranges overlap.  A slot is written by its first reader in data order
(`WorkingCopy::WritesSlot`: set-up values, energy pass, record, external capture, stopping scale, shift step), read by
all (folds, model pass, gradient; a site's cotangent is the sum of its readers' in point order).  Folds keep E1/E2,
A_y sigma and a_k A_0 in arrays (`FoldValues`), never in a slot's record.  The hand-over re-places a group at once.
Off the mesh a point is a group of its own.  Sites: natural_neon 373,449 -> 47,649, thick_target_resonance 74,634 ->
28,924, stopping_yield/models 179,438 -> 4,221.  `NodeTableCache` is still needed (re-placed blocks adopt other copies' rows).

### 2.7 Thread-count independence by construction

- The mesh and the node table are built single-threaded (or by a parallel loop that writes one node per index) from
  immutable inputs: the measurements, the quantized anchors, the settings. They are the same at every thread count.
- The evaluation's tasks are grid energies (or blocks of them). A task solves U(E_n) and writes only the observables
  O_n of that energy, into the slots the plan gives it; no task reads another's slots.
- A kernel row is applied by one task per data point, summing over its nodes in node order; the χ² is reduced in
  data order. No sum depends on which thread did a block or in which order blocks finished.
- Nothing is cached across evaluations beyond the stage table and the copy's own overflow table, both functions of
  (θ, ν) alone.

### 2.8 The planner of the shared mesh (branch 7 part 2, first milestone; 2026-11-12, `feature/upgrade-5-mesh-planner`)

`include/MeshPlanner.h`, `src/MeshPlanner.cpp`: `PlanMesh(anchors, supports, tails)`, a function of values (no engine
object): quantized anchors in E_x, one support [a_i, b_i] in E_x per folded point with its data set's settings (smooth
step, points per width p, core widths m), the tail constants of `AdaptiveIntegrationGrid`; out the nodes (ascending,
distinct) and per support the range of nodes its quadrature uses (the outer nodes of the two edge intervals: the edges
add no node). Not wired into the engine. The header's comment states the rules; in short:

- **Smooth nodes**: today's smooth step S = window / ceil(window / step) quantized down to a power of two,
  L = 2^floor(log2 S) MeV; the nodes are n·L (exact in floating point, a lattice of the nucleus with origin E_x = 0).
  The ladders are nested, so supports of different steps share nodes and every change of step is on a node of both.
- **Core nodes**: E_r ± k·Γ/p, k ≤ ceil(m p − 10⁻⁹) (today's nodes, bit for bit, when m p is an integer).
- **Tail nodes**: today's recurrence o_{k+1} = o_k + max(Γ/p, t o_k) from the core's edge, a function of the anchor and
  (p, m) alone, run while the step is below twice the class's largest ladder step.
- **Overlaps**: each lattice is a list of intervals; an interval is kept (both ends become nodes) unless one other
  lattice, active over all of it, is finer over all of it — finer than its length between two anchors' lattices or two
  ladders, finer than half its length between a ladder and an anchor's lattice. The half reproduces today's band where
  a point's grid keeps its smooth nodes beside a tail; without it the mesh never exceeded today's step either, but had
  up to 32 % fewer nodes than a point's own grid in 170 supports of `natural_neon`.

Why the step is never coarser than today's: at any E of a support, the lattice F that is finest there has its interval
around E kept (anything dominating it would be finer than F at E), and that interval is no longer than F's step at E,
which is ≤ today's step D(E) of every support covering E (L ≤ S, the same lattices). `MeshChangedRegion(before, after,
supports, tails)` returns the stretches outside which two anchor lists give the same nodes bit for bit: each changed
anchor's reach plus twice the largest ladder step.

**Tests** (`tests/upgrade/mesh_planner_test.cpp`, ctest `upgrade_mesh_planner`, 700,675 checks): against
`AdaptiveIntegrationGrid::GenerateGrid` on an isolated narrow level, two overlapping levels of different widths, a level
at a support's edge, a level outside in reach, no level, a 2-width convolution core, a zero-width support, and 200
random cases (seed fixed; one or mixed classes of settings): every mesh interval inside a support is no longer than
today's step D anywhere on it (D transcribed from `GenerateGrid`, and checked against its real grids: no real interval
exceeds 2.5 D); with one class of settings every node of today's lattices finer than the mesh's smooth step is a mesh
node, bit for bit; the range of a support is the tightest enclosing one. Shared: two overlapping data sets with the same
settings read the same nodes in the overlap; 7 identical supports plan the nodes and range of one. Order: 40 random
cases × 3 shuffles of supports, data-set numbers and anchors give the same bits; two disjoint groups planned together
and apart give the same nodes. Local: 100 random anchor changes (moved by quarter widths, two notches wider, moved
3 keV and narrower): the nodes outside `MeshChangedRegion` are the same bits. Shown failing with the smooth step one
notch too coarse on the ladder (122,979 failures).

Point by point the mesh interval can still be longer than the point's own interval today (on 1.7–3 % of a support's
length in the named cases): where today's grid interleaves its smooth lattice with a tail, its intervals are shorter
than its own step, at phases a shared lattice cannot reproduce. The bound of §2.4 is not affected (it is in terms of the
step), and the fold check below measures the error itself.

**Node counts** (`tests/upgrade/node_count --plan`, every project of the tree that sets up and has a fold; supports from
the set-up's sub-points, settings from `Unit::SubGridSpec`, anchors of the first folded point's entrance pair moved to
E_x and quantized there, taken before the widths are transformed, as the set-up saw them). Columns: units of the pass;
distinct E_x today; §1.2's estimate; planned mesh (planned nodes plus the exact energies of unfolded points); planned /
today; over-coverage (planned nodes with no node of today's grids, of any point, within the mesh's own spacing there);
nodes outside every support (outer nodes of edge intervals); cost of the ladder alone (∫ max 1/L over ∫ max 1/S on the
union of the supports); supports in which the mesh has fewer nodes than the point's own grid (inside plus the two outer
nodes); the smallest such ratio.

| project | units | distinct E_x today | estimate (§1.2) | planned | planned / today | over-coverage | outside | ladder cost | supports with fewer nodes | smallest ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| `13N` | 46121 | 2727 | 2786 | 2696 | 0.99 | 2 | 2 | 1.27 | 0 | 1.00 |
| `17O` | 19560 | 2021 | 1905 | 2000 | 0.99 | 95 | 2 | 1.35 | 0 | 1.00 |
| `18O_p_a_thm` | 14688 | 703 | 714 | 657 | 0.93 | 164 | 3 | 1.74 | 0 | 1.61 |
| `6Li_d` | 3162 | 3162 | 1182 | 1221 | 0.39 | 0 | 2 | 1.15 | 0 | 1.16 |
| `7Li_p_a` | 6807 | 4761 | 3415 | 3577 | 0.75 | 738 | 3 | 1.34 | 0 | 1.03 |
| `beam_profile_kernel` | 82264 | 3086 | 3161 | 3199 | 1.04 | 190 | 7 | 1.61 | 0 | 1.00 |
| `hybrid_potential` | 111976 | 2471 | 2650 | 2472 | 1.00 | 2 | 2 | 1.62 | 0 | 1.00 |
| `layer_yield/models/fit` | 60367 | 3133 | 2175 | 2112 | 0.67 | 7 | 2 | 1.47 | 0 | 1.00 |
| `layer_yield` | 87634 | 4137 | 2762 | 2851 | 0.69 | 9 | 4 | 1.61 | 0 | 1.00 |
| `natural_neon` | 292742 | 56696 | 25313 | 21070 | 0.37 | 699 | 17 | 1.47 | 8 | 0.99 |
| `rectangular_resolution/models/fit` | 29144 | 3189 | 2280 | 2207 | 0.69 | 4 | 2 | 1.35 | 5 | 1.00 |
| `shift_sqrt_subpoints/models/const` | 36841 | 2089 | 2125 | 2059 | 0.99 | 25 | 4 | 1.23 | 0 | 1.00 |
| `shift_sqrt_subpoints/models/dense` | 56028 | 8004 | 2961 | 4576 | 0.57 | 0 | 2 | 1.57 | 0 | 1.56 |
| `shift_sqrt_subpoints/models/narrow_const` | 44212 | 2158 | 2192 | 2127 | 0.99 | 27 | 4 | 1.23 | 0 | 1.00 |
| `shift_sqrt_subpoints/models/narrow_sqrt` | 44212 | 2158 | 2192 | 2127 | 0.99 | 27 | 4 | 1.23 | 0 | 1.00 |
| `shift_sqrt_subpoints/models/noshift` | 36841 | 2089 | 2125 | 2059 | 0.99 | 25 | 4 | 1.23 | 0 | 1.00 |
| `shift_sqrt_subpoints/models/shiftfn` | 36841 | 2089 | 2125 | 2059 | 0.99 | 25 | 4 | 1.23 | 0 | 1.00 |
| `shift_sqrt_subpoints/models/sqrt` | 36841 | 2089 | 2125 | 2059 | 0.99 | 25 | 4 | 1.23 | 0 | 1.00 |
| `stopping_yield/models/fit` | 176270 | 5014 | 2331 | 2130 | 0.42 | 3 | 2 | 1.34 | 0 | 1.01 |
| `stopping_yield` | 199089 | 7295 | 2552 | 2175 | 0.30 | 3 | 2 | 1.60 | 0 | 1.01 |
| `target_effect_ranges` | 5686 | 2164 | 1191 | 541 | 0.25 | 117 | 4 | 1.34 | 0 | 1.15 |
| `thick_target_resonance` | 48237 | 2148 | 2197 | 2068 | 0.96 | 19 | 4 | 1.23 | 0 | 1.01 |
| `units_scales/models/yield_story` | 37372 | 1782 | 1830 | 1751 | 0.98 | 2 | 2 | 1.57 | 0 | 1.00 |
| `upgrade/T2/T2` | 24319 | 5442 | 2782 | 2917 | 0.54 | 513 | 3 | 1.80 | 0 | 1.36 |
| `upgrade/T2/T2_true` | 23281 | 5359 | 2792 | 2905 | 0.54 | 495 | 3 | 1.79 | 0 | 1.36 |

No project's mesh exceeds 1.2 × today's distinct energies (largest 1.04, `beam_profile_kernel`: its grids already
coincide on the anchored lattice, and the ladder adds 1.61 × the smooth nodes there). The ladder costs 1.15–1.80 × the
smooth nodes; the mesh still saves where data sets overlap (T2 0.54, `natural_neon` 0.37, the stopping and layer
projects 0.30–0.69, `target_effect_ranges` 0.25, whose fallback and window points have uniform grids today). The
supports with fewer nodes (`natural_neon` 8, `rectangular_resolution` 5) have ratio ≥ 0.99: one node, the accounting of
today's edge nodes against the mesh's outer nodes. The shift projects are planned at their unshifted anchors (the
set-up's grids of a shifted data set are offset by its shift, `GridConfig::resonanceOffset`): shifted folded points
keep their own nodes until branch 8. The over-coverage of 7Li_p_a (738), `natural_neon` (699) and T2 (513) is the
ladder's finer smooth step where the windows of today's points were coarser.

**The fold check** (`tests/upgrade/mesh_fold_check`, no engine): a Breit–Wigner (one level; or two, the second 3 widths
above, twice as wide, amplitude 0.7 e^{2i}) folded with a Gaussian (σ = 10 keV, support ±5σ, 100 sub-points, p = 50,
m = 2) and with a rectangle [E0 − σ, E0] (50 sub-points, m = 20), 13 points E0 = E_r + k max(Γ, σ)/2, by the rule of
§2.4 on (a) each point's grid today and (b) the mesh planned for the 13 points together; reference: the same rule on
the union of both, every interval cut in 64. Largest relative error over the points; nodes distinct over the points.

| kernel | levels | Γ/σ | error today | error mesh | nodes today | nodes mesh |
|---|---|---|---|---|---|---|
| Gaussian | one | 0.01 | 2.04e-04 | 1.90e-04 | 610 | 562 |
| Gaussian | one | 0.03 | 2.61e-04 | 2.29e-04 | 564 | 516 |
| Gaussian | one | 0.1 | 3.87e-04 | 3.35e-04 | 518 | 472 |
| Gaussian | one | 0.3 | 4.99e-04 | 3.75e-04 | 472 | 426 |
| Gaussian | one | 1 | 3.50e-04 | 2.67e-04 | 392 | 362 |
| Gaussian | one | 3 | 8.43e-05 | 8.08e-05 | 690 | 490 |
| Gaussian | one | 10 | 5.62e-05 | 5.36e-05 | 874 | 719 |
| Gaussian | one | 30 | 7.25e-06 | 6.91e-06 | 1313 | 1355 |
| Gaussian | one | 100 | 6.65e-07 | 6.34e-07 | 1313 | 1365 |
| Gaussian | two | 0.01 | 9.30e-05 | 7.57e-05 | 715 | 749 |
| Gaussian | two | 0.03 | 1.35e-04 | 1.11e-04 | 674 | 660 |
| Gaussian | two | 0.1 | 2.45e-04 | 2.31e-04 | 624 | 582 |
| Gaussian | two | 0.3 | 3.47e-04 | 2.55e-04 | 566 | 533 |
| Gaussian | two | 1 | 3.85e-04 | 2.78e-04 | 482 | 496 |
| Gaussian | two | 3 | 1.04e-04 | 8.84e-05 | 895 | 803 |
| Gaussian | two | 10 | 5.13e-05 | 4.89e-05 | 1308 | 1333 |
| Gaussian | two | 30 | 6.63e-06 | 6.33e-06 | 1313 | 1356 |
| Gaussian | two | 100 | 6.09e-07 | 5.81e-07 | 1313 | 1365 |
| rectangle | one | 0.01 | 3.85e-04 | 1.43e-04 | 2402 | 2543 |
| rectangle | one | 0.03 | 1.62e-04 | 6.21e-05 | 2336 | 2471 |
| rectangle | one | 0.1 | 3.51e-05 | 1.31e-05 | 2172 | 2253 |
| rectangle | one | 0.3 | 3.05e-05 | 3.06e-05 | 1174 | 1160 |
| rectangle | one | 1 | 5.62e-05 | 3.20e-05 | 816 | 956 |
| rectangle | one | 3 | 1.61e-05 | 6.00e-06 | 663 | 1091 |
| rectangle | one | 10 | 2.50e-06 | 9.31e-07 | 663 | 1089 |
| rectangle | one | 30 | 2.94e-07 | 1.10e-07 | 663 | 1089 |
| rectangle | one | 100 | 2.66e-08 | 9.93e-09 | 663 | 1079 |
| rectangle | two | 0.01 | 4.26e-04 | 1.58e-04 | 3344 | 3489 |
| rectangle | two | 0.03 | 6.42e-05 | 3.16e-05 | 3235 | 3350 |
| rectangle | two | 0.1 | 2.97e-05 | 2.64e-05 | 2714 | 2717 |
| rectangle | two | 0.3 | 3.21e-05 | 3.24e-05 | 1174 | 1160 |
| rectangle | two | 1 | 4.22e-05 | 2.40e-05 | 1356 | 1774 |
| rectangle | two | 3 | 1.37e-05 | 5.08e-06 | 663 | 1090 |
| rectangle | two | 10 | 2.16e-06 | 8.05e-07 | 663 | 1089 |
| rectangle | two | 30 | 2.65e-07 | 9.86e-08 | 663 | 1089 |
| rectangle | two | 100 | 2.43e-08 | 9.05e-09 | 663 | 1079 |

The mesh's error is not larger than today's in 34 of 36 cases. The two others (rectangle, Γ/σ = 0.3, one and two
levels: 3.06 against 3.05 × 10⁻⁵ and 3.24 against 3.21 × 10⁻⁵) come from the edges: there the core lattice covers the
whole window, so both grids have the same interior nodes, and today evaluates σ at the support's edges while the mesh
interpolates the clipped edge intervals from the outer nodes; with the edges added to the mesh (`mesh_fold_check
--with-edges`) the errors are equal to the last digit.

**What the engine will hand the planner** (milestone 2): the supports from `WorkingCopy::SetUpTargetEffects`
(`src/TargetSetup.cpp`: `azure2::experiment::EffectSupport` of the effect's kernels, c.m. of the entrance pair, moved to
E_x as `AMatrixFunc::InEnergy` moves a unit's energy, THM by its kinematics); the settings from the target effect
(`NumSubPoints`, `GetResonanceWidthMultiplier`, `GetPointsPerWidth`, times `refineQuadrature`; recorded today per point
in `Unit::SubGridSpec`, `MakeSubGridSpec`); the anchors per nucleus from `AdaptiveIntegrationGrid::Anchors` with the frame
shift −(S + E_x) of a pair of that nucleus, at set-up from the input widths and in the hand-over from the observed
parameters of the current values (as `WorkingCopy::ReplaceSubPointGrids` builds them); the tail constants. The hand-over
needs from locality: `MeshChangedRegion` of the anchors before and after, the nodes inside it re-placed into the
overflow table, the kernel rows of the points whose support meets it rebuilt.

**Open, for the director or the PI.**
- The support's edges add no node: the edge interval is interpolated from a node outside the support. For a support that
  ends at a physical floor (a thick target at 1 keV, THM data below the entrance threshold) the outer node lies beyond
  the floor; milestone 2 needs either the edge as an exact node there or a rule for that node. It also moves the error of
  edge intervals by about 1 % of itself (above).
- Anchors are quantized in E_x, not in each entrance pair's c.m. frame as today: one lattice per nucleus. Numbers move by
  the quantization (class "quadrature" of §7.1).
- The ladder of powers of two costs up to 1.8 × the smooth nodes. A finer ladder (two nested families, 2^k and 3·2^k)
  would halve that, at the price of transitions between families; not needed while the mesh stays ≤ 1.04 × today.
- `MeshChangedRegion` widens by twice the largest ladder step of the nucleus; with a thick target whose smooth step is
  tens of keV the region re-placed after an anchor change is that much wider than the anchor's own reach.

## 3. Observables on the grid and the kernels

### 3.1 The three lines, concretely

```
O_n,r,q(θ)   = observable q of reaction r = (entrance α, exit α′) at grid energy E_n
m_i(θ, ν)    = c_i(ν) + Σ_n Σ_q K_i,(n,q)(ν) · O_n,r(i),q(θ)        (affine: c_i is |f_C|² for elastic data, else 0)
data term    = (m_i/J_i(ν) − n_k y_i^lab) / (n_k δy_i^lab)           (lab data, §4; c.m. data: J_i = 1)
```

The expression nodes of Phase 4 (`Sum`, `Ratio`, `Scale`, `Weighted`, `Transform`) stay above the kernels: they
combine the folded leaves of one point and are not energy kernels.

### 3.2 Observables (one value per grid energy and reaction)

| observable O_n | formula (energies c.m., σ in b) | replaces |
|---|---|---|
| integrated σ_r | (π/k_α²) Σ_J g_J Σ_{c∈α, c′∈α′} \|T^J_{cc′}\|², T from U and, for capture, the internal and external amplitudes | `GenMatrixFunc::CalculateBareCrossSection` (integrated branch), per unit |
| Legendre coefficients B_L,r | dσ/dΩ(θ) = Σ_{L≤L_max} B_L P_L(cos θ) (+ the Coulomb terms below for elastic); B_L from the pathway table's coupling coefficients | the angular sums of `CalculateBareCrossSection` (KL groups, interferences) and the per-unit `Unit::CalcLegendreP` (P_L moves to the kernel row) |
| nuclear partial amplitudes, elastic | A_{l s J}(E_n) = (U^J − 1) e^{2iω_l} / (2ik) per (l, s, J), complex; dσ/dΩ = \|f_C\|² + 2 Re Σ f_C*(θ) A P_l + Σ B_L^NN P_L | `Unit::CalcCoulombAmplitude`, `GenMatrixFunc::BuildCoulombInterference` |
| A_y numerator N_L | A_y dσ/dΩ = Σ_L N_L P_L¹(cos θ); the fold of A_y is (K·N)/(K·σ), as today's σ-weighted fold | `GenMatrixFunc::CalculateAmplitudeMatrix`, `PolarizationFunc`, `CalculateCaptureAnalyzingPower` |
| phase shifts δ_{Jls} | from U^J (eigenphase or diagonal, as the data set states); never folded | the phase-shift branch of `CalculateBareCrossSection` |
| total σ (optical theorem) | (2π/k²) Σ_J g_J Σ (1 − Re U_cc) | the total branch |
| β-delayed spectrum | f_β(Q − E) P_c(E) \|Σ_λ …\|² at E_n of the exit pair | the PType-20 branch of `CalculateBareCrossSection` |
| external capture | unit radial amplitudes per final state and multipole at E_n (`ECTable` already tabulates per term and energy), scaled by the ANC or γ in the evaluation | `Unit::CalculateECAmplitudes`, `GetECAmplitudeWithShift`, `ECAmplitudeCache` (per unit, interpolated when shifted) |
| THM σ_HOES | `THM.md` §1 at E_n, with the vertex tables per stage | `THMMatrixFunc::CalculateTHMCrossSection` per unit |

∂O/∂θ is today's adjoint (`AMatrixFunc::PointTBar`, `TBarAdjoint`, `LevelMatrixAdjoint`), taken once per (E_n, r)
with the cotangents of all points summed first (§8, gradient branch).

### 3.3 Kernels

Each kernel is a function of the node energies, the point's parameters and its own nuisances that returns the
weights w_in over the nodes of its support (the forward pass applies the same templated integrator to the node
values, so the arithmetic is today's; the weights are what the gradient reads — `TELinearForm` and
`Unit::TargetEffectWeights` already do this). Energies E c.m. of the entrance pair unless stated; ε the stopping
cross section per active nucleus (lab MeV cm²).

| kernel | formula | support | ∂ w / ∂(own parameter), analytic | replaces (`file:function`) |
|---|---|---|---|---|
| target layer / infinite / implanted | Y(E₀) = 10⁻²⁴ ∫_{E_out}^{E₀} ⟨σ⟩_g(E)/ε(E) dE, ⟨σ⟩_g a Gaussian average of width σ_eff² = σ_b² + σ_s²(E₀ − E); ∫_{E_out}^{E₀} dE/ε = N; thin target as average: Y/(10⁻²⁴ N) | [E_out − kσ_eff, E₀ + kσ_b] | ∂Y/∂N = 10⁻²⁴ ⟨σ⟩(E_out) (the back edge moves by −ε(E_out)); stopping scale s (ε → sε): ∂Y/∂s = (−Y + 10⁻²⁴ N ⟨σ⟩(E_out))/s; ∂/∂σ_b, ∂/∂c_strag from the Gaussian | `TargetFold.cpp: Unit::IntegrateKernelsImpl` (target branch), `GridQuadrature`, `TargetFallbackImpl`; `PointPass.cpp: Unit::CalculateThinTarget`; `TargetEffect.cpp: TargetThickness, MeanEnergyLoss, StragglingSigma, BohrStragglingSigma, GridExtent`; `TargetSetup.cpp: SetUpTargetEffects, AddTargetFallbackNodes`; `HandOver.cpp: FillDensitiesFromParams` (density by differences today) |
| Gaussian resolution | ∫ G(E − E₀; σ(E₀)) σ(E) dE / ∫ G dE over the support (normalised, as today) | E₀ ± r σ | ∂/∂E₀ (shift): w·(E − E₀)/σ² minus its mean; ∂/∂σ: w·((E − E₀)²/σ³ − 1/σ) minus its mean | `TargetFold.cpp` Gaussian branch; `TargetEffect::SigmaCM, GaussianKernel` |
| beam profile | Σ_i G W D σ / Σ_i G W (skewed Gaussians G, window W by erf, detailed balance D) | `BeamProfileSupport` | ∂/∂(shift of the window) from ∂W/∂a, ∂W/∂b (Gaussians) | `TargetFold.cpp` beam-profile branch; `TargetEffect::BeamProfileWeight` |
| UDR (time of flight) | w_i = ∫ R(τ) φ_i(τ) dτ, φ_i hat functions on the delay axis, Simpson exact | the table's τ range | no fitted parameter; ∂/∂E₀ through t(E₀) | `TargetFold.cpp` UDR branch; `TargetEffect::BuildUdrKernel, LoadUdrTable` |
| Doppler (free gas) | SAMMY free-gas kernel in √E, 4-point Gauss–Legendre pieces | `TargetEffect::DopplerWindow` | ∂/∂T exists (`Unit::DopplerTemperatureDerivative`) | `TargetFold.cpp: Unit::DopplerWeights` |
| rectangular resolution (branch 6, added by the recipe) | Y = (1/w) ∫_{E₀−w/2}^{E₀+w/2} σ dE, w lab × the pair's c.m./lab factor; format 3: `detector … energy_resolution rectangular width <w>` | [E₀ − w/2, E₀ + w/2] | ∂Y/∂w = [σ(E₀+w/2) + σ(E₀−w/2)]/(2w) − Y/w (`derivative`, tested; the fit takes central differences); w is a fit parameter (`width 14 keV fit`, kind `apparatus`, `apparatus_parameters.md`) | nothing: `experiment/src/KernelRectangular.cpp` |
| triangular resolution (drill of the apparatus-parameter recipe) | Y = ∫ t σ dE / ∫ t, t = 1 − \|E − E₀\|/h, h = f w/2 (w the base width, lab); format 3: `detector … energy_resolution triangular width <w> <unit> fixed\|fit …` | [E₀ − h, E₀ + h] | ∂Y/∂w = 2(R − Y)/w, R the rectangular average over the base (`derivative`, tested); w a fit parameter as the rectangular width | nothing: `experiment/src/KernelTriangular.cpp` |
| energy shift, folded data | the point's kernel centred at E₀(ν) = E + a + b√E + f(E; p) on fixed nodes | that kernel's | ∂w/∂ν = ∂w/∂E₀ · ∂E₀/∂ν | `ShiftFrames.cpp: Segment::UpdatePointEnergyWithShift` (rigid move of sub-points); `AZURECalc.cpp: FillShiftColumns` (central differences) |
| energy shift, unfolded data | exact node at E₀(ν) (§2.5); interpolated mode: linear or cubic weights | one node (or 2–4) | dO/dE · ∂E₀/∂ν | the same, and `WorkingCopy::FillEnergyShiftsFromParamsImpl` |
| detector attenuation Q_L, angular convolution | diagonal in L: row P_L(cos θ) Q_L | — | none fitted | `UnitLegendre.cpp: Unit::CalcLegendreP` (Q coefficients, `QTransformRule`) |
| electron screening | diagonal per node: f = exp(π η U_e / E) on σ before the fold | node | ∂f/∂U_e = (π η/E) f (exists) | `Unit::ScreeningFactor`, applied in `CalcEDependentValues` |
| THM fold | c.m. Gaussian, QF-window weights, OES S-fold (`THM.md` §2) | ±5σ | ∂/∂σ as Gaussian | the THM sub-point fold of `THMMatrixFunc` / `ThmNumerics` |
| SUM, RATIO, natural target, S factor, transmission | not energy kernels: expression nodes `Sum`, `Ratio`, `Weighted` (abundance, ∂ analytic), `Scaled`, `Transform` (Phase 4) | — | as in Phase 4 | nothing; they stay |

**Status (branch 6, `feature/upgrade-5-kernels`).** The kernels of the first five rows are rows of one table,
`Kernels()` in `experiment/src/Kernels.cpp` (`experiment/include/azure2/experiment/Kernels.h`): `target`
(`KernelTarget.cpp`: thin, thick, layer yield, with beam spread and straggling; the σ(Ē) fallback and the edge blend
it had are deleted, decision D3),
`gaussian`, `beam_profile`, `udr` (one file each) and the broadening `doppler` (`KernelDoppler.cpp`), each a fold
template on V = double (value) and V = `LinearForm` (the gradient's weights) and a support function (the interval the
set-up lays the grid over). `src/TargetFold.cpp` is deleted; the engine's fold is `FoldPoint` (`src/PointFold.cpp`):
it gathers the nodes and the point and makes one indirect call per point. Nodes and quadrature are unchanged, every
output is byte-identical. The derivative columns above are not yet written as kernel functions, except the
temperature derivative of `doppler` (`broadenDT`); a table kernel may carry `derivative` (∂/∂ its own parameter).

**Status (`feature/upgrade-5-apparatus-parameters`).** A table kernel's parameter whose row says `fittable` is a fit
parameter whenever the project writes it `fit` (the nuisance grammar the density and the stopping scale use): kind 12,
`apparatus`, one parameter per record and parameter, shared by every data set naming the record; its grid is laid out
at its layout value (`KernelParameters::layoutValues`); its derivative is central differences of the whole χ².  The
stopping scale is the first row of the same kind (`apparatus_parameters.md`).  No engine file names a kernel parameter.

The roadmap lists SUM, RATIO, natural targets and the S factor as kernel builders. They are already one table of
expression nodes with `value` and `backward`, tested alone; re-implementing them as kernels would write them twice.

## 4. Lab data with a fitted shift: the comparison in the frame of the data (M3 (b))

### 4.1 The rule

For a measurement whose values are lab differential (or any quantity that carries a lab → c.m. factor), with
J_i(ν) = dΩ_lab/dΩ_c.m. at the shifted energy E_0,i(ν) and the c.m. angle θ_c.m.,i(ν) of the lab angle there:

```
model in the lab:   m_i^lab(θ, ν) = m_i^c.m.(θ, ν) / J_i(ν)        (m^c.m. evaluated at E_0,i(ν), θ_c.m.,i(ν))
residual:           r_i = (m_i^lab − n y_i^lab) / (n δy_i^lab)        y^lab, δy^lab as read, never changed
χ² term:            Σ r_i²  or  rᵀ M⁻¹ r with a covariance (pulls whitened, as today)
```

### 4.2 What it moves: nothing but rounding (measured on the code, an identity)

Today `experiment::ConvertRow` (`experiment/src/Frames.cpp` 303–309, 343–345) multiplies **both** the value and the
error by the same J_i(ν): y^c.m. = J y^lab, δy^c.m. = J δy^lab; the residual is
`DataResidual(model, y^c.m., δy^c.m., n)` = (m − n y^c.m.)/(n δy^c.m.) (`include/ModelPass.h` 81,
`src/ModelPass.cpp` 362). Then

```
(m − n J y^lab) / (n J δy^lab) = (m/J − n y^lab) / (n δy^lab)
```

exactly: the pull of every point is the same in both frames. The weight 1/δy^c.m.(ν) does move with ν, but the datum
moves with it in the same ratio, so the χ² is not a model choice that differs between the frames. The covariance
terms work on pulls, and the relvector model is frozen in lab units (`DataCovariance.cpp` 269: m = y^lab + δy^lab r =
m^c.m./(n J)), so they are invariant too. The analytic gradient with respect to θ is invariant (J does not depend on
θ); the shift columns are central differences of the same residual (`FillShiftColumns`). SUM leaves are brought to the
base's frame by J_base/J_k (`Scale`) and RATIO data are dimensionless: in the lab, m/J_base = Σ f_k/J_k for both.

So M3 (b) moves **no χ² and no fit beyond floating-point rounding**: the two expressions differ by a few ulp per
residual (one product and one division in another order). Fits can end at different last digits through MIGRAD's
path; a calculation's outputs differ by ≤ 10⁻¹⁴ relative (inferred bound; the branch measures it). The projects
concerned are those with an active shift (fixed non-zero, free, group, √E or function) on lab-differential data:
by a search of the tree at least `tests/data_covariance/{dense,diag,fit,plain,rel}.azr`, `tests/angle_mapping`,
`tests/upgrade/T2` (segments 2 and 3), plus any old-format project with `shift`/`sqrtshift`/`shiftfn` on isDiff 1
data; the branch lists them from `--describe`.

What M3 (b) does change: the per-ν data frame loses y^c.m. and δy^c.m. (`FrameRow::yCM, dyCM`; `Segment::CMValue,
CMError`), the residual has one form for every data set, the written c.m. data columns become a display computed
as J·y^lab, and the lab-frame fold of §4.3 becomes expressible.

### 4.2a Status (branch `feature/upgrade-5-lab-frame`, 2026-11-11): done, and measured

The frame of a row holds kinematics only (`FrameRow`: E, θ_c.m., J); y and δy stay as read in the `Measurement`.
Every comparison is `DataResidual(f, J, y, δy, n)` = (f/J − n y)/(n δy) (`include/ModelPass.h`): the objective,
the per-segment χ², the profiled norms, the least-squares rows and their norm columns, both analytic gradients (the
model cotangent carries 1/J), the covariance pulls, the shift columns and the THM profiles.  The c.m. data columns
of the writers and of the API are J y and J δy, computed where they are read (byte-identical to before).
Measured against the same tree with the c.m. comparison (the branch's first commit, itself identical to the base):
- Projects that can move: those with a row whose J ≠ 1 (a probe of every project that loads without generated
  files): 13N, angle_mapping, channel_space, hybrid_potential, spin_orbit_potential, identical_entrance_reaction
  (c12c12, dd_pt), relativistic (he5_ang, li7pn_ang, pa_ang, each with _rel, li7pn_edge), sum_lab_frame (plain_data,
  sum_data, sum_free), T1, T1_start, T2, T2_true, T3, T3_start (`tools/classify.py` `LAB_FRAME_PROJECTS`).  Every
  other project: identical in three bases.
- Calculations (all projects × 3 bases, `tools/cls.sh`): every output file identical as printed, except
  `sum_lab_frame/models/sum_data` and `sum_free`, whose χ² is 3.99 × 10⁻¹⁷ (data equal to the model: χ² is rounding
  noise) and moves by 1–2 × 10⁻⁶ of itself, 7 × 10⁻²³ absolute.
- Model values: bit-identical (the model does not depend on the frame of the comparison).
- Per-point pulls, measured against the natural rounding scale ε (|f/J| + n|y|)/(n δy) of the two forms: at most
  1.3 ε (limit of §4.4: 8 ulp).  Relative to the pull itself the change is up to 2.7 × 10⁻¹² (T1) and 1.2 × 10⁻² for
  the pulls of 10⁻⁹ of `sum_lab_frame`: a pull that is a difference of two nearly equal numbers has no relative
  precision in either form.
- χ² at the start: ≤ 2 × 10⁻¹⁵ relative, except the 3.99 × 10⁻¹⁷ of `sum_lab_frame` (1.8 × 10⁻⁶).  χ² at the end
  point of the reference fit (restart from its `param.sav`): ≤ 2.5 × 10⁻¹⁵, except `sum_lab_frame/sum_free`
  (2 × 10⁻⁹ of 4 × 10⁻¹⁷).
- Fits at their ends (1 thread, Brune; MIGRAD unless named): end points move by rounding through MIGRAD's path, up to
  9 × 10⁻⁶ relative in a parameter (angle_mapping, a width), χ² up to 4.3 × 10⁻⁹ (spin_orbit_potential); the
  `--use-lm` fit of T2 is byte-identical.

### 4.3 A real difference it makes possible (decision D5)

Today a folded lab-differential point is folded at a fixed θ_c.m. (the point's) and compared through J at the point's
energy; physically, at a fixed lab angle both θ_c.m.(E) and J(E) vary across the fold window. In the lab frame the
kernel row can carry them per node: m_i^lab = Σ_n w_in Σ_L B_L(E_n) P_L(cos θ_c.m.(E_n; θ_lab)) / J(E_n; θ_lab).
This moves numbers by an amount of order ½ (Δ_window)² |d² ln(J P)/dE²| (inferred, small for keV windows; measured
when implemented).

### 4.4 The gate

Stated as: "every per-point pull of a calculation agrees with the reference to 8 ulp (relative 2 × 10⁻¹⁵); every
output number to 10⁻¹⁴ relative; every fit restarted from the reference fit's `param.sav` gives the reference χ² to
10⁻¹² relative at its first evaluation". Fits are not compared at their end points (their paths diverge by rounding).
Projects without lab-differential shifted data: bit-identical.

## 5. The yield of a layer of finite thickness

### 5.1 The record (format 3; `project_format3.md` §7, `data_kinds.md` §4.6)

```
target C13A
    kind           layer                                     # layer | infinite | implanted (profile) | sample | none
    active         13C                                       # the nucleus N counts
    composition    C 1                                       # stoichiometry (atoms per formula unit)
    areal_density  1.2e18 atoms/cm2 fit prior 1.2e18 15%     # N, active nuclei per cm^2
    stopping       evaluation scale 1 fit prior 1 4%         # or: formula "..." MeV*cm2/atom, or table "file" <units>
    straggling     bohr                                      # none | bohr | chu | yang | coefficient <c> keV^1/2
    phase          solid
dataset Y1
    quantity       yield 1/particle                          # or 1/uC (needs charge_state), counts/uC (efficiency)
    reaction       13C(p,g)14N
    target         C13A
    beam           B1                                        # spread sigma_b
    file           "data/y1.dat" columns E_lab value error
```

The quantity `yield` on `kind layer` was refused ("Phase 5 kernel", `project_format3.md` §7.2); branch 4 lifted the
refusal (status below §8). Units: E_lab (MeV), N (atoms cm⁻², active nuclei), ε (MeV cm² per active nucleus, lab), σ (b), Y
reactions per incident particle; per µC: Y × 10⁻⁶ C / (q e).

### 5.2 The formula, with the stopping module

```
Y(E₀) = 10⁻²⁴ ∫_{E_out}^{E₀} dE  ⟨σ⟩(E) / ε_eff(E)                    (lab E, σ at E_c.m. = E m_A/(m_a + m_A))
∫_{E_out}^{E₀} dE / ε_eff(E) = N                                       (E_out: the energy at the back of the layer)
⟨σ⟩(E) = ∫ G(E′ − E; σ_eff(E)) σ(E′) dE′ / ∫ G,   σ_eff² = σ_b² + σ_s²(E₀, E)
ε_eff = s · ε_compound(E) / x_active                                   (Bragg; s the stopping-scale nuisance)
σ_s²(x) = S(E(x))² ∫₀ˣ (dΩ²/dx′) / S(E(x′))² dx′                      (Tschalär propagation, Bohr/Chu/Yang dΩ²)
```

with the functions of `stopping/` (brief `AD_stopping_module_1.md`; names to follow that branch): ε(E) on a
`Material` and projectile; the energy after an areal density, E_out(E₀, N), and its inverse; their derivatives with
respect to the scale; the straggling σ_s with its propagation; the evaluation's relative uncertainty. The engine does
not integrate the range itself any more (`TargetEffect::TargetThickness`'s Newton–Simpson goes).

The quadrature is the target kernel of §3.3 on the grid nodes inside [E_out − kσ_eff, E₀ + kσ_b], with 1/ε
interpolated like σ (today's rule).

### 5.3 Limits (each a test)

- **Constant σ = σ₀** (any ε, no spread): Y = 10⁻²⁴ σ₀ N exactly, because ∫ dE/ε = N defines E_out. Detects an E_out
  inconsistent with the integrand's ε, a wrong 10⁻²⁴, a c.m./lab mix-up in ε.
- **Thin layer**, εN ≪ the scale of σ: Y → 10⁻²⁴ N σ(E₀ − Δ/2), Δ = E₀ − E_out; and the existing thin-target
  quantity is Y/(10⁻²⁴ N) (today's `σ` with a `layer` target, `data_kinds.md` eq. 2), computed by the same kernel.
- **Thick layer**, N → ∞ (E_out at the lower end of the grid): Y → today's thick-target yield (`kind infinite`) on the
  same grid, to the quadrature difference of the cut-off at the lowest node.
- **Isolated narrow resonance** (Breit–Wigner, ε, λ and partial widths constant across it): `data_kinds.md` eq. (7),
  Y = (λ_r² ωγ/(2π ε_r)) [arctan((E₀ − E_r)/(Γ/2)) − arctan((E₀ − Δ − E_r)/(Γ/2))] 10⁻²⁴, plateau λ_r² ωγ/(2ε_r).
- **Derivatives**: ∂Y/∂N = 10⁻²⁴ ⟨σ⟩(E_out) and ∂Y/∂s = (−Y + 10⁻²⁴ N ⟨σ⟩(E_out))/s against central differences;
  thin limit ∂Y/∂s → 0, thick limit ∂Y/∂s → −Y/s.
- **Per charge**: a ¹²C²⁺ beam gives half the yield per µC of ¹²C⁺.

The test project (`tests/layer_yield/`, synthetic, generator with the truth) holds one data set per limit.

## 6. How the stopping module is wired

### 6.1 What the target kernel calls

| need | `stopping/` function (as the AD brief specifies) | replaces |
|---|---|---|
| ε(E) per active nucleus | stopping cross section of the `Material` for the projectile, divided by the active atoms per formula unit, × 10⁻²¹ (eV/(10¹⁵ cm⁻²) → MeV cm²) | `TargetEffect::GetStoppingPowerEq()->Evaluate` (the user formula) |
| E_out(E₀, N) and N(E₀, E) | energy after an areal density; its inverse | `TargetEffect::TargetThickness` (Newton on Simpson) |
| ∂E_out/∂s, ∂E_out/∂N | derivatives with respect to the scale and N | finite differences of the density (`FillShiftColumns`, kind `Density`) |
| σ_s(E₀, E) | Bohr, Chu, Yang–O'Connor–Wang with Tschalär propagation; beam spread in quadrature | `TargetEffect::StragglingSigma` (c√ΔE), `BohrStragglingSigma` |
| mean energy loss | from E_out and the integrand | `TargetEffect::MeanEnergyLoss` |
| prior of the scale | the evaluation's relative uncertainty | — (new) |

### 6.2 The stopping-scale nuisance

One parameter s per (material, projectile species) of the project, multiplying ε_eff wherever that material is
crossed by that projectile; prior Gaussian with centre 1 and width the evaluation's stated relative uncertainty
(or the user's `prior` on the `stopping` field); fixed at 1 by default, fitted when the record says `fit`. It enters
the objective as one more prior term (`src/Objective.cpp`, the list of terms) and the parameter layout as a new kind
(`DataParameters.cpp`, `ParamIndexMap`). Its gradient: central differences of the whole χ² in branch 5 (as the
density); analytic through the kernel (§3.3) with the kernels branch. A user formula has no scale unless the record
gives one; an evaluation or a table has one, fixed at 1 unless the record says otherwise.

**Done in branch 5** (2026-11-11, `feature/upgrade-5-stopping-wire`): §6.1 for `evaluation` and `table` (formulas
unchanged: no output of a present project moves), §6.2, ERYA removed (§6.4 step 3). Not done: §6.4 step 2 (the old
integrators on the module's integrals). What calls what: `stopping.md` §14; the test with LunaSRIM's numbers:
`tests/stopping_yield`.

### 6.3 User formulas and ERYA

- **User formula** (`stopping formula "…" MeV*cm2/atom`): kept as an override, evaluated through `Equation` as today;
  the same kernel calls it in place of the evaluation. Every present project uses it, so the wiring changes no output.
- **User table** (`stopping table "file" <units>`): the module's tabulated override, with provenance.
- **ERYA is not called by the engine.** `external/erya` (with SRIM's `SRIM2013.xml`) is built under `USE_ERYA`
  (`CMakeLists.txt` 61, 375–388) and used only by `tests/reference/compound_stopping_reference_test.cpp`; the Qt GUI
  that used it is gone. No test project names it. Its removal is the CMake option, the directory and that test. The
  module's own stopping powers are SRIM-2013's since 2026-11-11 (ported from LunaSRIM, decision of the PI;
  `stopping.md`), and its tests compare them with LunaSRIM's values (`stopping/tests/data/srim_reference.csv`), so
  ERYA's SRIM file is no longer needed as a reference.

### 6.4 Order of removal

1. The kernel calls the module for `stopping evaluation` and `table`; the formula path stays (no output moves).
2. `TargetEffect::TargetThickness`, `MeanEnergyLoss`, `StragglingSigma`, `BohrStragglingSigma` call the module for
   formulas too (the module integrates any ε(E) it is given); the old integrators go. Outputs move by the difference
   of the two range integrations (Simpson-16 inside a Newton iteration stopped at 10⁻⁷ relative, against the module's
   adaptive quadrature; size inferred to be of that order, measured by the branch; class "quadrature").
3. ERYA, `USE_ERYA`, `SRIM2013.xml`, the reference test.

## 7. The gate of Phase 5

### 7.1 Classes of difference

Every changed number found by `tools/cls.sh` is put in one class; anything else is a failure to be explained.

| class | where allowed | criterion |
|---|---|---|
| **identical** | every project at every branch marked "bit-identical" in §8; at every branch, the projects that branch does not concern | `cls.sh`: identical |
| **rounding** | lab-frame comparison; gradient per energy | per-point pull to 8 ulp; output numbers ≤ 10⁻¹⁴ relative; fits: χ² at the reference `param.sav` ≤ 10⁻¹² relative |
| **quadrature** | the shared mesh; the shifted fold on fixed nodes; the range integration of the module | for each folded point, \|m_new − m_fine\| ≤ max(\|m_ref − m_fine\|, 10⁻¹⁰ \|m_ref\|), where m_fine is today's engine with every sub-point grid refined twofold (points per width ×2, base step ÷2, tail ratio ÷2): the new value must be at least as close to the converged value as today's |
| **physics** | a decision of §9 that changes a model (lab-frame fold per node, deleted windows / automatic tolerance) | listed per project with its size and the decision |

**The per-file rule for the shared mesh (decision of the director, 2026-11-12, to be shown to the PI).** The rule
above compares number by number; a mesh whose nodes lie elsewhere than today's, with an error of the same size, fails
it for about half the numbers by chance.  For the shared mesh a file is **quadrature** when its structure is the
reference's, no number that the refinement leaves unchanged has moved, and its largest relative distance to the fine
reference over its folded numbers, max |m_new − m_fine| / |m_fine| (m_ref ≠ m_fine), is not larger than **1.2 ×** the
same maximum of today's engine, **and** the root mean square of the same relative distances over all the file's folded
numbers, sqrt(mean (m_new − m_fine)² / m_fine²), is **not larger than 1.03 × today's**.  The 3 % is the resolution
of the comparison, not a tolerance for the engine: the k = 8 reference is itself 1/64 of today's distance from the
converged value, so the ratio of two distances to it is known to about 2/64.  Both are taken over every folded number
of the file, also those the new engine left in place.  Two further clauses, each a line of the tool: a file whose
largest new distance is at most **10⁻⁷** passes (both engines are converged there; 10⁻⁷ is the size by which the
Coulomb table moves a charged-channel number), and a project in the tool's table `NAMED` passes with its own factor
and is printed and counted as **named**, never as quadrature — two rows today: `hybrid_potential` (factor 1.25), for
the reason in the table below, and T2's `covariance.out` (factor 1.2, decision of the director, 2026-11-12): the
correlated χ² of the two lab-frame segments is dominated by one point (E 1.1228 MeV, pull +15) whose error shrank; the
old total was nearer the reference only because other points' errors cancelled part of it, and those shrank too (per
point the two segments are 3-4 times closer than before: largest 2.9e-4 / 3.0e-4 → 7.5e-5 / 1.2e-4).  Measured with a mesh that drops every other
node in the lower half of each point's range (13 folded projects of the quick set, Brune, against `tab-fine8` /
`tab-fine4`): the maximum alone fails 7 projects, the rms adds `13N`, `beam_profile_kernel`, `ratio_to_reference`
`folded` and `T2` (`T2_true` χ²); `thick_target_resonance` and `T2`'s own outputs pass both (the mesh there is finer
than today's grid by more than the dropped nodes cost).  With "no factor" on the rms (the first form of the rule) the engine
without the error failed in `13N` (rms 1.002 × and χ² 1.006 × today's: inside the resolution, hence the 3 %),
`layer_yield` `aa=3` (a fold exact to 3 × 10⁻⁸: 1.18 × and 1.44 ×, hence the 10⁻⁷ clause) and `hybrid_potential`
(1.21 ×).  A file above a bound is a failure to explain; the factors are not widened without a measured reason
written here.

`hybrid_potential`, the named project (measured by the agent, two rows re-run by the director; relative distance of
its 76 folded values to today's engine at k = 64, the largest and the rms agree to 0.3 %, the error has one sign):

| k | today's engine | shared mesh | mesh / today |
|---|---|---|---|
| 1 | 5.14 × 10⁻² | 6.17 × 10⁻² | 1.20 |
| 2 | 2.00 × 10⁻² | 1.71 × 10⁻² | 0.85 |
| 4 | 5.14 × 10⁻³ | 5.13 × 10⁻³ | 1.00 |
| 8 | 1.35 × 10⁻³ | 1.35 × 10⁻³ | 1.00 |
| 16 | 3.43 × 10⁻⁴ | 3.52 × 10⁻⁴ | 1.03 |

Both engines converge at second order to the same value; at the project's own settings both are 5–6 % from it (the
step is set by the resonance tails, 0.05 |E − E_r|; `points` has no effect: 10, 20, 40, 80 give the same bits), which
is before the asymptotic regime.  The finding for the owner is the 5 %, not the factor 1.2.  `tools/classify_quadrature.py` prints both classes per
file (the per-file one decides the summary and the exit status; the strict per-number count is printed beside it as
information) with the largest and rms distances of both engines.

**Measured on the merged engine** (director, 2026-11-12; whole tree, Brune basis, the shared mesh with shared sites
against the engine before the mesh; per project the largest value over its output files of the largest and of the rms
relative distance to the fine reference; `work/director/clsq_trial31.txt`).  Summary of the classifier: 577 files, 542
identical, 33 quadrature, 2 named, 0 failures (no-brune 550 / 517 / 31 / 2 / 0; use-park 577 / 542 / 33 / 2 / 0).

| project | largest, before | largest, mesh | rms, before | rms, mesh |
|---|---|---|---|---|
| `hybrid_potential` (named) | 5 × 10⁻² | 6 × 10⁻² | 5 × 10⁻² | 6 × 10⁻² |
| `ratio_to_reference`, folded | 1.1 × 10⁻² | 3.6 × 10⁻³ | 2.8 × 10⁻³ | 9.1 × 10⁻⁴ |
| T2 | 5.4 × 10⁻³ | 5.4 × 10⁻³ | 5.6 × 10⁻⁴ | 5.5 × 10⁻⁴ |
| `T2_true` | 4.8 × 10⁻³ | 4.8 × 10⁻³ | 5.3 × 10⁻⁴ | 5.2 × 10⁻⁴ |
| `stopping_yield`, fit | 4.2 × 10⁻³ | 2.0 × 10⁻³ | 8.2 × 10⁻⁴ | 3.8 × 10⁻⁴ |
| `thick_target_resonance` | 3.2 × 10⁻³ | 4.2 × 10⁻⁴ | 1.2 × 10⁻³ | 2.0 × 10⁻⁴ |
| `layer_yield`, fit | 2.1 × 10⁻³ | 9.9 × 10⁻⁴ | 9.6 × 10⁻⁴ | 4.4 × 10⁻⁴ |
| `shift_sqrt_subpoints`, unshifted | 1.3 × 10⁻³ | 6.7 × 10⁻⁴ | 6.2 × 10⁻⁴ | 3.4 × 10⁻⁴ |
| `layer_yield` | 8.5 × 10⁻⁴ | 4.9 × 10⁻⁴ | 2.0 × 10⁻⁴ | 8.3 × 10⁻⁵ |
| `natural_neon` | 7.8 × 10⁻⁴ | 2.9 × 10⁻⁴ | 1.1 × 10⁻⁴ | 5.1 × 10⁻⁵ |
| `target_effect_ranges` | 3.7 × 10⁻⁴ | 3.1 × 10⁻⁴ | 1.1 × 10⁻⁴ | 9.2 × 10⁻⁵ |
| `stopping_yield` | 2.6 × 10⁻⁴ | 1.1 × 10⁻⁴ | 1.1 × 10⁻⁴ | 5.6 × 10⁻⁵ |
| `rectangular_resolution`, fit | 1.4 × 10⁻⁴ | 5.3 × 10⁻⁵ | 5.4 × 10⁻⁵ | 2.6 × 10⁻⁵ |
| 13N | 1.4 × 10⁻⁴ | 1.4 × 10⁻⁴ | 3.5 × 10⁻⁵ | 3.5 × 10⁻⁵ |
| `beam_profile_kernel` | 8.5 × 10⁻⁵ | 4.7 × 10⁻⁵ | 3.5 × 10⁻⁵ | 2.1 × 10⁻⁵ |

The mesh is closer to the fine reference than the engine before it, or equal, in every project but the named one.
Where it gains most, two data sets or two segments fold over the same resonance and the mesh takes the finer of
their settings for both.  Fits (eight, 1 thread): the fitted parameters move by at most 0.14 of their fitted error
(`stopping_yield`, areal density, 7 × 10⁻⁵ relative) and, checked on that fit, towards the fit of the engine before
the mesh refined eightfold.

**The fine reference** is today's engine with k-fold refined grids (`AZURE_REFINE_QUADRATURE=k`): `tab-fine8` (k = 8)
wherever it has the file, `tab-fine4` (k = 4) otherwise (`beam_profile_kernel`, whose k = 8 run is too slow); the
classifier takes them as one argument, finest first (`tab-fine8,tab-fine4`), and prints the tag each file used.  The
error is of second order in k (measured by the director: today's default is about 64 times further from the converged
value than k = 8, k = 4 about 16 times), so k = 8 bounds the distance of today's engine to about 1.6 % of itself. A
number that k = 4 leaves unchanged to its last printed digit can still be folded (`natural_neon` `aa=2` line 723: k = 4
equal to today's, k = 8 one digit away): with k = 4 alone it counted as an unfolded number that moved.

### 7.2 What the director runs (minutes)

1. `tools/runall.sh` and `tools/cls.sh` against the reference tag, as today; on a bit-identical branch that is all.
2. On a quadrature branch, once per branch: a run of the *reference* engine with the refined grids (a flag the
   shared-mesh branch adds to the old code path in its first commit, before it replaces it: `--refine-quadrature=2`)
   gives the tag `tab-fine`; a script (`tools/classify_quadrature.py`, written by that branch) reads `tab-ref`,
   `tab-fine` and `tab-new` and classifies every changed number of a folded output by the criterion above. Unfolded
   projects (T1, T3, the shift projects without folds, 13N's unfolded data sets) must be **identical**.
3. Gradients = finite differences (`AZURE_GRAD_CHECK`) on T2, 13N and the layer-yield project; χ² smooth in a shift
   (a scan of 41 values of one free shift of T2 with no step above the second difference of its neighbours).

## 8. The order of work, as branches

| # | branch (`feature/upgrade-5-…`) | what it changes | files | deletes | gate | size | can run with |
|---|---|---|---|---|---|---|---|
| 1 (done) | `energy-tasks` | the plan's task is one compound energy (nucleus, E_x bits) with every unit at it: one level-matrix fill and solve per task, one T-matrix per (entrance, exit), the cross section per unit; sub-points are calculated in the energy pass, folds after it (the automatic-tolerance and fallback decisions read values already computed, so they choose the same) | `include/CalculationPlan.h`, `src/CalculationPlan.cpp`, `src/WorkingCopy.cpp` (`MapData`, `MapLeaves`, `CalculatePointsInParallel`), `src/PointPass.cpp`, `include/Unit.h`, `tests/upgrade/plan_test.cpp` | point mapping by (pairs, E) and `SharesCalculation`; the sub-point recursion of `Unit::Calculate`; probe copies | identical; a check mode (`AZURE_PLAN_CHECK`) compares the shared solve with a per-unit solve bit for bit on the whole tree; solves per evaluation reported (§1.2) | ~600 lines, 8–12 h | 3, 4, 5 |
| 2 (done, below) | `energy-tasks-gradient` | the analytic gradient per energy: cotangents of all units at E_n summed, one level-matrix adjoint per E_n | `src/AZUREGrad.cpp`, `src/AMatrixFunc.cpp` (adjoint entry points) | the grouped gradient by host (`BuildGradItems`), `ScopedPointMappingOff` | rounding (gradient sums reordered) | ~500, 8–12 h | after the pathway branch (it edits `AMatrixFunc`'s readers); with 3, 4, 5 |
| 3 | `lab-frame` | M3 (b): residual on y^lab, δy^lab and m/J(ν); frames keep kinematics only; the writer computes c.m. data for display | `experiment/.../Frames.h`, `experiment/src/Frames.cpp`, `include/Segment.h`, `src/Segment.cpp`, `src/ModelPass.cpp`, `src/AZURECalc.cpp` (`FitBarFn`), `src/DataOutputFiles.cpp` | `FrameRow::yCM, dyCM`, `CMValue/CMError` | rounding (§4.4) | ~250, 3–5 h | 1, 4, 5; **starts now** |
| 4 | `layer-yield` (**done**, below) | `yield` on a `layer` (formula stopping), the closed-form test project | `experiment/src/Quantities.cpp`, `experiment/src/Build*.cpp`, `src/TargetFold.cpp` (per-particle normalization), `src/TargetSetup.cpp`, `tests/layer_yield/`, `docs/source/…/data_formats.rst` | the refusal | identical (new quantity only) + the six tests of §5.3 | ~350 + test, 5–8 h | 1, 3; **starts now** |
| 5 | `stopping-wire` | §6: the kernel calls `stopping/`; `stopping evaluation` and `table`; the scale nuisance and its prior; ERYA removed | `src/TargetEffect.cpp`, `include/TargetEffect.h`, `src/HandOver.cpp`, `src/DataParameters.cpp`, `src/Objective.cpp`, `experiment/src/Validate.cpp`, `CMakeLists.txt`, `external/erya` (removed), `tests/reference/compound_stopping_reference_test.cpp` (moved) | `TargetThickness`'s integrator (step 2 of §6.4), ERYA | identical for formula projects (step 1); quadrature class for step 2 (its brief must name `external/erya`, which agents may not edit otherwise) | ~600, 8–12 h | after `stopping-1` is merged; with 1, 3 (after 4 or in coordination: both touch the target set-up) |
| 6 | `kernels` | the folds leave the unit: kernel functions of `experiment/` on (node energies, node values, kernel parameters), forward by the same templated integrator, gradient weights by the same code on `TELinearForm` | new `experiment/include/azure2/experiment/Kernels.h`, `experiment/src/Kernel*.cpp`; `src/TargetFold.cpp` (emptied, deleted), `src/AZUREGrad.cpp` (`GradTargetEffectAdjoint` calls the kernels) | `TargetFold.cpp`, the fold members of `Unit` | identical | ~900 moved, ~300 new, 6–10 h | after 1 and 4; with 2, 5 |
| 7 | `grid` | the shared mesh (§2.2–2.6); sub-points stop being units; the stage node table; kernel rows on mesh nodes; re-placement by anchors in the hand-over | new `experiment/.../EnergyGrid.h`, `experiment/src/EnergyGrid.cpp` (with `AdaptiveIntegrationGrid` moved in), `src/TargetSetup.cpp`, `src/WorkingCopyUnits.cpp`, `src/HandOver.cpp`, `include/Unit.h` | `ReplaceSubPointGrids`, `NodeTableCache.h`, sub-point units, `SubPointAccess.h`, per-unit sub-point tables; per D3 the automatic tolerance, windows, blend, fallback | quadrature (§7.1), unfolded projects identical | ~1,500, 16–24 h | after 6; alone (it touches the set-up) |
| 8 | `shift-kernel` | folded shifted data: the kernel moves on fixed nodes, ∂K/∂ν analytic; unfolded: exact node, ∂m/∂ν from dO/dE | `src/ShiftFrames.cpp`, `src/HandOver.cpp`, `src/AZURECalc.cpp` (`FillShiftColumns` for folded segments goes) | the rigid move of sub-points, shift tables of sub-points | quadrature; χ² smooth in shifts | ~400, 6–8 h | after 7; with 9 |
| 9 | `observables` | O_n = Legendre coefficients, elastic amplitudes, A_y numerators per (E_n, reaction); angle rows in the kernels; amplitude filters on dσ/dΩ (M8); lab-frame per-node fold if D5 | `src/GenMatrixFunc.cpp`, `src/UnitLegendre.cpp`, `src/PolarizationFunc.cpp`, `experiment/src/Kernel*.cpp` | per-unit angular sums, per-unit Legendre and Coulomb tables | rounding for the reordered angular sums; physics for D5 | ~1,200, 16–24 h | after 7 and after the pathway branch; with 8 |
| 10 | THM on the grid | σ_HOES at E_n, the THM kernel (`THM.md` §1–2) | `src/THMMatrixFunc.cpp`, `src/Thm*.cpp` | THM sub-point folds, THM-only mapping | quadrature | per `THM.md` | after 7 and 9 |

**Branch 4, done (2026-11-11, `feature/upgrade-5-layer-yield`).** The kernel is today's target integration
(`src/TargetFold.cpp`) not divided by N for an effect with `TargetEffect::IsLayerYield` (an effect of its own per data
set kind, so a σ and a yield on one layer share only the density, through the density group); E_out solves
∫ dE/ε = N with `GridQuadrature`'s rule on the point's own sub-points (`include/LayerExit.h`), not with
`TargetThickness`'s Newton–Simpson, which differs from it by the quadrature's error (up to 8 × 10⁻⁴ on the test's
coarsest grid). The grid reaches 2 % beyond the layer at its layout density; a grid that ends inside the layer gives
NaN, and the set-up refuses that case and a stopping that is not positive (`E_LAYER_EXIT_BELOW_MODEL`,
`E_LAYER_STOPPING`, `src/TargetLayer.cpp`). ∂Y/∂N stays central differences (the apparatus route of `AZURECalc.cpp`),
tested against 10⁻²⁴ σ(E_out). Where the stopping module plugs in: ε enters only as the sub-points' stopping (set in
`TargetSetup.cpp`'s set-up loops and `ReplaceSubPointGrids`) and is read by `GridQuadrature` and `LayerExitEnergy`
alike, so the scale s multiplies it there, once. Tests: `tests/layer_yield` (README).

**Branch 5, done (2026-11-11, `feature/upgrade-5-stopping-wire`).** A target says `stopping evaluation` (SRIM-2013 of
its `composition`, `phase`, `compound_correction`, for the projectile of the data set's entrance pair) or `stopping
table "file" <units>`; `TargetStopping` (`src/TargetStopping.cpp`) is the one place the stopping is computed, per
active nucleus ε_atom/f_active, and `TargetEffect::StoppingAt` multiplies the stopping scale once, where the sub-points
get their stopping and in the range integrals. The scale is parameter kind 12 (`stopping_scale.<material>.<projectile>`,
one per material and projectile, prior the record's or SRIM's stated accuracy), differenced on the whole χ² (central
differences, as the density). Every output of the present projects is identical; `tests/stopping_yield` compares the
plateau of a thick TiN layer with LunaSRIM's ε (8.3 × 10⁻⁵, the line shape), the width of a layer's step with
LunaSRIM's energy loss (2.4 × 10⁻⁷), a table of the library's values with the evaluation (2.1 × 10⁻¹⁰), and a fit of
two layers' N and the shared scale. pyazr binds the library (`pyazr.stopping`). ERYA is removed. Not done: §6.4 step
2; the analytic derivative of the scale.

**Branch 6, done (2026-11-11, `feature/upgrade-5-kernels`).** The folds of `src/TargetFold.cpp` (502 lines, deleted)
are kernels of `experiment/` listed in `Kernels()` (`experiment/src/Kernels.cpp`): `target`, `gaussian`,
`beam_profile`, `udr`, the broadening `doppler`, one file each, with the shared quadrature in `KernelQuadrature.h`;
`LayerExit.h`, `UdrKernel.h` and the linear form (`LinearForm.h`) moved to `experiment/`. A kernel row holds its
name, format-3 keyword, parameters, fold on double and on `LinearForm`, support and (optionally) a parameter
derivative; `TargetEffect` keeps the parameters as read and chooses its kernels once (`Kernels()`,
`KernelParametersFor`). The engine's fold is `src/PointFold.cpp` (one indirect call per point; the per-node work
stays inside the kernel, except a beam-spread formula, which is called through `EnergyFunction` as `SigmaCM` was).
Unit lost its ten fold members. The set-up takes the grid's interval from `EffectSupport` (the two copies of
`GridExtent`, the beam-profile, UDR and Doppler windows in `TargetSetup.cpp` are one call each). Every output of the
tree byte-identical. The proof of extension: the rectangular resolution (§3.3, last row) was added by the recipe in
`CLAUDE.md` with no engine file edited (the layer yield of branch 4 had touched 11 engine files).

**First steps before Phase 4's last branches are merged.** `lab-frame` (3) and `layer-yield` (4) touch no file of
the pathway branch (`CNuc*`, `GenMatrixFunc`, `AMatrixFunc` readers, `UnitExternalCapture`, `ExternalCaptureSetup`)
nor of the tree conversion (tests) nor of pyazr's file class; they can start at once. `lab-frame` touches
`src/Segment.cpp`, which the "engine reads format 3 only" branch also edits (the old reading path): a textual merge,
not a design conflict. `energy-tasks` (1) can start at once too if it leaves the external-capture lines of
`PointPass.cpp` (the per-sub-point `RecalcEDependentValues`) as they are; the pathway branch moves those amplitudes.
`stopping-wire` (5) waits for `stopping-1`. Three or four agents: {1, 3, 4} now, {5} when the stopping library
merges, {2, 6} next, {7}, then {8, 9}, then {10}.

**Expected speed (inferred from §1.2–1.3, to be measured by branch 1):** with one fill+solve per distinct energy and
one T-matrix per (energy, reaction), 13N's point pass goes from 258 ms to about 100 ms (fill+solve 77 → 5 ms, T-matrix
88 → 5 ms, cross sections unchanged), `thick_target_resonance`'s from 795 to about 230 ms, T2's from 70 to about
40 ms. Branch 9 then makes the per-unit cross section a short dot product.

**Branch 1 done (2026-11-11, `feature/upgrade-5-energy-tasks`), measured.** The point pass runs the plan's energy tasks
(`CalculationPlan::EnergyTasks`, `src/EnergyPass.cpp`): one fill and factorisation per (nucleus, route, E_x bits), one
T-matrix per group (entrance, exit, source of the external-capture amplitudes), the observable of every unit; then the
folds (`Unit::FoldSubPoints`). Every output of the tree is identical in the three bases. `AZURE_PLAN_CHECK=1` solved
every unit of the tree alone and compared it with its shared solve (three bases, 3,727,625 units): 0 differ, with the
key (nucleus, E_x bits) as designed -- but the T-matrix needs one more key than the pair: points and sub-points at one
energy may hold set-up external-capture amplitudes that differ in the last bits (found on
`tests/electron_screening/grad.azr`, 1,602 of 10,125 units, not on the tree's own projects; its gradient fit moved).
Per evaluation of a fit, 1 thread, median of `AZURE_PROFILE` (load 17–22 on 24 cores, ms; base = `b-trial18`):

| project | point pass, base → new | inferred | solves | T-matrices |
|---|---|---|---|---|
| T1 | 1.0 → 0.9 | — | 206 → 137 | 206 → 198 |
| T3 | 1.6 → 1.6 | — | 414 → 315 | 412 → 412 |
| T2 | 106.4 → 44.6 | 70 → 40 | 24,299 → 5,483 | 19,267 → 5,333 |
| 13N | 316.5 → 131.9 | 258 → 100 | 45,871 → 2,727 | 45,871 → 2,731 |
| `thick_target_resonance` | 1,102 → 264 | 795 → 230 | 48,237 → 2,148 | 48,237 → 15,036 |
| `natural_neon` | 3,207 → 1,718 | — | 291,838 → 56,696 | 292,742 → 131,199 |
| `18O_p_a_thm` | 23.5 → 23.4 | — | 14,688 → 703 | (THM) |

The cross sections are unchanged (now the largest stage: 13N 88 of 132 ms, `natural_neon` 896 of 1,718). THM gains
nothing: its fill was 7 of 23.5 ms, and the half-off-shell cross section, which reads each unit's channel rows, now
finds them out of cache (perf: `Unit::GetLoElement` 12 % of the run, from the fill's unit to the cross section's).
Whole runs (`perf stat` cycles, alternated with the base): fits of 8 MIGRAD calls 13N 1.95×, T2 1.66×, thick target
1.95×, T1 1.02×, T3 1.07×; calculations (6 runs each, set-up dominated) thick target 2.0×, `natural_neon` 1.34×, T2
1.11×, T1, T3, 13N, THM and `target_effect_ranges` within their scatter (±10–20 % at this load). Peak RSS within 2 %.
A 200-call T2 fit is not faster as a whole: at 1 thread 0.96× (its gradients' forward solves are still per host:
29.0 of 53 s), at 5 threads 0.89× (17.0 s against 14.8, two runs each; not explained, see HANDOFF).

**Branch 2 done (2026-11-11, `feature/upgrade-5-energy-tasks-gradient`), measured.** The analytic gradient, the least-squares rows and the model gradients run over the energy tasks (`src/GradEnergyTasks.cpp`):
the point pass, each point's cotangent back through its expression tree and fold to its units, then per task one fill
and factorisation, the T-matrix cotangents of each T-matrix group summed, one back-propagation of the merged bilinear
terms per J-group (rows: per row, or by the task's bilinear basis). The sums are formed in an order fixed by the tasks,
so a gradient and a fit give the same bits at any thread count. The 5-thread regression of branch 1 was the serial
rebuild of the tasks (a fit that moves a shift rebuilds them at nearly every pass: 1.15 s of a 15.6 s T2 fit at 5
threads); on flat arrays with parallel keys and sort it takes 0.37 s there (0.50 s over 234 passes at 3 threads).
Calculations identical in the whole tree; fit end points move by rounding (the gradient at the same parameters agrees
with the base to <= 4e-14 of its largest component; the natural_neon and screening fits stopped at the call limit
amplify that along their unconverged paths). Times (2026-11-12, machine shared with another user's 12 threads, load
10-13; min / median of 5 alternated runs for T2, median of 3 for the others; wall s):

| fit | b-trial19 (before branch 1) | b-trial21 (branch 1) | this branch |
|---|---|---|---|
| T2, 200 calls, 1 thread | 40.7 / 41.1 | 39.4 / 39.9 | 36.9 / 37.2 |
| T2, 200 calls, 3 threads | 18.9 / 20.6 | 20.6 / 20.7 | 18.4 / 19.1 |
| 13N `--use-gradient`, 200 calls | — | 50.7 | 26.9 |
| 13N `--use-lm` | — | 19.5 | 12.8 |
| thick target `--use-gradient` / `--use-lm` (no free parameter) | — | 13.8 / 7.9 | 5.2 / 4.6 |
| natural_neon `--use-gradient`, 30 calls | — | 219.6 | 105.7 |
| natural_neon `--use-lm`, 30 calls | — | 291.1 | 139.3 |

Peak RSS within 2 % of the base. Most of a T2 gradient is its finite-difference part (shifts, THM: about 1 s of the
1.1 s at 1 thread), which this branch does not touch. The 5-thread acceptance time was not measured (budget of 3
threads that night).

**Follow-up after the merge of branches 2, 5 and 6: delete the point mapping — done (2026-11-11, branch 7 part 1, `feature/upgrade-5-grid-nodes`).** The mapping no longer
saves anything: the energy tasks share the solves, and the gradient no longer uses it. Evidence: the whole tree with
`AZURE_NO_POINT_MAPPING=1` (tag `tab-egrad-nomap`, three bases) gives every output file identical to `tab-trial21a`
except the "Coulomb table: N lookups" line of `run.log` in `angle_mapping` and `sum_lab_frame/models/plain_data` (a
work counter; the director accepts that it moves); the T2 200-call fit at 1 thread took 38.4 / 38.3 s without the
mapping against 39.2 / 37.7 s with it (peak RSS 350.8 against 338.7 MB, +3.5 %); the `natural_neon` calculation
6.44 / 6.37 s against 6.46 / 6.34 s, RSS equal. What the follow-up deletes, by file:
- `src/WorkingCopy.cpp`: `MapData` becomes "build the plan" (`plan_.Build` with no mapping, one task per point);
  `SharesCalculation`, the shift-form helpers only it uses, `MapLeaves`, the switches `AZURE_NO_POINT_MAPPING`,
  `AZURE_NO_CROSS_SEGMENT_MAPPING`, `AZURE_NO_COMPONENT_MAPPING`; the `IsMapped`/`OwnsEDependentValues` conditions of
  the two set-up loops of `CalcEDependentValues` (every point computes its own values).
- `include/Unit.h`, `src/Unit.cpp`: `EnergyMap`, `IsMapped`, `GetMap`, `SetMap`, `ClearMapping`, `OwnsEDependentValues`,
  `SetOwnsEDependentValues`, `PlanMapped`, the members `energy_map_`, `own_e_dependent_`.
- `include/CalculationPlan.h`, `src/CalculationPlan.cpp`: `CalculationTask::mapped`, `MappedPoints`, the provider rule of
  `TaskUnits` (each unit is its own provider), the mapped clauses of `Check`; `src/EnergyPass.cpp`: the mapped points
  of a fold's ranges under `AZURE_PLAN_CHECK`.
- `src/PointPass.cpp`: the mapped points of `Unit::Calculate`, `MarkFailed`, `FoldSubPoints`; the `MappedPoints`
  parameter of `FoldThinTarget` (`src/TargetFold.cpp`, or its kernel after branch 6).
- `src/TargetSetup.cpp` (`ReplaceSubPointGrids`): `mappedSegment`, the "a mapped point keeps its host's grid" branch,
  the host/mapped grid invariant after re-placement.
- `src/ChannelTables.cpp` (`shareWithMappedOf`), `src/ExternalCaptureSetup.cpp` (the two `IsMapped` conditions of the
  component amplitudes and `intEC.dat`), `src/ShiftFrames.cpp` (the mapped host of a shifted point),
  `src/HandOver.cpp` (copying a mapped point's values; re-running `MapData` after a shift), `src/AZURECalc.cpp`
  (`MappingConfined` and the partners of the local shift columns), `src/AZURECalcMCMC.cpp`, `src/ModelPass.cpp`
  (`CalculatesItself`), `src/DataCheckFiles.cpp` (the mapping column of the data check file).
- Tests and documents: `tests/component_mapping/check.sh` (a, b: with/without), `tests/angle_mapping/check.sh`
  (a, b: the data check's mapping), `tests/shift_fn_groups/check.sh`, `tests/sum_lab_frame/check.sh`,
  `tests/pyazr/composite_gradient_test.py`, `tests/pyazr/failed_point_test.py` (the switches), the mapping sections and
  switches of `docs/source/developer/edata_and_target_effects.rst`, `.claude/skills/azure2-setup/references/data-definition.md`.
Verification: the full comparison (identical but the lookup lines), `AZURE_PLAN_CHECK=1` on the tree, the quick set.
Done as listed: every point is a task of its own (`CalculationPlan` holds one `PointRef` per task, `WorkingCopy::BuildPlan`
replaces `MapData`; the plan no longer changes after a shift, the energy tasks are rebuilt from the units as before);
`Unit` lost `EnergyMap`, `IsMapped`, `GetMap`, `SetMap`, `ClearMapping`, `OwnsEDependentValues`, `PlanMapped`,
`CopyEDependentValuesFrom` and three members; `Segment` lost `ShiftHost` and `UpdatePointGroupEnergyWithShift`; the
energy tasks lost their `providers` array (each unit is its own provider). The data check file keeps its mapping
column ("Not Mapped" on every row, the file's layout unchanged). The rewritten checks test what is still true:
`component_mapping` (a SUM extrapolation equals the sum of its terms; the components add no level-matrix solve),
`angle_mapping` (the shared solves of all segments equal each unit's own under `AZURE_PLAN_CHECK=1`; fewer energy tasks
than units). Evidence: the whole tree in three bases (tag `tab-nodes-m0` against `tab-trial23g`): 778 of 780, 751 of
753 and 778 of 780 files identical, the two others in each basis the "Coulomb table: N lookups" line of `run.log` in
`angle_mapping` (1516 → 1211 lookups, Brune) and `sum_lab_frame/models/plain_data` (110 → 230), as measured before
with the switch.

**Branch 7, part 1: the measuring tools for part 2 (2026-11-11, `feature/upgrade-5-grid-nodes`).**
- `AZURE_REFINE_QUADRATURE=k` (a diagnostic environment variable rather than the flag `--refine-quadrature` named in
  §7.2: diagnostics are `AZURE_*` variables, and `src/AZURE2.cpp` may not grow): every sub-point grid gets k times its
  points per width and k times its smooth intervals; the tail's step per distance follows the points per width
  (TailFactor ∝ 1/p), so it is halved at k = 2 with the constant `tailRatio` unchanged. Fallback nodes are not
  refined. Tested in `tests/layer_yield` (check 3b).
- `tools/classify_quadrature.py CMP REF FINE NEW`: §7.1's criterion on every changed number; numbers that the
  refinement leaves unchanged (unfolded) must not move.
- The thick layer against an independent integral: `tests/layer_yield` (K'), README there.
- Node counts (`tests/upgrade/node_count`, a tool built with the tests: every project of the tree that sets up, 24
  of them with a sub-point fold; set-up grids, menu 1; "mesh" is §1.2's estimate on today's grids, computed, not
  built):

| project | points | folded | units | distinct E_x today | shared mesh estimate | mesh / today |
|---|---|---|---|---|---|---|
| 13N | 593 | 37 | 46,121 | 2,727 | 2,786 | 1.02 |
| 17O | 23 | 23 | 19,560 | 2,021 | 1,905 | 0.94 |
| 18O_p_a_thm | 118 | 118 | 14,688 | 703 | 714 | 1.02 |
| 6Li_d | 62 | 62 | 3,162 | 3,162 | 1,182 | 0.37 |
| 7Li_p_a | 66 | 66 | 6,807 | 4,761 | 3,415 | 0.72 |
| beam_profile_kernel | 40 | 40 | 82,264 | 3,086 | 3,161 | 1.02 |
| hybrid_potential | 136 | 76 | 111,976 | 2,471 | 2,650 | 1.07 |
| layer_yield/models/fit | 41 | 41 | 60,367 | 3,133 | 2,175 | 0.69 |
| layer_yield/models/layer_yield | 116 | 87 | 87,634 | 4,137 | 2,762 | 0.67 |
| natural_neon | 6,271 | 3,107 | 292,742 | 56,696 | 25,313 | 0.45 |
| shift_sqrt_subpoints/models/{const,noshift,shiftfn,sqrt} | 28 | 28 | 36,841 | 2,089 | 2,125 | 1.02 |
| shift_sqrt_subpoints/models/dense | 28 | 28 | 56,028 | 8,004 | 2,961 | 0.37 |
| shift_sqrt_subpoints/models/narrow_{const,sqrt} | 28 | 28 | 44,212 | 2,158 | 2,192 | 1.02 |
| stopping_yield/models/fit | 93 | 93 | 176,270 | 5,014 | 2,331 | 0.46 |
| stopping_yield/models/stopping_yield | 175 | 175 | 199,089 | 7,295 | 2,552 | 0.35 |
| target_effect_ranges | 72 | 44 | 5,686 | 2,164 | 1,191 | 0.55 |
| thick_target_resonance | 56 | 56 | 48,237 | 2,148 | 2,197 | 1.02 |
| units_scales/models/yield_story | 33 | 33 | 37,372 | 1,782 | 1,830 | 1.03 |
| upgrade/T2 | 380 | 309 | 24,319 | 5,442 | 2,782 | 0.51 |
| upgrade/T2_true | 380 | 309 | 23,281 | 5,359 | 2,792 | 0.52 |

  No project's mesh exceeds 1.07 times today's distinct nodes (the stop condition of the brief was 1.5); where the
  grids already coincide on the anchored lattice the mesh adds the two edge nodes per point (≈ 1.02), elsewhere it
  halves the nodes (T2, `natural_neon`, the stopping and layer projects).
- The node table beside the sub-point units (`include/NodeTable.h`, check mode `AZURE_NODE_CHECK=1`): see the
  status of part 1 in §2.6.

## 9. Decisions for the PI, and risks

**D1. The grid inside a fit stage.** (a) Today's rule: re-placed in the hand-over when a quantized level anchor
changes, a function of θ alone (cost: the new nodes' channel functions, 15 ms of 93 in T2); (b) fixed through the
stage, wide enough to cover each narrow level's motion, re-planned between stages (no hand-over cost; a level that
moves further loses resolution, or a re-plan breaks MIGRAD's continuity). Recommendation (a).

**D2. Do the shared mesh (branch 7) at all?** (a) Yes, after the bit-identical branches; (b) stop at 1 + 6, the grid
being the union of today's per-point grids. Measured gain of (a) beyond (b): 1.9–2.2× fewer solves for T2 and
`natural_neon`, ≈ 1 for 13N, the thick target and THM; but (a) removes the sub-point units (98 % of the units: 142.7
MB per copy of `natural_neon`) and the grid machinery, and it is the only way the kernels become independent of
points. Recommendation (a), for memory and deletion, not for speed.

**D3. Automatic tolerance, fold windows, edge blending, the σ(Ē) fallback.** They exist to skip per-point sub-points;
on a shared mesh a node costs the same whoever uses it. (a) Delete them in branch 7 (the projects that use them move
by at most their own tolerance: in the tree `tests/target_effect_ranges`, and external projects such as 13C+α
ND2021r); (b) keep them as options. Recommendation (a), with the moved projects listed.

**D4. M3 (b).** It is an identity for the χ² (§4.2). (a) Implement it as the simplification it is (branch 3: one
residual form, no per-ν data frame), classified as rounding; (b) leave the c.m. comparison. Recommendation (a): it
deletes state and is what makes D5 expressible.

**D5. The fold of lab-differential data at a fixed lab angle.** (a) Per node θ_c.m.(E_n) and J(E_n) (exact); (b) at
the point's energy, as today. (a) moves folded lab-differential data by the second-order term of §4.3.
Recommendation (a), in branch 9, classified "physics" with the sizes listed.

**D6. The stopping scale.** (a) One per (material, projectile), fixed by default with its prior written; (b) one per
target; (c) fitted by default. Recommendation (a): the systematic is the evaluation's, shared by every target of
that material.

**D7. ERYA.** Not used by the engine (§6.3). Recommendation: remove it in branch 5 (or at once on its own).

**D8. Gradient bits.** Summing the cotangents of all units at one energy before one adjoint (branch 2) reorders sums;
gradients and fits move by rounding. Recommendation: accept, with the rounding class of §7.1.

**Risks, and the measurement that shows each early.**

- *Bit-equal energies across entrance pairs may not give bit-equal channel tables.* `ChannelKeys` forms every channel
  energy from E_x, so they should; the first commit of branch 1 compares the shared solve with a per-unit solve for
  every unit of the tree (`AZURE_PLAN_CHECK`). If one differs, the key becomes (entrance pair, E_x) and the T2 gain
  drops (measure: distinct (nucleus, entrance, E_x)).
- *Computing every sub-point where the automatic tolerance skipped most* can make branch 1 slower on such projects.
  Measure on `tests/target_effect_ranges` and the 13C+α project before and after.
- *Branch 7's mesh near many overlapping narrow resonances* takes the densest data set's step everywhere they overlap;
  the estimate of §1.2 says ≤ today's distinct-node count, but a project with very different `quadrature` settings
  per data set may grow. Measure the mesh size against the union of today's grids on every folded project first.
- *Shifted folded data before branch 8*: in branch 7 their nodes still move rigidly, so they keep their own nodes
  (today's own-segment mapping); only branch 8 puts them on the shared mesh.
- *The elastic kernel is affine and its rows depend on the node energy* (f_C(θ, E_n)): memory per (point, node, l);
  for `natural_neon` about today's sub-point count × l_max (inferred, a few tens of MB).
- *β-delayed and THM entrances* take E_x from the exit pair or the THM kinematics (`AMatrixFunc::InEnergy`): the
  grid key must use the same function; the check of the first risk covers it.

## 10. Where this note departs from the roadmap

1. **M3 (b)** was taken as a change of the fits ("a χ² weight that moves with a nuisance is a model choice nobody
   made", `data_model.md` M3). The weight moves, but the datum moves with it in the same ratio; the pull is the same
   in both frames (§4.2). The decision is a simplification that moves no fit beyond rounding.
2. **Most of the speed is bit-identical.** The roadmap's plan puts the speed in a new grid with a quadrature
   tolerance; the measurements put 5–22× of it in de-duplicating today's nodes, with no output change. That step
   (branch 1) is the start of Phase 5, and it is what Phase 6 calls "parallel over grid-energy blocks"; Phase 6 is
   left with the pools and copies.
3. **Energy shifts as interpolation weights** (Phase 5's kernel list) contradicts the exact-energies decision. For
   folded data the kernel should move on fixed nodes (an exact quadrature, smooth in ν, analytic derivative); for
   unfolded data the node moves; interpolation only in the opt-in mode.
4. **SUM, RATIO, natural targets and the S factor are not kernels**: Phase 4 made them expression nodes with
   `value` and `backward`; the kernel builders of Phase 5 are the energy-domain folds only.
5. **The quadrature gate** is better stated against today's own refinement (`tab-fine`) than as a tolerance per test:
   "at least as close to the converged answer as today's" needs no number chosen per project.
6. **ERYA's removal** is a CMake option and a reference test, not a migration: the engine never calls it.

## 11. Not established

- The speed of branches 1 and 7 is inferred from the counts and stage times; only branch 1 can measure it.
- The shared-mesh node counts are estimates (§1.2), not planned meshes.
- The size of the D5 change and of the rounding class of M3 (b) are inferred, not run.
- `shift_fn` was counted, not profiled; `natural_neon` was not profiled as a fit (a menu-1 calculation only).
- Whether the second channel-table pass of `thick_target_resonance` and T2 (rows asked again, none computed) is the
  component leaves or the external capture was not traced.

## 12. Director's review (2026-11-11)

Checked by the director, not taken from the report: the counts of §1.2 for 13N, T2 and `thick_target_resonance`
recounted from the raw dumps (45,872 / 2,728; 24,300 / 5,442; 48,239 / 2,149 solves / distinct (nucleus, E_x) printed
with 17 digits); the identity of §4.2 read in `experiment/src/Frames.cpp` (`ParticleValue` multiplies the value and the
error by the same Jacobian); the derivatives of §5.3 re-derived (∂Y/∂N, ∂Y/∂s). Not checked: the timings of §1.3,
the mesh estimates, everything §11 lists.

Status of the decisions of §9:
- D1 (the grid follows the level anchors, as today): taken; it keeps today's behaviour.
- D4 (M3 (b)): the PI had chosen (b) before it was known to be an identity for χ²; it stays, as a simplification.
- D8 (gradient sums reordered, rounding class): accepted; fits already depend on the thread count at that level.
- D2, D3, D5, D6, D7: **taken by the PI on 2026-11-11 ("yes to all")**, as recommended in §9:
  - D2: the shared mesh is built (branch 7).
  - D3: the automatic tolerance, the fold windows, the edge blending and the σ(Ē) fallback are deleted with branch 7.
    `tests/target_effect_ranges` is rewritten on the mesh's error bound; the converter turns the old options into a
    note of the conversion report, and projects outside the tree that use them (13C+α ND2021r) are converted with it.
  - D5: lab-differential data are folded with θ_c.m. and the Jacobian per node (branch 9, physics class: every moved
    project listed with the size).
  - D6: one stopping scale per (material, projectile), fixed by default, with its prior (branch 5).
  - D7: ERYA (`external/erya`, its CMake option and `tests/reference/compound_stopping_reference_test.cpp`'s use of
    it) is removed in branch 5, after the stopping library is merged and its comparison with ERYA is on record.

Corrections for the briefs:
- §5.1: the active fraction is the stoichiometric fraction times the isotopic abundance (`isotope` lines of the
  target), not the stoichiometry alone.
- §5.3, narrow resonance: the frames must be written out. With λ_r in the c.m. frame and ε_r in the lab frame per
  active nucleus the plateau is (λ_r²/2) ωγ (m_a + m_A)/(m_A ε_r); the test derives its closed form from the
  Breit–Wigner integral itself rather than quoting it.
- Branches 1, 4 and 3 are running (2026-11-11). The order of the others as in §8.
