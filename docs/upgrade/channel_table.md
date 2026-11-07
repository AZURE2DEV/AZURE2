# The per-stage channel table and the radial Coulomb waves of external capture

Every particle channel c = (α, l, s) enters the collision matrix at a compound-nucleus energy E_x through the
penetrability P_c, the shift function S_c, the hard-sphere phase factor e^{−iφ_c} and the Coulomb phase factor
e^{i(σ_l − σ_0)} at the channel energy E_c = (E_x − S_α) − e_α (MeV, c.m.; S_α the separation energy and e_α the
excitation energy of the pair α). These depend on (α, l, E_c) only, not on the channel spin s, the block J^π or any fit
parameter. This note fixes the table that holds them for the evaluation energies of a fit stage (ARCHITECTURE.md §1
decision 11, §3, §4; `coulomb.md` §2.3), the scratch that holds the Coulomb functions of the external-capture radial
integrals, what they replace, and how they are tested. It is step 3b (ii) and (ii′) of `cnuc_split.md` §3.

## 1. Where the repeated evaluations come from

The process-wide memo (`CoulFuncCache`) returned stored F, F′, G, G′ for any energy within 10⁻¹² MeV of a stored one.
Without it (branch `-nomemo`, c3dc22e5) the engine is 3–17× slower (`coulomb.md` §7). A `perf` profile of that build
(one thread, menu 1, `--no-coulomb-table`; `convolution_unify` `tgt_fixed` menu 3, two threads, first 240 s) shows
where the Coulomb functions are asked for again:

| project | share of the run in `coulsm` | where |
|---|---|---|
| T2 calc | 94 % | `EPoint::CalcEDependentValues` at setup (`EData::Initialize`), 95 % |
| 13N | 90 % | external capture at setup 77 % (node integrals of `ECTable` 44 %, the reduction √P at a_c 16 %), rows 16 % |
| natural_neon | 87 % | rows at setup 91 %, of which the isotopic component segments 60 % |
| `tgt_fixed` | 91 % | rows at setup of the target-integration sub-points, 93 % |

Nothing of it is in the evaluations: the four cases spend their extra time before the first χ². The repetitions are of
three kinds:

1. One (α, l, E_c) in several slots (J^π, channel) of a point: every channel spin s and every J^π block with the same
   l of the same pair repeats the same four numbers (`EDependentTables::Slot(jGroup, channel)`).
2. One E_x in several points: the component segments of a natural target repeat the energies of their data segment;
   target-integration and resolution sub-grids of neighbouring points share nodes; points of several data sets sit at
   the same energies.
3. The radial integrals of external capture: the terms (l_f, multipolarity, final state) of one entrance l_i at one
   node energy evaluate F, G at the same radii r, because GSL's `qagiu` places the same Gauss–Kronrod abscissae on the
   same sub-intervals; the reduction √P_{l_i}(a_c) is asked once per term and point.

Kinds 1 and 2 are evaluation energies of the stage; kind 3 is not (the radii are chosen by the adaptive quadrature).

## 2. The channel table (kinds 1 and 2)

**Definition.** A `ChannelTable` maps a key (α, l, E_c) to a row (P, S, e^{−iφ}, e^{i(σ_l − σ_0)}), dimensionless.
α is the engine's pair number, which fixes Z₁Z₂, the masses, the channel radius a_c, the hybrid potential and the
Coulomb route; E_c is the double formed by `EPoint`'s one expression, (E_x − S_α) − e_α, and the key compares it
exactly: two energies one ulp apart are two rows. No tolerance, so no row depends on which energy came first.

**Values.** One function, `ChannelTable::Compute(pair, l, E_c)`, gives a row, by the engine's rules unchanged: an
open channel (`ChannelFunc::IsClosed` false) takes P, S, e^{−iφ} from `CoulFunc` (physics/ `channel_functions`, the
Coulomb table or the hybrid model as configured) and e^{i(σ_l − σ_0)} = Π_{l′=1}^{l} (l′ + iη)/|l′ + iη|; a closed
channel (below, at or within the η > 100 band above threshold) takes S from `ChannelFunc::Shift` (the Whittaker branch)
and P = 0, both phase factors 1. The exact-threshold case of a scattering length (|E_c| ≤ 10⁻¹⁰ MeV) is not a row: it
stays in `EPoint`, as now. A row is a pure function of its key, so a value read from the table is bit for bit the
value `Compute` would give at that point.

**Who builds it, when.** `EData` builds it once per stage, before the threads that fill the points start:
`EData::CalcEDependentValues` (data segments) and `EData::InitializeComponentSegments` (component segments) first
collect the keys of every point and sub-point they will fill (`EPoint::ChannelKeys`, which walks the same J^π blocks
and channels and forms E_c with the same expression as `CalcEDependentValues`), sort them and remove duplicates, then
compute the rows in parallel, row i written only by the thread that holds i. The table is then published as a
`shared_ptr<const ChannelTable>` held by `EData` and passed to `EPoint::CalcEDependentValues`; nothing writes it
again. Rows depend only on their keys, so the table does not depend on the thread count.

**A new table, never a written one.** `Build(keys, previous)` copies the rows of the keys present in a previous table
of the same system (its pairs' Z, masses, radii, thresholds, hybrid tags and the Coulomb route compared at build) and
computes the others. The component segments' table is built so, from the data segments' one. A free energy shift or a
re-placed sub-grid moves the energies of one measurement between evaluations; the points refilled there do not read a
table in this step, they share their repeated (α, l, E_c) within the point (below). When the `EnergyGrid` of Phase 4
owns the evaluation energies, such a move builds `Build(new keys, current)` single-threaded between evaluations and
swaps the pointer, as ARCHITECTURE.md §4 states.

**A point without a table.** A point filled outside these two loops (an energy shift, a re-placed sub-grid, the
reaction-rate grid, a pyazr helper) computes its rows with `Compute`, sharing a repeated (α, l, E_c) among its own slots
(a list local to the call). The same function gives the same bits; only the time differs. The table counts its rows,
its lookups and the lookups it did not hold (`AZURE_PROFILE`), so a path that misses is visible.

**Readers.** `EPoint::CalcEDependentValues` takes its four numbers per slot from the row and forms the slot's values
as before: L⁰_c = S − B_c + iP (B_c is the block's), √P, and the two phase factors. The per-slot rows of
`EDependentTables` therefore stay in this step: they carry B_c, and 54 reads in `AMatrixFunc`, `RMatrixFunc`,
`THMMatrixFunc`, `GenMatrixFunc` and `EData` use them. Replacing them by an index into the table moves the
ownership of the per-point energies and is Phase 4 (`EnergyGrid`), as `coulomb.md` §2.3 said.

**Not in the table: `g_shiftDerivMemo`.** It holds S′ and S″ at the level energies E_λ in the Park basis. E_λ are fit
parameters θ, not evaluation energies: a per-stage table cannot hold them. It is an exact-argument memo, so it does
not make results depend on order or threads; it stays until `level_parameters` computes S, S′, S″ at E_λ once per
evaluation for all working copies (ARCHITECTURE.md §3), and is not changed here.

## 3. The radial Coulomb waves of external capture (kind 3)

The node integrals of `ECTable` (converged to 10⁻¹⁰), the direct integrals of `ECIntegral`, the reduction √P at
r = a_c and the entrance hard-sphere phase factor of each capture pathway evaluate F_l(η, kr), G_l(η, kr) at radii the
caller or the adaptive `qagiu` chooses. A thread-private scratch (`CoulFunc::SetRadialScratch`, set by these callers
only, not in the reaction-rate mode) holds them per block — Coulomb identity of the pair (Z₁, Z₂, μ, total mass if
relativistic, hybrid tag, GSL route), l, E_c — each block a hash map from the exact radius r to (F, F′, G, G′). At most
2²⁰ radii are held per thread; the scratch is emptied when full. It is not shared between threads and nothing else
reads it. Its values are those `CoulFunc` computes at the same exact arguments, so no result depends on what the
thread computed before, on the thread count or on the session, which is what the rule against caches (RULES §2.3)
guards.

A scratch of the last 16 blocks only (the first version) left 13N at 408 705 evaluations; the whole run needs 226 925
distinct (block, r); the bounded hash map reaches that number. The alternative, a fixed radial quadrature shared by
all terms, would change every external-capture integral at the 10⁻¹⁰ level of the node tolerance; it is left as an
open question.

## 4. What is deleted

- `CoulFuncCache` (`include/CoulFuncCache.h`, `src/CoulFuncCache.cpp`), its two reference tests, every
  `SetUseGlobalCache`, the `TransientScope`, the initialisation and clearing: the commits of `-nomemo` (c3dc22e5 and
  its docs, 37a3637b) re-cut on this branch.
- The memo-specific comments and the `cache_stats()` counters' meaning (they stay as keys, at 0, as on `-nomemo`).

## 5. Tests

- `channel_table_unit` (reference test): rows built from keys in two orders and with 1 and 8 threads are
  bit-identical, and equal to `Compute` at each key; `Find` returns no row for E_c one ulp from a key; a table built
  from a previous one equals one built from scratch; a previous table of a changed system (another radius) is not
  reused.
- Engine: outputs of the 179 `.azr` under `tests/` in three bases against the memo base (feature/upgrade 815d9b5f),
  ≤ 10⁻⁹ relative on observables or explained, and against `-nomemo`, where the rows are the same function at the same
  arguments and the outputs must be identical; the two pyazr THM byte checks five times; 1 vs 8 threads byte-identical
  on T1, T2, T3, 13N; T2 calc and fit200 (1, 8 threads), 13N, natural_neon and `convolution_unify` `tgt_fixed` no
  slower than the memo base on the same host.

Code: `include/ChannelTable.h`, `src/ChannelTable.cpp`, `EPoint::ChannelKeys`, `EPoint::CalcEDependentValues`,
`EData::CalcEDependentValues`, `EData::InitializeComponentSegments`, `CoulFunc::SetRadialScratch` and its callers in
`include/ECIntegral.h`, `src/ECTable.cpp`, `ECPathwayAmplitude` (`src/EPoint.cpp`).

## 6. As built (branch `feature/upgrade-1-channel-table`, 2026-11-07)

The table and the scratch as above, with two changes found by measurement: the scratch holds every block of the
thread up to 2²⁰ radii (16 blocks were not enough, §3), and it also serves the entrance phase factor of each capture
pathway (13N: 590 326 → 230 686 Coulomb evaluations). `EData::Initialize` releases the table after the component
segments, since nothing reads it until the energies change. `CoulFuncCache` is deleted (the `-nomemo` commits re-cut).

**Outputs.** 179 `.azr` × three bases, menu 1, one thread, every file compared number by number against the memo
base (feature/upgrade 815d9b5f, engine = 08516cfb): 486/689 (Brune), 458/662 (formal), 486/689 (Park) files
identical; observables ≤ 7.0 × 10⁻¹¹ relative (`energy_shift_sqrt`, `relativistic` he5 and he5_rel, T1 formal),
except the rounding residue of an exact zero (`ec_e1m1_phase` s_1, E1 at 0° and 180°, 1.5 × 10⁻²⁵ → 3.0 × 10⁻²⁵ beside
1.2 × 10⁻⁹ at 60°); 20 `intEC.dat` ≤ 4.0 × 10⁻¹¹ of the line's largest value; χ² of `energy_shift_sqrt` and
`shift_fn` in the 12th digit; every log without the memo's line. Against `-nomemo` (c3dc22e5) in the Brune basis:
689/689 files byte-identical, as it must be, since the table and the scratch return the values recomputation gives.
No golden and no `expected/` file is re-baselined.

**Checks** (GCC 11.5, GUI off, Condor): no compiler warnings; `channel_table_unit` passes; `tests/upgrade/check.sh`
PASS; `ctest -LE slow` 109/109; `tests/run_tests.sh` 132/132 (`shared_point_tables` (d) included); the pyazr THM byte
checks `pyazr_thm_profiled` and `pyazr_thm_fixed_angle` 5/5 each; 1 vs 8 threads byte-identical on T1, T2, T3 and
13N (all output files, `intEC.dat` included); ThreadSanitizer (clang 23, archer), section `upgrade` (T1–T3 calc,
extrapolation, rate, fit; T2 fit30), 8 threads: 0 reports.

**Speed** (host qa-xp-021, the memo base and this branch alternately, three repeats, `--no-coulomb-table`; mean and
range; rows in `tests/upgrade/bench/T2.tsv`):

| run | memo base | this branch | change |
|---|---|---|---|
| T2 calc, 1 thread | 5.21 s (4.77–5.54) | 4.86 s (4.50–5.18) | −7 % |
| T2 calc, 8 threads | 1.91 s (1.75–2.08) | 0.87 s (0.82–0.92) | −55 % |
| T2 fit200, 1 thread | 198.0 s (175.7–219.7) | 211.2 s (196.6–230.3) | **+7 %**, slower in each pair |
| T2 fit200, 8 threads | 118.5 s (111.8–130.7) | 105.7 s (102.8–108.3) | −11 % |
| 13N, 1 thread | 32.4 s (27.6–37.0) | 35.8 s (34.5–36.6) | +10 % (host noise ±15 %); a paired run on a quieter slot: 13.3 → 14.0 s, +5 % |
| natural_neon, 1 thread | 40.1 s (33.5–47.2) | 32.6 s (24.8–42.2) | −19 % |
| `convolution_unify` `tgt_fixed`, 2 threads | 97.2 s (87.8–105.9) | 97.3 s (90.8–102.3) | 0 % |

Peak RSS: T2 fit200 at 8 threads 481 → 425 MB; 13N 163 → 176 MB; natural_neon 1437 → 1484 MB; `tgt_fixed` 7668 →
7939 MB (+3.5 %).

**The gate is not met** for T2 fit200 at one thread (+7 %) and, within the noise, for 13N (+5 %). Where the time
goes (`perf`, one thread):

- T2 fit (60 calls): 72 % of the run in the Coulomb functions of `EData::FillEnergyShiftsFromParams` →
  `ESegment::UpdatePointEnergyWithShift` → `EPoint::RecalcEDependentValues` (memo base: 54 %). Each change of a free
  shift ν refills the rows of its points at the new energies, as decision 11 prescribes; the numerical gradient steps
  ν to ν ± h and back, and the memo found the energies of an earlier ν again. These points are filled outside the
  table (§2, "A point without a table") and share only within the point.
- 13N: about 20 % more samples in `coulsm` and 2 % in the scratch's hash lookups. The memo shared F, G at a_c
  between the point rows and the external-capture calls at the same channel energy (`Reduction`, the entrance phase
  factor), which the table and the scratch keep apart.

Proposed (next step): (1) the shift step sets the energies of a segment's points and sub-points first, builds
`Build(keys, previous)` from the segment's tables of its last two or three shift states (kept by the segment, not by
a process-wide memo), then refills the rows from it; (2) the external-capture callers read √P and e^{−iφ} at a_c from
the point's channel table where the key is present (the entrance and final pair of a capture are table pairs, with
E_c formed by the table's expression).

