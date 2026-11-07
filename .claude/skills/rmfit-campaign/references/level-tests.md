# Level tests in rmfit: levelscan, structure add, verdicts

This reference covers the two rmfit commands that test levels -- `levelscan` (a residual-driven screen of template
levels inside the full model) and `structure add` (a level add refitted against a same-budget control) -- and how to
read their output. Where to look for missing structure, and the width-floor recipe, belong to
`azure2-resonance-search`; the acceptance threshold τ and ρ are defined in `../SKILL.md` section 4.

## 1. Candidates

`candidates.json` (written before `structure add`):
`{"candidates": [{"energy": E_x (MeV), "jpi": "5/2+" or null to scan, "energy_fixed": true, "channels": [{"pair": n,
"L": l, "S": s, "gamma": Γ_c (eV, with sign)}, ...], "source": "...", "note": "..."}]}`. Sources: ENSDF/TUNL
(`pyazr.nds.fetch_levels`), the reaction's literature tables, and residual runs (`diagnose.residual_runs`) in the worst
segments. The channel width key is `gamma` (`gamma_eV`/`gamma_ev` are accepted too).

`jpi: null` with a pinned energy lets the angular data assign J^π: the window scan runs every J^π with a local polish
and sends the best to the full model. Record a disagreement with the compilation's J^π as an evaluation decision.

Candidate energies are rounded to 1 eV before they are written (a 16-digit energy in an `add_level` line corrupts the
model; `azure2-fit`).

## 2. levelscan

`rmfit levelscan [--ex E ...] [--top K] [--gamma keV|auto] [--scales ...] [--ex-offsets ...] [--jpi ...] [--refit N]`.

Candidate energies come from the incumbent's residual runs clustered in E_x across segments, plus dispersion-shaped
pairs (two adjacent clusters of opposite sign, level between them). At each energy one dummy level per J^π (all the
group's carrying channels, energy free) is appended; every (J^π, width partition, sign pattern, energy offset, scale)
template is evaluated frozen in the full model; the best `--gn-top` get one Gauss-Newton (GN) step on the level's
amplitudes and the norms; the best per J^π is reported (ledger round kind `levelscan`). `--refit N` sends the top N by
global gain to the ordinary level-add refit.

Options and what they are for:

| option | default | effect |
|---|---|---|
| `--grid M,P` | 5,20 | the screen variant's `<targetInt>` lines get resonanceWidthMultiplier M and pointsPerWidth P; ~6.7× cheaper than 20/50 for a 0.03 % χ² change; the base model is untouched, so verify/export/report stay converged; `none` keeps the model's own setting. Never go to M = 2 (0.2-0.7 % shift, larger than Δχ² = 1). Re-measure per project |
| `--differential-window keV` | 250 | re-rank candidates on differential segments within this half-width of the energy; 0 ranks on the total |
| `--scales` | 0.3,0.6,1.2,2.5,5,10 | trial widths as multiples of `--gamma` (×2 steps: reach matters more than density) |
| `--stage-scales a,b` / `--refine-top n` | 1.2,5 / 3 | two-stage screen: every partition × sign shape at two widths, then the best n shapes per J^π at every width and offset (~7× fewer evaluations); `none` restores the exhaustive stage |
| `--max-signs n` | | caps the random sign patterns of the equal-split family (least informative for ≥ 12 channels) |
| `--gamma auto [--gamma-bounds 10,1000]` | | template Γ = the cluster's measured run width, clamped; `--ex` candidates carry no runs and fall back to 10 keV |
| `--rank density` | | order clusters by χ² per MeV of structure, `score / max(width or span, 20 keV)`, so short secondary-γ runs are not outranked by long total-cross-section runs |
| `--seg-match REGEX` | | only segments whose data-file name matches make candidates |
| `--ex-range lo,hi` | | drop candidates outside before the `--top` cut (one pass with `--refit` replaces two-stage scans) |
| `--merge-kev` | 60 | drop residual candidates within this distance of a stronger one |
| `--gn-top`, `gn_per_jpi` | | GN slots; 2 per J^π are reserved. A row marked "(frozen, not GN-stepped)" in the top five means `--gn-top` is too small |

Template construction:
- A particle channel enters the templates only if its Wigner width Γ_W ≥ `dead_channel_ev` (1 eV); left-out channels
  are logged per J^π. Closed channels (an ANC, no Γ_W) are excluded from the partitions.
- Every template amplitude is γ = γ_W √(Γ_c/Γ_W) from the engine's Wigner amplitudes, clipped at the policy's θ² cap.
  Templates converted through the `.azr` (a Brune transform of a 1 eV dummy among keV-wide neighbours) are wrong by up
  to 100×; the screen logs the x0/Wigner ratio per group and flags it outside 0.3-3.
- Families: single channel, equal split, dominant, prior partitions (the group's existing levels, matched by
  (pair, ℓ, s)), and pairs among the top-4 channels by Γ_W plus the priors' channels (50/50, 70/30, 30/70, both signs).
  A level that lives in the product of two channels is invisible without the pair family.
- The response to a new level is quadratic in its amplitudes and saturates at resonance, so templates need a finite
  trial width; a linear matched filter is identically zero.

Cost: the pattern count dominates (0.9 s per frozen pattern on a 35k-point model). A per-J^π log line gives
"n channels, S shapes; N1 stage-1 + N2 refinement evaluations in T s". A few thousand patterns take 30-50 min per
energy. There is no energy-range option other than `--ex-range`.

## 3. Reading a screen

- **The global GN gain decides whether a level exists; the local (differential) gain discriminates J^π among
  candidates the global gain already flags.** A narrow neutron s-wave level is paid for by the neutron total and every
  (α,α) angle, not by the ±250 keV window. The local ranking counts only when its best gain is ≥ max(1,200, 10 % of the
  global recovery); for a level the angular data cannot see (e.g. no α width in (α,n) angular data) the global ranking
  decides.
- **A frozen template is a lower bound on a level's worth.** A null screen is not a reason to skip a pinned add when an
  independent peak in the data points at the energy: the window scan's local polish sees what a frozen screen plus one
  GN step cannot. The GN gain predicted released gains to 16 % for narrow levels; for broad levels it underestimates
  by ~3.5× (the neighbourhood re-shapes), so a broad candidate with a screen gain of ~2,000 deserves the refit.
- **Before trusting a null at an energy the data resolve**, grep the log for `Denominator less than zero`, then drop the
  hypothesised level in at plausible widths through the same code path (the screen's variant, the incumbent mapped on
  by key from `ev.x0`, dummies' widths zeroed) and check that the frozen gain is nonzero.
- **Before believing a "wrong J^π" flag**, run that drop-in with the level's own partition: a level that recovers 99 %
  through its own partition but not through the screen's templates is a grid gap.
- **A candidate on top of an existing same-J^π level** says that level's strength or partition is wrong there. Test it
  as a pinned add with the released polish, or test the existing level's width and partition (`relocate` keeps the
  partition and is the wrong tool). A one-level template cannot represent a second state on top of an existing one:
  the Brune transform is singular and the window scan returns +inf. Print Γ_W per channel before asking whether the
  existing level can absorb a candidate.
- **Leading globally, absent from the local top five, within 0.1 keV of a same-J^π level** -- stop before the refit.
- **A large removal** (hundreds of thousands of χ²) makes global distances meaningless; compare recoveries.
- **Clusters that converge on one energy** are one feature; screen it once.
- **Not every large residual run is a level**: runs are large in point count; a misfit spread over MeV is a model
  short of poles or a data question (band map and θ² audit, `../SKILL.md` 3b).

## 4. structure add

Stages: window J^π scan at pinned energy on a ±`half_width` (0.5 MeV) window model, local polish per J^π → full model,
pinned polish → optional local release → released global polish → same-budget control polish (and the same local
budget) of the incumbent with the level off → classify → adopt if accepted.

- **Energy.** The new level's energy is free in the released stage within ±max(50 keV, 3Γ), capped at `policy.e_window`
  (`cand.release_energy = False` opts out). "relocated" means it ended within 2 % of that window's edge.
  rmfit builds before 2026-09-19 (expires: 2027-02-01) wrote added levels with `levelFix = 1`; free such levels (token 4 of their
  channel lines) and polish.
- **Dead channels**: particle channels whose Γ_W at the level energy is below `dead_channel_ev` (1 eV) are zeroed and
  fixed before the released polish and left out of τ's k. A level can otherwise pass the θ² sanity with formal θ² ~0.6
  in channels with no physical width.
- **Omitted channels**: only the channels the candidate names are created; the group's others are written γ = 0 and
  fixed, so no width can move into them. `RMFIT_SEED_OMITTED=1` (or `seed_omitted=True`) seeds every open channel with
  Γ_W ≥ 10 eV at min(1 eV, 1 % Γ_W). After the window scan, seeded-omitted channels are reset to their seeds before
  the window amplitudes are copied to the full model.
- **Effective k**: the pattern's channels always count; a seeded channel counts only if the fit used it (partial width
  ≥ 1 % of the level's summed width, or θ² ≥ 0.05). The log prints k and τ for pattern / effective / all.
- **Local release** (`--local-release W --local-nfev N`; `round --local-release W` for signs): frees every level of any
  J^π within ±W MeV, all background widths and the norms, so a move that pays off only after a neighbour of another J^π
  re-partitions is not ranked low. The control gets the same local + global budget. Minutes per move: restrict it to
  residual-run energies.
- **Width floor**: `total_width_floor_windows` holds Γ_tot = Σ_open Γ_W,c θ_c² ≥ Γ_min for a level seeded in the window
  (one-sided penalty row, σ = 1e-4 relative). Use it with `RMFIT_SEED_OMITTED=1`. The per-channel
  `width_floor_windows` is only for a deliberate test of one channel's width. A floor never acts on a level whose
  amplitudes are all exactly 0 (the control).
- **Seeds**: a negative result for a level that could be α-dominated is not trustworthy from a neutron-heavy seed. Run
  both seeds (e.g. 2:1 n₀:α₀ and 2:1 α₀:n₀ in the lowest ℓ) and keep the better.
- **Window table**: candidate and control are scored on one fresh single session; the ±1/2/3Γ table of Δχ² (inside /
  outside, per data set inside ±3Γ) is logged and stored with `control_cid`.
- **Never quote the window-scan delta as the level's value**: a truncated window model exaggerates it by an order of
  magnitude; it only ranks J^π.
- `structure add` exits 0 on a failed scan: test the printed JSON `outcome`, never the exit code.
- `run_level_add` adopts an accepted structure as the campaign's base, so two adds in one campaign are sequential.
  Alternatives that must start from the same model each need their own campaign directory (copy rules in
  `policy-and-objective.md`).
- A candidate above the data range crashes the window scan (no segments in ±0.5 MeV). Find the true data ceiling by
  intersecting each segment's window with its file, and call `Campaign.run_level_add(cand, half_width=2.5)` for
  tail-effect candidates.
- Window fits that use one data set only are a screen for J^π, energy and width, never a start vector: their
  rearranged neighbours and background are rejected by the unseen data. With a channel's data absent from a window,
  its amplitudes run to the Wigner limit; fix or cap them.

## 5. Verdicts of a level test

| outcome | meaning | adopted |
|---|---|---|
| accept | passes τ, θ² caps, energy window, ρ | yes |
| marginal / worse | below τ / worse than control | no |
| unphysical | a width on the θ² cap (within 0.1 %), a total-width floor binding (within 1 %), or a failing Brune transform; width-bound pins (`width_bound_windows`) are not counted | no |
| collapse | max θ² < 1e-4 and summed physical width < 1 eV (the level vanished) | no |
| narrow | passes, but Γ below `min_level_width_ev` | no, recorded |
| broad | passes, but Γ above `max_level_width_ev` | no, recorded |
| splits | a same-J^π neighbour with free energy within 1 MeV, below `bg_energy`, unresolved (separation ≤ the larger width), lost ≥ 50 % of the new width and ≥ 5 % of its own | no, evaluator decides |
| relocated | energy ended at the window edge | no |
| redistributes | ρ > 0.5 | flag, not a veto |

Also logged: "N % of the summed data-set gains come from <data set>" (`top_dataset_share`); a share near 100 % in a
finely stepped data set is the lineshape signature. "PROVISIONAL: within 20 keV of a same-J^π level" marks a lineshape
hypothesis (resolution, a single point, the existing level's partition) for the evaluator. The acceptance criteria test
whether a level helps, not whether it is a level.

A collapse answers only "a narrow level is not wanted"; to test a hypothesis about a width, hold the width with a
floor. A level added for a narrow feature that is accepted but ends wider than 2.5× the feature's measured width has run
off into broad strength: keep it, re-screen the same feature on the new model, add again (at most three attempts). Two
consecutive splits or marginals on a feature are the stop signal.

## 6. Order of work

- **Broad structure first.** Fit levels with Γ ≳ 10 keV over the whole range (`--gamma 20 --scales 0.5,1,2,4`, or
  `--gamma auto`) before narrow alignment; record narrow candidates near same-J^π levels as PROVISIONAL and move on.
  Stripping the broad levels to "rebuild narrow" fails: every candidate grows into the missing broad strength.
- **Then narrow levels at their widths** -- either a per-channel `width_bound_windows` of ~2× the feature width
  (permanent on accept), or peak-weighted χ² (`point_weights`) -- then the underlying structure.
- **Never put a re-screen offset within 1.5 keV of an existing level** (checked against fitted energies from the
  ledger): AZURE2 merges same-J^π levels closer than 1 keV.
- **Two broad same-J^π poles sharing strength** is what a broad feature looks like in R-matrix, not a doublet to
  resolve. Several broad levels at the Wigner limit in one band say the model is short of poles there.
- **A level that ends every polish on the threshold guard** wants the other side: run two same-budget variants (where it
  is, and moved across by the same margin) and judge by verified objectives and by the data the threshold feeds.
- **Diagnosing an excess** (model above data) is not a level-scan job: zero each level's amplitudes in the offending
  pair one at a time and tabulate the per-data-set change; bound the levels whose removal helps most cheaply.
- **Removal screens cannot see a degenerate same-J^π pair** (zeroing one without letting the partner absorb is
  catastrophic by construction): test such a pair with a targeted remove-and-refit.
- **Recall test**: remove established levels one at a time and re-add them with the standard floored add; every level
  the method should find must come back (energy within ~10 keV, width within ~15 %).
