# rmfit policy, objective, covariances, sharding and Park mode

The campaign's `campaign.json` holds a bounds policy (`BoundsPolicy`: what is free and inside which box), a search
configuration (`SearchConfig`) and an objective (`ObjectiveSpec`: χ² or a robust loss, weights, priors, covariance
reference). This reference says what each knob does and the rules for using it. Pass long values with
`init --policy @policy.json` (also `--search`, `--objective`, `--flags`).

## 1. Bounds policy

All cap lookups (bounds, sanity, level adds, template clipping) go through `BoundsPolicy.cap_for(E)`, keyed on each
level's seed energy E (MeV) in the model file.

| key | meaning |
|---|---|
| `theta2_cap`, `bg_theta2_cap`, `bg_energy` | θ² cap for physical levels and for background levels (seed energy ≥ `bg_energy`) |
| `theta2_cap_windows` `[[e_min, e_max, cap], ...]` | overrides the cap for levels seeded in a window (first match wins); relax the cap only where needed, never globally (a global relaxation also frees every background pole and saturated channel) |
| `new_level_theta2_cap` | cap for the level under test only (a ±0.5 keV window at its seed energy, kept on accept); θ² = 3 is 3× Wigner: watch for absorbers |
| `width_bound_windows` `[[e_lo, e_hi, Γ_max (eV), pair?], ...]` | bounds each free partial width, θ_c² ≤ Γ_max/Γ_W,c; with a pair key only that pair's channels. The bound is per channel: two channel spins in one pair keep up to twice the bound |
| `width_floor_windows` `[[e_min, e_max, Γ_min (eV), pair], ...]` | one-channel floor (one-sided sign lock, at most 0.95 of the cap); use a ±0.5 keV window on the seed energy |
| `total_width_floor_windows` `[[e_min, e_max, Γ_min (eV)], ...]` | floor on the total open width (penalty row, `total_floor_sigma` 1e-4 relative; at 0.005 the floor sagged by a third) |
| `min_level_width_ev`, `max_level_width_ev` | width limits of the level-test verdicts `narrow` and `broad` (0 = off); they govern level tests only -- polishes can still widen existing levels |
| `threshold_margin_mev` | 0.003: free level energies stay on the seed's side of every particle threshold by this margin (the Brune transform diverges at a threshold; the rwa-space fit does not see it, so the failure otherwise appears at export) |
| `shifts_free` | energy shifts free in every step (default false) |
| `e_window` | maximum energy window of a released level |
| `window_min_half_mev` | minimum half-width of the ±mΓ table |

Rules:
- **Windows act by energy, not by level.** A level that moves next to a windowed level is clipped at the next start
  (+100k on 373k once). Shrink width-bound windows to well inside the distance to the nearest other level, and run the
  clip check before every submission.
- **Any bound that can push a parameter away from zero must exempt the "level off" state** (all amplitudes exactly 0).
- **A width cap alone does not prevent absorbers**: Γ_W in open neutron channels is several MeV. The absorber signature
  is channels on the θ² cap and a failing Brune transform, never the width alone. Regions with 100-600 keV structures
  need `max_level_width_ev` ~2e6; a 6 MeV "level" is a background pole.
- **Every J^π group needs a background pole.** Groups without one grow resonances into background; check pole coverage
  per group before a level search. One free-width pole per group, added and freed, is the first move when levels sit
  at their caps. Broad levels cannot be replaced by poles alone.

## 2. Seeding new amplitudes

- A zero amplitude is not an engine key and has zero derivative: it never moves. Seed only channels that are exactly
  zero, in the data's range.
- Seed in units of γ_W, never in keV: write a tiny physical seed so the key exists (1e-12 eV for an open channel, an ANC
  of 1e-6 fm^-1/2 for a closed one), then set each key to frac·γ_W (θ² = frac², 0.04 is plenty) in the engine and
  `bake` the vector. A 1 eV seed in a high-ℓ channel just above threshold (P ~ 1e-12) is a reduced width far above
  γ_W: the level's Brune transform fails ("Denominator less than zero") and its unseeded widths come back rescaled
  (×2.08 once, with a perfect round trip).
- After seeding: the seed job's log has no "Denominator less than zero" line; every unseeded width of the baked file
  matches the input (max relative change ~1e-3); the seed moves the objective by less than the polish's first step.
- Filter a `.sav` by the names a session reports (`snapshot`'s `parameter.name`), never by an index computed from the
  file: the engine's level index differs from file order.
- A channel that starts with zero strength must not enter with a free norm on its data (the norm falls to its floor).
- Sub-threshold (ANC) amplitudes: list zero entries of channels whose threshold lies a few hundred keV above a level and
  seed the nearest at +0.05 γ_W, with a same-budget control. They do not replace a level at the right energy.

## 3. Objective

- **Plain χ²** is the default and the objective once the evaluator has examined the error bars.
- **Robust losses** (`ObjectiveSpec`): reduced χ² per data set is unusable when data-set sizes span orders of magnitude
  (it hands the fit to small sets; `weight_power` 0.5 does not save it); a per-data-set Sivia scale abandons whole
  sets; Sivia with one global scale (`loss_scale: global`, the overall χ²/N) is the coherent "poorly described sets
  have underestimated errors" objective. The unscaled Sivia loss is degenerate when χ²/N ≫ 1. Run 150-iteration
  polishes under each candidate objective and read the per-file table before committing a campaign to one. Basins
  found under any objective are reusable as crossover donors under plain χ².
- **Switching objectives** needs its own polishes before any level test: compute per-data-set scales at a polished
  plain baseline, polish twice under the new objective, then test.
- **`point_weights`** `[[ex_lo, ex_hi, w], ...]`: χ² × w on rows whose E_x lies in a window (penalties unweighted).
  Deltas under weights are not comparable to τ; judge on a plain polish afterwards.
- **`param_priors`** `[{"key": "g[1/2+#2 p2 L1 S0.5]", "mean": 0.537, "sigma": 0.062}]`: a Gaussian row (x − mean)/σ on
  the free vector, i.e. on the observed (Brune) reduced-width amplitude, with the amplitude's sign; counted in the
  penalty. A prior on a key that is not free raises. AZURE2's own `<parameterSettings>` priors are not seen by rmfit.

## 4. Data covariances

- `covariance <file>`, `normgroup`, `shiftgroup`, `shiftfngroup`, free `sqrtshift`/`shiftfn` work in single and sharded
  evaluation. Workers run without the covariance and the evaluator whitens the assembled rows with the engine's
  algorithm (to 1e-14), so a group split across shards keeps its correlations. A pyazr without data covariances ignores
  the keyword; `SingleSession` compares its χ² with an engine session at its first evaluation and stops on mismatch.
- Keep covariance files under `data/`: structure, window and levelscan variants link only `data/`.
- Row weights scale the plain rows before whitening; the Sivia loss is refused with a covariance. Per-segment χ² are
  partition shares and can be negative; `ev.point_chi2(r)` / `ev.point_residuals(r)` give per-row values.
- **relvector terms** scale with a model. `init` saves the seed's model as `cov_reference.npz` and every evaluator of the
  campaign uses that one reference, so polish = verify = export; verify and export also log the self-consistent χ²
  (`Score.chi2_self`). To move the reference, start the next campaign from the export, or pass
  `init --cov-reference <npz>` to share one between a control and a variant. Never refreeze inside a least-squares
  run. Campaigns made by rmfit builds before 2026-10-25 (expires: 2027-02-01) have no reference: each step froze the
  model at its own start, so their objectives are not strictly comparable; trust the verify.
- **A relvector covariance plus a free normalization on the same rows is biased under either freeze** (data-unit
  reference: the norm runs up; model-following: down to the floor), and the more so the larger χ²/N. `init` refuses
  it. The t0 treatment (Ball et al., JHEP 05 (2010) 075): fix the norm at its nominal value (`varyNorm 0`), add one
  relvector term per former norm parameter (its prior × 1 over the rows it scaled) to the covariance file, keep the
  campaign reference as t0, refresh t0 by starting the next campaign from the export (one iteration usually suffices).
  To keep a norm free on purpose: `--objective '{"cov_free_norm_ok": ["<cov file>"]}'` (logs a warning). Theory and
  the engine-side rule: `azure2-uncertainty`.

## 5. Sharding and convolution

- Shards split segments by energy, composite sum segments included (cost × component count; `RMFIT_SPLIT_SUM=0`
  restores whole sums). Each free-shift segment gets its own shard while they are at most half the shards
  (`RMFIT_SHIFT_SOLO=0`). Shift and shift-function groups go to one shard and keep their keyword.
- Jacobian-cost balancing (`RMFIT_JAC_BALANCE`, on by default) times each shard once, fits time = R(α + βk) with k the
  shard's differenced shift columns, repacks greedily and caches the plan in `<workdir>/jac_balance.json` (a one-time
  start-up of a few minutes per new shards directory). Results are identical to 1e-15.
- Held shift columns are skipped in the Jacobian (`cols=`, pyazr `objective().jacobian(x, skip=)`).
  Engine builds before 2026-10-30 (expires: 2027-02-01) re-evaluate the whole shard for every differenced shift column; a self-contained
  segment's shift now re-evaluates only its own rows, forward only. A running campaign keeps the library it loaded:
  restart it, or run on a newer frozen build, to get this.
- rmfit re-implements `<targetInt>` Gaussian convolution for shards: one σ per data point, normalised kernel over ±3σ,
  the `sigmalab`, legacy-expression and constant σ forms, quoted energy windows with blend width (points outside every
  window are bare; model = b × convolved + (1 − b) × bare), and fixed-angle conv pieces. A line with an auto tolerance
  stays a plain piece. `RMFIT_CONV_PPS` sets the lattice density: 5 agrees with the engine's single session to 6e-6 on a
  large model, 3 is too coarse; fits with windowed lines run with 5. Differential segments are convolved only through
  these windowed pieces; check that a sharded objective is finite before submitting.
- When the engine's convolution changes, re-run a shard-vs-single comparison per convolved segment before trusting any
  sharded number. Job scripts for convolved models carry a guard on the fields of `TargetEffectSpec` (`sigma_lab`,
  `ranges`) so an old rmfit fails at once.
- A single-session Jacobian with free shifts can be very slow on a large model; `verify()` needs only the objective.
  Run such a call in the background with a cap.

## 6. Park mode (`"use_park": true`)

- rmfit forwards `--flags` unchanged, so a campaign can run with `--use-park` (or `--no-brune`); `bake()` writes the
  `#parametrization` tag (0/1/2) the CLI uses to convert parameter files between bases.
- The J > 0 wall (one row per level, min(J, 0)/1e-3, J = 1 − Σγ² dS/dE) is carried as residual rows counted in the
  penalty, and in Park mode the objective always takes the residual path.
  Results of rmfit builds before 2026-10-24 (expires: 2027-02-01) did not see the wall in the LM polish: check `m.park_norms(x)`. Broad levels sit at
  J ~ 0.2-0.4 on large models, and a ×1.3 amplitude step can put most of them past the wall.
- `--sav` seeds are applied by name with no basis conversion. Convert first: a CLI mode-1 run in the target mode with
  that `.sav` as external file writes a converted, tagged `output/param.par`; seed from that.
- `theta2_cap` bounds the fit amplitudes: the observed θ² in Park mode, the formal one in Brune mode -- the same number
  is a tighter bound in Park mode for broad levels.

## 7. Windowed fits and campaign copies

- Never refit levels of a frozen block against window data alone: the seam belongs to a joint fit. A MeV-wide level at
  a window's low edge is a block-mismatch absorber.
- A window model's background is not transferable: hand over resonances, not background poles.
- A pole far above the window (50 MeV) needs much larger reduced widths than one just above it for the same effect.
- Copying a campaign: `structures.azr` in the ledger is an absolute path. On the copy run
  `UPDATE structures SET azr = replace(azr, '<orig>/', '<copy>/')` on `ledger.sqlite`, use the real path of a copy
  reached through a symlink, and check `Campaign('.').incumbent()` resolves before submitting. An export in a copy is
  named after the copy's base; check the name before a downstream script looks for `<reaction>.azr`.
- A model variant written into a subdirectory needs the base model's `data/` linked beside it (segment paths are
  relative).
