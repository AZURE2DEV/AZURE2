# First fits: starting parameters before a global fit

A global fit started from poor parameters trades missing resonance strength for normalizations and stalls. This page
gives the procedure for building reasonable starting parameters for a new reaction before any global fit, the
channel-completeness rule for every fit range, and how to hand the result to an rmfit campaign. It is a provisional
procedure: re-evaluate it once the first reaction that used it has a converged evaluation, and change it by what worked.

## Procedure

1. **Fix the scale with one representative data set.** Fix its norms at 1 (the group's own data when they cover every
   channel); let every other set's norm vary. Sets with consistent shapes but a different scale then do not hold the
   fit back.
2. **Seed the lowest-l channel of every pair, for every level.** Give it a small nonzero amplitude, 0.01 γ_W (γ_W the
   Wigner amplitude, set in the engine's reduced-width units, not as a physical width), even where the literature lists
   no branch to that channel. A width of exactly 0 is fixed and never moves (SKILL.md, "Free and fixed").
3. **Every channel that opens inside the fit range is in the model, by default.** Whenever the fit range is extended
   (a new energy window, or a data set that reaches higher), list every partition open below the top of the range:
   each excited state of each residual nucleus (ENSDF) and the other light-particle partitions (d, t, ³He, ...), with
   thresholds E_x = S + E₂* (separation energy plus the excitation energy of the pair). Add a pair for each before
   fitting, without being asked: the lowest allowed l of each channel spin on every level; seeded per item 2 on levels
   at or within 0.5 MeV below the threshold and on background poles, inert (0, fixed) elsewhere. Add the new exits to
   every summed segment (total cross sections, total (n,p), (n,α), ...). Do not cap a data window at a threshold
   instead; a cap is a stop-gap, logged in the readme with the channel that forced it and listed for the evaluator. In
   one extension a ¹⁴N+n total cross section had been capped at E_n = 5.1 MeV because ¹⁴N(n,n₄) had no pair, which
   removed the strongest constraint, a norm-fixed total, from the top 1.2 MeV of the fit; seven pairs (n₄, α₃, d₀, p₁,
   n₅, n₆, p₂) were added by a builder script. Adding pairs: `azure2-setup`.
4. **Level by level, from low energy up.** Polish one level at a time (its energy and widths, plus the free norms)
   against the data around it, and see how well it can match.
5. **Then add levels** where a resonance is clearly present and no known level matches it in energy or width
   (`azure2-resonance-search` for the gated test).
6. **Only then a global fit**, staged as in SKILL.md, or as an rmfit campaign (below).

Why: global polishes from a poor start (χ²/N ≈ 200, new channels near zero) traded the missing strength for
normalizations. With every norm free, almost all fell to 0.38–0.68 within 150 evaluations; bounded to 0.7–1.5, eight
sat on a bound and the fit stalled at χ²/N = 12.7. With χ² ~ 10⁵, norm priors of a few percent cost nothing.

## Handing the start to an rmfit campaign

The campaign procedure is `rmfit-campaign`; the start is made in a new dated campaign directory under the reaction
directory (`$CAMPAIGN_ROOT/<reaction>/<M>-<D>-<YY>_rmfit_campaign/`), with `$RMFIT_HOME` on `PYTHONPATH`:

```bash
mkdir <dir> && cp <base>.azr <dir>/ && cp -r data <dir>/data && mkdir <dir>/output
python3 -m rmfit.cli <dir> init --base <dir>/<base>.azr [--sav <base>.sav] --nshards 24 \
    --flags '{"use_long_wavelength": false, "ignore_externals": true, "use_brune": true}' \
    --policy '{"theta2_cap": 1.0, "bg_theta2_cap": 3.0, "bg_energy": <Ex above the data>, "e_window": 0.3, "norm_range": [0.2, 5.0]}' \
    --search '{"top_k": 6, "pairs_per_group": 4, "stage1_nfev": 30, "stage2_nfev": 60, "union_nfev": 100}'
```

- `--flags` is passed to `pyazr.Session()`, so the basis is chosen there: `"use_brune": true` (Brune), `"use_park": true`
  (Park: the fit vector holds observed reduced-width amplitudes, Γ = 2Pγ²; `azure2-exchange`), `"use_brune": false`
  (standard Lane–Thomas). One basis per campaign: the ledger's vectors, `.npz` snapshots and `.sav` companions are all
  in it.
- **Objective options** (`--objective JSON`, stored as campaign.json `objective`; deBoer et al., Rev. Mod. Phys. 89
  (2017) 035007, section on fitting): `{"dataset_weight": "reduced", "dataset_key": "file", "loss": "sivia"}`.
  `dataset_weight: reduced` divides every data set's χ² by its number of points and rescales by N/K so the objective
  keeps a χ² magnitude (penalty rows unweighted); it stops a large set with small errors from dominating.
  `loss: sivia` replaces χ² by Sivia's broader likelihood, −2 ln[(1 − e^{−R²/2})/R²] − 2 ln 2 (about R²/2 for small
  residuals R, about 2 ln R² for outliers, which are down-weighted as if their error bars were inflated). It is applied
  to residuals standardized by each data set's reference χ²/N (`loss_scale: dataset`, computed by `init` from the
  plain baseline): without that, on a model with χ²/N of 20–300 every point sits in the logarithmic tail, the data lose
  about 10× their weight against the quadratic norm penalties, and the fit drifts (one polish took the plain χ²
  701k → 1,190k). Both methods are statistically incorrect but useful when error bars are underestimated. The
  objective is a property of a campaign: candidates verified under different objectives are not comparable, so a new
  objective means a new campaign directory (seeded from the previous export). `REPORT.md` and `segment_scores` keep
  the plain per-data-set χ² whatever the objective.
- `init` verifies the baseline (one full-model evaluation), registers the structure, writes `STATUS.md` and the batch
  job script. Re-running `init` keeps an existing configuration unless new values are given.
- **Policy.** `theta2_cap` is the Wigner-limit box for physical levels (θ² = γ²/γ²_W); `bg_theta2_cap` the box for
  background poles (levels with E_x ≥ `bg_energy`, whose energies stay fixed). 3.0 is a starting guess, not a physics
  default: check where the baseline's background amplitudes already sit (`init`'s cap audit); one ¹²C+α elastic model
  needed about 100. `e_window` is the energy freedom of a level per polish (MeV); `norm_range` multiplies the nominal
  normalization. `narrow_window` (e.g. `{"gamma_kev": 5, "factor": 3, "min_kev": 1}`) keeps a level whose free partial
  widths sum to less than `gamma_kev` within ±max(`min_kev`, `factor`·Γ) of its current energy instead of ±`e_window`,
  so a polish cannot carry a narrow level off its peak (a 9/2⁺ level once drifted 180 keV in one polish). The shard
  count is the core count; sharding is exact (χ² is additive over segments), so more shards only buy speed.
- **Cost check** (`campaign.log` after `init`): one sharded objective should take a few seconds and a Jacobian tens of
  seconds. If a segment with target-effect convolution dominates, rmfit folds it itself (`model.ConvSpec`); composite
  sum segments are allowed as convolution pieces.
