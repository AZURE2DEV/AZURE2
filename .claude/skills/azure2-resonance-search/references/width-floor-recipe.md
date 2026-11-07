# The width-floor recipe: a structure you can see but no search finds

Use this when the evaluator points at a broad structure that every screen has declared null. The failure it cures:
in a released polish a broad seed shrinks to Γ ~ eV while the neighbouring levels absorb the misfit, and the collapsed
candidate is then scored "nothing there". A floor on the candidate's total width keeps the broad hypothesis alive long
enough to be judged.

Worked case, for scale: a broad backward-angle (α,n₀) structure at E_x ≈ 9.4 MeV in ¹⁷O was null in template screens
(all J^π, Γ = 50–400 keV), in sign rounds and in a first local-release add (every J^π collapsed to Γ = 11–183 eV). With
the floor, a 1/2⁻ level at E_x = 9.348 MeV, Γ_n0 = 198 keV, Γ_α0 = 15 eV improved χ² by 10 950 against its control,
and held in a free refit (9.3474 MeV, 198.8 keV, 16 eV).

## 1. Localise and characterise the deficit before fitting

- Band map of χ² per 0.2 MeV of E_x by data set, from the incumbent's per-point residuals; then per 10 keV with the
  mean pull per band and per angle group.
- A broad deficit shows the model below the data by 2–3σ at every angle over ≥ 50 keV with one pull sign (in the
  worked case: the backward angles, cos θ_c.m. 0.35–0.71, over 400 keV). A narrow spike (two 10 keV bands, pull
  changing sign) is a line-shape or resolution question, not a level search.
- Write the anatomy into the work log first.

## 2. Width-floor level add, every J^π, local release, same-budget control

```bash
rmfit structure add --candidate i --local-release 0.5 --local-nfev 80 --pinned-nfev 60 --nfev 150
```

with, in the campaign policy,

```json
"total_width_floor_windows": [[E-0.0005, E+0.0005, 100000]]
```

(E_x in MeV, floor in eV) and `RMFIT_SEED_OMITTED=1` in the job environment, so that every open channel of the J^π
group exists and can take width.

- The floor acts on the total width over the level's open channels. The per-channel `width_floor_windows` forces the
  width into one named channel; use it only when the anatomy names the channel. A level add has only the channels the
  candidate names.
- Seeds: e.g. Γ_n0 = 150 keV and Γ_α0 = 30 keV in the lowest l.
- One campaign directory per J^π (an accepted level is adopted into its base, so groups must not share one); three
  groups per 16-slot node job.
- The floor is one-sided on the sign the amplitude starts with, never above 0.95 of the θ² cap, and skipped for an
  amplitude that is exactly 0 (the control's "level off" state).
- Check the control's start objective equals the incumbent's before reading any difference; a control that starts
  elsewhere is broken.
- Verdict: gain against the control, threshold τ from k (skill §6). An empty `at_bound` list at the end means the
  data, not the floor, chose the width.

## 3. Cross-checks that make the verdict believable

a. The same recipe on the other plausible hypotheses: widen the existing neighbour instead (a floor on its width). In
   the worked case a 3/2⁻ neighbour pinned at 400 and 500 keV was worse than the control by 2 349 and 9 938.
b. A data-only window multi-start on the complaining data set: ±0.5 MeV, normalizations fixed, random a₀ sign
   patterns, every J^π with a level, 250 evaluations each (about 1 h on one node). It is a screen that must agree on
   J^π and width; never a start vector. In the worked case it tied 1/2⁻ and 5/2⁻; the full data separated them
   (1/2⁻ gain 10 950, 5/2⁻ loss 14 965).
c. Attribution of the gain per data set and per angle (worked case: 72 % in one (α,n₀) set, all at backward angles,
   the symptom the evaluator named).

## 4. Full free refit of the accepted fit

Every parameter free, no floor, 400 evaluations. The level must hold on its own; that export is then the incumbent.

## 5. Propagate

Any model that freezes this block gets a new seed with the refreshed levels: copy energies and channel `gamma` fields
by (J^π, energy order), add the new levels fixed, keep the fixed flags, write a partial `.sav`, check the objective
frozen, recentre the policy windows, then refit. A splice cost is normal, because the free levels outside the block
were tuned to the old one.

## 6. Then the next region

Take it from the band map with the same recipe. Expect most J^π to be rejected quickly (high-spin groups by 10⁴–10⁵
in χ²: a 100 keV floor is above their Wigner limit). A candidate that passes on a data set other than the one aimed at
still counts.

Cost: 11 J^π take 4 node jobs of 6–8 h; the anatomy and the window screen are login-node work.
