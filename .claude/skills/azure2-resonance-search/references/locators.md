# Locators: measured performance and open improvements

Numbers that calibrate the three locators of the azure2-resonance-search skill (§1 residual ranking, §1a peak finder,
§1b blur map). All come from one dense spectrum, the 13C+α compound system (¹⁷O, E_x 4.1–13.4 MeV, 51 data sets,
a fitted model of 143 levels), and should be re-measured on a new reaction before they are quoted.

## Peak finder (§1a): recall against chance

Truth is the fitted level set; chance is the same candidate list shifted by 50–300 keV in E_x, scored at the same
tolerance. "Within" is the tolerance in E_x.

| level width Γ | within | found (any) | chance | found (≥ 2 data sets) | chance |
|---|---|---|---|---|---|
| 1–10 keV | 5 keV | 11/14 | 4.7 | 10/14 | 2.9 |
| 10–50 keV | 10 keV | 12/16 | 8.1 | 9/16 | 6.1 |
| 50–150 keV | 20 keV | 21/25 | 18.9 | 19/23 | 15.2 |
| ≥ 150 keV | 40–80 keV | all | ~all | all | ~all |

- Recall for Γ ≲ 50 keV is well above chance. For broad levels the chance rate is already near 100 %: a 300 keV level
  overlapping three others is not a peak.
- Precision cannot be judged against a fitted level set (the model is not the truth). The finder gave 293 candidates
  for 143 levels, many on the shoulders and interference lobes of real levels.
- Two level-free multi-set candidates (E_x 9.257 and 11.219 MeV) coincided with level-free residual regions of the
  §1 ranking (9.24–9.32 and 11.18–11.23 MeV). Agreement of independent locators is evidence.

## Blur validity (§1b)

rms of blur(full) − blur(broad-only) at D = 200 keV, per data set. Below ~20 % the blur leaves broad structure; above
~25 % the narrow levels carry the strength and a blurred deviation can be missing narrow area.

| data set | validity number | use |
|---|---|---|
| Cierjacks neutron total | 7 % | valid |
| Kelly | 7 % | valid |
| secondary γ sets | 1–10 % | valid |
| Si (α,α) | 12 % | valid |
| MANA and 5U (α,n₀) above E_x 10 MeV | 14–16 % | valid |
| Heil (α,α) | 17 % | valid |
| Brandenburg | 25 % | not valid |
| BandH | 27 % | not valid |
| ND 2021 (α,n₀), E_x 7–10 MeV | 43 % | not valid |

On this system the useful scale was D = 100 keV; at D = 200 keV the range E_x 4–10 MeV merged into one region. The
first deviation found was a forward-angle (α,n₀) deficit at E_x 10.0–10.8 MeV (ND 2021 at 0°: blurred data about twice the
blurred fit; MANA at 17°) on the same E_x as residual regions 10.05–10.46 MeV: a broad-level candidate.

## Global sign table

Σ fit / Σ(n·data) per data set was 0.86–0.90 for the (α,n₀) angular sets and 1.00 for the neutron total, with the
median point 1 % off. A fit with data-proportional errors undershoots the peaks; report this table with the blur map.

## Open improvements

- A broad-structure estimator on integrated data: multi-scale curvature, or local Breit–Wigner plus background fits.
- Lower bounds on J per candidate from the Legendre order of the residual (§2c of the skill).
- An automatic ENSDF cross-match of candidates (`pyazr.nds.fetch_levels`).
- A test on a sparse spectrum, where precision can be measured.
