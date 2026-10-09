# Calibration: what the cross-checks found (Stages 6–7)

Results of the full evaluations of the repository's THM examples (October 2026). Each is
a model average over full refits, quoted as mean ± stat ± model and classified by the
robustness rule. Use them to judge whether a new result is plausible and which axis to
expect to dominate.

## Key results

| reaction | quantity | value ± stat ± model | class | scale from |
|---|---|---|---|---|
| 7Li(p,α) | bare S(0) | 63.5 ± 3.1 ± 0.7 keV b | robust | direct data |
| | U_e (Li target) | 186 ± 44 ± 3 eV | robust; dataset-dependent | direct data |
| 6Li(d,α) | bare S(0) | 17.3 ± 0.9 ± 0.5 MeV b | robust | direct data |
| 15N(p,α₀) | bare S(0) | 70.7 ± 4.1 ± 0.2 MeV b | robust | direct data |
| | ωγ(312) | 0.709 ± 0.034 ± 0.009 keV | robust | direct data |
| 18O(p,α) | Γ_α1 (doublet) | 179.2 ± 2.5 ± 1.5 keV | robust | direct data |
| | Γ_p1 (doublet) | 5.17 ± 0.27 ± 0.09 keV | marginal | direct data |
| | ωγ(20)/ωγ(90) | (2.03 ± 0.20 ± 0.03) × 10⁻¹³ | robust | ratio, equal l |
| | ωγ(20) via LUNA ωγ(90) | (3.2 ± 0.5 ± 0.05) × 10⁻¹⁹ eV | robust | LUNA 90 keV |
| | ωγ(20) via the 144 keV anchor | (6.0 ± 0.8 ± 5.8) × 10⁻¹⁹ eV | model-dependent | LUNA 144 keV |
| | Γ_p(597) | 44(6) eV (DW); 0.6–9 eV (PW) | model-dependent | direct doublet |
| 17O(n,α) | ωγ(2⁺)/ωγ(3⁻) | 0.0435 ± 0.0053 ± 0.0082 | marginal | ratio |
| | ωγ(5⁻) | 93 ± 13 ± 29 meV | model-dependent | Wagemans 3⁻ |
| 19F(p,αγ) | ωγ(11) | (4.21 ± 0.64 ± 0.79) × 10⁻²⁹ eV | marginal | direct, 324 keV |
| | ωγ(213) | 0.0120 ± 0.0015 ± 0.0007 eV | marginal | direct |
| | ωγ(225) | (1.0 ± 0.4 ± 0.8) × 10⁻⁴ eV | model-dependent | direct |
| | rate/JUNA, T₉ = 0.1 | 0.81 ± 0.04 ± 0.01 | robust | direct |
| | rate/JUNA, T₉ = 0.01 | 0.21 ± 0.05 ± 0.03 | robust in the model; set by the undetermined Γ_α2(11) | direct |
| 12C+12C | E_r (7 of 8 levels), Γ (5 of 8) | σ_mod ≤ 5.3 keV | robust | – |
| | rate/CF88, T₉ = 0.5 / 1.2 | 33 / 1.51 (± 0.20/0.14 ± 0.69/0.43 dex, flat) | model-dependent | direct, 2.1–2.7 MeV |
| | Σωγ (E < 1.5 MeV) | 1.6 × 10⁻⁹ eV (± 0.31 ± 0.74 dex) | model-dependent | direct |

## Metrics of the quantities that test the reaction model

| quantity | σ_mod/σ_stat | σ_mod/\|m̄\| | leading ingredient | Δl | \|Δρ\|/d₀ | spectator |
|---|---|---|---|---|---|---|
| 15N S(0) | 0.06 | 0.3 % | radius | on shell | – | n |
| 7Li S(0) | 0.23 | 1.1 % | radius | on shell | – | d |
| 6Li S(0) | 0.56 | 2.9 % | radius | on shell | – | p |
| 18O Γ_p1 (doublet) | 0.33 (flat 5.2) | 1.7 % (27 %) | PW/DW vertex | 0 | – | n |
| 18O ωγ(20)/ωγ(90) | 0.15 | 1.5 % | none | 0 | 0.009 | n |
| 19F ωγ(11) | 1.23 (1.35) | 19 % (22 %) | PW/DW vertex | 0 | 0.23 | n |
| 17O ωγ(2⁺)/ωγ(3⁻) | 1.55 | 19 % | radius | 1 | 0.015 | p |
| 19F ωγ(213) | 0.47 (1.45) | 6 % (20 %) | vertex and acceptance | 1 | 0.08 | n |
| 18O ωγ(20) via 144 keV | 7.3 | 97 % | radius and vertex (zero of M₀) | 2 | 0.10 | n |
| 17O ωγ(5⁻) | 2.2 | 31 % | radius, then vertex | 2 | 0.03 | p |
| 19F ωγ(225) | 2.0 (8.9) | 80 % (190 %) | vertex, acceptance, doublet split | 3 | 0.07 | n |
| 12C rate, 0.5 GK | 3.5 | 0.69 dex | R(E) (×18), then radius (×10–17) | mixed | ≳ 1 | d |

## General findings

1. **Overlap with direct data makes the reaction model nearly irrelevant**
   (σ_mod/|m̄| ≤ 3 %). Flat weights raise the spread only by including variants the
   data reject.
2. **Equal l is not enough; the separation in ρ matters.** 18O 20/90 keV: 1.5 % at
   |Δρ|/d₀ = 0.009. 19F 11/324 keV: 19 % at 0.23.
3. **A different l brings in the radius and the vertex through the zeros of M_l.**
   The spread is 6–20 % for Δl = 1, 31–97 % for Δl = 2, and 80–190 % for Δl = 3. It is
   worst when the anchor sits near a zero of its own vertex.
4. **A charged spectator below its barrier: the transfer distortion can dominate.**
   For 12C+12C, R(E) leads and the radius follows; the data do not choose between plain
   and N_C × R models.
5. **The reading of the published points can matter more than any model ingredient**
   (7Li: 25–37 against 63.5 keV b).

Published THM uncertainties (3–46 %) contain no reaction-model spread. Every difference
from a published value was located in an analysis step (reading, area-to-strength
conversion, normalisation, extrapolation shape), never in the measured spectra.

## Case notes that recur

- **18O 2010 doublet.** On a fixed data treatment the vertex and radius move the doublet
  by less than stat in acceptable fits. The failing variants (pw with the 597 keV level,
  perlevel, onshell) move Γ_p1 to 9 keV, but only with direct norms of 1.2–1.5 and
  Δχ² ≥ 240. THM χ²/N stays ≥ 6 because the digitised cap-sized bars ⊕ 5 % are smaller
  than the scatter.
- **19F.** The full THM window (53 points) is incompatible with the direct 790/828 keV
  strengths in every vertex model (ωγ(828) 564–657 against 775(35); ωγ(790) → 0). The
  adopted window is E ≤ 0.45 MeV. The THM-only 213/225 doublet split is weak (ρ = −0.79).
- **17O.** The DW vertex agrees best with an independent direct set (Koehler & Graff
  0.2–20 keV: 0.88 ± 0.07, against 0.69–0.78 with PW). The THM spectrum has a minimum
  where K&G reported a 140 keV resonance, which argues against that resonance.
- **12C+12C.** The N_C × R fit reproduces the published THM rate. The plain fit, anchored
  on all direct data at 2.1–2.7 MeV, lies ×14 above it at 0.5 GK. The published
  normalisation window (α1 at 2.50–2.63 MeV) is the likely reason. Quote the rate with
  its window, its R(E) choice and ≥ 2 radii.

## Example pins

The CLI "Total Chi-Squared" (data + priors) of each example as stored. The data part is
in brackets where different; N is the number of points. The authoritative list is
*impl*, "Examples".

| example | file alone / N |
|---|---|
| `f19_pag_thm` | 79.408 (77.640) / 59 |
| `o18_lacognata2008` | 35.877 / 30 |
| `o18_lacognata2010` | 330.054 (329.692) / 161 |
| `c12c12_tumino2018` | 112.076 (111.651) / 247 |
| `o17_guardo2017_fit` | 12.996 / 23 |
| `n15_lacognata2007` | 224.425 (223.414) / 208 |
| `li7_tumino2006` | 448.575 (439.104) / 294 |
| `li6_pizzone2011` | 158.918 (156.991) / 132 |

Command: `printf '1\n\n\n7\n' | build/src/AZURE2 --no-gui --no-readline <name>.azr`, run
from the example directory, under `prlimit` and with `OMP_NUM_THREADS=1`. A different
value after a rebuild means the engine changed or the module is stale. Find out which
before fitting.
