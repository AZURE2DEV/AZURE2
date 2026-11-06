# Data provenance checklist (Stage 1)

Fill one record per THM dataset before any modelling. Keep the records with the
evaluation notes, outside the repository. Write the source in the `#` header of each
data file.

## Record template

| item | value | source (paper, section, table/figure) |
|---|---|---|
| reaction, Trojan horse a = x + s, spectator s | | |
| B_xs (MeV) | | |
| beam, beam energy (lab MeV), normal/inverse kinematics | | |
| spectator-momentum cut \|p_s\| (MeV/c), accepted angles | | |
| quantity tabulated: HOES yield / penetrability-corrected σ or S / fitted S | | |
| what was divided out: KF₃b, λ₃/λ₂, \|φ\|² only, La Cognata formula | | |
| penetrability correction: l values, radius, one constant or one per l | | |
| energy frame of the points (lab/c.m.), bin width | | |
| resolution: σ or FWHM, frame | | |
| errors: statistical, angular-integration, systematic; correlated? | | |
| subtracted background (formula, tabulated?) | | |
| normalisation window and reference data | | |
| channel radius, boundary value, vertex form of the authors' HOES fit | | |
| identical nuclei, odd-J levels | | |
| digitisation: method, error, cross-check | | |

## Telling HOES from on-shell-equivalent

Signs of a **HOES yield** (fit as a THM segment, code + 10):

- points below the entrance threshold (a penetrability-corrected σ vanishes there);
- the excitation function stays finite or rises toward threshold;
- narrow high-l levels appear as peaks of height comparable with broad low-l ones
  (17O: the f-wave 5⁻ at 0.6× the p-wave 3⁻; on shell it would be invisible);
- arbitrary units, and a fit "by the modified R-matrix".

Signs of an **on-shell-equivalent** set (fit as an ordinary segment with a free norm):

- "penetrability effects included / corrected", "normalised to direct data";
- absolute units (mb, MeV b); an EXFOR quantity such as SIG or SFC with DERIV;
- the points fall like e^{−2πη} toward threshold.

If the text is ambiguous, fit both readings, each with a free scale. Compare the THM
χ²/N, the direct norms and the residuals near threshold.

| case | on-shell reading | HOES reading | verdict |
|---|---|---|---|
| 7Li Tumino 2006 | THM χ²/N 3.2, direct norms 0.91–0.95 | χ²/N 9.8–12.8, norms 0.50–0.62, −10 to −14σ below 0.4 MeV | on-shell (paper text agrees) |
| 6Li Pizzone 2011 | χ² 102.8/62, free scale 1.11 | χ² 61.2/62 | on-shell by the paper text, although χ² prefers HOES |
| 15N La Cognata 2007 Table 3 | scale 0.96–1.01 | HOES × P₀e^{2πη}: THM χ² 85.5 vs 37.5 | on-shell |

The paper text decides, not the χ². A per-l renormalisation can leave a shape that a
free HOES model follows easily.

## Carrier and binding energy

- Read the carrier from the experiment section, not from earlier papers of the same
  group. 7Li(p,α) Tumino 2006 and 6Li(d,α) Pizzone 2011 both used ³He
  (B = 5.4935 MeV). The deuteron set-up (2H(⁷Li,αα)n) belongs to Lattuada 2001.
- Common values: d → p + n 2.224566 MeV; ³He → d + p 5.4935 MeV; ¹⁴N → ¹²C + d
  10.2723 MeV.
- Check the kinematics yourself (which body is beam and which target). The quasi-free
  energy and the reachable spectator momenta depend on it. THM Workspace >
  Experiments prints B(x+s) from the masses and the quasi-free energy.

## Digitisation

1. Prefer the PDF vector paths (exact). Otherwise use a tick-calibrated raster and check
   every major tick to ≤ 1 px.
2. Record the digitisation error, typically ~1 % or 0.03 dex. Cross-check against any
   independent digitisation or replot.
3. Error bars:
   - Read both ends. Bars asymmetric by > 30 % usually mean the marker or a curve
     overlaps one end; use the clean end.
   - A bar hidden by the marker is shorter than the marker radius. Cap it there
     (19F: ROOT starts a bar 0.0092 from the marker centre). Do not estimate it from a
     fraction of y without checking the trend of the visible bars.
   - Bars of cap size: add ~5 % in quadrature.
4. Look for a published table first. A table believed private may be open: the JUNA
   19F(p,αγ) S factors are Table I of the accepted manuscript. One digitised point
   was 17 % low against it.

## Error bars against scatter

- Fit a smooth model (one Gaussian per peak plus a polynomial) with the published bars.
  χ²/ν ≫ 1 means the bars are statistical only, or the angular-integration error is
  missing.
- 19F (Su 2025): χ²/ν ≈ 9 with the bars, against the quoted 2.0, which needs ~15 %
  errors. Handle it by scaling the statistical errors with √(χ²/ν) of the THM points
  and rescaling the AIC weights. Keep the file unchanged.
- χ²/N ≪ 1 (12C+12C, 0.24) means generous bars: AIC weights are not meaningful, so
  quote flat weights.

## Direct data

- EXFOR STATUS "CURVE" means read from a figure. DATA-ERR may cover only some points.
  Say in the header when per-point errors are assigned (for example "5 % assigned").
- Convert units once (mb → b, keV → MeV) and the frame once. AZURE2 data are lab
  energies, angles and cross sections. Check labels: an EN-CM column labelled MeV held
  keV (15N, C1788004).
- S → σ conversions: use the same masses as the model (integer against nuclear masses
  differ by up to 3.6 % at 78 keV for p+7Li).
- An identical exit pair (α+α) gets no (1+δ) factor in the model. Data that count both
  particles need ×½. A free norm of ~1.9 on one old dataset (Cassagnou) reflected a
  ×1.5 disagreement between datasets, not a double count. Keep such a norm free and do
  not anchor on it.

## Background

- Digitise the authors' subtracted background (for example the dashed line of
  ApJ 723 Fig. 4: 2186 E − 1109) and note it in the header.
- Never reuse a background fitted by one of your own variants for the others. In the
  old 18O 2010 example file, the background was a joint-fit line tuned to pw, 5.1 fm.
  It made the lowest point +200 %, three points negative, and the 6.1 fm variant
  χ² 2866 instead of 540.
- Let `background=linear` on the experiment line take up the residual smooth term.

## What to ask authors to publish (for the report, phrased as help)

1. The HOES yield after the division by KF·|φ|² and before any penetrability
   correction, with the per-l constants if a correction was applied.
2. Carrier, spectator, B, beam energy, momentum cut and accepted angles.
3. The KF convention, and whether |φ|² was divided out.
4. Resolution (σ or FWHM), frame and bin width.
5. Statistical and angular-integration errors separately, with their correlation.
6. Channel radius, boundary value and vertex form of the HOES fit, and the radius of
   the HOES-to-on-shell conversion.
7. R-matrix parameters with signs, B_c, radii and l values.
8. Coulomb term, acceptance, distortion and line shape used, if any.
9. Normalisation window and reference data.
10. The subtracted background, tabulated.
11. The treatment of identical nuclei and odd-J levels.
12. The fitted curve and the data, tabulated.

## What the papers we checked state

| paper | radius | boundary / vertex | KF | resolution |
|---|---|---|---|---|
| La Cognata PRL 2008, ApJ 708 (18O narrow) | not stated / symbol only | Breit–Wigner with M_i(E) | three-body, formula not written | 40 keV FWHM = σ 17 |
| La Cognata ApJ 723 (2010) (18O doublet) | 5.1 fm p, 5.7 fm α | B = S(E₁); fitted complex ratio L21, no M_l | not written | σ 17 keV |
| La Cognata 2006/07 (15N) | symbol only | PWIA + penetrability | "KF" | not quoted |
| Gulino 2013 / Guardo 2017 (17O) | 4.1 fm for the penetrability only | outgoing-wave (on-shell) form | Monte Carlo KF·\|φ\|² | σ 20/30 keV (fitted 22.2 ± 1.2) |
| Tumino 2006 (7Li) | none | PWIA + penetrability | eq. 4 (p_c³ misprint) | 80–120 keV, convention not stated |
| Pizzone 2011 (6Li) | none | PWIA | λ₃/λ₂, on-shell p_ax | none |
| Tumino 2018 + Reply (12C+12C) | not in open material (¹²C+d ≈ 3 fm) | per level (Tumino 2021 eq. 51) | not stated | σ 30 keV |
| Su 2025 (19F) | Woods–Saxon single-particle conversion, no R-matrix radius | PWIA | raw d³σ/\|φ\|² (`triple`) | not recorded here |

No primary analysis includes the Coulomb integral term C_l, and none tabulates its fitted
curve. The radius is therefore a variant (Stage 6), not an input.
