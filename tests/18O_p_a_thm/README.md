# 18O_p_a_thm — reproduction of a published THM R-matrix analysis

M. La Cognata, C. Spitaleri, A. M. Mukhamedzhanov, *Effect of high-energy
resonances on the 18O(p,α)15N reaction rate at AGB and post-AGB relevant
temperatures*, ApJ **723** (2010) 1512 — the 18O(p,α)15N THM S(E) across the
660 / 799 keV 1/2+ doublet of 19F, from 2H(18O,α15N)n.

This is the only case in `tests/` that is checked against a *published* THM
R-matrix curve rather than against AZURE2's own earlier output. It is where the
`vertex=` default of the `<thm>` block was decided (commit b351b90; see
`docs/source/theory/thm_implementation.rst`).

## What the project contains

* **Levels**: Table 3, row R, B_c = S_c(E_1) (formal parameters), two 1/2+
  levels; channels p+18O l=0, a = 5.1 fm (entrance, deuteron carrier,
  B = 2.224566 MeV in field 32) and α+15N l=1, a = 5.7 fm. AZURE2 takes Brune
  parameters, so the formal set was converted (Brune 2002: Ẽ is an eigenvalue of
  E_λδ + Σ_c γγᵀ(B_c − S_c(Ẽ)), γ̃ = eigenvector · γ) and entered as reduced-width
  amplitudes (`gammaIsRWA` = 1) at Ex = Ẽ + S_p (S_p = 7.99360, S_α = 4.01380 MeV,
  AME2020):

  | level | E_λ formal | γ_p | γ_α | Ẽ (Brune) | γ̃_p | γ̃_α |
  |---|---|---|---|---|---|---|
  | 1 | 0.609 | −0.407 | 0.212  | 0.60900 | −0.40700 | 0.21200 |
  | 2 | 0.845 |  0.495 | 0.0706 | 0.81220 |  0.44004 | 0.09664 |

  (MeV, MeV^½; E relative to p+18O.) B_p = S_0(0.609) = −1.06512,
  B_α = S_1(4.589) = −0.60476. The B = S(E_2) row converts to the same Brune set
  to < 0.1 %, and AZURE2's on-shell σ(p,α) from these inputs equals an
  independent formal R-matrix to 1e-5.
* **Brune on** (the CLI default) and **no `<thm>` block**, i.e. the default
  `vertex=constant`: the entrance vertex (B_c − 1) j_l − ρ j_l′ with the channel
  constant B_c = S_c(E_1) after the level sum. Since γᵀAγ is the same matrix in
  the Brune and formal representations, this equals a no-Brune run of the formal
  set (checked to 2e-11 in b351b90) — which is what the authors computed.
* **Resolution**: 17 keV Gaussian σ (c.m.; 40 keV FWHM in the companion papers),
  folded by `<targetInt>` on both segments; its σ is lab: 0.017952 MeV.
* **Segment 1** (`data/lc723_thm_points.dat`, 39 points): the THM S(E) points
  of Fig. 4, digitized, with their digitized error bars (mean of the two arms),
  turned into the HOES observable AZURE2 computes:
  σ_HOES ∝ (S − background) / (P_0(E) e^{2πη}).
* **Segment 2** (`data/lc723_band_mid.dat`, 79 points, 0.505–0.895 MeV): the
  mid-line of the published band, converted the same way, band half-width as
  error. `lc723_band.txt` holds the same band in S(E) with the background and
  the conversion factor; `check.sh` reads it.

Energies in the data files are proton lab energies; each segment's THM scale is
profiled out by the engine.

## The pins

`expected/chiSquared.out`: chi2 = 2130.08 (segment 1, 39 points) + 98.86
(segment 2, 79 points) = **2228.94**, a pin, not a fit. Segment 1 is large
because the digitized error bars are only 3–5 % while the points scatter
by more than that about the smooth curve, and the model has no freedom here.

`check.sh` encodes the verdict of the reproduction. It folds nothing itself:
it takes AZURE2's folded segment-2 curve, multiplies by P_0 e^{2πη}, fits one
normalization to (mid-line − background), adds the background, and measures the
relative deviation from the mid-line:

| run | rms 0.505–0.895 | rms 0.56–0.84 | max | in band | asserted |
|---|---|---|---|---|---|
| default, `vertex=constant` | 10.2 % | 5.2 % | 31 % | 0.71 | rms ≤ 12 %, peak ≤ 7 % |
| `vertex=perlevel` (default before 2026-09-25) | 26.8 % | 23.1 % | 55 % | 0.19 | peak ≥ 15 % |

The reproduction measured 10.1 / 5.1 % (no-Brune, Python folding of S) and
26.8 / 23.1 % (Brune per-level, AZURE2 folding) — folding σ in AZURE2 instead
of S in Python moves these by < 0.4 points.

Re-pinned 2026-09-28: the `<targetInt>` line lists both segments ("1,2"),
and AZURE2 converted the one shared effect's sigma lab -> c.m. once per listed
segment, so both segments were folded with 17 x 0.947 = 16.1 keV instead of
17 keV. Each listed segment now gets its own copy of the effect, converted
once. 2191.98 -> 2228.94 (segment 1 2095.59 -> 2130.08, segment 2 96.38 ->
98.86); the table above moved by less than 0.4 points (it was 10.2 / 5.2 %,
30 %, 0.72 and 26.8 / 23.1 %, 55 %, 0.15) and both verdicts stand.
`check.sh` also asserts that the "1,2" line gives output identical to one
line per segment.

## Verdict

With B_c = S(E_1) applied after the level sum, AZURE2 reproduces the published
band in the peak region at 5 % rms (77 % of it inside the band in the
reproduction). The per-level S_c(E_λ) vertex under Brune — the same on-shell
physics — puts the 799 keV peak about 40 % low: for interfering levels of one
Jπ it is not the vertex of any single R-matrix representation.

## Caveats (from the reproduction)

* **Fig. 4 is S(E), not the HOES cross section.** The conversion S = σ_HOES ·
  P_0 e^{2πη} (constant E/k² dropped) is the reading that makes the authors' own
  eq. (14) with their published φ21 fit best; the literal eq. (14) (extra
  factor E) would need a different φ21.
* **Only the band is published**, not the central curve: comparison is to the
  band mid-line, with the band half-width (±7–10 %) as the natural tolerance.
  Digitization error (~1 %) is negligible.
* **The authors' formula is not AZURE2's.** They replaced the vertex by a fitted
  complex ratio L21 = (γ1p/γ2p)·m21·e^{iφ21} (m21 = 1.09, φ21 = −0.32 rad) and
  added a linear background (S = 2186 E − 1108.7 MeV b, digitized, not given in
  the text). AZURE2 has no free complex vertex ratio; the agreement relies on the
  energy dependence of M_0(p) imitating it. The residual misfit is the low-energy
  flank, 0.50–0.57 MeV, +20–30 % (`coulombIntegral=1` removes most of it:
  5.4 % rms overall in the reproduction).
* **The paper is internally inconsistent for level 1**: its formal set gives
  Γ_p1 = 8.2 keV and Γ_1 = 190 keV against the quoted 11.1 ± 1.1 keV and
  199 ± 3 keV (Table 4). Level 2 agrees. The test uses the formal Table 3 set.

## Regenerating

The inputs come from an independent reproduction of the paper that is not part
of this repository (the digitized Fig. 4 in `digitize/`, the Coulomb functions
and the formal→Brune conversion in `model/rmat.py`, the level lines from
`azure/build_azr.py`); `THM_REPRO` names its directory:

    THM_REPRO=/path/to/reproduction python3 tests/18O_p_a_thm/make_inputs.py

(numpy, scipy, mpmath). Then re-pin `expected/chiSquared.out` from a run of the
project (menu 1, as `tests/run_tests.sh` does) and re-measure the table above
with `check.sh`. The test itself needs only bash and awk; the pair takes a
second.
