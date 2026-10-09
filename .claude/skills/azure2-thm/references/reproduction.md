# Reproducing a published THM analysis (Stage 7)

Reproduce before replacing. A reproduction separates what the engine does from how the
paper was read.

1. **Triage** (open sources only; do not bypass paywalls or captchas). Check whether
   the paper gives a full parameter set with signs, the resolution, the data, the
   fitted curve or only a band, the formula, the radii and B_c. Most papers lack at
   least one item. La Cognata et al. ApJ 723 (2010) was the only complete multilevel
   case.
2. **Collect** the formula (and its variants), formal against observed tables, radii,
   B, carrier, KF convention, σ, data windows, normalisation region and background.
3. **Digitise** the figure (`provenance-checklist.md`, "Digitisation"). Record the
   digitisation error.
4. **Compare.**
   - Allow exactly one free normalisation over a stated range.
   - Report rms and maximum in % (rms dex for log data), the peak-region values, and the
     fraction of points inside the published band.
   - Compute the floor first (for example a free Lorentzian per level and channel:
     12C+12C 0.8–1.9 % rms, 0.11 dex), so you know what "reproduced" can mean.
   - Evaluate the paper's formula independently in Python, to separate engine from
     interpretation (18O: AZURE2 against Python ≤ 0.04 %).
5. **Diagnose a non-reproduction**, in this order:
   1. Internal consistency: ΣΓ_c = Γ; formal against observed tables (ApJ 723: Γ_p1
      8.2 against 11.1 keV).
   2. Column permutations of width tables (all 24 for four columns). For 12C+12C the
      Γp0/Γp1 swap beat the next permutation by ×11 in χ², and per-level swap tests
      showed a transposition rather than scattered typos.
   3. Formula variants (a literal equation with an extra factor E against the
      OES-shaped reading; phases in radians or degrees).
   4. A missing or different background.
   5. Representation: Brune or formal B; `perlevel` or `constant` vertex.
   6. Frames (E_lab = 2E_cm for identical masses) and the (1+δ) factor for identical
      entrance nuclei.
6. **Verdict wording.** Write REPRODUCED, PARTIALLY REPRODUCED or NOT REPRODUCED,
   followed by the exact condition. Examples:
   - "only with the Γp0/Γp1 columns exchanged";
   - "with the constant vertex at B = S(E₁), 5 % peak rms, while the authors used a
     fitted complex ratio L21".

   Then give a numbered list of what the authors could publish to close the gap. Phrase
   it as help, never as criticism.

Reference reproductions:

- `tests/18O_p_a_thm`: ApJ 723 band, `vertex=constant` 10 % / 5 % rms (overall / peak),
  `perlevel` 23 %.
- 12C+12C (Tumino et al. 2018): the published table is reproduced only with its Γp0/Γp1
  columns exchanged. Integrating the published S* (four channels) reproduces their rate
  ratio to CF88 to 1–3 % at T₉ = 0.5–1.0. `examples/c12c12_tumino2018` is a joint
  refit with direct data, not the published table.
