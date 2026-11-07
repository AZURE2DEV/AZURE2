# Results from older builds

What builds before a given date computed differently, for anyone reusing archived outputs or running a frozen
binary. Each note expires; after that, the archive is assumed to be recomputed. The current behaviour is in the
Sphinx docs (`docs/source/getting_started/whats_new.rst` and the pages named below).

## External capture (M1)

- builds before 2026-10-23 (expires: 2027-02-01): M1 EC 6.25 % low in every mode; with `--no-long-wavelength` the M1
  EC of p+¹⁶O was × 0.31 (≈ 0 for N = Z); the `intEC.dat` signature ignored the switch.
- builds before 2026-10-25 (expires: 2027-02-01): wrong sign of the M1 spin term when the partner has half-integer
  spin (p+¹⁴N s-wave M1 × 3.3 or × 0.56 by channel spin; for two half-integer spins, p+¹³C, p+¹⁵N, p+⁷Li, the whole
  spin term had the wrong sign against the orbital term). Spin-0 targets were right. Refit capture fits on such systems
  (one p+¹⁴N model: χ² 43,130 → 43,082).
- builds before 2026-10-26 (expires: 2027-02-01): an M1 orbital term was added on channel-spin-flip pathways
  (s_f ≠ s_i), affecting M1 capture into p- or d-wave states of two-spin systems (p+¹⁴N, p+¹³C, p+¹⁵N, p+¹⁰B, p+⁶Li,
  p+²³Na); s-wave capture was unaffected (one p+¹⁴N model: χ² 43,082 → 43,065). Tests: `tests/ec_m1_exact`.
- A frozen old build reuses its own `intEC.*` files; the current build refuses them by signature and recomputes.

## Capture analyzing power and UPOS

- builds before 2026-10-26 (expires: 2027-02-01): the capture A_y of a target with spin had the channel-spin
  off-diagonal b_k terms sign-flipped (coupling order); spin-0 targets were right. Redo capture-A_y fits on targets with
  spin (`tests/capture_ay_coupling`). The capture route also returned 1.5 × the Madison A_y for a deuteron projectile;
  non-spin-½ projectiles are now refused.
- builds before 2026-10-26 (expires: 2027-02-01): the UPOS secondary-γ residual spin was taken from the decay index,
  wrong when a pair of another compound nucleus, a β-delayed pair or a pair outside the R-matrix J-groups is listed
  first (`tests/upos_pair_index`).
- builds before 2026-10-27 (expires: 2027-02-01): UPOS multipole mixing was skipped at any order where R_k^L'L'
  vanished (e.g. a 3/2 → 3/2 M1/E2 secondary γ at k = 2), and channels with l > 10 were truncated silently
  (`tests/mixing_guard_lorder`).

## Composite segments and frames

- builds before 2026-10-26 (expires: 2027-02-01): lab-frame differential SUM segments added c.m. terms of other exit
  pairs without J_base/J_term (−2.5 % to +6 % on one ¹³C(α,n₀+n₃) case; `tests/sum_lab_frame`).
- builds before 2026-10-27 (expires: 2027-02-01): RATIO denominators of another exit pair took the base segment's
  Jacobian (3 % on n₀/n₃ at 0° or 180° or with `relativistic=1`; 31 % for capture over elastic), and the lab output of a
  RATIO test segment in mode 3 was divided by J_base (17–20 % low) (`tests/ratio_frame`).
- builds before 2026-10-27 (expires: 2027-02-01): a SUM/RATIO component to another exit pair integrated its
  target-effect sub-points at the base pair's c.m. angle (0.5 % in one case; `tests/component_subpoint_angle`).
- builds before 2026-11-01 (expires: 2027-02-01): a SUM `Scaling` kept six decimals (1.5 × 10⁻⁷ became 0), a
  non-numeric scaling became 1.0, a zero RATIO denominator gave 0, and old-GUI component lists without angles were
  misread (some 2025 neutron-total lines lost components).
- builds before 2026-10-27 (expires: 2027-02-01): total capture written as exit key −1 with data type 0 read past the
  pair list (a crash, or meaningless numbers).

## Phase shifts, excitation energies, extrapolation grids

- builds before 2026-10-27 (expires: 2027-02-01): a phase-shift segment of a target with spin reported half the
  argument of the summed diagonal S elements (34° off in the test case; `tests/phase_channel_spin`).
- builds before 2026-10-28 (expires: 2027-02-01): the E_x output column lacked the entrance pair's E*
  (`tests/excitation_energy_column`).
- builds before 2026-10-26 (expires: 2027-02-01): a `param.sav` in mode 3 or a rate run shifted test grid k by data
  segment k's fitted shift; S "at E" from such files is at E + shift (`tests/extrap_sav_segments`).

## Rates

- builds before 2026-10-26 (expires: 2027-02-01): mode 5 splined σ against E (non-resonant charged-particle rates low by
  0.8 % at T9 = 0.01 for ¹²C(p,γ) and 1.2 % at T9 = 0.1 for ¹²C(α,γ)); set σ = 0 below E_c.m. = 1 keV (neutron 1/v
  rates 49 % low at T9 = 0.01, 14 % at 0.03, 2.8 % at 0.1, 0.6 % at 0.3); placed resonance nodes at the `.azr` levels,
  so a narrow level moved by more than ~64 widths by the parameter file could be lost; used k_B = 8.6171 × 10⁻²
  instead of 8.617333262 × 10⁻² MeV/GK (+0.1–0.3 % on narrow resonances). The current build agrees with an independent
  integral of the extrapolated σ to 10⁻⁴ (narrow level) and 2 × 10⁻⁵ (smooth) (`tests/reaction_rate`).

## pyazr

- builds before 2026-10-07 (expires: 2027-02-01): two sessions in one process shared the EC amplitude cache (wrong
  capture σ, up to 28 %, in the first session once a second existed). `angle-integrated-E1`/`-E2` segments returned the
  total from `calculate_rwa`.
- builds before 2026-10-08 (expires: 2027-02-01): free-vector entry points accepted a short vector silently.
- builds before 2026-11-03 (expires: 2027-02-01): the `.azr` was re-parsed on every call.
