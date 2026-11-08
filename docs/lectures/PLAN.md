# R-matrix lectures (plan)

Requested by the PI on 2026-11-08: a lecture series on R-matrix theory as AZURE2 implements it, in Markdown that can
be pasted into Notion (one `.md` per lecture, figures as PNG beside it), later also as slides (the `research-slides`
skill, from the same material). Each lecture says where in AZURE2 the quantity lives (file, function, upgrade part)
and how it is computed, with figures computed independently (mpmath) of the engine. Register: `docs/upgrade/STYLE.md`
(every symbol with its unit and frame, numbers not adjectives).

## Series

1. **The external region: Coulomb and Whittaker functions** — the radial equation outside a; F, G, η, ρ, turning
   point; P, S, φ (open channels); the Whittaker function W_{−η_b,l+½}(2κr) of bound channels and S below threshold,
   continuous through E = 0; the ANC ↔ γ conversion and N_f; the external-capture integrand and its peripherality;
   how AZURE2 computes them (`physics/coulomb`: coulsm/COULCC; `whittaker` via the integral representation of U,
   `docs/upgrade/coulomb.md` §1.4). **Status: figures done** (`01-external-region/figures/`, made by
   `scripts/lecture01_figures.py OUTDIR`, mpmath, ~2.5 min); **text to write** (`01-external-region/lecture.md`).
   Figures: (1) F_0, G_0 for η = 0 and 2 vs ρ, turning point; (2) P_l(E), p + ¹²C, l = 0–3, a = 4.11 fm;
   (3) S_l(E) of p + ¹⁶O, l = 0–2, −3 to +3 MeV through threshold, ¹⁷F 5/2⁺ and 1/2⁺ binding marked; (4) bound tails
   u² = C²W² of ¹⁷F 5/2⁺ (B = 0.600 MeV, l = 2, C² = 1.08 fm⁻¹) and 1/2⁺ (B = 0.105 MeV, l = 0, C² = 6490 fm⁻¹;
   Gagliardi et al., PRC 59, 1149 (1999)), and the fraction of the external norm beyond r (1/κ = 6.1 and 14.5 fm);
   (5) E1 external-capture integrand p + ¹⁶O → ¹⁷F at E = 0.3 MeV, l_i = 1, and its running integral.
2. **Channels, the level matrix and the collision matrix** — R and A matrix, boundary conditions B_c, the
   formal/Brune/Park bases (`level_parameters.md`), observed widths Γ = 2Pγ².
3. **External capture** — channel and non-channel capture, γ^ext, the `--ignore-externals` mask
   (`external_capture.md`).
4. **Observables** — cross sections, angular distributions, analysing powers, coupling coefficients.
5. **The apparatus fold and χ²** — target effects, resolution, normalizations, energy shifts.
6. **Fitting and uncertainties** — gradients (adjoint vs finite differences), bases in a fit, MCMC.
7. (later) THM in AZURE2.

Deliverables per lecture: `NN-topic/lecture.md` (Notion-ready: images referenced as `figures/*.png`), the figure
script in `scripts/`, then a slide deck.
