# T2 — p + ¹⁵N → ¹⁶O through an apparatus

T2 is a synthetic R-matrix project whose data were computed from known level parameters and known apparatus
parameters, then given noise. A fit of T2 must return those parameters within its own uncertainties. It exercises
the second and third lines of the calculation (`docs/upgrade/ARCHITECTURE.md` §1): the apparatus kernel K(ν) —
target integration with stopping and straggling, beam resolution, beam profile, Doppler broadening, energy-scale
shifts, composite observables, a natural target, the Trojan-horse (THM) fold — and the comparison with data through
normalizations, priors and a data covariance. It is the project on which the upgrade's apparatus kernels, priors and
THM observable are checked (`ARCHITECTURE.md` §9–§10), and the benchmark of `scripts/bench_upgrade.sh`.

## 1. The system

Compound nucleus ¹⁶O. Three particle pairs (atomic masses in u; separation energies of ¹⁶O in MeV; channel radii
a_c in fm):

| key | pair | J^π of the two particles | S (MeV) | a_c (fm) | note |
|---|---|---|---|---|---|
| 1 | p + ¹⁵N | ½⁺, ½⁻ | 12.1274 | 4.85 | entrance; carries B(d → p + n) = 2.224566 MeV for THM (columns 32–33) |
| 2 | α + ¹²C (g.s.) | 0⁺, 0⁺ | 7.16192 | 5.43 | α₀ |
| 3 | α + ¹²C (4.43891 MeV, 2⁺) | 0⁺, 2⁺ | 7.16192 | 5.43 | α₁; opens at E_x = 11.6008 MeV |

Masses: p 1.007825, ¹⁵N 15.000109, α 4.002603, ¹²C 12 (AME2020). Lab and c.m. energies of pair 1 are related by
E_lab = E_c.m. (m_p + m_15N)/m_15N = 1.06719 E_c.m. (non-relativistic, as the engine).

**Levels (the truth).** Brune parameters (C. R. Brune, Phys. Rev. C 66, 044611 (2002); the CLI default): the
observed excitation energy E_x (MeV) and the observed partial widths Γ_c (eV, at E_x), the sign of Γ being that of
the reduced-width amplitude γ_c (MeV^½). The two 1⁻ levels and the two 3⁻ levels are those of `tests/15N_p_a`
(fitted there to Schardt et al., Phys. Rev. 86, 527 (1952), and La Cognata et al., Phys. Rev. C 76, 065804 (2007));
the α₁ widths are invented for T2. Channels are written (s, l): channel spin and orbital momentum.

| level | J^π | E_x (MeV) | Γ_p (s=1, l=0 or 2) | Γ_α0 | Γ_α1 (s=2, l=1) | free in the fit (start value) |
|---|---|---|---|---|---|---|
| 1 | 1⁻ | 12.4208 | 496 eV (l=0) | 80.9 keV (l=1) | 50 eV | E_x (12.4258), Γ_α0 (75.0 keV) |
| 2 | 1⁻ | 13.1050 | 36.05 keV (l=0) | −689.6 keV (l=1) | 20.0 keV | Γ_p (40.0 keV), Γ_α1 (15.0 keV) |
| 3 | 1⁻ | 17.0 | 8305.81 keV (l=0) | −3369.465 keV (l=1) | 0 | background pole, fixed |
| 4 | 3⁻ | 13.0734 | 4.648 keV (l=2) | 43.94 keV (l=3) | 5.0 keV | — |
| 5 | 3⁻ | 13.2692 | 1.702 keV (l=2) | 12.92 keV (l=3) | 0 | E_x (13.2652), Γ_α0 (14.0 keV) |

The p-channel with the higher l of each level (l = 2 for 1⁻, l = 4 for 3⁻) has Γ = 0, fixed. The radii, the
background pole and the fixed widths are assumptions of the model, not determined by the data.

**Fixed widths at the start values.** What the fit holds fixed is the Brune reduced-width amplitude γ_c, while the
`.azr` file enters the observed width at the level energy, Γ_c = 2P_c γ_c² / (1 + Σ_c' γ_c'² dS_c'/dE) (P_c the
penetrability, S_c the shift function). For a level whose energy or one of whose widths starts off the truth, the
true Γ written at the start would fix the wrong γ: a 5 keV start offset of level 1 changes γ_p(1) from 0.46355 to
0.43667 MeV^½, and the fit can then not reach the truth. `make_t2.py` therefore finds, with the engine, the fixed
widths at the start values that give the true γ (to 10⁻⁶): Γ_p(1) = 547.90 eV, Γ_α1(1) = 53.309 eV,
Γ_α0(2) = −747.275 keV, Γ_p(5) = 1677.69 eV in `T2.azr` (recorded in `truth.json`). At the start only the free
parameters differ from the truth.

## 2. Observables and apparatus (13 segments, 330 points)

All data are lab-frame (the engine's input convention), E in MeV of the proton (lab), cross sections in b or b/sr.
The stopping cross section ε(E_lab) is in MeV cm² per atom of the entity the areal density N (atoms cm⁻²) counts;
x below is E_lab in MeV. Sub-point grids: segments 1, 4 and 5 use a resonance core of 10 widths with 20 points per
width (the engine's defaults are 20 and 50; with the defaults the thick target alone held 36,223 of 58,462 sub-points,
and its χ²/N moved by 7 × 10⁻³ relative between the two settings); generation and fit use the same grids, so the
setting does not bias the recovery. One evaluation has about 25,000 sub-points.

| seg | observable | apparatus (K) | comparison (ν, priors) | points, noise |
|---|---|---|---|---|
| 1 | σ(p,α₀), angle-integrated | thin TiN target: N = 5×10¹⁷ ¹⁵N cm⁻², ε = 1.6×10⁻²⁰ x^−0.45 per N atom (13.7 keV loss at 0.3 MeV); Gaussian beam spread σ_b = 1 keV (lab); Bohr straggling | norm n₁ (true 1.04; prior 5 % about 1); constant shift ΔE₁ (true +2.0 keV; prior 3 keV about 0) | 40, 0.20–1.40, 4 % |
| 2 | dσ/dΩ(p,α₀) at θ_lab = 90° | Gaussian resolution 2 keV (lab); detector attenuation Q_k = 1, 0.995, 0.98, 0.96, 0.93 (k = 0–4) | normgroup D (true 0.95; prior 8 %); shiftgroup G: ΔE (true −1.5 keV; prior 3 keV) plus shiftfn a₀(x − 0.8) (true a₀ = 0.0040; penalty 0.01); covariance `cov_diff.txt` | 30, 0.25–1.35, see below |
| 3 | dσ/dΩ(p,α₀) at θ_lab = 150° | as 2 | as 2 (same groups, same covariance) | 30, as 2 |
| 4 | thick-target yield of (p,α₀) per incident proton | N = 0 (infinitely thick: Y(E₀) = ∫₀^E₀ σ/ε dE), ε as segment 1; no beam spread | norm (true 0.98; prior 3 %); √E shift b√x (true b = 0.0015 MeV^½; penalty 0.003) | 30, 0.25–1.35, 3 % |
| 5 | σ(p,α₀) across level 5 | ¹⁵N implanted in Ta: N (true 2.2×10¹⁷; nominal 2.0×10¹⁷), ε = 6.0×10⁻²⁰ x^−0.30 per ¹⁵N atom (12.4 keV loss at the true N and 1.22 MeV); σ_b = 0.8 keV; straggling c_s = 0.6 keV^½ (c.m.) | fitted density `density_4` (prior 10 % about 2.0×10¹⁷; grids laid out to 3×10¹⁷) | 26, 1.17–1.27, 3 % |
| 6 | σ(p,α₀) across level 1 | ¹⁵N₂ gas, N = 10¹⁸ N cm⁻², ε = 0.9×10⁻²⁰ x^−0.5 per N atom (≈ 16 keV loss); free-gas Doppler broadening at T = 300 K (fixed); resolution σ = 0.004 x (lab, `sigmalab`) | norm (true 1.02; prior 5 %) | 25, 0.26–0.42, 4 % |
| 7 | σ(p,α₀) in reconstructed-energy windows | beam profile: one skewed Gaussian, ξ = 0.34 MeV, ω = 35 keV, α = 1.5 (lab); detector resolution s = 4 keV; windows of 15 keV from 0.270 to 0.420 MeV (data columns 5–6) | norm (true 0.97; prior from `prior_cov.txt`) | 10, 3 % |
| 8 | SUM σ(p,α₀) + σ(p,α₁) | Gaussian 2 keV | norm (true 1.03; prior from `prior_cov.txt`) | 30, 0.60–1.40, 4 % |
| 9 | RATIO σ(p,α₁)/σ(p,α₀) | none | norm fixed at 1 | 20, 0.70–1.40, 5 % |
| 10 | σ(p,α₀) on natural nitrogen | sample `natN`: the ¹⁶O term times the ¹⁵N abundance w (true 0.00364; nominal 0.0037, prior σ 0.0002) | norm fixed at 1 | 15, 0.25–0.45, 6 % |
| 11 | THM HOES, ²H(¹⁵N,α¹²C)n, run a | Gaussian resolution σ = 20 keV in the c.m. of p + ¹⁵N (entered as 21.344 keV lab, converted by m_15N/(m_p + m_15N)) | experiment[LC09]: one profiled scale shared with 12 | 30, E_c.m. 0.03–0.60, 7 % |
| 12 | THM HOES, run b | as 11 | as 11 | 28, E_c.m. 0.04–0.58, 7 % |
| 13 | THM on-shell-equivalent (OES) σ(p,α₀) | none (see §5) | free scale (true 1/1.35), no prior | 16, E_c.m. 0.05–0.50, 10 % |

**Covariance of segments 2–3** (`cov_diff.txt`, rows = the 60 points of the two files in order): C = D + v vᵀ,
D = diag((0.04 m_i)²), v_i = 0.03 m_i (E_i − 0.8)/0.55, m_i the noiseless model in the data's units — a 4 %
uncorrelated error and a 3 % efficiency tilt linear in energy and common to both angles.

**Prior covariance** (`prior_cov.txt`, `priorCovariance=` in `<config>`): a two-dimensional Gaussian on
(n₇, n₈), σ = (0.03, 0.04), correlation 0.5, its mean drawn once about the truth (an external normalization
measurement of both segments with one monitor).

**THM.** Segments 11–12 are half-off-energy-shell (HOES) excitation functions, observable code 10, computed with
the engine's defaults (`vertex=constant`, `kinematics=lacognata`, `entranceL=incoherent`; see
`docs/source/theory/thm_implementation.rst`), the p + ¹⁵N binding in the deuteron B = 2.224566 MeV, and the
three-body kinematics stated as in La Cognata et al. (2009): ¹⁵N on d at 60 MeV, spectator n. The two runs share
one arbitrary scale, profiled analytically (it is not a fit parameter); no background, no energy shift. Segment 13
holds on-shell-equivalent points (the THM result after division by the penetrability and normalization to direct
data), compared with the on-shell σ(p,α₀) with a free scale, as `THM.md` §2 requires: the two kinds are kept apart.

## 3. Generation and noise

`make_t2.py` contains the truth (the tables above). It writes `T2_true.azr` (every parameter fixed at its true
value) and computes the noiseless model at every data point with the engine, one run per segment (segments 2–3 and
11–12 together, since they share a norm, a shift or a THM scale), reading the c.m.
model from `AZUREOut_*.out` and returning it to the lab frame by the factor the engine applied to a placeholder
datum. Noise is Gaussian with standard deviation equal to the stated fraction of the noiseless model; segments 2–3
are drawn from the full covariance C. One random stream, `numpy.random.default_rng(20261106)`, is used in this
order: segments 2–3, then segments 1, 4, 5, …, 13, then the prior mean. Data values are the model divided by the
true normalization (the engine compares n·σ_data with σ_model). `T2.azr` is the fit project: the same data, the
free parameters at the start values of §1 and the nuisances at their nominal values.

Because the model is the engine's own, T2 tests that the parameters are identifiable through the apparatus and that
the fit returns them; it does not test that the engine's physics is right. That is the task of the identity checks
and `tests/reference` (`ARCHITECTURE.md` §10).

## 4. Recovery with today's engine

`check_recovery.py` fits `T2.azr` (menu 2 with `--use-lm`, see below) and compares each of the 19 free parameters p with its truth t by the pull (p − t)/δp, δp the fit's
parameter error from `param.sav`. Level parameters are compared in the units of `param.sav` (E_x in MeV, Brune
γ in MeV^½), their truth taken from a calculation of `T2_true.azr`. Tolerances: |pull| ≤ 3 and χ²/N ≤ 1.5 for the
committed data; |pull| ≤ 0.5 for noiseless (Asimov) data, where the optimum differs from the truth only by the pull
of priors centred on the nominal values.

Results, 2026-11-06, engine sources of commit 37978592 (`feature/upgrade`), 8 threads on an HTCondor slot:

| data | χ² (N = 330) | max \|pull\| | MIGRAD calls | fit wall time |
|---|---|---|---|---|
| committed (noisy) | 294.17 | 2.15 (`segment_2_energy_shift`) | 1082 | 585 s exact Coulomb; 72 s with a Coulomb table, same parameters to the 8 digits of `param.sav` |
| Asimov (noiseless) | 0.17 | 0.38 (`density_4`) | 1106 | 502 s |

At the truth the committed data give χ² = 311.6 for N = 330 plus 4.9 from the priors; the noiseless data give
3 × 10⁻¹⁰, which shows that generation and fit compute the same model.

Per parameter (committed data; start → fit ± δp, truth; E in MeV, γ in MeV^½, N in atoms cm⁻²):

| parameter | start | fit | truth | pull |
|---|---|---|---|---|
| E_x(1) `energy_1` | 12.4258 | 12.42093 ± 0.00018 | 12.4208 | +0.74 |
| γ_α0(1) `width_1_3` | 0.128808 | 0.13374 ± 0.00029 | 0.133815 | −0.28 |
| γ_p(2) `width_2_1` | 0.363719 | 0.3594 ± 0.0016 | 0.359443 | −0.01 |
| γ_α1(2) `width_2_4` | 0.490636 | 0.5890 ± 0.0068 | 0.589753 | −0.11 |
| E_x(5) `energy_5` | 13.2652 | 13.26954 ± 0.00018 | 13.2692 | +1.95 |
| γ_α0(5) `width_5_3` | 0.058009 | 0.05603 ± 0.00053 | 0.055691 | +0.65 |
| n₁ | 1 | 1.0380 ± 0.0077 | 1.04 | −0.26 |
| ΔE₁ (MeV) | 0 | 0.00220 ± 0.00026 | 0.0020 | +0.76 |
| n_D | 1 | 0.9537 ± 0.0063 | 0.95 | +0.59 |
| ΔE_G (MeV) | 0 | −0.00096 ± 0.00025 | −0.0015 | +2.15 |
| a₀ of shiftfn | 0 | 0.0031 ± 0.0006 | 0.0040 | −1.62 |
| n₄ | 1 | 0.9845 ± 0.0068 | 0.98 | +0.66 |
| b₄ (MeV^½) | 0 | 0.0019 ± 0.0007 | 0.0015 | +0.53 |
| N₅ `density_4` | 2.0 × 10¹⁷ | (2.083 ± 0.080) × 10¹⁷ | 2.2 × 10¹⁷ | −1.47 |
| n₆ | 1 | 1.0138 ± 0.0082 | 1.02 | −0.75 |
| n₇ | 1 | 0.9621 ± 0.0086 | 0.97 | −0.92 |
| n₈ | 1 | 1.0438 ± 0.0088 | 1.03 | +1.56 |
| w (¹⁵N in natN) | 0.0037 | 0.003700 ± 0.000055 | 0.00364 | +1.10 |
| n₁₃ (OES scale) | 1 | 0.726 ± 0.018 | 0.7407 | −0.82 |

The profiled THM scale of segments 11–12 is not a fit parameter and has no truth (the HOES data are in arbitrary
units); `thm_experiments.out` reports it.

**Which fitter.** For a project with THM segments the CLI switches on the analytic-gradient MIGRAD. On T2 that path
spent more than 3 min per gradient call at 8 threads and had not converged after 10 min. The THM segments have no
analytic residual Jacobian (the engine's `ComputeResidualJacobian` returns false), so `--use-lm` falls back to MIGRAD
with Minuit's numerical gradient, which converges (EDM 5 × 10⁻¹⁰). The recovery check, the golden fit and the
benchmark use `--use-lm` for this reason; an analytic THM Jacobian is a Phase 2/5 item.

**What the check would detect.** A kernel, prior or group that the fit evaluates differently from the generation
(the Asimov χ² at the truth would not be ~0, and pulls on the noiseless data would exceed 0.5), a parameter the
apparatus makes unidentifiable (δp large, or MIGRAD not converged), and a truth that lies outside the fit's parameter
space: before the fixed widths were matched (§1), the fit ended at χ² = 42.8 on the noiseless data, 2–11 δp from the
truth on 14 of the 19 parameters.

## 5. Not covered (TODO)

- **TODO (Phase 5, THM kernel): resolution fold of the OES points in S(E).** La Cognata et al. (2009) fold the
  astrophysical factor S(E) with a σ = 20 keV Gaussian; the engine's `<targetInt>` fold acts on σ(E), which at
  E_c.m. = 50 keV differs from the S fold by a factor of order 10² (`THM.md` §2). Segment 13 is therefore unfolded,
  in its generation and in the fit.
- **TODO: fitted Doppler temperature.** Segment 6 has T = 300 K fixed: the Doppler width of p + ¹⁵N at 0.3 MeV and
  300 K, Δ_D = (4 m_p E kT / m_15N)^½ ≈ 45 eV, is about 1800 times below the 81.4 keV total width of level 1, so the data cannot
  determine T. A fitted T needs a resonance a few hundred eV wide, which ¹⁶O does not offer here (covered by
  `tests/doppler_fgm`).
- **Natural target with a second compound nucleus.** Natural nitrogen is 99.636 % ¹⁴N, but ¹⁴N(p,α)¹¹C is closed
  below E_p ≈ 3.1 MeV (Q = −2.92 MeV), so segment 10 is one compound nucleus times its abundance. A sample sum over two compound
  nuclei is covered by `tests/natural_target`, not by T2.
- **Implanted-depth profile.** The engine has a uniform layer only; segment 5 is a uniform layer with the stopping of
  the host (Ta) per implanted atom. A non-uniform implantation profile has no engine support.
- Not in T2: SAMMY UDR resolution (neutron time of flight), transmission and self-shielded yield, electron
  screening, THM options beyond the defaults (spectator window, distortion, DW vertex, coherent background).

## 6. How to run (any machine)

```bash
export AZURE2_BINARY=/path/to/build/src/AZURE2         # a CLI build of this tree
python3 tests/upgrade/T2/make_t2.py                     # regenerate data and projects (numpy; ~60 s on 1 core)
python3 tests/upgrade/T2/check_recovery.py --threads 8  # recovery on the committed data
python3 tests/upgrade/T2/check_recovery.py --asimov     # recovery on noiseless data
python3 tests/upgrade/golden/T2/regen_compare.py        # compare with the golden outputs (tolerances inside)
scripts/bench_upgrade.sh                                 # time and peak RSS, appended to tests/upgrade/bench/T2.tsv
```

A calculation of `T2.azr` (menu 1) takes 2.4 s on one thread without a Coulomb table. Run AZURE2 by hand from a
copy of this directory: `printf '1\n\n\n7\n' | $AZURE2_BINARY --no-gui --no-readline T2.azr` (menu 1 =
calculate, 2 = fit).

Files: `make_t2.py` (generator, the truth), `truth.json` (the truth as written by it), `T2_true.azr`, `T2.azr`,
`data/` (one file per segment), `cov_diff.txt`, `prior_cov.txt`, `check_recovery.py`.
