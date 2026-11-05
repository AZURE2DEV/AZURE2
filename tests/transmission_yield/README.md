# transmission_yield

`transmission <n>` and `yield <n>` on a segment line (n in atoms/barn): the transmission
T = exp(-n σ) of a thick sample, and SAMMY's self-shielded capture or reaction yield
Y0 = (σ/σ_t)(1 − exp(−n σ_t)), σ_t the total cross section of the neutral entrance pair from
the optical theorem (a hidden `total 1` component of the segment). Under a sub-point target
effect the transform is applied at every sub-point **before** the kernel (SAMMY manual III.E.1):
a detector averages T(E), and ⟨e^{−nσ}⟩ ≠ e^{−n⟨σ⟩}.

The models are built at run time by `make_models.py`:
- from `tests/reaction_rate/models/n12C.azr`: n + ¹²C elastic and capture, with an s-wave 1/2⁺
  resonance at E_cm = 0.425 MeV and Γ_n = 34 keV. With n = 0.15 atoms/b, nσ_t ≈ 1 at the peak.
- from `tests/natural_neon/natural_neon.azr`: three compound nuclei, with natural sums of totals.

`analyze.py` does the comparisons.

| check | what |
|---|---|
| (a) | With no spread, T and Y0 match the formula applied to the plain σ_t and σ_γ segments, to 1e-9. |
| (b) | With a 15 keV (lab) Gaussian spread, T and Y0 match an independent integration, to 2e-4. The integration averages the transform of a 0.1 keV grid of unbroadened cross sections over the kernel (±5σ, normalised over the window). Transforming the *broadened* cross sections is off by 6 % (T) and 5 % (Y0) at the resonance: the check discriminates. |
| (c) | For a natural sum over three nuclei, T = exp(−n σ_nat), to 1e-9. The sample decomposition file notes that its total column is T. |
| (d) | The gradient check (`AZURE_GRAD_CHECK`) passes on a project with T and Y0 segments, with and without the spread, and every width free. The analytic gradient matches finite differences (to 1e-7 of the largest) and the Jacobian path (to 1e-10). Uses `--no-grid-replacement`, explained below. |
| (e) | Each of these is refused with an ERROR and produces no output: both keywords on one line, `yield` on the elastic total, `transmission` on a differential segment, a charged entrance pair, a natural sum under a sub-point spread, and `yield` on a sample segment. |
| (f) | On test segments (mode 3), the T and Y0 curves match the formula applied to the plain curves, to 1e-9. |
| (g) | A yield's σ_t is the total, not the elastic. The model is n + ²⁰Ne (natural_neon) with a free 1 keV width on every open reaction channel. A yield to α + ¹⁷O equals the formula with the `total 1` σ_t to 1e-9, at energies where the elastic differs from the total by more than 1e-3. |

The fitted sample thickness: `transmission n vary sigma` / `yield n vary sigma` makes n the fit parameter
`segment_<key>_thickness` (vary 0/1, sigma the Gaussian prior width in atoms/b, 0 = none), and
`thicknessgroup <label>` makes segments share the thickness of the group's first one.

| check | what |
|---|---|
| (h) | The data are the model's T (and Y0) at n = 0.18 atoms/b with 1 % errors, and every level is fixed. Starting from `transmission 0.15 1 0` (and `yield 0.15 1 0`), MIGRAD (`--use-gradient`) and `--use-lm` recover segment_1_thickness = 0.18 to 1e-4. |
| (i) | A `thicknessgroup` of the transmission and the yield: one parameter is fitted from both. The follower's fixed segment_2_thickness is written as the leader's value. |
| (j) | The gradient check runs with seven thicknesses free together with every width. The lines cover: no target effect, a Gaussian spread with a 0.01 atoms/b prior, Doppler alone (a thickness group), and Doppler plus the Gaussian. The analytic gradient matches finite differences to 1e-6 of the largest, and the Jacobian path to 1e-10. |
| (k) | A prior: `transmission 0.15 1 0.01` on 20 % data holds n between 0.15 and 0.18. chiSquared.out's `Total-Thickness-Chi-Squared` equals ((n − 0.15)/0.01)². |
| (l) | Each of these is refused with an ERROR naming the segment: vary 2, sigma < 0, sigma missing, n = 0, `thicknessgroup` without a label, and a group member without its own transmission or yield block. |
| (m) | With `REF_AZURE2` set to a build from before the fit form (7377779 or earlier, but with the transmission keyword), the plain `transmission n` / `yield n` outputs are byte-identical to it: AZUREOut and chiSquared.out of thin, conv, the test segments and the natural sum. |

An engine from c3d21f3 to 7377779 reads `transmission n` and silently ignores the fit form and
`thicknessgroup`: its keyword loop stops at the first token it does not know. On such a build the
fits of (h), (i) and (k) have no thickness parameter, and of (l) only n = 0 is refused.

Why (d) needs `--no-grid-replacement`: without it, a finite-difference step of the 1/2⁺ energy
or width re-places the sub-point grids. The difference quotient then differs from the derivative
at fixed grids by ~1e-2 on those two parameters. This happens for any sub-point segment, whether
or not it is transformed.

An engine older than November 2026 ignores both keywords silently and compares σ itself with
the data. On the build installed before the feature (123c747), 15 checks fail: (a), the
integration half of (b), (c), (f), (g) and all of (e). The "order matters" half of (b) and (d)
pass on any build: they compare plain cross sections with each other, and the gradient with
its own finite differences. A missing or crashed run is a failure, never a pass: with no
working binary all 19 checks fail.

    ./tests/transmission_yield/check.sh path/to/AZURE2      # ~40 s
    REF_AZURE2=path/to/reference/AZURE2 ./tests/transmission_yield/check.sh path/to/AZURE2   # with (m)
