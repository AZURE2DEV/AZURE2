# doppler_fgm

`doppler free_gas <T> K` on a target record: free-gas Doppler broadening of the data sets naming the target at the
effective temperature T (kelvin), SAMMY's free-gas model,

σ_D(E) = 1/(u√π v²) ∫₀^∞ [e^{−(v−w)²/u²} − e^{−(v+w)²/u²}] w² σ(w²) dw,  u² = kT m/M,

v = √E, w = √E′ in lab-equivalent energies. It is applied to the cross section of every sub-point **first**,
then any transmission/yield transform, then the beam's resolution kernel (SAMMY's order).

The models (format 3) are built at run time by `make_models.py` from `tests/reaction_rate/models/n12C.azr`, with the
1/2⁺ resonance (E_cm = 0.4247 MeV) narrowed from Γ_n = 34 keV to 200 eV, so that 300 K broadening (Doppler
width ≈ 60 eV) is a large effect. `analyze.py` integrates the free-gas formula independently on a 1 eV
grid of unbroadened cross sections (mode 3).

| check | what |
|---|---|
| (a) | Doppler alone: σ_t and σ_γ at the data points equal the independent integral (to 2e-4), and the step changes σ_t by more than 5 %. |
| (b) | Doppler + a 100 eV Gaussian beam spread: the Gaussian average of σ_D (to 2e-4); the Gaussian alone differs by more than 1 %. |
| (c) | transmission and yield: the transform of σ_D at the point (doppler alone), or the Gaussian average of the transform of σ_D (with the Gaussian), to 2e-4; transforming the broadened σ instead is more than 1 % off. |
| (d) | `AZURE_GRAD_CHECK` on a project with all three kinds: analytic gradient = finite differences (1e-6 of the largest) and the Jacobian path (1e-10), `--no-grid-replacement`. |
| (e) | a temperature of 0 K and none are refused with an ERROR and no output. |
| (f) | the fitted temperature, `doppler free_gas T K fit [prior T sigma] [grid_max Tmax]` (parameter `<target>.doppler`): data = the `dop` model's σ_t and σ_γ at 450 K (1 % errors, every level fixed), `free_gas 300 K fit grid_max 900`; MIGRAD (`--use-gradient`) and `--use-lm` recover T within 1 K. |
| (g) | `AZURE_GRAD_CHECK` on `grad_fit` (targets D0, D1 `fit grid_max 900` on the Doppler-only and the transmission/yield data sets, D2 `fit prior 300 20` on the Doppler + Gaussian transmission/yield ones), the three temperatures free with the widths: analytic = finite differences (1e-6 of the largest) and the Jacobian path (1e-10). The old model shared one temperature between the plain and the transmission lines; a cross section takes a target of kind none or layer and a transmission one of kind sample, so they are two targets now. |
| (h) | `fit prior 300 20` against the same values with 20 % errors: T lands between 300 K and the default layout temperature 363 K, and `Total-Doppler-Chi-Squared` = ((T − 300)/20)². |
| (i) | malformed fit forms are refused with an ERROR and no output: a prior without its width, a negative width, `grid_max` ≤ T, T < 0. The old line's `vary 2` has no format-3 form: the converter refuses it (`convert_constructs`). |
| (j) | only with `REF_AZURE2=<a reference build that reads format 3>`: the plain Doppler models' `AZUREOut` files and `chiSquared.out` are byte-identical to that build's. |

A fitted line lays its sub-point grids and every Doppler integral's window out once, for the layout temperature
(Tmax, else max(T + 3σ, 1.21 T), or 2 T without a prior), and the Gauss-Legendre pieces for the line's T; only
the kernel values move with the parameter, so the analytic derivative is exact. The fits recover 449.95 K rather
than 450.00 K: the data were computed with the fixed temperature's window and pieces, the fit with those of the fit
form, and the two quadratures differ by about 1e-6.

The data sets under the Gaussian use 1000 sub-points: the Gaussian's window is wider than the Doppler one and
the grid must still resolve the 200 eV resonance (400 sub-points leave 4e-4 at the σ_t dip, 1600 give 3e-5).

The old generator (old `<targetInt>` lines) is kept in `tests/convert/fixtures/tests/doppler_fgm/make_models.py` for
the pyazr tests of the old file class (`tests/pyazr/old_tree.py`).

    [REF_AZURE2=path/to/reference/AZURE2] ./tests/doppler_fgm/check.sh path/to/AZURE2      # ~1.5 min
