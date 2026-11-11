# transmission_yield

Quantities `transmission` and `neutron_yield` on a target of kind `sample` (its `areal_density` n in atoms/barn): the transmission
T = exp(-n σ) of a thick sample, and SAMMY's self-shielded capture or reaction yield
Y0 = (σ/σ_t)(1 − exp(−n σ_t)), σ_t the total cross section of the neutral entrance pair from
the optical theorem (a hidden `sigma_total` term of the data set). Under a sub-point target
effect the transform is applied at every sub-point **before** the kernel (SAMMY manual III.E.1):
a detector averages T(E), and ⟨e^{−nσ}⟩ ≠ e^{−n⟨σ⟩}.

The models (format 3) are built at run time by `make_models.py`:
- from `tests/reaction_rate/models/n12C.azr`: n + ¹²C elastic and capture, with an s-wave 1/2⁺
  resonance at E_cm = 0.425 MeV and Γ_n = 34 keV. With n = 0.15 atoms/b, nσ_t ≈ 1 at the peak.
- from `tests/natural_neon/natural_neon.azr`: three compound nuclei, with natural sums.

`analyze.py` does the comparisons.

| check | what |
|---|---|
| (a) | With no spread, T and Y0 match the formula applied to the plain σ_t and σ_γ segments, to 1e-9. |
| (b) | With a 15 keV (lab) Gaussian spread, T and Y0 match an independent integration, to 2e-4. The integration averages the transform of a 0.1 keV grid of unbroadened cross sections over the kernel (±5σ, normalised over the window). Transforming the *broadened* cross sections is off by 6 % (T) and 5 % (Y0) at the resonance: the check discriminates. |
| (c) | For a natural sum over three nuclei, T = exp(−n σ_nat), to 1e-9. The sample decomposition file notes that its total column is T. |
| (d) | The gradient check (`AZURE_GRAD_CHECK`) passes on a project with T and Y0 segments, with and without the spread, and every width free. The analytic gradient matches finite differences (to 1e-7 of the largest) and the Jacobian path (to 1e-10). Uses `--no-grid-replacement`, explained below. |
| (e) | Each of these is refused with an ERROR and produces no output: a yield of the elastic, a charged entrance pair, a natural sum under a sub-point spread, and a yield on a natural sample. (A transmission and a yield on one line, and a transmission of a differential cross section, have no format-3 form: the converter refuses the old lines, `convert_constructs`.) |
| (f) | On extrapolations (mode 3), the T and Y0 curves match the formula applied to the plain curves, to 1e-9. |
| (g) | A yield's σ_t is the total, not the elastic. The model is n + ²⁰Ne (natural_neon) with a free 1 keV width on every open reaction channel. A yield to α + ¹⁷O equals the formula with the `sigma_total` σ_t to 1e-9, at energies where the elastic differs from the total by more than 1e-3. |

The fitted sample thickness: `areal_density n atoms/b fit [prior n sigma]` on the sample makes n the fit parameter
`<target>.areal_density` (sigma the Gaussian prior width in atoms/b); data sets naming one sample share its
thickness (the old `thicknessgroup`).

| check | what |
|---|---|
| (h) | The data are the model's T (and Y0) at n = 0.18 atoms/b with 1 % errors, and every level is fixed. Starting from `0.15 atoms/b fit`, MIGRAD (`--use-gradient`) and `--use-lm` recover `Sd.areal_density` (`Sg.areal_density`) = 0.18 to 1e-4. |
| (i) | One sample for the transmission and the yield: one parameter `S.areal_density` is fitted from both; the second data set's row `S.areal_density[g]` in param.sav is written as its value. |
| (j) | The gradient check runs with seven thicknesses free together with every width. The lines cover: no target effect, a Gaussian spread with a 0.01 atoms/b prior, Doppler alone (one sample for a transmission and a yield), and Doppler plus the Gaussian. The analytic gradient matches finite differences to 1e-6 of the largest, and the Jacobian path to 1e-10. |
| (k) | A prior: `0.15 atoms/b fit prior 0.15 0.01` on 20 % data holds n between 0.15 and 0.18. chiSquared.out's `Total-Thickness-Chi-Squared` equals ((n − 0.15)/0.01)². |
| (l) | Each of these is refused with an ERROR naming the field: a negative or missing prior width, n = 0, and a cross section on the sample of a transmission. (The old `vary 2`, a `thicknessgroup` without its label and a group member without its own block are the converter's refusals, `convert_constructs`.) |
| (m) | With `REF_AZURE2` set to a reference build that reads format 3, the plain transmission / yield outputs are byte-identical to it: AZUREOut and chiSquared.out of thin, conv, the extrapolations and the natural sum. |

Why (d) needs `--no-grid-replacement`: without it, a finite-difference step of the 1/2⁺ energy
or width re-places the sub-point grids. The difference quotient then differs from the derivative
at fixed grids by ~1e-2 on those two parameters. This happens for any sub-point segment, whether
or not it is transformed.

A missing or crashed run is a failure, never a pass. The old generator (old segment lines) is kept in
`tests/convert/fixtures/tests/transmission_yield/make_models.py` for the pyazr tests of the old file class
(`tests/pyazr/old_tree.py`).

    ./tests/transmission_yield/check.sh path/to/AZURE2      # ~40 s
    REF_AZURE2=path/to/reference/AZURE2 ./tests/transmission_yield/check.sh path/to/AZURE2   # with (m)
