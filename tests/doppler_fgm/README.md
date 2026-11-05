# doppler_fgm

`doppler <T>` on a `<targetInt>` line: free-gas Doppler broadening of the listed segments at the effective
temperature T (kelvin), SAMMY's free-gas model,

σ_D(E) = 1/(u√π v²) ∫₀^∞ [e^{−(v−w)²/u²} − e^{−(v+w)²/u²}] w² σ(w²) dw,  u² = kT m/M,

v = √E, w = √E′ in lab-equivalent energies. It is applied to the cross section of every sub-point **first**,
then any `transmission`/`yield` transform, then the line's resolution kernel (SAMMY's order).

The models are built at run time by `make_models.py` from `tests/reaction_rate/models/n12C.azr`, with the
1/2⁺ resonance (E_cm = 0.4247 MeV) narrowed from Γ_n = 34 keV to 200 eV, so that 300 K broadening (Doppler
width ≈ 60 eV) is a large effect. `analyze.py` integrates the free-gas formula independently on a 1 eV
grid of unbroadened cross sections (mode 3).

| check | what |
|---|---|
| (a) | `doppler` alone: σ_t and σ_γ at the data points equal the independent integral (to 2e-4), and the step changes σ_t by more than 5 %. |
| (b) | `doppler` + a 100 eV Gaussian: the Gaussian average of σ_D (to 2e-4); the Gaussian alone differs by more than 1 %. |
| (c) | transmission and yield: the transform of σ_D at the point (doppler alone), or the Gaussian average of the transform of σ_D (with the Gaussian), to 2e-4; transforming the broadened σ instead is more than 1 % off. |
| (d) | `AZURE_GRAD_CHECK` on a project with all three kinds: analytic gradient = finite differences (1e-6 of the largest) and the Jacobian path (1e-10), `--no-grid-replacement`. |
| (e) | `doppler 0` and `doppler` without a temperature are refused with an ERROR and no output. |

The lines under the Gaussian use 1000 sub-points: the Gaussian's window is wider than the Doppler one and
the grid must still resolve the 200 eV resonance (400 sub-points leave 4e-4 at the σ_t dip, 1600 give 3e-5).

An engine older than November 2026 refuses the whole `<targetInt>` block ("Could not fill data object from
file"), so every check fails on it: it does not ignore the keyword.

    ./tests/doppler_fgm/check.sh path/to/AZURE2      # ~75 s
