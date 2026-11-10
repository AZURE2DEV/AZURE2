# Model data (`<modelData>`)

`check.sh` runs the engine on T2 at the generator's parameters (`tests/upgrade/T2/T2_true.azr`, data in
`tests/upgrade/T2/data`) with a `<modelData>` section appended, and compares the model value of each datum with a
number computed by hand from the widths written in `<levels>` (Brune parameters: the observed widths at the observed
energy, `tests/upgrade/T2/README.md` section 1).

Level 5: J^π = 3⁻, E_x = 13.2692 MeV; Γ_p = 1.702 keV (p + ¹⁵N, s = 1, l = 2; l = 4 has 0), Γ_α0 = 12.92 keV,
Γ_α1 = 0. Pair 1, p + ¹⁵N: spins ½ and ½, S = 12.1274 MeV, masses 1.007825 and 15.000109 u.

| datum | quantity | hand value |
|---|---|---|
| `wg_L5` | resonance strength | ω = 7/(2·2) = 7/4; ωγ = 1.75 × 1.702 × 12.92 / 14.622 = 2.6318028 keV |
| `Ex_L5` | E_x | 13.2692 MeV |
| `Ecm_L5` | E_c.m. of p + ¹⁵N | 13.2692 − 12.1274 = 1.1418 MeV |
| `Elab_L5` | E_lab of p + ¹⁵N | 1.1418 × 16.007934 / 15.000109 = 1.2185150815 MeV (non-relativistic, as the engine) |
| `G_L5` | total width | 1.702 + 12.92 + 0 = 14.622 keV |
| `Gp_L5` | partial width, p + ¹⁵N | 1702 eV |
| `BR_L5` | branching ratio Γ_α0/Γ | 12.92 / 14.622 = 88.360005 % |
| `a_p15N` | scattering length, p + ¹⁵N, s = 1 | the 1⁻ value of `scatteringLengths.out` (the same function) |
| `C_L1` | ANC of level 1 moved to E_x = 12.0 MeV (p + ¹⁵N closed) | 496 fm^−½, the value written |

Each value of the section is the hand value plus one error, so every residual is −1 and
`Total-Model-Data-Chi-Squared` is 7 plus the scattering length's (a₀ / 100 fm)². In the formal and Park bases the
values must be the same (observed quantities do not depend on the basis); the formal basis reaches them through the
formal → Brune iteration, which stops at 10⁻⁶ MeV in the energies, hence the relative tolerance 2 × 10⁻⁶ there.

The fit (segment 5 alone, the energy and α₀ width of level 5 free): a strength 3σ above the hand value (σ = 0.02 keV)
must pull Γ_α0 up, since ∂(ωγ)/∂Γ_α0 = ω Γ_p²/Γ² > 0 with Γ_p fixed. The refusals each change one thing.
