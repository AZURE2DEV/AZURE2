# layer_yield — the yield of a layer of finite thickness

`quantity yield` on a target of `kind layer` (docs/upgrade/energy_grid.md section 5, data_kinds.md section 4.6):
reactions per incident particle (or per µC) of a beam that crosses a layer of N active nuclei per cm². The engine's
kernel is the target integration of `src/TargetFold.cpp` on the point's sub-points, not divided by N, with the back of
the layer found on the same nodes (`include/LayerExit.h`).

## Formulas (frames and units written out)

Symbols: E the beam energy in the **lab** (MeV); E_0 at the surface, E_out at the back of the layer; E_cm = f E with
f = m_A/(m_a + m_A) (atomic masses in u, non-relativistic, as `Segment::CMEnergyFactor`); σ(E_cm) the cross section
in b, c.m.; ε(E) the stopping cross section per **active** nucleus, **lab**, MeV cm² (the target's `stopping
formula`, x = E/MeV); N the areal density of active nuclei, atoms/cm²; 10⁻²⁴ cm² per b.

```
(Y1)  Y(E_0) = 1e-24 ∫_{E_out}^{E_0} σ(f E) / ε(E) dE                   reactions per incident particle
(Y2)  ∫_{E_out}^{E_0} dE / ε(E) = N                                      defines E_out
(Y3)  Y_uC = Y × 1e-6 C / (q e),  e = 1.602176634e-19 C                  per µC of a beam of charge state q
```

The engine works in the c.m. frame: its sub-points carry E_cm and ε_cm = f ε(E_cm/f); since dE_cm/ε_cm = dE/ε, (Y1)
and (Y2) are the same integrals. The quadrature is `GridQuadrature`'s: two-point Gauss–Legendre on each interval
between sub-points, σ and ε interpolated linearly between the nodes; E_out solves (Y2) **with that rule on those
nodes** (bisection to the last bit inside the interval that holds it). With a beam spread the kernel averages σ/ε over
a Gaussian about each depth energy (today's thin-target kernel), which is ⟨σ⟩/ε to the curvature of 1/ε across it.

The active fraction. ε and N must count the same nuclei. With the user's formula (the only stopping of this
branch) ε is already per active nucleus: its stoichiometric part is in the formula (`experimental_effects.rst`,
"compound and mixed targets"). An `isotope` line on the target multiplies the yield by its abundance a (the
`Weighted` node of the model tree), which is the isotopic part of the active fraction: ε and N then count atoms of
the element, and ε_eff = ε/a, N_active = a N. The stopping module will divide its compound ε by the stoichiometric
fraction in the place named in `LayerExit.h` (with the stopping scale).

## Limits, each computed by `check.py` from `models/truth.json` (never from the engine)

**(C) Constant σ.** From (Y2): Y = 1e-24 σ_0 N exactly. On the engine's nodes the identity holds to rounding because
E_out comes from the same rule (`layer_kernel_test.cpp`: 1e-15). In the project the cross section is the s-wave hard
sphere of n + ¹²C (one J^π = 1/2⁺ channel, radius a = 0.1 fm, g = (2J+1)/((2s_n+1)(2s_C+1)) = 1):

```
(C1)  σ(E_cm) = (π/k²) g |1 − e^{−2ika}|² = 4π sin²(ka)/k²,   k² = 2 μ E_cm/(ħc)²,  μ = m_n m_C/(m_n + m_C) × 931.494 MeV
```

which varies as dlnσ/dlnE = −(ka)²/3 < 6e-6 across the layer. For σ linear in E, (Y1) is 1e-24 N σ(Ē) with the
weighted mean energy Ē = ∫E dE/ε / N. For ε = a0 x^a1 (b = 1 − a1, c = 2 − a1):

```
(C2)  E_out = (E_0^b − a0 b N)^(1/b),     Ē = (b/c)(E_0^c − E_out^c)/(E_0^b − E_out^b)
```

The tolerance is the slope of σ times the quadrature's offset of Ē (the weights 1/ε_lin are off by about (h/E)²/4,
h the node spacing) plus the output's 11 digits. Data sets `const` (ε = 1.2e-24/x, N = 1e20), `const_iso` (the same
target with `isotope 12C 0.5`: Y/2), `const_spread` (a 0.3 keV beam spread; 100 quadrature points: the Gaussian on
the depth grid limits it, measured 2.8e-8).

**(T) Thin layer**, Δ = E_0 − E_out much smaller than the scale of σ. Expanding σ about E_m = E_0 − Δ/2:

```
(T1)  Y = 1e-24 N [σ(E_m) + σ'(E_m)(Ē − E_m) + σ''(E_m) Δ²/24 + O(Δ³)]
```

`thin` (¹⁵N(p,α), ε = 4.5e-20 x^-0.6, N = 5e15: Δ ≈ 0.3 keV, E_0 = 0.60–0.80 MeV below the resonance) against
1e-24 N σ(E_m) with the remainder of (T1) estimated from the engine's bare σ at E_m and E_m ± 5 keV (`thin_bare`);
tolerance 1e-7 (the bare σ's 11 digits and the third order). And the existing thin-target quantity (data_kinds.md
eq. 2): `thin_sigma` (`quantity sigma` on the same target) is Y/(1e-24 N) to 1e-6 (two grids, the yield's 2 % longer,
and two E_out: the yield's from the quadrature, the thin target's from `TargetEffect::TargetThickness`).

**(K) Thick layer**, N large. (Y1) is additive in its lower limit, so with the infinite target (E_out = 1 keV) of the
same ε:

```
(K1)  Y_layer(E_0) = Y_∞(E_0) − Y_∞(E_out)
```

`thick` (N = 1.1e19: E_out = 0.375, 0.461, 0.540 MeV for E_0 = 1.00, 1.05, 1.10 MeV, the resonance inside) against `thick_inf` at E_0
and at E_out (C2); tolerance 1e-5: two adaptive grids across the resonance (3.2e-6; 1.5e-6 with the grids four times
finer). Y_layer/Y_∞(E_0) − 1 = −Y_∞(E_out)/Y_∞(E_0) = −4.3e-6 is printed.

**(R) Isolated narrow resonance.** One 1⁻ level of ¹⁶O at E_r = 0.9 MeV c.m. (E_x = 13.0274 MeV), proton s wave
Γ_p = 60 eV, α₀ p wave Γ_α = 40 eV (observed widths at E_r), no other level. The Breit–Wigner, in the c.m. frame:

```
(R1)  σ(E_cm) = (π/k²) ω Γ_p Γ_α / ((E_cm − E_r)² + Γ²/4),   Γ = Γ_p + Γ_α,   ω = (2J+1)/((2j_p+1)(2j_N+1)) = 3/4
```

with k, Γ_p, Γ_α constant across the scan (k = k_r at E_r) and ε constant (lab), (Y1) with E_cm = f E, dE = dE_cm/f:

```
(R2)  Y(E_0) = (1e-24 / (f ε)) (π/k_r²) ω Γ_p Γ_α ∫_{f E_out}^{f E_0} dE_cm / ((E_cm − E_r)² + Γ²/4)
(R3)         = 1e-24 (2π/k_r²) ωγ / (f ε) [ arctan(2(f E_0 − E_r)/Γ) − arctan(2(f (E_0 − Δ) − E_r)/Γ) ],   ωγ = ω Γ_p Γ_α / Γ
(R4)  plateau (Δ ≫ Γ): 1e-24 (λ_r²/2) ωγ (m_p + m_N)/(m_N ε),   λ_r = 2π/k_r (c.m.), ε lab per active nucleus
```

(R4) is (λ_r²/2) ωγ/ε_cm with ε_cm = f ε (data_kinds.md eq. 7 with its frames). Units: π/k² in fm² is 1e-2 b; widths
in MeV. `res` (ε = 5e-21 constant, N = 1e18: Δ = 5 keV lab; E_0 every 0.2 keV from E_r/f − 1.5 keV to + 6.5 keV)
against (R3); difference over the plateau, tolerance 3e-4. What limits it: the R-matrix line shape is not (R1) — k²,
the proton penetrability and the shift vary across the tails; the engine's bare σ (`res_bare`) differs from (R1) by
2.4e-4 within one width and 9.5e-3 at 30 widths (where σ is 1/3600 of the peak). Measured: 1.29e-4 of the plateau
(1.17e-4 with grids four times finer: the quadrature is 1e-5).

`deriv`: the same layer with N a fit parameter (`fit prior 1e18 10%`: the grid laid out for N0 + 3σ = 1.3 N0 and cut
at E_out), E_0 such that E_out lies on the resonance, against (R3) at N0: 1.5e-4 of the plateau, tolerance 3e-4.

**(D) Derivative with respect to N.** From (Y2), ∂E_out/∂N = −ε(E_out), so

```
(D1)  ∂Y/∂N = 1e-24 σ(f E_out)              (b per atoms/cm²: the cross section at the back of the layer)
```

The engine differentiates N by central differences (the other apparatus parameters' route, step 1e-6 (N + 1)).
`pyazr_layer_test.py` compares its Jacobian column with 1e-24 × the engine's bare σ at E_out: 1e-4 within ten widths
of the resonance, 3e-3 beyond, both the kernel's linear interpolation of σ between sub-points (its yield's exact
derivative is 1e-24 σ_lin(E_out)).

**(Q) Per unit charge.** `res_q1`, `res_q2`: the same layer in 1/uC with charge states 1 and 2 against (Y3) from `res`,
and their ratio 2, to the output's 11 digits.

## The fit

`models/fit.azr`: a scan of the resonance on a layer of N = 1.2e18 (Δ = 6 keV), 41 energies, data the arctan form
(R3) times 1 + 2 % Gaussian noise (seed 20261111), errors 2 %. Free: N (from 1e18) and Γ_p (from 45 eV). `check.sh`
requires both within two standard errors of the truth and χ²/ν within 3√(2/ν) of 1, and the same outputs at 1 and 2
threads (Minuit's call count left out: docs/upgrade/fit_threads.md section 3).

## Files

- `make_project.py` writes `models/layer_yield.azr`, `models/fit.azr`, `models/data/`, `models/truth.json`.
- `check.py TRUTH OUTPUT` the closed forms above against a menu-1 output; `check.sh AZURE2` runs it, the fit, the
  refusals and the pyazr test.
- `layer_kernel_test.cpp` (ctest `layer_kernel`) the kernel on hand-set sub-points: (C) to 1e-14, the weights' sum,
  E_out of a constant ε, a linear σ, NaN when the grid ends inside the layer.
- `pyazr_layer_test.py` (ctest `pyazr_layer_yield`) records built in Python, (R3), (D1), the CLI's values.
