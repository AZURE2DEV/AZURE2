# stopping_yield — the yield of a layer with SRIM-2013's stopping

A target whose record says `stopping evaluation` takes its stopping from the engine's `stopping/` library: SRIM-2013's
stopping powers, ported from LunaSRIM (docs/upgrade/stopping.md; **cite SRIM**: J. F. Ziegler, M. D. Ziegler and
J. P. Biersack, Nucl. Instrum. Methods Phys. Res. B 268 (2010) 1818). This test compares the engine's yields with numbers
that do not come from the engine: LunaSRIM's stopping and energy loss, computed once by `make_project.py` and kept in
`models/truth.json`.

## The convention (SetupSpecs.h `StoppingSpec`, TargetStopping.h)

The engine's stopping is per **active nucleus**, lab, MeV cm²:

```
(E1)  eps_active(E) = 1e-21 eps_atom(E) / f_active,      f_active = n_active / sum_i n_i
```

with eps_atom the material's stopping cross section per atom (eV/(10¹⁵ atoms/cm²), stopping/'s unit), n_i the atoms per
formula unit of the target's `composition`, n_active those of the active nucleus's element. Equivalently: the stopping
per formula unit divided by the active atoms per formula unit. The isotopic part of the active fraction is the
target's `isotope` line, which multiplies the yield (data_kinds.md section 4.6). `stopping_convention_test.cpp` (ctest
`stopping_convention`) checks (E1) on Ta₂O₅ with ¹⁸O active (divisor 5) and TiN with ¹⁴N active (divisor 1) against
2 eps_Ta + 5 eps_O and eps_Ti + eps_N put together from the library's per-atom values of each element in the solid, to
1e-14, and TiN at 0.5 MeV against the 12.00133921 eV/(10¹⁵ atoms/cm²) that `azure2-stopping` prints, to 1e-9.

## The projects

The narrow resonance of tests/layer_yield (p + ¹⁵N → ¹⁶O(1⁻) → α + ¹²C, E_r = 0.9 MeV c.m., Γ_p = 60 eV, Γ_α = 40 eV)
on layers of TiN (`composition Ti 1 N 1`, active ¹⁵N, solid, no straggling, no beam spread). The projectile's mass is
the pair's in `<levels>` (1.007825 u) in the engine and in LunaSRIM alike. f = m_N/(m_p + m_N) is E_cm/E_lab.

**(S1) Plateau of a thick layer** (`plateau`, target TP: N = 4 × 10¹⁸ ¹⁵N/cm², Δ ≈ 60 keV; E_0 = E_r/f + 2, 4, 8,
16 keV). The Breit–Wigner integral over the layer with k_r, the widths and eps constant across the resonance:

```
(S1)  Y(E_0) = 1e-24 (2 pi/k_r^2) omega gamma / (f eps_r) [atan(2(f E_0 - E_r)/Gamma) - atan(2(f E_out - E_r)/Gamma)]
```

whose plateau is 1e-24 (λ_r²/2) ωγ (m_p + m_N)/(m_N eps_r) (λ_r c.m., eps_r lab per active nucleus), eps_r =
eps_active(E_r/f) from LunaSRIM (1.634389 × 10⁻²⁰ MeV cm²), E_out from `lunasrim.straggling.energy_loss`. Tolerance
3 × 10⁻⁴; measured 8.3 × 10⁻⁵. What limits it: the R-matrix line shape is not the Breit–Wigner form (k², the proton
penetrability and the shift vary across the tails; tests/layer_yield measured the same 1.3 × 10⁻⁴ with a constant
stopping), not the stopping: the sub-point quadrature over the resonance is good to about 1.5 × 10⁻⁴ (60 points) and
was not tuned.

**(S2) Width of the step of a layer** (`step`, target TS: N = 5 × 10¹⁷, Δ ≈ 8.15 keV). The yield rises where E_0
reaches the resonance and falls where E_out does; the distance between the two half-height points is E_in − E_out of
the layer at the resonance. Both edges are scanned every 10 eV (lab) and five points in the middle. The engine's scan and
the Breit–Wigner scan (S1) with LunaSRIM's E_out at every energy are measured by the same rule (half the mean of the
middle points; linear interpolation), because the Lorentzian tails lower the middle of the step and move both points
by about 0.3 eV. Measured: 8.152753 keV against 8.152754 keV, 2.4 × 10⁻⁷; tolerance 2 × 10⁻⁵. (LunaSRIM's E_in − E_out at
E_r is 8.151921 keV; its RK4 integrator agrees with stopping/'s quadrature to 5 × 10⁻⁶, stopping.md section 9.1.)

**(S3) Evaluation against a table of the library's own values** (`plateau_table`, target TT: the layer of TP with
`stopping table "data/tin_p.tab" MeV eV/1e15atoms/cm2`). The table is the output of `azure2-stopping stopping
--ion-z 1 --ion-mass 1.007825 --material TiN --from 0.5 --to 1.5 --points 5000` (ten digits; its command is in its
header); `check.sh` checks that the committed file equals the tool's present output. Measured 2.1 × 10⁻¹⁰,
tolerance 10⁻⁸: the table's ten printed digits and its interpolation, linear in ln E – ln eps between points 2.2 × 10⁻⁴
apart in ln E; the two targets have the same sub-point grids, so the quadrature cancels.

**(S4) The fit** (`models/fit.azr`). Two TiN layers (F1, F2: N_true = 6 × 10¹⁷ and 9 × 10¹⁷, started at 5 × 10¹⁷ and
10¹⁸), each scanned across its whole step every 0.3 keV; the data are (S1) with eps → s eps, s_true = 1.03, and E_out
from LunaSRIM for the areal density s N (∫dE/(s eps) = N), times 1 + 1 % Gaussian noise (seed 20261111). The two targets
share one stopping scale, `stopping_scale.TiN.1H` (`stopping evaluation scale 1 fit prior 1 4%`). The plateau fixes s
(Y ∝ 1/s), each width s N_i. Measured: s = 1.02927 ± 0.00111 (pull −0.66), N1 = 6.0039 × 10¹⁷ ± 6.4 × 10¹⁴ (+0.61),
N2 = 9.0044 × 10¹⁷ ± 9.6 × 10¹⁴ (+0.46), χ²/ν = 71.5/90; the objective the fit minimised minus the data's χ² is the
prior ((s − 1)/0.04)² = 0.5354. Required: one scale parameter, every pull within 2, χ²/ν within 3√(2/ν) of 1, the prior
in the objective to 10⁻³. The fit runs at one thread (a fit at several threads is not reproducible to the last digit).

## Refusals (check.sh, section 5)

`E_FIELD_REQUIRED` (an evaluation without `composition`), `E_STOPPING_ACTIVE` (the active element not in the
composition), `E_STOPPING_PROJECTILE` (a ¹²C(n,n) data set on an evaluated target), `E_STOPPING_SCALE_CONFLICT` (two
TiN targets, one with `scale 1.02 fixed`), `E_STOPPING_RANGE` (a table from 0.93 MeV, above the back of the thick layer),
`E_STOPPING_DATA` (`AZURE2_SRIM_DATA_DIR` an empty directory: refused by `--validate` and by the run).

## Files

- `make_project.py` writes `models/` with LunaSRIM (any Python that imports `lunasrim` 1.0.0; no path of a machine), but
  `models/data/tin_p.tab`, the library's table (the command in its header).
- `check.py TRUTH OUTPUT` (S1)–(S3) against a menu-1 output; `check.sh AZURE2` the whole test (it takes
  `azure2-stopping` from `<build>/stopping/`).
- `stopping_convention_test.cpp` the convention (E1).
