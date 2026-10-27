# capture_ay_coupling

The capture analyzing power couples the channel spin light particle first, like the rest of AZURE2.

`CNuc::CalcCaptureAnalyzingPower` builds Seyler and Weller's `a_k` and `b_k` (Phys. Rev. C 20 (1979) 453,
Eqs. 20 and 21).  They couple the channel spin target first, `<j2 m2 j1 m1 | s nu>`; AZURE2 couples it
light particle first, `<j1 m1 j2 m2 | s nu>` (Lane and Thomas: `PolarizationFunc`, the M1 external-capture
spin operator, the secondary-gamma coupling).  The two differ by `(-1)^(j1 + j2 - s)` per channel-spin
amplitude, which on the product `R R'^*` of a `b_k` term is `(-1)^(s - s')`.  That factor was missing until
2026-10-26 (audit item 9): for a target with spin, every channel-spin off-diagonal term of the capture `A_y`
had the wrong sign, so fitting it was the same as flipping the sign of every width with odd `j1 + j2 - s`,
inconsistent with the particle `A_y` of the same widths.  A spin-0 target has one channel spin and was never
affected (`tests/13N_capture_ay`).

The same commit refuses an analyzing-power segment whose projectile is not spin-1/2: the particle route
returned 0 for any other spin, the capture route ran and returned `3 <S_y> / (x + 1)`, which is the Madison
`A_y` only for `x = 1/2` (1.5 times Ohlsen's Cartesian `A_y` for a deuteron).

## Models (`make_models.py`)

`twolevel`: n + a 1/2- target (channel spins 0 and 1) into a fictitious compound nucleus with a 0+ ground
state; two 1+ levels decaying by M1 to it, level A (E_cm 1.0 MeV) with a neutron width in channel spin 0
only, level B (1.3 MeV) in channel spin 1 only, both l = 1, gamma widths of 1 meV, no external capture.
The levels share no particle channel, so the level matrix is diagonal and the two capture amplitudes are

    T_s(E) = gamma_s gamma_gamma / D_lambda(E),
    D = E_lambda - E - gamma_s^2 [S(E) - S(E_lambda)] - i gamma_s^2 P(E)      (Brune parametrization)

with the closed-form neutral l = 1 functions S = -1/(1 + rho^2), P = rho^3/(1 + rho^2).  Every other
factor is common to the two amplitudes and drops out of `A_y`; between the levels the two have different
phases, so `A_y` reaches 0.19 and its sign is the coupling order.  `deuteron`: the same with a spin-1
projectile.

## Checks (`check.sh path/to/AZURE2`)

- (a) `verify.py`: the engine's `A_y` (mode 3, `--use-brune`) against a magnetic-substate calculation built
  from the engine's own reduced-width amplitudes (`parameters.out`): channel-spin amplitudes in the
  convention under which the unpolarized Seyler-Weller coefficients reproduce AZURE2's Blatt-Biedenharn
  cross section (02_angular_distributions.md of the 2026-10-25 audit), the projectile spin taken out of the
  channel spin light particle first, `A_y = 2 Im sum(A_up A_down^*) / sum(|A_up|^2 + |A_down|^2)`.
  65 points, tolerance 1e-5, and `A_y = 0` at 0 and 180 degrees.
- (b) the deuteron segment is refused with `ERROR: test segment #1: an analyzing power needs a spin-1/2
  projectile`.
- (c) `tests/13N_capture_ay` (spin-0 target) still gives its expected chi2.

dev 7d10253 (before the fix): (a) the engine's `A_y` is the exact negative of the truth (worst 0.37 against
it, 6.6e-7 against its negative); (b) the deuteron segment ran.  Fixed build: (a) 6.6e-7, (b) refused, (c) unchanged.
