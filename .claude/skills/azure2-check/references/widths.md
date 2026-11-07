# Wigner and Weisskopf limits in AZURE2

The dimensionless reduced width θ²_λc compares a level's width in channel c with the single-particle (Wigner) limit
of that channel at the channel radius a_c (fm). AZURE2 states the limit in two forms that are not interchangeable.
This file says which is which, how pyazr computes them, and how γ-ray strengths are judged.

## The two forms of the Wigner limit

1. **Width form (GUI, `ChannelWidth.wigner_width`), in eV.** Γ_W = 2 P_l(E_r) γ²_W, with P_l the penetrability of
   orbital angular momentum l at the resonance energy E_r (c.m.). It depends on energy. θ² = Γ/Γ_W, with Γ the
   observed partial width. Prefer this form: it can be checked one channel at a time in the GUI.
2. **Amplitude form (API, `Parameter.wigner_limit`, `ChannelWidth.wigner_gamma2`), in MeV.**
   γ²_W = 3ħ²/(2μ a_c²), with μ the reduced mass of the channel; already squared and independent of energy.
   θ²_formal = γ²_formal/γ²_W, with γ_formal the formal reduced-width amplitude (`g_int` in `parameters.out`,
   squared).

Never divide by `wigner_limit` twice. The two forms differ by the level-shift factor
J = 1/(1 + Σ_c γ_c² dS_c/dE): negligible for a narrow level, a factor 70 or more for a background pole at
E_x ≈ 30 MeV. AZURE2 reports observed partial widths, Γ_c = 2 P_c γ_c² J.

## Conventions of the unit

AZURE2's γ²_W is the Teichmann-Wigner sum-rule unit 3ħ²/(2μa²) (`pyazr.widths.teichmann_wigner()` returns it
unchanged). References that use ħ²/(μa²) (`lane_thomas_wigner()`, 2/3 of it) quote a θ² 1.5 times larger. State the
unit next to every θ² you report. Builds before 2026-10-27 (expires: 2027-02-01) had docstrings that called the
limit ħ²/(μa²); the value was always 3ħ²/(2μa²).

## How pyazr obtains Γ_W

Γ_W needs P_l(E_r), which the API does not expose. `wigner_widths` obtains it from AZURE2 itself: it transforms a
probe vector with every reduced-width amplitude set to a small ε, so the shift denominator is 1 and
Γ_c(ε) = 2 P_c ε². This costs one extra `TRANSFORM_RWA` at AZURE2's own radii and Coulomb functions, and reproduces
the GUI's Γ_W (³He+α, 7/2⁻ level at E_x = 4.57 MeV: Γ_W = 3.487 × 10⁵ eV, θ² = 0.494).

## γ-ray strengths

Strengths in Weisskopf units (W.u.) against the Weisskopf estimates, with E_γ in MeV and Γ in eV:

| Multipole | Γ_W.u. (eV) |
|---|---|
| E1 | 6.8 × 10⁻² A^{2/3} E_γ³ |
| E2 | 4.9 × 10⁻⁸ A^{4/3} E_γ⁵ |
| M1 | 2.1 × 10⁻² E_γ³ |

E1 to E4 and M1 to M4 are in `pyazr/widths.py`. A γ width is not proportional to γ² for a capture channel: the
external-capture term is linear in the amplitude. Convert with `transform_rwa`, never by squaring.
