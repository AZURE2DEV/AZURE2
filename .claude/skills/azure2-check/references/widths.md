# Wigner and Weisskopf limits in AZURE2

The dimensionless reduced width θ²_λc compares a level's width in channel c with the single-particle (Wigner) limit
of that channel at the channel radius a_c (fm). AZURE2 states the limit in two forms that are not interchangeable.
This file says which is which, how pyazr computes them, and how γ-ray strengths are judged.

## The two forms of the Wigner limit

1. **Width form (`ChannelWidth.wigner_width`), in eV.** Γ_W = 2 P_l(E_r) γ²_W, with P_l the penetrability of
   orbital angular momentum l at the resonance energy E_r (c.m.). It depends on energy. θ² = Γ/Γ_W, with Γ the
   observed partial width. Prefer this form: it can be checked one channel at a time.
2. **Amplitude form (API, `Parameter.wigner_limit`, `ChannelWidth.wigner_gamma2`), in MeV.**
   γ²_W = 3ħ²/(2μ a_c²), with μ the reduced mass of the channel; already squared and independent of energy.
   θ²_formal = γ²_formal/γ²_W, with γ_formal the formal reduced-width amplitude (`g_int` in `parameters.out`,
   squared).

Never divide by `wigner_limit` twice. The two forms differ by the level-shift factor
J = 1/(1 + Σ_c γ_c² dS_c/dE): negligible for a narrow level, a factor 70 or more for a background pole at
E_x ≈ 30 MeV. AZURE2 reports observed partial widths, Γ_c = 2 P_c γ_c² J.

## Conventions of the unit

AZURE2's γ²_W is the Teichmann-Wigner sum-rule unit 3ħ²/(2μa²): `pyazr.physics.wigner_limit(m1_u, m2_u, radius)`
gives it, and `pyazr.physics.single_particle_width(l, z1, z2, m1_u, m2_u, radius, e_cm)` the width form
(`docs/source/user_guide/physics_functions.rst`). References that use ħ²/(μa²) (2/3 of it) quote a θ² 1.5 times
larger. State the unit next to every θ² you report.

## How pyazr obtains Γ_W

`Session.wigner_widths(x)` is Γ_W = 2 P_l(E_r) γ²_W in eV for every free particle channel that is open at its
level: P_l the penetrability the engine formed the observed width with (the session's own channel functions, as
`Session.channel_functions(pair, e_cm, l)["P"]` gives them: AZURE2's radii and Coulomb functions, the hybrid model's
solution for a pair with a nuclear potential) and E_r the pair's c.m. energy at the level's physical energy (the
level's row of `Session.transform(x)`). The number is the same in Brune's, the formal and Park's parametrization.
A closed, a photon or a fixed channel has no entry there; `dimensionless_widths` gives `wigner_width` and `theta2`
of fixed open channels too, by the same definition, and its `level_energy`, `is_open` and `e_gamma` are taken at the
physical energy in every basis (`widths_check.py` prints them). Every row's `value` (the Γ of θ²), fixed rows
included, is the row of `Session.transform(x)` at the vector asked for: a fixed amplitude's observed width moves
with the level's other parameters. Example (T2 of `tests/upgrade`, α + ¹²C, l = 1, the
1⁻ level at E_x = 12.4258 MeV): Γ_W = 4.020 × 10⁶ eV.

A θ² of a level with fixed widths, or of a formal-basis session, computed with an older build is not this number:
`older-results.md`.

## γ-ray strengths

Strengths in Weisskopf units (W.u.) against the Weisskopf estimates, with E_γ in MeV and Γ in eV:

| Multipole | Γ_W.u. (eV) |
|---|---|
| E1 | 6.75 × 10⁻² A^{2/3} E_γ³ |
| E2 | 4.79 × 10⁻⁸ A^{4/3} E_γ⁵ |
| M1 | 2.07 × 10⁻² E_γ³ |

(The closed form for R = 1.2 A^{1/3} fm, Bohr and Mottelson vol. 1, eqs. 3C-38 and 3C-39; the textbooks round the
coefficients to 6.8 × 10⁻², 4.9 × 10⁻⁸, 2.1 × 10⁻².)

`Session.dimensionless_widths(x)` gives E1 to E4 and M1 to M4 (`ChannelWidth.weisskopf`, `.wu`; A is the compound
nucleus's mass number); one estimate is `pyazr.physics.weisskopf_width(e_gamma, multipolarity, mass_number,
radiation)`, in eV (the closed form for R = 1.2 A^{1/3} fm). A γ width is not proportional to γ² for a capture channel: the external-capture term is
linear in the amplitude. Convert with `Session.transform(x)`, never by squaring.

## The level norm J: saturated levels

In the Brune basis θ²/θ²_formal of an open channel is the level's J = 1/(1 + Σ_c γ_c² dS_c/dE) (Park, Phys. Rev. C
104 (2021) 064612), the same for every channel of the level; `scripts/widths_check.py` prints it per level (in the
Park basis `Session.park_norms(x)` gives it). 1 − J is the share of a bound state outside a_c (Thomas normalization),
by analogy for a resonance. J well below 0.5 marks a saturated level: its observed widths look ordinary while its
formal pole has run away (a ⁷Be evaluation: a 7/2⁻ level went from J = 0.91 to 0.05 when a channel opened, with
θ² = 0.99 and 0.45 that looked unremarkable). Background poles are exempt; say so. To bound J in a fit, a prior or
penalty must be stiff: at J ≈ 0.05 J is nearly flat in the widths, and a soft penalty leaves the saturated minimum in
place (azure2-fit). A level that sits on such a bound has a width set by the bound, not by the data.
