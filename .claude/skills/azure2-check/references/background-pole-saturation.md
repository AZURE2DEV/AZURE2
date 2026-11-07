# Saturated background poles: test and worked example

A background pole is saturated when its level-shift factor J = 1/(1 + Σ_c γ_c² dS_c/dE) is close to 0 and its
observed widths Γ_c = 2 P_c γ_c² J sit against their ceiling. χ² is then flat along the direction J → 0: the formal
amplitudes γ_c (MeV^½) are not constrained and the observed widths are set by the ceiling, not by the data. This
file gives the cap test and one worked example. The rule is in the parent `SKILL.md`.

## The cap test

One refit per cap k:

1. Bound each of the pole's amplitudes in `<parameterSettings>`:
   `Level N Channel k Width (eV) <start> -c c 0 0 0 rwa -1`, with c = k γ_W in MeV^½ and category `rwa`
   (`Level N Channel M` is `width_N_M` of `param.sav`: `docs/source/user_guide/fitting.rst`).
2. Start from the fit's `param.sav`, with the capped amplitudes clipped to 0.9 c.
3. Scan k = 100, 30, 10, 3. Record total χ², Δχ² against the unbounded fit, J, the amplitudes in units of γ_W, and
   the per-data-set χ² (which data pay for a tighter cap).

Reading: if the loosest cap costs nothing, prefer the capped solution, or report the pole as unconstrained. The k
at which χ² starts to rise measures how much background the data require. This is a different question from
holding all background poles to their Wigner limits, which can cost thousands in χ² on the same fit; one saturated
pole can usually be capped at no cost even when the background as a whole cannot be held to 1 γ_W.

## Worked example: 11B+α, 5/2⁺ background pole at E_x = 30 MeV

Unbounded fit: Brune amplitudes 743 (α), 344 and 256 (n, channel spin s = 1/2 and 3/2) and 52 (p) in units of
γ_W; Σ γ² dS/dE = 11763, J = 8.5 × 10⁻⁵.

| cap k | total χ² | Δχ² | J | amplitudes / γ_W (α, n s=1/2, n s=3/2, p) |
|---|---|---|---|---|
| none | 9925.60 | — | 8.5e-5 | 743, 344, 256, 52 |
| 100 | 9925.67 | +0.07 | 5.0e-3 | 97, 45, 33, 6.8 |
| 10 | 9938.12 | +12.5 | 0.28 | 10.0, 5.4, 3.5, 0.8 |
| 3 | 10014.62 | +89.0 | 0.59 | 3.0, 3.0, 1.7, 0.5 |
| 10 on n s=1/2 only | 9928.29 | +2.7 | 0.091 | 22, 10.0, 7.4, 1.5 |

The data require a strong 5/2⁺ background (the cost starts between k = 30 and 10, paid in the neutron-total and
(α,n) data), not a saturated one: a 60 times less saturated pole costs Δχ² = 0.07. Capping only the dominant channel
pulled the others down with it (α from 743 to 22 γ_W): the channels did not compensate each other, they moved
together along the J → 0 flat direction.

## Why a stored `.azr` can fail to transform

For a saturated pole the observed widths scale as 1/Σ γ² dS/dE, so a 1 % change in how dS/dE is computed moves them
by 1 %. An `.azr` that carries observed widths written by an older build (0.98 % above what a more accurate dS/dE
route gives for the same amplitudes) makes 1 − Σ negative, and AZURE2 prints
`Denominator less than zero while transforming`. The fit itself is unaffected: the `.azr` level values are only
starting values, and the fit is the `.azr` together with its `param.sav` (`azure2-fit`).
