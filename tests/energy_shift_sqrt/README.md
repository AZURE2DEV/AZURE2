# energy_shift_sqrt — the sqrt(E) term of a segment's energy shift

A segment line may end with `sqrtshift <b> <bError> <vary>`, which shifts the
segment's energies as E' = E + a + b*sqrt(E/MeV) (a is the ordinary constant
shift, b in MeV^1/2): the form an additive offset in an analyzing-magnet field
reading produces (E = kB^2, so dE ~ sqrt(E)), as opposed to the constant a.
Introduced 2026-10 after the 12C+a 15N(p,a1 gamma) 0-deg data of Bashkin
needed -15 keV at 3 MeV and nothing at 1.2 MeV.

Model: tests/identical_pp_res (p+p, E_cm 0.4-5 MeV, its own values as data).

Segments (chi2 = 24 per segment in the pristine model):

1. `xs_40.dat` with a = 0.004 MeV and b = 0.003 MeV^1/2, both fixed;
2. `xs_40_preshifted.dat`: the same yields written at E + 0.004 + 0.003*sqrt(E),
   no shift terms -- must give exactly segment 1's chi2 (the shift is applied
   to the original energies, so the two are the same data);
3. `ay_60.dat` with b = -0.002 (analyzing power, isDiff 7) and
4. `ay_60_preshifted.dat`, its twin;
5. `xs_90.dat`, untouched: still 24.

`expected/chiSquared.out` pins all five.  `check.sh` additionally requires
segments 1 = 2 and 3 = 4 to 1e-8, segment 1 != 24 (the shift is really
applied), and then fits a and b: `xs_40_mislabelled.dat` carries labels L
with L - 0.006 + 0.005*sqrt(L) = E exactly, every level fixed, a and b free
from 0 (penalty widths 1, i.e. unconstrained); the yields are the model's exact
values (the pristine 1 % offset divided out), so the truth is the chi2 = 0
minimum.  MIGRAD must return a = -0.006, b = 0.005 (|error| < 1e-4) and
chi2 = 0, once with its numerical gradient and once with `--use-gradient`
(AZURECalc::Gradient, where both shift terms are finite-differenced).
