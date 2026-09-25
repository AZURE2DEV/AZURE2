# identical_pp_mott — p+p Coulomb scattering is the spin-1/2 Mott formula

Two identical spin-1/2 particles, no nuclear interaction worth the name: the
only level is a zero-width 1G4 (l = 4), whose hard-sphere phase goes as
(ka)^9 and is below 1e-8 of the cross section here. What is left is the
Coulomb scattering of two identical fermions,

    dsigma/dOmega = (eta/2k)^2 [ csc^4(t/2) + sec^4(t/2)
                     - (1/2) 2 cos(eta ln tan^2(t/2)) / (sin^2(t/2) cos^2(t/2)) ],

the interference term weighted by (-1)^(2j)/(2j+1) = -1/2. It comes from the
exchange term entering channel spin s with (-1)^s: singlet (s = 0) symmetric,
triplet (s = 1) antisymmetric, averaged with weights 1/4 and 3/4.

## The data

`data/mott_<angle>.dat` hold **1.01 x the closed-form Mott cross section**
with a **1 % uncertainty**, at E_cm = 0.03-0.60 MeV (E_lab = 2 E_cm) and
theta_cm = 20, 45, 70, 90, 110, 135, 160 deg, computed with AZURE2's own
constants (include/Constants.h). A calculation that *is* Mott therefore
gives chi2 = 1 per point exactly, and the recorded chi2 = 140 over 140 points
says the deviation is below 1e-6 everywhere. The offset is there so that the
chi2 is dominated by a known number rather than by rounding noise, which
would not survive the suite's relative tolerance.

The spin-0 treatment (one exchange sign for every channel spin, interference
weight -1) gives chi2 = 3.8e5 on these data; at 90 deg it made the p+p
Coulomb cross section vanish.

The channels also exercise the identical-particle channel rule, L+S even:
1G4 (L = 4, S = 0) is allowed. The rule the code used before,
(-1)^(L+S) = -1 for fermions, is refused when the model is read.
