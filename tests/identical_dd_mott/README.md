# identical_dd_mott — d+d Coulomb scattering is the spin-1 Mott formula

The spin-1 counterpart of `identical_pp_mott`: two identical bosons of spin 1,
a single zero-width l = 4, S = 0 level (hard-sphere contribution below 1e-7
here), so the calculation must be

    dsigma/dOmega = (eta/2k)^2 [ csc^4(t/2) + sec^4(t/2)
                     + (1/3) 2 cos(eta ln tan^2(t/2)) / (sin^2(t/2) cos^2(t/2)) ].

Channel spins 0, 1, 2 carry exchange signs +, -, + with weights 1, 3, 5 out
of 9, so the interference weight is (1 - 3 + 5)/9 = +1/3 -- not the +1 of
spin-0 bosons that the code applied to every identical boson before.

Data: **1.01 x the closed form, 1 % uncertainty**, E_cm = 0.03-0.60 MeV
(E_lab = 2 E_cm), theta_cm = 20-160 deg; chi2 = 1 per point means Mott to
better than ~1e-6 (the largest single-point deviation here is the l = 4 hard
sphere at 0.6 MeV, a few 1e-7). The old treatment gives chi2 = 1.2e5.
