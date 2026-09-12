# beam_profile_kernel — the beam-profile (photodissociation) experimental effect

Two 12C(a,g0)16O differential-cm segments built from Haverson's (2026) eTPC
16O(g,a0) angular distributions (HIgS beams of 8.51 and 9.85 MeV, converted
by detailed balance), on the deBoer et al. (2017) 12C+alpha level scheme,
each folded with its `beamprofile` targetInt block: a skewed-Gaussian beam
profile (absolute, not point-centred), the 55 keV TPC energy resolution
window read from data-file columns 5-6, and the detailed-balance weight
(dbFlag = 1).  The 9.85 MeV beam covers the 0.6 keV wide 2+ at
Ex = 9.845 MeV, so the adaptive sub-point grid is exercised too.

## Where the data came from

**These are digitised points, not the author's numbers.** The thesis tabulates
only per-slice summaries (its table 8.24); the angular distributions exist in it
as raster figures (I.1-I.3), and `data/*.dat` were read off those panels: the
y-scale fixed by making each panel's bin sum match its tabulated corrected
counts, the x-scale by the bins tiling 0-180 degrees, and the uncertainties from
the drawn error bars. They were then converted to 12C(a,g0) by detailed balance,
f_db = mu c^2 E_cm / E_gamma^2 including the 16O recoil, which reproduces the
thesis' own factors to 0.06%.

So they carry a digitisation error of order a percent per bin on top of the
measurement, and they are not a substitute for the author's numeric data. They
are here because the kernel needs a realistic, awkwardly-shaped test case --
not as a reference dataset for the reaction. Cite the thesis, not this
directory:

> J. Haverson, *Photo-dissociation for studying the 12C(alpha,gamma)16O
> reaction*, PhD thesis (2026).

This pins (a) parsing of the keyword block and the optional data columns,
(b) the kernel normalisation and window shift, and (c) the lab -> c.m.
conversion of the kernel energies done once per effect.  The reference
chi2 was recorded from the run that was validated point by point against
an independent numpy implementation of the kernel
(R-matrix/12C+a_onefile/9-10-26_Haversen_test/kernel_reference.py).
