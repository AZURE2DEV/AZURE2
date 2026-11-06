# thm_rwa_grid — convergence of a convolution over RWA-given resonances

Runs `tests/17O` (THM segment, 21 keV Gaussian, levels whose widths are given
as reduced-width amplitudes via the `gammaIsRWA` field) at 150 and 600
integration points and requires the total chi^2 to agree to 1%.

Before the fix the adaptive grid read the amplitudes (MeV^1/2) as widths in eV,
underestimated a 46 eV resonance at E_cm = 74 keV as 0.7 eV, and the convolved
cross section overestimated that resonance's area by up to 10x in a way that
depended on the grid density: chi^2 13.8, 53 and 88 at 150, 600 and 2000 points.
After it: 96.69, 96.65, 96.88, and 1.6% agreement with an independent
fine-grid Gaussian fold everywhere.
