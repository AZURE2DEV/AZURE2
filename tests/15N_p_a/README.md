# 15N_p_a — a second compound nucleus: 16O, THM alongside direct data

15N(p,alpha)12C through the two interfering 1- levels of 16O — the CNO-cycle
leak reaction and a classic THM benchmark.

Model (built from scratch, not inherited from the 8Be evaluations):
pair 1 = p + 15N (S_p = 12.1274 MeV, a = 4.85 fm), pair 2 = alpha + 12C
(S_a = 7.16192 MeV, a = 5.43 fm); levels 1- 12.44, 1- 13.09, a fixed 1-
background at 17.0, and the two 3- levels 13.129 and 13.259 -- the second
is the narrow (Gamma ~ 21 keV) resonance Schardt resolves at
E_cm ~ 1.15 MeV. The file carries the *jointly fitted*
parameters: the fit ran against all three segments at once and the two 1-
alpha amplitudes emerge with opposite signs — the well-known destructive
interference between the 312 and 962 keV resonances.

Data:
- data/lacognata_a.dat, lacognata_b.dat — M. La Cognata et al., PRC 76,
  035801 (2007), EXFOR C17880031/2: the THM-derived on-shell cross section
  (two spectator-momentum selections), E_cm = 20-500 keV. EXFOR carries no
  point errors; 10% assigned. These are *derived* sigma, already
  penetrability-corrected, so unlike 6Li_d and 7Li_p_a they enter as
  ordinary segments with a free normalization and no penalty (the THM scale
  is arbitrary) rather than through the HOES machinery.
- data/schardt.dat — A. Schardt et al., Phys. Rev. 86, 527 (1952), EXFOR
  C0644002: 107 absolute points, E_cm = 0.10-1.50 MeV, covering both
  resonances; normalization uncertainty 15%.

The pinned chi2 = 4601.32/157 is deliberately *unconvolved*: the THM
excitation function carries an energy resolution that must be folded in
S-factor space (smearing sigma directly is wrong by many orders at the
Gamow tail, where the kernel leaks the far larger high-energy cross
section into the lowest points). The engine's <targetInt> convolution
works on sigma, so the folding lives in the analysis pipeline instead
(thm_15N/fit_n15.py), where the resolution fitted to sigma = 71.4 keV and
the joint chi2 is 810.6/157 with the Schardt block at 149.5/107 — the same
149.545 this pin shows, since no smearing touches the direct segment.
