# FRESCO benchmarks of physics/transfer

Reference numbers for `physics/tests/transfer_fresco_test.cpp` (what is compared and to what accuracy:
`docs/upgrade/transfer.md` section 9). One directory per case: `job.nin` (the FRESCO input of the reference run),
`fort4.dat` where the job reads a table on FRESCO's input file 4 (copy it to `fort.4` beside the input), and
`reference.txt` (the lines of FRESCO's output the test reads, for the reference run, a refined run and, for two
cases, further variants; the format is in the docstring of the script).

- Code: FRESCO 3.5, I. J. Thompson, github.com/I-Thompson/fresco commit c571c9e (GPL v3; not part of this
  repository), built with GNU Fortran 11.5.0, `-O2 -fallow-argument-mismatch -std=legacy`.
- Generated 2026-11-11 by `python3 physics/tests/data/make_fresco_benchmarks.py --fresco <path to fresco>`, which
  writes every file of this directory; each job runs as `fresco < job.nin > job.out` in a directory of its own.
- Constants: FRESCO's source fixes hbar c = 197.32705 MeV fm and u = 931.49432 MeV; every job sets
  `unitmass=1.0000004718 finec=137.03604742`, so that FRESCO's 2u/hbar^2 = 0.0478451 fm^-2 MeV^-1 and
  e^2 = 1.43996 MeV fm are those of `physics/Constants.h` (the test checks both from the output).
- Masses are mass numbers (u); binding energies and Q values are inputs.

| Case | What |
|---|---|
| `p15N_elastic` | p + 15N, E_lab = 5 MeV, Woods-Saxon + surface absorption + spin-orbit |
| `d15N_elastic` | d + 15N, E_lab = 10 MeV, spin 1 with spin-orbit |
| `li6c12_elastic` | 6Li + 12C, E_lab = 20 MeV, radii r0 (A_p^1/3 + A_t^1/3), strong Coulomb and absorption |
| `dn_finite_range` | 15N(d,n)16O(g.s.), E_lab = 10 MeV, finite range, post, no remnant; bound states n-p and p-15N (p1/2) |
| `dn_zero_range` | the same in zero range, D0 = -130.31 MeV fm^3/2, with recoil (`inh=2`) |
| `li6d_finite_range` | 12C(6Li,d)16O(g.s.), E_lab = 20 MeV, finite range, post, no remnant, the nuclear part of V_alpha-d as interaction (KRPOT; `fort4.dat` cancels the Coulomb term); bound states alpha-d (1 node) and alpha-12C (2 nodes) |
| `li6d_zero_range` | the same in zero range, D0 = 68.3056 MeV fm^3/2, `inh=2` |
| `li6d_finite_range_coulomb` | finite range with the whole V_alpha-d (Coulomb tail): FRESCO does not converge in its non-local range (the variants); data for the note, not compared by the test |
| `dn_bin_nonresonant_fr`, `dn_bin_narrow_fr`, `dn_bin_wide_fr` | 15N(d,n)[p + 15N] at E_lab = 10 MeV, finite range, post, no remnant, to a bin of p + 15N (kp 5, real, no spin-orbit) on [0, RSP = 40 fm]: l = 1 non-resonant (ISC = 2), l = 2 over a narrow resonance (ISC = 4), l = 0 inside a wide one (ISC = 4); NK = 1600; refined steps, rmatch 50, rnl 12 and NK 3200 as variants (`transfer_continuum_fresco_test`) |
| `dn_bin_*_zr` | the same in zero range (`inh=2`): FRESCO's zero range to a bin is not consistent with its finite range (`transfer.md` §15); data, not compared |
| `li6d_bin_wide_fr` | 12C(6Li,d)[alpha + 12C] at 20 MeV to a bin inside a wide d-wave resonance, RSP = 12 fm: FRESCO loses accuracy (ERWIN) at high L; data for `transfer.md` §15, not compared |
| `dn_prior`, `dn_prior_remnant`, `dn_post_remnant` | 15N(d,n)16O(g.s.) in the prior form, and with the full complex remnant (core-core n + 15N, kp 5): FRESCO does not converge in rnl (20, 26, 32 fm); data for `transfer.md` §17, not compared |

Generated 2026-11-12 (part 2), same command; FRESCO's bin listing is printed with `ipc=3` and read as the `bin*` lines.

Generated 2026-11-12 (part 3), same command; the `pol` lines are FRESCO's iT11 of the projectile (kqmax=1, Madison):

| Case | What |
|---|---|
| `dn_finite_range_so`, `dn_zero_range_so` | 15N(d,n)16O(g.s.) of `dn_*` with spin-orbit terms in d + 15N and n + 16O (`transfer_spin_orbit_test`) |
| `he3d_<gs|d52>_<finite|zero>_range[_so]` | 18O(3He,d)19F at E_lab = 25 MeV to the 1/2+ ground state (2s1/2) or the 5/2+ state at 0.197 MeV (1d5/2), finite range (post, no remnant, cutl = -3.2) and zero range (D0 = -152.53 MeV fm^3/2, inh=2), with or without spin-orbit in both distortions (`transfer_he3_fresco_test`) |
| `n14d_<gs|2p>_finite_range` | 12C(14N,d)24Mg at E_lab = 30 MeV to the 0+ ground state (12C + 12C, l = 0, six nodes) or the 2+ at 1.369 MeV (l = 2, five nodes), finite range, post, unsymmetrised; FRESCO with cutl = 0, cutr = -12, rnl = 16 (`transfer_n14_fresco_test`) |
