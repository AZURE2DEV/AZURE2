# natural_neon — a real natural-abundance project

Three compound nuclei in one file: n+²⁰Ne → ²¹Ne, n+²²Ne → ²³Ne and
n+²¹Ne → ²²Ne. The levels and parameters are those of the 2026-10-17 joint fit
of the neon total cross sections (R-matrix archive,
`17O+a/10-13-26_natNe_seg1/joint_fit.azr`; 90 free parameters, 6,328 points,
Brune parametrization). The natural abundances (0.9048 / 0.0925 / 0.0027) are
fixed in the `<abundances>` block.

`tests/natural_target` checks the same features on a small synthetic
p + ¹²C / ¹³C model. This suite runs them on a fitted model with real data, so
it fails if a change moves a real natural-target fit.

`run_tests.sh` runs a plain calculation with the CLI defaults and compares the
total and per-segment chi-squared against `expected/chiSquared.out`. The
normalizations are fixed at the joint-fit values, so every segment reproduces
the joint fit's own chi-squared.

`check.sh` runs the four test segments, with `<targetInt>` switched off,
because it folds test segments with the same numbers as the data segments it
names. It checks two things:
- the natural sum equals 0.9048 σ(²¹Ne) + 0.0925 σ(²³Ne) + 0.0027 σ(²²Ne) on
  one lab-energy grid, with each nucleus's total taken from its own test
  segment, to 1e-9;
- the decomposition file `AZUREOut_sample_seg=1.extrap` carries the three
  weighted terms, and they add up to the total.

About 55 s on two threads.

## Segments

| # | data (EXFOR) | sample | treatment | χ² (N) |
|---|---|---|---|---|
| 1 | Vaughn et al., Phys. Rev. 118, 683 (1960) (11106003) | natural Ne | natural sum over the three nuclei | 842.214 (76) |
| 2 | Almeida & Käppeler, ApJ 265, 417 (1983), low-energy run (21808) | natural Ne | natural sum | 548.689 (376) |
| 3 | Almeida & Käppeler 1983, low-energy run (21808) | ²¹Ne, 95.4 % | ²¹Ne only, folded σ = 1 keV | 443.244 (255) |
| 4 | Almeida & Käppeler 1983, low-energy run (21808) | ²²Ne, 99.8 % | ²²Ne only, folded σ = 1 keV | 267.78 (259) |
| 5 | Sikkema 1958 (20203) | ²²Ne, 98.1 % | ²²Ne only, folded 10.6 keV | 136.554 (35) |
| 6 | Sikkema 1958 (20261) | ²⁰Ne | ²⁰Ne only, folded 10.6 keV | 254.02 (35) |
| 7 | Junghans et al., EPJ Web Conf. 239, 01006 (2020) (23755005), E_lab 0.84–1.26 MeV | natural Ne | natural sum, folded with the energy-dependent nELBE resolution | 862.32 (272) |

The Almeida resolution is not tabulated; σ = 1 keV stands for the 1 % FWHM
assumed in the fit. Sikkema is folded with 25 keV FWHM.

Segment 7 is the joint fit's 0.8–1.2 MeV c.m. window of the Junghans total.
`data/Junghans_ntotal.dat` holds only that window's rows. The other Junghans
rows, the Almeida high-energy runs and the Hinterberger ²⁰Ne total are left
out to keep the run short.

The archive's engine (the feature-branch build that ran the fit) gave the same
chi-squared to the digits shown for segments 1–4 and 7. The two folded
Sikkema segments differ by 4e-5 (136.56) and 6e-4 (254.182).
