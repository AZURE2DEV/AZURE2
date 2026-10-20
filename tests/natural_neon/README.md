# natural_neon — a real natural-abundance project

Three compound nuclei in one file: n+²⁰Ne → ²¹Ne, n+²²Ne → ²³Ne and
n+²¹Ne → ²²Ne. The levels and parameters are those of refit1q, the model of
record of the neon total cross sections since 2026-10-19 (R-matrix archive,
`17O+a/10-13-26_natNe_seg1/refit1q_fit_resaved.azr`; 86 levels, 74 free parameters
and the Hinterberger energy shift, 6,328 points, Brune parametrization).
refit1q is refit1 polished with the two levels whose widths had gone to zero
removed (²¹Ne 5/2⁺ at 8.930 MeV and 3/2⁻ at 9.130 MeV). Earlier versions of
this suite were built from refit1 (2026-10-18) and from the joint fit
(`joint_fit.azr`, until 2026-10-18). The natural abundances (0.9048 / 0.0925 / 0.0027) are
fixed in the `<abundances>` block.

`tests/natural_target` checks the same features on a small synthetic
p + ¹²C / ¹³C model. This suite runs them on a fitted model with real data, so
it fails if a change moves a real natural-target fit.

`run_tests.sh` runs a plain calculation with the CLI defaults and compares the
total and per-segment chi-squared against `expected/chiSquared.out`. The
normalizations are fixed at the fitted values, so every segment reproduces
the fit's own chi-squared. The suite has no energy shift: refit1q's only
nonzero shift is on the Hinterberger data, which are not in the suite.

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
| 1 | Vaughn et al., Phys. Rev. 118, 683 (1960) (11106003) | natural Ne | natural sum over the three nuclei | 835.555 (76) |
| 2 | Almeida & Käppeler, ApJ 265, 417 (1983), low-energy run (21808) | natural Ne | natural sum | 559.24 (376) |
| 3 | Almeida & Käppeler 1983, low-energy run (21808) | ²¹Ne, 95.4 % | ²¹Ne only, folded σ = 1 keV | 433.353 (255) |
| 4 | Almeida & Käppeler 1983, low-energy run (21808) | ²²Ne, 99.8 % | ²²Ne only, folded σ = 1 keV | 295.392 (259) |
| 5 | Sikkema 1958 (20203) | ²²Ne, 98.1 % | ²²Ne only, folded 10.6 keV | 27.9777 (35) |
| 6 | Sikkema 1958 (20261) | ²⁰Ne | ²⁰Ne only, folded 10.6 keV | 246.304 (35) |
| 7 | Junghans et al., EPJ Web Conf. 239, 01006 (2020) (23755005), E_lab 0.84–1.26 MeV | natural Ne | natural sum, folded with the energy-dependent nELBE resolution | 947.456 (272) |

The Almeida resolution is not tabulated; σ = 1 keV stands for the 1 % FWHM
assumed in the fit. Sikkema is folded with 25 keV FWHM.

Segment 7 is the fit's 0.8–1.2 MeV c.m. window of the Junghans total.
`data/Junghans_ntotal.dat` holds only that window's rows. The other Junghans
rows, the Almeida high-energy runs and the Hinterberger ²⁰Ne total are left
out to keep the run short.

The levels are refit1q's, re-saved on 2026-10-19 from its `.sav` with dev
a163c06 (`refit1q_fit_resaved.azr`, written by pyazr `save_fit`). The fit's
own `refit1q_fit.azr` is not used: it gives the closed-channel strengths of
the ²¹Ne levels (the subthreshold 1/2⁺ level and the background poles) as
observed values, which the fit's build converted from its reduced widths with
the shift-function derivative of the time, a numerical derivative later found
to be off by up to 5e-3 for closed neutron channels (exact since dev
af3075e). Read back with the exact derivative they gave slightly different
reduced widths, and segments 1, 2 and 6 came out 2–4e-5 off the fit (835.523,
559.224 and 246.314 in the suite of that version). The `.sav` holds the
reduced widths themselves. The re-saved file reads back to them, so every
segment here equals the previous suite file evaluated with `refit1q_fit.sav`
(its segment rows dropped, since they number the full model's segments), digit
for digit. The energies
carry the `.sav`'s eight digits.

The folded segments depend on the target-integration grid, which is laid out
from the current levels on every evaluation since dev a163c06. A session of
the fit's build laid it out once from its start file's levels, so the fit's
own `output_refit1q/chiSquared.out` differs for those segments (Sikkema ²⁰Ne
246.119 there).
