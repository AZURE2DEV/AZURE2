# thm_options — the keys of the `<thm>` block

`check.sh` runs `tests/7Li_p_a` (THM segment, entrance pair 5 = 7Li+p with two
entrance l in one channel spin, 30 keV Gaussian) once per option, each with a
`<thm>` block appended to a temporary copy of the .azr, and compares the total
chi2 in `output/chiSquared.out`:

| block | chi2 |
|---|---|
| none / all defaults written out (`vertex=constant`) / comment only | 2138.51 (identical) |
| `vertex=perlevel`, `vertex=real` | 2181.41 (identical to each other; the pin of the default before 2026-09-25) |
| `entranceL=coherent` | 3196.17 |
| `entranceL=coherent` + `vertex=perlevel` | 2741.34 (AZURE2 before Sep 2026, whose pin was 2740.48; the 3e-4 is b6cc41b's dS/dE) |
| `kinematics=lacognata` / `triple` / `kf3body` / `lambda32` | 2138.51 (identical to none) / 2381.03 / 2625.34 / 5499.7 |
| `vertex=onshell` | 2070.84 (independent of the half-off-shell vertex choice) |
| `coulombIntegral=1` | 2112.47 |
| `spectatorEnergy=0` | identical to none |
| `spectatorEnergy=0.5` | 2195.75; `spectatorEnergy[5]=0.5` identical to it |
| `spectatorEnergy[3]=0.5` (not the entrance pair) | identical to none |

Two more cases are not `<thm>` options but the Trojan-horse binding energy:
the same models with the carrier the experiments actually used, 3He, whose
binding B(3He -> p + d) = 5.4935 MeV replaces the pinned projects' own
(d, 2.2246 MeV, in `7Li_p_a`; 6Li, 1.4735 MeV, in `6Li_d`) in field 32 of every
entrance channel line. See the *Data provenance* sections of those READMEs.

| project, binding | chi2 |
|---|---|
| `7Li_p_a`, B = 5.4935 (3He) | 1827.89 |
| `6Li_d`, B = 5.4935 (3He) | 483.474 |

These are regression pins of the same kind as the rest: the data are
penetrability-corrected and normalized to direct data, not raw half-off-shell
yields, so neither number is a physics benchmark.  (`run_tests.sh` takes the
first `.azr` of a project directory, so a second `.azr` beside each project
could not be pinned there; check.sh builds the variants on the fly.)

Section (i) is the energy-dependent weight `weight[k]=<file>` (a two-column
table E_cm, w multiplying the THM model of `<segmentsData>` line k, points
and folding sub-points, before the folding; log-linear interpolation). The
THM norm of `7Li_p_a` is free and so profiled out:

| weight on segment 1 | chi2 | norm (chiSquared.out) |
|---|---|---|
| none | 2138.51 (the none pin above) | 0.00493426 |
| w = 1 (table -1..20 MeV) | identical | identical |
| w = 2 (table 0..10 MeV), relative or absolute path | identical | doubled (0.00986851; the norm multiplies the data, n* = Smm/Smd) |
| w = 1 + E (11 rows, 0..10 MeV) | 868.963 | |
| `weightTest[2]` (a `<segmentsTest>` line, not read in a data run) | identical | |

`weightTest[3]` is checked in an extrapolation ("Calculate Segments Without
Data") with `<segmentsTest>` line 3 turned into a THM line (5 -> 4, isDiff
10, lab 0.5-3 MeV): with the two-row table (0, 1), (5, 4) every extrapolated
model value is the unweighted one times 4^(E_cm/5) to 1e-6.

With the table starting at 0 the 30 keV folding reaches below it
(E = -0.001 MeV): one `WARNING: <thm> weight ... outside its table` and the
end value. Refused at startup with `ERROR: <thm> weight[k]: ...` and a
non-zero exit: a missing file, energies not strictly increasing, w <= 0, three
columns, a data point outside the table (a table from 1 MeV), k beyond the
last `<segmentsData>` line, `weight[0]`, and a weight on a segment that is not
THM (a second, isDiff 0 copy of the segment line). A failure after the
configuration is read (these data-dependent checks) makes AZURE2 exit 255 as
a malformed block does; before 2026-09-26 AZURE2 exited 0 whenever the run
itself failed.

A misspelt key (`vertx=onshell`), an unknown value (`kinematics=kf2body`), a
negative spectator energy and an unterminated block must each print an
`ERROR: <thm> ...` line, make AZURE2 exit non-zero (255) and leave no
`chiSquared.out`.

Pins measured with the parameters of `tests/7Li_p_a` on 2026-09-25, re-measured
the same day after the default vertex moved from `perlevel` to `constant`
(commit b351b90); tolerance 1e-3 relative. The 3He-binding pins were measured
on 2026-09-26, after the shift-function fix of b6cc41b.

Re-pinned 2026-09-26: b6cc41b computes dS/dE without the old ~1e-9 noise,
which enters the Brune transformation of sub-threshold levels. That moved
every pin above by up to 3.2e-4 (2137.83 -> 2138.51 for no block). The
moves were inside the tolerance, but the pins are now exact again.
About 85 s with `OMP_NUM_THREADS=1`.
