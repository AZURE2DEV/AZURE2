# thm_options — the keys of the `<thm>` block

`check.sh` runs `tests/7Li_p_a` (THM segment, entrance pair 5 = 7Li+p with two
entrance l in one channel spin, 30 keV Gaussian) once per option, each with a
`<thm>` block appended to a temporary copy of the .azr, and compares the total
chi2 in `output/chiSquared.out`:

| block | chi2 |
|---|---|
| none / all defaults written out (`vertex=constant`) / comment only | 2138.48 (identical) |
| `vertex=perlevel`, `vertex=real` | 2181.45 (identical to each other; the default before 2026-09-25) |
| `entranceL=coherent` | 3196.77 |
| `entranceL=coherent` + `vertex=perlevel` | 2741.85 (AZURE2 before Sep 2026, whose pin was 2740.48; the difference is the two dS/dE fixes, b6cc41b and the step fix) |
| `kinematics=lacognata` / `triple` / `kf3body` / `lambda32` | 2138.48 (identical to none) / 2381.01 / 2625.33 / 5499.08 |
| `vertex=onshell` | 2070.75 (independent of the half-off-shell vertex choice) |
| `coulombIntegral=1` | 2112.45 |
| `spectatorEnergy=0` | identical to none |
| `spectatorEnergy=0.5` | 2195.69; `spectatorEnergy[5]=0.5` identical to it |
| `spectatorEnergy[3]=0.5` (not the entrance pair) | identical to none |

Two more cases are not `<thm>` options but the Trojan-horse binding energy:
the same models with the carrier the experiments actually used, 3He, whose
binding B(3He -> p + d) = 5.4935 MeV replaces the pinned projects' own
(d, 2.2246 MeV, in `7Li_p_a`; 6Li, 1.4735 MeV, in `6Li_d`) in field 32 of every
entrance channel line. See the *Data provenance* sections of those READMEs.

| project, binding | chi2 |
|---|---|
| `7Li_p_a`, B = 5.4935 (3He) | 1827.73 |
| `6Li_d`, B = 5.4935 (3He) | 482.892 (tolerance 5e-3, see below) |

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
| none | 2138.48 (the none pin above) | 0.00493262 |
| w = 1 (table -1..20 MeV) | identical | identical |
| w = 2 (table 0..10 MeV), relative or absolute path | identical | doubled (0.00986523; the norm multiplies the data, n* = Smm/Smd) |
| w = 1 + E (11 rows, 0..10 MeV) | 869.55 (tolerance 5e-3, see below) | |
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

Re-pinned 2026-09-28, after the dS/dE step fix (ChannelFunc::DerivativeStep;
the central difference of the shift function had a 1 eV step and was
round-off dominated). Every pin moved by less than 2e-4 except the three that
exposed the problem on GitHub CI, where the runners build for the x86-64
baseline and this machine's compiler for x86-64-v3 (FMA): 6Li_d with
B(3He) 483.474 -> 482.892 (1.2e-3; CI had 482.681), the w = 1 + E weight
868.963 -> 869.55 (6.8e-4; CI Linux/Windows 869.756, macOS 867.011) and
entranceL=coherent 3196.17 -> 3196.77 (1.9e-4). Several pins also still
predated 02ae97b. An FMA and a non-FMA build now agree on every case to the
printed digits, and perturbing every libm result by 1e-11 moves none.

Tolerances: 1e-3 relative, the suite default, except `TOL_WIDE` = 5e-3 for
6Li_d with B(3He) and the w = 1 + E weight. Those two are the most
sensitive to the platform: the 6Li_d model with the 3He binding has a level
at the limit of the Brune transformation (1 - sum Gamma dS/dE / 2P =
-4e-4), which amplifies any difference in dS/dE, and the ramp weight leans
on the lowest points, whose folds reach through threshold. Before the fix
CI saw 1.6e-3 (Linux/Windows) and 2.2e-3 (macOS) there; 5e-3 covers that
with margin until a macOS run shows how much of it the fix removed. The
identities (`== no block`, sorted == file order, w = 2 doubling the norm to
1e-5) stay exact.

Portability (GitHub CI runs this on Linux, macOS with BSD sed/awk, and
Windows MSYS2 with a native AZURE2.exe): the non-THM copy of the segment
line is written with awk to a temporary file (BSD sed has no GNU `sed -i`),
and the absolute weight path is written with `cygpath -m` where it exists,
since a /tmp/... path inside the .azr means nothing to the native binary
(AZURE2 recognises C:/... and C:\... as absolute since 2026-09-28).

