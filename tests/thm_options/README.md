# thm_options — the keys of the `<thm>` block

`check.sh` runs `tests/7Li_p_a` (THM segment, entrance pair 5 = 7Li+p with two
entrance l in one channel spin, 30 keV Gaussian) once per option, each with a
`<thm>` block appended to a temporary copy of the .azr, and compares the total
chi2 in `output/chiSquared.out`:

| block | chi2 |
|---|---|
| none / all defaults written out (`vertex=constant`) / comment only | 1753.15 (identical) |
| `vertex=perlevel`, `vertex=real` | 2180.69 (identical to each other; the pin of the default before 2026-09-25) |
| `entranceL=coherent` | 2706.00 |
| `entranceL=coherent` + `vertex=perlevel` | 2740.48 (legacy pin of AZURE2 before Sep 2026; measured 2740.55, rel 3e-5) |
| `kinematics=lacognata` / `triple` / `kf3body` / `lambda32` | 1753.15 (identical to none) / 1964.94 / 2191.45 / 5386.89 |
| `vertex=onshell` | 2070.17 (independent of the half-off-shell vertex choice) |
| `coulombIntegral=1` | 1641.34 |
| `spectatorEnergy=0` | identical to none |
| `spectatorEnergy=0.5` | 1713.61; `spectatorEnergy[5]=0.5` identical to it |
| `spectatorEnergy[3]=0.5` (not the entrance pair) | identical to none |

A misspelt key (`vertx=onshell`), an unknown value (`kinematics=kf2body`), a
negative spectator energy and an unterminated block must each print an
`ERROR: <thm> ...` line, make AZURE2 exit non-zero (255) and leave no
`chiSquared.out`.

Pins measured with the parameters of `tests/7Li_p_a` on 2026-09-25, re-measured
the same day after the default vertex moved from `perlevel` to `constant`
(commit b351b90); tolerance 1e-3 relative. About 55 s with `OMP_NUM_THREADS=2`.
