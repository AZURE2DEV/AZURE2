# thm_options — the keys of the `<thm>` block

`check.sh` runs `tests/7Li_p_a` (THM segment, entrance pair 5 = 7Li+p with two
entrance l in one channel spin, 30 keV Gaussian) once per option, each with a
`<thm>` block appended to a temporary copy of the .azr, and compares the total
chi2 in `output/chiSquared.out`:

| block | chi2 |
|---|---|
| none / all defaults written out / comment only | 2180.69 (identical) |
| `entranceL=coherent` | 2740.48 (legacy pin, rel 1e-3) |
| `kinematics=lacognata` / `triple` / `kf3body` / `lambda32` | 2180.69 / 2402.92 / 2632.81 / 5573.25 |
| `vertex=onshell` | 2070.17 |
| `coulombIntegral=1` | 2145.79 |
| `spectatorEnergy=0` | identical to none |
| `spectatorEnergy=0.5` | 2088.25; `spectatorEnergy[5]=0.5` identical to it |
| `spectatorEnergy[3]=0.5` (not the entrance pair) | identical to none |

A misspelt key (`vertx=onshell`), an unknown value (`kinematics=kf2body`), a
negative spectator energy and an unterminated block must each print an
`ERROR: <thm> ...` line, make AZURE2 exit non-zero (255) and leave no
`chiSquared.out`.

Pins measured with the parameters of `tests/7Li_p_a` on 2026-09-25; tolerance
1e-3 relative. About 45 s with `OMP_NUM_THREADS=2`.
