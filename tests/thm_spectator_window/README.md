# thm_spectator_window — averaging the THM vertex over the spectator momentum

`check.sh` runs `tests/18O_p_a_thm` (two 1/2+ levels of 19F, l = 0 entrance,
two THM segments with free norms, ~0.1 s a run) with `experiment[<name>] ...
ps=...` in a `<thm>` block appended to a temporary copy of the .azr
(docs/source/theory/thm_implementation.rst, "Spectator-momentum window").  Its
reaction, 2H(18O,α15N)n, has a deuteron Trojan horse: x = p, s = n,
μ_sx = 469.46 MeV, kinematics `beam=18O target=d spectator=n Ebeam=54`.

| case | expectation |
|---|---|
| (a) `ps=delta` | all output files byte-identical to the same line without the key; also without kinematics, and with `spectatorEnergy=0.4` |
| (b) `ps=hulthen:30-30`, `hulthen:29.9999-30.0001`, `gauss:60:29.9999-30.0001` | model = `spectatorEnergy` = 30²/2μ_sx at all 118 points to 1e-8 (measured: exact, 6e-10, 6e-10) |
| (c) `ps=table:` with constant weight on [20, 40] MeV/c | model = Simpson average of 41 single-`spectatorEnergy` runs (h = 0.5 MeV/c) to 1e-6 (measured 8e-8, the Simpson error); Hulthén [0, 40]: 16 nodes vs 32 to 1e-6 (measured 1e-10; 8 nodes 2e-8) |
| (d) a point p_s = 24 MeV/c vs the Hulthén window [0, 40], no folding | the vertex node at E ≈ 0.70 MeV (model/max = 2.4e-4) is filled (0.42 of the maximum) |
| (e) refusals | a window without kinematics; `hulthen:40-20`, a negative p, an unknown kind, b < a, FWHM 0; a missing, unordered, all-zero or one-row table; `psNodes=65`, `psNodes` without a window; `ps` given twice; a window with `spectatorEnergy` or `spectatorEnergy[1]`; `theta` still reserved |

`tests/pyazr/thm_spectator_window_test.py` checks the same through pyazr
(`AzrModel.set_thm_experiment(..., ps=..., psNodes=...)`, CLI == session,
`session.thm_vertex` against an independent evaluation, and a two-node window
equal to the weighted sum of two `spectatorEnergy` sessions to 1e-15).
