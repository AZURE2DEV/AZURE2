# thm_spectator_window — averaging the THM vertex over the spectator momentum

`check.sh` runs `tests/18O_p_a_thm` (two 1/2+ levels of 19F, l = 0 entrance,
two THM segments with free norms, ~0.1 s a run) with `experiment[<name>] ...
ps=...` in a `<thm>` block appended to a temporary copy of the .azr
(docs/source/theory/thm_implementation.rst, "Spectator-momentum window").  Its
reaction, 2H(18O,α15N)n, has a deuteron Trojan horse: x = p, s = n,
μ_sx = 469.46 MeV, kinematics `beam=18O target=d spectator=n Ebeam=54`; over
the data the kinematics reach |p_s| from |k_sF − β k_aA| (3-8 MeV/c) to
beyond 130 MeV/c.

At fixed E the spectator direction fixes q = |p_s|, and the three-body phase
space is d cos θ_cm = q dq/(β k_sF k_aA): the model is the mean over the
accepted directions with the weight |φ(q)|² q dq.  Until October 2026 the
weight was |φ|² p² dp on the whole window (the measure of events integrated
over E as well); the pins of (b)-(d) changed with it:

| case | old → new |
|---|---|
| (b) one node vs `spectatorEnergy` | exact → 4e-10 (q is recovered from cos θ_cm) |
| (c) flat table vs Simpson | weight 1 → weight q (a table is now \|φ\|², the measure is applied); 8e-8 → 1.6e-8 |
| (c) Hulthén [0, 40], 16 vs 32 nodes | 1e-10 → 1.3e-10 (8 nodes: 2e-8 → 1.6e-10) |
| (d) the vertex node at 0.70 MeV filled to | 0.42 → 0.48 of the maximum |

| case | expectation |
|---|---|
| (a) `ps=delta` | all output files byte-identical to the same line without the key; also without kinematics, and with `spectatorEnergy=0.4` |
| (b) `ps=hulthen:30-30`, `hulthen:29.9999-30.0001`, `gauss:60:29.9999-30.0001` | model = `spectatorEnergy` = 30²/2μ_sx at all 118 points to 1e-8 (measured 4e-10) |
| (c) `ps=table:` with constant \|φ\|² on [20, 40] MeV/c | model = (∫ q σ dq)/(∫ q dq) by Simpson over 41 single-`spectatorEnergy` runs (h = 0.5 MeV/c) to 1e-6 (measured 1.6e-8, the Simpson error); Hulthén [0, 40]: 16 nodes vs 32 to 1e-6 (measured 1.3e-10); the node tables (`ps_table`) linear in the nodes, < 4 kB per point at 32 |
| (d) a point p_s = 24 MeV/c vs the Hulthén window [0, 40], no folding | the vertex node at E ≈ 0.70 MeV (model/max = 2.4e-4) is filled (0.48 of the maximum) |
| (e) refusals | a window without kinematics; `hulthen:40-20`, a negative p, an unknown kind, b < a, FWHM 0; a missing, unordered, all-zero or one-row table; `psNodes=65`, `psNodes` without a window; `ps` given twice; a window with `spectatorEnergy` or `spectatorEnergy[1]`; a reversed `theta` window; a window out of reach (`hulthen:200-300`); a window with `spectatorAngle`; `psNodes` with `spectatorAngles`; `spectatorAngles` with neither a distribution nor a distortion |
| (f) one acceptance | `ps=hulthen:0-40` byte-identical to the same with `spectatorAngles=cm:0-180 spectatorAngleNodes=16`; `cm:160-180` with the cut changes the model (q ≤ 40 MeV/c needs θ_cm ≳ 147°) |

`tests/pyazr/thm_spectator_window_test.py` checks the same through pyazr
(`AzrModel.set_thm_experiment(..., ps=..., psNodes=...)`, CLI == session,
`session.thm_vertex` at every energy against an independent evaluation of
the three-body kinematics -- nodes in cos θ_cm, weights |φ(q)|², T_s, ρ to
1e-10 -- and one direction, `ps=hulthen:30-30`, equal to a `spectatorEnergy`
session to 1e-9).
