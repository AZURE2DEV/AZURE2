# target_auto_windows

Automatic target-integration windows from the level scheme (`"auto"` in the windows token, 2026-11-06).

Outside the energy windows of a thin-target integration a point is evaluated once, at the depth-averaged energy
E_bar of its target (`tests/target_window_fallback`). With `"auto"` the engine places the windows itself, per
segment from that segment's entrance pair (`EData.cpp` `TargetWindowsFor`):

- around every resonance narrow enough for sigma(E_bar) to miss its peak by more than the tolerance: for a
  Lorentzian of total width Gamma under a target of thickness Delta (energy loss plus twice the beam spread) the
  peak error is Delta^2/(3 Gamma^2), so a level is narrow when Gamma < Delta / sqrt(3 tol);
- around every particle threshold of the same compound nucleus (a target that straddles a threshold sees the
  cusp; 15 % at the 16O(6.05) n1 threshold in the 13C+alpha ND2021r data).

Each window reaches margin x max(Gamma, Delta) beyond the points whose target holds the resonance or threshold:
[E_R - w, E_R + Delta + w]. `"auto:<tol>:<margin>"` sets the numbers (defaults 0.01 and 3); explicit ranges in
the same token are added. A Note per segment lists levels, thresholds, windows and points integrated. The windows
are fixed for the session. Older engines read no range from `"auto"` and integrate everywhere.

Model: `tests/thick_target_resonance` (p+14N, 278 keV lab resonance), 57 points 0.240-0.520 MeV lab, one capture
pair, density 1e18, stopping a0/x^0.8.

| check | what |
|---|---|
| (a) | `"auto"`: the Note reports tolerance 0.01, margin 3, one narrow level, no threshold, some points left out (26 of 57 integrated) |
| (b) | output byte-identical to the line with the listed windows written out; the inside points are the full integration to every digit |
| (c) | every point within 1e-2 of the full integration (4.2e-3 at the window edge) |
| (d) | `"auto:0.001:1"` reported; fewer points integrated |
| (e) | a threshold (a pair whose heavy partner is excited by 0.42 MeV, zero ANC: sigma unchanged) gets a window |
| (f) | `"auto,0.48-0.50"`: the explicit range is added |
| (g) | data far from every level: no window, every point at E_bar |
| (h) | fitted density and a 10 keV blend: the gradient check stays on the analytic block, no disagreement |
| (i) | refused with an ERROR: `"auto"` on a thick target and on a line without a target integration |

Builds before the feature ignore `"auto"` and integrate every point: they fail 10 checks, all but (c).

On the 13C+alpha ND2021r model (16 lines, 11,459 points): full integration 400,616; hand-made windows (10-30)
401,435 with 48 % of the points integrated; `"auto"` 400,881 with 58 % (largest point shift 0.42 sigma);
`"auto:0.003"` 400,656 with 81 %; `"auto:0.001:5"` 400,621 with 99 %.

    ./tests/target_auto_windows/check.sh path/to/AZURE2      # ~5 s
