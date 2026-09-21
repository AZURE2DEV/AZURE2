# identical_pp_res — p+p with resonances: cross section and analyzing power

Identical spin-1/2 particles with nuclear scattering in every channel the
symmetry allows below l = 4: 1S0, 3P0, 3P1, 3P2 (mixed with 3F2) and 1D2,
all with sizeable widths, levels at E_cm = 1.0-2.5 MeV. It pins the two
routes that must symmetrize identically:

* the differential cross section (isDiff 4, 40/90/140 deg), built by the
  Blatt-Biedenharn route with the Coulomb amplitude taken per channel spin,
  f_C(t) + (-1)^s f_C(pi - t), and the nuclear terms x4 / interference x2;
* the vector analyzing power (isDiff 7, 60/120 deg), built from the
  channel-spin amplitude matrix, whose elastic pathways carry
  1 + (-1)^(l'+s') and whose Coulomb term is the same per-s amplitude.

The data are **this model's own values**: cross section x 1.01 with a 1 %
uncertainty, A_y + 0.01 with an absolute 0.01, E_cm = 0.4-5.0 MeV, so chi2 = 1
per point is the reference and any drift shows as a departure from 120. They
are a regression pin, not an independent check. The independent checks were
made when the change went in: the two routes agree (sigma / amplitude-matrix
spin sum constant in angle to 2e-9); A_y(pi - theta) = -A_y(theta), which
exchange symmetry demands for a polarized beam on an unpolarized identical
target (the Madison normal flips when k' -> -k'), holds to 4e-9 -- the
A_y behind the 120 deg data is the negative of that behind the 60 deg data;
and the
analytic parameter gradients of both observables match central differences.
The previous code gives chi2 = 2.3e5 on these data.
