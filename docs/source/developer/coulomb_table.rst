Universal Coulomb-function tables
=================================

Status: stage 1 of 4 (the table and its tool exist; the engine does not read
them yet). Design from the natural-neon work of October 2026
(``17O+a/10-13-26_natNe_seg1/coulomb_table_spec.md`` in the fit archive).

Why
---

A calculation spends most of its time on the channel functions. A ``perf``
profile of one cold evaluation of the natural-neon example (mode 1, four
threads, 10.5 min) gave about 40 % to the Coulomb functions of open channels
(Michel's routine, mostly for the charged α+¹⁷O channels), about 29 % to the
shift function of closed channels (a numerical derivative of GSL's Whittaker
function) and about 5 % to the R-matrix algebra. A data point pays this once;
the cost returns whenever its energy moves:

- a fitted energy shift moves every point of its segment on every evaluation;
- the adaptive sub-point grid of a ``<targetInt>`` line is placed once, from
  the level energies in the file, and never re-placed, because a new grid
  means new channel functions at every sub-point. A fit that moves narrow
  levels therefore minimizes against the start model's grid (natural neon
  refit1: 25,773.96 against 25,776.08 for the same parameters).

Closed neutral channels already use an exact recursion
(``ShftFunc::NeutralBoundShift``, 2026-10-18), which removed most of the 29 %
for the natural-neon example. Open neutral channels need nothing: Michel's
routine at η = 0 is already as fast as the Riccati-Bessel closed form. The
table is for charged channels.

Two variables make the table universal
--------------------------------------

For a channel with charges :math:`Z_1, Z_2`, reduced mass :math:`\mu`,
energy :math:`E` and radius :math:`a`, everything AZURE2 takes from the
Coulomb functions depends only on :math:`\eta = Z_1 Z_2 e^2 \mu / \hbar^2 k`,
:math:`\rho = k a` and :math:`l`:

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Quantity
     - Definition
     - Tabulated as
   * - :math:`P_l`
     - :math:`\rho / (F_l^2 + G_l^2)`
     - :math:`\ln P_l - \ln C_0^2(\eta)`, with :math:`C_0^2 = 2\pi\eta / (e^{2\pi\eta} - 1)`
   * - :math:`S_l`
     - :math:`\rho (F_l F_l' + G_l G_l') / (F_l^2 + G_l^2)`
     - :math:`S_l`
   * - :math:`\phi_l`
     - :math:`\arg(G_l + i F_l)`, continuous from :math:`\phi_l(0) = 0`
     - :math:`\phi_l - \rho`
   * - :math:`S_l` (closed)
     - :math:`z W'(z)/W(z)`, :math:`W = W_{-\eta_b, l+1/2}(z)`, :math:`z = 2\rho_b`
     - :math:`S_l + \rho_b`

So one table per :math:`l` over :math:`(\log_{10}\eta, \log_{10}\rho)` serves
every pair, charge, mass, radius and energy. A closed channel uses the same
table axes with :math:`\eta_b` and :math:`\rho_b = \kappa a`. :math:`F` and
:math:`G` are never tabulated: above the turning point they oscillate, while
the quantities above vary slowly. The phase is made continuous with
:math:`d\phi_l/d\rho = P_l/\rho` (the Wronskian is 1): at each node the branch
of :math:`\arg(G + iF)` nearest a Simpson integral of :math:`P_l` from the
previous node is taken.

A channel moves along a line of constant :math:`\eta\rho = 0.0344\,Z_1 Z_2\,
\mu[\mathrm{u}]\,a[\mathrm{fm}]`: about 0.8 for p+¹²C, 7 for ¹²C+α, 10 for
¹⁷O+α, 45 for ¹²C+¹²C and 120 for ¹⁶O+¹⁶O. The table is limited to
:math:`\eta\rho \le 150` (a build option); outside it, and outside the
:math:`\eta`, :math:`\rho` and :math:`l` ranges, the engine will call the exact
routine.

Reference routines
------------------

- Open channels: Michel's routine (``Coulomb_wave_functions_accurate``,
  AZURE2's default). On sampled points it agrees with ``mpmath`` at 30 digits
  to 1e-9 – 1e-12.
- Closed channels: :math:`S = -\rho_b + (l+1) - 2 a \rho_b\, U(a+1, b+1, z) /
  U(a, b, z)` with :math:`a = l + 1 + \eta_b`, :math:`b = 2l + 2`, from
  :math:`U' = -a\, U(a+1, b+1, z)` and GSL's scaled ``hyperg_U``. Exact, with
  no numerical derivative; at :math:`\eta_b = 0` it equals the neutral
  recursion to 1e-14.
- GSL's Coulomb functions (``--gsl-coul``) are **not** a usable reference:
  against Michel they are wrong by up to 3–6 % in :math:`P` and :math:`S` for
  :math:`0.1 \le \eta < 10` and :math:`1 \le \rho < 10` (0.03–0.3 % for
  :math:`\eta \ge 1`, :math:`\rho \ge 10`), and correct to ~1e-8 for
  :math:`\eta < 0.1` or :math:`\rho < 0.1`.

Interpolation
-------------

Tensor-product Lagrange interpolation through 8 nodes per axis (orders 6 and
bicubic Hermite are kept for comparison, ``--order``). It uses only the node
values, so the file stores nothing else. A lookup is a pure function of
:math:`(\eta, \rho, l)` and of the table file; it never depends on earlier
evaluations. That matters: ``CoulFuncCache`` once interpolated between energies
it had happened to cache, and was made exact-match because those
history-dependent values broke the adaptive quadrature of the reaction rate
(``CoulFuncCache.cpp``).

Stage-1 validation
------------------

``azure2-coulomb-table validate`` compares the table with the exact routines
on random points: uniform in :math:`(\log\eta, \log\rho)`, plus 5 % near the
turning point :math:`\rho_t = \eta + \sqrt{\eta^2 + l(l+1)}` and 5 % within two
cells of an edge. Errors: :math:`P` relative, :math:`S` absolute (relative
above :math:`|S| = 1`), :math:`\phi` absolute in radians; :math:`P` below
1e-250 (the exact value underflows) is skipped. Sweep of 2026-10-18,
:math:`l = 0..10`, :math:`\eta, \rho \in [10^{-3}, 100]`,
:math:`\eta\rho \le 150`, order 8, 400,000 points each:

.. list-table::
   :header-rows: 1
   :widths: 14 12 14 14 14 14 18

   * - Nodes per decade
     - File
     - :math:`P`
     - :math:`S`
     - :math:`\phi` (rad)
     - Closed :math:`S`
     - Same, :math:`\eta\rho \le 50`
   * - 60
     - 32 MB
     - 8.8e-7
     - 2.8e-6
     - 4.3e-7
     - 2.8e-8
     - :math:`S` 9.8e-7
   * - 80
     - 57 MB
     - 1.2e-7
     - 4.7e-7
     - 7.2e-8
     - 8.6e-9
     - :math:`S` 1.2e-7
   * - 100
     - 88 MB
     - 2.6e-8
     - 5.9e-8
     - 8.9e-9
     - 8.8e-9
     - :math:`S` 2.0e-8

The worst errors sit at the Coulomb turning point at :math:`\eta \approx 7`,
where the functions change over the Airy width while :math:`\ln P` changes by
about :math:`2\pi\eta`; they are spread evenly over :math:`l`. At the same
density, orders 4 and 6 are 5–30 times less accurate. 100 nodes per decade
meets the 1e-7 acceptance limit everywhere; the 100-per-decade table takes
5.3 min to build on four threads.

Timing in the validator, random access over the whole table: 1.5–3 µs per
lookup of the open and the closed values together, against 33–38 µs for the
exact routines (averaged over the domain; Michel's routine is slowest at large
:math:`\eta`). The engine's accesses are ordered in energy, so the real gain is
measured in stage 2.

The tool
--------

.. code-block:: console

   $ azure2-coulomb-table build table.bin --per-decade 100 --threads 4
   $ azure2-coulomb-table validate table.bin --n 400000 --threads 4
   $ azure2-coulomb-table info table.bin

``build`` takes ``--lmax`` (10), ``--eta`` and ``--rho`` ranges (1e-3 to 100),
``--per-decade`` (30) and ``--eta-rho-max`` (150). It writes ``FILE.new`` and
renames it over ``FILE``, so a job that has the old file open keeps reading
it. ``validate`` takes ``--n``, ``--seed``, ``--order``, ``--limit`` (1e-7) and
``--eta-rho-max`` (test only that part of the domain). The file is binary,
little-endian, with a header giving the format version, the axes, ``l_max``,
the :math:`\eta\rho` limit and the generator. ctest ``coulomb_table_small``
builds and validates a small table in a few seconds.

Remaining stages
----------------

2. Engine lookup behind ``--coulomb-table[=path]`` / ``--no-coulomb-table``
   (and the pyazr option), at the point where :math:`P`, :math:`S` and the
   hard-sphere phase are taken per channel (``EPoint::CalcEDependentValues``),
   default off; regression suites to 1e-6 in χ²; timing.
3. Default on; re-placing the target-integration grid as levels move; the
   refit1 fixed-point check.
4. User documentation and the CLI reference.

Decisions still open: where the table file lives and whether a copy ships;
the density (the table above) and :math:`l_{max}`; re-placing the grid every
evaluation or only when an energy moved; whether ``save_fit`` should write
fitted norms and shifts into the segment lines in the same change.
