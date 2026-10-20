Universal Coulomb-function tables
=================================

Status: complete (stages 1 to 4, October 2026). The table and its tool exist; the engine reads
the table by default (``--no-coulomb-table`` turns it off) and lays the
target-integration grids out again on every evaluation
(``--no-grid-replacement`` turns that off). Design from the natural-neon work of October 2026
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
``--per-decade`` (100) and ``--eta-rho-max`` (150). It writes ``FILE.new`` and
renames it over ``FILE``, so a job that has the old file open keeps reading
it. ``AZURE2 --build-coulomb-table[=FILE]`` and the GUI's **Configure → Build
Coulomb Table...** run the same build (``CoulombTable::BuildFile``) with the
default specification, so a ready-built AZURE2, which ships without the tool,
can make the table; the GUI reports progress from ``BuildSpec::rowsDone``,
cancels through ``BuildSpec::cancel``, and calls
``CoulombTable::ForgetMissing`` afterwards so the running process loads the
new file at its next calculation (ctest ``coulomb_table_build``). ``validate`` takes ``--n``, ``--seed``, ``--order``, ``--limit`` (1e-7) and
``--eta-rho-max`` (test only that part of the domain). The file is binary,
little-endian, with a header giving the format version, the axes, ``l_max``,
the :math:`\eta\rho` limit and the generator. ctest ``coulomb_table_small``
builds and validates a small table in a few seconds.

Stage 2: the engine reads the table
-----------------------------------

``--coulomb-table[=FILE]`` (CLI) or ``use_coulomb_table=True,
coulomb_table=...`` (pyazr) sets ``Config::useCoulombTable``. Below the
existing interfaces, so no caller changes:

- ``CoulFunc::Penetrability``, ``CoulFunc::PEShift`` and the new
  ``CoulFunc::HardSpherePhase`` (used by ``EPoint`` for the per-point phase)
  look up open charged channels, with a one-deep memo so a point's three
  requests cost one lookup;
- ``ShftFunc`` looks up closed charged channels;
- neutral pairs, pairs with a nuclear potential (hybrid model) and the
  external-capture integrals, which need :math:`F` and :math:`G` themselves at
  many radii, keep the exact routines;
- the table is loaded once per process (``CoulombTable::Shared``); a missing
  or unreadable file prints one warning and the run continues with the exact
  routines; at exit the run prints how many lookups the table served and how
  many fell outside it.

Validation, 2026-10-19 (the 100-per-decade table):

- With the table off, the natural-neon example is byte-identical to the build
  before this change; ctest 27/27, ``tests/run_tests.sh`` 29/29 with the
  table off and on.
- The natural-neon example, mode 1, four threads: 473 s without the table,
  88 s with it (20.5 million lookups, 7,695 outside the table), the same χ²
  and outputs to every printed digit.
- Several suites move with the table on by more than 1e-6: 13N 9.3e-4,
  β-delayed α 6.9e-4, thick-target resonance 3.5e-4, others 1e-5 – 9e-5. The
  table is not the cause. :math:`dS/dE` enters the conversions between
  observed and formal parameters, and the exact route takes it as a 1 eV
  central difference: of :math:`S` from a numerical radial derivative of the
  Whittaker function for closed channels (3e-4 – 4e-3 relative error for
  p+¹²C), and of Michel's :math:`S` for open ones (up to 3e-5 at a few MeV).
  The table's :math:`dS/dE` is within 1e-8 of the exact value. With the two
  comparison switches below the exact route reproduces the table at every
  printed digit in every suite without a nuclear potential:

  - ``AZURE_BOUND_SHIFT_U_RATIO=1``: closed charged channels from the exact
    U ratio (``CoulombTable::ExactBound``) instead of the numerical
    derivative of :math:`W`;
  - ``AZURE_SHIFT_DE_RICHARDSON=1``: open-channel :math:`dS/dE` by a
    Richardson-extrapolated central difference with a step of
    :math:`10^{-3} E`.

  Both became the default exact route in stage 3 (below).

Stage 3: default on, and the grid follows the levels
----------------------------------------------------

Three changes, 2026-10-19:

- **The table is on by default** (``Config::useCoulombTable``; pyazr
  ``use_coulomb_table=True``). Without a table file the run says once that it
  uses the exact routines; a file that was named, or exists and cannot be
  read, is a warning.
- **The exact route is accurate.** Closed charged channels take :math:`S`
  from the U ratio (``CoulombTable::ExactBound``) and open channels take
  :math:`dS/dE` by the Richardson-extrapolated difference, so the exact
  routines and the table agree at every printed digit.
  ``AZURE_BOUND_SHIFT_NUMERICAL=1`` and ``AZURE_SHIFT_DE_1EV=1`` restore the
  earlier routes for comparison.
- **Sub-point grids are re-placed** (``EData::ReplaceSubPointGrids``, called
  from ``FillEnergyShiftsFromParams`` on every evaluation;
  ``Config::replaceSubPointGrids``, CLI ``--no-grid-replacement``, pyazr
  ``replace_grids``). Each adaptive grid remembers how it was laid out
  (``EPoint::SubGridSpec``: window, node budget, step, width multiplier). The
  current fit values are converted to observed energies and widths on a
  private copy of the compound nucleus (``CNuc::TransformOut``, then
  ``GetTransformParams``), which is what a saved file would hold, and the same
  generator lays the grid out again. A point whose grid is unchanged (every
  node within 1e-12 relative) keeps its sub-points; the others get new ones,
  energised at once or, in a segment with an energy shift, by the shift step.
  The resonances are offset by the point's c.m. shift, so the dense part of
  the grid lands on each resonance after the sub-points move with the data.
  Segments that need external-capture amplitudes keep their grids, because
  those integrals are computed once at initialization. The four physical-
  parameter entry points (pyazr ``UpdateSegments`` and
  ``CalculateChi2Physical``, the two MCMC likelihoods) now fill the compound
  nucleus before the shifts, so the grids see the current levels.

The fixed-point check of the spec (natural neon, the objective at a fit's end
point in a session started from the fit's start file, A, and from its saved
file, B):

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - Model
     - A, re-placed
     - B, re-placed
     - A, placed once
     - B, placed once
   * - refit1
     - 25,776.076357
     - 25,776.076357
     - 25,773.959055
     - 25,776.084014
   * - refit1q
     - 25,693.609423
     - 25,693.609423
     - 25,691.859348
     - 25,694.137433

With the grid placed once the objective depended on the file the session
started from by 2.1–2.3; re-placed, the two agree to 1e-14. The objective a
fit reports is then a function of its parameters alone.

Cost: laying 29,322 grids out again takes 0.2–0.5 s per evaluation on four
threads when the levels have not moved; rebuilding them costs about as much
as a cold evaluation of those sub-points with the table (3,664 grids,
604,562 sub-points: 23 s, once, when a session starts from a file whose levels
are far from the point evaluated). The natural-neon example, mode 1: 85 s
re-placed, 82 s placed once, the same χ² to every printed digit. A fit moves
the levels a little at every step and rebuilds only the grids near the
resonances that moved.

The analytic gradient treats the grid as fixed for the step it differentiates;
the change of the quadrature nodes with the parameters is not part of it.

Stage 4: user documentation
---------------------------

The user guide describes the table under **Configure → Runtime Options**
(:doc:`../user_guide/configure_menu`) and the re-placed grids, with the two
adaptive-grid settings of the dialog, under **Experimental Effects**
(:doc:`../user_guide/experimental_effects`); the command-line options are in
:doc:`../reference/command_line`.

The table does not ship with AZURE2 (decided 2026-10-20): each user builds it
once with ``AZURE2 --build-coulomb-table`` or the GUI's **Configure → Build
Coulomb Table...**, about 5 minutes, and keeps it in
``~/.azure2/``. The file is 88 MB, and building it takes only minutes, so a
copy in the repository or the release is not worth its size.

Still open: whether ``save_fit`` should write fitted norms and shifts into the
segment lines.
