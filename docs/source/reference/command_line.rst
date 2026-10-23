Command-Line Usage
==================

While AZURE2 is primarily used through the GUI, it can also be executed from the
command line. This is useful for batch processing, remote execution on HPC
clusters, and automation.

Launching the GUI
-----------------

.. code-block:: bash

   AZURE2 [input_file.azr]

The Input file argument is optional. If provided, AZURE2 opens the GUI with the
specified project loaded.

Command-Line Mode
-----------------

.. code-block:: bash

   AZURE2 --no-gui [options] input_file.azr

The ``--no-gui`` flag launches AZURE2 in text-mode, with an interactive
interface similar to the original FORTRAN version. All calculation types
available in the GUI are accessible from the command line.

.. important::

   When running in command-line mode, the **Runtime Options** saved in the
   project file are **not** applied automatically. They must be specified as
   command-line flags each time.

Available Options
-----------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Flag
     - Description
   * - ``--help``
     - List all available command-line options. This is the only flag that
       can be used without specifying an Input file.
   * - ``--no-gui``
     - Run in command-line (text) mode.
   * - ``--no-readline``
     - Disable readline support for command-line input.
   * - ``--use-brune``
     - Use the Brune parameterization (equivalent to the GUI's
       "Use Brune formalism" option).  This is the default.
   * - ``--no-brune``
     - Use the standard Lane-Thomas parameterization (constant boundary
       conditions) instead of Brune's or Park's.
   * - ``--use-park``
     - Use Park's level-dependent boundary conditions (Phys. Rev. C 104,
       064612) in place of Brune's alternative level matrix.  The fit
       parameters are then the observed reduced width amplitudes,
       :math:`\Gamma_c = 2P_c\gamma_c^2`.  Same cross sections as the
       default; see the note below.
   * - ``--ignore-externals``
     - Ignore external width if internal width is zeroed (equivalent to the
       GUI option).
   * - ``--use-rmc``
     - Use the RMC capture formalism for neutron capture (equivalent to the
       GUI option).
   * - ``--gsl-coul``
     - Use GSL Coulomb functions (equivalent to the GUI option).  Not
       recommended: against the default routine, :math:`P` and :math:`S` are
       wrong by up to 3–6 % for :math:`0.1 \le \eta < 10`,
       :math:`1 \le \rho < 10` (see :doc:`../developer/coulomb_table`).
   * - ``--coulomb-table[=FILE]``
     - Take :math:`P`, :math:`S` and the hard-sphere phase of charged
       channels, open and closed, from the universal Coulomb table, with the
       exact routines outside it (the default).  ``FILE`` defaults to
       ``$AZURE2_COULOMB_TABLE``, else ``~/.azure2/coulomb_table.bin``; build
       it once with ``--build-coulomb-table`` (about 5 min on four
       threads, 88 MB).  Without a table file the run says so once and uses
       the exact routines.  The natural-neon example runs about four times faster.
       Pairs with a nuclear potential and neutral pairs are not tabulated
       (see :doc:`../developer/coulomb_table`).
   * - ``--no-coulomb-table``
     - Exact Coulomb functions everywhere.
   * - ``--build-coulomb-table[=FILE]``
     - Build the Coulomb table with the default specification on all OpenMP
       threads, write it to ``FILE`` (default as for ``--coulomb-table``),
       report progress in steps of 10 % and exit; no project file is needed.
       The same build as **Configure → Build Coulomb Table...** in the GUI.
   * - ``--no-long-wavelength``
     - Electric external capture without the long-wavelength approximation:
       the radial operator :math:`\bar e_L r^L` becomes the charge term of the
       extended Siegert theorem (Friar and Fallieros, Phys. Rev. C 29, 1645
       (1984)),
       :math:`L(2L+1)!!\,k_\gamma^{-L}\,[Z_1 g_L(k_\gamma r\,m_2/M) +
       (-1)^L Z_2 g_L(k_\gamma r\,m_1/M)]`,
       :math:`g_L(y) = \int_0^y j_L(t)/t\,dt`, which reduces to
       :math:`\bar e_L r^L` term by term as :math:`k_\gamma \to 0`.  The current
       tangential to :math:`r` is still left out.  Small for ordinary E1 and
       E2 (0.1--0.3 % for p + :sup:`16`\ O and :math:`\alpha` + :sup:`12`\ C
       E2), large where the long-wavelength E1 charge nearly vanishes
       (:math:`\alpha` + :sup:`12`\ C E1: 3.5--6 times the amplitude).  M1
       terms are the same in both modes.  Internal (resonant) photon channels
       are not affected.  pyazr: ``use_long_wavelength=False``.  Checked against
       an independent calculation in ``tests/ec_m1_exact``.
   * - ``--no-ec-table``
     - Compute every external-capture integral directly, at every data point
       and sub-point, and interpolate the amplitudes of shifted points
       linearly between the segment's own points (flat beyond its ends), as
       before the energy tables.  By default each capture term is tabulated in
       energy once, at nodes converged to 1e-10, and every point -- shifted or
       not, sub-point or rate energy -- is interpolated from the table
       (8-point Lagrange, 50 nodes per decade, error below 1e-9; see
       :doc:`../developer/ec_table`).  ``AZURE2_EC_TABLE_CHECK=1`` also computes
       every tabulated integral directly and reports the largest difference;
       ``AZURE2_EC_TABLE_NODES=N`` sets the node density.
   * - ``--no-grid-replacement``
     - Keep the adaptive sub-point grids of experimental effects
       (``<targetInt>``) where the input file's levels put them.  By default
       they are laid out again on every evaluation for the current levels
       and energy shifts, so a fit integrates each narrow resonance on a grid
       built for where the resonance is, and the result does not depend on
       which file the session started from.  Segments with external capture
       keep their grids either way.
   * - ``--no-transform``
     - Do not perform parameter transformations; input formal R-matrix
       parameters directly (equivalent to the GUI option).

Multiple options can be combined:

.. code-block:: bash

   AZURE2 --no-gui --use-brune --ignore-externals input_file.azr

Brune and Park parametrizations
-------------------------------

Both take the observed level energies as parameters and need no boundary
condition constants.  They differ only in how the reduced width amplitudes of a
level are normalized:

.. math::

   \gamma^{\rm Park}_{\lambda c} = \gamma^{\rm Brune}_{\lambda c}\,\sqrt{J_\lambda},
   \qquad
   J_\lambda = 1-\sum_c \left(\gamma^{\rm Park}_{\lambda c}\right)^2
   \left.\frac{dS_c}{dE}\right|_{E_\lambda}
   = \frac{1}{1+\sum_c \left(\gamma^{\rm Brune}_{\lambda c}\right)^2 (dS_c/dE)_{E_\lambda}} .

Park's level matrix is Brune's multiplied from both sides by
:math:`\mathrm{diag}(\sqrt{J_\lambda})`, which leaves the collision matrix
unchanged, so ``--use-park`` reproduces the default calculation
(``tests/park_formalism/check.sh``).  What changes is the parameter space the
minimizer works in:

* a Park amplitude is the observed width, :math:`\Gamma_c = 2P_c\gamma_c^2`
  (an ANC for a closed channel), with no factor depending on the other channels
  of the level.  Fixing, bounding or putting a prior on a width therefore acts
  on one parameter, and the parameters of a level are less correlated;
* every real Brune amplitude is a valid model, but Park amplitudes must satisfy
  :math:`J_\lambda > 0`; a level with :math:`J_\lambda \le 0` has no standard
  R-matrix counterpart (its widths exceed what the channel radii allow).  The
  fit objective therefore carries a penalty :math:`\sum_\lambda (J_\lambda/10^{-3})^2`
  over levels with :math:`J_\lambda < 0`, MCMC rejects such points, and a run
  that ends there prints a warning naming the level;
* analytic derivatives (``--use-lm``, the analytic-gradient minimizer, the
  covariance band, pyazr's ``chi2_and_grad`` / ``residual_jacobian``) are
  available in both modes;
* the Wigner-limit bounds act on the Park amplitudes, i.e. on the observed
  :math:`\theta^2`.

Input and output files are the same in both modes: energies, partial widths and
ANCs in ``<levels>`` and ``parameters.out``.  ``param.par`` / ``param.sav`` hold
the mode's own amplitudes; their first line, ``#parametrization`` (0 standard,
1 Brune, 2 Park), says which (a comment line, so ``numpy.loadtxt`` and other
positional readers are unaffected), and a file written in the other alternative mode
is converted on read (``gamma_Park = gamma_Brune sqrt(J)``), with a message.  A
file without the tag (written before it existed) is taken as Brune's.  The GUI
offers the mode as "Use Park parametrization" under Runtime Options.

Level-matrix solver and profiling
---------------------------------

At every energy each :math:`J^\pi` group is solved in *channel space* when it
has fewer coupled channels than levels: the bilinear forms of the level matrix
come from a :math:`c \times c` system for :math:`c` channels instead of the
:math:`n \times n` level matrix of :math:`n` levels (Woodbury identity, see
:doc:`../developer/level_matrix`).  The other groups, and the energies too close
to a pole of the formal R-matrix for the channel-space form to keep full
precision, use the level matrix as before.  Results agree with the level-matrix
solution to about :math:`10^{-12}` relative; the analytic gradient is computed
in the same way.  Two environment variables control and report this:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Variable
     - Effect
   * - ``AZURE_LEVEL_SPACE=1``
     - Solve every group with the level matrix (the route of earlier versions),
       for comparison.
   * - ``AZURE_PROFILE=1``
     - After every evaluation and gradient, write to standard error the time
       spent (summed over threads) in the matrix fill and solve, the T-matrix,
       the channel solves, the cross-section sums, the Coulomb functions and
       penetrabilities, and the stages of the gradient, with the number of
       channel-space and level-matrix solves and why a group fell back to the
       level matrix.  No cost when unset.

``AZURE_CHANNEL_VERIFY=1`` (slow) compares every channel-space solve with the
level matrix and reports the largest difference seen so far;
``AZURE_CHANNEL_DEBUG=1`` reports the formal representation of each group
(and why a group has none).

Examples
--------

Run a calculation with Brune formalism:

.. code-block:: bash

   AZURE2 --no-gui --use-brune my_analysis.azr

Show all available options:

.. code-block:: bash

   AZURE2 --help
