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

The answers to the prompts can be piped in, one per line, e.g.
``printf '1\n\n\n7\n' | AZURE2 --no-gui --no-readline project.azr`` (calculate
with data, parameters and external-capture integrals from the project).  If
the input ends at a prompt that needs an answer (the menu, the MINOS variance,
the reaction-rate and MCMC questions without a default), AZURE2 stops with
``ERROR: the input ended (end of file) at the prompt ...`` and exit status 1.
A prompt whose blank answer is its default (the parameter and integral file
names, the uncertainty-band questions, the MCMC spreads, overwriting
``samples.mcmc``) takes that default at end of input.

Exit status: 0 for a run that completed; 1 when the input ended at a prompt
that needs an answer; non-zero (-1, i.e. 255 on POSIX systems) when a run
fails -- a project AZURE2 refuses (``ERROR: ...``, e.g. a malformed ``<thm>``
block), a failed fill, initialization or parameter transformation, and an
MCMC run that is refused or fails. Scripts can therefore test the status
instead of parsing the log.

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
     - Use GSL Coulomb functions (equivalent to the GUI option).
   * - ``--no-transform``
     - Do not perform parameter transformations; input formal R-matrix
       parameters directly (equivalent to the GUI option).
   * - ``--no-long-wavelength``
     - Do not use the long-wavelength approximation for electric capture.
   * - ``--use-gradient``
     - Fit with MIGRAD and the analytic gradient (default: MIGRAD with
       numerical derivatives).
   * - ``--use-lm``
     - Fit with the Levenberg-Marquardt minimizer on the analytic Jacobian
       (falls back to MIGRAD).
   * - ``--use-gsl-lm``
     - Fit with GSL's trust-region least squares with geodesic acceleration
       (analytic Jacobian; falls back to MIGRAD).
   * - ``--covariance-band``
     - Compute the cross-section uncertainty band without asking (a fit saves
       its covariance to ``covariance.dat``; a calculation without data
       reuses it).
   * - ``--scale-covariance``
     - Scale the band's covariance to a reduced chi-squared of 1.

Flags can also be put in a file, one per line, named by the environment
variable ``AZURE_OPTIONS_FILE``; they are applied after those of the command
line. An unknown flag draws a ``WARNING`` and is ignored.

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
(``tests/park_formalism/check.sh``), THM (HOES) segments included: their
amplitude is bilinear in the reduced widths with level-diagonal factors only
(:doc:`../theory/thm_implementation`, "Brune and Park";
``tests/thm_park/check.sh``).  A channel entered as an amplitude
(``gammaIsRWA``, the 33rd field of a level line) holds Brune's amplitude in
either mode, so a file is the same model in both.  What changes is the
parameter space the minimizer works in:

* a Park amplitude is the observed width, :math:`\Gamma_c = 2P_c\gamma_c^2`
  (an ANC for a closed channel), with no factor depending on the other channels
  of the level.  Fixing, bounding or putting a prior on a width therefore acts
  on one parameter, and the parameters of a level are less correlated;
* every real Brune amplitude is a valid model, but Park amplitudes must satisfy
  :math:`J_\lambda > 0`; a level with :math:`J_\lambda \le 0` has no standard
  R-matrix counterpart (its widths exceed what the channel radii allow).  The
  fit objective therefore carries a penalty :math:`\sum_\lambda (J_\lambda/10^{-3})^2`
  over levels with :math:`J_\lambda < 0`, MCMC rejects such points (a chain
  started from such a point is reported, and refused if its whole starting
  ensemble is outside: :doc:`../user_guide/mcmc`), and a run that ends there
  prints a warning naming the level;
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

Examples
--------

Run a calculation with Brune formalism:

.. code-block:: bash

   AZURE2 --no-gui --use-brune my_analysis.azr

Show all available options:

.. code-block:: bash

   AZURE2 --help
