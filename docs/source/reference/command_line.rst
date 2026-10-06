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
       "Use Brune formalism" option).
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
  (an ANC for a closed channel), with no factor depending on the other channels;
* every real Brune amplitude is a valid model, but Park amplitudes must satisfy
  :math:`J_\lambda > 0`.  Nothing in the level matrix enforces that, and a
  level with :math:`J_\lambda \le 0` has no standard R-matrix counterpart;
* ``--use-park`` has no analytic derivatives (the fit falls back to numerical
  ones), and the Wigner-limit bounds act on the Park amplitudes.

Input and output files are the same in both modes: energies, partial widths and
ANCs in ``<levels>`` and ``parameters.out``.  Only ``param.par`` / ``param.sav``
hold the mode's own amplitudes, so a ``param.sav`` written in one mode must not
be read in the other.

Examples
--------

Run a calculation with Brune formalism:

.. code-block:: bash

   AZURE2 --no-gui --use-brune my_analysis.azr

Show all available options:

.. code-block:: bash

   AZURE2 --help
