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

Examples
--------

Run a calculation with Brune formalism:

.. code-block:: bash

   AZURE2 --no-gui --use-brune my_analysis.azr

Show all available options:

.. code-block:: bash

   AZURE2 --help
