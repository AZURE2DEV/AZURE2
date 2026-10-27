Configure Menu
==============

The **Configure** menu provides options for controlling the R-matrix formalism,
debug output, file directories, and runtime behavior.

Formalism
---------

The **Formalism** submenu allows choosing between two equivalent R-matrix
formalisms:

- **A-Matrix** (default) -- a level matrix formulation. Preferred when there are
  many channels and few levels.
- **R-Matrix** -- a channel matrix formulation. Preferred when there are many
  levels and few channels.

Both formalisms produce identical results; the choice is purely one of
computational efficiency.

Check Files
-----------

Selecting **Checks...** opens a dialog to control diagnostic output. AZURE2
can write detailed check files for debugging and verification. Each check file
can be set independently to one of three modes:

.. list-table::
   :widths: 20 80

   * - **None**
     - No output (default). Fastest execution.
   * - **Screen**
     - Print the check file contents to the terminal.
   * - **File**
     - Write to files with predefined names in the checks directory.

.. warning::

   Writing check files can drastically increase computation time. Keep them
   set to **None** unless debugging.

The check file directory is specified under **Directories** (see below).

Directories
-----------

Selecting **Directories...** opens a dialog where you can set:

- **Output Directory** -- where AZURE2 writes result files (cross sections,
  parameters, chi-squared, reaction rates, etc.).
- **Checks Directory** -- where debug check files are written.

Paths may be either absolute or relative to the Input file directory. These
directories **must exist** before running a calculation -- AZURE2 will not
create them automatically.

If no directories are specified, output files are written to the same directory
as the Input file.

Runtime Options
---------------

Selecting **Runtime Options...** opens a dialog with the following settings:

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - Option
     - Description
   * - **Use GSL Coulomb functions**
     - Use the GNU Scientific Library routines for Coulomb function calculations.
       Faster but potentially less accurate than the default method of N. Michel,
       *Computer Physics Communications* **176**, 232 (2007).
   * - **Use Brune formalism**
     - Use the Brune parameterization for R-matrix parameters
       (C.R. Brune, *Physical Review C* **66**, 044611, 2002). Recommended and
       enabled by default. More numerically stable than the classical approach.
       With it off, the physical parameters are transformed to formal ones per
       :math:`J^\pi` group; since 2026-10-27 that transformation keeps each
       level's identity whatever the order the levels are listed in (see the
       note in :doc:`levels_channels`).
   * - **Ignore external width if internal width is zeroed**
     - When enabled, the external gamma-ray width of a level is set to zero if
       no total gamma-ray width is specified.
   * - **Use RMC capture formalism**
     - Enable the Reich-Moore capture formalism. Currently only supported for
       (n, gamma) reactions. Do not use for other reaction types.
   * - **Do not perform parameter transformations**
     - Input R-matrix formal widths and pole energies directly, without
       transforming from physical parameters. Useful when starting from an older
       calculation. Note that formal parameters are radius and boundary-condition
       dependent.
   * - **Use Wigner limits for parameter limits**
     - Automatically apply Wigner limit bounds on fit parameters.
   * - **Use Hybrid Coulomb method**
     - Use a hybrid method for Coulomb function calculations that includes a
       nuclear potential (Woods-Saxon or Gaussian). This is the master switch;
       the Nuclear Potential tab then chooses which particle pairs it applies
       to, one at a time. Enables the Nuclear Potential
       tab.

.. note::

   The **Brune formalism** and **RMC capture formalism** are mutually exclusive.
   Enabling one will automatically disable the other.

.. note::

   When running AZURE2 from the command line (``--no-gui``), these runtime
   options are not applied from the saved configuration. They must be specified
   as command-line flags each time. See :doc:`/reference/command_line`.

Coulomb functions from a table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For charged-particle channels, open and closed, AZURE2 takes the
penetrability, the shift function and the hard-sphere phase from a
precomputed table when one is available. The table covers orbital angular
momenta 0 to 10 and the Sommerfeld parameter η and :math:`\rho = ka` from
10\ :sup:`−3` to 100 each, with :math:`\eta\rho \le 150`. Outside that range, for neutral channels, for
pairs with a nuclear potential (**Use Hybrid Coulomb method**) and for
external capture, AZURE2 computes the functions directly. The table agrees
with the direct calculation to better than 10\ :sup:`−7`, and a calculation
gives the same χ² to every printed digit with or without it. It is much
faster: the natural-neon example (:doc:`../examples/natural_neon`) takes about
1.5 minutes on four threads with the table and about 6 minutes without it.

The dialog has no switch for the table: the GUI uses it whenever it finds
the file. The file is ``~/.azure2/coulomb_table.bin``, or the file named by
the environment variable ``AZURE2_COULOMB_TABLE``. It is not part of the
download. Build it once, either with **Configure → Build Coulomb Table...**,
which shows the progress and can be cancelled, or from the command line:

.. code-block:: console

   $ AZURE2 --build-coulomb-table

This takes about 5 minutes on four threads and writes 88 MB.
``--build-coulomb-table=FILE`` writes another file. Calculations started after
the build use the table; there is no need to restart AZURE2. Without the
file, a run prints one line saying so and computes the functions directly;
the results are the same, only slower. A run that used the table reports, at
the end, how many values it took from it and how many fell outside it.

From the command line the table can be switched off with
``--no-coulomb-table`` or another file named with ``--coulomb-table=FILE``
(pyazr: ``use_coulomb_table=False``, ``coulomb_table="FILE"``); see
:doc:`/reference/command_line`. Within the table's range it takes the place
of both the default routine and **Use GSL Coulomb functions**. The design and
the validation are in :doc:`../developer/coulomb_table`.
