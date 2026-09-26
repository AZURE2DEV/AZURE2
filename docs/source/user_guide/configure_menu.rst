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

THM Options
-----------

Selecting **THM Options...** edits the optional ``<thm>`` block: the choices in
the Trojan Horse (half-off-energy-shell) observable of THM segments
(observable code 10 or more). They have no effect on other segments.

- **Entrance partial waves** (``entranceL``) -- incoherent (default) or
  coherent sum over the entrance orbital momenta.
- **Vertex boundary** (``vertex``) -- constant (default), perlevel or onshell
  boundary term in the transfer vertex.
- **Kinematic factors** (``kinematics``) -- lacognata (default), triple,
  kf3body or lambda32, matching how the HOES data were extracted.
- **External Coulomb term in the vertex** (``coulombIntegral``).
- **Spectator energy** (``spectatorEnergy``, MeV) for every THM entrance pair,
  and a table of per-pair values (``spectatorEnergy[<pair>]``).
- **Energy-dependent weight per segment** (``weight[<segment>]``,
  ``weightTest[<segment>]``) -- a two-column table w(E) multiplying the THM
  model of one segment, e.g. a Coulomb-distortion correction. The segment
  number is the line number in the Data (or Test) segments table, counting
  inactive lines. The **...** button picks the file; a file inside the project
  directory is stored relative to it.

Each control has a tooltip with the physics in one line; the full description
is in :doc:`../theory/thm_implementation`. Options left at their defaults are
not written, and a project whose options are all default has no ``<thm>``
block. If the dialog is not changed, the block is saved exactly as it was read.
