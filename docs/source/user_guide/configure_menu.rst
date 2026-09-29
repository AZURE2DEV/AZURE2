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

THM Workspace
-------------

Selecting **THM Workspace...** opens one window for everything that belongs to
the Trojan Horse (half-off-energy-shell) observable of THM segments
(observable code 10 or more): the optional ``<thm>`` block and the two THM
columns of the ``<levels>`` lines. The classic tabs do not show these; they
are edited here only. The window has three pages, and **Accept** checks them
with the rules AZURE2 applies at startup (a refusal shows AZURE2's own message
and the page it concerns) before anything is changed. **Cancel** changes
nothing.

If the project has no THM segment, the window says so and its pages are
disabled: tick *THM* on a segment in the Segments tab first.

**Model** -- the options of the ``<thm>`` block:

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

**Experiments** -- the ``experiment[<name>]`` lines: THM data segments
measured together (exit channels, angular bins or runs of one three-body
reaction) that share one profiled normalization and, optionally, a smooth
background. The table lists the experiments (name, segments, background,
reaction); **Add** and **Remove** are below it. The editor underneath shows the
selected experiment:

- **Name** -- letters, digits and ``_ - . +``.
- **Segments** -- a check list of the Data segments that can be added: THM
  segments with a free norm that are in no other experiment (the segments of an
  experiment share one norm, so a fixed one is refused).
- **Background** (``background``) -- none, const, linear or quadratic.
- **Three-body reaction** -- beam, target and spectator (a nuclide of the
  built-in table, or ``Z,A,mass`` with the nuclear mass in u) and the lab beam
  energy in MeV. The four are a unit: untick the group to write none of them.
  When they are complete, the page shows what AZURE2 prints for them: the
  binding energy B(x+s) of the Trojan horse and the quasi-free energy. The fit
  does not use them yet.

Keys of an experiment line that the page does not show are kept as written.

**Channels** -- the THM columns of the ``<levels>`` lines:

- **Binding energy** B (MeV, field 32) of the transferred particle in the
  Trojan horse, for each particle pair that is the entrance of a THM segment.
  It is a property of the pair and is written on every line of the pair.
- **Width input convention** (field 33) of each particle channel of those
  pairs: ticked, the width entered in the Levels tab is a reduced width
  amplitude (MeV\ :sup:`1/2`) rather than a partial width or ANC -- for a
  level known only from THM data, say. The value is not converted when the
  flag changes; re-enter it in the Levels tab, where it is labelled with its
  unit.

A pair with a binding energy, or a channel with the flag, is listed even if no
THM segment uses its pair, so that nothing in the file is out of reach.

Each control has a tooltip with the physics in one line; the full description
is in :doc:`../theory/thm_implementation`. Options left at their defaults are
not written, and a project whose options are all default and that has no
experiment has no ``<thm>`` block. The workspace writes nothing else: values
left as they were read are saved exactly as they were read, comments included,
and an edited experiment is written as one line in place of its first.
