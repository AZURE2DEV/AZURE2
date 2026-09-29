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
are edited here only. The window has three pages that edit and a fourth,
**Diagnostics**, that only shows. **Accept** checks the pages with the rules
AZURE2 applies at startup (a refusal shows AZURE2's own message and the page
it concerns) before anything is changed. **Cancel** changes nothing.
**Help** opens this section of the user guide. The pages carry no
explanatory text of their own: every field explains itself in its tooltip,
numbers are entered in spin boxes that show their unit, and values AZURE2
derives are shown as short label--value pairs next to the fields they come
from (their full text, as AZURE2 prints it, in their tooltips).

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
background. The list on the left shows the experiments (name and segments;
the background and the reaction in the tooltip of a row), with **Add** and
**Remove** below it. The editor on the right shows the selected experiment
in three sections -- *Experiment*, *Three-body reaction* and *Spectator
momentum window*:

- **Name** -- letters, digits and ``_ - . +``.
- **Segments** -- a check list of the Data segments that can be added: THM
  segments with a free norm that are in no other experiment (the segments of an
  experiment share one norm, so a fixed one is refused).
- **Background** (``background``) -- none, const, linear or quadratic.
- **Three-body reaction** -- beam, target and spectator (a nuclide of the
  built-in table, or ``Z,A,mass`` with the nuclear mass in u) and the lab beam
  energy in MeV. The four are a unit: untick the group to write none of them.
  When they are complete, the section shows what AZURE2 prints for them: the
  binding energy **B(x+s)** of the Trojan horse and the quasi-free energy
  **E_qf** (E(x+A) in its tooltip). The fit uses them only for the
  spectator's Coulomb line shape and momentum window, below.
- **Line shape: Coulomb** (``lineshape=on``) -- the last row of the
  three-body reaction section, offered once its four fields are filled. The
  charged spectator leaves in the Coulomb field of the resonance and of its
  decay products; this skews each resonance and shifts its peak (upwards for
  the usual sign of :math:`\zeta`). AZURE2 multiplies each level's exit
  amplitude by :math:`N_C` (Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020)
  233; Mukhamedzhanov, EPJA 58 (2022) 71). It needs the Brune
  parameterization, and the spectator must keep some energy at the highest
  data point; the page refuses otherwise, with AZURE2's words. When on, the
  row shows :math:`\zeta` of every exit pair at the lowest and highest point
  energy of the experiment's data -- the values AZURE2 writes to
  ``output/thm_experiments.out``. A neutral spectator gives
  :math:`\zeta = 0`, i.e. no change. Details:
  :doc:`../theory/thm_implementation`, "Coulomb line shape".

- **Spectator momentum window** (``ps``, ``psNodes``) -- its own section,
  offered once the four fields of the three-body reaction are filled. The off-shell x-A
  momentum of the entrance vertex depends on the spectator momentum
  :math:`p_s` (:math:`p_{xA}^2/2\mu_{xA} = E + B + p_s^2/2\mu_{sx}`), and THM
  data are averaged over the accepted :math:`p_s` window; AZURE2 then averages
  the HOES cross section over the window with the weight
  :math:`w(p) = |\phi(p)|^2 p^2`. This matters most near the nodes of the
  vertex. **Distribution**: *point* (the quasi-free :math:`p_s = 0`, the
  default; nothing is written), *Hulthén*, *Gaussian* or *table*. Only the
  fields of the chosen distribution are shown: **p_min** and **p_max** in
  MeV/c (Hulthén and Gaussian); the Hulthén **a** and **b** in
  fm\ :sup:`-1`, shown with the deuteron's standard values 0.2317 and 1.202
  and editable once **custom a, b** (next to the distribution) is ticked (e.g. an Eckart function for
  :sup:`3`\ He or :sup:`6`\ Li); the Gaussian **FWHM** of
  :math:`|\phi|^2` in MeV/c; the **table** file (two columns, :math:`p_s` and
  the event weight per unit :math:`p_s`; its range is the window), whose
  **...** button stores a file inside the project directory relative to it,
  as for weight tables. **Nodes** is the number of Gauss-Legendre
  nodes (1-64, default 16; ``psNodes`` is written only when it is not 16).
  The page writes ``ps=hulthen:pmin-pmax``, ``ps=hulthen:a,b:pmin-pmax``,
  ``ps=gauss:FWHM:pmin-pmax`` or ``ps=table:<file>``, keeping numbers as
  they were read or typed. With a window set, the section shows the mean
  spectator energy :math:`\langle T_s\rangle`; its tooltip has the rest of
  what AZURE2 prints at startup (the window, :math:`\mu_{sx}` and the range
  of :math:`T_s = p_s^2/2\mu_{sx}` over the nodes), computed with AZURE2's
  own code. A window cannot be combined with a
  non-zero **Spectator energy** (Model page) for the same entrance pair -- the
  window replaces it; AZURE2 refuses the pair, and so does the page. Unticking
  the three-body reaction drops the window with it. Details:
  :doc:`../theory/thm_implementation`, "Spectator-momentum window".

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

**Diagnostics** -- read-only plots of one THM data segment, computed by AZURE2
when **Compute** is pressed (never on an edit), for the parameters in the
Levels tab and the settings of the other pages as they stand -- including
changes not yet accepted. The engine runs on a temporary copy of the project,
in the background, while a busy bar runs; nothing is written into the project
or its output directory. One toolbar row holds the controls: the segment (its
experiment, if any, named in front), the J\ :sup:`π` group of the entrance
vertex panel and **Compute**. Below it, one status line gives the segment,
the number of points and the energy range; its tooltip has the details
(vertex option, B + T\ :sub:`s`, the spectator window's
:math:`\langle T_s\rangle`, line-shape levels not drawn). If the project or
the workspace changes after a computation, the status line says so; press
**Compute** again. The panels are framed cards of equal size in a grid (three
columns when the window is wide enough, else two) whose axes line up; the
energy axis is the data range of the segment (c.m. of the THM entrance pair),
and each plot carries the reaction in bold:

- **Entrance vertex** -- :math:`|M_l(E)|^2` at the quasi-free point, one curve
  per entrance orbital momentum of the J\ :sup:`π` group chosen in the
  toolbar, with the boundary the ``vertex`` option gives; dashed lines mark its
  nodes. A resonance close to a node is suppressed in the HOES cross section.
  With ``vertex=onshell`` or the Coulomb term the vertex is complex and has no
  nodes.
- **HOES and on-shell** -- the HOES cross section of the segment's channel
  (without resolution, weight or line shape; its arbitrary scale matched to
  the other curve) and the ordinary angle-integrated cross section of the same
  channel, on a logarithmic scale. Where the two differ in shape, the vertex
  and the kinematic factors are at work.
- **Line shape** (only for an experiment with the line shape on) --
  :math:`|N_C|^2` of the levels with a pole near the data (1 at the pole), and
  :math:`\zeta(E)` of the segment's exit pair.
- **Weight** (only when the segment has a weight table) -- :math:`w(E)` as the
  engine interpolates the table.
- **Spectator-momentum window** (only for an experiment with ``ps=``) --
  the event weight :math:`w(p_s) = |\phi(p_s)|^2 p_s^2` (a table: as given)
  over :math:`[p_\mathrm{min}, p_\mathrm{max}]`, normalized to unit area,
  with dots at the Gauss-Legendre nodes where AZURE2 evaluates the vertex.
  The **Entrance vertex** panel then also shows, dashed, the window average
  :math:`\langle |M_l|^2\rangle` that AZURE2 uses, next to the quasi-free
  :math:`|M_l|^2`: the window fills the vertex nodes. The status line's
  tooltip gives :math:`\langle T_s\rangle`. A table named relative to the project is found
  although the engine runs on a copy elsewhere.

The page needs a build with the engine API (``USE_API=ON``, the default),
as the Plot tab needs Qwt.

Each control has a tooltip with the physics in one line; the full description
is in :doc:`../theory/thm_implementation`. Options left at their defaults are
not written, and a project whose options are all default and that has no
experiment has no ``<thm>`` block. The workspace writes nothing else: values
left as they were read are saved exactly as they were read, comments included,
and an edited experiment is written as one line in place of its first.
