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
numbers are entered in spin boxes that show their unit (typed as AZURE2
reads them, ``1e-05`` included, and kept to the last digit), and values AZURE2
derives are shown as short label--value pairs next to the fields they come
from (their full text, as AZURE2 prints it, in their tooltips).

If the project has no THM segment, the window says so and its pages are
disabled: tick *THM* on a segment in the Segments tab first.

**Model** -- the options of the ``<thm>`` block:

- **Entrance partial waves** (``entranceL``) -- incoherent (default) or
  coherent sum over the entrance orbital momenta; coherent cannot be combined
  with an exit-angle window (**Experiments** page).
- **Vertex boundary** (``vertex``) -- constant (default), perlevel or onshell
  boundary term in the transfer vertex.
- **Kinematic factors** (``kinematics``) -- lacognata (default), triple,
  kf3body or lambda32, matching how the HOES data were extracted.
- **External Coulomb term in the vertex** (``coulombIntegral``); refused on
  Accept with an experiment whose a + A wave is distorted (distortion coulomb
  or optical without a plane a + A wave, or the DW vertex), which already
  contains that Coulomb force (:doc:`../theory/thm_implementation`, "Coulomb
  effects: what each option contains").
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
in four sections -- *Experiment*, *Three-body reaction*, *Spectator
momentum window* and *Distortion*:

- **Name** -- letters, digits and ``_ - . +``.
- **Segments** -- a check list of the Data segments that can be added: THM
  segments with a free norm that are in no other experiment (the segments of an
  experiment share one norm, so a fixed one is refused).
- **Background** (``background``) -- none, const, linear or quadratic.
- **Coherent background** (``cbackground``) -- a checkable section below
  *Experiment*; unchecked (the default) writes nothing. Checked, it shows a
  table with one row per term: **J**\ :sup:`π` (offered from the levels),
  **Exit** (the exit pair key of a segment of the experiment), **Channels**
  (*all*: every entrance :math:`(s,l)` and exit :math:`(s',l')` of the group,
  each with its own amplitude; or one, ``s,l,s',l'`` such as ``1/2,0,1/2,1``),
  **Form** (*const* or *linear*) and the start values **Re c0**, **Im c0**
  (and **Re c1**, **Im c1** for *linear*); a ticked value is fixed in the
  fit. **+** and **−** add and remove a term. The complex amplitude
  :math:`c_0 + c_1 E` times the entrance vertex is added to the resonant HOES
  amplitude before squaring: a THM-only background that interferes with the
  resonances (for non-resonant physics that is also in direct data, add a
  broad background level of the same J\ :sup:`π` on the Levels tab instead).
  Its free parameters appear on the **Fitting** tab (sub-tab *THM
  Background*), and a ``.sav`` loaded there writes their fitted values into
  the key (one explicit term per channel combination). The page refuses, in
  AZURE2's words, a malformed term, an exit pair of no segment of the
  experiment, a J\ :sup:`π` the levels do not have, and ``entranceL``
  *coherent*. Details: :doc:`../theory/thm_implementation`, "Coherent
  background".
- **Exit angle** (``theta``) -- *all (angle-integrated)*, the default
  (nothing written), or *window*, with :math:`\theta_\mathrm{min}` and
  :math:`\theta_\mathrm{max}` in degrees (0-180) beside it. With a window
  the model of every segment of the experiment is
  :math:`d\sigma/d\Omega` averaged over the window, where at a fixed angle
  the entrance partial waves and the J\ :sup:`π` groups interfere. The angle
  is the c.m. angle of particle 1 of the exit pair relative to particle 2,
  measured from the x-A direction (entrance particle 1 relative to 2) -- the
  :math:`\theta_\mathrm{cm}` of the Catania analyses when x is entrance
  particle 1; if the paper quotes the other exit particle's angle, give the
  supplementary window :math:`180-\theta_\mathrm{max}` to
  :math:`180-\theta_\mathrm{min}`. Equal ends are one angle (``theta=0-0``).
  The page writes ``theta=<min>-<max>``, the numbers as read or typed. A
  window cannot be combined with *coherent* **Entrance partial waves**
  (Model page): AZURE2 refuses it, and so does the page, in its words, as
  soon as the Experiments page is shown again; *window* is then not offered
  for an experiment that has none. AZURE2 also refuses a model with an exit
  orbital momentum above 40. Details: :doc:`../theory/thm_implementation`,
  "Fixed-angle observable".
- **Three-body reaction** -- beam, target and spectator (a nuclide of the
  built-in table, or ``Z,A,mass`` with the nuclear mass in u) and the lab beam
  energy in MeV. The four are a unit: untick the group to write none of them.
  When they are complete, the section shows what AZURE2 prints for them: the
  binding energy **B(x+s)** of the Trojan horse and the quasi-free energy
  **E_qf** (E(x+A) in its tooltip). The fit uses them only for the
  spectator's Coulomb line shape, momentum window and distortion factor,
  below. B(x+s) comes from the masses, while the THM vertex takes B from the
  entrance pair (field 32, **Channels** page); when the two differ by more
  than 1 keV, a warning row in this section says so, e.g. *B from masses
  5.493 MeV ≠ pair B 2.225 MeV (vertex uses the pair value)* (its tooltip:
  AZURE2's startup ``WARNING``), and the **Channels** page marks that pair's
  B with a warning icon. This is not a refusal: a made-up Trojan horse can
  stand in for the real one, but the kinematics then belong to it.
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

- **Spectator acceptance** (``ps``, ``psNodes``; ``spectatorAngles``,
  ``spectatorAngleNodes``) -- its own section,
  offered once the four fields of the three-body reaction are filled. The off-shell x-A
  momentum of the entrance vertex depends on the spectator momentum
  :math:`p_s` (:math:`p_{xA}^2/2\mu_{xA} = E + B + p_s^2/2\mu_{sx}`), and THM
  data are averaged over the accepted :math:`p_s` window. At fixed :math:`E`
  the spectator direction fixes :math:`p_s`, and AZURE2 averages the HOES
  cross section over the accepted directions with the fixed-:math:`E` event
  weight :math:`|\phi(p_s)|^2\,d\cos\theta_\mathrm{cm}`, i.e.
  :math:`|\phi|^2 p_s\,dp_s` over the :math:`p_s` the kinematics reach. This
  matters most near the nodes of the vertex. **Distribution**: *point* (the quasi-free :math:`p_s = 0`, the
  default; nothing is written), *Hulthén*, *Gaussian* or *table*. Only the
  fields of the chosen distribution are shown: **p_min** and **p_max** in
  MeV/c (Hulthén and Gaussian); the Hulthén **a** and **b** in
  fm\ :sup:`-1`, shown with the deuteron's standard values 0.2317 and 1.202
  and editable once **custom a, b** (next to the distribution) is ticked (e.g. an Eckart function for
  :sup:`3`\ He or :sup:`6`\ Li); the Gaussian **FWHM** of
  :math:`|\phi|^2` in MeV/c; the **table** file (two columns, :math:`p_s` and
  the momentum distribution :math:`|\phi(p_s)|^2`; its range is the window), whose
  **...** button stores a file inside the project directory relative to it,
  as for weight tables. **Nodes** is the number of Gauss-Legendre
  nodes in :math:`\cos\theta_\mathrm{cm}` (1-64, default 16; ``psNodes`` is
  written only when it is not 16; disabled with a direction window, which has
  its own).
  The page writes ``ps=hulthen:pmin-pmax``, ``ps=hulthen:a,b:pmin-pmax``,
  ``ps=gauss:FWHM:pmin-pmax`` or ``ps=table:<file>``, keeping numbers as
  they were read or typed. With a window set, the section shows the mean
  spectator energy :math:`\langle T_s\rangle` at the lowest and the highest
  point (the accepted directions, and so the nodes, follow :math:`E`); its
  tooltip has the rest of what AZURE2 prints at startup (the window,
  :math:`\mu_{sx}`), computed with AZURE2's own code. A window drops the
  Distortion section's single **Angle** (the window sets the directions). A window cannot be combined with a
  non-zero **Spectator energy** (Model page) for the same entrance pair -- the
  window replaces it; AZURE2 refuses the pair, and so does the page. Unticking
  the three-body reaction drops the window with it. Details:
  :doc:`../theory/thm_implementation`, "Spectator-momentum window".

  With a momentum window or a Coulomb or optical **Distortion** (below) the
  section also shows **Directions**: the spectator directions over which the
  vertex (with a momentum distribution), the distortion factor R(E), or the
  distorted-wave vertex, is averaged. *all inside the p_s window* / *one*
  (the default, nothing written: every direction inside the momentum window,
  or without one the single **Angle** of the Distortion section), *lab
  window* or *c.m. window* -- the polar angle of the spectator
  to the beam, **θ_min** and **θ_max** in degrees -- or *lab table* /
  *c.m. table*, a file of two columns, the angle and the acceptance (its
  range is the window; **...** stores it relative to the project directory);
  **Nodes** is the number of Gauss-Legendre nodes in cos θ_cm per c.m.
  interval (1-64, default 8; ``spectatorAngleNodes`` is written only when it
  is not 8). The page writes ``spectatorAngles=7-30``,
  ``spectatorAngles=cm:135-180`` or ``spectatorAngles=[cm:]table:<file>``;
  a window drops ``spectatorAngle`` and switches the Distortion section's
  **Angle** off (AZURE2 refuses both together), and going back to *one*, to no
  computed distortion or to no reaction drops it. At fixed E the direction
  fixes the spectator momentum, so a momentum window then only cuts it. The
  **R(E)** row shows the window's average. Details:
  :doc:`../theory/thm_implementation`, "Experimental acceptance".

- **Distortion** (``distortion``, ``opticalAA``, ``opticalSF``,
  ``spectatorAngle``, ``distortionRef``, ``distortionRatio``,
  ``boundState``) -- its own section. The PWA data reduction takes the
  transfer amplitude as constant; with the distortions of the a + A and
  s + F relative motion it varies with E, by orders of magnitude for a
  charged spectator below the s + F barrier (12C+12C). AZURE2 multiplies the
  model of every segment of the experiment, before the folding, by
  :math:`R(E) = \rho(E)/\rho(E_\mathrm{ref})` from a zero-range DWBA
  (Mukhamedzhanov & Pang, PRC 99 (2019) 064618). **Distortion**: *None* (the
  default; only this combo is shown), *Coulomb* (point-Coulomb waves in both
  channels), *Optical* (per channel) or *Table*. Coulomb and optical need the
  three-body reaction: they are offered once its four fields are filled, and
  unticking the reaction drops them. With Coulomb or optical the section
  shows **Angle** (``spectatorAngle``: *quasi-free*, the default; *lab*, an
  angle to the beam converted at every E; *c.m.*, fixed; degrees in the spin
  box beside it), **E_ref** (``distortionRef``, MeV, where R = 1; *auto* is
  the middle of the data -- only the scale, which the profiled norm absorbs),
  **Ratio** (``distortionRatio``: *DWBA/PWBA*, :math:`\rho = |M/M_{PW}|^2`,
  the default, or *DWBA*, :math:`\rho = |M|^2`, the papers' ratio) and
  **Bound state** (``boundState``: *Whittaker*, the default, or *Yukawa*
  tail of the s-x bound state, with an optional cut-off **r_min** in fm).
  Optical adds a row for the two channels, **a + A** (``opticalAA``) and
  **s + F** (``opticalSF``): *plane* (no distortion), *Coulomb* (the default),
  *global* or *Woods–Saxon*, whose **Edit…** button opens the ten parameters in a
  compact form -- real volume V, R, a; imaginary volume W, R_W, a_W;
  imaginary surface W_D, R_D, a_D (MeV and fm, depths > 0 attractive or
  absorptive, 0 switches a term off) and the Coulomb radius R_C (0: a point
  charge). *global* shows a second combo with the built-in global optical
  potentials (``ancai06``, ``daehnick80`` for deuterons, ``kd03`` for
  nucleons, ``bg71`` and ``liang09`` for t/3He, ``mcfadden66`` and
  ``avrigeanu94`` for alphas; the first one that describes the channel's
  light partner is chosen); its **Edit…** shows the reference, the validity
  range, the ten numbers the potential gives at the lowest and highest data
  energy (the s + F potential follows E_sF) and **Use outside the validity
  range** (``:extrapolate``), and the button's tooltip the depths at the two
  ends. A potential outside its range is refused, on the page as by AZURE2,
  unless that box is ticked (AZURE2 then warns). Details:
  :doc:`../theory/thm_implementation`, "Global optical potentials". *Table* shows the file (two columns, E and w, the ``weight[k]``
  format; every data point inside it), whose **...** button stores a file
  inside the project directory relative to it; a table needs no reaction.
  The **R(E)** row shows R at the lowest and highest point energy of the
  experiment's data (a table: its w there), computed with AZURE2's own code;
  its tooltip has the rest of what AZURE2 prints at startup
  (:math:`k_{aA}`, :math:`\eta_{aA}`, :math:`\kappa`, :math:`\eta_b`,
  :math:`\beta`, E_ref, and :math:`E_{sF}`, :math:`\eta_{sF}` and the
  spectator's c.m. angle at the two ends). Every control changes only its own
  key; the others stay as written. The page refuses what AZURE2 refuses, in
  its words: a malformed value, a key without its kind, a data point where
  the spectator has no energy left or a lab angle it cannot reach, a table
  that cannot be read or does not cover the data. Details:
  :doc:`../theory/thm_implementation`, "Distortion factor R(E)".

Keys of an experiment line that the page does not show are kept as written.

**Channels** -- the THM columns of the ``<levels>`` lines:

- **Binding energy** B (MeV, field 32) of the transferred particle in the
  Trojan horse, for each particle pair that is the entrance of a THM segment.
  It is a property of the pair and is written on every line of the pair.
  A warning icon inside the field marks a pair whose B differs from B(x+s)
  of an experiment's reaction (from the masses, **Experiments** page) by
  more than 1 keV; its tooltip names the experiment and both values.
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
- **Spectator-momentum window** (only for an experiment with ``ps=`` and the
  plane-wave vertex) -- the event weight per unit :math:`p_s` at the middle
  energy of the segment, :math:`A\,|\phi(p_s)|^2 p_s` (the fixed-:math:`E`
  measure; :math:`A` the acceptance of a direction window) over the accepted
  :math:`p_s` (inside :math:`[p_\mathrm{min}, p_\mathrm{max}]` and reachable
  at that energy), normalized to unit area, with dots at the nodes where
  AZURE2 evaluates the vertex there (they follow the energy).
  The **Entrance vertex** panel then also shows, dashed, the window average
  :math:`\langle |M_l|^2\rangle` that AZURE2 uses, next to the quasi-free
  :math:`|M_l|^2`: the window fills the vertex nodes. The status line's
  tooltip gives :math:`\langle T_s\rangle`. A table named relative to the project is found
  although the engine runs on a copy elsewhere.
- **Distortion R(E)** (only for an experiment with ``distortion=``) -- the
  factor that multiplies the model, as AZURE2 interpolates it, on a
  logarithmic scale when it spans decades, with :math:`|M|^2` (dashed; not
  for the DWBA ratio, where it is R itself) and :math:`|M_{PW}|^2` (dotted),
  each 1 at :math:`E_\mathrm{ref}`, which a vertical line marks: R is their
  ratio. A table: its :math:`w(E)`, found although the engine runs on a copy
  elsewhere. The status line's tooltip gives the settings and
  :math:`E_\mathrm{ref}`.
- **Angular distribution** (only for an experiment with a ``theta``
  window) -- :math:`d\sigma/d\Omega(\theta)` of the HOES observable over
  0-180° at one energy, the experiment's window shaded, and the window
  average -- the model of a point -- dashed across it (a dot for a single
  angle). The energy is the spin box **E** in the card's header, within the
  data; the first **Compute** takes the middle of the data, and a changed
  energy is used by the next **Compute**. AZURE2 has no call for the
  distribution at arbitrary angles, so the curve comes from a copy of the
  project run at that energy with 19 single-angle windows (every 10°) and
  the segment's experiment otherwise (reaction, line shape, momentum window);
  without resolution, weight, distortion or background, which at one energy
  only scale it. :math:`4\pi` times its average over 0-180° is the
  angle-integrated HOES cross section.

The page needs a build with the engine API (``USE_API=ON``, the default),
as the Plot tab needs Qwt.

Each control has a tooltip with the physics in one line; the full description
is in :doc:`../theory/thm_implementation`. Options left at their defaults are
not written, and a project whose options are all default and that has no
experiment has no ``<thm>`` block. The workspace writes nothing else: values
left as they were read are saved exactly as they were read, comments included,
and an edited experiment is written as one line in place of its first.
