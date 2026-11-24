Graphical Interface
===================

AZURE2 has two graphical front ends to one engine: the Qt setup program
``AZURESetup``, with the tabs users of earlier versions know, and the browser
page served by ``pyazr`` (:doc:`web_interface`). Both edit the same project
file, ``.azr`` in project format 3 (:doc:`../reference/data_formats`), and run
the same calculations as the ``AZURE2`` command-line program; neither holds
physics of its own. A run from ``AZURESetup`` writes the same output files,
byte for byte, as ``AZURE2 <command>`` on the saved project.

Building and starting the Qt program
------------------------------------

The Qt program needs Qt 6.4 or newer (Core, Gui, Widgets, Network, Svg; CI
builds it with the Qt 6.4 of Ubuntu 24.04) and is built only
when asked::

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_GUI=ON \
         -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
   cmake --build build --parallel
   build/gui/src/AZURESetup [project.azr]

Qt 6 asks for CMake 3.22 or newer; with an older CMake add
``-DQT_FORCE_MIN_CMAKE_VERSION_FOR_USING_QT=3.20`` (Qt's own switch; the build
of this program needs nothing newer). Without ``-DBUILD_GUI=ON`` nothing of the
Qt program is built and every other target is unchanged. The Plot tab is drawn
by the program itself (no Qwt is needed; ``-DUSE_QWT`` is gone).

The window
----------

While no project is open the window shows a start page: ``New Project`` and
``Open Project``, the recent projects, the example projects found on this
computer (``examples/`` and ``tests/upgrade/T1``-``T3`` of the checkout the
program was built from, an installed ``share/azure2``, or the directory named
by ``AZURE2_EXAMPLES``), and the switch ``Show expert options``. A click on a
project opens it; an ``.azr`` file dropped on the window is opened too.
``File > Close Project`` returns to the start page. What each tab is for is its
tab's tool tip.

Opening a project never shows a box: what there is to know about it (the
engine's check: errors and warnings; what the converter assumed for an older
project; a ``<levels>`` the Levels tab could not read) is one bar under the tool
bar with the number of notes and the first one, ``Details`` (a window that does
not block the main one, with ``Show Diagnostics``) and ``Dismiss``. The
Diagnostics panel opens by itself only for errors (a project that would not
run) or after ``Validate``; warnings are counted in the status bar, whose click
opens the panel.

**Standard and expert view.** ``View > Show Expert Options`` (also on the start
page and as ``Standard view`` / ``Expert view`` in the status bar; off at first,
kept between sessions) shows the options an occasional user does not need. In
standard view these are hidden: the tabs Nuclear Potential, THM and MCMC, the
sub-tabs Model Data, Per-channel radii, Gamma Decays and Isospin, the Kind
column of a table of one kind and the columns empty in every row, the gamma
multipolarity limits of a project without a photon pair, Export LaTeX, the call
limit, the covariance refreeze passes, the prior covariance file and the text
pages of ``<parameterSettings>`` and ``<mcmc>``, the commands MINOS and MCMC of
the Calculate tab's list and its options other than the parameter file, the
band, the rate's pairs and temperatures, MINOS's chi-squared rise and the
sampler's walkers and steps, the Physics tab's pole map (and its THM page
without THM data), and the Tools menu. The switch changes what is shown, never
the project or what the engine computes, and an expert option the project uses
(a record, a section, a value set) is shown in either view, marked as such in
its tool tip.

**Tables and forms.** A tab of records (data sets, extrapolations, targets,
beams, detectors, ...) shows the table above and the selected record's form
below, split by a handle: the table keeps at least five rows, ``Form`` closes
and opens the form (a double click on a row opens the record in a window of its
own), and the sizes you set are kept per table between sessions. The ``?``
beside the buttons says what the section is.

The tool bar holds the main actions, each with its shortcut in its tool tip:
New (Ctrl+N), Open (Ctrl+O), Save (Ctrl+S), Validate (F7), Calculate (F5), Fit
(F6), Stop (Shift+F5) and Plot (F8); the ``Run`` menu holds the same. Calculate
and Fit save the project and run that command with the options set in the
Calculate tab. The status bar shows the project file (with ``(modified)`` after
an edit; the title bar marks it too), whether the project is valid (``Valid``,
or the number of errors and warnings; a click opens the Diagnostics), and while
a run is going its progress and a Stop button.

``View > Appearance`` chooses the look: ``Follow the System`` (light or dark as
the desktop is; with a Qt older than 6.5 as the desktop's window colour is),
``Light`` or ``Dark``. The choice, the window's size and layout and the recent
projects are kept between sessions. ``AZURE2_THEME=native`` in the environment
leaves the platform's own style and colours.

The Overview page
-----------------

The first tab of an open project, shown when a project opens, says at a glance
what the evaluation is and where its fit stands. It only reads: the project as
it stands and the files of its output directory, through the engine's facade;
nothing is calculated again.

- **What to look at next**: plain sentences made from the numbers below, each a
  link to the tab or record concerned: the project does not validate, there is
  no run yet, the outputs are older than the project file or than a data file,
  the fit ended with a status other than ``converged``, data sets with
  chi-squared/N above 3 (the link opens the largest in the Plot tab),
  normalizations more than 2 sigma from their prior, parameters at a limit or
  without an error, active data sets the run left out. The buttons Validate,
  Calculate, Fit and Open Plot do what the tool bar's do.
- **The project**: the file, the compound nucleus, the pairs with their
  thresholds (E_x of the compound nucleus, MeV), the levels per J^pi as bars,
  the data sets and extrapolations active, the targets, beams and detectors and
  how many of them the active data sets use, the free parameters by kind, the
  basis (Brune, formal (Lane-Thomas) or Park) and A- or R-matrix, external
  capture, the minimizer chosen in Fitting Settings, the options set away from
  their default, and the validity (a click opens the Diagnostics).
- **The last run**: fit or calculation (inferred from the order in which a run
  writes its files: a fit writes ``fitStatus.out`` after the fitter and
  ``param.par`` at its start), when, a red mark when the project file or a data
  file is newer than the outputs, the fit's status as ``fitStatus.out`` gives it,
  the total chi-squared of the data (``chiSquared.out``'s
  ``Total-Chi-Squared``), the number of points N, the free parameters p and
  chi-squared/(N - p); then a table per data set (name, reaction, N,
  chi-squared, chi-squared/N with a bar, the norm the run used and its pull
  when the norm has a prior), sortable by a click on a header, the worst rows
  in bold and those above 3 in red. A click on a row opens that data set in the
  Plot tab. Without a run: ``No run yet`` and the buttons Calculate and Fit.
- **Charts**: chi-squared/N per data set (a dashed line at 1) and, during and
  after a fit run in this window, the best chi-squared against the evaluation
  (the engine writes no file of a fit's progress, so an earlier session's fit
  has no curve).

The page reads its files again after a run, a save, an edit and when it is
shown; the free parameters are counted by the engine's set-up of the saved
project in the background (seconds for a project with external capture).

Changes on disk
---------------

The window follows its files when another program changes them (a text
editor, pyazr, a script, a run on the command line):

- the project file changed and the window has no unsaved edit: it is read again
  at once, the tab and the selected records stay, and the status bar says
  ``Project reloaded: changed on disk``;
- the project file changed and the window has unsaved edits: nothing is
  overwritten; a bar under the tool bar offers ``Reload (discard mine)``,
  ``Keep mine`` (the window's version is written at the next Save) and ``Show
  differences`` (a line-by-line comparison of the two texts);
- new output files of a run (``chiSquared.out``, ``fitStatus.out``,
  ``AZUREOut_*``, ``param.sav``, ``samples.mcmc``): the Overview and the Plot
  tab read them (the plot keeps its check boxes and its zoom);
- a data file changed: the Overview's red mark names the data set.

The window's own Save is recognised by its content and causes no reload. A
change is read once the files have been quiet for a moment, so a file being
written is never read half-way; while a run of the window itself is going,
its outputs are read when it ends. ``View > Follow Changes on Disk`` turns
this off (kept between sessions).

Opening and saving a project
----------------------------

``File > Open`` reads a format-3 project. An older project (format 1 or 2) is
refused with the converter's command, and the program offers to convert it: the
converter of ``AZURE2 --convert`` writes a format-3 copy beside it (by default
``<name>.f3.azr``) and keeps the old file; ``File > Convert an Old Project``
does the same without opening. A project opened and saved without an edit is
written back byte for byte; an edit changes only the lines of what was edited
(a field of a record, a row of ``<levels>``, an entry of ``<config>``), with the
engine's own writers.

The engine checks the project as it stands: ``Validate`` (F7, also when a
project is opened) fills the Diagnostics panel below the tabs with the
diagnostics of ``AZURE2 --validate`` (severity, code, line, record and field,
message). The panel is hidden while there is nothing to show and opens by
itself when a validation finds an error or a warning; ``View > Diagnostics``
(Ctrl+Shift+D) shows or hides it. A double click on a diagnostic goes to its
tab and record. Before New, Open, Close or Quit the program asks whether to
save unsaved changes.

The tabs
--------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Tab
     - What it edits
   * - Overview
     - Nothing: the project and its last run at a glance (above).
   * - Particle Pairs
     - The pairs of the ``<levels>`` rows (:doc:`particle_pairs`); below them
       the per-channel radii (records ``channel_radius`` of ``<channelRadii>``).
   * - Levels and Channels
     - The levels and their channels (:doc:`levels_channels`); below them the
       gamma decays written as a total width and branchings (``<gammaDecays>``,
       :doc:`gamma_decay`) and the isospin-labelled levels (``<isospin>``).
   * - Segments
     - Segments From Data (the ``dataset``, ``normalization`` and
       ``energy_scale`` records of ``<data>``), Segments Without Data (the
       ``extrapolation`` records), Model Data (``<modelData>``)
       (:doc:`segments`). Spin observables, line windows, priors and nuisance
       parameters are fields of the records.
   * - Experimental Effects
     - Targets, beams and detectors (``<targets>``, :doc:`experimental_effects`),
       one list each: layers, isotopes and abundances, stopping and straggling,
       Doppler broadening, screening, fitted target parameters; beam spread and
       beam profile; resolution kernels, attenuation, UDR.
   * - Nuclear Potential
     - Hybrid-potential channels (``<potentials>``).
   * - THM
     - The Trojan horse experiments (``<thm>``): a form of the global options
       and one of each ``experiment[<name>]`` line, with the keys of the
       engine's own table (beam, target, spectator, beam energy, momentum
       distribution, distortion, the vertex model with ``vertexModel=ff``, the
       full transfer form factor, and its well, volume term and nodes, which
       are enabled once the vertex model is given; the background); a value the
       engine's reader refuses is named under the form.  ``Text`` shows the
       whole section; both edit the same lines.
   * - Fitting Settings
     - The fitter (the engine's list of fitters; MIGRAD unless another is
       chosen), the call limit and the covariance refreeze passes; the prior
       covariance file; the parameter table (below); the parameter settings
       (``<parameterSettings>``) and the sampler's block (``<mcmc>``) as text.
   * - Calculate
     - The command (calculate, fit, extrapolate, MINOS, reaction rate, and the
       sampler where the build has it) with its help in one line, and its
       options, each with its help beside it: the uncertainty band (a check
       box) and its scaling (enabled with the band), the parameter and capture
       files (a path and a chooser), the rate's pairs (drop-downs of the
       project's pairs; the exit also as -1, the total radiative capture) and
       temperatures (spin boxes in GK), MINOS's Δχ², the sampler's walkers and
       steps. A number left at its lowest value reads ``required`` or
       ``default ...`` and is not given to the program. The run's lines (warnings
       in amber, errors in red), and at the end the summary (below).
   * - MCMC
     - Shown only when the engine has the ``mcmc`` command. The Gaussian priors
       of the sampler (the ``<parameters>`` rows of ``<mcmc>``: Load Parameters
       lists the free parameters), and the run's settings: walkers, steps,
       initial spread, level-energy spread, threads, the sampled basis (rwa or
       physical) and a fresh start (overwrite ``samples.mcmc``). Run starts
       ``AZURE2 mcmc`` with them; the progress bar follows the steps and the
       acceptance fraction is shown as the program prints it. The samples are
       in ``output/samples.mcmc``; the burn-in is dropped when they are read,
       not by the sampler. ``Results`` reads the chain back (``Read Results``):
       the burn-in to leave out (steps), per parameter the mean, standard
       deviation and the 16, 50 and 84 % quantiles over the kept steps and
       walkers, the acceptance fraction (the share of walker moves that changed
       the position), and a plot of a parameter's trace (a curve per walker),
       its histogram, or the walkers' mean log posterior against the step.
   * - Plot
     - The data sets and extrapolations of the last run (below).
   * - Physics
     - The quantities between the level parameters and the observables: the
       physics functions, the project's channel functions, R-matrix and pole
       map, the THM parts, the stopping powers (below).

Records: tables and forms
-------------------------

The record tabs show each section (or each kind of ``<targets>``) as a table
with real columns, as the Segments tab of the earlier program did: a ``Use``
check box (a data set's or an extrapolation's ``active``), the name, then the
fields the engine's table names as a record's main ones: for a data set its
reaction, quantity, unit, data file, normalization and a ``Fit?`` box beside
it, the energy and angle windows, and its target, beam and detector; for a
target its kind, active nucleus, areal density, stopping and straggling; for a
beam its spread and profile. The whole record is the row's tool tip. The
columns sort by a click on their header, the filter line keeps the rows that
contain its text, and several rows can be selected (Shift, Ctrl).
F2, or a click on a selected cell, edits it in place: a list for a word, a unit
or the name of another record; the number of a normalization or of any value
with a ``fit`` box. ``Add...`` (Ins) opens the record's pop-up, as the earlier
program's segment dialogs: the kind, the name and every field a record of the
kind is usually given (for a data set the reaction, quantity, unit, data file,
columns, energy and angle windows, normalization, target, beam and detector;
for an extrapolation the reaction, quantity, unit, energy and angle grids and
the target, beam and detector), filled with the first choices of the project.
Under the data file its first lines are shown, with the number of rows and
columns of numbers, and the columns are made as many as the file has (check
their roles). ``OK`` checks the record with the engine: a refusal is shown
under its field and the pop-up stays; ``OK`` again writes it as it is (a field
left empty is written as ``?``, named by the diagnostics until it is given).
``Add another`` writes it and keeps the pop-up open with the same choices and
the next name, for a series of data files. A double click on a row (or
``Edit...``) opens the record in the same pop-up; the form below the table
stays. ``Duplicate`` (Ctrl+D),
``Remove`` (Del), ``Up`` and ``Down`` (Ctrl+Up, Ctrl+Down; only the record's
lines move in the file), ``Use`` and ``Leave out`` act on the selected records,
also from the context menu.

Below the table is the form of the chosen record. Nothing is typed that can be
chosen; every list comes from the engine's table of project format 3 and from
the project itself:

* a word of a closed list, a unit (of its dimension, or of the data set's
  quantity), the name of a target, beam, detector or shared normalization
  (with the target's epochs), and the quantity (by its title) are lists;
* a reaction is two lists, the entrance pair and the exit pair of the same
  compound nucleus (from the Particle Pairs tab), and the total capture where
  the project has a photon pair; a list beside them makes it a sum of reactions
  or a ratio of two: one row per term, its weight (empty: 1), its two pairs and,
  if asked, a fixed angle with its unit and frame (``Add term`` and ``-``); a
  term whose pairs the project has not is kept as written;
* a pair (of a channel radius, a potential, a gamma decay's final state) is a
  list of the project's pairs, still editable;
* a value that a fit may vary (a normalization, an areal density, a screening
  potential, a kernel's width) is its number, its unit, a ``Fit`` box and
  ``Prior...``: a Gaussian prior's centre and width (``%`` of the centre) and
  the largest value a fold's grid must span;
* a structured value (stopping, straggling, Doppler broadening, beam spread,
  beam profile, energy resolution, UDR, quadrature, ...) is a list of its forms
  (for the stopping: formula, evaluation, table), each explained in a line,
  and one labelled field per parameter with its unit; a formula is its
  expression and a field per parameter;
* a grid of an extrapolation is its kind (range, list, file, adaptive), its
  span, unit and frame, and a step or a number of points;
* a data file's columns are one list of roles per column, with the unit of an
  energy or angle column;
* a file is a line with ``Browse...``, relative to the project.

The quantity of a data set, an extrapolation or a model datum decides which of
its fields apply (``J``, ``l`` and ``s`` of a phase shift, the spin operator of
an analyzing power, the expansion and orders of Legendre coefficients, the
level of a resonance strength, ...): only those are shown, a ``*`` marks the
ones the quantity requires, and a tool tip says when one applies only under a
condition (a multipolarity only for a photon exit). Choosing another quantity
shows and hides them at once; a value in a field that disappears is kept in the
form until ``Apply``, which then leaves it out of the record and says so under
the form ("Dropped J: quantity 'sigma' does not use it"). In the ``Add...``
dialog a change of the quantity makes the rows again: the fields the new
quantity requires, its unit and columns chosen anew, the other values kept.

An optional field shows the value that applies when it is left out, greyed
(``yes (default)``), or as the line's placeholder (``default: 0``); nothing is
written for it.

A field that belongs with another is in its box (a target's composition, phase
and compound correction with the stopping; its epochs with the areal density;
the profile's truncation and detailed balance with the profile); a box is
closed unless the record sets it. The text button beside each field (an icon
with ``Ab``; its tool tip says so) shows the field as the project file writes
it, with its grammar; a value the lists cannot show opens as text. Each
field is checked by the engine's reader as it is edited, and its message is
shown under it. ``Apply`` checks the record against the whole project
(``--validate``), shows each refusal under its field and asks before writing a
record the engine refuses; choosing another record applies the form when the
engine accepts it, else the form stays. A record opened and applied without an
edit is written back byte for byte; an edit rewrites only the lines of the
fields that changed. ``Revert`` shows the record as the project holds it.

The form of a target, beam or detector says which data sets and
extrapolations use it (``Used by``); ``Assign to...`` checks them in a list.
An empty section shows a line that says how to add a record, and no empty form.

**What it does to a point.** Over the form of a target with a Doppler
broadening, a beam with an energy spread or a profile, and a detector with an
energy resolution, a curve per kernel shows what the record does to one point
of a data set that uses it, and its full width at half maximum under it, in keV
in the lab and in the c.m. frame. The curves are the engine's own kernels (the
free-gas Doppler broadening, the Gaussian spread, the rectangular and
triangular resolutions, the skewed-Gaussian profile) built as a run builds
them, from the form as it stands: a temperature, a width or a spread typed is
drawn at once, before ``Apply``. ``for`` chooses the data set (its entrance
pair's masses: the Doppler width depends on them) and ``at E_lab`` the beam
energy, by default the middle of that data set's energies. The Doppler
broadening of a free gas at temperature T is a Gaussian in velocity: at a lab
energy E much above kT its width is about 2.355 (2 E kT m/M)^1/2, m and M the
masses of the projectile and the target (0.17 keV for protons on 12C at
1.2 MeV and 300 K). A record no data set names yet is drawn for the project's
first particle pair at the energy typed, which the line under the curves says.

The Levels and Channels tab
---------------------------

The levels on the left, the channels of the selected level in the middle, the
channel limits and the details of the selected channel on the right; the gamma
decays and isospin records below. Opening a project selects its first level, so
its channels are on screen at once; the selection stays with the level through
every edit, filter, sort and undo.

* **Editing in place.** A double click, F2 or typing on a cell edits it: the
  energy E_x (MeV), the J^pi (J in a spin box of halves, the parity a +/-
  button), and in the channel table the value. The value column says what the
  number is and its unit: Γ [eV] for an open particle channel or a gamma width,
  C [fm^-1/2] (an ANC) for a closed channel, γ [MeV^1/2] for a row written as a
  reduced-width amplitude, μ [μ_N] and Q [b] for a static moment, B for a beta
  channel. A value that cannot be taken (not a number; a J^pi that has no
  channel with the pairs and limits; a level of the same J^pi within 1 keV of
  another, which the engine would join) is refused with the reason, at the cell
  and under the levels; the model is not changed. Enter commits and moves down.
  A changed J^pi gives the level the channels of its new J^pi; the channels it
  keeps keep their values.
* **Adding.** The bar under the levels: J, parity, E and ``Add`` (or Enter in
  the energy). The Insert key puts the cursor in the energy with the selected
  level's J^pi. ``+`` opens the dialog as before. ``Duplicate`` (Ctrl+D) copies
  the selected levels 10 keV above, with their widths and fix flags.
* **Removing.** Del (or ``-``, or the context menu) removes the selected levels
  with their channels after one confirmation that names them.
* **Undo.** Ctrl+Z undoes, Ctrl+Shift+Z (or Ctrl+Y) redoes, every change of
  this tab and of the Particle Pairs tab since the project was opened. Undone
  back to the opened state, the project is saved as it was opened, byte for
  byte, unless it was saved in between (a row the save rewrote keeps its content
  but may differ in its blanks).
* **Finding a level.** The search line (``3/2-`` or the digits of an energy),
  the J^pi list with the number of levels of each, ``Only included`` and the
  count of levels shown. A click on a header sorts; J^pi sorts by J, parity and
  energy, so the levels of one J^pi stand together, and the J^pi of the
  selected level is in bold. Ctrl+Right goes to the channels, Ctrl+Left back.
* **The context menu** (right click): add a level of this J^pi, duplicate,
  edit in the dialog, remove, include or leave out, fix or free the energies or
  the widths of the selected levels, show only this J^pi, add a background pole
  of this J^pi (a level 10 MeV above the J^pi's highest, energy fixed, widths 0:
  a starting point to edit). ``Fix All Widths`` and ``Fix All Energies`` act on
  the selected levels when several are selected, else on all.
* **J^pi without a level.** Under the limits, the J^pi the pairs allow up to the
  maximum l that no active level carries are listed, with what to do: with
  ``Hard-sphere scattering in every J^pi`` on they scatter as hard spheres;
  with it off they are left out of the calculation, and the text says to switch
  the option on or to add a level (a background pole is the usual choice). A
  click on a J^pi fills the add bar.
* **Channel details.** For the selected channel: the pair's particles, the
  separation energy, the channel radius (the channel's own from
  ``<channelRadii>`` when a record names it, the pair's beside it), and at the
  level's energy E_c = E_x - S - e, the penetrability P (below threshold the
  factor μa/(ħ²W²) that turns an ANC into a reduced width), the shift function S
  and dS/dE, the Wigner limit γ_W² = 3ħ²/(2μa²), the single-particle width
  Γ_sp = 2Pγ_W², the observed reduced width γ² = Γ/(2P) (C²/(2P) for an ANC)
  and θ² = γ²/γ_W², which the channel table also shows. These are the
  engine's channel functions (the same as at a data point); γ² is the
  thin-level value, not Brune's γ̃², which also needs the level's other widths
  (``parameters.out`` of a run gives that).

The Particle Pairs tab
----------------------

The headers carry their units. The excitation energies, the separation energy
and the channel radius are edited in place (a double click or F2; a value that
cannot be taken is refused with the reason); a double click elsewhere opens the
pair's dialog. Removing a pair asks first, naming the channels that go with it
and the data sets and extrapolations whose reaction names it (the engine
refuses them once the pair is gone; records that name a later pair by its key
then point one pair off). ``Duplicate`` (Ctrl+D) opens the dialog on a copy of
the selected pair, to change into another pair (an excited state, another
radius). Ctrl+Z undoes a pair's edit as in the Levels tab.

The pair's dialog: the pair type is a list (two particles, a photon pair of the
compound nucleus in a final level, a beta decay), and what a type does not use
is hidden (a photon pair has no light particle, separation energy or channel
radius); a look at another type and back gives the particle pair's values back.
Each particle is chosen by its element (a list of the mass table's symbols, or
the symbol typed) and its mass number A: Z and the atomic mass (AME2020) are
filled, and J and the parity where they are known without a spin table (n, p,
d, t, 3He, alpha, and 0+ for an even-even nucleus); otherwise the dialog says
to give them. J is a spin box of halves and the parity a +/- button. Every
value stays editable. ``Fill masses and S`` fills both masses and the
separation energy from the table.

The channels of a level (when it is added, or its J^pi or the limits change)
are the engine's rule, the one ``pyazr`` uses, applied to the pairs as the tab
holds them before a save.

The parameter table
-------------------

``Read Parameters`` (Fitting Settings) saves the project and lists every
parameter of a fit as the engine sets it up, under the names of ``param.sav``:
its kind, its start value in the fit's units (a width as its reduced-width
amplitude, MeV^1/2), whether it is fixed, its limits and its Gaussian prior.
The limits and the prior are the parameter's row of ``<parameterSettings>``;
for a level width they are partial widths in eV, which the engine converts;
empty limits are none. They apply to free parameters. Whether a parameter is
fixed, and its start value, are edited where they are written: the Levels tab,
a record's ``fixed``/``fit``. ``Show Fitted Values`` puts the values and errors
of the last fit's ``param.sav`` beside the start values. ``Take the Fitted
Values as the New Start`` copies ``param.sav`` to ``<project>.start.sav`` beside
the project and gives it to the next run as ``--parameter-file``: the run starts
from it, as ``AZURE2 <command> --parameter-file=<project>.start.sav`` does; the
project's own values are not changed (``pyazr``'s ``Session.save_fit`` writes
fitted values into a project).

The Plot tab
------------

The Plot tab shows what the last run wrote: for each data set the data with
their error bars and the calculated curve, for each extrapolation its curve,
read from the output files (``AZUREOut_*``, :doc:`../reference/output_files`)
and matched to the records by the engine. Several can be selected; each gets
its colour and marker in the legend, and a set of one point is a marker. The
x axis is the c.m. energy, the excitation energy or the c.m. angle; the y axis
the calculated column (the cross section in b, dσ/dΩ in b/sr, or the data set's
observable such as A_y or a phase shift) or the S factor (MeV b, for a cross
section). Either axis can be logarithmic (points that are not positive are then
left out). Where the run wrote the uncertainty band (``--covariance-band``) the
1σ band is shaded. The lower panel shows the residuals (data − calculation) /
error. The lists have a check box per curve and ``All`` / ``None``; a right
click on a curve chooses its colour and marker (``Default style`` returns it to
the program's). The mouse
wheel zooms about the cursor (with Ctrl the x axis only, with Shift the y axis
only); dragging a rectangle with the left button zooms to it; dragging with the
right or middle button, or with Shift and the left button, pans; a double click
or ``Reset Zoom`` shows the whole range. The cursor's coordinates are shown in
the plot's corner and above it. ``Export`` writes PNG, SVG or PDF (on white,
whatever the window's look), or CSV: the points drawn, per curve a block of
the data (x, y, error) and one of the calculation (x, y and the band's limits
where there is one). The tab is read again after
every run; a file older than the project (records added since the run) is
named in a note below the plot.

The Physics tab
---------------

The Physics tab draws the quantities between the level parameters and the
observables, each from the engine's own functions (the design:
``docs/upgrade/physics_tab.md``; the same calls the web page's Physics tab and
``pyazr.physics`` read). A list of pages stands on the left. Every page works
the same way: choose what to draw from the drop-downs, the swept argument and
its range (from, to, the number of points, log spacing), then ``Draw`` (this
curve alone) or ``Add curve`` (over the others: several l, radii, materials);
``Clear`` empties the plot. The plot has a legend, the cursor's read-out, log
or linear axes, zoom and pan as the Plot tab, and ``Export`` (PNG, SVG, PDF or
the points as CSV). A line under the plot says what was left out and why.
The computations run in worker threads: a busy bar and ``Cancel`` stand at the
bottom while one runs, and the window stays usable.

Functions
  Every function of the physics table (kinematics, Coulomb functions, P, S and
  the hard-sphere phase, the Wigner limit and single-particle width, the THM
  plane-wave vertex, ...): its arguments as fields with their units, started at
  the row's default or example values, its documentation above them, the
  swept argument as a drop-down. A point outside an argument's domain is a gap
  in the curve, and the line under the plot counts such points with their
  reason (``e_cm outside its domain (finite > 0)``). No project is needed.

Channels of the project
  The project's channel functions of a pair (drop-down of the particle pairs),
  an l and the channel radius (empty: the pair's): P, S, the hard-sphere phase,
  F, G and their derivatives at the channel radius, against the pair's c.m.
  energy, as the calculation uses them (the Coulomb table, a pair's nuclear
  potential).

R-matrix of the project
  One J^π group on an excitation-energy grid: the phase shifts and the
  inelasticities of the open particle channels, \|U_cc'\|, the Argand plot
  (Im U_cc against Re U_cc, the energy along the curve, a marker at each level),
  \|T_cc'\|, P and S per channel, \|A_λμ\| and log10 \|det A^-1\|; a check list
  chooses the channels or pairs drawn. The group's levels are dashed lines and,
  beside the plot, a table: J^π, level, the physical E_x and Γ (MeV), and per
  channel its observed width W (MeV) or its ANC (fm^-1/2), the numbers of
  ``parameters.out``. A click on a level centres the plot on it.

Pole map
  A map with the excitation energy along x and one row per J^π group (log10
  \|det A^-1\|, one calculation per group), per pair of open channels of one
  group (\|U_cc'\|) or per pair of its levels (\|A_λμ\|); the colour is the value,
  read on the colour bar (dark blue low, yellow high). The levels are white
  points on their group's row (or white lines on a map of one group). The
  cursor's read-out gives the energy, the row and the value; the wheel zooms in
  energy (Shift: the rows), a dragged rectangle zooms, the right button pans, a
  double click shows the whole map; a click on a row draws its cut below. This
  is a map on the real energy axis: a minimum of \|det A^-1\| is where the
  inverse level matrix is nearly singular, the energy of a narrow resonance.
  The complex-energy plane, where the poles of the collision matrix are the
  zeros of det A^-1 continued in E, is not available: the engine's Coulomb
  functions and level matrix take real energies only.

THM
  For a project with ``<thm>`` experiments: the entrance vertex (⟨\|M_l\|^2⟩
  over the spectator window, \|M_l\|^2 at p_s = 0), the Coulomb line shape and
  the distortion factor R(E) and its parts, against the c.m. energy of x + A.
  The distortion runs a DWBA per energy (about 17 ms each): ``Cancel`` stops it
  between groups of ten energies.

Stopping
  SRIM's stopping powers, ranges, energy loss and straggling (the stopping
  table, :doc:`experimental_effects`) of a material given by its
  formula and of a projectile given by Z and mass. SRIM's citation stands under
  the plot. No project is needed.

The project pages read the *saved* project through the engine's session (the
levels and the data set up as a calculation would: up to a few seconds); they
say so when no project is saved. ``Parameters`` above them chooses the levels'
values: the project's start values, the last fit (``param.sav`` of the output
directory, when a fit has written it) or a parameter file. After an edit, save
the project: the pages set the session up again when they are next shown, or at
``Re-read``. Pages for the kernel shapes of the apparatus and for a level /
interference decomposition are not there yet (the engine has no call for
them; the decomposition is a pyazr script, :doc:`calculate`).

Runs
----

``Save and Run`` (or Calculate and Fit in the tool bar) saves the project and
runs the chosen command in a worker thread: the window stays usable, the run's
lines appear as the program prints them, and the progress line shows the fit's
objective, beside a small curve of the best χ² against the evaluation. When the
run ends the summary gives the total χ² (and χ²/N over all points) and a table
of the data sets: the program's segment number, the data set's name, its
number of points N, its χ² and χ²/N; a click on a column's title sorts by it.
``Open Plot`` shows the Plot tab with the new curves. ``Stop`` ends a fit
at its best point; the outputs are written there and ``fitStatus.out`` says
``stopped`` (HESSE, MINOS and the band are not run after a stop). A calculation
without a fit runs to its end once its set-up is done.

The configuration menu
----------------------

``Configure > Formalism`` sets the A-matrix or R-matrix route; ``Checks``,
``Directories`` and ``Runtime Options`` set the entries of ``<config>``
(:doc:`configure_menu`), each with the old program's title and the engine's
description. The dialogs are drawn from the engine's table of ``<config>``, so
an entry the engine no longer has has no control: fold windows, the
adaptive-grid switch and GSL Coulomb functions are not offered.

The Tools menu
--------------

What lives in ``pyazr`` and not in the engine is started from ``Tools`` as
``python -m pyazr ...`` (:doc:`pyazr_command_line`), a separate process whose
output appears in the Calculate tab's log: Export ENDF-6 (``endf write``),
Import an ENDF-6 Evaluation as a New Project (``endf to-project``), ENDF-6
File Summary, Import an EXFOR Data Set (``import exfor``), ENSDF Levels of a
Nuclide (``import ensdf``), Build and Check the Coulomb Table
(``coulomb-table``), and Open the Web Page of This Project (``pyazr ui``,
served on 127.0.0.1 and opened in the browser; it stops with the window). The
Python interpreter is set in ``Configure > Directories`` (a setting of this
computer, ``python3`` by default); when it cannot import ``pyazr`` every item
is disabled and says why. The pyazr samplers (NUTS, emcee) and band jobs are
not commands of pyazr's command line: they are in the web page and in Python.

From the command line and Python
--------------------------------

Everything the Qt program does is also done by ``AZURE2 <command>``
(:doc:`../reference/command_line`) and by ``pyazr`` (:doc:`pyazr`). For the
program's internals see ``docs/upgrade/qt_gui.md`` in the repository.
