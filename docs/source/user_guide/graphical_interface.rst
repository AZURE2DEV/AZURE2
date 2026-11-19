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

The Qt program needs Qt 6 (Core, Gui, Widgets, Network, Svg) and is built only
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

The engine checks the project as it stands: ``Configure > Validate Project``
(F7) fills the Diagnostics panel below the tabs with the diagnostics of
``AZURE2 --validate`` (severity, code, line, record and field, message). A
double click on a diagnostic shows its record.

The tabs
--------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Tab
     - What it edits
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
     - Targets, beams and detectors (``<targets>``, :doc:`experimental_effects`):
       layers, isotopes and abundances, stopping and straggling, fitted target
       parameters, resolution kernels.
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
       sampler where the build has it) with its options: the uncertainty band
       and its scaling, the parameter and capture files, the rate's pairs and
       temperatures, the sampler's walkers and steps; the threads; the run's
       lines, the chi-squared of each data set and the files the run wrote.
       The rate's entrance and exit pairs are chosen from the project's pairs
       (the exit also as -1, the total radiative capture).
   * - MCMC
     - Shown only when the engine has the ``mcmc`` command. The Gaussian priors
       of the sampler (the ``<parameters>`` rows of ``<mcmc>``: Load Parameters
       lists the free parameters), and the run's settings: walkers, steps,
       initial spread, level-energy spread, threads, the sampled basis (rwa or
       physical) and a fresh start (overwrite ``samples.mcmc``). Run starts
       ``AZURE2 mcmc`` with them; the progress bar follows the steps and the
       acceptance fraction is shown as the program prints it. The samples are
       in ``output/samples.mcmc``; the burn-in is dropped when they are read
       (pyazr), not by the sampler.
   * - Plot
     - The data sets and extrapolations of the last run (below).

A record's form is drawn from the engine's table of project format 3: a closed
choice is a list, a unit one of its dimension's (or of the data set's quantity),
a file a chooser; the tool tip of each field gives its meaning and grammar.
Fields set rarely are behind ``Every field`` unless the record sets them.
``Apply`` writes the record; nothing is written before.

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
then point one pair off). The isotope fields of the pair's dialog complete
from the mass table. Ctrl+Z undoes a pair's edit as in the Levels tab.

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
error. Drag a rectangle to zoom; a double click or ``Reset Zoom`` shows the
whole range. ``Export`` writes PNG, SVG or PDF. The tab is read again after
every run; a file older than the project (records added since the run) is
named in a note below the plot.

Runs
----

``Save and Run`` saves the project and runs the chosen command in a worker
thread: the window stays usable, the run's lines appear as the program prints
them, and the progress line shows the fit's objective. ``Stop AZURE2`` ends a fit
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
