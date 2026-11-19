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

The Qt program needs Qt 6 (Core, Gui, Widgets, Network) and is built only when
asked::

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_GUI=ON \
         -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
   cmake --build build --parallel
   build/gui/src/AZURESetup [project.azr]

Qt 6 asks for CMake 3.22 or newer; with an older CMake add
``-DQT_FORCE_MIN_CMAKE_VERSION_FOR_USING_QT=3.20`` (Qt's own switch; the build
of this program needs nothing newer). Without ``-DBUILD_GUI=ON`` nothing of the
Qt program is built and every other target is unchanged. The Plot tab of earlier
versions needs Qwt and is not part of this version (``-DUSE_QWT`` is kept for it
and is off): plot the files in ``output/`` (:doc:`../reference/output_files`)
or use the browser page.

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
     - The Trojan horse experiments (``<thm>``), as text checked by the engine,
       including ``vertexModel=ff``, the full transfer form factor.
   * - Fitting Settings
     - The fitter (the engine's list of fitters), the call limit and the
       covariance refreeze passes; the prior covariance file; the parameter
       settings (``<parameterSettings>``) and the sampler's block (``<mcmc>``).
   * - Calculate
     - The command (calculate, fit, extrapolate, MINOS, reaction rate, and the
       sampler where the build has it) with its options: the uncertainty band
       and its scaling, the parameter and capture files, the rate's pairs and
       temperatures, the sampler's walkers and steps; the threads; the run's
       lines, the chi-squared of each data set and the files the run wrote.

A record's form is drawn from the engine's table of project format 3: a closed
choice is a list, a unit one of its dimension's (or of the data set's quantity),
a file a chooser; the tool tip of each field gives its meaning and grammar.
Fields set rarely are behind ``Every field`` unless the record sets them.
``Apply`` writes the record; nothing is written before.

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
description. Options the engine no longer has (fold windows, the adaptive-grid
switch, GSL Coulomb functions) have no control.

From the command line and Python
--------------------------------

Everything the Qt program does is also done by ``AZURE2 <command>``
(:doc:`../reference/command_line`) and by ``pyazr`` (:doc:`pyazr`). For the
program's internals see ``docs/upgrade/qt_gui.md`` in the repository.
