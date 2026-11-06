GUI Structure
=============

The setup utility (``AZURESetup``) is a Qt5 application linked into the
``AZURE2`` executable. It starts unless ``--no-gui`` is given. Its sources are
in ``gui/src`` and ``gui/include``. The Plot tab needs Qwt (``USE_QWT``), and
the whole GUI can be left out with ``-DBUILD_GUI=OFF``.

Main window and tabs
--------------------

``AZURESetup`` owns the main window, the menus and the ``.azr`` file as a
whole (``readFile`` / ``writeFile``). Each block of the file belongs to one tab
and one Qt model:

.. list-table::
   :widths: 30 35 35
   :header-rows: 1

   * - Tab
     - Classes
     - ``.azr`` block
   * - Particle Pairs
     - ``PairsTab``, ``PairsModel``, ``AddPairDialog``, ``MassCheckDialog``
     - ``<pairs>``
   * - Levels and Channels
     - ``LevelsTab``, ``LevelsModel``, ``ChannelsModel``, ``AddLevelDialog``
     - ``<levels>``
   * - Segments
     - ``SegmentsTab``, ``SegmentsDataModel``, ``SegmentsTestModel``,
       ``ScatteringLengthsModel``, ``AddSegDataDialog``, ``AddSegTestDialog``
     - ``<segmentsData>``, ``<segmentsTest>``, ``<scatteringLengths>``
   * - Experimental Effects
     - ``TargetIntTab``, ``TargetIntModel``, ``AddTargetIntDialog``
     - ``<targetInt>``
   * - Nuclear Potential
     - ``NuclearPotentialTab``
     - ``<potential>``
   * - Fitting Settings, Calculate, MCMC
     - ``FittingTab``, ``RunTab``, ``MCMCTab``
     - ``<config>``, ``<parameterSettings>``, ``<lastRun>``
   * - Plot
     - ``PlotTab``, ``AZUREPlot`` (Qwt)
     - reads ``output/``, not the ``.azr``

The sample abundances of natural targets (``<abundances>``, ``AbundancesModel``)
are shared by the Particle Pairs, Segments and Fitting Settings tabs.

Dialogs for data import are ``ExforDialog`` / ``ExforData`` (EXFOR; the same
Web-API as ``pyazr.nds``, so a parsing rule learned by one belongs in both) and
``IAEALevelDialog`` / ``IAEALevelData`` (ENSDF levels). The help text shown in
each tab is in ``InTabDocs.cpp``.

``MassTable`` reads the AME2020 atomic mass table that the GUI embeds as the
Qt resource ``:/ame/mass_1.mas20``. It is the file ``pyazr/data/mass_1.mas20``
that ``pyazr.masses`` reads, so both fill identical masses and separation
energies. A newer evaluation replaces that one file and the two version
strings named in ``pyazr/data/README.md``.

Running calculations
--------------------

``AZUREMainThread`` and ``AZUREMCMCThread`` run ``AZUREMain``, the same engine
entry point the command line uses, on a worker thread and report back with
signals. The GUI contains no physics of its own. Anything it displays is
computed by the engine.

Round trip
----------

Loading and saving a file without edits must reproduce it. A segment line's
keyword tail (``sqrtshift``, ``shiftfn``, ``normgroup``, ...) is kept verbatim
in ``extraTokens``, but only for the keywords listed in
``gui/include/SegmentKeywords.h``. A keyword the engine learns must be added
there. ``tests/consistency`` checks that this list matches the engine's and
pyazr's. The GUI tests in ``tests/gui`` run headless
(``QT_QPA_PLATFORM=offscreen``) on every CI platform.

Qt6
---

Homebrew disables ``qt@5`` and ``qwt-qt5`` on 2027-05-19, which sets the
deadline for porting to Qt6. New code should not use APIs that Qt6 removed or
moved (``QRegExp``, ``qrand``, ``QString::SkipEmptyParts``). New connections
should use the pointer-to-member form of ``connect``.
