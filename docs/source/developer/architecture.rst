Code Architecture
=================

AZURE2 is one C++ calculation engine with three front ends: the Qt setup
utility, the interactive command line, and the Python package ``pyazr``. All
three read the same ``.azr`` project file and must produce the same numbers.
This page maps where things live; :doc:`physics_core`, :doc:`level_matrix` and
:doc:`edata_and_target_effects` describe the physics inside the engine.

Repository layout
-----------------

.. list-table::
   :widths: 22 78
   :header-rows: 1

   * - Path
     - Contents
   * - ``src/``, ``include/``
     - The engine. Mostly one class per file pair (``CNuc.cpp`` / ``CNuc.h``).
   * - ``api/``
     - ``AZUREAPI``, a session object around the engine, and the pybind11
       module ``_azure2`` that exposes it to Python.
   * - ``pyazr/``
     - The Python package: the in-process driver, an ``.azr`` reader/writer,
       width transformations, uncertainty bands, and the EXFOR/ENSDF client.
   * - ``gui/``
     - The Qt5/Qwt setup utility (see :doc:`gui_structure`).
   * - ``numcmc/``
     - The affine-invariant ensemble MCMC sampler.
   * - ``external/``
     - Vendored upstream code, built in-tree and not modified: Minuit2, the
       Coulomb-function library ``coul``, and ``erya`` (SRIM stopping powers,
       pugixml).
   * - ``tests/``
     - The test suites (see :doc:`contributing`).
   * - ``docs/``
     - This site (Sphinx) and the C++ API reference (Doxygen).

A run, from the command line
----------------------------

1. ``src/AZURE2.cpp``, ``main``: handles ``--help`` and
   ``--build-coulomb-table``, starts the GUI unless ``--no-gui`` is given,
   parses the remaining flags into ``Config`` (``parseOptions``), reads the
   ``<config>`` block, and pre-scans the file for the segment pairs and
   external capture. It then shows the menu and runs ``AZUREMain`` for the
   chosen mode.
2. ``src/AZUREMain.cpp``, ``AZUREMain::operator()``:

   * ``CNuc::Fill`` reads ``<pairs>`` and ``<levels>`` (each level line through
     ``NucLine``). ``CNuc::Initialize`` builds the channels, the
     :math:`J^\pi` groups and the boundary conditions, and transforms the
     input parameters.
   * ``EData::Fill`` reads ``<segmentsData>`` / ``<segmentsTest>`` (through
     ``SegLine`` / ``ExtrapLine``) and the data files into ``ESegment`` and
     ``EPoint`` objects, plus ``<targetInt>`` (experimental effects).
     ``EData::Initialize`` sets up kinematics, external-capture integrals,
     Coulomb functions and the experimental-effect grids.
   * ``AZUREParams`` lays the fit parameters out for Minuit2.
   * ``AZURECalc`` is the objective: the total :math:`\chi^2` including the
     normalization, energy-shift, screening and abundance penalties.
     ``AZUREGrad`` supplies its analytic gradient and Jacobian, and
     ``AZURECalcMCMC`` the likelihood for ``numcmc``.
   * Minuit2's MIGRAD (or one of the least-squares minimizers) fits. Then
     ``EData::WriteOutputFiles`` writes ``output/``, and ``ReactionRate``
     integrates rates when asked.

3. The R-matrix itself: ``AMatrixFunc`` (the level matrix), ``GenMatrixFunc``,
   ``LevelMatrixSolver`` and ``ChannelSpace`` (the channel-space solution), and
   ``CoulFunc`` / ``CoulombTable`` / ``ShftFunc`` (penetrabilities and shift
   functions). External capture is in ``ECIntegral``, ``ECTable`` and
   ``ECAmplitudeCache``. Experimental effects are in ``TargetEffect``,
   ``Straggling``, ``AdaptiveIntegrationGrid`` and ``UdrKernel``.

The Python path
---------------

``pyazr.azure2`` creates an ``AZUREAPI`` session (``api/src/AZUREAPI.cpp``)
inside the Python process through ``_azure2``. The session holds a filled
``CNuc`` and ``EData``, and evaluates :math:`\chi^2`, gradients, residuals and
observables on parameter vectors without writing files. ``AZUREAPI`` has its
own entry points for the objective and its gradient. A change to what
``AZURECalc`` or ``AZUREGrad`` computes therefore has to be made there too,
and ``tests/pyazr`` compares the two.

One file, three readers
-----------------------

The ``.azr`` file is parsed independently by:

* the engine: ``Config.cpp`` (``<config>``), ``CNuc.cpp`` with ``NucLine.h``
  (``<levels>``), ``EData.cpp`` with ``SegLine.h`` / ``ExtrapLine.h``
  (segments and experimental effects), plus the pre-scans in ``AZURE2.cpp``;
* the GUI: ``AZURESetup.cpp`` (``readFile`` / ``writeFile``) and the tab
  models, with the segment keyword list in ``gui/include/SegmentKeywords.h``;
* pyazr: ``pyazr/azrfile.py`` (``AzrModel``) and ``pyazr/datasets.py``
  (``SegmentSet``).

A new keyword or field is therefore implemented three times. The GUI and pyazr
must also write the file back unchanged, including blocks they do not
understand. ``tests/gui/segments_roundtrip_test`` and
``tests/pyazr/roundtrip_test`` check the round trip, and ``tests/consistency``
checks that the three segment-keyword lists agree. The full list of places a
feature touches is in :doc:`contributing`.

Caches and their signatures
---------------------------

Expensive intermediate results are written next to the project and reused:
external-capture integrals (``output/intEC.dat``, ``intEC.extrap``), their
energy tables, and the Coulomb-function table (``~/.azure2/coulomb_table.bin``).
Each file carries a signature of everything that determines its contents (for
example ``EData::ECSignature``) and is recomputed when the signature does not
match. A change that alters what a cache would contain must change its
signature, or old files are silently reused with the old physics.
