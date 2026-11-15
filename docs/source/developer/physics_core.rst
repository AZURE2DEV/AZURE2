Physics Core
============

The library ``physics/`` (CMake targets ``azure2_physics`` and
``azure2_physics_transfer``) holds the physics of the R-matrix calculation as
functions of values: plain arrays and structs in physical units, no state, no
engine header. The engine, the binding and the web UI call the same
functions; every one of them is tested against a value from outside the code
(a closed form, an ``mpmath`` calculation, a published table, FRESCO).

Rules
-----

- It uses nothing else of AZURE2: no ``experiment/``, ``io/``, engine
  (``src/``), GSL or OpenMP. The target links nothing; the CTest test
  ``physics_standalone`` configures and builds ``physics/`` as a project of its
  own (``physics/tests/physics_standalone.cmake``); check 8 of the
  consistency test refuses any include but ``azure2/physics/`` headers and the
  C++ standard library.
- Public headers are in ``physics/include/azure2/physics/`` (included as
  ``"azure2/physics/X.h"``), sources in ``physics/src/``, tests in
  ``physics/tests/``. Everything is in the namespace ``azure2::physics``.
- One set of constants (``Constants.h``, CODATA 2018); the engine's
  ``include/Constants.h`` takes its values from there (check 10).
- Energies in MeV, lengths in fm, masses in u; each header states the frame
  (c.m. or lab) of what it takes.

What it contains
----------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Header
     - Contents
   * - ``System.h``, ``JpiBlock.h``
     - The reaction system: compound nuclei, partitions :math:`\alpha`,
       channels :math:`c = (\alpha, l, s)` and the :math:`J^\pi` blocks of
       levels (``azure2::physics::System``, ``JpiBlock``), plain values with
       0-based indices; spins stored doubled.
   * - ``Convention.h``, ``LevelParameters.h``
     - The bases of the level parameters (``LevelBasis``: formal, Brune,
       Park) and the coupling schemes, and the conversions between them:
       ``amplitudes_from_observed``, ``observed_from_amplitudes``,
       ``boundary_values``, ``formal_from_brune``, ``brune_from_formal``,
       ``park_from_brune``, ``brune_from_park``. Each works on one
       :math:`J^\pi` block of a ``System``.
   * - ``ChannelFunctions.h``, ``Derivatives.h``
     - Penetrability, shift function and hard-sphere phase of a channel from
       the Coulomb functions at :math:`\rho = k a_c`, their energy
       derivatives, :math:`\eta` and :math:`\rho` in non-relativistic or
       relativistic kinematics, the bound-state shift.
   * - ``LevelPairShift.h``
     - The shift function of one channel between two levels of one
       :math:`J^\pi` as the Brune and Park level matrices and the Brune
       :math:`\to` formal transformation hold it: the chord through
       :math:`(E_a, S_a)` and :math:`(E_b, S_b)`, its slope, its value at
       :math:`E = 0` and at the channel energy, and their derivatives by the
       level energies (``level_pair_shift``; header only, double and long
       double). Every place that writes these divided differences calls it
       (the level route, the R-matrix route, the adjoint, the derivative
       arrays, the transformed route, ``formal_from_brune``).
   * - ``CollisionMatrix.h``, ``LevelMatrixSolver.h``
     - The level matrix of one :math:`J^\pi` block at one energy from arrays
       (``LevelMatrixTables``; ``fill_level_matrix``: standard, Brune, Park,
       the Reich-Moore term), its factorization (``LevelMatrixSolver``), and
       the collision matrix: ``pathway_amplitude``, ``collision_element``,
       ``t_matrix_element`` per channel pair, the capture elements, and
       ``collision_matrix`` for a whole block (:math:`A`, :math:`U`,
       :math:`T`). :doc:`level_matrix`.
   * - ``ChannelSolve.h``
     - The channel route at one energy on the formal levels of a block
       (``channel_space::FormalBlock``, ``ChannelSolve``):
       ``solve_channel_space``, ``f_column``, ``apply_a``,
       ``channel_bilinear``, and the small LU routines.
   * - ``Coulomb.h``, ``Whittaker.h``
     - The Coulomb and Whittaker functions (``physics/coulomb/``, from
       pycoulomb, one directory per algorithm).
   * - ``AngularCoupling.h``, ``ObservableCoupling.h``
     - Clebsch-Gordan and Racah coefficients; the coupling coefficients of
       angular distributions, capture and analysing powers.
   * - ``Pathways.h``
     - The reaction pathways of one entrance partition.
   * - ``ExternalCapture.h``
     - Bound-state normalization and external reduced widths of capture.
   * - ``HybridPotential.h``
     - The nuclear potential of a hybrid particle pair and its radial
       equation.
   * - ``CoulombScattering.h``, ``AstrophysicalFactor.h``
     - The Rutherford and Mott cross sections; the S factor and its
       conversions.
   * - ``ResonanceQuantities.h``, ``Widths.h``
     - Resonance strength, total width, branching ratio; Wigner limit,
       single-particle width, photon penetrability.
   * - ``SymmetricEigen.h``
     - Small real symmetric (generalized) eigenproblems for the basis
       transformations.
   * - ``Elements.h``
     - Element symbols by :math:`Z`.

``azure2_physics_transfer`` (``physics/src/transfer``, headers in
``physics/include/azure2/physics/transfer/``) links ``physics/`` only:
optical potentials and the global optical sets (``OpticalPotential.h``,
``GlobalOptical.h``), bound states (``BoundState.h``), distorted waves and
continuum bins (``DistortedWaves.h``, ``Continuum.h``), the DWBA transfer
amplitude (``Transfer.h``) and the THM vertex (``ThmVertex.h``,
``ZeroRangeWeight.h``). It is checked against FRESCO
(``physics/tests/transfer_fresco_test.cpp``).

How the engine uses it
----------------------

The engine's compound-nucleus classes are the set-up view a ``System`` is
built from: ``BuildSystemOfCompound`` (``include/SystemOfCompound.h``) makes a
``SystemOfCompound``, the ``System`` with the map from its 0-based indices to the
engine's pairs, groups, levels and channels. ``CompoundSystem::Fill`` builds it once,
when the structure is complete (``CompoundSystem::SystemStructure``: the R-matrix
blocks; ``CompoundSystem::SystemAllGroups``: every group, for the reports), and every
reader takes that one; the structure does not change afterwards. ``CompoundSystem::TransformIn``,
``CompoundSystem::CalcBoundaryConditions`` and ``CompoundSystem::TransformOut`` call the
conversions of ``LevelParameters.h`` on it; ``CompoundSystem::BuildPathways`` builds
the ``PathwayTable`` from it. ``AMatrixFunc`` gathers a group's arrays at each
energy and calls ``CollisionMatrix.h``, ``LevelMatrixSolver.h`` and
``ChannelSolve.h`` (:doc:`level_matrix`). The channel functions reach the engine through
``ChannelFunc`` and the stage's ``ChannelTable``; the THM distortions
(``src/ThmDistortion.cpp``, ``src/ThmDwVertex.cpp``) call
``physics/src/transfer``.

The binding exposes the functions to Python with array arguments through one
table, ``BuildTable()`` in ``api/src/pybind_physics.cpp`` (``pyazr.physics``;
:doc:`extending`, "A physics function of the binding").

The ``System`` is built once at set-up and held by the ``CompoundSystem``, which is constant after it.

.. note:: **Planned.**
   The per-vector level values become arrays in ``JpiBlock``'s layout passed to these functions beside the
   ``System``.
