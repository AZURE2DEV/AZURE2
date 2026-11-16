The Data Model and the Folding
==============================

How the data of a project are held, how the points of a data set are
calculated on the shared energy mesh of their compound nucleus, and how an
apparatus effect (a target, a beam spread, a detector resolution) folds the
cross sections at the mesh's nodes into the value of a point. The level matrix
behind each node is :doc:`level_matrix`; who owns each object is
:doc:`state_and_ownership`.

The objects
-----------

``DataSet`` (``include/DataSet.h``)
   The data as the set-up leaves it, immutable: the data sets (``Segment``),
   the target effects (``TargetEffect``), the leaves of sums and ratios
   (``LeafRef``), the covariance groups (``DataCovariance::Group``), the
   meshes of the compound nuclei (``NucleusMesh``), the THM experiments
   (``ThmGroup``) and the parameter layout.

``Segment`` (``include/Segment.h``)
   One data set (or extrapolation): its measurement, its settings and the
   terms of its prediction. Term 0 is its own reaction; term *k* is its
   *k*-th leaf, a component of a ``SUM`` or ``RATIO`` (``TermRef``). Its
   prediction is an expression tree built once at set-up from the spec
   (``experiment/include/azure2/experiment/Prediction.h``,
   ``include/Predictions.h``).

``WorkingCopy`` (``include/WorkingCopy.h``)
   What an evaluation changes: one ``SegmentState`` per segment (norm, energy
   shifts, shift function, sample thickness, the terms' frames), the present
   values of the apparatus parameters (``TargetEffectValues``), the points and
   slots in two flat vectors, and the calculation plan.

``Site`` (``include/Site.h``)
   One calculation site: an energy and an angle at which the point pass
   solves the level matrix and writes an observable. It holds what depends on
   the point at that energy: the c.m., lab and compound energies, the angle
   and its Legendre rows, the channel rows, the Coulomb amplitude, the
   screening factor, the stopping, the external-capture amplitudes, a THM
   site's vertex pieces of its entrance channels (``ThmVertexRows``), and
   where its result record lies. The channel rows (:math:`L^0_c`,
   :math:`\sqrt{P_c}` and the two phase factors of every channel) depend on
   the compound energy alone (``Site::ChannelRowsKey``): one table per energy,
   computed by the first site there and held by every other site at that
   energy and by its node of the ``NodeTable``
   (``WorkingCopy::FillChannelRows``). The members are in three groups by
   who writes them: what the set-up fixes (the segment, the lab angle, the
   nominal energy, the settings, the set-up's external-capture amplitudes),
   what an evaluation's hand-over writes (``Site::Values``: the energies and
   the c.m. angle after a shift, the Legendre and channel rows, the Coulomb
   amplitudes, the screening energy, the stopping and the thickness, the row
   set), and the place of the result record.

``Unit`` (``include/Unit.h``)
   A data point: a ``Site`` with its measurement row, its target-effect fold
   over its slots and its result and model records. An unfolded point is
   calculated at its own site; a folded point reads the slots of its range,
   which are ``Site`` objects held by the working copy
   (``WorkingCopy::SlotAt``, ``Unit::SubPointFirst``,
   ``Unit::NumSubPoints``), in descending energy.

``azure2::Results`` (``include/Results.h``)
   The buffer of one evaluation. Each unit (point or slot) owns one record at
   a fixed offset: its observable (``kSigma``: a cross section in b, a
   differential cross section in b/sr, a phase shift in degrees, an
   analysing power, a Legendre coefficient, the THM cross section), the E1
   and E2 parts of a capture cross section, :math:`A_y`, :math:`A_0`, then the
   angular-distribution coefficients and the abundance-weighted terms of a
   sample. ``Results::model`` holds the model value of every data or test
   point.

The shared mesh
---------------

A folded point needs the cross section at the energies its kernel covers.
Instead of a grid per point, every folded point of one compound nucleus reads
one mesh of compound energies :math:`E_x` (MeV):

- **Planning.** ``PlanMesh`` (``include/MeshPlanner.h``) is a function of
  values: the anchors of the nucleus's levels (energy and width, quantized),
  one support per folded point (the interval of :math:`E_x` where its kernel is
  non-zero, with its data set's numerical settings) and the tail constants.
  It returns the nodes (``MeshPlan``) and, per support, the range of nodes its
  quadrature uses. Smooth nodes are multiples of a power of two, so that
  every support that asks for a node gets the same bits; near a level the
  nodes follow the level's width (a core of fixed pitch, then a geometric
  tail). The interval that contains an energy is never longer than the step a
  point's own adaptive grid (``AdaptiveIntegrationGrid``) would have there.
- **Set-up.** ``WorkingCopy::SetUpTargetEffects`` (``src/TargetSetup.cpp``)
  plans the mesh of each nucleus (``DataSet::Meshes``) and lays the slots out
  (``AddMeshSlotGroups``, ``src/MeshSlots.cpp``).
- **Moving levels and shifts.** The anchors move with the level parameters
  and a point's support with its energy shift. The hand-over
  (``WorkingCopy::ReplaceSubPointGrids``) plans again and re-places the slots
  only in the region that changed (``MeshChangedRegion``); a shift that keeps
  a point inside its planned nodes moves nothing (``include/ShiftedFold.h``).

Slots shared by points
----------------------

Where several points of one data set need the cross section at the same node
with everything else that enters it equal, it is calculated once
(``include/SlotGroup.h``). Such points form a sharing group: the group holds
one ``Site`` per node of the union of its points' ranges, in node order, and
each point's slots are a contiguous range of that block. The key of a group
(``SlotGroupKey``) holds the term's settings (the observable, the pairs, the
external-capture pathways; ``SameObservable``, ``include/TermSettings.h``) and
the target effect's stopping, screening and Legendre ``Q`` coefficients; not
the segment, so two data sets of one reaction share their sites, and not the
angle: a differential site holds its observables without an angle and its
readers bring their angular rows.

The angular rows of the readers of one block are a ``RowSet``
(``include/RowSet.h``), built once per angular identity
(``WorkingCopy::BuildRowSets``): :math:`P_L(\cos\theta)` with the ``Q``
coefficients, the Coulomb amplitudes of elastic scattering at each slot, the
harmonics of an analysing power. Data given at a lab angle are folded in the
lab: at each slot the c.m. angle and :math:`d\Omega_\mathrm{lab}/d\Omega_\mathrm{cm}`
follow from the slot's energy (``experiment::LabAngleAt``). The rows
contracted with the slots' records are the values a reader folds
(``Results::contracted``, ``WorkingCopy::ContractRowSets``).

A slot is written by the first point, in data order, whose range holds it
(``WorkingCopy::WritesSlot``), and read by every point of its group: the fold,
the model pass, the gradient (a slot's cotangent is the sum of its readers').
A slot whose solve fails fails every point that reads it.

The plan and the point pass
---------------------------

``CalculationPlan`` (``include/CalculationPlan.h``) has one task per point,
built at set-up (``WorkingCopy::BuildPlan``) and checked when built: every
point is calculated by exactly one task, which writes the fixed part of its
record and of the slots it writes, and nothing else. The ``NodeTable``
(``include/NodeTable.h``) lists every compound energy at which the pass solves
the level matrix -- an unfolded point or a slot -- once per (nucleus, route,
bits of :math:`E_x`).

``WorkingCopy::CalculatePointsInParallel`` (``src/EnergyPass.cpp``):

1. ``CalculationPlan::EnergyTasks`` groups the units of the work by node into
   energy tasks (``EnergyTaskList``); inside a task the units are grouped by
   (entrance, exit, source of the external-capture amplitudes) into T-matrix
   groups. The list is kept between passes and rebuilt only when a unit's
   energy or amplitude source has changed.
2. Each task runs ``SolveAtOneEnergy``: the level matrix filled and factored
   once, a T-matrix per group, the observable of each unit into its record.
   A THM site reads the level matrix's elements instead
   (``THMMatrixFunc``).
3. Each folded point is folded (``Unit::FoldSubPoints``,
   ``src/PointPass.cpp``).

With external capture, the slots' energy-dependent amplitudes are prepared
before the pass (``Unit::PrepareSubPoints``).

Kernels: the fold
-----------------

A target effect (``TargetEffect``, ``include/TargetEffect.h``) is built from
the target, beam and detector records a data set names
(``experiment/src/BuildEffects.cpp``). It selects its kernels from the table
``Kernels()`` (``experiment/src/Kernels.cpp``; ``TargetEffect::SelectKernels``)
and gives them their parameters (``TargetEffect::KernelParametersFor``): an
optional *broadening* applied to the node values, then an optional
*resolution* kernel that folds them into the point's value
(``EffectKernels``). The rows:

.. list-table::
   :widths: 18 82
   :header-rows: 1

   * - Kernel
     - What it does (file in ``experiment/src/``)
   * - ``target``
     - Integration over the depth of the target (thin, thick, layer yield),
       with the beam spread and the straggling as a Gaussian average at each
       depth (``KernelTarget.cpp``).
   * - ``gaussian``
     - Normalised Gaussian average about the point's energy, the width a
       constant or a formula of the beam record (``KernelGaussian.cpp``).
   * - ``beam_profile``
     - A weighted sum of skewed Gaussians: a broad gamma beam with an
       event-by-event reconstructed energy (``KernelBeamProfile.cpp``).
   * - ``udr``
     - A user-defined numerical resolution: a time-of-flight delay
       distribution per point from the detector's table (``KernelUdr.cpp``).
   * - ``doppler``
     - Free-gas Doppler broadening of the node values, before the resolution
       kernel (``KernelDoppler.cpp``).
   * - ``rectangular``
     - Average over a window of full width :math:`w` (``KernelRectangular.cpp``).
   * - ``triangular``
     - Average weighted by a triangle of base :math:`w`
       (``KernelTriangular.cpp``).

A kernel is plain functions of the point's nodes (``KernelNodes``: energies in
c.m. MeV, descending, and for a kernel that integrates over a target the
stopping at each node), the point (``KernelPoint``) and the effect's parameters
(``KernelParameters``). Its fold is a template on the value type: on
``double`` it gives the value, on ``LinearForm`` the exact weights
:math:`c_j` of :math:`y = \sum_j c_j \sigma_j`, from the same arithmetic. The
fold of an effect (``Fold``, ``FoldForm``) is the one dispatch: the broadening
at every node, then the resolution kernel. A point whose support is narrower
than ``MeshZeroWidth()`` has a single node and takes the kernels' thin limit.
On the shared mesh the two outer nodes of a point lie beyond its support and
serve only to interpolate: a kernel integrates over the support
(``KernelNodes::Low``, ``KernelNodes::High``), never over the end nodes' span.

In the engine, ``FoldPoint`` (``src/PointFold.cpp``) writes the folded value
into the point's own record from its slots; ``FoldValues`` folds given values
(the transform of a transmission or a yield, a Coulomb-reference leaf);
``FoldWeights`` gives the gradient the exact weights; the extent of the mesh
a point needs is ``EffectSupport``. A kernel's parameter with
``fittable = true`` becomes a fit parameter when the project says ``fit``
(``ApparatusRows()``, ``src/ApparatusParameters.cpp``). The stopping of a
target is computed in one place, ``TargetStopping::PerActiveNucleus``
(``src/TargetStopping.cpp``), lab MeV cm\ :sup:`2` per active nucleus.

The model pass and the residual
-------------------------------

After the point pass, ``azure2::EvaluateModel`` (``src/ModelPass.cpp``)
forms what every data or test point predicts from the leaf values in its
records: the base point's value (or its E1 or E2 part), times its abundance
in a sample, plus the scaled components of a ``SUM``, over the denominator of
a ``RATIO``, or the transmission or self-shielded yield of the weighted terms.
A leaf value is never overwritten by a model value. Values the R-matrix does
not compute (the Rutherford or Mott cross section of a ratio to Coulomb, a
tabulated reference) are leaves computed in ``PointInputs``
(``src/CoulombReference.cpp``, ``src/ReferenceTableLeaf.cpp``) and folded on
the point's own nodes by ``FoldReferenceAtPoint``.

The residual of a datum is written in one place (``Residual``,
``include/ModelPass.h``):

.. math::

   r = \frac{f/J - y\,n}{\delta y\,n},

with :math:`f` the model value (c.m. frame), :math:`J` the
:math:`d\Omega_\mathrm{lab}/d\Omega_\mathrm{cm}` of the row at its present
energy shift (1 for c.m., angle-integrated and ratio data), :math:`y` and
:math:`\delta y` the datum and its uncertainty as read, and :math:`n` the data
set's normalization, on the data side. Every path uses it: the objective,
the per-segment :math:`\chi^2`, the least-squares rows, the gradients, the
API, the covariance groups, THM.

Transmission and yield
~~~~~~~~~~~~~~~~~~~~~~

A transmission or a self-shielded yield is a nonlinear function of its cross
sections. Its points keep the plain cross sections (and share the solves of
other points at their energies); the model pass sums the terms at the point,
or under a target effect at every slot, transforms them, and folds the
transformed values once more with the broadening already applied
(``Broadening::kAlreadyApplied``). The gradient takes each term back through
:math:`f'` times the kernel weight (``FoldWeights``). Test:
``tests/transmission_yield``.

Failed points
~~~~~~~~~~~~~

A unit whose calculation throws (a singular matrix, a Coulomb function out of
range) is marked failed (``Site::MarkFailed``): its values are NaN, and a
warning naming the data set, the energy and the error is printed the first
time in a run.
A ``SUM`` or ``RATIO`` with a failed or non-finite term is NaN at that point,
and so is a ``RATIO`` whose denominator is 0. The :math:`\chi^2` is then NaN,
so a failed evaluation cannot pass for a result.

Checks
------

Environment switches that check the structure at run time
(:doc:`../reference/command_line`, "Environment switches"):
``AZURE_PLAN_CHECK=1`` checks that each step of the point pass writes only
its own slots and that each unit gets from its task's shared solve the bits
its own solve gives; ``AZURE_NODE_CHECK`` compares the node table with the
units' own energies and rows (``WorkingCopy::CheckNodeTable``,
``WorkingCopy::CheckReaderRows``); ``AZURE_POISON_RESULTS=1`` fills a buffer
with NaN before each evaluation, so that a value read before it is written
shows. Consistency checks 28 (one model pass and one residual) and 29 (the
calculation plan) hold the structure in the source.
