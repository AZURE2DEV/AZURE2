EData, Data Points, and Target Effects
=======================================

This page describes how experimental data is managed internally and how
target integration, beam convolution, and straggling effects are applied.
This is essential reading for anyone working on improvements to the
convolution routines.

Data Structure Overview
-----------------------

The data hierarchy from top to bottom:

``EData``
   Top-level container. Holds all segments, target effects, and manages the
   overall calculation loop. One ``EData`` object exists per calculation; it is
   cloned for thread-safe parallel fitting.

``ESegment``
   A calculation segment corresponding to a specific entrance/exit particle
   pair, energy range, angle range, and data type. Contains a vector of
   ``EPoint`` objects.

``EPoint``
   A single data point with energy, angle, cross section, and uncertainty.
   When target effects are active, an ``EPoint`` also holds a vector of
   **sub-points** (``integrationPoints_``) used for numerical integration.

``TargetEffect``
   Configuration for experimental effect corrections: Gaussian convolution
   parameters, stopping power equation, density, Q-coefficients, straggling,
   and the number of integration sub-points.

**Key files:**

- ``include/EData.h``, ``src/EData.cpp``
- ``include/ESegment.h``, ``src/ESegment.cpp``
- ``include/EPoint.h``, ``src/EPoint.cpp``
- ``include/TargetEffect.h``, ``src/TargetEffect.cpp``
- ``include/AdaptiveIntegrationGrid.h``, ``src/AdaptiveIntegrationGrid.cpp``

Data Loading and Initialization
-------------------------------

EData::Fill()
^^^^^^^^^^^^^

The ``EData::Fill()`` method (in ``src/EData.cpp``) orchestrates the loading
and initialization of all data:

1. **Read data files** -- for each data segment, parse the four-column data
   file and create ``EPoint`` objects for data points that fall within the
   segment's energy and angle ranges.

2. **Convert to CM frame** -- laboratory energies and angles are converted to
   center-of-mass frame values using the kinematics of the entrance/exit pair.

3. **Map points at a shared energy** -- points of the same reaction at the
   same energy, in one segment or in several, share one energy-dependent
   calculation: the first of them in data order (the host) calculates it and
   the others are mapped onto it via ``EnergyMap`` (see
   `Point mapping across segments`_ below).

4. **Create sub-points** -- for each data point that has an associated target
   effect, a grid of sub-points is generated for numerical integration
   (see below).

5. **Initialize sub-points** -- each sub-point's energy-dependent quantities
   (penetrabilities, Coulomb phases, shift functions) are pre-computed via
   ``EPoint::Initialize()``.

Sub-Point Creation
------------------

When a data point has a target effect, the code creates a grid of sub-points
spanning the energy range needed for integration. The sub-points are stored
inside the parent ``EPoint`` in the ``integrationPoints_`` vector.

Energy Range Determination
^^^^^^^^^^^^^^^^^^^^^^^^^^

The sub-point energy range depends on which effects are active:

**Target integration only:**

- From the beam energy (surface) down by the target thickness:

  .. math::

     E_\text{back} = E_\text{CM} - \Delta E_\text{target}

  where :math:`\Delta E_\text{target} = \epsilon(E) \cdot \rho` is the target
  thickness computed from the stopping power :math:`\epsilon(E)` and the
  areal density :math:`\rho`.
- Forward depth is zero (no energy above beam energy).

**Gaussian convolution only:**

- Symmetric range around the beam energy:

  .. math::

     E_\text{range} = E_\text{CM} \pm n_\sigma \cdot \sigma_b

  where :math:`\sigma_b` is the beam energy resolution (Gaussian sigma) and
  :math:`n_\sigma` is the ``convolutionRange`` parameter (typically 5).

**Target integration + convolution:**

- Backward: target thickness plus the convolution tail:

  .. math::

     E_\text{back} = E_\text{CM} - \Delta E_\text{target} - n_\sigma \cdot \sigma_b

- Forward: convolution tail above the beam energy:

  .. math::

     E_\text{forward} = E_\text{CM} + n_\sigma \cdot \sigma_b

**Straggling:**

When straggling is enabled, the backward range is further extended by the
straggling width at the back of the target:

.. math::

   \sigma_\text{straggling} = c_s \sqrt{\Delta E_\text{target}}

where :math:`c_s` is the straggling coefficient in keV\ :sup:`1/2` and
:math:`\Delta E_\text{target}` the energy loss in keV, **both in the centre of
mass**: the coefficient is a c.m. quantity and no frame factor is applied
(``TargetEffect::StragglingSigma``). With a coefficient of 0 the width is
Bohr's instead (``TargetEffect::BohrStragglingSigma``),

.. math::

   \sigma_\text{straggling} = \frac{m_2}{m_1 + m_2}\,
   \sqrt{4\pi Z_1^2 e^4 Z_2\, N x},

with :math:`N x` the areal density of target atoms traversed, the charges
those of the entrance pair (``TargetEffect::SetBohrCharges``, set beside
``ConvertSigmaToCM``), and :math:`e^2 = \alpha\hbar c`. For the grid extent
:math:`N x` is the line's density (a thick target: the thickness over the
stopping cross section at the surface).

Adaptive Grid Generation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   The grids are built in ``EData::Fill``, which runs *before* the input
   parameter transformation, so at that moment ``ALevel::GetGamma`` still
   returns the values read from the file (an observed partial width in eV for
   an open particle channel, an ANC for a closed one, :math:`\Gamma_\gamma`
   in eV for a photon channel), not reduced-width amplitudes.
   ``GridConfig::inputWidthsArePhysical`` (set from
   ``Config::TRANSFORM_PARAMETERS`` and ``CNuc::IsTransformedIn``) tells
   ``IdentifyResonances`` to sum those widths directly; treating them as
   amplitudes made every estimate degenerate to :math:`2P/(dS/dE)` and left
   narrow resonances unresolved (a 0.6 keV 2\ :sup:`+` in :sup:`12`\ C+α was
   reported 0.8 MeV wide).


Sub-points are placed on an **adaptive energy grid** generated by the
``AdaptiveIntegrationGrid`` class. The grid is denser near resonances and
coarser in smooth regions:

1. **Resonance detection** -- the generator scans all compound nucleus levels
   and identifies resonances that fall within the integration range. For each
   resonance, the total width :math:`\Gamma` is estimated from the reduced
   width amplitudes.

2. **Step size calculation** -- the local step size varies smoothly using a
   Gaussian falloff:

   .. math::

      \Delta E = \Delta E_\text{fine} \cdot e^{-d^2/2\sigma_r^2}
      + \Delta E_\text{base} \cdot (1 - e^{-d^2/2\sigma_r^2})

   where :math:`d` is the distance to the nearest resonance,
   :math:`\Delta E_\text{fine} = \Gamma / N_\text{ppw}` (points per width),
   and :math:`\sigma_r = \Gamma \cdot M / 2` (with :math:`M` being a multiplier).
   This ensures narrow resonances get a fine grid without discontinuities.

3. **Grid construction** -- starting from the highest energy, points are placed
   at adaptively determined intervals down to the lowest energy. The total
   number of points is bounded by ``maxPoints`` (the user-specified integration
   points count).

Each sub-point is an ``EPoint`` object with:

- The same angle as the parent point.
- Energy set to the grid point energy.
- Stopping power pre-computed at that energy (for target integration).

Calculation Flow for Target Effects
------------------------------------

During calculation (in ``AZURECalc``), the processing of a point with target
effects follows this order:

1. **Calculate each sub-point** -- the R-matrix cross section is computed at
   each sub-point energy by calling ``EPoint::Calculate()`` on each sub-point.
   This runs the full FillMatrices/InvertMatrices/CalculateTMatrix/
   CalculateCrossSection pipeline for each sub-point.

2. **Integrate** -- ``EPoint::IntegrateTargetEffect()`` is called on the parent
   point to numerically integrate the sub-point cross sections, yielding the
   effective yield that accounts for target and beam effects.

Integration Methods
-------------------

The integration is implemented in ``EPoint::IntegrateTargetEffect()``
(``src/EPoint.cpp``). All cases use **2-point Gauss-Legendre quadrature** on
each sub-interval, which is exact for cubic polynomials and handles
non-uniform (adaptive) grids correctly.

Cross section values between grid points are obtained by **linear
interpolation** from the two nearest sub-points.

Convolution Only
^^^^^^^^^^^^^^^^

For pure Gaussian convolution (no target integration), the yield at beam
energy :math:`E_0` is:

.. math::

   Y(E_0) = \int \sigma(E') \cdot g(E' - E_0) \, dE'

where :math:`g` is the Gaussian beam profile:

.. math::

   g(E' - E_0) = \frac{1}{\sqrt{2\pi}\,\sigma_b}
   \exp\!\left(-\frac{(E' - E_0)^2}{2\sigma_b^2}\right)

**Implementation:** The code loops over consecutive sub-point pairs
:math:`[E_i, E_{i+1}]`. For each interval, two Gauss-Legendre points are
evaluated. At each Gauss point, the cross section is linearly interpolated
from the bracketing sub-points, multiplied by the Gaussian convolution factor
(``TargetEffect::GetConvolutionFactor()``), and accumulated into the integral.

Target Integration Only
^^^^^^^^^^^^^^^^^^^^^^^^

For target integration without convolution, the yield is:

.. math::

   Y(E_0) = \int_{E_0 - \Delta}^{E_0} \frac{\sigma(E')}{\epsilon(E')} \, dE'

where :math:`\epsilon(E')` is the stopping cross section. The integral runs
from the beam energy at the target surface down to the energy at the back of
the target.

**Implementation:** Similar Gauss-Legendre quadrature over sub-intervals,
but the integrand is :math:`\sigma / \epsilon` (cross section divided by
stopping power). Only intervals within the target range
:math:`[E_\text{back}, E_\text{surface}]` contribute.

Target Integration + Convolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The full nested integral is:

.. math::

   Y(E_0) = \frac{1}{\rho} \int_{E_\text{back}}^{E_\text{surface}}
   \left[\int \frac{\sigma(E')}{\epsilon(E')} \cdot g(E' - E_d) \, dE'\right] dE_d

This is a **double integral**: the outer integral runs over target depth
(energy loss), and the inner integral convolves the cross section with the
beam energy distribution at each depth.

**Implementation:** Two nested Gauss-Legendre loops:

1. **Outer loop** (target depth): iterates over sub-point intervals within
   the target range :math:`[E_\text{back}, E_\text{surface}]`. At each Gauss
   point :math:`E_d`, the effective beam sigma is computed (including
   straggling if enabled).

2. **Inner loop** (convolution): for each depth :math:`E_d`, integrates
   :math:`\sigma(E')/\epsilon(E') \cdot g(E' - E_d)` over the range
   :math:`E_d \pm n_\sigma \cdot \sigma_\text{eff}`, where
   :math:`\sigma_\text{eff}` includes both beam resolution and straggling:

   .. math::

      \sigma_\text{eff} = \sqrt{\sigma_b^2 + \sigma_\text{straggling}^2}

   The straggling contribution grows with depth (c.m. energies and a c.m.
   coefficient):

   .. math::

      \sigma_\text{straggling} = c_s \sqrt{E_\text{surface} - E_d}

   or, with :math:`c_s = 0`, Bohr's width for the areal density between the
   surface and :math:`E_d`, :math:`N x = \int_{E_d}^{E_\text{surface}}
   dE/\epsilon(E)`, accumulated over the sub-points (trapezoid in
   :math:`1/\epsilon`; a node with a non-positive stopping cross section adds
   nothing).

Target Integration + Straggling (no beam convolution)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When straggling is enabled but Gaussian beam convolution is not, the code
uses the same nested double integral structure, but the convolution kernel
is purely from straggling (no beam sigma component). At the target surface
(:math:`\sigma_\text{straggling} = 0`), the inner integral reduces to a
direct :math:`\sigma/\epsilon` evaluation (no spreading).

Key Helper Functions
--------------------

``TargetEffect::GetConvolutionFactor(energy, centroid)``
   Returns the Gaussian weight :math:`g(E - E_0)` for a fixed sigma.

``TargetEffect::CalculateConvolutionFactor(energy, centroid, config)``
   Returns the Gaussian weight with an energy-dependent sigma (evaluated from
   a user-provided equation).

``TargetEffect::TargetThickness(energy, config)``
   Returns the energy loss :math:`\Delta` across the target at lab energy
   :math:`E`, the solution of :math:`\int_{E-\Delta}^{E} dE'/\epsilon(E') = N`
   (Newton from the first-order :math:`\epsilon(E)\,N`, Simpson on 16
   intervals; exact for a constant :math:`\epsilon`). Up to 2026-10-27 it
   returned the first-order product itself.

``TargetEffect::GetStoppingPowerEq()``
   Returns the ``Equation`` object for the parametrized stopping cross section
   (user-defined functional form with parameters ``a0``, ``a1``, ...).

``AdaptiveIntegrationGrid::GenerateGrid(startEnergy, endEnergy, compound)``
   Generates an adaptive energy grid with finer spacing near resonances.
   Returns a ``std::vector<double>`` of energies from high to low.

Summary: Data Flow Diagram
--------------------------

::

   Data File (lab frame)
       │
       ▼
   EData::Fill()
       ├── Create EPoint objects (convert to CM frame)
       ├── Map points at a shared energy (within and across segments)
       └── For target-effect points:
           ├── Compute target thickness from stopping power
           ├── Determine integration range (target + convolution + straggling)
           ├── Generate adaptive energy grid (AdaptiveIntegrationGrid)
           └── Create sub-point EPoints at grid energies
                   │
                   ▼
   AZURECalc calculation loop (per data point):
       ├── For each sub-point:
       │       └── EPoint::Calculate()
       │           ├── ClearMatrices()
       │           ├── FillMatrices()      ← build A⁻¹ or R matrix
       │           ├── InvertMatrices()    ← GSL LU decomposition
       │           ├── CalculateTMatrix()  ← collision matrix
       │           └── CalculateCrossSection() ← σ(E) at sub-point
       │
       └── EPoint::IntegrateTargetEffect()
           └── Gauss-Legendre quadrature over sub-point grid
               ├── Convolution only: ∫ σ(E') g(E'-E₀) dE'
               ├── Target only: ∫ σ(E')/ε(E') dE'
               └── Both: ∫∫ σ(E')/ε(E') g(E'-E_d) dE' dE_d
                           ▼
                   Final yield stored in parent EPoint

Beam-profile kernel
-------------------

``TargetEffect`` carries an optional beam-profile kernel (``IsBeamProfile``):
skewed-Gaussian components ``BeamProfileComponent {xi, omega, alpha, weight}``,
a detector resolution ``beamTpcSigma_``, a truncation and a
detailed-balance flag, parsed from the ``beamprofile`` keyword block at the
end of the ``targetInt`` line (older readers stop in front of the keyword).
``EData::ReadTargetEffectsFile`` converts its energies lab → c.m. once
(``ConvertBeamProfileToCM``) and sizes the sub-point grid from
``BeamProfileSupport`` intersected with the point's window ± 4 s.  The window
comes from ``DataLine`` extras (columns 5–6, ``EPoint::HasBinWindow``), is
converted with the point in ``ConvertLabEnergy`` and stored unshifted; the
branch in ``EPoint::IntegrateTargetEffect`` shifts it by the segment's current
c.m. energy shift, normalises numerically (numerator and denominator on the
same Gauss–Legendre points) and applies the detailed-balance weight from the
kinematics stored by ``EPoint::SetPhotoKinematics`` to the numerator only
(normalised at the point's unshifted energy; in the denominator it would
cancel that normalisation).
``TargetEffect::IsSubPointEffect`` is the single gate for every effect that
integrates over sub-points; ``EPoint::IntegrateTargetEffectComponents``
re-runs the combiner on the E1 and E2 components for isDiff 5/6 segments.

The convolution has a reference calculation of its own in
``tests/reference/beam_profile_reference_test.cpp``: it drives
``EPoint::IntegrateTargetEffect`` with hand-chosen sub-point cross sections and
compares against closed forms (a constant cross section must come back exactly)
and against a second, densely integrated implementation of the documented
formula. Unlike ``tests/beam_profile_kernel``, which pins a chi-squared this
code produced, it can tell you the answer is right rather than merely
unchanged. It needs no model and runs in milliseconds.

On the GUI side the block lives in ``TargetIntData`` (``gui/include/TargetIntModel.h``,
columns 22-26 of the table model) and is parsed and emitted by
``TargetIntTab::readFile``/``writeFile``.  The reader splits the block off at
its keyword *before* the quoted ranges token is handled, so the older numeric
chain parses exactly as it did; the writer emits it last and only when the
effect declares one, and formats the numbers with ``QString::number(v, 'g', 12)``
because ``QTextStream``'s default of six significant digits would round a
profile location of a few MeV on every save.  A block that declares more
components than it supplies is dropped rather than half-applied.

Energy ranges, blending and the automatic decision
--------------------------------------------------

``TargetEffect`` optionally carries lab-energy windows, a blend width and a
relative tolerance (three trailing tokens of the ``targetInt`` line).  During
``EData`` fill, a point whose blend weight is zero never receives the effect
number, so the rest of the machinery treats it as an ordinary point.  A
parent point with a fractional weight, or a positive tolerance, first
computes its bare value from a probe copy of itself that pretends not to
carry the effect; edge blending mixes the integrated and bare values with a
smoothstep weight, and the automatic decision probes the first, central and
last sub-point to estimate the size of the effect (curvature across the
window, plus the centroid shift under target integration), using the bare
value outright when the estimate is below the tolerance.  Points with mapped
observables skip only when every mapped observable passes its own estimate.

The combined branch of ``EPoint::IntegrateTargetEffect`` accepts either
convolution flavour: a fixed beam sigma, or the energy-dependent convolution
equation evaluated at each depth energy.  The latter combination used to fall
through to pure target integration on a grid sized for the convolution,
corrupting the yield.

Point mapping across segments
-----------------------------

Angular distributions are often stored one angle per segment (an EXFOR entry
per angle, say): sixteen one-angle segments of 12 000 points at 800 distinct
energies are not unusual.
Everything that depends on the energy alone -- the Coulomb functions,
penetrabilities and phases (``EPoint::CalcEDependentValues``), the level
matrix, its factorization and the T-matrix -- is the same for all the angles
at one energy, and only the observable (Legendre polynomials, Coulomb
amplitude, :math:`A_y`) differs.

``EData::MapData`` therefore groups the points of all segments (except
composite SUM/RATIO and total-capture segments) by entrance pair, exit pair and
energy, and maps each point onto the first point in data order (segment, then
point index) that shares its calculation (``SharesCalculation``): the same
original, laboratory and c.m. energies, the same target effect, blend weight
and sub-point grid (energies and stopping powers, point by point), and the
same kind of observable (cross section or :math:`A_y`, angular-distribution
coefficients or not).  The observable type may otherwise differ: an
angle-integrated and a differential point share the T-matrix.  The
normalization plays no part, since it scales the data, not the model.

Across segments the energy shifts must also agree *as functions of the fit
parameters* (``ShiftsAgree``): the same free parameters (a shift group's
shared constant shift, a shared shift function), no free :math:`\sqrt{E}`
coefficient (it is each segment's own), and the same fixed values.  Two
segments mapped onto each other are then at the same energy whatever a fit or
a finite-difference step does to the shifts; segments whose shifts only
coincide at the current values are not mapped.  Every change of a shift or of
a sub-point grid rebuilds the mapping (``EData::FillEnergyShiftsFromParams``).
The segment-local energy-shift derivatives of ``AZURECalc::Gradient`` and
``AZURECalc::ResidualJacobian`` stay local as long as a parameter's segments
are mapped only among themselves; they move all of them before calculating
any.

The host is calculated first by every loop over the whole data set, and
calculates the observables of the points mapped onto it
(``EPoint::Calculate``); ``AZURECalc::PrecalculatePoints`` hands out hosts
only.  The angle-independent products :math:`Z_1Z_2\,T_1T_2^*` of the
differential sum are formed once per T-matrix and summed by Legendre order
(``GenMatrixFunc::CalculateCrossSection``, from the flat interference table of
the decay, ``Decay::FlatInterferences``), as are the Coulomb-nuclear
interference terms of elastic scattering per channel spin and :math:`l`
(``BuildCoulombInterference``); every point that shares the T-matrix then
needs one term per order.  Particle exits other than UPOS take the same path.
The results are byte-identical with and without mapping.  The adjoint does the
same: the T-matrix cotangents of a differential point are, per pathway, a sum
over orders of coefficients taken once per T-matrix
(``AMatrixFunc::TBarProducts``).

Energy-grouped analytic gradient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The analytic gradient and the residual Jacobian (``AZUREGrad.cpp``) work on
the same groups (``BuildGradItems``): one work item per host with the points
mapped onto it, in parallel over items, so no two threads touch one point.  An
item solves the level matrix and T-matrix once, forms each member's observable
from them, and back-propagates:

- for the scalar gradient (``AccumulateEGammaGradient``), the members'
  T-matrix cotangents are summed and back-propagated once
  (``AMatrixFunc::PointTBar``, ``TBarAdjoint``; the back-propagation is linear
  in them);
- for the rows of the residual Jacobian and the model gradients
  (``GroupRows``), the derivatives of the group's level-matrix bilinear forms
  :math:`g_c^T A\, g_{c'}` are taken once, for the cotangents 1 and *i*
  (``AMatrixFunc::BilinearAdjoint``), and every member's row is a combination
  of them with its own angle-dependent coefficients -- unless back-propagating
  each member is cheaper (few members, many channel pairs) or the decay has
  external capture.

Electron screening (a ``screening`` block on a ``<targetInt>`` line) is a
factor :math:`f = e^{k U_e}` on each point's own observable, applied after the
bare one is formed from the shared T-matrix (``GenMatrixFunc::CalculateCrossSection``),
so points with different screening may share a calculation.  In the adjoint
``PointTBar`` scales the point's cotangent by :math:`f` and records
:math:`\partial Y/\partial U_e`; ``TBarAdjoint`` adds it, and the
bilinear-basis rows of ``GroupRows``, which do not back-propagate the
member, add it to the member's row directly.

The cotangent of the level matrix is a sum of rank-one terms,
:math:`\bar M = -\sum_t \beta_t\, \overline{A g_c}\,\overline{A g_{c'}}^T`,
formed from the solves :math:`A g_c` in :math:`O(n^2)` per term instead of the
two :math:`O(n^3)` products :math:`A^H \bar A A^H` (kept for external capture),
and its structural derivatives (standard, Brune and Park level matrices) are
a few row-times-matrix products (``AMatrixFunc::LevelMatrixAdjoint``).
On the channel-space route (:doc:`level_matrix`) the adjoint of a bilinear form
is taken in the formal representation of the level matrix instead, in
:math:`O(n c)` per term with no :math:`n \times n` matrix
(``AMatrixFunc::FormalBilinearAdjoint``); for the Brune and Park matrices the
scalar gradient sums the formal cotangents over all points and resolves them
into the level parameters once per pass (``ResolveFormalCotangents``), while
the Jacobian rows take the derivatives by free parameter below, with the solves
:math:`A g_c` of the Woodbury form.  Where the scalar gradient does take the level-matrix
adjoint (the level route, the pencil route), it defers it too: the level
matrix is :math:`K - EN - G\Lambda(E)G^T` with :math:`K, N` independent of the
energy, so only :math:`\sum_E \mathrm{Re}\,\bar M(E)` and
:math:`\sum_E E\,\mathrm{Re}\,\bar M(E)` are accumulated (``GradAccum::lmK``,
``lmN``) and contracted with the structural derivatives of :math:`K, N` once
per pass; the kernel term is differentiated at each energy from the rank-one
terms of :math:`\bar M`, in :math:`O(n c)` each.  The rows (residual Jacobian,
model gradients) need each point's derivative with respect to every free
parameter, so they go the other way: a level parameter changes
:math:`K - EN` by rank two in its level's row and column,
:math:`dK = e_l k^T + k e_l^T`, :math:`dN = e_l q^T + q e_l^T` (``ChannelSpace::
CurrentLevelDerivatives``, per parameter set and free mask), and with the
solves :math:`u = A g` its derivative is :math:`O(n)` per free parameter
(``AMatrixFunc::LevelAdjoint``, ``LevelBasis``; standard, Brune and Park
level matrices).
Points with components or sub-point target effects are still differentiated
one by one (``GradOnePoint``), with the same back-propagation.  For a SUM or
RATIO point without a sub-point target effect, ``GradOnePoint`` solves the
level matrix once at the base point and takes each term's T-matrix and
cotangent from it; the terms without external capture enter only through the
bilinear forms :math:`g_c^T A g_{c'}`, which do not depend on the exit pair,
so their cotangents are merged by channel pair and back-propagated once
(``AMatrixFunc::TBarAdjointMany``).

Segments with components (SUM, RATIO, total capture) are mapped as well, since
October 2026 (``EData::MapComponentSegments``): each component's point is
mapped onto the base segment's point at the same energy.  A component has its
own exit pair, so the host's calculation builds the component's T-matrix from
the solved level matrix (``GenMatrixFunc::RecalculateTMatrix``) before its
observable; the mapped points of a host come grouped by exit, so each exit's
T-matrix is built once.  A component keeps and re-energises its own
energy-dependent values and external-capture amplitudes
(``EPoint::OwnsEDependentValues``): its segment is shifted in parallel with the
base segment.  A family is not mapped with other segments, and the base's
points stay unmapped among themselves (the external-capture file layout
depends on which points are mapped).  On the 13C+alpha project of
``13C+a/10-22-26_merged_non16o_non1``, whose neutron total cross section sums
nine exit pairs over 13 129 points, one evaluation went from 13.1 s to 4.4 s
(136 136 level-matrix solutions to 23 203) and one residual Jacobian from
517 s to 141 s (four threads), every output identical and the Jacobian equal
to :math:`2 \times 10^{-13}`.  The gradient
equals the point-by-point one to rounding (about :math:`10^{-14}` of its
largest component on such a project).

On a 19-segment project of that kind (14 406 points, 694 free parameters, one
thread) one evaluation went from 7.0 s to 2.2 s, one residual Jacobian (the
gradient and G2 of ``--use-gradient``) from 51 s to 4.8 s, and 20 MIGRAD
calls from 355 s to 68 s; peak memory went down slightly.  Solving the level
matrix in channel space (:doc:`level_matrix`), the sums by Legendre order and
the gradient changes above then took them to 0.36 s, 1.4 s (1.0 s for the
adjoint gradient) and 17.5 s (two threads: 0.21 s, 0.75 s and 11 s).

Environment switches
^^^^^^^^^^^^^^^^^^^^

``AZURE_NO_CROSS_SEGMENT_MAPPING=1``
   Map points within each segment only, as before October 2026.  (A shifted
   capture segment is always mapped within itself only: its external-capture
   amplitudes at a shifted energy may be interpolated from its own points.)
``AZURE_NO_POINT_MAPPING=1``
   No point mapping at all: every point calculates its own energy-dependent
   values and T-matrix.
``AZURE_NO_COMPONENT_MAPPING=1``
   Every term of a SUM, RATIO or total-capture segment is its own calculation,
   in the evaluation and in the gradient, as before October 2026.
``AZURE_NO_GROUPED_GRADIENT=1``
   Differentiate every point on its own (``GradOnePoint``), as before.
``AZURE_GRAD_DENSE=1``
   Back-propagate through the level matrix with the dense
   :math:`A^H \bar A A^H` products.
``AZURE_LEVEL_SPACE=1``
   Solve every J-group with the level matrix (fill and LU), not in channel
   space (:doc:`level_matrix`).

None of them changes a result beyond rounding; they exist to check that and
to time the alternatives.
