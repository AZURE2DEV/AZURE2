Output Files
============

AZURE2 produces several output files in the configured output directory. An
important distinction: while all **input** is in the laboratory frame, all
quantities in **output** files are in the **center-of-mass frame**.

Parameter Files
---------------

param.par
^^^^^^^^^

Contains the initial formal R-matrix parameters (energies, reduced width
amplitudes, etc.) derived from the Input File. Primarily a check file with
limited direct use.

param.sav
^^^^^^^^^

Contains the best-fit formal R-matrix parameters after a fit is completed. This
file can be loaded back into AZURE2 to reproduce a fit or to use as starting
parameters for subsequent calculations (extrapolation, reaction rate, etc.).
In a project with an ``<abundances>`` block the isotopic abundances follow the
level parameters as ``abundance_<Z>_<A>`` (:doc:`data_formats`). An
experimental-effect line with a ``screening`` block adds its screening potential
(eV) as ``screening_<n>`` after the energy-shift parameters, and a fitted target
density (``densityfit``, one per line or ``densitygroup``) adds ``density_<n>``
(atoms/cm²) after those.

The ``segment_<k>_*`` lines (normalizations and energy shifts) belong to the
data segments of the fit that wrote the file. A run without data -- an
extrapolation (**Calculate Segments Without Data**) or a reaction rate --
numbers its own test grids 1, 2, ... and does not apply those lines (the run
log says so): a test segment's energies are the ones its line asks for. Until
October 2026 data segment *k*'s fitted shift was handed to test segment *k* by
name and moved its grid. ``pyazr``'s ``load_sav`` does the same in
extrapolation mode.

parameters.out
^^^^^^^^^^^^^^

Contains the **physical** (observable) parameters resulting from the fit. If
the user wishes to use these as new starting values, they must be manually
entered into the **Levels and Channels** tab.

relativistic.out
^^^^^^^^^^^^^^^^

Written with ``relativistic=1`` (:doc:`../user_guide/particle_pairs`). The
levels of ``parameters.out`` in the level-energy convention of the LANL code
EDA (ENDF ``KRL=1``): for each level its excitation energy :math:`E_x`, the
energy :math:`E_s = S_0 + T_0 + T_0^2/(2 m_0)` (:math:`T_0 = E_x - S_0` on the
reference pair of mass :math:`m_0`), the difference :math:`E_s - E_x` and the
factor :math:`dE_s/dE_x = 1 + T_0/m_0`; per particle channel the partial width
and reduced width amplitude, and the same times :math:`1 + T_0/m_0` and its
square root. The mapping is exact at each level energy; away from it the two
energy dependences differ at second order in :math:`T_0/m_0`.

normalizations.out
^^^^^^^^^^^^^^^^^^

One complete row per data segment whose normalization is fitted (a varied
normalization, or a member of a ``normgroup``): segment key, data file, angle
bounds and frame, energy bounds and frame, the fitted normalization and the
energy shift. It is a report only; the fitted normalizations are carried by
``param.sav`` and nothing reads this file back. (Before October 2026 every
other segment also left a row of six fields without key, file or angles.)

Cross Section Output
--------------------

AZUREOut_aa=\*_R=\*.out
^^^^^^^^^^^^^^^^^^^^^^^^

Output from **Calculate With Data** and **Fit With Data** modes. The filename
encodes the entrance (``aa``) and exit (``R``) particle pair indices.

Nine columns:

.. list-table::
   :widths: 10 90
   :header-rows: 1

   * - Col.
     - Description
   * - 1
     - Center-of-mass energy (MeV)
   * - 2
     - Excitation energy (MeV): :math:`E_{c.m.} + S + E^*` of the entrance pair
       (for a beta-delayed segment, of the exit pair); before October 2026 the
       pair's :math:`E^*` was left out, so the column was low by it for an
       excited-state pair
   * - 3
     - Center-of-mass angle (degrees)
   * - 4
     - Fit center-of-mass cross section (barns or barns/sr), or the fitted
       analyzing power (dimensionless) for an analyzing-power segment
   * - 5
     - Fit center-of-mass S-factor (MeV b or MeV b/sr)
   * - 6
     - Data center-of-mass cross section (barns or barns/sr), or the measured
       analyzing power (dimensionless) for an analyzing-power segment
   * - 7
     - Data center-of-mass cross section uncertainty (barns or barns/sr)
   * - 8
     - Data center-of-mass S-factor (MeV b or MeV b/sr)
   * - 9
     - Data center-of-mass S-factor uncertainty (MeV b or MeV b/sr)

When multiple segments share the same entrance and exit particle pairs, their
data are written to the same file in the order they appear in the **Segments**
tab, separated by a double blank line.

AZUREOut_aa=\*_R=\*.extrap
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Output from **Calculate Segments Without Data** mode. Same naming convention
as above. Five columns:

.. list-table::
   :widths: 10 90
   :header-rows: 1

   * - Col.
     - Description
   * - 1
     - Center-of-mass energy (MeV)
   * - 2
     - Excitation energy (MeV): :math:`E_{c.m.} + S + E^*` of the entrance pair
       (for a beta-delayed segment, of the exit pair); before October 2026 the
       pair's :math:`E^*` was left out, so the column was low by it for an
       excited-state pair
   * - 3
     - Center-of-mass angle (degrees)
   * - 4
     - Extrapolated center-of-mass cross section (barns or barns/sr), or the
       analyzing power (dimensionless) for an analyzing-power segment
   * - 5
     - Extrapolated center-of-mass S-factor (MeV b or MeV b/sr)

AZUREOut_sample_seg=\*.out / .extrap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One file per natural-target segment (``sample 1``, see
:doc:`../user_guide/segments`): the centre-of-mass energy, the weighted total
the segment is compared with, and one column per compound nucleus holding
that nucleus's abundance-weighted contribution, so the columns add up to the
total. ``.out`` for a data segment, ``.extrap`` for a test segment.

Uncertainty and Statistics
--------------------------

chiSquared.out
^^^^^^^^^^^^^^

One line per data segment, then a total::

    Segment#, Chi-Squared,  N,  Norm,  Norm-Chi-Squared
    1,823.88,17,1,0
    ...
    Total-Chi-Squared: 107456 Total-Norm-Chi-Squared: 0 Total-N: 415

``Chi-Squared`` and ``Total-Chi-Squared`` are the **data** term only;
``Norm-Chi-Squared`` is the separate penalty on a varied normalization, and
``N`` counts data points (not degrees of freedom). The quantity a fit actually
minimises is the sum of both — see :doc:`../user_guide/chi_squared`.

This file is the quickest scalar check that a run succeeded.

With a ``<scatteringLengths>`` block the totals line ends with
``Total-Scattering-Length-Chi-Squared``, the constraint term of
:doc:`../theory/scattering_length` (also part of the minimised sum).
The same holds for ``Total-Park-Chi-Squared`` and
``Total-Abundance-Chi-Squared``; the file written through pyazr's
``write_output_files`` carries all three evaluated at the parameters it was
asked for (before October 2026 the API left them at zero or stale).

scatteringLengths.out
^^^^^^^^^^^^^^^^^^^^^

Written after every calculation and fit (the run log only names the file): the
s-wave scattering length of every channel spin of every particle pair at the
parameters of the run::

     Pair  Particles     J^pi           x      Re a0 (fm)      Im a0 (fm)  note
        1  p+7Be         1-        2.0438    1.706777e+01   -0.000000e+00
        1  p+7Be         2-        2.0438   -5.018396e+00   -0.000000e+00
        2  p+7Be         0-        2.0438    2.583014e+01   -0.000000e+00  no level in this channel (hard sphere)
        2  p+7Be         1-        2.0438    9.603285e+00   -2.637272e-01

``x`` is :math:`(8Z_1Z_2e^2\mu a/\hbar^2)^{1/2}` (0 for a neutral pair).  A
nonzero imaginary part means a channel with a lower threshold is open; a
channel spin in which no level has a width gets the pure Coulomb (hard-sphere)
value and the note.  The measured lengths of the ``<scatteringLengths>`` block
follow, each with its model value and :math:`\chi^2`, and their total.

covariance.out
^^^^^^^^^^^^^^

Written with ``chiSquared.out`` when data segments carry a ``covariance``
block (:doc:`data_formats`): one line per covariance file, with the number of
fitted rows, the rows of the file, the segments, :math:`\chi^2 = r^T M^{-1} r`
and, for comparison, :math:`\sum r_i^2` (the same points with their
correlations ignored), and the number of ``relvector`` terms.

covarianceModel.out
^^^^^^^^^^^^^^^^^^^

For covariance files with ``relvector`` terms: the model :math:`m` they were
built with, frozen at the first evaluation of the run (and appended again at
every refreeze from pyazr), one line per fitted row: the row of the
covariance file, the segment key, the lab energy and angle, the data value,
:math:`m` and the uncertainty used, all in the units of the data file. With
it the :math:`\chi^2` of the run can be reproduced.

param.errors
^^^^^^^^^^^^

Contains the reduced width amplitudes and their asymmetric uncertainties from a
MINOS error analysis.

covariance_matrix.out
^^^^^^^^^^^^^^^^^^^^^

Contains the covariance and correlation matrices from a MINOS calculation,
providing a complete description of parameter correlations.

Other Files
-----------

intEC.dat
^^^^^^^^^

External capture integral values for data segments. This file can be reused to
speed up subsequent calculations, as long as:

- No calculation segments have been added or removed.
- No levels of a new :math:`J^\pi` have been added or removed.
- No channels have been added or removed.
- The channel radius has not changed.
- The integration grid has not changed: target thickness, energy straggling,
  convolution width and the adaptive-grid settings all move the sub-point
  energies the integrals are evaluated at.
- The Coulomb-function routine (``--gsl-coul``) and the hybrid potential are
  the same.
- Relativistic kinematics (``relativistic=1`` in ``<config>``) are on or off as
  before: the switch changes every channel's eta, rho and flux factor (it was
  missing from the signature before October 2026; a file written with it off is
  no longer read back by a run with it on).

Energies, widths, and ANCs of the R-matrix levels can be changed freely while
reusing this file; the energy of an external-capture final state cannot, since
it enters the integrals.

AZURE2 checks this itself.  Beside the file it writes ``intEC.dat.sig``, a
64-bit hash of every energy an integral is evaluated at and of every input the
integrals depend on.  When a saved file is offered for reuse, both its number of
amplitudes and its signature must match the calculation at hand; otherwise a
``WARNING`` says why and the integrals are recomputed.  A file with no
signature -- written by an older AZURE2, or copied without its ``.sig`` -- is
recomputed once, with a warning, and gets one.

With the external-capture energy tables (the default, see
:doc:`../developer/ec_table`) computing the integrals afresh costs seconds even
with resolution or target integration, so reusing the file saves little.  The
file stores every digit (17 significant), so a reused file gives the same
:math:`\chi^2` as the fresh calculation that wrote it.  Before 30 October 2026
it stored six, which moved segment :math:`\chi^2` values by up to about 1e-5
relative (the ``tests/13N`` total by 0.62); such files carry an older signature
and are recomputed once.

intEC.extrap
^^^^^^^^^^^^

Same as ``intEC.dat``, but for the calculation (extrapolation) segments; its
signature is ``intEC.extrap.sig``.

reactionrates.out
^^^^^^^^^^^^^^^^^

Contains temperatures (in GK) and calculated reaction rates
(in cm\ :sup:`3` mol\ :sup:`-1` s\ :sup:`-1`) from the **Calculate Reaction
Rate** mode.

cross_section_grid.dat
^^^^^^^^^^^^^^^^^^^^^^

The energy grid the reaction rate integrated: centre-of-mass energy (MeV), the
angle-integrated cross section (b) and the S-factor (MeV b) at every node (48 per
decade from 0.1 eV to 50 MeV, plus a cluster around every level of the parameters
being evaluated). The rate is the integral of a spline of ln σ against ln E
through these points, so this file is the first thing to plot when a rate looks
wrong.

samples.mcmc
^^^^^^^^^^^^

The MCMC chain, as CSV, one row per walker per step::

    step,walker,logprob,loglikelihood,logprior,param0,param1,...

Every accepted state appears exactly once, including the repeated states a
rejected proposal contributes — that repetition is how a Markov chain carries
probability mass, so the file must not be deduplicated. ``logprob`` equals
``loglikelihood + logprior`` exactly. See :doc:`../user_guide/mcmc` for how to
load it and what to check before using it.

walkers.mcmc
^^^^^^^^^^^^

The final position of every walker, written at the end of an MCMC run and when
one is stopped early. Its purpose is resuming: with it, a continued run picks
the ensemble up where it left off instead of re-scattering the walkers and
splicing a fresh burn-in into the middle of the chain.
