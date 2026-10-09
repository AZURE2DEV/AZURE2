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

``param.par`` and ``param.sav`` hold one row per parameter, ``name value
error``, after a first line ``#parametrization`` with the basis 0, 1 or 2
(standard, Brune, Park; :doc:`command_line`, "Brune and Park
parametrizations"), a comment line to positional readers such as
``numpy.loadtxt``. The rows are
matched **by name** when a file is read, by AZURE2 and by ``pyazr``: a name
the file lacks keeps its value from the project, and rows are appended as
parameters are added (``segment_N_energy_shift_sqrt``, the THM ``cbkg_*``),
so read them by name in your own scripts too.

param.fit
^^^^^^^^^

The parameters at the point being evaluated, written during a fit every 100
evaluations, in the format of ``param.sav`` (``#parametrization`` line
included): the state of an interrupted fit, which can be read back as a
parameter file in either mode. Files written before October 2026 lack the
line, and a run under ``--use-park`` takes their amplitudes as Brune's.

parameters.out
^^^^^^^^^^^^^^

Contains the **physical** (observable) parameters resulting from the fit. If
the user wishes to use these as new starting values, they must be manually
entered into the **Levels and Channels** tab.

normalizations.out
^^^^^^^^^^^^^^^^^^

Contains the fitted normalization factors for data segments where normalization
was varied. This file is automatically loaded when ``param.sav`` is selected.

shifts.out
^^^^^^^^^^

Written when a segment's energy shift or :math:`\sqrt{E}` shift coefficient
is varied: one row per such segment with its key, data file, angle and energy
bounds, norm, ``shift`` (MeV) and ``sqrt_shift`` (MeV\ :sup:`1/2`).

thm_experiments.out
^^^^^^^^^^^^^^^^^^^

Written by a calculation with data when the ``<thm>`` block defines THM
experiments (``experiment[<name>]`` lines): per experiment its segments,
background, number of points, chi-squared, status, the shared profiled norm
and the background coefficients with their uncertainties and covariance, and,
where the experiment has them, the coherent-background values (``cbkg``
rows), the line-shape :math:`\zeta` and :math:`E_{sF}`, the spectator window
(``ps_node``, ``<T_s>`` and ``ps_table`` rows), the distortion factor
(``distortion:`` and ``distortion_point`` rows) and the DW vertex
(``vertex:`` and ``dw_vertex_point`` rows). See
:doc:`../theory/thm_implementation`, "THM experiments" and the sections of
each option.

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
     - Excitation energy (MeV)
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
     - Excitation energy (MeV)
   * - 3
     - Center-of-mass angle (degrees)
   * - 4
     - Extrapolated center-of-mass cross section (barns or barns/sr), or the
       analyzing power (dimensionless) for an analyzing-power segment
   * - 5
     - Extrapolated center-of-mass S-factor (MeV b or MeV b/sr)

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
``N`` counts data points (not degrees of freedom). The ``Norm`` of a THM
segment is its profiled scale (shared by the segments of a THM experiment).
Under ``--use-park`` the total line ends with ``Total-Park-Chi-Squared``, the
:math:`J > 0` penalty. The quantity a fit actually
minimises is the sum of both — see :doc:`../user_guide/chi_squared`.

This file is the quickest scalar check that a run succeeded.

param.errors
^^^^^^^^^^^^

Contains the reduced width amplitudes and their asymmetric uncertainties from a
MINOS error analysis.

covariance_matrix.out
^^^^^^^^^^^^^^^^^^^^^

Contains the covariance and correlation matrices from a MINOS calculation,
providing a complete description of parameter correlations.

covariance.dat
^^^^^^^^^^^^^^

The parameter covariance a fit saves for the cross-section uncertainty band
(``--covariance-band``, or the band question of the CLI): the free level
energies, reduced widths and THM ``cbkg`` values, one row per free parameter.
A later calculation without data reads it for the band of an extrapolation;
``pyazr.bands.load_covariance`` reads it too.

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

intEC.extrap
^^^^^^^^^^^^

Same as ``intEC.dat``, but for the calculation (extrapolation) segments; its
signature is ``intEC.extrap.sig``.

reactionrates.out
^^^^^^^^^^^^^^^^^

Contains temperatures (in GK) and calculated reaction rates
(in cm\ :sup:`3` mol\ :sup:`-1` s\ :sup:`-1`) from the **Calculate Reaction
Rate** mode.

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
