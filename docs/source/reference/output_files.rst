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
A parameter is named after the record it belongs to (:doc:`format3`). A sample's
isotopic abundances follow the level parameters as ``<sample>.isotope.<isotope>``
(:doc:`data_formats`). A target with electron screening adds its screening
potential (eV) as ``<target>.screening`` after the energy-shift parameters, a
fitted areal density ``<target>.areal_density`` (atoms/cm²) after those, and a
fitted Doppler temperature ``<target>.doppler`` (K) after the densities. A
transmission or yield data set whose sample thickness is fitted adds
``<target>.areal_density`` (atoms/barn) after the Doppler temperatures; a later
member of a thickness group is written as ``<target>.areal_density[<data set>]``
with its leader's value.

The data-set lines (``<data set>.norm``, ``<data set>.shift.offset``, ...)
belong to the data sets of the fit that wrote the file. A run without data -- an
extrapolation (**Calculate Segments Without Data**) or a reaction rate --
names its extrapolations' parameters after their own records and does not apply
the data sets' lines (the run log says so): an extrapolation's energies are the
ones its record asks for. ``pyazr``'s ``Session(file, extrapolations=True,
parameter_file=path)`` does the same.

``param.par`` and ``param.sav`` hold one row per parameter, ``name value
error``, after a first line ``#parametrization`` with the basis 0, 1 or 2
(standard, Brune, Park; :doc:`command_line`, "Brune and Park
parametrizations"), a comment line to positional readers such as
``numpy.loadtxt``. The rows are
matched **by name** when a file is read, by AZURE2 and by ``pyazr``: a name
the file lacks keeps its value from the project, and rows are appended as
parameters are added (``<data set>.shift.sqrt``, the THM ``cbkg_*``),
so read them by name in your own scripts too.

**Kinds of parameter.** A name is a label for people and files: AZURE2 does
not read from it what a parameter is. Each parameter carries its kind and what
it belongs to, set where the parameter is made; a name in a file is matched
whole. pyazr shows them as ``Parameter.kind`` and the owner fields
(``segment_key``, ``effect``, ``sample``, ``thm_experiment``, ...), and
``pyazr.parameters.parameter_kinds()`` lists the kinds. In the order of the
vector:

.. list-table::
   :header-rows: 1
   :widths: 14 30 18 38

   * - kind
     - name
     - unit
     - what it is
   * - ``energy``
     - ``energy_<n>``
     - MeV
     - level energy of level *n* (counted over all J groups)
   * - ``width``
     - ``width_<n>_<c>``
     - MeV\ :sup:`1/2`
     - reduced-width amplitude of level *n* in channel *c*, in the run's basis
   * - ``abundance``
     - ``<sample>.isotope.<isotope>``
     - 1
     - isotopic abundance (``<abundances>``)
   * - ``norm``
     - ``<data set>.norm``; a norm group's ``<norm record>.value``
     - 1
     - normalization of a data set
   * - ``shift``
     - ``<data set>.shift.offset``; a shift group's ``<energy_scale record>.offset``
     - MeV (lab)
     - constant energy shift of a data set
   * - ``shift_sqrt``
     - ``<data set>.shift.sqrt``
     - MeV\ :sup:`1/2`
     - coefficient of the :math:`\sqrt{E}` shift of a data set
   * - ``shift_fn``
     - ``<energy_scale record>.a<j>``
     - as the function defines
     - parameter :math:`a_j` of an energy scale's shift function
   * - ``screening``
     - ``<target>.screening``
     - eV
     - electron-screening potential of a target
   * - ``density``
     - ``<target>.areal_density``
     - atoms/cm²
     - fitted areal density of a target (or of its density group)
   * - ``doppler``
     - ``<target>.doppler``
     - K
     - fitted Doppler temperature of a target
   * - ``thickness``
     - ``<target>.areal_density``
     - atoms/barn
     - fitted sample thickness of a transmission or yield data set
   * - ``cbkg``
     - ``cbkg_<experiment>_...``
     - HOES amplitude (:math:`c_1` per MeV)
     - Re or Im of a THM coherent-background coefficient
   * - ``apparatus``
     - ``stopping_scale.<material>.<projectile>``, ``<record>.<field>.<parameter>``
     - the parameter's (1 for a stopping scale, MeV for a resolution width)
     - a number of an apparatus record a fit varies, last in the vector, one per sharing key: the stopping scale
       :math:`s` of one material and projectile (every target of that material crossed by that projectile has the
       stopping :math:`s\,\epsilon`; the target's ``stopping ... scale`` field), or a table kernel's own parameter
       written as a nuisance (``detector D`` with ``energy_resolution rectangular width 14 keV fit``:
       ``D.energy_resolution.width``, one for every data set naming ``D``).  Its prior adds
       ``Total-Stopping-Scale-Chi-Squared`` or ``Total-Kernel-Parameter-Chi-Squared`` to ``chiSquared.out``.
   * - ``radiative``
     - ``<record>.total_width``, ``<record>.branching.<final state>``, ``<record>.mixing.<final state>``
     - eV (total width), 1 (branching, mixing ratio)
     - a number of a level's ``gamma_decay`` record (:doc:`../user_guide/gamma_decay`): the total radiative width
       :math:`\Gamma_\gamma`, the branching :math:`b_f` to a final state, or its mixing ratio :math:`\delta_f`
       (Krane-Steffen), at the position of one of the level's photon amplitudes (one number per photon channel);
       the level's photon amplitudes are computed from them at every evaluation.  ``parameters.out`` ends with a
       block ``RADIATIVE DECAYS`` that gives every :math:`b_f` (the remainder too), :math:`\delta_f` and the
       partial widths they imply.

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
``N`` counts data points (not degrees of freedom). The ``Norm`` of a THM
segment is its profiled scale (shared by the segments of a THM experiment).
Under ``--use-park`` the total line ends with ``Total-Park-Chi-Squared``, the
:math:`J > 0` penalty. The quantity a fit actually
minimises is the sum of both — see :doc:`../user_guide/chi_squared`.
The totals line also carries the other prior terms when they apply (energy
shifts, abundances, screening, densities, Doppler temperatures, thicknesses,
scattering lengths); with a ``priorCovariance=<file>`` key in ``<config>`` it
ends with ``Total-Prior-Covariance-Chi-Squared: <value>``, the multivariate
prior :math:`(\mathbf{p}-\boldsymbol{\mu})^T C^{-1} (\mathbf{p}-\boldsymbol{\mu})`
at the run's parameters, written with 15 significant digits (the older fields
keep their 6).

This file is the quickest scalar check that a run succeeded.  After a fit,
check ``fitStatus.out`` too: ``chiSquared.out`` is written whether or not the
minimizer converged.

fitStatus.out
^^^^^^^^^^^^^

Written by every MIGRAD or NLopt fit (and removed at the start of a fit, so it
never describes an earlier one)::

    minimizer   MIGRAD
    status      converged
    valid       1
    edm         6.06159e-08
    edm_goal    0.0002
    calls       489
    call_limit  none

``status`` is ``converged``, ``call_limit`` (stopped at the ``--max-calls``
limit) or ``invalid`` (MIGRAD's minimum is not valid: EDM above its goal or no
valid covariance); NLopt writes ``converged_<reason>``, ``call_limit`` or
``failed_<reason>`` and no EDM lines.  ``calls`` counts objective evaluations as
the minimizer counts them.  The Levenberg-Marquardt and GSL minimizers do not
write this file unless they fall back to MIGRAD.

With a ``<scatteringLengths>`` block the totals line ends with
``Total-Scattering-Length-Chi-Squared``, the constraint term of
:doc:`../theory/scattering_length` (also part of the minimised sum); with a
``<modelData>`` section (:doc:`data_formats`) it also carries
``Total-Model-Data-Chi-Squared``, the sum over its data.  The two are the
parts of one term of the objective, "model data".
The same holds for ``Total-Park-Chi-Squared`` and
``Total-Abundance-Chi-Squared``; the file written through pyazr's
``write_outputs`` carries all three evaluated at the parameters it was
asked for (before October 2026 the API left them at zero or stale).

modelData.out
^^^^^^^^^^^^^

Written after every calculation and fit of a project with a ``<modelData>``
section (and copied to the run log): one line per datum with its name, its
quantity, the model value at the parameters of the run, the value and error
of the datum in the unit it is written in, the residual
:math:`(q - d)/e` and its :math:`\chi^2` (two lines, ``Re`` and ``Im``, for a
scattering length with an imaginary part), and the total, which equals
``Total-Model-Data-Chi-Squared``.

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

covariance.dat
^^^^^^^^^^^^^^

Written by a fit (MIGRAD, ``--use-lm``, ``--use-gsl-lm``, MINOS) when the
cross-section uncertainty band is on (``--covariance-band``, or "y" at the band
prompt), and read back by a later run's band. It is the covariance of the fit's
**free level energies and reduced-width amplitudes** only -- the parameters a
calculated cross section depends on; normalizations, energy shifts and the other
nuisance parameters are left out -- in parameter-vector order (``param.sav``'s
order with every fixed and every non-R-matrix parameter dropped) and in
``param.sav``'s space (Minuit's external parameters: MeV and the reduced-width
amplitudes of the run's parametrization). MIGRAD's is HESSE's matrix, a
least-squares fit's :math:`(J^TJ + \text{priors})^{-1}`. With
``--scale-covariance`` it is multiplied by the reduced χ² when that is above 1.
The matrix is written one row per line; since 2026-11-05 it is preceded by
``#columns <name> ...`` (the column parameters' names) and, when scaled,
``#scaled <factor>`` -- comment lines that ``numpy.loadtxt``, the engine's own
reader and older builds skip. ``pyazr.prior_from_fit`` turns it and
``param.sav`` into a prior covariance file (:doc:`../user_guide/chi_squared`).

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
- Relativistic kinematics (``relativistic=1`` in ``<config>``) are on or off as
  before: the switch changes every channel's eta, rho and flux factor (it was
  missing from the signature before October 2026; a file written with it off is
  no longer read back by a run with it on).
- The Coulomb table is the same, or is off as before (``--no-coulomb-table``):
  the integrals take their Coulomb waves from it.  The table is recognised by
  what it holds (its header and a hash of its values), not by its file name or
  date, so a copy of the same table elsewhere reuses the file and a rebuilt
  table with other settings does not (before November 2026 the table was not in
  the signature at all).

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
Rate** mode. With the exit pair "Total capture" (key ``-1``) the rate is the total
radiative capture rate, summed over every capture pair of the entrance nucleus.

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
