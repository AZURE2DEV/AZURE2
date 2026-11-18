Fitting Tab
===========

.. note::

   This page describes the ``.azr`` content through the fields of the Qt setup
   utility of the released versions. This version has no setup utility;
   :doc:`graphical_interface` says how to do the same from the command line and
   with ``pyazr``.

The **Fitting** tab provides a centralized view and control over all parameters
that enter the R-matrix fit. It allows setting initial values, parameter limits,
and nuisance parameter flags.

Parameter Tables
----------------

Parameters are organized into sub-tabs, and one more, *THM Background*,
when THM is on (:doc:`configure_menu`, Runtime Options) and a THM experiment
has a coherent background (``cbackground=``, THM Workspace, Experiments
page): its free ``cbkg_*`` parameters, whose values edited here (or loaded
from a ``.sav``) go into the ``cbackground=`` key:

Level Parameters
^^^^^^^^^^^^^^^^

Displays one row for each level-channel combination with the following columns:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - **Name**
     - Identifier for the parameter (level energy or channel width), numbered
       as AZURE2 numbers its parameters: ``Level N`` is ``energy_N`` and
       ``Level N Channel M`` is ``width_N_M`` of ``param.sav``, counting the
       active levels in the order the project file lists them (by
       :math:`J^\pi`, then energy) and the channels within each
       :math:`J^\pi` group.
   * - **Value**
     - Current parameter value (editable).
   * - **Lower Limit**
     - Lower bound for the parameter during fitting (leave empty for no bound).
   * - **Upper Limit**
     - Upper bound for the parameter during fitting (leave empty for no bound).
   * - **Error**
     - Prior uncertainty for use as a nuisance parameter.
   * - **Fit Error**
     - Uncertainty from the fit result (read-only, populated after fitting).
   * - **Use as Nuisance**
     - Check to treat this parameter as a nuisance parameter with a Gaussian
       prior in the chi-squared calculation.

Abundance Parameters
^^^^^^^^^^^^^^^^^^^^

Shown only when the project carries an ``<abundances>`` block (a project
whose particle pairs form more than one compound nucleus, see
:doc:`particle_pairs`). One row per compound nucleus and sample, named
``abundance_<Z>_<A>`` (default sample) or ``abundance_<sample>_<Z>_<A>`` as
in ``param.sav``: **Value** is the isotopic abundance, **Error** the sigma of
its Gaussian prior (0 = none) and **Use as Nuisance** whether it is fitted at
all. The three are the same quantities as the Pairs tab's abundance table and
are written back there. In AZURE2's parameter vector the abundances sit
between the level parameters and the normalizations, sample by sample; one
that is not fitted takes no slot. The sum-to-one terms of the Pairs tab (one
per sample) are added by the engine and have no row here.

Normalization Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

Lists the normalization factors for each data segment that has **Vary Norm?**
enabled. These parameters scale the data yield during fitting.

The table has one more column, **Prior Centre**: empty, the prior on the
normalization is centred on its value in the Segments tab (the classic
behaviour); a number centres it there instead, whatever value the fit starts
from, so the Segments tab can hold the fitted normalization and the prior
stays on the experimental one (usually 1). It is saved as a ``prior_centre``
row of ``<parameterSettings>`` and survives *Load from .sav file*; see
:ref:`nominal-norm`.

A member of a normalization group (``normgroup <label>``, set in the
**Groups** field of the data segment dialog, :doc:`segments`) that is not the
group's first active segment has no normalization of its own: the engine uses
the first segment's. Its row is greyed and read-only, has no **Use as
Nuisance** check box, shows the first segment's value and reads, e.g.,
"tied to segment 1 (normgroup D)". To fit or fix the shared normalization,
change the first segment's row; to give the member its own, clear its group in
the Segments tab. The engine still carries a varied member's
parameter, ``<norm record>.value[<data set>]``, fixed at the shared value, so
parameter numbering does not depend on the groups.

Energy Shift Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

Lists the energy shift parameters for each data segment that has **Vary Energy
Shift?** enabled, with the same **Prior Centre** column. A segment's
:math:`\sqrt{E}` shift coefficient (the ``sqrtshift`` block in the input file,
see :doc:`segments`) is fitted by the engine under the name
``<data set>.shift.sqrt`` but is not listed here.

A member of a shift group (``shiftgroup <label>``) other than the group's
first active segment is shown tied to that segment in the same way as a
normalization-group member above ("tied to segment 1 (shiftgroup A)"). A
shift-function group (``shiftfngroup``) shares the function only, so its
members' constant shifts stay parameters of their own.

A shifted capture segment evaluates its external-capture amplitudes at the
shifted energies, from the same energy tables as every other point
(:doc:`../developer/ec_table`), so a shift gives the same cross sections as
data whose energies were moved by the same amount.  Before October 2026 (and
still with ``--no-ec-table``) the amplitudes of a shifted point were
interpolated linearly between the segment's own points and held constant
beyond its ends -- errors of several percent in the gaps between points and up
to orders of magnitude for low-energy points moved below the segment's lowest
energy.

Screening Parameters
^^^^^^^^^^^^^^^^^^^^

Shown only when an experimental-effect line carries electron screening (see
:doc:`experimental_effects`). One row per such line, named ``screening_<n>``
as in ``param.sav``, with :math:`n` the line's place among the *active*
experimental-effect lines: **Value** is the screening potential :math:`U_e` in
eV, **Error** the sigma of its Gaussian prior (0 = none) and **Use as
Nuisance** whether it is fitted. The three are the line's **U_e**, **Prior
Sigma** and **Fit U_e**, and an edit in either tab shows in the other. In
AZURE2's parameter vector the screening potentials come after the energy-shift
parameters; one that is not fitted takes no slot. **Load from .sav** sets the
value to the fitted :math:`U_e`, which is also the centre of its prior.

Target Density Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

Shown only when an experimental-effect line has a fitted target density (see
:doc:`experimental_effects`). One row per parameter, named ``density_<n>`` as in
``param.sav``: one per line with **Density Is a Fit Parameter** checked, a density
group counting once, under the number of its first line. **Value** is the density
in atoms/cm², **Error** the sigma of its prior in atoms/cm² (0 = none) and **Use as
Nuisance** whether it is fitted. They are the line's **Active Density**, **Prior
Sigma** and **Fit Density**, and an edit in either tab shows in the other; a new
value goes to every line of the group. A prior given on the line as a percentage
stays a percentage unless **Error** is edited here. In AZURE2's parameter vector
the densities come after the screening potentials; one that is not fitted takes no
slot. **Load from .sav** sets the density to the fitted value, which is also the
centre of its prior.

Doppler Temperature Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Shown only when an experimental-effect line has a fitted free-gas Doppler
temperature (``doppler T vary sigma [Tmax]``, see :doc:`experimental_effects`).
One row per such line, named ``doppler_<n>`` as in ``param.sav``, :math:`n` the
line's place among the active lines; a plain ``doppler T`` line has none.
**Value** is the temperature in K, **Error** the sigma of its prior in K (0 =
none) and **Use as Nuisance** whether it is fitted: the line's **Effective
Temperature T**, **Prior Sigma** and **Fit T**, synced both ways with the
Experimental Effects tab. In AZURE2's parameter vector the temperatures come after
the target densities; one that is not fitted takes no slot. **Load from .sav**
sets the temperature to the fitted value, which is also the centre of its prior.

Sample Thickness Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Shown only when a transmission or self-shielded-yield segment has the fit form
of its block (``transmission n vary sigma``, see :doc:`segments`). One row per
such active line, named ``<target>.areal_density`` as in ``param.sav``. **Value**
is :math:`n` in atoms/barn, **Error** the sigma of its prior (0 = none) and **Use
as Nuisance** whether it is fitted: the **Sample** row of the segment dialog
(n, Prior sigma, Fit), synced both ways with the Segments tab. A
``thicknessgroup`` follower's row is tied to its leader, as a normalization-group
follower's is: greyed, without a check box, showing the leader's value. As an
energy shift, a thickness parameter takes a slot in AZURE2's parameter vector
whether it is fitted or not. **Load from .sav** sets :math:`n` to the fitted value,
which is also the centre of its prior.

Prior Parameter Covariance
^^^^^^^^^^^^^^^^^^^^^^^^^^

The **Error** columns above are diagonal priors, one parameter at a time. A
multivariate Gaussian prior on several parameters at once -- an earlier
evaluation's parameters and their correlations, SAMMY's sequential update -- is a
file named under **Configure > Runtime Options** ("Prior Parameter Covariance",
``priorCovariance=<file>`` in ``<config>``); it has no table here. See
:doc:`chi_squared`, "Prior parameter covariance", for the file, the workflow with
``pyazr.prior_from_fit`` and how it combines with the priors of this tab (the file's prior replaces a parameter's
own), and :doc:`prior_builder` for the builder that writes the file from named parts (``pyazr.priors.build``, and
**Prior on the parameters…** in the Fit tab of the web page).

Controls
--------

Refresh from Current
   Rebuild the parameter table from the current GUI state (particle pairs,
   levels, segments).

Load from .sav File
   Import parameters from a ``param.sav`` file generated by a previous fit,
   fitted normalizations included (``normalizations.out`` is a report of them,
   not an input).

Populate Wigner Limits
   Automatically set upper and lower limits based on the Wigner limit
   (single-particle limit) for reduced width amplitudes. This provides
   physically motivated bounds for the fit.

Clear Limits
   Remove all upper and lower limit constraints from all parameters.

.. tip::

   Using the **Populate Wigner Limits** feature is a good starting point for
   constraining fits, especially for initial explorations of parameter space.
   Limits can then be adjusted manually as needed.

The Engine's Fitters
--------------------

The program and a pyazr session run the same fitters, from one table of the engine (``src/Fitters.cpp``;
``pyazr.kinds()["fitters"]``). The program takes the first row whose switch is given (the order below), MIGRAD
without one; a session names the row, ``Session.fit_with_engine(method=...)`` (:doc:`pyazr`).

.. list-table::
   :header-rows: 1
   :widths: 13 13 9 9 56

   * - Method
     - Switch
     - Call limit
     - Fallback
     - What it does
   * - ``nlopt``
     - ``--use-nlopt``
     - yes
     - (none)
     - NLopt's derivative-free Subplex (SBPLX) within the parameters' limits (unbounded: :math:`\pm10^{10}`); relative
       tolerances :math:`10^{-8}` in the objective and the parameters. No covariance, no MINOS. Only in builds
       configured with ``-DUSE_NLOPT=ON``; elsewhere the switch is refused (exit 2) with the switches of the build.
   * - ``lm``
     - ``--use-lm``
     - no
     - ``migrad``
     - Levenberg-Marquardt on the residual rows and their analytic Jacobian (priors as rows, limits by projection): at
       most 200 iterations, ending when the largest gradient component is
       below :math:`10^{-8}` or a step lowers the objective by less than :math:`10^{-9}` of it. Covariance
       :math:`(J^TJ + \text{priors})^{-1}`.
   * - ``migrad-gradient``
     - ``--use-gradient``
     - yes
     - (none)
     - Minuit2's MIGRAD with the analytic gradient of the objective and a Gauss-Newton second-derivative seed. HESSE
       covariance (the program runs HESSE with ``--covariance-band``), MINOS. The program also takes it for a project
       with THM data.
   * - ``migrad``
     - (none)
     - yes
     - (none)
     - Minuit2's MIGRAD with its own numerical gradient: the default. HESSE covariance, MINOS.

*Call limit*: the fitter honours ``--max-calls`` (``max_calls=`` of a session). A fitter that does not ends by its own
criterion and says so: a line ``Note: lm does not honour a call limit ...`` in the program's output, the same sentence
in the session result's ``notices``, and ``call_limit N (--max-calls) ignored`` in ``fitStatus.out``. *Fallback*: the
fitter that runs in its place for a model without an analytic Jacobian (no HESSE, no MINOS after it; its covariance is
MIGRAD's error matrix); the session result says so in ``fell_back_to_migrad`` and ``notices``, ``fitStatus.out`` in
``minimizer``. The program's usage text (``AZURE2 --help``) lists the switches of the fitters of the build from the
same table. ``gsl-lm`` (``--use-gsl-lm``, GSL's trust-region least squares on the same rows as ``lm``) was removed on
2026-11-18; its switch, and ``method="gsl-lm"`` of a session, are refused with the fitters to use instead (``lm``,
``migrad-gradient``), never run by another fitter in its place.

After a MIGRAD fit of a project whose data covariance has ``relvector`` terms, the model of those terms is re-frozen at
the fitted point and MIGRAD rerun up to ``--refreeze-covariance`` times. The program prints a progress line every
10 evaluations (``Iteration: N Chi-Squared: X``; LM one per iteration) and writes the intermediate files
every 100 evaluations (every 10 iterations for LM); a session receives the same steps as records
(``progress=``) and writes nothing.

How a Fit Ends
--------------

MIGRAD (the default minimizer) stops when the estimated distance to the minimum
(EDM) falls below its goal, :math:`2\times10^{-4}` times the error definition
(1 for a :math:`\chi^2` fit).  There is no limit on the number of function
calls unless the command-line option ``--max-calls=N`` sets one (see
:doc:`../reference/command_line`).  At the end of every fit the run says how
MIGRAD ended:

- ``MIGRAD converged: valid minimum, EDM ... (goal ...), N function calls``;
- ``WARNING: MIGRAD did NOT converge: stopped at the call limit ...`` or
  ``... invalid minimum ...`` (EDM above its goal, no valid covariance, HESSE
  failed).

The same is written to ``fitStatus.out`` in the output directory (see
:doc:`../reference/output_files`), by every fitter with the same keys.  A fit
that did not converge still writes ``param.sav``, ``chiSquared.out`` and the
other outputs, from the lowest objective the fit reached: a MIGRAD stopped at
its call limit ends at the lowest objective it evaluated, which can lie below
Minuit2's last state (Minuit2 returns its start when the limit falls while it
seeds the gradient); the run says so in a line, and the covariance is the one
of Minuit2's last state.  Rerun the fit from ``param.sav``.  Before 2026-10-31 MIGRAD and NLopt
stopped at a fixed 50,000 function calls without a message, which fits with
roughly a hundred or more free parameters could reach; their outputs looked
like a converged fit's.
