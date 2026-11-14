The reaction rate utility
=========================

The thermonuclear rate :math:`N_A\langle\sigma v\rangle` of a project, its uncertainty band and the Gamow window,
computed by the engine: the program's ``AZURE2 rate`` command (:doc:`calculate`, "Calculate Reaction Rate";
:doc:`../reference/command_line`), the module ``pyazr.rate``, and the **Rate** tab of the web page.  All three run
the same computation (``src/ReactionRate.cpp``); pyazr and the page add the parameters of a fit or a sampler, the
band and the comparison with a table.

What the engine integrates
--------------------------

.. math::

   N_A\langle\sigma v\rangle = N_A\left(\frac{8}{\pi\mu}\right)^{1/2}(kT)^{-3/2}
   \int \sigma(E)\,E\,e^{-E/kT}\,dE

with :math:`E` the c.m. energy of the entrance pair, :math:`\mu` its reduced mass, :math:`k` =
8.617333262 × 10⁻² MeV/GK (CODATA 2018), in cm³ mol⁻¹ s⁻¹ with :math:`T_9` in GK.

- **Range**: :math:`10^{-7}` to 50 MeV.  :math:`\sigma` is zero above 50 MeV and below the lowest energy where the
  computed :math:`\sigma` is positive (an underflowed Coulomb tail, a closed exit below its threshold).
- **Grid**: :math:`\sigma` is computed once on 48 points per decade and, around every level of the parameters being
  evaluated (their observed energies and widths: the particle width, or the total width under the R-matrix
  capture formalism), the level's energy and 41 points on each side from 1/32 to 64 widths; a
  Steffen spline of :math:`\ln\sigma` against :math:`\ln E` is integrated at every temperature (GSL's QAGP,
  relative accuracy :math:`10^{-7}`, the grid's nodes as breakpoints).  A narrow resonance is resolved by its
  cluster, wherever the parameters put it.  The grid is written to ``cross_section_grid.dat``.
- **Fallback**: a grid with fewer than three points of :math:`\sigma > 0` is integrated by a quadrature per
  temperature (the run says so).
- Below the lowest level and above the highest the integrand is the R-matrix :math:`\sigma` of the same
  parameters (background poles included); nothing is extrapolated otherwise.
- **Exit -1** is the total radiative capture: the sum over every capture pair of the entrance's compound nucleus.
- **Identical entrance particles**: :math:`\sigma` carries :math:`(1+\delta_{12})` and the rate excludes the
  :math:`1/(1+\delta_{12})` of the rate equation (the compilations' convention).

In Python
---------

.. code-block:: python

   from pyazr import rate
   r = rate.reaction_rate("tests/upgrade/T1/T1.azr", 1, 2, grid="table")        # 12C(p,g)13N
   r.T9, r.rate, r.units                                       # GK, cm^3 mol^-1 s^-1
   r = rate.reaction_rate(file, 1, 2, t_min=0.01, t_max=1, t_step=0.07)       # the program's grid
   r = rate.reaction_rate(file, 1, 2, temperatures=[0.1, 0.3], parameter_file="output/param.sav")

``reaction_rate(file, entrance, exit, values=None, *, grid | temperatures | t_min, t_max, t_step,
parameter_file=None, options=None)`` returns a ``Rate`` (``T9``, ``rate``, ``reaction``, ``source``, ``units``).
The temperatures: a named grid of ``TEMPERATURE_GRIDS`` (``table``: the 51 temperatures of the rate tables from
0.01 to 10 GK; ``log20``: 20 per decade), a list, or the program's grid :math:`T_{min} + i\,\Delta T`.  ``values``
are the parameters in the fit basis, by name (``dict(zip(names, x))`` of a session's fit vector; names the rate has
not, normalizations and shifts, are left out) or as a vector of every rate parameter.

``RateSession(file, entrance, exit)`` keeps the session for many rates (``.rate(values, ...)``, ``.names``,
``.values``); ``rate_band(session, values, covariance=C, names=...)`` propagates a covariance linearly in
:math:`\ln N_A\langle\sigma v\rangle` (central differences of the engine's rate, the steps of
``pyazr.bands.step_sizes``; the band is the rate times :math:`e^{\mp\sigma}`), ``rate_band(session, values,
draws=D, names=...)`` gives the 16th, 50th and 84th percentiles over the rates of drawn vectors (a sampler's chain
after its burn-in), each a ``RateBand`` (``low``, ``median``, ``high``, ``best``, ``how``, ``varied``).
``gamow_window(file, entrance, T9)`` gives the Gamow peak :math:`E_0 = (b\,kT/2)^{2/3}` and its width
:math:`\Delta = 4(E_0 kT/3)^{1/2}`, :math:`b = 2\pi\eta(E)E^{1/2}` from the engine's Sommerfeld parameter: the
window :math:`E_0 \mp \Delta/2` is where the data constrain a non-resonant charged-particle rate (NaN for a neutral
entrance).  ``write_rate_table(path, rate, band=None, comparison=None)`` writes the CSV with a header that names the
project, the parameters' source and the units; ``read_rate_table(path)`` reads a table to compare with (two or
three columns: :math:`T_9`, rate, optionally its uncertainty).  ``pyazr/examples/reaction_rate.py`` runs them.

A worked case
-------------

T1 (``tests/upgrade/T1``), :math:`^{12}`\ C(p,γ\ :sub:`0`), at :math:`T_9` = 0.1: :math:`N_A\langle\sigma v\rangle`
= 5.59562 × 10⁻⁷ cm³ mol⁻¹ s⁻¹, the Gamow peak at 84.7 keV, the window 53.5 to 115.9 keV.  pyazr's rates of T1 and
T3, and of the models of ``tests/reaction_rate`` on a grid, a list and with a parameter file, are the program's to
the six digits ``reactionrates.out`` holds (``tests/pyazr/reaction_rate_test.py``).

The rate at a fit's parameters is computed by a session set up for the rate, as the program's.  The parameters a
data session holds for the same ``.azr`` (its ``x0``) can differ from the rate session's own reading of the
``.azr`` in the last digits of the capture widths' amplitudes (1e-10 relative on T1), which moves the rate by
about 1e-9; a rate at a fit's or a sampler's point uses that point's values.

On the page
-----------

The **Rate** tab: the reaction (an entrance pair and an exit pair, or -1 for the total capture), the temperatures (a
named grid, or min, max and step in GK), the parameters (the working x, the project's values, or a finished fit or
sampler job's best point), the band (none; from that fit's covariance; from that sampler's draws, with their
number), and a comparison table of the project folder.  **Calculate** starts a job (kind ``rate``); its progress is
the number of rates computed.  At its end: the rate against :math:`T_9` on log-log axes with the band shaded and the
comparison dashed, the ratios to the median below, the Gamow window per temperature, and the table as CSV (also
written to ``output/rate_<job>.csv``), and the STARLIB and REACLIB files (below).  The page computes no rate, ratio or
window: they are the job's.

For rate libraries
------------------

``pyazr.rate_formats`` writes the rate in the formats of two rate libraries; the layouts are those of the formats'
own descriptions: REACLIB's from the JINA REACLIB Database's help page "REACLIB File Format" (the chapters of its
Table 1, the REACLIB 1 file: a header entry, then set entries of three lines of 74 characters,
``format(i1,4x,6a5,8x,a4,a1,a1,3x,1pe12.5)``, ``format(4e13.6)``, ``format(3e13.6)``), the set's rate
:math:`\exp(a_0 + \sum_{i=1}^{5} a_i T_9^{(2i-5)/3} + a_6\ln T_9)` of Cyburt et al. 2010 (ApJS 189, 240, eq. 1);
STARLIB's from Sallaska et al. 2013 (ApJS 207, 18), Appendix A "STARLIB Format" (a header line of 64 fields, then
:math:`T_9`, the recommended rate and its factor uncertainty on 60 temperatures from 1 MK to 10 GK).

- The species, the chapter and the Q value come from the project's pairs (``species_of(file, entrance, exit)``):
  :math:`Q` is the entrance pair's threshold minus the exit pair's (separation energies relative to the compound
  nucleus, plus excitation energies); exit -1 ends in the compound nucleus's ground state.  The labels and flags are
  the user's, among the formats' allowed values.
- **STARLIB** (``write_starlib(path, rate, band, species, label, flag)``): the rate must be computed on
  ``grid="starlib"`` (``STARLIB_T9``: the paper does not print the 60 temperatures; they are 0.001 to 0.009 GK and
  the 51 temperatures of the paper's rate tables from 0.010 to 10 GK) and with a band.  From a covariance band the
  recommended rate is the rate at the parameters and the factor uncertainty :math:`e^{\sigma}` of
  :math:`\ln N_A\langle\sigma v\rangle`; from draws the recommended rate is the median and the factor uncertainty
  :math:`(x_{84}/x_{16})^{1/2}` (the paper's :math:`x_{low} = e^{\mu-\sigma}`, :math:`x_{high} = e^{\mu+\sigma}`).
  Without a band nothing is written: no uncertainty is made up.  The paper gives no number format: the columns are
  written ``%.3e %.6e %.4e``.
- **REACLIB** (``fit_reaclib(rate, tolerance=0.05, max_sets=3)``, ``write_reaclib(path, fit, species, label,
  flag)``): one set is fitted to :math:`\ln N_A\langle\sigma v\rangle` by linear least squares, then two or three by
  nonlinear least squares of the logarithm of their sum (scipy), until the largest relative deviation over the
  computed temperatures in 0.01 to 10 GK is within the tolerance; the fit says how many sets it needed and its
  largest deviation, and is refused, with what was wrong, what is allowed and how to fix it, when the tolerance is
  not reached, when the rate is not computed over 0.01 to 10 GK (nothing is extrapolated) or when a rate there is
  not positive.  The numbers are written as ``%13.6e``.
- ``read_starlib(path)`` and ``read_reaclib(path)`` (with ``reaclib_rate(sets, T9)``) read a library's entry, e.g.
  to compare it with this rate (``read_rate_table`` takes a two- or three-column table).

T1's :math:`^{12}`\ C(p,γ\ :sub:`0`) on the table grid needs two sets (largest deviation 0.43 % at
:math:`T_9` = 1.25); a rate made of one known set is recovered to :math:`10^{-13}`
(``tests/pyazr/rate_formats_test.py``).
The Rate tab offers both files beside the CSV, with the labels and flags as form controls; a refusal is shown in its
words.
