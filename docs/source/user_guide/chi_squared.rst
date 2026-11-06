How the Chi-Squared Is Built
============================

Everything AZURE2 fits or samples is driven by one objective function. It is
worth knowing exactly what goes into it, because the number reported in
``chiSquared.out`` is *not* the whole objective, and the difference matters as
soon as normalizations are free.

The data term
-------------

For every point in every active segment, AZURE2 compares the calculated cross
section against the measured one, with the segment's normalization :math:`n`
applied to the *data* rather than to the theory:

.. math::

   \chi^2_{\text{data}} = \sum_{s}\ \sum_{i \in s}
   \frac{\bigl(\sigma^{\text{fit}}_i - n_s\,\sigma^{\text{data}}_i\bigr)^2}
        {\bigl(n_s\,\delta\sigma^{\text{data}}_i\bigr)^2}

Points with zero uncertainty are skipped rather than producing an infinity.
Because both the value and its error are scaled by :math:`n_s`, a normalization
that is free to move cannot reduce :math:`\chi^2` simply by shrinking the
uncertainties.

Correlated data: covariance matrices
------------------------------------

Data whose uncertainties are correlated (a run-by-run fit covariance, a
detector efficiency common to an excitation function, an analysis systematic
that moves several observables together) can carry their covariance matrix
:math:`C`: segments whose lines end with ``covariance <file>``
(:doc:`segments`; file format in :doc:`../reference/data_formats`) share it, and
their points enter as

.. math::

   \chi^2_{\text{cov}} = r^{T} M^{-1} r, \qquad
   r_i = \frac{\sigma^{\text{fit}}_i - n_s\,\sigma^{\text{data}}_i}{n_s\,\delta_i},
   \qquad M = D^{-1} C D^{-1},

with :math:`r` the same normalized residuals as above and
:math:`D = \mathrm{diag}(\delta_i)`. This is :math:`v^T C^{-1} v` with
:math:`v` the model-minus-data differences in the units of the data file, so
the value does not depend on :math:`\delta_i` (AZURE2 sets it to
:math:`\sqrt{C_{ii}}`), the matrix scales with the normalization exactly as the
error column does, and the lab-to-c.m. conversion acts on it as on the error
column. A diagonal :math:`C = \mathrm{diag}(\delta_i^2)` reproduces the
uncorrelated sum exactly. The covariance covers only the rows actually fitted:
points removed by the energy or angle limits drop out of the matrix (the
marginal distribution). Normalization and energy-shift penalties, normalization
groups and every fitter (MIGRAD with or without ``--use-gradient``,
``--use-lm``, ``--use-gsl-lm``, the MCMC likelihood) work unchanged; the
least-squares fitters see whitened residuals :math:`w = W r` with
:math:`W^T W = M^{-1}`, so :math:`\sum w_i^2 = \chi^2_{\text{cov}}`.

**Multiplicative uncertainties and Peelle's pertinent puzzle.** A fully
correlated *relative* uncertainty (an efficiency, a solid angle, an overall
scale) entered as :math:`(f\,d)(f\,d)^T`, built from the measured values
:math:`d`, biases a least-squares fit low (Peelle's pertinent puzzle). Two
remedies: fit the scale as a normalization with its penalty (the usual AZURE2
way, ``normgroup`` for one scale shared by several segments), or enter it as a
``relvector`` term, :math:`(f\,m)(f\,m)^T` with :math:`m` the **model** in the
data's units. AZURE2 freezes :math:`m` at the first evaluation of the run (the
starting parameters), so :math:`C` stays constant through a minimization, and
writes it to ``covarianceModel.out``. The iterative scheme (Chiba and Smith)
refits from the result: run AZURE2 again from the fitted ``param.sav``, or in
pyazr call ``refreeze_covariance()`` between fits.

``chiSquared.out`` gives a covariance segment the share
:math:`\sum_{i \in s} r_i (M^{-1} r)_i`; the shares add up to the group's
:math:`\chi^2` exactly, but a single share is not a goodness of fit of that
segment alone (it can even be negative). ``covariance.out`` lists each
covariance file's :math:`\chi^2` next to :math:`\sum r_i^2`, what the same
points would give with their correlations ignored.

The normalization penalty
-------------------------

A free normalization is not free of consequence. Each segment that has
**Vary Norm?** enabled and a non-zero **Norm Error** contributes

.. math::

   \chi^2_{\text{norm},s} =
   \left(\frac{n_s - n_s^{\text{nom}}}{n_s^{\text{nom}} \cdot \epsilon_s / 100}\right)^2

where :math:`n_s^{\text{nom}}` is the nominal normalization and
:math:`\epsilon_s` is the quoted systematic uncertainty **in percent** — the
value typed into the *Norm Error* column of the Segments tab. This is what
keeps a dataset's normalization near its experimental value instead of letting
it absorb every discrepancy in the model.

A segment with **Vary Norm?** enabled but no quoted error has no penalty, and
its normalization is genuinely unconstrained.

.. _nominal-norm:

The nominal normalization: the column, or an explicit prior centre
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A segment line has one normalization field, and AZURE2 reads it as
the *start value* of the fit (or the value a calculation uses) and, in a
classic file, also as the *centre* :math:`n_s^{\text{nom}}` of the penalty
above (which also scales its width, :math:`n_s^{\text{nom}}\epsilon_s/100`).
A fit leaves the column alone and writes the fitted normalizations to
``param.sav`` and ``normalizations.out``, so the prior stays where it was put.
Anything that copies fitted values into the column moves the prior with
them -- unless the centre is held explicitly. A row

.. code-block:: text

   <parameterSettings>
   ...
   segment_3_norm prior_centre 1
   segment_3_energy_shift prior_centre 0
   </parameterSettings>

centres the prior of segment 3's normalization on 1 (and of its energy shift
on 0), whatever the Data Norm. (Energy Shift) column holds; the column is then
only the start value. (The energy-shift prior is centred on 0 anyway since
25 October 2026, see below; on this branch an explicit
``segment_N_energy_shift prior_centre c`` row moves that centre to *c*.) ``3`` is the segment key, the line's position in
``<segmentsData>`` counting inactive lines, as in ``segment_3_norm`` of
``param.sav``. A normalization centre must be positive; a row naming no
active segment is reported and ignored. Without a row nothing changes: a file
that has none behaves exactly as before, and AZURE2 versions without the
feature skip these three-token rows (the other ``<parameterSettings>`` rows
have 9, 12 or 14) and centre the prior on the column.

With an explicit centre a file can carry the fitted normalizations *and* the
experimental priors, so that run on its own it is the best fit (data and
priors) and a fit started from it is pulled to the experimental values, not
to its own previous result:

- the Fitting tab shows the centre in the **Prior Centre** column of the
  Normalization and Energy Shifts tables (empty: the column value for a
  normalization, 0 for an energy shift) and keeps
  it when it loads a ``param.sav``; loading one still writes the fitted
  normalizations into the column, so in a project without explicit centres a
  fit started from the saved project is centred on them -- type the
  experimental value (usually 1) into **Prior Centre** first;
- pyazr's ``save_fit`` (default ``norms="fitted"``) writes the fitted norms
  and shifts into the columns and, for every normalization prior that thereby
  leaves its column, the centre the fit used as an explicit row (for a shift,
  only an explicit non-zero centre); the snapshot run on its own reproduces
  the fitted data :math:`\chi^2` and priors.
  ``save_fit(..., norms="nominal")`` keeps the columns and puts the fitted
  values only in the ``.sav`` (to be given as the parameter file: Run tab,
  *Parameters File*; CLI "use previous parameters").

The examples are stored this way: every direct normalization with a quoted
error (``n15_lacognata2007``, ``li7_tumino2006``, ``li6_pizzone2011``,
``o18_lacognata2010``, ``f19_pag_thm``, ``c12c12_tumino2018``) holds its
fitted value in the column and ``prior_centre 1``; the file alone gives the
best-fit total (the "Total Chi-Squared" on the CLI, data plus priors), listed
in the THM documentation ("Examples"). A fixed normalization pays its
penalty too when it sits away from an explicit centre (``f19_pag_thm`` and
``c12c12_tumino2018`` store fitted norms fixed), as it would after a
``param.sav`` moved it. A shift's penalty is measured from 0 (or from its
explicit centre), so writing fitted shifts into the column does not move it.

Segments that share one normalization (``normgroup``, :doc:`segments`)
contribute this term once, with the first member's :math:`n^{\text{nom}}` and
:math:`\epsilon`; the other members carry no normalization penalty.

The energy-shift penalty
------------------------

Segments with **Vary Energy Shift?** enabled contribute the analogous term,
this time with an absolute rather than a percentage error:

.. math::

   \chi^2_{\text{shift},s} =
   \left(\frac{\Delta E_s}{\delta(\Delta E_s)}\right)^2

measured from **zero**, i.e. from no shift: the value of :math:`\Delta E_s` in
the input file is only where the fit starts. (On this branch an explicit
``segment_N_energy_shift prior_centre c`` row, above, measures the constant
term from *c* instead.) A segment whose :math:`\sqrt{E}`
shift coefficient :math:`b_s` is free (the ``sqrtshift`` block, see
:doc:`segments`) adds the same kind of term, :math:`(b_s/\delta b_s)^2`, with
:math:`\delta b_s` in MeV\ :sup:`1/2`, and so does each free parameter
:math:`a_j` of an energy-shift function (the ``shiftfn`` block), once, on the
segment that declares the function: :math:`(a_j/\delta a_j)^2`. A free shift
term with no uncertainty has no penalty.

.. note::

   Before 25 October 2026 these terms were measured from the value in the input
   file, :math:`((\Delta E_s - \Delta E_s^{\text{file}})/\delta(\Delta E_s))^2`
   and likewise for :math:`b_s` and :math:`a_j`. A file whose shifts were seeded
   from an earlier fit therefore centred the penalty on that fit's shifts. Such
   files now give a larger :math:`\chi^2` (the shift part changes; the data part
   does not). The normalization penalty is unchanged.

Nuisance parameters
-------------------

Level energies and widths flagged **Use as Nuisance** in the Fitting tab add a
Gaussian term of the same shape, using the **Error** column as the width. This
is how prior knowledge of a resonance energy enters a least-squares fit. The
term is evaluated on the flagged parameter itself; before 2026-10-27 it was
evaluated on the wrong parameter whenever a fixed level parameter preceded the
flagged one in the parameter list (a fixed ground-state energy was enough),
so MIGRAD pulled another parameter toward the nominal value and a parameter
sitting at its nominal value still paid a penalty. The least-squares drivers
were not affected.

THM segments
------------

A THM (half-off-energy-shell) segment with a free normalization has an
arbitrary scale. It is not a fit parameter and pays no penalty: at every
evaluation the scale is profiled out in closed form,
:math:`n^* = S_{mm}/S_{md}` and
:math:`\chi^2 = S_{dd} - S_{md}^2/S_{mm}` with the weighted sums of model
:math:`\times` model, model :math:`\times` data and data :math:`\times`
data, and the segments of one THM experiment share it (with an optional
smooth background, profiled with it). The coherent THM background
(``cbkg_*``) is an ordinary parameter without a penalty. See
:doc:`../theory/thm_implementation`, "THM experiments" and "Normalization,
gradients and uncertainty bands".

Park's parametrization
----------------------

Under ``--use-park`` a level whose overlap :math:`J_\lambda` is negative adds
:math:`(J_\lambda/10^{-3})^2` to the objective (zero inside the physical
region), in every minimizer's cost and in the gradient. A calculation reports
it as ``Total-Park-Chi-Squared`` in ``chiSquared.out``. See
:doc:`../reference/command_line`, "Brune and Park parametrizations".

What ``chiSquared.out`` actually reports
----------------------------------------

The file separates the two contributions::

    Segment#, Chi-Squared,  N,  Norm,  Norm-Chi-Squared
    1,823.88,17,1,0
    2,24.335,20,1,0
    ...
    Total-Chi-Squared: 107456 Total-Norm-Chi-Squared: 0 Total-N: 415

- **Chi-Squared** (per segment) and **Total-Chi-Squared** are the *data* term
  only.
- **Norm-Chi-Squared** and **Total-Norm-Chi-Squared** are the normalization
  penalty.
- **N** is the number of points in the segment, and **Total-N** their sum.
- **Total-Shift-Chi-Squared**, at the end of the totals line whenever an
  energy shift is a fit parameter, is the sum of the energy-shift penalties
  (constant, :math:`\sqrt{E}` and shift-function terms).  ``shifts.out``
  breaks it down per shift function.  Files written before October 2026 lack
  it, and their totals fall short of the minimizer's own value by exactly this
  amount.
- **Total-Screening-Chi-Squared**, when an electron-screening potential has a
  prior (:doc:`experimental_effects`), is the sum of
  :math:`((U_e - U_e^{0})/\sigma)^2`.
- **Total-Density-Chi-Squared**, when a fitted target density has a prior
  (:doc:`experimental_effects`), is the sum of :math:`((N - N_0)/\sigma)^2`.
- **Total-Doppler-Chi-Squared**, when a fitted Doppler temperature has a prior
  (``doppler T vary sigma``, :doc:`experimental_effects`), is the sum of
  :math:`((T - T_0)/\sigma)^2`.
- **Total-Thickness-Chi-Squared**, when a fitted sample thickness has a prior
  (``transmission n vary sigma`` or ``yield n vary sigma``,
  :doc:`segments`), is the sum of :math:`((n - n_0)/\sigma)^2` (a
  ``thicknessgroup`` counts once, with its first segment's prior).
- **Total-Scattering-Length-Chi-Squared**, when the project has a
  ``<scatteringLengths>`` block, is the sum over its active lines of
  :math:`((a_0 - a_{\mathrm{exp}})/\sigma)^2` (real and, if given, imaginary
  part), with :math:`a_0` the s-wave scattering length at the pair's threshold
  (:doc:`../theory/scattering_length`).  ``scatteringLengths.out`` lists each
  line's model value and contribution.  The term is part of every fitter's
  objective and of the MCMC likelihood.
- **Total-Prior-Covariance-Chi-Squared**, whenever the ``<config>`` block names
  a prior covariance file (``priorCovariance=<file>``, see below), is
  :math:`(\mathbf{p}-\boldsymbol{\mu})^T C^{-1} (\mathbf{p}-\boldsymbol{\mu})`,
  written with 15 significant digits.

The quantity the minimizer actually descends is the sum of all of them, plus
the energy-shift and nuisance terms, which the file does not list. Quoting
``Total-Chi-Squared`` as "the χ²" of a fit with free normalizations understates
the objective. A calculation prints the whole objective, data and every prior,
on the console as ``Total Chi-Squared:`` (twelve significant digits); this is
the number the examples quote as "file alone"
(:doc:`../getting_started/examples`). Under ``--use-park`` the total line of
the file ends with ``Total-Park-Chi-Squared``.

.. note::

   ``Total-N`` counts data points, not degrees of freedom. Subtract the number
   of free parameters yourself before forming a reduced χ².

Prior parameter covariance
--------------------------

The priors above are diagonal: each acts on one parameter. A *prior parameter
covariance* is a multivariate Gaussian on a set of parameters, so a fit can start
from an earlier evaluation's parameters **and their correlations** -- the
sequential Bayesian update of SAMMY's "PUP" option. Fit A of data set 1 gives
parameters :math:`\boldsymbol{\mu}` and a covariance :math:`C`; fit B of data
set 2 then minimises

.. math::

   \chi^2_B = \chi^2_{\text{data 2}} + \ldots +
   (\mathbf{p}-\boldsymbol{\mu})^T C^{-1} (\mathbf{p}-\boldsymbol{\mu}),

which (for Gaussian errors and a model linear over the uncertainties) has the
minimum and the covariance of a joint fit of both data sets, without refitting
data set 1. The term acts on the parameters as the fit's parameter vector holds
them, the values ``param.sav`` lists: level energies in MeV, reduced-width
amplitudes of the run's parametrization (Brune by default; the
``#parametrization`` tag of ``param.sav`` says which), normalizations, energy
shifts and the other nuisance parameters in their own units.

**Switching it on.** One line in the ``<config>`` block
(:doc:`../reference/data_formats`), with a path relative to the ``.azr``'s
directory (an absolute path also works)::

    priorCovariance=prior_from_A.txt

The GUI has the field under **Configure > Runtime Options** ("Prior Parameter
Covariance", with a **Choose...** button); it is saved in the project, so the
command line and pyazr use it too. Removing the line switches it off.

**The file.** Plain text; ``#`` starts a comment. An optional
``#parametrization <0|1|2>`` line (``param.sav``'s tag: 0 standard, 1 Brune,
2 Park) before the names, then the parameter **names** on one line (as in
``param.sav``), their **means** on the next, and the :math:`n` rows of the
symmetric :math:`n \times n` **covariance** in the same order::

    #parametrization 1
    # widths of the 2- and 3- levels and the normalization of segment 4, from fit A
    width_2_1 width_3_1 segment_4_norm
    0.8206 0.8080 1.05
    5.0e-5 -2.1e-5 0
    -2.1e-5 2.1e-5 0
    0 0 4.0e-4

AZURE2 refuses the file, with an ``ERROR`` and before anything is written, when a
name is not a parameter of the model, the matrix is not symmetric (to
:math:`10^{-10}` of :math:`\sqrt{C_{ii}C_{jj}}`) or not positive definite (its
Cholesky factorization fails), the counts disagree (names, means, rows, entries
per row), or the tag names another parametrization than the run's: a Brune
prior's amplitudes are not Park's, and nothing is converted -- rerun the fit the
prior comes from in this run's parametrization (an untagged file is taken as it
is).

**Fixed parameters are marginalized.** A listed parameter that is fixed in the
run has its row and column dropped, with a Note naming it. That is the correct
prior on the others: the marginal of a multivariate Gaussian is the Gaussian of
the remaining block.

**Overlap with the diagonal priors.** A listed parameter can also carry a prior
of its own -- a normalization with an error on its segment line, a free energy
shift with an error, a screening, density, Doppler or thickness sigma, an
abundance sigma, a nuisance flag. A fit applies **both**, and the run prints a
Note listing such parameters. If the covariance already contains that
information (a normalization fitted in A with its own penalty, say), it is
counted twice: set the segment's error to 0 in B, or leave the parameter out of
the prior file.

**How it is computed.** The matrix is factored once, :math:`C = L L^T`; every
evaluation forms the whitened residuals :math:`\mathbf{r} = L^{-1}
(\mathbf{p}-\boldsymbol{\mu})` and adds :math:`\mathbf{r}\cdot\mathbf{r}`.
MIGRAD's analytic gradient adds :math:`2 C^{-1}(\mathbf{p}-\boldsymbol{\mu})`;
``--use-lm``, ``--use-gsl-lm`` and the Jacobian-path gradient append the
:math:`n` rows :math:`\mathbf{r}` with the constant Jacobian :math:`L^{-1}` to the
data rows, so the fitted covariance (``covariance.dat``, parameter errors)
includes the prior's information. Parameter limits act as they do for every
other term (the prior is evaluated on Minuit's external values). In the MCMC the
multivariate Gaussian, :math:`-\tfrac12\mathbf{r}\cdot\mathbf{r} -
\tfrac12\ln\det(2\pi C)`, is part of the log prior and **replaces** the diagonal
(automatic or user) prior of each of its parameters (a Note lists them); a chain
over physical parameters refuses a prior on a level energy or width, since the
prior is on the reduced-width amplitudes.

**Workflow: fit A, then B.**

1. Fit A with the uncertainty band on (``--covariance-band``, or "y" at the
   band prompt): the fit writes ``output/param.sav`` and ``output/covariance.dat``.
   ``covariance.dat`` holds the covariance of the **free level energies and
   reduced-width amplitudes** only (normalizations, shifts and the other nuisance
   parameters are not in it), in ``param.sav`` order and space; since
   2026-11-05 its first line ``#columns <names>`` names the columns. A fit with
   ``--scale-covariance`` (or "y" at the scaling prompt) multiplied it by the
   reduced χ² when that exceeded 1 (a ``#scaled`` line says by how much).
2. Write the prior file with pyazr (no engine needed)::

       from pyazr import prior_from_fit
       prior_from_fit("A/output/param.sav", "A/output/covariance.dat",
                      "B/prior_from_A.txt")                    # all of A's free level parameters
       prior_from_fit(..., names=["width_2_1", "width_3_1"])  # or a subset: the marginal block

   The file carries ``param.sav``'s ``#parametrization`` tag. For a
   ``covariance.dat`` written before 2026-11-05 (no ``#columns`` line), pass
   ``columns=[...]``, the free ``energy_*`` / ``width_*`` names of fit A in
   ``param.sav`` order -- ``param.sav`` alone does not say which were free.
3. In B's ``<config>``: ``priorCovariance=prior_from_A.txt``, and fit. B's level
   scheme must use the same parameter names (the same levels and channels in the
   same order); a name B does not have is refused. With B's data identical to
   A's, B returns A's parameters, the same data χ² and a prior term of ~0
   (``tests/prior_covariance``).

A prior from a fit of *other* data is the point of the method; a prior built from
the very data fit B also uses counts those data twice.

Consistency with the MCMC
-------------------------

The Bayesian sampler reaches the same objective by a different route, which is
worth understanding if you compare a fit against a posterior.

Its log-likelihood is built from the **data term alone**:

.. math::

   \ln \mathcal{L} = -\tfrac{1}{2}\,\chi^2_{\text{data}}

The normalization and energy-shift penalties reappear as *priors* instead:
AZURE2 derives a Gaussian prior for every varying normalization and energy
shift straight from the quoted experimental errors, with exactly the mean and
width implied by the penalties above. The resulting log-posterior therefore
matches the fit objective term for term, and a posterior mode coincides with a
χ² minimum. See :doc:`mcmc`.

The same distinction applies to :doc:`pyazr`: ``calculate_chi2_rwa`` and
``residual_jacobian`` return the **data** χ² only (``residual_jacobian`` plus the
rows of a prior parameter covariance, when the project has one). A least-squares fit driven
from Python that minimises those residuals alone will let the normalizations
drift to absorb every discrepancy, reaching a "better" χ² that AZURE2 itself
would never find. Append the penalty rows explicitly — the recipe is in the
pyazr chapter.
