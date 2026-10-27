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
the input file is only where the fit starts. A segment whose :math:`\sqrt{E}`
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
- **Total-Scattering-Length-Chi-Squared**, when the project has a
  ``<scatteringLengths>`` block, is the sum over its active lines of
  :math:`((a_0 - a_{\mathrm{exp}})/\sigma)^2` (real and, if given, imaginary
  part), with :math:`a_0` the s-wave scattering length at the pair's threshold
  (:doc:`../theory/scattering_length`).  ``scatteringLengths.out`` lists each
  line's model value and contribution.  The term is part of every fitter's
  objective and of the MCMC likelihood.

The quantity the minimizer actually descends is the sum of all of them. Quoting
``Total-Chi-Squared`` as "the χ²" of a fit with free normalizations understates
the objective.

.. note::

   ``Total-N`` counts data points, not degrees of freedom. Subtract the number
   of free parameters yourself before forming a reduced χ².

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
``residual_jacobian`` return the **data** χ² only. A least-squares fit driven
from Python that minimises those residuals alone will let the normalizations
drift to absorb every discrepancy, reaching a "better" χ² that AZURE2 itself
would never find. Append the penalty rows explicitly — the recipe is in the
pyazr chapter.
