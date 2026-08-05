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

The energy-shift penalty
------------------------

Segments with **Vary Energy Shift?** enabled contribute the analogous term,
this time with an absolute rather than a percentage error:

.. math::

   \chi^2_{\text{shift},s} =
   \left(\frac{\Delta E_s - \Delta E_s^{\text{nom}}}{\delta(\Delta E_s)}\right)^2

Nuisance parameters
-------------------

Level energies and widths flagged **Use as Nuisance** in the Fitting tab add a
Gaussian term of the same shape, using the **Error** column as the width. This
is how prior knowledge of a resonance energy enters a least-squares fit.

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
