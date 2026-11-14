MCMC Tab
========

.. note::

   This page describes the ``.azr`` content through the fields of the Qt setup
   utility of the released versions. This version has no setup utility;
   :doc:`graphical_interface` says how to do the same from the command line and
   with ``pyazr``.

The **MCMC** tab performs Bayesian parameter inference by sampling the
posterior directly, rather than locating a single best fit. It must be enabled
at compile time (``USE_MCMC``, on by default).

Where a Minuit fit returns one parameter set and a covariance matrix — a
Gaussian approximation valid near the minimum — sampling returns the shape of
the posterior itself. That matters for R-matrix work, where parameters are
often strongly correlated, bounded (a width cannot exceed the Wigner limit),
or simply not Gaussian.

The sampler
-----------

AZURE2 uses an **affine-invariant ensemble sampler** — the stretch move of
Goodman & Weare (2010), parallelised over two half-ensembles following
Foreman-Mackey et al. (2013). It is the same algorithm as ``emcee``'s
``StretchMove``, with :math:`a = 2`.

An ensemble of *walkers* explores the parameter space together, each proposing
moves along the line joining it to another walker. Because proposals are built
from the ensemble's own spread, the sampler adapts to correlated parameters
without being told about them — which is why it suits R-matrix posteriors,
where level energies and widths are rarely independent.

Two consequences follow from how the move works:

- The ensemble can never leave the affine hull of its starting positions, so
  **the walker count must exceed the number of free parameters** — comfortably.
  Two to four times is the usual advice, and ``emcee`` refuses to run below
  twice the dimension.
- Walkers are not independent chains. Statistics must be taken over the whole
  ensemble, not per walker.

Parameters and priors
---------------------

The **Parameters** sub-tab lists the level energies and widths, with the prior
you wish to impose on each:

.. list-table::
   :widths: 22 78
   :header-rows: 1

   * - Column
     - Description
   * - **Parameter**
     - Which level and channel this is — spin-parity, level energy, particle
       pair, :math:`L` and :math:`S`.
   * - **Current Value**
     - Starting value for the ensemble.
   * - **Prior Mean**
     - Centre of the Gaussian prior.
   * - **Prior Std**
     - Width of the Gaussian prior.
   * - **Use Gaussian Prior**
     - Impose the prior. Unchecked means an unbounded uniform prior.
   * - **Category**
     - ``level`` or ``level_rwa``; ``abundance`` for the isotopic abundances
       of a project with several compound nuclei (:doc:`particle_pairs`),
       whose prior is taken from the abundance table's value and sigma like a
       normalization's from the segment, and whose sum-to-one term enters the
       log-likelihood as it enters the chi-squared; ``cbkg`` for a THM coherent
       background.

Populate the table with **Load Physical Parameters** or **Load RWA
Parameters**. The latter reads ``param.par`` of the output directory by
parameter name, as AZURE2 reads a parameter file (a file in the other of
Brune's and Park's bases is converted).

.. note::

   **Normalizations and energy shifts are not listed, and do not need to be.**
   They are still sampled — AZURE2 builds their priors automatically from the
   errors quoted in the Segments tab (and, for a :math:`\sqrt{E}` shift
   coefficient, from its ``sqrtshift`` block), centred where the Minuit
   penalties are centred, prior centre rows included, so that they reproduce
   exactly the χ² penalties a Minuit fit would apply (see :doc:`chi_squared`).
   A normalization with no quoted error gets a uniform prior instead of an
   invented one, and the run log says so. The free norm of a THM segment is
   not a parameter at all: it is profiled at every evaluation.

   Level energies and widths carry no experimentally implied uncertainty, so
   those are yours to specify. Any left without a prior are sampled with an
   unbounded uniform prior; the run log lists them by name so nothing is
   unconstrained without you being told.

   A THM coherent background (``cbkg_*``) gets no automatic prior either: it
   is uniform unless you define one, and the run log counts them.

   The sampler tells the kinds apart by what each parameter is (its kind,
   recorded where AZURE2 creates it; :doc:`../reference/output_files`, "Kinds
   of parameter"), not by its name. A row of the ``<mcmc>`` ``<parameters>``
   block names its parameter by the whole name of ``param.par``
   (``energy_1``, ``width_1_2``, ``s01.norm``, ...), optionally as
   ``label (name)``, or by the short forms ``E<n>``, ``G<n><c>``, ``N<k>``,
   ``S<k>``, followed by the value, prior mean, prior sigma and the use flag
   (0 or 1). A row that names no parameter of the model stops the run with a
   message, rather than leaving its prior out.

   A prior parameter covariance (``priorCovariance=<file>`` in ``<config>``,
   :doc:`chi_squared`, "Prior parameter covariance") enters the log prior as the
   multivariate Gaussian :math:`-\tfrac12(\mathbf{p}-\boldsymbol{\mu})^T C^{-1}
   (\mathbf{p}-\boldsymbol{\mu}) - \tfrac12\ln\det(2\pi C)` and replaces the
   diagonal prior (automatic or from this table) of each parameter it lists; the
   run log names them. It is defined on the reduced-width amplitudes, so with
   **Use Reduced Widths** off a prior on a level energy or width is refused.

Limits
------

The sampler never leaves the fit's limits: the posterior is zero at a point
where any free parameter lies outside the limits a fit would give it. These
are the limits each kind of parameter carries (a normalization, an abundance,
a screening potential and a sample thickness are not negative; a fitted areal
density, Doppler temperature, stopping scale or kernel parameter stays inside
the range its grids are laid out for; a ``gamma_decay`` record's
:math:`\Gamma_\gamma \ge 0` and :math:`0 \le b_f \le 1`) and those written in
``<parameterSettings>`` (:doc:`chi_squared`). A parameter without a limit in
the fit has none in the sampler. A proposal outside is rejected without being
evaluated, and a walker of the starting ensemble drawn outside is drawn again
(the run log says how many).

With **Use Reduced Widths** off the sampler moves the observed widths and
energies, and the limits act on the reduced-width amplitudes and level
energies they are converted to: a ``<parameterSettings>`` limit on a width,
written in eV, is converted to an amplitude once at the start, as for a fit,
so at the start's level energy the sampler stops at the width written, and
where the level energy moves the edge in the width moves with the conversion.

A Gaussian prior of this table that is wider than a limit is cut by it: the
posterior is the prior times the likelihood inside the limits, zero outside,
so a prior centred near a limit does not describe the sampled distribution
there.

Sampling settings
-----------------

The program takes them as options of ``AZURE2 mcmc``
(:doc:`../reference/command_line`); the names below are the setup utility's.

**Number of Walkers** (``--walkers=N``)
   Size of the ensemble. Must exceed the number of free parameters — see
   above. Two to four times that is a reasonable starting point.

**Number of Steps** (``--steps=N``)
   Steps per walker. The total sample count is walkers × steps.

**Initial Chain Spread (%)** (``--parameter-spread=P``, default 5)
   Spread of the starting ensemble, as a percentage of each parameter's value.
   Applies to widths.

**Level Energy Spread (keV)** (``--energy-spread=E``, default 1)
   Spread of the starting ensemble for **level energies**, as an absolute
   energy, capped at 1 keV.

   .. warning::

      Level energies cannot be scattered by a percentage the way widths can.
      One percent of a 5 MeV level energy is 50 keV, which places walkers on
      completely different resonance structures; the ensemble never contracts
      and the chain does not converge. This is why the setting is separate and
      bounded.

   Normalizations and energy shifts are started inside their own prior width,
   since beginning outside it only wastes steps walking back in. A walker
   drawn outside the fit's limits is drawn again ("Limits").

**Number of Threads** (``--sampler-threads=N``)
   Walkers within a half-ensemble are evaluated in parallel. Results do not
   depend on the thread count: each walker carries its own generator, so a run
   is reproducible for a given seed whatever the threading.

**Fresh Start** (``--overwrite-samples``)
   Start a new chain, discarding any existing ``samples.mcmc``.

**Use Reduced Widths** (``--basis=rwa``; ``--basis=physical`` without)
   Sample reduced width amplitudes rather than physical widths. Usually
   preferable — it is the natural fit space, and avoids a transformation at
   every evaluation. Both modes sample the normalizations and energy shifts
   with the level parameters and apply them in the likelihood.

   .. note::

      Before 2026-10-26 the physical-parameter mode applied only the level
      parameters and the energy shifts: the sampled normalizations never
      reached the likelihood, so their posterior was their prior, the
      R-matrix posterior was taken at the starting normalizations, and
      ``param.mcmc`` recorded the starting values. A physical-mode chain with
      a free normalization from an older build should be rerun.

Monitoring a run
----------------

The **Progress** sub-tab reports the current step, the best walker's log
probability, log likelihood and log prior, and an estimated time remaining.

The run log additionally reports, at the end:

- the **acceptance fraction**. Healthy ensemble sampling sits roughly between
  0.2 and 0.5. Near zero means the chain is barely moving and the samples are
  not yet a usable posterior; a warning is printed below 0.05.
- any **walkers that started at zero posterior probability** and therefore
  cannot move — reduce the initial spread or check the priors.

Before the first step the starting ensemble is checked. Under Park's
parametrization (``--use-park``) a starting point where a level has
:math:`J = 1-\sum_c\gamma_c^2\,dS_c/dE \le 0` is reported, the level named
with its :math:`J`: the posterior is zero there (:math:`J<0` is rejected), so
walkers near it are rejected. If **no** walker of the starting ensemble has a
finite posterior probability — under Park, typically every walker of the ball
around such a point — the run is refused ("MCMC not started"): an ensemble
entirely outside the support carries no information on the posterior, and
the chain would be the start repeated with acceptance zero. No
``samples.mcmc`` is written, the GUI shows the report, and the command line
exits with a non-zero status (as it does for any MCMC run that fails).
Start from a fit made under ``--use-park`` (its penalty keeps :math:`J>0`), from
a Brune parameter file (converted on read) or with smaller widths. A start
from which part of the ensemble is inside the support runs, with the warning.

Output files
------------

``samples.mcmc``
   The chain, one row per walker per step::

       step,walker,logprob,loglikelihood,logprior,param0,param1,...

   Every accepted state appears exactly once, including the repeats a
   rejection contributes — that repetition *is* how a Markov chain encodes
   probability mass, so do not deduplicate it. ``logprob`` is exactly
   ``loglikelihood + logprior``.

   ``logprob`` is :math:`-\tfrac12\chi^2` of the fit's objective, term for
   term, plus the normalisation of the Gaussian priors, and the sampler
   rejects a point outside the support (Park's :math:`J \le 0`; the fit's
   limits of every parameter, "Limits" above).  ``logprior`` is the normalised Gaussian priors:
   the automatic ones, the nuisance priors, the prior covariance and those of
   this table.  ``loglikelihood`` is the rest of the objective: the data, the
   abundances' sum-to-one term, the model data, and the priors and the
   remainder's wall of a ``gamma_decay`` record (:doc:`gamma_decay`).  A
   point outside the support is never written.

   ``loglikelihood`` includes the measured scattering lengths of a
   ``<scatteringLengths>`` block (:doc:`../theory/scattering_length`).  For the
   posterior of a scattering length itself, evaluate
   ``Session.scattering_lengths`` on the chain's rows in pyazr.

``walkers.mcmc``
   The final position of every walker, so a later run continues the same
   ensemble rather than re-scattering it. Written on normal completion and on
   an early stop.

Resuming
--------

Leaving **Fresh Start** unchecked continues an existing run: AZURE2 counts the
steps already in ``samples.mcmc``, restores the ensemble from
``walkers.mcmc``, and appends. If the walker state is missing or does not match
the current model, the walkers restart from a fresh spread and the log says so
— the steps that follow are then a new burn-in rather than a continuation.

Analysing the chain
-------------------

The **Results** sub-tab shows the mean, standard deviation, median and 95%
credible interval per parameter, over every sample in the file. **Refresh
Results from File** recomputes them from an existing ``samples.mcmc``.

For anything beyond that — trace plots, corner plots, autocorrelation times —
read the CSV directly. It loads with one line of ``pandas`` or ``numpy``:

.. code-block:: python

   import numpy as np, pandas as pd

   df = pd.read_csv("output/samples.mcmc")
   nwalkers = df["walker"].nunique()
   nsteps   = df["step"].nunique()
   ndim     = len([c for c in df.columns if c.startswith("param")])

   # (steps, walkers, ndim), the shape emcee's own tooling expects
   chain = df.filter(like="param").to_numpy().reshape(nsteps, nwalkers, ndim)

   burn = nsteps // 5                 # discard burn-in, then flatten
   flat = chain[burn:].reshape(-1, ndim)

.. important::

   **Discard the burn-in.** The Results tab summarises the entire file,
   starting ensemble included. Early samples reflect where the walkers were
   put, not the posterior.

   **Check convergence before quoting anything.** AZURE2 reports the
   acceptance fraction but not an autocorrelation time, so it cannot tell you
   whether the chain is long enough. ``emcee.autocorr.integrated_time`` applied
   to the array above will, and a chain shorter than roughly 50 autocorrelation
   times should not be trusted.

Sampling from Python instead
----------------------------

Running the sampler through :doc:`pyazr` gives access to the wider ecosystem —
``emcee`` with its full range of moves and convergence diagnostics, or nested
sampling with ``dynesty``. ``s.objective().log_prob`` is the posterior this sampler takes (−½ × the
fit's objective, −∞ outside its support, without the priors' normalisation constant); see :doc:`pyazr`, "The
objective and the evaluation".  The worked example is ``pyazr/examples/fit_emcee.py``: ``log_prob``, one
session per pool worker, the autocorrelation times beside the chain's length and a rate per posterior sample.
