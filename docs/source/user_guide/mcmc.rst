MCMC Tab
========

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
     - ``level`` or ``level_rwa``.

Populate the table with **Load Physical Parameters** or **Load RWA
Parameters**.

.. note::

   **Normalizations and energy shifts are not listed, and do not need to be.**
   They are still sampled — AZURE2 builds their priors automatically from the
   errors quoted in the Segments tab, reproducing exactly the χ² penalties a
   Minuit fit would apply (see :doc:`chi_squared`). A normalization with no
   quoted error gets a uniform prior instead of an invented one, and the run
   log says so.

   Level energies and widths carry no experimentally implied uncertainty, so
   those are yours to specify. Any left without a prior are sampled with an
   unbounded uniform prior; the run log lists them by name so nothing is
   unconstrained without you being told.

Sampling settings
-----------------

**Number of Walkers**
   Size of the ensemble. Must exceed the number of free parameters — see
   above. Two to four times that is a reasonable starting point.

**Number of Steps**
   Steps per walker. The total sample count is walkers × steps.

**Initial Chain Spread (%)**
   Spread of the starting ensemble, as a percentage of each parameter's value.
   Applies to widths.

**Level Energy Spread (keV)**
   Spread of the starting ensemble for **level energies**, as an absolute
   energy, capped at 1 keV.

   .. warning::

      Level energies cannot be scattered by a percentage the way widths can.
      One percent of a 5 MeV level energy is 50 keV, which places walkers on
      completely different resonance structures; the ensemble never contracts
      and the chain does not converge. This is why the setting is separate and
      bounded.

   Normalizations and energy shifts are started inside their own prior width,
   since beginning outside it only wastes steps walking back in.

**Number of Threads**
   Walkers within a half-ensemble are evaluated in parallel. Results do not
   depend on the thread count: each walker carries its own generator, so a run
   is reproducible for a given seed whatever the threading.

**Fresh Start**
   Start a new chain, discarding any existing ``samples.mcmc``.

**Use Reduced Widths**
   Sample reduced width amplitudes rather than physical widths. Usually
   preferable — it is the natural fit space, and avoids a transformation at
   every evaluation.

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

Output files
------------

``samples.mcmc``
   The chain, one row per walker per step::

       step,walker,logprob,loglikelihood,logprior,param0,param1,...

   Every accepted state appears exactly once, including the repeats a
   rejection contributes — that repetition *is* how a Markov chain encodes
   probability mass, so do not deduplicate it. ``logprob`` is exactly
   ``loglikelihood + logprior``.

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
``emcee`` with its full range of moves and convergence diagnostics, ``zeus``,
or nested sampling with ``dynesty``. Worked examples ship with the package:
``pyazr/examples/fit_emcee.py`` and ``fit_zeus.py``.
