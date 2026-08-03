What AZURE2 Can Do
==================

A map of the whole package, with pointers to the chapter covering each part.
If you are looking for "can AZURE2 do X", start here.

The model
---------

AZURE2 implements **multi-channel, multi-level R-matrix theory** in the
formulation of Lane and Thomas. You describe a compound nucleus in terms of:

- **Particle pairs** — the reaction participants, with masses, charges, spins,
  separation and excitation energies, and a channel radius
  (:doc:`../user_guide/particle_pairs`).
- **Levels and channels** — resonances of the compound nucleus, each decaying
  through channels labelled by orbital angular momentum :math:`L`, channel
  spin :math:`S` and particle pair. Allowed channels are generated
  automatically from the pairs (:doc:`../user_guide/levels_channels`).
- **Segments** — what to compare against or compute: experimental datasets, or
  grids of energies and angles for prediction
  (:doc:`../user_guide/segments`).

Both the standard R-matrix parametrisation and the **alternative level matrix
of C. R. Brune** are supported, the latter with ``--use-brune``. Capture may
use the **Reich–Moore approximation** (``--use-rmc``), and external (direct)
capture is included.

Observables
-----------

- Angle-integrated and **differential cross sections**, in the lab or
  centre-of-mass frame
- **S-factors**
- **Angular distributions** (Legendre coefficients)
- **Phase shifts**
- **Total capture**, summed over final states
- **Reaction rates** as a function of temperature

All output is centre-of-mass, regardless of the frame the input used
(:doc:`../reference/output_files`).

Experimental effects
--------------------

Calculations can be corrected for the things that separate a measurement from
a point cross section (:doc:`../user_guide/experimental_effects`):

- **Target integration** — finite target thickness, with SRIM stopping powers
- **Beam energy convolution** — finite beam resolution and straggling
- **Detector geometry** — angular acceptance via Q-coefficients
- **Per-segment normalizations** and **energy shifts**, either fixed or
  treated as free parameters constrained by their quoted experimental
  uncertainties

Fitting
-------

Parameters are fitted by **least squares** against the data
(:doc:`../user_guide/fitting`). What enters the objective — the data term, the
normalization and energy-shift penalties, nuisance parameters — is set out in
:doc:`../user_guide/chi_squared`.

- **Minuit2** (MIGRAD) is the default minimizer
- A **Levenberg–Marquardt** minimizer using the analytic Jacobian is available
  with ``--use-lm``, and is usually much faster for width-dominated fits
- **MINOS** error analysis gives asymmetric uncertainties
- **Covariance-based uncertainty bands** on the calculated cross section
  (``--covariance-band``)
- Parameter **limits** and **Wigner-limit bounds**, which can be populated
  automatically

Bayesian inference
------------------

Instead of a single best fit, sample the posterior with an
**affine-invariant ensemble sampler** — the same algorithm as ``emcee``
(:doc:`../user_guide/mcmc`). Priors for normalizations and energy shifts are
built automatically from the quoted experimental errors; priors on level
energies and widths are yours to state. The chain is written as CSV for
analysis in any tool.

Scripting
---------

Everything above is reachable from Python through :doc:`../user_guide/pyazr`,
which runs headless AZURE2 processes and talks to them over a socket. That
covers custom minimizers, external samplers (``emcee``, ``zeus``, ``dynesty``),
parameter scans, cross-section decomposition into individual level and
interference contributions, dimensionless widths, and programmatic editing of
the model itself.

Interfaces
----------

.. list-table::
   :widths: 22 78
   :header-rows: 1

   * - Interface
     - Use it for
   * - **GUI**
     - Building and inspecting a model, running fits, plotting. The natural
       place to start.
   * - **Console** (``--no-gui``)
     - Scripted or remote runs, batch jobs, HPC
       (:doc:`../reference/command_line`).
   * - **Socket API** + ``pyazr``
     - Anything programmatic (:doc:`../user_guide/pyazr`).

What AZURE2 does not do
-----------------------

Stated plainly, so you do not go looking:

- It does not compute nuclear structure. Levels, spins and parities are input.
- It does not select a model for you. Which levels to include, and which
  channels to free, are physics decisions.
- The MCMC reports an acceptance fraction but no autocorrelation time, so it
  cannot tell you whether a chain has converged. Use ``emcee``'s diagnostics on
  the output file.
- Uncertainty bands are covariance-based and therefore Gaussian near the
  minimum. For a non-Gaussian posterior, sample it.
