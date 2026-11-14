Example scripts of pyazr
========================

The scripts of ``pyazr/examples/`` each tell one task from start to end, for an expert (or an agent working for one)
to copy.  Each runs from the repository's root on a project of ``tests/``, takes the project as its first argument,
writes files only where ``--out`` (or a printed temporary directory) says, and prints its numbers with their units;
its docstring says what it shows and every option.  ``tests/pyazr/examples_test.py`` runs every script by the command
below and pins one number of its output.  They are not installed with the package (they read the repository's test
projects); the optional packages they use are the extra ``examples`` (``pip install ".[examples]"``).  The projects
of ``examples/`` (THM evaluations) are listed in :doc:`../getting_started/examples`.

Each command is ``python3 pyazr/examples/`` followed by what the table shows.

Build and edit a project
------------------------

.. list-table::
   :widths: 24 44 32
   :header-rows: 1

   * - Script
     - Task
     - Command
   * - ``project_from_python.py``
     - Write a project record by record from Python: schema, a target, a beam, a data set refused and accepted,
       calculate (:doc:`project_from_python`).
     - ``project_from_python.py [WORKDIR]``
   * - ``print_scheme.py``
     - Read what the engine made of a project: pairs, levels, channels, widths, data sets.
     - ``print_scheme.py tests/upgrade/T1/T1.azr``
   * - ``edit_model.py``
     - Remove, add and switch off levels; remove and add data sets; calculate the edited copy.
     - ``edit_model.py tests/upgrade/T1/T1.azr --out DIR``
   * - ``exfor_fetch.py``
     - Fetch an EXFOR data set (or read a local copy), add it to a project, find its reference.
     - ``exfor_fetch.py tests/upgrade/T1/T1.azr --entry T0001003.1 --local tests/pyazr/data/exfor_master``
   * - ``covariance_builder.py``
     - Build a data set's covariance file from named parts and name it on the data set (:doc:`covariance_builder`).
     - ``covariance_builder.py tests/upgrade/T1/T1.azr el_lab090 --out DIR/T1_cov.azr``
   * - ``prior_builder.py``
     - Build a prior on the fit parameters from named parts and name it in the project (:doc:`prior_builder`).
     - ``prior_builder.py tests/upgrade/T1/T1.azr --out DIR/T1_prior.azr``
   * - ``nuisances.py``
     - Give data sets one shared normalization and one shared energy shift, each with a prior, and fit them.
     - ``nuisances.py tests/upgrade/T1/T1.azr el_lab090 el_lab125 el_lab160``
   * - ``natural_target.py``
     - A target of several isotopes: each nucleus's cross section weighted by its abundance.
     - ``natural_target.py tests/natural_target/models/M.azr``
   * - ``gamma_decay.py``
     - A level's radiative decay as Γγ, branchings and a remainder (a ``gamma_decay`` record), fitted.
     - ``gamma_decay.py tests/gamma_decay/models/record.azr``
   * - ``refusals.py``
     - What a refusal looks like and how a script reads it: the three parts, the diagnostics records.
     - ``refusals.py tests/upgrade/T3/T3.azr``

Calculate and compare
---------------------

.. list-table::
   :widths: 24 44 32
   :header-rows: 1

   * - Script
     - Task
     - Command
   * - ``per_dataset_chi2.py``
     - The χ² of each data set, worst first, and its whitened rows of the objective.
     - ``per_dataset_chi2.py tests/upgrade/T2/T2.azr``
   * - ``model_curve.py``
     - The model of every data set on a fine grid over its range (``Session.model_curve``), and the maximum a line
       through the points misses.
     - ``model_curve.py tests/upgrade/T2/T2.azr``
   * - ``decompose_cross_section.py``
     - Split a cross section into each level, its interference and the background; one level off or alone.
     - ``decompose_cross_section.py tests/upgrade/T1/T1.azr --level "1/2+ 2.365"``
   * - ``angular_distribution.py``
     - Legendre coefficients of an angular distribution at chosen energies.
     - ``angular_distribution.py tests/upgrade/T1/T1.azr``
   * - ``channel_radius_scan.py``
     - χ² against a channel radius, one edited copy and one session per radius.
     - ``channel_radius_scan.py tests/upgrade/T3/T3.azr --pair 2 --radii 3.5 4.0 4.5``
   * - ``nuclear_potential.py``
     - The hybrid nuclear potential of a pair: a depth scan, and why a fit vector belongs to its session.
     - ``nuclear_potential.py tests/hybrid_potential/hybrid_potential.azr``
   * - ``thm_dataset.py``
     - THM data sets of one experiment: the profiled norm they share, their χ² beside the direct data's.
     - ``thm_dataset.py tests/upgrade/T2/T2.azr``
   * - ``tabulate_cross_sections.py``
     - Tabulate σ(E) on an adaptive grid to a stated interpolation tolerance, for another code.
     - ``tabulate_cross_sections.py tests/upgrade/T1/T1.azr --out sigma.dat``

Fit
---

.. list-table::
   :widths: 24 44 32
   :header-rows: 1

   * - Script
     - Task
     - Command
   * - ``fit_engine.py``
     - Fit with the engine's MIGRAD, HESSE and MINOS in-process; save the fit as a project that reopens as the fit.
     - ``fit_engine.py tests/upgrade/T1/T1.azr --out DIR``
   * - ``fit_iminuit.py``
     - Fit with iminuit (MIGRAD, HESSE) on the engine's objective and its gradient.
     - ``fit_iminuit.py tests/upgrade/T3/T3.azr``
   * - ``fit_scipy.py``
     - Fit with scipy: L-BFGS-B with the gradient, least squares with the Jacobian.
     - ``fit_scipy.py tests/upgrade/T3/T3.azr``

Uncertainties
-------------

.. list-table::
   :widths: 24 44 32
   :header-rows: 1

   * - Script
     - Task
     - Command
   * - ``fit_emcee.py``
     - Sample the posterior with emcee: the objective's priors and limits, autocorrelation times, a rate per sample.
     - ``fit_emcee.py tests/upgrade/T3/T3.azr --steps 200 --rate 2 1``
   * - ``uncertainty_band.py``
     - The 1σ band of an S factor from a fit's covariance, on an extrapolation and on a new grid.
     - ``uncertainty_band.py tests/upgrade/T1/T1.azr``
   * - ``sensitivities.py``
     - Which parameters a data set constrains (d ln σ / d ln p); analytic against finite differences.
     - ``sensitivities.py tests/upgrade/T1/T1.azr --check``

Rates and S factors
-------------------

.. list-table::
   :widths: 24 44 32
   :header-rows: 1

   * - Script
     - Task
     - Command
   * - ``reaction_rate.py``
     - N\ :sub:`A`\ ⟨σv⟩ computed by the engine, the Gamow window, a band, REACLIB sets (:doc:`reaction_rate`).
     - ``reaction_rate.py tests/upgrade/T1/T1.azr --entrance 1 --exit 2``
   * - ``sfactor_extrapolation.py``
     - An S factor below the measured range on a c.m. grid, and S(0).
     - ``sfactor_extrapolation.py tests/upgrade/T1/T1.azr --entrance 1 --exit 2``

Evaluations
-----------

.. list-table::
   :widths: 24 44 32
   :header-rows: 1

   * - Script
     - Task
     - Command
   * - ``endf_evaluation.py``
     - An ENDF-6 evaluation beside a project: read it, draw it over the data, check the file (:doc:`endf`).
     - ``endf_evaluation.py``
   * - ``endf_start_model.py``
     - An evaluation's resonance parameters (LRF=7) as a new project, and a project written as an ENDF-6 file
       (:doc:`endf`).
     - ``endf_start_model.py``

Physics functions
-----------------

.. list-table::
   :widths: 24 44 32
   :header-rows: 1

   * - Script
     - Task
     - Command
   * - ``channel_functions.py``
     - Coulomb functions, penetrability, shift and hard-sphere phase of a channel over an energy grid.
     - ``channel_functions.py tests/upgrade/T1/T1.azr``
   * - ``ec_integrals.py``
     - The external-capture radial integrals, pathway by pathway.
     - ``ec_integrals.py tests/upgrade/T1/T1.azr``
   * - ``dimensionless_widths.py``
     - θ² of every particle width and Weisskopf units of every γ width.
     - ``dimensionless_widths.py tests/upgrade/T1/T1.azr``
   * - ``transform_widths.py``
     - Reduced-width amplitudes to physical widths with ``Session.transform``, beside 2γ²P; a chain mapped.
     - ``transform_widths.py tests/upgrade/T1/T1.azr``

Adding an example
-----------------

A new example is its file, its row in ``ROWS`` of ``tests/pyazr/examples_test.py`` and its row on this page
(:doc:`../developer/extending`, "An example script of pyazr"); consistency check 58 fails without either row.
