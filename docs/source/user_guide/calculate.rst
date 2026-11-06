Calculate Tab
=============

The **Calculate** tab controls the execution of AZURE2 calculations. Select a
calculation type from the drop-down menu and click **Save and Run** to execute.

Calculation Types
-----------------

Calculate With Data
^^^^^^^^^^^^^^^^^^^

Performs a single R-matrix calculation using the current parameters, comparing
against the experimental data defined in the **Segments** tab. No fitting is
performed.

This mode is useful for:

- Checking the quality of initial parameter values before fitting.
- Manually adjusting parameters (especially interference signs) until a
  reasonable starting point is reached.

Fit With Data
^^^^^^^^^^^^^

Performs an automated least-squares (:math:`\chi^2`) fit to the experimental
data using the MINUIT2 minimization package. The fit adjusts all unfixed
parameters to minimize the chi-squared.

During the fit:

- The :math:`\chi^2` value is printed every 10 evaluations.
- Output files (and ``param.fit``) are written every 100 evaluations, so
  progress can be monitored in the **Plot** tab.
- The final :math:`\chi^2` is printed to the output area and to
  ``chiSquared.out``.

.. important::

   Each fit overwrites all standard output files. If you are testing different
   fitting options, rename or move the output files before starting a new fit.

Extrapolate Without Data
^^^^^^^^^^^^^^^^^^^^^^^^

Performs a single calculation (no minimization) using the segments defined in the
**Segments Without Data** section. No experimental data is considered.

Common uses:

- Interpolate a cross section with finer energy spacing for plotting.
- Extrapolate to astrophysically relevant energies.
- Calculate the total cross section from a differential fit.
- Compute the cross section for the inverse reaction.
- Calculate a previously unobserved cross section using only level parameters.

Perform MINOS Error Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Runs the MINOS uncertainty estimation routine from MINUIT2 to calculate
asymmetric parameter uncertainties. Results are written to ``param.errors``
and ``covariance_matrix.out``.

.. tip::

   Running MINOS is also a useful test of fit robustness. If MINOS crashes, it
   often indicates that one or more parameters are unconstrained by the data.
   Try fixing different parameters to identify the problematic one(s).

Calculate Reaction Rate
^^^^^^^^^^^^^^^^^^^^^^^

Calculates the thermonuclear reaction rate by numerically integrating the cross
section over a Maxwell-Boltzmann energy distribution:

.. math::

   N_A \langle \sigma v \rangle = \left(\frac{8}{\pi\mu}\right)^{1/2}
   \frac{N_A}{(kT)^{3/2}} \int_0^\infty \sigma(E) \, E \, e^{-E/kT} \, dE

Additional fields appear for this mode:

**Entrance / Exit Pair Keys**
   Specify the particle pair keys for the reaction. For the **total radiative
   capture** rate choose "Total capture (all γ pairs)" as the exit pair (on the
   command line, exit pair ``-1``, the key a total-capture segment line carries):
   the rate then integrates the sum of the capture cross sections into every
   capture (γ, pair type 10) pair of the entrance pair's compound nucleus, e.g.
   capture to the ground state and to each bound excited state. The run log names
   the pairs summed, and ``cross_section_grid.dat`` holds the summed cross
   section. The entrance must be a particle pair. Builds before 2026-11-04 refuse
   the key; with them, add the rates of the individual capture pairs (the rate is
   linear in σ).

**Temperature Source**
   Choose between:

   - **Grid** -- specify minimum, maximum, and step temperatures (in GK).
   - **File** -- provide a file listing specific temperatures (single column, in GK).

Results are written to ``reactionrates.out``; the cross-section grid the rate
integrated is written to ``cross_section_grid.dat``.

.. warning::

   The rate integrates a spline of ln σ against ln E through a grid of 48 points per
   decade from 0.1 eV to 50 MeV (centre of mass) plus a cluster of points from 1/32
   to 64 widths around every level of the parameters being evaluated, so narrow
   resonances are resolved (a 22 eV level to 1e-4 in ``tests/reaction_rate``). The
   grid is written to ``cross_section_grid.dat``; a resonance outside 0.1 eV to
   50 MeV, or a level without widths, is not integrated. Builds before 2026-10-26
   splined σ itself (charged-particle rates about 1 % low at low temperature), set
   σ = 0 below 1 keV (neutron rates 49 % low at T9 = 0.01) and placed the clusters
   at the ``.azr`` level energies rather than at the parameters in use.

Minimizer Selection
-------------------

**Minuit2** (default)
   The standard MINUIT2 minimizer (MIGRAD) from CERN, with numerical
   derivatives.

**Minuit2 (analytic grad.)**
   MIGRAD with the analytic gradient (CLI ``--use-gradient``).

**Levenberg-Marquardt**
   Least squares on the analytic Jacobian (CLI ``--use-lm``); falls back to
   MIGRAD where the Jacobian is not available.

**GSL Trust-Region**
   GSL's trust-region least squares with geodesic acceleration (CLI
   ``--use-gsl-lm``); falls back to MIGRAD likewise.

**NLopt** (if built with ``USE_NLOPT``)
   Alternative minimizer with several algorithms:
   SBPLX, COBYLA, BOBYQA, NEWUOA, PRAXIS, Nelder-Mead.

The minimizer selection is only available for fitting modes.

**Chi-squared Variance**
   The target variance for the chi-squared minimization (default: 1.0). Only
   applies to fitting modes.

**Uncertainty**
   Compute the analytic cross-section uncertainty band from the fitted
   parameter covariance (CLI ``--covariance-band``): at the data energies for
   a fit, on the extrapolation grid for "Calculate Segments Without Data",
   from the covariance a previous fit saved (``covariance.dat``).
   **Scale covariance** inflates the covariance by the reduced chi-squared
   when it exceeds 1 (``--scale-covariance``).

**Set Wigner limits**
   Bound every reduced-width parameter by its Wigner (single-particle)
   limit during the fit.

Parameter and Integral Files
----------------------------

**Parameters File**
   - **Create New Parameters File** (default) -- use the parameters currently
     entered in the GUI.
   - **Use** -- load parameters from a previously saved file (e.g., ``param.sav``
     from a completed fit).

**External Capture Integrals File**
   - **Create New Integrals File** (default) -- compute the external capture
     integrals from scratch.
   - **Use** -- load a previously computed integrals file (``intEC.dat``) to
     speed up calculations.  Since the integrals come from energy tables
     (:doc:`../developer/ec_table`), a new file takes seconds even for
     convolved or target-integrated capture data, and the run log reports
     ``External capture table: N integrals, M computed at table nodes``.

   .. note::

      The integrals file can be reused as long as the calculation segments,
      channel radii, :math:`J^\pi` values and integration grids (target
      effects) have not changed. Level energies, widths, and ANCs can be
      changed freely. AZURE2 verifies this against the signature it keeps in
      ``intEC.dat.sig`` and recomputes, with a warning, a file that does not
      match (see :doc:`../reference/output_files`).

Execution Controls
------------------

**Save and Run**
   Saves the project and starts the selected calculation. Output is displayed
   in real time in the text area at the bottom.

**Stop AZURE2**
   Gracefully stops a running calculation. There may be a brief delay before
   the calculation concludes. Output files from the last completed iteration
   are preserved.
