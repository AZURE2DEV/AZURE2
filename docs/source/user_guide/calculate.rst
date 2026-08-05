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

- The :math:`\chi^2` values are updated every 1000 iterations.
- Output files are written periodically so progress can be monitored in the
  **Plot** tab.
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
   Specify the particle pair keys for the reaction.

**Temperature Source**
   Choose between:

   - **Grid** -- specify minimum, maximum, and step temperatures (in GK).
   - **File** -- provide a file listing specific temperatures (single column, in GK).

Results are written to ``reactionrates.dat``.

.. warning::

   Numerical integration is limited to broad resonance structures. If resonances
   are narrower than approximately 1 keV, the integration may fail. For narrow
   resonances, calculate the contribution separately using the narrow-resonance
   approximation.

Minimizer Selection
-------------------

**Minuit2** (default)
   The standard MINUIT2 minimizer from CERN.

**NLopt** (if available)
   Alternative minimizer with several algorithms:
   SBPLX, COBYLA, BOBYQA, NEWUOA, PRAXIS, Nelder-Mead.

The minimizer selection is only available for fitting modes.

**Chi-squared Variance**
   The target variance for the chi-squared minimization (default: 1.0). Only
   applies to fitting modes.

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
     speed up calculations.

   .. note::

      The integrals file can be reused as long as the calculation segments,
      channel radii, and :math:`J^\pi` values have not changed. Level energies,
      widths, and ANCs can be changed freely.

Execution Controls
------------------

**Save and Run**
   Saves the project and starts the selected calculation. Output is displayed
   in real time in the text area at the bottom.

**Stop AZURE2**
   Gracefully stops a running calculation. There may be a brief delay before
   the calculation concludes. Output files from the last completed iteration
   are preserved.
