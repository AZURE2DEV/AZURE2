Output Files
============

AZURE2 produces several output files in the configured output directory. An
important distinction: while all **input** is in the laboratory frame, all
quantities in **output** files are in the **center-of-mass frame**.

Parameter Files
---------------

param.par
^^^^^^^^^

Contains the initial formal R-matrix parameters (energies, reduced width
amplitudes, etc.) derived from the Input File. Primarily a check file with
limited direct use.

param.sav
^^^^^^^^^

Contains the best-fit formal R-matrix parameters after a fit is completed. This
file can be loaded back into AZURE2 to reproduce a fit or to use as starting
parameters for subsequent calculations (extrapolation, reaction rate, etc.).

parameters.out
^^^^^^^^^^^^^^

Contains the **physical** (observable) parameters resulting from the fit. If
the user wishes to use these as new starting values, they must be manually
entered into the **Levels and Channels** tab.

normalizations.out
^^^^^^^^^^^^^^^^^^

Contains the fitted normalization factors for data segments where normalization
was varied. This file is automatically loaded when ``param.sav`` is selected.

Cross Section Output
--------------------

AZUREOut_aa=\*_R=\*.out
^^^^^^^^^^^^^^^^^^^^^^^^

Output from **Calculate With Data** and **Fit With Data** modes. The filename
encodes the entrance (``aa``) and exit (``R``) particle pair indices.

Nine columns:

.. list-table::
   :widths: 10 90
   :header-rows: 1

   * - Col.
     - Description
   * - 1
     - Center-of-mass energy (MeV)
   * - 2
     - Excitation energy (MeV)
   * - 3
     - Center-of-mass angle (degrees)
   * - 4
     - Fit center-of-mass cross section (barns or barns/sr)
   * - 5
     - Fit center-of-mass S-factor (MeV b or MeV b/sr)
   * - 6
     - Data center-of-mass cross section (barns or barns/sr)
   * - 7
     - Data center-of-mass cross section uncertainty (barns or barns/sr)
   * - 8
     - Data center-of-mass S-factor (MeV b or MeV b/sr)
   * - 9
     - Data center-of-mass S-factor uncertainty (MeV b or MeV b/sr)

When multiple segments share the same entrance and exit particle pairs, their
data are written to the same file in the order they appear in the **Segments**
tab, separated by a double blank line.

AZUREOut_aa=\*_R=\*.extrap
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Output from **Calculate Segments Without Data** mode. Same naming convention
as above. Five columns:

.. list-table::
   :widths: 10 90
   :header-rows: 1

   * - Col.
     - Description
   * - 1
     - Center-of-mass energy (MeV)
   * - 2
     - Excitation energy (MeV)
   * - 3
     - Center-of-mass angle (degrees)
   * - 4
     - Extrapolated center-of-mass cross section (barns or barns/sr)
   * - 5
     - Extrapolated center-of-mass S-factor (MeV b or MeV b/sr)

Uncertainty and Statistics
--------------------------

chiSquared.out
^^^^^^^^^^^^^^

Contains the :math:`\chi^2/N` values for each data segment after a calculation
or fit. The final line gives the total :math:`\chi^2`.

param.errors
^^^^^^^^^^^^

Contains the reduced width amplitudes and their asymmetric uncertainties from a
MINOS error analysis.

covariance_matrix.out
^^^^^^^^^^^^^^^^^^^^^

Contains the covariance and correlation matrices from a MINOS calculation,
providing a complete description of parameter correlations.

Other Files
-----------

intEC.dat
^^^^^^^^^

External capture integral values for data segments. This file can be reused to
speed up subsequent calculations, as long as:

- No calculation segments have been added or removed.
- No levels of a new :math:`J^\pi` have been added or removed.
- No channels have been added or removed.
- The channel radius has not changed.

Level energies, widths, and ANCs can be changed freely while reusing this file.

intEC.extrap
^^^^^^^^^^^^

Same as ``intEC.dat``, but for the calculation (extrapolation) segments.

reactionrates.dat
^^^^^^^^^^^^^^^^^

Contains temperatures (in GK) and calculated reaction rates
(in cm\ :sup:`3` mol\ :sup:`-1` s\ :sup:`-1`) from the **Calculate Reaction
Rate** mode.

samples.mcmc
^^^^^^^^^^^^

Contains the MCMC chain samples from Bayesian parameter inference (if the MCMC
module is enabled).
