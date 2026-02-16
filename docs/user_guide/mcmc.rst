MCMC Tab
========

The **MCMC** tab provides Bayesian parameter inference via Markov Chain Monte
Carlo sampling. This feature must be enabled at compile time.

.. note::

   The MCMC module uses an affine-invariant ensemble sampler for efficient
   exploration of parameter space.

Parameters
----------

The **Parameters** sub-tab shows all fit parameters with the following columns:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - **Name**
     - Parameter identifier.
   * - **Value**
     - Starting value for the MCMC chain.
   * - **Prior Mean**
     - Center of the Gaussian prior (if enabled).
   * - **Prior Std**
     - Width of the Gaussian prior (if enabled).
   * - **Use Gaussian Prior**
     - Check to impose a Gaussian prior on this parameter.

Use the following buttons to populate the table:

- **Load Physical Parameters** -- load parameter values in physical units.
- **Load RWA Parameters** -- load parameter values as reduced width amplitudes.

Sampling Settings
-----------------

Configure the MCMC sampler:

**Number of Walkers**
   The number of independent walkers in the ensemble. Recommended:
   2--10 times the number of free parameters.

**Number of Steps**
   Total number of MCMC steps to perform.

**Chain Spread (%)**
   The initial spread of walkers around the starting values, as a percentage.

**Number of Threads**
   Number of parallel threads for the calculation.

**Fresh Start**
   When checked, start a new chain from scratch. When unchecked, resume from
   an existing ``samples.mcmc`` file.

**Use Reduced Widths**
   When checked, sample in reduced width amplitude space rather than physical
   parameter space.

Running MCMC
------------

Switch to the **Progress** sub-tab to monitor the sampling:

- **Progress bar** and **status label** show overall progress.
- **Current iteration**, **log probability**, **log likelihood**, and **log
  prior** are updated in real time.
- Runtime output is displayed in the text area.
- Click **Stop** to halt the sampling.

Results
-------

The **Results** sub-tab displays summary statistics after sampling completes:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Statistic
     - Description
   * - **Mean**
     - Mean value of the parameter across all samples.
   * - **Std Dev**
     - Standard deviation of the samples.
   * - **2.5%**
     - Lower bound of the 95% credible interval.
   * - **50%**
     - Median value.
   * - **97.5%**
     - Upper bound of the 95% credible interval.

Click **Refresh Results from File** to recalculate statistics from an existing
``samples.mcmc`` file (useful for post-processing or when resuming a previous
run).

Output
------

MCMC results are saved to the ``samples.mcmc`` file in the output directory.
This file can be loaded by external tools for further analysis (e.g., corner
plots, convergence diagnostics).
