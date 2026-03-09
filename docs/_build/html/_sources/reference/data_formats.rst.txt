Data Formats
============

Input Data Files
----------------

Experimental data files are plain text files with four whitespace-delimited
columns of real numbers (decimal or scientific notation):

.. list-table::
   :widths: 15 85
   :header-rows: 1

   * - Column
     - Description
   * - 1
     - Forward-kinematic laboratory frame **energy** (MeV)
   * - 2
     - Forward-kinematic laboratory frame **angle** (degrees)
   * - 3
     - Forward-kinematic laboratory frame **cross section** (barns or barns/sr)
   * - 4
     - Forward-kinematic laboratory frame **cross section uncertainty** (barns or barns/sr)

Example data file:

.. code-block:: text

   0.500  90.0  1.234e-03  5.6e-05
   0.600  90.0  2.345e-03  7.8e-05
   0.700  90.0  3.456e-03  9.0e-05
   0.800  90.0  4.567e-03  1.2e-04

Notes on Data Files
^^^^^^^^^^^^^^^^^^^

- **Angle column**: must always be present, even for angle-integrated data. In
  that case the angle is a dummy value and is not used in the calculation.
- **Sorting**: data may be entered in any order. However, for plotting purposes
  it is useful to sort by energy (for excitation curves) or by angle (for
  angular distributions).
- **Frame**: all input data must be in the **laboratory frame** with
  **forward kinematics** (light particle as projectile, heavy particle as
  target).
- **Units**: energies in MeV, angles in degrees, cross sections in barns
  (angle-integrated) or barns/sr (differential).
- **Format**: columns can be separated by spaces or tabs.

Phase Shift Data
^^^^^^^^^^^^^^^^

Phase shift data follows the same four-column format, but the third column
contains the phase shift value (in degrees) and the fourth column contains its
uncertainty. The convention is :math:`-90° < \theta_\text{lab} < 90°`.

Temperature Files
^^^^^^^^^^^^^^^^^

For reaction rate calculations, a temperature file is a single-column text file
listing temperatures in GK (gigakelvin):

.. code-block:: text

   0.01
   0.05
   0.10
   0.50
   1.00
   2.00
   5.00
   10.0

Project Files (.azr)
--------------------

The AZURE2 Input File (``.azr``) stores the complete project configuration in a
structured text format. While the format is not formally documented, it is a
plain text file that advanced users can modify directly for batch processing.

The file contains sections for:

- Configuration settings (formalism, directories, runtime options)
- Particle pair definitions
- Level and channel information
- Data segment definitions
- Test segment definitions
- Experimental effects configuration
- Fitting parameter settings
- MCMC settings (if enabled)

.. tip::

   The recommended practice is to always use the GUI to create and edit
   ``.azr`` files. The internal format may change between versions.
