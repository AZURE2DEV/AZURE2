Segments Tab
============

The **Segments** tab defines the energies and angles at which cross sections are
calculated. It is divided into two halves:

- **Segments From Data** (upper) -- for calculations that reference experimental
  data files.
- **Segments Without Data** (lower) -- for pure calculations and extrapolations.

Each segment is automatically assigned a numerical key that is referenced in the
**Experimental Effects** tab.

Segments From Data
------------------

These segments are used with the **Calculate Segments From Data** and **Fit
Segments From Data** calculation modes. Cross sections are computed only at the
energies and angles present in the experimental data file.

Creating a Data Segment
^^^^^^^^^^^^^^^^^^^^^^^

Click the **+** button in the lower-left corner of the upper frame. A dialog
appears with the following fields:

**Entrance and Exit Pair Keys**
   Select the entrance and exit particle pairs by their numerical keys (as
   assigned in the **Particle Pairs** tab). Values can be entered directly or
   set using the spinner controls.

**Energy Range**
   The **Low Energy** and **High Energy** fields select a range of data from the
   data file. To include all data, set the range to cover the entire file. Energies
   are in the laboratory frame, in MeV.

**Angle Range**
   For differential cross sections, specify the **Low Angle** and **High Angle**
   (in degrees, laboratory frame) to select a subset of the data. For angle-integrated
   or phase-shift data types, these fields are disabled.

**Data Type**
   Select from the drop-down menu:

   .. list-table::
      :widths: 30 70

      * - **Angle Integrated**
        - Angle-integrated cross section data.
      * - **Differential**
        - Differential cross section data at specific angles.
      * - **Phase Shift**
        - Phase shift data. Requires specifying the total angular momentum (*J*)
          and orbital angular momentum (*l*). The convention is
          :math:`-90° < \theta_\text{lab} < 90°`.
      * - **Total Capture**
        - Angle-integrated total capture cross section, summed over all gamma-ray
          transitions. Each significant gamma-ray cascade transition must be
          defined as a separate particle pair.
      * - **C.M. Differential**
        - Differential cross section data given in the center-of-mass frame.
      * - **Analyzing Power**
        - Vector analyzing power :math:`A_y` for a spin-1/2 projectile, in the
          Madison convention. Angles are centre-of-mass, as for **C.M.
          Differential**, and the data file carries
          ``E_lab  theta_cm  A_y  dA_y``. Because :math:`A_y` is a ratio,
          **Vary Norm?** is disabled for it -- a normalization factor has no
          meaning for a quantity that is already normalized. See
          :doc:`../theory/polarization_theory`.

**Data Normalization**
   A normalization factor applied to the data yield. Default is ``1.0``.

**Normalization Error (%)**
   The relative systematic uncertainty of the data set, entered as a percentage.
   This is included in the chi-squared calculation using the D'Agostini method.

**Vary Norm?**
   Check this box to allow the normalization factor to be varied as a free
   parameter during fitting.

**Energy Shift**
   An energy calibration shift applied to the data (in MeV).

**Vary Energy Shift?**
   Check this box to allow the energy shift to be varied during fitting.

**Data File**
   Path to the experimental data file. Use the **Choose...** button to browse,
   or enter the path directly. Paths can be absolute or relative to the Input
   file directory.

**Advanced Mode**
   When checked, enables advanced operations for combining multiple reaction
   channels:

   - **Sum** -- sum cross sections from multiple entrance/exit pair combinations.
   - **Ratio** -- compute the ratio of cross sections from different channels.

Editing and Managing Data Segments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **Double-click** a segment to edit it.
- Select a segment and click **-** to delete it.
- Use the **arrow buttons** on the right to reorder segments.
- Use **Check All** / **Uncheck All** to quickly enable or disable all segments.
- Use the **entrance** and **exit pair filter** dropdowns to show only segments
  for specific particle pairs.

.. tip::

   Data from different experiments should be placed in separate files. Multiple
   segments can reference the same data file with different energy or angle
   ranges to create separate excitation curves or angular distributions.

.. important::

   When including systematic uncertainties, carefully consider the segmentation
   scheme to avoid double-counting uncertainties across segments that share the
   same data.

Segments Without Data
---------------------

These segments define energy and angle ranges for pure calculations (no
experimental data), and are used with the **Calculate Segments Without Data**
mode. Typical uses include:

- Interpolating a cross section with finer energy steps than the data provides.
- Extrapolating to astrophysically relevant energies.
- Calculating a total cross section from a differential fit, or vice versa.
- Computing the cross section for the inverse reaction.
- Extracting angular distribution coefficients.

Creating a Test Segment
^^^^^^^^^^^^^^^^^^^^^^^

Click the **+** button in the lower half. The dialog is similar to the data
segment dialog but includes:

**Energy Range with Step Size**
   Specify **Low Energy**, **High Energy**, and **Energy Step** (all in MeV).
   For a single energy point, set low and high to the same value and step to zero.

**Angle Range with Step Size**
   Same as energy: specify **Low Angle**, **High Angle**, and **Angle Step**
   (in degrees).

**Data Types**
   In addition to the standard types, this mode supports:

   - **Angular Distribution Coefficient** -- output Legendre polynomial
     coefficients. Requires specifying the polynomial order.
   - **Analyzing Power** -- the vector analyzing power on the chosen energy and
     angle grid, in the centre-of-mass frame.

.. warning::

   Errors may occur if you specify an energy or angular range that is not
   kinematically allowed, or if the energy is too low and the penetrability
   becomes vanishingly small.

.. note::

   An analyzing power measured on a **thin** target should be given a segment
   with no target integration. :math:`A_y` averaged over a thick target is
   weighted by the cross section, and since Rutherford scattering diverges at
   low energy where :math:`A_y` is essentially zero, a thick target drives the
   average towards zero. This is physical, and it is explained in
   :doc:`../theory/polarization_implementation`.

   When an analyzing-power segment is plotted, the Plot tab switches the
   y-axis to a linear scale and to *Cross Section* automatically: a logarithmic
   axis cannot display a quantity that goes negative, and an S-factor
   conversion is meaningless for a ratio.
