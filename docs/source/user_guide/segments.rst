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
          Madison convention, for a **particle or a capture** exit channel.
          Angles are centre-of-mass, as for **C.M. Differential**, and the data
          file carries ``E_lab  theta_cm  A_y  dA_y``. Because :math:`A_y` is a
          ratio, **Vary Norm?** is disabled for it -- a normalization factor has
          no meaning for a quantity that is already normalized. See
          :doc:`../theory/polarization_theory`.

.. note::

   A measured **outgoing polarization** :math:`P` of :math:`A(a,b)B` is, by
   time reversal, the analyzing power of the inverse reaction
   :math:`B(b,a)A`, so it is entered as an **Analyzing Power** segment on the
   inverse channel (entrance :math:`b+B`, exit :math:`a+A`). Data published
   as the product :math:`P\,d\sigma/d\Omega` must be divided by the
   differential cross section first.

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

**sqrt(E) energy shift** (input file only)
   A second, optional shift term that grows with the beam energy. With it a
   segment's energies are shifted as

   .. math::

      E' = E + \Delta E + b\,\sqrt{E/\mathrm{MeV}}

   where :math:`\Delta E` is the constant **Energy Shift** above and :math:`b`
   (in MeV\ :sup:`1/2`) is the new coefficient. An additive offset in the
   field reading of an analyzing magnet (:math:`E = kB^2`) produces exactly
   this :math:`\sqrt{E}` dependence, whereas the constant term describes a
   fixed offset; a calibration that is right at low energy and drifts at high
   energy needs the second form. The two terms are independent: each has its
   own value, uncertainty and vary flag, and a segment may use either, both
   or neither.

   The term is written at the very end of the segment's ``<segmentsData>``
   line, after the data file and any composite/UPOS fields, as the keyword
   block ``sqrtshift <b> <bError> <vary>`` (``vary`` 1 or 0), e.g.::

      1 2 3 0 10 0 0 4 1 1 10 0 0.02 1 data/yield_0deg.dat 0 0 sqrtshift 0 0.01 1

   frees :math:`b` from 0 with a penalty width of 0.01 MeV\ :sup:`1/2`
   (:doc:`chi_squared`). A line without the block behaves exactly as before.
   The fit parameter is named ``segment_<key>_energy_shift_sqrt`` in
   ``param.sav``, and ``shifts.out`` reports it in its last column. The GUI
   preserves the block when it loads and saves a file but does not edit it.

   **Shared shift (shift groups).** A data set split into several segments
   (one per angle, say) usually has one beam-energy calibration, so its
   segments should share one shift rather than each fitting its own. Ending
   the lines with the keyword block ``shiftgroup <label>`` does that: all
   segments with the same label use the constant shift of the **first** of
   them in file order -- its value, uncertainty, vary flag and penalty -- and
   the other members' own shift value, uncertainty and vary flag are ignored,
   e.g.::

      1 2 1 0 3.567 157 161 4 1 1 10 0 0.02 1 data/pa0.dat 0 0 shiftgroup Frawley_pa0
      1 2 1 0 3.567 148 153 4 1 1 10 0 0.02 1 data/pa0.dat 0 0 shiftgroup Frawley_pa0

   fits one shift for both, with one 20 keV penalty. Each member keeps its
   ``segment_<key>_energy_shift`` entry in ``param.sav`` (fixed, written with
   the shared value), so parameter files stay compatible. The
   :math:`\sqrt{E}` term is not shared (use a shift function for that, below).
   The run log notes each group and its members. The GUI preserves the block
   but does not edit it.

**Energy-shift function** (input file only)
   Different calibration errors shift the energy scale in different ways: an
   offset in an analyzing magnet's field reading as :math:`\sqrt{E}`, an error
   in its calibration constant or in a time-of-flight path length as
   :math:`E`, a time-of-flight time-zero offset as :math:`E^{3/2}`. Instead of
   a fixed keyword per form, a segment may carry a general function,

   .. math::

      E' = E + \Delta E + f(x;\,a_0, a_1, \ldots),

   with :math:`x` the point's original laboratory energy in MeV and every
   :math:`a_j` a fit parameter with its own value, uncertainty and vary flag.
   The keyword block, at the end of the line like the others, is::

      shiftfn "<f(x)>" N v0 s0 vary0  v1 s1 vary1 ...

   with the expression in double quotes (blanks inside are ignored), ``N``
   the number of parameters and one value, uncertainty (penalty width) and
   vary flag per parameter. The expression is read like the ``<targetInt>``
   equations (:doc:`experimental_effects`): ``x``, parameters ``a0``, ``a1``,
   ..., numbers, ``+ - * / ^`` and ``sqrt``, ``exp``, ``ln``, ``log``,
   ``sin``, ``cos``, ``tan``, ``asin``, ``acos``, ``atan``.

   On a shift group the block goes on the group's **first** segment and the
   function applies to every member with one set of parameters -- one
   calibration for one experiment. Four Henderson :sup:`14`\ C(p,p) angles
   sharing one magnet-offset calibration::

      1 3 3 0 10  90  90 4 1 1 0 -0.0047 0.01 1 data/pp.dat 0 0 shiftgroup H shiftfn "a0*(sqrt(x)-sqrt(1.7))" 1 0 0.02 1
      1 3 3 0 10 125 125 4 1 1 0  0      0    0 data/pp.dat 0 0 shiftgroup H
      1 3 3 0 10 141 141 4 1 1 0  0      0    0 data/pp.dat 0 0 shiftgroup H
      1 3 3 0 10 165 165 4 1 1 0  0      0    0 data/pp.dat 0 0 shiftgroup H

   fit two parameters for the four segments: the shared constant shift
   (10 keV penalty) and the :math:`\sqrt{E}` coefficient (0.02
   MeV\ :sup:`1/2`).

   *Centre the function inside the data.* The constant :math:`\Delta E` and
   its penalty refer to the energy where :math:`f` vanishes. Written as
   ``a0*sqrt(x)``, that is :math:`E = 0`, far below the data, and the 10 keV
   penalty then constrains the shift at zero energy rather than any shift the
   data see; written as ``a0*(sqrt(x)-sqrt(1.7))`` it is the shift at
   1.7 MeV, inside the data. (The ``sqrtshift`` block cannot be centred; use a
   function where this matters.) Forms by mechanism, with :math:`E_0` a number
   you choose inside the data:

   .. list-table::
      :header-rows: 1

      * - Calibration error
        - Form
        - ``shiftfn``
      * - Analyzing-magnet field offset (:math:`E = kB^2`)
        - :math:`\propto \sqrt{E}`
        - ``"a0*(sqrt(x)-sqrt(E0))"``
      * - Magnet calibration constant :math:`k`
        - :math:`\propto E`
        - ``"a0*(x-E0)"`` (:math:`a_0 = \Delta k/k`)
      * - Time-of-flight path length :math:`L`
        - :math:`\propto E`
        - ``"a0*(x-E0)"`` (:math:`a_0 = 2\Delta L/L`)
      * - Time-of-flight time zero
        - :math:`\propto E^{3/2}`
        - ``"a0*(x^1.5-E0^1.5)"``

   Choose the form from the apparatus, not from :math:`\chi^2`.

   The parameters are named ``segment_<key>_energy_shift_fn<j>`` in
   ``param.sav`` (``<key>`` the segment declaring the function) and come after
   every other parameter, so files without a function keep their parameters
   exactly. Each free :math:`a_j` with a nonzero uncertainty adds a penalty
   (:doc:`chi_squared`). The function is evaluated at every point and every
   target-integration or convolution sub-point at its own original energy;
   the target thickness and stopping powers stay those of the unshifted
   energy, as for every shift (:doc:`experimental_effects`).

   AZURE2 stops at start-up for a follower that declares a function (put it
   on the group's first segment), for a function together with ``sqrtshift``
   on any member (two descriptions of one calibration; write the
   :math:`\sqrt{E}` term into the function), and for an expression that uses
   an undeclared parameter, leaves a declared one unused, or has a free term
   with no effect at the group's energies. It warns when the free terms --
   :math:`\Delta E` if varied, and each varied :math:`a_j` -- cannot be told
   apart over the group's energies (for example a :math:`\sqrt{x}` term on a
   single-energy angular distribution, where it only duplicates
   :math:`\Delta E`). The run log has one line per function with its
   members, parameters and that condition number, and ``shifts.out`` one
   ``#`` line per function with the shift at the group's lowest and highest
   energies and its penalty. The GUI preserves the block but does not edit
   it.

   A ``param.sav`` value sets an energy-shift parameter (constant,
   :math:`\sqrt{E}` or function) but never frees or fixes it: shifts follow
   the vary flag of their segment line, because their penalties do.

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
   y-axis to a linear scale and to *Cross Section* automatically: a
   logarithmic axis cannot display a quantity that goes negative, and it drops
   such points silently rather than complaining, so the curve would look
   merely sparse.
