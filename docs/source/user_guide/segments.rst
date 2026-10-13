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
   :math:`\sqrt{E}` term is not shared. The run log notes each group and its
   members. The GUI preserves the block but does not edit it.

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

Natural-abundance and compound targets (mixture segments)
---------------------------------------------------------

Neutron total cross sections, and other data taken on natural or compound
samples, are abundance-weighted sums over the isotopes present, each with its
own compound nucleus. An AZURE2 project describes one compound nucleus, so the
other isotopes enter a segment as fixed **tables**: the segment's model becomes

.. math::

   \sigma_\text{model}(E) = w_\text{self}\,\sigma_\text{R-matrix}(E)
   + \sum_i w_i\, T_i(E),

where :math:`w_\text{self}` is the abundance of this project's target isotope
(times its stoichiometry in a compound), and each :math:`T_i` is a tabulated
cross section of another isotope -- an ENDF/JEFF point-wise table, or the
extrapolation of another AZURE2 project -- with its own weight. For natural
neon fitted as :math:`^{20}\text{Ne}+n`, :math:`w_\text{self} = 0.9048` and the
:math:`^{22}\text{Ne}` and :math:`^{21}\text{Ne}` tables carry 0.0925 and
0.0027. The weights are plain numbers, not fit parameters.

The block is written at the end of the segment line (after any
``sqrtshift`` block), on a ``<segmentsData>`` or a ``<segmentsTest>`` line::

   ... data/Junghans_ntotal.dat 1 0 2 2 1 -999 2 3 -999 0 mixture 0.9048 2 "tables/22Ne_ntot.dat" 0.0925 "tables/21Ne_ntot.dat" 0.0027

i.e. ``mixture <w_self> <n> {"<file>" <w_i>} x n``. The example is a neutron
total written, as usual, as a **Sum** composite of the elastic exit and every
other exit pair; the mixture applies to the combined value. Each table is a
two-column text file, ``E_lab [MeV]  sigma [b]`` in strictly increasing lab
energy of the projectile (the axis the data file uses), ``#`` comments
allowed, interpolated linearly; its cross section is in the frame AZURE2
compares the segment in (barns for an angle-integrated or total cross section,
c.m. b/sr for a differential one). See :doc:`../reference/data_formats`.

Where the resolution function acts
   The tables are folded with the segment's **own** experimental effect --
   Gaussian beam spread, time-of-flight (UDR) kernel, beam profile or target
   integration from ``<targetInt>`` -- on the same sub-point grid and with the
   same kernel as the R-matrix part, so the whole sum is broadened together,
   as a transmission measurement sees it. The tables being
   parameter-independent, this is done once at initialization (and again when
   an energy shift moves the points), so a mixture costs nothing in the fit.
   A table must cover every energy the grid reaches, including the kernel's
   wings beyond the data range and the segment's nominal energy shift;
   otherwise the run stops and names the range needed. Give tables a grid
   fine enough to resolve the other isotope's resonances.

Outputs
   The fit column of ``AZUREOut_*.out`` / ``.extrap`` is the mixture total,
   the value compared with the data. A file
   ``AZUREOut_mixture_seg=<key>.out`` (or ``.extrap``) gives its
   decomposition per point: the energy, the total, :math:`w_\text{self}`
   times the R-matrix part, and each folded table term.

Fitting
   Derivatives scale by :math:`w_\text{self}` and the tables add nothing to
   them, in the analytic gradient and in the pyazr Jacobian alike. The usual
   normalization and energy-shift parameters apply to the whole sum.

Limits
   Supported on angle-integrated and differential cross sections (plain or
   **Sum** composites); refused on ratios, phase shifts, analyzing powers,
   total-capture and E1/E2-only segments, and on segments whose experimental
   effect uses an automatic-application tolerance. In the GUI the block is
   edited in the segment dialogs (Segments tab, add or double-click a row):
   tick **Natural-abundance or compound target**, set the abundance of this
   nucleus and add one row per table file with its weight; the sum of the
   weights is shown, and the dialog refuses the observables the engine refuses.
   Any other keyword block on the line (``sqrtshift``) is carried along
   unchanged. With pyazr use ``AzrModel.set_mixture``, ``set_test_mixture``,
   and ``pyazr.extrap_to_lab_table`` to turn another project's extrapolation
   into a table.

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
