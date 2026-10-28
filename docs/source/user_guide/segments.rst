Segments Tab
============

The **Segments** tab defines the energies and angles at which cross sections are
calculated. It is divided into two halves:

- **Segments From Data** (upper) -- for calculations that reference experimental
  data files, and measured scattering lengths (listed under the data segments).
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

**Entrance and Exit Pairs**
   Select the entrance and exit particle pairs from the drop-down lists. Each
   item starts with the pair's key (as assigned in the **Particle Pairs** tab);
   exit items read as the reaction with the chosen entrance pair. When the
   project has at least one gamma (capture) pair, the exit list ends with
   **Total capture (all γ pairs)**: the angle-integrated capture cross section
   summed over every gamma pair (see **Data Type** below).

**Energy Range**
   The **Low Energy** and **High Energy** fields select a range of data from the
   data file. To include all data, set the range to cover the entire file. Energies
   are in the laboratory frame, in MeV.

**Angle Range**
   For differential cross sections and analyzing powers, specify the **Low
   Angle** and **High Angle** (in degrees, in the frame of the data) to select a
   subset of the data. The box title names the frame: **C.M. Angle** for a
   C.M. differential cross section and an analyzing power, **Lab Angle**
   otherwise. For angle-integrated or phase-shift data types, these fields are
   disabled.

**Data Type**
   Select the observable from the drop-down menu. Two of them have a second
   drop-down next to them:

   .. list-table::
      :widths: 30 70

      * - **Angle Integrated**
        - Angle-integrated cross section data. With a gamma (capture) exit
          pair, a **Multipolarity** drop-down appears: **All** (default),
          **E1 only** or **E2 only**. The E1/E2 choices compare the data with
          only that multipolarity's part of the calculated capture cross
          section. With **Total capture (all γ pairs)** as the exit pair, the
          segment is the angle-integrated total capture cross section, summed
          over all gamma-ray transitions; each significant gamma-ray cascade
          transition must be defined as a separate particle pair. Total capture
          has no E1/E2 variant, so the multipolarity drop-down is hidden for it.
      * - **Differential**
        - Differential cross section data at specific angles. A **Frame**
          drop-down chooses **Lab** (default) or **C.M.**: the frame in which
          the data file gives its angles and cross sections. Both work for
          particle and capture exit pairs. Lab-frame capture data are converted
          to the c.m. frame with the photon aberration and solid-angle Jacobian
          (:doc:`../reference/data_formats`, "Frame"): a few percent and angle
          dependent for :sup:`12`\ C(:math:`\alpha,\gamma`) at a few MeV.
          Most low-energy capture measurements publish lab values.
      * - **Phase Shift**
        - Phase shift data. Requires specifying the total angular momentum (*J*)
          and orbital angular momentum (*l*), and, when the target has spin so
          that several channel spins *s* share that (*J*, *l*), the channel
          spin (the **s** field, written as ``spin <s>`` on the line): a phase
          shift belongs to one (*J*, *l*, *s*) partial wave, and an ambiguous
          segment is refused with the spins listed. The convention is
          :math:`-90° < \theta_\text{lab} < 90°`.
      * - **Analyzing Power**
        - Vector analyzing power :math:`A_y` for a spin-1/2 projectile (any
          other projectile spin is refused), in the Madison convention, for a
          **particle or a capture** exit channel.
          Angles are centre-of-mass, as for a C.M. differential cross section,
          and the data file carries ``E_lab  theta_cm  A_y  dA_y``. Because
          :math:`A_y` is a ratio, **Vary Norm?** is disabled for it -- a
          normalization factor has no meaning for a quantity that is already
          normalized. See :doc:`../theory/polarization_theory`.
      * - **Angular Coefficients**
        - Legendre coefficients of a reaction's centre-of-mass angular
          distribution, :math:`d\sigma/d\Omega_\text{cm} = A_0\,[1 + \sum_k
          a_k P_k(\cos\theta_\text{cm})]`, as data: the file carries
          ``E_lab  k  a_k  da_k`` with the Legendre order :math:`k` in the
          angle column, and the fit value of each point is the model's
          coefficient of that order -- the same coefficients a test segment of
          type **Angular Distribution** computes. Legendre orders are computed
          up to 20, so a particle channel with :math:`l > 10` truncates every
          differential cross section it enters (the run log warns for each
          such channel; angle-integrated values are complete). The trailing keyword
          ``absolute 1`` makes the data the absolute :math:`A_k` in b/sr
          (:math:`A_0 = \sigma/4\pi`). Target integration and resolution
          folding act on the absolute :math:`A_k`, which are linear in the
          cross section, and the ratio is formed afterwards. A normalised
          coefficient has no scale, so **Vary Norm?** is refused for it; an
          elastic pair (the expansion leaves the Coulomb amplitude out), total
          capture, composites and UPOS segments are refused. Energy shifts and
          the analytic gradient work as for any segment. Since October 2026
          (for the 11B(α,p)14C coefficients of Wang 1991 and Dayras 1976).

   While **Total capture (all γ pairs)** is the exit pair, only **Angle
   Integrated** can be chosen; the other observables are greyed out until a
   single exit pair is selected again.

   The choices map onto the observable code (``isDiff``) of the segment line in
   the ``.azr`` file, which is unchanged:

   .. list-table::
      :header-rows: 1
      :widths: 55 15 30

      * - Dialog choice
        - Code
        - Exit pair key written
      * - Angle Integrated (Multipolarity All, or a particle exit pair)
        - 0
        - the exit pair
      * - Differential, Frame Lab
        - 1
        - the exit pair
      * - Phase Shift
        - 2
        - the exit pair
      * - Angle Integrated, exit **Total capture (all γ pairs)**
        - 3
        - ``-1``
      * - Differential, Frame C.M.
        - 4
        - the exit pair
      * - Angle Integrated, gamma exit, Multipolarity E1 only
        - 5
        - the exit pair
      * - Angle Integrated, gamma exit, Multipolarity E2 only
        - 6
        - the exit pair
      * - Analyzing Power
        - 7
        - the exit pair
      * - Angular Coefficients
        - 9
        - the exit pair
      * - Analyzing Power
        - 7
        - the exit pair

   Editing a segment restores these controls from its code.

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

   **Shared normalization (normalization groups, input file only).** Segments
   that come from one measurement often share one absolute scale -- an
   excitation function and the angular distributions it was integrated from,
   or one data set split into per-angle segments. Ending their lines with the
   keyword block ``normgroup <label>`` makes all segments with the same label
   use the normalization of the **first** of them in file order -- its value,
   vary flag, uncertainty and penalty -- and the other members' own
   normalization, uncertainty and vary flag are ignored, e.g.::

      1 1 3 0 2.3 0 180 0 1 1 10 0.016 0.01 1 data/Dayras_ap.dat 0 0 normgroup Dayras
      1 1 3 1.572 1.572 0 180 4 1 1 0 0.030 0.01 1 data/Dayras_ang.dat 0 0 normgroup Dayras
      1 1 3 2.032 2.032 0 180 4 1 1 0 0.081 0.05 1 data/Dayras_ang.dat 0 0 normgroup Dayras

   fits one normalization for the three segments with one 10 % penalty.
   A member whose own line varies its normalization keeps its
   ``segment_<key>_norm`` entry in ``param.sav`` (fixed, written with the
   shared value), so parameter files and parameter numbering do not depend
   on the groups. ``chiSquared.out`` lists every member's (shared)
   normalization and the penalty on the first member only. The run log notes
   each group and its members. The keyword may be combined with
   ``shiftgroup``; writers put it before a ``sample`` block. The GUI preserves
   the block but does not edit it, and its Fitting tab still lists a member's
   own vary flag.

   **Correlated uncertainties (covariance matrices).** A
   line ending with the keyword block ``covariance <file>`` attaches its
   points to the data covariance matrix in ``<file>``; all segments naming
   the same file share it, and their points enter :math:`\chi^2` as
   :math:`r^T M^{-1} r` (:doc:`chi_squared`; file format in
   :doc:`../reference/data_formats`), e.g.::

      1 1 3 0 8 14.5 15.5 1 1 0 0 0 0 0 data/n0.dat 0 0 normgroup ALL covariance cov_n.txt
      1 1 3 0 8 29.5 30.5 1 1 0 0 0 0 0 data/n0.dat 0 0 normgroup ALL covariance cov_n.txt

   Rows of the matrix are the data rows of the files it lists, so one file
   may be split into several segments (here by angle window) and the energy
   and angle limits may change without a new matrix. The block may be
   combined with ``normgroup``, ``shiftgroup`` and the energy-shift blocks,
   and only data segments take it. In the GUI it is the **Covariance File**
   field of the data segment dialog (below).

**Covariance File**
   Optional. A data covariance file (the ``covariance <file>`` block above):
   every data segment that names the same file shares one covariance matrix.
   Type the path or use **Choose...**; the path is resolved like the data
   file's. Leave the field empty for uncorrelated errors (no block is
   written). Editing a segment shows the file of its line, and clearing the
   field removes the block.

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

   When the segments of one experiment need their own constant shifts (each
   angle digitized from a separate figure carries its own offset) but share
   one energy-dependent calibration, use ``shiftfngroup <label>`` instead of
   ``shiftgroup``: every segment with the label applies the function declared
   on the first of them, with that segment's parameters, and keeps its own
   constant :math:`\Delta E`, value, uncertainty and vary flag::

      1 3 3 0 10  90  90 4 1 1 0 -0.0047 0.01 1 data/pp.dat 0 0 shiftfngroup H shiftfn "a0*(sqrt(x)-sqrt(1.5))" 1 0 0.02 1
      1 3 3 0 10 125 125 4 1 1 0 -0.0009 0.01 1 data/pp.dat 0 0 shiftfngroup H
      1 3 3 0 10 141 141 4 1 1 0 -0.0117 0.01 1 data/pp.dat 0 0 shiftfngroup H
      1 3 3 0 10 165 165 4 1 1 0 -0.0183 0.01 1 data/pp.dat 0 0 shiftfngroup H

   fits four constants and one :math:`\sqrt{E}` coefficient; each
   constant is then that segment's shift at :math:`E_0` = 1.5 MeV. A segment
   may be in a ``shiftgroup`` or a ``shiftfngroup``, not both, and only the
   first member declares the ``shiftfn``. Members that are single-energy
   angular distributions add nothing to the coefficient: their own constants
   absorb it.

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
   :math:`(a_j/\delta a_j)^2`, measured from zero; the value on the line is
   only the starting value (:doc:`chi_squared`). The function is evaluated at every point and every
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

   Frames: a lab-frame differential **Sum** (Data Type "Differential", lab
   angles) is the sum of the terms' *lab* cross sections. Each term is
   evaluated in the c.m. frame of its own pair, so a component to another exit
   pair (a different Q value, or another compound nucleus of a natural target)
   carries its own lab-to-c.m. Jacobian; the engine applies the ratio of the
   base pair's Jacobian to the term's when it adds the term (since October
   2026; before, the c.m. values were added as they were, which for
   :math:`^{13}`\ C(α,n\ :sub:`1`\ +n\ :sub:`2`) at 6 MeV was −3 % at 0° and
   +6 % at 150°). A c.m. **Sum**, an angle-integrated one, a component of the
   same pair, or gamma terms of one entrance pair are unaffected. A target
   effect on a composite integrates each term at that term's own c.m. angle
   (its sub-points follow the component's conversion; before October 2026 a
   component to another pair, or a fixed-angle **Ratio** denominator, was
   integrated at the base pair's angle).

   A lab-frame differential **Ratio** is the ratio of the two terms' *lab*
   cross sections: each term is taken to the lab with the Jacobian of its own
   pair and exit type (photon or particle), so a denominator to another exit
   pair, at a fixed angle, or of the other kind (a capture yield over an
   elastic one) is compared in the right frame. Before October 2026 the
   denominator's Jacobian was computed with the base segment's pairs: right
   for a denominator of the same pair and for another particle pair at a
   general angle without ``relativistic=1``, but 3 % off for another exit pair
   with ``relativistic=1`` or at 0° / 180°, and 30 % off for a capture over an
   elastic ratio. A c.m. or angle-integrated **Ratio** carries no frame factor.
   The lab output (``kinflag.dat`` = 1) of a **Ratio** test segment prints the
   ratio itself.

**Unobserved Primary, Observed Secondary (UPOS)**
   For a particle exit pair :math:`b+B` whose residual nucleus :math:`B` is
   left in an excited state, the segment describes the **secondary gamma ray**
   of :math:`B \to C + \gamma` rather than the particle :math:`b`, which is not
   detected (Brune and deBoer, Phys. Rev. C **102**, 024628 (2020)). The angle
   is that of the gamma ray relative to the beam; the calculated quantity is
   :math:`d\sigma/d\Omega_\gamma` of Eq. (20) of that paper (100 % branching to
   :math:`C`; scale the data or the normalization for any other branching), and
   an angle-integrated UPOS segment is the cross section to :math:`B`.
   **Secondary Decay L** is the multipolarity :math:`L` of the gamma ray,
   **Final State J (Ic)** the spin of :math:`C`, and **Mixing Ratio (δ)**
   :math:`\delta = g_{L+1}/g_L` of an :math:`L`/:math:`L+1` mixture (0 for a
   pure multipole), in the convention of Eqs. (4)-(7) of the paper; the
   mixed radiation parameter :math:`(R_k^{LL} + 2\delta R_k^{LL'} + \delta^2
   R_k^{L'L'})/(1+\delta^2)` is applied at every order :math:`k` (before
   October 2026 it was skipped at an order where :math:`R_k^{L'L'}` happened
   to vanish, such as :math:`k = 2` of a :math:`3/2 \to 3/2` M1/E2 decay,
   which dropped the interference term and the normalization there). In the
   file the four fields follow the data file (and any composite fields) as
   ``isUPOS L Ic delta``, e.g. ``... data/n2g.dat 0 1 3 0 0`` for the E3 decay
   of :math:`^{16}`\ O(6.13, :math:`3^-`) to its :math:`0^+` ground state.

   When the exit pair has more than one channel spin (both :math:`b` and
   :math:`B` with nonzero spin), the gamma-ray angular distribution contains
   interference terms between the channel spins, so it depends on the
   relative signs of their reduced-width amplitudes; those signs are defined
   with the channel spin coupled light particle first,
   :math:`(I_b I_B m_b m_B | s'\nu')`, as for the analyzing powers. Eq. (22)
   of the paper as printed corresponds to the opposite order, which changes
   the sign of these terms only. Engines before the 2026-10 fix of the
   channel-spin pairing computed these interference terms incorrectly (see
   ``tests/secondary_gamma``).

   Give UPOS data in the centre-of-mass frame (data type differential c.m.):
   the gamma-ray angle in the rest frame of :math:`B` differs from the
   laboratory angle only by the aberration of order :math:`v_B/c`, while a
   laboratory-frame segment converts its angles with the kinematics of the
   particle :math:`b`, which does not apply to a gamma ray. An extrapolation
   (Segments Without Data) cannot be UPOS; to draw a UPOS curve, give a data
   segment a grid of pseudo-data.

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

Natural and Enriched Targets (Sample Segments)
----------------------------------------------

When the particle pairs form more than one compound nucleus
(:doc:`particle_pairs`), both segment dialogs offer a **Sample**
drop-down: *No*, the default sample, or any named sample of the Pairs tab.
A sample segment is the abundance-weighted sum: its own term, and each
component of a **Sum** composite, is multiplied by the isotopic abundance of
the compound nucleus its *entrance pair* forms in the chosen sample, as
entered in the Pairs tab's abundance table (or fitted, when the abundance is
a fit parameter). The row's reaction label is prefixed with the sample in brackets, *[default sample]* or
*[name]*. Data sets taken on different targets simply name different
samples.

For natural neon under a neutron beam, the total cross section is one Sum
composite with the components n+\ :sup:`20`\ Ne, n+\ :sup:`22`\ Ne and
n+\ :sup:`21`\ Ne (elastic exit pairs, one per nucleus) and the box ticked;
the engine then computes
:math:`0.9048\,\sigma(\mathrm{n}+{}^{20}\mathrm{Ne}) + 0.0925\,\sigma(\mathrm{n}+{}^{22}\mathrm{Ne}) + 0.0027\,\sigma(\mathrm{n}+{}^{21}\mathrm{Ne})`.
A plain (non-composite) sample segment is just its own cross section times
its nucleus's abundance, which is what an isotopically enriched target with a
known enrichment needs.

The box is only available for additive observables: angle-integrated cross
sections (including total capture and the E1/E2-only forms) and differential
cross sections in either frame, for data and test segments alike. It is cleared and disabled for a **Ratio**
composite, a phase shift, an analyzing power and a test segment of
angular-distribution coefficients. The entrance and exit pair of every term
must form the same compound nucleus. In the file the flag is the keyword
``sample 1`` (default sample) or ``sample <name>`` at the end of the
segment line (after a ``sqrtshift`` block, if any), and the engine refuses it
in a project without an ``<abundances>`` block or with a sample name the
block does not define. Each sample segment also gets an output file
``AZUREOut_sample_seg=<key>.out`` (``.extrap`` for a test segment) with the
total and one column per nucleus (:doc:`../reference/output_files`).

Neutron Total Cross Sections (``total 1``)
-------------------------------------------

A neutron total cross section can be given as a **Sum** composite of the
elastic channel and every reaction channel of the compound nucleus, but it is
cheaper, and cannot miss a channel, to take it from the optical theorem: the
keyword ``total 1`` at the end of the segment line makes the segment's term the
total cross section of its entrance pair,

.. math::

   \sigma_\mathrm{tot} = \frac{\pi}{k^2}\,\frac{1}{(2i_1+1)(2i_2+1)}
   \sum_J (2J+1) \sum_{c} 2\,\mathrm{Re}\,T_{cc}
   = \frac{2\pi}{k^2}\, g \sum_J (2J+1) \sum_c (1 - \mathrm{Re}\,U_{cc}),

summed over the entrance channels :math:`c`. Every channel of a
:math:`J^\pi` group is in its level matrix, so :math:`U` is unitary and
:math:`\sum_{c'} |T_{cc'}|^2 = 2\,\mathrm{Re}\,T_{cc}` over all exit
channels: the result equals the sum over every particle exit pair of the
compound nucleus (to :math:`10^{-11}` on the 13C+alpha neutron data), from one
T-matrix per point instead of one per exit pair. Photon channels have no width
in the level matrix and external capture is outside it, so capture is not
included; where it matters, keep the capture pair as a **Sum** component (a
component whose exit pair differs from its entrance pair stays what it is).

The segment must be angle-integrated (data type 0), its exit pair its entrance
pair, and the pair neutral (the total cross section of charged particles
diverges) and not of identical particles; the R-matrix-for-capture formalism
is refused. In a **Sum** composite every term whose exit pair is its entrance
pair is a total, so a natural-target total is one line, e.g.::

   1 2 2 0 2.33 0 180 0 1 0 5 0 0 0 data/ntot.dat 1 0 2 4 4 -999 6 6 -999 0 sample 1 total 1

for :math:`\sigma_\mathrm{tot}(\mathrm{n}+{}^{20}\mathrm{Ne})`,
:math:`\sigma_\mathrm{tot}(\mathrm{n}+{}^{22}\mathrm{Ne})` and
:math:`\sigma_\mathrm{tot}(\mathrm{n}+{}^{21}\mathrm{Ne})` weighted by the
abundances. The keyword applies to test segments (Segments Without Data) as
well. The GUI keeps it through a load and save but has no control for it yet.

.. warning::

   An engine older than October 2026 does not know the keyword and ignores
   it silently: it computes the elastic cross section alone.

Since October 2026 the terms of a **Sum**, **Ratio** or total-capture segment
also share one level-matrix solution per energy (each term still gets its own
T-matrix), in the calculation and in the analytic gradient; the results are
unchanged.

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
   The observable list and the **Total capture (all γ pairs)** exit item work
   as for data segments, with two differences: there is no **Multipolarity**
   drop-down (test segments have no E1/E2-only codes), and the list adds
   **Angular Distribution Coefficients**. In full:

   - **Angle Integrated** -- code 0; with the **Total capture (all γ pairs)**
     exit item, the total capture cross section, code 4 (exit key ``-1``).
   - **Differential** -- with **Frame** **Lab** (code 1) or **C.M.** (code 5).
   - **Phase Shift** -- code 2.
   - **Angular Distribution Coefficients** -- output Legendre polynomial
     coefficients, code 3. Requires specifying the polynomial order.
   - **Analyzing Power** -- the vector analyzing power on the chosen energy and
     angle grid, in the centre-of-mass frame, code 7.

.. warning::

   Errors may occur if you specify an energy or angular range that is not
   kinematically allowed, or if the energy is too low and the penetrability
   becomes vanishingly small.

Kinematic limits of the lab frame
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The lab-to-c.m. conversion of a differential segment (data or test, code 1)
covers the forward branch only, and the lab angle is always that of the light
ejectile (particle 1 of the exit pair):

- Where two c.m. angles reach the same lab angle (a reaction near the
  threshold of an endothermic exit, e.g. 7Li(p,n)7Be below about 2 MeV,
  or elastic scattering with a projectile heavier than the target), AZURE2
  models the forward c.m. angle. A lab angle beyond the kinematic maximum
  of that point's energy is refused: ``ERROR: ... beyond the kinematic
  maximum``. Before October 2026 such a point became NaN silently (one line
  in the log for a reaction, none for elastic scattering) and the run went
  on with ``nan`` in ``chiSquared.out`` and the extrapolation files.
  A lab-differential point below the threshold of its exit pair is refused
  for the same reason (no ejectile, nothing to convert).
- Recoil detection is not an option: data taken on the heavy recoil must be
  converted by hand (elastic scattering: :math:`\theta_{c.m.} = 180^\circ -
  2\theta_{recoil}`) and given as c.m. (code 5 / 4).
- Analyzing-power angles (code 7) are c.m.; lab-angle :math:`A_y` data need
  converting first (the value needs no Jacobian, the angle does).
- A UPOS secondary-gamma segment given as lab differential is converted with
  the kinematics of particle b (see UPOS below); the run says so with a
  warning.

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

Measured Scattering Lengths
---------------------------

AZURE2 computes the s-wave scattering length :math:`a_0` of every particle pair
at the pair's threshold after each calculation or fit and writes it to
``output/scatteringLengths.out``. Measured values enter the fit as the data
type **Scattering Length** of **Segments From Data**: each active one adds
:math:`((a_0 - \text{value})/\sigma)^2` to :math:`\chi^2`, so a scattering
length counts like one more data point. The definition and how the R-matrix
parameters determine :math:`a_0` are in :doc:`../theory/scattering_length`.

Click **+** under **Segments From Data** and choose the Data Type **Scattering
Length**. The dialog then shows **Pair** (the particle pair; gamma and beta pairs
have no channel spins), **Channel Spin**, the channel spin :math:`s`, which for
:math:`l = 0` is the :math:`J^\pi` of the channel, chosen from the values the
pair allows (:math:`|j_1 - j_2|` to :math:`j_1 + j_2`, only even ones for a pair
of identical particles), or **coherent**, the spin-weighted combination
:math:`\sum_s (2s+1)\,a_s / ((2j_1+1)(2j_2+1))` that a neutron measurement on
an unpolarized target usually gives (not defined for identical particles);
**Re a0** and its uncertainty in fm; and, optionally, **Im a0** and its
uncertainty for a pair with open reaction channels at its threshold. An
uncertainty of 0 leaves that part unconstrained; blank imaginary fields mean no
imaginary part. The other fields of a data segment (energies, angles, data file,
normalization, shifts) do not apply and are hidden.

The scattering lengths are listed as rows under the data segments, marked
**a**\ :math:`_0` instead of a row number: the reaction column names the pair,
the energy column reads *threshold*, the angle column the channel spin and the
data-file column the measured value. The check box, **Check All**, **Uncheck
All** and the pair filters act on them as on the segments; select one and click
**-** to remove it, or double-click it to edit it (its type cannot be changed to a
segment, nor a segment's to a scattering length). They are not segments: they
have no segment number, so adding or removing one does not renumber the
segments that **Experimental Effects** and the parameter names refer to.
Before October 2026 they had a table of their own at the bottom of the tab.

A row refers to its pair by the pair's key, the row number on the **Particle
Pairs** tab. When a pair is deleted or edited, rows whose pair no longer exists,
is no longer a particle pair, or no longer allows the row's channel spin are
not removed: they turn red, and the tooltip says why. They are saved as they
are, and the calculation stops with an error naming the line until they are
corrected. They are saved as a ``<scatteringLengths>`` block after
``</segmentsTest>``, and only when there is at least one; comments in a block
written by hand are not kept when the GUI saves the file.
