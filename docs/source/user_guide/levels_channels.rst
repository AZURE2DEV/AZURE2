Levels and Channels Tab
=======================

.. note::

   This page describes the ``.azr`` content through the fields of the Qt setup
   utility of the released versions. This version has no setup utility;
   :doc:`graphical_interface` says how to do the same from the command line and
   with ``pyazr``.

The **Levels and Channels** tab is used to enter the compound nucleus level
information. It uses the particle pair definitions from the **Particle Pairs**
tab to automatically calculate all allowed R-matrix channels.

Tab Layout
----------

The tab is divided into four frames:

1. **Compound Nucleus Levels** (left) -- the list of nuclear levels.
2. **Channels in Selected Level** (center) -- the allowed channels for the
   currently selected level.
3. **Channel Configuration** (top right) -- controls to limit the number of
   R-matrix channels.
4. **Channel Details** (bottom right) -- detailed information for the selected
   channel, including the width parameter input field.

Managing Levels
---------------

Adding a Level
^^^^^^^^^^^^^^

Click the **+** button in the lower-left corner. A dialog appears prompting for:

- **Compound Nucleus** -- only shown when the particle pairs form more than one
  (:doc:`particle_pairs`): which nucleus the level belongs to. The J\ :sup:`π`
  list follows the choice.
- **Excitation Energy** -- the excitation energy of the compound nucleus (in MeV).
- **Spin (J)** -- the total spin. Half-integer values should be given as decimals
  (e.g., ``1.5``).
- **Parity** -- select **+** or **-** from the drop-down menu.

After clicking **Accept**, the new level appears in the list and its allowed
channels are automatically calculated from the particle pair information.

Importing Levels from IAEA Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Click the **NDS** button to query the IAEA Nuclear Data Services for known level
schemes. AZURE2 automatically determines the compound nucleus from the defined
particle pairs (and asks which one when the pairs form several). The dialog
displays retrieved levels with their energies and spin-parities, allowing you
to select which levels to import.

Editing a Level
^^^^^^^^^^^^^^^

Select a level and double-click to modify it. To delete a level, select it
and click the **-** button.

Level Controls
^^^^^^^^^^^^^^

Include / Exclude
   The checkbox in the **Include?** column controls whether a level participates
   in the calculation. Uncheck to exclude a level without deleting it.

Fix Energy
   The checkbox in the **Fix?** column prevents the level energy from varying
   during a fit.

Other Buttons
   - **Fix All Widths** -- fix all channel width parameters at once.
   - **Fix All Energies** -- fix all level energies at once.
   - **Export LaTeX** -- export the level scheme as a LaTeX table for publication.

Levels of Several Compound Nuclei
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the particle pairs form more than one compound nucleus
(:doc:`particle_pairs`), every level belongs to one of them: the one chosen
in the dialog, or, for a level read from a file, the nucleus of its channels'
pair. The table gains a **Compound Nucleus** column and a drop-down above it
that shows one nucleus's levels at a time. A level's channels are built from
its own nucleus's pairs only, the missing-J\ :sup:`π` report is given per
nucleus, and two levels with the same J\ :sup:`π` and energy in different
nuclei are different levels. A level cannot be moved to another nucleus; add
it again there. The file keeps one block of level lines per nucleus; the
engine refuses a level whose channels name two nuclei, unless an
``isospin_level`` record declares it one isospin level (below).

Isospin Levels Across Compound Nuclei
-------------------------------------

Several compound nuclei can share isospin-labelled levels, as in Hale's
charge-independent R-matrix analysis of the A = 4 system (Tilley, Weller and
Hale, Nucl. Phys. A541, 1 (1992)). Such a level is one level ID whose
``<levels>`` lines span every nucleus it lives in (each line with that
nucleus's own excitation energy), declared by the records of ``<isospin>``
(:doc:`../reference/format3`)::

   <isospin>
   isospin_reference ref_4Li          # 4Li: energies measured from p+3He
       nucleus    4Li
       threshold  p+3He
   isospin_shift shift_4H             # 4H: its threshold pair and internal Coulomb shift
       nucleus    4H
       threshold  n+3H
       shift      -0.86 MeV fit prior -0.86 0.3
   isospin_shift shift_4He
       nucleus    4He
       threshold  p+3H
       shift      -0.64 MeV fit prior -0.64 0.3
   isospin_level T1                   # the level with ID 1 has T = 1
       level      1
       T          1
   isospin_level T0                   # the level with ID 2 has T = 0, soft-tied
       level      2
       T          0
       soft_energy    50 keV
       soft_amplitude 0.02
   isospin_mix dd_T1                  # level 1, channel d+d, l = 1, s = 1: an amplitude of its own
       level      1
       pair       d+2H
       l          1
       s          1
       amplitude  0.0 fit prior 0 0.1
   isospin_factor f                   # optional: in place of a computed Clebsch-Gordan factor
       level      1
       pair       p+3H
       factor     0.7071
   isospin_pair t                     # optional: the particle isospins of a pair
       pair       p+12C[15.11]
       t_light    1/2
       t_heavy    1
   </isospin>

- In each nucleus, a level's reduced-width amplitude in a channel of pair
  (a, b) is :math:`\langle t_a t_{3a}\, t_b t_{3b} | T T_3\rangle` times one
  amplitude per isospin channel (arrangement, l, s). Particle isospins default
  to :math:`t = |N-Z|/2`, :math:`t_3 = (Z-N)/2`; pairs whose particles match in
  mass number, spin, parity, isospin and state form one arrangement, and
  partner channels (same l, s and J) must share a channel radius (a
  ``channel_radius`` record may set it per l, s or J, the same for every pair
  of the arrangement). A particle's state is the rank of its excitation energy
  among the particles of the project with the same Z, A, spin and parity (the
  lowest listed is the first): mirror partners are ranked each in its own
  nucleus, so the n-th listed state of a spin and parity pairs with the n-th,
  and a nucleus that is its own mirror pairs each state with itself. The
  start-up report and ``isospin.out`` mark a state after the first, e.g.
  ``0+#2`` for :math:`^{12}`\ C(7.654) beside the 0+ ground state.
- The member in the reference nucleus holds the free parameters. The others
  are derived: energies
  :math:`E_m = E_\mathrm{ref} - T_\mathrm{ref} + T_m + \Delta_m - \Delta_\mathrm{ref}`
  (:math:`T` the threshold pair's excitation energy, :math:`\Delta_m` the
  nucleus's shift, the parameter ``<record>.shift`` of its ``isospin_shift``
  record; the reference has none), amplitudes the factor ratio times the free
  one. A channel with factor 0 is zero unless an ``isospin_mix`` record gives
  it an amplitude (the parameter ``<record>.amplitude``). Derived
  slots stay in the parameter list, fixed, and are written synced
  (``param.sav``); the two kinds of parameter follow the level parameters.
- The two soft widths relax one level's tie: its derived slots become free
  parameters with the penalty :math:`((x - x_\mathrm{tie})/\sigma)^2` each
  (forbidden channels centred on 0), reported with the priors of the shifts and
  mixes as ``Total-Isospin-Chi-Squared`` (``chiSquared.out``) and in
  ``isospin.out``; they start on their tie values. A soft level takes no
  ``isospin_mix`` record. Every fitter carries the penalty (the least-squares
  fitters as residual rows), and the sampler has it in its likelihood.
- The ties act on the fitted parameters of the run's basis. In the formal
  basis (``--no-brune``, Hale's formulation) these are the formal amplitudes
  and eigenenergies, and every J-group holding an isospin level uses
  :math:`B_c = -l_c`, so partners share their boundary values. In the Brune
  basis (the default) they are the Brune (on-resonance) amplitudes and the
  observed energies, and no boundary value enters. For a level alone in its
  J-group the Brune amplitudes equal the formal ones, so the amplitude ties are
  the same statement; with other levels of the same :math:`J^\pi` in the group
  (background poles) the two differ through level mixing. The energy relation
  differs in either case: an observed energy carries its nucleus's boundary
  (Thomas-Ehrman) shift,
  :math:`\tilde E = E_\lambda - \sum_c \gamma_c^2 (S_c(\tilde E) - B_c)`, which
  depends on the level, so one shift :math:`\Delta_m` per nucleus describes the
  formal energies better than the observed ones. In a Brune run give the
  energies a soft tie with a wide ``soft_energy``, or fix them at their measured
  values, unless one level per nucleus is tied. Park's basis (``--use-park``)
  and the RMC formalism are refused. The analytic gradient, the least-squares
  fitters and the sampler work through the ties (a derived slot's derivative
  passes to the free slots it follows); the sampler and a session take
  parameters in the fit's basis only, since the transformation from observed
  widths would mix the derived members into the free ones.
- Photon and beta channels of a J-group that holds an isospin level stay
  outside the ties: Hale's relations are between the particle reduced-width
  amplitudes of isospin partners, and an electromagnetic or weak amplitude is
  not charge independent. Each nucleus keeps its own gamma widths (and feeding
  amplitudes), free or fixed as its channel lines say; the start-up log counts
  them per member as "not tied". Capture to a tied bound state uses that
  member's own energy and particle amplitude, so in a tied model the
  asymptotic normalization of a derived member follows the reference nucleus's
  amplitude.
- The start-up log prints every factor; ``isospin.out`` lists the fitted
  (formal, or Brune) parameters of every member, and ``parameters.out`` labels
  the members.

``AZURE2 --convert`` turns the ``<isospin>`` block of the ``dev`` branch
(``reference``, ``shift``, ``level``, ``mix``, ``factor``, ``pairT``, ``soft``
lines) into these records, and ``AZURE2 --convert-params`` the names of its
parameter files (``isospin_shift_<Z>_<A>``, ``isospin_mix_<k>``) into the
records' (``shift_4H.shift``, ``mix_1.amplitude``).

Channels
--------

When a level is selected, its allowed channels appear in the **Channels in
Selected Level** frame. Each channel is characterized by:

- The particle pair it belongs to
- The channel spin (*S*)
- The orbital angular momentum (*L*)
- For gamma-ray channels: the radiative decay type

Selecting a channel displays its details in the **Channel Details** frame,
where the width parameter can be entered.

Width Parameters
^^^^^^^^^^^^^^^^

The width parameter entered at the bottom of the Channel Details frame
represents:

- **Partial width** for unbound levels (channels above the separation energy).
- **ANC** (Asymptotic Normalization Coefficient) for bound-state channels
  (automatically determined by the code).

.. tip::

   Relative interference signs between channels are set by changing the **sign**
   of the partial width or ANC.

Fix Width
   The **Fix?** checkbox in the channels list prevents the width parameter from
   varying during a fit. A width given as zero in the project file is held
   fixed in every parameter basis: also with ``--no-brune``, where its formal
   amplitude is a mixture of the other levels' and need not be zero, and with
   ``--no-ignore-externals``, where a zero photon width has the internal
   amplitude :math:`-\mathrm{Re}\,\gamma^{\rm ext}`. Give a non-zero start
   value to fit it. (Builds before 2026-11-08 tested the transformed amplitude,
   so in the formal basis rounding decided which of these were free.)

.. note::

   In the project file each level has one line per channel. The engine matches
   each line's width to its channel by (particle pair, *l*, *s*), so the lines of
   a hand-edited level may come in any order. A level that has no line for a
   channel which another level of its :math:`J^\pi` group lists gets a zero,
   fixed width for that channel. Both cases print a ``WARNING: level J = ...``
   line. Before 2026-10-17 the widths were stored in line order, so a reordered
   or shortened level silently got the wrong widths.

   Lines belong to one level when they carry the same level ID (the ninth
   column) within 1 keV, or the identical energy. Two levels of one
   :math:`J^\pi` closer than 1 keV with different IDs (a narrow level on top of
   a broad one) stay two levels, with a ``WARNING: level ID ... is within 1 keV``
   line. Before 2026-10-24 any lines within 1 keV were read as one level
   whatever their IDs: the later level's widths replaced the earlier one's and
   its parameters disappeared from the fit, with only a "lists channel ... twice"
   warning.

   The levels of a :math:`J^\pi` group may also be listed in any energy order
   (a background pole first, say). Before 2026-10-27 this mattered in the
   standard (non-Brune) parametrization: the transformation to formal
   parameters handed its ascending eigenvalues to the levels in file order, so
   the first-listed level carried the lowest formal energy, ``parameters.out``
   listed the group in ascending order rather than the file's, and a fit that
   fixed the pole's energy varied the low level instead. The cross sections
   were never affected (the set of levels is the same whatever the labels).
   The Brune parametrization, the default, never had the problem. Since
   2026-11-07 the boundary condition of the standard parametrization is set at
   the lowest level of the group, not the first listed, so the formal
   parameters themselves no longer depend on the listing order.

Channel Configuration
---------------------

The **Channel Configuration** frame (top right) controls which channels are
generated:

Maximum Orbital Momentum (L_max)
   The maximum orbital angular momentum quantum number. Default is 2.
   Must be set high enough for all relevant channels to appear.

Maximum Gamma Multipolarity
   The maximum electromagnetic multipolarity for gamma-ray channels. Default is 2.

Maximum Multipolarities per Decay
   The maximum number of multipolarities to consider per gamma-ray decay. Default
   is 2.

Include hard-sphere scattering for every :math:`J^\pi`
   On by default for a new project. AZURE2 builds one R-matrix per
   :math:`J^\pi`, and historically a :math:`J^\pi` entered the calculation only
   when a level carried it. A :math:`J^\pi` without a level was simply absent,
   so an elastic-scattering calculation lost the potential (hard-sphere)
   scattering of that partial wave -- a mistake nothing warned about. With this
   option every :math:`J^\pi` the particle pairs allow up to the maximum orbital
   momentum scatters, whether or not a level carries it: a :math:`J^\pi` with
   no level contributes its hard-sphere phase shifts and nothing else. The line
   under the box names the :math:`J^\pi` in that situation.

   Switching a level off then removes its resonance and keeps the partial wave,
   which is the natural decomposition into resonant and potential scattering;
   switching every level off gives the pure hard-sphere plus Coulomb cross
   section. To recover the old convention, in which a :math:`J^\pi` disappears
   with its levels, untick the box.

   The setting is stored in the project file (``includeAllJPi`` and
   ``maxOrbitalMomentum`` in the ``<config>`` block), so a command-line or
   ``pyazr`` run of the same file agrees with the GUI. A project saved before
   this option existed opens with it off and computes exactly as it always did.

.. warning::

   The user must carefully check that the channel configuration is set correctly.
   If the maximum orbital momentum is too low, channels may be missing, leading
   to incorrect results or crashes. The code requires **at least one channel**
   for each particle pair for each :math:`J^\pi`.

Dummy Levels
^^^^^^^^^^^^

Before the option above existed, the way to include the hard-sphere phase
shifts of a :math:`J^\pi` that had no level was a **dummy level**: a level at an
arbitrary (usually high) energy with that :math:`J^\pi`, its energy fixed and
all partial widths set to zero. Such levels are easily confused with background
poles, which do carry widths. With the option on, dummy levels are harmless (a
zero-width level adds nothing to the level sum, also at a data point exactly at
its energy since 2026-10-31; before, such a point lost its whole
:math:`J^\pi` group) and can be deleted; the
``tests/hard_sphere`` check of the source tree verifies that a model with a
dummy level per :math:`J^\pi` and the same model with the option on give
identical cross sections.

Gamma-Ray Calculations
----------------------

For gamma-ray calculations, the user should define a level corresponding to each
gamma-ray particle pair defined in the **Particle Pairs** tab. This allows
setting the ANCs and gamma-ray decay widths for each bound state. The energies
of these levels must match the excitation energies of the gamma-ray particle
pairs: the pair's final state is the level of its :math:`J^\pi` whose energy
lies within 1 keV of the pair's excitation energy. A pair with no such level
is refused when a segment uses it (warned about otherwise), and one whose
window holds two or more levels of that :math:`J^\pi` is refused as
ambiguous (before October 2026 the first one in the file was taken silently):
move the extra level by more than 1 keV, or give the pair the energy of the
level it decays to.

**A final state above a particle threshold is a narrow resonance.** A final
state of external capture that lies above the threshold (:math:`S` plus the
pair's excitations) of a particle pair with a channel in its :math:`J^\pi`
group is unbound in that channel. Its exterior is then the outgoing wave
:math:`O_l(k_f r)/O_l(k_f a)`, :math:`O = G + iF`, at :math:`E_f` above the
threshold, in place of the Whittaker function of a bound state; its strength
is the level's observed width in that channel, and the normalization
:math:`N_f` is continued as :math:`(\hbar^2/2\mu a)\,dS/dE`, which joins the
bound value at the threshold. The run says so in a ``NOTE``. This is right
while the width is small next to the photon energy and next to the distance
to the threshold. For a broad final state give the capture data set or
extrapolation the field ``line_window`` (the half-width :math:`W` of the
gamma-line window of the analysis, e.g. ``line_window 39 keV``): its cross
section is then integrated over the final state's line shape for
:math:`E_f` within :math:`E_R \pm W`, with the photon energy
:math:`E_x - E_x(E_f)` and the scattering wave of the final channel at each
:math:`E_f`; the line is centred on the final level's current energy, so a
fitted level energy moves it. Below the final state there is no photon
(:math:`E_\gamma = E_x - E_f \le 0`): the capture is zero there, for internal
and external capture alike.

The analytic derivatives do not run over the line shape, so a project with a
``line_window`` refuses everything that reads them. The program refuses, at
set-up and with no output, a fit with a fitter that reads the analytic
gradient or the Jacobian of the residual rows (``--use-gradient``,
``--use-lm``: the columns ``gradient`` and
``residual_rows`` of ``pyazr.kinds()["fitters"]``); fit with plain MIGRAD
(no fitter switch), whose gradient is numerical. A pyazr session raises
``pyazr.EngineError`` for ``Objective.value_and_grad``, ``Objective.grad``,
``Objective.jacobian``, ``Objective.residuals_and_jacobian``,
``Session.model_gradients``, the uncertainty bands made from them and
``fit_with_engine`` with those fitters; the values (``Objective(x)``,
``Objective.residuals``, ``Objective.log_prob``, ``Session.evaluate``) work.
A run that asks for the cross-section uncertainty band writes its outputs
without the band and says why.

One case is refused (``E_FINAL_STATE_UNBOUND``, by ``AZURE2 --validate``, by
every run and by pyazr): a final state with external capture at such a
threshold, or so close above the threshold of a charged pair that its
Sommerfeld parameter exceeds 100 (within about a keV for light nuclei), where
the outgoing wave cannot be evaluated. Lower the state below the threshold
(it is then bound) or move it higher.
