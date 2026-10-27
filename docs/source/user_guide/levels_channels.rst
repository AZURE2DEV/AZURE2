Levels and Channels Tab
=======================

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
engine refuses a level whose channels name two nuclei.

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
   varying during a fit. A width left at zero is automatically held fixed.

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
   The Brune parametrization, the default, never had the problem.

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
zero-width level adds nothing to the level sum) and can be deleted; the
``tests/hard_sphere`` check of the source tree verifies that a model with a
dummy level per :math:`J^\pi` and the same model with the option on give
identical cross sections.

Gamma-Ray Calculations
----------------------

For gamma-ray calculations, the user should define a level corresponding to each
gamma-ray particle pair defined in the **Particle Pairs** tab. This allows
setting the ANCs and gamma-ray decay widths for each bound state. The energies
of these levels must match the excitation energies of the gamma-ray particle
pairs.
