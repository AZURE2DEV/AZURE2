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
particle pairs. The dialog displays retrieved levels with their energies and
spin-parities, allowing you to select which levels to import.

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

.. warning::

   The user must carefully check that the channel configuration is set correctly.
   If the maximum orbital momentum is too low, channels may be missing, leading
   to incorrect results or crashes. The code requires **at least one channel**
   for each particle pair for each :math:`J^\pi`.

Dummy Levels
^^^^^^^^^^^^

Because the number of levels and angular momenta must be truncated, the code
only includes those :math:`J^\pi` values that are present in the defined levels.
To include all necessary hard-sphere phase shifts for scattering and external
capture calculations, the user may need to add **dummy levels**:

1. Create a level at an arbitrary energy with the missing :math:`J^\pi`.
2. Fix its energy.
3. Set all partial widths to zero.

The effects of including these dummy levels should be investigated for each
allowed spin-parity combination.

Gamma-Ray Calculations
----------------------

For gamma-ray calculations, the user should define a level corresponding to each
gamma-ray particle pair defined in the **Particle Pairs** tab. This allows
setting the ANCs and gamma-ray decay widths for each bound state. The energies
of these levels must match the excitation energies of the gamma-ray particle
pairs.
