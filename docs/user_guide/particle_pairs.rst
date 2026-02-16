Particle Pairs Tab
==================

The **Particle Pairs** tab is the first tab that should be filled out when
setting up a new calculation. It defines the reaction particle pairs -- the two
particles that fuse to form, or result from the decay of, the compound nucleus.

.. tip::

   Always start here. AZURE2 uses the particle pair information to automatically
   calculate allowed channels in the **Levels and Channels** tab and to populate
   options in the **Segments** tab.

Managing Particle Pairs
-----------------------

- Click the **+** button in the lower-left corner to add a new particle pair.
- Select a pair and click **-** to delete it.
- Double-click an existing pair to edit it.

Each particle pair is automatically assigned a numerical key (displayed on the
far left), which is referenced when creating segments.

.. important::

   The **first** particle pair must always be of type **(Particle, Particle)**.
   After the first, pairs may be defined in any order.

Particle Pair Types
-------------------

AZURE2 supports three types of particle pairs:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Type
     - Description
   * - **(Particle, Particle)**
     - Standard nuclear reaction pair (e.g., p + :sup:`14`\ N). Must be the
       first pair defined.
   * - **(Particle, Gamma)**
     - Radiative capture pair for gamma-ray transitions to bound states.
       Some fields are auto-populated and cannot be edited. Gamma-ray decays
       are limited to bound states.
   * - **(Beta Decay)**
     - Beta-delayed particle emission pair. Some fields are auto-populated.

Add Particle Pair Dialog
------------------------

When adding or editing a particle pair, a dialog window appears with the
following fields:

Particle Properties
^^^^^^^^^^^^^^^^^^^

For both the **light** and **heavy** particle:

Spin (J)
   The spin of the particle. Half-integer values should be entered as decimals
   (e.g., ``0.5`` for spin-1/2).

Parity
   The parity of the particle, selected from the drop-down menu: **+** or **-**.

Proton Number (Z)
   The number of protons in the particle.

Mass (M)
   The mass number of the particle.

Pair Properties
^^^^^^^^^^^^^^^

Excitation Energy
   The excitation energy of the heavy particle (in MeV). This is non-zero for
   transitions to excited states rather than the ground state.

Separation Energy
   The energy required to separate the compound system (in its ground state) into
   the constituent particle pair (in MeV).

Channel Radius
   The R-matrix channel radius (in fm). Different particle pairs need not share
   the same radius. A useful estimate is:

   .. math::

      R = R_0 \times (A_p^{1/3} + A_t^{1/3})

   .. warning::

      The channel radius is a model parameter, not a physical nuclear radius.
      Always test the sensitivity of your results to the channel radius by
      re-running the calculation with different values. The channel radius
      **cannot** currently be used as a fit parameter.

External Capture Multipolarities
   For **(Particle, Gamma)** pairs only. Select **E1** and/or **E2**
   multipolarities. The code automatically determines the allowed intrinsic and
   angular momentum combinations based on the defined resonances.
