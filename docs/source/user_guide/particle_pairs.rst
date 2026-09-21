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

Identical Particles
-------------------

A (Particle, Particle) pair whose two members have the same Z, mass, spin and
parity, with the heavy one in its ground state, is treated as a pair of
identical particles (:math:`\alpha+\alpha`, p+p, d+d, :math:`^3`\ He+\ :math:`^3`\ He,
:math:`^{12}`\ C+\ :math:`^{12}`\ C, ...). Nothing has to be switched on.

Allowed channels
   Exchange symmetry admits only channels with :math:`\ell + s` **even**, for
   bosons and fermions alike: 1S0, 3P\ :sub:`0,1,2`, 1D2, ... for p+p, and
   even :math:`\ell` only for a spin-0 pair. For a pair with spin, a channel
   with :math:`\ell + s` odd (3S1 in p+p, say) is refused when the model is
   read, with a message naming it; for a spin-0 pair it draws a warning.

Elastic scattering
   The differential cross section is symmetrized in each channel spin
   :math:`s`: the Coulomb amplitude is
   :math:`f_C(\theta) + (-1)^s f_C(\pi-\theta)`, and every nuclear pathway
   carries :math:`1 + (-1)^{\ell'+s'}` (a factor 2 on allowed channels). The
   spin-averaged Coulomb cross section is then the Mott formula for spin
   :math:`j`,

   .. math::

      \frac{d\sigma}{d\Omega} = \left(\frac{\eta}{2k}\right)^2 \left[
      \frac{1}{\sin^4\frac{\theta}{2}} + \frac{1}{\cos^4\frac{\theta}{2}}
      + \frac{(-1)^{2j}}{2j+1}\,
      \frac{2\cos\!\left(\eta \ln \tan^2\frac{\theta}{2}\right)}
           {\sin^2\frac{\theta}{2}\cos^2\frac{\theta}{2}} \right],

   with interference weight +1 for :math:`\alpha+\alpha`, :math:`-1/2` for
   p+p and :math:`+1/3` for d+d. The same symmetrization enters the
   polarization observables (isDiff 7 and 8), so for p+p
   :math:`A_y(\pi-\theta) = -A_y(\theta)`. The angle-integrated elastic cross
   section is reported per collision (half the integral over the full sphere)
   for any spin.

What is not symmetrized
   Reactions into or out of an identical pair (d+d :math:`\to` p+t, say) use
   only the channel rule above; no exchange factor is applied to their cross
   sections, the same convention as for spin-0 pairs. The vector analyzing
   power is defined for a spin-1/2 beam only (it is zero for d+d).
