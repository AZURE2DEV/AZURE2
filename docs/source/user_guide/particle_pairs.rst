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
   analyzing power (isDiff 7), so for p+p
   :math:`A_y(\pi-\theta) = -A_y(\theta)`. The angle-integrated elastic cross
   section is reported per collision (half the integral over the full sphere)
   for any spin.

Reactions out of an identical pair
   An identical pair in the **entrance** channel multiplies every cross
   section out of it by :math:`1+\delta_{12} = 2` -- reaction, capture and
   elastic, angle-integrated and differential:

   .. math::

      \sigma_{12\to34} = \frac{\pi}{k^2} \sum_J
      \frac{2J+1}{(2i_1+1)(2i_2+1)} \,(1+\delta_{12})
      \sum_{\ell s,\,\ell' s'} \left|T^J_{\ell' s',\ell s}\right|^2 ,

   the sum running over the channels symmetry allows; for a spin-0 pair this
   is the familiar :math:`(\pi/k^2)\sum_\ell (2\ell+1)[1+(-1)^\ell]\,T_\ell`
   of :math:`^{12}`\ C+\ :math:`^{12}`\ C fusion. It is the convention of the
   reciprocity theorem,
   :math:`w_{12} k_{12}^2 \sigma_{12\to34}/(1+\delta_{12}) =
   w_{34} k_{34}^2 \sigma_{34\to12}/(1+\delta_{34})` with
   :math:`w = (2i_1+1)(2i_2+1)`, and counts each reaction once. The
   differential cross section of a reaction (d(d,p)t, say) carries the same
   factor 2, so that its integral over :math:`4\pi` is the angle-integrated
   value; it is symmetric about 90 degrees. (Before September 2026 only the
   elastic cross section had the factor; reaction and capture cross sections
   out of an identical pair were a factor 2 low.)

Reactions into an identical pair
   An identical pair in the **exit** channel (:math:`^7`\ Li(p,\ :math:`\alpha`)\ :math:`\alpha`)
   takes no factor: the cross section counts reactions, not outgoing
   particles. Data that count both alphas of every event -- an angle-integrated
   yield summed over both particles, or a detector at :math:`\theta` that
   sees either alpha -- are twice the model and must be halved (or given a
   normalization of 2).

Reaction rates
   The rate is :math:`N_A\langle\sigma v\rangle` of this physical cross
   section, so for an identical entrance pair it includes
   :math:`1+\delta_{12}` through :math:`\sigma` and, like the rate
   compilations (Fowler, Caughlan and Zimmerman 1967; NACRE; REACLIB),
   **excludes** the pair-counting :math:`1/(1+\delta_{12})` of the rate
   equation :math:`r = n_1 n_2 \langle\sigma v\rangle / (1+\delta_{12})`.
   AZURE2 prints a reminder when it computes such a rate.

Trojan Horse (HOES) cross sections
   These are in arbitrary units, and every channel of a THM segment shares
   its entrance pair, so a :math:`1+\delta_{12}` would only rescale the fitted
   normalization; none is applied.

The vector analyzing power is defined for a spin-1/2 beam only (it is zero
for d+d).

Several Compound Nuclei in One Project
--------------------------------------

A natural or compound target gives an abundance-weighted sum of cross
sections of *different* compound nuclei: natural neon under a neutron beam is
0.9048 n+\ :sup:`20`\ Ne (:sup:`21`\ Ne), 0.0925 n+\ :sup:`22`\ Ne
(:sup:`23`\ Ne) and 0.0027 n+\ :sup:`21`\ Ne (:sup:`22`\ Ne). One project
can hold the level schemes of all of them.

The compound nucleus of a pair is inferred from its columns -- for a
(Particle, Particle) pair :math:`Z = Z_1 + Z_2` and :math:`A` the rounded sum
of the masses, for a (Particle, Gamma) pair the heavy particle itself, for a
(Beta Decay) pair the daughter -- so nothing has to be entered. Pairs with the
same :math:`(Z, A)` belong to one nucleus, and a project whose pairs all form
one nucleus behaves exactly as before. As soon as the pairs form more than one,
the tab shows

- a read-only **Compound Nucleus** column, and the pair drop-downs of the
  other tabs append the nucleus to each pair (``n + 22Ne -> 23Ne``);
- the **Isotopic abundances** table: one row per compound nucleus, with the
  abundance of the *target isotope* that forms it (the row is labelled with
  both, ``23Ne (22Ne+n)``), a **Fit?** box that makes the abundance a fit
  parameter named ``abundance_<Z>_<A>``, and a **Prior Sigma** that adds
  :math:`((w - w_0)/\sigma)^2` to the chi-squared when it is not zero. The
  **Default sum-to-one sigma** below the table adds :math:`((\sum w - 1)/s)^2`
  per sample; zero switches it off.
- the **Sample** row above the table. Data sets taken on targets of
  different composition (natural neon and an enriched ²²Ne target, say, or a
  compound target) each need their own abundance set: a *sample*. The table
  shows one sample at a time; **+ sample** adds a named one (one word),
  **- sample** removes it, and the sample's **own sum-to-one sigma** overrides
  the default (``0`` switches the term off for a stoichiometric sample whose
  weights do not sum to one, such as 2 and 3 for Al₂O₃). A named sample's
  parameters are ``abundance_<name>_<Z>_<A>``. Every segment picks its sample
  in the segment dialog (:doc:`segments`); the default sample is the one a
  plain *natural* segment uses.

Which segments are weighted by these abundances is decided segment by
segment (the **Sample** box, :doc:`segments`); the levels of each
nucleus are entered on the :doc:`levels_channels` tab, and the abundances
appear on the :doc:`fitting` and :doc:`mcmc` tabs. In the file the table is
the ``<abundances>`` block and the sum sigma the ``abundanceSumSigma`` key
(:doc:`../reference/data_formats`), named samples as ``<name>:``-prefixed
lines; a project without the block has no abundance parameters at all. The pyazr function
:func:`pyazr.merge_projects` builds such a project from several
single-nucleus ones.

.. note::

   The identity in the file is the **compound** nucleus (:sup:`23`\ Ne for
   n+\ :sup:`22`\ Ne), while the abundance is that of the target isotope
   (:sup:`22`\ Ne). Every label shows both; the parameter name uses the
   compound nucleus, ``abundance_10_23``.
