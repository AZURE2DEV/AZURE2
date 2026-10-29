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

What is not symmetrized
   Reactions into or out of an identical pair (d+d :math:`\to` p+t, say) use
   only the channel rule above; no exchange factor is applied to their cross
   sections, the same convention as for spin-0 pairs. The vector analyzing
   power is defined for a spin-1/2 beam only (it is zero for d+d).

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

Masses and separation energies
------------------------------

A pair carries both its two masses and its separation energy, and the two are
linked: for two particle pairs of the same compound nucleus,
:math:`S_b - S_a = (m_{1b} + m_{2b} - m_{1a} - m_{2a})\,c^2`. AZURE2 takes the
thresholds (and so the whole energy scale) from the separation energies and
uses the masses only in the kinematics (reduced masses, lab-to-c.m.
conversions), so a mismatch does not change a calculation by more than the
mass error itself (1e-4 u on a mass of 12 u is 1e-5). It is, however, how a
typo or a mixed mass convention shows up, and since October 2026 the engine
reports one when it reads the pairs:

.. code-block:: text

   **WARNING: separation energies and masses disagree: pair 2 (S = 18.8943 MeV,
   m1 + m2 = 8.025593633 u) and pair 1 (S = 17.2551 MeV, m1 + m2 = 8.023828466 u):
   the separation energies differ by 1.6392 MeV, the masses by 1.6442 MeV (-5.0 keV). ...

Each particle pair is compared with the first particle pair of its compound
nucleus; mismatches above 1 keV are listed (the first three, then a count). A
mismatch within 5 keV of a multiple of the electron mass adds a hint that atomic
and nuclear masses are mixed, typically a bare-proton mass (1.00728 u) next to
atomic masses for the rest. Gamma and beta pairs are not compared: their heavy
mass is only a label. A project that uses mass numbers as masses (12, 4, ...)
gets one note instead, since its masses carry no Q value. The cure is to take
the masses from one mass table (AME, atomic masses throughout) and the
separation energies from the same table: :math:`S = m_1 + m_2 - M` of the
compound nucleus.

Relativistic Kinematics
-----------------------

The **Relativistic kinematics** check box below the pairs table (``relativistic=1``
in the ``<config>`` block, so command-line and pyazr runs of the project use it
too) takes every channel's kinematics from the invariant mass :math:`s` instead
of the non-relativistic relations.  It follows the prescription of the LANL code
EDA (Haeck, Hale, Paris and Thompson, LA-UR-19-31084, the ENDF ``KRL=1`` flag).
For a pair of masses :math:`m_1, m_2` (:math:`m = m_1 + m_2`) at centre-of-mass
kinetic energy :math:`T = \sqrt{s} - m`:

.. math::

   s = m^2 + 2 m_2 E_\text{lab}, \qquad
   k^2 = \frac{(s - m^2)\,(s - (m_1 - m_2)^2)}{4 s}, \qquad
   \eta = \frac{Z_1 Z_2 \alpha\,(s - m_1^2 - m_2^2)}{2 k \sqrt{s}},

so the penetrabilities, shift functions, hard-sphere and Coulomb phases, the
:math:`\pi/k^2` factor and the bound-state decay constants (:math:`k^2 < 0`) all
use these, and a data point's :math:`T` comes from its lab energy through
:math:`s` (and, for a beta-delayed particle, from the decay at rest).  Particle
lab angles and differential cross sections are converted with the exact Lorentz
boost (:math:`\gamma`-ray data already were). Where two c.m. angles reach one
lab angle (just above the threshold of an endothermic exit, or a projectile
heavier than the target) the forward root is taken, and the lab angle has a
maximum, :math:`\sin\theta_\text{max} = \gamma_3^*/(\gamma\, g)` with
:math:`g = \beta/\beta_3^*`. That maximum is smaller than the
non-relativistic one, by 0.3 to 0.6° just above the 7Li(p,n) and 12C(α,p)
thresholds, and a point beyond it, or below the relativistic threshold, is
refused (before October 2026 the non-relativistic limits decided, so a point
between the two maxima was converted with the non-relativistic angle and
Jacobian without a warning). The conversions agree with the independent
reaction-kinematics calculator (det-lab) to 1e-7° in angle and 4e-9 in
:math:`d\Omega_\text{lab}/d\Omega_\text{cm}` for 7Li(p,n), t(d,n)α,
12C(α,p), p + 4He and 12C + p (inverse kinematics) up to 200 MeV.  At low energy everything reduces
to the familiar :math:`k^2 = 2\mu T`; threshold quantities (scattering lengths,
the Wigner limit, ANCs) keep the reduced mass.

What changes, and by how much.  Per channel the corrections are of order
:math:`T/m`: at :math:`E_\text{lab} = 20` MeV, :math:`k` of n+\ :sup:`4`\ He is
0.19 % larger, :math:`\eta` of p+\ :sup:`4`\ He 1.6 % larger, and c.m. angles and
cross sections move by up to 0.2° and 0.7 %.  The larger effect is between
partitions: the excitation energy of a point is :math:`T + S_c`, and with
:math:`T` from :math:`s` the energy scales of two pairs with a large Q value no
longer differ by a constant.  In :sup:`5`\ He (Q = 17.6 MeV) an n+\ :sup:`4`\ He
point at :math:`E_n = 20` MeV sits 27 keV away from where non-relativistic
kinematics puts it, relative to d+t.

Level energies.  The R-matrix denominators stay :math:`E_\lambda - E_x`: level
energies remain excitation energies, comparable with tabulations.  EDA instead
uses :math:`E(s) = (s - m_0^2)/(2 m_0)` of a reference pair; on the excitation
scale :math:`E_s = S_0 + T_0 + T_0^2/(2 m_0)`, :math:`T_0 = E_x - S_0`, which is
43 keV above :math:`E_x` at :math:`T_0 = 20` MeV in A = 5.  ``relativistic.out``
(:doc:`../reference/output_files`) lists every level in that convention, with
the widths scaled by :math:`dE_s/dE_x = 1 + T_0/m_0`, for comparison with EDA
parameters or ENDF ``KRL=1`` evaluations; ``relativisticReferencePair=<key>``
chooses the reference pair (default: the first particle pair).

Limits.  The hybrid nuclear-potential Coulomb functions (a ``<potential>``
block) solve a Schrödinger equation inside the channel radius and are refused
with ``relativistic=1``; they bend the scattering functions only, not the
bound-state tail, ANC conversion or NF integral of external capture
(:doc:`pyazr`, "What the potential does not touch").  Target effects (beam spread, target thickness, sub-point
grids) and energy shifts convert lab energy intervals to the c.m. frame with the
constant ratio :math:`m_2/m`; the exact ratio :math:`m_2/\sqrt{s}` differs by
:math:`T/m` of the interval (0.4 % of a few-keV shift at 20 MeV).  Photon
energies in capture channels do not include the recoil term.
