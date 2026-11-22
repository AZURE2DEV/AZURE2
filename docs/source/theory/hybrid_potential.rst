.. _theory-hybrid-potential:

The Hybrid Potential and its Spin–Orbit Term
============================================

A *hybrid* particle pair takes its channel functions not from the Coulomb functions :math:`F_l, G_l` but from the
solutions of a radial equation with a nuclear potential besides the Coulomb one.  This chapter states the equation,
the spin–orbit term and its convention, and what depends on the total angular momentum :math:`J` once the term is
present.  The design note is ``docs/upgrade/spin_orbit_potential.md``; the input is the ``potential`` record of
format 3 (:doc:`../reference/format3`), edited in pyazr as a record of the ``Project`` (:doc:`../user_guide/pyazr`).

The radial equation
-------------------

For a pair with reduced mass :math:`\mu` (MeV) and charge numbers :math:`Z_1, Z_2`, at the c.m. channel energy
:math:`E` (MeV), orbital angular momentum :math:`l` and total angular momentum :math:`J`,

.. math::

   u''(r) = \frac{2\mu}{\hbar^2}\,\bigl[V(r) - E\bigr]\,u(r),\qquad
   V(r) = V_C(r) + V_N(r) + \frac{\hbar^2 l(l+1)}{2\mu r^2} + v_{so}(r)\,\langle 2\,\mathbf{l}\cdot\mathbf{i}\rangle,

with :math:`V_C` the potential of a uniformly charged sphere of radius :math:`R_C` (``coulomb_radius``,
:math:`e^2 = \alpha\hbar c`), :math:`V_N` a Woods–Saxon :math:`-V_0/(1 + e^{(r-R)/a})` or a Gaussian
:math:`-V_0 e^{-(r/r_0)^2}` (:math:`V_0 > 0` attractive), and the spin–orbit term below.  The equation is integrated
inward by Numerov's method from the matching radius :math:`r_m` (15 fm), where the solution is the point-Coulomb
:math:`H^+_l = G_l + iF_l` (above threshold) or the Whittaker function :math:`W_{-\eta_b,l+1/2}(2\kappa r)` (below
it), to the channel radius :math:`a_c`.  The potential therefore acts between :math:`a_c` and :math:`r_m`; the
interior is the R-matrix's.  At :math:`a_c` the solution :math:`H` gives
:math:`S + iP = \rho H'/H` and the hard-sphere phase :math:`\phi = \arg H(a_c)`; below threshold the decaying
solution :math:`\tilde W` gives :math:`S = a\tilde W'/\tilde W` and the ANC-to-width factor.

The spin–orbit term
-------------------

.. math::

   v_{so}(r) = V_{so}\left(\frac{\hbar}{m_\pi c}\right)^2 \frac{1}{r}\frac{df}{dr},\qquad
   \left(\frac{\hbar}{m_\pi c}\right)^2 = 2.000\ \mathrm{fm}^2,\qquad
   f(r) = \frac{1}{1 + e^{(r - R_{so})/a_{so}}},

.. math::

   \langle 2\,\mathbf{l}\cdot\mathbf{i}\rangle = j(j+1) - l(l+1) - i(i+1),

with :math:`i` the spin of the particle the term acts on (``on light`` or ``on heavy``) and :math:`j` the coupled
angular momentum of :math:`\mathbf{l} + \mathbf{i}`.  This is FRESCO's convention (the projectile and target
spin–orbit form factors of its ``&POT`` namelist, with ``CONLS`` = 2.000 and the factor
:math:`j(j+1) - l(l+1) - s(s+1)`).  Since :math:`df/dr < 0`, :math:`V_{so} > 0` lowers :math:`j = l + i`.  For
:math:`i = 1/2`, :math:`\langle 2\,\mathbf{l}\cdot\mathbf{i}\rangle = \langle \mathbf{l}\cdot\boldsymbol{\sigma}\rangle`
is :math:`l` for :math:`j = l + 1/2` and :math:`-(l+1)` for :math:`j = l - 1/2`.

**The factor 2.**  The operator is :math:`2\,\mathbf{l}\cdot\mathbf{i}`, not :math:`\mathbf{l}\cdot\mathbf{i}`.  A
strength quoted with an :math:`\mathbf{l}\cdot\mathbf{s}` operator and no 2 corresponds to twice the :math:`V_{so}` of
this convention for the same physics; the record therefore requires the words ``operator 2l.s``, so that a strength
copied from elsewhere is checked against its convention.  Whether a published global set uses
:math:`(\hbar/m_\pi c)^2 = 2.000` fm\ :sup:`2` or :math:`(\hbar c/m_{\pi^\pm}c^2)^2 = 1.9989` fm\ :sup:`2` is a
difference of 5·10\ :sup:`-4` that the user must check against the source.

**One spinning particle.**  AZURE2 labels channels by the channel spin :math:`\mathbf{s} = \mathbf{i}_1 + \mathbf{i}_2`
with :math:`\mathbf{J} = \mathbf{s} + \mathbf{l}` (:doc:`conventions`).  When only one particle has spin,
:math:`s = i`, :math:`j = J`, and the term is diagonal in the channels: each :math:`(\alpha, l, J)` has its own scalar
equation.  When both particles have spin, a term on one of them couples the channel spins of one :math:`(l, J)`; that
coupled treatment is not implemented, and such a pair is refused (``E_SO_TWO_SPINS``), as is a term on a particle of
spin 0 (``E_SO_ZERO_SPIN``), a term with relativistic kinematics (``E_SO_RELATIVISTIC``) and on the entrance pair of
THM data (``E_SO_THM``).  The potential's numbers are not fit parameters.

What depends on J
-----------------

With a spin–orbit term every channel function of the pair depends on :math:`(l, J)`:

- :math:`P`, :math:`S` and :math:`e^{-i\phi}` of an open channel at a data point (the channel table holds a row per
  :math:`(\alpha, l, 2J, E_c)`; a pair without the term has one row per :math:`(\alpha, l, E_c)`, its :math:`2J` key
  being a placeholder);
- :math:`S`, :math:`S'`, :math:`S''` at the level energies (Brune and Park transformations, boundary values) and the
  closed channels;
- the ANC-to-reduced-width conversion of a bound state, through :math:`\tilde W_{lJ}(a_c)`;
- external capture: the entrance wave of :math:`(l_i, J_i)` and the bound state of :math:`(l_f, J_f)`.

The Coulomb phase :math:`e^{i(\sigma_l - \sigma_0)}` is the point-Coulomb one in every case.  In pyazr,
``channel_functions(pair, E, l, J=...)`` needs :math:`J` for a pair with the term (``E_SO_NEEDS_J``) and refuses it
for any other (``E_SO_NO_TERM``).

The interior of the potential as a channel's background
-------------------------------------------------------

.. note::

   A research option, off unless a ``potential`` record has a ``wave`` field with ``background interior``.  The
   design note, with the derivations, the tests and what is not done, is ``docs/upgrade/potential_background.md``.

The hybrid pair above uses the potential only outside the channel radius :math:`a_c`.  The hybrid R-matrix models of
Koonin, Tombrello and Fox and of Langanke and Koonin use its interior too: the potential describes the non-resonant
part of a partial wave, and levels are added only for what it does not produce.  With :math:`u_E` the regular
solution of the radial equation and

.. math::

   \ell(E) = \frac{a_c\,u_E'(a_c)}{u_E(a_c)}, \qquad R_\text{pot}(E) = \frac{1}{\ell(E) - B_\text{bg}},

:math:`R_\text{pot}` is the R-function of the potential in :math:`[0, a_c]` for the boundary constant
:math:`B_\text{bg}`: a sum of poles :math:`\gamma_p^2/(E_p - E)` at the eigenvalues of the potential in that box.
Added to the diagonal R-matrix element of the channel, :math:`R_{cc} = \sum_\lambda \gamma_{\lambda c}^2/(E_\lambda -
E) + R^0_c(E)`, it changes nothing in the level matrix or the collision matrix except the channel's functions at
:math:`a_c`, which become those of the outgoing wave :math:`H = G + iF` mapped by

.. math::

   \tilde H = H - R^0(\rho H' - B_\text{bg} H), \qquad \rho\tilde H' = \rho H' - B_\text{bg} R^0(\rho H' - B_\text{bg} H),

so that :math:`\tilde P = \rho/|\tilde H|^2`, :math:`\tilde\phi = \arg\tilde H` and :math:`\tilde S =
\operatorname{Re}\rho\tilde H'/\tilde H`.  With no level in the channel the collision matrix is then that of the
potential itself, whatever :math:`a_c` is.  Levels on top see :math:`\tilde S, \tilde P`: their observed (Brune)
energies and widths are those of a resonance on top of the background.

The record states it per partial wave::

    wave l <l> [J <J>] [depth <V> MeV] [width <r0> fm] background interior boundary <B> [subtract <n>]
               [capture interior [charge_e1 <q>] [charge_e2 <q>]]

*A state counted twice.*  A level that is also a state of the potential in that wave is in :math:`R_\text{pot}`
already.  ``subtract <n>`` removes the box state with :math:`n` nodes from the background,
:math:`R^0 = R_\text{pot} - \gamma_p^2/(E_p - E)`, and the level stands for it; or the depth is chosen so that the
potential has no state there.  The set-up warns of a level within one width of a resonance or bound state of the
potential in a wave with a background, and notes it when the wave has a ``subtract``.

*Capture* (``capture interior``).  External capture runs from :math:`a_c` outward; the background's part of the
entrance wave inside :math:`a_c` is :math:`D\,g_E(r)`, with :math:`D = \rho u' - B_\text{bg}u` the surface term of
the channel's wave at :math:`a_c` and :math:`g_E(r) = R_\text{pot}(E)\,u_E(r)/u_E(a_c)` less the subtracted box
states.  Its overlap with the final state inside :math:`a_c`,

.. math::

   K(E) = \int_0^{a_c} \varphi_f(r)\,w(r)\,g_E(r)\,dr, \qquad \varphi_f = u_{E_f}/u_{E_f}(a_c),

is added to both terms of the external capture (:math:`w` the multipole's weight).  A bound final state must be a
bound state of the potential in its partial wave within 1 keV of the level's energy (the run is refused otherwise,
with the depth that binds there); its size is the level's ANC.  With no levels the capture is the potential model's
direct capture and does not depend on :math:`a_c`.

- *The interior charge.*  The :math:`E_L` operator of the relative motion carries :math:`Z_1 (m_2/M)^L + Z_2
  (-m_1/M)^L`, which nearly vanishes for E1 between nuclei of equal :math:`Z/A` (:math:`-10^{-3}` for α +
  :sup:`12`\ C); the strength inside :math:`a_c` is then the interior's own.  ``charge_e1 q`` (``charge_e2 q``) adds
  :math:`e\,q\,r^L` to the operator of the interior capture of that wave, inside :math:`a_c` only.
- *An unbound final state* is a narrow resonance whose wave outside :math:`a_c` is the outgoing wave
  :math:`O_f(r)/O_f(a_c)`; inside it is continued in value by :math:`\varphi_f(r)` at :math:`E_f`.  There is no rule
  that the potential has a resonance there; the set-up prints the mismatch of the logarithmic derivatives.  For
  capture into a resonance of the potential itself, given as a level of the resonance's width, the cross section is
  the same at :math:`a_c` = 4.0, 5.35 and 6.5 fm to 0.4 % (α + :sup:`12`\ C, a state 5 % of its energy wide), against
  a factor of up to 4.6 without the interior.  The line shape of a broad final state (``line_window``) is refused
  with interior capture.

*The rate.*  The rate's energy grid has a cluster of points on each resonance of the entrance pair's potential in a
wave with a background, as it has on each level.

The numbers of the record (depths, :math:`B_\text{bg}`, the interior charge) are inputs, not fit parameters of the
engine; ``pyazr.record_fit`` fits them from outside (:doc:`../user_guide/pyazr`).

Checks
------

The term is checked against answers independent of the engine (``physics_spin_orbit``,
``tests/spin_orbit_potential``, ``tests/pyazr/spin_orbit_potential_test.py``): the angular factor against an
m-scheme matrix element built from Clebsch–Gordan coefficients (3.6·10\ :sup:`-15`); a square well with the Thomas
term of its sharp edge, a δ shell solved in closed form with Riccati–Bessel functions, against the solver on
Woods–Saxon edges of vanishing diffuseness (7.7·10\ :sup:`-4`); :math:`S, P, \phi` above and :math:`S`,
:math:`\ln\tilde W(a_c)` below threshold for p + :sup:`12`\ C against scipy's DOP853 with mpmath's Coulomb and
Whittaker functions (8.4·10\ :sup:`-8`); the phase shifts of a project against the same integration
(1.5·10\ :sup:`-9` rad); unitarity with the term; the external-capture integrals per :math:`J` against an
independent integration (6.5·10\ :sup:`-10` in the ratio with/without the term).  A pair without the term gives the
bits it gave before the term existed.
