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
