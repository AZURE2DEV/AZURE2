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
:math:`-V_0 e^{-(r/r_0)^2}` (:math:`V_0 > 0` attractive, :math:`V_0 < 0` repulsive), and the spin–orbit term below.  The equation is integrated
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

A potential that also acts inside the channel radius
----------------------------------------------------

With ``interior_boundary <b>`` on its record the potential acts on the whole radial range.  The R-matrix of the
channel is then the potential's own R-function plus the levels,

.. math::

   R(E) = R_{pot}(E) + \sum_\lambda \frac{\gamma_\lambda^2}{E_\lambda - E},\qquad
   R_{pot}(E) = \frac{1}{f(E) - b},\qquad f(E) = \left.\frac{a_c\,u'(a_c)}{u(a_c)}\right|_{u\ \mathrm{regular\ at}\ 0},

where :math:`u` is the solution of the radial equation above that is regular at the origin (integrated outward by
Numerov's method) and :math:`b` the boundary constant of the record.  A potential stands in this way for the broad
single-particle states and the smooth background that otherwise need background poles (Koonin, Tombrello and Fox,
Nucl. Phys. A 220 (1974) 221; Langanke and Koonin, Nucl. Phys. A 410 (1983) 334).  Folding :math:`R_{pot}` into the
channel functions leaves the collision matrix in its usual form with

.. math::

   S' + iP' = b + \frac{(L - b)(f - b)}{f - L},\qquad
   e^{-i\phi'} = e^{-i\phi}\,\frac{(f - L)^*}{|f - L|}\,\mathrm{sgn}(f - b),\qquad L = S + iP,

so every part of the calculation that takes :math:`P`, :math:`S` and :math:`\phi` (the level matrix, Brune's and
Park's parameters, the conversions of observed widths, external capture's entrance waves) takes the primed ones.
Without levels the channel's phase is the phase shift of the potential, at any channel radius.  The energy derivative
of :math:`f` is the interior norm, :math:`df/dE = -(2\mu a_c/\hbar^2)\int_0^{a_c} u^2\,dr / u(a_c)^2`.

Below threshold :math:`P = 0` and :math:`S' = b + (S - b)(f - b)/(f - S)`.  A bound level's amplitude at the surface
is :math:`\gamma\,(f - b)/(f - S)`, so that

.. math::

   C^2 = \frac{2\mu a_c}{\hbar^2 W(a_c)^2}\left(\frac{f - b}{f - S}\right)^2
         \frac{\gamma^2}{1 + \sum_c \gamma_c^2\,dS'_c/dE},

and :math:`(\hbar^2/2\mu a_c)\,dS'/dE` is the exterior norm integral :math:`\int_{a_c}^\infty W^2/W(a_c)^2\,dr` plus
the interior one of the potential's wave, :math:`((S - b)/(f - b))^2\int_0^{a_c}u^2\,dr/u(a_c)^2`, both times
:math:`((f - b)/(f - S))^2`.  An external-capture final state takes the same two factors.  In the limit of a large
reduced width at a bound state of the potential (:math:`f = S`) the ANC is the potential's single-particle ANC.  The
sign of an ANC stays the sign of its reduced width amplitude; where :math:`(f - b)/(f - S) < 0` the amplitude at the
surface has the opposite sign, which matters only for a level with two closed particle channels.

A level placed on a state the potential already produces describes it twice.  At set-up the bound states and
resonances of each such potential (:math:`f = S`; above threshold only zeros narrower than half their energy) are
found, and a level with a width in the channel that lies within one width of one, the larger of the level's observed
particle width and the state's, draws a warning; the calculation proceeds.

**Capture from the potential's interior.**  The interior wave function is a sum over the levels and over the
states of the potential, and the part the potential's states carry is closed:
:math:`\chi_c(r) = [u_c(r)/u_c(a_c)]\,R_{pot,c}\,(a_c\psi_c' - b_c\psi_c)` for a channel wave :math:`\psi_c` at the
surface.  The two terms of the external capture therefore continue inside the radius, when the entrance and the
final channel both have ``interior_boundary`` and the final state is bound:

* the non-resonant term's standing wave :math:`v(r)` continues as :math:`v(a_c)\,u(r)/u(a_c)`: without levels the
  term is the direct capture of the potential model, over the whole radial range;
* a level's channel wave :math:`O_c(r)/O_c(a_c)` continues as :math:`[(L_c - b_c)/(f_c - b_c)]\,u(r)/u(a_c)`
  (:math:`L = S` below threshold); the remainder of the level's surface amplitude is its own interior, the level's
  internal photon width;
* the bound final state :math:`C\,W(r)` continues as
  :math:`C\,W(a_c)\,[(S_f - b_f)/(f_f - b_f)]\,u_f(r)/u_f(a_c)`, the term already in its norm above.  At a bound
  state of the potential the factor is 1 and the final state is the potential's wave function.

The operator inside is the one outside.  ``interior_charge_e1 <q>`` (and ``interior_charge_e2``) on the record of
the entrance channel's :math:`l` adds :math:`q\,r^L` (:math:`q` in units of :math:`e`) inside the radius: for an
:math:`N = Z` pair the E1 operator of the relative motion nearly vanishes (:sup:`12`\ C + α:
:math:`Z_1 m_2/M - Z_2 m_1/M \approx -10^{-3}`) and the strength is the interior's (Koonin, Tombrello and Fox fit an
effective dipole strength).  Left out: the matrix element between the *level* part of the final state and the
potential part of the entrance wave; it vanishes when the final state is a state of the potential.  A final state
whose channel has no ``interior_boundary`` of its own :math:`l`, or that is unbound, gets no interior capture and a
warning at set-up.  These terms are integrated at every energy rather than interpolated from the energy tables of
the external capture: they carry the potential's resonances, and at a pole of :math:`R_{pot}` the tables' divisors
vanish or change sign.

Limits: a pair whose two particles both have spin is refused (``E_POT_INTERIOR_TWO_SPINS``; :math:`R_{pot}` would
couple the channel spins); no capture from the interior into an unbound final state; the potential's parameters and
the interior charge are not fit parameters.  A potential may be given per orbital angular
momentum with the record's ``l`` list.  The design note is ``docs/upgrade/interior_potential.md``.

Checks
------

The interior R-function is checked by ``physics_interior_potential`` (:math:`f`, the interior norm and the phase of
:sup:`12`\ C + α, l = 0 to 2, against scipy's DOP853 with mpmath's Coulomb functions, 2·10\ :sup:`-6`) and
``upgrade_interior_bound`` (the bound-level factor; the norm identity to 1.3·10\ :sup:`-5`; the single-particle ANC at
three bound states of the potential to 2.3·10\ :sup:`-5`; the phase without levels against the regular solution
matched at 14 fm, 2.3·10\ :sup:`-8` rad; the warning).  The capture from the interior is checked by
``upgrade_interior_capture`` (each term with the interior over the same term from outside alone, at three channel
radii and two boundary constants, above and below threshold, against scipy and mpmath: 2.7·10\ :sup:`-5`, and
9·10\ :sup:`-7` with the step at 0.01 fm; M1 between two states of one potential cancels to 6·10\ :sup:`-7`; the
interior charge) and ``interior_capture`` (the whole engine: the potential model's direct capture cross section is
the same at channel radii 4.0, 5.35 and 6.5 fm to 2.3·10\ :sup:`-4`, the tolerance of the direct integrals).

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
