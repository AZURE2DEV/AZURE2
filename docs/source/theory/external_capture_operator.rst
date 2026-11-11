The External-Capture Operator Beyond the Long-Wavelength Approximation
=======================================================================

External capture (EC) is the part of a radiative-capture amplitude that comes
from beyond the channel radius :math:`a`: the hard-sphere (or resonant-external)
scattering wave of the entrance pair, the Whittaker tail of the final bound
state with its asymptotic normalization coefficient (ANC), and an
electromagnetic operator between them.  By default AZURE2 uses the
long-wavelength approximation (LWA) for that operator.  ``--no-long-wavelength``
(pyazr ``use_long_wavelength=False``) replaces it by the exact operator of two
point clusters.  This chapter derives that operator and describes the checks
on it.  Internal (resonant) photon channels are not affected; their widths are
fit parameters.

Long-wavelength operator
------------------------

With the two clusters at :math:`c_1\mathbf r` and :math:`c_2\mathbf r`
(:math:`c_1 = m_2/M`, :math:`c_2 = -m_1/M`, :math:`\mathbf r` the relative
coordinate), the electric multipole operator in the LWA is

.. math::

   \mathcal O_{LM} = e\,\bar e_L\, r^L Y_{LM}(\hat r), \qquad
   \bar e_L = Z_1\left(\frac{m_2}{M}\right)^L + Z_2\left(-\frac{m_1}{M}\right)^L ,

the effective charge of Azuma et al., PRC 81, 045805 (2010), Eq. (24).  For
pairs with :math:`Z_1/m_1 \approx Z_2/m_2` (:math:`\alpha + {}^{12}\mathrm C`,
:math:`\alpha + {}^{10}\mathrm B`, ...) the E1 charge nearly vanishes
(:math:`-0.00098` for :math:`\alpha + {}^{12}\mathrm C`, from the mass defects
alone), and the first term that does not vanish is of order :math:`k_\gamma^2`.

The transverse electric field and its two parts
-----------------------------------------------

The photon's electric multipole field is
:math:`\mathbf V = k^{-1}\nabla\times[j_L(kr)\,\mathbf X_{LM}]`, with
:math:`k = E_\gamma/\hbar c`.  Split it into a gradient and a remainder
tangential to :math:`\mathbf r`:

.. math::

   \mathbf V = \nabla\Phi + \mathbf W, \qquad \mathbf W\cdot\hat r = 0,

   \Phi = \frac{i\sqrt{L(L+1)}}{k}\,Y_{LM}\,G_L(kr), \qquad
   G_L(x) = \int_0^x \frac{j_L(t)}{t}\,dt ,

   \mathbf W = i\,w_t(kr)\;r\nabla Y_{LM}, \qquad
   w_t(x) = \frac{(L+1)\,j_{L-1}(x) - L\,j_{L+1}(x)}{(2L+1)\sqrt{L(L+1)}}
            - \sqrt{L(L+1)}\,\frac{G_L(x)}{x} .

:math:`\Phi` is fixed by the radial component,
:math:`\hat r\cdot\mathbf V = i\sqrt{L(L+1)}\,j_L(kr)/(kr)\,Y_{LM}`.
:math:`w_t` is of order :math:`x^{L+1}`, so it vanishes in the LWA.

**Charge term.**  With current conservation, the coupling of the current to
:math:`\nabla\Phi` becomes :math:`\omega\int\rho\,\Phi`.  The radial operator
is then

.. math::

   \frac{L\,(2L+1)!!}{k^L}\sum_i Z_i s_i\, G_L(k|c_i| r), \qquad s_1 = 1,\ s_2 = (-1)^L .

This is the charge term of the extended Siegert theorem (J. L. Friar and
S. Fallieros, PRC 29, 1645 (1984)).  As :math:`k\to 0`,
:math:`G_L(x)\to x^L/(L(2L+1)!!)`, and the term reduces to
:math:`\bar e_L r^L` term by term.

**Convection-current term.**  The convection current of the relative motion,
:math:`\mathbf J = \sum_i e Z_i\,\mathbf p_i/m_i` (:math:`\mathbf p_1 = -\mathbf p_2 = \mathbf p`),
couples to :math:`\mathbf W`.  Take states :math:`u_i(r)/r\,Y_{l_i}` and
:math:`u_f(r)/r\,Y_{l_f}`.  Because :math:`\mathbf W` is tangential,
:math:`\nabla\cdot(f\,r\nabla Y_{LM}) = -L(L+1)fY_{LM}/r`, and
:math:`\int Y_{l_f}^*\,\nabla_\Omega Y_{LM}\cdot\nabla_\Omega Y_{l_i}` is the Gaunt
coefficient times :math:`\tfrac12[L(L+1)+l_i(l_i+1)-l_f(l_f+1)]`.  So

.. math::

   \langle f|\{\mathbf p, w\,r\nabla Y_{LM}\}|i\rangle
   = i\hbar\,[\,l_f(l_f+1)-l_i(l_i+1)\,]\,\langle Y_{l_f}|Y_{LM}|Y_{l_i}\rangle
     \int u_f\,u_i\,\frac{w}{r}\,dr .

This has the same angular factor as the charge term and the same two radial
functions, so it adds to the radial operator:

.. math::

   \frac{L\,(2L+1)!!}{k^L}\;
   \frac{-[\,l_f(l_f+1)-l_i(l_i+1)\,]}{\sqrt{L(L+1)}}
   \sum_i Z_i s_i\,\frac{\hbar c}{2 m_i c^2}\,\frac{w_t(k|c_i| r)}{r} .

For an N = Z pair, this term and the :math:`k^2 r^3` part of the charge term
are of the same order: the leading isoscalar E1 operator has both (Baye,
PRC 86, 034306 (2012)).  Spin-magnetization currents, which pairs with spin
would add, are not included.

What AZURE2 does
----------------

``--no-long-wavelength`` uses the charge term plus the convection-current term
for E1 and E2 external capture.  ``--no-ec-convection-current`` (pyazr
``ec_convection_current=False``) keeps the charge term alone, which was the
only exact form before 2026-10-23.  M1 external capture uses its
long-wavelength operator in either mode.  The ``intEC.dat`` signature records
both switches, so a saved file is recomputed when either one changes.

Size, for :math:`{}^{12}\mathrm C(\alpha,\gamma){}^{16}\mathrm O`
with :math:`a = 5.43` fm.  The E1 amplitude to the ground state is 3.5–4.7
times the LWA one with the charge term, and the convection current adds
another +10 %.  For the E1 amplitude to the 6.92 MeV :math:`2^+` state the
convection current contributes −5 % to −20 % (:math:`l = 1`) and +7 % to +28 %
(:math:`l = 3`) between 0.3 and 3 MeV.  On the archive's
:math:`{}^{12}\mathrm C + \alpha` model, at fixed parameters, the convection
current changes :math:`S_{\gamma_0}` by 0.1–0.4 % and the 6.13 MeV cascade by
up to 1 %.

Checks
------

* The current form of the gradient part equals the charge form when both
  states are exact eigenstates of one Hamiltonian.  This was checked for
  :math:`\alpha + {}^{12}\mathrm C` in a Woods–Saxon plus Coulomb potential, to
  1e-9 for each cluster.  It fixes every sign and factor of the derivation.
* ``tests/ec_m1_exact`` compares AZURE2's absolute external-capture cross
  sections with an independent golden-rule calculation (``truth.py``) for
  p + :sup:`16`\ O :math:`\to` :sup:`17`\ F(0.495) (E1 in all three forms,
  and M1) and for :math:`\alpha` + :sup:`12`\ C :math:`\to` :sup:`16`\ O(g.s.)
  (E1 + E2 in all three forms).  Since 2026-10-25 it also covers
  p + :sup:`14`\ N :math:`\to` :sup:`15`\ O by s-wave M1, with a moment on
  both clusters and a target of spin 1.

The M1 spin term :math:`g_1 \mathbf S_1 + g_2 \mathbf S_2` (the pair
columns g1 and g2) had a sign error before 2026-10-25.  A cluster's moment
term had the wrong sign whenever the other cluster had half-integer spin.  For
pairs such as p + :sup:`14`\ N the M1 external capture was wrong in size
(about a factor 0.5–3 depending on the channel spins).  For two clusters with
half-integer spin (p + :sup:`13`\ C, p + :sup:`15`\ N) the spin term had the
wrong sign against the orbital term.  Models with a spin-0 target, or with
only one non-zero g-factor and the other cluster of integer spin, were not
affected.
