Scattering Lengths
==================

AZURE2 computes the s-wave scattering length of every particle pair at the pair's
threshold, writes it after each calculation or fit, and can fit to measured
values.  This chapter gives the definition, how the R-matrix parameters
determine it, and what the code does.

Why a scattering length
-----------------------

At low energy the s-wave phase shift of a pair is summarized by two numbers,
the scattering length :math:`a_0` and the effective range :math:`r_0`.  They
are the quantities that ab initio and effective-field-theory calculations
predict and that low-energy capture depends on: in
:sup:`3`\ He(:math:`\alpha,\gamma`)\ :sup:`7`\ Be and
:sup:`7`\ Be(p,\ :math:`\gamma`)\ :sup:`8`\ B the s-wave external capture
interferes with the internal part through the same background level that sets
:math:`a_0`, so the scattering length constrains the extrapolation of the
S factor to zero energy (Paneru et al. 2019, 2024).  For neutrons the
scattering lengths are measured directly, by neutron interferometry or
transmission at thermal energies, often to a fraction of a percent, and are
among the most precise constraints a neutron-induced R-matrix analysis can use.

Definition
----------

Below any resonance the phase shift of a charged pair is dominated by the
Coulomb barrier, and the expansion that is analytic in energy is the
Coulomb-modified effective-range function

.. math::

   K(E) = C_0^2(\eta)\, k \cot\delta_0 + 2 k \eta\, h(\eta)
        = -\frac{1}{a_0} + \frac{1}{2} r_0 k^2 + \ldots,

with :math:`\eta` the Sommerfeld parameter,
:math:`C_0^2 = 2\pi\eta/(e^{2\pi\eta}-1)` and
:math:`h(\eta) = \mathrm{Re}\,\psi(i\eta) - \ln\eta`.  For a neutral pair
:math:`C_0 = 1`, the :math:`h` term is absent and :math:`k\cot\delta_0 =
-1/a_0 + \ldots`.  The sign convention is the usual one: a hard sphere of radius
:math:`R` has :math:`a_0 = R > 0`.

The R-matrix expression
-----------------------

At zero energy the radial equation outside the channel radius :math:`a` is the
Coulomb equation with :math:`k = 0`, whose solutions are modified Bessel
functions of :math:`x = (8 r / a_B)^{1/2}`, :math:`a_B = \hbar^2/(\mu Z_1 Z_2
e^2)`.  Normalized as the :math:`k \to 0` limits of :math:`F_0/(C_0 k)` and
:math:`C_0 G_0`,

.. math::

   f(r) = \frac{a_B x}{4} I_1(x), \qquad g(r) = x K_1(x),

so :math:`f \approx r` near the origin, :math:`g(0) = 1`, and their Wronskian is
:math:`-1`.  The zero-energy wave function is :math:`u = g - f/a_0`.  The
R-matrix ties its logarithmic derivative at :math:`r = a` to the level matrix:
in the notation of Brune (2002) the elastic element
:math:`M_{cc} = \tilde\gamma_c^T \tilde A\, \tilde\gamma_c` satisfies
:math:`M_{cc} = 1/(a u'/u - S_c)`, with :math:`S_c` the shift function at zero
energy.  Eliminating :math:`u` gives (Paneru et al. 2019, Eq. 7; Paneru et al.
2024, Eq. 26)

.. math::

   a_0 = -a\left[\frac{M_{cc}}{x^2 K_1^2(x)} - \frac{2 I_1(x)}{x^2 K_1(x)}\right],
   \qquad x = \left(\frac{8 Z_1 Z_2 e^2 \mu a}{\hbar^2}\right)^{1/2},

and for a neutral pair (:math:`x \to 0`) the familiar :math:`a_0 = a\,(1 - M_{cc})`.
:math:`M_{cc}` does not depend on the boundary condition, so the result is the
same in the standard, Brune and Park parametrizations.

The level matrix is evaluated exactly at the pair's threshold.  There the
entrance channel has zero penetrability and the shift function takes its
zero-energy limit,

.. math::

   S_l(0) = -l - \frac{x}{2}\,\frac{K_{2l}(x)}{K_{2l+1}(x)}
   \qquad (S_l(0) = -l \text{ for a neutral pair}),

which AZURE2 uses for any channel that sits within :math:`10^{-10}` MeV of its
threshold (the Whittaker function at zero binding, the general closed-channel
route, is ill-conditioned there).  Every other channel of the
:math:`J^\pi` group takes its usual values at that total energy: a closed
channel enters through its shift function, an open one -- a pair with a lower
threshold -- through :math:`S + iP`, and then :math:`M_{cc}` and :math:`a_0`
are complex, with :math:`\mathrm{Im}\,a_0 < 0` describing the flux lost to the
open channel.

Channels and combinations
-------------------------

For :math:`l = 0` the total spin equals the channel spin, :math:`J = s`, and the
parity is :math:`\pi_1\pi_2`, so a pair of spins :math:`j_1, j_2` has one
scattering length per channel spin :math:`s = |j_1 - j_2|, \ldots, j_1 + j_2`
(for p + :sup:`7`\ Be the :math:`a_{01}` and :math:`a_{02}` of the
literature).  An identical pair has only the channel spins that Bose or Fermi
symmetry allows, :math:`s` even.  A channel spin in which no level has a width
gets :math:`M_{cc} = 0`: the pure Coulomb (hard-sphere) value, flagged in the
output.

Neutron measurements usually give the *coherent* scattering length, the
spin-weighted average

.. math::

   a_{\mathrm{coh}} = \sum_s g_s\, a_s, \qquad g_s = \frac{2s+1}{(2j_1+1)(2j_2+1)},

which a constraint line selects with :math:`s = -1`.  Tables of neutron
scattering lengths (e.g. Sears) list *bound* lengths :math:`b`; the free
scattering length of the pair is :math:`a = b\,A/(A+1)` for a target of mass
number :math:`A`.

In the calculation
------------------

* **Output.**  ``output/scatteringLengths.out`` lists every
  channel of every particle pair (key, particles, :math:`J^\pi`, :math:`x`,
  :math:`\mathrm{Re}\,a_0`, :math:`\mathrm{Im}\,a_0`) at the parameters of the
  run, then the constraints with their :math:`\chi^2`.
* **Constraints.**  Each active line of the ``<scatteringLengths>`` block
  (:doc:`../reference/data_formats`) adds
  :math:`((\mathrm{Re}\,a_0 - a_{\mathrm{exp}})/\sigma)^2`, and the same for
  the imaginary part when its uncertainty is given, to the :math:`\chi^2` every
  fitter minimizes and to the MCMC likelihood.  ``chiSquared.out`` reports the
  sum as ``Total-Scattering-Length-Chi-Squared``.
* **Gradient.**  The derivative of :math:`a_0` with respect to the level
  energies and reduced widths is exact: it is the reverse mode of
  :math:`\tilde\gamma^T\tilde A\tilde\gamma` through the same level-matrix
  adjoint the cross sections use (``AMatrixFunc::BilinearAdjoint``, on the
  channel-space or the level route, :doc:`../developer/level_matrix`), so
  ``--use-gradient``, ``--use-lm`` and ``--use-gsl-lm`` treat a constraint like
  a data point (one residual row per constrained part).
* **pyazr.**  :meth:`azure2.scattering_lengths` returns the values at any
  parameter vector, so a posterior of :math:`a_0` follows from an MCMC chain;
  :class:`AzrModel` reads and writes the block.

Validation
----------

``tests/scattering_length`` checks the threshold evaluation four ways:

================================  =============  =============================
Case                              AZURE2 (fm)    Reference (fm)
================================  =============  =============================
:sup:`3`\ He + :sup:`4`\ He,      33.0996        33.10 (Paneru et al. 2024,
:math:`1/2^+`, Paneru 2024                       Table V parameters)
Table V
p + :sup:`7`\ Be, :math:`a_{01}`  17.0678        17.0678 (independent
                                                 calculation, mpmath)
p + :sup:`7`\ Be, :math:`a_{02}`  -5.0184        -5.0184 (same)
n + :sup:`16`\ O, one level       2.62136        :math:`a(1-\tilde\gamma^2/\tilde E)`
================================  =============  =============================

Independently of the threshold code, the same lengths follow from AZURE2's own
phase shifts at low energy by fitting :math:`K(E)` in :math:`k^2` and reading
:math:`-1/K(0)`, to within :math:`4\times 10^{-5}`.  The standard
parametrization gives the same values as Brune's to :math:`10^{-8}`, and the
analytic gradient of the constraint term agrees with finite differences to
:math:`4\times 10^{-10}` (standard, Brune) and :math:`4\times 10^{-8}` (Park).

The p + :sup:`7`\ Be values are those of the parameters printed in Table II of
Paneru et al. (2019); that paper quotes :math:`a_{01} = 17.34` and
:math:`a_{02} = -3.18` fm.  The first agrees within its uncertainty; the
second does not follow from the printed :math:`2^-` parameters (the
:math:`a_{02}` of a single level depends on nothing else), so the table entry or
the quoted value is not the one the fit used.

Not yet included: scattering volumes and effective ranges (:math:`l > 0`, and
:math:`r_0`).

References
----------

* S. N. Paneru et al., Phys. Rev. C **99**, 045807 (2019).
* S. N. Paneru et al., Phys. Rev. C **109**, 015802 (2024).
* C. R. Brune, Phys. Rev. C **66**, 044611 (2002).
* H. van Haeringen, *Charged-Particle Interactions* (Coulomb Press Leyden,
  1985), the Coulomb-modified effective-range expansion.
* V. F. Sears, Neutron News **3**, 26 (1992), neutron scattering lengths.
