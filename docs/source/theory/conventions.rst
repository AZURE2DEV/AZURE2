Angular-Momentum and Phase Conventions
======================================

This chapter states the angular-momentum coupling, the phase factors and the
sign conventions with which the present engine turns the level parameters of a
project into collision-matrix elements and observables. Each convention was
read from the code, and the file and line that fix it are given after the
equation. A reader needs these conventions to compare reduced-width amplitudes
with those of another code or of a structure calculation, to set the relative
signs of interfering levels by hand, and to know which parameter values are
equivalent. Where the engine departs from Lane and Thomas, Rev. Mod. Phys. 30
(1958) 257 (LT below), or from Descouvemont and Baye, Rep. Prog. Phys. 73
(2010) 036301 (DB below), the difference is stated. The last section describes
the ``Convention`` record of the upgrade's ``physics/`` library and the
conversions that use it. The engine does not use the record yet.

All energies are in MeV and lengths in fm. :math:`E` is the c.m. energy of the
entrance pair plus its threshold, i.e. the excitation energy of the compound
nucleus (``inEnergy`` in the code). A channel energy is
:math:`E_c = E - S_c - E^*_c` (MeV, c.m.), with :math:`S_c` the separation
energy and :math:`E^*_c` the excitation energy of the pair. Line numbers refer
to ``feature/upgrade`` at f8b7bf22.

Channels and coupling
---------------------

A partition :math:`\alpha` is a pair of particles: particle 1 is the light one
(spin :math:`j_1`, the ``j1`` column), particle 2 the heavy one (:math:`j_2`).
A particle channel is :math:`c = (\alpha, l, s)`, coupled as

.. math::
   :label: conv-coupling

   \mathbf s = \mathbf j_1 + \mathbf j_2, \qquad
   \mathbf J = \mathbf s + \mathbf l, \qquad
   |(j_1 j_2) s\,\nu\rangle = \sum_{m_1 m_2} (j_1 m_1\, j_2 m_2 | s\,\nu)\,|j_1 m_1\rangle |j_2 m_2\rangle ,

with :math:`s` running from :math:`|j_1 - j_2|` to :math:`j_1 + j_2`
(``CNuc.cpp:1440``, ``PolarizationFunc.cpp:50-59``).

- **Channel spin, light particle first.** The engine couples particle 1 first.
  This is the order of LT Sec. III.2a, which the code quotes
  (``PolarizationFunc.cpp:333-343``; ``ParticleCG``, lines 21-38). The same
  order enters the M1 spin operator of external capture
  (``ECIntegral.cpp:436-452``) and the secondary-:math:`\gamma` (UPOS)
  coefficients (``CNuc.cpp:1918-1937``, whose comment records that Eq. (22) of
  Brune and deBoer, PRC 102, 024628 (2020), as printed, couples the heavy
  particle first). Coupling the other way multiplies every amplitude of
  channel spin :math:`s` by :math:`(-1)^{j_1+j_2-s}`.
  ``AZURE2_SWAP_COUPLING_ORDER=1`` reverses the order in the analyzing power
  only, as a convention-invariance test.
- **Orbital coupling, channel spin first.** The spin-amplitude route of the
  analyzing powers couples :math:`(s\,\nu\; l\,0 | J\,\nu)`
  (``PolarizationFunc.cpp:115-127``). The unpolarized cross section uses the
  Blatt–Biedenharn coefficients (Eq. :eq:`conv-z`). The two routes give
  cross sections whose ratio is constant in angle to :math:`10^{-9}`
  (:doc:`polarization_implementation`, Gate 1). That check fixes the relative
  orbital–spin phase of the two routes; it does not compare them with a third
  convention.
- **Photon channels.** For a :math:`\gamma` + nucleus pair the channel's
  :math:`l` is the multipolarity :math:`L`. The type is electric when
  :math:`\pi_\lambda \pi_f = (-1)^L` and magnetic otherwise
  (``AChannel.cpp:17-22``).
- **Clebsch–Gordan and Racah coefficients.** They follow Condon and Shortley:
  :math:`(j_1 m_1 j_2 m_2|JM) = (-1)^{j_1-j_2+M}\sqrt{2J+1}\,
  \begin{pmatrix} j_1 & j_2 & J\\ m_1 & m_2 & -M\end{pmatrix}`
  (``AngCoeff.cpp:9-21``). The Racah coefficient is
  :math:`W(abcd;ef) = (-1)^{a+b+c+d}\{a\,b\,e;\,d\,c\,f\}`
  (``AngCoeff.cpp:23-34``). The 3j, 6j and 9j symbols are GSL's.
  :math:`Y_l^m` carries the Condon–Shortley phase (``AngCoeff.cpp:52-69``).
  :math:`P_l^1` does not, following Seyler and Weller, PRC 20 (1979) 453
  (``AngCoeff.cpp:44-50``).

No factor :math:`i^l` appears in the particle channels. It is absent from the
channel functions, from :math:`U` and from the angular coefficients, which are
real (Eq. :eq:`conv-z`). An explicit :math:`i^{l}` enters in one place only:
the electric external-capture amplitude (Eq. :eq:`conv-ec-e`). The spin-amplitude
route multiplies every pathway by :math:`i\sqrt{2l+1}`, the :math:`i` common to
all of them (``PolarizationFunc.cpp:125-126``).

Channel functions
-----------------

For a particle channel with :math:`E_c > 0`, at the channel radius
:math:`a_c` (fm, one per partition), with :math:`\rho = k_c a_c`,
:math:`\eta` the Sommerfeld parameter and :math:`F_l, G_l` the regular and
irregular Coulomb functions,

.. math::
   :label: conv-ps

   P_c = \frac{\rho}{F_l^2 + G_l^2}, \qquad
   S_c = \rho\,\frac{F_l F_l' + G_l G_l'}{F_l^2 + G_l^2}, \qquad
   L^0_c = S_c - B_c + i P_c

(``CoulFunc.cpp:252-291``; ``EPoint.cpp:1496-1500``). Below threshold,
:math:`P_c = 0` and :math:`S_c = \rho\,W'/W` with the Whittaker function
:math:`W_{-\eta,l+1/2}(2\kappa r)` (``ShftFunc.cpp:61-78``). Just above a
charged-particle threshold, where :math:`P_c` vanishes in double precision
(:math:`\eta > 100`), :math:`S_c` is continued through :math:`E_c = 0`
(``ChannelFunc.cpp:183-194``).

The phases are

.. math::
   :label: conv-phases

   e^{-i\phi_c} = \frac{G_l - iF_l}{\sqrt{F_l^2 + G_l^2}}, \qquad
   e^{i\omega_c} = \prod_{n=1}^{l} \frac{n + i\eta}{\sqrt{n^2+\eta^2}}
   = e^{i(\sigma_l - \sigma_0)}, \qquad
   \Omega_c = e^{i(\omega_c - \phi_c)}

(``CoulFunc.h:39-45``; ``EPoint.cpp:1505-1508``).

- **Hard-sphere phase.** :math:`\phi_c = \arctan(F_l/G_l)` at
  :math:`\rho`, as in LT. The hard-sphere phase *shift* is :math:`-\phi_c`;
  pyazr's ``coulomb_functions`` returns this value as ``delta_hs``
  (``AZUREAPI.cpp:1871``).
- **Coulomb-phase reference.** :math:`\omega_c = \sigma_l - \sigma_0`, as in
  LT. All amplitudes are therefore relative to the common factor
  :math:`e^{2i\sigma_0}`, which drops out of every observable of one pair.
- **Photon channels.** :math:`L^0_c = 0`, :math:`\Omega_c = 1` and
  :math:`P_c^{1/2} = (E_\gamma/\hbar c)^{L+1/2}`, with
  :math:`E_\gamma = E - E_f` (MeV) and no recoil correction
  (``EPoint.cpp:1522-1528``). With ``--use-rmc``, :math:`P_c = 1`.
- **β-decay channels.** :math:`P_c^{1/2}` is the square root of the
  integrated Fermi function (``EPoint.cpp:1529-1538``).

Boundary condition
^^^^^^^^^^^^^^^^^^

In the standard (``--no-brune``) parametrization, every particle channel of a
:math:`J^\pi` group has the constant boundary condition

.. math::
   :label: conv-bc

   B_c = S_c(E_1),

where :math:`E_1` is the energy of the **first level listed** in the file for
that :math:`J^\pi` that is an R-matrix level (``boundary_values``,
``physics/src/LevelParameters.cpp:305-324``, called by ``CNuc::CalcBoundaryConditions``). A group without levels has
:math:`B_c = 0`, which then has no effect on the result. Observables do not
depend on :math:`B_c` as long as the transformation of the parameters
(below) is exact. The formal parameters do depend on :math:`B_c`, so they
depend on the order of the levels in the file. In the Brune and Park
parametrizations :math:`B_c` cancels from the level matrix
(``AMatrixFunc.cpp:302-311``). LT leave :math:`B_c` free. DB usually take
:math:`B_c = 0` in their examples; the R-matrix is defined the same way.

Reduced widths and their signs
------------------------------

The level matrix is built from real reduced-width amplitudes
:math:`\gamma_{\lambda c}`. Their unit is MeV\ :sup:`1/2` for particle
channels and MeV\ :sup:`1/2` fm\ :sup:`L+1/2` for photon channels. In the
standard parametrization,

.. math::
   :label: conv-amatrix

   (A^{-1})_{\lambda\mu} = (E_\lambda - E)\,\delta_{\lambda\mu}
   - \sum_c \gamma_{\lambda c}\gamma_{\mu c}\,L^0_c ,

as in LT (``AMatrixFunc.cpp:270-322``). The Brune form replaces
:math:`B_c` by the level-dependent shift terms :math:`S_c(\lambda,\mu)` of
:doc:`../developer/level_matrix` (lines 302-311), and the Park form uses
Park's overlap :math:`J` (lines 294-301 and 313-317).

Which basis the numbers in ``<levels>`` belong to is set on the command line
(:doc:`../reference/command_line`). The conversion at start-up is
``CNuc::TransformIn``, which calls ``level_parameters``
(``physics/src/LevelParameters.cpp``; ``docs/upgrade/level_parameters.md``). For the Brune basis (the
default) it is

.. math::
   :label: conv-brune

   \tilde\gamma_{\lambda c}^2 =
   \frac{|\Gamma_{\lambda c}|}{2P_c(E_\lambda)}\;
   \frac{1 + \sum_{c'\in\mathrm{rwa}} \tilde\gamma_{\lambda c'}^2 S'_{c'}(E_\lambda)}
        {1 - \sum_{c'\notin\mathrm{rwa}} |\Gamma_{\lambda c'}|\,S'_{c'}(E_\lambda)/(2P_{c'}(E_\lambda))},
   \qquad
   \operatorname{sign}\tilde\gamma_{\lambda c} = \operatorname{sign}\Gamma_{\lambda c}

(``amplitudes_from_observed``, ``LevelParameters.cpp:103-242``). Here :math:`\Gamma_{\lambda c}` is the signed
input (eV; the code divides it by :math:`10^6` to MeV), :math:`S' = dS/dE` (MeV\ :sup:`-1`), the
sums run over particle channels only, and "rwa" denotes the channels flagged
``gammaIsRWA`` (``input_is_rwa`` in the API's parameter records). The input of
such a channel is :math:`\tilde\gamma_{\lambda c}` itself, in MeV\ :sup:`1/2`,
with its sign.

- **Open particle channel.** The input is a signed observed partial width
  (eV). The sign of the input becomes the sign of the amplitude.
- **Bound particle channel** (:math:`E_\lambda` below the threshold of the
  channel). The input is a signed ANC :math:`C` (fm\ :sup:`-1/2`). In
  Eq. :eq:`conv-brune` it is converted with
  :math:`|\Gamma|/(2P) \to \hbar^2 C^2 W_l^2(2\kappa a)/(2\mu a)`
  (``LevelParameters.cpp:161-176``). The sign of :math:`C` becomes the sign of
  :math:`\tilde\gamma`.
- **Photon channel.** The input is a signed :math:`\Gamma_\gamma` (eV), with
  :math:`P = (E_\gamma/\hbar c)^{2L+1}` and the same particle-channel factor
  of Eq. :eq:`conv-brune` (``LevelParameters.cpp:178-198``). When external capture is
  on, the input is the total amplitude: the internal amplitude is
  :math:`\gamma^{\rm int} = \operatorname{sign}(\Gamma)\,[\tilde\gamma^2 -
  (\mathrm{Im}\,\gamma^{\rm ext})^2]^{1/2} - \mathrm{Re}\,\gamma^{\rm ext}`
  (``subtract_external_capture``, ``LevelParameters.cpp:328-352``), with :math:`\gamma^{\rm ext}` from
  ``CNuc::CalcExternalWidth`` (``CNuc.cpp:3196-3278``).
- **Static moments.** A photon channel of a level at the energy and
  :math:`J^\pi` of the final state carries a static moment (M1 or E2 only):
  a penetrability factor of its own applies (:math:`3J/[4(J+1)\mu_N^2]` for
  M1), and the sign is reversed for half-integer :math:`J`
  (``LevelParameters.cpp:185-195``).
- **β-decay channels** (F, GT). The input is the amplitude itself.
- **Park basis** (``--use-park``). The amplitudes are the observed ones,
  :math:`\Gamma = 2P\gamma^2`. The same sign rule applies; a ``gammaIsRWA``
  input is read as Brune's amplitude and multiplied by :math:`\sqrt{J}`
  (``LevelParameters.cpp:210-240``).
- **Formal basis** (``--no-brune``). The Brune amplitudes are transformed to
  formal ones by the generalized eigenproblem :math:`N a = E M a` of Brune,
  PRC 66, 044611 (2002) (``formal_from_brune``, ``LevelParameters.cpp:356-434``; the
  engine passes GSL's eigensolver). GSL fixes the sign of each
  eigenvector, so the overall sign of a formal level's amplitude vector is
  arbitrary. The relative signs of the channels of one level are preserved,
  and observables are even in the overall sign of a level
  (:math:`\gamma^T A \gamma` is invariant under it). pyazr's ``save_fit``
  check accordingly compares amplitudes "up to a common sign of the level"
  (:doc:`../user_guide/pyazr`).
- **No transformation** (``--no-transform``). The ``<levels>`` values are
  used as :math:`E_\lambda` and :math:`\gamma_{\lambda c}` directly.

On output (``CNuc::TransformOut``) the reported width carries the sign of the
amplitude, and for photon channels with external capture the sign of the real
part of the total amplitude (``CNuc.cpp:2483-2499``).

Collision matrix and T-matrix
-----------------------------

For particle and photon channels :math:`c, c'` of one :math:`J^\pi`,

.. math::
   :label: conv-u

   U_{cc'} = \Omega_c\Omega_{c'}\Big[\delta_{cc'} + 2i\,P_c^{1/2}P_{c'}^{1/2}
   \sum_{\lambda\mu}\gamma_{\lambda c}A_{\lambda\mu}\gamma_{\mu c'}\Big],
   \qquad
   T_{cc'} = e^{2i\omega_c}\delta_{cc'} - U_{cc'}

(``AMatrixFunc.cpp:540-556``; the R-matrix route ``RMatrixFunc.cpp:365-387``
uses the same phases). This is LT's level-matrix form of :math:`U` with LT's
:math:`T`. Let :math:`S^{\rm nuc}` be the nuclear collision matrix of the
textbook amplitude
:math:`f = f_C + (2ik)^{-1}\sum (2l+1)e^{2i\sigma_l}(S^{\rm nuc}_l - 1)P_l`
(spinless notation). AZURE2's :math:`U` is
:math:`e^{2i\omega_c}S^{\rm nuc}`. The phase shift a segment reports is
:math:`\delta = \tfrac12\arg(U_{cc}e^{-2i\omega_c})` in degrees, in
:math:`(-90°, 90°]`, and in :math:`[0°, 180°)` for identical particles
(``GenMatrixFunc.cpp:517-566``). It includes the hard-sphere part
:math:`-\phi_c`.

With ``--use-rmc`` (Reich–Moore), the photon channels add
:math:`-i\sum_\gamma\gamma_{\lambda\gamma}^2` to the diagonal of
:math:`A^{-1}` (``AMatrixFunc.cpp:290-293``), the capture T-matrix is
:math:`-i\cdot 2iP_c^{1/2}(\gamma^TA\gamma)P_{c'}^{1/2}` without the phases
:math:`\Omega` (line 549), and the capture cross section is the unitarity
deficit of the particle channels (``GenMatrixFunc.cpp:54-101``). Without it, photon channels have :math:`L^0 = 0`:
**radiative widths do not enter the level matrix.** A level's total width
there is the sum of its particle widths. The photon widths also stay out of
the Brune factor of Eq. :eq:`conv-brune`.

Photon channels
^^^^^^^^^^^^^^^

The internal capture amplitude is Eq. :eq:`conv-u` with :math:`c'` a photon
channel:

.. math::
   :label: conv-int

   T^{\rm int}_{c\gamma} = -2i\,\Omega_c\,P_c^{1/2}
   \left(\frac{E_\gamma}{\hbar c}\right)^{L+1/2}
   \sum_{\lambda\mu}\gamma_{\lambda c}A_{\lambda\mu}\gamma_{\mu\gamma}.

The external-capture amplitude to a final level :math:`f` with particle
channel :math:`c_f = (\alpha, l_f, s)` is (``EPoint.cpp:1946-2005``,
``ECIntegral.cpp:316-459``, ``AMatrixFunc.cpp:559-583``)

.. math::
   :label: conv-ext

   T^{\rm ext}_{c\gamma} = \gamma_{f c_f}\,\mathcal N_f^{1/2}
   \frac{(2\mu a/\hbar^2)^{1/2}}{W_{l_f}(2\kappa_f a)}\;
   \Omega_c\left(\frac{E_\gamma}{\hbar c}\right)^{L+1/2}
   \mathcal A_{L}\;\mathcal I ,

with the radial overlap, for a non-channel-capture pathway,

.. math::

   \mathcal I = \frac{(\mu c^2/2E_c)^{1/4}}{(\hbar c)^{1/2}}
   \int_a^\infty \big[\cos\phi_c\,F_l(kr) - \sin\phi_c\,G_l(kr)\big]\,
   W_{l_f}(2\kappa_f r)\,r^{L}\,dr

(:math:`r^0` for M1; the exact operator of
:doc:`external_capture_operator` replaces the effective charge by a function of
:math:`r`). The bracket is the hard-sphere wave function, which vanishes at
:math:`r = a`. :math:`\mathcal N_f` is the bound-state normalization
(``SqrtNFFactor``), so the prefactor is the ANC of the final level, with the
sign of :math:`\gamma_{fc_f}`. The angular factors are

.. math::
   :label: conv-ec-e

   \mathcal A_{EL} = -i\,\bar e_L\,
   \frac{[8(2L+1)(L+1)/L]^{1/2}}{(2L+1)!!}\;
   i^{\,l + L - l_f}\,(l\,0\,L\,0|l_f\,0)\,\hat l\,\hat J_f\,
   W(L\,l_f\,J\,s;\,l\,J_f),
   \qquad
   \bar e_L = (\alpha\hbar c)^{1/2}\Big[Z_1\Big(\frac{m_2}{M}\Big)^L
   + Z_2\Big(-\frac{m_1}{M}\Big)^L\Big],

.. math::
   :label: conv-ec-m

   \mathcal A_{M1} = +i\,\frac{\alpha^{1/2}(\hbar c)^{3/2}}{2m_pc^2}
   \Big(\frac{16}{3}\Big)^{1/2}\hat J_f\,
   \big[\delta_{s s_f}\,\bar\mu\,\{l(l+1)(2l+1)\}^{1/2} W(1\,l\,J\,s;\,l\,J_f)
   + \text{spin term}\big],

with :math:`\hat x = (2x+1)^{1/2}`, :math:`J` the spin of the capturing
:math:`J^\pi` and :math:`\bar\mu = \mu\,m_p(Z_1/m_1^2 + Z_2/m_2^2)` (masses
in u). The spin term is :math:`\langle (j_1j_2)s_f\|g_1\mathbf S_1 +
g_2\mathbf S_2\|(j_1j_2)s\rangle` in the coupling order of
Eq. :eq:`conv-coupling` (``ECIntegral.cpp:436-458``).

- **E/M multipole phases.** The electric amplitude carries
  :math:`-i\,i^{\,l+L-l_f}` and the magnetic one :math:`+i`. Their relative
  phase is pinned against an m-scheme calculation by ``tests/ec_e1m1_phase``.
- **Sign of external relative to internal capture.** The two amplitudes add
  coherently in every T-matrix sum with the same :math:`(J, l, L)`
  (``GenMatrixFunc.cpp:136-218``). Their relative sign is that of
  :math:`\gamma_{fc_f}` (the final ANC) times
  :math:`\mathcal A_L\mathcal I` in Eq. :eq:`conv-ext`, against that of
  :math:`\gamma_{\lambda c}\gamma_{\mu\gamma}` times :math:`-2i\Omega_c` in
  Eq. :eq:`conv-int`. A positive :math:`\Gamma_\gamma` therefore means
  constructive or destructive interference with external capture depending
  on the ANC sign and the angular factors. No test compares this relative sign
  with an independent calculation (see "Not established" below).
- **Channel capture** (``isChannelCapture``) replaces the hard-sphere wave in
  :math:`\mathcal I` by :math:`-\tfrac{i}{2}P_c^{1/2}e^{-i\phi_c}(G_l + iF_l)`
  and multiplies the pathway by :math:`2iP^{1/2}\gamma^TA\gamma` of the
  particle pathway
  (``ECIntegral.cpp:361-366``; ``AMatrixFunc.cpp:573-582``).

Observables
-----------

Let :math:`k` be the entrance wave number (fm\ :sup:`-1`) and
:math:`g_I = [(2j_1+1)(2j_2+1)]^{-1}`. Cross sections are formed in
fm\ :sup:`2` and divided by 100 to give b.

Angle-integrated cross section (``GenMatrixFunc.cpp:128-257``):

.. math::
   :label: conv-sigma

   \sigma_{\alpha\to\alpha'} = \frac{\pi}{k^2}\,g_I\,(1+\delta_{12})
   \sum_{s s'}\sum_{J l l'} (2J+1)\,
   \Big|\sum_{\rm pathways} T\Big|^2 ,

where the inner sum adds internal and external amplitudes with the same
:math:`(s, s', J, l, l')`; for capture :math:`l'` is the multipolarity.
:math:`(1+\delta_{12}) = 2` for identical particles in the entrance pair and
1 otherwise. The total cross section uses the optical theorem,
:math:`\sigma_{\rm tot} = (\pi/k^2)\,g_I\sum (2J+1)\,2\,\mathrm{Re}\,T_{cc}`
(lines 103-123).

Angular distributions (Blatt–Biedenharn; ``CNuc.cpp:1833-1982``,
``GenMatrixFunc.cpp:261-426``):

.. math::
   :label: conv-ad

   \frac{d\sigma}{d\Omega} = \frac{g_I}{k^2}\sum_{s s'}\sum_L
   \sum_{\rm pairs} z_L\,\mathrm{Re}\big[T_1 T_2^*\big]\,P_L(\cos\theta),

.. math::
   :label: conv-z

   Z(l_1J_1l_2J_2; sL) = \hat l_1\hat l_2\hat J_1\hat J_2\,
   (l_1 0\, l_2 0|L 0)\,W(l_1J_1l_2J_2; sL),

with :math:`z_L = \tfrac14(-1)^{s'-s}Z(l_1J_1l_2J_2;sL)\,Z(l_1'J_1l_2'J_2;s'L)`
for a particle exit pair. For a photon exit,
:math:`z_L = \tfrac14(-1)^{1+s-j_f}Z(l_1J_1l_2J_2;sL)\,\bar Z`, with
:math:`\bar Z = \hat L_1\hat L_2\hat J_1\hat J_2\,(L_1 1\, L_2\, {-1}|L 0)
\,W(L_1J_1L_2J_2; j_fL)`, where :math:`j_f` is the final spin. The orders
allowed are those with :math:`l_1+l_2+L` even and
:math:`l'_1+l'_2+w_1+w_2+L` even, where :math:`w = 1` for E and 0 for M and
particle exits (``CNuc.cpp:1903-1916``). Angular-distribution segments report
:math:`a_k` in :math:`d\sigma/d\Omega = (\sigma/4\pi)[1+\sum_k a_kP_k]`
(``GenMatrixFunc.cpp:427-438``). Secondary-γ (UPOS) coefficients follow Brune
and deBoer (2020) Eq. (22) in the coupling order of Eq. :eq:`conv-coupling`
(``CNuc.cpp:1918-1937``).

Elastic scattering adds the Coulomb amplitude and its interference
(``EPoint.cpp:1890-1934``, ``GenMatrixFunc.cpp:439-516``):

.. math::
   :label: conv-elastic

   \frac{d\sigma}{d\Omega} = \frac{\pi}{k^2}|C_\alpha(\theta)|^2
   + \frac{g_I}{k^2}\sum z_L\,\mathrm{Re}[T_1T_2^*]P_L
   + \frac{\pi}{k^2}\,\frac{1}{\sqrt\pi}\,
   \mathrm{Re}\Big[i\,C_\alpha(\theta)\sum_{J s l} g_I(2J+1)\,T^*_{cc}\,P_l(\cos\theta)\Big],

.. math::
   :label: conv-coulomb

   C_\alpha(\theta) = \frac{\eta}{2\sqrt\pi}\,\csc^2\frac\theta2\,
   \exp\Big(-2i\eta\ln\sin\frac\theta2\Big).

The physical Coulomb amplitude, relative to :math:`e^{2i\sigma_0}`, is
:math:`f_C = -\eta/(2k\sin^2\tfrac\theta2)\,
e^{-i\eta\ln\sin^2(\theta/2)} = -(\sqrt\pi/k)\,C_\alpha`. With :math:`T` of
Eq. :eq:`conv-u`, the last term of Eq. :eq:`conv-elastic` equals
:math:`2\,\mathrm{Re}[f_C^* f_N]`, so the interference has the physical sign.
``tests/coulomb_nuclear_sign`` checks this against an mpmath partial-wave sum
for p + :sup:`12`\ C hard-sphere scattering: 30 points agree to
:math:`8\times10^{-11}` on both matrix routes and in both parametrizations.
Reversing the interference term changes the cross section by 7–74 % on that
grid.

**Identical particles.** Channels with :math:`l + s` odd are refused for
pairs with spin and warned about for spin 0 (``CNuc.cpp:671-728, 947-980``).
In elastic scattering the nuclear term is multiplied by 4 and the
interference by 2. The Coulomb amplitude of channel spin :math:`s` is
:math:`C_\alpha(\theta) + (-1)^s C_\alpha(\pi-\theta)`; for spin-0 bosons it
is the symmetric sum (``EPoint.cpp:616-624, 1903-1925``;
``GenMatrixFunc.cpp:439-514``). A reaction with an identical entrance pair
takes the factor 2 of Eq. :eq:`conv-sigma`. A reaction with an identical exit
pair takes no factor: it counts events, not particles. The sign :math:`+1` for
integer and :math:`-1` for half-integer :math:`j` is set in ``PPair.cpp:60-72``.

Comparison with Lane–Thomas and Descouvemont–Baye
-------------------------------------------------

The comparison below is with the published forms as AZURE2's sources cite
them (LT Sec. III.2a is quoted in ``PolarizationFunc.cpp``). The paper texts
were not re-read for this page.

- **The same as LT:** channel spin :math:`\mathbf s = \mathbf j_1 + \mathbf j_2`
  with Condon–Shortley coefficients; :math:`P_c`, :math:`S_c`,
  :math:`\phi_c`, :math:`\omega_c = \sigma_l-\sigma_0` and
  :math:`\Omega_c`; the level matrix Eq. :eq:`conv-amatrix`; :math:`U` and
  :math:`T = e^{2i\omega}-U` of Eq. :eq:`conv-u`; real :math:`\gamma_{\lambda c}`
  with :math:`\Gamma = 2P\gamma^2` for formal widths.
- **Different from LT:**
  (i) the default parameters are Brune's (and optionally Park's), not LT's
  formal ones. The input is a signed observed width or ANC, and the sign of
  the input defines the sign of the amplitude;
  (ii) :math:`B_c = S_c(E_1)` at the first listed level (Eq. :eq:`conv-bc`),
  not a free constant;
  (iii) photon channels have :math:`L^0 = 0` and do not enter :math:`A`
  except with ``--use-rmc``;
  (iv) external capture is added outside the R-matrix (Eq. :eq:`conv-ext`),
  with its own multipole phases.
- **Relative to DB:** the R-matrix :math:`R_{cc'} = \sum_\lambda
  \gamma_{\lambda c}\gamma_{\lambda c'}/(E_\lambda - E)` has the same form.
  If DB's collision matrix is the nuclear :math:`S^{\rm nuc}` of the textbook
  amplitude (not checked), AZURE2's differs from it by the Coulomb phase,
  :math:`U^{\rm AZURE2} = e^{2i(\sigma_l-\sigma_0)}S^{\rm nuc}`
  (Eq. :eq:`conv-u`). AZURE2's spin-amplitude route couples in the order
  :math:`[[\phi_1\otimes\phi_2]^s\otimes Y_l]^J` and its channel functions
  carry no :math:`i^l`. Whether DB use the same order and phase was not
  checked.

Convention record
-----------------

.. note::

   Implemented in the upgrade's ``physics/`` library
   (``physics/include/azure2/physics/Convention.h``,
   ``LevelParameters.h``; ``docs/upgrade/level_parameters.md``). The engine
   does not use it yet: it computes in the single convention described above,
   and nothing in a project file states a convention.

A ``Convention`` record in the ``System`` states what a project's level
parameters mean. ``level_parameters`` converts from it to the internal
convention above. On the :math:`\gamma` columns of a :math:`J^\pi` block the
conversion is an orthogonal block matrix of 6j symbols and signs; it does not
depend on the parameters, so it is its own derivative. The fields, their
values and the internal value:

.. list-table::
   :widths: 25 35 40
   :header-rows: 1

   * - Field
     - Values
     - Internal value (this page)
   * - ``coupling``
     - ``LS`` (channel spin), ``JJ``
     - ``LS``
   * - ``spin_order``
     - ``LightFirst``, ``HeavyFirst`` (jj: the particle whose spin couples
       with :math:`l`)
     - ``LightFirst``: :math:`(j_1 m_1 j_2 m_2|s\nu)`
   * - ``orbital_order``
     - ``SpinFirst``, ``OrbitalFirst``
     - ``SpinFirst``: :math:`(s\nu\,l m|JM)`
   * - ``orbital_phase``
     - ``None``, ``IPowL``, ``IPowMinusL``
     - ``None`` on particle channels
   * - ``electric_sign``, ``magnetic_sign``
     - :math:`\pm 1` on the E or M photon amplitudes
     - :math:`+1` (the phases of Eqs. :eq:`conv-ec-e`, :eq:`conv-ec-m`)
   * - ``coulomb_phase``
     - ``SigmaLMinusSigma0``, ``SigmaL``
     - ``SigmaLMinusSigma0`` (no change of :math:`\gamma`)
   * - ``basis``
     - ``Formal``, ``Brune``, ``Park``
     - ``Brune``
   * - ``boundary``
     - ``FirstListedLevel``, ``LowestLevel``
     - ``FirstListedLevel`` (formal basis only)
   * - ``formal_sign``
     - ``AsComputed``, ``LargestPositive``
     - ``AsComputed`` (the eigensolver's sign)

The proposals of this page (B\ :sub:`c` from the lowest level; the largest
:math:`|\gamma|` of a formal level positive) are the options ``LowestLevel``
and ``LargestPositive``. They are off by default and change no result of the
engine. The width sign (the sign of the input :math:`\Gamma`, ANC or
amplitude becomes the sign of the amplitude) is not an option: it is the
definition of a signed width. A width written in another code's convention
of signs is converted through the amplitudes (``electric_sign``,
``magnetic_sign``, ``orbital_phase``).

For jj coupling the recoupling coefficient is

.. math::

   \langle ((j_1 j_2)s, l)J \,|\, ((j_x l)j, j_y)J \rangle
   = (-1)^{4j_x + 3j_y + 3l - s - j}\,\hat s\,\hat j\,
   \begin{Bmatrix} j_x & j_y & s\\ J & l & j \end{Bmatrix}

with a further :math:`(-1)^{j_1+j_2-s}` when :math:`x` is the heavy particle
and :math:`(-1)^{l+j_x-j}` for the order :math:`(l\,j_x)j`; the
:math:`i^l` convention is represented by the sign :math:`(-1)^{\lfloor
l/2\rfloor}` (``docs/upgrade/level_parameters.md`` §5).
``physics_angular_coupling`` checks every coefficient against an m-scheme
overlap of Clebsch–Gordan products. ``physics_level_parameters_from_cnuc``
runs T1–T3 in six conventions: the formal parameters agree with those of the
internal convention to :math:`1.5\times10^{-15}` relative. This is a test at
the parameter level; the engine reads no convention yet. A channel missing
from a jj ↔ (l, s) recoupling is reported: T3's :math:`2^+` block lists only
:math:`s = 2` of the two channel spins its p + :sup:`7`\ Li and
n + :sup:`7`\ Be :math:`l = 1` channels allow, and has no jj form.

Not established
---------------

- **Sign of external relative to internal capture.** The formulas above fix it
  in the code, and ``tests/ec_e1m1_phase`` pins E against M within external
  capture. No test compares the sign of the internal–external interference
  (for a given sign of :math:`\Gamma_\gamma` and of the ANC) with an
  independent potential-model or m-scheme calculation.
- **The orbital–spin order** is established relative to the Blatt–Biedenharn
  route only (Gate 1 of :doc:`polarization_implementation`). No independent
  m-scheme calculation of an unpolarized cross section with mixed :math:`l`
  checks it.
- **The comparisons with LT and DB** rest on the formulas the code cites and
  on the standard forms of the level-matrix theory, not on a new reading of the
  papers (no PDF reader was available on the machine where this page was
  written). Equation and section numbers of LT and DB are therefore not
  given, except LT Sec. III.2a, which the code quotes.
