Trojan Horse (HOES) Observable
==============================

A segment whose observable code is 10 or more (``isDiff >= 10``) is a Trojan
Horse Method (THM) segment: its data are the half-off-energy-shell (HOES)
excitation function :math:`A(x,c)B` extracted from a three-body reaction
:math:`a + A \to c + B + s` with :math:`a = x + s`. The entrance pair carries
the binding energy :math:`B_{xs}` of the transferred particle in the optional
32nd field of its channel lines (the optional 33rd is the ``gammaIsRWA``
flag, which can only follow it). ``THMMatrixFunc::CalculateTHMCrossSection``
evaluates

.. math::

   \sigma_\mathrm{HOES}(E) \propto \sum_J (2J+1) \sum_{f} K(E)\, 2P_f(E)
   \sum_{s,l} \Bigl| \sum_{\lambda\lambda'} \gamma_{\lambda f}
   A_{\lambda\lambda'} V_{\lambda'}^{(s,l)} \Bigr|^2 ,
   \qquad
   V_{\lambda}^{(s,l)} = \sum_{c \in (s,l)} \gamma_{\lambda c}\, M_{l}(p, L_c),

with the half-off-shell momentum :math:`p = \sqrt{2\mu_{xA}(E + B_{xs})}/\hbar`
and the transfer form factor
:math:`M_l = (L_c - 1)\, j_l(pa) - pa\, j_l'(pa) \; [+\, C_l(E)]`
(Tumino et al., ARNPS 71 (2021) eqs. 50–51; Tribble et al., RPP 77 (2014)
106901 eqs. 2.76, 2.79). :math:`K(E)` collects the kinematic factors. The
choices in this expression are set by the ``<thm>`` block.

Options (``<thm>`` block)
-------------------------

The block is optional and may appear anywhere in the ``.azr`` file; the engine
uses the first line that starts with ``<thm>`` (``Config::ReadThmBlock``). One
``key=value`` per line, ``#`` starts a comment. A file without the block, or
with every key at its default, gives identical results
(``tests/thm_options``).

.. code-block:: text

   <thm>
   entranceL=incoherent     # coherent | incoherent
   vertex=constant          # constant | perlevel (alias real) | onshell
   kinematics=lacognata     # lacognata | triple | kf3body | lambda32
   coulombIntegral=0        # 0 | 1
   spectatorEnergy=0        # MeV, every THM entrance pair
   spectatorEnergy[5]=0.4   # MeV, entrance pair 5 only (overrides the above)
   </thm>

An unknown key, an unknown value, a negative spectator energy or a missing
``</thm>`` prints ``ERROR: <thm> ...``, and AZURE2 exits with a non-zero status
before any calculation. The GUI and ``pyazr.AzrModel`` have no editor for the
block but carry it through a save unchanged; the GUI writes it after
``</targetInt>``.

``entranceL`` — interference between entrance orbital momenta
   ``incoherent`` (default) sums the entrance partial waves of one channel spin
   :math:`s` incoherently in :math:`l`; ``coherent`` adds them in one amplitude,
   which is what AZURE2 did before September 2026 and what mrmpy does.

   With the quantization axis along :math:`\hat p_{xA}` the transferred-particle
   plane wave contributes only :math:`m_l = 0`. The HOES observable is
   integrated over the exit direction and summed over spin projections, and

   .. math::

      \sum_{M, m_s} \langle s\, m_s\, l\, 0 | J M \rangle
                    \langle s\, m_s\, l'\, 0 | J M \rangle
      = \delta_{ll'}\, \frac{2J+1}{2l+1},

   so the :math:`l \ne l'` cross terms vanish. The :math:`1/(2l+1)` cancels
   against :math:`|Y_{l0}(\hat z)|^2 = (2l+1)/4\pi` of the plane-wave
   expansion, leaving the :math:`l`-independent weight :math:`2J+1`. The
   coherent sum is exact only in the limit of the exit direction fixed along
   :math:`\hat p_{xA}`; it is kept for comparison with older fits.
   ``tests/7Li_p_a`` (two entrance :math:`l` in one channel spin): 2137.83
   incoherent, 3195.37 coherent, same parameters (2180.69 / 2740.48 with
   ``vertex=perlevel``, the pins before 2026-09-25).

``vertex`` — boundary term at the entrance vertex
   The transfer form factor :math:`M_l = (L_c - 1)\, j_l(pa) - pa\, j_l'(pa)`
   needs a boundary value :math:`L_c`:

   ``constant`` (default since 2026-09-25)
      a constant :math:`B_c = S_c(E_1)`, the shift at the lowest-energy level
      of the :math:`J^\pi` group, applied after the level sum, in either
      formalism. This is the vertex of the formal R-matrix at that
      :math:`B_c` — La Cognata's working formula (Tumino 2021 eq. 51) as the
      authors evaluate it. :math:`B_c` is taken from the lowest level, not
      from AZURE2's channel boundary condition (the first level read), so the
      result does not depend on the order of levels in the file; the vertex
      :math:`B_c` is free because :math:`\gamma^T A\gamma` does not depend on
      the R-matrix boundary.
   ``perlevel`` (alias ``real``; the default until 2026-09-25)
      the per-level shift :math:`S_c(E_\lambda)` of each level under the Brune
      parametrization (as mrmpy, ``vertex_boundary="per_level"``); without
      Brune it falls back to AZURE2's channel boundary condition (first
      level read).
   ``onshell``
      the logarithmic derivative of the on-shell outgoing wave,
      :math:`L_c(E) = S_c(E) + i P_c(E)` (Tribble 2014 eq. 2.76;
      Mukhamedzhanov et al., PRC 96 024623 (2017) eq. 28), independent of any
      boundary constant, which makes the vertex complex. It differs from the
      real choices where :math:`P_c` is not small next to :math:`S_c - L_c`,
      i.e. near and above the entrance barrier.

   *Why* ``constant``. The level matrix sandwiched between the widths,
   :math:`\gamma^T A \gamma`, is the same matrix in the Brune and in the
   formal representation — that is why every on-shell observable agrees
   between them. The HOES amplitude
   :math:`\sum_{\lambda\lambda'} \gamma_{\lambda f} A_{\lambda\lambda'}
   \gamma_{\lambda' c} M_l(p, L_c)` is of that form only if :math:`M_l` does
   not depend on the level: then :math:`M_l(p, B_c)` factors out of the level
   sum and the result is representation-independent, equal to a no-Brune run
   of the formal parameters (checked to 2e-11). With ``perlevel`` the factor
   :math:`M_l(p, S_c(E_\lambda))` sits inside the sum, and once two levels of
   one :math:`J^\pi` interfere the amplitude is not that of any single
   R-matrix representation: the same on-shell physics gives a different HOES
   curve. For a :math:`J^\pi` group with one level the two coincide
   (``tests/17O`` is unchanged).

   *Evidence*. ``tests/18O_p_a_thm`` reproduces La Cognata, Spitaleri &
   Mukhamedzhanov, ApJ 723 (2010) 1512 — 18O(p,α)15N, the interfering
   660 / 799 keV 1/2\ :sup:`+` doublet, their Table 3 formal parameters with
   :math:`B = S(E_1)`, 17 keV folding. Against the mid-line of their published
   band (one normalization, 0.505–0.895 MeV), ``constant`` gives 10 % rms
   overall and 5 % in the peak region 0.56–0.84 MeV; ``perlevel`` gives 27 % /
   23 %, with the 799 keV peak about 40 % low; ``onshell`` 31 % / 25 %. The
   test's ``check.sh`` asserts the first two.

   Pins that moved with the default (same parameters): ``tests/7Li_p_a``
   2180.69 → 2137.83, ``tests/6Li_d`` 1076.06 → 682.099; ``vertex=perlevel``
   reproduces the old values. ``tests/7Li_p_a`` with ``vertex=onshell``:
   2070.17.

``kinematics`` — factors left in the model by the data reduction
   HOES data are the measured triple-differential cross section divided by a
   kinematic factor (KF) and by the momentum distribution
   :math:`|\phi(p_{xs})|^2`. Whatever the KF did not remove must be in
   :math:`K(E)`. The choice is therefore fixed by what the experimenters
   divided by, which the paper states (usually next to the KF formula):

   ============== ============================= =========================================================
   value          :math:`K(E)`                  data were ...
   ============== ============================= =========================================================
   ``lacognata``  :math:`k_f/\mu_f`             fitted with La Cognata's working formula (Tumino 2021
                                                eq. 50); default, the historical behaviour
   ``triple``     1                             raw :math:`d^3\sigma / |\phi|^2`, no KF
                                                (Mukhamedzhanov 2017 eq. 34)
   ``kf3body``    :math:`1/(\mu_f k_f)`         divided by the full three-body phase-space KF
                                                (Typel & Baur, Ann. Phys. 305 (2003) eq. 16;
                                                Tumino 2021 eq. 15)
   ``lambda32``   :math:`1/k_i` (on shell)      divided by :math:`\mathrm{KF} = \lambda_3/\lambda_2`
                                                (Pizzone et al., PRC 83 045801 (2011);
                                                Tumino 2021 eq. 22)
   ============== ============================= =========================================================

   :math:`k_f`, :math:`\mu_f` are the exit-channel momentum and reduced mass,
   :math:`k_i` the on-shell entrance momentum at :math:`|E|`. The factors
   change the energy dependence, not only the scale, so the wrong choice
   distorts the relative heights of resonances. ``tests/7Li_p_a``:
   2137.83 / 2380.33 / 2624.63 / 5499.45 for the four values in table order
   (2180.69 / 2402.92 / 2632.81 / 5573.25 with ``vertex=perlevel``).

``coulombIntegral`` — external Coulomb term of the vertex
   ``1`` adds

   .. math::

      C_l(E) = 2\eta k \int_a^\infty dr\, \frac{O_l(kr)}{O_l(ka)}\, j_l(pr)

   to :math:`M_l` (Tribble 2014 eq. 2.79; Mukhamedzhanov 2017 eqs. 26–27;
   Typel & Baur eq. A.4, which reduces to the same bracket). For
   :math:`E < 0`, :math:`O_l` is replaced by the Whittaker function
   :math:`W_{-\eta, l+1/2}(2\kappa r)`. Without it (default ``0``) the vertex
   is the surface term alone, which neglects the Coulomb interaction of
   :math:`x` and :math:`A` outside the channel radius; the term matters for
   charged entrance pairs at energies well below the barrier. It is skipped
   for :math:`Z_1 Z_2 = 0` and costs about a factor 3–4 in run time.
   ``tests/7Li_p_a``: 2111.79 (2145.79 with ``vertex=perlevel``).
   Implementation and numerics:
   ``ThmCoulombTerm`` (``src/ThmFunc.cpp``), checked against
   ``tests/reference/thm_coulomb_term_reference.py``.

``spectatorEnergy`` — spectator motion in the half-off-shell momentum
   The quasi-free picture puts the spectator :math:`s` at rest, so
   :math:`p^2/2\mu_{xA} = E + B_{xs}`. A non-zero value :math:`T_s` (MeV,
   :math:`\ge 0`) replaces this by :math:`E + B_{xs} + T_s` with
   :math:`T_s = \langle p_{sx}^2 \rangle / 2\mu_{sx}` (Typel & Baur eq. 11).
   ``spectatorEnergy=`` applies to every THM entrance pair;
   ``spectatorEnergy[k]=`` to the entrance pair with key ``k`` (the pair
   number that the segment's entrance key refers to) and overrides the global
   value for it. Default 0. ``tests/7Li_p_a``: 2195.06 at 0.5 MeV (2088.25 with
   ``vertex=perlevel``).

Which ``kinematics=`` to use
   Read the data paper's definition of the extracted quantity. Division by the
   full three-body KF → ``kf3body``; by :math:`\lambda_3/\lambda_2` →
   ``lambda32``; triple cross section over :math:`|\phi|^2` only → ``triple``;
   data analysed with La Cognata's modified R-matrix formula, or the
   convention unknown → ``lacognata``. A THM segment always has a free
   arbitrary normalization, so only the energy dependence of :math:`K(E)`
   matters.
