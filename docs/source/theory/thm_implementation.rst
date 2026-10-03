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
   weight[1]=R_E.dat        # w(E) multiplying the model of <segmentsData> line 1
   weightTest[2]=R_E.dat    # the same for <segmentsTest> line 2
   experiment[E1] segments=1,2 background=linear   # see "THM experiments" below
   experiment[E2] segments=3 beam=14N target=12C spectator=d Ebeam=30 distortion=coulomb
   experiment[E3] segments=4 theta=50-70   # dsigma/dOmega averaged over theta_cm = 50-70 deg
   experiment[E4] segments=5 cbackground=1/2+:2=0.1,-0.05   # interfering THM-only background
   </thm>

An unknown key, an unknown value, a negative spectator energy, a missing
``</thm>`` or a bad weight (below) prints ``ERROR: <thm> ...``, and AZURE2
exits with a non-zero status before any calculation; ``pyazr`` refuses the
project the same way. ``pyazr.AzrModel`` reads and edits the block with the
GUI's rules (below): ``thm_options()`` returns the options as a dict,
``set_thm_option(key, value)`` and ``clear_thm_option(key)`` change one,
``set_thm_weight(segment, path, test=False)`` / ``clear_thm_weight`` a weight
table (the segment must exist and be THM; the table is read with the engine's
rules, a relative path from the directory of the .azr). A value the engine
would refuse raises ``ValueError`` and leaves the model unchanged; a project
whose block the engine refuses raises on any of them. Without an edit the
block is carried through a save unchanged.

The GUI edits the block under *Configure > THM Workspace...*, page *Model*
(:doc:`../user_guide/configure_menu`, "THM Workspace"): the four global
choices, the spectator energy, a table of per-pair spectator energies and a
table of per-segment weight files (with a file chooser; a file inside the
project directory is stored relative to it). Every control has a one-line
tooltip. The dialog applies the engine's rules (same keys and values,
spectator energies :math:`\ge 0`, every weight file read with the engine's
own reader) before it accepts; whether a weighted segment exists, is THM and
has its points inside the table is checked by the engine at startup. The GUI
writes the block after ``</targetInt>``. A project whose options are not
changed in the dialog keeps its block byte for byte
(``tests/gui/thm_block_test``). When they are changed, only non-default keys
are written: comment and blank lines of the block stay where they were, a
line whose value did not change is kept verbatim, a changed one is rewritten
keeping its indentation and inline comment, a key that returns to its default
is removed, and new keys are appended. Options that are all default remove
the block, comments included, since an empty block and no block are the same
to the engine (``tests/gui/thm_options_dialog_test``). A block with a line the
engine would refuse is not opened in the editor (it is kept as it is).

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
   coherent sum is not the fixed-angle limit either: along
   :math:`\hat p_{xA}` the partial waves add with Clebsch-Gordan weights
   :math:`\sqrt{2l+1}\,\langle s\,\nu\,l\,0|J\,\nu\rangle` and the
   :math:`J^\pi` groups interfere ("Fixed-angle observable" below, which
   computes that observable with ``theta=``). ``coherent`` is kept for
   comparison with older fits, and cannot be combined with a ``theta``
   window.
   ``tests/7Li_p_a`` (two entrance :math:`l` in one channel spin): 2138.48
   incoherent, 3196.77 coherent, same parameters (2181.45 / 2741.85 with
   ``vertex=perlevel``; 2180.69 / 2740.48 were the pins before 2026-09-25).
   The ``tests/7Li_p_a`` values in this section are those after the
   shift-function fix b6cc41b, which moved them by up to 3.2e-4, and the
   dS/dE step fix of 2026-09-28 (up to 1.9e-4), unless marked otherwise.

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
   reproduces the old values. Since b6cc41b and the dS/dE step fix these
   are 2181.45 → 2138.48 and 682.147. ``tests/7Li_p_a`` with
   ``vertex=onshell``: 2070.75.

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
   2138.48 / 2381.01 / 2625.33 / 5499.08 for the four values in table order
   (2181.41 / 2403.67 / 2633.59 / 5573.51 with ``vertex=perlevel``, before
   the dS/dE step fix).

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
   ``tests/7Li_p_a``: 2112.47 (2146.48 with ``vertex=perlevel``).
   ``coulombIntegral=1`` is refused together with an experiment whose
   :math:`a + A` wave is distorted (``distortion=coulomb``, ``optical``
   without ``opticalAA=plane``, or ``vertexModel=dw``), which contains the
   same :math:`x`-:math:`A` Coulomb force ("Coulomb effects: what each option
   contains" below). Implementation and numerics:
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
   value for it. Default 0. ``tests/7Li_p_a``: 2195.75 at 0.5 MeV (2088.91 with
   ``vertex=perlevel``). A single :math:`T_s` stands for the whole accepted
   range of spectator momenta; ``ps=`` on a THM experiment line averages over
   that range instead ("Spectator-momentum window" below). The two exclude
   each other for one entrance pair.

``weight`` — energy-dependent weight of one segment's model
   ``weight[k]=<file>`` multiplies the THM model cross section of segment
   ``k`` by :math:`w(E)` at every point and at every sub-point of the
   experimental-effect (resolution) folding, i.e. before the folding, with
   :math:`E` the c.m. energy of the THM entrance pair at that (sub-)point.
   ``k`` is the segment's line number in ``<segmentsData>``, counting inactive
   lines — the numbering of ``segment_<k>_norm``, of the first column of
   ``normalizations.out`` and of ``chiSquared.out``, and of the output files
   (``ESegment::GetSegmentKey``). ``weightTest[k]=<file>`` does the same for
   line ``k`` of ``<segmentsTest>`` (extrapolation of a THM segment with the
   weighted model; ``<segmentsTest>`` lines are numbered on their own, from 1).
   A relative path is taken from the directory of the ``.azr`` file; a path
   cannot contain ``#`` (it starts a comment).

   The file has two columns, :math:`E` (MeV) and :math:`w(E) > 0`, one row per
   line, at least two rows, strictly increasing in :math:`E`; ``#`` starts a
   comment. Between rows :math:`w` is interpolated log-linearly (linear in
   :math:`E`, linear in :math:`\ln w`: exponential pieces, which suit factors
   that vary by orders of magnitude). Every point of the segment must lie
   inside the table, or AZURE2 stops at startup; folding sub-points (and points
   moved by a fitted energy shift) beyond either end get the end value, and the
   first such evaluation prints one ``WARNING``.

   Errors (``ERROR: <thm> ...``, exit non-zero): a file that cannot be read,
   a line that is not two numbers, :math:`w \le 0`, energies not strictly
   increasing, fewer than two rows, a ``k`` beyond the last line of the block,
   a weight on a segment that is not THM (``isDiff < 10``), a point outside
   the table. A weight on an inactive (or unusable) line is ignored with a
   ``WARNING``.

   With a free (profiled) THM normalization a constant :math:`w` changes
   nothing but the norm: :math:`\chi^2` is identical, and since the norm
   multiplies the *data* (:math:`n^* = S_{mm}/S_{md}`, ``ESegment::
   ProfileNormChiSquared``), ``w ≡ 2`` *doubles* the norm written to
   ``chiSquared.out`` and ``normalizations.out``. Only the energy dependence
   of :math:`w` matters, so the choice of :math:`E_\mathrm{norm}` below is
   immaterial. ``tests/thm_options/check.sh`` pins both and an energy ramp.

   *Physics: the Coulomb-distortion factor.* The plane-wave (PWA) THM analysis
   takes the transfer amplitude of :math:`a + A \to s + F^*` as constant over
   the measured range (Mukhamedzhanov et al., J. Phys. G 35 (2008) 014016:
   :math:`M^\mathrm{DW}` "practically constant on the interval of a few hundred
   keV"), so the HOES excitation function is the measured yield divided by
   :math:`\mathrm{KF}\,|\phi_a(p_{sx})|^2` alone. When the projectile or the
   spectator moves below a Coulomb barrier this fails: the distorted-wave
   amplitude :math:`M(E)` varies with :math:`E` through the initial and the
   final (spectator–residual) Coulomb interaction. Mukhamedzhanov & Pang,
   PRC 99 (2019) 064618, eqs. 20–24 (DWBA, FRESCO), and Mukhamedzhanov,
   arXiv:2609.04498, eqs. 22–30 (zero-range prior amplitude, Nordsieck
   integral), find for 12C(14N,α 20Ne)d and 12C(16O,α 20Ne)α a factor that
   varies by 10\ :sup:`2`–10\ :sup:`5` across the measured range — mostly the
   exit-channel Coulomb penetrability of the spectator. The PWA-extracted
   quantity is then

   .. math::

      \sigma^\mathrm{PWA}(E) = R(E)\, \sigma^\mathrm{HOES}(E), \qquad
      R(E) = \frac{|M(E)|^2}{|M(E_\mathrm{norm})|^2}
      \quad \text{(arXiv:2609.04498 eqs. 29–30)},

   so the published :math:`S^*(E)` has to be *divided* by :math:`R(E)` — or,
   equivalently and without touching the data, the HOES model that is fitted
   to them *multiplied* by :math:`R(E)`. That is what ``weight`` does, with
   :math:`w(E) = R(E)`. The definition that avoids every convention is
   :math:`w(E) = S^*_\mathrm{PWA}(E) / S^*_\mathrm{corrected}(E)` up to a
   constant: the published curves do not all define "R" the same way round
   (1806.08828 eq. 24 writes a ratio of DWBA cross sections at
   :math:`E_\mathrm{norm}` over :math:`E`, and its figures show the factor that
   multiplies :math:`S^*`). For 12C+12C the correction *lowers* the
   low-energy :math:`S^*`, so :math:`w` grows toward low :math:`E`, by about
   10\ :sup:`2` from 2.5 to 0.8 MeV. For light, weakly charged systems
   (7Li(p,α), 6Li(d,α), 17O(n,α)) PWA and DWBA agree to about 10 % and no
   weight is needed. The same hook can carry any other correction that
   multiplies the HOES cross section; Mukhamedzhanov's line-shape factor
   :math:`N_C` for charged spectators is built in, per level and inside the
   coherent sum (``lineshape=on``, "Coulomb line shape" below), which a weight
   of the summed cross section cannot do.

   *Where the table comes from.* Since September 2026 AZURE2 computes
   :math:`R(E)` itself in the zero-range DWBA of Mukhamedzhanov's papers,
   with point-Coulomb or optical-potential distortions, for the segments of
   a THM experiment that has the kinematics (``distortion=coulomb|optical``,
   "Distortion factor R(E)" below); ``distortion=table:<file>`` applies a
   table of this format to all segments of an experiment. A ``weight[k]``
   table remains the hook for anything else (a finite-range DWBA or CDCC,
   e.g. FRESCO, with the experiment's angular acceptance), made outside, for
   instance

   * from the published curves (1806.08828 Fig. 10; 2609.04498 Figs. 7–9),
     digitized, taking care of which way round they are defined;
   * from a DWBA of the transfer reaction at each :math:`E`:
     :math:`w(E) = d\sigma^\mathrm{DW}(E) / d\sigma^\mathrm{PWA}(E)` up to a
     constant;
   * as an order-of-magnitude estimate, from the Coulomb penetrability of the
     spectator in the exit channel: it leaves the residual nucleus with
     :math:`E_{sF} = Q - E` (12C(14N,α 20Ne)d:
     :math:`E_{d+{}^{24}\mathrm{Mg}} = 3.573\,\mathrm{MeV} - E`), and
     :math:`w(E) \propto P_0(E_{sF})` at :math:`r \approx 5` fm. For 12C+12C
     this reproduces the 2026 curve to within a factor 2 (the factor that
     multiplies :math:`S^*`, :math:`1/w` in the scale of their curve, is
     4.4×10\ :sup:`-3` at 0.8 MeV against their 3.7–4.2×10\ :sup:`-3`).

Which ``kinematics=`` to use
   Read the data paper's definition of the extracted quantity. Division by the
   full three-body KF → ``kf3body``; by :math:`\lambda_3/\lambda_2` →
   ``lambda32``; triple cross section over :math:`|\phi|^2` only → ``triple``;
   data analysed with La Cognata's modified R-matrix formula, or the
   convention unknown → ``lacognata``. A THM segment always has a free
   arbitrary normalization, so only the energy dependence of :math:`K(E)`
   matters.

THM experiments
---------------

Segments measured in one experiment -- several exit channels, angular bins or
runs of the same three-body reaction -- carry the *same* arbitrary THM scale.
A line

.. code-block:: text

   experiment[<name>] key=value key=value ...

in the ``<thm>`` block groups them. Keys may be spread over several lines of
the same name (they merge; a key given twice is an error). ``<name>`` is
letters, digits and ``_ - . +``.

``segments=`` (required)
   ``<segmentsData>`` line numbers, counting inactive lines (the numbering of
   ``weight[k]``, of ``chiSquared.out`` and of ``normalizations.out``), as a
   comma list with ranges: ``1,2,5-7``. Each must be a THM segment with a free
   norm, and in no other experiment. An inactive line is left out with a
   ``WARNING``.
``background=none|const|linear|quadratic``
   A smooth background added to the model (default ``none``), below.
``beam=``, ``target=``, ``spectator=``, ``Ebeam=``
   The three-body reaction: nuclides from the built-in table (n p d t 3He 4He
   6Li 7Li 9Be 10B 11B 12C 13C 14N 15N 16O 17O 18O 19F 20Ne 23Na 24Mg; AME2020
   atomic masses minus the electrons plus their binding, i.e. nuclear masses in
   u) or ``Z,A,mass`` (nuclear mass in u), and the lab beam energy in MeV. All
   four or none. Only ``lineshape=on`` and a ``ps`` window use them: AZURE2 checks that one of beam and
   target is a nucleus of the segments' entrance pair and the other is the
   second nucleus plus the spectator (the Trojan horse :math:`a = x + s`), and
   prints :math:`B_{xs}` and the quasi-free energy
   :math:`E_{xA} - B_{xs}` (spectator at rest in the lab when the Trojan horse
   is the target, moving with the beam velocity when it is the beam).
``lineshape=on|off``
   The Coulomb line-shape factor of the spectator (default ``off``; needs the
   four kinematics keys), below.
``ps=delta|hulthen:...|gauss:...|table:<file>``, ``psNodes=``
   The spectator-momentum window: the HOES cross section is averaged over
   the spectator directions whose :math:`|p_s|` lies in it, with the weight
   :math:`|\phi(p_s)|^2\,d\cos\theta_\mathrm{cm}` (default ``delta``, the
   quasi-free point; the other forms need the four kinematics keys), and its
   number of Gauss-Legendre nodes in :math:`\cos\theta_\mathrm{cm}` (1-64,
   default 16), below.
``distortion=none|coulomb|optical|table:<file>``, ``opticalAA=``, ``opticalSF=``, ``spectatorAngle=``, ``distortionRef=``, ``distortionRatio=``, ``boundState=``
   The distortion factor :math:`R(E)` multiplying the model of every
   segment (default ``none``; ``coulomb`` and ``optical`` need the four
   kinematics keys), below ("Distortion factor R(E)"). With a ``ps`` window
   ``distortionRatio=dw`` is refused ("Coulomb effects" below).
``theta=all|<thmin>-<thmax>``
   The angular window of the exit pair (degrees): the model of every segment
   is the HOES :math:`d\sigma/d\Omega` averaged over it instead of the
   angle-integrated cross section (default ``all``), below ("Fixed-angle
   observable").
``cbackground=<term>[;<term>...]``
   A coherent (interfering) THM-only background amplitude per J\ :sup:`π`,
   entrance and exit channel, whose real and imaginary parts are fit
   parameters (default: none), below ("Coherent background").
``vertexModel=pw|dw``
   The entrance vertex of every segment: the plane-wave :math:`M_l`
   (default ``pw``) or the surface term of the prior-form DWBA built from the
   experiment's distorted waves (``dw``; needs ``distortion=coulomb|optical``,
   and then replaces :math:`R(E)`), below ("Distorted-wave entrance vertex").

Anything else -- an unknown key or nuclide, a malformed value, partial
kinematics, a segment that does not exist, is not THM, has a fixed norm or is
in two experiments, too few points for the linear parameters -- prints
``ERROR: <thm> experiment[<name>]: ...`` and AZURE2 exits with a non-zero
status. A file without experiment lines, and an experiment of one segment
without background, give results identical to before, bit for bit
(``tests/thm_experiment``).

*Shared norm.* All segments of an experiment share one profiled norm: the
sums :math:`S_{mm}, S_{md}, S_{dd}` of the previous section run over the points
of all of them, :math:`n^* = S_{mm}/S_{md}` and
:math:`\chi^2 = S_{dd} - S_{md}^2/S_{mm}`. Every segment of the experiment
carries :math:`n^*` in ``chiSquared.out`` and ``normalizations.out``, and its
own points' share of the :math:`\chi^2`. Segments in no experiment keep their
own profile.

*Background.* With ``background=const|linear|quadratic`` the curve compared
with the data is

.. math::

   f_i = n^{-1}\,\bigl(m_i + b_0 + b_1 E_i + b_2 E_i^2\bigr),

:math:`m_i` the THM model at point :math:`i` (resolution-folded, weighted),
:math:`E_i` the point's c.m. energy of the THM entrance pair (MeV), and
:math:`b_k` in the units of the model (the data scaled by :math:`n`, the
scale of the output files). The background is added to the *folded* model: it
is smooth on the scale of the resolution, so folding it would change nothing
but the cost; for the same reason it is not multiplied by ``weight[k]``.
:math:`s = 1/n` and :math:`a_k = s\,b_k` are linear parameters and are
eliminated together by weighted linear least squares over all points of the
experiment: with the design matrix :math:`A = [m, 1, E, E^2]` (as many
background columns as terms), :math:`\tilde A = A/e`, :math:`\tilde y = d/e`,

.. math::

   c = (s, a_0, \dots) = G^{-1}\tilde A^T\tilde y, \qquad
   G = \tilde A^T\tilde A, \qquad
   \operatorname{cov}(c) = G^{-1},

solved by Cholesky on the Jacobi-scaled :math:`G` (``SolveThmProfile`` in
``src/ThmExperiment.cpp``). Points with :math:`e_i = 0` carry no weight. The
output file writes :math:`m_i + b(E_i)` as the fitted curve next to the data
scaled by :math:`n^*`. The background may come out negative: it is a
phenomenological term (non-quasi-free or sequential contributions, an
imperfect subtraction), not a positive-definite cross section, and a sign
constraint would make the profile non-linear; judge a negative background on
the plot. :math:`n^* = 1/s` must be positive.

*Degenerate cases.* A singular :math:`G` (the model a combination of the
background terms, or collinear background columns), or :math:`s \le 0` (no
positive overlap of model and data), leave the norm at 1, as the single-segment
profile does, and profile the background alone; if its own normal matrix is
singular too, nothing is profiled. ``thm_experiments.out`` says which in its
``status`` line. Fewer points (with an error) than linear parameters plus one
is refused at startup.

*Output.* ``output/thm_experiments.out`` lists, per experiment, its segments,
background, points, :math:`\chi^2`, status, the norm and the :math:`b_k` with
their uncertainties and covariance, transformed from :math:`\operatorname{cov}(c)`
by :math:`n = 1/s`, :math:`b_k = a_k/s`. They are the uncertainties of the
profile at fixed R-matrix parameters, not scaled by :math:`\chi^2/\nu`.
``pyazr``: ``session.thm_background(name)`` (norm, ``b``, ``cov``, ``sigma``,
``chi2``, ``status``) and ``session.thm_experiments()``; ``segment_norms()``
gives the shared :math:`n^*` for each segment. ``AzrModel.thm_experiments()``,
``set_thm_experiment(name, segments, background=..., beam=..., target=...,
spectator=..., Ebeam=...)`` (replaces the record with one line) and
``clear_thm_experiment(name)`` edit the lines with the engine's rules. In the
GUI, *Configure > THM Workspace...*, page *Experiments*, edits them
(:doc:`../user_guide/configure_menu`, "THM Workspace"), including
``lineshape``; keys it does not show, such as ``ps`` and ``psNodes``, are kept
as written.

*Derivatives.* Where the model Jacobian :math:`J_m = \partial m/\partial p`
is used (MIGRAD's THM gradient, ``pyazr``'s ``residual_jacobian`` and
``chi2_and_grad``, the covariance band), the dependence of the profiled
:math:`c^*(p)` is included exactly. With :math:`\rho = \tilde y - \tilde A c`
and :math:`\dot{\tilde A} = [J_m/e, 0, \dots]` (only the model column
depends on :math:`p`), differentiating the normal equations gives
(Golub & Pereyra, SIAM J. Numer. Anal. 10 (1973) 413)

.. math::

   \dot c = G^{-1}\bigl(\dot{\tilde A}^T\rho - \tilde A^T\dot{\tilde A}\,c\bigr),
   \qquad
   \dot r = \dot{\tilde A}\,c + \tilde A\,\dot c
          = P_\perp \dot{\tilde A}\,c + \tilde A G^{-1}\dot{\tilde A}^T\rho,

:math:`P_\perp = 1 - \tilde A G^{-1}\tilde A^T` the projector onto the
orthogonal complement of the columns of :math:`\tilde A`. The first term alone
(Kaufman's variable-projection approximation) is exact only at zero residual;
the second is kept. For the model column alone this is the :math:`\partial
s/\partial p` below. ``ThmProfileDerivative``; checked against central
differences in ``tests/pyazr/thm_experiment_test.py`` (all six parameters of
``tests/18O_p_a_thm``, linear background, rel. :math:`\le 5\times 10^{-6}`)
and, for the band, in ``tests/thm_band/check.sh``. The band row of a point is
the derivative of the curve the output shows, :math:`(s J_m + A\dot c)/s`.

Coherent background
-------------------

``background=`` adds a smooth term to the *cross section*: whatever it stands
for, it does not interfere with the resonances. Two kinds of non-resonant
physics can interfere, and AZURE2 treats them differently.

*Non-resonant x + A → b + B in the same partial wave.* Direct coupling of the
entrance and exit channels and the tails of distant levels belong to the
two-body reaction itself. They are in the direct data as well as in the HOES
amplitude, with the same reduced widths, and the R matrix describes them with
a **background pole**: a broad level of the same J\ :sup:`π` far above the
data, with widths in the entrance and exit channels. In the HOES amplitude
:math:`\sum_{\lambda\lambda'}\gamma_{\lambda f}A_{\lambda\lambda'}V_{\lambda'}`
such a level enters through the level matrix exactly as the resonances do, so
it acts on THM and direct segments consistently and keeps the S matrix unitary
(its phase relative to the resonances is fixed by the widths, not free). It
needs no option: add the level on the Levels tab (or ``<levels>``) and free
its widths. This is the term Mukhamedzhanov et al. drop from the THM amplitude
when they "neglect the direct coupling between the initial x + A and final
b + B channels" (J. Phys. G 35 (2008) 014016, p. 2) and the non-resonant
:math:`V^N_{bB}\,G\,V_{sx}` term of Mukhamedzhanov, Kadyrov & Pang, EPJA 56
(2020) 233, eqs. (48)–(49), which "does not produce the resonance peak ... and
can be treated as a background".

*Non-quasi-free mechanisms.* Sequential decay through other two-body
resonances, direct breakup and other three-body processes feed the same final
state :math:`c + B + s` at the same kinematics but do not go through the
x + A system. They are absent from direct data. The THM analyses remove them by
kinematic cuts (spectator momentum :math:`p_s \lesssim 30`–50 MeV/c, the loci
of the relative energies, angular correlations) and fit what is left with an
incoherent polynomial (Tribble et al., RPP 77 (2014) 106901, §4.2–4.3 and
Fig. 31: three Gaussians plus a first-order polynomial for 18O(p,α); La
Cognata et al. 2010: a linear background subtracted; Spitaleri et al.,
PRC 95 (2017) 035801, eq. (11): an incoherent sum with a non-resonant
polynomial). None of these analyses models an interference between a non-QF
and the QF amplitude, and none estimates it; ``background=`` is that
incoherent term. A residue of a non-QF amplitude in the same final state can,
however, interfere: projected on the partial waves of the b + B relative
motion at the quasi-free kinematics it has a component in every J\ :sup:`π`
and exit channel. ``cbackground=`` models that component.

``cbackground=<term>[;<term>...]`` on an experiment line, a term being

.. code-block:: text

   <J><+|->:<exit pair key>[:<s>,<l>,<s'>,<l'>][:const|:linear][=<Re c0>,<Im c0>[,<Re c1>,<Im c1>]]

adds, for each combination of an entrance channel :math:`(s,l)` and an exit
channel :math:`(s',l')` of the J\ :sup:`π` group (all of them, or the one
named), a complex amplitude to the HOES amplitude of that combination before
it is squared:

.. math::

   x^J_{(s'l'),(sl)} \to \sqrt{K\,2P_{c'}}\;e^{i(\omega_{c'}-\phi_{c'})}
   \Bigl[\sum_{\lambda\lambda'}\gamma_{\lambda c'}N_\lambda A_{\lambda\lambda'}V^{(s,l)}_{\lambda'}
   + c(E)\,M_l(p, B_c)\Bigr],
   \qquad c(E) = c_0 + c_1 E,

:math:`M_l` the vertex of the entrance channel without a width (the boundary
:math:`B_c` of ``vertex``: :math:`S_c` at the lowest level for ``constant``,
:math:`L_c(E)` for ``onshell``, the channel boundary condition for
``perlevel``; the node of a spectator window; the two components of the DW
vertex), :math:`E` the c.m. energy of the entrance pair in MeV. :math:`c` takes
the place of :math:`\sum\gamma_{\lambda c'}A_{\lambda\lambda'}\gamma_{\lambda' c}`:
dimensionless, it is the THM-only analogue of a level-matrix element, and it
carries the threshold behaviour of the exit channel (:math:`\sqrt{P_{c'}}`) and
the energy dependence of the vertex as a background pole would. It has no line
shape (:math:`N_\lambda` belongs to a level). For the angle-integrated
observable the result is
:math:`(2J+1)\,K\,2P_{c'}\,|M_l|^2\,|\sum\gamma A\gamma + c|^2` for one level
and one channel; at a fixed angle (``theta=``) the background interferes across
J\ :sup:`π` groups like the resonant amplitudes. One complex :math:`c` per
combination: without the channels a term gives every combination its own
amplitude (all starting at the given values); ``const`` (default) has
:math:`c_1 = 0`. Values default to 0; a value followed by ``f`` is fixed.

*Fit parameters.* :math:`c` enters the model quadratically, so unlike the
incoherent background it cannot be profiled linearly: Re and Im of
:math:`c_0` (and :math:`c_1`) are ordinary fit parameters, the last block of
the parameter vector after the energy shifts (none without the key), named
``cbkg_<experiment>_<J><π>_<exit>_<s>,<l>,<s'>,<l'>_re0`` (``im0``, ``re1``,
``im1``). They are in ``param.par`` / ``param.sav`` (an external parameter file
sets them by name), MIGRAD fits them with the THM part of its gradient (central
differences of the THM χ² only), the CLI writes their values to
``output/thm_experiments.out`` (``cbkg`` lines), and ``<parameterSettings>`` can
limit them by name. Their derivatives in ``ComputeTHMRows`` are central
differences of the HOES model, exact up to round-off since the model is
quadratic in them; covariance bands keep to the R-matrix parameters, as they do
for norms. ``pyazr``: kind ``"cbkg"`` (``parameters.cbkg``, with the J group,
entrance channel and exit pair), ``residual_jacobian`` and ``chi2_and_grad``
carry their columns, ``thm_experiments()[name]["cbkg"]`` their values;
``save_fit`` writes the fitted values back into ``cbackground=`` (one explicit
term per combination); ``AzrModel.set_thm_experiment(..., cbackground=...)``
and ``set_thm_cbackground(name, value)`` edit the key with the engine's rules.
The GUI: *THM Workspace*, page *Experiments*, section *Coherent background*;
the free parameters appear on the Fitting tab (tab *THM Background*), and
loading a ``.sav`` there writes their values into the key.

*Refused* (``ERROR: <thm> experiment[<name>]: cbackground ...``): a malformed
term, a J\ :sup:`π` the model has no level of, an exit pair that no segment of
the experiment has, a J\ :sup:`π` group that does not couple the entrance pair
to that exit pair or has no such channels, a combination given twice,
segments with different entrance pairs, and ``entranceL=coherent`` (its
buckets merge the l of a channel spin). Without the key nothing changes, and a
background fixed at zero gives the files of none.

*Which one to use.* For non-resonant physics of the x + A system --
anything that is also in direct data -- use a background pole: it constrains
THM and direct data together and preserves unitarity. Use ``cbackground`` only
for what direct data cannot contain, as a test of whether a THM-direct tension
can be a non-QF interference. Read a good fit with care: :math:`c` is
strongly correlated with the resonance parameters (a constant :math:`c`
shifts and skews a peak much as a nearby level does), the quadratic form has
in general two solutions of equal χ² that differ in the phase of :math:`c`,
and the incoherent ``background=`` and the coherent term compete for the same
smooth part of the spectrum. A background that the data require only through
the interference is evidence of a missing amplitude, not of its origin.

*Validation.* ``tests/thm_coherent_background/check.sh``: (a, b) without the
key, and with a background fixed at zero (angle-integrated, ``theta`` window,
``ps`` window with a linear incoherent background, DW vertex), the output is
byte-identical; (c) a toy model with one 1/2\ :sup:`+` level and neutral
:math:`l = 0` channels against the closed form
:math:`2\,(k_f/\mu_f)\,2P_f\,|M_0(\gamma_f\gamma_c A + c(E))|^2`, constant and
linear :math:`c`, to :math:`4\times10^{-9}` (the formula's own difference
quotient), and :math:`m(c)/m(0) = |\gamma_f\gamma_cA + c|^2/|\gamma_f\gamma_cA|^2`
to :math:`1.4\times10^{-10}`; a 0–180° window times :math:`4\pi` equals the
angle-integrated model with the background to :math:`6\times10^{-11}`; (d)
synthetic 18O(p,α) data made with :math:`c_0 = 0.3 - 0.2i` are fitted by
MIGRAD from :math:`c_0 = 0` back to :math:`0.29999897 - 0.19999902i`; (f) the
refusals. ``tests/pyazr/thm_coherent_background_test.py``: ``AzrModel``, CLI ==
session, the Jacobian of all ten free parameters (six R-matrix, four ``cbkg``
of a linear term) against central differences (rel. about :math:`10^{-9}` for
the ``cbkg`` columns; the test asks for :math:`10^{-4}`), ``save_fit`` round trip, a fixed zero background bit
for bit equal to none.

*Effect study* (October 2026; joint least-squares fits with ``pyazr``'s
``residual_jacobian``, the fitter in the user's script; not part of the
tests).

18O(p,α)15N, ``examples/o18_lacognata2010`` (THM 39 points, La Cognata
2010, with the linear S-factor background of an earlier joint fit,
:math:`-5152 + 7002\,E` MeV b, subtracted -- not the authors' dashed line
of their Fig. 4, which lies above it below 0.84 MeV; Mak 1978, 32; Amsel 1967, 90;
the two 1/2\ :sup:`+` levels' energies and four widths and the two direct
norms free). χ² per data set:

======================================================  =========  ========  =========
THM model                                               THM (39)   Mak (32)  Amsel (90)
======================================================  =========  ========  =========
no background                                           410.2      52.8      119.3
``background=linear`` (incoherent)                      410.9      54.2      117.2
``cbackground=1/2+:2`` (const)                          324.1      50.8      125.2
``cbackground=1/2+:2:linear``                           306.0      55.3      117.6
``background=linear cbackground=1/2+:2``                291.9      46.0      133.7
THM alone, no background (no direct data)               198.1      --        --
THM alone, ``cbackground=1/2+:2``                       192.2      --        --
======================================================  =========  ========  =========

The constant term converges to :math:`c_0 = -0.019 + 0.069i` from six
starting points (0, ±0.3, ±0.3i, 1 + i), one minimum. Of the
212 units of χ² that the direct data cost the THM shape (410 jointly
against 198 alone) the coherent background recovers 86 (constant), 104
(linear) and 118 (with the incoherent term too, at +8 on the direct data);
the THM alone barely needs it (198 → 192). It does not resolve the tension:
:math:`\chi^2/N` of the THM stays 7.5–8.3 with the direct data, and 5.1
without them, i.e. the published points scatter beyond their errors about
any smooth model.

19F(p,α\ :sub:`2`)16O, the full 53-point THM window of Su et al. 2025
(``examples/f19_pag_thm`` model, stage J2 parameters -- 22 R-matrix
parameters free, JUNA and Spyrou with 10 % norm priors, the direct ωγ and
Γ as penalty rows, the −9.17 keV THM shift; start ``J2_s9``):

==================================  =========  ==========  ===========  =========  ==============  ==============
model                               THM (53)   JUNA (20)   Spyrou (9)   penalty    ωγ(828) (eV)    ωγ(564) (eV)
==================================  =========  ==========  ===========  =========  ==============  ==============
no coherent term                    203.3      13.8        3.6          343.0      183             9.6
``cbackground=2-:6``                145.8      14.1        3.2          240.3      315             3.0
``cbackground=1-:6;2-:6``           147.2      14.1        3.2          238.5      318             3.0
==================================  =========  ==========  ===========  =========  ==============  ==============

(direct: ωγ(828) = 775(35) eV, ωγ(564) = 48(7) eV, ωγ(790) = 17(5) eV, which
stays at about 0.) Of the eight 2\ :sup:`-` amplitudes only
:math:`(s,l) = (1,1) \to (s',l') = (3,2)` moves, to
:math:`c_0 = -0.27 - 0.82i`; the 1\ :sup:`-` ones stay near zero
(:math:`|c_0| \approx 0.03`). The interference lowers the THM χ² by 58 and
the penalty by 103, but the full window still does not fit with the direct
strengths: ωγ(828) rises from 183 to 315 eV, not to 775 eV, and ωγ(564)
falls further. The tension there is the l-dependence of the vertex
(``thm_19F/RESULTS.md``, section 4), which a smooth amplitude in one
partial wave cannot undo. In both cases the coherent background is a
useful diagnostic -- it says how much of a THM-direct mismatch a smooth
interfering amplitude can absorb -- and in neither does it make the THM and
the direct data consistent.

Coulomb line shape
------------------

In :math:`A + a(x+s) \to s + F^*(x+A) \to s + b + B` the spectator :math:`s`
leaves in the Coulomb field of the resonance :math:`F^*`, and after
:math:`F^*` has decayed in that of :math:`b` and :math:`B`. Mukhamedzhanov,
Kadyrov & Pang, EPJA 56 (2020) 233 (arXiv:2007.13331), eqs. (55)-(57), derive
the resonant THM amplitude with these three-body Coulomb interactions: the
Breit-Wigner factor :math:`(E_0 - E - i\Gamma/2)^{-1}` becomes the branch
point :math:`(E_0 - E - i\Gamma/2)^{-1-i\zeta}` times Coulomb factors that do
not depend on the level,

.. math::

   N_C = \frac{\Gamma(1+i\eta_{bs})\Gamma(1+i\eta_{Bs})}{\Gamma(1+i[\eta_{bs}+\eta_{Bs}])}
         F(-i\eta_{Bs}, -i\eta_{bs}, 1; -1)\,(-\gamma_{(0)})^{i\eta_{bs}}(-\nu)^{i\eta_{Bs}}
         \Bigl(E_0 - E - i\frac{\Gamma}{2}\Bigr)^{-i\zeta},
   \qquad \zeta = \eta_{sb} + \eta_{sB} - \eta_0

(eqs. 56-57; :math:`\eta_{ij} = Z_iZ_j\alpha\mu_{ij}/k_{ij}`, :math:`\eta_0`
the :math:`s`-:math:`F^*` parameter of the intermediate state), and for a
narrow level

.. math::

   |N_C|^2 = \frac{\sinh[\pi(\eta_{sb}+\eta_{sB})]}{\sinh(\pi\eta_{sb})\sinh(\pi\eta_{sB})}
             \frac{\pi\eta_{sb}\eta_{sB}}{\eta_{sb}+\eta_{sB}}
             \frac{\pi\eta_\zeta}{\sinh(\pi\eta_\zeta)}\,|F|^2\,
             \exp\Bigl[2\zeta\arctan\frac{2(E_0 - E)}{\Gamma}\Bigr]

(2020 eq. 62 = Mukhamedzhanov, EPJA 58 (2022) 71, eq. 32). Both papers then
set :math:`|N_C| = 1` in their fits. ``lineshape=on`` on an experiment line
applies it:

*zeta.* :math:`\eta_{sb}` depends on the direction of :math:`b`, which the
HOES cross section integrates over; the papers' own tractable limit (2020 p.
15 case 2, 2022 eq. 34) is :math:`m_B \gg m_s, m_b`, so that
:math:`k_{sB} \approx k_{sF}`, with :math:`|\eta_{sb}| \ll 1` dropped. AZURE2
takes that limit:

.. math::

   \zeta(E) = \eta_{sB} - \eta_0
            = \frac{Z_s\,\alpha\,(Z_B\,\mu_{sB} - Z_F\,\mu_{sF})}{\hbar c\,k_{sF}(E)},
   \qquad
   k_{sF} = \frac{\sqrt{2\mu_{sF}E_{sF}}}{\hbar c}, \quad
   E_{sF} = E_{aA} - B_{xs} - E,

:math:`b` the lighter and :math:`B` the heavier nucleus of the segment's exit
pair, :math:`F = x + A` (:math:`Z_F = Z_b + Z_B`, :math:`m_F = m_x + m_A` from
the entrance pair), :math:`E_{aA}` the beam-target c.m. energy (non-
relativistic, from ``Ebeam`` and the table masses), :math:`B_{xs}` from the
masses, :math:`E` the :math:`x + A` c.m. energy from its threshold (energy
conservation in the three-body final state; all energies of the exit channel
move with :math:`E`, so :math:`E_0 - E_{bB} = E_\lambda - E`). Since
:math:`Z_B < Z_F`, :math:`\zeta < 0` whenever :math:`\mu_{sB} \approx
\mu_{sF}`: after :math:`F^*` decays the spectator is repelled by :math:`Z_B`
only, gains less energy from the field than it would have, and the balance
goes to :math:`b + B` -- the peaks move *up* in :math:`E` (the
post-collision-interaction shift of atomic physics, ref. [18] of the 2020
paper). :math:`\eta_0` is evaluated at the point's :math:`k_{sF}(E)` rather
than at the pole momentum :math:`k_0`: the two agree at the pole (narrow
resonance), and it keeps :math:`\zeta` common to all levels at one energy, so
that the arbitrary energy unit of the complex power is a common phase.
Taking the :math:`s`-:math:`F` Sommerfeld parameter alone would be
:math:`\zeta = -\eta_0` (eq. 57 with :math:`\eta_{sb} = \eta_{sB} = 0`,
no final-state interaction); the final-state :math:`\eta_{sB}` cancels most
of it, leaving :math:`|\zeta| \approx (Z_b/Z_F)\,\eta_0` (for
12C(14N,d): 0.18 :math:`\eta_0` into :math:`\alpha + {}^{20}`\ Ne, 0.09
:math:`\eta_0` into :math:`p + {}^{23}`\ Na).

*Per level, coherently.* Level :math:`\lambda`'s exit amplitude -- the level
that decays to :math:`b + B` -- is multiplied by

.. math::

   N_{C,\lambda}(E) = e^{\pi\zeta/2}\,
     \Bigl(E_\lambda - E - i\frac{\Gamma_\lambda}{2}\Bigr)^{-i\zeta}
     \quad (\text{MeV}), \qquad
   |N_{C,\lambda}|^2 = \exp\Bigl[2\zeta\arctan\frac{2(E_\lambda - E)}{\Gamma_\lambda}\Bigr],

inside the coherent sum over levels,
:math:`\sum_\lambda N_{C,\lambda}\gamma_{\lambda c'}\sum_\mu A_{\lambda\mu}v_\mu`
(``THMMatrixFunc``). For an isolated level the HOES cross section is multiplied
by eq. (62)'s exponential, normalized to 1 at the pole; the level-independent
factors of eq. (56) and the prefactor of eq. (62) are smooth in :math:`E` and
are left to the arbitrary THM normalization. The pole is the level's observed
energy and total width at the current parameters, Brune (required:
``lineshape=on`` with the formal parameterization is refused):
:math:`\Gamma_\lambda = \sum_c 2\gamma_c^2P_c(E_\lambda) /
(1 + \sum_c\gamma_c^2\,dS_c/dE)` over the open particle channels plus the
radiative widths, as ``parameters.out`` writes them. The papers derive
:math:`N_C` for one isolated narrow level. For overlapping levels the poles of
the level matrix are not exactly the Brune :math:`(E_\lambda, \Gamma_\lambda)`,
and the factor on the exit index is an approximation; for a broad level
(:math:`\Gamma \gtrsim` the scale on which :math:`\zeta` varies) the
narrow-resonance form of eq. (62) is itself an approximation.

*Checks and output.* Every data point must leave :math:`E_{sF} > 0` (refused
otherwise; a sub-point of a folding grid beyond the limit takes
:math:`E_{sF} = 1` keV with a ``WARNING``). A neutral spectator gives
:math:`\zeta = 0` and results identical to ``lineshape=off`` bit for bit.
``thm_experiments.out`` adds, per line-shape experiment, :math:`E_{sF}`,
:math:`\eta_0` and, per exit pair, :math:`\zeta` and the size of the neglected
:math:`\eta_{sb}` (averaged over the :math:`b` direction,
:math:`\langle 1/|v_s - v_b|\rangle = 1/\max(v_s, v_b)` in the :math:`F^*`
frame) at the lowest and highest point energy. ``pyazr``:
``session.thm_lineshape(name, energies)`` returns :math:`E_{sF}`,
:math:`\eta_0` and, per exit pair, :math:`\zeta`, :math:`\eta_{sb}` and per
level the pole and :math:`|N_C|^2`; ``AzrModel.set_thm_experiment(...,
lineshape=True)``. Extrapolation segments (``<segmentsTest>``) are in no
experiment and carry no line shape. ``tests/thm_lineshape`` and
``tests/pyazr/thm_lineshape_test.py``.

*Size.* ``tests/18O_p_a_thm`` with a made-up 18O(3He,α15N)d at 115 MeV
(:math:`E_{sF} \approx 10.4` MeV): :math:`\zeta \approx -0.14`, the peak of the
0.61 MeV level moves up by 13 keV and the model changes by 0.79-1.44 across
it -- but the neglected :math:`\eta_{sb} \approx 0.13` is as large as
:math:`\zeta`. For 2H(18O,α15N)n the spectator is a neutron: :math:`N_C = 1`.
For 12C(14N,α20Ne / p23Na)d at 30 MeV (Tumino et al. 2018; ``examples/
c12c12_tumino2018`` at its parameters, the four THM segments one experiment,
as the example's own line ``experiment[E1] segments=1-4 beam=14N
target=12C spectator=d Ebeam=30`` has them: ``lineshape=on`` and
``distortion=`` can be switched on there, also in the GUI's THM workspace;
with both off the example gives :math:`\chi^2` 111.651, THM 61.02):
:math:`E_{sF}` = 2.76-0.88 MeV over :math:`E` = 0.82-2.69 MeV,
:math:`\eta_0` = 1.55-2.74, :math:`\zeta` = -0.28 to -0.49 into
:math:`\alpha_{0,1}` and -0.13 to -0.24 into :math:`p_{0,1}`; the folded
model changes by a factor 0.51-3.6 (α) and 0.73-2.0 (p) across the range,
rising with :math:`E` (most levels lie below a given energy, and above its
pole a level is raised by up to :math:`e^{\pi|\zeta|}`), the low-energy
peaks move up by 0-15 keV, and the THM :math:`\chi^2` at the published
parameters goes from 61 to 479 (shared norm profiled): a fit with
:math:`N_C` would move the levels. The neglected :math:`\eta_{sb}` is
0.26-0.29 into the α channels (half of :math:`|\zeta|`) and 0.07-0.10 into
the p channels: the case-2 limit is marginal for α.

Spectator-momentum window
-------------------------

*Kinematics.* In :math:`A + a \to s + F^*`, :math:`a = (s\,x)`, the
transferred particle :math:`x` is virtual, and energy-momentum conservation at
the vertices :math:`a \to s + x` and :math:`x + A \to F` gives the off-shell
:math:`x`-:math:`A` momentum

.. math::

   \frac{p_{xA}^2}{2\mu_{xA}} = E + B_{xs} + \frac{p_{s}^2}{2\mu_{sx}}

(Mukhamedzhanov et al., PRC 96 (2017) 024623, eqs. 29-32; Typel & Baur, Ann.
Phys. 305 (2003) 228, eq. 11), :math:`p_s = p_{sx}` the :math:`s`-:math:`x`
relative momentum in :math:`a`, which is the spectator momentum in the rest
frame of :math:`a` (2017 eq. 30). The quasi-free point is :math:`p_s = 0`.
Experiments accept events with :math:`|p_s|` up to some tens of MeV/c (for
the deuteron typically 30-40 MeV/c, below :math:`\kappa_{sx}\hbar c = 45.7`
MeV/c; Tribble et al., RPP 77 (2014) 106901, sec. 4.1), and across that
window :math:`\rho = p_{xA} a` changes -- for a deuteron by up to
:math:`T_s = p_s^2/2\mu_{sx} = 1.7` MeV at 40 MeV/c on top of
:math:`E + B = 2.2`-3 MeV. The vertex
:math:`M_l(\rho) = (B_c - 1)j_l(\rho) - \rho j_l'(\rho)\,(+\,C_l)` sees this,
and most strongly near its nodes.

*The measure at fixed* :math:`E`. In :math:`A + a \to s + F^*`,
:math:`F^* \to b + B`, the c.m. energy :math:`E` of :math:`x + A` fixes the
:math:`s`-:math:`F` energy, :math:`E_{sF} = E_{aA} - B_{xs} - E`, so
:math:`k_{sF}` has a fixed length, and the spectator momentum in :math:`a`
(2017 eqs. 29-30, :math:`\mathbf q = \mathbf p_{sx}`, with :math:`\beta =
m_s/m_a`)

.. math::

   q^2 = k_{sF}^2 + \beta^2 k_{aA}^2 - 2\beta k_{sF}k_{aA}\,x, \qquad
   x = \hat k_{sF}\cdot\hat k_{aA},

is fixed by the direction of :math:`\mathbf k_{sF}`: at fixed :math:`E`,
:math:`|p_s|` and the spectator direction are one variable, and only
:math:`|k_{sF} - \beta k_{aA}| \le q \le k_{sF} + \beta k_{aA}` is reached.
The derivation of the weight of a datum, step by step:

1. *Cross section.* The PWIA triple differential cross section in the c.m.
   variables (2017 eq. 33; Typel & Baur eqs. 13-14) is

   .. math::

      \frac{d^3\sigma}{d\Omega_{bB}\,d\Omega_{sF}\,dE} = K(E)\,\overline{|\mathcal M|^2},
      \qquad K(E) = \frac{\mu_{aA}\mu_{sF}\mu_{bB}\,k_{sF}k_{bB}}{(2\pi)^5 k_{aA}}\quad(\hbar = 1),

   and :math:`K` depends on :math:`E` only. Integrated over
   :math:`\Omega_{bB}` (2017 eq. 34),

   .. math::

      \frac{d^2\sigma}{d\Omega_{sF}\,dE} \propto |\phi_a(q)|^2 \sum_l P_l^{-1} S_l(E)\,
      |\tilde M_l(p_{xA})|^2 \equiv |\phi_a(q)|^2\,\sigma(E; q),

   per unit solid angle of :math:`\mathbf k_{sF}`, with :math:`p_{xA}^2/2\mu_{xA}
   = E + B_{xs} + q^2/2\mu_{sx}` (2017 eq. 31). The phase space at fixed
   :math:`E` is :math:`d\Omega_{sF}`, uniform.
2. *The lab form and KF.* The measured variables are, e.g., :math:`(E_b,
   \Omega_b, \Omega_B)`, and :math:`d^3\sigma/dE_b\,d\Omega_b\,d\Omega_B = KF\,
   |\phi_a(q)|^2\,d\sigma^\mathrm{HOES}/d\Omega` (Tribble 2014 eq. 4.1;
   Spitaleri et al. 2017 eq. 3). KF "contains the final-state phase-space
   factor": it is :math:`K` times the Jacobian of
   :math:`(E, \Omega_{sF}, \Omega_{bB}) \to (E_b, \Omega_b, \Omega_B)`
   (Typel & Baur eqs. 15-16; Tumino et al., ARNPS 71 (2021) 345, eqs.
   15-16), so :math:`KF\,dE_b\,d\Omega_b\,d\Omega_B = K(E)\,dE\,d\Omega_{sF}\,
   d\Omega_{bB}`. KF varies across the acceptance in the lab variables only
   because the lab volume element does; in the c.m. variables it is
   :math:`K(E)`. The same holds for the ratio :math:`\lambda_3/\lambda_2` of
   Pizzone et al., PRC 83 (2011) 045801 eqs. 3-10 (:math:`\lambda_2` depends
   on :math:`E` alone).
3. *The data.* The HOES datum of an energy bin is the yield of the bin
   divided by the Monte Carlo integral of :math:`KF\,|\phi_a|^2` over the
   same bin and acceptance (Spitaleri 2017 eq. 15 and text: "the product
   :math:`KF\cdot|\Phi|^2_\mathrm{exp}` is calculated by using a Monte Carlo
   simulation"; Pizzone 2011 p. 045801-5). With the detection probability
   :math:`\varepsilon`,

   .. math::

      \bar\sigma_\mathrm{bin} = \frac{\int_\mathrm{bin} \varepsilon\,KF\,|\phi_a|^2\,
      \frac{d\sigma^\mathrm{HOES}}{d\Omega}\,dE_b\,d\Omega_b\,d\Omega_B}
      {\int_\mathrm{bin} \varepsilon\,KF\,|\phi_a|^2\,dE_b\,d\Omega_b\,d\Omega_B}
      = \frac{\int dE\,K(E) \int d\Omega_{sF}\,A\,|\phi_a(q)|^2\,\sigma(E; q)}
             {\int dE\,K(E) \int d\Omega_{sF}\,A\,|\phi_a(q)|^2},

   :math:`A(\theta)` the acceptance of a spectator direction integrated over
   its azimuth and over the :math:`b`-:math:`B` direction. It is a ratio of
   integrals: the division by :math:`|\phi_a|^2` is not event by event, and
   :math:`|\phi_a|^2` stays in the weight (the datum is the
   :math:`|\phi_a|^2`-weighted mean; only a cross section constant across the
   accepted :math:`q` makes it drop out). KF does not appear: it is the
   Jacobian. AZURE2 evaluates the model at fixed :math:`E` and folds it with
   the resolution, which takes care of the bin; :math:`K(E)` varies by
   :math:`O(\Delta E/E_{sF})` across one.
4. *At fixed* :math:`E`, with :math:`d\Omega_{sF} = d\cos\theta_\mathrm{cm}\,
   d\varphi` (:math:`\theta_\mathrm{cm}` the angle of :math:`\mathbf k_{sF}`
   to the beam, :math:`x = \pm\cos\theta_\mathrm{cm}`) and
   :math:`d\cos\theta_\mathrm{cm} = q\,dq/(\beta k_{sF}k_{aA})`:

   .. math::

      \bar\sigma(E) = \frac{\int A\,|\phi_a(q)|^2\,\sigma(E; q)\;q\,dq}
                           {\int A\,|\phi_a(q)|^2\;q\,dq}
      \qquad (\text{plane-wave vertex}),

   over the accepted :math:`q`: inside :math:`[p_\mathrm{min},
   p_\mathrm{max}]`, reached at :math:`E`, and in the directions of the
   acceptance. The weight is :math:`|\phi|^2 q\,dq`, not :math:`|\phi|^2
   q^2\,dq`.

   With the DW vertex and :math:`R(E)` the amplitude depends on
   :math:`\mathbf k_{sF}` through :math:`x` and the distorted waves, not on
   :math:`q` alone, and :math:`|\phi_a|^2` is replaced by the plane-wave limit
   :math:`|\tilde\varphi(q)|^2` of the vertex's own bound state (the DW
   source is normalized by :math:`4\pi\tilde\varphi(q)`): the same measure,

   .. math::

      \langle |M_l|^2 \rangle = \frac{\int A\,|\tilde\varphi(q)|^2\,
      |M_l(x)|^2\,dx}{\int A\,|\tilde\varphi(q)|^2\,dx}, \qquad
      \bar R \propto \frac{\int A\,|M(x)|^2\,dx}{\int A\,|M_\mathrm{PW}(x)|^2\,dx},

   ("Experimental acceptance" below). Both assume that the :math:`|\phi|^2`
   of the Monte Carlo is the model's (the ``ps`` distribution, or
   :math:`\tilde\varphi` for the DW vertex); otherwise the datum carries the
   smooth factor :math:`\int A|\phi_\mathrm{model}|^2/\int
   A|\phi_\mathrm{MC}|^2`.

*The isotropic measure* :math:`|\phi|^2 p^2\,dp` *(AZURE2 before October
2026).* :math:`p^2\,dp\,d\Omega_p` is :math:`d^3q = d^3k_{sF} = \mu_{sF}k_{sF}\,
dE_{sF}\,d\Omega_{sF}`: the measure of events integrated over :math:`E` as
well, e.g. the :math:`|p_s|` spectrum of all coincidences. A HOES datum is
at one :math:`E` bin, where only the sphere :math:`|\mathbf k_{sF}| =
\mathrm{const}` is left; the isotropic weight also put nodes at :math:`q`
the kinematics do not reach at that :math:`E`. No HOES data convention
corresponds to it, so it was removed rather than kept as an option. (Dividing
each event by its own :math:`KF\,|\phi|^2` and normalizing by the accepted
phase space, which no extraction we know of does, would give
:math:`A\,K(E)/KF\,d\Omega_{sF}`: the fixed-:math:`E` measure again, not
:math:`p^2\,dp`.) The experimental momentum distribution
:math:`|\Phi(p_s)|^2_\mathrm{exp} \propto Y/KF` in narrow :math:`E` and angle
windows (Spitaleri 2017 eq. 9; Pizzone 2011 Fig. 8) is :math:`|\phi|^2`
itself -- what a ``ps`` table holds.

*Implementation.* With ``ps=`` on the experiment line every THM point of the
experiment (and every sub-point of its resolution folding) takes the
accepted directions at its own :math:`E` -- Gauss-Legendre nodes in
:math:`\cos\theta_\mathrm{cm}` on every accepted interval, the directions
``ThmDistortion::AngleNodes`` gives :math:`R(E)` and the DW vertex --
evaluates the entrance vertex (the surface term and the Coulomb term
:math:`C_l`, which depends on :math:`p_{xA}` too) at :math:`T_k =
q_k^2/2\mu_{sx}` added to :math:`E + B` (what ``spectatorEnergy`` does with
one value), and sums the HOES cross section of the nodes with the normalized
weights :math:`w_k \propto \omega_k A_k |\phi(q_k)|^2`, :math:`\omega_k` the
Gauss-Legendre weights in :math:`\cos\theta_\mathrm{cm}`. The level matrix
is common to the nodes; only the vertex changes. A ``ps`` window alone
accepts every direction whose :math:`q` lies in it (``psNodes`` nodes); with
``spectatorAngles=`` ("Experimental acceptance") the window is the
intersection, with ``spectatorAngleNodes`` nodes per c.m. interval: one
acceptance for the plane-wave vertex, :math:`R(E)` and the DW vertex. A
folding sub-point beyond the reach of the window (no accepted direction)
takes the nodes of the nearest data point; a data point without one is
refused.

*Syntax* (momenta in MeV/c, :math:`0 \le p_{\min} \le p_{\max}`):

``ps=delta``
   The quasi-free point (default): the vertex at :math:`E + B`, or at
   :math:`E + B + T_s` with ``spectatorEnergy``. Byte-identical to no key.
``ps=hulthen:pmin-pmax``, ``ps=hulthen:a,b:pmin-pmax``
   :math:`\phi(p) \propto (a^2 + q^2)^{-1} - (b^2 + q^2)^{-1}`, :math:`q =
   p/\hbar c`, the Hulthén function of the deuteron with :math:`a = 0.2317`,
   :math:`b = 1.202` fm\ :sup:`-1` (Tribble 2014 eq. 4.4), or other
   :math:`0 < a < b` -- the same form is the Eckart function used for
   :sup:`3`\ He and :sup:`6`\ Li (Tribble 2014 sec. 4.1), with :math:`a`
   usually :math:`\sqrt{2\mu_{sx}B_{xs}}/\hbar c`.
``ps=gauss:FWHM:pmin-pmax``
   :math:`|\phi(p)|^2 = \exp(-4\ln 2\,p^2/\mathrm{FWHM}^2)` (the FWHM of
   :math:`|\phi|^2` along a line through :math:`p = 0`, as momentum
   distributions are quoted).
``ps=table:<file>``
   The momentum distribution :math:`|\phi(p)|^2` itself (any scale; e.g. the
   experimental :math:`|\Phi(p_s)|^2_\mathrm{exp} \propto Y/KF` of Spitaleri
   2017 eq. 9) -- not a :math:`|p_s|` spectrum of accepted events, which
   contains the measure and the acceptance (divide it by :math:`p` and by the
   acceptance first). Until October 2026 the table was the event weight per
   unit :math:`p`; a table written for that convention changes meaning. Two
   columns, :math:`p` (MeV/c, strictly increasing, :math:`\ge 0`) and
   :math:`|\phi|^2 \ge 0`, ``#`` comments, at least two rows, some positive
   value; linear between rows. The window is the table's range. A relative path is taken from the directory of the
   ``.azr``, exactly as for ``weight[k]=`` files (the same resolution in
   ``Config::ReadThmBlock``: a leading ``/`` or ``\`` or a drive letter is
   absolute), so a copy of the project elsewhere must make it absolute or
   copy the table; the name cannot contain blanks or ``#``.
``psNodes=N``
   Gauss-Legendre nodes in :math:`\cos\theta_\mathrm{cm}` (1-64, default 16)
   on the directions inside the window, at every energy; with
   ``spectatorAngles=`` its ``spectatorAngleNodes`` instead (``psNodes`` is
   then refused). With :math:`p_{\min} = p_{\max}` the window is one
   direction of weight 1 (the vertex at :math:`T_s = p^2/2\mu_{sx}`). A table
   with kinks is integrated to :math:`O(h^2)` only; raise ``psNodes`` or
   smooth it.

:math:`\mu_{sx}` and the reachable :math:`q` come from the experiment's
kinematics (beam, target, spectator, :math:`x` = Trojan horse minus
spectator; built-in table masses), so a window needs ``beam``, ``target``,
``spectator`` and ``Ebeam``. Refused, with ``ERROR: <thm> experiment[...]``:
a window without kinematics, a malformed ``ps`` or ``psNodes``, ``psNodes``
without a window or with ``spectatorAngles``, an unreadable or invalid table,
a window together with a non-zero ``spectatorEnergy`` (global or for the
experiment's entrance pair) or with ``spectatorAngle`` (one direction), and a
data point at which no direction of the window is accepted (the window
beyond the reachable :math:`q`). The window is printed at startup (the
distribution, the cut, the nodes, :math:`\mu_{sx}`, and at the lowest and
the highest point the number of nodes, their :math:`|p_s|` range and
:math:`\langle T_s \rangle`) and written to ``thm_experiments.out``
(``ps_node`` rows at those two energies: :math:`E`,
:math:`\theta_\mathrm{cm}`, :math:`p_k`, :math:`w_k`, :math:`T_k`; a
``<T_s>`` row after each). Extrapolation segments (``<segmentsTest>``) are
in no experiment and keep the quasi-free vertex. Cost: the vertex and the
exit-channel sums once per node (the level matrix once). Memory: every point
and folding sub-point stores the vertex pieces of its entrance-pair channels
at its own nodes and their weights (16 bytes per channel and node, 8 per
node, shared with its mapped points); the total is the ``ps_table`` row of
``thm_experiments.out`` (points with sub-points, bytes). Until October 2026
the pieces of every channel were kept, about 117 MB per node for the 19F
model of ``examples/f19_pag_thm`` with its full window, and 16 nodes did not
fit in 2 GB; now that model needs 481 MB in all with 16 nodes in a pyazr
session (398 MB without a window; the CLI 636 and 520 MB, getrusage).

*pyazr.* ``AzrModel.set_thm_experiment(..., ps="hulthen:0-40", psNodes=16)``;
``session.thm_vertex(name, energies)`` returns, per energy (lists over
:math:`E` of arrays over its nodes), the nodes ``p_s``, ``weights``, ``T_s``,
``theta_cm`` and ``rho``, and per entrance channel
and level the window average :math:`\langle |M_l|^2 \rangle` and the
quasi-free :math:`|M_l(p_s = 0)|^2` with the boundary the vertex uses for
that level (for a ``delta`` experiment one node, :math:`p_s = 0`,
:math:`T` = its ``spectatorEnergy``).

*Checks* (``tests/thm_spectator_window``, ``tests/pyazr/thm_spectator_window_test.py``,
on ``tests/18O_p_a_thm`` with its deuteron Trojan horse, :math:`\mu_{sx} =
469.46` MeV, 2H(18O,α15N)n at 54 MeV, reachable :math:`|p_s|` from 3-8 to
beyond 130 MeV/c): ``ps=delta`` byte-identical to no key; a window shrinking
to a point equals ``spectatorEnergy`` :math:`= p^2/2\mu_{sx}` (one node and a
:math:`2\times 10^{-4}` MeV/c window: :math:`4\times 10^{-10}`); a flat
:math:`|\phi|^2` table on [20, 40] MeV/c equals :math:`\int q\,\sigma\,dq/
\int q\,dq` by Simpson's rule over 41 single-``spectatorEnergy`` runs to
:math:`1.6\times 10^{-8}` (the Simpson error); one direction
(``hulthen:30-30``) equals a ``spectatorEnergy`` session to
:math:`10^{-9}`; the Hulthén window [0, 40] with 16 nodes agrees with 32 to
:math:`1.3\times 10^{-10}` (8 nodes: :math:`1.6\times 10^{-10}`); ``ps``
alone is byte-identical to ``spectatorAngles=cm:0-180`` with as many nodes;
``thm_vertex`` at every energy against an independent evaluation of the
three-body kinematics (nodes in :math:`\cos\theta_\mathrm{cm}`, weights
:math:`\omega_k|\phi(q_k)|^2`, :math:`T_k`, :math:`\rho`: :math:`10^{-10}`)
and of :math:`\langle|M_0|^2\rangle`; refusals (out of reach, with
``spectatorAngle``, ``psNodes`` with ``spectatorAngles``, and the others
above). ``tests/reference/thm_spectator_angles_test``: a ``ps`` window alone
has the measure of every direction inside the cut (brute force,
:math:`4\times 10^{-6}`). Without ``ps`` and ``spectatorAngles`` every output
file is byte-identical to the previous build (``tests/18O_p_a_thm`` in eight
configurations -- no experiment, kinematics, R(E) Coulomb and optical at
qf, lab and c.m. angles, DW, ``theta``, ``spectatorEnergy`` -- and the
7Li, 12C+12C, 19F and 18O(p,α) 2010 examples).

*Size* (models at fixed parameters, not refitted; scratch studies, not in
the repository).

- A vertex node. With the spectator at a fixed :math:`p_s = 24` MeV/c,
  ``tests/18O_p_a_thm`` (no folding) puts :math:`M_0 = 0` at :math:`E \approx
  0.70` MeV: the model there drops to :math:`2.4\times 10^{-4}` of its
  maximum. The Hulthén window [0, 40] fills it to 0.48 of the maximum
  (:math:`5\times 10^{4}` times the point value; 0.42 with the isotropic
  measure). Across the whole data range the window changes the (folded) model
  by up to 67 % (73 %) against the quasi-free point.
- :sup:`7`\ Li(p,α): the HOES model of ``examples/li7_tumino2006`` as it
  was before its THM points were refitted as on-shell-equivalent cross
  sections (October 2026). It used :math:`B = 2.2246` MeV and is placed
  here in a deuteron Trojan horse, 2H(7Li,αα)n at 19 MeV, with Hulthén
  [0, 40]. That is the set-up of Lattuada et al. (2001), not of these data.
  The fixed-:math:`E` kinematics bound what such a window can mean: at
  19 MeV :math:`E_{sF} > 0` needs :math:`E < 2.01` MeV, and above 1.92 MeV
  the smallest reachable :math:`|p_s|` is above 40 MeV/c (43-68 MeV/c at
  1.92 MeV), so the window is refused there. The THM points reach 6.9 MeV,
  which the isotropic measure did not notice.

  On the 18 THM points below 1.83 MeV (the segment cut at
  :math:`E_\mathrm{lab}` = 2.15 MeV):

  * the THM model rises by 1.12-1.17 below 1 MeV and by 1.18-1.31 up to
    1.82 MeV;
  * :math:`\langle T_s\rangle` goes from 0.69 MeV at 0.08 MeV to 1.61 MeV
    at 1.82 MeV, where only 38-40 MeV/c is reachable;
  * the THM :math:`\chi^2` (norm profiled) goes 341 → 310.

  Refitted jointly with the direct data, the window raises the model's
  bare :math:`S(0)` from 45.9 to 50.8 keV b. On the scale of the direct
  data (:math:`S(0)` divided by the fitted low-energy norm) it goes only
  from 61.5 to 63.0 keV b.

  The points were measured differently. Tumino et al. (2006) used a
  33 MeV :sup:`3`\ He beam on :sup:`7`\ Li, with the spectator d and
  :math:`B = 5.4935` MeV, which puts the quasi-free point at
  :math:`E = 4.15` MeV. With Hulthén :math:`a = 0.42`, :math:`b = 1.2`
  fm\ :sup:`-1` and the paper's [0, 30] MeV/c, the window is reachable at
  every point from 0.08 to 6.9 MeV. The smallest reachable :math:`|p_s|` is
  2.3 MeV/c near 4 MeV and 28 MeV/c at 0.08 MeV; :math:`\langle T_s\rangle`
  is 0.34-0.67 MeV. A joint HOES refit with that window changes the total
  :math:`\chi^2` from 1105 to 1063 and the direct-scale :math:`S(0)` from
  66.5 to 66.4 keV b.

  An earlier version of this paragraph took a 33 MeV :sup:`7`\ Li beam on
  :sup:`3`\ He (:math:`E_{qf} = -1.35` MeV), the inverse of the experiment,
  and found the window out of reach above 1.62 MeV. Its numbers describe
  that set-up only: on the 15 points below 1.51 MeV, 0.98-1.16 below 1 MeV,
  1.17-1.18 above, and :math:`\chi^2` 218 → 176. With the isotropic measure
  on all points the deuteron set-up gave 1.14-1.29 and :math:`\chi^2`
  847 → 933, and that inverse set-up 0.91-1.15 and 1154 → 1227.
- 12C+12C (``examples/c12c12_tumino2018``, 12C(14N,α/p)d at 30 MeV, the four
  THM segments one experiment with a free norm; :math:`B = 10.27` MeV,
  :math:`\mu_{sx} = 1606` MeV; the d-12C distribution in 14N as an Eckart
  function with :math:`a = \sqrt{2\mu_{sx}B}/\hbar c = 0.918` fm\ :sup:`-1`,
  :math:`b = 2.5` fm\ :sup:`-1`): with [0, 40] MeV/c the kinematics reach
  only 38.9-40 MeV/c at the lowest point (0.82 MeV; :math:`\langle
  T_s\rangle` = 0.48 MeV) and 4.5-40 MeV/c at the highest (2.69 MeV, 0.25
  MeV), so :math:`\rho \approx 13.4`-14.7 moves by up to 0.25, a fair part of
  the oscillation of :math:`j_l` for :math:`l \le 8`. The model changes by
  0.19-2.26 level by level (p\ :sub:`0` 0.21-1.10, p\ :sub:`1` 0.27-2.26,
  α\ :sub:`0` 0.19-1.60, α\ :sub:`1` 0.26-1.18), and the THM :math:`\chi^2`
  at the published parameters goes 61 → 2863 (1578 with the isotropic
  measure); with the 80 MeV/c that the setup accepts ("Experimental
  acceptance") 7406. ``thm_vertex``: the :math:`l = 4` vertex has a node at
  :math:`E = 2.02` MeV (:math:`|M_4|^2 = 8\times 10^{-6}` against a maximum of
  0.68), which the window fills to :math:`4.3\times 10^{-2}` (4.8\ :math:`\times
  10^{-2}` isotropic). Published 12C+12C THM fits assume the quasi-free
  vertex; with the window the parameters would move.

Distortion factor R(E)
----------------------

*Physics.* The PWA analysis of a THM experiment takes the transfer amplitude
of :math:`a + A \to s + F^*` as constant, so the HOES excitation function is
the yield divided by the kinematic factor and the momentum distribution
:math:`|\phi_a(p_{sx})|^2` alone ("weight" above). With the distortions of
the initial (:math:`a + A`) and final (:math:`s + F`) relative motion the
amplitude varies with :math:`E` -- for a charged spectator below the
:math:`s + F` barrier by orders of magnitude, through the penetrability of
the spectator in the exit channel (Mukhamedzhanov & Pang, PRC 99 (2019)
064618; Mukhamedzhanov, arXiv:2609.04498). ``distortion=coulomb|optical`` on
an experiment line computes that variation and multiplies the model of every
segment of the experiment by it, before the folding, like ``weight[k]``.

*Amplitude.* Zero-range prior form (2019 eqs. 20-24; 2026 eqs. 22-28): with
the :math:`x`-:math:`A` vertex of zero range, :math:`\mathbf r_{sF} =
\mathbf r_{sx} \equiv \mathbf r` and :math:`\mathbf r_{aA} = \beta\mathbf
r`, :math:`\beta = m_s/m_a`, and

.. math::

   M(E) = \int d^3r\, \chi^{(-)*}_{\mathbf k_{sF}}(\mathbf r)\,
          \phi_{sx}(r)\, \chi^{(+)}_{\mathbf k_{aA}}(\beta\mathbf r)
        = \frac{4\pi}{k_{sF}\,\beta k_{aA}} \sum_l (2l+1)\,
          e^{i(\sigma_l^{aA} + \sigma_l^{sF})} P_l(x)
          \int_{r_\mathrm{min}}^\infty dr\, \phi_{sx}(r)\,
          u_l^{sF}(k_{sF} r)\, u_l^{aA}(\beta k_{aA} r),

with :math:`\chi^{(-)*}_{\mathbf k} = \chi^{(+)}_{-\mathbf k}` (2026 eq. 23),
:math:`x = \hat k_{sF}\cdot\hat k_{aA}`, :math:`\sigma_l` the Coulomb
phases and :math:`u_l` the regular radial waves normalized to :math:`u_l \to
F_l + T_l H_l^+` (:math:`u_l = F_l` for point Coulomb). The :math:`s`-:math:`x`
bound state is taken in :math:`l_{sx} = 0` from its tail: the Whittaker
function :math:`\phi = W_{-\eta_b,1/2}(2\kappa r)/r` (default,
``boundState=whittaker``, :math:`\eta_b = Z_sZ_x\alpha\mu_{sx}/\hbar\kappa`)
or the Yukawa :math:`e^{-\kappa r}/r` (``yukawa``), :math:`\kappa =
\sqrt{2\mu_{sx}B_{xs}}/\hbar`, zero below :math:`r_\mathrm{min}`
(``boundState=whittaker:3``, fm; 2026 uses its tail for :math:`r \ge 3` fm).
The normalization of :math:`\phi` cancels in :math:`R`. The plane-wave limit
is the Fourier transform of :math:`\phi` at the spectator momentum in
:math:`a`, :math:`\mathbf q = \mathbf k_{sF} - \beta\mathbf k_{aA}`,

.. math::

   M_\mathrm{PW}(E) = 4\pi\int dr\, r^2 j_0(qr)\,\phi_{sx}(r),

the momentum distribution by which the PWA data reduction divides
(:math:`4\pi/(\kappa^2 + q^2)` for the Yukawa tail).

*The factor.*

.. math::

   R(E) = \frac{\rho(E)}{\rho(E_\mathrm{ref})}, \qquad
   \rho = \frac{|M|^2}{|M_\mathrm{PW}|^2}\ (\texttt{distortionRatio=dwpw},
   \text{default})
   \quad\text{or}\quad \rho = |M|^2\ (\texttt{dw}),

and the model is multiplied by :math:`R`: the PWA-extracted quantity is
:math:`R\,\sigma^\mathrm{HOES}`, so multiplying the model by :math:`R` is the
same as dividing the published :math:`S^*` by :math:`R` (2026 eqs. 29-30; the
2019 eq. 24 and both papers' figures show :math:`1/R`, the factor that
multiplies :math:`S^*`, with :math:`E_\mathrm{norm} = 2.664` MeV). ``dwpw`` is
the correction to data that were divided by :math:`|\phi_a(p_{sx})|^2` at
the event's momentum: without distortion :math:`R \equiv 1`. ``dw`` is the
papers' ratio, which also carries the energy dependence of
:math:`|\phi_a(q(E))|^2` at a fixed angle (they drop :math:`\phi_a` as
energy independent). :math:`E_\mathrm{ref}` (``distortionRef=``, MeV,
default the middle of the experiment's data range) only sets the scale, which
the profiled norm absorbs: :math:`\chi^2` does not depend on it.

*Kinematics* (from ``beam``, ``target``, ``spectator``, ``Ebeam``, nuclear
masses of the built-in table): :math:`E_{aA}` the non-relativistic beam-target
c.m. energy, :math:`E_{sF} = E_{aA} - B_{xs} - E` (as for the line shape; a
data point with :math:`E_{sF} \le 0` is refused), :math:`F = x + A`,
:math:`B_{xs}` from the masses. If the masses and field 32 of the entrance
pair disagree on :math:`B_{xs}` by more than 1 keV AZURE2 prints a
``WARNING`` (the vertex uses field 32, these kinematics the masses).
``spectatorAngle=`` sets :math:`x`:

``qf`` (default)
   :math:`\hat k_{sF} = \hat k_{aA}`, the direction of least :math:`q`:
   :math:`\theta_s = 0` in the c.m. (and in the lab) when the Trojan horse
   is the beam, :math:`180^\circ` when it is the target.
``<deg>``
   a lab angle of the spectator to the beam, converted to the c.m. at every
   :math:`E` (non-relativistic, forward branch :math:`\theta_\mathrm{cm} =
   \theta_\mathrm{lab} + \arcsin(\gamma\sin\theta_\mathrm{lab})`,
   :math:`\gamma = V_\mathrm{cm}/v_s`); an angle the spectator cannot reach
   at a data point is refused (on the grid beyond the data, the largest
   reachable angle is used).
``cm:<deg>``
   a c.m. angle to the beam, fixed.

*Distorted waves.* ``distortion=coulomb``: point-Coulomb waves in both
channels. ``distortion=optical``: per channel, ``opticalAA=`` (:math:`a + A`)
and ``opticalSF=`` (:math:`s + F`), each ``coulomb`` (default), ``plane`` (no
distortion at all, as the 2026 Figs. 5-6 switch one channel off) or ten
numbers ``V,R,a,W,RW,aW,WD,RD,aD,RC``, or the name of a built-in global
optical potential (``ancai06``, ``kd03``, ... -- "Global optical potentials"
below):

.. math::

   U(r) = V_C(r) - V f(r; R, a) - iW f(r; R_W, a_W)
          - 4iW_D \frac{e^{(r - R_D)/a_D}}{(1 + e^{(r - R_D)/a_D})^2},
   \qquad f = \frac{1}{1 + e^{(r - R)/a}},

depths in MeV (positive: attractive and absorptive), radii in fm (not
reduced radii), :math:`V_C` a uniform sphere of radius ``RC`` or, for
``RC=0``, a point charge. A depth 0 switches its term off. A global
potential is evaluated at the channel's energy: fixed for :math:`a + A`, at
:math:`E_{sF}(E)` of every tabulated energy for :math:`s + F`.

*Numerics* (``src/ThmDistortion.cpp``). The radial equations are integrated
outward from the origin by Numerov's method (complex for an optical
potential), started from the power series, with step :math:`h = 0.02` fm for
:math:`r` (smaller if a local wave number exceeds :math:`0.1/h`) and
:math:`\beta h` for :math:`a + A`, so that :math:`u^{aA}` is read at
:math:`\beta r` on its own grid; beyond the potential and the turning point
each wave is matched at two points a quarter wavelength apart to
:math:`c_1F_l + c_2G_l` from AZURE2's Coulomb functions (COUL, as CoulFunc)
and divided by :math:`c_1 - ic_2`. The radial integrals are Simpson sums from
:math:`r_\mathrm{min}` to :math:`r_\mathrm{min} + 50/\kappa`: the bound state
makes them absolutely convergent, so no complex rotation or damping is needed
(unlike the external Coulomb term :math:`C_l`); the size of the cut is
checked and a ``WARNING`` printed above :math:`10^{-8}`. The :math:`a + A`
waves are computed once, up to the :math:`l` where :math:`(2l+1)\int|\phi
u_l^{aA}|` has fallen below :math:`10^{-17}` of its maximum; at each energy
the sum stops once the bound on the remaining terms is below
:math:`10^{-13}|M|` (12C(14N,d): :math:`l \le 23`; a deuteron Trojan horse,
:math:`\kappa = 0.23` fm\ :sup:`-1`: :math:`l \approx 50`-80). :math:`\ln R`
is tabulated at 10 keV steps from 0.5 MeV below the lowest to 0.5 MeV above
the highest point (short of the spectator threshold) and interpolated by
cubic Lagrange (error :math:`< 10^{-6}`); points beyond take the end value
with one ``WARNING``. Cost: about 1 s at startup for 12C(14N,d), nothing per
evaluation. With ``dwpw`` and a cutoff :math:`r_\mathrm{min} > 0`,
:math:`M_\mathrm{PW}` can pass through zero at large :math:`q`; AZURE2 warns,
and ``dw`` avoids the singularity.

*Output.* The startup summary gives :math:`k_{aA}`, :math:`\eta_{aA}`,
:math:`\kappa`, :math:`\eta_b`, :math:`\beta`, :math:`E_\mathrm{ref}`, the
grid, the highest :math:`l` and :math:`R` at the ends of the data;
``thm_experiments.out`` a ``distortion:`` line and ``distortion_point`` rows
(:math:`E`, :math:`E_{sF}`, :math:`\eta_{sF}`, :math:`\theta_\mathrm{cm}`,
:math:`|M|^2`, :math:`|M_\mathrm{PW}|^2`, :math:`R`, :math:`l_\mathrm{max}`)
at the lowest point, :math:`E_\mathrm{ref}` and the highest point.
``pyazr``: ``session.thm_distortion(name, energies)``;
``AzrModel.set_thm_experiment(..., distortion="optical",
opticalAA=[V, R, a, W, RW, aW, WD, RD, aD, RC], opticalSF="coulomb",
spectatorAngle=8, distortionRef=2.664, distortionRatio="dw",
boundState="whittaker:3")``. ``distortion=table:<file>`` (the ``weight[k]``
format and interpolation, relative to the ``.azr``, every point inside the
table) multiplies all segments of the experiment by its :math:`w(E)`; it needs
no kinematics. A segment that also has ``weight[k]`` gets both (warned).
Extrapolation segments are in no experiment and carry no distortion.

*Checks.* ``tests/reference/thm_distortion_test`` (ctest ``thm_distortion``):
:math:`|M|^2` and :math:`M_\mathrm{PW}` against an independent evaluation
(``tests/reference/thm_distortion_reference.py``: mpmath Coulomb functions and
adaptive quadrature; scipy's DOP853 for a complex Woods-Saxon in both
channels) -- point Coulomb to :math:`10^{-8}`, the optical case to
:math:`5\times 10^{-6}`, for 12C(14N,d) forward and at 90°, with a 3 fm
cutoff, and a made-up 18O(3He,d) with the Trojan horse as target; plane
waves in both channels equal :math:`4\pi/(\kappa^2 + q^2)` to
:math:`2\times 10^{-9}` at 0, 70 and 180°, so :math:`R \equiv 1`; nuclear
part off (ten zeros) equals point Coulomb exactly. ``tests/thm_distortion``
(CLI): ``distortion=none`` byte-identical to no key; a neutral spectator with
a plane :math:`a + A` wave gives :math:`R = 1` to :math:`5\times 10^{-11}`; the
model ratio equals the directly evaluated :math:`R`; the 12C+12C factors
below; refusals. ``tests/pyazr/thm_distortion_test.py``: ``AzrModel``, CLI ==
session, ``thm_distortion``, a table.

*12C+12C: the published curves.* 12C(14N,α/p)d at 30 MeV (Tumino et al.,
Nature 557 (2018) 687): :math:`E_{aA} = 13.84` MeV, :math:`B = 10.272` MeV,
:math:`E_{sF} = 3.573 - E`, :math:`\eta_{aA} = 4.52`, :math:`\eta_{sF} =
1.55`-2.74 over :math:`E` = 0.8-2.66 MeV, :math:`\kappa = 0.921`
fm\ :sup:`-1`, :math:`\eta_b = 0.39`, :math:`\beta = 1/7`. Against the
factors that multiply :math:`S^*`, :math:`1/R` normalized at 2.664 MeV,
digitized from the vector figures (digitisation 0.03 dex), over 0.8-2.55 MeV:

================================================ ==================== ===================
AZURE2 settings                                  2019 Fig. 10          2026 Fig. 9
                                                 (FRESCO, Coulomb)     (forward angles)
================================================ ==================== ===================
``coulomb`` (qf, whittaker, dwpw), Ebeam 30      0.026 / 0.046 dex     0.12 / 0.16
same, ``boundState=yukawa:3``                    0.007 / 0.015         0.13 / 0.21
qf, whittaker, ``dw``, Ebeam 30                  0.040 / 0.076         0.11 / 0.14
qf, whittaker, ``dw``, Ebeam 30.11               0.16 / 0.23           0.007 / 0.023
qf, whittaker, dwpw, Ebeam 30.11                 0.14 / 0.20           0.012 / 0.019
================================================ ==================== ===================

(rms / max of :math:`\log_{10}` differences.) The zero-range Coulomb DWBA at
the quasi-free angle reproduces the 2019 FRESCO curve within the digitisation
with the default settings (:math:`1/R = 2.6\times 10^{-3}` at 0.8 MeV,
:math:`3.3\times 10^{-3}` at 1.0, 0.30 at 2.5, against their 2.3, 3.0 and
:math:`302\times 10^{-3}`). The 2026 forward-angle curve is
steeper near the normalization point; at the nominal kinematics AZURE2 lies
0.11 dex (a factor 1.3) below it at 0.8-1.5 MeV, and it is reproduced to 0.007
dex if :math:`E_{sF}` is 50 keV higher (``Ebeam=30.11``, :math:`E_{aA}` + 51
keV), so the two papers differ mainly in the :math:`s`-:math:`F` energy their
kinematics assign, where the barrier penetrability is steep (:math:`E_{sF}
\approx 0.9` MeV at the normalization point). Angles other than forward do
worse for 2019 (0.09-0.4 dex); against 2026 at the nominal kinematics
``cm:90`` with ``dw`` comes to 0.030 dex, which we take for a coincidence
(the paper computes 0° and 8° and finds them nearly equal, as AZURE2 does:
8° in the lab gives 0.083 dex). Which value of :math:`E_{sF}` the 2026
calculation used is not stated.

*Size* (models at the published parameters, not refitted; the four THM
segments of ``examples/c12c12_tumino2018`` one experiment with a free norm,
30 keV folding). ``distortion=coulomb`` (defaults): :math:`R` = 5.4 at 0.82
MeV to 0.012 at 2.69 MeV (:math:`E_\mathrm{ref}` = 1.75 MeV), the folded model
changes by a factor 0.012-5.4 across the data (a factor 440), and the THM
:math:`\chi^2` goes 61 → 3247 (192 points); with the 2026-like ``dw``,
Ebeam 30.11: 61 → 2999. Published 12C+12C :math:`S^*` would move by 2-3
orders of magnitude below 1.5 MeV. For a light system with a neutron
spectator, 2H(18O,α15N)n at 54 MeV (``examples/o18_lacognata2010``, only the
:math:`d + {}^{18}`\ O Coulomb wave, :math:`\eta_{aA} = 0.73`, is
distorted): :math:`R` = 0.93-1.08 over 0.51-0.89 MeV and the THM
:math:`\chi^2` 414 → 397 (``dw``: 398); for the made-up charged 18O(3He,d)
at 115 MeV (:math:`E_{sF} \approx 10` MeV, above the barrier): :math:`R` =
0.998-1.002 at the quasi-free angle, 0.99-1.01 at 30° in the lab.

*Scope and limits.* Zero range at the :math:`x`-:math:`A` vertex, the
:math:`s`-:math:`x` bound state in :math:`l_{sx} = 0` and in its asymptotic
form, one spectator angle or, with ``spectatorAngles``, the acceptance in
its polar angle ("Experimental acceptance" below), no
spin-orbit term, no three-body (line-shape) Coulomb
effects -- those are ``lineshape=on``, which multiplies independently and
does not overlap with :math:`R` ("Coulomb effects: what each option
contains" below).
Whether a DWBA vertex should replace the PWA one at all is the question the
papers debate; the factor lets a fit show what it would change. The
:math:`l`-dependent vertex of the same amplitude, of which :math:`R(E)` is the
zero-range limit, is ``vertexModel=dw`` (next section).

Distorted-wave entrance vertex
------------------------------

*Why.* The plane-wave vertex :math:`M_l = (B - 1)\, j_l(pa) - pa\, j_l'(pa)`
sets how much each entrance partial wave contributes, and the ratio between
partial waves is the largest model dependence of a THM analysis that
converts peak areas into strengths. For 19F(p,αγ)16O (19F + d at 55 MeV)
:math:`|M_1(213\,\mathrm{keV})/M_0(324\,\mathrm{keV})|^2` is 1.10, 3.00 and
19.7 at :math:`a_p` = 4.1, 5.1 and 6.1 fm, and the THM strength of the
213 keV resonance follows it. The distortion factor :math:`R(E)` of the
previous section multiplies every partial wave by the same number and cannot
change such a ratio. ``vertexModel=dw`` on an experiment line replaces
:math:`M_l` by the vertex computed from the distorted waves of :math:`a + A`
and :math:`s + F`, with the experiment's distortion settings
(``distortion=coulomb|optical``, ``opticalAA``, ``opticalSF``,
``spectatorAngle``, ``boundState``). ``vertexModel=pw`` (the default) is the
plane-wave vertex of the first section, byte for byte.

*Amplitude.* Three structureless particles :math:`s`, :math:`x`, :math:`A`;
the Trojan horse :math:`a = (s\,x)` is bound in :math:`l_{sx} = 0` with wave
function :math:`\varphi(r_{sx})`, and :math:`F = x + A`. The prior-form DWBA
amplitude of :math:`a + A \to s + F^*` is (Mukhamedzhanov, PRC 84 (2011)
044616, sec. III; Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020) 233,
eqs. 19-20)

.. math::

   M = \bigl\langle \chi^{(-)}_{sF}\, \Upsilon_{xA} \bigm| U_{sA} + V_{xA} - U_{aA}
       \bigm| \varphi\, \chi^{(+)}_{aA} \bigr\rangle ,

:math:`\Upsilon_{xA}` the :math:`x`-:math:`A` scattering state with
:math:`F^*`. Split the :math:`\mathbf r_{xA}` integral at :math:`r_{xA} = a`.
Inside, the operator is :math:`[V_{xA} + U_{sF}] + (U_{sA} + V_{sx} -
U_{sF}) - [V_{sx} + U_{aA}]`, and the bracketed terms are the Hamiltonians of
the final and the initial channel at the same total energy. Green's theorem
in :math:`\mathbf r_{xA}` (the kinetic energy is :math:`T_{xA} + T_{sF}`,
and the :math:`\mathbf r_{sF}` integral runs over all space) turns the
internal prior amplitude into the internal post amplitude plus a surface
term (2011 eqs. 36-38; 2020 eqs. 28-31). The internal post amplitude is
dropped, as it is for the plane-wave vertex (2020: "we disregard
:math:`M_\mathrm{int}^{DW(post)}`"), which leaves

.. math::

   M \simeq M_S + M^{\mathrm{prior}}_{\mathrm{ext}} .

In the Jacobi coordinates :math:`\mathbf r \equiv \mathbf r_{xA}` and
:math:`\mathbf u \equiv \mathbf r_{sx}`,

.. math::

   \mathbf r_{sF} = \alpha\,\mathbf r + \mathbf u, \qquad
   \mathbf r_{aA} = \mathbf r + \beta\,\mathbf u, \qquad
   \alpha = \frac{m_A}{m_F}, \quad \beta = \frac{m_s}{m_a}

(2020 eq. 34), and :math:`d^3r_{sF} = d^3u` at fixed :math:`\mathbf r`. The
*source* of the vertex is the initial state projected on the final
spectator wave at fixed :math:`\mathbf r_{xA}`:

.. math::

   S(\mathbf r) = \int d^3r_{sF}\; \chi^{(-)*}_{\mathbf k_{sF}}(\mathbf r_{sF})\,
                  \varphi(r_{sx})\, \chi^{(+)}_{\mathbf k_{aA}}(\mathbf r_{aA})
                = \int d^3u\; \varphi(u)\,
                  \chi^{(-)*}_{\mathbf k_{sF}}(\alpha\mathbf r + \mathbf u)\,
                  \chi^{(+)}_{\mathbf k_{aA}}(\mathbf r + \beta\mathbf u) .

The surface term is the Wronskian of :math:`\Upsilon` and :math:`S` on the
sphere :math:`r = a` (2020 eqs. 31-32; the derivative is taken at fixed
:math:`\mathbf r_{sF}`, which is how :math:`S` is defined). Outside,
:math:`\Upsilon` of channel :math:`c = (s, l, J)` is the R-matrix external
wave :math:`O_l(kr)/O_l(ka)` times the level-matrix combination
:math:`\sum_{\lambda\lambda'} \gamma_{\lambda f} A_{\lambda\lambda'}
\gamma_{\lambda' c}` of the first section. With

.. math::

   S_{lm}(r) = \int d\Omega_r\, Y^*_{lm}(\hat r)\, S(\mathbf r),

the surface term of channel :math:`c` and projection :math:`m` is that
combination times

.. math::

   V_{lm}(B) = (B - 1)\, S_{lm}(a) - a\, S'_{lm}(a),
   \qquad B = a\,\frac{O_l'(ka)}{O_l(ka)} ,

and the external prior term adds

.. math::

   E_{lm} = \frac{2\mu_{xA}}{\hbar^2}\int_a^\infty dr\, r\, \frac{O_l(kr)}{O_l(ka)}
            \bigl[(V^C_{xA} + U_{sA} - U_{aA})\, S\bigr]_{lm}(r)

(:math:`[\dots]_{lm}` the projection of the operator acting inside the
:math:`\mathbf u` integral). Replacing the on-shell logarithmic derivative
:math:`B` by a boundary constant of the R-matrix (``vertex=constant``,
``perlevel``) is the same step as for the plane-wave vertex: it concerns
:math:`\Upsilon`, not the source, so ``vertex=`` works as before, and
``onshell`` gives the complex :math:`B = S_c + iP_c`.

*Plane-wave limit.* With plane waves in both channels, the substitution
:math:`\mathbf u \to \mathbf r_{sF} - \alpha\mathbf r` and
:math:`\alpha\beta + m_x (m_a + m_A)/(m_a m_F) = 1` give

.. math::

   S_\mathrm{PW}(\mathbf r) = \tilde\varphi(q)\, e^{i\mathbf p\cdot\mathbf r},
   \qquad
   \mathbf q = \mathbf k_{sF} - \beta\,\mathbf k_{aA}, \quad
   \mathbf p = \mathbf k_{aA} - \alpha\,\mathbf k_{sF}, \quad
   \tilde\varphi(q) = \int d^3u\, e^{-i\mathbf q\cdot\mathbf u}\varphi(u) :

:math:`\mathbf p` is the off-shell :math:`x`-:math:`A` momentum and
:math:`\mathbf q` the :math:`s`-:math:`x` momentum of 2017 eqs. 29-30 (2020
eq. 36). Since :math:`\alpha/\mu_{xA} = \beta/\mu_{sx} = 1/m_x`,

.. math::

   \frac{p^2}{2\mu_{xA}} - \frac{q^2}{2\mu_{sx}}
   = \frac{k_{aA}^2}{2\mu_{aA}} - \frac{k_{sF}^2}{2\mu_{sF}} = E + B_{xs} ,

the spectator relation of the ``ps`` window (2017 eq. 31); with nuclear
masses (:math:`m_a < m_s + m_x`) it holds to :math:`O(B_{xs}/m_a c^2)`, about
:math:`10^{-3}` for a deuteron. The plane-wave expansion
:math:`e^{i\mathbf p\cdot\mathbf r} = 4\pi\sum_{lm} i^l j_l(pr)\,
Y^*_{lm}(\hat p)Y_{lm}(\hat r)` gives

.. math::

   S^\mathrm{PW}_{lm}(r) = 4\pi\,\tilde\varphi(q)\, i^l\, Y^*_{lm}(\hat p)\, j_l(pr),
   \qquad
   V^\mathrm{PW}_{lm}(B) = 4\pi\,\tilde\varphi(q)\, i^l\, Y^*_{lm}(\hat p)\,
   \bigl[(B - 1)\, j_l(pa) - pa\, j_l'(pa)\bigr],

exactly the plane-wave vertex :math:`M_l(p)`. The operator of the external
term is then :math:`V^C_{xA}` alone (no :math:`s`-:math:`A` or
:math:`a`-:math:`A` interaction), :math:`[V^C_{xA} S]_{lm} = (\hbar^2/2\mu)
(2\eta k/r)\, S^\mathrm{PW}_{lm}`, and :math:`E_{lm}` becomes
:math:`4\pi\tilde\varphi\, i^l Y^*_{lm}(\hat p)\, C_l(E)` with the Coulomb
term of ``coulombIntegral`` (Tribble 2014 eq. 2.79; Typel & Baur eq. A.4). The
distorted-wave vertex contains the present one as its limit, :math:`M_l` and
:math:`M_l + C_l`, up to the factor :math:`4\pi\tilde\varphi(q)\, i^l\,
Y^*_{lm}(\hat p)`, whose :math:`l` dependence drops out of the observable
(next paragraph).

*Orbital projections and the observable.* A plane-wave source has only
:math:`m = 0` along :math:`\hat p`; a distorted one does not, since it
depends on the directions of both :math:`\mathbf k_{aA}` and
:math:`\mathbf k_{sF}`. The angle-integrated HOES observable sums over the
exit direction and all spin projections, and
:math:`\sum_{M m_s}\langle s\,m_s\,l\,m|J\,M\rangle
\langle s\,m_s\,l'\,m'|J\,M\rangle = \delta_{ll'}\delta_{mm'}(2J+1)/(2l+1)`
makes it incoherent in :math:`l` *and* :math:`m`. The vertex squared of the
first section is therefore replaced by

.. math::

   |M_l|^2 \;\to\; \frac{4\pi}{2l+1}\sum_{m=-l}^{l}
   \Bigl|\frac{V_{lm}(B)}{4\pi\tilde\varphi(q)}\Bigr|^2 ,

which does not depend on a quantization axis and equals :math:`|M_l(p)|^2`
for plane waves, since :math:`\sum_m |Y_{lm}(\hat p)|^2 = (2l+1)/4\pi`.
:math:`V_{lm}` is linear in :math:`B`: with :math:`s_m = S_{lm}(a)/4\pi
\tilde\varphi` and :math:`d_m = a S'_{lm}(a)/4\pi\tilde\varphi`,

.. math::

   \frac{4\pi}{2l+1}\sum_m |c_1 s_m + c_2 d_m|^2 = c^\dagger G\, c, \qquad
   G = \frac{4\pi}{2l+1}\begin{pmatrix} \sum|s_m|^2 & \sum s_m^* d_m \\
                                   \sum d_m^* s_m & \sum |d_m|^2\end{pmatrix},

:math:`c = (B - 1, -1)`. With a level-dependent :math:`B` (``perlevel``) the
level sum sits inside the modulus; writing :math:`G = L^\dagger L`
(Cholesky, pivoted on the larger diagonal) the amplitude splits into two
components :math:`M^{(k)}(B) = L_{k1}(B - 1) - L_{k2}`, each summed over the
levels coherently and squared, and the two are added: the vertex of an
entrance channel becomes two incoherent buckets. For plane waves :math:`G`
has rank one, the second component vanishes and the first is :math:`M_l`
up to a phase. At the ``qf`` direction (below) only :math:`m = 0` enters and
:math:`G` has rank one with distortion too.

*Finite range in* :math:`s`-:math:`x`. The source keeps the full
:math:`\mathbf r_{sx}` dependence, with the bound state of the distortion
factor (``boundState=whittaker|yukawa[:rmin]``: the Whittaker or Yukawa tail,
:math:`l_{sx} = 0`). A zero-range :math:`\varphi` would be a poor
approximation here: in the prior form :math:`\varphi` itself enters, not
:math:`V_{sx}\varphi`, and for a deuteron (:math:`1/\kappa = 4.3` fm) it
spreads :math:`\mathbf r_{aA}` by :math:`\beta/\kappa \approx 2` fm and
:math:`\mathbf r_{sF}` by 4 fm around a surface point at 4-6 fm, which is
where the distortions act. The plane-wave limit holds for any
:math:`\varphi`, and with finite range the source at the centre is
:math:`S(\mathbf 0) = M(E)` of the distortion factor exactly, which ties the
two together.

*Normalization, and the relation to* :math:`R(E)`. The data are divided by
:math:`|\varphi(p_s)|^2` at the event's momentum, so the vertex is divided by
:math:`4\pi\tilde\varphi(q)` of its own kinematics (the ``dwpw`` convention
of :math:`R`); without distortion it is the plane-wave vertex at
:math:`p(E)`. The distortion factor is the same source at the centre of
:math:`F`: :math:`R_\mathrm{dwpw}(E) \propto |S(\mathbf 0)/S_\mathrm{PW}(\mathbf
0)|^2`. If the distortion were constant across the surface,
:math:`S(\mathbf r) = [S(\mathbf 0)/\tilde\varphi]\, S_\mathrm{PW}(\mathbf r)`
for :math:`r \le a`, the DW vertex would be :math:`R^{1/2}` times the
plane-wave one for every :math:`l`, and the model would be :math:`R(E)` times
the plane-wave model. That factorization is the step of 2020 eqs. 35-39,
where the off-shell :math:`\mathbf p_{xA}` of the distorted waves' momentum
distribution is replaced by :math:`\mathbf k_{xA} = \mathbf k_{aA} -
\alpha\mathbf k_{sF}` and the surface term becomes the zero-range amplitude
times the plane-wave off-shell factor :math:`W_l`: :math:`R(E)` is the
:math:`r_{xA} \to 0`, :math:`l`-independent limit of the DW vertex.
``vertexModel=dw`` evaluates eq. 32 without that replacement, and the
distortion of the source across the surface differs between partial waves,
most visibly through the local momentum of the Trojan horse in the
:math:`a + A` potential. The two carry the same overall energy dependence
(12C+12C below), so they are never applied together: with
``vertexModel=dw`` the experiment's distortion settings build the vertex
and :math:`R(E)` is not applied; ``distortionRatio`` and ``distortionRef``,
which only concern :math:`R`, are refused.

*The external term is not computed.* Outside :math:`a` the DW operator is the
three-body remnant :math:`V^C_{xA} + U_{sA} - U_{aA}`. For point charges
(2020 eq. 26) it is :math:`e^2 Z_A (Z_x/r_{xA} + Z_s/r_{sA} - Z_a/r_{aA})`,
whose monopole cancels (:math:`Z_x + Z_s = Z_a`): the :math:`x`-:math:`A`
Coulomb force that :math:`C_l` describes in the plane-wave limit is mostly
already in the :math:`a + A` Coulomb wave, and what remains falls off like a
dipole, :math:`1/r^2`. Its nuclear part needs the :math:`s`-:math:`A`
optical potential, which is not an input here, and the source would be
needed on :math:`a \le r < \infty`, where the number of partial waves grows
with :math:`r`. The DW vertex is therefore the surface term alone: the
counterpart of ``coulombIntegral=0`` (the default, and the setting of every
example), to which it reduces without distortion. ``coulombIntegral=1``
together with ``vertexModel=dw`` is refused. 2011 and 2020 describe
:math:`M^\mathrm{prior}_\mathrm{ext}` as small and neglect it "in some cases
with a reasonable choice of the channel radius".

*Kinematics of a node.* Without a ``ps`` window the vertex is taken at the
spectator direction of ``spectatorAngle`` (default ``qf``,
:math:`\hat k_{sF} = \hat k_{aA}`), as :math:`R(E)` is; :math:`q = |\mathbf
k_{sF} - \beta\mathbf k_{aA}|` is then small but not zero, and the vertex
without distortion is the plane-wave one with ``spectatorEnergy``
:math:`= q^2/2\mu_{sx}` at that energy. With a ``ps`` window (or
``spectatorAngles``) the vertex is averaged over the accepted directions
("Experimental acceptance"): :math:`q` and the direction are one variable,
:math:`\cos\theta = (k_{sF}^2 + \beta^2 k_{aA}^2 - q^2)/(2\beta k_{sF}
k_{aA})`, only :math:`|k_{sF} - \beta k_{aA}| \le q \le k_{sF} + \beta
k_{aA}` is reached at a given energy, and the nodes are Gauss-Legendre in
:math:`\cos\theta_\mathrm{cm}` on the directions whose :math:`q` is inside
the cut, weight :math:`|\tilde\varphi(q)|^2\,d\cos\theta_\mathrm{cm}` (the
fixed-:math:`E` measure, "Spectator-momentum window"); without distortion
this is the plane-wave window with :math:`\tilde\varphi` for :math:`\phi`.
(Until October 2026 the nodes were Gauss-Legendre in :math:`q` with the
weight :math:`|\phi|^2 q^2`.) ``spectatorAngle`` together with a window is
refused (the window sets the direction), and a data point whose window is
out of reach is refused at startup.

*Numerics* (``src/ThmDwVertex.cpp``). With the distorted waves in partial
waves (the conventions of the previous section,
:math:`\chi^{(+)}_{\mathbf k}(\mathbf R) = (4\pi/kR)\sum_L i^L e^{i\sigma_L}
u_L(kR) \sum_M Y^*_{LM}(\hat k) Y_{LM}(\hat R)` and
:math:`\chi^{(-)*}_{\mathbf k} = \chi^{(+)}_{-\mathbf k}`), rotational
invariance reduces the six-dimensional surface integral to reduced
amplitudes that do not depend on the directions of :math:`\mathbf k_{aA}`
and :math:`\mathbf k_{sF}`:

.. math::

   S_{lm}(r) = \frac{(4\pi)^2}{k_{sF} k_{aA}} \sum_{L_s L_a} (-i)^{L_s} i^{L_a}
   e^{i(\sigma^{sF}_{L_s} + \sigma^{aA}_{L_a})}\, h^l_{L_s L_a}(r)\,
   \bigl\{Y_{L_s}(\hat k_{sF}) \otimes Y_{L_a}(\hat k_{aA})\bigr\}^*_{lm},

.. math::

   h^l_{L_s L_a}(r) = \sqrt{\frac{4\pi}{2l+1}}\; 2\pi \int u^2 du\,
   \varphi(u) \int_{-1}^{1} d\cos\theta_u\; f_{L_s}(R_s)\, f_{L_a}(R_a)\,
   K^l_{L_s L_a}(\theta_s, \theta_a),

with :math:`f_L(R) = u_L(kR)/R`, :math:`\mathbf r` along :math:`z` and
:math:`\mathbf u` in the :math:`xz` plane, :math:`\mathbf R_s = \alpha r\hat
z + \mathbf u` and :math:`\mathbf R_a = r\hat z + \beta\mathbf u` at polar
angles :math:`\theta_s`, :math:`\theta_a` (the azimuth of :math:`\mathbf u`
integrates to :math:`2\pi`; only :math:`L_s + L_a + l` even contributes).
:math:`K = \sum_M \langle L_s\,M\,L_a\,{-M}|l\,0\rangle\,\bar Y_{L_s
M}(\theta_s)\, \bar Y_{L_a, -M}(\theta_a)`, :math:`\bar Y_{LM}(\theta) =
Y_{LM}(\theta, 0)`, is the zero component of a bipolar harmonic; rotating
:math:`\hat R_s` onto :math:`z` leaves

.. math::

   K^l_{L_s L_a} = \sqrt{\frac{2L_s+1}{2l+1}} \sum_{|\nu| \le l} (-1)^\nu
   \langle L_s\,0\,L_a\,\nu|l\,\nu\rangle\, \bar Y_{l\nu}(\theta_s)\,
   \bar Y_{L_a\nu}(\theta_a - \theta_s),

:math:`O(l)` terms instead of :math:`O(L)`. :math:`S'_{lm}(a)` is the same
integral with the :math:`r` derivative of the integrand at fixed
:math:`\mathbf u`. In the frame :math:`\hat k_{aA} = \hat z` the bipolar
harmonic is :math:`\langle L_s\,m\,L_a\,0|l\,m\rangle\, \bar
Y_{L_s m}(\theta_{sa}) \sqrt{(2L_a+1)/4\pi}`, so every node of a window (an
angle :math:`\theta_{sa}` between :math:`\mathbf k_{sF}` and
:math:`\mathbf k_{aA}`) costs only this sum. The geometry and
:math:`f_{L_a}` do not depend on :math:`E`; per energy only
:math:`f_{L_s}` changes, and batches of energies share one pass over the
quadrature points. The waves are the radial waves of the distortion factor
(Numerov from the origin, matched to COUL; a plane wave is the
Riccati-Bessel function), tabulated at 0.02 fm (less where a local wave
number exceeds 5 fm\ :sup:`-1`) and read by six-point Lagrange
interpolation (value and derivative). The :math:`u` integral runs from
:math:`r_\mathrm{min}` to :math:`r_\mathrm{min} + 26/\kappa` (the Yukawa tail
beyond is below :math:`10^{-10}`) in 16-point Gauss-Legendre panels of at
most 4 radians of :math:`(k_{sF} + \beta k_{aA})u` and 5 fm;
:math:`\cos\theta_u` on :math:`1.2\,L_\mathrm{max} + 20` nodes;
:math:`L \le kR_\mathrm{max} + 14` in each channel, :math:`R_\mathrm{max} =
a + \beta u_\mathrm{max}` and :math:`\alpha a + u_\mathrm{max}`. The
Clebsch-Gordan coefficients come from the Racah sum in long double
(checked against GSL, and to :math:`10^{-12}` in orthonormality at
:math:`L = 60`). :math:`G` is tabulated per :math:`l` and node at 20 keV
steps from 0.3 MeV below the lowest point to 0.3 MeV above the highest (short
of the spectator threshold) and interpolated by cubic Lagrange (the node
weights of a window likewise, renormalized); points beyond take the end value
with one ``WARNING``. Cost at startup, two threads: 6 s for 19F + d
(:math:`\kappa = 0.23` fm\ :sup:`-1`, :math:`l \le 4`, 58 energies,
:math:`L \le 61`), 4 s for 12C(14N,d) (:math:`l \le 8`, 125 energies); per
evaluation the tables cost nothing measurable (a window averages the Gram
matrix of the directions, one node: the same time per evaluation as without
it). Memory: the 19F session peaks at about 400 MB with or without the
window.

*Syntax, output, pyazr.* ``vertexModel=pw|dw`` on the experiment line.
``dw`` needs ``distortion=coulomb`` or ``optical`` (and so the kinematics);
refused, with ``ERROR: <thm> experiment[...]``: any other ``vertexModel``,
``dw`` without a computed distortion, with ``distortionRef``,
``distortionRatio``, ``theta`` (the fixed-angle sum carries only
:math:`m_l = 0` about :math:`\hat p_{xA}`), ``spectatorAngle`` with a
``ps`` window, with ``coulombIntegral=1``, ``entranceL=coherent`` or a
``spectatorEnergy`` for the entrance pair, and a window out of reach at a
data point. The startup summary gives :math:`\alpha`, :math:`\beta`,
:math:`k_{aA}`, :math:`\eta_{aA}`, :math:`\kappa`, the entrance :math:`l`,
the partial waves, the quadrature, the grid and the time;
``thm_experiments.out`` a ``vertex:`` line and ``dw_vertex_point`` rows
(:math:`E`, :math:`q`, :math:`pa`, :math:`l`, :math:`G_{11}`,
:math:`G_{22}`, :math:`G_{12}`) at the lowest, middle and highest point.
``AzrModel.set_thm_experiment(..., distortion="optical",
opticalAA=[...], vertexModel="dw")`` checks the same rules (also
``set_thm_option`` against an existing dw experiment);
``session.thm_vertex(name, energies)`` returns ``model = "dw"``, ``M2`` (the
DW vertex over the nodes the engine uses), ``M2_qf`` (at the
``spectatorAngle`` direction), ``M2_pw`` (the plane-wave vertex at the same
:math:`p`), the nodes per energy (``dw_q``, ``dw_weights``) and
``dw_q_delta``, ``dw_p_delta``. The GUI keeps the key as written; its
angular-distribution diagnostic says it is not available with ``dw``.

*Checks.* ``tests/reference/thm_dw_vertex_test`` (ctest ``thm_dw_vertex``):
(b) plane waves in both channels give the plane-wave Gram matrix at
:math:`p = |\mathbf k_{aA} - \alpha\mathbf k_{sF}|` to :math:`10^{-9}` and
the ratios :math:`|M_l|^2/|M_0|^2` for real and complex :math:`B` to
:math:`10^{-8}`, for 19F(d,n) at the quasi-free direction (:math:`l \le 4`),
12C(14N,d) at 40° in the c.m. (:math:`l \le 8`, every :math:`m` enters) and
a ``ps`` window: the directions with :math:`q` in the cut, their
:math:`q` and weights :math:`\omega_k|\tilde\varphi(q_k)|^2` (against an
independent radial quadrature of :math:`\tilde\varphi`, :math:`10^{-9}`),
and the averaged :math:`G` against the weighted plane-wave Gram matrices;
(c) point Coulomb in :math:`d` + 19F: :math:`G` for :math:`l = 0` (324 keV)
and :math:`l = 1` (213 keV) at 5.136 fm against
``thm_dw_vertex_reference.py``, a direct three-dimensional quadrature of
:math:`S(\mathbf r)` with the closed-form Coulomb wave
:math:`e^{-\pi\eta/2}\Gamma(1+i\eta)e^{ikz}\,{}_1F_1(-i\eta, 1, ik(R-z))`
(mpmath), projected on :math:`Y_{l0}` -- no partial waves, Numerov or tables
shared -- to :math:`2\times 10^{-7}` of the largest entry (the reference is
converged to :math:`3\times 10^{-8}`); the Cholesky components reproduce
:math:`c^\dagger G c`; interpolation between grid nodes to :math:`10^{-6}`.
``tests/thm_dw_vertex/check.sh`` (CLI, ``tests/18O_p_a_thm``, 2H(18O,α15N)n
at 54 MeV): ``vertexModel=pw`` byte-identical to no key, with and without
``distortion=coulomb``; plane waves: the ``dw_vertex_point`` rows are
:math:`j_0^2`, :math:`(\rho j_0')^2`, :math:`j_0\rho j_0'` to
:math:`10^{-9}` and the model is the plane-wave one within 5 % (3.5 % found; the
kinematic :math:`p(E)` against :math:`q = 0`, largest near the node of
:math:`M_0`); Coulomb: no :math:`R(E)`, a window changes the model; eleven
refusals. ``tests/pyazr/thm_dw_vertex_test.py``: ``AzrModel`` and its
refusals, CLI == session, ``thm_vertex`` (plane waves: ``M2_qf`` ==
``M2_pw`` to :math:`10^{-8}`), and without folding the model ratio dw/pw
equals the vertex ratio of ``thm_vertex`` at every point (:math:`10^{-8}`),
also with a window.

*Size: the radius dependence of the* :math:`l` *ratios*. 19F(p,αγ)16O,
``examples/f19_pag_thm`` (19F + d at 55 MeV, neutron spectator,
``vertex=constant``: :math:`B` is :math:`S_l` at the lowest level of the
:math:`J^\pi` group, the −448 keV 1\ :sup:`+` for :math:`l = 0` and the 213
keV 2\ :sup:`−` for :math:`l = 1`), :math:`a_p` of the p + 19F pairs varied,
nothing refitted. Optical potentials: :math:`d` + 19F with the An & Cai (PRC
73 (2006) 054605) global parameters at :math:`E_d` = 5.8 MeV, n + 20Ne
Koning-Delaroche at 2.6 MeV (outside its mass range, :math:`A \ge 24`), both
written out at one energy (``opticalAA=92.57,3.066,
0.753,1.467,3.581,0.588,10.65,3.711,0.696,3.477``, ``opticalSF=53.4,3.13,
0.675,0.39,3.13,0.675,6.86,3.53,0.541,0``); window: Hulthén 0-50 MeV/c (the
:math:`p_s` cut of Su et al., PRL 135 (2025) 182701).

======================== ====================== ====================== ======================
vertex                   :math:`|M_1(213)/`     :math:`|M_1(828)/`     :math:`|M_0(11)/`
                         :math:`M_0(324)|^2`    :math:`M_0(324)|^2`    :math:`M_0(324)|^2`
                         (4.1 / 5.1 / 6.1 fm)   (4.1 / 5.1 / 6.1 fm)   (4.1 / 5.1 / 6.1 fm)
======================== ====================== ====================== ======================
pw, :math:`q = 0`        1.10 / 3.00 / 19.7     1.16 / 2.83 / 15.7     1.24 / 1.55 / 3.17
pw, window               2.23 / 8.07 / 6.49     2.21 / 7.04 / 4.56     1.26 / 1.57 / 0.86
dw, plane waves (qf)     1.11 / 3.01 / 19.9     1.16 / 2.87 / 16.4     1.21 / 1.47 / 2.81
dw, Coulomb (qf)         0.171 / 0.235 / 0.341  0.242 / 0.319 / 0.433  0.99 / 1.04 / 1.10
dw, Coulomb, window      0.36 / 0.52 / 0.81     0.43 / 0.60 / 0.88     0.98 / 1.03 / 1.10
dw, optical :math:`aA`   0.44 / 0.56 / 0.70     0.51 / 0.62 / 0.74     1.06 / 1.09 / 1.14
dw, optical (qf)         0.33 / 0.41 / 0.52     0.40 / 0.48 / 0.58     1.04 / 1.08 / 1.13
dw, optical, window      0.89 / 0.94 / 1.06     0.91 / 0.93 / 1.00     1.05 / 1.07 / 1.09
======================== ====================== ====================== ======================

The plane-wave ratios change by factors 18, 13 and 2.6 over the three radii;
with the distortion they change by 1.6-2.5 (quasi-free) and by 1.04-1.19
with the optical potentials and the window (the window rows: the accepted
directions with the fixed-:math:`E` measure; before October 2026 the
isotropic one gave pw 2.94 / 11.9 / 3.78, dw Coulomb 0.42 / 0.63 / 1.04, dw
optical 1.29 / 1.24 / 1.32 for the first column). The DW source is the Trojan horse
decelerated by the 19F Coulomb field where it breaks up: at 5-6 fm the
deuteron has about 3 of its 5.3 MeV, and :math:`p` (a difference of
:math:`k_{aA}` and :math:`\alpha k_{sF}`) roughly halves, which moves the
:math:`l = 1` to :math:`l = 0` ratio down by an order of magnitude; the
nuclear attraction at the surface gives part of the momentum back. The
:math:`l` ratios themselves are therefore model dependent at the level of the
optical potentials and of the window (a factor 1.5-2 between Coulomb and
optical waves at the quasi-free direction, 1.1-2.5 with the window, and
1.7-2.7 between the quasi-free direction and the window with optical waves), but
they hardly depend on the channel radius any more. The PWIA stripping conversion of Su et al.
corresponds to :math:`|M_1(213)/M_0(324)|^2 \approx 2.1` and
:math:`|M_1(828)/M_0(324)|^2 \approx 0.66` in these units (the plane-wave
vertex at 4.7 fm, and a conversion factor 4.3 times ours at 5.1 fm).

*Size: 12C+12C against* :math:`R(E)`. ``examples/c12c12_tumino2018``,
12C(14N,α/p)d at 30 MeV, point Coulomb in both channels, quasi-free: the
ratio :math:`F_l(E) = |M_l^\mathrm{dw}|^2/|M_l^\mathrm{pw}(p)|^2`,
normalized at 1.75 MeV, against the distortion factor :math:`R(E)` (dwpw,
same :math:`E_\mathrm{ref}`) over 0.85-2.65 MeV. The sum over the entrance
channels, :math:`\sum|M^\mathrm{dw}|^2/\sum|M^\mathrm{pw}|^2`, follows
:math:`R` to 0.06 dex rms (0.10 dex at most) while both fall by a factor 330:
the energy dependence of the vertex is the :math:`s + F` Coulomb barrier of
the deuteron, common to all :math:`l`, and :math:`R(E)` is indeed its
:math:`l`-summed limit. The single partial waves deviate from :math:`R` by
0.03 (:math:`l = 8`) to 0.43 dex rms (:math:`l = 2`, 0.15 for :math:`l = 0`,
0.32 for :math:`l = 6`); where the plane-wave vertex has a node (:math:`l =
4` at 1.95 MeV) the DW vertex has none, and the ratio is not defined there.

*Size: refits of 19F(p,αγ)16O.* The joint fit of ``examples/f19_pag_thm``
(THM of Su et al., JUNA and Spyrou direct data, and penalty rows for 17
direct :math:`\omega\gamma` and :math:`\Gamma`, all :math:`a_p` = 5.136 fm,
the THM energy shifted by −9.17 keV lab) refitted with each vertex, in the
adopted THM window (:math:`E \le 0.45` MeV, 28 points) and the full one (53
points); optical potentials and window as above. Strengths in eV; every fit
stopped at its evaluation limit (``least_squares``, 20-60 Jacobians, from
the plane-wave fit or a previous stage), so differences of a few units of
:math:`\chi^2` are not significant.

.. list-table::
   :header-rows: 1

   * - fit
     - THM χ²/N
     - JUNA, Spyrou
     - penalty (17)
     - ωγ(11) [1e-29]
     - ωγ(213)
     - ωγ(226)
     - ωγ(828)
   * - adopted window, pw
     - 63.8 / 28
     - 12.8, 5.2
     - 7.5
     - 3.68
     - 0.0116
     - 9.3e-5
     - 775
   * - adopted, pw + window
     - 61.3 / 28
     - 16.0, 3.4
     - 36.1
     - 3.68
     - 0.0058
     - 1e-8
     - 775
   * - adopted, dw Coulomb
     - 166.4 / 28
     - 14.6, 6.7
     - 102.0
     - 4.32
     - 0.0135
     - 4.9e-3
     - 776
   * - adopted, dw optical
     - 114.8 / 28
     - 14.2, 6.6
     - 60.7
     - 4.53
     - 0.0136
     - 4.0e-3
     - 776
   * - adopted, dw optical + window
     - 81.3 / 28
     - 13.0, 5.8
     - 3.8
     - 5.53
     - 0.0127
     - 4.0e-4
     - 775
   * - full window, pw
     - 204.0 / 53
     - 13.8, 3.6
     - 342.3
     - 5.04
     - 0.0104
     - 2.3e-4
     - 184
   * - full, pw + window
     - 217.4 / 53
     - 23.7, 3.1
     - 508.8
     - 5.01
     - 0.0049
     - 2e-5
     - 67
   * - full, dw Coulomb
     - 242.1 / 53
     - 19.8, 7.0
     - 122.4
     - 4.63
     - 0.0130
     - 5.0e-3
     - 727
   * - full, dw optical
     - 199.4 / 53
     - 18.0, 6.8
     - 83.0
     - 5.06
     - 0.0130
     - 4.3e-3
     - 667
   * - full, dw Coulomb + window
     - 170.6 / 53
     - 13.5, 5.6
     - 32.7
     - 6.70
     - 0.0118
     - 1.9e-3
     - 655
   * - full, dw optical + window
     - 204.1 / 53
     - 12.2, 5.2
     - 62.3
     - 7.22
     - 0.0118
     - 5.4e-4
     - 564
   * - direct (penalty rows)
     - 
     - 
     - 
     - (7.5 ± 3.0) NACRE
     - 0.0126(13)
     - 0.0011(4)
     - 775(35)

In the full window the plane-wave vertex cannot hold the 828 keV (l = 1)
strength at its direct value: the fit trades 184 eV (67 eV with the window)
against the THM peak shape. With the DW vertex the same fit keeps 667-727 eV
at the quasi-free direction (564 eV with optical waves and the window, 655 eV
with Coulomb waves and the window), the penalty rows fall from 342 to 83-122
(33-62 with the window) and the total :math:`\chi^2` from 566 to 309
(optical); the THM :math:`\chi^2` improves little (199-242 at the quasi-free
direction, 171-204 with the window, against 204), the 790 keV (l = 2)
strength stays at 0-10 eV against 17(5) in every vertex model, and with
optical waves and the window the 564 keV strength drops to 34 eV against
48(7). So the DW vertex removes most of the l = 1 conflict but the full
window still does not fit together with the direct strengths; the adopted
window stays. Inside it the DW vertex at the quasi-free direction is worse
than the plane-wave one (THM :math:`\chi^2` 115-166 against 64), mainly
because it lifts the 225 keV (l = 3) strength to 4-5 × 10\ :sup:`-3` eV
against 1.1(4) × 10\ :sup:`-3`; with the window it keeps the 213 keV
strength at its direct value with small penalties (3.8 against 7.5) but is
still worse in the THM :math:`\chi^2` (81 against 64).
:math:`\omega\gamma(11)` moves from 3.7 to 4.3-4.5 × 10\ :sup:`-29` eV
(quasi-free DW, at the upper end of the plane-wave radius range 2.3-4.4)
and to 5.5 with the window; the window with plane waves moves the 213 keV
strength to half its direct value (−5.2σ). The rows "with the window" are
refits (October 2026) with the fixed-:math:`E` measure of the ``ps`` window
and ``spectatorAngles=cm:135-180`` (equal to ``ps`` alone, "Experimental
acceptance"), each from the earlier best fit; the full-window dw optical fit
moved by 0.1 in :math:`\chi^2` over a further 60 evaluations. The fits made
before October 2026 with the isotropic measure :math:`|\phi|^2 p^2\,dp`
gave 76.7 / 223.3 (dw optical, adopted / full window) and 57.8 / 278.5 (pw)
for the THM :math:`\chi^2`, ωγ(11) 5.39 and ωγ(828) 493 eV (dw optical).

*Scope and limits.* The surface term of the prior-form DWBA with the
internal post amplitude neglected (as for the plane-wave vertex) and without
the external prior term (``coulombIntegral=1`` refused); the angle-integrated
observable only (no ``theta``); the :math:`s`-:math:`x` bound state in
:math:`l_{sx} = 0` and in its asymptotic form; one spectator direction or
the acceptance in the spectator's polar angle
(``spectatorAngles`` and ``ps``, "Experimental acceptance"); no
spin-orbit term; the optical potentials are an input whose choice now
carries the :math:`l` dependence. With a cutoff :math:`r_\mathrm{min} > 0`
the normalization :math:`\tilde\varphi(q)` can pass through zero at large
:math:`q` (as :math:`M_\mathrm{PW}` of the distortion factor can), where the
normalized vertex is singular; an exact zero on the grid is refused. ``vertexModel=pw`` against ``dw`` (with
Coulomb and with optical waves) belongs in the model-dependence protocol of
any analysis that converts peak areas of different :math:`l` into strengths.

Global optical potentials
-------------------------

*Why.* With ``vertexModel=dw`` the ratios between entrance partial waves
hardly depend on the channel radius any more, but they depend on the distorted
waves (previous section: a factor 1.5-3 between Coulomb and optical waves).
To scan that dependence with standard choices instead of typing ten numbers
per channel, ``opticalAA=`` and ``opticalSF=`` also take the name of a
built-in global parametrisation (``src/ThmOptical.cpp``), for the distortion
factor and the DW vertex alike:

================ ========= ======================================================= =========== ====================
name             for       reference                                               target A    :math:`E_\mathrm{lab}`
                                                                                               (MeV)
================ ========= ======================================================= =========== ====================
``ancai06``      d         An & Cai, PRC 73 (2006) 054605                          12-238      0-183
``daehnick80``   d         Daehnick, Childs & Vrcelj, PRC 21 (1980) 2253           27-238      11.8-90
``kd03``         n, p      Koning & Delaroche, NPA 713 (2003) 231, global set      24-209      0.001-200
``bg71``         t, 3He    Becchetti & Greenlees (1971; RIPL-3 7100, 8100)         40-208      1-40
``liang09``      3He       Liang, Li & Cai, J. Phys. G 36 (2009) 085104            9-208       0-270
``mcfadden66``   4He       McFadden & Satchler, NPA 84 (1966) 177                  16-208      1-25
``avrigeanu94``  4He       Avrigeanu, Hodgson & Avrigeanu, PRC 49 (1994) 2136      16-250      1-73
================ ========= ======================================================= =========== ====================

The validity ranges are the papers' (as summarised in RIPL-3, Capote et al.,
Nucl. Data Sheets 110 (2009) 3107, and for the deuteron sets in Chin. Phys. C, doi:10.1088/1674-1137/acb2bc
(2023)); a lower energy limit that a paper does not
state is 0. Heavy ions have no entry: there is no standard Woods-Saxon
global set for them (the São Paulo potential is a double folding; the
Akyüz-Winther Woods-Saxon parametrisation of the proximity potential exists
in several versions and has no imaginary part), and the imaginary part near
the barrier is a modelling choice -- 14N + 12C and the like take the ten
numbers.

*Which nucleus, which energy.* The projectile is the partner of the channel
the model is made for (:math:`a` or :math:`A` in :math:`a + A`, :math:`s` or
:math:`F = x + A` in :math:`s + F`), the other is the target, whose
:math:`A`, :math:`Z` and :math:`N` set the depths and radii; a model that
describes neither partner is refused. The energy is the projectile's lab
energy on the target at rest at the channel's c.m. energy,
:math:`E_\mathrm{lab} = E_\mathrm{cm}(m_p + m_t)/m_t` (non-relativistic, as
all THM kinematics in AZURE2): :math:`E_{aA}` for :math:`a + A`, fixed by
the beam, and :math:`E_{sF} = E_{aA} - B_{xs} - E` for :math:`s + F`, which
changes across the data. The :math:`s + F` potential is therefore evaluated
anew at every energy at which a wave is computed -- each node of the
:math:`\ln R` grid (10 keV), of the DW vertex grid (20 keV), and every direct
evaluation (``thm_distortion``, the output rows) -- and the radial step takes
the deepest potential over :math:`0 < E_{sF} \le E_{sF}(E_\mathrm{low})`. For
19F + d (neutron spectator) the n + 20Ne potential goes from
:math:`E_n` = 3.30 MeV at the lowest THM point to 2.73 MeV at the highest
(:math:`W_D` 7.05 → 6.90 MeV); for 12C(14N,d), d + 24Mg from 3.0 to 1.0 MeV.

*Conventions* (mapped onto ``V,R,a,W,RW,aW,WD,RD,aD,RC``). Radii are
:math:`r_0 A_t^{1/3}` in fm (target mass number only, as in all these
papers); the Coulomb term is a uniform sphere of radius
:math:`r_C A_t^{1/3}` (none for neutrons). The surface absorption of every
model is the derivative form :math:`-4a_DW_D\,df/dr = 4W_De^x/(1+e^x)^2`
(Koning & Delaroche, An & Cai, Daehnick, Liang), exactly AZURE2's,
so :math:`W_D` maps one to one; Becchetti-Greenlees, McFadden-Satchler and
Avrigeanu have volume absorption only. Dropped or changed: the spin-orbit
terms (real and, for KD, imaginary) -- the distorted waves have no spin;
negative imaginary depths are set to 0 (Liang's volume term below 23 MeV, on
light targets the surface term carries the absorption); KD's relativistic
kinematics (its energy argument is the lab energy as here, its wave numbers
are not used); KD's Coulomb correction of the proton real depth,
:math:`V_C v_1(v_2 - 2v_3(E - E_F) + 3v_4(E - E_F)^2)` with
:math:`V_C = 1.73 Z/r_C A^{1/3}`, is kept. ``daehnick80`` is the global set
as the FRONT front end of TWOFNR (J. A. Tostevin, Surrey) codes it --
:math:`V = 88.5 - 0.26E + 0.88ZA^{-1/3}`, :math:`a = 0.709 + 0.0017E`,
:math:`W_V = (12.2 + 0.026E)(1 - e^{-(E/100)^2})`, :math:`W_D = (12.2 +
0.026E)e^{-(E/100)^2}`, :math:`a_I = 0.53 + 0.07A^{1/3} - 0.04\sum_i
e^{-((M_i - N)/2)^2}` over the neutron magic numbers -- with
non-relativistic kinematics; the paper's relativistic variant is not
offered (we could not check its table labels against the paper itself).

*Validity.* At startup the target mass and the projectile's lab energy over
the data (both ends for :math:`s + F`) are compared with the range above;
outside, the experiment is refused with the range in the message. Written as
``<name>:extrapolate`` the potential is used there anyway and AZURE2 prints a
``WARNING`` with what lies outside. Refusing by default keeps an
extrapolation visible in the input: THM channels are light and slow, and
several of them are outside every global set -- n + 20Ne is below the KD
mass range (:math:`A \ge 24`) and d + 19F below Daehnick's in mass and energy
(:math:`E_d` = 5.8 MeV), while An & Cai covers d + 19F and d + 24Mg
(12C(14N,d), :math:`E_d` = 1-3 MeV). The grid beyond the data (up to
0.5 MeV for :math:`R`, 0.3 MeV for the vertex) is evaluated without a check.

*Syntax, output, pyazr, GUI.* ``opticalAA=ancai06 opticalSF=kd03:extrapolate``
with ``distortion=optical``; a misspelt name or option is refused with the
list of names. The startup summary (and ``thm_experiments.out``) gives each
global channel with its parameters, for :math:`s + F` at both ends of the
data, e.g. ``s + F: kd03:extrapolate (n on 20Ne) Woods-Saxon V=53.3201 ...
at E_lab = 2.72825 MeV (E = 0.45 MeV) to V=53.0977 ... at E_lab = 3.29548
MeV (E = -0.09 MeV)``. ``AzrModel.set_thm_experiment(..., opticalAA="ancai06",
opticalSF="kd03:extrapolate")`` checks the name, the projectile, the target
mass range and the :math:`a + A` energy as the engine does (the :math:`s + F`
energy needs the data and is the engine's). The GUI's optical row offers
*plane*, *Coulomb*, *global* (with a combo of the names; the first one that
describes the channel's light partner is chosen) and *Woods–Saxon*; for a
global potential **Edit…** shows the reference, the validity range, the ten
numbers at the two ends of the data, and the *:extrapolate* switch, and the
button's tooltip the depths at the ends.

*Checks.* ``tests/reference/thm_optical_test`` (ctest ``thm_optical``): (a)
the ten numbers of every model at 23 (projectile, target, energy) points --
including d + 19F at 5.83 MeV, d + 24Mg at 2 MeV and n + 20Ne at 2.6 MeV --
against ``tests/reference/thm_optical_reference.py``, which reads the RIPL-3
library file ``om-parameter-u.dat`` and evaluates its entries by the
library's own coefficient formulas (standard and Koning forms: 6200, 2405,
5405, 7100, 8100, 9100, 9600) and transcribes FRONT21 for Daehnick and Liang
(RIPL's Daehnick entries 6112-6116 are polynomial refits of the imaginary
part for single nuclei; their real depth and :math:`a_I` agree with the
formulas to 4 digits for 40Ca, 90Zr, 120Sn and 208Pb, the absorptive depths
to 1 %), to :math:`10^{-9}`; the An & Cai d + 19F numbers reproduce the ten
numbers used before for ``examples/f19_pag_thm``. (b) Elastic scattering,
d + 40Ca at 56 MeV with ``daehnick80`` and no spin-orbit term: the ratio to
Rutherford from the engine's waves (``ThmDistortion::Wave``, Numerov + COUL,
which now also returns :math:`T_l`) agrees with scipy's DOP853 and mpmath's
Coulomb functions to :math:`1.3\times 10^{-4}` at 12 angles, and with the
measured ratio (Hatanaka et al. 1980, EXFOR E0682-022) within 0.16 dex at
15-78° (0.066 dex rms) -- the potential is used as published and the
solver handles it. (c) In the distortion factor the global :math:`s + F`
amplitude at :math:`E` equals the one with the ten numbers written out at
:math:`E_{sF}(E)` (identical); the refusals and the ``:extrapolate`` warning;
the parser. ``tests/thm_distortion`` (CLI): ``opticalAA=ancai06`` equals its
ten numbers at :math:`E_{aA}` in the model, ``kd03:extrapolate`` for n + 19F
is warned and reported at both ends, four refusals.
``tests/pyazr/thm_distortion_test.py`` and ``tests/gui/thm_workspace_test``:
the same checks in ``AzrModel`` and on the page, round trips byte for byte.

*Size: the spread of standard choices, 19F(p,αγ)16O.* ``examples/f19_pag_thm``
(19F + d at 55 MeV, :math:`a_p` = 5.136 fm, ``vertex=constant``, as in the
radius table of the previous section; the example before October 2026,
THM :math:`\chi^2` 63.9), ``vertexModel=dw``; :math:`l` ratios
from the whole THM range (so that 828 keV lies on the vertex grid), the THM
:math:`\chi^2` of the adopted window (28 points) at the example's parameters
(fitted with the plane-wave vertex; profiled norm and linear background,
nothing refitted). "custom" are the ten numbers used before (An & Cai and
Koning-Delaroche written out at :math:`E_d` = 5.83 MeV and :math:`E_n` = 2.6
MeV); the global :math:`s + F` potential instead follows :math:`E_n` = 3.30-2.73
MeV.

============================================ ================================ ===============================
:math:`a + A` / :math:`s + F`                quasi-free                       window (Hulthén 0-50 MeV/c)
                                             :math:`|M_1(213)/M_0(324)|^2`,   same, :math:`\chi^2`
                                             :math:`|M_1(828)/M_0(324)|^2`,
                                             :math:`|M_0(11)/M_0(324)|^2`,
                                             :math:`\chi^2`
============================================ ================================ ===============================
plane-wave vertex                            3.00, 2.83, 1.55; 63.9           11.9, 10.1; 1525
Coulomb / (neutral: plane)                   0.235, 0.333, 1.04; 680          0.634, 0.725; 479
custom / custom                              0.407, 0.486, 1.083; 611         1.242, 1.179; 176
``ancai06`` / ``kd03:extrapolate``           0.406, 0.487, 1.080; 611         1.242, 1.183; 177
``daehnick80:extrapolate`` / ``kd03:extr.``  0.402, 0.487, 1.079; 614         1.243, 1.194; 180
``ancai06`` / plane                          0.562, 0.632, 1.095; 541
============================================ ================================ ===============================

The two deuteron sets give the same ratios to 1 % (Daehnick far below its
range, :math:`E_d` = 5.8 MeV against 11.8 MeV, and below its mass range), and
the energy-dependent KD potential the same as its value at 2.6 MeV to
0.3 %. The steps that matter are elsewhere: point Coulomb to an optical
potential in :math:`d` + 19F (a factor 1.7), the n + 20Ne absorption (1.4),
and the momentum window (3). With standard choices the optical-model spread
of the 19F :math:`l` ratios is about 1 %, well below the window and vertex
model dependences.

*Size: 12C+12C.* ``examples/c12c12_tumino2018``, 12C(14N,α/p)d at 30 MeV,
the four THM segments one experiment with free norm (192 points), the
quasi-free direction, ``vertexModel=dw``, the THM :math:`\chi^2` at the
published parameters (61 with the plane-wave vertex). 14N + 12C has no
global entry: the :math:`a + A` wave is point Coulomb throughout, and the
d + 24Mg potential is varied (:math:`E_d` = 1.0-3.0 MeV; ``ancai06`` in
range, ``daehnick80`` extrapolated in mass and energy). Ratios at
:math:`E` = 1.5 MeV, and the energy dependence of the vertex,
:math:`\sum_l|M_l|^2` at 0.9 MeV over that at 2.6 MeV:

=============================== ====================================== ================== ============
:math:`s + F` (d + 24Mg)        :math:`|M_l/M_0|^2`, l = 2, 4, 6, 8    0.9 / 2.6 MeV      :math:`\chi^2`
=============================== ====================================== ================== ============
plane-wave vertex               0.726, 0.138, 0.331, 1.67              1.25               61
point Coulomb                   1.18, 1.09, 3.96, 4.75                 383                4713
``ancai06``                     0.059, 0.423, 0.420, 0.053             258                4598
``daehnick80:extrapolate``      0.058, 0.426, 0.434, 0.061             269                4624
plane (no s + F distortion)     0.316, 0.075, 0.695, 0.578             1.08               9396
=============================== ====================================== ================== ============

The deuteron absorption in d + 24Mg reshapes the partial-wave pattern
(:math:`|M_2/M_0|^2` from 1.18 to 0.06, :math:`|M_8/M_0|^2` from 4.8 to
0.06) and softens the barrier energy dependence by a third, but the two
standard deuteron sets agree to 3 % for :math:`l \le 6` (15 % for
:math:`l = 8`) and 4 % in the energy dependence. As for 19F, the spread among
standard global potentials is small next to the choice between Coulomb and
optical waves; for 12C+12C the remaining open input is the 14N + 12C
potential, which no global set covers and which stays the user's ten
numbers.

Experimental acceptance
-----------------------

*Why.* Without an acceptance the distortion factor :math:`R(E)` and the
distorted-wave vertex are evaluated at one spectator direction
(``spectatorAngle``: the quasi-free direction or one lab or c.m. angle). The
data integrate over the detector acceptance: ranges of the angles of the two
detected particles, and so of the spectator direction, together with the
accepted :math:`|p_s|` window. ``spectatorAngles=`` gives the accepted
spectator directions; with ``ps=`` (the :math:`|p_s|` cut and the momentum
distribution) it is one acceptance, over which the plane-wave vertex,
:math:`R(E)` and the DW vertex are averaged.

*Kinematics: one variable at fixed* :math:`E`. In :math:`a + A \to s + c + C`
through :math:`F^* \to c + C`, the energy :math:`E` of :math:`x + A` fixes
:math:`E_{sF} = E_{aA} - B_{xs} - E`, so :math:`k_{sF}` has a fixed length.
What remains at fixed :math:`E` is the direction of :math:`\mathbf k_{sF}`
(two angles) and the :math:`c`-:math:`C` direction (two angles), which the
angle-integrated HOES observable integrates. The DW vertex and
:math:`R(E)` depend on the directions only through
:math:`x = \hat k_{sF}\cdot\hat k_{aA}` (rotational invariance: the reduced
amplitudes of the previous sections and the partial-wave sum of :math:`R`
do not depend on the frame), and :math:`\hat k_{aA}` is the beam axis. So
the azimuth of the spectator does not enter the vertex at all; it enters
only the acceptance, which is therefore taken as integrated over the
azimuth at each polar angle (a full ring, or a pair of coplanar detectors,
gives the same vertex; only the acceptance per polar angle matters). The
direction also fixes both momenta of the vertex,

.. math::

   q^2 = k_{sF}^2 + \beta^2 k_{aA}^2 - 2\beta k_{sF}k_{aA}\,x, \qquad
   p^2 = k_{aA}^2 + \alpha^2 k_{sF}^2 - 2\alpha k_{aA}k_{sF}\,x,

:math:`q = |p_s|` the spectator momentum in :math:`a` and :math:`p` the
off-shell :math:`x`-:math:`A` momentum: at fixed :math:`E`, :math:`|p_s|`
and the spectator angle are the same variable, and a two-dimensional
:math:`(|p_s|, \theta_s)` quadrature would be degenerate. A ``ps`` window
(the cut of the data reduction) is then a cut on :math:`x`: the accepted
directions are those of the angle window whose :math:`q` lies in
:math:`[p_\mathrm{min}, p_\mathrm{max}]` (a ``ps`` window alone: every
direction whose :math:`q` lies in it). Its distribution (Hulthén, Gauss,
table) weights the plane-wave vertex; :math:`R(E)` and the DW vertex carry
:math:`|\tilde\varphi(q)|^2` of their own bound state instead.

*What is averaged.* At fixed :math:`E` the three-body phase space is
:math:`d\Omega_s` (the length of :math:`k_{sF}` is fixed), i.e.
:math:`d\cos\theta_\mathrm{cm}\,d\varphi`, :math:`\theta_\mathrm{cm}` the
polar angle of :math:`\mathbf k_{sF}` (the spectator in the :math:`a + A`
c.m.) to the beam. With an acceptance :math:`A(\theta)` the yield of the
energy bin is :math:`\int A\,|M(\theta)|^2\,d\cos\theta_\mathrm{cm}` (times
the HOES cross section and the kinematic factor, taken as constant across
the window as for ``ps``), and the PWA reduction divides it by the Monte
Carlo of :math:`|\phi_a(q)|^2` over the same events. Events of different
directions are different final states and add incoherently, so the
cross section is averaged, not the amplitude:

.. math::

   \bar R(E) \propto \frac{\int A\,|M|^2\,d\cos\theta_\mathrm{cm}}
                          {\int A\,|M_\mathrm{PW}|^2\,d\cos\theta_\mathrm{cm}}
   \quad(\texttt{dwpw}), \qquad
   \frac{\int A\,|M|^2\,d\cos\theta_\mathrm{cm}}{\int A\,d\cos\theta_\mathrm{cm}}
   \quad(\texttt{dw}),

normalized at :math:`E_\mathrm{ref}` as before; for the DW vertex, normalized
by :math:`4\pi\tilde\varphi(q)` of its own direction, the event weight is
:math:`A\,|\tilde\varphi(q)|^2\,d\cos\theta_\mathrm{cm}`:

.. math::

   \langle |M_l|^2 \rangle = \frac{\int A\,|\tilde\varphi(q)|^2\,c^\dagger G(\theta)\,c\,d\cos\theta_\mathrm{cm}}
                                  {\int A\,|\tilde\varphi(q)|^2\,d\cos\theta_\mathrm{cm}}
                           = c^\dagger \bar G\, c .

The HOES model is linear in :math:`G`: summed over the Cholesky components,
:math:`|X a_k + Y d_k|^2` gives :math:`|X|^2G_{11} + |Y|^2G_{22} + 2\,\mathrm{Re}
(X^*Y G_{12})` whatever the level sums :math:`X`, :math:`Y`, also with a
level-dependent boundary and with :math:`N_C`. So the window average is the
vertex of the averaged Gram matrix :math:`\bar G`, one node: the window costs
nothing per evaluation. Without distortion both averages reduce to the
plane-wave quantities (:math:`R \equiv 1`; the DW vertex with plane waves
is the plane-wave vertex at the :math:`p` of each direction).

*Relation to the* ``ps`` *window.* The measure is the one derived in
"Spectator-momentum window": :math:`d\cos\theta_\mathrm{cm} = q\,dq/(\beta
k_{sF}k_{aA})` at fixed :math:`E`. Until October 2026 the ``ps`` window
averaged instead with :math:`|\phi|^2 q^2\,dq` (an isotropic
:math:`\mathbf p_s`, the measure of events integrated over :math:`E`), the
plane-wave vertex on the whole :math:`[p_\mathrm{min}, p_\mathrm{max}]` and
the DW vertex on its reachable part; for 19F (below) the adopted-window THM
:math:`\chi^2` of the DW vertex was 176 with ``ps=hulthen:0-50`` that way
and 296 with ``spectatorAngles=cm:0-180`` and the same cut, and a ``ps``
table of weight :math:`|\phi_\mathrm{H}|^2 p` reproduced 302: the measure,
not :math:`\phi`, made the difference. Now a ``ps`` window alone is
``spectatorAngles=cm:0-180`` with its cut (``tests/thm_spectator_window``:
byte-identical with as many nodes), for every vertex model. For the
plane-wave vertex, ``spectatorAngles`` with a ``ps`` distribution averages
:math:`|M_l(p_{xA}(q))|^2` with :math:`A\,|\phi(q)|^2\,d\cos\theta_\mathrm{cm}`;
without a distribution (``ps=delta``) and without a computed distortion it
has nothing to average and is refused.

*Lab windows.* The spectator's lab angle follows from
:math:`\tan\theta_\mathrm{lab} = \sin\theta_\mathrm{cm}/(\cos\theta_\mathrm{cm} +
\gamma)`, :math:`\gamma = V_\mathrm{cm}/v_s` with :math:`v_s = \hbar
k_{sF}/m_s` the spectator's c.m. speed. For :math:`\gamma < 1` it is
monotonic, and :math:`[\theta_1, \theta_2]` maps to one c.m. interval
:math:`[f(\theta_1), f(\theta_2)]`, :math:`f(\theta) = \theta +
\arcsin(\gamma\sin\theta)`. For :math:`\gamma > 1` (a spectator slower in the
c.m. than the c.m. itself: 12C(14N,d) at 30 MeV above :math:`E = 2.2`
MeV, 2H(18O,α15N)n at 54 MeV above 0.63 MeV, 2H(19F,α16O)n at 55 MeV above
0.54 MeV) the lab angle reaches at most
:math:`\theta_m = \arcsin(1/\gamma)` and every lab angle below it has two c.m.
angles, the forward :math:`f(\theta)` and the backward :math:`\theta + \pi -
\arcsin(\gamma\sin\theta)`: the window maps to two intervals, both accepted
(a detector sees both; the backward one has the larger :math:`q` when the
Trojan horse is the beam, the smaller when it is the target). Each interval
gets ``spectatorAngleNodes`` Gauss-Legendre nodes in
:math:`\cos\theta_\mathrm{cm}`; whether a second interval exists is decided
at the top of the energy grid (an empty one has weight 0). A window of zero
width is its direction; on two branches the limit of a shrinking window,
weight :math:`|d\cos\theta_\mathrm{cm}/d\theta_\mathrm{lab}|` per branch
(:math:`(\gamma \pm 1)^2` at 0°). ``spectatorAngle=<lab>`` takes the
forward branch only.

*Syntax* (angles in degrees, polar angles of the spectator to the beam,
:math:`0 \le \theta_\mathrm{min} \le \theta_\mathrm{max} \le 180`):

``spectatorAngles=thmin-thmax``
   a uniform lab window, converted to c.m. intervals at every energy.
``spectatorAngles=cm:thmin-thmax``
   a uniform c.m. window (:math:`\theta_\mathrm{cm}` of :math:`\mathbf
   k_{sF}`; the quasi-free direction is 0° when the Trojan horse is the
   beam, 180° when it is the target).
``spectatorAngles=table:<file>``, ``spectatorAngles=cm:table:<file>``
   the acceptance :math:`A(\theta)` of a spectator emitted in that direction
   (the azimuth-integrated detection probability, any scale), in the lab or
   c.m. angle: two columns, angle (strictly increasing, 0-180) and
   :math:`A \ge 0`, ``#`` comments, at least two rows, some positive value;
   linear between rows; the window is its range; relative to the ``.azr``
   as ``ps`` tables. Not the angular distribution of accepted events (that
   contains the cross section and :math:`|\phi|^2` already). A table with
   kinks is integrated to :math:`O(h^2)` only.
``spectatorAngleNodes=N``
   Gauss-Legendre nodes per c.m. interval (1-64, default 8).

It needs ``distortion=coulomb`` or ``optical``, or a ``ps`` distribution
(plane-wave vertex). Refused, with ``ERROR: <thm> experiment[...]``: with
neither, together with ``spectatorAngle``, ``spectatorAngleNodes`` without a
window, ``psNodes`` with a window (its nodes are ``spectatorAngleNodes``),
malformed values, an unreadable or invalid table, and a data point at which
no direction of the window is accepted (for a lab window beyond
:math:`\theta_m`, or nothing inside the ``ps`` cut). Grid energies beyond the data without an accepted
direction take the nearest value. ``distortionRatio=dw`` with a ``ps``
window stays refused (:math:`|\phi|^2` twice).

*Numerics and cost.* :math:`R`: the radial integrals once per energy; only
:math:`P_l(x)` changes between directions, so the window costs
:math:`O(\text{nodes}\times l_\mathrm{max})` per grid energy on top of
them (19F, optical waves, one calculation of the adopted window: 14.2 s at the
quasi-free direction, 14.3 s with 8 nodes). The partial-wave sum stops when the bound on the rest is
below :math:`10^{-13}` of the smallest node amplitude. DW vertex: the Gram
matrix of every direction from the reduced amplitudes the vertex has
anyway, :math:`\bar G` per grid energy; tables :math:`n_E \times n_l \times
4` as without a window (plus the directions for the report), the same
evaluation time. Memory, 19F full data window (53 points, 30 keV folding;
CLI, getrusage): peak RSS 603 MB with :math:`R(E)`, ``ps=hulthen:0-50`` (8
nodes per point with ``spectatorAngles=cm:135-180``), 649 MB with the
``ps`` window alone (16 nodes; the plane-wave vertex tables dominate; 636 MB
without :math:`R`, 520 MB without a window), and 522 MB with
``vertexModel=dw``, the ``ps`` window and the angle window (523 MB without
the angle window). Eight nodes are
converged for the cases below: 16 nodes change no :math:`\chi^2` beyond
:math:`10^{-8}`.

*Output, pyazr.* The startup summary and the ``distortion:`` / ``vertex:``
lines of ``thm_experiments.out`` name the window, the nodes, the branches
and the cut; the ``distortion_point`` rows give
:math:`\langle|M|^2\rangle` and :math:`\langle|M_\mathrm{PW}|^2\rangle`
and the acceptance-weighted mean :math:`\theta_\mathrm{cm}`.
``AzrModel.set_thm_experiment(..., spectatorAngles="cm:135-180" | (7, 30),
spectatorAngleNodes=8)``; ``session.thm_distortion`` returns the averages;
``session.thm_vertex`` the averaged DW vertex (one node, ``dw_q`` its mean
:math:`q`) and the directions at the nearest grid energy
(``angle_theta_cm``, ``angle_q``, ``angle_weights``). In the GUI the
Experiments page shows the window under the momentum window ("Spectator
acceptance").

*Checks.* ``tests/reference/thm_spectator_angles_test`` (ctest
``thm_spectator_angles``): the weights of a uniform window sum to the
accepted solid angle against a brute-force sum over :math:`4\times 10^6`
bins (c.m. and lab windows, one and two branches, the whole sphere = 2 to
:math:`10^{-15}`, a ``ps`` cut, acceptance ramps in the c.m. and lab angle;
:math:`\le 5\times 10^{-6}`, the bins), every node inside the window, and
the zero-width branch weights against a shrinking window (:math:`2\times
10^{-6}`) and :math:`(\gamma-1)^2/(\gamma+1)^2` at 0°.
``tests/thm_spectator_angles/check.sh`` (``tests/18O_p_a_thm``):
``cm:25-25`` and the lab ``12-12`` (one branch, Ebeam 100 MeV) equal
``spectatorAngle=cm:25`` and ``12`` for :math:`R` and the DW model to all
printed digits, a :math:`10^{-4}`-degree window around 25° likewise; the
uniform ``cm:120-180`` equals Simpson's rule over 21 single-direction runs
with Richardson extrapolation -- :math:`\langle|M|^2\rangle`,
:math:`\langle|M_\mathrm{PW}|^2\rangle` to :math:`3\times 10^{-6}`, the
no-folding DW model weighted by :math:`(\kappa^2 + q^2)^{-2}` to
:math:`4\times 10^{-6}` (41 runs: :math:`6\times 10^{-8}` and
:math:`5\times 10^{-8}`); the lab ``0-180`` over two branches equals
``cm:0-180`` to the printed digits; with the DW vertex a ``ps`` window
covering every :math:`q` and a flat table change nothing; twelve refusals.
``tests/pyazr/thm_spectator_angles_test.py``: ``AzrModel``, CLI == session
for :math:`R` and the DW vertex, ``thm_distortion``, and the reported
:math:`\bar G` is the one the model uses (:math:`2\times 10^{-9}`). Without
the key every output file is byte-identical to the previous build (nine
configurations of ``tests/18O_p_a_thm``: no experiment, R(E) Coulomb and
optical at qf, lab and c.m. angles, DW with and without ``ps``, the
plane-wave ``ps`` window, ``theta``; at the time of that change -- the
``ps`` window itself moved to the fixed-:math:`E` measure in October 2026,
"Spectator-momentum window").

*The windows of two experiments.* 19F(p,αγ)16O, 2H(19F,α16O)n at 55 MeV
(Su et al., PRL 135 (2025) 182701, supplement secs. III-V): position-sensitive
detectors at 5.6-12.4° (A, B: 16O) and 15.7-31.6°, 16.5-32.4° (C, D: α) on
both sides of the beam, coincidences A-C and B-D, :math:`|p_s| < 50` MeV/c; a
Monte Carlo of that geometry (non-relativistic, isotropic :math:`F^*` decay,
out-of-plane acceptances of 0.7-1.6° assumed -- the result does not depend
on them) puts the neutron at :math:`\theta_\mathrm{cm} = 135`-180° (90 %:
139-175°, Hulthén-weighted 142-176°) at every :math:`E` = 0.1-0.9 MeV. The
cut alone requires :math:`\theta_\mathrm{cm} \ge 131.7°`
(:math:`m_n V_\mathrm{cm} = 67` MeV/c), so here the momentum cut, not the
detectors, sets the window: ``spectatorAngles=cm:135-180`` with
``ps=hulthen:0-50``. 12C(14N,α20Ne/p23Na)d at 30 MeV (Tumino et al., Nature
557 (2018) 687; setup in Tumino et al., Nuovo Cimento 42 C (2019) 55): the
spectator deuteron is detected, in coincidence with the α or p, by
telescopes on both sides of the beam covering 7-30°; :math:`|p_s|` up to
about 80 MeV/c from the phase space. That is ``spectatorAngles=7-30`` (lab)
directly, :math:`\theta_\mathrm{cm} \approx 12`-68° (90 %: 15-56°); its
backward branch (:math:`q > 100` MeV/c) is removed by the 80 MeV/c cut and
changes little anyway (DW :math:`\chi^2` 5670.0 with, 5669.9 without the
cut).

*Size* (models at the parameters of the examples, not refitted).
19F (``examples/f19_pag_thm`` before October 2026 -- THM :math:`\chi^2` 63.9
in the adopted window; the refitted example has 54.2 --, :math:`a_p` = 5.136 fm, ``vertex=constant``;
optical potentials of the DW section): the :math:`l` ratios of the DW vertex
and the THM :math:`\chi^2` of the adopted window (28 points) and of all 53
points:

====================================== ==================== ==================== =================== =========== ========
vertex / distortion                    :math:`|M_1(213)/`   :math:`|M_1(828)/`   :math:`|M_0(11)/`   χ² adopted  χ² full
                                       :math:`M_0(324)|^2`  :math:`M_0(324)|^2`  :math:`M_0(324)|^2`
====================================== ==================== ==================== =================== =========== ========
pw                                     3.00                 2.83                 1.55                63.9        3101
pw, R(E) optical, qf                                                                                 87.0        3106
pw, R(E) optical, cm:135-180                                                                         76.5        3106
pw, ps 0-50                            8.07                 7.04                 1.57                1114        3112
pw, R(E) optical, ps 0-50                                                                            1139        (a)
dw optical, qf                         0.41                 0.49                 1.08                611         2203
dw optical, cm:135-180                 0.95                 0.96                 1.04                281         3043
dw optical, ps 0-50                    0.94                 0.93                 1.07                296         3038
dw optical, cm:135-180 + ps cut        0.94                 0.93                 1.07                296         3038
dw optical, cm:140-176 + ps cut        0.93                 0.94                 1.04                315         3027
dw Coulomb, qf                         0.24                 0.33                 1.04                680         1861
dw Coulomb, ps 0-50                    0.52                 0.60                 1.03                536         2775
dw Coulomb, cm:135-180 + ps cut        0.52                 0.61                 1.03                536         2775
====================================== ==================== ==================== =================== =========== ========

(a) The profile of norm and linear background is degenerate there.
``ps 0-50`` alone is the cut on every direction, which here is
:math:`\theta_\mathrm{cm} \ge 131.7°`: the same directions as
``cm:135-180`` with the cut (the :math:`\chi^2` agree to :math:`10^{-9}`,
16 against 8 nodes); before October 2026 the ``ps`` rows used the isotropic
measure (dw optical 1.24, 1.18, 1.08 and :math:`\chi^2` 176 / 3082; dw
Coulomb 0.63, 0.73, 1.03 and 479 / 2895; pw 1525 adopted). With the
acceptance the :math:`l = 1` to :math:`l = 0` ratio of the DW vertex more
than doubles against the quasi-free direction (0.41 → 0.94 optical, 0.24 →
0.52 Coulomb); narrowing the window to the central 90 % changes it by 1 %.
With the plane-wave vertex the window raises it from 3.0 to 8.1 at 5.136 fm
and the adopted-window :math:`\chi^2` at the published parameters from 64 to
1114. For :math:`R(E)` (neutral spectator,
only the a + A wave distorted) the window moves :math:`\chi^2` from 87 to 77
(64 without :math:`R`). 12C+12C (``examples/c12c12_tumino2018``, the four THM
segments one experiment, 30 keV folding, point Coulomb): :math:`R` = 5.44 at
0.82 MeV to 0.0115 at 2.69 MeV at the quasi-free direction, 4.62 to 0.0215
over the lab window 7-30° (5.15 to 0.0120 over ``cm:12-68``): the window
halves the slope's range (a factor 215 instead of 470) and the THM
:math:`\chi^2` goes 3247 → 2907 (``cm:12-68``: 3133; 61 without :math:`R`).
The DW vertex ratios :math:`|M_l/M_0|^2` at :math:`E` = 1.0, 1.5, 2.0, 2.5
MeV:

=============================== ======================= ======================= =======================
vertex                          :math:`l = 2`           :math:`l = 4`           :math:`l = 6`
=============================== ======================= ======================= =======================
pw                              0.86 0.73 0.55 0.25     0.36 0.14 0.0005 0.58   0.05 0.33 1.31 6.0
dw, qf                          1.01 1.18 1.40 1.68     0.98 1.09 1.24 1.47     3.4 4.0 4.6 5.4
dw, lab 7-30 + cut 80 MeV/c     0.93 1.01 1.11 1.19     0.62 0.62 0.61 0.58     2.9 3.0 3.0 2.8
dw, cm:12-68                    0.83 0.91 1.02 1.15     0.41 0.43 0.46 0.53     2.3 2.4 2.5 2.6
=============================== ======================= ======================= =======================

The acceptance lowers :math:`l = 4` against :math:`l = 0` by a factor 1.6-2.5
and flattens its energy dependence; the THM :math:`\chi^2` with the DW
vertex goes 4713 (qf) → 5670 (lab 7-30) and 6625 (``cm:12-68``), all far
from the plane-wave 61 at the published parameters: the DW vertex needs
refitted strengths, and the direction window is part of its model
dependence, of the size of the choice between Coulomb and optical waves.

*Scope and limits.* The acceptance in the spectator's polar angle (the
vertex does not depend on its azimuth; the table is azimuth-integrated); the
:math:`c`-:math:`C` direction integrated as in the angle-integrated HOES
observable (an acceptance that correlates the two is beyond it); the
kinematic factor constant across the window; one acceptance for every
energy (a measured :math:`A(\theta, E)` would need a table per energy bin).

Fixed-angle observable
----------------------

A THM measurement rarely covers the full solid angle of the exit pair. The
three-body yield is binned in the c.m. angle of the two detected particles,
:math:`\theta_\mathrm{cm} = \arccos(\hat k_{xA}\cdot\hat k_{bB})` (Tribble et
al., RPP 77 (2014) 106901, section 4.3), and either integrated over the whole
range after an angular-distribution fit (La Cognata et al., arXiv:0909.4716,
18O(p,α)15N) or kept in a window: :math:`\theta_\mathrm{cm}` = 50–70° for
7Li(p,α)α (Tumino et al., EPJA 27 s01 (2006) 243), 95–110° for
6Li(n,α)3H (Gulino et al. 2010, Tribble Fig. 32). A window observable is
:math:`d\sigma/d\Omega` averaged over the window, not :math:`\sigma/4\pi`:
at a fixed angle the entrance partial waves of one channel spin and the
:math:`J^\pi` groups interfere, and only the :math:`4\pi` integral removes
the cross terms. ``theta=<thmin>-<thmax>`` on an experiment line computes it
for every segment of the experiment (``ThmAngular.h``,
``THMMatrixFunc::CalculateTHMCrossSection``).

*Angle.* :math:`\theta` is the c.m. angle of particle 1 of the segment's exit
pair relative to particle 2, measured from the relative momentum of particle
1 of the entrance pair relative to particle 2 -- the convention of AZURE2's
ordinary differential segments (``EPoint::ConvertLabAngle``). When the
transferred particle :math:`x` is particle 1 of the entrance pair and the
angle of the papers refers to particle 1 of the exit pair, it is their
:math:`\theta_\mathrm{cm}`; otherwise give the supplementary window
:math:`180 - \theta_\mathrm{max}` to :math:`180 - \theta_\mathrm{min}`. The
direction of :math:`\vec k_{xA}` is that of the off-shell :math:`x`
(:math:`\vec p_x = \vec p_a - \vec p_s` in the quasi-free picture), which
the experiments reconstruct event by event. Degrees,
:math:`0 \le \theta_\mathrm{min} \le \theta_\mathrm{max} \le 180`;
:math:`\theta_\mathrm{min} = \theta_\mathrm{max}` is a single angle. ``theta=all``
(the default) is the angle-integrated observable of the previous sections,
bit for bit.

*Amplitude.* With the quantization axis along :math:`\hat p_{xA}` the plane
wave of :math:`x` contributes only :math:`m_l = 0`, with
:math:`Y_{l0}(\hat z)\propto\sqrt{2l+1}` -- exactly the structure of an
ordinary unpolarized beam along :math:`z`, whose Coulomb or plane wave also
has :math:`m_l = 0` only. The HOES amplitude for entrance channel spin
projection :math:`\nu` and exit :math:`\nu'` is therefore the ordinary
reaction amplitude (Lane & Thomas, RMP 30 (1958) 257, section VIII) with the
T-matrix element of each :math:`J^\pi` group replaced by the HOES partial
amplitude:

.. math::

   F_{\nu\nu'}(\theta) = \sum_{J l l'} \sqrt{2l+1}\,
   \langle s\,\nu\,l\,0|J\,\nu\rangle
   \langle s'\,\nu'\,l'\,m'|J\,\nu\rangle\,
   x^J_{(s'l'),(sl)}\, Y_{l'}^{m'}(\theta, 0),
   \qquad m' = \nu - \nu',

.. math::

   x^J_{(s'l'),(sl)} = \sqrt{K(E)\,2P_{c'}}\; e^{i(\omega_{c'} - \phi_{c'})}
   \sum_{\lambda\lambda'} \gamma_{\lambda c'} N_\lambda A_{\lambda\lambda'}
   V_{\lambda'}^{(s,l)} ,

:math:`c' = (s', l')` an exit channel of the group,
:math:`V^{(s,l)}` the vertex of the first section, :math:`N_\lambda` the line
shape (1 without ``lineshape=on``), :math:`e^{i\omega}` and
:math:`e^{-i\phi}` the Coulomb and hard-sphere phases of the exit channel as
in AZURE2's T matrix (``AMatrixFunc::CalculateTMatrix``). The ordinary
entrance factor :math:`e^{i(\omega_c - \phi_c)}\sqrt{P_c}` -- the incident
Coulomb wave's :math:`e^{i\sigma_l}` and the surface value of :math:`F_l` --
is replaced by the vertex :math:`M_l` of the plane wave, which has neither
(both waves carry the same :math:`i^l` of the partial-wave expansion, which
AZURE2's phase convention absorbs). This is the structure of Tribble et al.
eqs. (2.75) and (2.79): the on-shell S-matrix element with its entrance
factor :math:`P_l^{-1/2} e^{i\delta^\mathrm{hs}_l}` taken out and the
plane-wave vertex put in, times the spherical harmonics of
:math:`\hat k_{bB}` and :math:`\hat p_{xA}`. The exit phases matter only here: the
angle-integrated sum is incoherent in the exit channels.

*Legendre form.* Summing :math:`|F|^2` over the spin projections with the
Racah algebra gives the Blatt-Biedenharn form (RMP 24 (1952) 258) that
AZURE2 uses for ordinary angular distributions (``CNuc::CalcAngularDists``,
the same :math:`\bar Z` and the same :math:`(-1)^{s'-s}/4`):

.. math::

   \sum_{\nu\nu'} |F_{\nu\nu'}(\theta)|^2 = \frac{1}{\pi}\sum_L b_L P_L(\cos\theta),
   \qquad
   b_L = \sum_{ij} \frac{(-1)^{s'-s}}{4}\,
   \bar Z(l_i J_i l_j J_j; s L)\,\bar Z(l'_i J_i l'_j J_j; s' L)\,
   \mathrm{Re}\,x_i x_j^*,

.. math::

   \bar Z(l_1 J_1 l_2 J_2; s L) = \sqrt{(2l_1+1)(2l_2+1)(2J_1+1)(2J_2+1)}\,
   \langle l_1\,0\,l_2\,0|L\,0\rangle\, W(l_1 J_1 l_2 J_2; s L),

the double sum over the partial amplitudes of all :math:`J^\pi` groups with
the same entrance channel spin :math:`s` and the same exit :math:`s'` (the
spins are unpolarized and unobserved: different channel spins add
incoherently). The Decay/KGroup objects of ordinary segments are not reused
-- the THM segment has no T matrix -- but the algebra is theirs, and the
direct sum above equals it (checked to :math:`10^{-16}`, below). Since
:math:`4 b_0 = \sum_J (2J+1)\sum |x|^2`, the :math:`4\pi` integral of the
observable is the angle-integrated HOES cross section of the first section
(the :math:`1/(2l+1)` argument of ``entranceL`` is the :math:`L = 0` term).
The model of a point is

.. math::

   \Bigl\langle \frac{d\sigma}{d\Omega} \Bigr\rangle
   = \frac{1}{\pi}\sum_L b_L\,
     \frac{\int_{\cos\theta_\mathrm{max}}^{\cos\theta_\mathrm{min}} P_L(u)\,du}
          {\cos\theta_\mathrm{min} - \cos\theta_\mathrm{max}},

the average over the window's solid angle; 0–180 gives
:math:`\sigma_\mathrm{HOES}/4\pi` (with the exact :math:`\pi`; AZURE2's
constant ``pi`` is 3.14159265). The averages of :math:`P_L` are taken once
at startup by 48-point Gauss-Legendre in :math:`\cos\theta` (nodes by
Newton's method), exact for :math:`L \le 95`; orders up to 80 are kept, so a
model with an exit :math:`l' > 40` is refused. The geometric coefficients
depend only on the list of partial waves and are cached per thread; the cost
is not measurable (7Li(p,α): 2.9 s either way).

*Limit* :math:`\theta \to 0`. Along the axis :math:`Y_{l'}^{m'}(0) =
\delta_{m'0}\sqrt{(2l'+1)/4\pi}`, so :math:`\nu' = \nu` and

.. math::

   \frac{d\sigma}{d\Omega}(0) = \frac{1}{4\pi}\sum_{s s'\nu}
   \Bigl|\sum_{J l l'}\sqrt{(2l+1)(2l'+1)}\,
   \langle s\,\nu\,l\,0|J\,\nu\rangle\langle s'\,\nu\,l'\,0|J\,\nu\rangle\,
   x^J_{(s'l'),(sl)}\Bigr|^2 ,

coherent in :math:`l` *and* :math:`J`, with Clebsch-Gordan weights. For an
:math:`\alpha + \alpha` exit (:math:`s' = 0`, :math:`l' = J`) only
:math:`\nu = 0` survives, and for 7Li + p the channel spin :math:`s = 2`
drops out entirely at 0° (:math:`\langle 2\,0\,l\,0|2\,0\rangle = 0` for odd
:math:`l`). This is not ``entranceL=coherent`` (equal weights, :math:`J^\pi`
groups incoherent), which is therefore neither the angle-integrated nor the
fixed-angle observable; ``tests/7Li_p_a``: :math:`\chi^2` 1747.26 at
``theta=0-0``, 3196.77 with ``entranceL=coherent``, 2138.48
angle-integrated (same parameters, profiled norm).

*Identical particles.* A channel of two identical nuclei is admitted only
with :math:`l + s` even (``CNuc``), so the exit waves of an
:math:`\alpha + \alpha` pair all have even :math:`l'` and the distribution
has even :math:`L` only: it is symmetric about 90°, and a window and its
mirror give the same model. No factor is applied, as in AZURE2's
differential branch for an identical exit pair and in the angle-integrated
THM observable: the model is per event, its :math:`4\pi` integral is the
angle-integrated one, and a detector that counts either particle sees twice
it -- a constant, absorbed by the norm. An identical entrance pair (12C+12C)
also gives even :math:`L` only (both amplitudes of a :math:`b_L` term have
the same :math:`s`, so :math:`l_i + l_j` is even).

*Everything else combines unchanged.* A spectator-momentum window averages
the observable over its nodes with their weights (each node its own set of
:math:`x`); the line shape sits inside :math:`x`; ``weight[k]`` and the
distortion factor multiply the result before the resolution folding; the
folding, the shared norm and the background of the experiment act on the
model as before. The distortion factor is computed at one spectator angle and
does not depend on :math:`\theta`. ``entranceL=coherent`` with a window is
refused (the window computes the :math:`l` interference exactly).

Errors (``ERROR: <thm> experiment[<name>]: ...``, exit non-zero): a value
that is not ``all`` or two numbers ``a-b`` with
:math:`0 \le a \le b \le 180` (a single number, a reversed or negative window,
beyond 180°), ``theta`` twice, a window together with
``entranceL=coherent``. ``pyazr``: ``AzrModel.set_thm_experiment(...,
theta="50-70")`` (or ``(50, 70)``, ``"all"``) with the same rules; the
session's residuals, Jacobian and output files follow. The GUI keeps the key
as written.

*Validation.* ``tests/reference/thm_fixed_angle_test`` (ctest
``thm_fixed_angle``) takes the partial amplitudes of a toy two-level model
(7Li+p→α+α-like: :math:`J = 2^+` with two levels and :math:`(s,l)` = (1,1),
(1,3), (2,1), (2,3), plus :math:`0^+`; and a spin-1/2 case with
:math:`1/2^\pm`, :math:`3/2^-` and odd :math:`L`) from
``thm_fixed_angle_reference.py``, which does the literal M-sum above with
sympy's exact Clebsch-Gordan coefficients, mpmath's spherical harmonics and
quadrature: :math:`b_L`, :math:`d\sigma/d\Omega` at 0, 37, 90, 143, 180° and
five windows agree to :math:`5\times 10^{-16}`; 0–180 equals
:math:`\sum(2J+1)|x|^2/4\pi` to :math:`3\times 10^{-16}`;
:math:`\theta = 0` equals the axis sum (:math:`\nu' = \nu`) and windows
:math:`0`–:math:`t` approach it as :math:`t^2`; the identical-boson case has
no odd :math:`b_L` and is symmetric, the spin-1/2 case is not.
``tests/thm_fixed_angle/check.sh`` (CLI): ``theta=all`` byte-identical to no
key (and the model to no block); 0–180 times :math:`4\pi` equals the
angle-integrated model to :math:`4\times 10^{-11}` (the output file's 11
digits) for 7Li(p,α) and 17O(n,α), with the same :math:`\chi^2`; the
one-group 18O(p,α) model (:math:`1/2^+`, isotropic) gives the same shape for
any window, also with a ``ps`` window and a linear background; windows
0–0.5° and 0–1° approach 0–0 as :math:`t^2` (ratio 3.998–4.003);
30–60° equals 120–150° for 7Li(p,α)α and 10–40° equals 140–170° for
6Li(d,α)α, while 17O(n,α)14C (1\ :sup:`-`, 2\ :sup:`+`, 3\ :sup:`-`,
5\ :sup:`-`) differs by up to 27 % between 20–60° and 120–160°; refusals.
``tests/pyazr/thm_fixed_angle_test.py``: ``AzrModel``, CLI == session,
:math:`4\pi` model(0–180) == model to :math:`2\times 10^{-15}` in double
precision, residuals and Jacobian unchanged, symmetry.

*Size* (``tests/7Li_p_a``, Paneru parameters, not refitted, 30 keV folding,
profiled norm). The Tumino 2006 window ``theta=50-70`` changes the shape of
the THM model by 14 % rms and up to 45 % over the 66 points (both scaled to
the same mean): relative to the angle-integrated curve it is 16 % higher at
0.08 MeV, 14 % lower across the 2.6 MeV peak (0.86 at 2.1–2.7 MeV) and 25 %
higher across the 5 MeV one (1.25–1.27 at 4.3–5.1 MeV); :math:`\chi^2`
2138.48 → 2349.86. (``theta=30-60``: 1765.64; ``theta=0-0``: 1747.26.)
The data are penetrability-corrected and normalized to direct data (see
``tests/7Li_p_a/README.md``), so this is the size of the effect, not a
statement about the fit. 18O(p,α)15N (``tests/18O_p_a_thm``): the data are
integrated over the whole angular range (La Cognata et al.,
arXiv:0909.4716), and with one :math:`1/2^+` group the model is isotropic
anyway -- no change for any window.

Coulomb effects: what each option contains
------------------------------------------

Four options describe Coulomb physics of the same three-body reaction
:math:`a(x+s) + A \to s + F^* \to s + b + B` in factorised approximations:
the Coulomb term :math:`C_l` (``coulombIntegral=1``), the line shape
:math:`N_C` (``lineshape=on``), the distortion factor :math:`R(E)`
(``distortion=``) and the DW vertex (``vertexModel=dw``). Which combinations
count an interaction twice follows from where each one sits in the amplitude
of Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020) 233 (arXiv:2007.13331),
Appendix A. The resonant part of the Green's function of the
:math:`s + F` system is decomposed over the :math:`s`-:math:`F` Coulomb
scattering states (eqs. 129, 136),

.. math::

   M \propto \int \frac{d\mathbf k_{sF}}{(2\pi)^3}\,
   \frac{\langle \Phi^{C(-)}_{f} | \Psi^{C(-)}_{\mathbf k_{sF}} \rangle\;
         M^{(tr)}(\mathbf k_{sF}, \mathbf k_{aA})}{k_R^2 - k_{sF}^2},
   \qquad
   M^{(tr)} = \bigl\langle \Psi^{C(-)}_{\mathbf k_{sF}} \tilde\varphi_{R(xA)}
   \bigm| \Delta V \bigm| \varphi_a \Psi^{C(+)}_{\mathbf k_{aA}} \bigr\rangle

(eqs. 138-140). The intermediate :math:`s`-:math:`F` Coulomb state appears
twice, once on each side of the completeness relation: as the bra of the
*formation* amplitude :math:`M^{(tr)}` (eqs. 23-42), and as the ket whose
overlap with the final three-body Coulomb state :math:`\Phi^{C(-)}_f`
(:math:`s`-:math:`b` and :math:`s`-:math:`B` after :math:`F^*` has decayed)
gives the factor :math:`N_C(E_0 - E - i\Gamma/2)^{-1-i\zeta}` of eqs. 55-57.
The result is the product :math:`M^{(tr)} N_C/(E_0 - E - i\Gamma/2)`
(eq. 55); Mukhamedzhanov, EPJA 58 (2022) 71, eq. 18 with eqs. 24-26
multiplies :math:`|N_C|^2` by the zero-range DWBA cross section in the same
way, with the :math:`s + F` wave at :math:`k_{sF} = \sqrt{2\mu_{sF}E_{sF}}`
of the event. In AZURE2:

.. list-table::
   :header-rows: 1
   :widths: 16 34 34 16

   * - option
     - interaction, region
     - where in the amplitude
     - :math:`l`, level
   * - vertex :math:`M_l` (``pw``)
     - :math:`x`-:math:`A` at the channel radius; plane :math:`a + A`,
       :math:`s + F`
     - surface term of :math:`M^{(tr)}` with plane waves (2020 eqs. 37-41)
     - per :math:`l`; common to the levels
   * - ``coulombIntegral=1``
     - :math:`x`-:math:`A` Coulomb outside :math:`a`, plane :math:`a + A`
     - external prior term of :math:`M^{(tr)}` with :math:`U_{aA} = U_{sA} =
       0` (eq. 33; Tribble eq. 2.79)
     - per :math:`l`
   * - ``distortion`` (:math:`R`)
     - :math:`a + A` (Coulomb of :math:`Z_a = Z_x + Z_s`, optional nuclear)
       and :math:`s + F` (:math:`Z_sZ_F`) in the formation region
     - :math:`M^{(tr)}` in zero range, :math:`r_{xA} \to 0` (eqs. 35-39),
       at the running :math:`E_{sF}`
     - common to all :math:`l` and levels
   * - ``vertexModel=dw``
     - the same distortions, finite range in :math:`s`-:math:`x`, at the
       channel radius
     - surface term of :math:`M^{(tr)}` (eq. 32) at the running
       :math:`E_{sF}`; no external term
     - per :math:`l`; common to the levels
   * - ``lineshape=on`` (:math:`N_C`)
     - :math:`s`-:math:`F^*` beyond the formation region while :math:`F^*`
       lives; :math:`s`-:math:`b`, :math:`s`-:math:`B` after its decay
       (post-collision interaction)
     - the ket side: overlap of the intermediate with the final three-body
       Coulomb state (eqs. 55-62)
     - per level (its pole)

The window ``ps`` (spectator momentum) and ``theta`` (exit angle) carry no
Coulomb interaction; ``weight[k]`` and ``distortion=table`` are whatever the
table holds.

*(a) DW vertex and* :math:`N_C` *do not overlap.* The DW vertex is
:math:`M^{(tr)}` (its surface term); :math:`N_C` is the other side of the
completeness relation. Two limits separate them. When the decay does not
change the Coulomb field the spectator sees (:math:`Z_b = 0`, :math:`m_B
\gg m_b`: :math:`\eta_{sB} = \eta_0`, :math:`\zeta = 0`), the final state
is the intermediate :math:`s`-:math:`F` Coulomb wave,
:math:`\langle\Phi^{C(-)}_{\mathbf p}|\Psi^{C(-)}_{\mathbf k}\rangle =
(2\pi)^3\delta(\mathbf p - \mathbf k)`, the integral collapses to
:math:`M^{(tr)}(\mathbf p)/(k_R^2 - p^2)` -- the vertex with its full
:math:`s + F` Coulomb distortion at the running energy, as ``dw`` and
:math:`R` evaluate it -- and :math:`N_C = 1`. Conversely, with the
:math:`s + F` wave of the vertex switched off (``opticalSF=plane``)
:math:`N_C` is unchanged. The :math:`-\eta_0` in :math:`\zeta = \eta_{sb} +
\eta_{sB} - \eta_0` is therefore not a subtraction of the vertex's
:math:`s`-:math:`F` distortion: it is the long-range Coulomb phase of the
intermediate state, which ends when :math:`F^*` decays, replaced by that of
the products. Both factors use the same :math:`k_{sF}(E)` (``lineshape``
takes :math:`\eta_0` at the point's :math:`k_{sF}`, the DW vertex and
:math:`R` their :math:`s + F` wave there), so they are consistent with each
other. What is left over is a cross term: with :math:`\zeta \ne 0` the
integral samples :math:`M^{(tr)}` between the pole and the final momentum,
an error of order :math:`\delta E\, \partial_E \ln|M^{(tr)}|^2`, with
:math:`\delta E` the post-collision energy transfer. 2020 eq. 55 takes
:math:`M^{(tr)}` at the pole :math:`k_{0(sF)}` instead, as does the
per-pole formation amplitude of Lei (arXiv:2605.16890); AZURE2's running
energy is exact at :math:`\zeta = 0` and agrees with the pole value at each
peak of a narrow level. For 12C(14N,d) the peaks move by up to 15 keV and
:math:`\partial_E \ln R \approx 3.2` MeV\ :sup:`-1`, so the cross term is at
most about 5 % of the vertex at a peak -- not implemented. ``lineshape=on``
with ``vertexModel=dw`` (or with :math:`R`) is allowed.

*(b) The external terms.* In the plane-wave limit the external prior term
of the DW vertex is :math:`C_l` ("Distorted-wave entrance vertex", plane-wave
limit). With a distorted :math:`a + A` wave the operator outside the
radius is :math:`V^C_{xA} + V^C_{sA} - U^C_{aA} = e^2 Z_A(Z_x/r_{xA} +
Z_s/r_{sA} - Z_a/r_{aA})` (2020 eq. 26), whose monopole cancels; with
:math:`\mathbf r_{sA} = \mathbf r + \mathbf u`, :math:`\mathbf r_{aA} =
\mathbf r + \beta\mathbf u`, what remains is the dipole
:math:`-e^2Z_A(Z_s m_x - Z_x m_s)/m_a\; \hat r\cdot\mathbf u/r^2`, and in
the zero range of :math:`R` (:math:`\mathbf u \to 0`) nothing. The
:math:`x`-:math:`A` Coulomb force that :math:`C_l` adds is the :math:`x`
part of the :math:`a`-:math:`A` force that generated the :math:`a + A`
Coulomb wave; adding :math:`C_l` to a model built on that wave counts it
twice. Hence ``coulombIntegral=1`` is refused with ``vertexModel=dw`` (as
before) *and* with :math:`R(E)` whenever its :math:`a + A` wave is
distorted (``distortion=coulomb``, or ``optical`` with ``opticalAA``
``coulomb`` or a potential); with ``opticalAA=plane`` :math:`C_l` is the
external term of that plane wave and is allowed. ``distortion=table``
with ``coulombIntegral=1`` is warned (the table's content is not known).
The zero range of :math:`R` does not carry the *energy dependence* of
:math:`C_l` (the :math:`x`-:math:`A` outgoing wave outside :math:`a` at
the energy :math:`E`): the two are different approximations of the same
force, of which one may be chosen. For 12C(14N,d) the dipole vanishes too
(:math:`Z_s m_x = Z_x m_s`, both :math:`N = Z`).

*(c)* :math:`R(E)` *and* :math:`N_C` *do not overlap*, for the reason of
(a): :math:`R` is the zero-range, :math:`l`-summed limit of the same
:math:`M^{(tr)}`. 2022 eqs. 18 and 40 contain both. Without folding the
model with both is the product of the two factors exactly
(``tests/thm_coulomb_consistency``: :math:`7\times 10^{-11}`).

*(d) Other combinations.* :math:`R` and ``dw`` are the same amplitude and
are never applied together (``vertexModel=dw`` does not apply :math:`R`).
``distortionRatio=dw`` makes :math:`R = |M|^2`, which carries the momentum
distribution :math:`|\tilde\varphi(q)|^2` at the spectator direction; a
``ps`` window already averages the model with the event weight
:math:`|\varphi(p_s)|^2\,d\cos\theta_\mathrm{cm}` of data divided by
:math:`|\varphi|^2`, so the pair is refused (``dwpw`` divides it out and is
allowed). With a ``ps`` window :math:`R` is averaged over the same accepted
directions as the plane-wave vertex ("Experimental acceptance"); the model
is then the product of the two averages rather than the average of the
product, an approximation of the size of the variation of :math:`R` across
the window (12C(14N,d): :math:`R` at 20° in the c.m. differs from the
quasi-free one by 0.993-1.003 over the data; ``vertexModel=dw`` evaluates
the distortion in every direction). ``weight[k]`` on a segment with :math:`R` or the
DW vertex is warned, as before. ``theta`` does not touch the Coulomb
factors.

*Allowed combinations.*

.. list-table::
   :header-rows: 1
   :widths: 22 18 20 22 18

   * -
     - ``coulombIntegral=1``
     - ``lineshape=on``
     - ``distortion`` (:math:`R`)
     - ``vertexModel=dw``
   * - ``coulombIntegral=1``
     - --
     - allowed
     - refused (``coulomb``; ``optical`` unless ``opticalAA=plane``);
       ``table`` warned
     - refused
   * - ``lineshape=on``
     - allowed
     - --
     - allowed
     - allowed
   * - ``distortion`` (:math:`R`)
     - see above
     - allowed
     - --
     - :math:`R` not applied
   * - ``ps`` window
     - allowed
     - allowed
     - ``distortionRatio=dw`` refused; ``dwpw`` allowed
     - allowed (per node)

The rules are in ``CheckThmExperiments`` and ``CheckThmCoulombConsistency``
(``src/ThmExperiment.cpp``), called by ``Config::ReadThmBlock`` and by the
GUI's THM workspace; ``pyazr.AzrModel`` applies them in
``set_thm_experiment`` and ``set_thm_option``. Files without these
combinations give results identical to before, byte for byte.

*Size* (models at the published parameters, not refitted; THM :math:`\chi^2`
with the shared norm profiled). 12C(14N,α/p)d at 30 MeV
(``examples/c12c12_tumino2018``, the four THM segments one experiment, 30 keV
folding; plane-wave vertex :math:`\chi^2` 61.0):

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - model
     - THM :math:`\chi^2`
     - factor on the model (per channel)
   * - plane-wave vertex
     - 61.0
     - --
   * - :math:`C_l`
     - 66.1
     - 0.29-0.42 (shape 1.09-1.45)
   * - :math:`N_C`
     - 479
     - shape 2.0-6.9
   * - :math:`R` (coulomb, qf, dwpw)
     - 3247
     - 0.012-5.4 (shape 181-440)
   * - :math:`R` + :math:`N_C`
     - 2861
     - :math:`N_C` shape 2.0-6.7 on :math:`R`
   * - :math:`R` + :math:`C_l` (now refused)
     - 3224
     - :math:`C_l` again 0.29-0.42 (shape 1.09-1.42)
   * - dw (Coulomb)
     - 4713
     - --
   * - dw + :math:`N_C`
     - 5120
     - :math:`N_C` shape 1.9-8.2 on dw
   * - :math:`C_l` + :math:`N_C`
     - 529
     - --
   * - :math:`R` + window (Eckart, 0-40 MeV/c)
     - 5746
     - --
   * - same, ``distortionRatio=dw`` (now refused)
     - 5647
     - :math:`|\tilde\varphi(q)|^2` again, 0.952-1.020
   * - :math:`R` at ``spectatorAngle=cm:20``
     - 3231
     - 0.993-1.003 against qf

The double-counted pieces that are now refused are small next to :math:`R`
here: :math:`C_l` changes the shape by up to 45 % and the second
:math:`|\tilde\varphi|^2` by 7 %. :math:`N_C` multiplies the :math:`R` model
as it does the plane-wave one (the shape factors differ by the 30 keV folding
of the steep :math:`R`); on the DW vertex its effect differs (1.9-8.2)
because the DW vertex changes the relative weights of the entrance partial
waves in the coherent level sum, not because of an overlap. 19F(p,αγ)16O
from 19F(d,n) at 55 MeV (``examples/f19_pag_thm`` as it was before October
2026, THM :math:`\chi^2` 63.9; the current example, with the published JUNA
table and the 0.0092 cap on the hidden THM bars, has 54.2):
the spectator is a neutron, :math:`N_C = 1` (byte-identical);
:math:`C_l` (:math:`p` + 19F, :math:`Z_xZ_A = 9`, deep below the barrier)
changes the model by 0.26-1.0 (:math:`\chi^2` 133); :math:`R` (only the
:math:`d` + 19F wave is distorted) by 0.94-1.16 (:math:`\chi^2` 76.5); both
together, now refused, 204, the :math:`C_l` factor 0.28-1.0 on :math:`R`:
for this deep sub-barrier proton the double count would have been the
largest effect of the model.

Normalization, gradients and uncertainty bands
----------------------------------------------

A THM segment with a free norm has its arbitrary scale profiled out
(``ESegment::IsProfiledNorm``): :math:`n^* = S_{mm}/S_{md}` and
:math:`\chi^2 = S_{dd} - S_{md}^2/S_{mm}`, no penalty, no fit parameter (the
segments of a THM experiment share one, with an optional background: see
"THM experiments"). The
analytic adjoint differentiates the T-matrix observable, not the HOES one, so
every derivative of a THM point is taken by central differences of the HOES
model, :math:`J_m`, with the dependence of the profiled scale
:math:`s = 1/n^*` added analytically (``ComputeTHMRows`` in ``AZUREGrad.h``;
used by the MIGRAD gradient, ``pyazr``'s Jacobian and gradients, and the CLI
band).

The cross-section band (``--covariance-band``) of a THM point is that of the
model *as it lies against the data*. The output file shows the model
:math:`m` next to the data scaled by :math:`n^*`, and the parameters move
both, so the band is taken for :math:`q(p) = m(p)\,n^*(p_0)/n^*(p)`:

.. math::

   \delta q = \sqrt{g\,\Sigma\,g^T}, \qquad
   g = J_m + m\,\frac{\partial s/\partial p}{s}

(:math:`\partial s/\partial p = 0` for a fixed THM norm). A change of the
parameters that only rescales :math:`m` is absorbed by :math:`n^*` and gives
no band, as it should: the THM scale is arbitrary.
``tests/thm_band/check.sh`` checks the band against finite differences of the
CLI's own output on ``tests/18O_p_a_thm`` with six free parameters.

Model averaging
---------------

The THM options above are not parameters the data determine: the channel
radius, ``vertexModel`` (pw or dw), ``vertex``, the optical potentials and the
``ps`` window each define a different model, and a strength or an S factor
read off a THM fit moves with them (the 19F :math:`l` ratios by factors 1.1-18
in "Distorted-wave entrance vertex"). The model-dependence protocol fits
every variant on its own and quotes the spread. ``pyazr.modelavg`` turns the
fitted variants into one number per quantity; ``scripts/thm_model_average.py``
produces the variants (one project per variant through ``AzrModel``, fitted
one after the other with scipy's ``least_squares`` on ``residuals`` and
``residual_jacobian`` plus the norm and shift penalty rows; pyazr itself does
no fitting).

*Weights.* For variant :math:`i` with minimised :math:`\chi^2_i`, :math:`N`
points and :math:`k_i` free parameters (the THM norms and backgrounds that
are profiled count: they are eliminated in closed form, but they are fitted),

.. math::

   w_i = \frac{\pi_i\, e^{-\Delta_i/2}}{\sum_j \pi_j\, e^{-\Delta_j/2}}, \qquad
   \Delta_i = \mathrm{IC}_i - \min_j \mathrm{IC}_j, \qquad
   \mathrm{IC} = \chi^2/s + 2k \;\;\text{(Akaike)}

with the alternatives :math:`\chi^2/s + k\ln N` (BIC), :math:`\chi^2/s`
alone, and flat weights; :math:`\pi_i` is an optional prior weight (e.g. 0.5 for
the extreme radii of a scan). Akaike weights estimate the relative expected
Kullback-Leibler distance of each model from the truth and need no nesting:
a 4.1 fm and a 6.1 fm model, or a plane- and a distorted-wave vertex, are not
special cases of one another, so a likelihood-ratio test does not apply,
while :math:`\Delta\mathrm{AIC}` does. :math:`\Delta = 2` gives a weight ratio
0.37, :math:`\Delta = 10` gives 0.007. Variations of the radius, the vertex
or the potentials do not change :math:`k`, so for such a grid Akaike, BIC and
:math:`\chi^2`-only weights are identical; they differ only when a variant adds
parameters (a background, a coherent background term, a free energy shift).
The weights compare models of the same data: variants with different
:math:`N` (another energy window) are not comparable this way, and the module
warns.

*Spreads.* For a quantity with value :math:`m_i` and statistical variance
:math:`s_i^2` in variant :math:`i`,

.. math::

   \bar m = \sum_i w_i m_i, \qquad
   S^2 = \sum_i w_i s_i^2, \qquad
   M^2 = \sum_i w_i (m_i - \bar m)^2, \qquad
   T^2 = S^2 + M^2 .

:math:`S` is the statistical error of a typical variant (the weighted mean of
the variances, not reduced by the number of variants: they all fit the same
data); :math:`M` is the model spread, the scatter of the best-fit values
between models at their weights; :math:`T` is the total of the law of total
variance. They are kept separate in the output because they answer different
questions: :math:`S` shrinks with more or better data, :math:`M` does not,
and it is :math:`M` that says how much the reaction theory still matters. The
averaged covariance matrices are :math:`\sum_i w_i C_i` and
:math:`\sum_i w_i (m_i-\bar m)(m_i-\bar m)^T`. The driver takes
:math:`C_i = D (J^TJ)^{-1} D^T` at each variant's end point, with :math:`J` the
residual Jacobian (data and penalty rows) and :math:`D` the derivative of the
physical quantities (Brune-transformed energies and widths, strengths) with
respect to the reduced-width amplitudes, by central differences of the
transform; :math:`(J^TJ)^{-1}` is taken from the SVD of :math:`J` with its
columns scaled to unit norm, since the columns of one model differ by many
orders of magnitude and the product :math:`J^TJ` loses the poorly constrained
directions to round-off. It is not scaled by :math:`\chi^2/\nu`, it is a
linearization (meaningless for a width the data do not constrain, whose
error comes out larger than any physical value), and a fit stopped at its
evaluation limit gives only an approximate one.

*Over-confident weights.* :math:`\Delta\chi^2` between variants scales as
:math:`1/\sigma^2` of the data. With underestimated errors (THM points with
statistical errors only, digitised errors, :math:`\chi^2/\nu` of 2-3) a
difference that is really a few units becomes tens, and the Akaike weights put
everything on one variant: the model spread then vanishes for the wrong
reason. Rescaling every :math:`\chi^2` by :math:`s = \max(1, \chi^2/\nu)` of
the best variant (``rescale="best"``, ``--rescale best``; the PDG scale factor)
is the minimum correction, and we recommend it whenever the best variant has
:math:`\chi^2/\nu` well above 1. Even then the weights only say which model
reproduces the *shape* of these data best, which for THM includes everything
the HOES model leaves out (resolution, background, digitisation); quote the
flat-weight range next to the weighted one, and treat a variant that takes
all the weight with suspicion rather than as a selection.

*Names and the averaged model.* Quantities are matched across variants by
name: :func:`pyazr.modelavg.parameter_label` names a level energy
``E[<Jπ>#<n>]`` and a width ``G[<Jπ>#<n>;p<pair key>;L<l>;S<s>]`` (``n`` the
engine's level number in the :math:`J^\pi` group, the pair key as the file
writes it), which depends on the level scheme only. ``write_averaged_azr``
writes the averaged energies and widths into a copy of a template ``.azr``
through ``AzrModel``'s setters (a name the template does not have is skipped
with a warning); the radii and THM settings stay the template's. Physical
(observed) widths are averaged, which is meaningful across radii; reduced-width
amplitudes are not, and channels entered as amplitudes (``gammaIsRWA``) carry
them. A width carries the sign of its amplitude; where the sign differs
between variants of appreciable weight the mean is not meaningful (the driver
lists such quantities in its summary). The averaged file is a representative
model that opens in the GUI, not a fit: its :math:`\chi^2` is not any
variant's.

*Size: 19F(p,αγ)16O.* ``examples/f19_pag_thm`` (THM in the adopted window
:math:`E \le 0.45` MeV with JUNA and Spyrou, 25 free widths, 57 points, no
penalty rows for the direct strengths), :math:`a_p` of the three p + 19F
pairs 4.1 / 5.1 / 6.1 fm × ``vertexModel`` pw / dw (``ancai06/kd03:extrapolate``,
quasi-free, no window), each refitted from the example (40 evaluations; the
5.1 fm dw fit stopped on ``xtol`` after 17). :math:`\chi^2` 73.9 / 70.5 / 71.3
(pw) and 87.6 / 91.4 / 84.2 (dw), :math:`\chi^2/\nu` 2.4-3.2. Akaike weights:
pw 0.10 / 0.54 / 0.36, dw ≤ 0.001 (:math:`\Delta\mathrm{AIC}` 14-21); with
``rescale="best"`` (:math:`s = 2.43`) pw 0.20 / 0.41 / 0.35 and dw 0.006-0.024.
In eV (Akaike, unscaled; mean, stat, model; range over the six variants):

====================  ==============  ===========  ===========  ==================
quantity              mean            stat         model        range
====================  ==============  ===========  ===========  ==================
ωγ(213)               0.0065          0.060        0.0037       0.0027-0.020
ωγ(324)               21.4            335          8.1          12.2-32.2
ωγ(11) [1e-29]        2.69            28           0.23         2.48-4.25
====================  ==============  ===========  ===========  ==================

With flat weights ωγ(213) = 0.0138 (model 0.0065) and ωγ(11) = 3.10 (model
0.57) × 10\ :sup:`-29`: the dw variants carry the larger 213 keV strengths
(0.018-0.020 against 0.003-0.016 for pw), and the Akaike weights, which give
them nothing, decide the quoted value. The linearised statistical errors are
larger than the values: without the direct-strength rows the THM scale is
free and the absolute widths rest on the direct data below 0.35 MeV, so in
this example the model spread is the informative number, and the statistical
one says that the data alone do not fix the strengths. (These numbers are
from a rerun with every variant in its own process, after the
``vertex=constant`` memo fix of October 2026: the first run, all six
variants in one process, had the 5.1 and 6.1 fm fits differ by up to
4 × 10\ :sup:`-5` in :math:`\chi^2` and 4 × 10\ :sup:`-6` in the strengths, and
the statistical error of ωγ(213) was 0.059.)

*The driver.* ``scripts/thm_model_average.py <project.azr> --out <dir>``
with the axes ``--radius-pairs K.. --radii R..``, ``--vertex-model pw dw``,
``--vertex constant perlevel onshell``, ``--optical AA/SF ..`` (global
potential names, e.g. ``ancai06/kd03:extrapolate``, or ``coulomb``; applied to
the dw variants), ``--ps delta hulthen:0-50``, ``--lineshape on off`` (the spectator's
:math:`N_C`), ``--distortion none coulomb optical`` (:math:`R(E)` of the pw
variants; ``optical`` takes the potentials of an ``AA/SF`` ``--optical``
value or the project's; the axis collapses for dw variants, whose distortion
is the vertex's own), or a JSON ``--spec``; the grid is their product. ``--strength NAME=Jπ@E`` with ``--strength-in`` /
``--strength-out`` (file pair keys) adds
:math:`\omega\gamma = \frac{2J+1}{(2j_1+1)(2j_2+1)}\Gamma_\mathrm{in}\Gamma_\mathrm{out}/\Gamma`
(open channels, eV) as a derived quantity with its propagated error;
``--weights``, ``--rescale``, ``--prior radius=6.1:0.5``.
``--penalty-hook file.py[:func]`` adds signed residual rows ``func(session,
x)`` (direct strengths against their measured values, priors; with their
Jacobian, or by central differences of ``func``) to the fit, the
:math:`\chi^2`, the weights and the covariance; ``--x-scale`` sets
``least_squares``' ``x_scale`` (``jac`` by default; 1 when such rows have
column norms many orders of magnitude apart from the data's -- on the 19F
model with 17 strength rows ``jac`` rejected every step). Each variant is
fitted in a fresh Python process that hands its result back as JSON, so no
engine state can pass from one variant to the next (a ``vertex=constant``
memo keyed by object address did, before October 2026); ``--in-process``
fits them in the driver's process. A variant AzrModel
or the engine refuses (``vertexModel=dw`` without ``distortion=``, an
optical potential outside its range) is recorded with the reason and
skipped. ``--dry-run`` lists the grid and the refusals AzrModel already
knows, without numpy or the engine. Output: ``variants.csv`` (one row per
variant), ``variants.json``, ``average.json``/``.csv``, ``summary.txt``,
``averaged/<project>_avg.azr`` beside a copy of the data, and ``work/`` with
each variant's project and fitted snapshot. One engine session at a time, by
design: the driver itself holds none. ``tests/pyazr/model_average_test.py`` checks the weights and spreads
against closed forms and the writer; ``tests/pyazr/thm_model_average_test.py``
runs the driver on ``tests/18O_p_a_thm``.
