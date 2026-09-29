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
   coherent sum is exact only in the limit of the exit direction fixed along
   :math:`\hat p_{xA}`; it is kept for comparison with older fits.
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

   *Where the table comes from.* AZURE2 does not compute :math:`R(E)`: it
   needs the spectator and the residual system (masses, charges, the bound
   state of :math:`a = x + s`, the kinematics and the angular acceptance of
   the experiment) and a DWBA/CDCC (e.g. FRESCO) or Nordsieck-integral
   calculation, none of which is part of an R-matrix project. Make the table
   outside and name it here, for instance

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
   The spectator-momentum window over which the HOES cross section is
   averaged (default ``delta``, the quasi-free point; the other forms need the
   four kinematics keys), and its number of Gauss-Legendre nodes (1-64,
   default 16), below.
``theta=``, ``distortion=``
   Reserved for later versions (angular acceptance, distortion):
   refused with ``not implemented yet``, so that no file can rely on them
   silently.

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
c12c12_tumino2018`` at its parameters, the four THM segments one experiment):
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

*What is averaged.* At one :math:`p_s` the THM double differential cross
section is

.. math::

   \frac{d^2\sigma}{d\Omega_{sF}\,dE} \propto |\phi_a(p_s)|^2 \sum_l
   \Bigl|\dots M_l\bigl(p_{xA}(p_s)\bigr)\Bigr|^2

(2017 eq. 34): :math:`\phi_a` multiplies the whole cross section and the
vertex is evaluated at that event's :math:`p_{xA}`. Events with different
:math:`p_s` are different final states (the spectator's momentum is measured),
so they add incoherently: the *cross section* is averaged, not the amplitude.
The HOES "data" are then obtained (2017 eq. 35; Tribble 2014 eq. 4.1) by
dividing the yield of an energy bin by the product of the kinematic factor
and :math:`|\phi_a|^2`, computed by a Monte Carlo simulation of the setup over
the same bin and :math:`p_s` window (Pizzone et al., PRC 83 (2011) 045801;
Spitaleri et al., PRC 95 (2017) 035801, eqs. 3-4 with :math:`|p_n| \le 30`
MeV/c), i.e. :math:`N(E)/\int_W |\phi_a|^2\,d^3p_s`. The
yield itself is :math:`\int_W |\phi_a(p_s)|^2 \sigma(E; p_s)\,d^3p_s`, so the
datum at :math:`E` is the event-weighted mean

.. math::

   \bar\sigma(E) = \frac{\int_{p_{\min}}^{p_{\max}} w(p)\,\sigma(E; p)\,dp}
                        {\int_{p_{\min}}^{p_{\max}} w(p)\,dp},
   \qquad w(p) = |\phi_a(p)|^2 p^2

for an isotropic window :math:`p_{\min} \le |p_s| \le p_{\max}` (the kinematic
factor taken as constant across one bin). This is what AZURE2 computes: with
``ps=`` on the experiment line, every THM point of the experiment (and every
sub-point of its resolution folding) evaluates the entrance vertex -- the
surface term and the Coulomb term :math:`C_l`, which depends on :math:`p_{xA}`
too -- at the nodes :math:`p_k` of a Gauss-Legendre rule on
:math:`[p_{\min}, p_{\max}]` with :math:`T_k = p_k^2/2\mu_{sx}` added to
:math:`E + B` (exactly what ``spectatorEnergy`` does with one value), and sums
the HOES cross section of each node with the normalized weights
:math:`w_k \propto \omega_k w(p_k)`, :math:`\sum_k w_k = 1`, :math:`\omega_k`
the Gauss-Legendre weights. The level matrix is common to the nodes; only the
vertex changes. If the data were instead divided by :math:`|\phi_a|^2`
event by event, the weight is the phase space :math:`p^2` alone -- use a
table.

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
   The event weight itself, :math:`w(p)` per unit :math:`p` -- i.e.
   :math:`|\phi|^2 p^2`, or :math:`p^2` alone, or a measured :math:`|p_s|`
   distribution of the accepted events -- two columns, :math:`p` (MeV/c,
   strictly increasing, :math:`\ge 0`) and :math:`w \ge 0`, ``#`` comments, at
   least two rows, some positive weight; linear between rows. The window is
   the table's range. A relative path is taken from the directory of the
   ``.azr``, exactly as for ``weight[k]=`` files (the same resolution in
   ``Config::ReadThmBlock``: a leading ``/`` or ``\`` or a drive letter is
   absolute), so a copy of the project elsewhere must make it absolute or
   copy the table; the name cannot contain blanks or ``#``.
``psNodes=N``
   Gauss-Legendre nodes (1-64, default 16) on the whole window. With
   :math:`p_{\min} = p_{\max}` the window is one node of weight 1. A table
   with kinks is integrated to :math:`O(h^2)` only; raise ``psNodes`` or
   smooth it.

:math:`\mu_{sx}` comes from the experiment's kinematics (spectator and the
nucleus :math:`x` = Trojan horse minus spectator; built-in table masses), so a
window needs ``beam``, ``target``, ``spectator`` and ``Ebeam`` (``Ebeam`` is
not used by the window). Refused, with ``ERROR: <thm> experiment[...]``: a
window without kinematics, a malformed ``ps`` or ``psNodes``, ``psNodes``
without a window, an unreadable or invalid table, and a window together with
a non-zero ``spectatorEnergy`` (global or for the experiment's entrance pair).
The window is printed at startup (nodes, :math:`\mu_{sx}`, the range of
:math:`T_s` and :math:`\langle T_s \rangle`) and written to
``thm_experiments.out`` (one ``ps_node`` row per node: :math:`p_k`,
:math:`w_k`, :math:`T_k`). Extrapolation segments (``<segmentsTest>``) are in
no experiment and keep the quasi-free vertex. Cost: the vertex and the
exit-channel sums once per node (the level matrix once); the 7Li and 12C+12C
examples below run in about 1.8 and 1.1 times the time without a window.

*pyazr.* ``AzrModel.set_thm_experiment(..., ps="hulthen:0-40", psNodes=16)``;
``session.thm_vertex(name, energies)`` returns the nodes, weights,
:math:`T_k`, :math:`\rho` at every node and energy, and per entrance channel
and level the window average :math:`\langle |M_l|^2 \rangle` and the
quasi-free :math:`|M_l(p_s = 0)|^2` with the boundary the vertex uses for
that level (for a ``delta`` experiment one node, :math:`p_s = 0`,
:math:`T` = its ``spectatorEnergy``).

*Checks* (``tests/thm_spectator_window``, ``tests/pyazr/thm_spectator_window_test.py``,
on ``tests/18O_p_a_thm`` with its deuteron Trojan horse, :math:`\mu_{sx} =
469.46` MeV): ``ps=delta`` byte-identical to no key; a window shrinking to a
point equals ``spectatorEnergy`` :math:`= p^2/2\mu_{sx}` (one node: exactly; a
:math:`2\times 10^{-4}` MeV/c window: :math:`6\times 10^{-10}`); a flat table on
[20, 40] MeV/c equals the Simpson average of 41 single-``spectatorEnergy``
runs to :math:`8\times 10^{-8}` (the Simpson error); two nodes equal
:math:`w_1 m(T_1) + w_2 m(T_2)` from two sessions to :math:`10^{-15}`; the
Hulthén window [0, 40] with 16 nodes agrees with 32 to :math:`10^{-10}` (8
nodes: :math:`2\times 10^{-8}`); ``thm_vertex`` against an independent
evaluation of nodes, weights, :math:`\rho` and :math:`\langle|M_0|^2\rangle`.

*Size* (models at fixed parameters, not refitted; scratch studies, not in
the repository).

- A vertex node. With the spectator at a fixed :math:`p_s = 24` MeV/c,
  ``tests/18O_p_a_thm`` (no folding) puts :math:`M_0 = 0` at :math:`E \approx
  0.70` MeV: the model there drops to :math:`2.4\times 10^{-4}` of its
  maximum. The Hulthén window [0, 40] fills it to 0.42 of the maximum
  (:math:`5\times 10^{4}` times the point value). Across the whole data
  range the window changes the (folded) model by up to 73 % against the
  quasi-free point.
- :sup:`7`\ Li(p,α) (``examples/li7_tumino2006``, deuteron Trojan horse as
  its :math:`B = 2.2246` MeV, Hulthén [0, 40]): the THM model rises by
  1.14-1.18 below 1 MeV and by 1.15-1.29 up to 3 MeV (0.92 at 4.5 MeV), so
  the shape changes by ~3 % below 1 MeV and ~40 % over 0.08-6.9 MeV; the THM
  :math:`\chi^2` (norm profiled) goes 847 → 933. Replacing the window by its
  mean, ``spectatorEnergy`` :math:`= \langle T_s\rangle = 0.88` MeV, gives the
  window's model to 1 % (:math:`\chi^2` 942): here the window acts as a shift
  of :math:`E + B`. With a :sup:`3`\ He Trojan horse (:math:`B = 5.4935` MeV,
  Eckart :math:`a = 0.42`, :math:`b = 1.2` fm\ :sup:`-1`,
  :math:`\mu_{sx} = 625.4` MeV, [0, 40]): 0.91-1.15, :math:`\chi^2` 1154 →
  1227.
- 12C+12C (``examples/c12c12_tumino2018``, 12C(14N,α/p)d at 30 MeV, the four
  THM segments one experiment with a free norm; :math:`B = 10.27` MeV,
  :math:`\mu_{sx} = 1606` MeV; the d-12C distribution in 14N as an Eckart
  function with :math:`a = \sqrt{2\mu_{sx}B}/\hbar c = 0.918` fm\ :sup:`-1`,
  :math:`b = 2.5` fm\ :sup:`-1`, nearly flat over [0, 40] MeV/c, so the
  weight is close to :math:`p^2`; :math:`T_s \le 0.49` MeV): :math:`\rho
  \approx 13.4`-14.7 moves by up to 0.25, a fair part of the oscillation of
  :math:`j_l` for :math:`l \le 8`. The model changes by 0.32-1.73 level by
  level (α\ :sub:`0` 0.34-1.12, p\ :sub:`1` 0.43-1.73), and the THM
  :math:`\chi^2` at the published parameters goes 61 → 1578; with
  ``spectatorEnergy`` = :math:`\langle T_s\rangle = 0.295` MeV instead it is
  1707, and the window differs from that mean by up to 1.73: near the nodes
  the average matters, not only the shift. ``thm_vertex``: the :math:`l = 4`
  vertex has a node at :math:`E = 2.02` MeV (:math:`|M_4|^2 = 8\times
  10^{-6}` against a maximum of 0.68), which the window fills to
  :math:`4.8\times 10^{-2}`. Published 12C+12C THM fits assume the
  quasi-free vertex; with the window the parameters would move.

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
