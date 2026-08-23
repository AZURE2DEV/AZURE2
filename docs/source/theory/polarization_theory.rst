Polarization Observables in R-Matrix Theory
===========================================

This chapter develops the theory behind the vector analyzing power as AZURE2
computes it. It assumes the R-matrix formalism of Lane and Thomas but nothing
about spin algebra beyond angular-momentum coupling.

Why an analyzing power is worth measuring
-----------------------------------------

A differential cross section is a sum of squared amplitudes. Squaring destroys
phase, and phase is where much of the physics lives: whether two levels
interfere constructively or destructively, whether a resonance is a
:math:`3/2^-` or a :math:`5/2^+`, whether a broad background amplitude has the
sign that a direct-capture calculation predicts. A cross section constrains
these things only weakly, because a fit can usually trade a phase against a
width and land on the same curve.

A polarization observable is different. It is built from an *interference*
between amplitudes that differ in the spin projection of the projectile, and it
appears in the numerator of a ratio whose denominator is the cross section. The
normalization cancels, the overall phase cancels, and what survives is precisely
the relative phase. Where a cross section has a smooth maximum, an analyzing
power can swing from :math:`+1` to :math:`-1` over a few tens of keV.

This is why the elastic scattering of polarized protons from :sup:`12`\ C is used
as a beam polarimeter, and why the analyzing power of :sup:`9`\ Be(p,d) or
:sup:`9`\ Be(p,\ :math:`\alpha`) adds information that no amount of cross-section
data supplies. For astrophysics the payoff is indirect but real: a spin and
parity assignment that a polarization measurement settles is one fewer degree
of freedom in the extrapolation of an S-factor to stellar energies.

The spin structure of the scattering amplitude
----------------------------------------------

For spinless particles a single complex function :math:`f(\theta)` carries all
the information, and :math:`d\sigma/d\Omega = |f|^2`. With spin, the amplitude
becomes a matrix in spin space. Write the entrance channel spin as :math:`s`
with projection :math:`\nu`, and the exit channel spin as :math:`s'` with
projection :math:`\nu'`. The object we need is

.. math::

   M_{s'\nu' s\nu}(\theta),

the amplitude for the system to enter with spin projection :math:`\nu` and leave
with :math:`\nu'`. The channel spin is the vector sum of the two intrinsic
spins, :math:`\mathbf{s} = \mathbf{j}_1 + \mathbf{j}_2`, so for a spin-1/2
projectile on a spin-0 target :math:`s = 1/2` and :math:`\nu = \pm 1/2` is just
the projectile's own spin projection.

Every observable follows from :math:`M` by a density-matrix construction. If the
beam is described by a spin density matrix :math:`\rho_{\text{in}}`, the outgoing
state is

.. math::

   \rho_{\text{out}} = M \rho_{\text{in}} M^{\dagger},

and any observable is a trace of :math:`\rho_{\text{out}}` against the operator
that the detector measures. The unpolarized cross section is the trace with the
identity,

.. math::

   \frac{d\sigma}{d\Omega} = \frac{1}{2s+1}\sum_{s'\nu'}\sum_{\nu}
                             \left| M_{s'\nu' s\nu} \right|^2,

averaging over entrance projections and summing over exit ones. The vector
analyzing power is the trace with a Pauli matrix. Once :math:`M` exists there is
no further angular-momentum algebra to do; this is the reason for organizing the
calculation around it.

Seyler's expression for the transition matrix
---------------------------------------------

The step that makes this practical for an R-matrix code is that the amplitude
can be written directly in terms of the collision matrix :math:`U` — the very
object an R-matrix code already builds. The expression is Lane and Thomas's
[LaneThomas1958]_; we follow the form in which Seyler [Seyler1969]_ quotes it
as his Eq. (4), which is the one AZURE2's channel-spin bookkeeping matches
term for term:

.. math::
   :label: seyler

   M_{s'\nu' s\nu}(\theta) = \frac{\sqrt{\pi}}{k}\Bigg[
     -\,C(\theta)\,\delta_{ss'}\delta_{\nu\nu'}
     + i \sum_{J l l'} \sqrt{2l+1}\;
       (s\,l\,\nu\,0 | J\,\nu)\;
       (s'\,l'\,\nu'\,\nu-\nu' | J\,\nu)\\
       \times\; e^{i(\omega_l + \omega_{l'})}
       \left( \delta_{ss'}\delta_{ll'} - U^J_{s'l' sl} \right)
       Y_{l'}^{\nu-\nu'}(\theta,0)
   \Bigg].

Each piece earns its place:

*The Coulomb term* :math:`C(\theta)` is the Rutherford amplitude. It is diagonal
in channel spin and its projection, since a pure Coulomb field does not touch
spin, and it exists only when the entrance and exit pairs are the same.

*The first Clebsch–Gordan coefficient* :math:`(s\,l\,\nu\,0|J\nu)` has the
entrance orbital projection fixed at zero. This is the choice of quantization
axis: the beam travels along :math:`\hat z`, so :math:`\mathbf{l}\cdot\hat z = 0`
for the incoming wave, and the total projection is therefore :math:`\nu`, the
channel-spin projection alone.

*The second Clebsch–Gordan coefficient* :math:`(s'\,l'\,\nu'\,\nu-\nu'|J\nu)`
enforces conservation of the total projection: whatever the exit channel spin
takes as :math:`\nu'`, the exit orbital motion must carry
:math:`\mu = \nu - \nu'`.

*The spherical harmonic* :math:`Y_{l'}^{\nu-\nu'}` is the angular function of
that outgoing orbital motion, and it is the structural reason this observable
needed new machinery rather than a new formula on top of the old.

An unpolarized cross section sums over spin projections, and only the
:math:`\mu = \nu - \nu' = 0` terms survive that sum. Its angular dependence is
therefore carried entirely by :math:`Y_{l'}^{0} \propto P_{l'}(\cos\theta)`,
which is why an R-matrix code written for cross sections — AZURE2 included —
computes Legendre polynomials and stores nothing else.

A non-zero analyzing power comes precisely from the spin-flip amplitudes,
:math:`\nu \neq \nu'`, hence :math:`\mu \neq 0`, hence
:math:`Y_{l'}^{\mu \neq 0}`, which needs the *associated* Legendre functions
:math:`P_{l'}^{\mu}`. Those terms are not merely absent from a
Legendre-only code; they are not representable in it. This is why the
implementation adds an angular basis rather than post-processing the existing
one — see :doc:`polarization_implementation`.

*The bracket* :math:`\delta_{ss'}\delta_{ll'} - U^J_{s'l'sl}`, dressed with the
Coulomb phases :math:`e^{i(\omega_l+\omega_{l'})}`, is the transition matrix.
Section :doc:`polarization_implementation` shows that AZURE2 already forms
exactly this quantity, phases included, so nothing has to be reconstructed.

The vector analyzing power
--------------------------

The restriction that matters is on the **projectile**: the vector analyzing
power is a spin-1/2 beam observable, so :math:`j_1 = 1/2` throughout. The
*target* spin is unrestricted, and the target is taken unpolarized, so its
projection is traced over.

That trace is the one subtlety. :math:`A_y` is defined with respect to the
polarization of the projectile alone, but :eq:`seyler` delivers amplitudes in
the **channel-spin** basis, in which projectile and target spin are already
coupled. The entrance index must therefore be decomposed before a Pauli matrix
can act on the projectile:

.. math::
   :label: decompose

   M_{\text{out};\,m_1 m_2} =
     \sum_s \left( j_1\, j_2\, m_1\, m_2 \,\middle|\, s,\, m_1{+}m_2 \right)
     M_{\text{out};\, s,\, m_1+m_2}.

The coupling order and phases follow Lane and Thomas — :math:`\mathbf{s} =
\mathbf{I}_1 + \mathbf{I}_2` with Condon–Shortley coefficients — since that is
the convention the rest of AZURE2 is built on. Nothing in the unpolarized cross
section fixes this order, because it adds channel spins incoherently and is
blind to it; the analyzing power is not.

For a spin-0 target the sum collapses to a single term, the channel spin *is*
the projectile's spin, and :math:`M` is a plain :math:`2\times 2` matrix in the
projectile projection — the case of
:math:`{}^{12}\mathrm{C}(\vec p, p)`. For a spin-1/2 target such as
:sup:`15`\ N the channel spins are 0 and 1 and never 1/2, so a code that looks
for a channel spin of 1/2 finds nothing and reports zero for a perfectly well
defined, non-zero observable.

Take the quantization axis for the polarization along

.. math::

   \hat{\mathbf{n}} = \frac{\mathbf{k}_{\text{in}} \times \mathbf{k}_{\text{out}}}
                           {|\mathbf{k}_{\text{in}} \times \mathbf{k}_{\text{out}}|},

the normal to the scattering plane. This is the Madison convention, and it is
the only sensible choice: parity conservation forbids a vector polarization
along any other direction, so :math:`A_x = A_z = 0` identically and
:math:`A_y` is the whole of the vector analyzing power.

The observable is a ratio of traces: the outgoing density matrix against
:math:`\sigma_y` acting in the projectile's spin space, over the same trace
against the identity,

.. math::

   A_y = \frac{\mathrm{Tr}\left( M \sigma_y M^{\dagger} \right)}
              {\mathrm{Tr}\left( M M^{\dagger} \right)},

both traces running over the exit configuration and over the target projection
:math:`m_2`, which is where the unpolarized target is averaged away.

Writing this out with :math:`(\sigma_y)_{+-} = -i` and
:math:`(\sigma_y)_{-+} = +i` gives the form the code evaluates:

.. math::
   :label: ay

   A_y = \frac{2 \displaystyle\sum_{s'\nu'}
                \mathrm{Im}\!\left[\,
                  M_{s'\nu',\,s\,+\frac{1}{2}}\;
                  M^{*}_{s'\nu',\,s\,-\frac{1}{2}} \right]}
             {\displaystyle\sum_{s'\nu'}
                \left( \left|M_{s'\nu',\,s\,+\frac{1}{2}}\right|^2
                     + \left|M_{s'\nu',\,s\,-\frac{1}{2}}\right|^2 \right)}.

The structure is worth pausing on. The numerator is an interference between the
two entrance projections; the denominator is the cross section. If the two
projections scatter identically the numerator vanishes. If they scatter
differently but *in phase*, the imaginary part vanishes and :math:`A_y` is still
zero. **A non-zero analyzing power requires both spin dependence and a relative
phase**, which is exactly why it is a sharp probe of interference and why it
peaks near resonances, where phases move fast.

Four places where the analyzing power must vanish
--------------------------------------------------

These are the checks that tell you an implementation is right, and each one
fails differently when it is wrong.

**At** :math:`\theta = 0` **and** :math:`\theta = 180^\circ`. The normal
:math:`\hat{\mathbf{n}}` is undefined when :math:`\mathbf{k}_{\text{out}}` is
parallel or antiparallel to :math:`\mathbf{k}_{\text{in}}`, so there is no
direction for the polarization to point along. Formally, the spherical harmonic
:math:`Y_{l'}^{\mu}` vanishes at :math:`\theta = 0` for every
:math:`\mu \neq 0`, so no spin-flip amplitude survives. Every curve in Baumann's
fig. 1 goes to zero at both ends.

**In the pure Coulomb limit.** Rutherford scattering is spin-independent: the
Coulomb term in :eq:`seyler` is diagonal in :math:`\nu`, so both projections
acquire the same amplitude and the same phase. At energies well below any
resonance :math:`A_y \to 0`. This is a useful trap — checking an implementation
only at low energy will show a perfectly correct zero and prove nothing.

**Without spin-orbit splitting.** If the nuclear interaction did not distinguish
:math:`j = l + 1/2` from :math:`j = l - 1/2`, the sum over :math:`J` in
:eq:`seyler` would collapse and the spin-flip amplitude would cancel. The
analyzing power in :sup:`12`\ C(:math:`\vec p`,p) is large precisely because the
:math:`p_{3/2}` and :math:`d_{5/2}` resonances are spin-orbit partners of states
that lie elsewhere.

**When only one partial wave contributes.** A single resonance with no
background produces a common phase, which cancels in the ratio. The structure in
:math:`A_y` between 1.6 and 1.8 MeV comes from the :math:`3/2^-` and
:math:`5/2^+` resonances interfering *with each other* and with the non-resonant
:math:`s_{1/2}` and :math:`p_{1/2}` phases.

The bound :math:`|A_y| \le 1` follows from the Cauchy–Schwarz inequality applied
to :eq:`ay` and is the cheapest sanity check available; a value outside
:math:`[-1,1]` means the amplitude matrix is not a consistent set of amplitudes.

The sign, and why it cannot be checked internally
--------------------------------------------------

Every check listed above is invariant under :math:`A_y \to -A_y`. The bound
holds, the zeros stay zeros, the Coulomb limit is still zero. The overall sign
depends on a chain of conventions — the direction assigned to
:math:`\hat{\mathbf{n}}`, the Condon–Shortley phase in the spherical harmonics,
the ordering of the coupling in the Clebsch–Gordan coefficients — and getting
any one of them backwards flips it without disturbing anything else.

The sign therefore has to be fixed against measurement. The comparison used here
is described in the next chapter: Baumann *et al.* mark four points where
:math:`|A_y|` reaches unity, three positive and one negative, and reproducing
that pattern — not merely the positions — is what pins the convention down.

Analyzing power through a target of finite thickness
-----------------------------------------------------

One subtlety has no analogue in cross-section work and is easy to get wrong.

A cross section measured on a target of finite thickness is an integral: the
beam loses energy as it goes, and the yield is

.. math::

   Y = \int_{E_{\text{back}}}^{E_{\text{surface}}}
        \frac{\sigma(E)}{\varepsilon(E)}\, dE

with :math:`\varepsilon` the stopping power. An analyzing power is not an
integral of this kind, because it is a *ratio*. Averaging :math:`A_y(E)` along
the target with equal weight has no physical meaning. What the experiment
returns is the ratio of the polarized and unpolarized yields, so the correct
average is weighted by the cross section:

.. math::
   :label: aytarget

   \langle A_y \rangle =
      \frac{\displaystyle\int A_y(E)\,\sigma(E)\,dE}
           {\displaystyle\int \sigma(E)\,dE}.

The depths where the reaction is likely contribute more, as they must.

The consequence is quantitative and severe. In the :sup:`12`\ C + p system below
2 MeV the cross section is dominated by Rutherford scattering, which diverges as
:math:`E^{-2}` at low energy, while :math:`A_y` there is essentially zero. A
thick target therefore dilutes the analyzing power towards zero: for the gas
target in the ``tests/13N`` evaluation, whose integration grid runs from 0.42 to
1.54 MeV, the weighted average of a resonant :math:`A_y \approx 0.8` comes out
around :math:`10^{-6}`. This is physics, not a numerical artefact, and it is why
analyzing-power measurements use thin targets — Baumann's were 85 nm of
:sup:`12`\ C, about 3.1 keV of energy loss at 1.7 MeV.

Scope
-----

What is implemented covers the vector analyzing power :math:`A_y` for a
spin-1/2 projectile on a target of **any** spin, in **both** particle and
capture exit channels — the two by different routes, since the photon exit has
no amplitude matrix of the form of :eq:`seyler`.

*Capture channels* use Seyler and Weller [SeylerWeller1979]_, who give the
Legendre coefficients of the angular distribution directly in the channel-spin
representation — the representation AZURE2 already works in, so their
:math:`R` is the T-matrix element the code already forms. Writing

.. math::
   :label: swcapture

   \sigma(\theta,\phi) = N \sum_k \Bigl[ a_k P_k(\cos\theta)
                       + b_k P_k^1(\cos\theta)\, p_y \Bigr]

their Eqs. (20) and (21) give :math:`a_k` and :math:`b_k` as sums over pairs of
reaction pathways :math:`t = \{p L b l s\}` weighted by
:math:`\mathrm{Re}\,R R'^*` and :math:`\mathrm{Re}\,(i R R'^*)`, so that

.. math::

   A_y(\theta) = \frac{\sum_k b_k P_k^1(\cos\theta)}
                      {\sum_k a_k P_k(\cos\theta)} .

The structural point is that :math:`a_k` requires :math:`s = s'` while
:math:`b_k` does not. The channel-spin off-diagonal terms are exactly the
information a polarization measurement adds and a cross section cannot carry,
and they are why the pathway pairs have to be enumerated across channel-spin
groups rather than within one.

Their Eq. (12) carries a :math:`(-1)^M` phase that several earlier treatments
drop — so their :math:`P_L^M` is the associated Legendre function *without* the
Condon-Shortley phase, and for :math:`M = 1` that is a sign. Two further
conventions had to be settled numerically against their worked example; see
:doc:`polarization_implementation`.

One extension remains unimplemented and should not be assumed to work:
*tensor observables* (:math:`T_{20}`, :math:`T_{22}`, and the rest) need a
spin-1 projectile and the corresponding rank-2 operators. For particle channels
the amplitude matrix built here is general enough in principle — it carries all
channel spins and projections — but the observable side would have to be
written; for capture, Seyler and Weller's Eqs. (22)–(25) give the
:math:`c_k`, :math:`d_k` and :math:`e_k` coefficients that would be needed, in
the same notation as the :math:`b_k` already coded.

.. [LaneThomas1958] A. M. Lane and R. G. Thomas, *Reviews of Modern Physics*
   **30** (1958) 257. The formalism AZURE2 is built on; the source of
   :eq:`seyler` and of the channel-spin coupling convention used throughout.

.. [Seyler1969] R. G. Seyler, *Nuclear Physics* **A124** (1969) 253. Quotes
   the Lane–Thomas amplitude as its Eq. (4) before specialising to spin-1/2
   on spin-1; only that general equation is used here.

.. [SeylerWeller1979] R. G. Seyler and H. R. Weller,
   *Physical Review C* **20** (1979) 453.

.. [Baumann1992] R. Baumann *et al.*, *Nuclear Physics* **A542** (1992) 53.
