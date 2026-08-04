Computing the Analyzing Power in AZURE2
=======================================

This chapter is the companion to :doc:`polarization_theory`. It describes what
was actually built, how each step was tested, what went wrong, and how far the
result can be trusted. Everything here is checkable against the source.

What the code already had, and the one thing it did not
--------------------------------------------------------

Seyler's Eq. (4) needs five ingredients. AZURE2 had four of them:

* **the collision matrix** :math:`U`, from the A- or R-matrix inversion;
* **Coulomb and hard-sphere phases**, ``EPoint::GetExpCoulombPhase`` and the
  hard-sphere phases stored per channel;
* **the Coulomb amplitude** :math:`C(\theta)`, ``EPoint::GetCoulombAmplitude``;
* **Clebsch–Gordan and Racah coefficients**, ``AngCoeff::ClebGord`` and
  ``AngCoeff::Racah``;
* and the whole code is written in the **channel-spin representation** already,
  so :math:`s` and :math:`s'` are available per channel rather than needing to be
  recoupled out of a :math:`jj` scheme.

The missing ingredient was the angular functions. A search of ``src/`` and
``include/`` turned up exactly one angular function, ``EPoint::GetLegendreP(L)``
— no associated Legendre functions and no spherical harmonics anywhere. As
:doc:`polarization_theory` explains, this is not an oversight but a direct
consequence of what the code had been asked to compute: unpolarized cross
sections need only :math:`\mu = 0`.

GSL was already a dependency and provides ``gsl_sf_legendre_sphPlm``, the
normalized associated Legendre function, so the addition was small
(``src/AngCoeff.cpp:37``):

.. code-block:: cpp

   complex AngCoeff::SphericalHarmonic(int l, int m, double theta, double phi) {
     if (l < 0 || std::abs(m) > l) return complex(0.0, 0.0);
     const int am = std::abs(m);
     const double norm = gsl_sf_legendre_sphPlm(l, am, std::cos(theta));
     complex y = norm * complex(std::cos(am * phi), std::sin(am * phi));
     if (m < 0) { y = std::conj(y); if (am % 2) y = -y; }
     return y;
   }

The negative-\ :math:`m` branch applies
:math:`Y_l^{-m} = (-1)^m \left(Y_l^{m}\right)^{*}`. GSL's ``sphPlm`` already
carries the Condon–Shortley phase, so it is not applied again here — applying it
twice is a mistake that leaves :math:`|A_y|` untouched and flips its sign, which
is precisely the class of error that internal checks cannot catch.

The step that made the rest easy
---------------------------------

The most useful thing found during this work was not written; it was already
there. Look at how ``AMatrixFunc`` forms its T-matrix element
(``src/AMatrixFunc.cpp:408``):

.. code-block:: cpp

   complex uphase = coulombPhaseEn * hardSpherePhaseEn
                  * coulombPhaseEx * hardSpherePhaseEx;
   complex umatrix = 2.0 * complex(0.0,1.0) * sqrtPenEn * sqrtPenEx
                   * this->GetUBilinear(jNum, chNum, chpNum);
   complex tphase  = coulombPhaseEn * coulombPhaseEn;
   if (chNum == chpNum) tmatrix = tphase - uphase * (1.0 + umatrix);
   else                 tmatrix = -uphase * umatrix;

Since :math:`U = \texttt{uphase}\,(1 + \texttt{umatrix})` and ``tphase`` is
:math:`e^{2i\omega_l}`, this is

.. math::

   \texttt{tmatrix} = e^{i(\omega_l + \omega_{l'})}
                      \left( \delta_{ss'}\delta_{ll'} - U^J_{s'l'sl} \right),

which is the bracket of Eq. :eq:`seyler` **with its exponential already
applied**. No phase has to be reconstructed, no convention has to be guessed,
and — most importantly — the polarization code is guaranteed to see the same
:math:`U` that the cross-section code sees, including every boundary-condition
and Brune-parametrization subtlety. Had this not been the case the sensible
approach would have been to rebuild :math:`U` from scratch, which would have
doubled the surface area for disagreement.

The amplitude matrix
--------------------

``Polarization::AmplitudeMatrix`` (``include/PolarizationFunc.h``) holds one
complex number per :math:`(s,\nu,s',\nu')` and is filled by accumulating one
reaction pathway at a time. The constructor enumerates the channel spins
available in the entrance and exit pairs, from :math:`|j_1 - j_2|` to
:math:`j_1 + j_2`, and allocates a slot for every projection pair.

``AddPathway`` is Eq. :eq:`seyler` transcribed:

.. code-block:: cpp

   for (double v = -s; v <= s + 1.e-6; v += 1.0) {
     const double cg1 = AngCoeff::ClebGord(s, (double)l, jValue, v, 0.0, v);
     if (std::fabs(cg1) < 1.e-12) continue;
     for (double vp = -sp; vp <= sp + 1.e-6; vp += 1.0) {
       const double mu = v - vp;
       if (std::fabs(mu) > lp + 1.e-6) continue;
       const double cg2 = AngCoeff::ClebGord(sp, (double)lp, jValue, vp, mu, v);
       if (std::fabs(cg2) < 1.e-12) continue;
       const complex y = AngCoeff::SphericalHarmonic(lp, (int)std::lround(mu), theta_);
       At(s, v, sp, vp) += complex(0.0,1.0) * std::sqrt(2.0*l+1.0)
                           * cg1 * cg2 * tMatrixElement * y;
     }
   }

The loops run over half-integers by stepping in units of one from
:math:`-s`, with a tolerance on the upper bound, and comparisons of spins use a
tolerance rather than equality — a spin of :math:`1/2` is exactly representable
but :math:`\nu - \nu'` accumulated in floating point is not reliably so.

``AddCoulomb`` adds :math:`-C(\theta)` on the diagonal in
:math:`(s,\nu)`, and only when the entrance and exit pairs coincide.
``AnalyzingPowerAy`` evaluates Eq. :eq:`ay` directly.

Validating it, in the order the checks were applied
----------------------------------------------------

The strategy was to test the amplitude matrix against something already trusted
*before* asking it for a new observable, so that a disagreement could only come
from the new code.

**Gate 1 — the unpolarized cross section, checked as angle independence.**
:math:`M` must reproduce what ``GenMatrixFunc::CalculateCrossSection`` produces
by the Blatt–Biedenharn route. The two differ by kinematic factors that are
tedious to match exactly, so the test used was sharper and needed no bookkeeping:
at fixed energy, the *ratio* of the two must be constant in angle. Any error in
the coupling order, or in the choice of :math:`l` against :math:`l'` inside the
spherical harmonic, gives an angle-dependent ratio. The ratio came out constant
to between :math:`10^{-9}` and :math:`10^{-11}`, which is the level of the
accumulated floating-point noise.

This gate is the one that matters. It tests the Clebsch–Gordan ordering, the
spherical-harmonic indices, the phase convention and the pathway enumeration all
at once, against a code path that has been in production for years.

**Gate 2 — the bounds and the zeros.** :math:`|A_y| \le 1` everywhere;
:math:`A_y \to 0` at :math:`\theta \to 0`; :math:`A_y \to 0` far below the
resonances where scattering is pure Coulomb.

A false negative was encountered here and is worth recording. On first inspection
:math:`A_y` appeared to be identically zero, and the natural conclusion was that
the spin-flip amplitudes were not being built. In fact only low-energy points had
been examined, where :math:`A_y = 0` is the *correct* answer. Scanning the full
energy range showed the expected structure. A test whose passing value and whose
failure value are both zero proves nothing.

**Gate 3 — the sign, against data.** As :doc:`polarization_theory` explains, no
internal check can fix the overall sign. Baumann *et al.* [Baumann1992]_ mark
four points in their fig. 4 where :math:`|A_y|` reaches unity. Scanning the
``tests/13N`` model over the same window reproduces all four, including the
pattern of three positive extrema and one negative:

.. list-table::
   :header-rows: 1
   :widths: 30 30 20

   * - Baumann fig. 4
     - computed
     - :math:`\Delta E`
   * - 1.665 MeV, 97°, +
     - 1.670, 100°, +0.998
     - 5 keV
   * - 1.725 MeV, 120°, +
     - 1.740, 125°, +0.998
     - 15 keV
   * - 1.750 MeV, 80°, −
     - 1.750, 80°, −0.992
     - 0
   * - 1.770 MeV, 145°, +
     - 1.775, 145°, +0.998
     - 5 keV

The agreement in position is within the precision of reading crosses off a
contour plot. The agreement in sign pattern is what fixes the convention.

That this works at all is a non-trivial statement about the model rather than
just the code: the ``tests/13N`` evaluation carries a :math:`3/2^-` level at
:math:`E_x = 3.503` MeV and a :math:`5/2^+` at 3.545 MeV, fitted to cross-section
and capture data with no knowledge of any polarization measurement, and
Baumann's table 1 gives 3.499 and 3.546 MeV from their phase-shift analysis of
:math:`A_y`. Two independent routes to the same two states.

Making it an observable
-----------------------

Computing :math:`A_y` at a point is one thing; letting a user *ask* for it is
another. The observable code 7 was added to the segment machinery, and two
defects surfaced that are worth describing because both were silent.

**A valid extrapolation line was rejected outright.** ``ExtrapLine`` reads its
ten fields and then calls ``getline`` to look for optional advanced-segment
data. For a plain segment those ten fields are the whole line, so the stream is
already at end-of-file and ``getline`` sets ``failbit``. ``EData::MakePoints``
tests ``rdstate()`` and returns :math:`-1` on it — which aborts the parse of
``<segmentsTest>`` entirely. One plain segment therefore discarded every test
segment in the file, with no message. It had gone unnoticed because the GUI
always writes the trailing field; only a hand-written or generated ``.azr``
trips it. The fix is to read only when something is left.

**Target-effect sub-points never learned what they were.** A segment carrying
target integration is not evaluated at its own points. ``EPoint::Calculate``
builds sub-points spanning the target thickness, evaluates those, and combines
them into a yield. Sub-points are copied from the parent *before* the observable
was stamped on it, so they never became analyzing-power points and the
:math:`A_y` branch never ran for them. The symptom was that of three declared
:math:`A_y` segments, exactly one — the one that happened to have no target
effect — returned an analyzing power, and the other two returned cross sections.

Diagnosing this took longer than it should have because two plausible
explanations had to be eliminated first: the stored points were verified to
carry the flag correctly after ``MakePoints``, and pointer identity showed that
the points being *calculated* were different objects entirely, with a null parent
``EData``. Only then did the sub-point mechanism become the obvious candidate.
``AddSubPoint`` now propagates the observable (``src/EPoint.cpp:1577``).

Averaging a ratio over a target
--------------------------------

Fixing the flag exposed the physics point of Eq. :eq:`aytarget`: the existing
integrator computes a yield, and running it on :math:`A_y` as though
:math:`A_y` were a cross section is meaningless.

``EPoint::IntegrateTargetEffectForObservable`` (``src/EPoint.cpp:1498``) runs the
existing yield integrator twice — once on :math:`\sigma`, once on the product
:math:`\sigma A_y` — and divides:

.. code-block:: cpp

   for (int i=1;i<=n;i++) sigma[i-1] = GetSubPoint(i)->GetFitCrossSection();
   IntegrateTargetEffect(configure);
   const double denominator = GetFitCrossSection();

   for (int i=1;i<=n;i++)
     GetSubPoint(i)->SetFitCrossSection(sigma[i-1] *
                                        GetSubPoint(i)->GetAnalyzingPower());
   IntegrateTargetEffect(configure);
   const double numerator = GetFitCrossSection();

   for (int i=1;i<=n;i++) GetSubPoint(i)->SetFitCrossSection(sigma[i-1]);
   SetFitCrossSection(std::fabs(denominator) > 0.0 ? numerator/denominator : 0.0);

Reusing the integrator rather than writing a second one is deliberate: the
quadrature, the straggling model and the energy-loss treatment are then
identical to the cross-section path by construction, and there is only one of
them to maintain. The cost is one extra integration per analyzing-power point,
which is negligible against the R-matrix evaluation itself.

This also required :math:`A_y` to be carried *beside* the cross section rather
than replacing it, since :math:`\sigma` is the weight. Only a non-sub-point
substitutes :math:`A_y` into the reported value, which keeps output files,
:math:`\chi^2` and plotting free of special cases.

The dilution this produces is dramatic and is a genuine prediction, not a bug.
For the thick gas target in ``tests/13N`` the sub-point grid runs from 0.42 to
1.54 MeV; Rutherford scattering makes :math:`\sigma` diverge at the low end where
:math:`A_y \approx 0`, and a resonant :math:`A_y \approx 0.8` averages down to
about :math:`10^{-6}`. Comparisons against thin-target data must use segments
with no target integration.

.. warning::

   ``<targetInt>`` assigns effects to segments **by index**, and
   ``<segmentsTest>`` shares that numbering with ``<segmentsData>``. Test
   segments 1 and 2 therefore silently inherit whatever the data segments 1 and 2
   declare. This is not specific to analyzing powers, but it is much more visible
   here, because the wrong answer is a factor of :math:`10^{6}` rather than a few
   percent.

Getting data to test against
-----------------------------

Baumann *et al.* publish no table of :math:`A_y`. Their results are contour
plots, plus six angular distributions in fig. 1 at
:math:`E_p = 1.618, 1.658, 1.708, 1.738, 1.758` and 1.779 MeV. Those six panels
are the only numerically recoverable data in the paper, so they were digitized.

The procedure (``tests/13N/digitize_fig1.py``) is:

#. Render page 3 at 600 dpi with ``pdftoppm``.
#. Locate the six panel frames by finding long horizontal and vertical dark runs.
#. Calibrate the axes on the *printed tick labels* rather than on the frame — the
   frame does not coincide with the axis limits. The x labels 0 and 150 degrees
   give 6.827 px/degree; the y ticks turn out to be spaced by 0.2 in
   :math:`A_y`, not the 0.25 one would guess from the 1.0 / 0. / −1.0 labels,
   and this was checked by locating the label glyphs themselves.
#. Remove the thin dashed :math:`A_y = 0` rule with a vertical morphological
   opening — the curve is about 12 px thick and survives, the rule is not.
#. Keep every elongated connected component, discarding compact ones, which
   drops the ``1.618 MeV`` label glyphs while keeping a curve that has
   fragmented.
#. Trace the curve column by column with a continuity constraint, and sample at
   10-degree intervals from 40° to 160°.

Two honest caveats attach to the result. First, what is traced is the *fitted
curve* of the paper's phase-shift analysis, which the plotted points lie on;
these are the published analysis sampled on a grid, not the raw measurements.
Second, the quoted uncertainty of 0.04 is **digitization error, not the
experiment's**. Baumann quote statistical errors below 0.004 — ten times smaller.
The 0.04 is roughly the drawn line width, 12 px against 466 px per unit of
:math:`A_y`, plus the calibration residual. Points where the trace could not
follow a steep section were dropped rather than interpolated, leaving 71 of a
possible 78.

Angles are centre-of-mass, as printed on the figure axis, which is what
observable 7 expects; energies are laboratory.

Fitting it
----------

The analyzing-power data live in ``tests/13N`` as segments 11--16, one per
energy, beside the capture and scattering data of the same compound nucleus.
Fitting :math:`A_y` *alone* -- every level fixed except the excitation energy and
proton width of the two resonances the energy range is sensitive to -- takes
:math:`\chi^2` from 236.8 to 85.2, or 1.20 per point and 1.27 per degree of
freedom (``tests/13N/README.md`` gives the recipe):

.. list-table::
   :header-rows: 1
   :widths: 26 22 26 22

   * - parameter
     - fitted
     - Baumann table 1
     - start
   * - :math:`3/2^-`  :math:`E_x`
     - 3.4954 MeV
     - 3.499 MeV
     - 3.5032 MeV
   * - :math:`3/2^-`  :math:`\Gamma_p`
     - 50.5 keV
     - 57 keV
     - 55.2 keV
   * - :math:`5/2^+`  :math:`E_x`
     - 3.5430 MeV
     - 3.546 MeV
     - 3.5453 MeV
   * - :math:`5/2^+`  :math:`\Gamma_p`
     - 53.6 keV
     - 50 keV
     - 49.0 keV

Excitation energies agree to 3–4 keV and widths to about 10%. That is what
digitized data can support, and no more should be claimed from it: with
uncertainties an order of magnitude looser than the measurement, the fit cannot
be expected to recover the published parameters more tightly. The point of the
exercise is that an R-matrix model fitted to :math:`A_y` alone lands on the same
two resonances, from the correct side, with sensible widths.

Differentiating it
------------------

The analyzing power has an exact adjoint, which matters because without one a
fit either falls back to numerical derivatives or -- worse -- uses the wrong
ones. Before this was written, an analyzing-power point fell through
``AMatrixFunc::PointAdjoint`` into the differential-cross-section branch, since
:math:`A_y` segments are flagged differential. It silently returned
:math:`\partial\sigma/\partial p` where :math:`\partial A_y/\partial p` was
wanted.

One structural fact makes the derivative easy. ``AddPathway`` is called in a
loop over exactly the ``(k, m)`` K-group and M-group indices that
``PointAdjoint`` uses for its T-matrix cotangents ``tBar``, so the forward and
reverse passes are the same loop and no index mapping is needed. And since
:math:`M` is *linear* in :math:`T`, the coefficient that ``AddPathway``
multiplies :math:`T` by **is** the derivative:

.. math::

   M_{s'\nu' s\nu} = \sum_{k,m} c^{(k,m)}_{s'\nu' s\nu}\, T_{k,m}
                     + \text{(Coulomb)},
   \qquad
   \frac{\partial A_y}{\partial T^{*}_{k,m}}
     = \sum_{s'\nu' s\nu} \overline{c^{(k,m)}_{s'\nu' s\nu}}\;
       \frac{\partial A_y}{\partial M^{*}_{s'\nu' s\nu}} .

The Coulomb amplitude depends only on energy, so it drops out.

With :math:`u_i = M_{s'\nu',+1/2}` and :math:`d_i = M_{s'\nu',-1/2}`, and
:math:`A_y = N/D` as in :eq:`ay`, the Wirtinger derivatives are

.. math::

   \frac{\partial N}{\partial u_i^{*}} = i\, d_i, \quad
   \frac{\partial N}{\partial d_i^{*}} = -i\, u_i, \quad
   \frac{\partial D}{\partial u_i^{*}} = u_i, \quad
   \frac{\partial D}{\partial d_i^{*}} = d_i,

so that

.. math::

   \frac{\partial A_y}{\partial u_i^{*}} = \frac{i\, d_i}{D}
       - \frac{N}{D^{2}}\, u_i ,
   \qquad
   \frac{\partial A_y}{\partial d_i^{*}} = -\frac{i\, u_i}{D}
       - \frac{N}{D^{2}}\, d_i .

``AmplitudeMatrix::AnalyzingPowerBar`` returns twice these, matching the
cotangent convention the rest of ``AMatrixFunc`` uses -- it finally takes
``Re(conj(bar) * dz/dp)`` with no further factor, so the bar must carry the
:math:`2` that a real function of a complex variable requires.
``PathwayAdjoint`` then walks the ``AddPathway`` loop and contracts. Everything
downstream -- T to :math:`U`, :math:`U` to the level matrix, level matrix to
:math:`E_\lambda` and :math:`\gamma` -- is the machinery that already existed.

Verified against central differences on all 71 analyzing-power rows of
``tests/13N`` and all 14 free parameters. Agreement is at :math:`10^{-9}` to
:math:`10^{-10}` for twelve of them. The two level-energy parameters, whose
derivatives are the largest, show :math:`4\times10^{-6}` at a step of
:math:`10^{-5}` -- and that residue falls as :math:`h^2`
(:math:`1.5\times10^{-5}`, :math:`1.5\times10^{-7}`, :math:`1.7\times10^{-9}`
for :math:`h = 10^{-4}, 10^{-5}, 10^{-6}`), which is finite-difference
truncation rather than an error in the adjoint. The capture-channel widths give
identically zero, correctly: they cannot affect elastic scattering.

A point whose :math:`A_y` is averaged over a target is deliberately not
supported. That quantity is a ratio of two integrals (Eq. :eq:`aytarget`), so
its derivative needs the quotient rule across sub-points and both integrals
differentiated. Rather than approximate it, such a point reports itself
unsupported, which makes the whole Jacobian unavailable and returns the fit to
numerical derivatives.

Using it
--------

In the GUI, choose *Analyzing Power* from the **Data Type** menu when adding a
data or test segment. It behaves like a centre-of-mass differential segment:
the angle fields are enabled, **Vary Norm?** is disabled because a
normalization factor means nothing for a ratio, and the Plot tab switches the
y-axis to linear and to *Cross Section* when such a segment is drawn, since a
logarithmic axis cannot show a quantity that goes negative.

By hand, the observable is code 7. Angles are centre-of-mass and energies
laboratory, as for any differential segment:

.. code-block:: text

   <segmentsData>
   1  1  1  1.5  1.9  40  160  7  1  0  5  0  0.005  0  data/baumann_ay.dat 0 0
   </segmentsData>

The data file carries ``E_lab  theta_cm  A_y  dA_y``. From Python:

.. code-block:: python

   from pyazr import azure2
   with azure2("13N.azr") as azr:
       ay = azr.calculate_analyzing_power_rwa(azr.params_rwa)

Because :math:`A_y` is reported in place of the cross section, everything
downstream — :math:`\chi^2`, output files, plotting, fitting, MCMC — works
without modification.

Two things to keep in mind. A normalization factor is meaningless for a ratio, so
leave ``varyNorm`` at 0. And use segments without target integration when
comparing against thin-target data, for the reason given above.

What is not done
----------------

* **Tensor observables** need a spin-1 projectile and rank-2 operators. The
  amplitude matrix carries all channel spins and projections already, so the
  missing part is the observable side, but it is not written.
* **Capture channels** need Seyler and Weller rather than Seyler.
* **Analytic derivatives through a target integration.** :math:`A_y` itself now
  has an exact adjoint (see below), but a point whose :math:`A_y` is averaged
  over a target is a ratio of two integrals and is not differentiated
  analytically. Such a point returns *unsupported*, which makes the whole
  Jacobian unavailable and falls the fit back to numerical derivatives -- coarse,
  but never wrong.
