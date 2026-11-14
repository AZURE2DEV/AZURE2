.. _example-ratio-to-rutherford:

Worked Example: α + :sup:`12`\ C Elastic Scattering as a Ratio to Rutherford
=============================================================================

Elastic angular distributions of charged particles are often published as the
ratio of the measured differential cross section to the Rutherford cross
section,

.. math::

   R(E,\theta) = \frac{d\sigma/d\Omega\,(E,\theta)}{d\sigma_R/d\Omega\,(E,\theta)},
   \qquad
   \frac{d\sigma_R}{d\Omega} = \frac{\eta^2}{4k^2\sin^4(\theta/2)}
   = \left(\frac{Z_1 Z_2 e^2}{4E}\right)^2 \frac{1}{\sin^4(\theta/2)},

with :math:`E` the c.m. energy (MeV), :math:`\theta` the c.m. angle,
:math:`k` the wave number (fm\ :sup:`-1`), :math:`\eta` the Sommerfeld
parameter and :math:`e^2 = \alpha\hbar c = 1.43996` MeV fm (the right-hand form
is the non-relativistic one).  The charge and solid angle of the experiment
cancel in the ratio, which is why such data are published this way.  About a
fifth of the light-nucleus elastic angular distributions in EXFOR are ratios
to Rutherford (code ``,DA,,RTH``).

In project format 3 such data are their own quantity, ``ratio_to_rutherford``
(``ratio_to_mott`` for identical particles).  The engine divides its elastic
cross section by its own Rutherford cross section at the same c.m. energy and
angle, with the same masses, constants and kinematics, and folds both with
the same apparatus kernel.  Multiplying the data by a Rutherford cross section
beforehand, the only route before November 2026, gives the same χ² only
without an energy shift and without a fold: a fitted shift moves
:math:`\sigma_R(E)` with the point in the ratio form but not in data that carry
:math:`\sigma_R(E)` frozen in (:math:`\sigma_R \propto E^{-2}`, so a 1 %
energy shift is a 2 % change), and a beam spread or a target thickness
averages numerator and denominator separately: the model is
:math:`\langle d\sigma/d\Omega\rangle / \langle d\sigma_R/d\Omega\rangle`, not
:math:`\langle R\rangle`.

The project
-----------

The example is the user story of ``docs/upgrade/project_format3.md``
section 9.2: α + :sup:`12`\ C elastic scattering at five lab angles (30°,
60°, 90°, 120°, 150°), 13 energies from 3.0 to 6.0 MeV (lab) in one file, with
a 3 % common normalization uncertainty.  The data are synthetic (the model
below times a fixed ±3 % pattern, divided by 1.02, with 4 % errors); the files
are ``tests/ratio_to_rutherford/models/aC_story.azr`` and
``models/data/aC_rth.dat``.  The levels are three :sup:`16`\ O states above the
α threshold (7.16192 MeV): 1\ :sup:`−` at 9.585 MeV (Γ\ :sub:`α` = 420 keV),
2\ :sup:`+` at 11.52 MeV (71 keV) and 3\ :sup:`−` at 11.60 MeV (800 keV), with
a channel radius of 5 fm and every other partial wave up to l = 4 as
hard-sphere scattering (``includeAllJPi=1``, ``maxOrbitalMomentum=4``).

The data set::

   <data>
   dataset aC_RTH
       reaction       12C(a,a)12C
       quantity       ratio_to_rutherford
       unit           1
       file           data/aC_rth.dat
       columns        E_lab[MeV] theta_lab[deg] value error
       norm           1 fit prior 1 3%
   </data>

The data file has the columns the ``columns`` field names: the lab energy
(MeV), the lab angle (deg), the ratio and its uncertainty.  The lab angle is
converted to the c.m. angle with α + :sup:`12`\ C elastic kinematics when the
file is read; the ratio is not converted, because it is the same at
corresponding lab and c.m. angles (numerator and denominator are differential
in the same solid angle).  An extrapolation of the same quantity on a finer
grid::

   <extrapolations>
   extrapolation aC_RTH_curve
       reaction       12C(a,a)12C
       quantity       ratio_to_rutherford
       unit           1
       energies       range 3.0 6.0 MeV lab step 0.05
       angles         list 30 60 90 120 150 deg lab
   </extrapolations>

What the engine understood (``AZURE2 --describe aC_story.azr``)::

   dataset aC_RTH (data set 1, line 28)
     ratio_to_rutherford: differential cross section, lab angle of 12C(a,a)12C
     divided by the Rutherford cross section at the same c.m. energy and angle, folded alike; the value is
     the same in both frames (not converted at load)
       12C(a,a)12C: entrance pair 1 (a+12C), exit pair 1 (a+12C)
     ...
     parameters:
       aC_RTH.norm (norm, fit)
       aC_RTH.shift.offset (energy shift, fixed)
       aC_RTH.shift.sqrt (energy shift sqrt, fixed)

and the model expression (``AZURE2 --print-predictions``,
``pyazr.predictions``): ``data aC_RTH: ratio(dsdo_cm(1->1), rutherford(1->1))``, the c.m.
elastic cross section of pair 1 over the Rutherford cross section of the same
pair at the same point.  The five angles share one normalization parameter,
``aC_RTH.norm``, with its 3 % prior.  If the five angles were five files, five
data sets would name one ``normalization`` record.

A calculation (``calculate``) gives χ² = 38.98 for the 65 points at the starting
normalization 1; a fit of the three energies, three widths and the
normalization (``fit --use-gradient``) ends at χ² = 18.52 with the
normalization at 1.0176 (normalization penalty 0.35), against the 1.02 the
data were made with.  The output file ``AZUREOut_aa=1_R=1.out`` holds the
model ratio and the data ratio in its columns 4 and 6; ``.extrap`` the
ratio on the extrapolation grid.

Identical particles: ratio to Mott
----------------------------------

For two identical particles of spin :math:`s` (α + α, :sup:`12`\ C +
:sup:`12`\ C, p + p) the Coulomb amplitude is symmetrised, and the reference
is the Mott cross section

.. math::

   \frac{d\sigma_M}{d\Omega} = \frac{\eta^2}{4k^2}\left[\frac{1}{\sin^4(\theta/2)}
   + \frac{1}{\cos^4(\theta/2)} + \frac{2(-1)^{2s}}{2s+1}\,
   \frac{\cos\!\left(\eta\ln\tan^2(\theta/2)\right)}{\sin^2(\theta/2)\cos^2(\theta/2)}\right],

symmetric about 90°.  The engine's elastic cross section of identical
particles is already the symmetrised one (for α + α,
:math:`|f_C(\theta)+f_C(\pi-\theta)+2f_N(\theta)|^2`, checked against an
independent partial-wave sum in ``tests/ratio_to_rutherford``), so such data
are ``quantity ratio_to_mott`` with ``reaction 4He(a,a)4He``.
``ratio_to_rutherford`` is refused there (``E_REACTION_NOT_ALLOWED``, the
message names ``ratio_to_mott``), and ``ratio_to_mott`` is refused for
different particles.  An EXFOR ``RTH`` entry of identical particles does not
say which reference was used; the entry's text decides (``pyazr.sources.nds``
asks, below).

What is refused
---------------

* A neutral pair (no Coulomb scattering) and an inelastic group (the
  lab-to-c.m. Jacobians of the reaction and of elastic Rutherford scattering
  differ): give ``dsigma_domega``.
* Folds the Coulomb reference does not follow, ``E_NOT_BUILT``: a
  Doppler-broadened target, electron screening, a detector's attenuation
  coefficients (they attenuate the Legendre terms of the nuclear amplitude;
  :math:`d\sigma_R/d\Omega \propto \theta^{-4}` has no finite Legendre
  series).  A beam spread, a target layer and a beam profile are folded: the
  Rutherford cross section is evaluated at the
  point's own sub-points and integrated by the point's own kernel.

From EXFOR
----------

``pyazr.sources.nds.fetch_exfor(id).format3_quantity()`` maps ``,DA,,RTH`` and
``PAR,DA,,RTH`` of an elastic reaction to ``ratio_to_rutherford`` and
``,DA,,MOT`` of identical particles to ``ratio_to_mott``, and refuses, with
the reason, ``RTH`` of identical particles (pass
``quantity="ratio_to_mott"`` once the entry says the reference is Mott),
``MOT`` of different particles (the relativistic point-nucleus cross section
of the EXFOR dictionary) and an inelastic group.
``ExforData.to_dataset(data_dir, name, norm=...)`` writes the ratios as they
are and returns the ``dataset`` record for ``<data>``.  ``to_azr`` refuses a
ratio to Rutherford.

Checks
------

``tests/ratio_to_rutherford/check.sh`` and ``tests/pyazr/ratio_to_rutherford_test.py``:
with every width zero the ratio equals :math:`|f_C+f_{HS}|^2/|f_C|^2` of an
mpmath partial-wave sum to 4×10\ :sup:`−10` (α + :sup:`12`\ C, c.m. and lab
angles; α + α with the Mott denominator); without shift or fold the χ² of ratio
data equals, to 9×10\ :sup:`−15`, that of the same data multiplied by the
Rutherford cross section and given as ``dsigma_domega``; with a 5 keV shift
the two differ by the amount the independent script computes (0.400 on
χ² ≈ 12.5); with a 0.2 MeV (lab) Gaussian beam spread the model is the
ratio of the two folds to 2×10\ :sup:`−6`, while the fold of the ratio differs
by up to 7×10\ :sup:`−3`; the gradient equals central differences.
