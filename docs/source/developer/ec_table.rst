External-capture energy tables
==============================

Status: on by default since October 2026 (``--no-ec-table`` turns it off).
Code: ``include/ECTable.h``, ``src/ECTable.cpp``; callers
``EPoint::CalculateECAmplitudes`` and ``EPoint::GetECAmplitudeWithShift``.
Measurements: ``12C+a_onefile/10-20-26_EC_table_timing`` in the fit archive.

Why
---

An external-capture amplitude needs a radial integral from the channel radius
to infinity of a Coulomb function of the entrance channel against the bound
state's Whittaker function (``ECIntegral``, GSL ``qagiu``, about 1 ms).  It was
computed directly for every data point and every sub-point, once at start-up:

- resolution folding or target integration of capture data multiplies the
  number of integrals by the number of sub-points.  On the
  :sup:`12`\ C(α,γ) model, a 12.7 keV Gaussian resolution on the capture data
  took the integrals from 30,558 to 894,361 and the start-up from 15 s to
  243 s, 95 % of it in the integrals;
- the reaction rate integrates afresh at every energy of its grid (2418 for
  :sup:`12`\ C(α,γ\ :sub:`0`), 28 s, about 88 % of it in the Coulomb functions
  of the radial integrand);
- an energy shift moves points away from the energies their integrals were
  computed at.  The amplitudes were then interpolated linearly between the
  energies of the segment's own points and held flat beyond its ends.  On the
  :sup:`12`\ C(α,γ) cascades a shift of −30 keV put the cross section 9–15 %
  off at the segment ends and 1–8 % off in the large gaps between points; on
  the 13N test a −30 keV shift of the LUNA data put the lowest points off by
  a factor of 57 (the amplitude falls steeply and was held constant).

Yet every one of those integrals is a sample of one smooth function of energy
per capture term.

What is tabulated
-----------------

A *term* is everything the integral depends on except the energy: the final
pair (charges, masses, spins, g-factors, channel radius, threshold), the
entrance and final orbital and spin quantum numbers, the total spins, the
multipolarity and radiation type, channel capture or not, the final level's
energy, and the options that change the integrand (long-wavelength
approximation, GSL Coulomb functions, rate mode, the Coulomb table, the hybrid
potential's tag).  The key holds these values, not positions in a model, so
several models or pyazr sessions share one process-wide table and an entry can
never belong to a different calculation.

The energy variable is :math:`u = \ln E_\mathrm{out}`, with
:math:`E_\mathrm{out}` the c.m. energy above the final pair's threshold; nodes
are uniform in :math:`u`, 50 per decade by default.  Before tabulation the
steep Coulomb behaviour is divided out, exactly and cheaply (the Coulomb table
answers both at the channel radius):

- an ordinary term is divided by :math:`\sqrt{P_{l_i}(E_\mathrm{out})}` of the
  final pair -- the integral falls like the Coulomb penetration;
- a channel-capture term is divided by the final pair's hard-sphere phase
  factor, which it carries and which turns with energy.

The amplitude at any energy is then the 8-point Lagrange interpolant of the
stencil around :math:`u`, times the factor at that energy, times the exact
entrance-channel factors (Coulomb and hard-sphere phases, photon-energy power)
that ``EPoint`` always computed.  Sub-threshold entrance channels
(:math:`E_\mathrm{out} \le 0`, the Whittaker-Whittaker branch) are computed
directly.

Nodes are computed on demand, the first time an energy within reach of them is
asked for, and kept for the life of the process.  Any number of OpenMP threads
may ask at once: a thread claims the missing nodes of its stencil under the
term's lock, computes them outside it, and waits for nodes another thread is
computing (``std::shared_future``).  A node is converged to a relative
tolerance of 1e-10 (``ECIntegral::SetRelativeTolerance``); where the quadrature
cannot reach that it is computed at the direct integrals' 1e-4, and a node that
fails altogether sends the energies that need it to the direct integral.

Shifted points
--------------

``EPoint::GetECAmplitudeWithShift`` returns the stored amplitude while the
point is at the energy it was computed at; at any other energy it evaluates the
pathway at the point's current energy (``ECPathwayAmplitude``, the same code as
start-up) and keeps the result until the energy moves again.  The final level's
energy is the one in force at the first shifted evaluation, so shifted and
unshifted segments stay consistent where an MCMC moves that level.  With
``--no-ec-table`` the old linear interpolation through ``ECAmplitudeCache`` is
used.

Accuracy
--------

``AZURE2_EC_TABLE_CHECK=1`` computes every tabulated integral a second time,
directly and converged to 1e-10, and the run log reports the largest relative
difference and, per term, where it occurs.  On the 13N test model
(8 terms, :math:`E_\mathrm{out}` = 0.068–2.34 MeV):

=====================  =================================
nodes per decade        largest relative difference
=====================  =================================
25                      3.1e-9
50 (default)            1.6e-10
100                     8.4e-11
=====================  =================================

from 50 nodes per decade on, the difference is the convergence of the reference
integrals themselves.  Against the direct integrals of ``--no-ec-table``
(converged only to 1e-4 relative, ``qagiu``'s setting) the amplitudes differ by
1e-4 of each term's largest value -- that is the old integrals' scatter, which
the table no longer has; where a term passes near zero this is up to 0.9 % of
the small amplitude there, and 1e-6 of the cross section.

Cost
----

13N test model with its target integration (369,696 integrals): 1400 node
integrals and 11 s for the whole run on two threads, against 177 s with ``--no-ec-table``
(the same total :math:`\chi^2`, 336,454, and every segment to the printed digits).
The 12C(α,γ) measurements are in the fit archive's
``12C+a_onefile/10-20-26_EC_table_timing/readme``.

Tests
-----

``tests/ec_shift/check.sh``: shifted capture segments against the same data
moved in energy (identical), the shift's effect, the table against converged
integrals (1e-7) and against ``--no-ec-table`` (1e-3).  Every physics suite of
``tests/run_tests.sh`` with external capture now runs on the table.
