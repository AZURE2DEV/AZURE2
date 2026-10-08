Experimental Effects Tab
========================

The **Experimental Effects** tab is used to apply corrections for experimental
effects such as beam energy loss in targets, beam energy resolution, detector
geometry, and energy straggling.

.. warning::

   The target integration and beam resolution convolution routines implemented
   in AZURE2 are basic and may not cover all experimental situations. The
   developers strongly recommend evaluating these routines on a case-by-case
   basis. Modifications to the source code may be necessary for specific
   experimental setups.

.. warning::

   There are known issues when using both the target convolution and target
   integration routines simultaneously. Exercise extreme caution if combining
   these options.

Overview
--------

The experimental effects are modeled as:

.. math::

   F(E_0) = \int_{E_0 - \Delta}^{E_0} \frac{\sigma(E')}{\epsilon(E')}
   \int_{-\infty}^{+\infty} g(E - E_0) \, dE' \, dE

where :math:`\sigma(E')` is the true cross section, :math:`g(E' - E)` is a
spreading function representing the beam energy distribution, and
:math:`\epsilon(E')` is the stopping cross section.

The spreading function is a Gaussian:

.. math::

   g(E - E_0) = \frac{1}{\sqrt{2\pi}\,\sigma_b}
   \exp\left(-\frac{(E - E_0)^2}{2\sigma_b^2}\right)

Managing Experimental Effects
-----------------------------

- Click **+** to create a new experimental effects entry.
- Select an entry and click **-** to delete it.
- Double-click to edit.

.. note::

   Experimental effects entries apply to both data segments and calculation
   segments simultaneously. Remember to enable or disable them appropriately
   depending on the calculation being performed.

Add Experimental Effect Dialog
------------------------------

Associated Segments
^^^^^^^^^^^^^^^^^^^

The **Segments List** field specifies which calculation segments (from the
**Segments** tab) this experimental effect applies to. Enter segment numbers
using:

- Comma-separated values: ``3,4,5,7,8,9``
- Ranges: ``3-9``
- Combinations: ``3,6,7-14``

Integration Points
^^^^^^^^^^^^^^^^^^^

The number of points used for numerical integration when computing energy
convolution or target integration. The required number depends on how rapidly
the cross section changes with energy. Adjust using the spinner or enter a
value directly.

Gaussian Energy Convolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Check **Include Gaussian Convolution** to convolve the calculated cross section
with a Gaussian beam energy distribution.

**Sigma** (MeV)
   The standard deviation :math:`\sigma_b` of the Gaussian convolution function
   (not the FWHM; FWHM :math:`= 2.355\,\sigma_b`). Although beam resolution is
   typically of order keV, the value must be entered in **MeV** (e.g., ``0.001``
   for 1 keV), in the laboratory frame; AZURE2 converts it to the centre of mass.

Target Integration
^^^^^^^^^^^^^^^^^^

Check **Include Target Integration** to account for beam energy loss in the
target. The yield is computed as

.. math::

   Y(E_0) = \frac{1}{N} \int_{E_0 - \Delta}^{E_0} \frac{\sigma(E')}{\epsilon(E')}\, dE',
   \qquad \Delta = \epsilon(E_0)\, N ,

where :math:`N` is the areal density entered below and :math:`\epsilon(E)` the
stopping cross section entered below.

**Active Density** (atoms/cm\ :sup:`2`)
   The areal density of the active target material (the nuclei producing the
   reactions of interest in a mixed-material target).

**Stopping Cross Section**
   The effective stopping cross section must be entered as a continuous function
   of energy using a parameterized equation:

   - The variable ``y`` represents the stopping cross section, in
     MeV cm\ :sup:`2` per *atom* (or per whatever entity the density counts --
     see below). Tabulated values in eV per 10\ :sup:`15` atoms/cm\ :sup:`2`
     (SRIM, Ziegler) are converted by dividing by ``1e21``.
   - The variable ``x`` represents the laboratory energy in MeV.
   - Parameters are labeled ``a0``, ``a1``, ``a2``, etc.

   **Example** -- a second-order polynomial with 3 parameters::

      y = a0 + a1*x + a2*x^2

   Set the **Number of Parameters** to ``3`` and enter the values for ``a0``,
   ``a1``, and ``a2`` in the table.

   AZURE2 also provides tools to look up stopping powers by element or compound
   formula (**Fetch from ERYA**).

**Compound and mixed targets: one stopping power is enough, but it must
be normalized to the same thing the density counts.**

   The stopping cross section and the areal density enter the calculation
   only through their product (the energy loss :math:`\Delta = \epsilon N`)
   and through :math:`1/(N\epsilon)` in the integrand, so any pair
   :math:`(\epsilon, N)` that refers to the *same entity* gives the same
   yield. You do not need a separate "compound" and "effective" stopping
   power; you need the one that matches your :math:`N`. For a target
   A\ :sub:`x`\ B\ :sub:`y` in which A is the reacting nucleus, Bragg's rule
   gives the stopping cross section per molecule
   :math:`\epsilon_{\rm mol} = x\,\epsilon_A + y\,\epsilon_B`, and the three
   consistent choices are

   ============================  =================================================  ==========================
   entity                        stopping cross section :math:`\epsilon`            density :math:`N`
   ============================  =================================================  ==========================
   active atom (recommended)     :math:`\epsilon_A + (y/x)\,\epsilon_B`              active atoms / cm\ :sup:`2`
   molecule                      :math:`x\,\epsilon_A + y\,\epsilon_B`               molecules / cm\ :sup:`2`
   average atom                  :math:`(x\,\epsilon_A + y\,\epsilon_B)/(x+y)`       *all* atoms / cm\ :sup:`2`
   ============================  =================================================  ==========================

   The mistakes are the mixed pairings: the pure-element :math:`\epsilon_A`
   with the active density (the energy loss is missing the
   :math:`(y/x)\,\epsilon_B` term -- for TiN or Ta\ :sub:`2`\ O\ :sub:`5` a
   large fraction of it), or a per-average-atom value with the active
   density (the loss is too small by :math:`(x+y)/x`).

   **Fetch from ERYA** with a compound formula does this for you: choose the
   reacting nucleus in the **Active Element** box next to the formula (it is
   filled from the formula; the first element is preselected) and the fetched
   equation is the per-active-atom sum, ready to pair with the Active
   Density. The **Calculate ΔE** preview uses the same normalization.

   .. warning::

      Before September 2026 the ERYA compound fetch returned the
      *average-atom* value (each element weighted by its stoichiometry divided
      by the *total* stoichiometry) with no way to choose otherwise. A
      ``<targetInt>`` entry written by an older GUI for a compound target is
      therefore consistent only with the areal density of **all** atoms; with
      the active density it makes the target too thin by :math:`(x+y)/x`.
      Re-fetch it, or multiply its equation by that factor.

   If the experiment reports the energy loss in the target directly (for
   example "the target thickness corresponded to 1.5 keV at the resonance"),
   that number is more reliable than a stoichiometry and a nominal thickness:
   choose :math:`N` so that :math:`\epsilon(E_0)\,N` reproduces the quoted
   loss at the energy of interest.

Restricting an Effect to Energy Ranges
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default an experimental effect covers every point of the segments it lists.
Three optional controls in the *Add/Edit Experimental Effect* dialog refine
this without splitting the data into artificial segments:

- **Apply in Energy Ranges** -- a comma-separated list of laboratory-energy
  windows, for example ``0.42-0.61,1.20-1.35`` (MeV).  Points outside every
  window are computed as ordinary points; leave the field empty to cover the
  whole segment.
- **Blend Width** -- the width (MeV) of a smooth transition at each window
  edge.  With a hard edge (``0``) the modelled curve can show a small step
  where the convolution switches off; a positive width ramps continuously
  between the convolved and the unconvolved curve.
- **Auto Tolerance** -- a relative tolerance that makes the decision
  automatic: at each point the code estimates how much the effect would change
  the observable and skips the integration where the change is below the
  tolerance.  Any discontinuity this introduces is bounded by the tolerance,
  and smooth regions stop paying for integration they do not need.  It can be
  combined with explicit ranges or used on its own.

In the ``.azr`` file these appear as optional tokens at the end of the
``targetInt`` line, for example ``"1.95-2.55" 0.12 0.002``; files that do not
use them are written exactly as before, and remain readable by older versions
of AZURE2.

Beam-profile kernel (photodissociation in a broad γ beam)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Gaussian convolution above is centred on each data point.  A
measurement in a broad, asymmetric photon beam whose reaction energy is
reconstructed event by event (a TPC at HIγS, for example) averages the cross
section differently: the beam profile is an *absolute* energy distribution,
and each data point is the slice of that beam that the detector reconstructed
inside an energy window :math:`[a, b]` with a Gaussian resolution
:math:`s`.  For such data the model is

.. math::

   \langle\sigma\rangle = \frac{\int G(E)\,W(E)\,D(E)\,\sigma(E)\,dE}
                                 {\int G(E)\,W(E)\,D(E)\,dE},
   \qquad
   W(E) = \tfrac{1}{2}\left[\operatorname{erf}\frac{b-E}{s\sqrt2}
                              - \operatorname{erf}\frac{a-E}{s\sqrt2}\right]

with :math:`G(E)` the beam profile, a weighted sum of skewed Gaussians

.. math::

   G(E\,|\,\xi,\omega,\alpha) = \frac{1}{\omega\sqrt{2\pi}}
   \exp\left[-\frac{(E-\xi)^2}{2\omega^2}\right]
   \left[1+\operatorname{erf}\frac{\alpha\,(E-\xi)}{\omega\sqrt2}\right],

and :math:`D(E)` the detailed-balance factor of the inverse reaction
relative to its value at the point's own energy (``dbFlag`` = 1), so that a
capture cross section is averaged the way the photodissociation measurement
averaged it.  The formalism follows Haverson (2026), appendix A.

The effect is written as a trailing block of the ``targetInt`` line, after
the optional straggling and energy-range tokens::

   beamprofile N  xi_1 omega_1 alpha_1 w_1  ...  xi_N omega_N alpha_N w_N  s  nCut  dbFlag

All energies are **laboratory energies of the entrance channel** in MeV,
like every other energy AZURE2 reads (a photon-beam profile has to be shifted
by the Q-value and converted first).  ``nCut`` > 0 zeroes each component
outside its mean ± ``nCut`` standard deviations (0 = full profile).  The
window :math:`[a, b]` of each point is read from optional **columns 5 and 6
of the data file** (lab MeV); a point without them is averaged over the
whole beam.  A segment energy shift moves the window with the point.  The
sub-point grid covers the beam profile, narrowed to the window plus four
resolution widths, on the adaptive grid described for the other effects.
The kernel applies to any observable of the segment, including the E1- and
E2-only capture components, and is supported by the analytic gradient.
When the beam covers a resonance much narrower than the base sub-point step,
the adaptive grid anchors a fine lattice on it; give the effect a resonance-width
multiplier of at least 20 (the tenth field before ``beamprofile`` in the line,
``resonance_width_multiplier`` in pyazr) so that the lattice also covers the
Lorentzian tails -- with 10 widths, linear interpolation of the tails biased a
slice dominated by a 0.6 keV 2\ :sup:`+` by 18 %, with 20 widths the result
agreed with a dense reference integration to 0.04 %.

In the GUI the effect is set up under **Include Beam Profile** in the
*Add/Edit Experimental Effect* dialog: a component count, a table of
``xi``/``omega``/``alpha``/weight rows, the detector resolution sigma, the
truncation, and the **Weight by detailed balance** box.  The tab reads and
writes the tokens, so a project can be loaded, edited and saved without
losing them.  ``pyazr.AzrModel.add_target_effect`` writes the same line
from a script.

User-defined numerical resolution function (SAMMY UDR)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Time-of-flight facilities increasingly publish their resolution function
as a *numerical* distribution of flight-time delays, obtained from a
simulation of the neutron source and the beam line, rather than as an
analytic parametrization: n_TOF at CERN distributes the response of each
experimental area this way, in the format of SAMMY's **user-defined
resolution function** (UDR).  AZURE2 reads the same files.

The file holds, for each of a set of neutron energies, the probability
density :math:`R(\tau;E)` that a neutron of energy :math:`E` is detected
with a delay :math:`\tau` (microseconds) relative to its true flight time
over the nominal path.  A positive delay -- moderation in the target, a
longer path in the detector -- makes the neutron look slower than it is,
so the true energies :math:`E'` contributing to a data point at nominal
energy :math:`E_0` are those with

.. math::

   t(E') = t(E_0) - \tau, \qquad t(E) = \frac{L}{v(E)},

:math:`L` the nominal flight path and :math:`v` the (relativistic) velocity.
The fitted value of the point is the exact integral of the piecewise-linear
tabulated kernel against the cross section interpolated linearly between
the sub-points of the integration grid,

.. math::

   \langle\sigma\rangle(E_0) = \frac{\sum_i w_i\,\sigma(E_i)}{\sum_i w_i},
   \qquad
   w_i = \int R(\tau;E_0)\,\phi_i(\tau)\,d\tau ,

with :math:`\phi_i` the hat function of sub-point :math:`i` on the delay
axis -- the scheme of SAMMY's ``Udr_Resb``.  Before use the tabulated
function at the point's energy is obtained by linear interpolation between
the two bracketing tabulated energies, re-centred on its centroid (so that
the mean delay, which is part of the nominal flight path, is removed --
SAMMY does the same; it can be switched off), and convolved with a Gaussian
**burst width** and a rectangular **channel width** in time, both given in
nanoseconds.

File format (SAMMY Table III C5.1): any descriptive text, a line containing
at least five hyphens, then blocks of an energy line (eV; anything after the
number is ignored) followed by ``time  density`` pairs in microseconds and
per microsecond, blocks separated by blank lines.  The normalisation of the
tabulated density is irrelevant; the kernel is normalised numerically.  Only
the energies within a decade of the segment's range are read, so the 146 MB
n_TOF files load in a few seconds.

The effect is written as a trailing block of the ``targetInt`` line, after
the optional ranges and beam-profile tokens::

   udr "<file>"  L  burstFWHM  channel  centred

with ``L`` in metres, the widths in ns (0 = none) and ``centred`` 1 or 0.
A *negative* channel value ``-n`` stands for ``n`` bins per decade of
energy, the logarithmic binning of n_TOF data, and gives each point a
channel of :math:`t\,\ln 10/(2n)` at its own flight time (SAMMY manual
eq. III C3 a.14).
The flight path must be the one the experiment used to turn flight times
into the energies of the data file.

Two points about the n_TOF files in particular (n_TOF-PUB-2021-001).  The
proton burst is already folded into them -- the Transport Code that makes
the histogram smooths in time with the nominal 7 ns r.m.s. pulse, and the
tabulated delay distribution keeps a 6-7 ns r.m.s. width up to the highest
energies, where moderation contributes nothing -- so ``burstFWHM`` should
be 0 for them, and a separate file exists per pulse type (7 ns high
intensity, 14 ns low intensity).  Check any new file the same way: the
width of its highest-energy blocks is the burst it contains.  And the EAR2
file is specific to one experiment (flight path, sample size and alignment
change its shape), produced for that measurement with the Transport Code
and ``RF2sammy``; the n_TOF validation of EAR2 does not let SAMMY re-centre
the function but fits the flight path instead, which here is ``centred 0``
with ``L`` the length the data's time-to-energy conversion used.  The sub-point grid of each point covers
the kernel's support in true energy, on the adaptive grid described for the
other effects, so the usual advice on ``resonance_width_multiplier`` and
``points_per_width`` applies.  The kernel applies to any observable of the
segment and is supported by the analytic gradient.

In the GUI, check **Include User-Defined Resolution Function** in the
*Add/Edit Experimental Effect* dialog and give the file (``Browse...``),
the flight path, the burst and channel widths and the centring choice.
``pyazr.AzrModel.add_target_effect(..., udr_file=, flight_path=,
udr_burst=, udr_channel=, udr_centred=)`` writes the same line from a
script.

Straggling
^^^^^^^^^^

Check **Include Straggling** to account for energy straggling of beam particles
in the target. Enter the straggling coefficient in the provided field.

Attenuation Coefficients (Q-Coefficients)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Attenuation coefficients correct for the finite solid angle of detectors in
close geometry, following the method of M. E. Rose, *Physical Review* **91**,
610 (1953).

The angular distribution is corrected as:

.. math::

   W(\theta) = \sum_{i=0}^{\infty} a_i \, Q_i \, P_i(\cos\theta)

where :math:`Q_i` are the attenuation coefficients.

Set the number of coefficients using the spinner and enter the :math:`Q_i`
values in the table (default value is 1.0 for each).
