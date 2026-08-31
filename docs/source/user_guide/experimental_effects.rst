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
   The full width at half maximum of the Gaussian convolution function. Although
   beam resolution is typically of order keV, the value must be entered in **MeV**
   (e.g., ``0.001`` for 1 keV).

Target Integration
^^^^^^^^^^^^^^^^^^

Check **Include Target Integration** to account for beam energy loss in the
target.

**Active Density** (atoms/cm\ :sup:`2`)
   The areal density of the active target material (the nuclei producing the
   reactions of interest in a mixed-material target).

**Stopping Cross Section**
   The effective stopping cross section must be entered as a continuous function
   of energy using a parameterized equation:

   - The variable ``y`` represents the stopping cross section.
   - The variable ``x`` represents the energy.
   - Parameters are labeled ``a0``, ``a1``, ``a2``, etc.

   **Example** -- a second-order polynomial with 3 parameters::

      y = a0 + a1*x + a2*x^2

   Set the **Number of Parameters** to ``3`` and enter the values for ``a0``,
   ``a1``, and ``a2`` in the table.

   AZURE2 also provides tools to look up stopping powers by element or compound
   formula.

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
