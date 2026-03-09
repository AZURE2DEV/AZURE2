Troubleshooting
===============

While AZURE2 is fairly stable, errors can occur in certain situations. This
page describes the most common issues and their solutions.

Common Issues
-------------

Kinematic Limits
^^^^^^^^^^^^^^^^

**Symptom**: GSL error or hard crash during a pure calculation or extrapolation.

**Cause**: An energy or angular range was specified that is not kinematically
allowed, especially when extrapolating to very low energies.

**Solution**: Check the energy and angle ranges in your segments. Ensure they
are physically meaningful for the reaction being studied.

Threshold Effects
^^^^^^^^^^^^^^^^^

**Symptom**: GSL errors when calculating near the threshold of a reaction
channel.

**Cause**: When multiple entrance channels are defined, the code attempts to
calculate the penetrability for each channel. Very close to a threshold, the
penetrability can become extremely small, causing numerical issues.

**Solution**: Exclude data points that are very close to the threshold. Note
that this issue is exacerbated by target integration, convolution, and reaction
rate routines that create sub-points near the threshold.

Unconstrained Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom**: Crashes during fitting, especially with the Brune parameterization
disabled.

**Cause**: One or more fit parameters are very poorly constrained by the data,
causing them to be varied to unphysical values. The parameter transformation
code then fails.

**Solution**:

1. Use the **Brune parameterization** (recommended and more numerically stable).
2. Fix suspected unconstrained parameters and re-run the fit to identify the
   problematic one(s).
3. Add parameter limits in the **Fitting** tab.

MINOS Crashes
^^^^^^^^^^^^^

**Symptom**: The MINOS error analysis routine crashes.

**Cause**: MINOS is very sensitive to unconstrained parameters. A crash
usually indicates that the fit is not robust.

**Solution**: Running MINOS is itself a good test of fit quality. If it crashes,
fix different parameters to isolate the unconstrained ones. A successful MINOS
run is a strong indicator of a robust fit.

Reaction Rate Integration Failures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom**: Crashes or incorrect results during reaction rate calculation.

**Cause**: The numerical integration encounters a threshold or a very narrow
resonance that cannot be resolved by the adaptive integration routine.

**Solution**: For narrow resonances (total width less than approximately 1 keV),
calculate the reaction rate contribution separately using the narrow-resonance
approximation. For threshold issues, create a fine energy extrapolation and
perform the numerical integration externally (e.g., Simpson's rule).

Hanging Calculations
^^^^^^^^^^^^^^^^^^^^

**Symptom**: AZURE2 runs indefinitely without producing output or crashing.

**Cause**: Incorrectly formatted input data files.

**Solution**: Verify that data files are properly formatted with exactly four
whitespace-delimited columns of numerical data. Check for extra blank lines,
non-numeric characters, or missing columns.

Missing Channels
^^^^^^^^^^^^^^^^

**Symptom**: Unexpected physics results or crashes, especially for scattering
or external capture calculations.

**Cause**: The maximum orbital momentum (L_max) in the **Levels and Channels**
tab is set too low, resulting in missing channels.

**Solution**: Increase L_max and verify that at least one channel exists for
each particle pair for each :math:`J^\pi`. Consider adding dummy levels to
include all necessary hard-sphere phase shifts.

Getting Help
------------

For questions and support, contact the developers at:

   azure@nd.edu
