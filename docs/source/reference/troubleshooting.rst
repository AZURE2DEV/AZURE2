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

No Standard R-Matrix Corresponds to the Levels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom**: With the Brune parameterization off (``--no-brune``, pyazr
``use_brune=False``) the run stops at the parameter transformation: "No
standard R-matrix (constant boundary conditions) corresponds to the levels of
J^pi = ...", then "Calculation was aborted." (pyazr: ``AZURE2Error`` on
opening the session).

**Cause**: The standard formalism needs the input (observed) energies and
widths transformed to formal ones, a generalized eigenvalue problem
:math:`N a = E M a` with
:math:`M = 1 - \sum_c \gamma_c \gamma'_c\, dS_c/dE` (Brune, Phys. Rev. C 66,
044611).  It has a real solution only when :math:`M` is positive definite;
widths too large for the channel radii -- typically broad background poles of
a Brune fit -- make it indefinite, and no standard R-matrix reproduces those
levels.  (Before 2026-11-01 the failure went unnoticed: the run crashed with
"z is not finite in log_Gamma".)  The warnings "Denominator less than zero
while transforming" printed just before name the levels whose own
:math:`1 - \sum_c \gamma_c^2\, dS_c/dE` is negative.

**Solution**: Use the Brune parameterization (the default), or reduce the
widths of the levels listed (or enlarge the channel radii).

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

**Symptom**: A reaction rate that disagrees with an integral of the extrapolated
cross section, or that changes when a level is moved by a parameter file.

**Cause**: The rate tabulates the cross section once on an energy grid (0.1 eV to
50 MeV in the centre of mass, 48 points per decade, plus a cluster of points from
1/32 to 64 widths around every level of the parameters being evaluated), splines
ln σ against ln E and integrates each temperature over that spline. A resonance
the parameters have moved outside the energy range is not integrated at all, and
a level with no width in any channel is not treated as a resonance.

**Solution**: Look at ``output/cross_section_grid.dat`` (E_cm, σ, S-factor), which
is exactly the grid the rate integrated, and at ``parameters.out`` for the observed
widths; compare with a fine mode-3 extrapolation integrated externally (the
regression test ``tests/reaction_rate`` does this with its ``truth.py``). Older
builds (before 2026-10-26) splined σ itself and set σ = 0 below 1 keV, which made
charged-particle rates about 1 % low at low temperature and neutron (1/v) rates
49 % low at T9 = 0.01, and placed the resonance clusters at the .azr level
energies rather than at the parameters being evaluated.

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
each particle pair for each :math:`J^\pi`. For elastic scattering, tick
**Include hard-sphere scattering for every** :math:`J^\pi` on the **Levels and
Channels** tab (on by default in new projects; off in projects saved before it
existed): every :math:`J^\pi` up to L_max then scatters whether or not a level
carries it, and the line under the box lists the :math:`J^\pi` that have none.
Older projects did this with fixed zero-width dummy levels, which still work.

File Dialogs on macOS
^^^^^^^^^^^^^^^^^^^^^

**Symptom**: on a Mac, a "Choose..." button (parameter file, data file,
directory) or File > Open does not let you pick a file. This was reported with
the native macOS file dialog and builds before November 2026.

**Cause**: the native macOS file panel did not work with the Qt 5 build of
the macOS application. Since November 2026 AZURE2 uses Qt's own file dialog
for every file and directory choice on macOS. It looks different from the
Finder dialog (no Finder sidebar or iCloud shortcuts) but works the same way.

**Solution**: use a current build. To try the native dialog anyway, start
AZURE2 from a terminal with ``AZURE2_NATIVE_DIALOGS=1`` in the environment:
``AZURE2_NATIVE_DIALOGS=1 /Applications/AZURE2.app/Contents/MacOS/AZURE2``
(``open -a`` does not pass the variable on). ``AZURE2_NATIVE_DIALOGS=0``
gives Qt's dialog on Linux and Windows as well.

Getting Help
------------

For questions and support, contact the developers at:

   azure@nd.edu
