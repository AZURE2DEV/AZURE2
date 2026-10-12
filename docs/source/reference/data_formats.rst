Data Formats
============

Input Data Files
----------------

Experimental data files are plain text files with four whitespace-delimited
columns of real numbers (decimal or scientific notation):

.. list-table::
   :widths: 15 85
   :header-rows: 1

   * - Column
     - Description
   * - 1
     - Forward-kinematic laboratory frame **energy** (MeV)
   * - 2
     - Forward-kinematic laboratory frame **angle** (degrees)
   * - 3
     - Forward-kinematic laboratory frame **cross section** (barns or barns/sr)
   * - 4
     - Forward-kinematic laboratory frame **cross section uncertainty** (barns or barns/sr)

Example data file:

.. code-block:: text

   0.500  90.0  1.234e-03  5.6e-05
   0.600  90.0  2.345e-03  7.8e-05
   0.700  90.0  3.456e-03  9.0e-05
   0.800  90.0  4.567e-03  1.2e-04

Notes on Data Files
^^^^^^^^^^^^^^^^^^^


- Blank lines, and lines whose first non-blank character is ``#``, are
  comments and are ignored, so a data file may carry a header describing its
  provenance and units. Any other line must have at least the four numeric
  columns; a line that does not parse is reported with its file, line number
  and text, and the segment is dropped with the usual *Could Not Fill Segment*
  warning rather than silently skipped -- a mangled row is a corrupt file, not
  a comment.
- **Angle column**: must always be present, even for angle-integrated data. In
  that case the angle is a dummy value and is not used in the calculation.
- **Sorting**: data may be entered in any order. However, for plotting purposes
  it is useful to sort by energy (for excitation curves) or by angle (for
  angular distributions).
- **Frame**: all input data must be in the **laboratory frame** with
  **forward kinematics** (light particle as projectile, heavy particle as
  target).
- **Units**: energies in MeV, angles in degrees, cross sections in barns
  (angle-integrated) or barns/sr (differential).
- **Format**: columns can be separated by spaces or tabs.

Phase Shift Data
^^^^^^^^^^^^^^^^

Phase shift data follows the same four-column format, but the third column
contains the phase shift value (in degrees) and the fourth column contains its
uncertainty. The convention is :math:`-90° < \theta_\text{lab} < 90°`.

Mixture Table Files
^^^^^^^^^^^^^^^^^^^

A segment on a natural-abundance or compound target can add the cross sections
of the other isotopes as fixed tables (the ``mixture`` block of the segment
line; see :doc:`../user_guide/segments`). A table is a two-column text file:

.. code-block:: text

   # 22Ne(n,tot) ENDF/B-VIII.0, retrieved 2026-10-12
   # E_lab[MeV]  sigma[b]
   1.0000000000e-01  2.8541000000e+00
   1.0100000000e-01  2.8536000000e+00
   ...

- **Energy**: the lab energy of the projectile in MeV, as in a data file,
  strictly increasing. AZURE2 interpolates linearly between rows, so the grid
  must resolve the resonances of the tabulated isotope; it must also extend
  beyond the data range by the width of the segment's resolution function and
  by its energy shift, or the run stops and reports the range needed.
- **Cross section**: in the frame AZURE2 compares the segment in -- barns for
  angle-integrated and total cross sections, c.m. barns/sr for a differential
  segment.
- Lines starting with ``#`` and blank lines are ignored; at least two rows are
  required.

``pyazr.extrap_to_lab_table`` writes such a table from another AZURE2
project's ``AZUREOut_*.extrap`` (run with its ``<targetInt>`` lines inactive,
so the table is a bare cross section), converting its c.m. energies to the lab
with that project's entrance-pair masses. ``pyazr.mixture.endf_to_table``
writes one from the point-wise MF3 section of an ENDF-6 file, and
``endf_total_to_table`` does the same for the total cross section of a file
whose MF2 holds a Breit-Wigner (SLBW/MLBW) resolved-resonance range,
reconstructing that range at 0 K (``breit_wigner_total``, ENDF-102 App. D)
and joining it to the point-wise part above; Reich-Moore and R-matrix-limited
ranges need a processing code (NJOY, PREPRO).

Temperature Files
^^^^^^^^^^^^^^^^^

For reaction rate calculations, a temperature file is a single-column text file
listing temperatures in GK (gigakelvin):

.. code-block:: text

   0.01
   0.05
   0.10
   0.50
   1.00
   2.00
   5.00
   10.0

Project Files (.azr)
--------------------

The AZURE2 Input File (``.azr``) stores the complete project configuration in a
structured text format. While the format is not formally documented, it is a
plain text file that advanced users can modify directly for batch processing.

The file contains sections for:

- Configuration settings (formalism, directories, runtime options)
- Particle pair definitions
- Level and channel information
- Data segment definitions
- Test segment definitions
- Experimental effects configuration
- Fitting parameter settings
- MCMC settings (if enabled)

The ``<config>`` block holds a fixed sequence of lines (formalism, output and
checks directories, the check-file switches) followed by optional ``key=value``
lines, which a reader that does not know a key skips:

``includeAllJPi=1``
   Hard-sphere scattering in every :math:`J^\pi` the particle pairs allow,
   whether or not a level carries it (see :doc:`../user_guide/levels_channels`).
   Absent or ``0``: only the :math:`J^\pi` of active levels enter the calculation,
   as in projects saved before the key existed.

``maxOrbitalMomentum=N``
   The orbital-momentum limit of that enumeration, i.e. the **Maximum Orbital
   Momentum** of the Levels and Channels tab. Absent: the largest :math:`l` of
   any particle channel in the ``<levels>`` block.

A ``<targetInt>`` line may end with keyword blocks (``beamprofile ...``,
``udr ...``, ``sigmalab``). ``sigmalab`` marks an energy-dependent beam-spread
expression as laboratory in, laboratory out (see
:doc:`../user_guide/experimental_effects`); an expression without it is read in
the pre-October-2026 centre-of-mass convention.

.. tip::

   The recommended practice is to always use the GUI to create and edit
   ``.azr`` files. The internal format may change between versions.
