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

``abundanceSumSigma=s``
   The default sum-to-one width: every sample of the ``<abundances>`` block
   (below) without a ``sum`` line of its own adds
   :math:`((\sum_k w_k - 1)/s)^2` to the chi-squared, :math:`w_k` its
   abundances; absent or ``0`` is off.

Several compound nuclei.  The compound nucleus of a pair is inferred from its
columns (particle pair :math:`Z = z_1 + z_2`, :math:`A = \mathrm{round}(m_1 +
m_2)`; gamma pair the heavy particle; beta pair the daughter), and a project
may hold levels of several. An optional ``<abundances>`` block right after
``</levels>`` has one line per compound nucleus::

   <abundances>
   10 21 0.9048 0 0
   10 23 0.0925 1 0.01
   10 22 0.0027 0 0
   </abundances>

``Z A value vary sigma``: the isotopic abundance (of the target isotope that
forms the compound nucleus), whether it is a fit parameter (then named
``abundance_<Z>_<A>``, emitted right after the level parameters) and the sigma
of a Gaussian prior on it (0 = none). Without the block every abundance is 1,
fixed, and no parameter exists. A nucleus not formed by any pair draws a
warning.

**Samples.** Data sets taken on targets of different composition (a natural
one and an enriched one, or a compound such as Al\ :sub:`2`\ O\ :sub:`3`)
need several abundance sets in one project. A line prefixed with
``<name>:`` belongs to the named sample; the lines without a prefix form the
unnamed (default) sample; ``[<name>:] sum <s>`` gives a sample its own
sum-to-one width (``0`` off, absent: ``abundanceSumSigma``)::

   <abundances>
   10 21 0.9048 0 0
   10 23 0.0925 0 0
   10 22 0.0027 0 0
   enr22Ne: 10 23 0.99 1 0.005
   enr22Ne: 10 21 0.01 0 0
   enr22Ne: sum 0.01
   </abundances>

A nucleus without a line in a sample is 1, fixed. Each sample's abundances
are parameters, named ``abundance_<name>_<Z>_<A>`` for a named sample,
emitted sample by sample (the unnamed one first, in file order) right after
the level parameters. A sample name is one token without ``:`` or ``#``.
Stoichiometric weights need not sum to one: give that sample ``sum 0``.

A ``<segmentsData>`` or ``<segmentsTest>`` line may end with the keyword
blocks ``sqrtshift <b> <bError> <vary>``, ``natural <1|name>`` and
``shiftgroup <label>`` (in any order; writers put ``natural`` last).
``natural 1`` multiplies the segment's own term and each composite component
by the abundance of the term's entrance pair's compound nucleus in the
unnamed sample (or the only sample, if every sample is named);
``natural <name>`` uses the named sample. It is allowed for
angle-integrated, differential and capture observables (data ``isDiff`` 0,
1, 3, 4, 5, 6; test 0, 1, 4, 5), never for a ratio, only in a project with
an ``<abundances>`` block, and a name the block does not define is refused.

.. note::

   Compatibility: an engine older than this feature ignores the
   ``<abundances>`` block and the key, but reads ``natural 1`` as positional
   fields and misparses the line; a GUI older than it drops the block on
   save. A single-nucleus project is unaffected either way.

A ``<targetInt>`` line may end with keyword blocks (``beamprofile ...``,
``udr ...``, ``sigmalab``). ``sigmalab`` marks an energy-dependent beam-spread
expression as laboratory in, laboratory out (see
:doc:`../user_guide/experimental_effects`); an expression without it is read in
the pre-October-2026 centre-of-mass convention.

.. tip::

   The recommended practice is to always use the GUI to create and edit
   ``.azr`` files. The internal format may change between versions.
