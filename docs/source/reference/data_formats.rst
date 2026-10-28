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
- **Frame**: energies are always **laboratory** energies in **forward
  kinematics** (light particle as projectile, heavy particle as target).
  Angles and cross sections are laboratory values too, except for the
  **C.M.** differential observable (code 4 in data segments, 5 in test
  segments) and the analyzing power (code 7), which take centre-of-mass angles
  and cross sections.  Both frames work for particle and for capture (gamma)
  exit pairs.  A lab-frame capture point is converted with the photon
  aberration, :math:`\cos\theta_{cm} = (\cos\theta - \beta)/(1 - \beta\cos\theta)`,
  and the solid-angle Jacobian
  :math:`(1 - \beta\cos\theta)^2/(1 - \beta^2)`, where :math:`\beta` is the
  c.m. velocity.  Many low-energy capture measurements publish lab angles and
  cross sections without the transformation; with :math:`\beta \approx 0.01`
  (:sup:`12`\ C(:math:`\alpha,\gamma`) at a few MeV) the difference is a few
  percent and angle dependent, so enter such data as **Lab**.
- **Units**: energies in MeV, angles in degrees, cross sections in barns
  (angle-integrated) or barns/sr (differential).
- **Format**: columns can be separated by spaces or tabs.

Data Covariance Files
---------------------

A data segment whose line ends with ``covariance <file>`` (see
:doc:`../user_guide/segments` and :doc:`../user_guide/chi_squared`) takes the
covariance of its points from ``<file>``, a plain-text file with ``#``
comments and free-form numbers (any whitespace, line breaks anywhere)::

   # neutron groups of one experiment
   file data/n0.dat 1173        # rows 1..1173: the data rows of data/n0.dat
   file data/n1.dat 475         # rows 1174..1648
   block 3                      # a 3 x 3 symmetric block ...
     12 40 1180                 #   ... on these rows (1-based, any order)
     1e-9 2e-10 1e-10  2e-10 3e-9 0  1e-10 0 2e-9
   diagonal                     # N numbers, added to the diagonal
     ...
   vector                       # N numbers v: adds v v^T (repeatable)
     ...
   relvector                    # N fractions f: adds (f m)(f m)^T, m = the model
     ...
   matrix                       # N x N numbers: adds a full symmetric matrix
     ...

- The ``file`` lines come first. Each names a data file (path relative to the
  directory AZURE2 runs in, as on the segment line) and its number of data
  rows, counted as AZURE2 reads them: lines that are neither blank nor
  comments. The rows of the matrix are these data rows, file after file; a
  count that does not match the file is an error. ``size N`` may follow as
  a check.
- The terms (``block``, ``diagonal``, ``vector``, ``relvector``,
  ``matrix``) add up to :math:`C`, in the units of the files' value column
  (lab frame, as in the data file: (b/sr)\ :sup:`2` for differential cross
  sections). Blocks may overlap and need not be contiguous; ``block`` and
  ``matrix`` must be symmetric (to :math:`10^{-6}` relative).
- ``relvector`` multiplies its fractions by the model :math:`m` (in the
  data's units), frozen at the first evaluation of the run and written to
  ``covarianceModel.out``; see :doc:`../user_guide/chi_squared`.
- Every point of a segment naming the file must be a row of it (its data
  file listed), and no row may be fitted by two segments (overlapping
  windows on one file): both are errors. Rows not fitted drop out.
- AZURE2 never forms :math:`C`. It splits it into :math:`B` (``diagonal``,
  ``block`` and ``matrix`` terms, factored by connected components) and the
  rank-one ``vector``/``relvector`` terms, and whitens with
  :math:`W = (I + V T V^T) L_B^{-1}`, :math:`V = L_B^{-1} U`, at a cost of
  :math:`O(N k)` per evaluation for :math:`k` rank-one terms (about a
  millisecond for :math:`N = 14\,000`, :math:`k = 31`). :math:`B` must
  therefore be positive definite by itself: every fitted row needs a
  diagonal or block variance of its own. A ``matrix`` term makes one dense
  component of all rows.
- The data file's uncertainty column is replaced by :math:`\sqrt{C_{ii}}`
  (with the data value standing in for :math:`m` in ``relvector`` terms); the
  run log warns if it differed by more than 1 %. :math:`\chi^2` does not
  depend on it.

Phase Shift Data
^^^^^^^^^^^^^^^^

Phase shift data follows the same four-column format, but the third column
contains the phase shift value (in degrees) and the fourth column contains its
uncertainty. The convention is :math:`-90° < \theta_\text{lab} < 90°`.

Angular-Coefficient Data
^^^^^^^^^^^^^^^^^^^^^^^^

A segment of type **Angular Coefficients** (``isDiff`` 9) reads the same
four columns as ``E_lab  k  a_k  da_k``: the lab energy, the Legendre order
:math:`k` (an integer, 0-20, in the angle column), the coefficient
:math:`a_k = A_k/A_0` of :math:`d\sigma/d\Omega_\text{cm} = A_0\,[1 +
\sum_k a_k P_k(\cos\theta_\text{cm})]` and its uncertainty. Several orders
at one energy are separate rows. With the segment keyword ``absolute 1`` the
third column is the absolute :math:`A_k` in b/sr (:math:`A_0 = \sigma/4\pi`).
A normalised coefficient carries no normalization (the segment's norm must stay
fixed at 1).

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

``relativistic=1``
   Relativistic two-body kinematics: channel wave numbers, Sommerfeld
   parameters, lab <-> c.m. energies, angles and cross sections from the
   invariant mass (:doc:`../user_guide/particle_pairs`, "Relativistic
   Kinematics"). Absent or ``0``: non-relativistic, as before.

``relativisticReferencePair=<key>``
   The pair whose :math:`E(s)` ``relativistic.out`` uses; absent or ``0``: the
   first particle pair of each compound nucleus.

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

Total capture is the exit key ``-1`` together with the total-capture data
type (3 on a data line, 4 on a test line); the reader expands such a segment
into one term per capture pair. The key without the type is refused (before
October 2026 the engine read past the start of its pair list for such a line
and could crash); the type with a real exit key has always meant total
capture, the exit key being ignored, and now draws a warning saying so.

A ``<segmentsData>`` line may end with the keyword blocks
``sqrtshift <b> <bError> <vary>``, ``sample <1|name>``, ``shiftgroup <label>``,
``shiftfngroup <label>``, ``normgroup <label>``, ``covariance <file>``,
``shiftfn "<f(x)>" N ...``, ``total <0|1>`` and ``spin <s>`` (in any order;
writers put ``sample`` last). A ``<segmentsTest>`` line honours ``sample``,
``total`` and ``spin``; the shift, norm and shift-function blocks have no
meaning without data and are skipped there (before October 2026 they stopped
the keyword reader, so a ``sample`` or ``total`` written after them was lost),
and ``covariance`` on a test line is refused (:doc:`../user_guide/segments`).
The three shift fields of a data line (``shift shiftError vary``) are
optional: a line may carry one shift value or none before the data file, and
the reader tells a value from a file name by the whole token, so a file whose
name starts with digits (``12C_p.dat``) is a file (before October 2026 it was
read as a 12 MeV shift and the file ``C_p.dat``).

The ``<parameterSettings>`` block (written by the Fitting tab) has one line per
parameter: a label of several words, then ``value lower upper error fit_error
nuisance category minuit_index``. The numeric fields accept ``nan`` and
``inf`` (any case, signed): a fit error written from a singular covariance can
be either, and before October 2026 the parser threw on them and the
command-line run ended before any output. The fit error is a report, not an
input; a limit that is not finite is ignored (no limit).
``natural`` is the keyword's older name and is still read.
``spin <s>`` on a phase-shift line (``isDiff`` 2) names the channel spin of
the partial wave: a phase shift is one (J, l, s) wave, and the engine refuses a
phase-shift segment whose (J, l) has several channel spins in the level scheme
unless the line says which (before October 2026 it summed the S-matrix
elements of all of them and reported half the argument, which is not a phase
shift of anything; an older engine ignores the keyword). A single-spin (J, l),
such as any spin-0 target, needs no keyword.
``total 1`` makes every term whose exit pair is its entrance pair the total
cross section of that (neutral) pair, by the optical theorem
(:doc:`../user_guide/segments`, "Neutron Total Cross Sections"); an older
engine ignores it and computes the elastic cross section.
``sample 1`` multiplies the segment's own term and each composite component
by the abundance of the term's entrance pair's compound nucleus in the
unnamed sample (or the only sample, if every sample is named);
``sample <name>`` uses the named sample. It is allowed for
angle-integrated, differential and capture observables (data ``isDiff`` 0,
1, 3, 4, 5, 6; test 0, 1, 4, 5), never for a ratio, only in a project with
an ``<abundances>`` block, and a name the block does not define is refused.

.. note::

   Compatibility: an engine older than this feature ignores the
   ``<abundances>`` block and the key, but reads ``sample 1`` as positional
   fields and misparses the line; a GUI older than it drops the block on
   save. A single-nucleus project is unaffected either way.

**Measured scattering lengths.** An optional ``<scatteringLengths>`` block
(anywhere in the file; the GUI writes it after ``</segmentsTest>``) lists
measured s-wave scattering lengths, in fm, that the fit must reproduce
(:doc:`../theory/scattering_length`)::

   <scatteringLengths>
   # active pairKey s value sigma [imValue imSigma]
   1 1 1 17.34 1.2
   1 1 2 -3.18 0.5
   1 3 -1 3.59 0.02
   1 2 1 9.0 0.5 -0.3 0.05
   </scatteringLengths>

``active pairKey s value sigma``: whether the line is used, the particle pair
(its key), the channel spin :math:`s` (:math:`J = s`, parity
:math:`\pi_1\pi_2` for :math:`l = 0`), and the real part of the scattering
length with its uncertainty.  ``s = -1`` selects the coherent combination
:math:`\sum_s (2s+1)\,a_s / ((2j_1+1)(2j_2+1))` of a neutron-type
measurement (not for an identical pair).  Two optional numbers constrain the
imaginary part, for a pair whose threshold lies above an open channel.  Each
active line adds :math:`((\mathrm{Re}\,a_0 - \mathrm{value})/\mathrm{sigma})^2`
(and the same for the imaginary part when its sigma is nonzero) to the
chi-squared; ``sigma = 0`` leaves that part unconstrained.  ``#`` starts a
comment.  An unknown pair key, a photon pair, a channel spin the pair does not
have, or a malformed line stops the run with a message naming the line.
Neutron tables list *bound* lengths :math:`b`; enter the free length
:math:`b\,A/(A+1)`.

.. note::

   Compatibility: an engine older than this feature ignores the block (the
   constraints then do nothing); a GUI older than it drops the block on save.

A ``<targetInt>`` line may end with keyword blocks (``beamprofile ...``,
``udr ...``, ``sigmalab``, ``screening Ue vary sigma``). ``screening`` is
electron screening of the listed segments: the screening potential in eV, 1 to
make it a fit parameter (``screening_<n>``), and the width in eV of a Gaussian
prior about the value given (0 = none); see
:doc:`../user_guide/experimental_effects`. An engine older than October 2026
refuses a line with this keyword, and an older GUI drops it on save. ``sigmalab`` marks an energy-dependent beam-spread
expression as laboratory in, laboratory out (see
:doc:`../user_guide/experimental_effects`); an expression without it is read in
the pre-October-2026 centre-of-mass convention.

.. tip::

   The recommended practice is to always use the GUI to create and edit
   ``.azr`` files. The internal format may change between versions.
