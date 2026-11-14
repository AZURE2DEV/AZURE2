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
  A point at a lab angle that is folded over energy (a target, a beam spread,
  a detector resolution) is folded in the laboratory: at every energy
  :math:`E` of its fold the c.m. angle :math:`\theta_{cm}(E;\theta_{lab})` and
  the Jacobian :math:`J(E;\theta_{lab}) = d\Omega_{lab}/d\Omega_{cm}` are those
  of its lab angle at :math:`E`, so the model is the kernel's average of the
  lab cross section :math:`\sigma_{cm}(E,\theta_{cm}(E))/J(E)`.  Until
  2026-11-12 every energy of the fold took the point's own c.m. angle and
  Jacobian; the difference grows with the window's width and with the slope
  of the angle and the cross section across it (measured: up to
  :math:`2.6\times10^{-5}` relative for a 2 keV beam spread on
  :sup:`15`\ N(p,\ :math:`\alpha`) at 150° lab; :math:`2.2\times10^{-2}`
  for a 300 keV window 0.4 MeV above the :sup:`7`\ Li(p,n) threshold at 30°
  lab, ``tests/lab_frame_fold``).  An analysing power is given at c.m. angles
  only.
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
  project file's directory, as on the segment line) and its number of data
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
uncertainty. The segment line names the elastic partial wave :math:`(J, l)`,
and its channel spin :math:`s` where needed (``spin <s>``, below). The quantity
is the nuclear phase shift of that wave,

.. math::

   \delta = \tfrac12 \arg\big(U_{cc}\,e^{-2i\omega_c}\big),

measured from the Coulomb phase :math:`\sigma_l` and including the
hard-sphere part :math:`-\phi_c` (:math:`U`, :math:`\omega_c` and
:math:`\phi_c` as in :doc:`../theory/conventions`). The engine reports it in
:math:`-90° < \delta \le 90°`, and in :math:`0° \le \delta < 180°` when the
two particles of the pair are identical (``GenMatrixFunc.cpp:553-564``). The
residual is the plain difference of data and model, not reduced modulo
180°, so the data must be given in the same range.

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

In format 3 (``quantity legendre``, ``columns E_lab[MeV] order value error``)
the field ``expansion`` names the normalisation, for data sets and
extrapolations alike (EXFOR's codes in brackets), with
:math:`d\sigma/d\Omega_\text{cm} = \sum_L A_L P_L(\cos\theta_\text{cm})`:

==============  ==============================  ===========  ==============
``expansion``   coefficient of order :math:`L`  unit         EXFOR
==============  ==============================  ===========  ==============
``a_over_a0``   :math:`A_L/A_0`                 ``1``        ``LEG/RS``
``absolute``    :math:`A_L`                     b/sr         ``LEG``
``rsl``         :math:`A_L/((2L+1)A_0)`         ``1``        ``LEG/RSL``
``2l2``         :math:`2A_L/(2L+1)`             b/sr         ``LEG/2L2``
``4pi``         :math:`4\pi A_L`                b            ``LEG/4PI``
==============  ==============================  ===========  ==============

The normalised forms (``a_over_a0``, ``rsl``) keep their norm fixed; the
engine computes :math:`A_L/A_0` or :math:`A_L` and multiplies by the factor of
:math:`L` the table implies, so a fold, an energy shift and the gradient act
on all five alike.  An extrapolation writes :math:`L = 0 \ldots K` of its
``orders 0 K`` in its expansion (the ``.acoeff`` file and pyazr's
``evaluate(x).data[i].extra["legendre"]``).  ``tests/legendre_expansions`` checks every
order against a Gauss--Legendre projection of the engine's own
:math:`d\sigma/d\Omega`.

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
- Fitting parameter settings (``<parameterSettings>``: limits, errors and
  nuisance flags, and optional ``segment_N_norm prior_centre c`` rows that
  hold the centre of a normalization or energy-shift prior apart from the
  segment's own value; see :ref:`nominal-norm`)
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

``priorCovariance=<file>``
   A multivariate Gaussian prior :math:`(\mathbf{p}-\boldsymbol{\mu})^T C^{-1}
   (\mathbf{p}-\boldsymbol{\mu})` on the parameters the file lists
   (:doc:`../user_guide/chi_squared`, "Prior parameter covariance"). The path is
   relative to the ``.azr``'s directory (or absolute); everything after the first
   ``#`` on the line is a comment. Absent or empty: none. The file is plain text
   with ``#`` comments: an optional ``#parametrization <0|1|2>`` line (as the first
   line of ``param.sav``), the parameter names (``param.sav``'s), their means, and
   the :math:`n` rows of the symmetric, positive-definite :math:`n\times n`
   covariance::

       #parametrization 1
       width_2_1 segment_4_norm
       0.8206 1.05
       5.0e-5 1.0e-6
       1.0e-6 4.0e-4

   Refused with an ``ERROR`` (nothing is written): an unknown name, a duplicate
   name, a non-number, a count mismatch, a non-symmetric (1e-10 relative) or
   non-positive-definite matrix, a tag other than the run's parametrization, a
   file that cannot be read. Parameters fixed in the run are marginalized. Only a
   run with data uses it (an extrapolation or a rate ignores the key). Engines
   older than 2026-11-05 ignore the key silently.

``use_brune=0|1`` and the ten other options of the calculation
   ``use_park``, ``transform``, ``ignore_externals``, ``use_long_wavelength``,
   ``ec_convection_current``, ``use_rmc``, ``use_gsl_coul``,
   ``use_coulomb_table``, ``use_ec_table``, ``replace_grids``, under the names
   of pyazr's keywords: a
   flag on the command line or a keyword of a session wins over the key, the
   key over the default; absent: the default, as before
   (:ref:`options-as-config-keys`, with the flags, the defaults and the pairs
   that are refused).

The ``<levels>`` block has one line per channel of a level, 31
whitespace-separated columns (read by ``azure2::io::LevelRow``, ``io/include/azure2/io/LevelRows.h``), and
optionally one, two or three more:

.. list-table::
   :widths: 12 18 70
   :header-rows: 1

   * - Column
     - Name
     - Meaning
   * - 1-4
     - ``J pi E fix``
     - the level's spin, parity, excitation energy (MeV) and fixed flag
   * - 5-6
     - ``aa ir``
     - entrance key (written as ``1``, not read) and the channel's pair key
   * - 7-8
     - ``2s 2l``
     - twice the channel spin and twice the orbital angular momentum
   * - 9-11
     - ``id active chfix``
     - level ID, level active flag, width fixed flag
   * - 12
     - ``gamma``
     - partial width (eV) of an open channel, ANC (fm\ :sup:`-1/2`) of a closed one
   * - 13-16
     - ``j1 pi1 j2 pi2``
     - spin and parity of the light (1) and heavy (2) particle **in the state
       the pair has them in**
   * - 17
     - ``e2``
     - excitation energy of the heavy particle (MeV); of the final level for a
       gamma pair
   * - 18-21
     - ``m1 m2 z1 z2``
     - **ground-state** atomic masses (u) and charges
   * - 22-23
     - ``Sin S``
     - the entrance pair's separation energy (written, not read) and this pair's
       ground-state separation energy (MeV)
   * - 24-26
     - ``j3 pi3 e3``
     - unused (written as ``0 0 0.0``)
   * - 27-28
     - ``type R``
     - pair type (0 particle, 10 gamma, 20 beta) and channel radius (fm)
   * - 29-30
     - ``g1 g2``
     - g-factors of the two particles
   * - 31
     - ``ecmask``
     - external-capture multipolarity mask
   * - 32-34
     - ``B rwa e1``
     - **optional**, told apart by their number: one trailing field is
       ``e1``; two are ``B rwa``; three are ``B rwa e1``. ``B`` is the THM
       binding energy of the transferred particle (MeV, per pair), ``rwa``
       the width flag of the channel (1: column 12 is a reduced-width
       amplitude in MeV\ :sup:`1/2`; 0: a physical width or ANC), and ``e1``
       the excitation energy of the light particle (MeV), since November
       2026. Absent means 0. The pair opens at :math:`S + e_1 + e_2`.

A pair's columns (13-21, 23, 27-31, ``B`` and ``e1``) are read from the first
line that names it. ``e1`` is written **only when it is nonzero**. pyazr keeps
the layout of each line (31 or 33 columns) and appends ``e1`` to it; the GUI
always writes ``B`` and ``rwa`` (33 columns, both 0 for a project without THM
content) and appends ``e1`` as column 34. A project without excited light
particles therefore keeps its 31 columns through pyazr, and has 33 after a GUI
save. An engine or GUI from before November 2026 does not know ``e1``: the
engine ignores it silently (the pair then opens at :math:`S + e_2`), and the
GUI drops it when it saves the project. The engine refuses a nonzero ``e1`` on a gamma or beta pair, a
negative one, and two different values on lines of one pair. Two particles of
a particle pair are identical when Z, mass, spin, parity, g-factor and
excitation (``e1 = e2``) are all equal (:doc:`../user_guide/particle_pairs`).

A pair that no level names yet is held by a line of its own, a *pair row*:
inactive (``active`` 0) with level ID 0, J 0+, E 0, l 0, s 0 and width 0, the
pair's columns as on any line. The engine reads it as any inactive line (its
pair belongs to the system only with ``includeAllJPi``); pyazr's
``Project.levels()`` lists no level for it, ``Project.add_pair`` writes it and
``Project.add_level`` removes it once a level names the pair.

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
``shiftfn "<f(x)>" N ...``, ``total <0|1>``, ``transmission <n> [<vary> <sigma>]``,
``yield <n> [<vary> <sigma>]``, ``thicknessgroup <label>`` and ``spin <s>`` (in any order;
writers put ``sample`` last). A ``<segmentsTest>`` line honours ``sample``,
``total``, ``transmission``, ``yield`` and ``spin``; the shift, norm, thickness-group and shift-function blocks
have no meaning without data and are skipped there, as is the fit form of a ``transmission`` / ``yield`` block
(a test segment uses its :math:`n`) (before October 2026 they stopped
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
``transmission <n>`` makes the segment the transmission
:math:`e^{-n\sigma}` of a sample :math:`n` atoms/barn thick, and ``yield <n>``
the self-shielded yield :math:`(\sigma/\sigma_\mathrm{tot})(1 -
e^{-n\sigma_\mathrm{tot}})`; the data file then holds transmissions or yields,
and under a sub-point target effect the transform comes before the kernel
(:doc:`../user_guide/segments`, "Sample Transmission and Self-Shielded
Yield"). An engine older than November 2026 ignores both keywords silently.
The fit form ``transmission <n> <vary> <sigma>`` (``vary`` 0 or 1, ``sigma`` the
Gaussian prior width in atoms/barn, 0 = none) makes :math:`n` the fit parameter
``segment_<key>_thickness``, and ``thicknessgroup <label>`` makes the segments
with the label share the thickness of the first of them (each member needs its
own ``transmission`` or ``yield`` block); a malformed fit form is refused. The
engines of early November 2026 without the fit form read :math:`n` and ignore
the rest of the block and ``thicknessgroup`` silently.
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

**Model data.** An optional ``<modelData>`` section holds literature values of
level properties -- a resonance strength, a resonance energy, a total or
partial width, a branching ratio, an ANC -- and scattering lengths, each
compared with a function of the level parameters rather than with the model
at a beam energy (``docs/upgrade/project_format3.md`` section 3.10).  It is
the first section written in the records of project format 3
(:doc:`format3`): a header ``datum <name>`` and one ``<key> <value>`` line per
field, ``#`` comments anywhere.  T2's level 5 (``tests/upgrade/T2``), with
the values ``tests/model_data`` checks::

   <modelData>
   datum wg_L5
       quantity   resonance_strength
       level      5                  # the level ID of <levels> (column 9)
       J          3-                 # required: the level's J^pi
       reaction   15N(p,a)12C        # entrance and exit pairs
       value      2.63 keV
       error      0.15 keV
       source     "author, journal volume page (year), table 2"
   datum Er_L5
       quantity   resonance_energy
       level      5
       J          3-
       frame      E_lab              # E_x, or E_cm / E_lab of a pair
       pair       p+15N
       value      1.2185 MeV
       error      0.3 keV
   </modelData>

A level is named by its ID (column 9 of ``<levels>``) and its :math:`J^\pi`,
which must agree; a pair by its particles and states (``p+15N``,
``a+12C[4.439]``: an excitation energy matches to the digits written,
``g+16O[6.13]`` a capture pair).  The fields of each quantity, its units and
its EXFOR codes are the table of :doc:`format3`.  Every datum adds
:math:`((q - d)/e)^2` to the chi-squared, with :math:`q` evaluated at every
step of a fit from the **observed** energy and widths of the level (those
``parameters.out`` reports, in the Brune, formal and Park bases alike):

* ``resonance_strength``: :math:`\omega\gamma = \omega\Gamma_a\Gamma_b/\Gamma`,
  :math:`\omega = (2J+1)(1+\delta_{aA})/((2j_a+1)(2j_A+1))`, with
  :math:`\Gamma_a`, :math:`\Gamma_b` the observed widths summed over the
  channels of the entrance and exit pairs and :math:`\Gamma` over all
  channels (eV, meV, ueV, keV).  :math:`\Gamma` includes the radiative widths, which the level matrix leaves out,
  so a calculated yield of the same level has the strength :math:`\omega\Gamma_a\Gamma_b/\sum\Gamma_\text{particle}`,
  larger by :math:`1+\Gamma_\gamma/\sum\Gamma_\text{particle}` (:doc:`../user_guide/gamma_decay`);
* ``resonance_energy``: :math:`E_x` of the level, or :math:`E_\mathrm{c.m.} = E_x - S - E^*`
  and :math:`E_\mathrm{lab}` of the pair (``relativistic=1`` honoured);
* ``total_width``, ``partial_width`` (``pair``), ``branching_ratio``
  (:math:`\Gamma_b/\Gamma`, ``reaction``);
* ``anc``: :math:`|C|` of a closed channel (``pair``, ``l``, ``s``; 0 if a fit
  moves the level above the threshold);
* ``scattering_length``: as a ``<scatteringLengths>`` line (``pair``, ``s``
  or ``coherent``, ``value``/``error``, optional ``imag``/``imag_error``).

The engine refuses, with the code in the message: a level ID no active level
has (``E_UNKNOWN_LEVEL``), a :math:`J^\pi` that is not the level's
(``E_LEVEL_JPI``), a reaction or pair no pair has (``E_UNKNOWN_PAIR``) or two
pairs have (``E_AMBIGUOUS_PAIR``), a pair that is no channel of the level
(``E_NOT_A_CHANNEL``), a strength of a level below its entrance threshold
(``E_OUT_OF_RANGE``: give an ``anc``), an ANC of an open channel, two data of
one name, and every error of the format-3 table.  It warns when eq. (8) and
the area of the resonance differ (``docs/upgrade/data_kinds.md`` section 4.7):
:math:`\Gamma/E_r > 1\,\%` (``W_BROAD_LEVEL``) or another level of the same
:math:`J^\pi` within :math:`10\,\Gamma` (``W_NEARBY_LEVEL``), from the widths
written in ``<levels>``; a non-resonant amplitude in the same channel, which
also changes the area, is not detected.  The table of the model values is
written to ``modelData.out`` and the run log, and the term's total to
``chiSquared.out`` (:doc:`output_files`); pyazr gives the same values as
records, ``pyazr.model_data(model)``.  The level quantities' gradient is a
central difference over the parameters of the level's :math:`J^\pi` group;
the scattering lengths' is exact.

A ``<targetInt>`` line may end with keyword blocks (``beamprofile ...``,
``udr ...``, ``sigmalab``, ``screening Ue vary sigma``). ``screening`` is
electron screening of the listed segments: the screening potential in eV, 1 to
make it a fit parameter (``screening_<n>``), and the width in eV of a Gaussian
prior about the value given (0 = none); see
:doc:`../user_guide/experimental_effects`. An engine older than October 2026
refuses a line with this keyword, and an older GUI drops it on save. ``sigmalab`` marks an energy-dependent beam-spread
expression as laboratory in, laboratory out (see
:doc:`../user_guide/experimental_effects`); an expression without it is read in
the pre-October-2026 centre-of-mass convention. ``densityfit vary sigma [max]``
and ``densitygroup <label>`` make a thin target's density a fit parameter
(``density_<n>``) with an optional prior (atoms/cm², or ``2.3%`` of the line's
density), the largest density its sub-point grids span, and a label shared by
lines that see one target; see :doc:`../user_guide/experimental_effects`. Engines
older than the feature refuse a line with these keywords; GUIs older than it drop
them on save. ``doppler T`` is free-gas Doppler broadening of the listed
segments at the effective temperature T (kelvin), applied to the cross section
before any ``transmission``/``yield`` transform and before the line's resolution
kernel; a line may carry it alone. Engines older than November 2026 refuse a line
with it; older GUIs drop it on save. ``doppler T vary sigma [Tmax]`` makes T the fit
parameter ``doppler_<n>`` (vary 0 or 1, a prior of width sigma kelvin about T, 0 =
none, and the largest temperature the sub-point grids span, by default
max(T + 3 sigma, 1.21 T), or 2 T without a prior); see
:doc:`../user_guide/experimental_effects`. Engines older than November 2026 refuse
it; a build that knows only ``doppler T`` ignores the numbers after T (a fixed
temperature, no parameter) and its GUI misreads the line.

.. tip::

   The recommended practice is to create and edit ``.azr`` files with
   ``pyazr.Project`` (:doc:`../user_guide/pyazr`), whose load-and-save round trip
   (every byte kept) ``tests/pyazr/project_test`` checks. The internal format may change
   between versions.

Project Files in Format 3
-------------------------

A project in format 3 (``format=3`` in ``<config>``; the design is
``docs/upgrade/project_format3.md``, the table of every record kind, field,
quantity, unit and diagnostic is :doc:`format3`) keeps ``<config>``,
``<levels>``, ``<parameterSettings>``, ``<mcmc>`` and ``<thm>`` and replaces
``<segmentsData>``, ``<segmentsTest>``, ``<targetInt>``, ``<abundances>`` and
``<scatteringLengths>`` by named records: ``dataset``, ``normalization`` and
``energy_scale`` in ``<data>``, ``extrapolation`` in ``<extrapolations>``,
``target``, ``beam`` and ``detector`` in ``<targets>``, ``datum`` in
``<modelData>``.  A record is a header line ``<kind> <name>`` and one
``<key> <value>`` line per field; ``#`` starts a comment.  Records refer to one
another by name, in any order.  An excerpt of T2 (``tests/upgrade/T2``, written
by ``AZURE2 --convert T2.azr T2_f3.azr``; the whole file is
``experiment/tests/data/T2_format3.azr``)::

   <config>
   format=3
   </config>
   <targets>
   target tgt_s01_int_thin                         # <targetInt> line 1
       kind           layer
       active         15N
       areal_density  5e17 atoms/cm2 fixed
       stopping       formula "a0*x^a1" a0=1.6e-20 a1=-0.45 MeV*cm2/atom
       straggling     bohr
   beam beam_s01_int_thin
       spread         0.001 MeV lab
   </targets>
   <data>
   normalization D                                 # shared by s02 and s03
       value          1 fit prior 1 8%
   energy_scale G
       offset         0 MeV fit prior 0 0.003
       function       "a0*(x-0.8)" MeV
       param          a0 0 fit prior 0 0.01
   dataset s01_int_thin
       reaction       15N(p,a)12C
       quantity       sigma
       unit           b
       file           data/s01_int_thin.dat
       columns        E_lab[MeV] _ value error
       energy_window  0.199999 1.400001 MeV
       norm           1 fit prior 1 5%
       shift          offset 0 MeV fit prior 0 0.003
       target         tgt_s01_int_thin
       beam           beam_s01_int_thin
       quadrature     points 20 width_multiplier 10 per_width 20
   dataset s02_diff090
       reaction       15N(p,a)12C
       quantity       dsigma_domega
       unit           b/sr
       file           data/s02_diff090.dat
       columns        E_lab[MeV] theta_lab[deg] value error
       energy_window  0.249999 1.350001 MeV
       angle_window   90 90 deg
       norm           D
       shift          G
   </data>

The engine reads such a file in one pass: the records are read and checked
against the table, the reactions are resolved against the pairs of
``<levels>`` (a label names a pair by its particles and states; where it names
two, ``entrance_pair`` / ``exit_pair`` give the key), the model data against
its levels, and the data files are opened **relative to the working
directory** (as for an old project).  Any error refuses the project with the
same diagnostics ``AZURE2 --validate`` prints.  ``AZURE2 --describe FILE``
(``pyazr.describe``) prints what was understood: per data set the observable,
the pair keys of each reaction term, the data file's resolved path, the
effects linked to it and its fit parameters.

Where the engine reading a format-3 file behaves differently from the old
engine on the same project (``AZURE2 --convert`` reports the rest):

* A data file that cannot be read refuses a run with data (``E_DATA_FILE``);
  the old engine skipped that data set with a warning and fitted the others.
  A run without data (extrapolations, menu 3, and pyazr calls that only
  extrapolate) does not read the data files: one that cannot be read is the
  warning ``W_DATA_FILE_NOT_READ``.  ``--validate`` and ``--describe`` require
  them.
* A member of a varied shared ``normalization`` has a fixed parameter of its own
  (``D.value[<data set>]``) whether or not its old line varied it, and the
  start and step of a follower's fixed parameters are the record's.
* Scattering lengths are model data: ``chiSquared.out`` labels them
  ``Model-Data`` and ``modelData.out`` lists them.
* An abundance-weighted output names its sample by its target.

Ratio to Rutherford and to Mott.  ``quantity ratio_to_rutherford`` (an elastic
reaction of two different charged particles) and ``ratio_to_mott`` (identical
particles) are :math:`d\sigma/d\Omega` over the Rutherford or Mott cross
section at the same c.m. energy and angle (``unit 1``; lab or c.m. angles,
the angle converted at load, the value not, since it is the same in both
frames).  The engine computes the reference with its own :math:`\eta` and
:math:`k` (the pair's masses, the run's kinematics; Mott with the
interference weight :math:`(-1)^{2s}/(2s+1)` of the particles' spin) and folds
it with the data set's beam spread, target layer or beam profile on the same
sub-points as the cross section: the model is the ratio of the two folds.  A
fitted energy shift moves both.  Refused: ``ratio_to_rutherford`` for
identical particles and ``ratio_to_mott`` for different ones, a neutral pair,
an inelastic group, and a Doppler-broadened target, electron screening and
detector attenuation coefficients on such a data set (``E_NOT_BUILT``).  The worked example is
:doc:`../examples/ratio_to_rutherford`.

An energy shift (``shift``) is a calibration of the beam energy: a folded point
measured at :math:`E` with the shift :math:`s` is folded as the same point
without a shift at :math:`E + s`, the target's stopping taken at the shifted
energies (:doc:`../user_guide/experimental_effects`).

Parameter names follow the records: a data set's own are
``<data set>.norm``, ``<data set>.shift.offset`` and ``<data set>.shift.sqrt``
(every data set has the two shift parameters, fixed at 0 without a shift); a
shared record's belong to the first active data set of the group
(``D.value``, ``G.offset``, ``G.a0``, ``<target>.areal_density``) and a later
member keeps a fixed parameter of its own, ``D.value[s03_diff150]``; the
effects of a target are ``<target>.areal_density``, ``<target>.doppler``,
``<target>.screening``, an isotope's abundance ``<target>.isotope.<nuclide>``.
``AZURE2 --convert`` reports the map from the old names
(``segment_2_norm``) to these in its ``--json-report`` (``parameters``) and
renames the rows of ``<parameterSettings>`` by it.

Units.  Every number may be written in any unit of its dimension in the unit
table (:doc:`format3`); the engine converts it to its own unit (MeV, deg, b,
b/sr, MeV b, atoms/cm\ :sup:`2`, MeV cm\ :sup:`2`/atom, eV for U\ :sub:`e`,
keV\ :sup:`1/2` for a straggling coefficient) with the table's factor, in one
function: a factor that is a power of ten divides by the exact integer, so a
value written as ``12.5 mb`` is the same double as ``0.0125 b`` (a value exact in
the written unit gives the correctly rounded engine value).  The fields a record
holds (windows, shifts, grids, a target's ``areal_density`` -- atoms/b or
atoms/cm\ :sup:`2`, for a layer and for a neutron sample --, U\ :sub:`e`, beam
spread and profile, detector resolution, UDR times and flight path, a
``stopping`` or ``spread`` formula, whose value is then multiplied by the unit's
factor) are converted when the project is read.  The columns of a data file are
converted row by row where the points are made: the energy column
(``E_lab[keV]``, ``E_cm[MeV]`` -- a c.m. energy is taken to the lab energy with
the entrance pair's ``<levels>`` masses, as a ``cm`` grid --, and the
``window_lo``, ``window_hi`` columns with it), the angle column (``[rad]``), and
the value and error columns (the data set's ``unit``).  The ``energy_window`` is in
the frame and unit of the energy column.  Outputs (``AZUREOut*``,
``normalizations.out``, ...) list the data in the engine's units.  An
extrapolation in another unit than the engine's is computed in that unit (the
model times the inverse factor: ``unit mb`` lists millibarn).  Refused
(``E_NOT_BUILT``): a value unit other than the engine's on a data set with a
``covariance`` file (the file's absolute entries are in the engine unit and are
shared by data sets), an ``E_cm`` column with ``window_lo``/``window_hi`` (the
windows of a beam-profile kernel are lab energies).  Mass-thickness units
(µg/cm\ :sup:`2`) are not in the table, so they are refused as unknown units.

S factor.  ``quantity sfactor`` (units ``MeV*b``, ``keV*b``, ``eV*b``; an ``E_lab``
or ``E_cm`` column; a charged entrance pair) compares the data with
:math:`S(E) = E\,\sigma(E)\,e^{2\pi\eta(E)}` of the model, with :math:`E` the
c.m. energy of each point as the model evaluates it -- after a fitted energy
shift, which therefore moves the factor -- and :math:`\eta = Z_1 Z_2 \alpha
\sqrt{\mu c^2/2E}` (non-relativistic, whatever ``relativistic`` says) with the
pair's ``<levels>`` masses.  The data are not converted to a cross section.  With
a target or beam fold the model is :math:`E e^{2\pi\eta(E)}\langle\sigma\rangle` at
the point's energy.  An extrapolation of ``quantity sfactor`` lists S, and an
adaptive grid of it refines on S: on the Coulomb slope below a resonance, where
the cross section falls by orders of magnitude and S is smooth, the relative test
on S needs a fraction of the points the test on the cross section takes (10
against 125 for 12C(p,γ) at 0.05 -- 0.25 MeV, ``tests/units_scales``); with a
narrow peak in the range both are governed by its floor (10\ :sup:`-3` of the
largest value).  Refused: an uncharged entrance
(``E_REACTION_NOT_ALLOWED``) and a THM data set (``E_NOT_BUILT``).  In the
``AZUREOut`` files of such a data set the cross-section columns hold
:math:`\sigma = S/(E e^{2\pi\eta})` and the S-factor columns S.

Modified S factor.  ``quantity modified_sfactor`` with the field ``shape <g>
1/MeV`` (or ``1/keV``) compares the data with the modified S factor of
heavy-ion fusion, :math:`S^*(E) = E\,\sigma(E)\,e^{2\pi\eta(E) + gE}`
(Patterson, Winkler and Zaidins, ApJ 157 (1969) 367), in the units of
``quantity sfactor``, with the same :math:`E` and :math:`\eta`; everything said
above of the S factor holds for it (energy shifts, folds, adaptive grids).  The
constant :math:`g` is given, never computed: 0.46 /MeV for
:sup:`12`\ C + :sup:`12`\ C, in general
:math:`0.122\sqrt{\mu R^3/(Z_1 Z_2)}` /MeV with :math:`\mu` in u and the
radius :math:`R` in fm.  With ``shape 0 1/MeV`` the model is the S factor bit
for bit.  Refused: no ``shape`` (``E_FIELD_REQUIRED``), a ``shape`` that is no
inverse energy (``E_UNIT_NOT_ALLOWED``), an uncharged entrance
(``E_REACTION_NOT_ALLOWED``), a THM data set (``E_NOT_BUILT``).  In the
``AZUREOut`` files the cross-section columns hold :math:`\sigma` and the
S-factor columns :math:`S^*` (``tests/modified_sfactor``).

Ratio to a reference cross section.  ``quantity ratio_to_reference`` (unit
``1``; an ``E_lab`` or ``E_cm`` column) with ``reference <file> [unit <u>]``
compares the data with :math:`\langle\sigma\rangle/\langle\sigma_\text{ref}\rangle`:
the file has two columns, the energy in the data set's energy role, frame and
unit (an extrapolation: its grid's) and the reference cross section (``b``
unless ``unit`` names another cross-section unit); between rows it is
interpolated linearly in :math:`\ln\sigma` against :math:`E`, and never
extrapolated.  With a target or beam fold the reference is folded on the
point's own sub-points, as the Rutherford reference is.  Refused: no
``reference`` (``E_FIELD_REQUIRED``), a file that is not a table (fewer than two
rows, energies not increasing), a value :math:`\le 0`, a data row or grid point
outside the table (``E_REFERENCE_TABLE``).  A fold's sub-point or a shifted
energy outside the table makes that point's model NaN: give a table that covers
the fold.  ``tests/ratio_to_reference``.

Yields.  ``quantity yield`` on a ``target`` of kind ``infinite`` compares the
data with the thick-target yield the engine folds,
:math:`Y(E_0) = 10^{-24}\int_{1\,\mathrm{keV}}^{E_0}\langle\sigma\rangle/\varepsilon\,dE`
(reactions per incident particle; ``1/particle``, or ``1/particle/sr`` with a lab
angle).  On a ``target`` of kind ``layer`` the integral stops at the back of the
layer, :math:`Y(E_0) = 10^{-24}\int_{E_{out}}^{E_0}\langle\sigma\rangle/\varepsilon\,dE`
with :math:`\int_{E_{out}}^{E_0}dE/\varepsilon = N`: :math:`E_0` and :math:`E_{out}`
lab energies (MeV), :math:`\sigma` in b (c.m., at :math:`E_{cm} = E\,m_A/(m_a+m_A)`),
:math:`\varepsilon` the target's stopping formula (lab MeV cm\ :sup:`2` per active
nucleus, the nucleus :math:`N` counts) and :math:`N` its ``areal_density``
(atoms/cm\ :sup:`2`; fitted when the field says ``fit``).  The engine finds
:math:`E_{out}` with the quadrature of the yield itself, so a constant cross section
gives :math:`10^{-24}\sigma N`; a narrow resonance gives the arctan step whose
plateau is :math:`10^{-24}(\lambda_r^2/2)\,\omega\gamma\,(m_a+m_A)/(m_A\varepsilon)`
(``tests/layer_yield``, README).  The same target with ``quantity sigma`` is the
target-averaged cross section, the yield divided by :math:`10^{-24}N`.  An
``isotope`` line multiplies the yield by its abundance (then :math:`\varepsilon` and
:math:`N` count atoms of the element).  Per unit of beam charge (``1/uC``, ``1/mC``, ``1/C``, with
``charge_state q``) the model is multiplied by the beam particles per unit charge,
:math:`10^{-6}\,\mathrm{C}/(q e)` for 1/µC; counts per unit charge
(``counts/uC``, ...) need ``efficiency``, a fixed number that multiplies the model
as well (a fitted efficiency is refused: it cannot be told apart from the
normalization; give its uncertainty as the norm's prior); ``arbitrary`` needs a
free norm without prior.  A per-charge unit without ``charge_state`` is
``E_FIELD_REQUIRED``.  The set-up refuses a stopping that is not positive and
finite at an energy the layer's integration reaches (``E_LAYER_STOPPING``) and a
layer the beam does not leave above the lowest energy the model is evaluated at, 1
keV c.m. (``E_LAYER_EXIT_BELOW_MODEL``: use ``kind infinite``).  The worked example is ``tests/units_scales``
(``models/yield_story.azr``).

The stopping of a target.  ``stopping formula "<expr>" <unit>`` is the stopping per
active nucleus (lab, x = E_lab/MeV).  ``stopping evaluation`` takes SRIM-2013's
stopping powers (AZURE2's ``stopping/`` library, ported from LunaSRIM; every result
that uses them must cite SRIM: J. F. Ziegler, M. D. Ziegler and J. P. Biersack,
Nucl. Instrum. Methods Phys. Res. B 268 (2010) 1818) for the target's
``composition`` (element symbols and atoms per formula unit, required), its
``phase`` (``solid`` or ``gas``; left out: gas for a material of one element that is
a gas in its natural state, solid otherwise, so a gaseous compound must say
``gas``) and ``compound_correction`` (SRIM's factor of Bragg's rule, left out 1), for
the projectile of the entrance pair of each data set that names the target (its
mass from ``<levels>``).  ``stopping table "file" MeV|keV eV/1e15atoms/cm2|MeV*cm2/g``
reads a measured stopping of the composition, two columns (lab energy, stopping per
atom of the material), ``#`` comments.  For both, the stopping per active nucleus is
the material's per atom divided by the atom fraction of the active nucleus's element:
for Ta\ :sub:`2`\ O\ :sub:`5` with :sup:`18`\ O active the stopping per formula
unit over 5.  ``scale <s> fixed|fit [prior <s> <width>]`` on the ``stopping`` field
is the stopping-scale nuisance (:math:`\varepsilon \to s\,\varepsilon`): one
parameter ``stopping_scale.<material>.<projectile>`` per material and projectile,
shared by every target of that material (targets that share it write it alike); left
out it is fixed at 1 (a formula has none); ``fit`` without a prior takes SRIM's
stated accuracy for the projectile class (H 4.0 %, He 3.9 %, Li 4.8 %, heavier
5.8 %).  The run's log, ``--validate`` and ``--describe`` print where each target's
stopping comes from, with the citation.  Refusals: ``E_STOPPING_ACTIVE`` (the active
element not in the composition), ``E_STOPPING_PROJECTILE`` (a neutral projectile),
``E_STOPPING_SCALE_CONFLICT``, ``E_STOPPING_DATA`` (the build has no SRIM data files
and neither ``--srim-data-dir=DIR`` nor ``AZURE2_SRIM_DATA_DIR`` names a directory
with them; a table that cannot be read),
``E_STOPPING_RANGE`` (an energy outside a table).  The test is
``tests/stopping_yield`` (README: the numbers against LunaSRIM).

.. _extrapolation-grids:

Extrapolation Grids
^^^^^^^^^^^^^^^^^^^

The ``energies`` field of an ``extrapolation`` record (format 3) says where the
model is computed.  Each form names its unit and frame; ``lab`` energies are
those of the light particle of the entrance pair on the heavy one at rest, as in
the data files, and ``cm`` energies are entrance-channel c.m. energies, taken to
lab energies with the pair's masses from ``<levels>`` (``E_lab = E_cm (m1 +
m2)/m2``, or the relativistic inverse with ``relativistic=1``) so that the
engine computes, and the output lists, the c.m. energies asked for (to rounding,
10⁻¹⁶ relative).  Any energy (angle) unit of the table may be written; the values
are converted to MeV (deg).  The same forms describe an angle grid
(``angles``; an angle grid is not adaptive).

``range lo hi <unit> <frame> step s``
   ``E_k = lo + k s`` for ``k = 0 … ⌊(hi − lo)/s + 10⁻⁹⌋``::

      energies   range 0.02 0.8 MeV cm step 0.002

``range lo hi <unit> <frame> points n [log]``
   ``n`` points from ``lo`` to ``hi``, evenly spaced or, with ``log``, evenly
   spaced in ``ln E`` (``lo > 0``)::

      energies   range 20 800 keV cm points 200 log

``list x1 x2 … <unit> <frame>``
   The energies as written, in that order::

      energies   list 0.025 0.05 0.1 0.2 MeV cm

``file "path" <unit> <frame>``
   One energy per line; ``#`` starts a comment, blank lines are skipped.  The
   path is relative to the project's directory, as a data file's;
   ``AZURE2 --describe`` prints the path it opened and the number of energies::

      energies   file "grids/gamow_window.txt" keV cm

``adaptive lo hi <unit> <frame> tolerance t max_points n start_step s``
   A grid made by the extrapolation run (menu 3) at the run's parameters, for
   a narrow resonance that a uniform grid would need a tiny step everywhere
   for::

      energies   adaptive 0.2 0.4 MeV cm tolerance 1e-3 max_points 2000 start_step 0.01

   The run starts from the uniform grid of ``start_step`` over ``[lo, hi]``
   (``hi`` included) and adds, as seeds, ``E_λ + kΓ_λ`` for ``k = 0, ±1/2,
   ±1, ±2, ±4`` of every level of the entrance pair's compound nucleus inside
   the range (the observed energy and total width at the run's parameters).
   Then, sweep by sweep, the model is computed at the midpoint ``m`` of every
   interval ``[a, b]`` made in the sweep before (all intervals in the first),
   and an interval is split -- ``m`` joins the grid, its two halves are tested
   in the next sweep -- when the model differs from the straight line between
   its ends by more than the tolerance:
   ``|f(m) − (f(a) + f(b))/2| > t · max(|f(m)|, 10⁻³ max|f|)``, the maximum over
   the values computed so far.  Below 10⁻³ of its largest value a curve is thus
   refined to the absolute deviation ``t · 10⁻³ max|f|``, not relative to
   itself.  The midpoint of an interval that is not split is dropped, so every
   interval of the final grid has been tested.  The refinement ends when no
   interval is split, or when a sweep's splits would make the grid longer than
   ``max_points``: then the splits with the largest deviations are made up to
   that count and the run warns with the largest deviation left and where.
   ``f`` is the value the point pass computes, before the output's frame and
   S-factor conversions: the c.m. cross section (or differential cross
   section, phase shift, analyzing power), and for ``quantity sfactor`` the S
   factor.  With an angle grid the energies are
   refined at every angle at once: an energy joins the grid, at all angles,
   when any angle needs it (one grid for all angles).  The run log reports the
   count (``Adaptive grid of extrapolation #1: 251 energies (...)``); the
   output files list the grid's points as for any other grid, and so does the
   data check file (``checks/data.chk``, the Data Check line of
   ``<config>``): it is written after the refinement.

   The refinement sees what the start grid and the level seeds bracket.  A
   structure that is not at a level -- an interference minimum, a threshold
   -- and is narrower than ``start_step`` can lie inside an interval whose
   ends and midpoint are on a straight line, and is then not found: choose
   ``start_step`` to resolve such structures.

   The grid depends on the parameters: a run with other parameters can give
   other energies.  Which intervals are split is decided once per sweep from
   the computed values, so the grid does not depend on the number of threads.
   A single level of width 1 eV at 300 keV (``tests/extrapolation_grids``)
   needs fewer than 400 points for its peak and its integral to 10⁻³; a uniform
   grid of the same count misses the peak by more than half.  pyazr sessions
   (``pyazr.Session``) evaluate at many parameter vectors and refuse an adaptive
   grid; ``pyazr.refine_grid`` runs the same refinement on any Python function,
   and ``pyazr.extrapolation_grid`` gives the energies of the other forms.
