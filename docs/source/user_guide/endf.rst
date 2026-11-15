ENDF-6 evaluations
==================

``pyazr.endf`` reads an evaluated file in the ENDF-6 format, says what it holds, gives its cross sections and
angular distributions at the points of a data set of your project, and checks the file's format. It is written
for evaluators: the words are those of the formats manual, ENDF-102 (2023 edition, BNL), and every rule the module
applies is cited by the manual's printed page in its source and in its messages. Each call has a command of
``python -m pyazr`` beside it (:doc:`pyazr_command_line`); ``--json`` gives the same answer as data.

``pyazr/examples/endf_evaluation.py`` does the four things below on a hand-written p + 12C file and the project T1
of the tests, or on a file and a project you name; ``pyazr/examples/endf_start_model.py`` starts a project from the
ENDF/B-VIII.1 3He + 4He extract and writes T1 as an evaluation.

Read a file and see what it holds
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Python
     - command line
   * - ``m = pyazr.endf.read("h-002_He_004.endf")``
     - ``python -m pyazr endf summary h-002_He_004.endf``
   * - ``print(m.summary())``
     -

The summary names the material (MAT), the target and the projectile (from NSUB), the library line, the lab, the
evaluation date and the authors, LRP in words, and every section: its file and reaction (``(z,p0)``, ``(z,a1)``,
...), its lines, and what it holds -- for MF=3 the number of energies, their range and QM, QI; for MF=6 the frame
(LCT) and each product's law, with LANG or LTP and the incident energies; for MF=2 the format of each range and,
for R-Matrix Limited (LRF=7), KRM, IFG, KRL and the numbers of pairs and spin groups. The descriptive text of
MF=1/MT=451 follows. A file that holds several materials needs ``mat=`` (``--mat``).

``m.cross_section(600)`` is MF=3/MT=600: ``xs(energies_eV)`` gives barns at lab energies, with the section's own
interpolation laws (INT 1-6, the charged-particle law INT=6 included); ``xs.qm``, ``xs.qi`` (eV) and ``xs.table``
(the TAB1 record). ``m.distribution(2)`` is MF=6/MT=2 with its products. A value outside a table is NaN, never zero
and never extrapolated.

Draw an evaluation over your data
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Python
     - command line
   * - ``c = m.curves(s.evaluate(x).data[i], session=s)``
     - ``python -m pyazr endf curves FILE PROJECT --record NAME``

``c.values`` are the file's values at the data set's points (``c.energies``: MeV, c.m.; ``c.angles``: deg, c.m.),
in the data set's quantity and the engine's unit (b, b/sr, MeV b, or a ratio); ``c.coverage`` is False where the
file has no value (below its first energy, for example); ``c.source`` names the sections used; ``c.notes`` what to
know. ``print(c)`` prints the points. Without a session, a ``pyazr.endf.CurveRequest`` gives the quantity,
the points and the two pairs by hand.

Which section gives which quantity:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - the data set
     - the file
   * - elastic ``dsigma_domega``
     - MF=6/MT=2 LAW=5 (LTP 1, 2, 12, 14, 15; identical particles too)
   * - elastic ``ratio_to_rutherford``, ``ratio_to_mott``
     - the same, divided by the engine's Rutherford or Mott cross section
   * - ``sigma``, ``sfactor`` of a particle exit
     - MF=3/MT of the exit particle's family (600+k protons, 800+k alphas, ...), k from the file's
       QM - QI = the residual's excitation
   * - ``dsigma_domega`` of a particle exit
     - MF=3 x MF=6 LAW=2 (Legendre or tabulated): sigma y p / 2 pi
   * - ``sigma`` of capture to every final state
     - MF=3/MT=102
   * - ``sigma``, ``sfactor``, ``dsigma_domega`` of capture to one final state
     - MF=3/MT=102 x the yield of its primary photon in MF=6/MT=102 (found by its energy, QI minus the final
       state's excitation), and its LAW=2 distribution for ``dsigma_domega``

Anything else (an analyzing power, a phase shift, an excited light particle) is
refused with what was asked, what is allowed and what to do (``pyazr.FormatRefused``), never drawn as zero.

The conversions are the engine's, not a second copy in Python: the lab energy of the projectile is the engine's
lab energy of each point (or its ``cm_to_lab_energy`` when the file's projectile is the pair's second particle);
the Coulomb terms of LAW=5 (eta, the Rutherford and Mott cross sections, eqs. 6.9-6.12 of the manual) are the rows
of :doc:`physics_functions` on the file's masses (AWI, AWR). Between the incident energies of MF=6 the module
interpolates the value at fixed angle by the section's law. For LAW=5 with a tabulated distribution (LTP=12, 14,
15) it reads the file as the manual says (p. 146): sigma_NI(E) of MF=3 times P_NI(mu, E) of MF=6, each
interpolated under its own law, plus the Coulomb part at the energy asked. ``m.curves(data, session=s,
elastic_reading="product")`` interpolates sigma_NI(E_i) P_NI(mu, E_i) of the tabulated energies instead; where
the two readings differ by more than 1e-3 the curve's notes say so (it happens between two energies where
sigma_NI changes sign). For LTP=1 and 2 the nuclear and interference part sigma - sigma_c is interpolated.

Check a file's format
---------------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Python
     - command line
   * - ``print(pyazr.endf.check_format(path))``
     - ``python -m pyazr endf check FILE`` (its first part)

The module's own check needs no outside program. It reads every line and the records of MF=1/MT=451, MF=2/MT=151
(LRF=7), MF=3 and MF=6, and reports each finding with its line, its section and the manual's page:

- **errors**: a control number that is not a number; a line without MAT, MF, MT; MAT, MF, MT out of order; a
  section without SEND, a file without FEND, a material without MEND, no TEND; a count (NPL, NR, NP, NZ) larger than
  the lines of the record, or lines left after the last record; the directory of MF=1/MT=451 listing a section that
  is not there, missing one that is, or giving another number of lines; an interpolation law outside 1-6 in MF=3;
  energies that decrease; characters outside ASCII 32-126.
- **warnings**: characters beyond column 80 (undefined, p. 30); wrong sequence numbers; non-zero fields in an end
  record; fields after a record's last value that are neither blank nor zero (a count smaller than what was
  written: the reader cannot see it otherwise); the HEAD's ZA, AWR different from MF=1; LRP=2 (for derived files
  only, p. 42); counts above the limits of Appendix G (p. 391, strict for primary evaluations).
- **notes**: a file without sequence numbers (optional, declared obsolete in 2013).

Run the ENDF checking codes
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Python
     - command line
   * - ``pyazr.endf.checkers()``
     - (listed in the answer of ``endf check``)
   * - ``r = pyazr.endf.check(path)``; ``print(r)``
     - ``python -m pyazr endf check FILE [--program checkr]``
   * - ``c = pyazr.endf.standardize(path)``; ``print(c.summary())``
     - ``python -m pyazr endf standardize FILE``
   * - ``c.accept(to="new.endf")``
     - ``python -m pyazr endf standardize FILE --write new.endf``

pyazr installs none of the checking codes; it finds them. ``checkers()`` looks in the directories of
``$AZURE2_ENDF_UTILITIES`` (CHECKR, FIZCON, PSYCHE, INTER, STANEF of the ENDF utility codes, NNDC and IAEA) and
``$AZURE2_PREPRO`` (PREPRO's FIXUP, DICTIN, ENDF2C, LINEAR, RECENT, D.E. Cullen), then on ``PATH``. Do not put
PREPRO's directory on ``PATH``: it holds a program named ``activate``.

``check(path)`` runs the module's own format check first (``r.own``, named as its own), then CHECKR, FIZCON and
PSYCHE, each in a fresh directory of its own, on a copy of the file, with its input deck (the input file, the message
file, ``Y`` for the whole file, ``DONE``) and a time limit. Each report is read into findings: MAT, MF, MT, the
code's record or sequence number, the code's own text, and for messages pyazr knows what they mean for an evaluator
(``finding.meaning``). The raw report and the terminal text are kept (``r.reports[i].report``). A crash, a time-out
or a report that ends early is a finding of the run. The exit status of the codes says nothing about findings: they
exit 0 with errors in the report.

What the codes say about charged-particle files, and what it means:

- *LRP=2 ONLY ALLOWED IN DERIVED FILES* (CHECKR): the rule of p. 42 of the manual. The ENDF/B-VIII.1
  charged-particle files with R-matrix parameters use LRP=2 for a full MF=3 beside MF=2.
- *NE = 400 OUT OF RANGE 1 - 200* followed by *SECTION CANNOT BE CHECKED* (CHECKR): CHECKR's own limit of 200
  incident energies; the manual allows 2000 in MF=6 (Appendix G). The section is not checked, not wrong.
- *AP NOT IN RANGE* in MF=2 (PSYCHE): PSYCHE reads the record after the energy range as the spin and radius of the
  older formats; in R-Matrix Limited that record holds IFG, KRM, NJS, KRL. Not a finding about the file.
- *INVALID ANGULAR REPRESENTATION LTP=15* (CHECKR 8.24): the 2023 manual defines LTP=15 (p. 147); CHECKR does not
  know it yet.
- *AWI SHOULD BE SET TO ...* (FIZCON): AWI is the projectile's nuclear mass in neutron masses (p. 24), not its atomic
  mass.
- *CHECK NORMALIZATION* (FIZCON): a tabulated distribution whose integral over the cosine is not 1.

What "correct" does and does not do
-----------------------------------

``standardize(path)`` runs a code that rewrites a file -- STANEF (tape label, directory, sequence numbers, numbers
in the classic form), or PREPRO's FIXUP, DICTIN, ENDF2C -- on a copy, and returns a candidate. Nothing replaces the
original: ``c.diff()`` lists per section what changed (bookkeeping, numbers written another way with every value
equal, numbers changed: how many, the largest relative change and how many digits a number lost), and
``c.accept(to=...)`` writes the candidate where you say (by default beside the original as
``<name>.<program>.endf``; over an existing file only with ``overwrite=True``).

Read the difference before keeping a candidate. STANEF writes every number in the classic form of 7 significant
digits, MF=2 included: on a file written with 9-10 digits (as FUDGE writes them) it rounds the resonance
parameters, which moves an R-matrix fit's chi-squared (about one unit in 10 000 on the projects measured in the
design note). The built STANEF 8.10 also rewrote a number written without a decimal point (``35070343``, an
energy in eV of the ENDF/B-VIII.1 3He + 4He file) as ``3.507034+3``, 10 000 times smaller: the difference shows it
as a relative change of 1. No program changes a parameter to make a check pass, and pyazr runs no program you did not
ask for.

Trying the codes on your machine::

   export AZURE2_ENDF_UTILITIES=/path/to/endf-utility-codes/bin AZURE2_PREPRO=/path/to/prepro/bin
   python -m pyazr endf check h-002_He_004.endf
   python -m pyazr endf standardize h-002_He_004.endf

Start a project from an evaluation
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Python
     - command line
   * - ``r = pyazr.endf.read(path).resonances(isotope=0, energy_range=0)``
     - ``python -m pyazr endf summary FILE`` (what MF=2 holds)
   * - ``project, report = pyazr.endf.to_project(r, "new.azr", material=m)``
     - ``python -m pyazr endf to-project FILE NEW.azr [--separation-energy S] [--json]``

``Material.resonances()`` reads MF=2/MT=151 in the R-Matrix Limited format (LRF=7) into a record in the file's own
units and frame (ER in eV, the laboratory energy of the incident particle; GAM in eV^1/2; masses in neutron masses;
radii in 1e-12 cm); ``print(r)`` lists the pairs, spin groups, channels and levels. ``to_project`` makes a new
format-3 project of it (never an existing file) in the Brune basis, every level and width fixed:

- the incident pair is the one with MT=2 (p. 71); its separation energy from the compound nucleus comes from
  AME2020 (the file holds no mass of the compound nucleus), or ``separation_energy=`` in MeV; every other pair takes
  the file's Q value, and the residual's excitation from MF=3's QM - QI of its MT (p. 105) when the file has that
  section. Masses are the file's; the report compares each with AME2020 (the 3He + 4He file of ENDF/B-VIII.1 gives
  6Li 21.7 keV below AME2020, and the engine warns that separation energies and masses disagree);
- ER and GAM go to E_x and the c.m. scale with E_lab / E_cm = (m_A + m_B)/m_B of the incident pair (pp. 66,
  355-356); the file's boundary conditions BND are kept: the engine turns the formal levels with those B_c into
  Brune energies, widths (eV) and ANCs (fm^-1/2), the observed form a project holds;
- a spin group without levels is not written (AZURE2 computes the hard-sphere phases itself); a level above EH is
  marked as background in the report.

On the 3He + 4He extract of the tests the project reproduces the file's own MF=3/MT=600 (4He(3He,p)6Li) to 1.1e-6
above 7.2 MeV (median 7e-8); the two points just above the p + 6Li threshold do not follow from the file's MF=2.

A file with IFG=0 holds widths at the level energy instead of amplitudes: ``to_project`` converts them as the
manual says (p. 73, eq. D.45), gamma = sign(GAM) sqrt(abs(GAM) / 2P), with P the engine's penetrability of the
channel at the level, at abs(ER) for a negative ER (the channel's c.m. energy abs(ER) m_B/(m_A + m_B) + Q,
eq. D.82); a width on a channel closed there cannot be converted and is refused.

Refused, each with what, what is allowed and how to fix it: a range that is not LRF=7 or has NRO=1, KRM other than
4 (Reich-Moore and Breit-Wigner files are not read as a start model; a Reich-Moore file's eliminated widths would
need a conversion to the Brune basis that is not built), KRL=1 (no equations in the manual), a photon pair (with
KRM=4 its channels enter the level matrix; AZURE2 keeps radiative widths out of the total width), a fission pair,
KBK > 0 or KPS > 0 (background R-matrix, tabulated phases), APE different from APT or two radii in one pair (AZURE2
has one radius per pair), SHF=0 or PNT=-1 on a particle channel. The manual does not say in which order the
channel spin couples the two spins; the module reads
SCH as AZURE2's (light particle first). No cross section or angular distribution an ENDF-6 file holds depends on
that order: exchanging the two spins changes each channel by a sign, and unpolarized cross sections sum over the
channel spins incoherently (measured on the ENDF/B-VIII.1 3He + 4He file: every amplitude of one channel spin of p +
6Li flipped, sigma and dsigma/dOmega of 4He(3He,p)6Li unchanged to the last bit). Analysing powers do depend on it:
a project made from an evaluation and fitted to analysing-power data rests on AZURE2's order, which the conversion
report says when a pair has two non-zero spins.

Write an evaluation
-------------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Python
     - command line
   * - ``record, losses = pyazr.endf.from_session(s, x)``
     - (the MF=2 that would be written, and what it cannot carry)
   * - ``w = pyazr.endf.write(s, "new.endf", x, header={"MAT": 625, "AUTH": "..."})``
     - ``python -m pyazr endf write PROJECT NEW.endf [--header h.json] [--low E --high E] [--json]``

``write`` writes the session's model at the fit vector ``x``:

- **MF=1/MT=451** (pp. 38-42): LRP=2, the projectile's AWI (its nuclear mass in neutron masses, p. 24), EMAX, NSUB;
  the three standard text lines and HSUB from ``header`` (MAT, NLIB, NVER, LREL, NMOD, ALAB, EDATE, AUTH, REF, DDATE,
  RDATE, LIBRARY, TEXT: the evaluator's fields; MAT defaults to 9999); a text that says what wrote the file, that
  MF=3 holds the engine's cross sections and that MF=2 reproduces them, the photon widths of the levels, and every
  loss; the directory.
- **MF=2/MT=151**: the particle channels in the formal basis (LRF=7, KRM=4, IFG=1), with the engine's boundary
  values B_c = S_c(E_1) of each spin group written as BND, ER and GAM on the laboratory scale of the incident pair.
  Photon channels are not in MF=2 (the PI's decision): capture is MF=3/MT=102.
- **MF=3** of every particle exit (MT 600 + k, 650 + k, ..., 800 + k) and of the total capture (MT=102, internal
  and external capture), on one grid of laboratory energies made by ``pyazr.refine_grid``, the engine's refinement
  of an adaptive extrapolation grid: an interval is split while its midpoint differs from linear interpolation by
  more than ``tolerance`` (default 1e-3) times the value, or times 1e-3 of the largest value where the cross section
  is smaller than that. The refinement starts from 41 even points, every threshold of an exit and every level energy
  of the project with the two energies half its width away (so a level narrower than the spacing is found), and
  the quarter points of every interval are checked by the same rule afterwards (a midpoint can lie on the chord of
  a steep flank). INT=2. On T1, T2 and 13N linear interpolation in the written file gives the engine's value on
  20001 points per section to 1.0e-3 by that rule.
- Before the file is written, MF=2 is read back into a project and its particle reactions computed on the grid: they
  must equal MF=3 to 1e-5 (T2: 8e-8), else the write is refused.
- **MF=6** (pp. 133-147), in the c.m. (LCT=2), on the energy grid of MF=3 (the grid's rule also covers the
  distributions): each particle exit as LAW=2 (the engine's Legendre coefficients, normalised as eq. 6.8) with the
  residual as LAW=4; capture to each final state as its primary photon (LAW=2, its energy at zero incident energy in
  AWP, the yield sigma_f/sigma) and the recoil (LAW=4). The cascade below the final states is not given (AZURE2 does
  not hold it). No MF=6 section has more than 200 energies and no LAW=5 panel more than 201 cosines (the limits of
  CHECKR, which other codes share; the PI's choice): where the tolerance would need more, the section is thinned to
  the limit by the same rule and MF=1, ``w.accuracy`` and the end of ``w.losses`` say the tolerance it reaches.
  ``mf6_limits="manual"`` writes to the format's own limits instead (2,000 energies, Appendix G of the manual).
- **Charged-particle elastic scattering** (MF=6/MT=2 LAW=5, pp. 145-147), in one of two forms
  (``elastic_form``; the default is ``pyazr.endf.writer.ELASTIC_FORM``):

  - ``"amplitudes"`` (the default): the nuclear amplitude expansion, LTP=1, with 1.0 in MF=3/MT=2 (p. 145). At
    each energy the file holds the engine's coefficients b_l and a_l of eqs. 6.13 (distinguishable particles) and
    6.14 (identical particles), from which the manual's formula gives the engine's elastic cross section at every
    angle (measured to 1e-11 on T1, T2, 13N, p + p and d + d). The coefficients come from the collision matrix
    of the engine's evaluation (the point output ``nuclear_amplitudes`` of an elastic data set). Between two
    energies a reader interpolates sigma_e - sigma_c at a fixed angle; the energies are chosen by the grid's rule
    on the cross section that interpolation gives at 17 angles. Measured within CHECKR's limits: T1 200 energies,
    tolerance stated 6.4e-3, largest deviation 4.3e-3 between the written energies (median 9e-5); T2 117 energies,
    1.1e-3 stated, 1.8e-4 measured; 13N 200 energies, 6.0e-3 stated, 3.4e-3 measured. With
    ``mf6_limits="manual"`` T1 and 13N take 491 and 465 energies and reach 1.1e-3.
  - ``"tabulated"``: LTP=12, P_NI tabulated in mu (each energy on the cosines it needs) and sigma_NI in MF=3/MT=2
    (p. 150), from mu = -1 (0 for identical particles) to ``elastic_max_cosine`` (0.98; at smaller angles the file
    gives the Coulomb term alone). Elastic scattering has energies of its own, chosen so that a code reading the
    file as the manual says (sigma_NI of MF=3 times P_NI of MF=6, each interpolated) gets the cross section to the
    tolerance: wherever sigma_NI changes sign, the two written energies around the change are at most 1e-5 x E
    apart, and between them that reading is wrong by up to the nuclear part of the cross section. MF=1 and the
    answer (``w.elastic``) list these brackets, the number of energies and the tolerance reached. Within CHECKR's
    limits the 200 energies cannot follow sigma_NI near its zeros: on T1 and 13N the section, read as the manual
    says, is not usable there (tolerance stated 5.5 and 1.3; median error 1e-3), on T2 it holds 0.12.

  On T1, T2 and 13N the file's dsigma/dOmega equals the engine's at the data's points and between the tabulated
  energies and angles to 3e-3 (LAW=2 and photons) and 5e-3 (elastic scattering) of the cross section, or to the
  tolerance the file states where CHECKR's limits thinned a section.

- **MF=32** (pp. 253-259), when asked (``write(..., covariance="hesse")``): the covariance of MF=2's parameters
  (ER and GAM as MF=2 lists them, eV and eV^1/2) in the general form LCOMP=1, one block. It is HESSE (iminuit) of
  the fit's objective in a session of the formal basis opened at ``x`` (the fit written as a project and opened with
  ``use_brune=False``); the norms and other parameters of the fit and the photon channels are marginalised. A
  resonance whose parameters were all fixed is not listed; a fixed parameter of a listed resonance has a zero row.
  The formal session holds the project's fixed parameters fixed in the formal basis, which is not the fit's
  constraint where a level energy or width is fixed (on T1 the two covariances differ by up to 3.9 sigma_i sigma_j;
  ``docs/upgrade/endf.md`` 13.13). ``Material.covariance()`` reads MF=32 back: names ``"<J^pi>#<k>.ER"`` and
  ``"<J^pi>#<k>.GAM<c>"``, values, matrix, ``sigma``, ``correlation``. The fit's covariance in its own basis
  (``EngineFitResult.covariance``) is refused as ``covariance=``: MF=32 holds the formal parameters.

Listed in MF=1 and in the answer (``w.losses``, ``w.not_written``), each only when the project has it: every
section a reaction of the project needs and the file does not hold, with the reason; the photon channels (their
widths are listed in the text); the data sets with their records; the Trojan-horse experiments; beta-decay pairs;
an excited light particle. Refused: a session in the formal basis, Reich-Moore
(``use_rmc``), relativistic kinematics, a nuclear potential, several compound nuclei, an exit pair outside ENDF's
reaction families, a session without an exit. A cross section far below 1e-3 of its largest value in the energy
range is good to the absolute tolerance only: give a lower range (``energy_range_MeV``) for sub-barrier capture.

Checking an evaluation's resonance parameters
---------------------------------------------

``pyazr.endf.check_resonances(m)`` (or ``pyazr.endf.check(path, resonances=True)``, ``pyazr endf check FILE
--resonances``) does for a file with MF=2 in LRF=7 what the checking codes do not (they count the records of
LRF=7): MF=2 read into a project and the engine's cross section of every particle exit at the energies of its MF=3
section, against MF=3 (per section the largest and median relative difference and the points off by more than
1e-4, and points where MF=3 lies orders of magnitude above what MF=2 gives); each nuclide's mass against AME2020,
each pair's Q against AME2020's and against the Q the file's own masses give; every flag the MF=1 text names
(LRF, KRM, IFG, LTP, LRP) against the record that holds it. Each finding has a severity, the section and the
manual's page, as the format check's. On the ENDF/B-VIII.1 3He + 4He file it finds the first point of MF=3/MT=600
37 orders of magnitude above its own MF=2, the 6Li mass 21.7 keV below AME2020 while the Q value agrees with
AME2020 (the two from different tables), and a text that says LTP=2 where the record holds LTP=1.

What the checking codes say about a written file
------------------------------------------------

Measured with CHECKR 8.24, FIZCON 8.19 and PSYCHE 8.11 on the files written from T1, T2 and 13N:

- *LRP=2 ONLY ALLOWED IN DERIVED FILES* (CHECKR): expected, the PI's decision (the library's charged-particle LRF=7
  files have the same).
- *MATERIAL NUMBER SHOULD BE 625* (CHECKR, warning, for p + 12C): the default MAT 9999; the MAT is the evaluator's (``header``).
- *NE = ... OUT OF RANGE 1 - 200* in MF=6, then *SECTION CANNOT BE CHECKED* and *OUT OF SEQUENCE* in the next
  section (CHECKR), and *NL = ... OUT OF RANGE 1 - 201* for a LAW=5 panel: CHECKR applies File 5's limit of 200
  incident energies and File 4's 201 angular points to MF=6; the manual's limit for MF=6 is 2,000 energies
  (Appendix G, p. 391). A file written with the default limits has neither; one written with
  ``mf6_limits="manual"`` may. The library's files have the first.
- *DISCRETE 2-BODY LAW NOT PERMITTED FOR MT= 102* (FIZCON): the manual allows LAW=2 for the primary photons of
  MT=102 (p. 143); FIZCON's list of reactions omits it.
- *CHECK NORMALIZATION= 0.998889* (FIZCON, a LAW=5 LTP=12 panel where sigma_NI is close to zero and P_NI large; 81
  of the 200 panels of T1 within CHECKR's limits): FIZCON sums in single precision, and where sigma_NI is near zero
  its sign and size are decided by the cosines a panel keeps. The amplitude form has no such panel.
- *AP NOT IN RANGE* and *SPIN OF ISOTOPE ... DISAGREES WITH ...* (PSYCHE): PSYCHE's reading of LRF=7 (above).
- *THE CALCULATED Q 6.67452E+06 DISSAGREES WITH THE GIVEN Q 1.94349E+06* for MT=102 (PSYCHE, note): PSYCHE's Q
  comes from its own mass formula (its function ``AMASS``: a liquid-drop formula with shell terms; a table only for
  n, p, d, t, 3He and alpha), which gives 6.67 MeV for 12C + p -> 13N; the file's QM = QI = 1.94349 MeV is the
  AME2020 proton separation energy of 13N. PSYCHE's number, not the file.
- *THERMAL SIGMA* (negative) and *RESONANCE INTEGRAL* (PSYCHE, notes): PSYCHE's neutron quantities, computed
  below the first energy of the table (0.12 MeV), where the file has no value.

On the page
-----------

The web interface's **ENDF** tab (:doc:`web_interface`) makes the same calls: the files of the project folder and
what each holds (``summary()``), every data set with the files' values, the data and the model and the ratios
(``curves(data, session=s)``), the codes found (``checkers()``), **Check** (``check(path, codes=..., resonances=True)``) and
**Correct** (``standardize(path, program=...)``, then ``accept(to=...)`` or, after a second question,
``accept(overwrite=True)``) as jobs.  **New project from a file** shows ``read(path).resonances(isotope=,
energy_range=)`` (or its refusal) before anything is written, runs ``to_project(record, path,
separation_energy=..., material=m)`` as a job into a directory of its own and shows the report (``print(conversion)``
in Python), then writes the project as a new file of the folder; the page serves one project, so **Open … here**
replaces the open one (``pyazr ui new.azr`` opens a page of its own).  **Export this project** shows
``from_session(s, x)[1]`` (the losses and the photon widths) before the run and runs ``write(s, path, x,
header={...}, energy_range_MeV=..., tolerance=..., number_form=...)`` as a job; the file then appears among the tab's
files, ready for **Check**.  The model is drawn on the grid of ``Session.model_curve`` and each file's values on the
same grid (``curves(curve, session=s)``: a curve is taken as a data set is).

What reading refuses and what it reports
----------------------------------------

``pyazr.endf.read`` refuses (``pyazr.FormatError``, with file:line:columns, what was found, what is allowed and
the page) what makes a record unreadable: a control number that is not an integer, a line too short to hold MAT,
MF, MT, a number it cannot read (``1.0D+3``: the letter D is not allowed, p. 31), a count larger than the record's
lines, lines left after the last record of a section it is asked for. It reads, and leaves to the check, what a
reader can live with: columns beyond 80, wrong or missing sequence numbers, a directory that does not match.
