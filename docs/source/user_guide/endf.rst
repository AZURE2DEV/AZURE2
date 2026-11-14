ENDF-6 evaluations
==================

``pyazr.endf`` reads an evaluated file in the ENDF-6 format, says what it holds, gives its cross sections and
angular distributions at the points of a data set of your project, and checks the file's format. It is written
for evaluators: the words are those of the formats manual, ENDF-102 (2023 edition, BNL), and every rule the module
applies is cited by the manual's printed page in its source and in its messages. Each call has a command of
``python -m pyazr`` beside it (:doc:`pyazr_command_line`); ``--json`` gives the same answer as data.

``pyazr/examples/endf_evaluation.py`` does the four things below on a hand-written p + 12C file and the project T1
of the tests, or on a file and a project you name.

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

Anything else (an analyzing power, a phase shift, capture to one final state, an excited light particle) is
refused with what was asked, what is allowed and what to do (``pyazr.FormatRefused``), never drawn as zero.

The conversions are the engine's, not a second copy in Python: the lab energy of the projectile is the engine's
lab energy of each point (or its ``cm_to_lab_energy`` when the file's projectile is the pair's second particle);
the Coulomb terms of LAW=5 (eta, the Rutherford and Mott cross sections, eqs. 6.9-6.12 of the manual) are the rows
of :doc:`physics_functions` on the file's masses (AWI, AWR). Between the incident energies of MF=6 the module
interpolates the value at fixed angle by the section's law; for LAW=5 it interpolates the nuclear and interference
part, sigma - sigma_c, as the manual tabulates it for LTP=12-15, and adds the Coulomb part at the energy asked.

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

What reading refuses and what it reports
----------------------------------------

``pyazr.endf.read`` refuses (``pyazr.FormatError``, with file:line:columns, what was found, what is allowed and
the page) what makes a record unreadable: a control number that is not an integer, a line too short to hold MAT,
MF, MT, a number it cannot read (``1.0D+3``: the letter D is not allowed, p. 31), a count larger than the record's
lines, lines left after the last record of a section it is asked for. It reads, and leaves to the check, what a
reader can live with: columns beyond 80, wrong or missing sequence numbers, a directory that does not match.
