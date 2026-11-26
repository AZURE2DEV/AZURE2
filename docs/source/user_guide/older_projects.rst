Coming from an older AZURE2
===========================

This AZURE2 reads project format 3 only (:doc:`../reference/format3`): every
data set, extrapolation, target, beam, detector, shared normalization, energy
scale, model datum and hybrid potential is a named record.  A project written by
an earlier AZURE2 or by the Qt GUI (format 1, or format 2 with ``format=2`` in
``<config>``) keeps its data in ``<segmentsData>``, ``<segmentsTest>``,
``<targetInt>``, ``<abundances>`` and ``<scatteringLengths>``, numbered by line.
Such a file is converted once; the original is kept.

What happens with an old file
-----------------------------

Every command that reads a project -- a calculation or fit, ``--validate``,
``--describe``, ``--print-predictions``, a pyazr session or
``pyazr.Project.open`` -- refuses a file that is not format 3 before reading
anything of it.  The command line exits with status 3 and ends with::

   ERROR: T2.azr is a format-1 project (<segmentsData> at line 53). This AZURE2 reads format 3 only.
          Convert it:  AZURE2 --convert T2.azr T2_f3.azr     (or pyazr.convert_azr("T2.azr", "T2_f3.azr"))
          The converter keeps T2.azr unchanged and reports every change it makes.

pyazr raises ``pyazr.FormatError`` with the same text.  A file that says
``format=3`` but still holds an old section is refused the same way (convert
the original instead), and so is a file without ``<config>`` (an empty file,
or not a project).

Converting a project
--------------------

.. code-block:: bash

   AZURE2 --convert T2.azr T2_f3.azr --json-report T2_f3.json

or, from Python, ``pyazr.convert_azr("T2.azr", "T2_f3.azr",
report="T2_f3.json")``.  The converter reads the data files beside the old
project (to check the columns) but never writes one; it writes the new project,
a converted copy ``<name>.f3<ext>`` of a prior-covariance file, and the report,
and overwrites nothing.  The new file is validated before it is written.  The
results of the converted project are those of the old one: every converted test
project of the repository gives the old engine's outputs file for file
(``docs/upgrade/format3_tree.md``).

A line the old grammar cannot read
----------------------------------

Old files were often edited by hand or written by other versions of AZURE2.  The
converter reads, each with a note in the report saying what it assumed, the
variants whose meaning is unambiguous: a UTF-8 byte-order mark before the first
tag, Windows line ends, tabs between the fields, a comment line starting with
``#`` in ``<segmentsData>``, ``<segmentsTest>`` or ``<targetInt>`` (skipped), a
``<targetInt>`` segment list or equation with blanks inside its quotes
(``"8, 9"``), an unquoted data-file name with a blank when that file exists
beside the project, and a ``<targetInt>`` line that ends after its empty
stopping-power equation (its parameter count can only be 0).

Any other line it cannot read stops the conversion, and the refusal says where
the line is and what could be meant, for example::

   ERROR: --convert: p.azr is not converted: data segment #3 ends with 'data/meyer_84.dat 0 0', ...
   line 69 of the file, in <segmentsData>:
       1 1 1 0 5 84.3 84.3 1 1 0 5 0 0.005 data/meyer_84.dat 0 0
     expected: isActive entranceKey exitKey minE maxE minAngle maxAngle isDiff [J l, for isDiff 2] norm varyNorm
               normError [energyShift [shiftError varyShift]] dataFile [keyword blocks]
     possible readings (none is applied by itself; take one, or edit the line, in a copy of the project and convert the copy):
       1) the energy shift's vary flag is missing: 0 (the shift fixed):
          1 1 1 0 5 84.3 84.3 1 1 0 5 0 0.005 0 data/meyer_84.dat 0 0
       2) the normalization's error is missing: 0 (no prior); the next three numbers are the energy shift, its error and its vary flag:
          1 1 1 0 5 84.3 84.3 1 1 0 0 5 0 0.005 data/meyer_84.dat 0 0

The readings the converter proposes are a missing normalization error or shift
vary flag written as 0, a missing ``isDiff`` of a test line as 0 or 1, decimal
commas read as points, and an unquoted file name with blanks quoted.  None is
taken without you: a normalization, an energy shift or a data type changes the
fit.  A data line without a data file, or with a number such as ``0,005`` where
the file's name stands, is refused (it has no reading the converter can
propose).

The Qt program (:doc:`graphical_interface`) shows such a line in a window: the
file and line number, the line itself in an editable field, the fields expected
and one button per reading.  ``Convert again`` writes a copy of the old project
beside it with the line as edited (``<name>.repaired.azr``; the old file is
kept) and converts the copy; the lines repaired and the converter's notes are
listed in the window's notes bar.

Files of every AZURE2 version
-----------------------------

The converter reads the forms the project file has had since before 2020; each
form below converts with the meaning the program of its era gave it, or is
refused with the reason (the timeline and the fixtures:
``tests/convert/fixtures/history/README.md``).

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Form (era)
     - What the converter does
   * - ``<levels>`` rows of 30 columns and an ``<externalCapture>`` section
       (before 2020)
     - Converts: the external-capture mask of each active line set on its
       exit pair, as the 2020 setup program imported such files (``minJ`` and
       ``maxJ`` were not used by that import).  An inactive line is refused.
   * - The 2020 file (``origin/qt5``): eleven ``<config>`` lines with full
       paths, ``<lastRun>``, data lines whose file name is the rest of the line
     - Converts; also with Windows line ends, tabs, blank lines, empty
       sections, paths with blanks or backslashes, an unclosed ``<lastRun>``
       at the end of the file.
   * - Target lines with convolution coefficients but no equation
       (2025-01-26 to 01-30)
     - Converts: sigma(E) = c0 + c1 E + ... written as the equation.
   * - Energy shifts (2025-08), composite (component) segments with two,
       three or four numbers per component, secondary decays, scalings
       (2025-08 to 2026-05)
     - Converts, the layout read as the setup program of the era read it.
   * - A number as wide as its column run into the next (written by the
       setup program between 2026-09-25 and 10-07)
     - Refused: two readings; a data, test or target line shows the split as
       a proposed reading.
   * - Polarization times cross section, isDiff 8 (2026-09-17 to 09-28);
       dev's ``mixture`` block (2026-10-12)
     - Refused: format 3 has no quantity or record for them; the reason says
       what to do instead.
   * - dev's spin observables written as type 9 with an operator (until
       2026-11-11)
     - Converts as type 10.
   * - External capture into a final state that is no level of
       ``<levels>`` (run by AZURE2 before 2026-10-26, which added a hidden
       level with a free width)
     - Refused at conversion: add the level to a copy of the old project (or
       switch that pair's external capture off) and convert the copy.
   * - Data files not at the path the project names (another computer, a
       moved folder)
     - Converts with a note per file; the Qt program lists them in one note
       with **Locate Data Files...**, which finds each file by its name in a
       folder you choose.

The conversion and the reading of a project do not depend on the language
settings of the computer: a decimal-comma locale (German, Italian, Polish, ...)
gives the same files as any other.

What the report says
--------------------

The report (``--json-report``, or the dict ``convert_azr`` returns) lists:

* **notes**: meanings the old file left implicit and the new one states, each
  also as a ``# CONVERTED:`` comment above the field (a thick target's numbers
  per incident particle, a Gaussian convolution as a beam spread, an effect
  that also folded the extrapolation with the same line number), and values the
  old engine did not use, which are dropped (a group member's own norm or
  shift).
  One note is about a value that is dropped with ``<lastRun>`` (the Qt GUI's
  run settings): when the old program's check box "use GSL Coulomb functions"
  was ticked (bit 12 of the first number of ``<lastRun>``), the report says that
  the old program computed the project with GSL's Coulomb functions, that this
  program has its own only, and that results inside the Coulomb barrier can
  differ from the old ones by up to about 5 % in a penetrability
  (``docs/upgrade/coulomb.md``, sections 8.3 and 9.1).  The note is in the report
  only; the converted project is the same with or without it;
* **questions**: where the old meaning cannot be written without a choice (a
  fitted temperature shared by lines that become two targets).  Answer them with ``--answers A.json``
  (``{"<record>.<field>": value}``); without answers the conversion stops, or,
  with ``--allow-questions``, writes ``?`` (the file then does not validate
  until it is answered);
* **parameters**: the old parameter names and their new ones
  (``segment_3_norm`` -> ``s03.norm``, ``segment_2_energy_shift`` ->
  ``s02.shift.offset``, ``density_1`` -> ``<target>.areal_density``).  A data
  set is named after its data file; the level parameters (``energy_<k>``,
  ``width_<k>_<c>``) keep their names.

Two pairs with the same particles and states (``12C+12C`` as pair 1 for THM
data and pair 6 for direct data) are told apart by their ``<levels>`` keys: the
converter writes ``entrance_pair`` and ``exit_pair`` on a data set or an
extrapolation whose reaction label names two pairs.  On a sum or a ratio the
keys apply to every term whose label names two pairs, so the old line's terms
must take the same pair there; a line whose terms take two such pairs (one term
from pair 1, another from pair 6) has no format-3 form and is refused, with the
advice to split it into one data set per pair.

The converter refuses a line the old engine refused (with the old engine's
reason), so a file that did not run before does not run after.  It also refuses
what it would otherwise drop or copy through unread, since the converted project
would then compute as if it were absent:

* a section it does not convert: a section of a newer AZURE2 (its ``dev``
  branch) that this build does not compute yet is named as such (the ``dev``
  branch's ``<channelRadii>`` lines are converted, one ``channel_radius``
  record a line, and its ``<isospin>`` lines into ``isospin_*`` records,
  :doc:`levels_channels`); any other unknown section is refused as one that
  no AZURE2 the converter knows wrote.  Keep running such a project with the
  program that wrote it, or convert a copy without the section if its numbers
  may change;
* a data or extrapolation line of a newer program's type: a type followed by a
  spin operator (``b:y``, ``t:y``, ``b:y*t:y``, ...) is a general spin
  observable of the ``dev`` branch, which this build does not compute.

The engine refuses the same sections in a format-3 file (``E_UNKNOWN_SECTION``),
and a ``<config>`` key it does not have (``E_UNREAD_LINE``).

Parameter files
---------------

A ``param.sav``, ``param.par`` or ``param.fit`` written for the old project
names its data-set and target parameters by the old names.  Given as the start
of a run of the converted project, such a file applies only its level
parameters, and the run says so (``W_PARAMETER_FILE_NAMES``, or
``E_PARAMETER_FILE_NO_MATCH`` when no name matches).  Rename it with the
conversion's map:

.. code-block:: bash

   AZURE2 --convert-params param.sav param_f3.sav --project T2.azr

``--project`` names the *old* project (the map is that of its conversion);
``pyazr.convert_params("param.sav", "param_f3.sav", "T2.azr")`` does the same.
Every other byte of the file is kept.  The converted T2 test project calculated
from its old fit's renamed ``param.sav`` gives that fit's χ²
(``tests/convert_params``).

The Qt GUI still reads and writes the old format; convert what it writes
before running it with this AZURE2.
