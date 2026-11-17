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
