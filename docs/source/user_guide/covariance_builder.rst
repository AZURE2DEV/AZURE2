The covariance builder
======================

A data set whose record ends with ``covariance <file>`` enters :math:`\chi^2` with the covariance of its points
(:doc:`chi_squared`, "Correlated data: covariance matrices"; the file format in :doc:`../reference/data_formats`).
The file is a sum of terms the engine reads; the builder writes it from named *parts*, checks it as the engine will,
and shows what the engine will make of it.  It is the module ``pyazr.covariance`` and, on the web page, the
**Covariance…** button of a data set in the Data tab.

What the engine does with the file
----------------------------------

- :math:`\chi^2 = v^T C^{-1} v` over the points the data sets fit, :math:`v` the model minus the data in the
  engine's unit of the quantity (b, b/sr, ...: the value column after the data set's column units).  Rows outside a
  data set's window are in the file (its ``file`` line counts every data row) and drop out.
- The ``diagonal``, ``block`` and ``matrix`` terms on the fitted rows must be positive definite by themselves: every
  fitted row needs a variance of its own.  A fully correlated part (``vector``, ``relvector``) cannot stand alone.
- A ``relvector`` term scales with the **model**, frozen at the first evaluation of a run (``covarianceModel.out``),
  which avoids Peelle's pertinent puzzle; a ``vector`` term is absolute (made here from the data values).
- The error column of a fitted row is replaced by :math:`\sqrt{C_{ii}}`, the data value standing in for the model;
  the run warns when the two differ by more than 1 %.  :math:`\chi^2` does not depend on the column.
- A data set's normalization :math:`n` scales the data and the covariance together (as it scales the error column):
  a normalization fitted with a prior of :math:`\epsilon` % and a fully correlated scale part of :math:`\epsilon` %
  are the same uncertainty counted twice.  Keep one; the builder says so when both are there.

The parts
---------

=================  ==============================================================  ==============
kind               what it is                                                      term written
=================  ==============================================================  ==============
``own``            the points' own, uncorrelated errors: the error column            ``diagonal``
                   (``source=column``), a percentage of the value
                   (``source=percent``), or what the column holds beyond the
                   other parts (``source=remainder``: the column is the total)
``scale``          a common scale uncertainty of ``percent`` %, fully correlated     ``relvector``
                   over every row (``basis=model``) or of the data values            or ``vector``
                   (``basis=data``)
``range``          the same on the rows with ``e_min`` ≤ E_lab ≤ ``e_max`` (MeV)    as ``scale``
``length``         ``percent`` % of the value, correlated as                        ``block``
                   :math:`\exp(-|E_i-E_j|/L)` (``shape=exponential``) or
                   :math:`\exp(-(E_i-E_j)^2/2L^2)` (``gaussian``), ``length``
                   :math:`L` in MeV (lab), optionally on ``e_min`` to ``e_max``
``correlation``    a given correlation matrix (``file``: N × N over every row       ``block``
                   covered) with the error column or ``percent`` % of the value
=================  ==============================================================  ==============

The builder refuses a part it cannot make (an unknown field, a remainder below zero, a correlation file of the
wrong size or not a correlation matrix), and a matrix the engine would refuse: no own variance on a fitted row, or a
part that is not positive semi-definite (named, even where the other parts would hide it).  The parts are written on
the file's ``#pyazr.covariance`` comment line, so the builder reopens with them (``read_parts``); a file without the
line was written by hand.

In Python
---------

.. code-block:: python

   from pyazr import Project, covariance
   p = Project.open("tests/data_covariance/plain.azr")       # the data in tests/angle_mapping/data
   c = covariance.build(p, ["p12c_60"], [{"kind": "own"},
                                         {"kind": "scale", "percent": 5},
                                         {"kind": "length", "percent": 3, "length": 0.2}])
   print(c.refusal, c.checks["positive_definite"], c.notes)
   print(c.total_error[:3], c.error_column[:3])                # sqrt(C_ii) beside the file's errors
   c.write("cov_p12c_60.txt")
   p.set("p12c_60", covariance="cov_p12c_60.txt")

``build(project, datasets, parts)`` returns a ``Covariance``: ``rows`` (every row covered, as the engine reads the
data files: ``_azure2.data_rows``), ``terms``, ``variance`` and ``total_error`` per row, ``matrix`` and
``correlation`` on the fitted rows, ``checks`` (``positive_definite``, ``own_variance``, ``smallest_eigenvalue``,
``sigma_mismatch`` and ``worst_mismatch`` as the engine counts them), ``refusal`` with ``refusal_parts``, ``notes``;
``write(path)``, ``text()``, ``summary()`` (plain data).  ``rows_of(project, datasets)`` gives the rows alone.
``pyazr/examples/covariance_builder.py`` runs the case below.

A worked case
-------------

``p12c_60`` of ``tests/data_covariance/plain.azr`` (61 rows of :math:`d\sigma/d\Omega` at 60°, 0.40 to 2.20 MeV,
errors 3 %), with its own errors, a 5 % model scale and 3 % correlated over 0.2 MeV: the first row's total error is
0.304028 b/sr against the column's 0.139092 (all 61 rows differ by more than 1 %, the worst by 54.25 %, which the
engine's warning repeats), and the engine's :math:`\chi^2` of the data set with the written file, 9.38942514115,
equals :math:`v^T C^{-1} v` formed in Python from a session without it (``tests/pyazr/covariance_builder_test.py``).
The data set's normalization is fitted with a 5 % prior: the builder's notes say that the 5 % scale part counts it a
second time.

On the page
-----------

**Covariance…** on a data set opens the builder (``#/covariance/<name>``): the parts as rows (kind, its fields with
units and docs, **Remove**, **+ part**), the data sets the covariance joins, **Build**, the total error per point
beside the file's (a plot and a table), the correlation on the fitted rows as a heat map, the refusal with the parts
it names, the notes, and **Write and use** with a proposed file name (``cov_<name>.txt``).  Writing sets the
``covariance`` field of the data sets and writes the file as one edit: **Undo** removes the file (or puts back the one
it replaced) and the field together.
