The prior builder
=================

A prior on the fit parameters -- a multivariate Gaussian over several of them, ``priorCovariance=<file>`` in
``<config>`` (:doc:`chi_squared`, "Prior parameter covariance") -- is a file of names, centres and a covariance.
The builder writes it from named *parts* over the project's parameter table, checks it as the engine will, and says
what the engine will do with it.  It is ``pyazr.priors.build`` and, on the web page, **Prior on the parameters…**
in the Fit tab.  It has a page of its own beside
:doc:`covariance_builder` (the covariance of a data set's points): the two files hold different things (data rows
against fit parameters) and are named in different places (a data set's field against a ``<config>`` key).

What the engine does with the file
----------------------------------

- It adds :math:`(\mathbf{p}-\boldsymbol{\mu})^T C^{-1}(\mathbf{p}-\boldsymbol{\mu})` to the objective, over the
  listed parameters as the fit vector holds them (the ``param.sav`` space: level energies in MeV, reduced-width
  amplitudes of the run's parametrization, normalizations and the other nuisance parameters in their own units).
- A listed parameter the run fixes is **marginalized**: its row and column are dropped, with a Note.
- A listed parameter with a prior **of its own** (a normalization's error, a shift's error, a nuisance's sigma) has
  that prior **replaced**: its one-parameter term is not added, and the run's Note says so.  Measured on the p+p
  model of ``tests/norm_groups``: with a 10 % prior of its own on a normalization and a 5 % prior on it in the file,
  the objective is the data :math:`\chi^2` plus the file's term only (``tests/pyazr/prior_builder_test.py``).
- The file's ``#parametrization`` tag must be the run's (0 standard, 1 Brune, 2 Park): amplitudes are not
  converted.  The builder writes the run's tag.
- The engine refuses a name the model does not have, a matrix that is not symmetric (to :math:`10^{-10}` of
  :math:`\sqrt{C_{ii}C_{jj}}`) or not positive definite (a Cholesky pivot not above :math:`10^{-13}` of its
  diagonal entry), and a file whose counts disagree.

The parts
---------

===============  ===========================================================================================
kind             what it gives
===============  ===========================================================================================
``fit``          an earlier fit: the centres of its ``param.sav`` and the covariance of its
                 ``covariance.dat`` (its free level energies and widths; written by a fit with
                 ``--covariance-band``) for the ``names`` chosen (all its columns when none): the numbers
                 ``prior_from_fit`` writes.  Its tag must be the run's.
``width``        a standard deviation per parameter, uncorrelated: ``sigma`` in the parameter's unit or
                 ``percent`` of the centre; the centre each parameter's start value, or ``centre`` given.
``correlation``  correlation coefficients between parameters of the other parts (``pairs`` of
                 ``[name, name, rho]``, or ``names`` and a ``matrix``); the sigmas stay.  A pair inside one
                 ``fit`` or ``chain`` part is refused (that part has its correlation already).
``chain``        the sample mean and covariance of a sampler's chain after its burn-in: ``chain_<id>.npy`` of
                 a sampler job of the web page, with ``chain_<id>.json`` beside it (the names and the burn-in;
                 ``discard`` overrides it).  A chain does not record its basis: it is taken as the run's.
===============  ===========================================================================================

Every parameter is in one part.  The names, kinds, units, fixed marks and limits come from the session's
parameter table (``Session.parameters``), not from a list of the builder's own.

**Checks.**  Each is a refusal (``pyazr.ParameterError`` with what was wrong, what is allowed and how to fix it) or
a note:

- refused: a name the parameter table does not have (the close names are said), a parameter in two parts, a part
  that cannot be made (an unknown kind or field, ``sigma`` and ``percent`` both or neither, a percentage of a
  centre 0, a correlation of magnitude 1 or of a parameter in no part, a pair twice, a fit's name its covariance
  lacks, a chain without its JSON or with fewer draws than parameters), a fit of another parametrization;
- refused when written (``refusal`` is set, the record still shown): a matrix that is not symmetric or not
  positive definite, with the smallest eigenvalue, the condition number and the parts concerned;
- notes: a parameter the run fixes (the engine marginalizes it), a parameter with a prior of its own (the engine
  replaces it), a centre outside the parameter's limits, a fit without a tag or with a scaled covariance, a chain's
  basis.

The parts are written on the file's ``#pyazr.prior`` comment line, so the builder reopens with them
(``read_parts``); a file without the line was written by hand or by ``prior_from_fit``.

In Python
---------

.. code-block:: python

   from pyazr import Project, Session, priors
   p = Project.open("project.azr")
   with Session(p.path) as s:
       prior = priors.build(s.parameters, [
           {"kind": "width", "names": ["xs_40_s08.norm", "xs_90_s09.norm"], "percent": 5},
           {"kind": "width", "names": ["width_1_1"], "sigma": 0.05, "centre": 0.6},
           {"kind": "correlation", "pairs": [["xs_40_s08.norm", "width_1_1", 0.4]]},
       ], directory=p.directory)
   print(prior)                    # names, centres, sigmas, the correlation matrix, the checks, the notes
   prior.use(p, "prior.txt")       # writes prior.txt beside the project and sets priorCovariance=prior.txt
   p.save()

``build(parameters, parts, directory=".")`` returns a ``Prior``: ``names``, ``mean`` (the centres), ``cov``,
``sigma``, ``correlation``, ``units``, ``kinds``, ``fixed``, ``bounds``, ``basis``, ``parts``, ``part_of``,
``checks`` (``symmetric``, ``positive_definite``, ``smallest_eigenvalue``, ``largest_eigenvalue``,
``condition_number``), ``refusal`` (what, allowed, fix) with ``refusal_parts``, ``notes``; ``text()``,
``write(path)``, ``use(project, path)`` and ``summary()`` (plain data).  ``PART_KINDS`` is the table of kinds and
their fields.  ``pyazr/examples/prior_builder.py`` runs the case below on a project.

The file's grammar has one definition, the engine's reader; ``read_prior`` reads as it does and refuses what it
refuses (a tag other than 0, 1 or 2 or without a value, a name twice, ``nan``, ``inf``, ``1_0``, ``1d0``, a
non-ASCII space between tokens, a matrix not symmetric or not positive definite), so a file pyazr reads is a file
the engine reads.  What only a run can check -- a name the model lacks, the tag of another basis -- the builder
checks against the parameter table.

A worked case
-------------

The p+p model of ``tests/norm_groups`` (each data set its own free normalization, the first with a 10 % prior of its
own): 5 % on both normalizations around 1, 0.05 MeV\ :sup:`1/2` on ``width_1_1`` around 0.6, correlations 0.4
(first normalization, width) and -0.2 (the two normalizations).  The condition number is 2.618; the note says the
first normalization's own prior is replaced.  At the normalizations 1.07 and 0.96 and the width 0.63 the engine's
term is 2.25, :math:`d^T C^{-1} d` formed by hand (to :math:`10^{-10}`), and the objective is the data
:math:`\chi^2` plus that term, without the 10 % prior (``tests/pyazr/prior_builder_test.py``).

On the page
-----------

**Prior on the parameters…**, above the Fit tab's parameter table, opens the builder (``#/prior``; one click in the
first view as in the expert's).  It shows:

- the file the project names (``priorCovariance``), whether the builder wrote it (it then opens with its parts) or
  it was written by hand, and **Remove the prior** (the key taken out of ``<config>``; the file stays; Undo puts the
  key back);
- the parameter table with a check box per row (name, kind, unit, x, free or fixed, limits, a prior of its own; a
  filter and "free only");
- the parts as rows: a kind chosen from the list and **+ part** (a part ``width`` takes the parameters checked),
  each field as its kind states it: a file chosen among the project folder's files (``param.sav``,
  ``covariance.dat``, ``chain_<id>.npy``), the checked parameters taken into a names field ("take the checked"), a
  number box with its unit, a pair of parameters chosen from two lists with its coefficient;
- **Build**: the engine's refusal (what, allowed, fix, the parts concerned) or the check (eigenvalues, condition
  number), the notes, the centres and sigmas, the correlation matrix as a heat map (a table without plotly.js), the
  file's text;
- **Write and use** with the proposed file name (``prior_parameters.txt``): the file and the key as one edit; Undo
  removes the file (or puts back the one it replaced) and the key together.  Writing again under the same name is
  an edit too.

The page computes nothing of the prior: the routes call ``pyazr.priors.build`` on the session's parameter table.
