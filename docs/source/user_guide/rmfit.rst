rmfit — Global Fitting Campaigns
================================

``rmfit`` (directory ``rmfit/`` in this repository) runs a global R-matrix fit
with the in-process engine (:doc:`pyazr`) as a series of resumable, checked
steps. A script executes the steps on a cluster node or a workstation. A
person, or an LLM agent such as Claude Code, supervises and resumes them from
plain-text state files. It is built for models too large to fit by hand:
hundreds of free parameters, tens of thousands of data points.

What it does:

- **sharded evaluation**: the data sets are split across single-thread engine
  sessions. The chi-squared adds up over data sets, so one objective or
  Jacobian of a large model costs seconds. A data set the engine folds (it
  names a target, a beam or a detector) or one that shares an energy scale is
  never split;
- **bounded polishes** in reduced-width-amplitude space with the engine's
  residual rows and their Jacobian (:doc:`pyazr`, ``Objective.residuals`` and
  ``Objective.jacobian``), inside Wigner-limit bounds and energy windows;
- **interference-sign rounds** and a population stage (sign-perturbed restarts
  and crossover between them) to move between basins;
- **structure searches**: level add, remove and relocate, each with a stated
  hypothesis, and a residual-driven level search;
- a declared objective, ``spec.json``. Data treatment (normalizations,
  excluded rows, which sets are fitted) is the evaluator's. A question the spec
  does not answer goes into a decision queue (``DECISIONS.md``) and the run
  waits for it;
- pre-fit audits, a post-fit review that judges the fit point by point
  (``rmfit review``), and posterior sampling with automatic convergence gates
  as the last stage (``rmfit uq``);
- a ledger (SQLite) of every candidate, a ``STATUS.md`` a new session can read
  to continue, and an engine pin: a campaign records the AZURE2 build it
  started with and refuses another.

Two backends
------------

The same ``rmfit`` package runs on two versions of pyazr. ``rmfit/backend.py``
reads the version of the pyazr that ``AZURE2_ROOT`` provides, once, and a
process runs on that one backend:

- **pyazr 3** (this branch): every project is in format 3; an older file is
  converted first (``AZURE2 --convert OLD NEW``, :doc:`older_projects`). rmfit
  reads and edits a project only through ``pyazr.Project``
  (``rmfit/azr3.py``), the data rows and the priors are the engine's
  (``Session.objective()``), a fit is written with ``Session.save_fit``, and
  ``rmfit uq`` samples with ``rmfit/sampling.py``;
- **pyazr 2** (the ``dev`` branch and its ``.azr`` text files): the code paths
  of ``dev``'s rmfit. ``rmfit/mixture.py`` runs on pyazr 2 only; on pyazr 3
  several compound nuclei in one project take its place.

The search, the ledger, the driver, the review criteria and the parameter keys
are one code for both. The design note is ``docs/upgrade/rmfit_two_backends.md``.

Requirements and setup
----------------------

rmfit needs a built ``pyazr`` (the ``_azure2`` module; see :doc:`pyazr`),
NumPy and SciPy. The review figures, the level search and the sampler also use
Matplotlib. rmfit is not part of the pip package: run it from the checkout,
with the directory that contains ``rmfit/`` on ``PYTHONPATH``. It imports
``pyazr`` from ``AZURE2_ROOT``: the ``AZURE2_ROOT`` environment variable, else
the AZURE2 checkout ``rmfit/`` sits in, else ``~/AZURE2``.

.. code-block:: bash

   export PYTHONPATH=/path/to/AZURE2
   python3 -m rmfit.cli campaign_dir init --base model.azr --nshards 8
   python3 -m rmfit.cli campaign_dir status

``init`` checks the starting model (one full evaluation), records the engine
and writes a job script for the campaign. Steps that need the engine run as
jobs, not on a login node.

The engine pin is the MD5 sum of the compiled module in
``$AZURE2_ROOT/pyazr/`` (``python3 -m rmfit.cli campaign_dir engine`` shows
the pinned and the imported engine), and every job script checks that file
before Python starts. A build that writes the module elsewhere
(``-DPYAZR_OUTPUT_DIR``, found through ``PYAZR_ENGINE_DIR``) gives a campaign
no pin; build with the module in ``pyazr/`` for a pinned campaign.

Job scripts and site settings
-----------------------------

Settings that belong to a site or a user, not to a campaign, live in
``~/.rmfit/site.json`` (or the file named by ``$RMFIT_SITE``). An environment
variable overrides the file. The scheduler is one of:

- ``"sge"`` (the default): csh scripts, with the queue from ``sge_queue`` and
  the parallel environment from ``sge_pe``;
- ``"slurm"``: bash scripts with ``#SBATCH`` lines from ``slurm_partition``,
  ``slurm_account``, ``slurm_time`` and ``slurm_extra``;
- ``"shell"``: bash scripts with no scheduler, for a workstation or a node
  you already hold. They log to a file, refuse a second run on the same
  campaign, and run the drive's legs one after another.

``setup`` lists the lines every script runs before Python (modules, a Python
with NumPy and SciPy), written in the script's shell. ``mail`` is the job mail
address. For example:

.. code-block:: json

   {"scheduler": "slurm", "mail": "you@example.org", "slurm_partition": "compute",
    "slurm_time": "2-00:00:00", "setup": ["module load python/3.12 gsl"]}

The SGE templates have run on a production cluster. The Slurm and plain-shell
templates are checked by the tests (bash, a simulated job environment); the
Slurm ones have not yet run under a real Slurm scheduler.

Tests
-----

Each script of ``rmfit/tests/`` is a ctest, ``rmfit_<name>`` with the label
``rmfit``. A script that is the same on both backends is ``test_<name>.py``;
where the two differ, ``test_<name>.py`` is the pyazr 2 version and
``test_p3_<name>.py`` the pyazr 3 one, and each is skipped on the other
backend. The engine-free scripts take seconds and run with
``ctest -LE slow``; the scripts that start engine sessions on the
``tests/13N`` project take minutes to tens of minutes and carry the label
``slow``:

.. code-block:: bash

   ctest --test-dir build -L rmfit -LE slow     # seconds
   ctest --test-dir build -R '^rmfit_' -L slow  # the engine scripts
   python3 rmfit/tests/test_site_mail.py        # one script by hand

Further reading
---------------

- ``rmfit/README.md``: the module table and the command line;
- ``rmfit/DESIGN.md``: principles, the closed loop (``rmfit drive``), UQ,
  benchmarks and open items;
- the Claude Code skill ``rmfit-campaign`` (``rmfit/.claude/skills/rmfit-campaign/``,
  which lives with rmfit): the process guide (how to start, drive, judge and resume a campaign).
