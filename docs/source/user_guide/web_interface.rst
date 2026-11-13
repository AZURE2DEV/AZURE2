The web interface
=================

The web interface of AZURE2 is a page in your browser, served by a small server on your machine that holds one
project.  It replaces the Qt setup utility of earlier versions.  The page shows the project as the engine reads it
and sends every edit to the engine, which validates it first; nothing in the page decides whether a value is valid,
and every number it shows (the χ², the model values) is the engine's.  Its API, which scripts and AI agents use
directly, is described in :doc:`web_ui_api`.

Starting it
-----------

.. code-block:: bash

   pip install pyazr[web]            # fastapi and uvicorn
   python -m pyazr.ui project.azr    # prints http://127.0.0.1:<port>/?token=<token>

Open the printed URL in a browser.  The token in it is exchanged for a cookie on the first load, and the address
bar then shows the page without it; keep the URL to open the page again (it is also written to
``.<name>.ui/url`` beside the project, readable by you only).  Stop the server with Ctrl-C.  Edits are kept in
``.<name>.ui.azr`` beside the project until you press **Save** (Ctrl+S): only Save writes your file.

**On a remote machine** (a cluster's login node, a workstation reached over SSH) the server binds to the loopback
address only.  Start it there with a fixed port and forward that port:

.. code-block:: bash

   ssh -L 28400:127.0.0.1:28400 user@remote     # on your laptop
   python -m pyazr.ui project.azr --port 28400  # on the remote machine, in that session

and open the printed URL (``http://127.0.0.1:28400/?token=...``) in the browser of your laptop.  Every other user of
the remote machine can reach its loopback ports: the token is what keeps them out of your project, so do not share
the URL.

What the page shows
-------------------

The header shows the project's name, whether it has unsaved edits, the revision, **Undo**, **Redo**, **Save**, and
the total χ² of the current evaluation with its number of points.  The line at the bottom shows the engine's last
message.  The keys are listed by ``?``: ``1`` … ``7`` switch tabs, Ctrl+Z and Ctrl+Shift+Z undo and redo,
``/`` filters the current table, the arrow keys and Enter move through a table, Ctrl+Enter applies a form, Esc
reverts it, ``t`` switches between the system's theme, light and dark.  Each view has its own address
(``#/data/s02_diff090``), so a link or a reload brings you back to it.

**Project.**  The file and its working copy; the problems the engine's validator reports, each with a link to the
record it names; the ``<config>`` entries, edited as text (an empty value removes a key); the session's options
(Brune, Park, the Coulomb table and the others) with the values in use and the program's flags, read-only for now;
the history of edits; the project's text, read-only, with **Replace text** for an edit of the whole; the engine's
version, build and data files.

**Nucleus.**  The particle pairs (radius, masses, separation energy and excitations can be edited; choose a row);
the levels; for the chosen level its energy, its fixed and active flags and its channels, whose widths are edited in
place (Enter sends a width; the status line then shows the χ² before and after, and Undo takes the edit back).  A
width is written as the file has it (Γ in eV, or an ANC in fm^-1/2 for a bound channel).  The γ-decay and potential
records are edited with the same forms as the records of the Data tab.

**Data.**  The records of every section of the project file, one group per section (data sets, targets,
extrapolations, model data); the list shows the program's χ² per data set, its points and normalization.  Choosing
a record opens its form, generated from the engine's schema, with the engine's description of the record beside it
and the engine's messages under the fields they concern; **Apply** (Ctrl+Enter) sends only the fields you changed.
A data set also shows the first rows of its data file and a preview plot of the data (points with errors) and the
model at the working parameters, against energy or angle; **Export data (CSV)** saves the numbers drawn with their
units.  **+ new**, **Duplicate**, **Rename** and **Delete** add, copy, rename (with every reference) and remove
records; each is one entry of the history.

**Fit.**  On the left the job form, drawn from the server's job table: the kind (fit, MINOS), the method (iminuit's
MIGRAD with HESSE, scipy's L-BFGS-B; the engine's fitters, the sampler, bands and the rate are listed greyed with the
pyazr call they wait for), the method's settings with their defaults, the start point (the working parameters, the
values of the project text, or the best point of a job) and **Run**.  A job runs in its own process on the text of
the current revision with one thread of the server's budget; while it runs, evaluations wait (the header shows
χ² —).  On the right the job: its calls, the objective (χ² + priors) at the last call and the best so far, the
time, a plot of the objective against the call (log scale), the last lines of its log, **Stop** (it ends before its
next objective call and keeps its best point) and **Kill** (at once; the best point reported so far is kept).  A
reload finds a running job again.  When it has ended: the recipe's own words (converged or not), the objective at
the start and at the best point, the χ² and N at the best point, **Accept into project** (the best point written
into the project text as one history entry; Undo takes it back; the χ² of the project as written then stands beside
the fit's, and a difference is said, not hidden), and after an iminuit fit **minos (iminuit)…**, which fills the
form for a MINOS job from that fit.  Below, the table of every parameter: name, kind (its description on hover),
free or fixed, the unit of x, the working x (a free parameter's value is edited here: Enter sets it, the header's χ²
follows; **Reset x** goes back to the text's values), the physical value, bounds, prior, and for the job shown its
best point, errors and MINOS intervals; **compare with start** adds the job's start point.  Filter by kind (the
kinds the engine lists), by "free only" or by text (``/``).  A fixed value is a value of the project text: it is
edited in its record (a link where the parameter's name begins with a record's) or level.  After a fit with HESSE,
the covariance and the correlation as tables and CSV.

**Results.**  Everything at the working parameters (the vector of the last edit or Undo), from one answer of the
server (``GET /api/evaluation``).  Tick the data sets and extrapolations to show (grouped as on the Data tab, with
**all** and **none**); choose x (c.m. energy, excitation energy, c.m. angle; **log x**) and y (as measured, or the
S factor where the data set is a cross section with a charged entrance; **log y**).  Each series gets a plot of its
own (**together** draws them in one plot per unit), with the data times the normalization in use and their errors,
the model, and under it, on the same x axis, the engine's residuals (m − n·d)/(n·σ), whose squares are the χ² of the
fit.  The model of a data set is drawn at its points (a resonance narrower than their spacing is not drawn); an
extrapolation is drawn on its grid, so an extrapolation with an ``adaptive`` grid over a data set's range shows the
model between its points.  **levels** marks the levels of the project on the excitation-energy axis.  Below:
angular distributions at a chosen energy, the χ² of every data set as the program's ``chiSquared.out`` splits it
(with the sum of the squares of the residuals drawn), the model data against their values, and **Write the
program's output files**, which writes ``AZUREOut_*``, ``chiSquared.out`` and the rest into the project's output
directory and lists them with their first rows.  Each plot exports its numbers as CSV (17 significant digits) and
the figure as SVG or PNG.  The view is also a link: ``#/results?x=e_x&y=s_factor&logy=1&levels=1&show=a,b``.

Not there yet
-------------

The tabs **Physics** and **THM** say which step brings them; the ``<thm>`` section is edited through the API. In the
Fit tab the engine's own fitters (MIGRAD, LM, their MINOS), the emcee sampler and bands wait for pyazr calls, and a
fixed parameter is not edited from the table. Adding or removing a particle pair, generating channels from a J^π,
the session's options as editable keys, the file browser, EXFOR and ENSDF imports, the Wigner limits and physical
widths of the channels, the lab energy as an x axis, the model on a dense grid between data points, the reaction
rate and uncertainty bands, and live updates pushed by the server (the page asks every three seconds whether another
client changed the project) come later. The page was checked in a headless Chrome; plotly.js 4.1.2 (MIT licence) is
part of the package, so the page loads nothing from the network, and without it the plots appear as tables of their
numbers.