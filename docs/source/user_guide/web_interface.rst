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
the URL.  The forward must end at the loopback address **of the machine the server runs on**
(``127.0.0.1:28400`` as seen from that machine): a forward that reaches the machine by its network name from a
gateway (``28400:remote:28400`` in a session on the gateway) finds nothing listening.  Through a gateway, either
log in to the remote machine with the gateway as jump host (``ssh -J user@gateway -L
28400:127.0.0.1:28400 user@remote``) or forward in two steps.

**If the browser shows "the token is needed"**, the address was opened without its ``?token=...`` part (a bookmark
of the page, a browser restarted since) or with the token of another server: paste the token into the field of that
page, or open the printed URL again.  Several projects served on one machine keep separate cookies (one per port).

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

**Fit.**  On the left the job form, drawn from the server's job table: the kind (fit, MINOS), the method, the
method's settings with their defaults, the start point (the working parameters, the values of the project text, or
the best point of a job) and **Run**.  The methods are first the engine's own fitters, one per row of the engine's
fitter table (``pyazr.kinds()["fitters"]``, the fits of the command-line program): ``engine_migrad`` (the program's
default fit and the form's), ``engine_gsl-lm``, ``engine_lm``, ``engine_migrad-gradient``; under the method the
engine's description of the fitter and what it gives (the program's switch, a covariance, whether MINOS can follow);
settings ``max_calls`` (the program's ``--max-calls``) and, where the fitter allows it, ``minos``.  Then the recipes
on the objective: iminuit's MIGRAD with HESSE and scipy's L-BFGS-B; the sampler (emcee); the bands; the rate is
listed greyed with the pyazr call it waits for.  A job runs in its own process on the text of
the current revision with one thread of the server's budget; while it runs, evaluations wait (the header shows
χ² —).  On the right the job: its calls, the objective (χ² + priors) at the last call and the best so far, the
time, a plot of the objective against the call (log scale), the last lines of its log, **Stop** (it ends before its
next objective call and keeps its best point) and **Kill** (at once; the best point reported so far is kept).  A
reload finds a running job again.  When it has ended: the recipe's own words (converged or not), the objective at
the start and at the best point, the χ² and N at the best point, **Accept into project** (the best point written
into the project text as one history entry; Undo takes it back; the χ² of the project as written then stands beside
the fit's, and a difference is said, not hidden), and after a fit that MINOS can follow (the engine's MIGRAD, with or
without the gradient, or iminuit's) **minos (…)…**, which fills the form for a MINOS job from that fit's best point.
An engine fit ends where the same fit of ``Session.fit_with_engine`` and of the program ends; its result says the
engine's status in its own words (``converged``, ``call_limit``, ``stopped``, ...), its errors are those the program
writes into ``param.sav``, its covariance HESSE's after MIGRAD and (JᵀJ + priors)⁻¹ after LM and GSL-LM.  While it
runs, the panel shows the stage (fit, then HESSE and MINOS, whose calls are counted on) and, for LM, its damping;
**Stop** ends it at the best point it has seen.  Where LM or GSL-LM cannot run (a data set without an analytic
Jacobian) the engine runs MIGRAD in its place: the result shows the engine's notice and the fitter that ran.
Below, the table of every parameter: name, kind (its description on hover), free or fixed, the unit of x,
the working x (a free parameter's value is edited here: Enter sets it, the header's χ²
follows; **Reset x** goes back to the text's values), the physical value, bounds, prior, and for the job shown its
best point, errors and MINOS intervals; **compare with start** adds the job's start point.  Filter by kind (the
kinds the engine lists), by "free only" or by text (``/``).  A fixed value is a value of the project text: it is
edited in its record (a link where the parameter's name begins with a record's) or level.  After a fit with HESSE,
the covariance and the correlation as tables and CSV.

After a fit, **sample (emcee)…** fills the form for the sampler from the fit's best point (walkers, stored steps,
the start walkers' relative spread, thin, seed, the burn-in left out of the summary).  It runs emcee on the
log-probability −χ²/2 of the objective (zero outside the fit's limits; a start point outside them is refused, naming
the parameter) in one process; the panel shows the step, the acceptance and a running autocorrelation time τ, and
**Stop** ends it between steps.  The chain is written to ``output/chain_<id>.npy`` (stored steps × walkers ×
parameters) with a JSON summary beside it; the page shows, per parameter, the mean, standard deviation, the 2.5,
16, 50, 84 and 97.5 % quantiles, τ, the effective draws and the split R̂, and the draws' covariance.  A sampler
started from a sampler job's best point continues its chain, also a sampler that was killed or interrupted by a
restart of the server (from its chain as last saved, every 50 stored steps): **sample (emcee)…** is offered beside
it.  A chain shorter than about 50 τ is not quoted
(:doc:`mcmc`).  **bands (covariance)…** after a fit with a covariance computes the 1σ band of the model of every data
set and extrapolation at the fit's best point (linear propagation over the level parameters, ``pyazr.bands``);
**bands (draws)…** after a sampler evaluates the model at draws of its chain and takes the 16th, 50th and 84th
percentiles at each point.  The result says which parameters the band varied (the level parameters for the linear
band, every free parameter for the band from draws) and which were held at the start point.

**Results.**  Everything at the working parameters (the vector of the last edit or Undo), from one answer of the
server (``GET /api/evaluation``).  Tick the data sets and extrapolations to show (grouped as on the Data tab, with
**all** and **none**); choose x (c.m. energy, excitation energy, c.m. angle; **log x**) and y (as measured, or the
S factor where the data set is a cross section with a charged entrance; **log y**).  Each series gets a plot of its
own (**together** draws them in one plot per unit), with the data times the normalization in use and their errors,
the model, and under it, on the same x axis, the engine's residuals (m − n·d)/(n·σ), whose squares are the χ² of the
fit.  The model of a data set is drawn at its points (a resonance narrower than their spacing is not drawn); an
extrapolation is drawn on its grid, so an extrapolation with an ``adaptive`` grid over a data set's range shows the
model between its points.  **levels** marks the levels of the project on the excitation-energy axis.  **band**
shades the band of the newest band job of the Fit tab around each series; the head says which job and project text
it belongs to and which parameters it varied, and its legend names the job and the point it is drawn around (the
model at that job's start point: accept the point to see the model line inside it).  When the project's text is a
job's accepted best point and the band drawn is around another point, the head offers the band at the accepted
point.  The CSV holds the band's edges.  Below:
angular distributions at a chosen energy, the χ² of every data set as the program's ``chiSquared.out`` splits it
(with the sum of the squares of the residuals drawn), the model data against their values, and **Write the
program's output files**, which writes ``AZUREOut_*``, ``chiSquared.out`` and the rest into the project's output
directory and lists them with their first rows.  Each plot exports its numbers as CSV (17 significant digits) and
the figure as SVG or PNG.  The view is also a link: ``#/results?x=e_x&y=s_factor&logy=1&levels=1&show=a,b``.

Physics
-------

The **Physics** tab evaluates the engine's own functions (kinematics, Coulomb functions, channel functions, widths,
resonance quantities, Coulomb scattering, the THM vertex; the stopping powers, ranges, energy losses and straggling
of the stopping library; and under "this project" the session's channel functions of a pair) on a range, from the
same tables as ``pyazr.physics`` and ``pyazr.stopping``.  A stopping plot shows SRIM's citation under it, and its
CSV and Python text carry it: every use of these stopping powers cites SRIM.  Pick a function, choose the argument to plot against,
edit the form (several values of an argument, separated by blanks, give one line each) and press **Plot**; every
function plots at once with typical values.  The status line says how many elements were computed and why the others
were not (they are gaps in the plot, never zeros); the CSV holds the arrays as returned, one line per element; **Copy
as Python** gives the ``pyazr.physics`` call that makes the same arrays.  A view is a link:
``#/physics/coulomb_waves?from=-2&to=20&spacing=linear&v.l=0 1 2 3``.

THM
---

The **THM** tab is listed when the project has THM content (a ``<thm>`` section, or a data set that names a THM
experiment).  **Model options**: a form of the global keys of ``<thm>`` (vertex, kinematics, entranceL,
coulombIntegral, spectatorEnergy), each with its choices and the engine's default; a key the form does not list
(``spectatorEnergy[2]``, one pair's spectator energy) is set under "other key".  **Experiments**: the table of the
experiments with the data sets that name them and their keys as written; a row opens the form of its keys (beam,
target, spectator, Ebeam, background, ps window, distortion, optical potentials, vertexModel, coherent background,
...) with "data sets that join", and **Remove**; **+ experiment** adds one.  Every change is one edit of the project:
**Undo** takes it back, and what the engine refuses is shown in its words.  **Channels**: the THM pairs with their
binding energy B and their channels, with whether a width is an amplitude; these are read from the session and
cannot be edited here yet.  **Diagnostics**: the session's THM calls at the working parameters: the experiments'
profiled norms and backgrounds (with their σ and covariance, and the coherent background's parameters), the entrance
vertex per channel and level, the Coulomb line shape per level, and the distortion factor R(E), each with its form,
its numbers, its plots (an energy outside the reach of the spectator window is a gap), CSV and Copy as Python, as in
the Physics tab.  A view is a link: ``#/thm/LC09?fn=thm_vertex``.

Not there yet
-------------

In the Fit tab a fixed parameter is not edited from the table.  The THM binding energy and amplitude flag are not
edited from the THM tab.
Adding or removing a particle pair, generating channels from a J^π,
the session's options as editable keys, the file browser, EXFOR and ENSDF imports, the Wigner limits and physical
widths of the channels, the lab energy as an x axis, the model on a dense grid between data points, the reaction
rate, and live updates pushed by the server (the page asks every three seconds whether another
client changed the project) come later. The page was checked in a headless Chrome; plotly.js 4.1.2 (MIT licence) is
part of the package, so the page loads nothing from the network, and without it the plots appear as tables of their
numbers.