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

The first view and the expert's
-------------------------------

The page opens in its **first view**: the tabs follow the old setup program's path and words, **Project**,
**Particle Pairs**, **Levels and Channels**, **Segments** (Segments From Data, Segments Without Data, Model
data), **Experimental Effects** (targets, beams, detectors), **Fit** and **Plot**.  Each screen holds what its task
needs: a field the engine's table marks as rarely set (the source, EXFOR and DOI of a data set, its pair keys,
THM experiment, reference cross section, UPOS, a target's phase, Doppler broadening or screening, ...) waits under
**More fields** unless the record sets it, and every field's meaning is written under its title.  The Project tab
says what to do next (add the pairs, then the levels, then the data, then Calculate).  The γ decays and potentials,
the session's options, the history and the project's text are one click away under **More**.

Tick **expert** (top right; the key ``e``) for the expert's view: the tabs Project, Nucleus, Data, Fit, Results,
Physics and THM with every field, as before.  The choice is kept in this browser; ``?depth=expert`` or
``?depth=first`` in an address opens one view in that depth.  The two depths show the same project: every
address means the same thing in both.

**Calculate** (in the header on every tab; the key ``c``) asks the engine for the current text and opens the
Plot tab with every data set, its model and its residuals, the χ² per data set and in total, and what changed
since the last Calculate: Δχ² in total and per data set, the data set that moved most first (choose a row to go
to its plot), and the plots in that order.  The status line also gives the change of χ² since the last Calculate
after every edit.  The page computes no χ²: each number is the engine's, and a change is the difference of two.

Six tasks, step by step
~~~~~~~~~~~~~~~~~~~~~~~

*Add a particle pair.*  **Particle Pairs** › **+ pair**: the kind (particle, photon, beta), each partner's J,
parity, Z, mass (u) and g, the separation energy, the radius and the rest › **Add pair**.  The engine checks every
number and says what it refuses.

*Add a level and set a width free.*  **Levels and Channels** › **+ level**: J, parity, E (MeV), l max ›
**Channels** (the channels that J^π has with the project's pairs) › tick those you want › **Add level**.  Choose
the level, untick a channel's **fixed** and type its width (Enter sends it).

*Add a data set with its file and a normalisation.*  **Segments** › **+ new**: the reaction, the quantity, the
unit, the data file (the first rows of the file are shown once it is set), the columns, the normalisation (a value
of its own with a mode and a prior, or a shared one) › **Apply**.

*Add a target effect.*  **Experimental Effects** › **+ new**: a target (its kind, the active nucleus, the
composition, the areal density, stopping and straggling), a beam or a detector › **Apply**; then **Segments** ›
the data set › **Target** (or Beam, Detector) › **Apply**.

*Run a fit and look at it.*  **Fit** › the method › **Run**; follow its progress; **Accept into project** ›
**Calculate**: every data set against the new model, with what the fit changed.

*Give a data set a covariance.*  **Segments** › the data set › **Covariance file**: a file of the project folder
(``data_formats.rst``, covariance files) › **Apply**.

What the page shows
-------------------

The header carries the original program's icon (the red and green spheres with the blue "(A)"; it is also the
browser tab's icon) and shows the project's name, whether it has unsaved edits, the revision, **Undo**, **Redo**,
**Save**, and the total χ² of the current evaluation with its number of points.  **About** gives the version of the
program and of pyazr, the reference every publication that uses AZURE2's results cites (R.E. Azuma et al.,
*Physical Review C* **81**, 045805 (2010)) and the link to this documentation.  The page's colours are the icon's:
its blue is the accent, its red marks an error and its green a success, and each of these states also carries a
symbol (✕, ⚠, ✓), so that nobody has to tell red from green.  The line at the bottom shows the engine's last
message.  The keys are listed by ``?``: ``1`` … ``7`` switch tabs, Ctrl+Z and Ctrl+Shift+Z undo and redo,
``/`` filters the current table, the arrow keys and Enter move through a table, Ctrl+Enter applies a form, Esc
reverts it, ``t`` switches between the system's theme, light and dark.  Each view has its own address
(``#/data/s02_diff090``), so a link or a reload brings you back to it.

**Project.**  The file and its working copy; the problems the engine's validator reports, each with a link to the
record it names; the **Run settings**: every ``<config>`` entry the engine knows, in groups (the old program's
Runtime Options first: the formalism, Brune, Park, the parameter transformation, external widths, RMC, GSL Coulomb
functions; then physics, fit, directories, checks, numerics), each with its control (a switch, a choice, a number
box), its key, what it does, what the engine does when the project does not set it, and whether the project sets
it; a change is applied at once and Undo takes it back; **×** removes a key; a search box finds a setting by any
word; the rarely set ones wait under **More settings**, and **Advanced** has the entries as text (any key); the
session's options (Brune, Park, the Coulomb table and the others) with the values in use and the program's flags;
the history of edits; the project's text, read-only, with **Replace text** for an edit of the whole; the engine's
version, build and data files.

**Nucleus.**  The particle pairs: choose a row to edit its radius, masses, separation energy, excitations and THM
binding energy B (the fields are those the server's pair edit takes), or to remove it with **Remove pair** (refused
while a level has a channel of the pair, in the engine's words; then **Remove pair … and its channels** removes
those channel rows too, and with them their fit parameters).  **+ pair** adds one: its kind (particle, photon,
beta), its partners (J, parity, Z, mass in u, g), separation energy, radius and the rest; a pair no level names yet
is kept in the file by a row of its own.  The levels; **+ level** adds one: J (0.5, 1, ...), parity, E and whether
its energy is fixed, then **Channels** lists the channels that J^π has with the project's pairs (the engine's rule,
up to ``l max``), each added with width 0, fixed (with every one checked, a J^π the project has takes its group's
channels).  For the chosen level its energy, its fixed and active flags and its channels, whose widths are edited in
place (Enter sends a width; the status line then shows the χ² before and after, and Undo takes the edit back), with
whether each width is an amplitude.  A width is written as the file has it (Γ in eV, an ANC in fm^-1/2 for a bound
channel, or an amplitude in MeV^1/2 where the row says so; ticking the amplitude box does not convert the number).  The
γ-decay and potential
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
engine's description of the fitter and what it gives (the program's switch, a covariance, whether MINOS can follow,
whether it honours a call limit, the fitter it falls back to);
settings ``max_calls`` (the program's ``--max-calls``) and, where the fitter allows it, ``minos``.  Then the recipes
on the objective: iminuit's MIGRAD with HESSE and scipy's L-BFGS-B; the sampler (emcee); the bands; the reaction rate (drawn by the
Rate tab).  A job runs in its own process on the text of
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
Jacobian) the engine runs MIGRAD in its place: the result shows the engine's notice and the fitter that ran.  The
result shows the fit's outcome in the words every fitter shares: ``converged`` in a green box, anything else
(``call_limit``, ...) in a warning box saying the best point is not a minimum; and the engine's notices (a call
limit the fitter ignored, a call limit reached).  After the engine's MINOS the parameter table's column **MINOS
validity** says for each side whether the crossing was found or why not (``at_limit``, ``call_limit``,
``new_minimum``).
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

Rate and covariance
-------------------

The **Rate** tab computes the reaction rate of a pair of pairs as the program's ``rate`` command does, at the working
x or at a fit's or a sampler's best point, with the band of that fit's covariance or that sampler's draws, the Gamow
window per temperature, a comparison table and the table as CSV (:doc:`reaction_rate`).  **Covariance…** on a data set
of the Data tab opens the covariance builder (``#/covariance/<name>``): the data set's covariance from named parts,
checked as the engine will read it, written and named on the data set as one edit (:doc:`covariance_builder`).

ENDF evaluations
----------------

The **ENDF** tab (in the expert's view; ``#/endf``) puts an evaluated file beside the project, with
``pyazr.endf`` (:doc:`endf`). It is named after the format because that is the word an evaluator searches for, and
both of its screens start from an ENDF-6 file.

- **Over the data** (``#/endf?files=a.endf,lib/b.endf``): the ENDF-6 files of the project folder (and two levels
  below) with a box each; a path inside the project's root can be typed. Each file shown has a colour, and a line
  that opens what it holds (MAT, target and projectile, library, lab, date, authors, LRP, EMAX, the sections in the
  manual's words and the header text). For every data set, the files' values at its points in its quantity and unit
  are drawn over the data (times its normalization) and the model at the working x, with the ratios of each file
  to the data and to the model below; where a file has nothing for a data set (another reaction, a quantity an ENDF
  file does not hold), its refusal says what, what is allowed and what to do, in the data set's place. Which files
  are shown is part of the address, not of the project.
- **Check and correct** (``#/endf/check?file=a.endf``): the checking codes pyazr knows, each found or not (found
  through ``AZURE2_ENDF_UTILITIES``, ``AZURE2_PREPRO`` and ``PATH`` of the server), with what it does and what it
  does to numbers. **Check** starts a job: pyazr's own format check (always; each finding with its line and the
  page of ENDF-102), then each code ticked, on a copy, with a time limit; each code's findings are shown apart, with
  the meaning pyazr gives a known message (the known false alarms: LRP=2 in the charged-particle files, PSYCHE's
  "AP NOT IN RANGE" in R-Matrix Limited, CHECKR and LTP=15), the counts its report states and the raw report one
  click away. **Correct** runs a code that rewrites (STANEF, FIXUP, DICTIN, ENDF2C) on a copy, after saying what it
  does to numbers; the candidate's difference per section (unchanged, bookkeeping, form, numbers with the digits
  lost, restructured) is shown, then **Keep as a new file** (beside the original unless named) or **Replace**, which
  asks a second time. Neither changes the project text.

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
...) with "data sets that join", and **Remove**; **+ experiment** adds one.  The forms are the engine's table of the
keys it reads, so a key the engine gains appears without a change of the page; the beam, target and spectator boxes
offer the nuclei the engine knows by name and the optical-potential boxes its global sets (``pyazr.kinds()``).  Every
change is one edit of the project:
**Undo** takes it back, and what the engine refuses is shown in its words.  **Channels**: the THM pairs with their
binding energy B and their channels, with whether a width is an amplitude, as the session reads them; a row opens
the pair or the level in the Nucleus tab, where both are edited.  **Diagnostics**: the session's THM calls at the
working parameters: the experiments'
profiled norms and backgrounds (with their σ and covariance, and the coherent background's parameters), the entrance
vertex per channel and level, the Coulomb line shape per level, and the distortion factor R(E), each with its form,
its numbers, its plots (an energy outside the reach of the spectator window is a gap), CSV and Copy as Python, as in
the Physics tab.  A view is a link: ``#/thm/LC09?fn=thm_vertex``.

Not there yet
-------------

In the Fit tab a fixed parameter is not edited from the table.  The THM binding energy and amplitude flag are
edited in the Nucleus tab, not in the THM tab.  The session's options as editable keys, the file browser, EXFOR and
ENSDF imports, the Wigner limits and physical
widths of the channels, the lab energy as an x axis, the model on a dense grid between data points, and live updates pushed by the server (the page asks every three seconds whether another
client changed the project) come later. The page was checked in a headless Chrome; plotly.js 4.1.2 (MIT licence) is
part of the package, so the page loads nothing from the network, and without it the plots appear as tables of their
numbers.