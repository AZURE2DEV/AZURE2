The web UI's API
================

The web UI of AZURE2 is a server on your machine that holds one project and answers JSON requests; the page in
the browser is one client of it, a script or an AI agent is another.  This page describes the server's API (the page
itself: :doc:`web_interface`; the design is ``docs/upgrade/web_ui.md``).

Starting the server
-------------------

.. code-block:: bash

   pip install pyazr[web]            # fastapi and uvicorn
   python -m pyazr.ui project.azr    # prints http://127.0.0.1:<port>/?token=<token>

Options: ``--port N`` (default: a free port above 20000), ``--threads N`` (the compute budget, default 1: one
computation at a time), ``--token T`` (a fixed token instead of a random one), ``--root DIR`` (the directory the
server may read and write; default the project's).  The server binds to a loopback address only (``--host`` refuses
anything else; for remote use, forward the port with SSH) and writes its URL also to ``.<name>.ui/url`` beside the
project (readable by you only).  Stop it with Ctrl-C: its worker processes end with it.

Every request needs the token: the header ``Authorization: Bearer <token>``, or ``?token=<token>`` in the URL (which
sets a cookie).  Without it the answer is 401 with ``error``, ``message`` and ``fix`` in JSON (a browser that asks for
the page gets a page with a field for the token instead).  ``GET /api/describe`` lists every route with its
documentation and the pyazr call behind it; ``GET /openapi.json`` is the machine-readable form.

What the server keeps
---------------------

- **The project text is the state.**  Every edit goes through ``pyazr.Project`` and is validated by the engine
  before it is accepted; a refused edit answers 422 with the engine's diagnostics and changes nothing.  A refusal's
  body has ``error`` (the class), ``message`` (the whole text) and the message's parts ``what``, ``allowed``, ``fix``
  as pyazr's errors carry them (``null`` where the message is the engine's own text: then ``message`` is all there
  is).  Every
  accepted edit gives a new revision number, and every edit names the revision it was made on (``"rev"`` in the
  body, or an ``If-Match`` header): an edit on an old revision answers 412, so two clients cannot overwrite each
  other.  ``POST /api/project/undo`` and ``/redo`` step through the history of texts.
- **Your file is written only by** ``POST /api/project/save``.  Until then the current text is kept in
  ``.<name>.ui.azr`` beside the project; a server started after a crash offers it back (``crash_copy`` in
  ``GET /api/project``, ``POST /api/project/restore``).
- **The engine runs in worker processes.**  Evaluations use a session of the current text, built again when the
  text changes; the free parameters' working values ``x`` (``GET``/``PATCH /api/x``, in the fit basis) move without
  a rebuild.  The χ² per data set is the program's (the rows of its ``chiSquared.out``, six significant digits).
- **The model as a curve**: ``GET /api/evaluation/curve/<name>`` is ``Session.model_curve`` at the working x: per
  grid point ``e_cm``, ``e_lab``, ``e_x``, ``angle_cm``, ``model`` (the engine's unit, as ``model`` of
  ``/api/evaluation/data/<name>``) and ``at_data`` (one of the data set's points), the ``branches`` (``start``,
  ``stop``, ``fixed``: ``angle_lab``, ``angle_cm`` or ``energy_lab``, ``value``), ``variable`` (``energy`` or
  ``angle``) and ``s_factor`` where offered.  The curves of every data set are made in one session of the
  extrapolations (their folds on one mesh, as the data's), kept until the text or an energy shift changes them.  A
  data set that is not one function of one variable, or whose columns give each point its own window, answers 422
  with ``what``, ``allowed``, ``fix``.
- **Jobs** (``POST /api/jobs``; the table at ``GET /api/jobs/kinds``): the engine's fitters (``"method":
  "engine_migrad"``, one per row of ``pyazr.kinds()["fitters"]``, through ``Session.fit_with_engine``) and their
  MINOS, recipe fits (scipy L-BFGS-B, iminuit MIGRAD + HESSE) and iminuit's MINOS, the emcee sampler
  (``"kind": "sample"``; its chain in ``output/chain_<id>.npy``) and bands (``"kind": "bands"``, ``"method":
  "covariance"`` or ``"draws"``, started with ``"x0": {"job": <a fit's or a sampler's id>}``, whose record the job
  takes as its source), each in a process of their own,
  with progress (``GET /api/jobs/<id>``), a cooperative stop (``/stop``), a kill (``/kill``), a checkpoint of the
  best point, and ``/accept`` to write the best point into the project text.  With a budget of one thread, an
  evaluation asked while a job runs answers 409 with the job's id.
- **Functions on arrays** (the Physics tab's routes; no project needed): ``GET /api/physics`` and ``GET
  /api/stopping`` list the rows of the engine's physics and stopping tables (``pyazr.kinds()["physics"]``,
  ``["stopping"]``: arguments with unit, domain, example and sweep, returned fields with units; every stopping row
  carries SRIM's citation, owed by every use).  ``POST /api/physics/<fn>`` and ``POST /api/stopping/<fn>`` take
  ``{"arguments": {name: value}}``: a value is a number, a nested list of numbers (broadcast with numpy's rules), a
  choice by name or a list of names, a text (a stopping row's ``formula``), or an axis ``{"from", "to", "points",
  "spacing"}`` (``linear``, ``log``, ``integer``) that the server makes with numpy.  The answer: the returned fields
  as lists, their units, a status per element with the words of each status (an element outside the domain is
  ``"nan"``, never a refusal), the arguments as used, ``python`` (the pyazr call that makes the same arrays; run it
  to check) and, for stopping, ``citation``; at most 20000 elements per request.  ``GET /api/session/functions`` and
  ``POST /api/session/functions/<fn>`` do the same for the session's own functions, the rows of
  ``Session.describe()["functions"]`` (``channel_functions`` of a pair of the project, and the THM diagnostics
  ``thm_experiments``, ``thm_background``, ``thm_lineshape``, ``thm_vertex``, ``thm_distortion`` with an
  ``experiment`` argument), in the session worker at the working parameters (the fit vector ``params`` is not an
  argument of the route); an argument with ``"required": false`` may be left out.  An answer adds ``lines`` (a field
  over channels or levels: one list per line, labelled by the fields its row names), ``facts`` (the call's scalars
  with units) and, for the THM profiles, ``tables``.  The session's own answer of each THM call, as plain data, is
  ``POST
  /api/session/thm_<call>`` (``experiment``, ``energies``).
- **Pairs and levels** (the Nucleus tab's routes): ``GET /api/project/pairs``; ``POST /api/project/pairs``
  (``Project.add_pair``: ``kind``, ``light`` and ``heavy`` as ``{J, parity, Z, mass, g}``, ``separation_energy``,
  ``radius``, ...; the answer's ``result`` is the new key), ``PATCH /api/project/pairs/<key>`` (``set_pair``, with
  ``binding_energy`` the THM B), ``DELETE /api/project/pairs/<key>`` (``remove_pair``; refused while a level names
  the pair unless ``"remove_channels": true``); ``GET /api/project/allowed_channels?J=&parity=`` (the channels a J^π
  has with the project's pairs: ``[{pair, l, s}]``; ``l_max`` and the photon limits as ``Project.allowed_channels``
  takes them); ``POST /api/project/levels`` (``channels`` a list or ``"allowed"``); ``PATCH
  /api/project/levels/<label>/channels/<pair>/<l>/<s>`` (``width``, ``fixed``, ``width_is_amplitude``).
- **What a fit returns**: an engine fit job's ``result`` has ``outcome`` (the words every fitter shares:
  ``converged``, ``call_limit``, ...), ``notices`` (the engine's sentences), and after MINOS per parameter ``lower``,
  ``upper``, ``valid`` and ``sides`` (``{lower, upper}``, each ``{valid, at_limit, call_limit, new_minimum}``); the
  rows of ``GET /api/jobs/kinds`` carry the fitter table's columns among their ``facts`` (``honours a call limit``,
  ``falls back to``).
- **THM** (the THM tab's routes): ``GET /api/project/thm`` (``null`` without a ``<thm>`` section), ``PATCH
  /api/project/thm`` (``{"options": {key: value or null}}``), ``PUT /api/project/thm/experiments/<name>``
  (``{"keys": {...}, "datasets": [names that join]}``) and ``DELETE`` (refused while a data set names it); ``GET
  /api/thm/keys`` describes the keys as record kinds of a schema (value type, choices, doc): the engine's table of
  the keys its parser reads, ``pyazr.kinds()["thm"]``.
- **Files and references**: ``GET /api/files`` lists the project directory's files (``path`` relative to the
  project, as a record's field names it; ``size``; ``kind``: ``project``, ``table`` with ``columns``, ``text``,
  ``binary``; not ``output/``, ``checks/`` or hidden files, never outside the root); ``GET
  /api/files/<path>?rows=`` gives a file's first rows.  ``GET /api/project/records/<name>/preview?field=&value=`` is
  what the engine makes of the record with one field set (a data file, say) on a copy of the project: its
  description (``rows_read``, ``rows_kept``, ``kept_E_lab_MeV``, ...) or the refusal; the project is not changed.
  ``GET /api/project/references`` lists every reference of a record to another (``record``, ``field``, ``name``,
  ``epoch``; ``Project.references``).  ``GET /api/kinds`` carries the words of a pair's fields (``pairs``:
  ``fields``, ``partner_fields``, ``arguments``, ``kinds``) and of a data set's description (``description``).
- Routes that need a pyazr call that does not exist yet (the level matrix, ...) answer 501 naming that call.
- ``GET /api/covariance/{name}``, ``POST /api/covariance/build``, ``POST /api/covariance/write`` serve the covariance
  builder (``pyazr.covariance``; the write is one history entry with its file); ``GET /api/rate/grids`` the rate's
  temperature grids; the rate itself is the job ``rate``/``engine`` (``pyazr.rate``).
- **ENDF-6 files** (the ENDF tab, ``pyazr.endf``): ``GET /api/endf/files`` lists the project folder's ENDF-6 files;
  ``POST /api/endf/summary`` (``{"path": ...}``) is ``read(path).summary()``; ``POST /api/endf/curves``
  (``{"files": [...]}``) gives per data set the files' values at its points at the working x
  (``Material.curves(data, session=s)``) beside the data times its normalization and the model, with the ratios
  ``ratio_to_data`` (and its error from the data's) and ``ratio_to_model``, or the refusal's ``what``, ``allowed``,
  ``fix``; ``GET /api/endf/checkers`` lists the codes found and every code pyazr knows with what it does to
  numbers. The check is the job ``endf_check``/``codes`` (``pyazr.endf.check``; settings ``path``, ``codes``,
  ``timeout``), the rewrite the job ``endf_standardize``/``program`` (``pyazr.endf.standardize``; settings ``path``,
  ``program``, ``mf2``, ``timeout``), whose candidate ``POST /api/endf/accept`` (``{"job_id", "to", "overwrite"}``)
  keeps; these jobs build no session. A path outside the root is refused (403) before anything reads it.  The
  curves route also gives, per data set, its model curve (``curve``) and each file's values on its grid
  (``on_curve``).
- **A new project from an ENDF-6 file, and the project as one**: ``POST /api/endf/resonances`` (``{"path",
  "isotope", "energy_range"}``) is what MF=2 holds (``Material.resonances``: ``text``, ``el_eV``, ``eh_eV``,
  ``entrance``, ``pairs``, ``groups``) or ``refused``, nothing written; the job ``endf_import``/``to_project``
  (settings ``path``, ``isotope``, ``energy_range``, ``separation_energy``: ``ame2020`` or MeV, ``name``) runs
  ``pyazr.endf.to_project`` in a directory of its own (``result.project``: ``report``, ``report_text``, ``text``,
  ``name``); ``POST /api/endf/keep_project`` (``{"job_id", "to"}``) writes it as a new ``.azr`` of the folder (never
  over a file; the answer's ``open_path`` is what ``POST /api/project/open`` takes).  ``GET /api/endf/export`` is
  what a file of the current text at the working x would lose (``pyazr.endf.from_session``: ``losses``,
  ``photon_widths``, ``not_written``, or ``refused``), the data's ``energy_range_MeV`` and a ``name`` not in use;
  the job ``endf_write``/``writer`` (settings ``name``, ``mat``, ``authors``, ``lab``, ``date``, ``reference``,
  ``library``, ``e_low``, ``e_high`` in MeV lab, ``tolerance``, ``number_form``) runs ``pyazr.endf.write`` into a
  new file of the folder (``result.written``: ``path``, ``sections``, ``grids``, ``mf2_check``, ``losses``,
  ``not_written``).
- ``GET /api/describe`` carries ``first_view``: the columns of the lists of data sets (``description``) and of
  pairs (``pairs``) the page's first view shows (``pyazr/ui/first_view.py``).
- ``GET /api/prior``, ``POST /api/prior/build`` (``{"parts": [...]}``), ``POST /api/prior/write`` (``{"parts":
  [...], "path": ..., "rev": ...}``) serve the prior builder (``pyazr.priors.build`` on the session's parameter
  table; the write sets ``priorCovariance`` and writes the file as one history entry); removing the prior is ``PATCH
  /api/project/config`` with ``{"entries": {"priorCovariance": null}}``.

A session from Python
---------------------

With the server running on a copy of the test project ``tests/upgrade/T2`` (``url`` and ``token`` as it printed
them), widen the first channel of the first level, compare the χ², and undo:

.. code-block:: python

   import httpx
   api = httpx.Client(base_url=url, headers={"Authorization": f"Bearer {token}"}, timeout=600)
   before = api.post("/api/evaluate").json()
   level = api.get("/api/project/levels").json()[0]
   ch = level["channels"][0]
   rev = api.get("/api/project").json()["rev"]
   api.patch(f"/api/project/levels/{level['label']}/channels/{ch['pair']}/{ch['l']}/{ch['s']}",
             json={"width": 2 * ch["width"], "rev": rev}).raise_for_status()
   after = api.post("/api/evaluate").json()
   print(before["objective"], "->", after["objective"], [(d["name"], d["chi2"]) for d in after["datasets"]][:3])
   api.post("/api/project/undo", json={"rev": rev + 1}).raise_for_status()

``tests/pyazr/ui/server_test.py`` runs this block against a server process.
