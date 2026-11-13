The web UI's API
================

The web UI of AZURE2 is a server on your machine that holds one project and answers JSON requests; the page in
the browser is one client of it, a script or an AI agent is another.  This page describes the server as it stands
today: the API without the page (the design is ``docs/upgrade/web_ui.md``).

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

Every request needs the token: the header ``Authorization: Bearer <token>``, or ``?token=<token>`` in the URL
(which sets a cookie).  ``GET /api/describe`` lists every route with its documentation and the pyazr call behind
it; ``GET /openapi.json`` is the machine-readable form.

What the server keeps
---------------------

- **The project text is the state.**  Every edit goes through :class:`pyazr.Project` and is validated by the engine
  before it is accepted; a refused edit answers 422 with the engine's diagnostics and changes nothing.  Every
  accepted edit gives a new revision number, and every edit names the revision it was made on (``"rev"`` in the
  body, or an ``If-Match`` header): an edit on an old revision answers 412, so two clients cannot overwrite each
  other.  ``POST /api/project/undo`` and ``/redo`` step through the history of texts.
- **Your file is written only by** ``POST /api/project/save``.  Until then the current text is kept in
  ``.<name>.ui.azr`` beside the project; a server started after a crash offers it back (``crash_copy`` in
  ``GET /api/project``, ``POST /api/project/restore``).
- **The engine runs in worker processes.**  Evaluations use a session of the current text, built again when the
  text changes; the free parameters' working values ``x`` (``GET``/``PATCH /api/x``, in the fit basis) move without
  a rebuild.  The χ² per data set is the program's (the rows of its ``chiSquared.out``, six significant digits).
- **Jobs** (``POST /api/jobs``): recipe fits (scipy L-BFGS-B, iminuit MIGRAD + HESSE) in a process of their own,
  with progress (``GET /api/jobs/<id>``), a cooperative stop (``/stop``), a kill (``/kill``), a checkpoint of the
  best point, and ``/accept`` to write the best point into the project text.  With a budget of one thread, an
  evaluation asked while a job runs answers 409 with the job's id.
- Routes that need a pyazr call that does not exist yet (the engine's own fitters, the reaction rate, the physics
  table, ...) answer 501 naming that call.

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
