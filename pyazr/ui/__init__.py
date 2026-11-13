"""The web UI of AZURE2: its server and API (docs/upgrade/web_ui.md).  Not imported by ``import pyazr``.

Start it with ``python -m pyazr.ui project.azr [--port N] [--threads N] [--token T]``: it prints the URL with its
token; ``GET /api/describe`` lists every route.  From Python: :func:`pyazr.ui.server.main` and
:func:`pyazr.ui.server.create_app`.  It needs the optional packages fastapi and uvicorn (``pip install pyazr[web]``).

Modules: ``server`` (the app, :class:`~pyazr.ui.server.Server`, ``main``), ``routes`` (the route table), ``state``
(the project text and its history), ``workers`` (worker processes, the compute budget), ``jobs`` (the job table),
``engine`` (the adapter to today's pyazr session), ``errors`` (exceptions to HTTP statuses).
"""
