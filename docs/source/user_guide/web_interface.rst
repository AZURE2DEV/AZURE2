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

Not there yet
-------------

The tabs **Fit**, **Results**, **Physics** and **THM** say which step brings them.  Until then a fit runs through
the API (``POST /api/jobs``) or pyazr, and the ``<thm>`` section is edited through the API.  Adding or removing a
particle pair, generating channels from a J^π, the session's options as editable keys, the file browser, EXFOR and
ENSDF imports, the Wigner limits and physical widths of the channels, plots of whole fits and extrapolations, and
live updates pushed by the server (the page asks every three seconds whether another client changed the project)
come later.  The page was checked in a headless Chrome; plotly.js 4.1.2 (MIT licence) is part of the package, so
the page loads nothing from the network, and without it the plots appear as tables of their numbers.
