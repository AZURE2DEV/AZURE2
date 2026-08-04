pyazr — the Python Interface
============================

``pyazr`` drives AZURE2 from Python. It spawns headless AZURE2 processes and
talks to them over the socket API, so the full R-matrix engine — the same code
the GUI uses — becomes callable from a script.

This is the route to anything the GUI does not offer directly: custom
minimizers, external samplers, parameter scans, systematic studies over many
model variants, and publication figures built from the model rather than from
exported files.

Installation
------------

.. code-block:: bash

   pip install -e .          # from the repository root
   pip install -e ".[all]"   # plus what the examples need

Core dependencies are NumPy, mpmath and SciPy. mpmath is imported when
``pyazr.transform`` loads, so it is required rather than optional; the samplers
and plotting libraries the examples use are the extras.

AZURE2 itself is C++ and is not installed by pip. ``pyazr`` finds the binary
via ``$AZURE2_BINARY``, then ``build/src/AZURE2``, then ``$PATH``, and needs it
built with ``USE_API=ON`` (the default).

First steps
-----------

.. code-block:: python

   from pyazr import azure2
   import numpy as np

   with azure2("13N.azr") as azr:
       best = np.asarray(azr.params_rwa, float)   # free parameters
       chi2 = np.sum(azr.calculate_chi2_rwa(best))
       print(chi2)

The context manager shuts the subprocess down; without it, call ``close()``.

.. important::

   **One AZURE2 session per process.** Opening a second ``azure2()`` in the
   same interpreter desynchronises the parameter metadata. To sweep model
   variants, loop in the shell — one process per variant — or build a single
   temporary ``.azr`` carrying everything the run needs.

   **Run from the directory containing the** ``.azr``. The file stores its
   ``output/``, ``checks/`` and data paths relative to itself. ``pyazr`` does
   this for you via ``cwd=``, which defaults to the ``.azr``'s directory.

Two parameter conventions
-------------------------

AZURE2 exposes its parameters in two spaces, and mixing them silently produces
wrong answers.

**Reduced width amplitudes** (``*_rwa`` methods) are the natural fit space and
the only one with analytic derivatives. Default to these.

.. code-block:: python

   azr.params_rwa               # the free parameter vector
   azr.calculate_chi2_rwa(x)    # chi-squared per segment
   azr.calculate_rwa(x)         # cross section per segment
   azr.residual_jacobian(x)     # residuals + analytic Jacobian

**Physical parameters** (level energies in MeV, partial widths in eV) are what
``parameters.out`` reports and what the GUI shows.

.. code-block:: python

   azr.params                   # physical vector
   azr.calculate(x)             # cross sections from physical parameters
   azr.transform_rwa(x)         # rwa -> physical

Both vectors contain only the **free** parameters, in ``.azr`` order.

.. warning::

   The ``gamma`` field in a ``.azr`` ``<levels>`` block is **not** a reduced
   width amplitude — it holds the physical value: Γ in eV for an open particle
   channel, an ANC for a closed one, Γ\ :sub:`γ` in eV for a photon channel.
   The two differ by factors of 10² to 10⁷, and a file written with the wrong
   one loads without complaint and is wrong. Convert with ``transform_rwa``.

Inspecting a model
------------------

Every parameter carries what it physically *is*, so you never have to parse a
name like ``width_1_2``:

.. code-block:: python

   print(azr.level_scheme)          # pairs -> J-groups -> levels -> channels
   print(azr.parameters.table())    # every parameter with its metadata
   print(azr.datasets.table())      # per segment: file, reaction, E range, norm error
   print(azr.pairs)                 # masses, charges, spins, separation energies

   for w in azr.parameters.widths.free:
       print(w.name, w.jpi, "L =", w.L, "S =", w.S, "pair", w.pair)

Filtered views (``.free``, ``.energies``, ``.widths``, ``.norms``,
``.shifts``) and lookups (``.by_level(...)``, ``.by_name(...)``,
``.by_physical_level()``) are all available. A ``LevelKey`` prints as
``5/2-#2@6.588MeV``; ``(jgroup, level)`` is its identity, since AZURE2 restarts
level numbering inside every J-group.

Chi-squared and derivatives
---------------------------

.. code-block:: python

   chi2   = np.sum(azr.calculate_chi2_rwa(x))
   val, g = azr.chi2_and_grad(x)          # analytic gradient
   r, J   = azr.residual_jacobian(x)      # sum(r**2) == chi2

``residual_jacobian`` computes the whole Jacobian for roughly the cost of two
forward evaluations, which makes Gauss-Newton and Levenberg-Marquardt fits
cheap. It also yields the per-segment χ² the API does not expose directly, by
slicing the residuals by segment length.

.. warning::

   These return the **data** χ² only. AZURE2's own fit objective additionally
   penalises every free normalization and energy shift. Minimise the residuals
   alone and the normalizations will drift to absorb every discrepancy,
   reaching a χ² that AZURE2 would never find. Append the penalty rows — see
   :doc:`chi_squared` and ``pyazr/examples/`` for the recipe.

Editing the model
-----------------

``AzrModel`` parses the ``<levels>`` block and re-emits everything else
verbatim, so the original file is never modified:

.. code-block:: python

   from pyazr import AzrModel

   mdl = AzrModel.from_file("7Be.azr")
   mdl.remove_level(jpi="1/2+", energy=20)
   mdl.add_level(J=1.5, parity=+1, energy=8.6,
                 channels=[dict(pair=1, L=2, S=0.5, gamma=1000.0, fixed=False)])
   path = mdl.write("_variant.azr")

Because AZURE2 reads its model from the file, an edited scheme takes effect by
launching a fresh instance from the new file.

One rule catches people out: **all levels of one J**\ :sup:`π` **share a
channel set**. Adding a level to an existing J\ :sup:`π` clones that group's
channel structure; only a brand-new J\ :sup:`π` lets your channel list define
the group.

Also available: ``set_channel_radius``, ``set_segment_norm``,
``set_segment_active``, ``set_segment_datafile``, ``set_extrapolations`` and
``apply_fit``, which writes a fit result back into a ``.azr`` you can reopen in
the GUI.

Decomposing a cross section
---------------------------

At a fitted point, without refitting, each number is a component's raw
contribution:

.. code-block:: python

   azr.extrap_mode()
   full = azr.calculate_rwa(best)
   off  = azr.calculate_rwa(azr.without_level(best, jpi="5/2-", energy=6.588))
   only = azr.calculate_rwa(azr.only_level(best,   jpi="5/2-", energy=6.588))

``full - off`` is everything the level does, resonance plus interference;
``only`` minus the non-resonant background is the bare resonance. The
difference between them is pure interference, which is block-diagonal in
J\ :sup:`π` — only same-J\ :sup:`π` levels interfere, and that is a useful
check on any decomposition.

Angular distributions
---------------------

AZURE2 writes an angular distribution as

.. math:: W(\theta) = \sum_k a_k P_k(\cos\theta)

and computes the :math:`a_k` **only** for segments declared as angular
distributions -- ``observable="angular-distribution"`` with an ``order``.
Everything else yields empty arrays.

For the grids a model already declares, on a live instance:

.. code-block:: python

   dists = azr.calculate_angular_dists_rwa(x)
   # one entry per segment; each a list with one array of coefficients per point

For arbitrary energies, which is usually what is wanted:

.. code-block:: python

   from pyazr import angular_distribution

   e_cm, coeffs = angular_distribution("model.azr", [0.05, 0.1, 0.2],
                                       entrance=1, exit=2, order=4)

``coeffs[i, k]`` is :math:`a_k` at energy *i*. Input energies are **lab** by
default (``lab=False`` for centre-of-mass); the returned energies are always
centre-of-mass. Energies AZURE2 could not evaluate come back as ``NaN`` rather
than being dropped, so the rows always line up with the input.

Each call writes a temporary model requesting exactly those energies and starts
one AZURE2 process, so pass every energy in a single call rather than looping.

Reading the result: :math:`a_0` is the normalisation, so :math:`a_0 = 1` with
everything else zero means isotropic. That is the correct answer for a
resonance formed in an s-wave -- the compound nucleus has no preferred
direction -- and is what :sup:`3`\ H+d gives at low energy. Anisotropy shows
as the higher orders departing from zero.

Worked example: ``pyazr/examples/angular_distribution.py``.

Analyzing power
---------------

The vector analyzing power :math:`A_y` is observable code 7. Declare a segment
with ``observable="analyzing-power"`` and it is reported in place of the cross
section, so :math:`\chi^2`, fitting and plotting need no special handling:

.. code-block:: python

   ay = azr.calculate_analyzing_power_rwa(azr.params_rwa)

Data files carry ``E_lab  theta_cm  A_y  dA_y`` -- angles are centre-of-mass,
unlike an ordinary differential segment. Leave ``vary_norm`` off: a
normalization factor is meaningless for a ratio.

Two traps are worth knowing about. Comparing against thin-target data requires
segments with **no target integration**, because :math:`A_y` averaged over a
thick target is weighted by the cross section and Rutherford scattering drives
that average towards zero. And the analytic-derivative paths
(``residual_jacobian``, ``chi2_and_grad``) do not know about :math:`A_y` yet, so
do not use them on analyzing-power segments.

See :doc:`../theory/polarization_theory` for the formalism,
:doc:`../theory/polarization_implementation` for the implementation, and
``tests/13N`` (segments 11--16) for a worked comparison against measured data.

Cleaning up stray processes
---------------------------

A session stops the instances it started, and an :mod:`atexit` hook covers a
script that forgets to ``close()``. What nothing can cover is the interpreter
dying outright -- ``kill -9``, or a notebook kernel restarting -- which leaves
AZURE2 processes running with nobody to talk to them.

.. code-block:: bash

   python -m pyazr.cleanup          # list what would be reaped
   python -m pyazr.cleanup --kill   # reap it
   pyazr-cleanup --kill             # same, when pip-installed

Only instances whose **parent process has died** are candidates, so a
concurrently running script's instances are never touched, however many it has.
``pyazr.find_orphans()`` and ``pyazr.kill_orphans()`` are the API equivalents.

Dimensionless widths
--------------------

.. code-block:: python

   t = azr.dimensionless_widths(best)
   print(t.particles.table())         # theta^2 per particle channel
   print(t.photons.nonzero.table())   # Weisskopf units per gamma channel
   [c for c in t.particles if c.theta2 and c.theta2 > 1]   # unphysical

θ² > 1 exceeds the Wigner limit and is unphysical; γ-ray strengths in Weisskopf
units are the corresponding sanity check for capture.

Examples
--------

Worked scripts ship in ``pyazr/examples/``:

.. list-table::
   :widths: 32 68
   :header-rows: 1

   * - Script
     - What it shows
   * - ``print_scheme.py``
     - Reading a model's level scheme and dataset provenance.
   * - ``edit_scheme.py``
     - Adding and removing levels, writing a new ``.azr``.
   * - ``deactivate_level.py``
     - Switching a resonance off without removing it.
   * - ``transform_widths.py``
     - Reduced width amplitudes to physical partial widths.
   * - ``dimensionless_widths.py``
     - θ² and Weisskopf units for a whole fit.
   * - ``save_fit_to_azr.py``
     - Writing a fit back into a ``.azr``, verifying it round-trips.
   * - ``uncertainty_band.py``
     - Cross-section uncertainty bands from a saved fit covariance.
   * - ``fit_emcee.py``
     - MCMC sampling with ``emcee``, one AZURE2 instance per pool worker.
   * - ``fit_zeus.py``
     - The same with ``zeus``.
   * - ``angular_distribution.py``
     - Legendre coefficients at chosen energies, with an optional plot.

Running several instances
-------------------------

``nprocs=N`` spawns N independent AZURE2 instances; ``proc=i`` selects one.
This is how a sampler evaluates one χ² per walker in parallel:

.. code-block:: python

   azr = azure2("13N.azr", nprocs=8)
   chi2 = np.sum(azr.calculate_chi2_rwa(theta, proc=worker_index))

Ports are assigned by the OS, so parallel sessions never collide.

Data mode and extrapolation mode
--------------------------------

``data_mode()`` (the default) evaluates the ``<segmentsData>`` segments — this
is what χ² uses. ``extrap_mode()`` switches to the ``<segmentsTest>`` grids for
predictions on arbitrary energies and angles. Both re-initialise every
instance, so switch sparingly and batch the work.

Segment indexing differs between them, and getting it wrong misaligns results
silently: in data mode segment *i* is ``azr.datasets[i]``; in extrapolation
mode **only active test segments are returned**, so segment *i* is
``azr.extrapolations.active[i]``.

.. note::

   Delete ``output/intEC.extrap`` whenever the ``<segmentsTest>`` grid changes.
   AZURE2 caches external-capture integrals there and silently reuses them on a
   different grid, which corrupts capture cross sections. It is safe to delete;
   it only costs time to rebuild.

Frames and units
----------------

**Input is lab frame** — the energies and angles in a ``.azr`` and in data
files. **All output and every API result is centre-of-mass.** This includes
extrapolation grids: the energies you *set* are lab, the energies you *get
back* are c.m. Convert with ``E_lab = E_cm · (m_beam + m_target) / m_target``.

Plot different reaction channels against **excitation energy**
(``calculate_excitation_energy``); it is the only axis shared by all entrance
channels.
