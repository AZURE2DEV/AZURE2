pyazr — the Python Interface
============================

``pyazr`` drives AZURE2 from Python in-process: the R-matrix engine is compiled
into a pybind11 extension module (``_azure2``) and a session is a real C++
``AZUREAPI`` object living in the interpreter. There are no subprocesses, no
sockets and no port bookkeeping — the full R-matrix engine, the same code the
GUI uses, becomes callable from a script.

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

AZURE2 itself is C++ and is not installed by pip. The engine is the pybind11
module ``_azure2``, built by CMake with ``USE_API=ON`` (the default) into the
``pyazr/`` package directory; ``import pyazr`` loads it from there (or from the
pip-installed package data).

First steps
-----------

.. code-block:: python

   from pyazr import azure2
   import numpy as np

   with azure2("13N.azr") as azr:
       best = np.asarray(azr.params_rwa, float)   # free parameters
       chi2 = np.sum(azr.calculate_chi2_rwa(best))
       print(chi2)

The context manager releases the engine — several MB of compound nucleus and
data — rather than waiting for the garbage collector; without it, call
``close()``. A closed session raises if used again.

.. note::

   **One session per** ``azure2()`` **object, and any number of them.** Each
   ``azure2()`` owns an independent engine, so several can coexist in one
   interpreter — useful for sweeping model variants. They do not interfere:
   each enters its own directory for the duration of a call and leaves again,
   so your own working directory is never changed.

   **But drive them from one thread.** The engine is not reentrant and keeps
   process-wide state. For parallelism give each *process* its own session
   (see `Running in parallel`_).

   **The** ``.azr`` **stores its** ``output/``, ``checks/`` **and data paths
   relative to itself**, which is what ``cwd=`` handles; it defaults to the
   ``.azr``'s own directory, so paths you pass ``pyazr`` are resolved from
   wherever *you* are.

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

``residual_jacobian`` computes the level-energy, reduced-width and
normalization columns for roughly the cost of two forward evaluations, which
makes Gauss-Newton and Levenberg-Marquardt fits cheap. It also yields the
per-segment χ² the API does not expose directly, by slicing the residuals by
segment length.

Energy-shift columns are finite-differenced instead, at two extra residual
evaluations each. A shift translates the energy axis of a whole segment, so
what is wanted is :math:`\partial\,\mathrm{model}/\partial E`, and AZURE2
applies a shift by rebuilding every energy-dependent quantity of the affected
points --- level matrix, penetrabilities, Coulomb and hard-sphere phases,
external-capture amplitudes, and the target-integration sub-point grid. A fit
with many free shifts is dominated by their columns.

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
opening a fresh ``azure2()`` session on the new file.

One rule catches people out: **all levels of one J**\ :sup:`π` **share a
channel set**. Adding a level to an existing J\ :sup:`π` clones that group's
channel structure; only a brand-new J\ :sup:`π` lets your channel list define
the group.

Also available: ``set_channel_radius``, ``set_segment_norm``,
``set_segment_active``, ``set_segment_datafile``, ``add_data_segment``,
``remove_data_segments``, ``clear_data_segments``, ``set_extrapolations`` and
``apply_fit``, which writes a fit result back into a ``.azr`` you can reopen in
the GUI.

.. warning::

   Adding or removing data segments changes which energies AZURE2 evaluates.
   The external-capture integrals are cached in ``output/intEC.dat`` /
   ``output/intEC.extrap`` keyed on the *grid*, not on the segment selection,
   and AZURE2 silently reuses a stale file.  After any data edit, delete those
   caches (or give the edited model its own output directory) so the integrals
   are recomputed; ``azure2.recalculate_external_capture()`` forces it inside a
   live session.  See ``pyazr/examples/edit_model.py``.

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

Fetching data from EXFOR and NDS
--------------------------------

``pyazr.nds`` wraps the IAEA web services for experimental nuclear data and
converts the results straight into AZURE2 form (network required).  EXFOR
holds the measured reaction data; LiveChart/ENSDF the evaluated level schemes.

.. code-block:: python

   from pyazr import nds

   # search EXFOR: target, reaction (proj,exit) and quantity filters
   hits = nds.search_exfor(target="C-13", reaction="p,g", quantity="SIG")
   # -> [ExforDataset O2599004 6-C-13(P,G)7-N-14,,SIG,,SFC n=31 ...]

   # fetch one dataset, convert to an AZURE2 data file + segment kwargs
   data = nds.fetch_exfor("O2599004")                 # S-factor, B*KEV
   kw = data.to_azr("run/data", entrance=1, exit=2,
                    observable="total-capture")       # lab E, barns
   AzrModel.from_file("13N.azr").add_data_segment(**kw).write("13N_new.azr")

   # level scheme of the compound nucleus, from ENSDF
   for lv in nds.fetch_levels("14n"):
       print(lv.energy_mev, lv.jp, lv.half_life)

   # the paper behind the dataset, resolved to a DOI
   ref = nds.reference("O2599004")
   doi = nds.resolve_doi(ref)   # -> 10.1103/physrevlett.131.162701

Frames and units are handled automatically: ``EN-CM`` energies are converted
to lab, ``B*KEV`` S-factors to barns via the Sommerfeld factor, ``NB/SR``
differentials to b/sr.  Pick ``observable`` to match the data's angle frame
(EXFOR ``ANG-CM`` → ``"differential-cm"`` or ``"analyzing-power"``; lab
``ANG`` → ``"differential"``).  See the ``nds-explorer`` skill for the full
quantity/unit reference and failure modes.

One case cannot be decided automatically. Ratio-to-Rutherford data
(``,,RTH``) is dimensionless with an angle column, which is exactly what an
analyzing power looks like — and ``x4get`` drops the quantity suffix that
distinguishes them, so only the *search* result knows. ``to_azr`` warns and
passes such values through; give it ``rutherford=True`` to multiply by the
Coulomb cross section, or ``rutherford=False`` once you have checked.

The module needs nothing beyond NumPy and the standard library, and doubles as
a command-line tool::

   python -m pyazr.nds search --target C-13 --reaction p,g --quantity SIG
   python -m pyazr.nds download O2599004 -o data/skowronski.dat
   python -m pyazr.nds reference O2599004

The Qt setup utility has its own EXFOR dialog, backed by
``gui/src/ExforData.cpp``. The two are independent implementations of the same
Web-API, so a parsing rule learned by either belongs in both.

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

Each call writes a temporary model requesting exactly those energies and opens
one ``azure2()`` session on it, so pass every energy in a single call rather
than looping.

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

One trap is worth knowing about: comparing against thin-target data requires
segments with **no target integration**, because :math:`A_y` averaged over a
thick target is weighted by the cross section and Rutherford scattering drives
that average towards zero.

``residual_jacobian`` and ``chi2_and_grad`` do differentiate :math:`A_y`
exactly. The one exception is an analyzing-power point that also carries target
integration: that is a ratio of two integrals, is not differentiated
analytically, and makes the analytic Jacobian unavailable for the whole fit
rather than returning something approximate.

See :doc:`../theory/polarization_theory` for the formalism,
:doc:`../theory/polarization_implementation` for the implementation, and
``tests/13N`` (segments 11--16) for a worked comparison against measured data.

Cleaning up
-----------

There are no stray processes to reap: a session runs in-process, so when the
interpreter exits, the engine goes with it. No orphans, no ``pyazr.cleanup``,
no port collisions.

There is still memory. A session holds its compound nucleus and data — several
MB for a small model, more for a real one — so in a loop over model variants
use the context manager (or ``close()``) rather than leaving it to the garbage
collector, which may otherwise hold two models at once.

Dimensionless widths
--------------------

.. code-block:: python

   t = azr.dimensionless_widths(best)
   print(t.particles.table())         # theta^2 per particle channel
   print(t.photons.nonzero.table())   # Weisskopf units per gamma channel
   [c for c in t.particles if c.theta2 and c.theta2 > 1]   # unphysical

θ² > 1 exceeds the Wigner limit and is unphysical; γ-ray strengths in Weisskopf
units are the corresponding sanity check for capture.

The external region, and the caches
-----------------------------------

The quantities that describe everything outside the channel radius can be asked
for directly. They are what the penetrabilities in the level matrix are built
from, what sets the hard-sphere phase, and what the external-capture integrals
integrate.

.. code-block:: python

   c = azr.coulomb_functions(pair=1, energies=E, L=0)   # radius=0 -> channel radius
   c["F"], c["G"], c["P"], c["S"], c["delta_hs"]

   paths = azr.ec_integrals(pair=1, energies=E)         # one entry per EC pathway
   paths[0]["li"], paths[0]["lf"], paths[0]["radiation"], paths[0]["value"]

   azr.cache_stats()   # queries, hits, hit_rate, entries, keys, disabled_keys

The Coulomb functions follow the run's own configuration, so the same call
returns the accurate routine's values, GSL's (``--gsl-coul``), or the Numerov
solution through a nuclear potential (the hybrid model of the
``<potential>`` block). Comparing them is how one sees what those options do to
the external region.

The hybrid nuclear potential, per particle pair
----------------------------------------------

A nuclear potential belongs to a **particle pair**: it bends the radial wave
functions of that channel and no other. Each pair therefore carries its own
setting, and ``pair=0`` is the default that a pair without one inherits.

.. code-block:: python

   azr.nuclear_potentials()          # {0: default, 1: ..., 2: ...}, resolved
   azr.nuclear_potential(1)          # one pair; .own says if it is its own

   azr.set_nuclear_potential(1, type="WoodsSaxon", enabled=True,
                             V0=20.0, R=3.6, a=0.6)
   azr.set_nuclear_potential(2, type="Gaussian", enabled=True, V0=60.0, r0=4.0)
   azr.set_nuclear_potential(2, enabled=False)   # off for pair 2 alone
   azr.clear_nuclear_potential(1)                # follow the default again

Anything left at ``None`` keeps what the pair already resolves to, so switching
one pair off is just ``set_nuclear_potential(2, enabled=False)``. Setting
several pairs is cheaper with ``reinitialize=False`` on all but the last: the
potential is read when ``CoulFunc`` is constructed, deep inside a calculation,
so it only reaches the model when the model is rebuilt.

Two things that will bite:

- **Re-read the parameter vector afterwards.** The potential moves the
  penetrabilities and shift functions, and those are what map physical widths
  to reduced-width amplitudes --- so ``params_rwa`` is re-derived by the
  rebuild, and a vector captured beforehand no longer describes the same model.
  Feed the old one back and the chi-squared is quietly wrong. This is the same
  caveat as :meth:`~pyazr.azure2.set_channel_radius`, for the same reason.
- **A fit made without the potential is not a fit made with it.** Refit before
  reading anything off the model.

The equivalent in a ``.azr`` is the ``<potential>`` block, where keys before the
first ``pair=`` are the default and a ``pair=`` line opens a section that starts
from it:

.. code-block:: text

   <potential>
   useHybridPotential=1
   potentialType=0        # 0 = Woods-Saxon, 1 = Gaussian
   V0=40
   R=3.6
   a=0.6
   pair=1                 # pair 1 alone, starting from the default above
   V0=20
   pair=2
   useHybridPotential=0
   </potential>

A file with no ``pair=`` line reads exactly as it did when the model was
global. The GUI's Nuclear Potential tab edits one pair at a time and writes
the same block, and ``--no-gui`` reads it through the same parser, so the three
cannot drift. Worked example: ``pyazr/examples/nuclear_potential.py``;
``tests/hybrid_potential`` pins the behaviour.

External-capture integrals are the most expensive part of a capture
calculation, which is why the Coulomb functions they need are memoized.
``cache_stats`` makes that visible: asking for the same integrals twice on
``16O(p,gamma)17F`` takes 14.9 s and then 0.44 s, with the hit rate rising from
82% to 91%. ``disabled_keys`` counts the memos that have given up because too
few of their entries were being asked for twice --- which is what a *varying*
energy shift produces, since it moves every point energy at every iteration.

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
   * - ``edit_model.py``
     - The full edit loop: add/remove resonances *and* data segments,
       recalculate the external-capture integrals when the data change, and
       save the edited file.
   * - ``exfor_fetch.py``
     - Pulling data from the IAEA EXFOR/NDS web services (cross sections,
       S-factors, differentials, analyzing powers, levels) and dropping it
       into a model as new segments — with the DOI of the paper behind the
       dataset.
   * - ``deactivate_level.py``
     - Switching a resonance off without removing it.
   * - ``transform_widths.py``
     - Reduced width amplitudes to physical partial widths.
   * - ``dimensionless_widths.py``
     - θ² and Weisskopf units for a whole fit.
   * - ``coulomb_functions.py``
     - Coulomb functions, penetrability and hard-sphere phase over an
       energy grid.
   * - ``ec_integrals.py``
     - External-capture integrals per pathway, and what caching them buys.
   * - ``save_fit_to_azr.py``
     - Writing a fit back into a ``.azr``, verifying it round-trips.
   * - ``uncertainty_band.py``
     - Cross-section uncertainty bands from a saved fit covariance.
   * - ``sensitivities.py``
     - Which parameters a dataset constrains, as d ln σ / d ln p; and the
       analytic sensitivities checked against finite differences.
   * - ``fit_emcee.py``
     - MCMC sampling with ``emcee``, one in-process engine per pool worker.
   * - ``fit_zeus.py``
     - The same with ``zeus``.
   * - ``angular_distribution.py``
     - Legendre coefficients at chosen energies, with an optional plot.
   * - ``per_dataset_chi2.py``
     - Slicing the residual vector to get the χ² of each experiment.
   * - ``decompose_cross_section.py``
     - Separating resonance, interference and background contributions at
       the fitted parameters.
   * - ``sfactor_extrapolation.py``
     - An S factor below the measured range, and S(0).
   * - ``channel_radius_scan.py``
     - χ² against the channel radius, one in-process session per radius.
   * - ``reaction_rate.py``
     - N\ :sub:`A`\ ⟨σv⟩ by integrating the extrapolated cross section
       over a Maxwell--Boltzmann distribution.

Running in parallel
-------------------

There are no AZURE2 instances to pool — every ``azure2()`` object is already an
independent engine in-process. But the engine is not thread-safe, so the unit
of parallelism is the *process*: give each worker its own session by
constructing it at module scope in the worker module. That works whichever way
the pool starts its workers — under ``spawn`` each re-imports the module and
builds its own engine, under ``fork`` each inherits a copy — and each then
evaluates one χ² per walker:

.. code-block:: python

   # worker module
   azr = azure2("13N.azr")

   def log_prob(theta):
       return -0.5 * np.sum(azr.calculate_chi2_rwa(theta)) + prior(theta)

See ``pyazr/examples/fit_emcee.py`` and ``fit_zeus.py`` for the full recipe.

The same applies inside ``pyazr`` itself: ``sensitivities``,
``uncertainty_bands`` and ``extrapolation_bands`` take ``nprocs=N``, which
spreads the finite-difference columns over ``N`` worker processes, each with
its own engine and its own output directory (seeded from yours, so the
external-capture integrals are reused rather than rebuilt).

It is not free, and it is not always a win. Every worker pays one model
initialisation before it evaluates a single column, so the columns have to be
worth more than that startup — which in practice means hundreds of free
R-matrix parameters. On the 14-parameter ``tests/13N`` model it *loses* at any
``nprocs``, both in data mode (0.8 s serial) and on an 1845-point
extrapolation grid (53 s serial, 56 s over four workers). It is also
irrelevant to the default ``method="analytic"``, which is one call whatever the
parameter count. Reach for it only when the analytic path is unavailable and
the model is large, and measure rather than assume.

Data mode and extrapolation mode
--------------------------------

``data_mode()`` (the default) evaluates the ``<segmentsData>`` segments — this
is what χ² uses. ``extrap_mode()`` switches to the ``<segmentsTest>`` grids for
predictions on arbitrary energies and angles. Both re-initialise the engine, so
switch sparingly and batch the work.

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
