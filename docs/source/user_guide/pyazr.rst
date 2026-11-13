pyazr — the Python Interface
============================

``pyazr`` drives AZURE2 from Python in-process: the R-matrix engine is compiled
into a pybind11 extension module (``_azure2``) and a session is a real C++
``AZUREAPI`` object living in the interpreter. There are no subprocesses, no
sockets and no port bookkeeping — the full R-matrix engine, the same code the
command-line program runs, becomes callable from a script.

This is the route to anything the command-line menu does not offer directly: custom
minimizers, external samplers, parameter scans, systematic studies over many
model variants, and publication figures built from the model rather than from
exported files.

Installation
------------

``pip install .`` from the repository's root builds the engine and installs
``pyazr`` with it, the commands ``pyazr`` and ``AZURE2`` and the data files;
the Coulomb table is fetched or built once with ``pyazr coulomb-table fetch``
or ``build``. Python 3.10 or later, NumPy 1.24 or later. What the wheel holds,
the SRIM citation and the developers' in-tree build:
:doc:`../getting_started/installation` ("The Python package").

First steps
-----------

.. code-block:: python

   from pyazr import Session
   import numpy as np

   with Session("13N.azr") as azr:
       best = np.asarray(azr.params_rwa, float)   # free parameters
       chi2 = np.sum(azr.calculate_chi2_rwa(best))
       print(chi2)

The context manager releases the engine — several MB of compound nucleus and
data — rather than waiting for the garbage collector; without it, call
``close()``. A closed session raises if used again.

.. note::

   **One session per** ``Session()`` **object, and any number of them.** Each
   ``Session()`` owns an independent engine, so several can coexist in one
   interpreter — useful for sweeping model variants. They do not interfere:
   each reads and writes the model's relative paths in its own directory
   (``cwd``, by default the ``.azr`` file's), and the process's working
   directory is neither used nor changed.

   **Two sessions may run on two threads at once; calls on one session run
   one after the other** (each holds the session's lock).  A session owns its
   model, results, scratch and diagnostics' counters; the evaluations release
   the GIL, so two sessions' χ² or Jacobians run side by side and give the
   bytes each gives alone (``tests/pyazr/concurrent_sessions_test.py``).  What a process still
   shares: immutable tables (the Coulomb tables, SRIM's stopping data, the
   bytes of input files read with the cache, keyed by path and file stamp),
   the environment switches (read once, when the first session is made: a
   variable set later changes nothing), GSL's error handler, and the state of
   an MCMC run, which is process-wide until the sampler's run object exists
   (run one MCMC per process).  A session's set-up keeps the GIL.  For more
   cores than one evaluation uses, give each *process* its own session (see
   `Running in parallel`_).

   **The** ``.azr`` **stores its** ``output/``, ``checks/`` **and data paths
   relative to itself**, which is what ``cwd=`` handles; it defaults to the
   ``.azr``'s own directory, so paths you pass ``pyazr`` are resolved from
   wherever *you* are.  ``pyazr.validate(path)`` and ``pyazr.describe(path)``
   read a project's data files by the same rule, relative to the ``.azr``'s
   directory (before November 2026: relative to the working directory).

Session options
---------------

``Session(file, *, cwd=None, data_mode=True, use_brune=True,
ignore_externals=True, transform=True, use_long_wavelength=True,
use_gsl_coul=False, use_rmc=False, use_park=False, ...)`` takes the runtime
options of the CLI (:doc:`../reference/command_line`) as keywords.  All but
``cwd``, ``data_mode`` and ``coulomb_table`` are also keys of the project's
``<config>`` under the same name (:ref:`options-as-config-keys`): a keyword
given wins over the key, the key over the default the signature shows, and a
keyword left out is not given.  ``session.options`` holds the values in
effect, ``session.option_sources`` where each came from (``"default"``,
``"project"``, ``"keyword"``); ``Project.set_config(use_brune=False)`` writes
``use_brune=0``.  Keywords that cannot hold together are refused as the
program refuses the flags (``AZURE2Error``, ``E_CONFIG_OPTIONS_CONTRADICT``,
with what is allowed and the fix): ``use_rmc=True`` needs ``use_brune=False``,
``use_park=True`` with ``use_brune=False`` is refused, and
``ec_convection_current=False`` needs ``use_long_wavelength=False``.
``pyazr.kinds()["options"]`` lists them all with type,
default, ``<config>`` key, flag and meaning. ``use_park=True`` is ``--use-park``:
``params_rwa`` then holds Park's amplitudes, ``park_norms()`` gives every
level's overlap :math:`J`, and ``basis`` is 0, 1 or 2 (standard, Brune, Park:
the ``#parametrization`` tag of ``param.sav``). With ``transform=False``
(``--no-transform``) nothing is converted: the numbers of the file are taken
as they stand as formal level energies and reduced-width amplitudes, so a
file that holds physical widths (eV) and energies is a different model.

Asking the package about itself
-------------------------------

Three calls answer, as plain data (dicts and lists that ``json.dumps`` takes), what an agent or a script would
otherwise look up in the sources:

.. code-block:: python

   import pyazr
   pyazr.describe()        # version, engine build (git, MCMC, NLopt, OpenMP), threads,
                           # the Coulomb table in use (path and content identity, or "exact"),
                           # the stopping library and its citation, the mass table
   k = pyazr.kinds()
   k["records"]            # record kinds of format 3 and their fields (type, grammar, units, doc)
   k["data"]               # data kinds (quantities): units, engine unit, energy frames, angle columns
   k["parameters"]         # parameter kinds (Parameter.kind): unit, owner fields, doc
   k["options"]            # session options: type, default, the program's flag, doc
   k["physics"]            # the physics functions the compiled module binds (group, examples, sweeps)
   k["stopping"]           # the stopping functions as rows of the same shape, each with SRIM's citation
   k["errors"]             # the error classes below, with their bases
   s = pyazr.Session("model.azr")
   s.describe()            # pairs, levels, data sets (kind, points, unit, frame), extrapolations,
                           # parameters (name, kind, unit, value, free), options in effect, output directory

``kinds()`` is read from the engine's own tables (the schema of ``AZURE2 --schema``, the parameter-kind table, the
binding's option table, the names the binding registers): a data kind or a parameter kind added to the engine appears
there without an edit of Python (``tests/pyazr/describe_test.py`` compares it with the program).
``pyazr.describe(path)`` (with a path) is, as before, the engine's description of a project (``AZURE2 --describe``),
and ``help(pyazr)`` names the objects and shows a session of ten lines.

Errors
------

Every error pyazr raises on purpose derives from ``pyazr.Error``; each class also derives from the builtin exception
the call raised before, so an existing ``except ValueError`` or ``except KeyError`` keeps catching it.

========================  ========================================  ================================================
class                     also a                                    raised for
========================  ========================================  ================================================
``FormatError``           ``ValueError``                            a project text or record edit the engine refuses
``ParameterError``        ``ValueError``                            a parameter vector, name or value that is wrong
``ParameterKeyError``     ``ParameterError``, ``KeyError``          an unknown parameter, level, pair, data set
``ParameterIndexError``   ``ParameterError``, ``IndexError``        an unknown pair key (``coulomb_functions``)
``ParameterVectorError``  ``ParameterError``, ``AZURE2Error``       a free vector of the wrong length, or a NaN
``EngineError``           ``RuntimeError``                          a set-up or evaluation the engine refuses;
                                                                    the class ``AZURE2Error``, which keeps its name
``DomainError``           ``ValueError``                            a physics argument outside its domain
``DomainZeroDivision``    ``DomainError``, ``ZeroDivisionError``    the same, where the formula divided by zero
``SessionClosed``         ``RuntimeError``                          a call on a closed session
``MissingFile``           ``FileNotFoundError``                     a project, directory or ``param.sav`` missing
========================  ========================================  ================================================

A message says what was wrong with the value given, what is allowed, and which call shows or fixes it, for example::

   pyazr.ParameterVectorError: calculate_objective_rwa: entry 0 of the free parameter vector (width_1_1) is nan;
   every entry must be a finite number.  Session.params_rwa is such a vector (the start values);
   Session.parameters.free names its entries in order.

A notice of the engine that a session meets -- an external-capture amplitude file that cannot be written because the
output directory does not exist, a normalization or shift file that cannot be written -- is a ``pyazr.EngineWarning``
naming the file and the consequence (it was a line on standard output before).  Make it an error with
``warnings.simplefilter("error", pyazr.EngineWarning)``.

Documented zeros stay zeros: ``coulomb_functions`` below threshold (E <= 0) gives the shift function alone, and an
energy at which the Coulomb routine fails gives zeros; ``whittaker_wave`` gives ``W`` 0 or inf outside the double
range (``mantissa`` and ``log10_scale`` keep the value).

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
``parameters.out`` reports.

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
   print(azr.datasets.table())      # per data set: record, file, reaction, E range
   print(azr.pairs)                 # masses, charges, spins, separation energies

   for w in azr.parameters.widths.free:
       print(w.name, w.jpi, "L =", w.L, "S =", w.S, "pair", w.pair)

A pair prints as ``Pair(#7 key=8 ...)``: ``#7`` is the engine's internal pair
number, which ``Parameter.pair`` carries, and ``key=8`` the pair key of the
``.azr`` file, which the engine's messages and ``isospin.out`` use. The two
differ whenever the file numbers its pairs in another order than the engine
counts them, so match a parameter to a reaction through ``azr.pairs`` rather
than by reading the key as the number.

``Pair.excitation`` is the pair's **total** excitation, light plus heavy
particle, so ``Pair.threshold`` = ``sep_energy + excitation`` is where its
channels open. ``Pair.light_excitation`` is the light particle's share (the
optional 32nd column of a ``<levels>`` line, November 2026) and
``Pair.heavy_excitation`` the heavy one's; ``Pair.label`` marks an excited
particle (``'2+*(4.439 MeV) + 0+'``). In a :class:`pyazr.Project` the pair's
``light_excitation`` and ``heavy_excitation`` are fields of
:meth:`~pyazr.Project.pairs`; ``set_pair(key, light_excitation=E)`` writes the
32nd column on every row of the pair (a row without the column gets it only
for a value other than 0), ``heavy_excitation=`` the heavy one's
(:doc:`../reference/data_formats`).

``azr.parameters`` is the table of the fit parameters (a :class:`pyazr.Parameters`; its class
``ParameterSet`` adds the views of earlier releases), one row per parameter in the order of the
engine's vector. Each row is the engine's record of the parameter (``Session.parameter_records()``),
written where the parameter is created, never parsed from its name, which is only a label:
``name``, ``kind``, ``unit``, ``start`` (in the fit vector's basis and ``unit``), ``value`` (as the
project states it: a partial width in eV, an ANC, an energy in MeV), ``fixed``, ``free_index`` (its
entry of the fit vector), ``bounds``, and the fields of its kind (a level's ``level_label``, ``J``,
``parity``, ``level_energy``; a channel's ``L``, ``S``, ``pair``, ``pair_key``, ``radiation_type``,
``wigner_limit``, ``width_is_amplitude``; a data set's ``segment_key``; a nuisance's
``prior_centre`` and ``prior_sigma``; ...). ``pyazr.kinds()["parameter_fields"]`` lists every field
with its unit, the kinds it describes and its doc; a field of another kind reads ``None``.

.. code-block:: python

   P = azr.parameters
   x = P.x0                                   # a copy of the start of the fit vector
   w = P.of_level("3/2- 3.5032").of_kind("width")
   x2 = P.put(x, w.free, 2 * P.take(x, w.free))   # the level's free widths doubled
   P.write("fit.sav", x2)                     # a parameter file, by name
   x3 = P.read("fit.sav")                     # and back
   print(P.free.table())                      # name, kind, unit, start, free, bounds, owner

The fit vector is the only vector: ``x0``, ``put``, ``take``, ``read`` and ``write`` work on it, and a
selection (``free``, ``of_kind(...)``, ``of_level(...)``, a slice) is again a table with the same
``x0``. ``names``, ``kinds``, ``units``, ``free_mask`` and ``bounds`` are its columns (``bounds`` the
limits the engine creates the parameter with, ``(-inf, inf)`` where it sets none). The views of
earlier releases (``.energies``, ``.widths``, ``.norms``, ``.shifts``, ``.sqrt_shifts``,
``.fn_shifts``, ``.abundances``, ``.screenings``, ``.densities``, ``.dopplers``, ``.thicknesses``,
``.cbkg``, ``.fixed``, ``.by_level(...)``, ``.by_physical_level()``) are ``of_kind`` and the level
selections under their old names. ``pyazr.snapshot.write_values(azr, project, x)`` writes the values
of ``x`` into a :class:`pyazr.Project` of the session's file in memory (what ``save_fit`` writes, without saving).

``pyazr.kinds()["parameters"]`` is the engine's table of kinds (key, unit, owner fields, flags), in
the order of the parameter vector (the abundances follow the level block):

===============  ===================================================
``kind``         parameters
===============  ===================================================
``energy``       level energies (``energy_*``)
``width``        reduced-width amplitudes (``width_*``)
``norm``         data-set normalizations
``shift``        data-set energy shifts (the constant term)
``shift_sqrt``   :math:`\sqrt{E}` shift coefficients
``shift_fn``     shift-function parameters a0, a1, ... (``fn_term`` = j)
``abundance``    isotopic abundances (after the level block;
                 ``sample_name`` its sample)
``screening``    electron-screening potentials
``density``      fitted target densities
``doppler``      fitted Doppler temperatures
``thickness``    fitted sample thicknesses
``cbkg``         THM coherent background (``cbkg_*``);
                 ``thm_experiment``, ``thm_jpi``, ``thm_exit_key``,
                 ``thm_channels`` and ``thm_term`` (0 Re c0, 1 Im
                 c0, 2 Re c1, 3 Im c1) say which coefficient
``apparatus``    a number of an apparatus record (a stopping scale, a
                 kernel's parameter)
``radiative``    a number of a level's ``gamma_decay`` record
                 (``decay_record``, ``decay_part``, ``final_state``)
===============  ===================================================

A THM segment's free norm is profiled, not a parameter, and has no entry. A ``LevelKey`` prints as
``5/2-#2@6.588MeV``; ``(jgroup, level)`` is its identity, since AZURE2 restarts
level numbering inside every J-group.

Chi-squared and derivatives
---------------------------

.. code-block:: python

   chi2   = np.sum(azr.calculate_chi2_rwa(x))
   val, g = azr.chi2_and_grad(x)          # analytic gradient
   r, J   = azr.residual_jacobian(x)      # sum(r**2) == chi2

Segments with a data covariance (``covariance <file>``, :doc:`chi_squared`)
return whitened residuals :math:`w = W r` and their Jacobian, so
``sum(r**2) == chi2`` still holds and Gauss-Newton fits need no change; the
whitened rows mix the points of the covariance group. With ``relvector``
terms the model in them is frozen at the first evaluation;
``azr.refreeze_covariance()`` makes the next evaluation freeze it again at its
parameters (fit, refreeze, evaluate at the fit, fit again).

``residual_jacobian`` computes the level-energy, reduced-width and
normalization columns for roughly the cost of two forward evaluations, which
makes Gauss-Newton and Levenberg-Marquardt fits cheap. It also yields the
per-segment χ² the API does not expose directly, by slicing the residuals by
segment length.

Energy-shift columns are finite-differenced instead, from two forward residual
evaluations each (no adjoint). A shift whose segments are self-contained (no
components, not a total-capture sum, no point mapping onto other segments) is
differenced on those segments' points alone; any other shift, and a fitted
target density, re-evaluates every point, in one parallel pass per side. (Until
30 October 2026 each side re-solved the whole Jacobian, adjoints included, and
re-evaluated every point whatever the shift.) A shift translates the energy axis of a whole segment, so
what is wanted is :math:`\partial\,\mathrm{model}/\partial E`, and AZURE2
applies a shift by rebuilding every energy-dependent quantity of the affected
points --- level matrix, penetrabilities, Coulomb and hard-sphere phases,
external-capture amplitudes, and the target-integration sub-point grid. A fit
with many free shifts is dominated by their columns. This covers the
:math:`\sqrt{E}` coefficients and the shift-function parameters
(``kind == "shift_fn"``, :attr:`Parameter.fn_term` = j) as well;
``penalties()`` and ``objective()`` include their penalties.

A fit that varies only some of the free parameters can say which of these
columns it does not need: ``residual_jacobian(x, skip=[...])`` (free-column
indices) leaves their finite differences out and returns those columns as zero;
every other column is unchanged. rmfit's ``polish`` passes the columns its
bounds hold (since 30 October 2026, before which a held shift cost its two
evaluations on every Jacobian all the same).

A fitted target density (``areal_density ... fit [prior ...]`` on a ``target``
record, shared by every data set that names the target,
:doc:`experimental_effects`) is a parameter ``<target>.areal_density`` of kind
``"density"`` (``ParameterSet.densities``; :attr:`Parameter.effect` is its target,
``prior_centre`` / ``prior_sigma`` its prior and ``layout_max`` the largest
density its grids span). Its column is finite-differenced the same way, and
``penalties()["density"]`` is its prior, which ``objective()`` adds.  The
field is read and written as any record field: ``project["d1"]["areal_density"]``,
``project.set("d1", areal_density=Nuisance(5e17, "fit", unit="atoms/cm2",
prior=(5e17, "5%")))``.

A fitted free-gas Doppler temperature (``doppler free_gas T K fit [prior T0
sigma] [grid_max Tmax]`` on a ``target`` record, :doc:`experimental_effects`) is
a parameter ``<target>.doppler`` of kind ``"doppler"``
(``ParameterSet.dopplers``; :attr:`Parameter.effect` is its target,
``prior_centre`` / ``prior_sigma`` its prior and ``layout_max`` the layout
temperature its grids span). Unlike a density its gradient entry and Jacobian
column are analytic, through the Doppler weights. ``penalties()["doppler"]`` is
its prior, which ``objective()`` adds.  The field is edited as any record field
(``project.set("D0", doppler="free_gas 300 K fit grid_max 900")``).

A fitted sample thickness (``areal_density n atoms/b fit [prior n sigma]`` on a
``target`` record of kind ``sample``, named by ``transmission`` or
``neutron_yield`` data sets, :doc:`segments`) is a parameter
``<target>.areal_density`` of kind ``"thickness"`` (``ParameterSet.thicknesses``;
:attr:`Parameter.segment_key` is the first data set on the target,
``prior_centre`` / ``prior_sigma`` its prior, ``sample_kind`` ``"transmission"``
or ``"yield"``); data sets that name one target share the parameter. Its
gradient entry and Jacobian column are analytic. ``penalties()["thickness"]`` is
its prior, which ``objective()`` adds.

A prior parameter covariance (``priorCovariance=<file>`` in ``<config>``,
:doc:`chi_squared`, "Prior parameter covariance") is set up by the engine when the
session starts (a refused file makes ``Session()`` raise; listed parameters fixed in
the model are marginalized). ``penalties()["prior_covariance"]`` is its term
:math:`(\mathbf{p}-\boldsymbol{\mu})^T C^{-1}(\mathbf{p}-\boldsymbol{\mu})`, which
``objective()`` adds (it is not part of ``calculate_chi2_rwa``);
``prior_covariance_residuals(x)`` returns the whitened residuals
:math:`L^{-1}(\mathbf{p}-\boldsymbol{\mu})` (:math:`C = LL^T`) in the order of
``prior_covariance_names``, and ``residual_jacobian`` appends them as its last rows,
with the constant Jacobian :math:`L^{-1}` on the prior's columns -- so
:math:`\sum r^2` of ``residual_jacobian`` is ``calculate_chi2_rwa`` plus the prior
term. ``Project.config()["priorCovariance"]`` reads the key and
``Project.set_config(priorCovariance=...)`` sets it (``None`` removes it).
:mod:`pyazr.priors` needs no engine: ``prior_from_fit(sav_path, covariance_path,
out_path, names=None, columns=None)`` writes the file from a fit's ``param.sav``
(means, ``#parametrization`` tag) and ``covariance.dat`` (the free level energies and
widths; ``names`` picks a subset, i.e. the marginal block; ``columns`` names the
columns of a ``covariance.dat`` from before its ``#columns`` header),
``write_prior(path, names, mean, cov, basis=None)`` any prior (checked symmetric and
positive definite), ``read_prior`` and ``read_covariance_dat`` read them back.

.. warning::

   These return the **data** χ² only. AZURE2's own fit objective additionally
   penalises every free normalization and energy shift. Minimise the residuals
   alone and the normalizations will drift to absorb every discrepancy,
   reaching a χ² that AZURE2 would never find. Append the penalty rows — see
   :doc:`chi_squared` and ``pyazr/examples/`` for the recipe.

The objective and the evaluation
--------------------------------

Two objects hold what a fit and a calculation need, every number from the engine.

``s.fit_objective()`` returns a :class:`pyazr.Objective`: the objective the program's fit minimises (the data χ²
and every prior and constraint) as one callable on the free vector ``x``, in the basis the session fits in.
``f(x)`` is ``s.objective(x)`` bit for bit; ``f.grad(x)`` and ``f.value_and_grad(x)`` its gradient;
``f.residuals(x)`` and ``f.jacobian(x)`` the objective as least-squares rows, ``sum(f.residuals(x)**2) == f(x)``
(the data rows, then one row per prior and constraint); ``f.log_prob(x)`` is ``-0.5 * f(x)``, ``-inf`` where the
program's sampler puts zero probability (outside ``f.bounds``, the fit's limits of a parameter of any kind, a
``<parameterSettings>`` block's among them, and where Park's overlap J of a level is negative); ``f.x0``,
``f.names``, ``f.bounds`` and ``f.errordef`` (1) for the minimiser.  Calls from several threads on one session run one after the
other.

``s.evaluate(x)`` returns a :class:`pyazr.Evaluation` from one engine evaluation: ``objective``; ``chi2``, the data
χ² (``Total-Chi-Squared`` of ``chiSquared.out``); ``terms``, each term of the engine's objective by its name;
``data[i]`` per data set with ``model``, ``residuals``, ``chi2`` (its row of ``chiSquared.out``), ``norm``,
``energies``, ``angles``, ``values``, ``errors`` and ``extra`` (the lab energy and angle, the S-factor conversion,
the E1 and E2 parts of a capture, Legendre coefficients: ``pyazr.kinds()["data_outputs"]``); ``model_data``;
``table()``.  ``ev.objective - ev.chi2`` are the priors: on the fitted T2 the objective is 297.86, the data χ² 293.94.

.. code-block:: python

   s = pyazr.Session("T2.azr")
   f = s.fit_objective()

   # scipy: the gradient, or the rows and their Jacobian (pyazr/examples/fit_scipy.py)
   r = scipy.optimize.minimize(f.value_and_grad, f.x0, jac=True, bounds=f.bounds, method="L-BFGS-B")
   q = scipy.optimize.least_squares(f.residuals, f.x0, jac=f.jacobian, x_scale="jac")

   # iminuit: MIGRAD and HESSE (pyazr/examples/fit_iminuit.py)
   m = iminuit.Minuit(f, f.x0, grad=f.grad, name=f.names)
   m.limits = f.bounds
   m.migrad(); m.hesse()

   # emcee: log_prob; with a pool, one session per worker (pyazr/examples/fit_emcee.py)
   sampler = emcee.EnsembleSampler(nwalkers, len(f.x0), f.log_prob)

   ev = s.evaluate(np.array(m.values))
   print(ev.terms)
   print(ev.table(0))

The parameters of a project differ in scale by many orders (an areal density in atoms/cm² beside a width in
MeV\ :sup:`1/2`): ``fit_scipy.py`` runs L-BFGS-B on ``x`` divided by ``|x0|``.  The three recipes are tested
(``tests/pyazr/fit_recipes_test.py``): on T2, iminuit's MIGRAD and both scipy recipes reach the program's own
MIGRAD minimum (iminuit to 1e-6; the scipy recipes go 3e-6 below it, where the program's
MIGRAD stops early).  For walls at every bound in a sampler, wrap ``log_prob`` in one line:
``lambda x: f.log_prob(x) if all(lo <= v <= hi for v, (lo, hi) in zip(x, f.bounds)) else -np.inf``.

Editing the model
-----------------

pyazr reads project format 3 only: a session or :meth:`pyazr.Project.open` on a
format-1 or format-2 file raises :class:`pyazr.FormatError` with the engine's
message, which names the conversion (``AZURE2 --convert OLD NEW``, or
``pyazr.convert_azr(old, new)`` from Python).

:class:`pyazr.Project` edits a project the way the engine reads it: records by
name and field, the ``<levels>`` rows by level and channel, ``<config>`` and
``<thm>`` by key.  The engine owns every grammar (it makes each edited line and
refuses what a run would refuse); every edit is validated, a new error raises
:class:`~pyazr.FormatError` and leaves the project unchanged, and the lines an
edit does not touch keep their bytes.  Nothing is written until ``save``:

.. code-block:: python

   from pyazr import Project, Nuisance

   p = Project.open("13N.azr")
   p.set("meyer_84", norm=Nuisance(1.0, "fit", prior=(1.0, "5%")))   # a record field
   p.set_level(("3/2-", 3.5), energy=3.5021)                          # a level by J^pi and E_x (MeV)
   lv = p.level(("1/2+", 2.365))
   c = lv.channels[0]
   p.set_channel(lv, c.pair, c.l, c.s, width=35000.0, fixed=False)    # a channel by pair key, l, s
   p.set_pair(1, radius=4.5)                                          # every row of pair 1 (fm)
   p.remove_level(("1/2+", 20))                                       # a background pole
   p.add_level(1.5, +1, 8.6, [dict(pair=1, l=2, s=0.5, width=1000.0)])
   p.add_dataset("roughton", reaction=p.reaction(1, -1), quantity="sigma", unit="b",
                 file="data/roughton.dat", columns="E_lab[MeV] _ value error", norm="1 fixed")
   p.set_config(output_dir="output_variant/")
   path = p.save("_variant.azr")                                      # the original stays as it was
   with p.session() as m:                                             # a calculation of the edited project
       print(m.calculate_chi2_rwa(m.params_rwa))

A level is named by its label (``"3/2- 3.5032"``, as :meth:`~pyazr.Project.levels`
gives it), by ``(J^pi, E_x)`` within 1 keV, or by the :class:`~pyazr.ProjectLevel`
itself; ``levels()`` gives each level's rows, its compound nucleus and the
engine's (group, level) indices, ``pairs()`` the particle pairs (masses, charges,
spins, radius, excitations).

Because AZURE2 reads its model from the file, an edited scheme takes effect by
opening a fresh session on the saved file.

One rule catches people out: **all levels of one J**\ :sup:`π` **share a
channel set**. Adding a level to an existing J\ :sup:`π` clones that group's
channel structure (channels not given get width 0, fixed); only a brand-new
J\ :sup:`π` lets your channel list define the group.

The other edits: ``add`` / ``set`` / ``remove`` / ``rename`` of any record (data
sets, extrapolations, normalizations, energy scales, targets, beams, detectors,
model data, potentials; ``rename`` follows every reference and the parameter
names of ``<parameterSettings>``); ``config()`` / ``set_config(**entries)``;
``thm()``, ``set_thm(**options)``, ``set_thm_experiment(name, **keys)`` and
``remove_thm_experiment(name)`` for the ``<thm>`` section (a data set joins an
experiment by its own field ``thm <experiment>``).  The THM keys are the
engine's (:doc:`../theory/thm_implementation`): ``cbackground="1/2+:2"`` adds a
coherent THM-only background amplitude whose Re and Im are fit parameters of
kind ``"cbkg"`` (``parameters.cbkg``, written back into the key by
``save_fit``); ``ps="hulthen:0-40"`` averages the HOES cross section over the
accepted spectator directions; ``distortion="coulomb"`` multiplies the model by
the zero-range DWBA distortion factor :math:`R(E)`, with ``opticalAA``,
``opticalSF``, ``spectatorAngle``, ``distortionRef``, ``distortionRatio`` and
``boundState``; ``theta="50-70"`` turns the model into the fixed-angle
observable; ``distortion="optical", vertexModel="dw"`` replaces the plane-wave
entrance vertex by the DWBA surface term; ``spectatorAngles="cm:135-180"``
averages over the accepted spectator directions.  A session reports what they
do: ``thm_background(name)``, ``thm_lineshape(name, energies)``,
``thm_vertex(name, energies, strict=False)`` (its ``reached`` entry is False at
an energy the window or the DW vertex grid does not reach; such energies raise a
``UserWarning``, with ``strict=True`` a ``ValueError``) and
``thm_distortion(name, energies)``.  A ``<thm>`` edit is checked as the engine
reads a format-3 block; a refusal raises :class:`~pyazr.FormatError` with the
engine's message and leaves the project unchanged.

To snapshot a fit, use :meth:`~pyazr.Session.save_fit`:

.. code-block:: python

   azr, sav = m.save_fit("7Be_fit.azr")        # or save_fit(path, x_best)

``save_fit(path, x=None, param_sav=True, verify=True, norms="fitted",
close_session=False)`` converts the fit to the physical values a
``<levels>`` row holds and writes them into the rows they come from (fixed
widths included, at the physical value the fit had: a fixed Brune amplitude
keeps its amplitude, not its width; a level by the engine's (group, level), a
channel by its pair key, l and s -- :func:`pyazr.snapshot.place_level_values`),
writes the fitted free nuisances (normalizations, energy-scale terms, fitted
densities, temperatures, thicknesses, abundances) into the records they come
from, writes the companion ``param.sav`` carrying every parameter, and then
reopens what it wrote and checks it against the fit (``verify``: the
physical values transformed back, the reduced-width amplitudes the file reads
for the particle and beta-decay channels of each level against the fit's, up
to a common sign of the level, and with ``norms="fitted"`` the data
:math:`\chi^2` and the prior terms of the snapshot run alone, to
:math:`10^{-6}`) --- removing both files and raising if they disagree.  (The
amplitude check is there because the physical one cannot see a value two
parameters map onto: a negative beta-decay feeding amplitude came back
positive from the physical transform until 2026-10-12.)  A coherent THM
background is written back into its experiment's ``cbackground=`` key.  A
format-3 nuisance states its prior with its own centre (``1.0 fit prior 1.0
5%``), so the snapshot's fitted value leaves the prior where the fit had it, and
a fit started from the snapshot is pulled to the same centres, not to the fitted
values; ``save_fit(path, x, norms="nominal")`` keeps the loaded file's values in
the fields and the fitted ones in the ``.sav`` only (run with it as the parameter
file to reproduce the fit; :meth:`~pyazr.Session.load_sav` matches its
names to the run's as AZURE2 does: names the project lacks warn with
``W_PARAMETER_FILE_NAMES``, a file with no name of the run is refused with
``E_PARAMETER_FILE_NO_MATCH``).  See
:ref:`nominal-norm`.  The check opens a second engine beside the fit's; as
the last call of a fit, ``save_fit(path, x, close_session=True)`` closes the
fit's session first, so that only one is in memory (a 19F THM session is
about 220 MB, 390 MB with the verifying one beside it).  ``path`` is always
explicit; nothing is written in place.  Note a ``.azr`` names its data files
relative to itself, so a snapshot runs from the directory the original did.

.. warning::

   Adding or removing data segments changes which energies AZURE2 evaluates.
   The external-capture integrals are cached in ``output/intEC.dat`` /
   ``output/intEC.extrap`` for the energies of one grid.  AZURE2 records a
   signature of that grid beside each file (``.sig``) and recomputes, with a
   warning, a cache that no longer matches -- but only at session start.
   After any data edit, delete those caches (or give the edited model its own
   output directory) so the integrals are recomputed;
   ``Session.recalculate_external_capture()`` forces it inside a live session.
   See ``pyazr/examples/edit_model.py``.

Parameter and output files
--------------------------

Parameter files (``param.par``, ``param.sav``, ``param.fit``) are read by
**name**, as AZURE2 reads them: ``Session.read_sav(path)`` returns
``{name: value}`` (the ``#parametrization`` line and comments skipped),
``full_rwa_from_sav(path)`` applies a file to the session's full vector
(a name the file lacks keeps the session's value; a file with no name in
common is refused), ``update_rwa_params_from_sav()`` and
``update_sav_from_rwa_params(x)`` (writes ``param.sav.new``, rows updated in
place, missing free rows appended) do the same with the run's own
``param.sav``, and ``pyazr.bands.best_fit_params`` loads a fit for the
bands. Do not read these files by row position: rows are added as the code
grows (the :math:`\sqrt{E}` coefficients, the tag line), and a positional
read of a file from another layout is silently off by one or more rows.

``write_output_files(params=None)`` writes the CLI's output files
(``AZUREOut_*``, ``chiSquared.out``, ...) for a parameter vector, into
``output_dir``: the output directory of the project's ``<config>``, relative
to ``cwd``, as the engine resolves it (not necessarily ``output/``).

``pyazr.bands`` (``sensitivities``, ``uncertainty_bands``,
``extrapolation_bands``) spans the columns the CLI's ``covariance.dat`` spans:
level energies, reduced widths and THM ``cbkg`` values (kinds 0, 1 and 9);
norms and shifts only rescale or move the data and carry no band.

Masses and separation energies
------------------------------

``pyazr.masses`` reads the AME2020 atomic mass table that ships with pyazr
(``pyazr/data/mass_1.mas20``). Particles are
named by isotope (``"12C"``, ``"C12"``, ``"c-12"``) or shorthand (``"p"``,
``"n"``, ``"d"``, ``"t"``, ``"a"``, ``"helion"``):

.. code-block:: python

   from pyazr import masses, Project

   masses.atomic_mass("p")               # 1.007825031898 u, atomic
   masses.separation_energy("p", "12C")  # 1.943490064 MeV, from the mass excesses
   masses.pair_masses("a", "12C")        # Z1, M1, Z2, M2, S and the compound (8, 16)

   p = Project.open("13N.azr")
   for row in masses.mass_check(p):      # every pair against the table
       print(row.describe())
   masses.fill_masses(p)                 # masses and S of every pair from the table
   masses.fill_masses(p, pairs=[2], separation=False)   # masses of pair 2 only
   p.save("13N_ame.azr")

``fill_masses`` identifies each particle by its charge and its mass rounded
to a mass number, so mass numbers and a bare-proton mass are recognised. It
sets the masses and separation energy of the chosen pairs through
``Project.set_pair`` (the engine writes them on every row of the pair, and the
entrance-S copy on the rows whose entrance pair changed), and returns the list
of changes. A gamma pair gets the compound nucleus's mass and keeps its
S of 0. A beta pair gets its parent's mass. Moving S moves the pair's
threshold, so compare the fit before and after.

Several compound nuclei (natural targets)
-----------------------------------------

A natural or compound target is an abundance-weighted sum of cross sections
of different compound nuclei (see the Particle Pairs page of the user guide).
In format 3 the target is a ``target`` record with one ``isotope`` field per
target nucleus (``isotope 12C 0.6 fit prior 0.6 0.02``) and an ``isotope_sum``
field; a data set on it names the sum of its reactions (``12C(p,p)12C +
13C(p,p)13C``), each term weighted by its isotope's abundance.
``merge_projects`` builds one such project from single-nucleus ones: the
engine merges their ``<levels>`` rows (each member's pair keys and level IDs
moved past those before it), every record is kept (a name an earlier member
uses is renamed ``<name>_<k>`` with its references), data-file paths are
rewritten relative to the merged file, ``isotopes`` adds the target record, and
a *shared* group of data sets (the same natural-target data set present in
every member) becomes one record whose reaction is the sum of the members':

.. code-block:: python

   from pyazr import merge_projects, Session

   report = merge_projects(
       ["20Ne/17O+a.azr", "22Ne/22Ne+n.azr", "21Ne/21Ne+n.azr"], "natNe.azr",
       isotopes={"20Ne": 0.9048, "22Ne": 0.0925, "21Ne": 0.0027},
       shared=[{0: "vaughn", 1: "vaughn", 2: "vaughn"}],   # one dict per common data set: member -> record name
   )
   print(report.pair_map, report.record_map, report.warnings)

   m = Session("natNe.azr")
   m.nuclei                     # the compound nuclei and their abundances
   m.parameters.abundances      # the abundance parameters (<target>.isotope.<nucleus>), kind "abundance"
   m.penalties()["abundance"]   # the prior and sum-to-one terms at the current point

Member 0's ``<config>``, ``<thm>``, ``<parameterSettings>`` and ``<mcmc>`` are
kept, the others' are not (a warning says so).  On the file side the isotopes
are record fields (``project["natural"]["isotope"]``, ``project.set("natural",
isotope=[...])``) and :attr:`ProjectLevel.nucleus` gives each level's compound
nucleus; ``Pair.nucleus`` and ``Parameter.nucleus`` identify the nucleus at run
time; ``save_fit`` writes a fitted abundance back into its isotope field.

Scattering lengths
------------------

``scattering_lengths(params)`` returns the s-wave scattering length of every
channel spin of every particle pair at a free RWA vector
(:doc:`../theory/scattering_length`), so the posterior of :math:`a_0` follows
from an MCMC chain sample by sample:

.. code-block:: python

   m = Session("7Be+p.azr")
   for v in m.scattering_lengths():          # current parameters
       print(v.pair, v.jpi, v.a0, v.has_levels)
   a02 = [next(v.a0.real for v in m.scattering_lengths(x) if v.pair == 1 and v.jpi == "2-")
          for x in chain]                     # chain: free RWA vectors

``a0`` is complex when a lower channel is open at the pair's threshold.  The
measured lengths -- ``datum`` records of quantity ``scattering_length`` in
``<modelData>`` -- are already inside
``calculate_chi2_rwa`` (and its gradient and ``residual_jacobian``, one row per
constrained part after the data rows); ``penalties()["scattering_length"]``
reports their sum, and ``scattering_length_chi2(x)`` the same from the
model-data term.  On the file side they are records: ``add_datum``,
``set`` and ``remove`` as for any record.

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

   # fetch one dataset; an S factor becomes a format-3 data set of quantity sfactor
   data = nds.fetch_exfor("O2599004")                 # S-factor, B*KEV
   record = data.to_sfactor_dataset("run/data", "S_luna", norm="1 fit prior 1 5%")
   # 'dataset S_luna\n    reaction 13C(p,g)14N\n    quantity sfactor\n    unit MeV*b ...'
   # a cross section: an AZURE2 data file, then a data set record
   xs = nds.fetch_exfor(entry_id)                      # an entry in B, MB, ... (not SFC)
   kw = xs.to_azr("run/data", entrance=1, exit=2)
   p = Project.open("run/13N.azr")
   p.add_dataset("xs", reaction=p.reaction(1, -1), quantity="sigma", unit="b", norm="1 fixed",
                 file=os.path.relpath(kw["data_file"], p.directory), columns="E_lab[MeV] _ value error")
   p.save()

   # without the network, or when the IAEA answers with a challenge page:
   # the same dataset from a local copy of the EXFOR master files
   data = nds.fetch_exfor("O2599004", source="local")

   # level scheme of the compound nucleus, from ENSDF (IAEA LiveChart, or the
   # NNDC ENSDF API when LiveChart cannot be reached; source="nndc" forces it)
   for lv in nds.fetch_levels("14n"):
       print(lv.energy_mev, lv.jp, lv.half_life)

   # the paper behind the dataset, resolved to a DOI
   ref = nds.reference("O2599004")
   doi = nds.resolve_doi(ref)   # -> 10.1103/physrevlett.131.162701

Frames and units are handled automatically: ``EN-CM`` energies are converted
to lab, ``NB/SR`` differentials to b/sr.  S-factor data (``B*KEV``, an ``SFC``
code) are format-3 quantity ``sfactor``: ``to_sfactor_dataset`` writes them as S
(MeV b, the energy as EXFOR gives it, ``E_cm`` for ``EN-CM``), and the engine
multiplies its model by :math:`E e^{2\pi\eta}` at each point with the masses of
``<levels>``; ``to_azr`` refuses them (until November 2026 it converted them to
barns with mass numbers, which moves the cross section by +6.7 % for
p + :sup:`12`\ C at 100 keV).  ``sfactor_to_cross_section(S, E_cm_kev, z1, z2, m1_u,
m2_u)`` converts by hand with the engine's function and the atomic masses
(``_azure2.sfactor_conversion`` gives the factor).  Pick ``observable`` to match the data's angle frame
(EXFOR ``ANG-CM`` → ``"differential-cm"`` or ``"analyzing-power"``; lab
``ANG`` → ``"differential"``).  See the ``azure2-setup`` skill for the full
quantity/unit reference and failure modes.

Ratio-to-Rutherford data (``,DA,,RTH``; ``,DA,,MOT`` for identical
particles) are a quantity of project format 3, ``ratio_to_rutherford``
(``ratio_to_mott``), which the engine divides by its own Coulomb cross section
at the point's energy, folded as the data are
(:doc:`../examples/ratio_to_rutherford`).  ``format3_quantity()`` names the
quantity of an entry's code, and ``to_dataset`` writes the ratios as they are
and returns the ``dataset`` record::

   data = nds.fetch_exfor("T0001004", source="local")   # 12C(p,el),,DA,,RTH
   data.format3_quantity()                               # 'ratio_to_rutherford'
   record = data.to_dataset("run/data", "pC_rth", norm="1 fit prior 1 3%")
   # 'dataset pC_rth\n    reaction 12C(p,p)12C\n    quantity ratio_to_rutherford ...'

``to_azr`` refuses such data: until November 2026 it multiplied them by a
Rutherford cross section, which an energy shift did not move and a fold did not
average, with mass numbers.  ``RTH`` of identical particles is refused until
the entry's text says the reference (pass ``quantity="ratio_to_mott"``), and
``MOT`` of different particles (the relativistic point-nucleus cross section)
is refused.  A dimensionless differential whose reaction code was dropped
(``x4get`` drops the suffix that says ``RTH``) is taken for an analyzing power
with a warning; ``rutherford=False`` says it is one.

The module needs nothing beyond NumPy and the standard library, and doubles as
a command-line tool::

   python -m pyazr.nds search --target C-13 --reaction p,g --quantity SIG
   python -m pyazr.nds download O2599004 -o data/skowronski.dat
   python -m pyazr.nds reference O2599004

The IAEA services are sometimes unreachable, or answer every request with a
Cloudflare challenge page. ``fetch_exfor`` and ``reference`` then read the
entry from a local copy of the EXFOR library, with a warning. The copy is a
clone of `IAEA-NDS/exfor_master <https://github.com/IAEA-NDS/exfor_master>`_
(about 1.5 GB; ``git clone --depth 1``, refreshed with ``git pull``) at
``$AZURE2_EXFOR_MASTER`` or ``~/.azure2/exfor_master``.
``nds.exfor_master_path()`` says which copy is used. ``source="web"`` or
``source="local"`` uses only one of the two, and ``python -m pyazr.nds
download <ID> --source local`` does the same from the command line. The local
copy gives the original units only (``plus=0``). A subentry holding several
datasets needs the pointer, as in ``"C1010003.2"``. ``search_exfor`` needs the
Web-API.

The Qt setup utility of the released versions has its own EXFOR reader
(``gui/src/ExforData.cpp`` on ``dev``), an independent implementation of the
same Web-API.

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
one ``Session()`` session on it, so pass every energy in a single call rather
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

Angular-distribution coefficients
---------------------------------

Legendre coefficients of a reaction's c.m. angular distribution are a data set
of quantity ``legendre``: data rows ``E_lab  k  a_k  da_k`` with the order
:math:`k` in the angle column; ``expansion a_over_a0`` (unit ``1``) gives
:math:`a_k = A_k/A_0` of :math:`d\sigma/d\Omega_\text{cm} = A_0[1 + \sum_k a_k P_k]`
(the coefficients a ``legendre`` extrapolation computes), ``expansion absolute``
(unit ``b/sr``) the absolute :math:`A_k`. The model's coefficient of each row's
order comes back from :meth:`~pyazr.Session.calculate_rwa` in the data slot, so
:math:`\chi^2`, fitting and the analytic gradient need no special handling;
target integration folds the absolute :math:`A_k` before the ratio. Normalised
coefficients cannot carry a free norm (the engine refuses it); an elastic
reaction, total capture and sums are refused.

.. code-block:: python

   p.add_dataset("wang_ak", reaction=p.reaction(2, 1), quantity="legendre", unit="1", expansion="a_over_a0",
                 file="data/wang_ak.dat", columns="E_lab[MeV] order value error", energy_window="0.8 1.9 MeV",
                 norm="1 fixed")

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
units are the corresponding sanity check for capture. AZURE2's Wigner limit
(``Parameter.wigner_limit``, ``ChannelWidth.wigner_gamma2``) is the
Teichmann-Wigner sum-rule unit 3ħ²/(2μa²); a θ² quoted against the other
common unit, ħ²/(μa²) (``lane_thomas_wigner()``), is 1.5 times larger.

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

   azr.cache_stats()   # input_files_read, input_files_served, input_files_held (and six
                       # Coulomb-memo counters of earlier builds, now always 0)

The engine's physics functions need no session: ``pyazr.physics`` has one function per row of the binding's physics
table (kinematics, Coulomb functions and phase, penetrability, shift function and its energy derivative, hard-sphere
phase, Wigner limit, single-particle width, photon penetrability, resonance quantities, Coulomb cross sections, S-factor
conversions), each on broadcast numpy arrays with its arguments, units and domain in its ``help()``; an element outside
the domain raises ``pyazr.DomainError``, or with ``errors="mask"`` gives NaN and a status per element.  The page
:doc:`physics_functions` lists them; ``pyazr.kinds()["physics"]`` gives the same rows as data.  Where the session's
own setting matters (its Coulomb table, GSL, a nuclear potential), ``azr.channel_functions(pair, e_cm, l)`` gives F,
G, their derivatives, P, S and the hard-sphere phase shift of one of its pairs.

The bound-channel Whittaker function needs no session: ``pyazr.whittaker_wave(l, eta_b, rho)`` returns
:math:`W_{-\eta_b,l+1/2}(2\rho)` (unit ANC, :math:`W \sim e^{-\rho}(2\rho)^{-\eta_b}`) on broadcast arrays of the
dimensionless :math:`l`, :math:`\eta_b = Z_1Z_2e^2\mu/(\hbar^2\kappa)` and :math:`\rho = \kappa r`, as
``mantissa`` and ``log10_scale`` (the value leaves the double range just below a charged threshold), ``log10_W``,
``W`` and ``rho_dlog`` (:math:`\rho W'/W`, the shift function of a closed channel at :math:`\rho = \kappa a`). It
is the engine's GSL-free function (physics/coulomb), good to :math:`10^{-12}` against mpmath for
:math:`l \le 10`, :math:`\eta_b \le 1000`, :math:`10^{-4} \le \rho \le 200`; an invalid argument raises
``ValueError`` naming the element. ``pyazr.whittaker`` (pair, radius and binding energy in MeV) is the older Python
copy of the engine's GSL path.

.. code-block:: python

   w = pyazr.whittaker_wave(1, eta_b=2.3, rho=np.linspace(0.1, 10, 100))
   w["W"], w["rho_dlog"]

Stopping powers need no session either: ``pyazr.stopping`` is AZURE2's stopping library, SRIM-2013's stopping
powers ported from LunaSRIM (cite SRIM: ``pyazr.stopping.CITATION``), the same numbers a target with ``stopping
evaluation`` uses.  Energies are lab MeV of the projectile, stopping cross sections eV/(10\ :sup:`15` atoms/cm\ :sup:`2`)
per atom of the material, areal densities atoms/cm\ :sup:`2` of the material, widths MeV:

.. code-block:: python

   from pyazr import stopping
   ta2o5 = stopping.Material("Ta2O5")                    # solid; phase="gas" for a gaseous compound
   eps = stopping.stopping(ta2o5, "p", [0.2, 0.5, 1.0])  # dict: electronic, nuclear, total, rel_uncertainty
   per_o = eps["total"] / ta2o5.atom_fraction(8)         # per O atom: the engine's stopping per active 18O
   stopping.energy_after(ta2o5, "p", [1.0], [1e18])      # MeV after 1e18 atoms/cm2
   stopping.straggling(ta2o5, "a", [2.0], [1e18], beam_sigma=[0.001])["total"]
   stopping.range_of(ta2o5, (6, 12.0), [10.0])           # a 12C ion: (Z1, mass in u)
   print(stopping.provenance(ta2o5, "p"))                # SRIM-2013, the citation, the data per element

The same functions are rows of a table (``pyazr.stopping.table()``, also ``pyazr.kinds()["stopping"]``) in the shape
of the physics rows, each with SRIM's citation; ``stopping.evaluate(name, formula, **arguments)`` evaluates a row on
broadcast arrays with a status per element (NaN where the status is not 0), bit for bit the functions above:

.. code-block:: python

   result, status = stopping.evaluate("stopping_power", "Ta2O5", z1=1, mass_u=1.00727646688,
                                      energy=np.geomspace(0.001, 10, 200))
   result, status = stopping.evaluate("energy_after", "Ta2O5", z1=1, mass_u=1.00727646688, energy_in=0.1,
                                      areal_density=[1e17, 1e21])   # status [0, -2]: beyond the range

A rebuild of the model (``set_channel_radius``, a nuclear-potential change, a mode switch) reads no input
file that has not changed on disk (the ``input_files_*`` counters, for the whole process) and does not
print the setup ``Note:`` lines again: a session prints each once. Before
AZURE2 dev 2026-11-04 a likelihood that sampled a channel radius grew the memo by about 100 kB per new
radius without bound, re-read the ``.azr`` and every data file at each call, and printed every Note at
each call.

The Coulomb functions follow the run's own configuration, so the same call
returns the accurate routine's values, GSL's (``--gsl-coul``), or the Numerov
solution through a nuclear potential (the hybrid model of the
``<potential>`` block). Comparing them is how one sees what those options do to
the external region.

The hybrid values come from a Numerov integration of H⁺ = G + iF inward from
the 15 fm matching radius on a grid of about 0.0275 fm, laid so that the
requested radius is a grid node; the functions and their derivatives are good
to about 10\ :sup:`−8` there (with the potential switched to zero depth, P and S
agree with the Coulomb functions to 6·10\ :sup:`−9` and 2·10\ :sup:`−8` for p + :sup:`12`\ C).
The Coulomb term of the potential is that of a uniformly charged sphere of radius
:math:`R_C` (the record's ``coulomb_radius``; 3.5 fm in an old ``<potential>`` block) with e\ :sup:`2` = αħc, the constants of the Coulomb functions it is matched to;
builds before 2026-11-10 used e\ :sup:`2` = 1.4399764 MeV fm, which put a 10\ :sup:`−5`
difference into P and S that was not the integration's.
Before October 2026 the values were interpolated linearly between nodes and the
derivatives were the slope of the bracketing cell, a first-order quantity: the
shift function of a hybrid pair was off by 4·10\ :sup:`−3` (l = 0) to
4·10\ :sup:`−2` (l = 3) for p + :sup:`12`\ C at 5.4 fm, depending on where the
radius fell between nodes. A potential of vanishing depth now reproduces the
plain Coulomb functions (``tests/pyazr/hybrid_derivative_test.py``).

Below threshold the potential acts as well. A closed channel or a bound state of a
hybrid pair uses the decaying solution :math:`\tilde W_l(r)` of the same radial
equation at :math:`E = -B`, equal to the Whittaker function
:math:`W_{-\eta_b,l+1/2}(2\kappa r)` beyond the 15 fm matching radius and
integrated inward from there by the same Numerov grid. It gives the shift function
:math:`S = a\tilde W'(a)/\tilde W(a)` (level-parameter transformations, boundary
values, closed channels), the ANC-to-reduced-width conversion
:math:`\mu a/(\hbar^2 \tilde W(a)^2)` (the ANC keeps its meaning: the bound state
is :math:`C\,W_{-\eta_b,l+1/2}(2\kappa r)` at large r), and the final-state tail
and NF integral of external capture. At threshold the zero-energy decaying solution
:math:`\sqrt{r}K_{2l+1}(\sqrt{8cr})` (:math:`r^{-l}` without Coulomb) is the outer
solution, closer to threshold than :math:`\eta_b = 1000` S is interpolated linearly
to it, and in the band :math:`\eta > 100` just above threshold
:math:`S(E) = 2S(0) - S(-E)`, as for the Coulomb functions; S is then continuous
through threshold (to within :math:`|dS/dE|\,|E|`, about :math:`10^{-3}` at 1 keV for
p + :sup:`12`\ C). Against an independent integration (scipy DOP853 with mpmath's
Whittaker functions) S agrees to :math:`10^{-7}`, :math:`\ln\tilde W(a)` to
:math:`10^{-9}` and :math:`\int_a^\infty \tilde W^2\,dr` to :math:`3\cdot10^{-8}`
(``physics_hybrid_potential``). Builds before 2026-11-10 used the undistorted
Whittaker function below threshold, so S of a hybrid pair jumped at threshold (by
about 1 for p + :sup:`12`\ C with V₀ = 20 MeV) and a bound level was transformed
with the Coulomb S while the scattering used the distorted one.
``coulomb_functions`` returns this S at :math:`E \le 0` (channel radius only).

The hybrid nuclear potential, per particle pair
-----------------------------------------------

A nuclear potential belongs to a **particle pair**: it bends the radial wave
functions of that channel and no other. Each pair therefore carries its own
setting, and ``pair=0`` is the default that a pair without one inherits.

.. code-block:: python

   azr.nuclear_potentials()          # {0: default, 1: ..., 2: ...}, resolved
   azr.nuclear_potential(1)          # one pair; .own says if it is its own

   azr.set_nuclear_potential(1, type="WoodsSaxon", enabled=True,
                             V0=20.0, R=3.6, a=0.6)
   azr.set_nuclear_potential(2, type="Gaussian", enabled=True, V0=60.0, r0=4.0,
                             coulomb_radius=4.0)
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
  caveat as :meth:`~pyazr.Session.set_channel_radius`, for the same reason.
- **A fit made without the potential is not a fit made with it.** Refit before
  reading anything off the model.

``coulomb_radius`` (fm, default 3.5) is the radius :math:`R_C` of the uniformly
charged sphere of the potential's Coulomb term; ``nuclear_potential(p).coulomb_radius``
reads it back. ``matching_radius`` (fm) is where the solution inside the
potential is matched to the Coulomb functions and ``step`` (fm) the radial step
of the integration (15 fm and 0.0275 fm unless the record states others);
``set_nuclear_potential(1, matching_radius=12.0, step=0.02)`` sets them and
``nuclear_potential(p).matching_radius`` / ``.step`` read them back.

In a project of format 3 a potential is a ``potential`` record of the
``<potentials>`` section, one per particle pair; a pair without a record has
none, and there is no default record:

.. code-block:: text

   <potentials>
   potential p12C
       pair            p+12C
       coulomb_radius  3.5 fm
       central         woods_saxon depth 20 MeV radius 3.6 fm diffuseness 0.6 fm
       matching_radius 15 fm
       step            0.0275 fm
   </potentials>

``central`` is ``woods_saxon depth radius diffuseness`` or ``gaussian depth
width``; the depth is positive for an attractive well. The pair is named by its
label, with ``pair_key <key>`` only where two pairs of ``<levels>`` have that label.
``coulomb_radius`` is required; ``matching_radius`` and ``step`` default to 15 fm
and 0.0275 fm (the run log notes it). Lengths and energies may be written in any
unit of the unit table (``keV``, ``m``). The fields, their diagnostics and the
refusals (a photon or beta pair, two records for one pair, ``relativistic=1``,
the entrance pair of THM data) are in :doc:`../reference/format3`.
A pair with one spinning particle can have a spin–orbit term (FRESCO's
convention, operator :math:`2\,\mathbf{l}\cdot\mathbf{i}`; :doc:`../theory/hybrid_potential`):
the record's ``spin_orbit woods_saxon depth 6 MeV radius 3.6 fm diffuseness 0.6 fm on light operator 2l.s``,
or in pyazr

.. code-block:: python

   azr.set_nuclear_potential(1, spin_orbit={"depth": 6.0, "radius": 3.6,
                                            "diffuseness": 0.6, "on": "light"})
   azr.nuclear_potential(1).spin_orbit        # the dict, or None
   azr.set_nuclear_potential(1, spin_orbit=False)   # remove it
   azr.coulomb_functions(1, energies, L=1, J=1.5)   # J required for such a pair

The channel functions of such a pair depend on J: ``coulomb_functions`` needs
``J`` for it (``E_SO_NEEDS_J``) and refuses one for any other pair
(``E_SO_NO_TERM``).  A term on a pair whose two particles both have spin
(``E_SO_TWO_SPINS``), on a spin-0 particle (``E_SO_ZERO_SPIN``) or on the
default ``pair=0`` is refused with a ``ValueError``.

``save_fit`` of a format-3 project writes the session's potentials as these
records when ``set_nuclear_potential`` changed them (the converter's writer: a
session default becomes one record per particle pair it applies to).  An older
project's ``<potential>`` block is not written: ``save_fit`` of a format-1/2
project whose potential the session changed raises ``ValueError`` with the
conversion command (``pyazr.convert_azr(old, new)``).

An older project has the ``<potential>`` block instead, where keys before the
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

Its Coulomb radius is 3.5 fm. ``AZURE2 --convert OLD NEW`` writes one
``potential`` record per particle pair the block gives a potential (the default
expanded, ``coulomb_radius 3.5 fm`` written out; a section with
``useHybridPotential=0`` writes nothing and is noted), and a format-3 project
that still holds the block is refused (``E_OLD_POTENTIAL_BLOCK``). The engine,
pyazr and the converter read the block through one parser, so they cannot drift.
Worked example: ``pyazr/examples/nuclear_potential.py``; ``tests/hybrid_potential``
pins the behaviour.

The potential belongs to the session. Several sessions in one process keep
their own default and per-pair settings: a second session's ``<potential>``
block, or ``set_nuclear_potential`` on it, leaves the first session computing
with the settings it was opened with. (Before October 2026 the manager behind
the settings was one process-wide object, so the last session to read a
``<potential>`` block silently reconfigured every other live session;
``tests/pyazr/potential_sessions_test.py`` pins the isolation.)

The Coulomb functions are computed at every call: there is no process-wide
memo (``docs/upgrade/coulomb.md``). Earlier builds kept one, which took
energies within :math:`10^{-12}` MeV as equal, so that a second session in one
process could differ from a fresh process in the last digits; now every session
gives the bits of a fresh process. In its place the channel functions of the
evaluation energies are computed once per pair, :math:`l` and exact channel
energy, in a table built before the point loops, and the radial Coulomb waves of
the external-capture integrals once per exact radius in a per-thread scratch
(``docs/upgrade/channel_table.md``); neither depends on what was computed
before. A fit that moves a free energy shift keeps the tables of the
measurement's last three shift states, so the gradient's steps back to an
earlier shift do not recompute them. The Coulomb table (``--build-coulomb-table``)
remains the way to make P, S and the phases fast.

A project as records (format 3)
-------------------------------

:class:`pyazr.Project` builds and edits a format-3 project as named records
(``docs/upgrade/project_format3.md`` section 4); the recipe from schema to
calculation is :doc:`project_from_python`.

.. code-block:: python

   from pyazr import Project, Nuisance, FormatError
   p = Project.open("T2.azr")            # a format-1/2 file: FormatError with the conversion command
   p["s02_diff090"].norm                 # 'D': the field's text; p["s02_diff090"]["norm"] too
   p["s01_int_thin"].set(norm="1.05 fit prior 1 5%")   # replaces that one line
   p.add_normalization("D2", value=Nuisance(1, "fit", prior=(1, "8%")))
   p.rename("D2", "D_late")              # every reference and parameter name follows
   p.remove("D_late")                    # refused while a record names it
   p.validate()                          # [] or the engine's diagnostics (AZURE2 --validate --json)
   print(p.describe("s09_ratio"))        # AZURE2 --describe of that record
   p.save("T2_edited.azr")
   with p.session() as m:                # pyazr.Session on the saved project
       print(sum(m.segment_chi2()))

``Project.new(path, nuclear_part=file)`` starts from another format-3 file's
``<config>``, ``<levels>`` and other sections without records.  There is one
``add_<kind>`` per record kind of the engine's table (``add_dataset``,
``add_normalization``, ``add_energy_scale``, ``add_extrapolation``,
``add_target``, ``add_beam``, ``add_detector``, ``add_datum``,
``add_potential``), generated at import from ``pyazr.schema()``:
``help(Project.add_dataset)`` lists the fields with their help text and units.
A value is the file's text, a number, a ``(value, "unit")`` tuple, a
:class:`pyazr.Nuisance` or a list (a field given more than once); ``#`` starts
the field's comment, as in the file.  Every edit is checked by the engine's
validator on the project as it would stand and refused with
:class:`pyazr.FormatError` (``.diagnostics``, ``.codes``: the engine's) when
it adds an error; the project is then unchanged.  The project is the engine's
lossless document: opened and saved without an edit it is written byte for
byte, and an edit changes only the lines of the fields it names.

``save_fit(path, x, norms="fitted")`` of a format-3 project writes the free
norms, shift terms (offset and square-root) and shift-function parameters into
the nuisances of their records (a shared normalization's or energy scale's
own), each in the unit it is written in, and the target isotopes' abundances;
the priors keep their stated centres, so the snapshot run alone gives the fit's
chi-squared and priors.

Model averaging
---------------

A result that depends on choices the data cannot make (channel radius,
plane- or distorted-wave THM vertex, optical potentials, spectator window) is
quoted from several fitted *variants*. ``pyazr.modelavg`` does the
bookkeeping of that average and nothing else: the fits are yours.

.. code-block:: python

   from pyazr import Variant, model_average, parameter_label, write_averaged_azr

   # after fitting one variant: physical values keyed by a file-stable name
   names  = [parameter_label(p, azr.pairs) for p in free_params]  # "G[2-#1;p1;L1;S1]"
   v = Variant("r5.1_dw", {"radius": 5.1, "vertexModel": "dw"},
               chi2=chi2, npoints=N, nfree=k, values=dict(zip(names, phys)),
               covariance=(names, cov_phys),            # or errors={name: sigma}
               derived={"wg(213)": (wg, wg_err)}, prior=1.0)

   avg = model_average(variants, weights="aic", rescale=None)  # "bic", "chi2", "flat"; "best"
   avg.weights                      # one row per variant: chi2, N, k, dIC, prior, weight
   a = avg.derived["wg(213)"]       # a.mean, a.stat, a.model, a.total
   print(avg.summary())
   write_averaged_azr("adopted.azr", "averaged.azr", avg)    # means into a copy

The weights are :math:`w_i \propto \pi_i\,e^{-\mathrm{IC}_i/2}` with
:math:`\mathrm{IC} = \chi^2 + 2k` (Akaike, the default), :math:`\chi^2 + k\ln N`
(BIC), :math:`\chi^2` alone, or 1 (flat), and :math:`\pi_i` an optional prior
per variant. ``rescale="best"`` divides every :math:`\chi^2` by
:math:`\max(1, \chi^2/\nu)` of the best variant first. Each quantity gets the
weighted mean, the **statistical** spread (square root of the weighted mean of
the variances), the **model** spread (square root of the weighted variance of
the means) and their quadrature sum, all kept separate. ``avg.covariance()``
gives the averaged matrices. Variants are matched by name; a quantity missing
from some variants is averaged over those that have it (``coverage`` is their
summed weight).

``parameter_label`` names a level energy ``E[<Jπ>#<n>]`` and a width
``G[<Jπ>#<n>;p<pair key>;L<l>;S<s>]``, with ``n`` the engine's level number
within the :math:`J^\pi` group, so the name depends on the ``.azr`` and not on
the radius or the THM options. ``write_averaged_azr`` sets the averaged energies
and widths of a template into its rows through :class:`pyazr.Project` (names it does not find are
skipped with a warning) and keeps everything else, so the result opens in the
GUI. It is a representative model, not a fit.

The module is pure Python (no numpy). ``scripts/thm_model_average.py`` drives a
whole THM grid: one project per variant along the axes ``--radius-pairs`` /
``--radii``, ``--vertex-model``, ``--vertex``, ``--optical``, ``--ps``,
``--lineshape`` and ``--distortion`` (or a JSON ``--spec``), each fitted in a
fresh Python process (``--in-process`` to keep them in one), one after the
other, with scipy's ``least_squares`` on ``residuals``/``residual_jacobian``
plus the penalty rows. ``--penalty-hook file.py:func`` adds rows of your own
(measured strengths, priors) to the fit, the :math:`\chi^2`, the weights and
the covariance; ``--x-scale 1`` replaces ``least_squares``' default
``x_scale="jac"`` when those rows make it reject every step;
``--strength`` and ``--derived-hook`` add derived quantities. Each fit ends
with ``save_fit(..., close_session=True)``, so one engine is in memory at a
time; refused variants are skipped with the reason, and ``--dry-run`` lists
the grid without numpy or the engine. The theory, every option and the
caveats are in :doc:`../theory/thm_implementation`, "Model averaging".

The transform to physical parameters (``transform_rwa``, ``transform_physical``)
works at the level energies, which a fit or a sampler moves at every call, and
evaluates the external part of each radiative width there. It shares Coulomb
functions within the call but leaves none in the memo afterwards, so a sampler
that transforms at every step (physical priors, a fixed total gamma width) keeps
a constant footprint. Before October 2026 every such call added a few thousand
entries for good: about 4 GB per worker after 17,000 evaluations of a
12C(α,γ) chain.

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
   * - ``fit_scipy.py``
     - ``scipy.optimize.minimize`` with the gradient and ``least_squares`` with the Jacobian on ``Objective``.
   * - ``fit_iminuit.py``
     - MIGRAD and HESSE of iminuit on ``Objective``.
   * - ``fit_emcee.py``
     - MCMC sampling with ``emcee`` on ``Objective.log_prob``, one session per pool worker.
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

There are no AZURE2 instances to pool — every ``Session()`` object is already an
independent engine in-process. Two sessions can evaluate on two threads at
once (above), but a sampler's walkers share one session, whose calls run one
after the other, so the simple unit of parallelism is the *process*:
give each worker its own session by
constructing it at module scope in the worker module. That works whichever way
the pool starts its workers — under ``spawn`` each re-imports the module and
builds its own engine, under ``fork`` each inherits a copy — and each then
evaluates one χ² per walker:

.. code-block:: python

   # worker module
   azr = Session("13N.azr")

   f = azr.fit_objective()

   def log_prob(theta):
       return f.log_prob(theta)

See ``pyazr/examples/fit_emcee.py`` for the full recipe (``fit_zeus.py`` is the older pattern with
``calculate_chi2_rwa`` and priors of its own).

Inside one session the engine uses threads: every evaluation (``calculate_chi2_rwa``,
``calculate_chi2``, ``calculate_rwa``, ``calculate``, ``calculate_all_rwa``) spreads
its data points over ``OMP_NUM_THREADS`` threads, as the command-line χ² does, and so
does ``residual_jacobian``. (Before 2026-10-30 the evaluations ran their points one at a
time on one thread.) With several worker processes, keep workers × ``OMP_NUM_THREADS``
at or below the cores you have; ``AZURE_SERIAL_CHI2=1``, set before the first
session is made (the environment is read once), restores the one-at-a-time
evaluation. The results agree with it to rounding (a few units in the last place).

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

``data_mode()`` (the default) evaluates the active ``dataset`` records — this
is what χ² uses. ``extrap_mode()`` switches to the ``extrapolation`` grids for
predictions on arbitrary energies and angles. Both re-initialise the engine, so
switch sparingly and batch the work.

Segment indexing differs between them, and getting it wrong misaligns results
silently: in data mode segment *i* is ``azr.datasets[i]``; in extrapolation
mode **only active extrapolations are returned**, so segment *i* is
``azr.extrapolations.active[i]``.

.. note::

   Delete ``output/intEC.extrap`` whenever an extrapolation grid changes.
   AZURE2 caches external-capture integrals there; a session that starts on a
   changed grid detects it from ``intEC.extrap.sig`` and recomputes, but a
   running session does not re-check. It is safe to delete; it only costs time
   to rebuild.

Extrapolation grids
-------------------

``pyazr.extrapolation_grid(text, masses=None, relativistic=False,
project_dir="", angles=False)`` gives the energies of a grid written as in a
format-3 ``extrapolation`` record (:ref:`extrapolation-grids`) -- the engine's
own reader, so the values are the ones the record gives -- as lab energies in
MeV; a c.m. grid needs the entrance pair's masses in u::

   E = pyazr.extrapolation_grid("range 20 800 keV cm points 200 log", masses=(1.00728, 12.0))

``pyazr.refine_grid(start, model, tolerance, max_points)`` runs the refinement
of an ``adaptive`` grid on any Python function ``model(xs)`` (one value, or one
sequence of values, per abscissa) and returns the grid, the values and whether
``max_points`` stopped it.  A session (``pyazr.Session``) in extrapolation mode
refuses a project with an ``adaptive`` grid: the session evaluates at many
parameter vectors, and the grid is made by the extrapolation run (menu 3) at
one.

Frames and units
----------------

**Input is lab frame** — the energies and angles in a ``.azr`` and in data
files. **All output and every API result is centre-of-mass.** This includes
extrapolation grids: the energies you *set* are lab, the energies you *get
back* are c.m. Convert with ``E_lab = E_cm · (m_beam + m_target) / m_target``.

Plot different reaction channels against **excitation energy**
(``calculate_excitation_energy``); it is the only axis shared by all entrance
channels.

The web UI's server
-------------------

``python -m pyazr.ui project.azr`` serves a project to the web UI's API (JSON over HTTP on 127.0.0.1, with a
token); it needs the optional packages of ``pip install pyazr[web]``.  See :doc:`web_ui_api`.
