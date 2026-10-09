What's new on the thm branch
============================

The ``thm`` branch adds the Trojan Horse Method to AZURE2 and carries a set of
fixes and features for classic projects; all of ``dev`` up to 6f3228e is
merged into it. ``MERGE_NOTES.md`` at the top of the repository lists, with
commits and re-pinned tests, every difference a classic project (no
``<thm>`` block, no THM segment) sees against ``dev``. This page is the
summary.

Trojan Horse Method (opt-in)
----------------------------

Nothing below acts on a project without THM content. In the GUI it is
switched on by *Configure > Runtime Options > Use Trojan Horse Method
(THM)*, and by itself for a project that has THM content; the setting is not
stored in the file. The engine and ``pyazr`` need no switch.

- **The HOES observable** of a THM segment (observable code + 10), with the
  binding energy of the Trojan horse in field 32 of the entrance pair's level
  lines and channel widths optionally entered as reduced-width amplitudes
  (field 33). Options in the optional ``<thm>`` block: incoherent or coherent
  entrance partial waves, the vertex boundary (``constant`` by default,
  representation-invariant), the kinematic factor of the data reduction, the
  external Coulomb term, a spectator energy, and energy-dependent weight
  tables (:doc:`../theory/thm_implementation`).
- **THM experiments**: segments of one measurement share one profiled
  normalization and an optional smooth background
  (``output/thm_experiments.out``).
- **Physics of the three-body reaction**, per experiment: a coherent THM-only
  background with fitted complex amplitudes (``cbkg`` parameters), the
  Coulomb line shape :math:`N_C` of a charged spectator, the
  spectator-momentum window and the angular acceptance (with the measure at
  fixed :math:`E`), the distortion factor :math:`R(E)` from a zero-range
  DWBA, a distorted-wave entrance vertex, built-in global optical potentials,
  and a fixed-angle observable. Combinations that would count a Coulomb
  interaction twice are refused.
- **Brune and Park**: THM segments under ``--use-park`` give Brune's model;
  an amplitude in the file is Brune's in both modes.
- **Fits and errors**: the profiled scale enters the Jacobian, the MIGRAD
  gradient and the covariance band exactly; ``prior_centre`` rows keep norm
  priors on their experimental value when a file stores fitted norms.
- **GUI**: *Configure > THM Workspace...* with the pages Model, Experiments,
  Channels and Diagnostics (vertex, HOES and on-shell shapes, line shape,
  weight, window, :math:`R(E)` and angular distribution, computed by the
  engine on a copy); the *THM* tick in the segment dialogs; the *THM
  Background* sub-tab and the *Prior Centre* column of the Fitting tab.
- **pyazr**: ``AzrModel`` edits the ``<thm>`` block with the engine's rules
  (``set_thm_option``, ``set_thm_experiment``, ...); a session reports
  ``thm_vertex``, ``thm_distortion``, ``thm_lineshape``, ``thm_background``;
  ``pyazr.modelavg`` and ``scripts/thm_model_average.py`` fit and average a
  grid of model variants, each in a fresh process.
- **Examples**: nine THM projects, most with their direct data
  (:doc:`examples`).

For classic projects
--------------------

*Results that move* (each re-pinned and explained in ``MERGE_NOTES.md``,
section A): the shift function below threshold and one helper for
:math:`S`, :math:`P` and :math:`dS/dE` at a channel energy (no NaN at or near
threshold); a larger, bounded step for :math:`dS/dE`; the adaptive
integration grid (geometric tails, quantised anchors, refreshed during a
fit); :math:`1 + \delta_{12} = 2` for every cross section out of an identical
entrance pair (dev applied it to elastic scattering only); a
``<targetInt>`` line shared by several segments converted once per segment;
nuisance priors and fixed parameters honoured by every minimizer.

*From dev, merged*: the :math:`\sqrt{E}` energy-shift term per segment
(``sqrtshift`` block), Park's parametrization (``--use-park``) and
``--no-brune``, observable code 8 (polarization times cross section)
removed, the energy-dependent convolution window from the c.m. energy, the
resident-memory fix of long fits, ``param.fit`` at the evaluated point, and
the ``#parametrization`` line of parameter files with readers that match
rows by name.

*Fixed on this branch* (and present on ``dev``): end of input at a CLI
prompt that needs an answer stops the run (exit 1) instead of looping; MCMC
classifies the parameters by name; an MCMC start under Park with no valid
walker is refused, and a failed MCMC run exits non-zero; the Park
:math:`J > 0` penalty is part of every minimizer's cost; ``param.fit``
carries the ``#parametrization`` line, so a fit under ``--use-park`` resumes
from it.

*Behaviour*: the CLI exits non-zero when a run fails; a malformed
``<potential>`` block is refused; the console total is printed to twelve
digits and includes every prior; ``--help`` lists exactly the flags the
command line takes. The GUI saves doubles with round-trip
precision, keeps level lines' fields 32 and 33, and renumbers what refers to
a data segment by number when segments are moved, added or deleted.
``pyazr``'s ``save_fit`` writes the fitted norms and the matching
``prior_centre`` rows and verifies the snapshot.

*Speed*: a Coulomb-function memo that gives up on nearby energies keeps the
exact energies it has seen, so a second session in one process (``pyazr``,
the GUI's diagnostics) starts much faster; results are unchanged.

*Tests*: a time limit on every engine run (``tests/lib/guard.sh``), shared
check helpers, the exit-77 skip convention of the ``pyazr`` tests, and CI that
fails when those tests are skipped (:doc:`../developer/contributing`).
