Radiative Widths as Total Width and Branchings
==============================================

A level that decays by several γ-ray transitions is usually known from experiment as a total radiative width
:math:`\Gamma_\gamma` (or a resonance strength), the branchings :math:`b_f` to its final states and, for a transition
with two multipolarities, the mixing ratio :math:`\delta_f`. In ``<levels>`` the same level carries one partial width
(eV) per final state and multipolarity. A ``gamma_decay`` record in a section ``<gammaDecays>`` lets you write, fix,
fit and read the level's radiative widths in the first form; the engine computes the partial widths from it at every
step of a fit and is otherwise unchanged. The design note is ``docs/upgrade/gamma_parametrisation.md``.

Definitions
-----------

For the final states :math:`f` of the level's photon channels, each with one multipolarity :math:`L` or two
(:math:`L` and :math:`L+1`):

.. math::

   \Gamma_{f,L} = \sigma_f\,\frac{\Gamma_\gamma b_f}{1+\delta_f^2}, \qquad
   \Gamma_{f,L+1} = -\sigma_f\,\frac{\Gamma_\gamma b_f\,\delta_f|\delta_f|}{1+\delta_f^2}, \qquad
   \sum_f b_f = 1 .

* :math:`\Gamma_\gamma` (eV) is the sum of the level's radiative partial widths. As the radiative widths themselves
  (:doc:`../theory/conventions`), it is **not** part of the resonance's total width: the total width in the level
  matrix is the sum of the particle widths.
* :math:`b_f` is the branching of the γ decay to :math:`f` (not of the level's total width). One final state, the
  ``remainder``, takes :math:`1-\sum_{\text{others}} b_f`; it is not a fit parameter.
* :math:`\delta_f` is the mixing ratio in the convention of Krane and Steffen (Phys. Rev. C 2, 724 (1970); also
  ENSDF's), :math:`\delta^2 = \Gamma_{L+1}/\Gamma_L`. Its sign was fixed against the engine's γ angular distributions:
  partial widths of the same sign are :math:`\delta < 0`.
* :math:`\sigma_f = \pm 1` is the interference sign of the final state (``sign``), as the sign of a width in
  ``<levels>``: it matters only for the interference with other levels' amplitudes to the same final state.

Writing the record
------------------

The 259.7 keV resonance of :sup:`14`\ N(p,γ):sup:`15`\ O (level ID 5, :math:`J^\pi = 1/2^+`) of
``tests/gamma_decay/models/record.azr``:

.. code-block:: text

   <gammaDecays>
   gamma_decay level5
       level        5
       J            1/2+
       total_width  0.0385 eV fit prior 0.040 0.004
       branching    g+15O 0.0150 fit
       branching    g+15O[6.7931] 0.230 fit prior 0.23 0.02
       branching    g+15O[6.1763] 0.584 fit
       branching    g+15O[6.8594] 0.0014 fixed
       remainder    g+15O[5.183] prior 0.17 0.02
       mixing       g+15O[6.7931] -0.12 fit
       sign         g+15O[6.1763] -
   </gammaDecays>

``level`` and ``J`` name the level as a ``<modelData>`` datum does. A final state is its photon pair, written as a
reaction label writes it (``g+15O[6.7931]``; ``g+15O`` is the ground state). Every number is
``<value> fixed|fit [prior <centre> <width>]``; a measured branching (ENSDF, a paper) is a prior. A final state of the
level's channels that the record does not name has :math:`b_f = 0`, fixed; a final state with two multipolarities
and no ``mixing`` is pure :math:`L`. The level's photon lines stay in ``<levels>`` (they define the channels) with
width 0; their fixed flags are not read.

The engine refuses, with a code: a photon width of the level written in ``<levels>`` as well
(``E_GAMMA_MIXTURE``: write 0 there or delete the record); fixed branchings that do not sum to 1, or a fitted one
without a ``remainder`` (``E_BRANCHING_SUM``); branchings that leave a negative remainder (``E_BRANCHING_REMAINDER``); a
mixing ratio for a final state with one multipolarity (``E_MIXING_ONE_MULTIPOLE``). At set-up it refuses the formal
basis (``--no-brune``) and a level with external capture on one of the record's channels (write such a level with
partial widths).

Fitting and reading
-------------------

The record's numbers are fit parameters of the kind ``radiative``, named ``level5.total_width``,
``level5.branching.g+15O[6.7931]``, ``level5.mixing.g+15O[6.7931]``; they appear in ``param.par``, ``param.sav``
and pyazr's ``Session.parameters`` (``Parameter.decay_record``, ``decay_part``, ``final_state``). A branching is kept
in :math:`[0, 1]` and :math:`\Gamma_\gamma \ge 0`. A mixing ratio may start at 0, but angle-integrated data depend on
:math:`\delta^2` alone, so they do not determine its sign and :math:`\delta = 0` is a stationary point there: start away
from 0 or fit an angular distribution as well. When the free branchings sum above 1 the remainder is 0 and a wall
in the objective brings the sum back. ``parameters.out`` ends with a block that gives :math:`\Gamma_\gamma`, every
:math:`b_f` (the remainder included), :math:`\delta_f` and the partial widths, with a note for a branching at 0.

From Python, an agent reads the record, changes one branching and refits:

.. code-block:: python

   from pyazr import Project, Nuisance, azure2

   project = Project.open("record.azr")
   print(project["level5"]["branching"])            # the lines as written
   project.set("level5", total_width=Nuisance(0.0385, "fit", unit="eV"))
   lines = [("g+15O[6.7931] 0.25 fit" if l.startswith("g+15O[6.7931]") else l)
            for l in project["level5"]["branching"]]
   project.set("level5", branching=lines)
   project.save("record_changed.azr")
   with azure2("record_changed.azr") as s:
       print([(p.name, p.value) for p in s.parameters if p.kind == "radiative"])

A new record is ``project.add("gamma_decay", "level5", level=5, J="1/2+", total_width=..., branching=[...],
remainder="g+15O[5.183]")``.

Resonance strengths
-------------------

Because the radiative widths stay out of the total width of the level matrix, the strength of a narrow resonance in a
calculated yield is :math:`\omega\Gamma_\gamma\Gamma_p/\sum\Gamma_\text{particle}`. The ``resonance_strength`` datum
of ``<modelData>`` (:doc:`../reference/data_formats`) divides by :math:`\Gamma = \sum\Gamma_\text{particle} +
\Gamma_\gamma` instead, so the two differ by the factor :math:`1+\Gamma_\gamma/\sum\Gamma_\text{particle}`: 4×10⁻⁵ for the
259.7 keV level of :sup:`14`\ N(p,γ) (:math:`\Gamma_p = 975` eV, :math:`\Gamma_\gamma = 0.0385` eV), 2 for
:math:`\Gamma_p = \Gamma_\gamma`.
