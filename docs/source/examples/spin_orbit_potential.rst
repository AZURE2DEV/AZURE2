.. _example-spin-orbit-potential:

Worked Example: p + :sup:`12`\ C with a Spin–Orbit Potential
==============================================================

This example adds a spin–orbit term to the hybrid potential of p + :sup:`12`\ C and shows what it does to the
:math:`p_{1/2}` and :math:`p_{3/2}` waves.  The files are in ``tests/spin_orbit_potential``; the theory is
:doc:`../theory/hybrid_potential`.

The record
----------

.. code-block:: text

   <potentials>
   potential p12C
       pair           p+12C
       coulomb_radius 3.5 fm
       central        woods_saxon depth 20 MeV radius 3.6 fm diffuseness 0.6 fm
       spin_orbit     woods_saxon depth 6 MeV radius 3.6 fm diffuseness 0.6 fm on light operator 2l.s
   </potentials>

The proton (``on light``) is the only particle with spin, so the term is diagonal: :math:`\langle 2\,\mathbf{l}
\cdot\mathbf{i}\rangle = +1` for :math:`p_{3/2}` and :math:`-2` for :math:`p_{1/2}`.  ``AZURE2 --describe`` prints
the term with its convention.

The phase shifts
----------------

``models/phase.azr`` has one level of zero width in each of the :math:`1/2^-`, :math:`3/2^-` and :math:`1/2^+`
groups, so the phase shifts are the potential's, :math:`\delta_{lJ} = -\phi_{lJ}`:

.. code-block:: text

   extrapolation p1half
       reaction       12C(p,p)12C
       quantity       phase_shift
       unit           deg
       energies       list 0.5 1.0 2.0 3.0 MeV cm
       J              1/2
       l              1
       s              1/2

At :math:`E_{cm}` = 3 MeV the engine gives :math:`\delta(p_{1/2}) = -7.238^\circ`, :math:`\delta(p_{3/2}) =
-5.611^\circ` and :math:`\delta(s_{1/2}) = -27.997^\circ` (the s wave has no spin–orbit term); an independent
integration (``physics/tests/data/make_spin_orbit_reference.py``) agrees to 10\ :sup:`-9` rad.

In pyazr
--------

.. code-block:: python

   m = pyazr.azure2("hybrid_potential.azr")
   m.set_nuclear_potential(1, spin_orbit={"depth": 6.0, "radius": 3.6, "diffuseness": 0.6, "on": "light"})
   p12 = m.coulomb_functions(1, [0.5, 1.0, 3.0], L=1, J=0.5)
   p32 = m.coulomb_functions(1, [0.5, 1.0, 3.0], L=1, J=1.5)

At :math:`E_{cm}` = 0.5 MeV and :math:`a_c` = 3.6 fm the shift function is :math:`-0.7518` (:math:`p_{1/2}`) and
:math:`-0.2428` (:math:`p_{3/2}`), the penetrability 3.583·10\ :sup:`-3` and 4.418·10\ :sup:`-3`.  The elastic
:math:`\chi^2` of the Meyer data at 144° moves from 2259.87 to 2117.49 with the term; the level parameters are
re-derived with the new :math:`S` and :math:`P`, as for any change of the potential.
