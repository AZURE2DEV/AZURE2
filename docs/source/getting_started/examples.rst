Examples
========

The projects in ``examples/`` are Trojan Horse (THM) evaluations, each with
the direct data it is fitted with. A project folder holds ``<name>.azr`` (the
best fit), ``data/`` (every file names its source and what its points are in
``#`` header lines) and empty ``checks/`` and ``output/`` directories. Open
the ``.azr`` in the GUI (THM is switched on by itself for these projects:
:doc:`../user_guide/configure_menu`, "Runtime Options") or run it from the
command line. The THM options are described in
:doc:`../theory/thm_implementation`.

"File alone" is the CLI calculation of the file as stored,

.. code-block:: bash

   cd examples/<name>
   printf '1\n\n\n7\n' | AZURE2 --no-gui --no-readline <name>.azr

(menu 1, calculate with data; blank answers take the parameters and the
external-capture integrals from the project; 7 exits). The number given is
the ``Total Chi-Squared`` line, data plus norm priors; the data part
(``Total-Chi-Squared`` of ``output/chiSquared.out``) and the priors are in
brackets where they differ. Every direct normalization with a quoted error
holds its fitted value in the norm column and an explicit ``prior_centre 1``
row (:ref:`nominal-norm`), so the file alone reproduces the fit and a refit
is pulled to 1. Reproduced with the ``thm`` branch in October 2026 (one run
of each: 11 s for 19F, 25 s for 12C+12C, the others a few seconds).

THM (HOES) fits
---------------

.. list-table::
   :header-rows: 1
   :widths: 18 27 33 22

   * - example
     - reaction and data
     - what it shows
     - file alone / N
   * - ``o18_lacognata2008``
     - 18O(p,α)15N; THM only: La Cognata et al., PRL 101 (2008) 152501,
       Fig. 2, with their fitted linear background subtracted
     - the three narrow levels below 0.2 MeV in a HOES spectrum; the
       starting point for strengths from peak areas against an anchor
     - 35.877 / 30
   * - ``o18_lacognata2010``
     - 18O(p,α)15N from 2H(18O,α15N)n at 54 MeV; THM: La Cognata,
       Spitaleri & Mukhamedzhanov, ApJ 723 (2010) 1512, Fig. 4, the
       authors' background subtracted; direct: Mak 1978, Amsel 1967
       (differential)
     - the interfering 1/2\ :sup:`+` doublet; the 3/2\ :sup:`−` 597.6 keV
       level at the LUNA widths, held fixed; an experiment line with a
       profiled linear background and the distorted-wave vertex with global
       optical potentials (``ancai06``, ``kd03:extrapolate``)
     - 330.054 (329.692 + 0.362) / 161; THM 247.14 / 39
   * - ``o17_guardo2017_fit``
     - 17O(n,α)14C from 2H(17O,α14C)p; Guardo et al., PRC 95 (2017) 025807
     - genuine HOES yields: points below the neutron threshold, the f-wave
       5\ :sup:`−` peak next to the p-wave 3\ :sup:`−`; a neutron entrance
       channel with a charged spectator
     - 12.996 / 23
   * - ``f19_pag_thm``
     - 19F(p,αγ)16O from 2H(19F,α16O)n at 55 MeV; THM: Su et al., PRL 135
       (2025) 182701, Fig. 1, :math:`E \le 0.45` MeV; direct: JUNA (Zhang
       et al., PRC 106 (2022) 055803, Table I), Spyrou et al. (2000)
     - HOES and direct data together; ``kinematics=triple``; a profiled
       linear background; a fixed THM energy shift of −9.17 keV (lab); the
       driver example of "Model averaging"
     - 79.408 (77.640 + 1.768) / 59; THM 54.17 / 28; direct norms stored
       fixed
   * - ``c12c12_tumino2018``
     - 12C(12C,α)20Ne and 12C(12C,p)23Na (α\ :sub:`0,1`, p\ :sub:`0,1`)
       from 12C(14N,α20Ne / p23Na)d at 30 MeV; THM: Tumino et al., Nature
       557 (2018) 687, Fig. 1; direct: Spillane 2007 (γ-ray production), Tan 2024,
       Nippert 2025, Jiang 2018
     - four exit channels as one experiment with one profiled norm
       (``experiment[E1] segments=1-4`` with the three-body kinematics, so
       the line shape, :math:`R(E)` and the DW vertex can be switched on);
       an identical-boson entrance pair
     - 112.076 (111.651 + 0.425) / 247; THM 61.02 / 192; direct norms
       stored fixed
   * - ``o17_guardo``
     - 17O(n,α)14C, the Guardo et al. 2017 points as first added to the
       repository
     - the original input: channel widths entered as reduced-width
       amplitudes (field 33) and no ``<thm>`` block (all defaults); its
       ``output/`` holds the files it was committed with
     - 67.594 / 23

On-shell-equivalent THM data
----------------------------

These papers publish THM points that are already corrected for the
penetrability and normalized to direct data. They are fitted as ordinary
on-shell cross sections with a free norm, so no THM option acts on them
(:doc:`../theory/thm_implementation`, "HOES yields or on-shell-equivalent
points?").

.. list-table::
   :header-rows: 1
   :widths: 18 37 23 22

   * - example
     - reaction and data
     - what it shows
     - file alone / N
   * - ``n15_lacognata2007``
     - 15N(p,α\ :sub:`0`)12C; THM: La Cognata et al., PRC 76 (2007) 065804,
       Table 3 (an on-shell S factor); direct: Redder 1982, Schardt 1952
     - an S-factor table converted to a cross section; direct norms with
       their quoted errors
     - 224.425 (223.414 + 1.012) / 208
   * - ``li7_tumino2006``
     - 7Li(p,α)4He; THM: Tumino et al., EPJA 27 s01 (2006) 243, at
       :math:`\theta_\mathrm{cm}` = 60°; direct: Rolfs & Kavanagh 1986,
       Cruz 2005/2008, Cassagnou 1962, Mani 1964
     - an identical exit pair (data that count both α halved); a free norm
       with no prior (Cassagnou)
     - 448.575 (439.104 + 9.471) / 294
   * - ``li6_pizzone2011``
     - 6Li(d,α)4He; THM: Pizzone et al., PRC 83 (2011) 045801; direct:
       Engstler 1992, Elwyn 1977, McClenahan 1975, Jeronymo 1962
     - the same reading for a second Li reaction
     - 158.918 (156.991 + 1.927) / 132

``examples/run_gui.sh`` and ``examples/run_mcmc.sh`` start the GUI and an
MCMC script in the Docker image of ``packaging/docker`` (README.md,
"Containers").
