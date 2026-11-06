# AZURE2 upgrade

AZURE2 computes cross sections, angular distributions, polarisation observables and astrophysical S factors from a
multi-level, multi-channel R-matrix description of the compound nucleus, and fits that description to measured data
corrected for the experimental apparatus. Every result it produces passes through the same chain: level parameters
determine the collision matrix at each energy; the collision matrix determines the observables; the apparatus
(target thickness, stopping and straggling, resolution, Doppler broadening) folds those observables into the
quantities actually measured; and the χ² compares them with the data. The present engine implements this chain
correctly but not visibly: one datum, its channel tables, its computed results and the apparatus corrections share a
single object, and parallel evaluation copies and re-shares that object. About 55 % of the roughly 126 fixes made
on `dev` between August and November 2026 touched those objects, and at least 21 repaired defects introduced by
earlier performance changes; the physics itself was rarely at fault.

The upgrade rebuilds the engine so that its structure follows the chain itself. Each stage computes one quantity and
depends only on the stages before it: `physics/` (levels, channels, collision matrix, observables), `experiment/`
(measurements, materials, stopping, the apparatus kernel), `inference/` (parameters, the objective, the fitter) and
`io/` (projects and nuclear-data libraries). The calculation is written in three lines (`ARCHITECTURE.md` §1), and a
new reaction, observable or apparatus effect enters at one stated place and is tested against synthetic projects
whose true parameters are known. The rebuild is staged: every phase must reproduce the results of the present
engine, or explain each difference, before the code it replaces is removed. It also adds what the present code
lacks: a Trojan-horse (THM) observable with DWBA vertices, AZURE2's own stopping and straggling evaluation, Coulomb
functions from pycoulomb, pyazr 3.0 and a local web interface in place of Qt.

Read in this order:

1. `HANDOFF.md` — where things stand, what is open, what is next.
2. `ARCHITECTURE.md` — the target design.
3. `ROADMAP.md` — phases 0–9 with tasks and gates.
4. `THM.md` — the Trojan horse method in the new structure.
5. `RULES.md` — rules for code, tests and skills.
6. `STYLE.md` — how documentation is written.
7. `figures/structure.html` — a summary page with the two figures (open locally in a browser). Where it and the `.md` files differ, the `.md` files win.

For an agent continuing this work: start from `HANDOFF.md`, follow `ROADMAP.md` phase by phase, obey `RULES.md` and
the "Upgrade" section of the top-level `CLAUDE.md`, and update `HANDOFF.md` and the roadmap ticks at the end of every
session.
