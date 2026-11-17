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
`io/` (projects and nuclear-data libraries). (That was the plan; what was built is in `FOR_THE_OWNER.md` section 5:
`inference/` holds the parameter vector as plain data, the priors and the sums of the objective's term list, while
the terms and the fitters are engine code (`inference.md`); the arithmetic of the collision matrix and of the observables is in `physics/`,
their analytic gradient and the walk over a site's pathways are engine classes.) The calculation is written in three
lines (`ARCHITECTURE.md` §1), and a
new reaction, observable or apparatus effect enters at one stated place and is tested against synthetic projects
whose true parameters are known. The rebuild is staged: every phase must reproduce the results of the present
engine, or explain each difference, before the code it replaces is removed. It also adds what the present code
lacks: a Trojan-horse (THM) observable with DWBA vertices, AZURE2's own stopping and straggling evaluation, Coulomb
functions from pycoulomb, pyazr 3.0 and a local web interface in place of Qt.

Read in this order:

1. `FOR_THE_OWNER.md` — what the branch is, how to try it, where its numbers differ, how to add to it, what is not
   finished (section 10), and the reading order for an AI agent (section 11).
2. `dev_features.md` — where each feature of `dev` lands in this structure.
3. `ROADMAP.md` — phases 0–9 as planned and ticked; its block "State at the hand-over" is the state per box.
4. `ARCHITECTURE.md` — the target design (a plan: where the code differs, `FOR_THE_OWNER.md` section 10 says so).
5. `THM.md` — the Trojan horse method in the new structure.
6. `RULES.md` — rules for code, tests and skills. `STYLE.md` — how documentation is written.
7. `HANDOFF.md` — the log of every session, newest first (more than 12,000 lines): search it, do not read it whole.
8. `figures/structure.html` — a summary page with the two figures (open locally in a browser). Where it and the
   `.md` files differ, the `.md` files win.

For an agent continuing this work: follow `FOR_THE_OWNER.md` section 11, obey `RULES.md` and the "Upgrade" section
of the top-level `CLAUDE.md`, and add an entry at the top of `HANDOFF.md` at the end of every session.
