# AZURE2 upgrade

A staged rebuild of the AZURE2 engine into a structure shaped like the physics: fast C++, one-way module
dependencies, one evaluation chain shared by every reaction and apparatus, a simple THM with DWBA vertices, AZURE2's
own stopping and straggling, Coulomb functions from pycoulomb, pyazr 3.0 and a local web UI that replaces Qt.

Read in this order:

1. `HANDOFF.md` — where things stand, what is open, what is next.
2. `ARCHITECTURE.md` — the target design.
3. `ROADMAP.md` — phases 0–9 with tasks and gates.
4. `THM.md` — the Trojan horse method in the new structure.
5. `RULES.md` — rules for code, tests and skills.
6. `figures/structure.html` — a summary page with the two figures (open locally in a browser). Where it and the `.md` files differ, the `.md` files win.

For an agent continuing this work: start from `HANDOFF.md`, follow `ROADMAP.md` phase by phase, obey `RULES.md` and
the "Upgrade" section of the top-level `CLAUDE.md`, and update `HANDOFF.md` and the roadmap ticks at the end of every
session.
