# The AZURE2 upgrade: where it stands, and how to add to it

A short account for the owner of AZURE2, 2026-11-09. The details are in `HANDOFF.md`, `ROADMAP.md` and the design
notes of `docs/upgrade/`.

## Why

AZURE2's physics is right; its structure made it fragile. One object held a datum, its channel tables, its computed
values and the apparatus corrections, and parallel evaluation copied and shared it. Of the roughly 126 fixes on
`dev` between August and November 2026, about 55 % touched those objects and at least 21 repaired defects that an
earlier performance change had introduced. Every feature also had to be written several times: the χ² terms in four
code paths, the `.azr` grammar in the engine, in pyazr and in the GUI, the masses in two tables. A feature added to
one copy and not the others gives a wrong number that looks right, which is the worst bug a code used for
evaluations can have.

The upgrade (branch `feature/upgrade`) rebuilds the engine so that its structure follows the calculation — levels →
collision matrix → observables → apparatus → χ² — with one implementation of each formula and format.

## What is done (Phases 0–3 of 9)

| Phase | What it gives | Checked by |
|---|---|---|
| 0. Measure and protect | Synthetic projects T1–T3 with known true parameters; golden outputs; a speed benchmark; nightly thread and memory sanitizers; the rules in `CLAUDE.md`; ten skills | — |
| 1. Physics core | `physics/`, a standalone library (no engine headers, no GSL): Coulomb and Whittaker functions from pycoulomb, level-parameter conventions and transformations (Brune, formal, Park), external capture, coupling coefficients, reaction pathways. The Coulomb memo is replaced by a per-stage table: results no longer depend on the thread count | mpmath and closed forms; outputs bit-identical or every difference classified |
| 2. One objective | One ordered list of χ² terms used by the fit, the gradient, least squares, the MCMC and the Python API; pyazr computes no prior itself | χ² bit-identical where no decision applied; decisions D1–D7 recorded |
| 3. One project model | `io/`: the `.azr` read once, losslessly, by the engine and pyazr alike; a format version, `AZURE2 --convert`, unchanged round trips of every file; one mass table | 179 projects × 3 bases identical |

After every phase: the 179 project files of the test suite in three bases (Brune, formal, Park) compared file by
file with the previous engine, the 133 regression projects, about 130 unit and pyazr tests, and 1 against 8 threads.

**Speed** (T2, 340 points with every apparatus kernel; this machine, same build options; `tests/upgrade/bench/`):

| T2, mean of three runs | before Phase 1 (93f8cbc5) | after Phase 3 | change |
|---|---|---|---|
| calculation, 1 thread | 0.93 s | 0.90 s | −3 % (within the spread) |
| calculation, 8 threads | 0.42 s | 0.24 s | −43 % |
| fit, 200 calls, 1 thread | 63.2 s | 56.2 s | −11 % |
| fit, 200 calls, 8 threads | 36.4 s | 30.7 s | −16 % |
| fit peak memory, 8 threads | 488 MB | 643 MB | +32 % |

The memory comes from tables kept for the numerical derivative with respect to energy shifts; they go when that
derivative becomes analytic. Speed was not the aim of Phases 1–3 (structure was); Phases 5 and 6 are where it is
won.

## What users notice

- Old `.azr` files load as before: they are converted in memory, with a note in the log.
  `AZURE2 --convert old.azr new.azr` (or `pyazr.convert_azr`) writes the new format, which drops sections nothing
  reads any more.
- pyazr needs the compiled engine module to read a project: there is one reader, the engine's.
- Results are unchanged except where a decision was taken and recorded (for example D1: the MCMC now profiles THM
  normalisations, which corrected a diverging chain). Each change is listed with its size in the design notes.
- No Qt GUI on this branch: a local web interface replaces it (Phase 9), checked against a 201-item parity list
  taken from the Qt code.
- The Coulomb table is rebuilt once (format 2).

## How to add a feature now

`CLAUDE.md` (section "Recipes") lists, for each kind of addition — a segment keyword, a `<config>` key or section, a
χ² term, a physics function, a command-line flag — every place to touch. Each recipe has a check in
`tests/consistency/repo_consistency_test.py` that fails when it is bypassed:

- a second parser of the project file, or a section read outside `io/` → fails (check 12);
- a χ² term added outside `src/Objective.cpp` → fails (check 11);
- a copy of a physical constant with another value → fails (check 10);
- a segment keyword known to the engine and not to pyazr, a flag missing from `--help` or the docs, an unregistered
  test → fail (checks 1–6).

`scripts/check.sh --quick` runs these in under a second, then the build and the unit tests; `scripts/check.sh` adds
the 133 regression projects. An AI assistant working in the repository reads `CLAUDE.md` first, so it is given the
recipes and is stopped by the checks. **Ask of a new feature: does `scripts/check.sh` pass, and which test would
fail if the feature were wrong?**

## What comes next

Phases 4–9 (`ROADMAP.md`): data as arrays converted to c.m. at load (4); one energy grid with the apparatus as
kernels, AZURE2's own stopping, and the THM observable (5); evaluation without copies (6); DWBA vertices for THM (7);
pyazr 3.0 with wheels (8); the web interface, and removal of Qt, ERYA and numcmc (9). The first new feature to go
through the recipes is extrapolation segments on a user-given or adaptive energy grid.

**`dev` and the upgrade.** `dev` has 18 commits since the upgrade branched. They are merged when the reconstruction
is complete, not before. Most are documentation; the target-effect windows, the UDR kernels and test fixes port as
they are; three changes (the heavy-ion Coulomb fallback, the coincident-level limits, pyazr NUTS sampling) are redone
in the new structure; two are obsolete (the assessment is in `HANDOFF.md`).
Until then, new physics is best written on the upgrade branch, following the recipes.

**Known open problem.** A fit's path still depends on the number of threads: on T2 the χ² after 200 MIGRAD calls
differs in the 9th digit between 1, 2 and 4 threads, which an unconverged fit amplifies (691.6 against 595.9 at 1
and 8 threads). Calculations are identical at any thread count. It was already present when Phase 1 began, and is traced
in `HANDOFF.md`; Phase 6 (evaluation without copies) removes the suspected cause if it is not found sooner.

**Two values for the maintainers to decide**, since changing them moves results: the electron mass in the
β-decay Fermi function (CODATA 2002) and the proton mass in `ECIntegral.cpp` (older than CODATA 2018).
