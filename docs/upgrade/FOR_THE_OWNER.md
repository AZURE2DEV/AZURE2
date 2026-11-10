# The AZURE2 upgrade: where it stands, and how to add to it

A short account for the owner of AZURE2, 2026-11-10. The details are in `HANDOFF.md`, `ROADMAP.md` and the design
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

## What is done (Phases 0–3 of 9; Phase 4 begun)

| Phase | What it gives | Checked by |
|---|---|---|
| 0. Measure and protect | Synthetic projects T1–T3 with known true parameters; golden outputs; a speed benchmark; nightly thread and memory sanitizers; the rules in `CLAUDE.md`; ten skills | — |
| 1. Physics core | `physics/`, a standalone library (no engine headers, no GSL): Coulomb and Whittaker functions from pycoulomb, level-parameter conventions and transformations (Brune, formal, Park), external capture, coupling coefficients, reaction pathways. The Coulomb memo is replaced by a per-stage table: calculations no longer depend on the thread count | mpmath and closed forms; outputs bit-identical or every difference classified |
| 2. One objective | One ordered list of χ² terms used by the fit, the gradient, least squares, the MCMC and the Python API; pyazr computes no prior itself | χ² bit-identical where no decision applied; decisions D1–D7 recorded |
| 3. One project model | `io/`: the `.azr` read once, losslessly, by the engine and pyazr alike; a format version, `AZURE2 --convert`, unchanged round trips of every file; one mass table | 179 projects × 3 bases identical |

After every phase: the 179 project files of the test suite in three bases (Brune, formal, Park) compared file by
file with the previous engine, the 133 regression projects, about 130 unit and pyazr tests, and calculations at 1 against 8 threads.

**Speed** (T2, 340 points with every apparatus kernel; this machine, same build options; `tests/upgrade/bench/`):

| T2, mean of three runs | before Phase 1 (93f8cbc5) | after Phase 3 | change |
|---|---|---|---|
| calculation, 1 thread | 0.93 s | 0.90 s | −3 % (within the spread) |
| calculation, 8 threads | 0.42 s | 0.24 s | −43 % |
| fit, 200 calls, 1 thread | 63.2 s | 56.2 s | −11 % |
| fit, 200 calls, 8 threads | 36.4 s | 30.7 s | −16 % (see the open problem below) |
| fit peak memory, 8 threads | 488 MB | 643 MB | +32 % |

The memory comes from tables kept for the numerical derivative with respect to energy shifts; they go when that
derivative becomes analytic. Speed was not the aim of Phases 1–3 (structure was); Phases 5 and 6 are where it is
won.

## What users notice

- Today old `.azr` files still load (converted in memory, with a note in the log), and
  `AZURE2 --convert old.azr new.azr` (or `pyazr.convert_azr`) writes the current format. When format 3 lands
  (Phase 4, below) the engine reads only that; the converter turns any older file into it.
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

**Phase 4, begun: data are data** (`data_model.md`, `data_kinds.md`). Two things:

1. *Structure.* A datum becomes an immutable record; what an evaluation computes goes to a buffer the caller owns;
   what the model predicts for a data set is a small expression tree (a cross section, a sum, a ratio, a scale) with
   one forward and one adjoint implementation. That composition was written in about 20 loops and twice more in
   the gradient, and the residual formula 13 times; since 2026-11-10 both are one function each (the model pass,
   `src/ModelPass.cpp`; a check fails on a second copy). The gradient's hand-written adjoint is the next to go.
2. *What a user can fit.* A data set will state in words what it is — quantity, unit, frame, reference, target,
   normalisation — instead of an integer code and column positions:
   `quantity dsigma/dOmega_ratio reference rutherford frame lab …`, `quantity yield target C13_a per charge q=1 …`.
   Targets get names and are defined once. Model-level data have their own section, `<modelData>` (done:
   resonance strength ωγ, resonance energy, total and partial width, branching ratio, ANC, scattering length — values
   from the literature compared with a function of the level parameters, as one χ² term, `modelData.out`,
   `pyazr.model_data`). First new kinds of beam data: ratio to Rutherford or Mott, S factor, thick-target yields per charge, ωγ; the order
   follows what EXFOR holds for light nuclei (10,132 charged-particle data sets counted, `data_kinds.md` §3).

**Decision taken with Jakub (2026-11-09): no backward compatibility inside the engine.** The engine reads one
grammar (format 3); `AZURE2 --convert` and `pyazr.convert_azr` turn any old file into it, and an old file given
directly is refused with a message naming the converter. The gate is that every converted test project gives outputs
identical to the original on the previous engine. This supersedes "old files load as before" above once format 3 lands.

**For programs and AI agents:** `AZURE2 --schema` (every record kind, field, quantity, unit and diagnostic code as
JSON; done, with `pyazr.schema()` and a reference page generated from it), `AZURE2 --validate FILE [--json]` (all
problems of a file at once, each with a code and a line; done, with `pyazr.validate`) and `--describe` (how the engine
understood each data set: quantity, reactions, linked effects, parameters; in work). One table of quantities in the
engine feeds the validation, Python, the GUI and the documentation: a new kind of data is one row of that table, a
tree builder and a test.

**The converter exists** (`convert/`, `AZURE2 --convert OLD NEW --to 3`, `pyazr.convert_azr`): of the 189 project
files of the test tree, 186 convert to format 3 and validate without a question, and the 3 refused are files the old
engine refuses too. It reports the map from the old parameter names (`segment_5_norm`) to the new ones, which carry
the data set's name. The engine's set-up no longer takes file lines: it takes plain records ("specs": data set,
reaction, target effect, abundances), and `Schema::Build` makes those from a validated format-3 file. A test compares
them, member by member and bit by bit, with what the old file gives: 183 of the 186 converted files are equal and 3
are recorded differences (files that exist to test a refusal). What is left before the engine runs a format-3 file
directly is the wiring of that path into the program and the comparison of the outputs (in work).

**State and extension points** (`ARCHITECTURE.md` §9). The only singletons are constant tables; whatever an
evaluation writes is owned by its caller; an extension is one explicit table row, never a self-registering class.
Done so far in Phase 4: the values an evaluation computes live in one buffer (`Results`) that the caller owns and
passes in, instead of in the data points (a debug mode poisons the buffer, so a value read before it is written fails
a test), and what the objective computes (profiled normalisations, χ² per data set, penalty totals) lives there too,
so an evaluation writes nothing into the data; the model expressions (sum, ratio, scale, sample transform) and the lab → c.m. frames have one
implementation each in a standalone library `experiment/`, with unit tests against closed forms; what a fit parameter
is (a level energy, a width, a normalisation, an abundance …) is a record given where the parameter is created, and
no code recognises a parameter by the spelling of its name any more (`parameter_kinds.md`; a check fails on a new
name test), which is what lets a data set's parameters carry the user's names.

**`CNuc` is not yet split completely.** Phase 1 moved its physics into `physics/` and left the class as the owner of
reading, the parameter vector and printing: 2,907 lines. Since then its printing, its consistency checks, its
messages and the `<levels>` grammar (now in `io/`) have moved to files of their own: `CNuc.cpp` is 1,954 lines. What
remains is the parameter vector and the set-up from the level rows (`cnuc_split.md`, "Remaining steps").

**No source file may grow.** `tests/consistency/file_sizes.tsv` lists the 18 hand-written files above 1,000 lines
with their sizes; the consistency test fails when a listed file grows or another file passes 1,000 lines, and the
numbers are only lowered. New code therefore goes into small files with one job, which is also what an AI assistant
reads and changes most reliably.

Then Phases 5–9 (`ROADMAP.md`): one energy grid with the apparatus as kernels, AZURE2's own stopping, and the THM
observable (5); evaluation without copies (6); DWBA vertices for THM (7); pyazr 3.0 with wheels (8); the web
interface, and removal of Qt, ERYA and numcmc (9).

**An idea outside AZURE2 (Jakub):** an experiment logbook — yields, target scans, on-line analysis — as a separate
small Python package on the same `experiment/` and stopping libraries, exporting named targets and data sets to a
project. AZURE2 itself stays an R-matrix code. Not before Phase 5 (it needs the stopping module).

**`dev` and the upgrade.** `dev` has 18 commits since the upgrade branched. They are merged when the reconstruction
is complete, not before. Most are documentation; the target-effect windows, the UDR kernels and test fixes port as
they are; three changes (the heavy-ion Coulomb fallback, the coincident-level limits, pyazr NUTS sampling) are redone
in the new structure; two are obsolete (the assessment is in `HANDOFF.md`).
Until then, new physics is best written on the upgrade branch, following the recipes.

**Fits and the number of threads.** A fit's path depended on the number of threads (T2 after 200 MIGRAD calls: χ²
691.6 at 1 thread, 595.9 at 8). Two causes were found (`fit_threads.md`). The first is fixed: a working copy kept a
sub-point grid that lay within 10⁻¹² of the new one, so its χ² depended, in the last bits, on what it had evaluated
before. The second is in the vendored Minuit2: its call counter is not atomic, so with several threads a fit limited
by `--max-calls` can run one more iteration. Two atomic increments remove it (shown on a scratch build); the patch
touches `external/` and waits for the maintainers' decision. Calculations are identical at any thread count.

**Two values for the maintainers to decide**, since changing them moves results: the electron mass in the
β-decay Fermi function (CODATA 2002) and the proton mass in `ECIntegral.cpp` (older than CODATA 2018).
