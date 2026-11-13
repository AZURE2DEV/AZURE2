# AZURE2 upgrade: rules

§1 and §2 are in force on `feature/upgrade` now (they are copied into the "Upgrade" section of `CLAUDE.md` on this
branch). §3 replaces the "Skills" section of `CLAUDE.md` when Phase 0 rewrites the skills.

## 1. Shape

1. Every class or function is a physical or experimental object or stage, and lives in `physics/`, `experiment/`,
   `inference/` or `io/`. Plumbing (threads, scratch, profiling) goes in `support/`.
2. Modules use only the modules below them (`ARCHITECTURE.md` §2). Each is its own CMake library target; a wrong
   dependency must fail to compile, never be worked around.
3. Physics stages are functions on arrays (Jπ blocks hold E[λ], γ[λ, c]); no objects that point at each other.
4. New physics enters at one extension point (`ARCHITECTURE.md` §9), changes one symbol of the three lines, and
   appears in one of the test projects T1–T3. No other place takes new behaviour.
5. One implementation of each formula and file format, in C++ (Coulomb functions: pycoulomb). A Python copy needs a
   test against the core. Library readers are the exception: they live in pyazr `sources/`, run at project setup,
   and the C++ core reads only the project.
6. The core computes in one stated angular-momentum and phase convention; conversions happen only in
   `level_parameters`.
7. Every physics function is exposed in the binding with array inputs, so the web UI can plot it.
8. Front ends (CLI, pyazr, web UI) own no physics and parse no files.
9. No network access in the engine. Every imported number carries its provenance (library, record, version, checksum).

## 2. Performance and threads

1. Setup state (System, measurements, grid, tables, kernels) is read-only after the fit stage starts.
2. Code in a parallel region writes only to its own output slice or per-task scratch. Shared objects are `const`.
3. No caches beyond the per-stage tables. A new cache needs a written reason, a complete key, and two tests: changing
   any input misses; a long run stays bounded in memory.
4. A perf change states time and peak memory before and after on the T2 benchmark (1 and 8 threads) and passes TSan
   and ASan at 8 threads.
5. Results do not depend on the thread count (fixed-order reductions).

## 3. Tests and skills

1. Every feature appears in T1, T2 or T3, and has an identity check (no data) wherever one exists.
2. A reference value changes only with an independent reason written in the commit (unchanged from `CLAUDE.md`).
3. Skills teach how to use AZURE2; the Sphinx docs say how it behaves; `CLAUDE.md` says how to change it.
4. A `SKILL.md` has at most 500 lines (hard limit 800); detail goes in `references/*.md` of at most 300 lines.
5. No commit hashes, job IDs, session narratives, quotes or "fixed on" notes in skills. A version note only as
   "builds before vX.Y / YYYY-MM-DD" with `expires: YYYY-MM-DD` at most 90 days ahead.
6. One owner skill per task; others link to it. Repo and user-level skills never share a name or description.
7. A description is at most 600 characters, names its triggers, and says when to use a sibling skill instead.
8. No personal absolute paths (`/groups/...`, `/Users/...`); use environment variables.
9. `tools/check_skills.py` enforces 4, 5, 6, 8, 10, 11 and that docs links resolve, in pre-commit and CI.
10. No skill names a removed pyazr name (the column "earlier name" of the migration table of
    `docs/source/user_guide/pyazr.rst`), except in a "builds before" version note.
11. A script of a skill is a tool: argument parser, `--help`, a docstring saying what it reads and writes, exit codes,
    pyazr's public names only. Every `scripts/*.py` is run by the skill's `scripts/selftest.py` (PASS/FAIL lines,
    exit 1 on a failure, 77 without the engine, under 120 s at one thread; longer: `selftest_slow.py`, label
    `slow`), registered in ctest by `tests/skills/CMakeLists.txt` as `skill_<name>`.

The skill set (Phase 0):

| Skill | Use it to | Replaces |
|---|---|---|
| azure2-setup | build a project: ENSDF levels, AME masses, EXFOR data, conventions, channel radii, first checks | azr-project-builder, nds-explorer, parts of azure2-eval |
| azure2-calculate | observables at any E and θ, extrapolation, S(0), rates, decomposition | parts of azure2-eval |
| azure2-fit | free/fixed, staged fits, nuisances, priors, multiple starts, convergence, editing levels and background poles | parts of azure2-eval, azure2-legacy-params, azure2-drop-dataset |
| azure2-check | decide whether a result can be trusted: pulls, residual structure, χ²/ν, limits, Wigner limits, radius and background-pole sensitivity, grid convergence, thread independence, frames and units | scattered "golden rules" |
| azure2-uncertainty | Hessian, profile likelihood, sampling on the objective, bands on observables and rates, systematic axes and their combination | data-covariance, parts of azure2-eval |
| azure2-plot | publication figures: excitation functions with pulls, angular distributions, S factor with bands and decomposition, rate ratios | plotting notes |
| azure2-targets | materials, solid/implanted/gas targets, stopping and straggling, resolution, Doppler, stopping-scale nuisance | target notes |
| azure2-thm | THM projects: HOES and OES data, vertex models, DWBA, kernel, normalization | THM notes outside the repo |
| azure2-resonance-search | gated search for missing or misassigned levels | repo resonance-search (renamed) |
| azure2-exchange | conventions, GNDS/ENDF export, rmfit and SAMMY parameters | ferdinand-gnds |

`rmfit-campaign` stays separate, trimmed to its procedure (its 1,700-line refinement log goes; ~15 rules lifted).
General user-level skills (uncertainty-audit, prc-figures, data-auditor, physics-referee) are linked, not copied.
The user-level `~/.claude/skills/azure2` (LUNA-Skills) must narrow its description so it does not load next to the
repo skills.
