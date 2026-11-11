# The AZURE2 upgrade: where it stands, and how to add to it

A short account for the owner of AZURE2, 2026-11-11. The details are in `HANDOFF.md` (newest first), `ROADMAP.md`
and the design notes of `docs/upgrade/`. Nothing of this is on `dev` or `main`; it lives on the local branch
`feature/upgrade` and is merged with `dev` only when the reconstruction is complete.

## Why

AZURE2's physics is right; its structure made it fragile. One object held a datum, its channel tables, its computed
values and the apparatus corrections, and parallel evaluation copied and shared it. Of the roughly 126 fixes on
`dev` between August and November 2026, about 55 % touched those objects and at least 21 repaired defects that an
earlier performance change had introduced. Every feature also had to be written several times: the χ² terms in four
code paths, the `.azr` grammar in the engine, in pyazr and in the GUI, the masses in two tables. A feature added to
one copy and not the others gives a wrong number that looks right, which is the worst bug a code used for
evaluations can have.

The upgrade rebuilds the engine so that its structure follows the calculation — levels → collision matrix →
observables → apparatus → χ² — with one implementation of each formula and each grammar, state that is owned by
whoever asks for a calculation, and a project file that says in words what each data set is.

## What is done

| Phase | What it gives | State |
|---|---|---|
| 0. Measure and protect | Synthetic projects T1–T3 with known true parameters; golden outputs; a speed benchmark; the rules in `CLAUDE.md`; skills for AI assistants | done |
| 1. Physics core | `physics/`, a standalone library (no engine headers, no GSL): Coulomb and Whittaker functions, level-parameter conventions and transformations (Brune, formal, Park), external capture, coupling coefficients, reaction pathways. Calculations no longer depend on the thread count | done |
| 2. One objective | One ordered list of χ² terms used by the fit, the gradient, least squares, the MCMC and the Python API | done |
| 3. One project model | `io/`: the project file read once, losslessly, by the engine and pyazr alike; one mass table | done |
| 4. Data are data | Named records with units; new kinds of data; the data classes replaced; `CNuc` reduced | most of it done, see below |
| 5–9 | Energy grid and apparatus kernels, own stopping, finite targets for yields; evaluation without copies; DWBA for THM; pyazr 3.0; web interface | not begun |

**Stopping.** AZURE2 has its own stopping and straggling library, `stopping/` (materials by formula, energy loss,
exit energy and range of a layer, straggling propagated through the target, each value with an uncertainty and its
source; `docs/upgrade/stopping.md`). A target says `stopping evaluation` and the engine takes the stopping of its
composition from it, per active nucleus, with one stopping-scale nuisance per material and projectile (fixed at 1, or
fitted with SRIM's stated accuracy as its prior); `pyazr.stopping` gives the same numbers in Python. Its stopping powers are
SRIM-2013's for any ion in any element, solid or gas, ported to C++ from the LunaSRIM package (they agree with it to
rounding) with SRIM's two data files, and every output cites SRIM (decision of 2026-11-11). One point is yours:
those two files stay under the terms of SRIM's authors (LunaSRIM's README reads them as free for non-commercial use),
so whether the public repository ships them or asks the user to fetch them is for you to say; the library works
either way (`stopping/data/srim/README.md`).

Phase 7 has begun as a separate library that the engine does not use yet: optical potentials, bound states,
distorted waves and the DWBA transfer to a bound state, which agree with FRESCO to better than 0.2 %
(`docs/upgrade/transfer.md`).

**How every step is checked.** Each merge is verified as a whole: the 183 project files of the test suite run in
three bases (Brune, formal, Park) and every output file is compared byte for byte with the previous head (726 + 699
+ 726 files); 182 unit and pyazr tests; the 145 regression projects with their own checks; calculations at 1
against 6 threads; and the "converted tree": every old project converted to the new format must give the old
outputs (483 files identical, 38 with a listed, explained difference such as a parameter's new name, 0 unexplained).
A number that moves is either removed or recorded with its cause and size.

## Phase 4 in plain words

**The project file (format 3).** A data set is a named record that says what it is:

```
dataset Y_curve
    reaction       15N(p,a)12C
    quantity       yield
    unit           counts/uC
    charge_state   1
    file           data/yield_curve.dat
    columns        E_lab[keV] _ value error
    norm           1 fit prior 1 3%
    target         N15T
```

Targets, beams, detectors, normalisations, energy scales, potentials and extrapolations are records too, defined
once and referred to by name. Every number has a unit (a data file may be in keV or MeV, degrees or radians, b, mb or
μb, lab or c.m. energies). Fit parameters are named after the user's records, not `segment_5_norm`.
`AZURE2 --convert old.azr new.azr --to 3` (or `pyazr.convert_azr`) turns any old file into this and reports every
change it makes; of the 197 project files of the test tree 185 convert without a question, and the rest are format 3
already or are files that exist to test a refusal.

**What a user can fit now that could not be fitted before:**

- ratio to Rutherford and to Mott; the S factor as a quantity (the engine converts with the pair's masses);
- yields per charge or per particle, linked to a named target, thick or a layer of finite thickness (a resonance scan
  through a target of known or fitted areal density; the stopping a formula, a table, or SRIM-2013's evaluation of the
  target's composition with a fitted stopping scale);
- values from the literature compared with functions of the level parameters, as χ² terms of their own
  (`<modelData>`): resonance strength ωγ, resonance energy, total and partial width, branching ratio, ANC,
  scattering length;
- extrapolations on a list of energies, a file, a number of points, a logarithmic or an adaptive grid, in the lab or
  the c.m. frame;
- a nuclear (hybrid) potential as a record of the project, with a spin–orbit term for a pair with one spinning
  particle.

The order follows what EXFOR holds for light nuclei (`data_kinds.md`). A new kind of data is one row of the table of
quantities, a small expression (cross section, sum, ratio, scale) and a test.

**For programs and AI agents.** `AZURE2 --schema` (every record kind, field, quantity, unit and diagnostic code, as
JSON), `--validate FILE [--json]` (all problems at once, each with a code, a line and a fix), `--describe FILE` (how
the engine understood each data set: reactions, rows read and kept, linked effects, parameter names),
`--print-predictions FILE` (the expression each data set is compared with). In Python the same records:

```python
p = pyazr.Project.open("project.azr")
p.add_dataset("Y_curve", reaction="15N(p,a)12C", quantity="yield", unit="counts/uC", ...)
p["Y_curve"].set(norm=pyazr.Nuisance(1.0, "fit", prior=(1, "3%")))
p.validate()
with p.session() as m:              # saves, and returns the calculation object
    chi2 = sum(m.segment_chi2())
```

pyazr spells no field, quantity or unit itself: the calls and their help text are generated from the engine's schema,
and every edit is checked by the engine's one validator, so a wrong unit or a reference to a missing target is
refused at once with the engine's message. An edit changes only the lines it touches.

**Structure.**

- What an evaluation computes lives in a buffer the caller owns, not in the data points; the model's prediction for
  a data set is one small expression tree with one forward and one adjoint implementation (it was written in about
  20 loops, twice more in the gradient, and the residual formula 13 times).
- The three data classes (`EData`, `ESegment`, `EPoint`: 10,665 lines in three files) are **replaced**, not patched.
  What the set-up fixes is one immutable `DataSet`, shared by every thread through a const pointer (the compiler
  refuses a write after set-up); what an evaluation changes is a small flat `WorkingCopy`, copied with `=` (11.9 MB
  instead of 18.6 MB for the largest test project); a point is one entry of a flat array (`Unit`, 480 bytes; it was
  608), which still carries the apparatus fold as member functions until Phase 5. The code is in files with one job
  each (the largest has 1,788 lines; it was 5,457), and the old class names are kept out by the consistency test.
- `CNuc.cpp` went from 2,907 to 1,548 lines: printing, checks, messages, the `<levels>` grammar and the reaction
  pathways are files of their own. The pathways of a run are one immutable table built at set-up and shared by all
  threads, read by every observable and by external capture; the six classes of the old pathway trees are deleted;
  the capture analyzing-power terms are built at set-up too (a lazily built table behind a mutex is gone).
- State rule (`ARCHITECTURE.md` §9): the only singletons are constant tables; whatever an evaluation writes is owned
  by its caller; an extension is one explicit table row, never a self-registering class.
- No source file may grow: `tests/consistency/file_sizes.tsv` lists the 15 hand-written files above 1,000 lines; the
  consistency test fails when one grows or another passes 1,000 lines. New code goes into small files with one job,
  which is also what an AI assistant reads and changes most reliably.

## What users notice

- Old `.azr` files still load in the engine today; pyazr (3.0.0.dev0) already refuses them with a message that
  names the converter. Every project of the test tree is format 3, so the next branch makes the engine refuse an old
  file in the same words and deletes its old reading path. Decision taken with Jakub: no backward compatibility
  inside the engine, a converter instead.
- pyazr needs the compiled engine module to read a project: there is one reader, the engine's.
- No Qt GUI on this branch: a local web interface built on the same records and schema replaces it (Phase 9),
  checked against a 201-item parity list taken from the Qt code.
- Results are unchanged except where a decision was taken and recorded. Those that change a number or refuse a
  project that used to run:
  - the MCMC profiles THM normalisations as the fit does (it corrected a diverging chain);
  - the gradient of an E1-only or E2-only sum with a component was wrong (0.29 and 1.0 of the largest element
    against central differences) and is right now: fits with `--use-gradient` of such segments change;
  - `azure2.ec_integrals` in the formal basis used the wrong energy; it gives other (correct) numbers;
  - a project with the same final state on two photon pairs was read without a message and gave a wrong capture
    cross section for the first pair (a factor 2.2 in the test); it is now refused, naming both pairs;
  - a point whose exit pair no pathway of its entrance reaches is refused (the last decay was read before);
  - external capture with a beta pair or a pair of another nucleus numbered below the photon pair read the wrong
    pathways (a crash in the test); it now finds the photon pair's pathways by the pair;
  - in format 3 a data run with an unreadable data file is refused (the old engine skipped the segment with a
    warning); an extrapolation run does not need the data files.
  - data given in the lab frame are compared with the model in the lab frame (the model divided by
    dΩ_lab/dΩ_c.m., the datum and its error as read), as the PI decided: χ² and fits are the same to rounding.
  - two resolution kernels on one data set (two of: a target integration, a beam spread without a target, a beam
    profile, a UDR) are refused; before, one of them was not applied, without a message (a target integration
    with a beam spread is one kernel and works as before);
  - a THM data set that names an experiment the `<thm>` section does not define is refused at reading
    (`E_THM_EXPERIMENT_UNKNOWN`);
  - a format-3 `<thm>` section that still lists its members the old way (`segments=`, `weight[k]=`) is refused at
    reading (`E_THM_OLD_LIST`): in format 3 the data sets name their experiment;
  - a laboratory-frame differential point at the kinematic limit (dΩ_lab/dΩ_c.m. ≤ 10⁻⁶) is refused at set-up
    (`E_LAB_JACOBIAN_AT_LIMIT`); without the refusal such a point gave χ² of 10¹⁷.

## Speed

Measured before Phase 4 (T2, 340 points with every apparatus kernel; `tests/upgrade/bench/`):

| T2, mean of three runs | before Phase 1 | after Phase 3 | change |
|---|---|---|---|
| calculation, 1 thread | 0.93 s | 0.90 s | −3 % (within the spread) |
| calculation, 8 threads | 0.42 s | 0.24 s | −43 % |
| fit, 200 calls, 1 thread | 63.2 s | 56.2 s | −11 % |
| fit, 200 calls, 8 threads | 36.4 s | 30.7 s | −16 % |
| fit peak memory, 8 threads | 488 MB | 643 MB | +32 % |

Each Phase 4 branch was measured against its own start and stayed within about ±2 % in CPU cycles. The whole of
Phase 4 has not yet been measured against the Phase 3 head in one sitting. Speed is not the aim of Phases 1–4
(structure is); Phases 5 and 6 are where it is won.

Phase 5 so far (2026-11-11; 1 thread on a machine shared with other work, so read the ratios, not the seconds):

| | before | now |
|---|---|---|
| `natural_neon`, one calculation: time (against the code at the branch point with `dev`) | 41.2 s | 8.0 s |
| `natural_neon`, one calculation: peak memory (same two builds) | 4.52 GB | 1.46 GB |
| 13N fit with `--use-gradient`, 200 calls (against the build before the gradient branch) | 50.7 s | 26.9 s |
| 13N fit with `--use-lm` | 19.5 s | 12.8 s |
| `natural_neon` fit with `--use-gradient`, 30 calls | 219.6 s | 105.7 s |
| `natural_neon` fit with `--use-lm` | 291.1 s | 139.3 s |
| T2 fit, plain MIGRAD, 200 calls | 39.9 s | 37.2 s |

The first two rows were measured by the director (both runs one after the other, χ²/N of the first segment 10.9941
in both); the fits are the gradient agent's series, read by the director from its timing files and not repeated. The
gain is one level-matrix solve per distinct compound energy instead of one per data point or sub-point. A plain
MIGRAD fit of a project without target integration gains little, as T2 shows.

## How to add a feature now

`CLAUDE.md` (section "Recipes") lists, for each kind of addition — a field of a record, a quantity, a unit, a new
kind of data set, a datum that depends only on the level parameters, a kind of fit parameter, a χ² term, a physics
function, a `<config>` key, a command-line flag, a diagnostic code — every place to touch. Each recipe has a check
in `tests/consistency/repo_consistency_test.py` that fails when it is bypassed: a second parser of the project file,
a χ² term outside `src/Objective.cpp`, a second copy of the model pass or of a physical constant, a flag missing
from `--help` or the docs, an unregistered test. `scripts/check.sh --quick` runs these in under a second, then the
build and the unit tests; `scripts/check.sh` adds the regression projects. An AI assistant working in the repository
reads `CLAUDE.md` first, so it is given the recipes and is stopped by the checks. **Ask of a new feature: does
`scripts/check.sh` pass, and which test would fail if the feature were wrong?**

## What is in work and what comes next

In work now (2026-11-11, night): the sub-points of a folded data point stop being objects of their own and become
entries of one node table per compound nucleus (the old point mapping deleted; every output must stay the same: the
structural half of the shared energy grid); one generic kind of fit parameter for the parameters of apparatus
effects, so that a kernel's width or a target's stopping scale is fitted without new engine code; the engine
refusing old files, with its old reading path deleted (the close of Phase 4).

**Do the checks still test what they say? (merged 2026-11-11, evening.)** A regression check that edits a project's
text can go on passing after the conversion to format 3 because its pattern no longer matches anything and it
compares a run with itself. All 126 checks were run with logging stand-ins for the text tools and the engine
(`tools/edit_audit.sh`): 1,658 edit places, 40 that change nothing, 36 of them by design and 4 vacuous, in two
checks; a third check tested nothing because the engine ignores a `<parameterSettings>` row that names no
parameter. The three are rewritten and each fails with an injected error. A guard keeps it so
(`tests/lib/edit_guard.sh`: an edit that changes nothing fails the check; a rule of the consistency test; a recipe
in `CLAUDE.md`). Two defects of the test harness were found in the same batch: 33 pyazr tests ended as passed when
the engine module could not be imported, and one was registered so that it never ran in a separate build
directory; both fixed.

**A crash that a full verification did not find.** The head merged at 19:25 on 2026-11-11 crashed in a
least-squares fit (`--use-lm`) or a covariance band of a project with no free parameter. The whole verification
had passed: no test ran that case. The agent that wrote the gradient found it in its timing series and fixed it
with a test (`upgrade_energy_task_gradient`); the fix was merged the same evening. The lesson is the one of the
audit above: the verification compares what the tests cover. The director then wrote the missing check
(`tests/fit_modes`: the four minimizer routes, each with the uncertainty band, on projects with no free parameter,
with level parameters, with level and target parameters, with target parameters only). On its first run it found a
second crash, **older than the upgrade**: a MIGRAD fit with `--covariance-band` of a project without free
parameters ends in a segmentation fault inside Minuit2 (HESSE on a state with no variable parameter; the code is
the same on `dev`). The engine now skips HESSE there. The check also showed two things about the band a user should
know: it propagates the covariance of level energies, widths and the THM coherent background only, so **a fitted
target parameter (areal density, stopping scale) does not widen the band of a yield**, and a fit whose only free
parameters are such ones gets no band. The message said "this fit produced no usable covariance"; it now says
which parameters the band takes. Whether fitted target parameters should enter the band is a question for you.

Merged on 2026-11-11: the data classes replaced (`DataSet`, `WorkingCopy`, `Segment`, `Unit`); **all 183
convertible project files in format 3**, the THM projects included; validation that reports what a run would
refuse, and a refused project that ends with its refusal as the last line; pyazr's functions on a format-3 project,
with edits of levels, `<config>` and `<thm>` done by the engine; the old pathway trees deleted (external capture
reads the pathway table; a crash on one ordering of the particle pairs fixed); the design note of Phase 5
(`energy_grid.md`) with its decisions taken, and its first branches: the level matrix solved once per compound
energy instead of once per data point or sub-point (a fit evaluation of 13N 132 ms instead of 317, of the
thick-target project 264 ms instead of 1,102, every output unchanged); the analytic gradient on the same per-energy
tasks, with its sums in a fixed order, so that a fit gives the same numbers at any thread count (see "For the
maintainers to decide"; the old gradient routes deleted); the yield of a layer of finite thickness as a
quantity (its accuracy over a narrow resonance on a thick layer is 1.5e-4, the present quadrature's); laboratory
data compared in the laboratory frame (last-digit rounding in 23 projects); SRIM-2013's stopping powers in
`stopping/`, ported from LunaSRIM (agreement with LunaSRIM to the printing precision), and the engine taking its
stopping from there: a target says `stopping evaluation`, with one stopping scale per material and projectile as a
nuisance parameter with a prior (checked by the director on Si₃N₄ against an own integral with LunaSRIM: 8e-5 on a
thick layer, 2e-4 of the plateau at the back edge of a thin one); ERYA removed; the apparatus effects as kernel
functions in one table (next paragraph).

**What an apparatus effect costs now (Phase 5, branch 6).** The folds (target, beam spread, beam profile, UDR,
Doppler) are no longer branches inside the engine's point class: each is a kernel of `experiment/` in a file of its
own, listed in one table (`experiment/src/Kernels.cpp`) with its format-3 keyword, its parameters and the interval its
grid must span; the engine calls whatever the table holds, and the gradient reads the same code. Every output is
unchanged. As the test of that claim a new effect was added by the written recipe (`CLAUDE.md`, "A new apparatus
kernel"): a rectangular energy resolution of full width w (`energy_resolution rectangular width 20 keV` on a
detector). It took one new kernel file (69 lines), one line in the table, one in the CMake list, tests (cases in
`experiment/tests/kernels_test.cpp`, a fixture line, a test project `tests/rectangular_resolution/`), the regenerated
schema and three documentation rows -- and no engine file. The yield of a finite layer, added a week earlier before
this table existed, had touched 11 engine files. **The limits of that count** (director's review, 2026-11-11):
(1) it holds for a resolution kernel that stands alone and has scalar parameters; a kernel that needs a new datum of
the point (as Doppler needs the temperature) also edits the two parameter blocks of `Kernels.h` and the engine
function that fills them (`TargetEffect::KernelParametersFor`); (2) kernels do not compose yet: the new kernel on a
data set that also has a target, a beam spread, a beam profile or a UDR is refused, so it cannot fold a
target-integrated yield until the shared grid (Phase 5, branch 7) chains kernels; (3) its width is declared fittable
in its row but cannot be fitted: a parameter a fit varies still needs its own kind of fit parameter, and a new kind
is not cheap -- the stopping scale of branch 5 edited 20 existing files. One generic kind that reads the kernel
table's rows is the planned remedy.

pyazr reads format 3 only since 3.0.0.dev0 (2026-11-11): `pyazr/azrfile.py` (3,756 lines of a second grammar in
Python) and `datasets.py` are deleted, every pyazr test runs on the tree's format-3 projects, and a session or
`Project.open` on an old file says how to convert it.  Then, to close Phase 4: the engine refuses old files and its
old reading path is deleted; the rest of `CNuc` (the parameter vector and the set-up from the level rows).

Phases 5–9 (`ROADMAP.md`): one energy grid with the apparatus as kernels, AZURE2's own stopping module, yields on
finite layers, the comparison of energy-shifted data in the lab frame, and the THM observable (5); evaluation without
copies (6); DWBA vertices for THM (7); pyazr 3.0 with wheels (8); the web interface, and removal of Qt (ERYA is gone) and
numcmc (9).

**An idea outside AZURE2 (Jakub):** an experiment logbook — yields, target scans, on-line analysis — as a separate
small Python package on the same `experiment/` and stopping libraries, exporting named targets and data sets to a
project. AZURE2 itself stays an R-matrix code. Not before Phase 5 (it needs the stopping module).

**`dev` and the upgrade.** `dev` has 18 commits since the upgrade branched. Most are documentation; the
target-effect windows, the UDR kernels and test fixes port as they are; three changes (the heavy-ion Coulomb
fallback, the coincident-level limits, pyazr NUTS sampling) are redone in the new structure; two are obsolete (the
assessment is in `HANDOFF.md`). Until the merge, new physics is best written on the upgrade branch, following the
recipes.

## For the maintainers to decide

- **Fits and the number of threads.** Calculations are identical at any thread count. Since the gradient branch
  (2026-11-11) so are the results of a fit: six fits (T2 with the analytic gradient, MIGRAD and Levenberg-Marquardt,
  13N, a fitted stopping scale, a fitted layer) give the same parameters, covariance, χ² and data output at 1 and 3
  threads and from run to run (before, two runs of one build at 5 threads differed in 8 files). What still varies
  is book-keeping: the number of calls printed in `fitStatus.out` (205, 204, 203 for one fit), the order of the
  "Iteration" lines in the log, and the intermediate snapshot `param.fit`. The cause is in the vendored Minuit2: its
  numerical-gradient loop calls the objective from several threads and its call counter is not atomic, so a fit
  limited by `--max-calls` could in principle stop one iteration apart. Two atomic increments remove it (shown on a
  scratch build, `fit_threads.md`); the patch touches `external/` and is not applied (the PI's decision); the engine
  could count the calls itself instead.
- **Two constants.** The electron mass in the β-decay Fermi function (CODATA 2002) and the proton mass of the M1
  capture integral (older than CODATA 2018) are kept at their values so that no result moves; they are named
  constants in `include/Constants.h` (`electronMassOfFermiFunction`, `protonMassOfM1Capture`) and can be changed in
  one place.
- **The refusals listed under "What users notice"**: each makes a project that ran before stop with a message.
