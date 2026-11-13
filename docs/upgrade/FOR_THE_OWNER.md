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

Phase 7 is a separate library, which the THM distortions of the engine now call (below): optical potentials, bound states,
distorted waves and the DWBA transfer to a bound state, which agree with FRESCO to better than 0.2 %
(`docs/upgrade/transfer.md`). It now also computes the transfer into the unbound x + A system (FRESCO's continuum
bins, within 0.5 % of FRESCO for the deuteron), and shows that inside the channel radius this amplitude is the
level's vertex times the resonance's energy dependence, which is what the THM fit assumes. The distorted waves of a
transfer may now carry the spin-orbit term of published optical potentials (it also gives the beam's analysing
power), and the two other Trojan horses of THM experiments, ³He and ¹⁴N, agree with FRESCO to better than 0.1 % for
transfers to bound states (FRESCO's own step refinement moves its ³He results by up to 0.2 %; the ¹⁴N cases use test
potentials and do not symmetrise the two ¹²C). Transfers of ⁶Li to continuum bins are not shown to 1 %: a second code,
Frescox, has the numerics of FRESCO there (`transfer.md` §23).

The engine's THM distortions are now this library's. The distortion factor R(E) (`distortion=coulomb|optical`) no
longer has a solver of its own: the waves, the Trojan-horse tail, the sum over partial waves and its plane-wave limit
come from `physics/transfer`, and so does the distorted-wave vertex (`vertexModel=dw`): the projection of its
source on the channel's harmonics is one library function (`SourceProjection.h`, `transfer.md` §25) with its own
tests (a closed form with plane waves, the library's `EntranceVertex` built from it, the engine's former numbers);
the engine keeps only the energy grid, the directions and the 2 × 2 Gram matrix, and the vertex moved by less than
3 × 10⁻¹¹ when the quadrature moved (the speed of its check +7 %). When the waves moved, the numbers moved by less
than one part in 10⁸ (|M|² at most 2.9 × 10⁻⁹, the vertex tables 6.7 × 10⁻¹⁰); no project of the test suite uses a distortion, so no output file changed. Refining the old numerics
shows that they were converged to about 5 × 10⁻⁷ for an optical potential (step and matching radius each 3 × 10⁻⁶,
partly cancelling), which the new code shares. The seven global optical sets (`opticalAA=ancai06`, ...) are one table
of the library: adding a set is one row in `physics/src/transfer/GlobalOptical.cpp` and its test points (recipe
"Add a global optical set" in `CLAUDE.md`); their spin-orbit terms are still left out.

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

The engine's physics functions (penetrability, shift, Coulomb and Whittaker functions, kinematics, Wigner limit, ...)
are `pyazr.physics.<name>` on numpy arrays, each a row of one table in the binding (`api/src/pybind_physics.cpp`)
that states its arguments, units and domain; a new one is a row and a test.

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

- **Old `.azr` files are refused.** The engine and pyazr read project format 3 only (decision taken with Jakub: no
  backward compatibility inside the engine, a converter instead). Every command that reads a project -- a run,
  `--validate`, `--describe`, `--print-predictions`, a pyazr session -- refuses a format-1 or format-2 file before
  reading anything of it, exit code 3, with one message:

  ```
  ERROR: T2.azr is a format-1 project (<segmentsData> at line 53). This AZURE2 reads format 3 only.
         Convert it:  AZURE2 --convert T2.azr T2_f3.azr     (or pyazr.convert_azr("T2.azr", "T2_f3.azr"))
         The converter keeps T2.azr unchanged and reports every change it makes.
  ```

  The converter (`AZURE2 --convert OLD NEW`, `--to 3` no longer needed; the format-2 output is gone) is the only
  reader of an old file; `AZURE2 --convert-params old.sav new.sav --project OLD.azr` (`pyazr.convert_params`)
  renames a parameter file of the old project to the conversion's names. User guide: "Coming from an older
  AZURE2" (`docs/source/user_guide/older_projects.rst`).
- A `<parameterSettings>` row that names no parameter of the project is ignored with a warning
  (`W_PARAMETER_SETTING_UNKNOWN`, its line and name); it was ignored without a word.
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
  - a SUM or RATIO whose capture reaction appears only as a later term (no data set or extrapolation of the
    project has that capture as its own reaction) ran from the command line **without external capture and without
    a message**: in the test the capture term was 0.54 to 2.9 times the right one, depending on the energy. The
    switch looked at the first term of each reaction only; the same is in `dev`. Every term counts now
    (`tests/component_mapping`, case (d)). pyazr always had external capture on and was not affected;
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
  - a point given at a lab angle with a resolution or a target effect (a target, a beam spread, a detector
    resolution) is folded at the lab angle (decision D5, 2026-11-12, observables branch): inside the window the c.m.
    angle and dΩ_lab/dΩ_c.m. follow the energy, where before every energy of the window took the point's own c.m.
    angle and Jacobian.  The largest move in the test tree is T2's (15N(p,α) at 150° lab, a 2 keV beam spread):
    2.6 × 10⁻⁵ of a model value, 1/1900 of the datum's error.  The move grows with the window's width times the slope
    of the cross section across it (and the slope of the angle): for a 300 keV resolution 0.4 MeV above the 7Li(p,n)
    threshold at 30° lab the old fold was 2.2 % off (`tests/lab_frame_fold`, checked against a brute-force integral of
    the unfolded cross section).  Elastic scattering moves only by rounding (its c.m. angle at a fixed lab angle does
    not depend on the energy without relativistic kinematics).  Where the window, or a fitted energy shift, reaches
    an energy at which the lab angle is beyond the kinematic maximum, the point fails with a warning that names the
    cause and the remedy (the model is NaN and χ² is 10¹⁰⁰); the file is not refused.
  - the `<config>` key `useAdaptiveGrid` is gone (the PI, 2026-11-12). With it at 0 a folded point (a target, a beam
    spread or profile, a detector resolution) was integrated on a uniform grid of its own, which misses a resonance
    narrower than its step; now every folded point is integrated on the mesh, which places nodes on every level and is
    never coarser than the uniform grid of the same number of points. A file that still has the line is refused with
    a message that says to delete it (`E_CONFIG_KEY_REMOVED`); the converter drops it from an old file and, when it
    was 0, says what changes. The reaction rate always takes the tabulated grid (the key at 0 chose a per-temperature
    quadrature, which `tests/reaction_rate` had recorded as close to its reference as the tabulated grid; it stays
    only as the fallback of a grid with fewer than three points of σ > 0). The check of the PI's condition, a 9.7 eV level of 14N(p,γ)
    in a 3.2 × 10¹⁹ at/cm² layer (`tests/shift_sqrt_subpoints/models/narrow_sqrt.azr`), first yield:

    | grid | default | 8 times the points | 64 times |
    |---|---|---|---|
    | mesh (the key at 1, now the only grid) | 2.18022e-6 | 2.17984e-6 | 2.17983e-6 |
    | uniform (the key at 0) | 2.06e-8 | 6.00e-8 | 5.09e-7 |

    Projects with the key at 1 give the same numbers as before, bit for bit. The ones with it at 0 move to the
    converged fold: in the tree, `examples/c12c12_tumino2018` (12C+12C THM data with a beam profile) moves by up to
    57 % in a folded value, χ² 111.661 → 117.516 (117.516 with 8 times the points), and the test model of the
    uniform grid (`shift_sqrt_subpoints/models/dense`, now deleted) by 9.4 × 10⁻⁷.

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

The shared energy mesh (merged 2026-11-12; one calculation, 1 thread, best of three runs, time and peak memory):

| project | before the mesh | with the mesh and shared sites |
|---|---|---|
| `natural_neon` | 6.2 s, 1390 MB | 1.3 s, 340 MB |
| `thick_target_resonance` | 0.91 s, 388 MB | 0.59 s, 272 MB |
| T2 | 0.26 s, 138 MB | 0.26 s, 134 MB |
| 13N | 8.3 s, 225 MB | 7.7 s, 122 MB |
| `stopping_yield` | 0.68 s, 333 MB | 0.11 s, 105 MB |
| `layer_yield` | 0.30 s, 200 MB | 0.10 s, 113 MB |
| `hybrid_potential` | 1.7 s, 338 MB | 1.1 s, 106 MB |

Against the code at the branch point with `dev`, `natural_neon` went from 41 s and 4.5 GB to 1.3 s and 0.34 GB. Fits:
the gradient fit of the stopping-yield test 43 s → 8 s, of 13N 27 s → 18 s, of T2 36 s → 34 s. The gain comes from
two things together: one mesh of energies per compound nucleus, so that data sets folded over the same energies use
the same level-matrix solves, and one stored calculation site per mesh node for all the points of a data set that
read it. The mesh without the shared sites was slower than the code before it (7.0 s and 1.7 GB on `natural_neon`);
it was never merged alone.

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

### Measured: what three new kinds of data cost (the extension drill, 2026-11-12)

The claim "easy to extend" was put to a test. An AI assistant that had never seen the engine was given only
`CLAUDE.md` and asked to add three kinds of data that did not exist, of rising difficulty, and to log every place
where the recipe was enough, silent or wrong. The three features are merged; each has a test project with an
independent answer (`truth.py`) and was shown failing with an injected error.

| What was added | `experiment/` (the table of data kinds) | `physics/` | engine (`src/`, `include/`) | API | Python | time | recipe steps: enough / silent / wrong |
|---|---|---|---|---|---|---|---|
| A scale of a computed observable: the modified S factor S* = σ E exp(2πη + gE) of heavy-ion fusion (`quantity modified_sfactor`, `shape <g> 1/MeV`) | 8 files, +49 −11 | 3 files, +20 −2 | 5 files, +27 −16 | +15 | 0 | 19 min | 2 / 4 / 1 |
| A form the table declared and the engine refused: the Legendre expansions `rsl`, `2l2`, `4pi`, and `absolute` for extrapolations | 8 files, +133 −12 | 0 | 4 files, +22 −4 | +3 −3 | 0 | 18 min | 1 / 4 / 0 |
| A value the R-matrix does not compute: a cross section relative to a tabulated standard (`quantity ratio_to_reference`, `reference <file> unit <u>`) | 12 files, +221 −4 | 0 | 7 files, +68 −7 (40 of them one new leaf file) | 0 | 0 | 17 min | 1 / 3 / 1 |

(Hand-written lines; tests, generated schema and docs apart. Python needed no line in any lane: a new quantity is
written and read through the generic record interface.)

What this says, without rounding up:

- **Another normalisation of something already computed is a table row.** No engine file.
- **Another scale of an observable is about 50 lines and two or three engine files.** Not zero: the factor is applied
  in the model pass and written in the output columns, and the S-factor conversion is still written in two places
  (`src/CoulombReference.cpp` and the column writer of `src/DataOutputFiles.cpp`).
- **Another value from outside the R-matrix is about 150 lines**, and it is the dearest lane for reasons of design,
  not of the recipe: the leaf's value is chosen by a hand-written conditional in `src/ModelPass.cpp`; its data travel
  in `ReactionSpec`; a folded point's nodes can be reached only through the point object; a refusal that depends on
  how far a fold reaches cannot yet be made at set-up (a fold or shift that leaves the reference table gives a
  not-a-number point, documented, instead of a refusal); and every output column sits in one file. The node table
  (merged the same night) removes the third; the shared energy mesh, whose supports are known at set-up, makes the fourth
  possible; the first is a table waiting to be written.
- **The recipes were silent more often than they were enough** (11 silent, 2 wrong, 4 enough over the three lanes).
  Every silent or wrong step is now a sentence of the recipe or a check of the consistency test (check 35 is one).
  The checks that stopped real mistakes during the drill: `experiment_quantities`, `convert_spec_equivalence`,
  `experiment_prediction`, `experiment_build`.

Open from the drill: the three non-default Legendre definitions follow `docs/upgrade/data_kinds.md` and were not
checked against the LEXFOR text itself; differential data relative to a reference are not built; pyazr's
`calculate_sfactor` knows only the plain S factor.

## No hidden state (four steps, 2026-11-12 and 13; not finished)

The aim: every number the engine prints depends only on the project, the command line and a short list of
environment switches (`docs/upgrade/evaluation.md` section 2), not on what ran before in the same process or on how
many threads ran it.  What is done: no per-thread cache is left and no "print this warning once per process" flag.
What the engine keeps between calls for speed now lives in objects that a run owns: the scratch of the energy tasks
(in the evaluation's result buffer), the run's per-thread scratch of the external-capture integrals and the THM
vertex (in the run's `Config`), and the warnings printed once per run.  `tests/pyazr/session_isolation_test.py` runs
three projects in one process one after another, in both orders, at 1 and 3 threads, and checks that each gives the
bytes of a fresh process.  The second step removed the global
pointer to the run's settings and the amplitude cache the Python binding swapped per call: each particle pair is built
with the run's settings of its functions (relativistic or not, which Coulomb table), and the external-capture term
table, Park's shift-derivative memo and the `--no-ec-table` amplitudes belong to the run.  Two sessions may now
evaluate at the same moment on two Python threads: `tests/pyazr/concurrent_sessions_test.py` (slow, about 11 minutes)
runs pairs of sessions at once 20 times, and two different Coulomb table files one after another, and each gives the
bytes of a fresh process; the engine before this step failed it.  A project's relative paths (data, output,
checks, tables, intEC.dat, files named at a prompt) are now the project file's directory's, or the directory a pyazr
session names, never the process's working directory, and pyazr no longer changes it: two sessions in two directories
run side by side, and the command-line program started from anywhere reads and writes the files it would read and
write started in the project's directory (`tests/run_from_elsewhere`).  This is a change of behaviour of the
command-line program: before, relative paths were the working directory's.  A fit and the MCMC make a copy of the
working data only when no idle one exists (one copy less in memory), and the signature of `intEC.dat` names the
Coulomb table by its content, so a file written with another table is recomputed, not reused.  The fourth step: every
environment variable the engine reads is one row of one table (`support/src/Environment.cpp`, 39 switches, each with
its kind: changes numbers, changes cost, prints, test hook), read once when the process starts, and the manual's list
(`command_line.rst`, "Environment switches") is checked against it; no other code reads the environment.  The
switches that can move a number are copied into the run's settings (`Config::environment`), so making one of them an
option is one line; `AZURE2 --describe` ends with the switches that are set.  The counters of the diagnostics
(profiling, the Coulomb table's lookups, the gradient checks) belong to the run.  A variable set inside a running
process (pyazr's `os.environ`) is no longer seen: set it before the process, or the first session, starts.  To add a
switch: the recipe "Add a diagnostic switch" in `CLAUDE.md`.  What is not done: the fitter's and the sampler's own
members, the MCMC's progress globals.  Every one is a row of
`tests/consistency/allowed_state.tsv` with the step that removes it (53 rows; 15 more are permanent: constant
tables, the Coulomb table's content identity, locks of output streams, and the bytes of input files, kept per
process under a key of path, size and modification time, which `cache_stats()` counts: say if you want them per
session).  What you
must not do: add a `static` or `thread_local` variable that keeps a value, a global pointer, or a `mutable` member; a
consistency check refuses them.

## Radiative widths as a total width and branchings (merged 2026-11-13)

A level that decays by several γ rays can now be written as its total radiative width, its branchings to the final
states and the mixing ratios, instead of one partial width per final state and multipolarity: a record
`gamma_decay` in a section `<gammaDecays>` of the project file, which a fit varies and reads directly, with measured
branchings as priors. The engine is not changed: before every calculation the record's numbers are turned into the
partial widths and then into the amplitudes the engine has always used, by the engine's own conversion; a project
written with the record and the same project with the partial widths give the same numbers. The mixing ratio is
Krane and Steffen's (ENSDF's); its sign was fixed by comparing the engine's γ angular distributions with their
formula. The total radiative width stays out of the resonance's total width, as the radiative widths always have;
the user guide says what that means for a resonance-strength datum (a difference of Γ_γ/Γ_particle). A level whose
transitions have external capture (¹⁴N(p,γ)¹⁵O into the ground state and 6.79 MeV, ¹²C(p,γ), ¹⁷O(p,γ)) can be written
so too: the record's widths are the observed ones, external part included, and the engine recomputes the external
part at every evaluation, so they stay what the record says while a fit moves a proton width, the level energy or an
ANC. Not covered: the formal basis (`--no-brune`), refused with a message. From the record's merge until the fix of
the same day, the paths that take parameters in observed form ignored the record's numbers: pyazr's
`calculate_chi2` and `update_segments` given a new Γ_γ or branching computed the project as written, and the MCMC in
observed parameters (`Use Reduced Width Amplitudes: no`) had a likelihood that did not move with them; the grid
placement at narrow resonances read the branchings as widths in eV, but only into a sum no output used. Fixed and
tested (`tests/gamma_decay/check.sh`, `tests/pyazr/gamma_decay_test.py`); fits in the CLI and the MCMC in amplitudes
were never affected.
The ensemble sampler (MCMC, both bases) did not have the record's priors and the remainder's wall, and let a
branching leave [0, 1] and Γ_γ go negative (samples reached Γ_γ = −6.6×10¹³); it now samples −½ × the fit's
objective, every term from the same table, and rejects a point outside the record's limits (tested,
`tests/reference/sampler_objective_test.cpp`). The sampler's run (its file, step, best walker, callbacks and stop
request) now belongs to the sampler object and the call that runs it, so two samplers in one process share nothing.

## The centre of a prior is the number written (2026-11-13)

A fitted number of a format-3 project is written `<value> fit prior <centre> <width>`. For norms and the constant
energy shift the program used the centre written; for the other kinds it did not. The abundance of an isotope, the
`isotope_sum`, the screening potential, the Doppler temperature, a transmission sample's thickness, the sqrt(E) shift
and the shift function's parameters measured the prior from the start value (or from 0, or from 1 for the sum)
without a word; the areal density of a target, the stopping scale and a kernel's parameter (the rectangular width)
refused a centre different from the value. Since format 3 (2026-11-10) this cost the prior to anyone who continued a
fit from its fitted values written into the project, which is what pyazr's `save_fit` writes: the fitted value with
the original prior. A fit continued so had no prior on those numbers, and `save_fit`, which checks its own result,
refused it after any fit with a free abundance or a free areal density. The program's own way of continuing a fit (its
`param.sav` read back) kept the centre of the project and was not affected. Now the centre written is the centre for
every kind, one rule (`parameter_kinds.md` section 4): a start away from the centre pays its prior from the first
evaluation. No project of the tree had a centre different from its value where it mattered, so no number of the tree
moved. Found with it and mended: `save_fit` did not write a fitted areal density, Doppler temperature, screening
potential or sample thickness back into the project at all (it kept the old value), which also made it refuse its own
result. The old AZURE2 (`origin/dev`) has no centre to drop: its `<abundances>` line is `Z A value vary sigma` and the
prior is centred on the value on the line, so a user who copies fitted abundances into that line moves the prior with
them; through `param.sav` it was right. The grids and the limits of a fitted density, Doppler temperature, stopping
scale or kernel parameter are laid out from the centre too (when the start lies inside the range so laid out), so a
fit written back with `save_fit` and computed again gives the fit's χ²; before, the Doppler temperature's moved by
6e-6 and every layout moved with the start.

## What is in work and what comes next

Merged in the night of 2026-11-12, each after the whole verification: the
engine refusing old files, with its old reading path deleted (the close of Phase 4 except the rest of `CNuc`); the
nodes of a folded data point as entries of a node table (the old point mapping and the sub-point objects are
deleted: a data point is a `Unit`, each of its nodes a `Site`, every output identical, the fit as fast as before
and natural_neon's peak memory 2.5 % lower); the extension drill (previous section); and the planner of the shared
mesh as a pure function with its tests (`include/MeshPlanner.h`): one mesh per compound nucleus, planned from the
kernels' supports, the levels and each data set's numerical settings, so that two data sets folded over the same
resonance share their level-matrix solves. On the projects of the tree the planned mesh has at most 1.04 times
today's distinct energies and far fewer where data sets overlap (T2 0.54, natural_neon 0.37).

**The engine is on the shared mesh (merged 2026-11-12, morning). This is the first step of the upgrade that changes
numbers on purpose.** Folded yields (target integration, beam resolution, any apparatus kernel) move by about the
quadrature error they had; unfolded numbers do not move at all. What was measured before merging:

- *How far the old default was from a converged fold.* The old engine was run with its quadrature refined eightfold
  (and fourfold where eightfold is too slow). Largest relative distance of the old default to that reference: 5e-3
  on T2, 3e-3 on the thick-target project, 1e-4 on 13N, and 5 % on `hybrid_potential`, a project whose own settings
  are too coarse in both engines (a finding about that project's input, not about the mesh).
- *That every changed number is at least as good.* Of 577 compared output files (Brune basis) 542 are identical, 33
  moved and are as close to the refined reference as before or closer (thick-target project 3.2e-3 → 4.2e-4,
  `natural_neon` 7.8e-4 → 2.9e-4, the stopping-yield fit 4.2e-3 → 2.0e-3; T2 and 13N unchanged in size), 2 belong to
  `hybrid_potential` (6 % against 5 %), none failed. The table per project is in `docs/upgrade/energy_grid.md` §7.1.
  **The acceptance rule is ours (confirmed by Jakub on 2026-11-12) and should be confirmed by you**: a file passes when its largest distance to the
  refined reference is at most 1.2 times the old engine's and its root-mean-square distance at most 1.03 times
  (the 3 % is the resolution of the comparison).
- *What it does to fits.* Eight test fits: fitted parameters move by at most 0.14 of their fitted error (an areal
  density, 7e-5 relative), and, checked on that fit, towards the result of the refined engine.
- *Stored expectations.* Six regression checks compared with numbers of the old quadrature; each new number was shown
  to be the closer one before the expectation was regenerated. One check runs with a wider tolerance (1e-5 instead of
  1e-6) until a data set with a fitted energy shift is on the mesh too (the next branch).

**Trojan-horse data on the shared mesh (merged 2026-11-12).** A THM data point is now folded like
any other: its half-off-shell cross section is computed once per mesh energy and shared by every point and data set of
the same experiment that needs it, and at an energy where direct data are folded too the level matrix is solved once
for both. The THM cross section at an energy depends on the energy, the reaction and the experiment's settings only,
never on the data point, which is what makes this possible without new physics. The THM numbers move by the quadrature
error they had and come closer to the converged fold in every file (18O(p,α) THM: largest distance 4.2e-4 → 2.7e-4);
T2's direct data move a little too, because the THM points join the same mesh. One stored number changed: an option kept
for old projects (`kinematics=lambda32`) divides by the entrance momentum, which vanishes at threshold, so its χ² rests on
the quadrature there (5605 → 5562; the converged value is about 5507). The gradient of a THM fit is now ten times
cheaper (18O(p,α): 88 → 9 ms), with the same agreement with finite differences as before. One defect of the engine
was found on the way and is fixed: with the analytic gradient (`--use-gradient`), the component of a fitted energy
shift on a THM data set with a profiled normalisation could hold only the derivative of the shift's prior (a fit of
the 18O(p,α) THM data with the level energy and the shift free ended "invalid" with a wrong error of the energy; it
converges now). The `<config>` key `useAdaptiveGrid`, which when off gave a folded point a uniform grid of its own
(the last path beside the mesh), was removed afterwards (below, "What users notice").

**Four mechanisms deleted (merged 2026-11-12).** The automatic tolerance, the fold windows, the edge blending and
the σ(Ē) fallback existed only to make a point's private grid affordable. On the shared mesh a node costs the same
whoever uses it, so they are gone: a folded data set names its effect and, optionally, the mesh's settings, and every
point of it is folded with the whole kernel. In numbers: 44 engine, format and converter files, 874 lines removed
against 161 added; the set-up of target effects 968 → 762 lines; `fold_windows` is no longer a field of the format (a
file that writes it is refused with a message that names the field and says to remove it; an old file converts with a
note that names the line). Every project of the test tree keeps its bytes except the one that used a window
(`beam_spread_ay`: χ² 3838.66 → 3760.95, one A_y point by 1.5 of its error): its new numbers are, bit for bit,
what the old engine gives with the window widened to all energies, so the move is the points the window had left
without the beam spread. **This is a behaviour change you should know of**: a window was a restriction of the physics,
not an approximation, and an analysis that relied on it gets the full fold. Speed and memory of the measured projects
are unchanged.

**Energy shifts on the fixed mesh (merged 2026-11-12).** A folded data set with a fitted energy shift (constant,
sqrt(E), a shift function, a follower) used to keep nodes of its own that moved with the shift, so every change of the
shift recomputed the R-matrix at all of them, and a cache with locks existed for those nodes. Now the nodes stay and
the kernel moves: the mesh is planned for the data set's energies widened by the reach of the shift (three prior
widths), and a shift inside the reach costs a fold and no solve. In numbers, on the synthetic project T2: the
level-matrix solves of one gradient 158,929 → 48,526, peak memory 300 → 259 MB, one calculation 0.23 → 0.28 s (the
reach adds nodes). Only the projects with a folded shifted data set change (T2, its truth file, the five models of
`shift_sqrt_subpoints`), each toward the old engine with an eightfold refined quadrature. **Three things you should
know of.** (1) The fitted shifts of T2 move by 0.13 of their errors, onto the refined engine's fit (within 0.01 σ).
(2) When a fit walks a shift beyond the planned reach, the mesh is planned again and χ² steps by the quadrature
difference there (measured 1.3e-2 in 10,500); a prior on the shift sets the reach. The old engine had steps of that
size inside the range a fit visits. (3) A thick- or thin-target fold of a shifted point read the stopping power at
the node's energy minus the shift, as the old rigid move did, so that numbers stayed; the PI decided on 2026-11-12 to
read it at the true energy (the paragraph below). Two defects of the first version were found in the review and fixed with tests before
the merge (a stopping formula called at negative energy; an error stop for shifted data with external capture).

**An energy shift is the true beam energy (branch of 2026-11-13).** A folded data set with an energy shift took the
stopping power, the target thickness and the depth window at the measured (unshifted) beam energy and only moved them
with the shift, so a yield on a resonance was off by (d ln ε/d ln E) × shift/E, 0.9 % for a 10 keV shift at 270 keV;
now a shifted point is folded exactly as the same point without a shift at the shifted energy.  Your files change only
where a data set has both an energy shift and a target integration with an energy-dependent stopping: in the test tree
T2's truth file (0.21 % at most, χ² 330.43 → 330.84) and the `shift_sqrt_subpoints` models (−0.89 to −0.97 %); a
fitted shift on such data moves accordingly.
Not moved to the shifted energy, and yours to decide: the beam-profile kernel's window and detailed-balance norm, a
user-defined resolution kernel's table, and the 5 σ support of a Gaussian beam spread whose width depends on energy
(about 1e-7). None of them reads the stopping; whether they should follow depends on what a shift means for those
data (a calibration of the beam energy, or of the reconstructed energy).

**Six follow-ups (merged 2026-11-12).** What a user sees of them: a yield on a target layer thinner than 1 eV is
σ(E)·N (the file was refused before, with a misleading message; with an energy shift on the data set it is still a
failed point, in work); a command-line flag that takes a value and gets none or a bad one is an error with the usage
line and exit code 2, where it was a warning; the Jacobian's shift column of a summed data set with a shift of its own
moves every term (it moved the first only). No number of the test tree moved.

**Observables without the angle (merged 2026-11-12).** The engine stores, per mesh energy and reaction, the
observables without the angle (integrated cross section, Legendre coefficients, elastic amplitudes, the analysing
power's coefficients: one table of kinds), and a data point's angle enters only in its own row. Twenty angles at one
energy share one sum over pathways, and two experiments on the same reaction share everything up to the row. Numbers
moved by rounding only, with one named physics change: a point given at a lab angle is folded in the lab (the list
of behaviour changes above; T2 moves by 2.6 × 10⁻⁵). Measured: the project with the most angles per energy
(`beam_profile_kernel`) runs in 16 s instead of 5 min 27 s and in 134 MB instead of 913 MB at one thread, and the T2
fits are about 30 % faster. A new observable is six named places (`CLAUDE.md`, the recipe), a new ratio observable
needs no fold code.

**A second batch of small corrections (merged 2026-11-12).** Behaviour you may meet: a target layer thinner than 1 eV
is calculated with an energy shift on its data too; a fit with the analytic gradient or least squares and a free
energy shift on a thin layer no longer ends the program; an unknown command-line option is an error (exit code 2),
`--use-nlopt` in a build without NLopt included; an analysing power on a β-delayed or photon entrance is refused
when the file is read (it was calculated as zero); a project whose energy mesh would need more than a million nodes
(a wide prior on the energy shift of a very thin layer) is refused with a message that names the data set and what to
change, where it ran out of memory before. No number of the test tree moved. Known and left for a later branch: on
layers below 1 eV the fitted error of an energy shift is too small (χ² is a staircase in the shift there), and the
cost of a shifted data set grows with the width of the shift's prior over the width of the layer.

In work now: the diagnostic switches read once into one record (the last step of "no hidden state"; no number may
change); the distorted-wave vertex's quadrature moved into the transfer library; the `gamma_decay` record in the
places that handle observed level parameters (the defect named in the section on radiative widths). Merged on
2026-11-13: the `gamma_decay` record on a level whose transitions have external capture; the third step of "no
hidden state" (section above); the stopping of a data set
with an energy shift is taken at the shifted beam energy (numbers of such data sets on a target change, by the
logarithmic slope of the stopping times shift over energy); the total radiative width, the gamma branchings and the
mixing ratios of a level as fit parameters (previous section); the Trojan-horse distortion factor and the waves of
the distorted-wave vertex computed by the transfer library, with the global optical potentials as one table (the
unit tests' quantities moved by at most 1e-8; no project of the tree has a distortion). Merged in the earlier batch: one generic kind of fit parameter for the
parameters of apparatus effects (see "How to add a feature now").

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
target-integrated yield until the shared grid (Phase 5, branch 7) chains kernels; (3) a fitted parameter of an apparatus record
cost a kind of its own -- the stopping scale of branch 5 edited 20 existing engine files.

**What a fitted apparatus parameter costs now** (`docs/upgrade/apparatus_parameters.md`). One parameter kind,
`apparatus`, carries every number of a target, beam or detector record that a fit varies: the stopping scale and a
kernel's own parameter.  A kernel parameter declared fittable in its table row is written in the project with the
words every fitted number already uses (`width 14 keV fit prior 14 2`), is one parameter for every data set that
names the record, and goes through one engine file for its slot in the vector, its hand-over, its prior and its
outputs.  The test of that claim: the rectangular width is fitted (pseudo-data made at 20 keV: 20.026 ± 0.013 keV;
without noise 19.999998 keV), the stopping scale moved onto the kind with every output of its fit unchanged, and a
triangular resolution with a fitted width was added by the written recipe: one kernel file (68 lines), two lines in
the table, one CMake line, tests, the regenerated schema, three documentation rows -- **0 engine files, about 6
minutes** (before: 20 engine files).  The limits: the fit takes central differences of the whole χ² for it (the
kernel's analytic derivative is tested but not used yet); a number that reaches the calculation other than as a
kernel's value is a new row of the kind's table (one function file, one table line); the density, Doppler
temperature and sample thickness keep their own kinds (moving them would reorder every project's parameters).

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
UDR kernels and test fixes port as they are; three changes (the heavy-ion Coulomb
fallback, the coincident-level limits, pyazr NUTS sampling) are redone in the new structure; four are obsolete, among
them the two on automatic target-effect windows, which the upgrade deletes (the
assessment is in `HANDOFF.md`). Until the merge, new physics is best written on the upgrade branch, following the
recipes.

## The web interface: the page (2026-11-13)

`python -m pyazr.ui project.azr` prints a URL; the browser shows the project. The page has a header (project, unsaved
mark, revision, Undo, Redo, Save, the engine's total χ²) and the tabs Project (file, the validator's problems with
links, `<config>`, the session's options read-only, history, text, engine), Nucleus (pairs, levels and channels: a
width typed in changes the χ², Undo takes it back; γ-decay and potential records) and Data (every section of the file,
the χ² per data set, one form per record generated from the engine's schema with the engine's description beside it,
the first rows of the data file, a preview plot of data and model). Fit, Results, Physics and THM are placeholders
naming the step that brings them. Nothing is hand-written per record kind except which fields share a line; a new
value type of the schema needs one widget. plotly.js 4.1.2 (the cartesian bundle, MIT licence) is vendored under
`pyazr/ui/static/vendor/`: the page loads nothing from the network, and without plotly the plots are tables. It was
checked with node's test runner and a headless Chrome on this machine; no person has used it in a real browser yet.
User guide: "The web interface".

The Fit tab (2026-11-13, second step of the page) runs the recipe fits (iminuit's MIGRAD with HESSE, scipy's
L-BFGS-B) and iminuit's MINOS as jobs from a form drawn from the server's job table, follows a running job (calls,
objective, its plot against the call, Stop, Kill), writes an accepted fit into the project as one Undo-able edit and
shows the χ² of the project as written beside the fit's; the table of every parameter edits the working x of the free
ones. The engine's own fitters and the sampler are listed with the pyazr call they wait for.

The Results tab (2026-11-13): for every data set and extrapolation, the data times the normalization in use, the
model and the engine's residuals under it, against c.m. energy, excitation energy or angle, as measured or as S
factor; the χ² table equal to the program's `chiSquared.out`, model data, angular distributions at one energy, and
the program's output files written on request. The numbers equal the program's `AZUREOut_*` files (tested on T1,
T2); the model of a data set is drawn at its points, a dense grid waits for pyazr.

## For the maintainers to decide

- **Fits and the number of threads.** Calculations are identical at any thread count. Since the gradient branch
  (2026-11-11) so are the results of a fit: six fits (T2 with the analytic gradient, MIGRAD and Levenberg-Marquardt,
  13N, a fitted stopping scale, a fitted layer) give the same parameters, covariance, χ² and data output at 1 and 3
  threads and from run to run (before, two runs of one build at 5 threads differed in 8 files). What still varies
  is book-keeping: the number of calls printed in `fitStatus.out` (205, 204, 203 for one fit), the order of the
  "Iteration" lines in the log, and the intermediate snapshot `param.fit`. The cause is in the vendored Minuit2: its
  numerical-gradient loop calls the objective from several threads and its call counter is not atomic, so a fit
  limited by `--max-calls` could in principle stop one iteration apart. Two atomic increments remove it (shown on a
  scratch build, `fit_threads.md`); the patch touches `external/` and is not applied (the PI's decision).
  Minuit2's call counter is not thread-safe under its own OpenMP gradient loop; not patched; AZURE2 counts itself
  (2026-11-13, branch `feature/upgrade-8-gamma-branchings`): `fitStatus.out` and the log give the objective's own
  atomic count of MIGRAD's calls, equal to Minuit2's at 1 thread and the same at any thread count; `--max-calls`
  still acts on Minuit2's counter, so under threads that limit is approximate by a few calls.
- **Two constants.** The electron mass in the β-decay Fermi function (CODATA 2002) and the proton mass of the M1
  capture integral (older than CODATA 2018) are kept at their values so that no result moves; they are named
  constants in `include/Constants.h` (`electronMassOfFermiFunction`, `protonMassOfM1Capture`) and can be changed in
  one place.
- **The refusals listed under "What users notice"**: each makes a project that ran before stop with a message.
