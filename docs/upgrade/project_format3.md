# Project format 3: experimental data, targets and model data — design

A fit compares measured numbers with model values. For that comparison the engine must know, for every data set,
what the numbers are (quantity, unit, frame of each column), which reaction they belong to, what apparatus they passed
through (target, beam, detector), and which nuisances (normalization, energy scale, areal density) they carry. Today
this is written as positional columns, integer codes and lists of segment indices (§1). This note designs the
replacement, "format 3": one self-describing record per data set, target, beam, detector, shared nuisance,
extrapolation and model datum, with named fields, and with every closed choice taken from one table that the engine
exports. It gives the grammar (§3), the same records as pyazr calls, GUI forms and agent operations (§4–§6), the
table itself (§7), the converter from today's files (§8), and the decisions left to the director (§9). It changes no
code. The engine reads format 3 only; old files go through `AZURE2 --convert` (maintainer, 2026-11-09,
`data_model.md` §5).

**What is built (step 3b, 2026-11-10).** The table, `--schema`, `--validate` (§6.1, 6.2, 7); the converter (§8.1–8.3);
`Schema::Build` (records → the engine's set-up specs, model data included; §8 construct by construct, checked by
`convert_spec_equivalence` over every converting test file); the engine runs a format-3 file
(`Config::ReadConfigFile`: records → Validate → `<levels>` tables → Build → set-up; §7.3) and `--describe` (§6.3, a
plain listing of the specs, not yet the expression tree of S10); parameter names from the records (§8.4's renaming:
the converter's `parameters` map). Not built: §8.4's refusal of old files and `--convert-params` (step 3d), pyazr's
model on records (step 3c), units other than the engine's (refused, `E_NOT_BUILT`), the quantities and grids marked
"a later step" in the table, `Measurement` (S7: the fixed follower parameters of a shared normalization remain).
Differences from the old engine on the same project (the gate of §8.5, classified; `reference/data_formats.rst`): an
unreadable data file refuses a run with data (`E_DATA_FILE`; the old engine skipped the data set with a warning), a run
without data does not read them (`W_DATA_FILE_NOT_READ`); a follower of a varied shared normalization always has its
fixed parameter `<record>.value[<data set>]`; scattering lengths are reported as model data; abundance-weighted
outputs name the sample by its target; an angle-integrated data set reports the angle window 0–180°.

Inputs: `data_model.md` (Phase 4: `Measurement`, model expression trees, decisions M1–M13) and `data_kinds.md` (what
the engine accepts, EXFOR, the physics of each kind, K1–K9). Line numbers refer to `feature/upgrade-4-data` at
`eb79d95a`. **[V]** read in the code or a file; **[I]** inferred, not run. A small survey script (1 thread, in the
session scratchpad, not committed) counted constructs in the 189 `.azr` files under `tests/` and `examples/`; its
counts are quoted as "survey".

## 1. What today's grammar does, and why it must go

| construct | where | what is wrong with it for a user, a GUI or an agent |
|---|---|---|
| segment line: 8 positional numbers, then `J l` only when isDiff = 2, then norm, then a shift block of 3, 1 or 0 tokens decided by whether tokens look like numbers | `include/SegLine.h:104-106, 146-174` [V] | the meaning of a token depends on the tokens before it; a file name that parses as a number was once read as a shift (`SegLine.h:146-149`) |
| observable as an integer `isDiff` whose meaning differs between data (`ESegment.cpp:38-83`) and test lines (`ESegment.cpp:166-196`): 3, 4, 5 mean different things | `pyazr/datasets.py:77-98`, `pyazr/azrfile.py:2563-2568` copy both tables [V] | a number, not a word; two tables; pyazr holds its own copies |
| composite block `1 op n e x θ [s]` with three layouts told apart by counting numbers | `include/ComponentList.h:30-88` [V] | reactions as pair keys; a short list silently becomes a plain segment with a warning (`ComponentList.h:71-75`) |
| UPOS block `isUPOS L Ic δ` after the composite block | `SegLine.h:199-210` [V] | positional; a missing field silently resets to 0 (`SegLine.h:205-208`) |
| trailing keywords; an unknown one ends the line silently | `SegLine.h:382-384` [V] | a typo drops every later keyword without a message; test lines skip data-only keywords silently (`include/ExtrapLine.h:104-131`) |
| `<targetInt>` line: 9 positional fields, a parameter count, then optional fields detected by peeking at the next character | `src/TargetEffect.cpp:61-173` [V] | a field is optional only by position; the target has no name |
| `<targetInt>` segment list `"3,6,7-14"`, by segment index | `src/EData.cpp:1500-1525` [V] | adding or moving a segment line silently re-targets the effect; the GUI renumbers lists on every move (`qt_parity.md` §6.1) |
| a test segment k gets the effect of a `<targetInt>` line naming data segment k | `EData::MakePoints` (`src/EData.cpp:688`, effects at `:958`) [V]; 32 projects of the survey have both test segments and an active `<targetInt>` line | a positional coincidence that changes extrapolations; `tests/natural_neon/README.md` switches `<targetInt>` off in `check.sh` for this reason |
| fit parameters named by position: `segment_<k>_norm`, `density_<n>` (n = active line), `screening_<n>`, `doppler_<n>`, `abundance_<sample>_<Z>_<A>` | `src/CNuc.cpp:535-537`, `src/AZURE2.cpp:1171-1212` [V] | `prior_cov.txt` of T2 names `segment_7_norm segment_8_norm` (`tests/upgrade/T2/prior_cov.txt:2`): a priors file that breaks when a segment is inserted |
| code that recognises a parameter's kind by its name prefix | `src/AZUREParams.cpp:9-15, 123-125`, `src/AZURECalc.cpp:1431`, `include/ParameterLabel.h:27-48`, `src/ParameterLimitsManager.cpp:245, 304` [V] | with user-chosen names (format 3) this cannot work; §8.4 |
| group sharing by label, "first segment in file order owns the value, the others' values are ignored" | `SegLine.h:227-259` (shiftgroup, normgroup, shiftfngroup), `TargetEffect.h:95-96` (densitygroup) [V] | order carries meaning; ignored numbers stay in the file and look live |
| `<abundances>` per **compound nucleus** (Z, A of 16O for p+15N) | `data_formats.rst` (abundances block), `tests/natural_neon/natural_neon.azr` [V] | a natural target is described by its isotopes, not by the compound nuclei they form |

The grammar is read by four hand-written routines (`SegLine`, `ExtrapLine`, `ComponentList`, `TargetEffect`), each
with its own tokenizer, plus the `<thm>`, `<abundances>` and `<scatteringLengths>` readers. Every new keyword so
far was code in one of them and a line in pyazr's `_KEYWORDS` (`pyazr/datasets.py:120-122`).

## 2. Principles

What "intuitive in the GUI, in Python and for an agent" means here, as rules a reviewer can check:

1. **One record per thing, with a name.** A data set, a target, a beam, a detector, a shared normalization, a shared
   energy scale, an extrapolation and a model datum are each one record `kind name` followed by named fields. Nothing
   is identified by its position in a list or in a line; adding or removing a record never changes the meaning of
   another. References are by name.
2. **One spelling everywhere.** A field is called the same in the file (`charge_state 1`), in the pyazr call
   (`charge_state=1`), in the GUI form (field id `charge_state`, label from the table), in the JSON schema, in
   diagnostics, and in parameter names (`run1.norm`).
3. **One table.** Every closed choice — record kinds, fields, quantities, units with factors, frames, column roles,
   norm modes, target kinds, stopping and straggling models, Legendre expansions — is a row of one table compiled into
   the engine (§7). The engine validates files with it; the binding exports it; pyazr validates arguments and builds
   its signatures from it; the GUI builds its drop-downs and conditional fields from it; the reference docs are
   generated from it; agents read it as JSON. A front end that spells a quantity name itself fails a consistency check
   (§7.4). A new quantity is one row, a tree builder and a test.
4. **Reactions by what they are.** `15N(p,a)12C`, `15N(p,a)12C[4.439]`, `12C(p,g)13N[total]`; sums and ratios
   written as `+` and `/` of reactions. Pair keys never appear in format 3.
5. **Units stated, never implied.** Every dimensioned number carries its unit token; every column of a data file has
   a declared role with its frame (`E_lab`, `theta_cm`). There is no default unit and no default frame.
6. **Loud.** A required field that is missing, an unknown field, a field that does not belong to the quantity, a
   duplicated field, a reference to an unknown name: each is an error naming file, line, record and field. Nothing is
   ignored. A default that can change a number does not exist; the defaults that remain (quadrature settings,
   `active yes`) are listed in the table and written to the log when used.
7. **The engine explains itself.** `--describe` prints, for every record, the expression tree, units and the
   parameters it creates, before anything is computed (§6).

Not principles but consequences: the file stays the lossless `io::Document` of Phase 3 (`docs/upgrade/project.md` §2;
`io/include/azure2/io/Document.h:20-80`); a pyazr or GUI edit replaces only the lines of the fields it changes.

## 3. The format

### 3.1 Records, not keyword lines (decision F1)

A record is a header line `<kind> <name>` followed by field lines `<key> <value tokens>`, inside a `<section>`. It
ends at the next header line or at the section's closing tag. Indentation is not significant (writers indent four
spaces). `#` starts a comment anywhere outside quotes. Tokens are separated by blanks; a token with blanks is quoted
(`"data/run 1.dat"`). Keys are lower case with underscores; record kinds and field keys are disjoint sets, so a line is
a header or a field by its first token alone.
**[Correction, step 2 (2026-11-10): the note's own records break both rules. `norm`, `target`, `beam` and `detector` are record kinds and
also data-set fields, and `J` is an upper-case key. Implemented: kinds and keys are disjoint *within a section*
(checked by `experiment_schema`), each section is read with its own kinds, and the shared-normalization record is
renamed `normalization` (it lives in `<data>` beside the field `norm`); keys are `[A-Za-z][A-Za-z0-9_]*`.]**

Why records and not today's one-line keyword lists:

| criterion | one line per data set, keyword blocks (`data_kinds.md` §5.1, K1) | a record, one field per line (chosen) |
|---|---|---|
| hand editing | T2's segment 2 is already 186 characters (`tests/upgrade/T2/T2.azr:55`); with quantity, unit, frame, target and beam blocks it doubles | a field per line; comments per field |
| diffs | any change rewrites the whole line | the diff is the changed field |
| lossless document, local edits | an edit rewrites a line that holds 15 fields | an edit replaces one line; `Document::ReplaceLines` already does this (`Document.h:58-63`) |
| error messages | "line 55, token 23" | "T2.azr:58: dataset s02_diff090, field angle_window: …" |
| parser | arity of every keyword must be known to find the next one (today's silent stop, `SegLine.h:382-384`) | split at line ends; each field's value grammar is a table entry; one routine (§7.3) |
| logbook tool, spreadsheet, agent | must reproduce token order | writes `key value` lines in any order; a CSV row maps to a record by column headers = field keys (§6.5) |

TOML or YAML for the whole file was considered (F1 (b)): standard parsers exist, but none keeps comments and layout for
in-place edits from C++ without a second document model, and the `<levels>` and `<thm>` sections are not in that
form. A record block inside today's sections keeps one document model and one tag rule.

### 3.2 Sections

| section | record kinds | replaces |
|---|---|---|
| `<data>` | `dataset`, `normalization` (§3.1 correction), `energy_scale` | `<segmentsData>`, the group keywords |
| `<extrapolations>` | `extrapolation` | `<segmentsTest>` |
| `<targets>` | `target`, `beam`, `detector` | `<targetInt>`, `<abundances>`, `transmission`/`yield` sample thicknesses |
| `<modelData>` | `datum` | `<scatteringLengths>`; new: ωγ, E_r, Γ, branching, ANC (M10) |
| `<config>` | `format=3` and the existing keys | `format=2` (`io/src/Convert.cpp:30-38`) |

`<levels>`, `<thm>` (without `segments=`, §3.9), `<parameterSettings>` and `<mcmc>` keep their grammar in this
note. Names are resolved after the whole file is read, so records may appear in any order and in any of the sections
listed for their kind.

### 3.3 Values

- **Names**: `[A-Za-z][A-Za-z0-9_-]*`, unique across the project (one namespace for all records). `.` and `@` are
  reserved: `.` joins a record and a field in a parameter name (`run1.norm`), `@` names a target epoch (`C13A@run2`).
- **Numbers** are copied as written: the reader keeps the token text, and a writer never re-prints a number it did not
  change (the bit-identical gate depends on it, §8.5).
- **A number with a unit**: `0.2 1.4 MeV` — the unit token follows the numbers it applies to. Unit tokens are rows of
  the unit table (§7.2); `µb` is an alias of `ub`.
- **A nuisance** (normalization, shift term, areal density, temperature, U_e, abundance, efficiency, stopping scale):

  ```
  <value> [<unit>] fixed|fit|profiled [prior <centre> <width>[%]] [grid_max <value>]
  ```

  **[Correction, step 2 (2026-11-10): T2's THM data sets write `norm profiled` without a value (§3.9, §3.11). Implemented:
  `<value> [<unit>] fixed|fit [prior …] [grid_max …] | profiled`, `profiled` alone.]**

  The mode is required. `prior` always states its centre (today's `prior_centre` rows of `<parameterSettings>` become
  the centre, `data_formats.rst:189`); a width with `%` is relative to the centre. `fit` without `prior` is a free
  parameter. `profiled` is allowed only where the table says (THM norms). `grid_max` is the largest value the
  sub-point grids must span (today's `densityfit … max`, `TargetEffect.cpp:232-240`, and `doppler … Tmax`).
- **A reaction expression**: `term { + term }` or `term / term`, `term := [<number> *] <reaction> [@ <angle> <unit>
  <frame>]` (§3.5).
- **A formula**: a quoted expression in `x` with named parameters, `formula "a0*x^a1" a0=1.6e-20 a1=-0.45` (read by
  today's `Equation` class).

### 3.4 `dataset`

```
dataset <name>
    reaction       <reaction expression>                    required
    quantity       <quantity>                                required (§7 rows)
    unit           <unit of the value column>               required
    file           <path>                                    required
    columns        <role> <role> ...                         required (§3.4.1)
    energy_window  <lo> <hi> <unit>                          optional: rows outside are not used
    angle_window   <lo> <hi> <unit>                          optional; only for quantities with an angle column
    norm           <nuisance> | <name of a norm record>      required
    shift          offset <nuisance> | sqrt <nuisance> | <name of an energy_scale>   optional, repeatable; terms add
    covariance     <path>                                    optional (format: data_formats.rst:68-120)
    target         <target name>[@<epoch>]                   per quantity (§7)
    beam           <beam name>                               optional
    detector       <detector name>                           optional
    quadrature     points <n> width_multiplier <w> per_width <p>        points required with a fold; w, p default 20, 50 (logged)
    fold_windows   <lo>-<hi> ... <unit> | auto[:tol[:margin]]  blend <w> <unit>  auto_tolerance <t>   optional
    active         yes|no                                    optional, default yes
    source, exfor, doi                                       optional text, never read by the physics
    + quantity-specific fields: J l s (phase_shift), expansion (legendre), multipole (sigma, sfactor),
      secondary_gamma (UPOS), charge_state, efficiency (yield), thm (THM data), reference frame of a fixed angle
```

#### 3.4.1 Columns

`columns` names the role of each column of the data file, in order. A column's frame is part of its role:

| role | meaning | unit |
|---|---|---|
| `E_lab` | beam energy of the light particle on the heavy one at rest (forward kinematics) | in the role: `E_lab[MeV]` |
| `E_cm` | entrance-channel c.m. energy | `E_cm[keV]` … |
| `E_decay` | lab energy of the detected particle of a decay at rest (β-delayed) | |
| `theta_lab`, `theta_cm` | detection angle | `theta_cm[deg]`, `[rad]` |
| `order` | Legendre order k (integer) | — |
| `value`, `error` | the datum and its 1σ uncertainty, in `unit` | `unit` |
| `error_percent` | the uncertainty in per cent of the value | % |
| `window_lo`, `window_hi` | reconstructed-energy window of a beam-profile kernel (today columns 5–6, `qt_parity.md` §7) | as the energy column |
| `_` | a column that is not read | — |

**[Correction, step 3a (2026-11-10): the role `E_x` is removed. No quantity row of §7.2 takes an excitation-energy
column and the old engine reads none (data energies are E_lab, `data_formats.rst`); `frame E_x` of a `resonance_energy`
datum (§3.10) is a field value, not a column role. `columns E_x[MeV] …` is now `E_UNKNOWN_ROLE`.]**

The unit of the energy and angle columns is written in the role: `E_lab[MeV]`, `E_cm[keV]`, `theta_cm[deg]`.
The converter writes `E_lab[MeV]` and `theta_lab[deg]` (today's fixed units, `data_formats.rst:7-31`). An angle-
integrated quantity has no angle role, so the dummy angle column of today's files is `_`.

The value's frame follows the angle's: dσ/dΩ with `theta_cm` is a c.m. cross section; a lab value at a c.m. angle is
refused. This replaces a single `frame` field (`data_kinds.md` §5.1) because the commonest differential layout in
today's files has E in the lab and θ, dσ/dΩ in the c.m. (isDiff 4: 65 data and 16 test lines in the survey), which one
frame cannot describe (decision F5).

A data file may instead start with a header line naming the roles (`# columns: E_lab[keV] value error`) and the record
then says `columns from_file` — for files written by a logbook tool or a spreadsheet (§6.5).

### 3.5 Reactions

A reaction is `<target nucleus>(<projectile>,<ejectile>)<residual>[<state>]`:

- nuclei as `<A><Symbol>`; light particles `n p d t h a g` (h = ³He, g = γ); `e-`/`b-`, `b+` for β pairs [I: the label of
  a β-delayed entrance pair (type 20) needs the `<levels>` grammar's pair records and is left to that note];
- the residual state: absent or `[0]` for the ground state, `[<E_x in MeV>]` for an excited state, matched to a pair
  whose heavy-particle excitation agrees to the digits written (`[4.439]` matches 4.43891 MeV; ambiguity is an
  error naming the candidates); `[total]` for the sum over every γ pair of the compound nucleus (today exit −1,
  isDiff 3 data / 4 test);
- a light particle's excitation `a[e1]` analogously (`data_formats.rst`, column 34 `e1`);
- the reaction must name exactly one entrance pair and one exit pair of the project; otherwise the error lists the
  project's pair labels.
- **[Step 3a, director's decisions D2, D3 (2026-11-10):** the target nucleus may carry a state, `16O[2.143](n,n)16O[2.143]`
  — an entrance pair whose heavy particle is excited — matched to the pair's heavy-particle excitation exactly as the
  residual state is; `[total]` is not a target state. An elastic reaction has the same state on both sides. Where a
  label names two pairs of the project (the same particles and states, e.g. `examples/c12c12_tumino2018`: pairs 1 and
  6, ¹²C + ¹²C, which differ in their channel lines and in the THM binding-energy column), the data set or
  extrapolation writes `entrance_pair <key>` and/or `exit_pair <key>`, the pair key of `<levels>` (column 6), a
  written label as the level ID of F13 (a). They are accepted only where the label alone is ambiguous and required
  there (the reader checks both against `<levels>`: a redundant key is a second source of truth); they apply to a
  single reaction only — a sum or a ratio with an ambiguous term waits for named pairs in `<levels>`. The rule "a THM
  data set takes the pair carrying the THM columns" was not taken: the old engine uses the key written, the two pairs
  of that file carry different channels (pair 1 has eight levels pair 6 lacks), and the THM columns are read only by the
  THM vertex (`EDataThm.cpp`, `EPoint.cpp`); the rule would hold in that file but is no meaning the old engine has.]**

Ordinal labels (`a1`, `p2`) are refused: an ordinal counts the states the project happens to include, which is a
position. `a0`, `p0`, `n0`, `g0` are accepted as aliases of the ground state.

Composite expressions (today SUM, RATIO, `include/ComponentList.h`):

```
reaction  15N(p,a)12C + 15N(p,a)12C[4.439]                 # SUM, terms in the order written
reaction  15N(p,a)12C[4.439] / 15N(p,a)12C                 # RATIO, one numerator, one denominator
reaction  12C(p,p)12C / 12C(p,p)12C @ 90 deg cm            # denominator at a fixed angle
reaction  20Ne(n,n)20Ne + 0.5*20Ne(n,a)17O                 # a scaled term
```

**[Correction, step 2 (2026-11-10): the third example is a ratio of differential cross sections, quantity `dsigma_domega` in unit `1`, which
the `dsigma_domega` row of §7.2 did not allow. Implemented: that row takes `1` with a ratio reaction, as `sigma` does.]**

No parentheses and no ratio of sums: the tree has `Sum`, `Scale`, `Ratio(base, t1)` (`data_model.md` §2.3), and the
grammar writes exactly those. The first term is the base term whose frame the others are converted to (today's rule,
`ESegment.cpp:1487-1493` via `data_model.md` §1.5).

### 3.6 Shared nuisances: `norm` and `energy_scale` records (decision F6)

```
normalization D
    value      1 fit prior 1 8%

energy_scale G
    offset     0 MeV fit prior 0 0.003
    sqrt       0 MeV^1/2 fixed                               # E' = E + b sqrt(E), b in MeV^1/2 and E in MeV
    function   "a0*(x-0.8)" MeV                              # x = E_lab of the point (MeV); the result in MeV
    param      a0 0 fit prior 0 0.01                          # repeatable, one per parameter
```

**[Correction, step 2 (2026-11-10): the record kind is `normalization D` (§3.1 correction); the function states the unit of its result,
`function "a0*(x-0.8)" MeV` (principle 5); with b in MeV^1/2 the sqrt term is E' = E + b sqrt(E), E in MeV (b sqrt(E/MeV)
would make b an energy; `include/SegLine.h:214` has the same comment).]** **[Step 3a: the example above is written
with these corrections; the comment of `SegLine.h` (now `convert/include/azure2/convert/SegLine.h`) says E' = E + a +
b sqrt(E) with E in MeV and b in MeV^1/2.]**

A data set says `norm D` or `shift G`. The value lives in one place; no data set carries a copy that is ignored.
`shiftgroup G` (shared offset and function) is `shift G` on each member; `shiftfngroup F` (shared function, own
offset) is `shift F` plus the member's own `shift offset …`. The engine sums the terms in a fixed order (offset, sqrt,
function), not in file order, so E' is computed as today (`SegLine.h:211-259`).

### 3.7 `target`, `beam`, `detector`: who owns which effect (decision F7)

Every apparatus effect of today's `<targetInt>` (`src/TargetEffect.cpp:61-297`) belongs to the object that physically
causes it:

| effect | today | owner in format 3 | field |
|---|---|---|---|
| areal density N, stopping ε(E), straggling, density fit and group | positional fields, `densityfit`, `densitygroup` | `target` | `areal_density`, `stopping`, `straggling` |
| thick target (N = 0, `TargetEffect.cpp:334`) | density 0 | `target` | `kind infinite` |
| isotopic abundances of a natural or enriched target | `<abundances>` per compound nucleus, `sample <name>` | `target` | `isotope` (repeatable), `isotope_sum` |
| neutron sample thickness n (atoms/b) | `transmission n`, `yield n`, `thicknessgroup` | `target` (`kind sample`) | `areal_density … atoms/b` |
| free-gas Doppler broadening of target atoms | `doppler T [vary σ Tmax]` | `target` | `doppler free_gas <nuisance K>` |
| electron screening (depends on the target's electrons) | `screening Ue vary σ` | `target` | `screening <nuisance eV>` |
| beam energy spread, constant or σ(E) | `isConvolution σ`, `isConvCoefficients "expr"`, `sigmalab` | `beam` | `spread` |
| beam energy distribution of skewed Gaussians | `beamprofile` | `beam` | `profile` (repeatable), `truncation`, `detailed_balance` |
| detector energy resolution of a reconstructed energy (TPC σ, THM c.m. resolution) | the `beamprofile` TPC σ; a Gaussian for THM | `detector` | `energy_resolution` |
| finite-solid-angle attenuation Q_k | `isQCoefficients n Q…` | `detector` | `attenuation` |
| time-of-flight resolution (SAMMY UDR) | `udr` | `detector` | `udr` |
| sub-point grid: points, resonance width multiplier, points per width | positional | `dataset` / `extrapolation` | `quadrature` |
| lab-energy windows, auto windows, blend width, auto tolerance | quoted range token and two numbers | `dataset` / `extrapolation` | `fold_windows` |

A target referenced by several data sets is one set of parameters (one N, one U_e, one T): that is today's per-line
sharing (`TargetEffect.h:57-63`, "a screening parameter belongs to the line, shared by its copies") and today's
`densitygroup`, without the line. The quadrature belongs to the data set because it is a numerical setting of a fold
applied to that data set's points, not a property of the target.

```
target C13A
    kind           layer                                     # layer | infinite | sample | none
    active         13C                                       # the nucleus the areal density counts
    composition    C 1                                       # element stoichiometry; needed by 'stopping evaluation'
    areal_density  1.2e18 atoms/cm2 fit prior 1.2e18 15%
    epoch          run2 areal_density 1.0e18 atoms/cm2 fit prior 1.2e18 15%
    stopping       formula "a0*x^a1" a0=1.6e-20 a1=-0.45 MeV*cm2/atom   # per active atom, x = E_lab/MeV
    #  or          evaluation scale 1 fit prior 1 4%           (Phase 5 stopping, ARCHITECTURE.md §7)
    straggling     bohr                                      # none | bohr | coefficient <c> keV^1/2 cm
    phase          solid                                     # solid | gas; only with 'stopping evaluation'
    #              [Correction, step 2: refused here, since this target's stopping is a formula; fixture without it]
    isotope        13C 0.99 fixed                            # repeatable; abundance weights of Weighted
    doppler        free_gas 300 K fixed
    screening      1000 eV fit prior 1000 200

beam B1
    spread         1 keV lab                                 # constant sigma (lab), or
    #  spread      formula "a0*x" a0=0.004 lab               # sigma(E_lab) (today's sigmalab), or 'cm' (legacy)
    profile        skewed_gaussian xi 0.34 MeV omega 0.035 MeV alpha 1.5 weight 1
    truncation     0
    detailed_balance no

detector D90
    attenuation    1 0.995 0.98 0.96 0.93                    # Q_0 .. Q_4, lab
    energy_resolution 4 keV lab
    udr            file "res.udr" flight_path 20 m burst 5 ns channel 1 ns centred yes
```

**Epochs.** A target that changes during a campaign (carbon build-up, loss) has one `areal_density` per epoch; the
stoichiometry, stopping and isotopes are shared. A data set names `target C13A@run2`. The parameter is
`C13A@run2.areal_density` (`data_kinds.md` §4.6, K5).

### 3.8 `extrapolation`

The same fields as a data set except `file`, `columns`, `norm`, `shift`, `covariance`, `*_window` (refused with a
message, where `ExtrapLine.h:104-131` skips them today), and with a grid:

```
extrapolation S_low
    reaction   12C(p,g)13N
    quantity   sfactor
    unit       keV*b
    energies   range 0.02 0.8 MeV cm step 0.002
    #          range 0.02 0.8 MeV cm points 200 log
    #          list 0.025 0.05 0.1 0.2 MeV cm
    #          file "grid.txt" MeV lab
    #          adaptive 0.02 0.8 MeV cm tolerance 1e-3 max_points 2000 start_step 0.01
    angles     none                                          # or range 0 180 deg cm step 5, list …
    target     ...                                           # only when written: no positional inheritance
```

- `range lo hi step s`: E_k = lo + k·s for k = 0 … ⌊(hi − lo)/s + 10⁻⁹⌋. Today's loop accumulates `energy += eStep`
  (`src/EData.cpp:787-822`); the two differ (decision F10, §8.5).
- `adaptive`: start from the uniform grid of `start_step`, seed with E_λ ± k Γ_λ of every level in the range, and bisect
  every interval where the model differs from linear interpolation between its ends by more than `tolerance`
  (relative), until no interval is refined or `max_points` is reached (then a warning). The grid is a function of the
  current parameters, computed in the extrapolation run (extrapolations never enter the χ²); `AdaptiveIntegrationGrid`
  is reused where it fits (HANDOFF 2026-11-10, "Requested feature").
- The frame of the grid is stated (`lab` or `cm`); `cm` replaces pyazr's own conversion (`pyazr/azrfile.py:1522-1540`).

**[Implemented (brief K, 2026-11-10).** One reader of the grid forms, `experiment/Grids.h` (`ReadGridValues`), used by
`Build` and by `pyazr.extrapolation_grid`; a `cm` energy grid is taken to lab energies by `Frames`' `BeamLabFromCM`
(the inverse of the engine's map) with the `<levels>` masses of the entrance pair; a `file` grid is read relative to
the project as a data file. The refinement is a pure function, `experiment/AdaptiveGrid.h` (`RefineGrid`, given the
model as a callable; `pyazr.refine_grid` on a Python function). Its rule as built: an interval is split when
|f(m) − (f(a)+f(b))/2| > t·max(|f(m)|, 10⁻³ max|f|); the midpoint of an interval not split is dropped; seeds
E_λ + kΓ_λ, k = 0, ±1/2, ±1, ±2, ±4 (observed values at the run's parameters); `max_points` makes the splits with the
largest deviations and warns with the largest left. The extrapolation run (`src/AdaptiveExtrapolation.cpp`, called
from `AZUREMain` before the data set is initialized) computes each sweep's batch with the point pass on an
extrapolation data set made by `EData::MakePoints` (the parameter layout the run's), then makes the run's data set on
the refined grids; with an angle grid the energies are refined at all angles at once (one union grid). A pyazr session
refuses an adaptive grid. `AdaptiveIntegrationGrid` (the sub-point grid of a fold, refined from a level list and a
tolerance on the folded integral) is not reused: its rule and its inputs are those of an integral, not of a curve.
User page: `data_formats.rst`, "Extrapolation Grids".]**

### 3.9 THM data

A THM data set is the same record with `thm <experiment>` (today isDiff + 10, `ESegment.cpp:38-43`), naming an
`experiment[<name>]` of `<thm>`, or `thm global` for one without an experiment. `segments=` disappears from `<thm>`
(today `src/ThmExperiment.cpp:77-87` parses an index list); the experiment's members are the data sets that name it.
`norm profiled` is the norm mode of THM data sets in an experiment. Target integration on a THM data set stays refused
(`src/EData.cpp:1504-1515`).

### 3.10 `datum` (model data, no beam energy)

```
datum wg_L5
    quantity   resonance_strength
    level      5                       # the level ID of <levels> (column 9), a label written in the file
    J          3-                      # required: must equal the level's J^pi (catches a wrong ID)
    reaction   15N(p,a)12C             # entrance and exit of omega*gamma
    value      2.63 keV
    error      0.15 keV
    source     "author, journal volume page (year), table 2"
```

Quantities: `resonance_strength`, `resonance_energy` (`frame E_x | E_cm | E_lab`), `total_width`, `partial_width`
(`pair`), `branching_ratio` (`reaction`), `anc` (`pair`, `l`, `s`; closed channels only), `scattering_length`
(`pair`, `s`, `value` and `error` for Re, optional `imag`, `imag_error`), all ((q(θ) − d)/e)² in one
`ModelDataTerm` (`data_model.md` §4, M10). ωγ is eq. (8) of `data_kinds.md` §4.7 with observed widths at E_r; a
warning at load when Γ/E_r > 1 % or another level of the same J^π lies within 10 Γ; refused for E_r < 0 (K4).

### 3.11 T2 in format 3, beside the old text

Old (`tests/upgrade/T2/T2.azr:49-93`, sections relevant here):

```
<abundances>
natN: 8 16 0.0037 1 0.0002
natN: sum 0
</abundances>
<segmentsData>
1 1 2 0.199999 1.400001 0 180 0 1 1 5 0 0.003 1 data/s01_int_thin.dat 0 0
1 1 2 0.249999 1.350001 90 90 1 1 1 8 0 0.003 1 data/s02_diff090.dat 0 0 normgroup D shiftgroup G shiftfn "a0*(x-0.8)" 1 0 0.01 1 covariance cov_diff.txt
1 1 2 0.249999 1.350001 150 150 1 1 0 0 0 0 0 data/s03_diff150.dat 0 0 normgroup D shiftgroup G covariance cov_diff.txt
1 1 2 0.249999 1.350001 0 180 0 1 1 3 0 0 0 data/s04_thick.dat 0 0 sqrtshift 0 0.003 1
1 1 2 1.169999 1.270001 0 180 0 1 0 0 0 0 0 data/s05_implanted.dat 0 0
1 1 2 0.259999 0.420001 0 180 0 1 1 5 0 0 0 data/s06_gas.dat 0 0
1 1 2 0.277499 0.412501 0 180 0 1 1 0 0 0 0 data/s07_profile.dat 0 0
1 1 2 0.599999 1.400001 0 180 0 1 1 0 0 0 0 data/s08_sum.dat 1 0 1 1 3 -999
1 1 3 0.699999 1.400001 0 180 0 1 0 0 0 0 0 data/s09_ratio.dat 1 1 1 1 2 -999
1 1 2 0.249999 0.450001 0 180 0 1 0 0 0 0 0 data/s10_natural.dat 0 0 sample natN
1 1 2 0.032015 0.640314 0 180 10 1 1 0 0 0 0 data/s11_thm_hoes_a.dat 0 0
1 1 2 0.042687 0.618970 0 180 10 1 1 0 0 0 0 data/s12_thm_hoes_b.dat 0 0
1 1 2 0.053358 0.533595 0 180 0 1 1 0 0 0 0 data/s13_thm_oes.dat 0 0
</segmentsData>
<targetInt>
1 "1" 20 1 0.001 1 5e17 "a0*x^a1" 2 1.6e-20 -0.45 0 0 0 "" 0 1 0 10 20
1 "2,3" 20 1 0.002 0 0 "" 0 1 5 1 0.995 0.98 0.96 0.93 0 "" 0
1 "4" 40 0 0 1 0 "a0*x^a1" 2 1.6e-20 -0.45 0 0 0 "" 0 0 0 10 20
1 "5" 40 1 0.0008 1 2e+17 "a0*x^a1" 2 6.0e-20 -0.30 0 0 0 "" 0 1 0.6 10 20 densityfit 1 10% 3e+17
1 "6" 40 0 0 1 1e18 "a0*x^a1" 2 0.9e-20 -0.50 0 0 1 "a0*x" 1 0.004 0 sigmalab doppler 300
1 "7" 60 0 0 0 0 "" 0 0 0 0 "" 0 0 0 20 50 beamprofile 1 0.34 0.035 1.5 1.0 0.004 0 0
1 "8" 20 1 0.002 0 0 "" 0 0 0 0 "" 0
1 "11,12" 30 1 0.021344 0 0 "" 0 0 0 0 "" 0
</targetInt>
<thm>
experiment[LC09] segments=11,12 beam=15N target=d spectator=n Ebeam=60
</thm>
```

**[Step 3a: the block below is the note's draft; what the converter writes is `convert/tests/fixtures/T2.f3.azr`
and, as records, `experiment/tests/data/T2_format3.azr`. They differ from the draft in three places, each a defect of
the draft: s02, s03 and s08 share `beam_targetInt2` (and s02, s03 `det_targetInt2`; F7: identical beams are one
record, shared records are not named after one data set); the excited state is `[4.43891]`.]**

New (what `AZURE2 --convert T2.azr T2_f3.azr` writes; names from the data-file stems; pairs: 1 = p + ¹⁵N, 2 = α +
¹²C, 3 = α + ¹²C(4.43891); `<targetInt>` lines without straggling fields take today's grid defaults 20/50, written
out):

```
<targets>
target tgt_s01_int_thin                         # <targetInt> line 1
    kind           layer
    active         15N
    areal_density  5e17 atoms/cm2 fixed
    stopping       formula "a0*x^a1" a0=1.6e-20 a1=-0.45 MeV*cm2/atom
    straggling     bohr
beam beam_s01_int_thin
    spread         0.001 MeV lab
beam beam_s02_diff090                           # <targetInt> line 2 (data sets s02, s03)
    spread         0.002 MeV lab
detector det_s02_diff090
    attenuation    1 0.995 0.98 0.96 0.93
target tgt_s04_thick                            # line 3
    kind           infinite
    active         15N
    stopping       formula "a0*x^a1" a0=1.6e-20 a1=-0.45 MeV*cm2/atom
    straggling     none
target tgt_s05_implanted                        # line 4
    kind           layer
    active         15N
    areal_density  2e+17 atoms/cm2 fit prior 2e+17 10% grid_max 3e+17
    stopping       formula "a0*x^a1" a0=6.0e-20 a1=-0.30 MeV*cm2/atom
    straggling     coefficient 0.6 keV^1/2 cm
beam beam_s05_implanted
    spread         0.0008 MeV lab
target tgt_s06_gas                              # line 5
    kind           layer
    active         15N
    areal_density  1e18 atoms/cm2 fixed
    stopping       formula "a0*x^a1" a0=0.9e-20 a1=-0.50 MeV*cm2/atom
    straggling     none
    doppler        free_gas 300 K fixed
beam beam_s06_gas
    spread         formula "a0*x" a0=0.004 lab          # [Correction, step 2: ... a0=0.004 MeV lab, the unit of sigma]
beam beam_s07_profile                           # line 6
    profile        skewed_gaussian xi 0.34 MeV omega 0.035 MeV alpha 1.5 weight 1.0
    truncation     0
    detailed_balance no
detector det_s07_profile
    energy_resolution 0.004 MeV lab
beam beam_s08_sum                               # line 7
    spread         0.002 MeV lab
beam beam_LC09                                  # line 8 (s11, s12)
    spread         0.021344 MeV lab
target natN                                     # <abundances> sample natN; compound 16O = p + 15N
    kind           none
    isotope        15N 0.0037 fit prior 0.0037 0.0002
    isotope_sum    none                         # "natN: sum 0"
</targets>

<data>
norm D                                          # normgroup D, value of s02 (s03's 1 0 0 was ignored)
    value          1 fit prior 1 8%
energy_scale G                                  # shiftgroup G, from s02 (s03's 0 0 0 was ignored)
    offset         0 MeV fit prior 0 0.003
    function       "a0*(x-0.8)"                 # [Correction, step 2: "a0*(x-0.8)" MeV; the record is 'normalization D']
    param          a0 0 fit prior 0 0.01

dataset s01_int_thin
    reaction       15N(p,a)12C
    quantity       sigma
    unit           b
    file           data/s01_int_thin.dat
    columns        E_lab[MeV] _ value error
    energy_window  0.199999 1.400001 MeV
    norm           1 fit prior 1 5%
    shift          offset 0 MeV fit prior 0 0.003
    target         tgt_s01_int_thin
    beam           beam_s01_int_thin
    quadrature     points 20 width_multiplier 10 per_width 20
dataset s02_diff090
    reaction       15N(p,a)12C
    quantity       dsigma_domega
    unit           b/sr
    file           data/s02_diff090.dat
    columns        E_lab[MeV] theta_lab[deg] value error
    energy_window  0.249999 1.350001 MeV
    angle_window   90 90 deg
    norm           D
    shift          G
    covariance     cov_diff.txt
    beam           beam_s02_diff090
    detector       det_s02_diff090
    quadrature     points 20
dataset s03_diff150                             # as s02 with angle_window 150 150 deg, file data/s03_diff150.dat
    ...
dataset s04_thick
    reaction       15N(p,a)12C
    quantity       yield
    unit           1/particle                   # CONVERTED: the old engine compared these numbers per incident particle
    file           data/s04_thick.dat
    columns        E_lab[MeV] _ value error
    energy_window  0.249999 1.350001 MeV
    norm           1 fit prior 1 3%
    shift          sqrt 0 MeV^1/2 fit prior 0 0.003
    target         tgt_s04_thick
    quadrature     points 40 width_multiplier 10 per_width 20
dataset s05_implanted                           # sigma, b; norm 1 fixed; target tgt_s05_implanted; beam beam_s05_implanted;
    ...                                         # quadrature points 40 width_multiplier 10 per_width 20
dataset s06_gas                                 # sigma, b; norm 1 fit prior 1 5%; target tgt_s06_gas; beam beam_s06_gas;
    ...                                         # quadrature points 40
dataset s07_profile
    reaction       15N(p,a)12C
    quantity       sigma
    unit           b
    file           data/s07_profile.dat
    columns        E_lab[MeV] _ value error window_lo window_hi
    energy_window  0.277499 0.412501 MeV
    norm           1 fit                        # its prior is in prior_cov.f3.txt (s07_profile.norm)
    beam           beam_s07_profile
    detector       det_s07_profile
    quadrature     points 60 width_multiplier 20 per_width 50
dataset s08_sum
    reaction       15N(p,a)12C + 15N(p,a)12C[4.439]
    quantity       sigma
    unit           b
    file           data/s08_sum.dat
    columns        E_lab[MeV] _ value error
    energy_window  0.599999 1.400001 MeV
    norm           1 fit
    beam           beam_s08_sum
    quadrature     points 20
dataset s09_ratio
    reaction       15N(p,a)12C[4.439] / 15N(p,a)12C
    quantity       sigma
    unit           1
    file           data/s09_ratio.dat
    columns        E_lab[MeV] _ value error
    energy_window  0.699999 1.400001 MeV
    norm           1 fixed
dataset s10_natural
    reaction       15N(p,a)12C
    quantity       sigma
    unit           b
    file           data/s10_natural.dat
    columns        E_lab[MeV] _ value error
    energy_window  0.249999 0.450001 MeV
    norm           1 fixed
    target         natN
dataset s11_thm_hoes_a
    reaction       15N(p,a)12C
    quantity       sigma
    unit           arbitrary
    thm            LC09
    file           data/s11_thm_hoes_a.dat
    columns        E_lab[MeV] _ value error
    energy_window  0.032015 0.640314 MeV
    norm           profiled
    beam           beam_LC09
    quadrature     points 30
dataset s12_thm_hoes_b                          # as s11, file data/s12_thm_hoes_b.dat, window 0.042687 0.618970 MeV
    ...
dataset s13_thm_oes
    reaction       15N(p,a)12C
    quantity       sigma
    unit           b
    file           data/s13_thm_oes.dat
    columns        E_lab[MeV] _ value error
    energy_window  0.053358 0.533595 MeV
    norm           1 fit
</data>

<thm>
experiment[LC09] beam=15N target=d spectator=n Ebeam=60
</thm>
```

and `<config>` gets `priorCovariance=prior_cov.f3.txt` (a converted copy whose name row reads
`s07_profile.norm s08_sum.norm`) and `format=3`. The `...` lines above abbreviate fields identical to the record
before; the converter writes every field. Parameter names after conversion: `s01_int_thin.norm`,
`s01_int_thin.shift.offset`, `D.value`, `G.offset`, `G.a0`, `s04_thick.shift.sqrt`, `tgt_s05_implanted.areal_density`,
`natN.isotope.15N`, … (today `segment_1_norm`, `segment_2_norm`, `segment_2_energy_shift`, `density_4`,
`abundance_natN_8_16`).

What the old text left implicit and the new text states: the frame and unit of every column; that segment 4's numbers
are reactions per incident proton; which beam spread is a beam property and which (line 8) is the THM energy
resolution entered as a lab beam spread (the converter cannot tell, §8.3); that s03's own norm and shift are not used;
that `natN` is a target of ¹⁵N atoms, not a sample of ¹⁶O compound nuclei. T2 has no UDR and no screening; their
conversion is shown in §8.1 with lines from `tests/electron_screening` and `tests/beam_profile_kernel`.

## 4. The same records in Python

### 4.1 Calls

pyazr 3 holds the project as the engine's document (`_azure2.Document`, `project.md` §5) and edits records through
the binding. One call per record kind; the keyword arguments are the field keys of the table, the values the same
text as in the file, or Python values that the binding turns into that text:

```python
from pyazr import Project
p = Project.open("T2.azr")                     # refuses a format-1/2 file: "run pyazr.convert_azr(...)"

p.add_target(name, *, kind, **fields)          # target
p.add_beam(name, **fields)                     # beam
p.add_detector(name, **fields)                 # detector
p.add_norm(name, *, value)                     # shared normalization
p.add_energy_scale(name, **fields)             # shared energy scale
p.add_dataset(name, *, reaction, quantity, unit, file, columns, norm, **fields)
p.add_extrapolation(name, *, reaction, quantity, unit, energies, **fields)
p.add_datum(name, *, quantity, value, error, **fields)

p.records()                                    # {name: Record}, every kind
p["s02_diff090"].norm                          # 'D'  (the field text)
p["s02_diff090"].set(norm="1 fit prior 1 8%")  # edits one line of the document
p.rename("s02_diff090", "diff90")              # renames the record and every reference to it (and parameter names)
p.remove("s03_diff150")                        # refused while another record references it, naming the referrers
p.validate()                                   # list[Diagnostic], every problem (§6.2)
p.describe("s09_ratio")                        # the tree as text and as a dict (§6.3)
p.save("T2.azr")
```

Required keyword-only arguments are the fields the table marks required for every row of the kind; quantity-specific
required fields are checked when the record is validated, at the call. The signature and the docstring of each
`add_*` are generated at import from `_azure2.schema()` (`inspect.Signature`), so `help(p.add_dataset)` lists the
fields, and an IDE completes them; pyazr spells no field name itself.

Structured values serialise to the file text: `norm=pyazr.Nuisance(1.0, "fit", prior=(1.0, "5%"))` ≡
`norm="1 fit prior 1 5%"`; `areal_density=(5e17, "atoms/cm2")` ≡ `"5e17 atoms/cm2 fixed"` only if the mode is given:
a tuple without a mode is refused (the mode is required in the file too).

### 4.2 What is validated and what the user sees

Each call sends the record to `_azure2.validate_record(kind, name, fields, project)`, the engine's validator (§7.3),
and raises `pyazr.FormatError` carrying the list of diagnostics; the message is the engine's text:

```
pyazr.FormatError: dataset 'Y_curve':
  field 'charge_state' is missing: unit 'counts/uC' is per charge, and the number of beam particles per µC
    depends on the charge state q (integer >= 1). Add charge_state=<q>.
  field 'efficiency' is missing: unit 'counts/uC' counts detected photons; give the detection efficiency times the
    branching as a nuisance, e.g. efficiency="0.0182 fit prior 0.0182 3%".
```

References to other records are checked at the call against the project as it stands (add the target before the data
set that uses it); `Project.open` resolves them order-independently.

## 5. The same records in the GUI

The web UI (`ARCHITECTURE.md` §11) gets the table from `/api/schema` (= `pyazr.schema()`) and builds every form
from it. A form field's id is the field key, its label and help text come from the table, its widget from the value
type:

| value type | widget |
|---|---|
| enum (quantity, unit, kind, mode, frame, straggling model, expansion, multipole) | drop-down filled from the table; options that the table refuses for the current reaction are shown greyed with the refusal message as tooltip |
| reaction expression | entrance-pair and exit-pair drop-downs from the project's pair labels (`²⁷Al(p,γ)` style, `qt_parity.md` §4), exit list includes `[total]`; "+ term" for a sum, "÷" for a ratio; scale and fixed angle per term |
| number with unit | number box + unit drop-down (units of the field's dimension) |
| nuisance | value, unit, mode radio (fixed / fit / profiled where allowed), prior centre and width (absolute or %), grid max |
| reference | drop-down of records of the right kind, "+ New …" opens that record's form |
| columns | the first rows of the file, one role drop-down above each column, roles limited to those the quantity accepts |
| formula | text box with parameter table and a plot of the formula over the data's energy range |

Fields appear only where the table says they apply: choosing `quantity yield` shows `target` (required), `charge_state`
and `efficiency` only for per-charge and counts units; `phase_shift` shows J, l, s; `legendre` shows `expansion` and
sets the column role `order`; `ratio_to_rutherford` appears only for an elastic, charged, non-identical reaction.

**Preview.** Beside the form, the engine's description of the record (`_azure2.describe`, §6.3), updated on each
edit: the tree in words ("ratio of dσ/dΩ(c.m.) of ¹²C(α,α)¹²C to the Rutherford cross section, lab angles converted
with α + ¹²C elastic kinematics; folded over target TaO (layer, N fixed) and beam B1; normalization fitted, prior 1 ±
3 %"), the units of each column after reading, the number of rows read and excluded with the reason, the parameters
the record creates, and a plot of the data with the current model. The text is produced by the engine, so the GUI,
pyazr's `describe()` and `AZURE2 --describe` say the same.

Tables. The Data tab lists records by name with reaction, quantity, unit, target, norm; references are links. Deleting
a record that is referenced is refused with the list of referrers (no renumbering exists to do).

## 6. The same records for an agent

An AI agent is the third front end. It needs to write a valid project without reading prose, to learn every problem
at once in a form it can act on, and to check what a record means before computing.

### 6.1 Schema

`AZURE2 --schema` prints, and `pyazr.schema()` returns, one JSON document generated from the table (§7):

```json
{"schema_version": "3.0", "engine_version": "2.2.0",
 "sections": {"data": ["dataset", "normalization", "energy_scale"], "targets": ["target", "beam", "detector"],
              "extrapolations": ["extrapolation"], "modelData": ["datum"]},
 "value_types": {"nuisance": {"grammar": "<value> [<unit>] fixed|fit|profiled [prior <centre> <width>[%]] [grid_max <value>]"}, "...": {}},
 "records": {"dataset": {"fields": {
     "reaction": {"type": "reaction_expr", "required": true, "doc": "..."},
     "unit": {"type": "enum", "values_from": "quantity.units", "required": true},
     "charge_state": {"type": "integer", "min": 1, "applies_when": {"unit_dimension": "per_charge"}, "required_when": "applies"},
     "...": {}}}},
 "quantities": {"ratio_to_rutherford": {"title": "...", "units": {"1": 1.0}, "columns": {"energy": ["E_lab", "E_cm"], "angle": ["theta_lab", "theta_cm"]},
     "requires": [], "allows": ["target", "beam", "detector", "quadrature", "fold_windows"], "reaction": "elastic, charged, non-identical",
     "tree": "Ratio(Fold[dsdo_cm(r)], Fold[rutherford(r)])", "exfor": [",DA,,RTH", "PAR,DA,,RTH"],
     "refusals": [{"code": "RTH_IDENTICAL", "message": "..."}]}},
 "units": {"mb": {"dimension": "cross_section", "factor": 1e-3, "engine_unit": "b"}},
 "diagnostic_codes": {"E_FIELD_REQUIRED": "...", "E_UNKNOWN_FIELD": "...", "ASK_CHARGE_STATE": "..."}}
```

Everything a front end needs is in it: record kinds, fields with types, units and allowed values, required and
conditional fields, which fields apply to which quantity, refusals, diagnostic codes. It is versioned
(`schema_version`); a golden copy is committed and a change needs a version bump (check 21).
**[Correction, step 2 (2026-11-10): `engine_version` is left out of the document: the committed copy would change with every engine version
without a change of the table, which check 21 cannot tell apart; the engine version is the binary's own.]**

### 6.2 Validate

`AZURE2 --validate project.azr [--json]` and `Project.validate()` read the project and its data files (rows, columns,
kinematic ranges), compute no observable, and return **every** problem, not the first:

```json
[{"severity": "error", "code": "E_FIELD_REQUIRED", "file": "13C.azr", "line": 41, "column": 1,
  "section": "data", "record_kind": "dataset", "record": "Y_curve", "field": "charge_state",
  "found": null, "allowed": "integer >= 1",
  "message": "unit counts/uC is per charge: the data set needs charge_state (the beam's charge state q)",
  "fix": {"action": "insert_field", "after_line": 40, "text": "    charge_state  ?"}},
 {"severity": "error", "code": "E_UNKNOWN_FIELD", "file": "13C.azr", "line": 44, "record": "Y_curve",
  "field": "taget", "allowed": ["target", "..."], "message": "unknown field 'taget'",
  "fix": {"action": "replace_key", "line": 44, "text": "target"}}]
```

A `fix` changes syntax only (a misspelled key, a missing required field as a `?` placeholder, a unit alias). It never
supplies a number: `?` is itself refused by the reader ("field charge_state is '?': a value is needed"), so an
unanswered question cannot run (decision F18). The exit code is 0 with no error, 2 with errors. The text form prints
the same records one per line, `file:line: record.field: message`.

The converter's questions use the same structure with `"severity": "ask"` (§8.3), and `--answers answers.json`
(`{"Y_curve.charge_state": 1}`) supplies them; the GUI shows them as a dialog, pyazr raises `pyazr.NeedsInput` with
`.questions`.

### 6.3 Describe

`AZURE2 --describe project.azr [--json] [--record NAME]` (`data_model.md` S9, `--print-predictions`) prints for each
record its tree, the units and frames of its columns after reading, the counts of rows read and excluded, the
apparatus applied in order, and the fit parameters it creates:

```
s09_ratio: Ratio( sigma(15N(p,a)12C[4.439]), sigma(15N(p,a)12C) ) * J(den)/J(num)
   value 1 (dimensionless), E_lab MeV; 20 rows, 0 excluded; no apparatus; norm fixed 1; parameters: none
```

The JSON gives the tree as nested `{"node": "Ratio", "children": [...]}` objects with leaf selectors (`kind`, `entrance`,
`exit`, `filter`, `options`), so an agent can compare what it meant with what the engine will compute.

### 6.4 Format properties that make an agent's edits safe

- One record per block, one unique name, greppable: `grep -n '^dataset\|s02_diff090' T2.azr` finds the record and
  every reference to it.
- No meaning by position or order: records in any order; fields in any order; a field given twice is an error.
- No cross-references by index: targets, norms, scales, THM experiments, priors and parameter names refer to names.
- Local edits: adding a data set adds a block; it renumbers nothing and changes no other record's meaning, parameter
  name or prior. Removing a referenced record is refused, naming the referrers.
- Parameter names derive from record names (`<record>.<field>[.<sub>]`), so a prior file or `<parameterSettings>` row
  written today still means the same parameter after any edit that does not rename that record; `rename` renames
  them everywhere in the project and in the project's prior-covariance file.

### 6.5 Writing a project from a logbook tool or a spreadsheet

A `target` and a `dataset` are self-describing records, so a tool outside AZURE2 writes them without knowing the
rest of the project: `pyazr.records_from_csv(path, kind)` reads a sheet whose header row holds field keys (one row per
record, `name` column required) and returns records for `Project.add`; the same validator runs. Data files may carry
their column roles in a header line (`columns from_file`). The logbook the maintainer considers would emit exactly
these records; it needs no AZURE2 code beyond the schema.

## 7. The table

### 7.1 Columns of the quantity table

| column | content |
|---|---|
| `name` | the field value (`ratio_to_rutherford`); lower case, underscores |
| `aliases` | accepted alternatives, normalised by the reader and written back as `name` (`dsdo`, `rth`) |
| `title`, `doc` | GUI label and help; the docs table row |
| `record_kinds` | `dataset`, `extrapolation`, `datum` |
| `units` | allowed units of the value column, each with its factor to the engine unit and its dimension |
| `columns` | allowed roles per axis (energy, angle, order, window) and which are required |
| `reaction` | constraint on the reaction expression (elastic only, charged, non-identical, neutral entrance, γ exit, composites allowed, ratio only) |
| `requires`, `allows` | fields required by and allowed for this quantity (others are refused, not ignored) |
| `conditional` | fields required when another field has a value (per-charge unit ⇒ `charge_state`) |
| `norm_modes` | allowed modes; `prior` allowed or not (arbitrary units: `fit` without prior only) |
| `tree` | the builder in `experiment/` and its template, printed by `--describe` |
| `exfor` | the SF5–SF8 codes (with modifiers) that map to the row, used by the EXFOR importer and its tests |
| `refusals` | conditions with a code and a message |

### 7.2 Rows

Engine units: MeV, deg, b, b/sr, MeV·b, atoms/cm², MeV cm²/atom, K, eV (widths, U_e), fm. Unit factors: b 1, mb 10⁻³,
ub 10⁻⁶, nb 10⁻⁹ (also /sr); MeV·b 1, keV·b 10⁻³, eV·b 10⁻⁶; eV 10⁻⁶, keV 10⁻³ MeV; atoms/b = 10²⁴ atoms/cm²; per-charge
units 1/uC, 1/mC, 1/C, with the model scaled by 10⁻⁶ C/(q e) for 1/uC etc., e = 1.602176634 × 10⁻¹⁹ C from
`physics/include/azure2/physics/Constants.h`. Data sets:

| name | units | columns (energy; angle) | requires / allows | tree | EXFOR | refusals |
|---|---|---|---|---|---|---|
| `sigma` | b family; `1` with a `/` reaction; `arbitrary` | E_lab, E_cm, E_decay; none | allows `multipole E1|E2`, `secondary_gamma`, `target` (thin layer: target-averaged σ, `data_kinds.md` eq. 2), `beam`, `detector`, `thm` | σ(r) [filter], Sum/Scale/Ratio from the reaction, Weighted if the target has isotopes, folds inside the leaf (M4) | `,SIG`, `PAR,SIG`, `,SIG,,TTA`; `/`, `+` combinations | charged-particle elastic σ ("infinite; give dsigma_domega or ratio_to_rutherford", LEXFOR p. S.2); `multipole` with a particle exit; target `kind infinite` ("use quantity yield") |
| `sigma_nuclear` [step 3a, D1] | b family; `1` with a `/` reaction | E_lab, E_cm; none | every term elastic with a charged entrance; sums and ratios as for `sigma` (a natural target's isotopes, `tests/natural_target`); allows `target`, `beam`, `detector` | the angle-integrated elastic cross section of the nuclear amplitude alone, π/k² Σ g \|T\|² (`GenMatrixFunc`); for a charged pair it leaves out the Coulomb amplitude and its interference and is not a measurable cross section: for tests, diagnostics and comparison with phase-shift formulas | — | inelastic; neutral entrance (there `sigma` is the same integral and the one name); in a `dataset` the warning `W_NOT_MEASURABLE` |
| `sigma_total` | b family | E_lab, E_cm; none | reaction elastic, neutral entrance | σ_tot(a→a) by the optical theorem (`GM` 103–119) | `,SIG` with SF3 `TOT` | charged entrance |
| `dsigma_domega` | b/sr family | E_lab, E_cm; theta_lab, theta_cm (value in the angle's frame) | allows `secondary_gamma`, target, beam, detector, thm | dσ/dΩ_cm(r, θ), lab data converted at load | `,DA`, `PAR,DA` | `[total]` exit; `,DA,,REL` without `norm fit` |
| `sfactor` | MeV·b, keV·b, eV·b | E_cm, E_lab; none | allows `multipole` | Scale_SF(E)∘σ, S = Eσe^{2πη} with the engine's η and atomic masses (`data_kinds.md` §4.1) | `,SIG,,SFC`, `PAR,SIG,,SFC`, `EP/PAR,SIG,,SFC` | neutral entrance |
| `ratio_to_rutherford` | 1 | E_lab, E_cm; theta_lab, theta_cm | reaction elastic, charged, not identical | Ratio(Fold[dσ/dΩ], Fold[σ_R]); new leaf `rutherford` (S14) | `,DA,,RTH`, `PAR,DA,,RTH` (elastic) | identical particles ("use ratio_to_mott", K2); neutral; inelastic ("not defined: give dsigma_domega") |
| `ratio_to_mott` | 1 | as above | reaction elastic, identical particles | Ratio(Fold[dσ/dΩ], Fold[σ_Mott]) | `,DA,,MOT` (identical-particle reading only, `data_kinds.md` §3.2 finding 2) | non-identical |
| `legendre` | 1 for normalised expansions, b/sr for `absolute`, b for `4pi` | E_lab, E_cm; `order` | requires `expansion a_over_a0|absolute|rsl|2l2|4pi`; norm `fixed` for normalised expansions (today's rule, `data_formats.rst:143-154`) | linear maps of the c.m. coefficient leaves A_L | `,DA,,LEG`, `LEG/RS`, `LEG/RSL`, `LEG/2L2`, `LEG/4PI` with `NUMBER-CM` | lab-frame coefficients; `RSD`, `RS0`, `1K2`, `LEG/RTH` (K3) |
| `analyzing_power` | 1 | E_lab, E_cm; theta_lab, theta_cm | spin-½ projectile | A_y(r, θ), folded as A_yσ over σ | `,POL/DA,,ANA` | other projectile spins; tensor analysing powers (`data_kinds.md` §5.2 item 7) |
| `phase_shift` | deg, rad | E_lab, E_cm; none | requires `J`, `l`, `s`; elastic | δ_{J l s} (`data_formats.rst:122-141`) | — | non-elastic; (J, l, s) not allowed by the pair |
| `yield` | 1/particle, 1/particle/sr (with an angle column), 1/uC, 1/mC, 1/C, counts/uC, arbitrary | E_lab; none or theta_lab | requires `target` (`layer` or `infinite`); per-charge units require `charge_state`; counts units require `efficiency`; `arbitrary` requires `norm … fit` without prior | Fold_target∘σ (or dσ/dΩ) × N·10⁻²⁴ for a layer, ∫σ/ε for infinite (`data_kinds.md` eqs. 3, 6), × 10⁻⁶ C/(q e) for per-µC, × ε for counts | `,MLT,,TT`, `,PY,,TT`, `TT/CH`, `PAR,MLT,G,TT/CH` | target `kind none` or `sample`; `TTY` activation yields (K9); per-charge without `charge_state` |
| `transmission` | 1 | E_lab; none | target `kind sample`; neutral entrance | Transform_T(n; Sum) (`data_model.md` §3) | `,TRN` [I] | charged entrance |
| `neutron_yield` | 1 | E_lab; none | target `kind sample`; neutral entrance | Transform_Y(n; …) (SAMMY self-shielded yield; today's `yield n`) | — | charged entrance |
| `beta_delayed_spectrum` | arbitrary | E_decay; none | β entrance pair; norm `fit` | the leaves with EnergyIn = DecayLab | — | non-β entrance |

Model data (`datum`):

| name | units | requires | value | EXFOR | refusals |
|---|---|---|---|---|---|
| `resonance_strength` | eV, meV, ueV, keV | `level`, `J`, `reaction` | ωγ, eq. (8) of `data_kinds.md` §4.7 | `,WID/STR`, `PAR,WID/STR` | E_r < 0 (use `anc`); `WID/STR,,RG` without the spins (importer) |
| `resonance_energy` | MeV, keV, eV | `level`, `J`, `frame E_x|E_cm|E_lab`, `pair` for E_cm and E_lab | observed E_λ in the stated frame | `,EN` | — |
| `total_width`, `partial_width` | eV, keV, MeV | `level`, `J`; `pair` for partial | observed Γ at E_r | `,WID`, `PAR,WID` | `,WID,,G`, `,WID,,2G` unless the importer removes g |
| `branching_ratio` | 1, % | `level`, `J`, `reaction` | Γ_b/Γ | — | — |
| `anc` | fm^-1/2 | `level`, `J`, `pair`, `l`, `s` | C from γ̃ through the Whittaker function at a_c | — | open channel |
| `scattering_length` | fm | `pair`, `s` (or `coherent`); optional `imag`, `imag_error` | Re a₀, Im a₀ (`include/ScatteringLength.h`) | `,AMP`, `COH,AMP` | charged pair |

Other closed lists in the same table: record kinds and their fields; column roles; frames (`lab`, `cm`, `E_x`);
norm modes; target kinds (`layer`, `infinite`, `sample`, `none`); stopping models (`formula`, `evaluation`);
straggling models (`none`, `bohr`, `coefficient`); beam-profile shapes (`skewed_gaussian`); multipoles (`E1`, `E2`,
later `M1` …, M8); Legendre expansions; grid kinds (`range`, `list`, `file`, `adaptive`); quadrature defaults
(width multiplier 20, points per width 50, the `AdaptiveIntegrationGrid::GridConfig` defaults,
`TargetEffect.cpp:24-28`, logged when used; `points` has no default and is required on a data set with a fold, as every
`<targetInt>` line states it today, `TargetEffect.cpp:61`).

### 7.3 Where it lives and how it is read

- `io/` (no engine headers, `project.md` §2) gets one generic `RecordReader`: it splits a section into records and
  fields with their line and column, handles comments and quotes, and knows no field. It returns plain
  `Record{kind, name, line, fields[{key, tokens, line, col}]}`.
- `experiment/` (`data_model.md` M1, M13) holds the table: `experiment/src/Schema.cpp` (record kinds, fields, value
  types), `Quantities.cpp` (the rows of §7.2), `Units.cpp`. `Schema::Validate(records)` is the one table-driven routine:
  for each field it looks up the key, dispatches on the field's value type (a closed set of about twelve: number,
  integer, number-with-unit, nuisance, enum, reference, reaction expression, path, formula with parameters, list of
  numbers, column roles, grid), and collects diagnostics; `Schema::Build` then makes the `Measurement`, the
  `ModelDatum` and the apparatus specs, the tree from the quantity row's builder.
- A new field is a row; a new quantity is a row plus a tree builder plus a test; a new value type (rare) is code in
  the dispatcher plus a test.
- Element symbols for reaction labels: one table in `physics/` (Z → symbol), checked against `pyazr/data/mass_1.mas20`.
- Binding (`api/src/pybind_schema.cpp`): `_azure2.schema()` (the JSON of §6.1 as a dict), `validate_text`,
  `validate_record`, `describe`, `convert`. The CLI's `--schema`, `--validate`, `--describe`, `--convert` call the same
  functions.
- Docs: `tools/gen_schema_docs.py` writes `docs/source/reference/format3.rst` (the tables of this section) from the
  schema; the generated file is committed.

### 7.4 Consistency checks (continuing `data_model.md` §7, checks 13–16)

- **17.** Field keys, quantity names and unit tokens are spelled in `experiment/src/{Schema,Quantities,Units}.cpp`
  only: a grep of `src/`, `api/`, `pyazr/`, the web UI and `io/` for any of the table's strings (taken from the
  built schema) fails, except generated files and tests.
- **18.** `docs/source/reference/format3.rst` equals the output of `gen_schema_docs.py`.
- **19.** Every quantity row has a tree builder, a unit test named after it in `experiment/tests/`, at least one
  test-project data set, and, if it lists EXFOR codes, an importer test on one real entry.
- **20.** Every `.azr` under `tests/` and `examples/` is format 3 and validates with zero errors; `--convert` of the
  committed old fixtures reproduces the committed converted files byte for byte; converting a format-3 file is
  refused with "already format 3".
- **21.** The schema JSON equals the committed `tests/schema/schema_v3.json` unless `schema_version` changed.
- **22.** pyazr's `_KEYWORDS`, `_DATA_CODE`, `_OBSERVABLE`, `_EXTRAP_OBSERVABLE` are gone (today
  `pyazr/datasets.py:77-122`, `pyazr/azrfile.py:2563-2568`); check 6 (engine and pyazr keyword sets equal) is
  replaced by 17.
- **23.** No code recognises a parameter's kind by its name (prefix matches on `segment_`, `density_`, `screening_`,
  `doppler_`, `abundance_`, §1): the parameter vector carries a kind attribute (`ARCHITECTURE.md` §3, `inference/`
  `ParameterVector`).
- **24.** Every diagnostic emitted by the validator has a code listed in the schema, and every code has a test that
  triggers it.

## 8. The converter

`AZURE2 --convert old.azr new.azr [--json-report r.json] [--answers a.json]`, `pyazr.convert_azr(old, new)`. It runs
`io::ConvertToFormat2` first (`io/src/Convert.cpp:40-120`), then the format-3 step, then parses its own output with the
format-3 validator and refuses to write a file that does not validate (except `?` placeholders of answered-later
questions, written only with `--allow-questions`). It lives in a target `azure2_convert` linking `io`, `experiment` and
`physics` (it needs pair labels and the table). It reads the data files (to decide windows and column counts) and
never writes them. It never overwrites its input or a file it did not create; side files get a `.f3` suffix.

### 8.1 Construct by construct

| old | new | note |
|---|---|---|
| data line, active flag 0 | `active no` | the record keeps its place in output numbering (§8.5) |
| entrance, exit keys | `reaction` label | exit −1 → `[total]`; ambiguous label (two pairs with the same partition and state) → ask |
| Emin Emax | `energy_window … MeV`, tokens verbatim | |
| θmin θmax | `angle_window … deg` for quantities with an angle; dropped for angle-integrated data [settled, step 3a: `ESegment::IsInSegment` tests the angle window of differential data (isDiff 1, 4, 7) only, so the old engine never filtered angle-integrated rows by angle; no question] | |
| data isDiff 0, 1, 2, 3, 4, 5, 6, 7, 9 | `sigma`; `dsigma_domega` + `theta_lab`; `phase_shift` + J, l (+ s); `sigma` + `[total]`; `dsigma_domega` + `theta_cm`; `sigma` + `multipole E1`; `… E2`; `analyzing_power` + `theta_cm`; `legendre` + `order` + `expansion a_over_a0` (`absolute` with `absolute 1`) | survey: 170, 257, 7, 17, 65, 0, 0, 99, 10 lines |
| test isDiff 0, 1, 2, 3 + K, 4, 5, 7 | `sigma`; `dsigma_domega` lab; `phase_shift`; `legendre expansion a_over_a0 orders 0-K` [settled, step 3a: `GenMatrixFunc` divides each coefficient by the angle-integrated cross section (`angularCoeff[k] … 4 / angleIntegratedXS`) and `EData` writes a_0 … a_K, i.e. a_k = A_k/A_0, not the absolute A_k]; `sigma [total]`; `dsigma_domega` c.m.; `analyzing_power` | survey: 112, 74, 17, 1, 3, 16, 2 |
| isDiff + 10 | `thm <experiment>` from `<thm>` `segments=`, or `thm global` | 19 data, 6 test lines |
| isDiff 8, 18 and other refused codes | refused, as the old engine refused them (`EData.cpp:397`) | |
| `spin s` | `s`; absent and unique from the pair spins → written; absent and not unique → the old engine refused the line | |
| norm, varyNorm, normError | `norm <n> fixed|fit [prior <c> <e>%]`, c = the `prior_centre` row if any, else n | vary 0 with an error: `fixed prior …` kept (a fixed parameter pays no prior, D5) |
| shift, shiftError, vary (3, 1 or 0 tokens) | `shift offset <s> MeV fixed|fit [prior <c> <e>]`; omitted when 0 and fixed | |
| file | `file` (quoted if it has blanks) | |
| composite `1 0 …` / `1 1 …`, all three layouts | `reaction` sum `+` / ratio `/`; scaling −999 → none, else `c*`; angle −999 → none, else `@ θ deg <frame of the line>` | survey: 29 data, 54 test composites; the short-list warning case (`ComponentList.h:71-75`) → plain reaction and a note |
| UPOS `1 L Ic δ` | `secondary_gamma L <L> Ic <Ic> delta <δ>` | |
| `sqrtshift b e v` | `shift sqrt <b> MeV^1/2 …` | 14 |
| `shiftfn` on one segment | `energy_scale <ds>_scale` + `shift <ds>_scale` | 13 |
| `shiftgroup L` | `energy_scale L` from the first member (offset, sqrt, function); members `shift L`; the other members' shift values → note "not used by the old engine; dropped" | 20 |
| `shiftfngroup L` | `energy_scale L` with the function; members keep `shift offset` | 2 |
| `normgroup L` | `norm L` record from the first member; members `norm L`; others' values → note | 6 |
| `covariance f` | `covariance f` | 34 |
| `sample s` / `natural 1` | `target` of kind `none` (or the data set's existing target) with `isotope` lines from `<abundances>` | 19 + 5 |
| `total 1` | `sigma_total` | 2 |
| `transmission n [v σ]`, `yield n [v σ]`, `thicknessgroup` | `transmission` / `neutron_yield` + target `kind sample`, `areal_density n atoms/b …`; a group → one target | |
| `<segmentsTest>` line k | `extrapolation x<k>`, grid `range lo hi MeV lab step s` (§8.5) | 232 lines |
| `<targetInt>` line | one `target` (if it integrates, screens, or Dopplers), one `beam` (if it convolves or has a profile), one `detector` (Q, UDR, TPC σ); `quadrature`, `fold_windows` on each listed data set **and on the extrapolation with the same index**, which the old engine also folded (`EData.cpp:958`) — written explicitly, with a note | 140 lines, 35 inactive (dropped with a note), 64 multi-segment |
| a listed index that is not a segment | dropped with a note (the old engine skipped it, `EData.cpp:1502`) | |
| density 0 with integration | `kind infinite` | 21 |
| `densityfit v σ [max]`, `densitygroup L` | `areal_density … fit prior … grid_max …`; a group → one target with the leader's density; a member's own density → note (ignored by the old engine; `tests/target_density/fit_group.azr` has 1e18 and an unused 2e18) | members that differ in stopping or straggling → ask |
| straggling flag, coefficient (0 = Bohr) | `straggling none|bohr|coefficient c keV^1/2 cm` | |
| constant σ; expression with `sigmalab`; expression without | `beam spread σ MeV lab`; `spread formula … lab`; `spread formula … cm` (the pre-2026-10 convention, `TargetEffect.h:46-49`) | every Gaussian becomes a beam spread: the kernel is the same, the converter cannot tell a detector resolution |
| `beamprofile n {ξ ω α w} tpc trunc db` | `beam profile …` ×n, `truncation`, `detailed_balance`; `detector energy_resolution tpc MeV lab`; columns 5–6 → `window_lo window_hi` | 4 |
| Q coefficients | `detector attenuation` | |
| `udr "f" L burst channel centred` | `detector udr file … flight_path L m burst … ns channel … ns centred yes|no` | |
| `screening Ue v σ` (`1 "6,7" 200 1 0.002 … screening 1000 1 200`, `tests/electron_screening`) | `target … screening 1000 eV fit prior 1000 200` | 13 |
| `doppler T [v σ Tmax]` | `target … doppler free_gas T K …` | 2 |
| ranges, auto windows, blend, tolerance | `fold_windows` | |
| `<abundances>` `Z A w v σ` per compound nucleus | `isotope <target nucleus> w …`, the target nucleus = the heavy particle of the entrance pairs whose compound nucleus is (Z, A) among the data sets using the sample; several different target nuclei for one compound nucleus → ask; `sum σ` and `abundanceSumSigma=` → `isotope_sum prior 1 σ` | `natural_neon`: 21Ne → 20Ne 0.9048, 23Ne → 22Ne 0.0925, 22Ne → 21Ne 0.0027 |
| `<scatteringLengths>` `a pair s Re σ [Im σ]` | `datum` `scattering_length` | 4 files |
| `<thm>` `segments=`, `weight[k]` | data-set field `thm <exp>`; `thm_weight` on the extrapolation | 10 files |
| `<parameterSettings>` rows `segment_<k>_*`, `density_<n>`, `screening_<n>`, `doppler_<n>`, `abundance_*` | renamed; `prior_centre` rows folded into the fields and removed | level rows unchanged |
| `<mcmc>` `<parameters>` priors | renamed | 31 files |
| `priorCovariance=f` | `f.f3` with renamed names; key updated | 1 file |

`natural_neon` segment 1 (`tests/natural_neon/natural_neon.azr`, line `1 2 2 0 2.33 … data/Vaughn_ntot.dat 1 0 6 2 1
-999 2 3 -999 4 4 -999 4 5 -999 6 6 -999 6 7 -999 0 natural 1`) becomes

```
dataset Vaughn_ntot
    reaction   20Ne(n,n)20Ne + 20Ne(n,a)17O + 20Ne(n,n)20Ne[1.634] + 22Ne(n,n)22Ne + 22Ne(n,n)22Ne[1.275] + 21Ne(n,n)21Ne + 21Ne(n,n)21Ne[0.351]
    quantity   sigma
    unit       b
    target     natNe
    ...
```

and `tests/13N` segment 1 (`artemov.dat`, a thick-target line `"1,2"` with density 0) becomes `quantity yield unit
1/particle target tgt_artemov` (kind `infinite`, stopping formula of eight parameters) with the CONVERTED note: the
file's numbers (4 × 10⁻¹⁴ at 0.19 MeV) were compared as reactions per incident proton, which is what the record now
says.

**[Corrections, step 3a (2026-11-10), from the converter (`convert/`) and the old engine's code:**
- **`shiftgroup L`**: the group shares the constant shift and the leader's `shiftfn` only (`EData::ResolveShiftGroups`,
  `ResolveShiftFunctions`); each member keeps its own `sqrtshift` (a sqrt term and a function on one member are refused
  by the old engine). `energy_scale L` holds the leader's offset and function; the members write `shift L` and their
  own `shift sqrt …`.
- **Shift priors** are measured from zero (`ESegment::AddShiftPenalty`): `shift offset <s> MeV fit prior <c> <σ>` has
  c = 0 unless a `prior_centre` row gives one; `shift sqrt` and the `param`s of a function `prior 0 <σ>`. A shift of 0
  that is not fitted is omitted. A norm keeps its prior on `fixed` (as in the table above); a shift only when fitted.
- **THM norms**: every THM data set with a fitted norm has it profiled (`ESegment::IsProfiledNorm`), `thm global` too
  (§3.9 and the validator now allow it); its unit is `arbitrary`. A THM test line is `thm global` (experiments list
  data lines only).
- **Records shared by several data sets** (F7) are named after their `<targetInt>` line (`beam_targetInt2`), after a
  THM experiment when all users are its data sets (`beam_LC09`), after the data set (`tgt_s01_int_thin`) only when they
  belong to one; identical beams and detectors (and targets without a free parameter) are one record. Data sets are
  named from their file names, extrapolations `x<k>`, model data `scattering_length_<n>`; a name taken already gets
  `_2`, `_3`.
- **Residual states** are written as `<levels>` writes the pair's excitation (`[4.43891]`); §3.11's T2 and the fixture
  `experiment/tests/data/T2_format3.azr` follow (the comparison is the ctest `convert_fixtures`).
- **A phase shift without `spin s`** takes the one channel spin of its (J, l) among the entrance pair's active
  `<levels>` channels (and the hard-sphere channels with `includeAllJPi=1`), as `EData::Initialize` requires; several:
  refused, as the old engine refused it.
- **Test grids** of one point (Emin = Emax, or a step of 0: the old loop's `break`) are written as `list`.
- **Numbers** are written as the old line wrote them when a token of the line reads back to the same double (the
  longest such token), else as the shortest decimal that reads back to it.
- **Side files** are `<stem>.f3<extension>` (`prior_cov.f3.txt`, as in §3.11).
- **Also taken by the table** (they are read by the old engine and were refused, `experiment/` commit of step 3a):
  legendre data folded with a target, beam or detector; beam-profile window columns on `dsigma_domega`; a fixed norm
  on `beta_delayed_spectrum`; `scattering_length` of a charged pair (Coulomb-modified).]**

### 8.2 What the converter writes without asking, and the note it leaves

The converter preserves the meaning the old engine gave a file; it does not reinterpret it. Where that meaning was
implicit and is now explicit, it writes it with a `# CONVERTED:` comment on the field and a line in the report:

- the unit of a thick-target (density 0) data set: `yield` in `1/particle` (`data_kinds.md` §2.3);
- the unit of a thin-target data set: `sigma` in `b`, target-averaged;
- every Gaussian convolution as a beam spread;
- values the old engine ignored (group members' norms and shifts, a density-group member's density, inactive
  `<targetInt>` lines, indices that name no segment): dropped and listed;
- the extrapolations that took a data line's effect by index: the effect written on them.

Old files cannot hold per-charge yields (no field for q exists, `data_kinds.md` §2.3), so a charge state is never a
conversion question; it is one for the EXFOR importer (§6, story d).

### 8.3 What it asks

Asked (severity `ask`, §6.2), with the answer key in brackets:

- a reaction label that does not identify one pair (two pairs with the same particles and state) [`<ds>.reaction`];
- an angle window that excludes rows of an angle-integrated data set [`<ds>.angle_window`];
- a density group whose members differ in stopping, straggling or Doppler [`<group>.split`];
- a sample whose compound nucleus is reached from two different target nuclei in the data sets that use it
  [`<sample>.isotope.<Z>_<A>`];
- a data set with both a `transmission`/`yield` sample thickness and a target-integration line (two areal densities)
  [`<ds>.target`];
- a data set with a `sample` whose `<targetInt>` target is shared with data sets without the sample [`<ds>.target`].

The survey suggested that none of the 189 test files raises a question. **[Step 3a, measured with the converter
(`tools/convert_tree.sh`): 18 files raise questions, of three kinds the note did not foresee — (i) an angle-integrated
cross section of the elastic scattering of a charged pair (16 files: `tests/natural_target`, `tests/hard_sphere`,
`tests/coulomb_nuclear_sign`, `tests/identical_entrance_reaction`), which the old engine computed as the integral of
the nuclear amplitude alone, π/k² Σ g |T|², and which §7.2 refuses as `sigma` without offering another quantity;
(ii) two pairs with the same particles and states (`examples/c12c12_tumino2018`: pairs 1 and 6, ¹²C + ¹²C, differing
only in the THM columns of `<levels>`), which no reaction label can tell apart; (iii) an entrance pair with an excited
target (`tests/thm_coherent_background/toy`: n + ¹⁶O(2.143)), which the grammar of §3.5 cannot write. Each is a
decision for the format (a quantity, pair names, a target state), not for the converter.]** **[Decided (director, 2026-11-10) and implemented in step 3a:
(i) quantity `sigma_nuclear` (§7.2), which the converter writes with a `# CONVERTED:` comment; (ii) `entrance_pair` /
`exit_pair` (§3.5); (iii) the target state (§3.5). The tree now converts without a question (§8.5).]**

### 8.4 The engine refuses older files

A project without `format=3`, or with any of `<segmentsData>`, `<segmentsTest>`, `<targetInt>`, `<abundances>`,
`<scatteringLengths>`, is refused before anything is read:

```
ERROR: T2.azr is a format-2 project (<segmentsData> at line 53). This AZURE2 reads format 3 only.
       Convert it:  AZURE2 --convert T2.azr T2_f3.azr     (or pyazr.convert_azr("T2.azr", "T2_f3.azr"))
       The converter keeps T2.azr unchanged and reports every change it makes.
```

Exit code 3, distinct from a parse error. The same for `param.sav`-style files given with `--param`: the converter
renames their segment and effect rows by the same map (`--convert-params old.sav new.sav --project T2_f3.azr`).

Consequence for the code (check 23): parameter names are user names, so every place that tests a name prefix
(`src/AZUREParams.cpp:9-15, 98, 123-125`, `src/AZURECalc.cpp:1431`, `include/ParameterLabel.h:27-48`,
`src/ParameterLimitsManager.cpp:245, 304`, the labels of `src/AZURE2.cpp:1171-1212`) must read the parameter's kind
from the parameter vector instead. This is a step before the reader switch.

### 8.5 Converting the test suite and the gate

1. Convert every `.azr` under `tests/` and `examples/` (189 files; the 179 run by `tools/runall.sh` are those with a
   project directory) with the converter at one thread into a copy of the tree; the report of each conversion is kept.
2. `AZURE2 --validate --json` on every converted file: zero errors.
3. Run the engine before format 3 on the originals and the format-3 engine on the converted files, menu 1, the three
   bases (`tools/runall.sh` with `JOBS` within the thread budget), and compare with `tools/cls.sh`
   (`classify.py`; `ARCHITECTURE.md` section 10). Required: every output of every data set bit-identical, `chiSquared.out` and `normalizations.out`
   included. This holds only if (a) the reader keeps number tokens verbatim and parses them with the old function
   (`strtod`), (b) prior widths are computed by the old expressions in the old order (norm: n₀σ/100 inside the penalty,
   `src/AZURECalc.cpp:901`; density: σ/100·N₀ at read, `TargetEffect.cpp:355`), (c) data sets keep their file order,
   which is their output number (`src/EData.cpp:2542`), and (d) every `<targetInt>` copy per segment and its per-line
   parameters map to the same kernels and parameters.
4. Extrapolation outputs: identical except where the grid rule changes (F10). Measured on the 232 test lines of the
   survey: with E_k = lo + k·s, 53 grids gain or lose their last point and 135 more move nodes by at most 3.9 × 10⁻¹²
   relative. These outputs are compared at matching energies and classified "grid rule", with the list of the 53.
5. `tests/run_tests.sh` 133/133 and `ctest -LE slow` on the converted tree, then pyazr's round-trip tests (an
   unedited converted file writes back byte for byte).
**[Step 3a, director's D4 (2026-11-10): the converter converts what is written and does not copy the old engine's
semantic refusals; they are the format-3 reader's build checks. `tools/convert_tree.sh` runs the old engine (menu 1,
one thread, as `tools/run_one.sh`) on every original and lists the files it refuses that the converter converts. On
the tree of step 3a (186 of 189 files convert and validate; the 3 refused are refused by the old engine too), the
format-3 engine of 3b must refuse each of these five with the reason given:
- `tests/angular_coefficients/models/bad_elastic.azr` — "segment #1 (angular coefficients): an elastic distribution
  includes the Coulomb amplitude, which the Legendre expansion leaves out";
- `tests/ec_final_state_listed/models/missing_group.azr` — "external capture into pair 3 needs its final state,
  J^pi = 0.5+ at E_x = 0.4953 MeV, as a level of the project …; that J^pi group does not exist";
- `tests/ec_final_state_listed/models/missing_level.azr` — the same final state, "none is listed within 1 keV in that
  J^pi group";
- `tests/natural_target/models/bad_level.azr` — "level ID 1 (J = 1+, E = 8.2506 MeV) has channels of two compound
  nuclei, Z=7 A=13 and Z=7 A=14";
- `tests/relativistic/models/li7pn_edge_rel.azr` — "segment #1: lab angle 44.9 deg at E_lab = 1.9 MeV: beyond the
  kinematic maximum 44.5984 deg (relativistic kinematics …)".
Menu 1 only: a refusal that only an extrapolation run (menu 3) or a data file written by a test's `check.sh` reaches is
not in this list (the `check.sh` of `tests/natural_target`, `tests/angular_coefficients`, `tests/electron_screening`
and the others with `bad_*` files say which runs they expect to fail).]**

**The gate as a tool** (`tools/format3_gate.sh`, `tools/format3_gate.py`): items 1, 3 and 4 with exact grids, menu 1,

```
WORK=<scratch> REF=<tag of the original tree's run> JOBS=3 tools/format3_gate.sh <engine> <source tree> <tag>
```

converts the source tree's `tests/` with `tools/convert_tree.sh` and `--exact-grids` (without the old engine's pass)
into `$WORK/gate-<tag>/converted`, makes it runnable (`format3_gate.py tree`: each original whose conversion
validates, or fails only with `E_DATA_FILE`, is replaced by its conversion under its own name; the others are listed
in `tree.txt`), runs it in the three bases with `runall.sh` as `tab-<tag>` (one thread per run, `JOBS` runs at once)
and compares every output file with those of the reference tag `REF` (`format3_gate.py compare`; `param.par` after
each conversion's parameter-name map; `run.log` not compared). A file is identical, or its difference consists of
classes each listed for that file and project in the table `CLASSES` of `format3_gate.py` (class, file, projects,
why), or it differs. The counts per basis are printed, the list in `$WORK/gate-<tag>/gate.txt`; the exit status is 1
when a file differs or is missing, so a new difference fails the gate until a row with its reason is added to the
table. On the head of the fourth wave (64b8e7ae): of 510 (483 formal) files per basis 472 (445) identical, 38
classified, 0 differ, 0 missing.

6. The converted files replace the originals in one commit; the originals move to `tests/convert/fixtures/` with their
   committed conversions (check 20), so the converter is tested on every construct for as long as it exists.

## 9. Worked user stories

Each story is given as file, Python, GUI and agent. Pair labels assume a project with the pairs named.

### 9.1 (a) ¹³C(p,γ) thick-target yield curve in counts per µC, charge state 1, shared target, thickness fitted

Physics: a ¹³C layer on a backing, its areal density N fitted; a yield curve across the 1.75 MeV resonance region
measured as counts of the primary γ to the ¹⁴N ground state per µC of H⁺ beam, and a second data set (a narrow-
resonance scan used to profile the same target) sharing N. Per-charge data need q; counts need the efficiency × γ
branching ε.

| | (a) |
|---|---|
| file | `target C13A` `kind layer`, `active 13C`, `composition C 1`, `areal_density 1.2e18 atoms/cm2 fit prior 1.2e18 15%`, `stopping evaluation scale 1 fit prior 1 4%`, `straggling bohr` [correction, step 2: and `phase solid`, required with `stopping evaluation`]; `beam pbeam` `spread 0.3 keV lab`; `dataset Y_curve` `reaction 13C(p,g)14N`, `quantity yield`, `unit counts/uC`, `charge_state 1`, `efficiency 0.0182 fit prior 0.0182 3%`, `file data/yield_curve.dat`, `columns E_lab[keV] value error`, `norm 1 fixed`, `target C13A`, `beam pbeam` [correction, step 2: and `quadrature points <n>`, required on a folded data set, §7.2]; `dataset Y_scan` the same with its file and `efficiency 0.0175 fit prior 0.0175 3%` |
| Python | `p.add_target("C13A", kind="layer", active="13C", composition="C 1", areal_density="1.2e18 atoms/cm2 fit prior 1.2e18 15%", stopping="evaluation scale 1 fit prior 1 4%", straggling="bohr")`; `p.add_beam("pbeam", spread="0.3 keV lab")`; `p.add_dataset("Y_curve", reaction="13C(p,g)14N", quantity="yield", unit="counts/uC", charge_state=1, efficiency="0.0182 fit prior 0.0182 3%", file="data/yield_curve.dat", columns="E_lab[keV] value error", norm="1 fixed", target="C13A", beam="pbeam")`; without `charge_state` the error of §4.2 |
| GUI | Data → "+ Data set": reaction p + ¹³C → γ + ¹⁴N (g.s.); quantity "Yield"; unit drop-down → "counts/µC" makes "Charge state q" and "Efficiency × branching" appear; Target "+ New target…" opens the target form (kind "Layer" shows active nucleus, composition, areal density with Fit and prior, stopping, straggling, phase); the second data set picks C13A from the drop-down. Preview: "yield of ¹³C(p,γ)¹⁴N(g.s.) over target C13A (layer of ¹³C, N fitted, prior 15 %, stopping evaluation with scale fitted), beam spread 0.3 keV, × 10⁻⁶ C/(1·e) × ε (fitted, prior 3 %) → counts/µC; parameters C13A.areal_density, C13A.stopping.scale, Y_curve.efficiency" |
| agent | reads `quantities.yield` in `--schema` (requires target; `counts/uC` ⇒ `charge_state`, `efficiency`); writes the three records; `--validate --json` → `[]`; `--describe --record Y_curve --json` → `Scale(1e-6/(q e)) · Scale(ε) · Fold_layer(C13A)[sigma(13C(p,g)14N)]` |

Note: a yield per particle on a finite layer, N·10⁻²⁴⟨σ⟩, with N fitted in the prefactor as well as in the energy
window, is a Phase 5 kernel; today's thin target divides by N (`data_kinds.md` eq. 2) and today's density fit moves only
the window (`TargetEffect.h:89-101`).

### 9.2 (b) α + ¹²C elastic scattering as ratio to Rutherford, lab frame, five angles, 3 % common normalization

| | (b) |
|---|---|
| file | `dataset aC_RTH` `reaction 12C(a,a)12C`, `quantity ratio_to_rutherford`, `unit 1`, `file data/aC_rth.dat`, `columns E_lab[MeV] theta_lab[deg] value error` (the five angles are rows of one file), `norm 1 fit prior 1 3%` |
| Python | `p.add_dataset("aC_RTH", reaction="12C(a,a)12C", quantity="ratio_to_rutherford", unit="1", file="data/aC_rth.dat", columns="E_lab[MeV] theta_lab[deg] value error", norm="1 fit prior 1 3%")`. With `reaction="4He(a,a)4He"`: `FormatError: dataset 'aC_RTH', field 'quantity': 4He(a,a)4He has identical particles; the Coulomb reference is Mott-symmetrised: use quantity 'ratio_to_mott'` |
| GUI | reaction α + ¹²C → α + ¹²C; quantity list shows "Ratio to Rutherford" (enabled: elastic, charged, non-identical); column roles from the file preview; the form reports "5 angles: 30, 60, 90, 120, 150 deg (lab)". Preview: "ratio of dσ/dΩ(c.m.) of ¹²C(α,α)¹²C to the Rutherford cross section at the same c.m. angle; lab angles converted with α + ¹²C elastic kinematics; no apparatus; normalization fitted, prior 1 ± 3 % (one parameter for all five angles: aC_RTH.norm)" |
| agent | `--describe` returns `Ratio(dsdo_cm(12C(a,a)12C), rutherford(12C(a,a)12C))`, value unit 1, angles lab → c.m. |

If the five angles are five files, five data sets name one `normalization aC` record (§3.1 correction) (`value 1 fit prior 1 3%`): one parameter,
one prior.

### 9.3 (c) ωγ and E_r of one level from the literature

For T2's level 5 (3⁻ at E_x = 13.2652 MeV in `T2.azr`, the start value of the fit; the truth of
`tests/upgrade/T2/truth.json` is 13.2692 MeV; Γ_p = 1.702 keV, Γ_α0 = 12.92 keV, Γ_α1 = 0): ω = 7/4 for p + ¹⁵N (spins ½, ½), so
ωγ = 1.75 × 1.702 × 12.92 / 14.62 = 2.63 keV (computed for this note from `tests/upgrade/T2/README.md` §1).

| | (c) |
|---|---|
| file | `<modelData>`: `datum Er_L5` `quantity resonance_energy`, `level 5`, `J 3-`, `frame E_x`, `value 13.2652 MeV`, `error 0.3 keV`, `source "…"`; `datum wg_L5` `quantity resonance_strength`, `level 5`, `J 3-`, `reaction 15N(p,a)12C`, `value 2.63 keV`, `error 0.15 keV`, `source "…"` |
| Python | `p.add_datum("wg_L5", quantity="resonance_strength", level=5, J="3-", reaction="15N(p,a)12C", value="2.63 keV", error="0.15 keV", source="…")`; with `J="2+"`: `FormatError: datum 'wg_L5', field 'J': level 5 is 3- (E_x = 13.2652 MeV); the datum says 2+` |
| GUI | Model data → "+ Datum": quantity; level drop-down "5: 3⁻ 13.2652 MeV" (fills J, shown read-only); reaction from the level's channels; value and error with unit drop-downs. Preview: "ωγ = ω Γ_p Γ_α0/Γ of level 5 with observed widths at E_r, ω = 7/4; model now 2.63 keV; warning: Γ/E_r = 1.3 % > 1 %, the parameter form (8) and the area of the resonance differ (data_kinds.md §4.7)" |
| agent | `--validate --json` returns the Γ/E_r warning with code `W_BROAD_LEVEL`; `--describe` gives the formula and the parameters it depends on |

### 9.4 (d) The same data from EXFOR

`pyazr import exfor <dataset-id> --project p.azr [--name aC_RTH]` (pyazr `sources/`, `ARCHITECTURE.md` §3) writes the
data file with a provenance header and adds the record; the GUI's EXFOR dialog calls the same function.

Filled by the importer, each with its source recorded in the file header and in `exfor`, `doi`, `source`:
`reaction` (the forward form of a tautology, `data_kinds.md` §3.1, matched to the project's pair labels); `quantity`
(SF5–SF8 against the table's `exfor` column: `,DA,,RTH` → `ratio_to_rutherford`); `unit` (Dict. 25 factor, not
prefix matching); `columns` with frames from the headings (`EN` → `E_lab`, `EN-CM` → `E_cm`, `ANG` → `theta_lab`,
`ANG-CM` → `theta_cm`, `NUMBER-CM` → `order`); `error` as the quadrature sum of `ERR-S` and partials flagged U;
`norm` (`REL`/`ARB-UNITS` → `fit` without prior; else `fit prior 1 <F-flagged partials or ERR-SYS>%`); a covariance
file for P/C-flagged partials (`data_kinds.md` §6).

Asked (`ask` diagnostics; GUI dialog; `pyazr.NeedsInput`), never guessed:
- the residual state when `E-LVL`/`E` matches no pair of the project, or several;
- `RTH` for identical particles: Rutherford or Mott (K2);
- for `TT/CH` yields: `charge_state` (EXFOR does not code q), the `target` (thickness, composition, stopping are not
  in EXFOR), and `efficiency` if the data are counts;
- an entry without `ERR-ANALYS`, or with `ERR-T` only: whether there is a normalization uncertainty;
- a `DEP` data set whose parent is in the project; a `SPSDD` (superseded) entry;
- `MOT` whose text means the relativistic point-nucleus cross section (refused).

Refused with the table's message: every code in `data_kinds.md` §5.4.

## 10. Decisions for the director

**Taken (director, 2026-11-09).** F1–F20 as recommended, with these additions:

- *F7.* The converter writes identical `beam`, `detector` and `target` records once and lets the data sets share the
  name; a record is named after a data set only when it belongs to one.
- *F9.* `efficiency` and `norm` multiply: the reader refuses a data set in which both are free without a prior
  (degenerate), with a message naming the two fields.
- *F10.* The gate is run twice. With `--convert --exact-grids` the converter writes every extrapolation grid as the
  `list` of today's accumulated energies: outputs must be bit-identical to the originals on the previous engine. With
  the default (`range`), the differences must be confined to extrapolation rows and are classified (53 grids with
  another point count, 135 moved by ≤ 3.9 × 10⁻¹² relative).
- *F11.* `experiment/` may use `io/`'s plain record structs (`azure2/io/Record.h`); `io/` includes nothing of
  `experiment/`. The rule "experiment/ uses physics/ only" (`data_model.md` M13) is relaxed by exactly this header.
- *F13.* Level names wait for the `<levels>` note, written when step 7 of the `CNuc` split takes the reader out of
  `CNuc` (`cnuc_split.md`).
- *Order of work.* (1) parameter kinds as an attribute of the parameter vector, no recognition by name prefix
  (bit-identical; prerequisite of user-chosen names, check 23); (2) the generic record reader in `io/`, the table in
  `experiment/`, `--schema`, with unit tests, not wired; (3) the format-3 reader producing the specs the engine
  consumes, the converter holding the four old parsers, `--validate`, `--describe`, the two gate runs; the old
  parsers then leave the engine; (4) pyazr and the web forms on the table (their own copies deleted, checks 17–24);
  (5) new quantities, one row each, in the order of `data_kinds.md` §5.3.


**F1. Record syntax.** (a) Key–value record blocks inside today's sections; (b) TOML for the whole file; (c) keyword
blocks on one line per data set (`data_kinds.md` K1). Recommendation (a): one field per line gives local diffs,
line-precise errors and local edits with the existing `io::Document`, and the parser has no arity problem (§3.1). K1
is withdrawn: its reason (old files keep working) no longer holds.

**F2. Record end.** (a) The next header line or the section end; (b) an explicit `end`. Recommendation (a): headers
and field keys are disjoint sets, so a record's extent is unambiguous; `end` adds a line per record and a way to forget
one.

**F3. Names and parameter names.** (a) One namespace, `[A-Za-z][A-Za-z0-9_-]*`, parameters `<record>.<field>[.<sub>]`;
(b) per-kind namespaces. Recommendation (a): a name then means one thing in a grep, a diagnostic and a prior file.

**F4. Residual states.** (a) By excitation energy in brackets, matched to the digits written, `[total]` for the γ sum;
(b) ordinals `a1`. Recommendation (a): an ordinal counts the states in the project, i.e. a position.

**F5. Frames.** (a) Column roles carrying frame and unit (`E_lab[MeV] theta_cm[deg]`), `columns` required; (b) one
`frame` field for the record (`data_kinds.md` §5.1). Recommendation (a): (b) cannot write isDiff 4 (E lab, θ and value
c.m.), the second commonest differential layout.

**F6. Shared nuisances.** (a) Named `norm` and `energy_scale` records referenced by data sets; (b) group labels with a
first-member owner (today). Recommendation (a): no order dependence and no ignored copies.

**F7. Apparatus ownership.** As in §3.7: target (N, stopping, straggling, isotopes, Doppler, screening), beam (spread,
profile), detector (Q, energy resolution, UDR), data set (quadrature, fold windows). Alternative: one `apparatus`
record per data set as today's line. Recommendation: the split, because a target's parameters are shared by
construction and the GUI form of each object is short.

**F8. Charge state.** (a) On the data set beside the unit (K6); (b) on the beam. Recommendation (a): it matters only
for per-charge units, and a beam record shared by runs of different charge states would be wrong silently.

**F9. Counts.** (a) `counts/uC` with a required `efficiency` nuisance (ε × branching, its own parameter and prior);
(b) refuse counts; the user divides by ε and adds its uncertainty to the norm prior. Recommendation (a): the story
asked for counts, and the efficiency then appears in the file with its uncertainty instead of inside a norm prior.

**F10. Extrapolation grid.** (a) E_k = lo + k·s (exact, the last point kept); (b) today's accumulation. Recommendation
(a): it is what a user means; the gate classifies the 53 + 135 changed grids (§8.5 item 4). Extrapolations do not enter
the χ², so no fit moves.

**F11. Where the table lives.** (a) `experiment/`, with a generic `RecordReader` in `io/`; (b) all in `io/`.
Recommendation (a): the table builds trees and needs physics constants (e for per-charge units); `io/` stays dependent
on nothing (M13, check 9).

**F12. Converter strictness.** (a) Write the old meaning explicitly with `# CONVERTED` notes, ask only where the old
meaning is undefined (§8.3), `--strict` turns notes into questions; (b) ask about every implicit unit. Recommendation
(a): the gate requires unattended conversion of 189 files, and the notes say exactly what was assumed.

**F13. Level references in `<modelData>`.** (a) The level ID of `<levels>` plus a required J^π; (b) match by J^π and
energy; (c) level names, when `<levels>` gets a format-3 note. Recommendation (a) now, (c) later: an energy moves in a
fit; the ID is a written label, and J^π catches a wrong one.

**F14. Prior centres.** Always written in the field (`prior <centre> <width>`); `prior_centre` rows go. Alternative:
centre = value unless stated. Recommendation: always written, since a fit that writes back its value would otherwise
move the prior.

**F15. Side files.** The converter writes `.f3` copies of a prior-covariance file and of `param.sav` files it is given,
never overwriting. Alternative: inline priors in a `<priors>` section. Recommendation: copies now; a section can come
later.

**F16. pyazr validation timing.** (a) At each call, references against the project as it stands; (b) only at
`validate()`/`save()`. Recommendation (a), loud and early; `Project.open` resolves out-of-order files.

**F17. Schema export.** (a) AZURE2's own JSON (§6.1), versioned; (b) a JSON Schema (draft 2020-12) of an equivalent
JSON project. Recommendation (a): the project is not JSON, and the table's content (units with factors, trees,
refusals, EXFOR codes) does not fit JSON Schema without extensions.

**F18. Fixes in diagnostics.** (a) Syntax fixes only, `?` placeholders for values, `?` refused by the reader;
(b) suggested values. Recommendation (a): a suggested number is a default by another name.

**F19. Negative or zero errors in data files.** `tests/13N/data/baumann_1618.dat` has negative errors (10 % of a
negative A_y). (a) Refuse zero; accept negative with a warning and use |δy| (the χ² is unchanged, the gate holds);
(b) refuse both. Recommendation (a).

**F20. Extrapolation records.** (a) A separate kind refusing data-only fields; (b) a data set without a file.
Recommendation (a): today's silent skipping of data-only keywords on test lines (`ExtrapLine.h:104-131`) becomes an
error instead of a mode.

**Deliberately not allowed:** positional references of any kind; a default unit, frame or norm mode; an unknown,
duplicated or inapplicable field (refused, never ignored); a value the engine would ignore; per-charge units without
`charge_state`; a prior on an arbitrary-unit norm; `ratio_to_rutherford` for neutral or identical pairs; lab-frame
Legendre coefficients; ratios of sums and parentheses in reactions (until a tree node needs them); a free `predict`
expression (M9 as decided); target folds on THM data; inline data points in the project (data stay in files; model data
are not points); include files; reading format 1 or 2 in the engine.

## 11. What this note contradicts in `data_model.md` and `data_kinds.md`

- `data_kinds.md` §5.1 and K1 (keyword blocks on the segment line, old lines readable, `<targets>` as a one-line
  table): replaced by records (F1), because compatibility is no longer required and the line form fails every
  criterion of §3.1.
- `data_kinds.md` §5.1, "the frame applies to E, θ and value together": replaced by column roles (F5).
- `data_kinds.md` §5.1, "the stoichiometry gives both ε_eff and the abundance weights": kept as the target's two fields
  `composition` (stopping) and `isotope` (weights) because today's stopping is a user formula per active atom
  (`experimental_effects.rst:263-320`) that knows no stoichiometry; with `stopping evaluation` both come from one record.
- `data_kinds.md` §5.2 items 4, 5: `quantity ratio reference rutherford` and `legendre norm rs|…` become the quantities
  `ratio_to_rutherford`/`ratio_to_mott` and the field `expansion`: `norm` is the normalization nuisance, and one
  drop-down is simpler than two linked ones.
- `data_kinds.md` §5.2 item 11 (thick-target yields only): `yield` on a finite `layer` target is added; it is a Phase 5
  kernel.
- `data_model.md` §2.2 `Measurement::key` "segment number in the file": the record name; the number survives only as
  the output order.
- `data_model.md` §3 text form `predict "<expr>"`: not in format 3 (as M9 was decided); its reaction part becomes the
  reaction expression (§3.5).
- `data_model.md` M12 "keep both isDiff tables, map both in one place (`SegmentKinds`)": both tables live in the
  converter only (the maintainer's override already says so); `yield` is the charged-particle yield and the SAMMY
  yield is `neutron_yield`.
- `ARCHITECTURE.md` §3 "today's segment types are read by `io/` into these trees, so old projects load unchanged":
  superseded by the converter.

## 12. Not settled

- The `<levels>` grammar (pairs and levels as named records) is outside this note. Reaction labels are built from the
  pair columns of `<levels>` (`data_formats.rst`, levels table); β-pair and light-excitation labels are sketched only.
- [Settled, step 3a: `ESegment::Fill` does not drop angle-integrated rows by the angle window; test isDiff 3 is
  a_k = A_k/A_0, k = 0 … K (§8.1).]
- The finite-layer yield per particle with N in the prefactor (story a) and the stopping evaluation are Phase 5 work;
  the format reserves the fields.
- `,TRN` as the EXFOR code for transmission was not checked against Dict. 236.
- [Run, step 3a: 18 test files raise converter questions (§8.3).]
- The adaptive extrapolation grid's seeding and tolerance are a proposal (HANDOFF 2026-11-10); its test (a grid that
  reproduces a narrow resonance's area to the stated tolerance) is to be written with it.
- Where the web UI's own code lives was not checked; check 17 assumes it is under `pyazr/`.
