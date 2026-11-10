# Segments, tail keywords and the other `.azr` blocks

This reference says how a data or test segment line is read, what each tail keyword changes in the computed
observable, and how `<targetInt>` and parameter-file entries are matched to segments. The behaviour of each keyword
is specified in `docs/source/reference/data_formats.rst` and `docs/source/user_guide/segments.rst`; fitting with
norms and shifts is azure2-fit; covariance files are azure2-uncertainty.

## 1. Blocks of a project file

- `<config>`: the formalism flag (A-matrix or R-matrix), `output/` and `checks/` directories, check switches, then
  optional `key=value` lines (`includeAllJPi`, `maxOrbitalMomentum`, `abundanceSumSigma`, `relativistic`, ...).
  A-matrix versus R-matrix is a computational choice only: both give the same results. The A-matrix (level matrix)
  is faster with many channels and few levels, the R-matrix (channel matrix) with many levels and few channels.
- `<levels>`: references/azr-levels.md.
- `<segmentsData>`, `<segmentsTest>`: §2 and §3.
- `<targetInt>`: target and resolution effects (azure2-targets), matched by key (§4).
- `<parameterSettings>`: limits, errors, nuisance priors, category, Minuit index. "Level N Energy (MeV)" and
  "Level N Channel M Width (eV)" are the engine's `energy_N` and `width_N_M` of the same file (active levels J-group
  by J-group, M the channel within the J-group); a block written with the line `# level numbering: engine` uses that
  numbering. Older Qt GUI saves placed limits and priors on other levels (azure2-fit).
- `<mcmc>`: walkers, steps, threads.

Prefer `AzrModel` to hand edits.

## 2. The segment line

```
isActive entranceKey exitKey minE maxE minA maxA isDiff [phaseJ phaseL] dataNorm varyNorm dataNormError
[energyShift energyShiftError varyEnergyShift] dataFile [tail]
tail = isAdvanced [operationType nComp {entrance exit angle [scaling]}...] isUPOS [L Ic delta] [keyword blocks]
```

Parser: `convert/include/azure2/convert/SegLine.h`. Energies E_lab (MeV), angles θ_lab (deg) unless the observable code says c.m.
`isDiff + 10` marks a THM (HOES) segment (azure2-thm).

- **Formatting.** The engine tokenizes on any whitespace, but the Qt GUI's segment editor (released builds) expects
  each field at its own column offset and mis-displays a tab-joined line. Match the padding of an unmodified neighbour and diff the two.
- **Old layouts.** The three shift fields before the file are optional; the reader decides the layout by whole
  tokens (`tests/segline_formats`). A misspelt keyword (`normgruop`) is ignored by engine and pyazr by design (older
  engines must ignore newer keywords): when a keyword seems to do nothing, read the segment summary in the run log.
- **Read the tokens after the data-file name before modelling, plotting or extrapolating any segment.** They change
  what the segment computes. Audit one-liner (prints every segment with a non-trivial tail):

  ```
  awk '/<segmentsData>/{f=1;n=0;next}/<\/segmentsData>/{f=0}f{n++;s=$0;sub(/.*\.dat/,"",s);if(s!~/^ *0 +0 *$/)print n": "$0}' model.azr
  ```

### Composite and UPOS tails

- `isAdvanced = 1` makes a composite. `operationType` 0 = SUM of the listed components (the segment's own exit is
  included), 1 = RATIO (dimensionless; the lab→c.m. cross-section conversion is skipped). `nComp = -1` marks the
  newer format: the next integer is the count and each component carries a fourth `scaling` token. A component
  angle <= −900 means unset.
- `isUPOS` is read whether or not the segment is advanced. 1 = unobserved primary, observed secondary: the
  observable is the secondary γ ray from the state the exit pair is left in (e.g. the 4.44 MeV γ of ¹²C(2⁺) in
  ¹⁵N(p,α₁γ)), not the primary cross section. `L` is the γ multipolarity, `Ic` the spin of the state it feeds,
  `delta` its mixing ratio.
- Reading a tail: `file 0 0` plain; `file 0 1 2 0 0` UPOS, E2 to a 0⁺ state, no mixing; `file 1 0 3 ...` a sum of
  three components.
- **An extrapolation cannot reproduce a UPOS segment.** A test line with the same entrance, exit and angle computes
  the primary cross section, which differed from a UPOS γ yield by factors 0.1–3.5 point to point. To draw a UPOS
  curve, copy the segment line unchanged, point it at a pseudo-data grid file (`E 0 1000 100` rows), run mode 1 with
  the `.sav`, and read column 4 of that segment's `AZUREOut` block; check it against the real segment at the data
  energies (agreement 1e-5).

### Keyword blocks

Each is documented in `docs/source/reference/data_formats.rst`; what to remember when reading a model:

| keyword | effect | pyazr / test |
|---|---|---|
| `sqrtshift b bErr vary` | E' = E + ΔE + b·(E/MeV)^½, b in MeV^½, penalty (b/bErr)²; models an offset in an analysing-magnet field reading (E = kB²) | `Segment.energy_shift_sqrt`; `tests/energy_shift_sqrt` |
| `shiftfn "<f(x)>" N v0 s0 vary0 ...` | adds f(x; a₀..a_{N−1}), x the original E_lab (MeV), to the shift; parameters `segment_<key>_energy_shift_fn<j>` | `Segment.shift_fn`; `tests/shift_fn` |
| `shiftgroup <label>` | all members use the constant shift (value, error, vary, penalty) of the first member in file order | `tests/shift_groups` |
| `shiftfngroup <label>` | all members apply the first member's `shiftfn`, each keeps its own constant shift | `Segment.shift_fn_group`; `tests/shift_fn_groups` |
| `normgroup <label>` | all members use the first member's norm (value, vary, error, one penalty); use for an excitation function and the angular distributions it was integrated from | `Segment.norm_group`; `tests/norm_groups` |
| `covariance <file>` | correlated data uncertainties; χ² = rᵀM⁻¹r over the segments naming the file | `Segment.covariance`; `tests/data_covariance`; azure2-uncertainty |
| `sample 1` / `sample <name>` | natural-abundance target (references/data-definition.md §3) | `Segment.natural_sample`; `tests/natural_target` |
| `total 1` | optical-theorem neutron total (references/data-definition.md §2) | `docs/source/user_guide/segments.rst` |

Start-up errors: a segment in both a `shiftgroup` and a `shiftfngroup`; a first `shiftfngroup` member without a
`shiftfn`; a later member declaring one; `shiftfn` together with `sqrtshift` on a group member. A `.sav` value never
frees or fixes a shift parameter: shifts follow the `.azr` vary flag. Shift penalties are measured from zero
(azure2-fit).

Covariance files: rows are the data rows of the `file <path> <count>` header, not the points left after cuts.
`relvector` terms are scaled by the model frozen at the run's first evaluation (`covarianceModel.out`): a fit that
starts far from its minimum minimizes a χ² built on the start's model. The CLI re-freezes after MIGRAD and
`--refreeze-covariance[=N]` repeats until the objectives agree; with older binaries score by mode 1 and refit from
the result's `param.sav` (`docs/source/user_guide/chi_squared.rst`). The data-scaled `vector` term for a
normalization is Peelle's puzzle: use norms (`normgroup`) or `relvector` instead.

## 3. Test segments and extrapolation grids

`<segmentsTest>` lines define extrapolation grids. A test line skips the data-only keywords (`sqrtshift`,
`shiftgroup`, `shiftfngroup`, `normgroup`, `absolute`, `shiftfn`) and still reads `sample`, `total` and `spin`
written after them. Grids, modes and the `intEC.extrap` cache are azure2-calculate.

## 4. Keys: how `<targetInt>` and parameter files find a segment

- A `<targetInt>` entry is matched to a segment purely by its numeric key column, whatever section the line lives
  in (`EData::SetUpTargetEffects`). `<segmentsData>` keys run 1..N in file order; `<segmentsTest>` keys are a
  separate 1..M counter. Appending a test line (by hand, `add_extrapolation` or `set_extrapolations`) can land on a
  key an unrelated data segment claims for a target effect and inherit its convolution. Robust layout: mirror every
  `<segmentsData>` line into `<segmentsTest>` one to one, in the same order, inactive by default (remapping the
  `isDiff` codes between the sections), and switch on only the lines needed.
- The same collision reaches the parameter file. A mode-3 run given a fit's `param.sav` applies
  `segment_<k>_energy_shift` and `_norm` of data segment k to test segment k, because parameters are matched by
  name: an (α,n) Legendre-coefficient test segment came out 5.8 keV (lab) high. For extrapolations pass a parameter
  file with every `segment_*` line removed, or check the output energy grid against the requested one.
- Deleting a `<segmentsData>` line renumbers the keys that `<targetInt>` and `param.sav` use; set a line inactive
  instead (azure2-fit).

## 5. The data file (`.dat`)

Four whitespace-delimited columns, lab frame, forward kinematics: E_lab (MeV), θ_lab (deg), cross section (b, or
b/sr if differential), uncertainty. The angle column is required even for angle-integrated data (a dummy value).
Phase-shift data carry the phase (deg) in column 3. Analysing-power data carry A_y (dimensionless, signed) in
column 3, and their angle column is **c.m.**, not lab. Data covariance files: `docs/source/reference/data_formats.rst`
("Data Covariance Files").
