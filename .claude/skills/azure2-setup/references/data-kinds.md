# Data sets of format 3: what each kind measures, and how to write it

A data set is a `dataset` record of `<data>`: a reaction, a quantity with its unit, a data file with the role of each
column, windows, a normalization, and the target, beam or detector that fold it. The grammar of every field and the
engine's refusals are in `docs/source/reference/format3.rst` (generated from the table the engine reads; `pyazr.kinds()
["data"]` gives the same as data); what each kind computes is in `docs/source/user_guide/segments.rst`. This file is
the chooser and the pitfalls.

## 1. The record

```text
dataset meyer_84
    reaction       12C(p,g)13N          # a sum: "12C(p,g)13N + 12C(p,g)13N[2.365]"; a ratio: "A / B"
    quantity       sigma
    unit           b
    file           data/meyer.dat       # relative to the project file
    columns        E_lab[MeV] value error_percent
    energy_window  0.1 2.5 MeV
    norm           1 fit prior 1 5%     # or: fixed, profiled, or the name of a shared normalization record
    shift          offset 0 MeV fit prior 0 0.002    # or the name of an energy_scale record
    covariance     data/meyer_cov.txt   # optional: correlated errors (azure2-uncertainty)
    target         C_foil               # optional: folds the model (azure2-targets)
    doi            10.1016/...          # provenance: source, exfor, doi
```

From Python: `Project.add("dataset", name, **fields)` (or `add_dataset`), `Project.set(name, field=value)`,
`Project.set(name, field=None)` removes a field; every edit is checked by the engine and refused with a
`FormatError` naming the field. `AZURE2 --describe model.azr --record meyer_84` shows what the engine made of it:
the observable, the reaction tree with the pair keys it resolved to, the rows read and kept inside the windows.

- **Columns carry their frame and unit**: `E_lab[MeV]`, `E_cm[keV]`, `theta_lab[deg]`, `theta_cm[deg]`, `value`,
  `error`, `error_percent`, `order` (Legendre), `window_lo`/`window_hi` (a beam profile). A file of c.m. energies is
  read as such; the value's frame is the angle's (a dσ/dΩ with `theta_cm` is a c.m. value). Every output is c.m.
- **An angle-integrated file needs no angle column**; a differential one needs one, or an `angle_window` that
  selects rows of a file holding several angles.
- **Parameter names follow the record**: `<data set>.norm`, `<data set>.shift`, `<target>.areal_density`, ...
  (`s.parameters.table()`). A parameter file of a project before its conversion is renamed with
  `AZURE2 --convert-params`. Renaming a record (`Project.rename`) renames its parameters, so a `param.sav` written
  before the rename no longer matches: the program warns `W_PARAMETER_FILE_NAMES`.
- **Sharing**: a `normalization` record named by several data sets is one norm for all (an excitation function and
  the angular distributions integrated from it); an `energy_scale` record (`offset`, `sqrt`, `function` with `param`s)
  is one calibration for all data sets that name it in `shift`.
- **`active no`** keeps a data set in the file but out of the calculation. Parameters are named by record, so
  removing one renames nothing; the numbers of `chiSquared.out` (a data set's position, the `key` of `--describe`)
  do move.

## 2. Which quantity

| quantity | measures | note |
|---|---|---|
| `sigma` | σ(E), angle-integrated | capture: `multipole E1` or `E2` for one part; on a thin layer the layer average |
| `sigma_total` | neutral-projectile total σ from the optical theorem | cannot leave an open channel out; capture not in it |
| `sigma_nuclear` | the nuclear elastic integral π/k² Σ g\|T\|² | not measurable for a charged pair: diagnostics |
| `dsigma_domega` | dσ/dΩ at the row's angle | frame of the value = frame of the angle column |
| `ratio_to_rutherford`, `ratio_to_mott` | dσ/dΩ over Rutherford (Mott for identical particles), same c.m. angle | |
| `ratio_to_reference` | σ over a tabulated reference σ (a monitor) | `reference <file>` in the data set's energy role |
| `sfactor`, `modified_sfactor` | S(E) = Eσ exp(2πη); S* with exp(gE), `shape g` | the engine's η and atomic masses |
| `legendre` | c.m. Legendre coefficients | `expansion a_over_a0|absolute|rsl|2l2|4pi`, an `order` column |
| `analyzing_power` | A_y of a spin-1/2 projectile | folded as A_y σ over σ; no norm |
| `phase_shift` | δ_{J l s} (deg) | `J`, `l`, `s` fields |
| `yield` | reactions per particle (or per charge) on a layer or an infinite target | azure2-targets |
| `transmission`, `neutron_yield` | T = exp(−nσ) and the self-shielded yield of a sample | azure2-targets |
| `beta_delayed_spectrum` | particles after β decay, arbitrary units | β pair (type 20) |
| model data (`datum` records of `<modelData>`) | `resonance_strength`, `resonance_energy`, `total_width`, `partial_width`, `branching_ratio`, `anc`, `scattering_length` | a measured property of a level or a pair, a term of the objective |

A secondary γ ray of a state the exit pair is left in (the 4.44 MeV γ of ¹²C after (p,α₁)) is the data set's
`secondary_gamma L <l> Ic <J> delta <δ>` field; it is not the primary angular distribution, and an extrapolation of
the primary reaction does not reproduce it.

## 3. Totals, sums and ratios

- **A total is summed over every channel it measures.** A 4π (α,n) total yield written as (α,n₀) alone looks like
  missing physics above the next threshold (one case: Δχ² = 3.35 × 10⁴ once the missing exits were added). Write
  the sum in the reaction (`22Ne(a,n0)25Mg + 22Ne(a,n1)25Mg`); for a neutron total use `sigma_total`, which cannot
  leave a channel out. Before refitting, evaluate the corrected definition at the current parameters: the right
  definition usually wins already. Two total measurements that track each other across a threshold with no step
  measure the same quantity.
- **Two data sets on one data file are two data sets**, each compared with its own reaction: a total written as several
  data sets with different exits on one file compares partial cross sections with the total. The converter of old
  projects writes them as they were; look for one file named by several data sets with different reactions (`awk
  '$1=="file"{c[$2]++} END{for(f in c) if(c[f]>1) print f}' model.azr`) and replace them by one data set with the sum.
- **A ratio** (`reaction A / B`, unit `1`) takes each term with its own lab→c.m. factor; it has no normalization
  meaning (fix the norm at 1).
- `AZURE2 --print-predictions model.azr` prints the expression tree each data set evaluates (`sum(...)`,
  `ratio(...)`, `abundance(...)`): read it for every composite.

## 4. Several compound nuclei: natural targets

- A natural target is a `target` record with one `isotope` field per target nucleus (`isotope 20Ne 0.9048 fixed`,
  `fit prior ...` to fit it) and `isotope_sum`; a data set on it names the sum of its reactions, each term weighted
  by its nucleus's abundance. The nucleus of a pair is inferred from its masses and charges.
- Build it with `pyazr.merge_projects(projects, out, isotopes={...}, shared=[{0: "vaughn", 1: "vaughn"}])` rather than
  by hand: the engine merges the `<levels>` rows (pair keys and level IDs moved past the earlier members'), keeps
  every record (renaming collisions) and adds the target; member 0's `<config>`, `<thm>`, `<parameterSettings>` and
  `<mcmc>` are kept. `pyazr.rst`, "Several compound nuclei".
- Cross-check a merge: each nucleus alone (the other abundances 0) reproduces its own project's χ² to 1e-10
  (`tests/natural_target`). Reference model: `tests/natural_neon`.

## 5. Provenance of the data

Record `source`, `exfor` and `doi` on every data set (fetching and converting from EXFOR: `nds-exfor.md`). Check the
paper for what was measured (a ratio, a yield, a total, a secondary γ ray), the frame of the energies and angles, the
error model (statistical, common systematic, shape; their covariance file: azure2-uncertainty) and whether a point
is a derived value (a paper's own S(0) or extrapolation in EXFOR is not a measurement: never fit it with the data it
came from).
