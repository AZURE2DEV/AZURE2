# THM data sets and the `<thm>` block: mechanics

How a THM measurement enters an AZURE2 project and how to edit it. Physics, equations and pinned values:
`docs/source/theory/thm_implementation.rst`. The THM path is the modified R-matrix, half-off-shell (HOES) one: the
reaction is induced on a nucleon or cluster x bound in a carrier a = (x s), so it proceeds below the Coulomb barrier
without penetrability suppression; the entrance channel is off the energy shell, the exit channel on it.

## THM data sets

- **A data set is THM when it names an experiment** (`thm <experiment>` or `thm global`); its quantity is the ordinary
  one (`sigma`, `dsigma_domega`, also `analyzing_power`), so THM composes with every observable.
  `AZURE2 --print-predictions` prints such a term with `thm_` before it; `Session.thm_experiments()` lists the
  experiments with their data sets.

- **The analytic gradient does not cover THM.** The adjoint differentiates the on-shell T-matrix observable, not the
  HOES one. THM data sets are left out of the analytic accumulation; their contribution to the energy, width and
  normalization derivatives comes from central differences (`AZURECalc::Gradient`, guarded by `haveTHM`). All other
  data sets stay analytic.

## Widths given as reduced-width amplitudes

A channel line can carry `gammaIsRWA` (field 33). Its `gamma` column is then an amplitude γ in MeV^½, not a partial
width Γ in eV or an ANC in fm^−½: useful for a level known only from a THM measurement. The convention is kept on
output (`CNuc::GetTransformParams` returns the amplitude unchanged), so a fit written back never flips it.

```python
ch.gamma_is_rwa            # AzrChannel: is this channel's gamma an amplitude?
w.width_is_amplitude       # Parameter: the same, from the live API
p.binding_energy           # pair: THM binding energy B_xs (MeV), 0 if not THM
```

Set both as Project edits, not by editing the row text: `Project.set_pair(key, binding_energy=B)` writes B on every
row of the pair (a particle pair, B >= 0); `Project.set_channel(level, pair, l, s, width=g, width_is_amplitude=True)`
sets the flag of one channel row with its amplitude g in MeV^½ (a particle channel; the value is not converted, so
give it in the new convention in the same call). `Project.pairs()[k].binding_energy` and
`Project.level(label).channel(pair, l, s).amplitude` read them back without a session.

- **Never convert an rwa-declared width by hand.** `Session.save_fit` respects the flag; hand-written
  amplitude-to-width arithmetic corrupts the file by orders of magnitude.
- **`Parameter.wigner_limit` is γ²_W in MeV**, the limit on the reduced width squared. It matches neither convention
  `value` can be in. Form θ² with `Session.dimensionless_widths(x)` (azure2-check, `widths_check.py`), not by dividing
  by it.

## The `<thm>` block

Optional, anywhere in the `.azr`, one `key=value` per line. An unknown key or value, or a missing `</thm>`, makes
AZURE2 exit with `ERROR: <thm> ...` before computing anything.

| key | values | meaning |
|---|---|---|
| `entranceL=` | `incoherent` (default), `coherent` | sum over entrance l |
| `vertex=` | `constant` (default), `perlevel` (alias `real`), `onshell` | boundary in the plane-wave vertex |
| `kinematics=` | `lacognata` (default), `triple`, `kf3body`, `lambda32` | factor K(E) left in the data |
| `coulombIntegral=` | `0`, `1` | x–A Coulomb term outside a |
| `spectatorEnergy=`, `spectatorEnergy[<pair>]=` | MeV | mean spectator energy T_s |
| `experiment[<name>]` | `beam= target= spectator= Ebeam= background=linear ...` | shared profiled norm, background |

- `kinematics=`: choose from what the data paper divided by. Full three-body KF → `kf3body`; λ₃/λ₂ → `lambda32`; only
  |φ|² → `triple`; La Cognata's formula or unknown → `lacognata`. The wrong one distorts relative resonance heights.
- `vertex=constant` puts B_c = S_c(E₁) of the lowest level of each J^π in the vertex. It is independent of level order
  and identical with Brune on or off. `perlevel` puts each level's S_c(E_λ) in the vertex (as the mrmpy code does) and
  depends on the representation once levels of one J^π interfere; use it only to reproduce numbers made with it.
  `tests/18O_p_a_thm` reproduces La Cognata et al., ApJ 723 (2010) 1512 with `constant` (5 % rms at the peaks) and not
  with `perlevel` (23 %). Regression of all keys: `tests/thm_options`.
- A weight w(E) multiplying the THM model of one data set is that data set's `thm_weight <file>` field.
- `experiment[<name>]`: the data sets that name it (`thm <name>`) share one profiled norm plus a profiled background;
  results in `output/thm_experiments.out` and `session.thm_background(name)`; test `tests/thm_experiment`.
  Per-experiment keys: `references/model-choices.md`.

### Editing it from Python

`pyazr.Project` (format 3) applies the engine's rules and keeps comments and untouched lines:

```python
p = Project.open("17O.azr")
p.thm()                                    # the global options and the experiments
p.set_thm(kinematics="kf3body")            # None removes an option
p.set("o17_guardo", thm_weight="R_of_E.dat")   # a THM data set's own weight table (path relative to the project)
p.set_thm_experiment("E", background="...")    # set keys of an experiment's line; None removes one
p.set("o17_guardo", thm="E")               # a data set joins an experiment by its own field
p.remove_thm_experiment("E")
```

A refusal raises `pyazr.FormatError` with the engine's message. By hand, append `"\n<thm>\n...\n</thm>\n"`: an `.azr`
may lack a final newline, and `<thm>` must start its own line.

## Reference projects

| test | system | what it pins |
|---|---|---|
| `tests/17O` | ¹⁷O(n,α)¹⁴C, Guardo/Sergi et al., 23 points, one THM data set | χ² = 96.69; neutron entrance, rwa widths |
| `tests/7Li_p_a` | ⁷Li(p,α) | χ² = 2138.51; charged entrance, two entrance l |
| `tests/6Li_d` | ⁶Li(d,α) | χ² = 682.133 |
| `tests/18O_p_a_thm` | ¹⁸O(p,α) | χ² = 2191.98 and `check.sh` against the published band; two interfering levels of one J^π |
| `tests/thm_rwa_grid`, `tests/thm_options` | | `check.sh` |

They exist because THM shares `Segment`, `Unit`, `AZURECalc` and the gradient code with the on-shell path and
diverges only at the observable.
