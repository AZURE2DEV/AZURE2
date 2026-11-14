# pyazr examples

The scripts in `pyazr/examples/` run on the projects of the repository (`tests/upgrade/T1`, `T2`, `T3`, ...), each
from the repository's root by the command its docstring states; the full list by task, with the commands, is
`docs/source/user_guide/pyazr_examples.rst` (`tests/pyazr/examples_test.py` runs every one; consistency check 58
keeps the page, the test, the directory and the names of this list one set).  Owner skill of
each task in the last column.

| Topic | Script | What it shows | Owner |
|---|---|---|---|
| model and scheme | `print_scheme.py` | the level scheme of a model | azure2-setup |
| | `edit_model.py` | add/remove/switch off levels and data sets, a prior and limits on level parameters, write a new `.azr` | azure2-fit |
| | `project_from_python.py` | a project built from records on another project's nuclear part | azure2-setup |
| | `refusals.py` | what a refusal looks like and how a script reads its three parts | this skill |
| widths | `transform_widths.py` | reduced-width amplitudes → partial widths (with the shift-function denominator) | azure2-exchange |
| | `dimensionless_widths.py` | widths in units of the Wigner limit | azure2-check |
| observables | `angular_distribution.py` | Legendre coefficients at given energies (`pyazr.angular_distribution`) | this skill |
| | `sfactor_extrapolation.py` | S factor below the data and S(0) | this skill |
| | `decompose_cross_section.py` | level, interference and external-capture parts; a level off or alone at run time | this skill |
| | `reaction_rate.py` | N_A⟨σv⟩ by the engine (`pyazr.rate`), Gamow window, band from a covariance | this skill |
| | `tabulate_cross_sections.py` | adaptive σ(E) tables (`pyazr.tabulate`) | this skill |
| fitting and UQ | `fit_engine.py` | the engine's MIGRAD/HESSE/MINOS, `Session.save_fit` and the read-back | azure2-fit |
| | `fit_scipy.py`, `fit_iminuit.py` | the objective under scipy and iminuit | azure2-fit |
| | `fit_emcee.py` | ensemble MCMC on the objective, a rate per sample | azure2-uncertainty |
| | `nuisances.py`, `covariance_builder.py`, `prior_builder.py` | norm priors and shared shifts; a data covariance; a prior covariance | azure2-fit |
| | `per_dataset_chi2.py` | χ² per data set (several segments) | azure2-check |
| | `sensitivities.py`, `uncertainty_band.py` | d(observable)/d(parameter); covariance bands | azure2-uncertainty |
| model internals | `channel_functions.py` | F, G, P and the hard-sphere phase | this skill |
| | `ec_integrals.py` | external-capture radial integrals and their cache | this skill |
| | `channel_radius_scan.py` | χ² against the channel radius | azure2-check |
| | `nuclear_potential.py` | the hybrid Coulomb + Woods–Saxon model per pair | this skill |
| data | `exfor_fetch.py` | `pyazr.sources.nds` end to end | azure2-setup |
| | `natural_target.py`, `gamma_decay.py`, `thm_dataset.py` | a natural target; a `gamma_decay` record; a THM data set | azure2-targets, azure2-thm |
| evaluations | `endf_evaluation.py`, `endf_start_model.py` | an ENDF-6 file read and checked; a start model from one | azure2-exchange |

Notes for running them:

- Set `OMP_NUM_THREADS` before importing numpy (the samplers set it to 1 per worker).
- A script takes its project as the first argument and writes only where `--out` (or `--plot`, `--csv`) says, or
  into a temporary directory it prints; nothing appears beside the project.
- `angular_distribution` writes a temporary `.azr` and starts an engine per call; pass all energies in one call.
