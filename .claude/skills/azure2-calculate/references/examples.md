# pyazr examples

The scripts in `pyazr/examples/` run on the models shipped beside them. Run `ls pyazr/examples` rather than trusting
this list to stay complete. Owner skill of each task in the last column.

| Topic | Script | What it shows | Owner |
|---|---|---|---|
| model and scheme | `print_scheme.py` | the level scheme of a model | azure2-setup |
| | `edit_scheme.py`, `edit_model.py` | add/remove levels and segments, write a new `.azr` | azure2-fit |
| | `deactivate_level.py` | switch a level off, file-level and at run time | azure2-fit |
| widths | `transform_widths.py` | reduced-width amplitudes → partial widths (with the shift-function denominator) | azure2-exchange |
| | `dimensionless_widths.py` | widths in units of the Wigner limit | azure2-check |
| observables | `angular_distribution.py` | Legendre coefficients at given energies (`pyazr.angular_distribution`) | this skill |
| | `sfactor_extrapolation.py` | S factor below the data and S(0) | this skill |
| | `decompose_cross_section.py` | level, interference and external-capture parts | this skill |
| | `reaction_rate.py` | N_A⟨σv⟩ by the engine (`pyazr.rate`), Gamow window, band from a covariance | this skill |
| | `tabulate_cross_sections.py` | adaptive σ(E) tables (`pyazr.tabulate`) | this skill |
| fitting and UQ | `fit_emcee.py`, `fit_zeus.py` | ensemble MCMC on the objective | azure2-uncertainty |
| | `per_dataset_chi2.py` | χ² per data set (several segments) | azure2-check |
| | `sensitivities.py`, `uncertainty_band.py` | d(observable)/d(parameter); covariance bands | azure2-uncertainty |
| snapshots | `save_fit_to_azr.py` | `Session.save_fit`: a fit baked into a `.azr` + `.sav` | azure2-fit |
| model internals | `channel_functions.py` | F, G, P and the hard-sphere phase | this skill |
| | `ec_integrals.py` | external-capture radial integrals and their cache | this skill |
| | `channel_radius_scan.py` | χ² against the channel radius | azure2-check |
| | `nuclear_potential.py` | the hybrid Coulomb + Woods–Saxon model per pair | this skill |
| data | `exfor_fetch.py` | `pyazr.sources.nds` end to end | azure2-setup |

Notes for running them:

- Set `OMP_NUM_THREADS` before importing numpy (the samplers set it to 1 per worker).
- Scripts that write output use the model's `output/` directory; run them on a copy of a production model.
- `angular_distribution` writes a temporary `.azr` and starts an engine per call; pass all energies in one call.
