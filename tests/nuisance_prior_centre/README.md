# nuisance_prior_centre — the centre of a nuisance parameter's prior is the number written

`<value> fit prior <centre> <width>` on a nuisance parameter: the penalty is ((p − centre)/width)², with the centre
as written and the start at the value written, for every kind that can carry a prior (norm, the three shift terms,
abundance, `isotope_sum`, screening, areal density, Doppler temperature, sample thickness, stopping scale, a
kernel's parameter; `docs/upgrade/parameter_kinds.md` section 4). Until 2026-11-13 the sqrt(E) shift, the shift
function, the abundance, `isotope_sum`, screening, Doppler temperature and thickness measured the penalty from the
start (or from 0 / 1), and the areal density, the stopping scale and a kernel's parameter refused a centre away from
the value.

`check.sh path/to/AZURE2` moves each line of copies of projects of the tree so that value ≠ centre and compares the
prior terms of `chiSquared.out` with numbers put together by hand (the header of `check.sh` lists them); with
`PYAZR_ENGINE_DIR` set, `write_back.py` checks that a fitted screening potential (written in keV), Doppler
temperature, areal density and sample thickness are written back into their records in the unit written, the prior
untouched. No committed outputs: every number is in `check.sh`.

The older rows of `tests/prior_centre` (the `prior_centre` rows of `<parameterSettings>` read from old lines) are
another path to the same centres of a norm and a shift.

`continued.sh path/to/AZURE2` (run by `check.sh` as its check (h)): a fit continued from its fitted values computes
the fit it continues. For a fitted areal density (`tests/target_density/fit_prior.azr`), Doppler temperature
(`tests/doppler_fgm`'s `fit_prior`), stopping scale (`tests/stopping_yield/models/fit.azr`, densities fixed) and a
kernel's width (`tests/rectangular_resolution/models/fit.azr` with a prior), it fits, writes the fitted value into
the record (the value token only), and checks that the written project at the fit's parameters gives the fit's
`chiSquared.out`, that a second fit stays at the minimum, and that a start outside the range laid out from the centre
is laid out from the start with one line saying so. With `PYAZR_ENGINE_DIR`, `continued.py` adds `save_fit`'s own
check and compares the layout values the sessions report (`parameter_kinds.md` section 4, "The layout of a fitted
value"). On the build before the rule, the Doppler temperature's χ² differed by 6e-6 and every layout moved.
