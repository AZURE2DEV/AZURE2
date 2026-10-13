# 16N beta-delayed alpha emission

The smallest model with a beta-decay pair: 16N (2-) beta-decays to 1- states
of 16O that break up into 12C + alpha. Three 1- levels feed the alpha channel:
the bound 7.117 MeV state (fixed ANC), the 9.58 MeV resonance and a 20 MeV
background pole whose beta amplitude is **negative**. Pair key 1 is
alpha + 12C, pair key 2 the beta pair (pair type 20, Gamow-Teller channel);
the level lines were taken from the 12C+alpha evaluation, with the amplitudes
at one of its fitted points.

`data/bda.dat` is the model's own spectrum on a 0.6-3.0 MeV grid with 5 %
pseudo-noise from a fixed seed (`make_data.py`), so the chi-squared of the
shipped levels is of order the number of points and the relative comparison
in `run_tests.sh` is meaningful. Flipping the sign of the background pole's
beta amplitude alone changes it by orders of magnitude.

What it guards (`tests/pyazr/beta_sign_test.py`): `CNuc::GetTransformParams`
returned the absolute value of a beta amplitude until 2026-10-12, so every
consumer of the physical vector (the GUI's level write-back, the MCMC's
physical parameters, pyazr's `save_fit`) wrote a negative amplitude back
positive, and a snapshot of the 12C+alpha fit re-read with its beta-delayed
alpha chi-squared 750 -> 29724. `save_fit(verify=True)` compared physical
values, which both sides got wrong alike; it now also compares the re-read
amplitudes of the particle and beta channels of each level with the fit's.
