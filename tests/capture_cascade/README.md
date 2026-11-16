# capture_cascade — the secondary gamma ray of a capture final state (secondary_gamma on a capture reaction)

A data set with `secondary_gamma` whose exit is a capture pair is the secondary gamma ray of the capture final state
B after an unobserved primary capture photon: A(a, gamma_1)B, B -> C + gamma_2, gamma_2 detected. An unobserved photon
summed over its direction and both helicities leaves no interference between multipoles and orients B as a spin-0
particle of orbital angular momentum L would, so Eq. (22) of Brune and deBoer, PRC 102, 024628 (2020) applies with
l' = L, I_b = 0 and s' = I_B = J_f (the extension named in that paper's conclusions). Internal and external capture
amplitudes of one multipole interfere. Direct feeding of B only.

**Models** (`models/`, dev's, written in format 3 by `AZURE2 --convert` with 0 notes), each with four data sets of
one capture pair: the cascade distribution, the angle-integrated cross section, the primary photon's distribution
and an angle-integrated cascade data set.

- `single_e1`, `multipoles`, `twoj_plus`, `twoj_minus`: resonant capture in 13C + alpha (one level and one
  multipole; two primary multipoles of one level with a mixed secondary transition; two levels of different J through
  one multipole, with either relative sign);
- `ec_s`, `ec_p`: external capture alone on a 1/2- target, into an s state (isotropic only as a coherent sum over J)
  and a d -> p single-particle transition;
- `ecres_both`: a resonance and external capture together, at three energies across the resonance;
- `grad`: a fit with a free level and a free ANC, for the gradient.

**Truth** (`expected/`): W(theta) = 4 pi (dsigma/dOmega)/sigma from sums over magnetic substates that use neither the
cascade coefficient nor R_k (dev's `make_models.py`, which also wrote the models); for the resonant and the
external-capture models the PRIMARY distribution of the same amplitudes is checked too, which ties the conventions to
code that existed before. `verify.py` compares (1e-7; 2e-11 to 8e-11 measured) and checks that the angle-integrated
cascade data set equals the capture cross section (1e-12). `check.sh` adds the analytic gradient against finite
differences (`AZURE_GRAD_CHECK`, worst 9e-9 of the largest component).

**Origin.** dev's `tests/capture_cascade` (2026-11-15). This engine's values equal dev's on every model (1e-10; two
points of a primary distribution that vanish at 0 and 180 degrees differ in their rounding noise, 1e-23 against
3e-9).

**Run.** `./tests/capture_cascade/check.sh build/src/AZURE2` (15 s; `verify.py` needs numpy and is skipped without).
