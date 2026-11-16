# isospin_levels

Isospin-labelled levels across compound nuclei: the records of `<isospin>` (`isospin_reference`, `isospin_shift`,
`isospin_level`, `isospin_mix`, `isospin_factor`, `isospin_pair`; `include/IsospinTies.h`). dev's test of 2026-10-17
(`tests/isospin_levels` there), re-written for format 3.

`check.sh path/to/AZURE2` runs 58 checks (56 in a build without the sampler); its header lists them, (a) to (k). What it shows, in short: the factors are
Hale's A = 4 relations; a derived member is exactly the stand-alone level with its parameters (formal and Brune basis,
with and without `relativistic=1`); the refusals; four fitters recover two shifts and a mix from pseudo-data; the
analytic gradient and the Jacobian path agree with finite differences through hard and soft ties; photon channels stay
outside the ties; excited particles are ranked.

## Models

`models/*.azr` are dev's models converted by `AZURE2 --convert` (the `<isospin>` block becomes records; the `<levels>`
rows are unchanged):

| file | what |
|---|---|
| `iso.azr` | the toy A = 4 model: p+3He (4Li, the reference), n+t (4H), p+t, n+3He, d+d (4He); a 2- T=1 level in the three nuclei, a 2- T=0 level in 4He; shifts -0.86 and -0.64 MeV |
| `iso_mix.azr` | the same with a fixed isospin-forbidden d+d amplitude 0.05 MeV^1/2 on the T=1 level: the truth of the recovery fits |
| `iso_fit.azr` | every level parameter fixed; the shifts (from -0.5, -0.4) and the mix (from 0) free; data sets `nt.dat`, `ptn.dat`, `ddp.dat`, `ddn.dat`, which `check.sh` writes from the truth run |
| `iso_free.azr` | the same data, every reference parameter free, priors on the shifts (0.3 MeV) and the mix (0.1): the gradient checks |
| `iso_soft.azr` | `iso_free.azr` without the mix, the T=1 level soft-tied (0.05 MeV, 0.02 MeV^1/2) |
| `standalone_4H.azr` | 4H alone with one level; `check.sh` writes the 4H member's energy and amplitude into its `<levels>` row (run with `--no-transform`) |
| `bad_undeclared.azr`, `bad_radius.azr`, `bad_mix.azr`, `bad_T.azr` | refused: a level ID across nuclei without a record; n+t at 5.0 fm; a mix on an allowed channel; T=0 with a line in 4Li |
| `ecA.azr`, `ecM.azr` | dev's `tests/two_nuclei_ec` models: 12C(p,g0)13N alone, and with 12C(n,g0)13C |
| `ecA_iso.azr` | `ecA.azr` with every level in a one-nucleus block |
| `ecM_tied.azr` | `ecM.azr` with the ground states and the 1/2+ states of 13N and 13C tied |
| `ecM_exc1.azr`, `ecM_exc2.azr` | `ecM_tied.azr` with n + 12C*(7.654) in 13C, and also p + 12C*(7.654) in 13N |

## Differences from dev's test

- The parameters are named by their records (`shift_4H.shift`, `mix_1.amplitude`), not `isospin_shift_1_4` and
  `isospin_mix_1`.
- dev refused a soft tie with the least-squares fitters; here the soft ties are residual rows of the objective, and
  the check compares the fitter's soft penalty with MIGRAD's.
- A cross section out of the d+d pair is twice dev's (this branch's 1 + delta_12 of an identical entrance pair,
  `docs/upgrade/FOR_THE_OWNER.md`); the checks compare runs of one build, so none depends on it.
