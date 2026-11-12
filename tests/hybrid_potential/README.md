# hybrid_potential — the hybrid Coulomb model, per particle pair

A nuclear potential belongs to a **particle pair**: it bends the radial wave
functions of that channel and no other. This project exercises that, and the
precedence rules that go with it.

It is `13N` cut down to the three Meyer elastic-scattering segments
(`pair1 -> pair1`), so there is no capture and no external-capture integral to
compute — the run is quick and the only thing moving the numbers is the
Coulomb functions themselves.

## The potential record

The project is format 3: the potential of a pair is one `potential` record of `<potentials>`
(`docs/upgrade/project_format3.md` section 3.10a):

```
<potentials>
potential p12C
    pair            p+12C
    coulomb_radius  3.5 fm
    central         woods_saxon depth 20 MeV radius 3.6 fm diffuseness 0.6 fm
    matching_radius 15 fm
    step            0.0275 fm
</potentials>
```

A pair without a record has no potential (the photon pair here). The old file
(`tests/convert/fixtures/tests/hybrid_potential/hybrid_potential.azr`) states it as
an old `<potential>` block: a default (`V0=40`) that a pair without a section of its
own inherits, a `pair=1` section overriding it with `V0=20`, and `pair=2` with
`useHybridPotential=0`. The converter writes one record per particle pair from it
(`convert_spec_equivalence` compares the potential of every particle pair).

## What the recorded chi-squared proves

In the old file the default is `V0 = 40` and pair 1 overrides it with `V0 = 20`.
Pair 1 is the only particle pair, so the recorded number is the `V0 = 20` answer.
That single number discriminated three ways:

| if this were broken | chi-squared would be |
|---|---|
| — (correct) | 3.6598e6 |
| the per-pair override ignored, default used | the `V0 = 40` answer |
| the hybrid model not reaching the CLI at all | 4.3549e6 |

## Tolerance

This project carries its own `tolerance` file (2e-2) instead of the suite's
default 1e-3, because the hybrid model's Numerov integration is far more
floating-point sensitive than the plain Coulomb path -- the same compiler at
`-O0` instead of `-O3` moves segment 2 by 1.1e-3, and the CI runners spread to
2.2e-3, while `13N` moves by 3e-6 under the same comparison. The reasoning and
the measurements are written out in `tolerance` itself.

## check_per_pair.py

One project yields one number, which cannot show that the *per-pair* rules are
right. `check_per_pair.py` covers those directly: that pair 1's record is its
setting (there is no default in format 3), that switching a pair
off returns exactly the no-hybrid answer, that setting a potential through the
pyazr API equals loading the same potential from a file, and that the Coulomb
memo is keyed on the potential so a change to one is actually seen.

It also checks that the project as committed gives the total of
`expected/chiSquared.out` within `tolerance`, read from those files. It runs in
a temporary copy of this folder, in CTest as `pyazr_hybrid_per_pair` (about
20 s), or by hand with the compiled engine:

```bash
python3 tests/hybrid_potential/check_per_pair.py
```

## Below threshold

The bound ¹³N ground state (E_x = 0, 1.944 MeV below the p + ¹²C threshold, l = 1) is transformed with the shift
function of the decaying solution of the radial equation **with** the potential (`physics/HybridPotential.h`,
`DecayingSolution`), so S is continuous through threshold; its values against an independent integration are in
`physics/tests/data/hybrid_bound.csv` (`physics_hybrid_potential`, `tests/pyazr/hybrid_threshold_test.py`).

## One thing that will bite

The potential changes the penetrabilities and shift functions, and those are
what map physical widths to reduced-width amplitudes. So a change **re-derives
`params_rwa`**, and a vector captured before it no longer describes the model —
feed the old one back and the chi-squared is quietly wrong. Re-read
`m.params_rwa` after the call. Check 7 pins exactly this.

It is the same caveat as a channel-radius change, and for the same reason.

Re-pinned 2026-09-26: the shift-function fix b6cc41b removed ~1e-9 noise from S(E) that ShftFunc::EnergyDerivative amplified into dS/dE, which enters the Brune transformation of sub-threshold levels. One segment moved by 4.9e-5. It had already drifted by 6e-5 before b6cc41b, below the tolerance and never re-pinned.

Re-pinned 2026-09-26: the adaptive grid was reworked (geometric tails beyond each resonance's core, quantized anchors, only the finest lattice kept where lattices overlap; tests/thm_narrow_fold). The target-integration segments 1 and 2 integrate Rutherford-dominated elastic yields down to 1 keV. That integral is not converged at 50 points per width in either version. At 100 points per width the total is 2.4723e6 (new grid) and 2.4725e6 (old grid); at 50, the new grid gives 2.46612e6 and the old one 2.46215e6. Segments 1 and 2 moved by 1.4e-3 and 2.6e-3, and the total 2.46215e6 -> 2.46612e6 (+1.6e-3). The sub-point count is unchanged (1219 and 1821 per point).

Re-pinned 2026-09-28 (segment 3 only, 5e-5, inside the tolerance): the dS/dE step fix (ChannelFunc::DerivativeStep, commit ec416de). 2239.58 -> 2239.69.

Re-pinned 2026-11-12: the shared mesh (one energy mesh per compound nucleus, `docs/upgrade/energy_grid.md` §2.2, §7.1). Segments 1 and 2 2.09385e6 / 370034 -> 2.09231e6 / 369535, total 2.46614e6 -> 2.46411e6 (-8e-4). Today's engine refined eightfold gives 2.10131e6 / 372458 (2.47603e6): this project is the named exception of §7.1 (both engines are 5-6 % from the converged value at its own settings; the mesh is 1.2 x further at k = 1, equal from k = 4).
