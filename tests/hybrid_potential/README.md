# hybrid_potential — the hybrid Coulomb model, per particle pair

A nuclear potential belongs to a **particle pair**: it bends the radial wave
functions of that channel and no other. This project exercises that, and the
precedence rules that go with it.

It is `13N` cut down to the three Meyer elastic-scattering segments
(`pair1 -> pair1`), so there is no capture and no external-capture integral to
compute — the run is quick and the only thing moving the numbers is the
Coulomb functions themselves.

## The `<potential>` block

```
<potential>
useHybridPotential=1     # the master switch
useAdaptiveGrid=1
potentialType=0          # 0 = Woods-Saxon, 1 = Gaussian
V0=40                    # the default: what a pair with no setting inherits
R=3.6
a=0.6
pair=1                   # ... and everything below applies to pair 1 alone
useHybridPotential=1
potentialType=0
V0=20
R=3.6
a=0.6
pair=2
useHybridPotential=0     # the photon pair opts out
</potential>
```

Keys before the first `pair=` are the default. A `pair=` line opens a section
that starts from that default and only states what differs. A file with no
`pair=` line is read exactly as it was when the model was global, which is what
keeps every existing project working.

The same parser serves the GUI, `--no-gui` and pyazr
(`Config::ReadPotentialBlock`), so the three cannot drift.

## What the recorded chi-squared proves

The default is `V0 = 40` and pair 1 overrides it with `V0 = 20`. Pair 1 is the
only particle pair, so the recorded number is the `V0 = 20` answer. That single
number discriminates three ways:

| if this were broken | chi-squared would be |
|---|---|
| — (correct) | 3.6598e6 |
| the per-pair override ignored, default used | the `V0 = 40` answer |
| the hybrid model not reaching the CLI at all | 4.3549e6 |

## check_per_pair.py

One project yields one number, which cannot show that the *per-pair* rules are
right. `check_per_pair.py` covers those directly: that a pair's own setting
beats the default, that a pair without one inherits it, that switching a pair
off returns exactly the no-hybrid answer, that setting a potential through the
pyazr API equals loading the same potential from a file, and that the Coulomb
memo is keyed on the potential so a change to one is actually seen.

Run it from this directory:

```bash
python3 check_per_pair.py
```

## One thing that will bite

The potential changes the penetrabilities and shift functions, and those are
what map physical widths to reduced-width amplitudes. So a change **re-derives
`params_rwa`**, and a vector captured before it no longer describes the model —
feed the old one back and the chi-squared is quietly wrong. Re-read
`m.params_rwa` after the call. Check 7 pins exactly this.

It is the same caveat as a channel-radius change, and for the same reason.
