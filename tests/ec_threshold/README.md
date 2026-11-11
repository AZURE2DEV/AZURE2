# ec_threshold — external capture across a charged-particle threshold

A 16O(n,g)17O(g.s.) model where the capture runs partly through the alpha + 13C channel (channel capture). That
channel opens at E_cm(n) = 2.21556 MeV. The data grid is dense across that threshold. `make_models.py` (deleted; see "How the models were made") writes the
model; `check.sh` runs it with the external-capture energy tables (default) and with `--no-ec-table`, and checks:

- (a) every cross section is finite and positive;
- (b) the excitation function is smooth: each point within 1e-5 of the line through its neighbours;
- (c) the two routes agree to 1e-6.

**What it caught (fixed 2026-10-24).** Within about 1 keV above an open charged channel's threshold, the Sommerfeld
parameter exceeds about 110. The penetrability then underflows while G overflows, and the channel-capture integral
sqrt(P) int (G + iF) W came out as 0, as an overflow, as NaN, or (from eta ~ 300) finite but about 30 % low.

In this model the old build gave 1e177 b with direct integrals and 1e305 b with the energy tables. In a fold of
16O(n,g) data with a 200-point resolution kernel (13C+a, Holt), one sub-point 0.2 keV above the threshold made two
data points NaN.

`ECIntegral::CalculateDirect` now interpolates channel capture linearly between two values: its threshold value
(the closed-channel Whittaker-ratio form at -E_c/100) and the open form at E_c, the energy where eta = 100. This is
continuous at both ends. Direct capture from such a channel is suppressed by exp(-pi eta) < exp(-300), so it is set
to zero.

Two related changes: the external-capture cache now refuses non-finite amplitudes (it used to refuse only NaN), and
`EffectiveCharge` evaluates int_0^y j_L(t)/t dt in closed form (E1, E2) with a series at small y. The adaptive
quadrature it replaces threw GSL roundoff errors on the long Whittaker tails near a threshold.

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show 3e9b5298:tests/ec_threshold/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Model for tests/ec_threshold: 16O(n,g)17O(g.s.) across the alpha + 13C threshold of 17O.

Pairs: n + 16O (key 1, entrance), alpha + 13C (key 2; opens at E_cm(n) = 6.35869 - 4.14313 = 2.21556 MeV), photon to the
17O 5/2+ ground state (key 3, E1 external capture).  The g.s. has ANCs in both particle channels; a 5/2- level at
Ex = 9 MeV carries n (l = 3), alpha (l = 2) and an E1 photon width, so the capture includes channel capture through
the alpha + 13C channel, whose Coulomb functions over- and underflow within ~1 keV above its threshold.
The data grid is dense across the threshold (E_cm 2.195-2.235 MeV, 0.25 keV steps, plus points 0.01-0.5 keV around it).

Run from this directory: python3 make_models.py   (writes models/)
```
