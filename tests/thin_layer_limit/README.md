# thin_layer_limit — a target that loses less than 1 eV

A folded point's nodes span its kernel's support.  A support narrower than `MeshZeroWidth()` (1e-6 MeV, `MeshPlanner.h`)
is a single node at its upper edge, the point's own energy E0 for a target without a beam spread.  The fold of a single
node is the thin limit of the effect's kernels (`Kernel::thinLimit`, applied in `Compose` of
`experiment/src/Kernels.cpp`): the target-averaged cross section σ(E0) of a thin target, 1e-24 σ(E0) N for a yield on a
layer (reactions per incident particle; σ in b, N in atoms/cm²).  Before, every kernel returned 0 for one node, and a
yield on such a layer was refused at the set-up as a layer "the beam does not leave".

`models/thin_limit.azr`: the 15N(p,α)12C model of `tests/layer_yield` (one broad 1/2+ level), a layer of 15N with
ε = 4.5e-22 x^-0.6 lab MeV cm² per nucleus (x = E/MeV), three energies E_lab = 0.6, 0.7, 0.8 MeV, as a yield on the
layer (`thin`), as a target-averaged cross section (`thin_sigma`) and bare at E0 and 1 keV either side (`bare`).
`check.sh` replaces its areal density (1e14 in the file) by 33 values from 1e13 to 1e17 atoms/cm² (8 per decade); the energy loss ε N is 1 eV
near N = 1.7e15.

The tolerance.  Above the boundary the depth integral gives σ at the depth-averaged energy, E0 − T/2 for a thin layer
(T = ε N the energy loss), so Y/(1e-24 N) − σ(E0) ≈ −(T/2) dσ/dE; the check allows T |dσ/dE| (a factor 2 of margin)
plus 1e-9 (the output's 11 digits).  Below it the single node is E0 itself, and both quantities equal σ(E0) to 1e-9.
The step between the two neighbouring densities on either side of the boundary is then below the same tolerance: the
yield over N is continuous through 1 eV (measured: 6e-6, 8e-6 and 9e-6 at the three energies, half the bound).

**With an energy shift** (`ShiftedFold.h`, the paragraph on single nodes).  The check also runs the model at 1e14 with a
shift on the two folded data sets (0, ±2 keV, ±20 keV fixed; ±20 keV varied without a prior; +20 keV and −2 keV varied
with a 10 eV prior) and through 1e14–1e16 with +20 keV fixed and varied, against the bare data set calculated at the
shifted energies (its data file shifted, the data set not).  A fixed shift keeps the single node at the shifted point's
own energy: σ(E0 + s) to 1e-11.  A varied shift without a prior has a planned support 7 times the layer's (its reach is
3 widths each side), still one node below about 2.4e14 at 0.6 MeV, its upper edge: σ at up to 7 T above E0 + s, 2e-6 off
at 1e14 (3 T|σ′|/σ).  A 10 eV prior gives a planned support of many nodes and the depth integral (3.2e-7 = T/2 |σ′|/σ).
Before, every shifted single-node point but s = 0 fixed was NaN.

Not in the check: a prior of keV on this layer.  The planned support (± 3 prior widths) keeps the smooth step of the
layer's own grid (its width over 20 points, quantized: about 3e-9 MeV), so a 3 keV prior asks for millions of nodes
(28 GB after 3 minutes, stopped; 0.1 keV: 10.6 s and 2.4 GB).  And at −20 keV with a 10 eV prior the target-averaged
σ is 6.2e-6 off (the yield 3.2e-7, as expected; +20 keV, −2 keV and 1e16 agree): not explained, listed in the
follow-ups' report.
