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
