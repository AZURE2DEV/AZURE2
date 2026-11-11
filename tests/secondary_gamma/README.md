# tests/secondary_gamma

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show 3093e074:tests/secondary_gamma/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Build the secondary-gamma (UPOS) test models and the expected angular distributions.

13C + alpha in 17O.  Exit pairs: n + 16O(6.13, 3-) (key 5; channel spins 5/2 and 7/2, E3 to the 0+ ground state)
and alpha + 13C(3.68, 3/2-) (key 6; channel spin 3/2, M1/E2 to the 1/2- ground state).  Single-level models are
evaluated exactly at the level energy, where the exit T-matrix elements of one l' are proportional to
sign(Gamma) sqrt|Gamma| with the same energy factor and phase, so the distribution follows from the widths alone
(truth.py).  Every model also has an angle-integrated segment of the same reaction, which normalises the UPOS cross
section: W(theta) = 4 pi (dsigma/dOmega_gamma) / sigma.

  spin2_plus / spin2_minus / spin2_zero : 5/2+, l' = 1 in s' = 5/2 (Gamma 20 keV) and 7/2 (+10, -10, 0 keV)
  spin2_multil                          : as plus, with l' = 3 in both channel spins as well (-5 and +8 keV)
  mixing                                : 5/2+ -> alpha + 13C(3/2-), l' = 1 and 3; M1/E2, delta = -0.094
  twolevel                              : 5/2+ and 3/2+ -> alpha + 13C(3/2-) (one channel spin), J1 != J2 interference;
                                          expected = the pre-fix engine (dev 6c6932d), unchanged by the fix
  grad                                  : 5/2+ and 3/2+ -> n + 16O(3-), free energies and widths, UPOS data at
                                          three energies -- for the AZURE_GRAD_CHECK gradient test

Run from this directory: python3 make_models.py  (needs numpy, scipy and sympy; writes models/ and expected/).
```
