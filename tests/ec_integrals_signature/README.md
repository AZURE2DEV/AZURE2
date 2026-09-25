# ec_integrals_signature — a saved intEC file is not reused on a different grid

`output/intEC.dat` holds the external-capture amplitudes at the sub-point
energies of one integration grid. Energy straggling, a target thickness, the
adaptive-grid settings move those energies without changing how many there
are; a channel radius, a final-state energy, `--gsl-coul` or the hybrid
potential change the integrand at the same energies. The count check that
guarded reuse could see none of this, so a stale file was read back silently
and gave a wrong chi-squared that looked like physics.

AZURE2 now writes `intEC.dat.sig` beside the file: a 64-bit FNV-1a hash of
every energy an amplitude is evaluated at, the structure of each capture
block, the pairs involved and the Coulomb-function options
(`EData::ECSignature`). A file offered for reuse must match in count *and*
signature, otherwise a `WARNING` says why and the integrals are recomputed.

`check.sh` runs `tests/13N` cut down to one gas-target capture segment
(Artemov, target integration with external capture; sub-point density lowered
to keep it short), and checks that

1. a fresh run writes the signature;
2. the file is reused, without a warning, on the model it was built for, and
   gives the fresh chi-squared to the file's precision (intEC.dat stores six
   significant digits, so a reused file moves chi2 in the fifth or sixth);
3. with the 12C+p channel radius changed from 3.6 to 4.0 fm -- same number of
   amplitudes, so the count alone would let it through -- the file is
   rejected with the signature warning, and the chi-squared equals a fresh
   calculation's;
4. with straggling switched on the answer is the fresh one, whether the grid
   moved (recomputed) or not (reused). In this model it does not move: the
   straggling kernel changes, the sub-point energies do not;
5. a file with no signature (written before signatures existed) is
   recomputed, with a warning;
6. `--gsl-coul` invalidates a file computed with the default routine.

`run_tests.sh` runs it after the chi-squared comparisons; by hand:

```
./tests/ec_integrals_signature/check.sh path/to/AZURE2
```

It works in a temporary directory and leaves the tree untouched.
