# ec_relativistic_signature

Relativistic kinematics are part of the external-capture integral signature and of the in-process
node table key (audit 2026-10-25, checklist t3-25).

With `relativistic=1` in `<config>` every channel's eta, rho and flux factor change, so the
external-capture amplitudes change too (0.1 to 1 % through the penetrabilities). Up to dev b80b857
the switch was neither in the signature written beside `intEC.dat` (`EData::ECSignature`) nor in the
`ECTable` key (`TermKey`): a relativistic run handed a file written without the switch read it back
silently, the pattern of the long-wavelength signature bug fixed in 41a468b, and a pyazr session
evaluating a model and then its relativistic variant reused the nodes computed for the other.

The signature now adds "relativistic" when the switch is on, so every non-relativistic file keeps
its signature (nothing in the archive is recomputed), while a relativistic run no longer accepts a
non-relativistic file; the table key carries the switch as a field.

## Checks

`check.sh path/to/AZURE2`, on `tests/13N` cut down to its gas-target capture segment:

- (a) the relativistic chi2 differs from the non-relativistic one by more than 1e-4 (teeth);
- (b) the relativistic run given the non-relativistic `intEC.dat` warns that the file was computed
  for different inputs, recomputes, and equals the fresh relativistic chi2;
- (c) the non-relativistic run given its own file reuses it without a warning.

`tests/pyazr/ec_relativistic_test.py` (ctest `pyazr_ec_relativistic`): in one process, the
relativistic objective evaluated after the non-relativistic one equals that of a fresh process to
1e-9, and the two objectives differ by more than 1e-4.

dev b80b857: (b) already held, because the signature hashes the c.m. energy of every sub-point and
the switch moves those energies (the reviewer's "read back silently" needs a grid whose energies
coincide, which the field added now covers). The regression the old build does fail is the pyazr one:
its relativistic objective after a non-relativistic evaluation reuses the other's nodes.
