# thm_experiment — THM experiments: shared norm and background

`check.sh` runs `tests/18O_p_a_thm` (two THM segments with free norms, all
parameters fixed, ~0.1 s a run) with `experiment[<name>] key=value ...` lines
in a `<thm>` block appended to a temporary copy of the .azr
(docs/source/theory/thm_implementation.rst, "THM experiments"):

| case | expectation |
|---|---|
| (a) `experiment[A] segments=1` + `experiment[B] segments=2` (B with kinematics) | `chiSquared.out`, `normalizations.out`, `AZUREOut_aa=1_R=2.out` byte-identical to the run without the block; the kinematics summary prints B(d) = 2.2246 MeV |
| (b) `experiment[A] segments=1,2` | n\* and chi2 in `thm_experiments.out` equal sum S_mm / sum S_md and sum S_dd - (sum S_md)^2 / sum S_mm computed here from the models of run (a) and the data files, to 1e-8 (measured ~1e-11; the models are read at 10 digits); both segments carry n\* (2.96005e-05, chi2 2241.17 against 2228.94 with a norm each) |
| (c) `background=linear` on data made as s m + a0 + a1 E | n = 1/s, b0 = a0/s, b1 = a1/s recovered to 1e-8, chi2 ~ 1e-17, the written curve (model + b(E)) equals the scaled data; the same data without the background give chi2 > 1 |
| (d) reserved keys `theta`, `distortion` (`lineshape`: tests/thm_lineshape; `ps`: tests/thm_spectator_window); unknown key; unknown nuclide; bad `Z,A,mass`; partial kinematics; bad background, segment list or `Ebeam`; no `segments=`; a key given twice; a segment in two experiments, not existing, not THM or with a fixed norm; kinematics that do not make the entrance pair | `ERROR: <thm> experiment[...]` and a non-zero exit |

The derivatives through the profile are checked in
`tests/pyazr/thm_experiment_test.py` (residual Jacobian and chi2 gradient
against central differences) and in `tests/thm_band/check.sh` (the CLI band
of a THM experiment with a linear background against finite differences of
the output). About 1.5 s.
