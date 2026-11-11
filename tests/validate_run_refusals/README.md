# validate_run_refusals

`AZURE2 --validate` predicts what a run refuses: each variant of two small format-3 projects is validated and then
run (menu 1), and the check wants the diagnostic code where the run refuses (`check.sh` lists the cases).

- `models/rth.f3`: alpha + 12C ratio to Rutherford on the hard sphere, from `tests/ratio_to_rutherford` (data
  `models/data/rth.dat`, its `aC_rth_cm.dat`): rows that do not parse (`E_DATA_ROW`), a covariance file missing or
  with a variance 0 (`E_COVARIANCE_FILE`), an uncertainty 0 (`W_DATA_ERROR`, a warning: the chi^2 skips the point).
- `models/capture.f3`: 16O(p,g)17F(0.495 MeV) with E1 external capture, `tests/ec_final_state_listed/models/listed.azr`
  converted to format 3 (data `models/data/ecm.dat`): the final-state level removed, doubled within 1 keV
  (`E_FINAL_STATE_NOT_LEVEL`), or shared by a second photon pair (`E_FINAL_STATE_SHARED`).

The models end in `.f3`, not `.azr`, so that `tools/runall.sh` does not run them as projects; `check.sh` copies them.
