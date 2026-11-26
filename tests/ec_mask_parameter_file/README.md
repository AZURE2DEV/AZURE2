# ec_mask_parameter_file — the `--ignore-externals` mask follows a parameter file

`check.sh` (run by `run_tests.sh`) checks, on the 13N reference model (`../13N`, 12C(p,γ)) and the E1 width of its
3/2⁻ level at 3.5032 MeV (`width_5_2`), that a parameter file is input for the mask of
`docs/upgrade/external_capture.md` section 3:

1. **A 0 in the parameter file is a 0 in `<levels>`.** The project with `width_5_2` exactly 0 in the file gives the
   capture cross sections of the project whose `<levels>` holds 0 there (limit 1e-12 with the same file, 1e-6 against
   the run without a file, whose other parameters have all their digits), and the run says so in a note.
2. **A non-zero value in the file is not masked.** With 0 in `<levels>` and the non-zero amplitude in the file the
   result is that of the project with the width in both (limit 1e-12), with a note; no note when file and `<levels>`
   agree.
3. **The channel matters** (the cross sections move by more than 1e-3 when it is switched off), so 1 and 2 prove
   something.

Before 2026-11-25 only `<levels>` decided the mask: a project whose `<levels>` held a photon width that the fit's
`param.sav` held as 0 kept that channel's external capture (checks 1 and 2 fail on such a build).  Found on a
12C(p,γ) project of the archive, whose capture angular distributions differed by 1.25e-3 from the engine before the
restructuring.
