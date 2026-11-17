# two_nuclei_ec — external capture in the second compound nucleus of a project

dev's `tests/two_nuclei_ec` (dev `6e918ed1`), brought here on 2026-11-15: the three models converted with
`AZURE2 --convert` (no note, no question), the script rewritten to the `extrapolate` command.

`models/A.azr` is 12C(p,γ0)13N and `models/B.azr` 12C(n,γ0)13C, each with an external-capture ground state;
`models/M.azr` holds both. The two nuclei do not couple, so each capture extrapolation of the merged project must
equal its standalone project (1e-10; measured: identical). dev died here with a segmentation fault until its repair;
this branch finds the capture exit by the pair's number (`tests/ec_exit_position`) and was right before this test.

Run three times: the default (Brune), `--no-brune` and `--no-ec-table` (dev's script ran the default twice, once
under the name `--use-brune`). Six checks, a few seconds.
