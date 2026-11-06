# excitation_energy_column

The excitation-energy column of the output files includes the pair's
excitation energy (audit 2026-10-25, 01_kinematics_frames B3; 2026-10-28).

`EPoint::ConvertLabEnergy`, `ConvertExcitationEnergy` and `ConvertDecayEnergy`
set the point's excitation energy to E_cm + S, while the physics uses
E_cm + S + E* of the pair (`CalcEDependentValues`). The E_x column of
`AZUREOut_*` and the API's excitation energies were therefore low by E* for a
segment whose entrance pair (exit pair, for beta-delayed data) is an
excited-state pair. Output only: no cross section or chi-squared depended on
it.

## Checks (`check.sh path/to/AZURE2`)

`tests/7Li_p_ay`, mode 3, two angle-integrated test segments:

- (a) pair 3 (n + 7Be*, S = 18.898 MeV, E* = 0.429 MeV): column 2 minus
  column 1 of `AZUREOut_aa=3_R=3.extrap` is 19.327 MeV at every energy;
- (b) pair 5 (p + 7Li, E* = 0): 17.25 MeV, unchanged.

dev b1ad0cd: (a) gives 18.898 (0.429 MeV low).
