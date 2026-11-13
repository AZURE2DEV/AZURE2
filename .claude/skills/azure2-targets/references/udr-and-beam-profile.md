# Beam-profile and UDR kernels: how to set them up

The equations of both kernels are in `docs/source/user_guide/experimental_effects.rst` ("Beam-profile kernel" and
"User-defined numerical resolution function"). This page holds the set-up recipe and the checks.

## Beam profile (broad γ beam, event-by-event energy reconstruction)

The point-centred Gaussian is the wrong kernel for HIγS TPC data: the beam profile G(E) is an absolute
distribution and each datum is an energy slice [a, b] of it seen through the detector resolution s.

1. One data set per nominal beam energy, each naming its own `beam` record (the profile) and the `detector` record of
   the TPC resolution (`Project.add_beam`, `add_detector`, `set`; fields in `docs/source/reference/format3.rst`).
2. The per-point window [a, b] goes in two columns of the data file with the roles `window_lo[MeV]` and
   `window_hi[MeV]` (lab energies of the entrance channel; a photon-beam profile is shifted by the Q value and
   converted first). The detailed-balance weight of the profile averages a capture σ the way the inverse
   measurement did.
3. Give the data set a resonance-width multiplier of at least 20 when a resonance narrower than the base sub-point
   step lies under the beam (10 widths biased a slice dominated by a 0.6 keV 2⁺ by 18 %; 20 agreed with a dense
   reference to 4 × 10⁻⁴).
6. Check: `tests/beam_profile_kernel/` compares the engine with an independent numpy kernel.

## UDR (tabulated TOF resolution, SAMMY format)

The `udr` field of a `detector` record: `file <path> flight_path <L> m burst <B> ns channel <C> ns centred yes|no`.
A negative channel value −n means n bins per decade (n_TOF: 10 below 300 keV, 100 above). Only the blocks within a
decade of the data set's energies are read, so the 150 MB n_TOF EAR1 file loads in seconds.

Checks before a real fit:

- **Right file for the station.** An EAR1 file does not describe EAR2 data. EAR2 files are made per experiment
  (flight path, sample size, alignment) with the Transport Code and `RF2sammy`; obtain the file the authors' SAMMY
  fit used.
- **Burst already folded in** (n_TOF files): `burst 0 ns`. The delay rms of the highest-energy blocks (6–7 ns in the
  EAR1 HI file) is the burst the file contains; HI (7 ns) and LI (14 ns) pulses have separate files.
- **Centring.** `centred yes` re-centres each interpolated kernel on its centroid (SAMMY's default). If the data's
  TOF → E conversion used a fitted flight path L₀ that already absorbs the mean moderation delay, use `centred no`
  with that L₀.
- **Nominal energy.** The kernel is evaluated at the point's nominal energy (as SAMMY does), not at each
  sub-point's true energy.
- **One channel width per line.** Split a data set where the binning changes if both parts matter.
- **Validate on the cluster, not the login node**: a synthetic Gaussian-in-time UDR against a numpy fold of the same
  data set agreed to 7.6 × 10⁻⁴ relative; repeat that check for a new file. `tests/reference/udr_reference_test.cpp`
  pins moments, centring, widths in quadrature, interpolation and kinematics.
