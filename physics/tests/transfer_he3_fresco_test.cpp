// transfer_he3_fresco (docs/upgrade/transfer.md sections 20 and 22): the 3He Trojan horse to a bound state against
// FRESCO 3.5, 18O(3He,d)19F at E_lab = 25 MeV, to the ground state 1/2+ (2s1/2) and to the 5/2+ state at 0.197 MeV
// (1d5/2), each in finite range (post, no remnant) and zero range (D0 = -152.53 MeV fm^3/2, recoil), with and
// without spin-orbit terms in the 3He + 18O and d + 19F distortions: eight cases he3d_<state>_<range>[_so].  Per
// case: constants, kinematics, potentials, bound states (3He = d + p in a Woods-Saxon well, 19F = 18O + p), then
// dsigma/dOmega from 0 to 90 deg above 1 % of its maximum to 1 % and the 3He's iT11 to 1e-3 (absolute), with
// FRESCO's reference against its refined run (CheckTransfer, CheckAnalysingPower).
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/transfer/Transfer.h"
#include "transfer_fresco_common.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;
using namespace fresco_reference;

namespace {

OpticalPotential He3O18(bool spin_orbit) {  // 3He + 18O (he3d_* kp 1)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 18, 0, 1.3);
  p.real_volume = {155.0, 1.2, 0.72};
  p.imaginary_volume = {10.0, 1.4, 0.88};
  if (spin_orbit) p.real_spin_orbit = {2.5, 1.2, 0.72};
  return p;
}
OpticalPotential D19F(bool spin_orbit) {  // d + 19F (he3d_* kp 2)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 19, 0, 1.3);
  p.real_volume = {88.0, 1.17, 0.75};
  p.imaginary_volume = {2.0, 1.32, 0.65};
  p.imaginary_surface = {12.0, 1.32, 0.65};
  if (spin_orbit) p.real_spin_orbit = {6.0, 1.07, 0.66};
  return p;
}
OpticalPotential PD() {  // p-d binding of 3He (he3d_* kp 3)
  OpticalPotential p = Potential(RadiusConvention::Absolute, 1, 0, 1.5);
  p.real_volume = {60.0, 1.5, 0.5};
  return p;
}
OpticalPotential P18O() {  // p-18O binding (he3d_* kp 4)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 18, 0, 1.25);
  p.real_volume = {50.0, 1.25, 0.65};
  p.real_spin_orbit = {6.0, 1.25, 0.65};
  return p;
}

void TestCase(bool excited, bool spin_orbit) {
  const std::string tag = std::string(excited ? "d52" : "gs"), so = spin_orbit ? "_so" : "";
  const std::string fr_name = "he3d_" + tag + "_finite_range" + so, zr_name = "he3d_" + tag + "_zero_range" + so;
  const auto runs = ReadCase(fr_name);
  if (runs.empty()) return;
  const std::vector<Line> &ref = runs.at("reference");
  const double E_cm = 25.0 * 18.0 / 21.0, Q = 7.9937 - 5.4935 - (excited ? 0.1972 : 0.0);
  const ScatteringRequest in = Scattering(2, 8, 3.0, 18.0, 1, E_cm, 0.0, 36);
  const ScatteringRequest out = Scattering(1, 9, 2.0, 19.0, 2, E_cm, Q, 36);
  CheckConstants(fr_name, ref);
  CheckKinematics(fr_name, ref, "3He", in);
  CheckKinematics(fr_name, ref, "d", out);
  CheckPotentials(fr_name, ref, {{1, He3O18(spin_orbit)}, {2, D19F(spin_orbit)}, {3, PD()}, {4, P18O()}});
  TransferReaction r;
  r.mass_a_u = 3.0, r.mass_A_u = 18.0, r.mass_s_u = 2.0, r.mass_F_u = 19.0, r.two_s_x = 1;
  const BoundState he3 = SolveBoundState(PD(), Bound(1, 1, 2.0, 1.0, 0, 1, 1, 5.4935, 0));
  r.projectile.radial = RadialOverlapOf(he3);
  r.projectile.two_J_core = 2, r.projectile.two_J_composite = 1;
  const BoundState f19 = excited ? SolveBoundState(P18O(), Bound(1, 8, 18.0, 1.0, 2, 5, 1, 7.7965, 0))
                                 : SolveBoundState(P18O(), Bound(1, 8, 18.0, 1.0, 0, 1, 1, 7.9937, 1));
  r.residual.radial = RadialOverlapOf(f19);
  r.residual.two_J_core = 0, r.residual.two_J_composite = excited ? 5 : 1;
  // The p + d ANC: ours 1.631716 at h = 0.05, 0.025, 0.0125 fm; FRESCO's plateau 1.63168 at 8, 10, 12 fm and at
  // hcm = 0.025, 0.0125 fm: 2.2e-5 apart, not resolved (transfer.md section 20); it moves sigma by < 5e-5.
  CheckBound(fr_name, ref, 1, PD(), he3, true, 3e-5);
  CheckBound(fr_name, ref, 2, P18O(), f19, false);
  const DistortedWaves wi = SolveDistortedWaves(He3O18(spin_orbit), in),
                       wf = SolveDistortedWaves(D19F(spin_orbit), out);
  TransferGrid g;
  g.step_fm = 0.1;
  g.r_max_fm = 25.0;
  g.vertex_range_fm = 16.0;  // V_pd u_3He keeps its Coulomb tail: 10 fm leaves 2e-3 (transfer.md section 20)
  g.angle_nodes = 48;
  const TransferObservables fr = FiniteRangeObservables(r, wi, wf, g, Angles());
  CheckTransfer(fr_name, fr.cross_section, 0.01);
  CheckAnalysingPower(fr_name, fr.analysing_power_iT11, 1e-3);
  const TransferObservables zr = ZeroRangeObservables(r, -152.53, wi, wf, Angles());
  CheckTransfer(zr_name, zr.cross_section, 0.01);
  CheckAnalysingPower(zr_name, zr.analysing_power_iT11, 1e-3);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s <tests/data/fresco>\n", argv[0]);
    return 2;
  }
  data_dir = argv[1];
  for (bool excited : {false, true})
    for (bool spin_orbit : {false, true}) TestCase(excited, spin_orbit);
  if (failures) {
    std::printf("transfer_he3_fresco: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_he3_fresco: all checks passed\n");
  return 0;
}
