// transfer_n14_fresco (docs/upgrade/transfer.md sections 21 and 22): the 14N Trojan horse to a bound 12C + 12C
// cluster state against FRESCO 3.5, 12C(14N,d)24Mg at E_lab = 30 MeV, finite range (post, no remnant, the two 12C
// not symmetrised), to the 0+ ground state (l = 0, six nodes) and the 2+ at 1.369 MeV (l = 2, five nodes): cases
// n14d_<gs|2p>_finite_range.  Constants, kinematics, potentials, bound states (14N = 12C + d, 24Mg = 12C + 12C),
// dsigma/dOmega from 0 to 90 deg above 1 % of its maximum to 1 %, FRESCO's reference against its refined run.
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/transfer/Transfer.h"
#include "transfer_fresco_common.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;
using namespace fresco_reference;

namespace {

// The potentials of the engine's tests/reference/thm_distortion_test.cpp ("optical"), absolute radii.
OpticalPotential N14C12() {  // 14N + 12C (n14d_* kp 1)
  OpticalPotential p = Potential(RadiusConvention::Absolute, 1, 0, 5.5);
  p.real_volume = {50.0, 5.5, 0.6};
  p.imaginary_volume = {20.0, 5.6, 0.6};
  return p;
}
OpticalPotential D24Mg() {  // d + 24Mg (n14d_* kp 2)
  OpticalPotential p = Potential(RadiusConvention::Absolute, 1, 0, 3.8);
  p.real_volume = {90.0, 3.4, 0.75};
  p.imaginary_surface = {10.0, 3.8, 0.65};
  return p;
}
OpticalPotential C12D() {  // 12C-d binding of 14N (n14d_* kp 3)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 12, 0, 1.25);
  p.real_volume = {100.0, 1.25, 0.65};
  return p;
}
OpticalPotential C12C12() {  // 12C-12C binding of 24Mg (n14d_* kp 4)
  OpticalPotential p = Potential(RadiusConvention::BothMasses, 12, 12, 0.9);
  p.real_volume = {200.0, 0.9, 0.7};
  return p;
}

void TestCase(bool excited) {
  const std::string name = std::string("n14d_") + (excited ? "2p" : "gs") + "_finite_range";
  const auto runs = ReadCase(name);
  if (runs.empty()) return;
  const std::vector<Line> &ref = runs.at("reference");
  const double E_cm = 30.0 * 12.0 / 26.0, Q = 13.933 - 10.272 - (excited ? 1.3687 : 0.0);
  const ScatteringRequest in = Scattering(7, 6, 14.0, 12.0, 2, E_cm, 0.0, 39);
  const ScatteringRequest out = Scattering(1, 12, 2.0, 24.0, 2, E_cm, Q, 39);
  CheckConstants(name, ref);
  CheckKinematics(name, ref, "14N", in);
  CheckKinematics(name, ref, "d", out);
  CheckPotentials(name, ref, {{1, N14C12()}, {2, D24Mg()}, {3, C12D()}, {4, C12C12()}});
  TransferReaction r;
  r.mass_a_u = 14.0, r.mass_A_u = 12.0, r.mass_s_u = 2.0, r.mass_F_u = 24.0, r.two_s_x = 0;
  const BoundState n14 = SolveBoundState(C12D(), Bound(6, 1, 2.0, 12.0, 0, 0, 0, 10.272, 1));
  r.projectile.radial = RadialOverlapOf(n14);
  r.projectile.two_J_core = 2, r.projectile.two_J_composite = 2;
  // 2n + l = 12 oscillator quanta of 12C-12C relative motion (transfer.md section 21)
  const BoundState mg = excited ? SolveBoundState(C12C12(), Bound(6, 6, 12.0, 12.0, 2, 4, 0, 13.933 - 1.3687, 5))
                                : SolveBoundState(C12C12(), Bound(6, 6, 12.0, 12.0, 0, 0, 0, 13.933, 6));
  r.residual.radial = RadialOverlapOf(mg);
  r.residual.two_J_core = 0, r.residual.two_J_composite = excited ? 4 : 0;
  CheckBound(name, ref, 1, C12D(), n14, true);
  // FRESCO's u/W of the 12C-12C state still rises at 12 fm (g.s.: 5.17985, 5.25601, 5.26051 x 1e9 at 8, 10, 12 fm),
  // and ours moves by 8e-5 between h = 0.05 and 0.0125 fm (5.26038, 5.26081 x 1e9): 1e-4 for this ANC.
  CheckBound(name, ref, 2, C12C12(), mg, false, 1e-4);
  const DistortedWaves wi = SolveDistortedWaves(N14C12(), in), wf = SolveDistortedWaves(D24Mg(), out);
  TransferGrid g;
  g.step_fm = 0.1;
  g.r_max_fm = 25.0;
  g.vertex_range_fm = 10.0;
  g.angle_nodes = 48;
  CheckTransfer(name, FiniteRangeCrossSection(r, wi, wf, g, Angles()), 0.01);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s <tests/data/fresco>\n", argv[0]);
    return 2;
  }
  data_dir = argv[1];
  TestCase(false);
  TestCase(true);
  if (failures) {
    std::printf("transfer_n14_fresco: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_n14_fresco: all checks passed\n");
  return 0;
}
