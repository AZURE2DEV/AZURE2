// transfer_spin_orbit (docs/upgrade/transfer.md sections 19-20): the transfer amplitude with spin-orbit terms in the
// distorted waves, for 15N(d,n)16O(g.s.) at E_lab = 10 MeV (the reaction of part 1's dn_* cases):
//
//  1. Reduction: waves flagged as having a spin-orbit term of strength 1e-300 MeV (every (l, j) solved, all equal)
//     give the cross section and iT11 of waves without one (one per l) to 1e-12, finite and zero range: the sums
//     over j of Eq. (19.5) are then the identity (19.7), and the kernel's pairs of waves (Eq. (6.4)) are indexed
//     right.
//  2. Reflection in the reaction plane (Eq. (19.9)), with spin-orbit in both channels: T(-mu, -M_F; -m, -M_A) =
//     (-1)^(l_a + l_F) (-1)^(sum of J - M) T(mu, M_F; m, M_A) to 1e-10 of the largest |T|, and Tr(T S_x T+) = 0;
//     neither is an input of the code (it never uses the symmetry).
//  3. FRESCO (dn_finite_range_so, dn_zero_range_so): potentials with their spin-orbit terms, bound states,
//     dsigma/dOmega from 0 to 90 deg above 1 % of its maximum to 1 %, the deuteron's iT11 to 1e-3 (absolute), and
//     FRESCO's reference against its refined run (0.2 %, 1e-3).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/transfer/Transfer.h"
#include "azure2/physics/transfer/TransferAmplitudes.h"
#include "azure2/physics/transfer/TransferKernel.h"
#include "transfer_fresco_common.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;
using namespace fresco_reference;

namespace {

OpticalPotential N16OSpinOrbit() {  // n + 16O with spin-orbit (dn_*_so kp 2)
  OpticalPotential p = N16O();
  p.real_spin_orbit = {5.5, 1.2, 0.65};
  return p;
}

TransferReaction DN(BoundState &np, BoundState &p15n) {
  TransferReaction r;
  r.mass_a_u = 2.0, r.mass_A_u = 15.0, r.mass_s_u = 1.0, r.mass_F_u = 16.0, r.two_s_x = 1;
  np = SolveBoundState(NP(), Bound(0, 1, 1.0, 1.0, 0, 1, 1, 2.2246, 0));
  p15n = SolveBoundState(P15NBound(), Bound(1, 7, 15.0, 1.0, 1, 1, 1, 12.1274, 0));
  r.projectile.radial = RadialOverlapOf(np);
  r.projectile.two_J_core = 1, r.projectile.two_J_composite = 2;
  r.residual.radial = RadialOverlapOf(p15n);
  r.residual.two_J_core = 1, r.residual.two_J_composite = 0;
  return r;
}

ScatteringRequest Entrance(int two_s) { return Scattering(1, 7, 2.0, 15.0, two_s, 10.0 * 15.0 / 17.0, 0.0, 30); }
ScatteringRequest Exit(int two_s) { return Scattering(0, 8, 1.0, 16.0, two_s, 10.0 * 15.0 / 17.0, 9.9028, 30); }

TransferGrid Grid() {
  TransferGrid g;
  g.step_fm = 0.1;
  g.r_max_fm = 25.0;
  g.vertex_range_fm = 10.0;
  g.angle_nodes = 48;
  return g;
}

// ---------------------------------------------------------------- 1. reduction
void TestReduction() {
  BoundState np, p15n;
  const TransferReaction r = DN(np, p15n);
  OpticalPotential d_flagged = D15N(false), n_flagged = N16O();
  d_flagged.real_spin_orbit = {1e-300, 1.07, 0.66};  // a spin-orbit term that changes no wave: every j solved
  n_flagged.real_spin_orbit = {1e-300, 1.2, 0.65};
  const DistortedWaves wi = SolveDistortedWaves(D15N(false), Entrance(2)), wf = SolveDistortedWaves(N16O(), Exit(1));
  const DistortedWaves wi_j = SolveDistortedWaves(d_flagged, Entrance(2)), wf_j = SolveDistortedWaves(n_flagged, Exit(1));
  Check(wi_j.spin_orbit && wf_j.spin_orbit && wi_j.waves.size() > wi.waves.size(), "reduction: waves per j solved");
  const std::vector<double> th = Angles();
  const TransferObservables fr = FiniteRangeObservables(r, wi, wf, Grid(), th),
                            fr_j = FiniteRangeObservables(r, wi_j, wf_j, Grid(), th);
  const TransferObservables zr = ZeroRangeObservables(r, -130.31, wi, wf, th),
                            zr_j = ZeroRangeObservables(r, -130.31, wi_j, wf_j, th);
  double worst = 0.0, worst_a = 0.0;
  for (std::size_t i = 0; i < th.size(); ++i) {
    worst = std::max({worst, std::fabs(fr_j.cross_section[i] / fr.cross_section[i] - 1.0),
                      std::fabs(zr_j.cross_section[i] / zr.cross_section[i] - 1.0)});
    worst_a = std::max({worst_a, std::fabs(fr_j.analysing_power_iT11[i] - fr.analysing_power_iT11[i]),
                        std::fabs(zr_j.analysing_power_iT11[i] - zr.analysing_power_iT11[i])});
  }
  char buf[200];
  std::snprintf(buf, sizeof buf, "reduction without spin-orbit: cross sections %.1e (relative), iT11 %.1e", worst,
                worst_a);
  Check(worst < 1e-12 && worst_a < 1e-12, buf);
  std::printf("%s\n", buf);
}

// ---------------------------------------------------------------- 2. reflection in the reaction plane
double Sx(int two_s, int two_m, int two_mt) {  // <m | S_x | m~>
  const double s = 0.5 * two_s, mt = 0.5 * two_mt;
  if (two_m == two_mt + 2) return 0.5 * std::sqrt(s * (s + 1.0) - mt * (mt + 1.0));
  if (two_m == two_mt - 2) return 0.5 * std::sqrt(s * (s + 1.0) - mt * (mt - 1.0));
  return 0.0;
}

void TestReflection() {
  BoundState np, p15n;
  const TransferReaction r = DN(np, p15n);
  const DistortedWaves wi = SolveDistortedWaves(D15N(true), Entrance(2)),
                       wf = SolveDistortedWaves(N16OSpinOrbit(), Exit(1));
  const double alpha = 15.0 / 16.0, beta = 0.5, jac = 1.0 - alpha * beta;
  const KernelIntegrals K = FiniteRangeIntegrals(r, wi, wf, Grid());
  const AngularCoefficients c =
      AngularCoefficientsOf(r, wi, wf, K, 16.0 * M_PI * M_PI / (wi.k_fm * wf.k_fm * std::pow(jac, 3)));
  const int eta = (r.projectile.radial.l + r.residual.radial.l) % 2 ? -1 : 1;
  double worst = 0.0, worst_x = 0.0;
  for (double th : {10.0, 35.0, 80.0}) {
    const AmplitudeTable t = AmplitudesAt(r, c, th);
    double big = 0.0, sum = 0.0;
    for (const auto &v : t.T) big = std::max(big, std::abs(v)), sum += std::norm(v);
    std::complex<double> trace_x = 0.0;
    for (int mu = -t.two_Js; mu <= t.two_Js; mu += 2)
      for (int MF = -t.two_JF; MF <= t.two_JF; MF += 2)
        for (int m = -t.two_Ja; m <= t.two_Ja; m += 2)
          for (int MA = -t.two_JA; MA <= t.two_JA; MA += 2) {
            const int phase = ((t.two_Js - mu) + (t.two_JF - MF) + (t.two_Ja - m) + (t.two_JA - MA)) / 2;
            const double sign = eta * ((phase % 2) ? -1.0 : 1.0);
            const std::complex<double> a = t.T[t.Index(mu, MF, m, MA)];
            worst = std::max(worst, std::abs(a - sign * t.T[t.Index(-mu, -MF, -m, -MA)]) / big);
            for (int mt = -t.two_Ja; mt <= t.two_Ja; mt += 2)
              trace_x += a * Sx(t.two_Ja, m, mt) * std::conj(t.T[t.Index(mu, MF, mt, MA)]);
          }
    worst_x = std::max(worst_x, std::abs(trace_x) / sum);
  }
  char buf[200];
  std::snprintf(buf, sizeof buf, "reflection in the reaction plane: %.1e of max |T|; Tr(T S_x T+)/Tr(T T+) %.1e", worst,
                worst_x);
  Check(worst < 1e-10 && worst_x < 1e-10, buf);
  std::printf("%s\n", buf);
}

// ---------------------------------------------------------------- 3. FRESCO
void TestFresco() {
  const auto runs = ReadCase("dn_finite_range_so");
  if (runs.empty()) return;
  const std::vector<Line> &ref = runs.at("reference");
  CheckConstants("dn_finite_range_so", ref);
  CheckKinematics("dn_finite_range_so", ref, "d", Entrance(2));
  CheckKinematics("dn_finite_range_so", ref, "n", Exit(1));
  CheckPotentials("dn_finite_range_so", ref, {{1, D15N(true)}, {2, N16OSpinOrbit()}, {3, NP()}, {4, P15NBound()}});
  BoundState np, p15n;
  const TransferReaction r = DN(np, p15n);
  CheckBound("dn_finite_range_so", ref, 1, NP(), np, true);
  CheckBound("dn_finite_range_so", ref, 2, P15NBound(), p15n, false);
  const DistortedWaves wi = SolveDistortedWaves(D15N(true), Entrance(2)),
                       wf = SolveDistortedWaves(N16OSpinOrbit(), Exit(1));
  const TransferObservables fr = FiniteRangeObservables(r, wi, wf, Grid(), Angles());
  CheckTransfer("dn_finite_range_so", fr.cross_section, 0.01);
  CheckAnalysingPower("dn_finite_range_so", fr.analysing_power_iT11, 1e-3);
  const TransferObservables zr = ZeroRangeObservables(r, -130.31, wi, wf, Angles());
  CheckTransfer("dn_zero_range_so", zr.cross_section, 0.01);
  CheckAnalysingPower("dn_zero_range_so", zr.analysing_power_iT11, 1e-3);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s <tests/data/fresco>\n", argv[0]);
    return 2;
  }
  data_dir = argv[1];
  TestReduction();
  TestReflection();
  TestFresco();
  if (failures) {
    std::printf("transfer_spin_orbit: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_spin_orbit: all checks passed\n");
  return 0;
}
