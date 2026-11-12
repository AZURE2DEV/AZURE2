// transfer_continuum_transfer (docs/upgrade/transfer.md sections 15-17): the transfer to the x + A continuum and the
// forms of the interaction, against answers that do not need FRESCO.  15N(d,n) at E_lab = 10 MeV with the inputs of
// tests/data/fresco/dn_* (transfer_fresco_common.h).
//
//  1. The narrow-bin limit (section 15): d2sigma/dOmega dE to p + 15N at E = 2 MeV (l = 1, from the scattering state
//     in the delta(E - E') normalisation) against sigma_bin/dE of bins of width dE = 0.05 and 0.025 MeV around E
//     (Simpson's k rule, the support of both cut at 40 fm): 1e-3 at 0.025 MeV, and the error falls as dE^2 (ratio of
//     the two errors 4 +- 0.3).
//  2. The interior vertex and the level function (section 16): for p + 15N, l = 2, V = 53 MeV (one isolated
//     resonance) the vertex int_0^a w V u_E dr (w = r j_2(pr), delta(E - E') normalisation) against A_lambda(E) W_lambda
//     (A_lambda the projection of u_E on the level function, W_lambda the vertex of ThmVertex.h): equal at E_r
//     (1e-10), and the numbers of the table of section 16 at E_r -+ Gamma for a = 5, 6, 7 fm to 1e-3 relative.
//  3. Post = prior with both remnants (section 17), 15N(d,n)16O(g.s.): to 2e-3 at 0-90 deg (rho <= 40 fm).
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/coulomb/Coulomb.h"
#include "azure2/physics/transfer/RadialNumerics.h"
#include "azure2/physics/transfer/ThmVertex.h"
#include "azure2/physics/transfer/Transfer.h"
#include "transfer_fresco_common.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;
using namespace fresco_reference;

namespace {

const double kEcmDeuteron = 10.0 * 15.0 / 17.0;  // MeV

OpticalPotential P15NContinuum(double depth) {  // p + 15N, kp 5 of tests/data/fresco/dn_bin_*
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 15, 0, 1.25);
  p.real_volume = {depth, 1.25, 0.65};
  return p;
}

TransferReaction DeuteronToContinuum() {
  TransferReaction r;
  r.mass_a_u = 2.0, r.mass_A_u = 15.0, r.mass_s_u = 1.0, r.mass_F_u = 16.0, r.two_s_x = 1;
  BoundStateRequest bq = Bound(0, 1, 1.0, 1.0, 0, 1, 1, 2.2246, 0);
  r.projectile.radial = RadialOverlapOf(SolveBoundState(NP(), bq));
  r.projectile.two_J_core = 1, r.projectile.two_J_composite = 2;
  r.residual.two_J_core = 1, r.residual.two_J_composite = 2;
  return r;
}

// ---------------------------------------------------------------- 1. narrow-bin limit
void TestNarrowBinLimit() {
  TransferReaction r = DeuteronToContinuum();
  ContinuumRequest cq;
  cq.pair = {1, 7, 15.0 / 16.0};
  cq.l = 1, cq.two_j = 1, cq.two_s = 1, cq.step_fm = 0.05, cq.r_max_fm = 40.0;
  const double E = 2.0;
  ScatteringRequest in = Scattering(1, 7, 2.0, 15.0, 0, kEcmDeuteron, 0.0, 30);
  ScatteringRequest out = Scattering(0, 8, 1.0, 16.0, 0, kEcmDeuteron, -2.2246 - E, 30);
  in.matching_radius_fm = out.matching_radius_fm = 52.0;
  const DistortedWaves wi = SolveDistortedWaves(D15N(false), in), wf = SolveDistortedWaves(N16O(), out);
  const std::vector<double> th = {1e-3, 20.0, 40.0, 60.0, 90.0};
  TransferGrid g;
  g.r_max_fm = 50.0;
  r.residual.kind = OverlapKind::ScatteringState;
  r.residual.radial = RadialOverlapOf(SolveScatteringState(P15NContinuum(50.0), cq, E), cq.pair.mu_u);
  const std::vector<double> direct = FiniteRangeCrossSection(r, wi, wf, g, th);
  double err[2] = {0.0, 0.0};
  const double widths[2] = {0.05, 0.025};
  for (int k = 0; k < 2; ++k) {
    r.residual.kind = OverlapKind::Bin;
    r.residual.radial =
        RadialOverlapOf(MakeBin(P15NContinuum(50.0), cq, E, widths[k], 21, BinWeight::Flat, BinRule::Simpson));
    const std::vector<double> xs = FiniteRangeCrossSection(r, wi, wf, g, th);
    for (std::size_t i = 0; i < th.size(); ++i) err[k] = std::max(err[k], std::fabs(xs[i] / widths[k] / direct[i] - 1.0));
  }
  char buf[200];
  std::snprintf(buf, sizeof buf, "narrow-bin limit: |sigma_bin/dE / (d2sigma/dOmega dE) - 1| %.2e (dE = 0.05 MeV), "
                "%.2e (0.025 MeV), ratio %.2f; d2sigma/dOmega dE(0 deg) = %.4f mb/(sr MeV)", err[0], err[1],
                err[0] / err[1], direct[0]);
  Check(err[1] < 1e-3 && std::fabs(err[0] / err[1] - 4.0) < 0.3, buf);
  std::printf("%s\n", buf);
}

// ---------------------------------------------------------------- 2. interior vertex and level function
struct VertexRow {
  double a, rem_low, rem_high, ext_centre;  // the table of transfer.md section 16
};

void TestLevelFunction() {
  const OpticalPotential p = P15NContinuum(53.0);
  ContinuumRequest q;
  q.pair = {1, 7, 15.0 / 16.0};
  q.l = 2, q.two_j = 5, q.two_s = 1, q.step_fm = 0.01, q.r_max_fm = 40.0;
  const double mu = q.pair.mu_u * kAtomicMassUnit_MeV, c = 2.0 * mu / (kHbarC_MeV_fm * kHbarC_MeV_fm);
  auto delta = [&](double E) { return SolveScatteringState(p, q, E).phase_shift_rad; };
  double lo = 1.5, hi = 2.3;  // delta = pi/2: the principal value jumps from +pi/2 to -pi/2
  for (int it = 0; it < 60; ++it) (delta(0.5 * (lo + hi)) > 0.0 ? lo : hi) = 0.5 * (lo + hi);
  const double Er = 0.5 * (lo + hi), dh = 1e-4;
  const double up = delta(Er + dh), down = delta(Er - dh);  // across the jump of the principal value
  const double G = 2.0 / ((up - down + (up < down ? std::acos(-1.0) : 0.0)) / (2.0 * dh));
  char buf[300];
  std::snprintf(buf, sizeof buf, "level: E_r %.6f MeV (table 1.901482), Gamma %.6f MeV (0.095190)", Er, G);
  Check(std::fabs(Er - 1.901482) < 2e-6 && std::fabs(G - 0.095190) < 2e-6, buf);
  std::printf("%s\n", buf);
  const ScatteringState sr = SolveScatteringState(p, q, Er);
  const VertexRow table[3] = {{5.0, -1.262e-03, 1.281e-03, 8.756e-02},
                              {6.0, -6.466e-03, 6.440e-03, 2.059e-02},
                              {7.0, -1.873e-02, 1.815e-02, 4.708e-03}};
  for (const VertexRow &row : table) {
    const int na = static_cast<int>(std::lround(row.a / q.step_fm));
    LevelFunction lv;
    lv.step_fm = q.step_fm, lv.energy_MeV = Er, lv.mu_MeV = mu;
    std::vector<double> sq(na + 1);
    for (int i = 0; i <= na; ++i) {
      lv.u.push_back(sr.u[i]);
      lv.potential_MeV.push_back(sr.potential_MeV[i]);
      sq[i] = sr.u[i] * sr.u[i];
    }
    const double norm = std::sqrt(Simpson(sq, q.step_fm));
    for (double &x : lv.u) x /= norm;
    const double B = BoundaryValue(lv);
    double rem[3], ext = 0.0;
    const double energies[3] = {Er - G, Er, Er + G};
    for (int e = 0; e < 3; ++e) {
      const ScatteringState s = SolveScatteringState(p, q, energies[e]);
      const double fE = NormalisationFactor(s, q.pair.mu_u, ContinuumNormalisation::DeltaE);
      const double eps = energies[e] + 2.2246, pm = std::sqrt(c * eps);
      std::vector<std::complex<double>> w(na + 1);
      std::vector<double> fin(na + 1), proj(na + 1), fout(s.u.size() - na);
      for (int i = 0; i <= na; ++i) {
        const double r = i * q.step_fm;
        const double F = r > 0.0 ? coulomb::coulomb(2, 0.0, pm * r, coulomb::Method::coulfg).F / pm : 0.0;
        w[i] = F;
        fin[i] = F * s.potential_MeV[i] * s.u[i] * fE;
        proj[i] = lv.u[i] * s.u[i] * fE;
      }
      for (std::size_t i = na; i < s.u.size(); ++i) {
        const double r = i * q.step_fm, F = coulomb::coulomb(2, 0.0, pm * r, coulomb::Method::coulfg).F / pm;
        fout[i - na] = F * (s.potential_MeV[i] - s.coulomb_MeV[i]) * s.u[i] * fE;
      }
      const double dwa = coulomb::coulomb(2, 0.0, pm * row.a, coulomb::Method::coulfg).Fp;
      const double Wl = EntranceVertex(lv, w, dwa, eps, B).volume.real();
      const double Wint = Simpson(fin, q.step_fm), A = Simpson(proj, q.step_fm);
      rem[e] = (Wint - A * Wl) / Wint;
      if (e == 1) ext = Simpson(fout, q.step_fm) / Wint;
    }
    std::snprintf(buf, sizeof buf, "level, a = %.0f fm: remainder %.3e, %.1e, %.3e (E_r -+ Gamma), exterior/interior "
                  "%.3e at E_r", row.a, rem[0], rem[1], rem[2], ext);
    Check(std::fabs(rem[1]) < 1e-10 && std::fabs(rem[0] / row.rem_low - 1.0) < 1e-3 &&
              std::fabs(rem[2] / row.rem_high - 1.0) < 1e-3 && std::fabs(ext / row.ext_centre - 1.0) < 1e-3, buf);
    std::printf("%s\n", buf);
  }
}

// ---------------------------------------------------------------- 3. post = prior with the remnants
void TestPostPrior() {
  TransferReaction r;
  r.mass_a_u = 2.0, r.mass_A_u = 15.0, r.mass_s_u = 1.0, r.mass_F_u = 16.0, r.two_s_x = 1;
  BoundStateRequest bq = Bound(0, 1, 1.0, 1.0, 0, 1, 1, 2.2246, 0);
  bq.r_max_fm = 50.0;  // the deuteron tail: the prior-form kernel needs u_a, not V u_a
  r.projectile.radial = RadialOverlapOf(SolveBoundState(NP(), bq));
  r.projectile.two_J_core = 1, r.projectile.two_J_composite = 2;
  r.residual.radial = RadialOverlapOf(SolveBoundState(P15NBound(), Bound(1, 7, 15.0, 1.0, 1, 1, 1, 12.1274, 0)));
  r.residual.two_J_core = 1, r.residual.two_J_composite = 0;
  ScatteringRequest in = Scattering(1, 7, 2.0, 15.0, 0, kEcmDeuteron, 0.0, 30);
  ScatteringRequest out = Scattering(0, 8, 1.0, 16.0, 0, kEcmDeuteron, 9.9028, 30);
  in.matching_radius_fm = out.matching_radius_fm = 50.0;
  const DistortedWaves wi = SolveDistortedWaves(D15N(false), in), wf = SolveDistortedWaves(N16O(), out);
  std::vector<double> th;
  for (int i = 0; i <= 18; ++i) th.push_back(i == 0 ? 1e-3 : 5.0 * i);
  TransferGrid g;
  g.r_max_fm = 45.0;
  g.vertex_range_fm = 40.0;
  r.interaction = PostFormWithRemnant(N15Core(), 0, N16O(), 0);
  const std::vector<double> post = FiniteRangeCrossSection(r, wi, wf, g, th);
  r.interaction = PriorFormWithRemnant(N15Core(), 0, D15N(false), 7);
  const std::vector<double> prior = FiniteRangeCrossSection(r, wi, wf, g, th);
  double worst = 0.0, at = 0.0, peak = *std::max_element(post.begin(), post.end());
  for (std::size_t i = 0; i < th.size(); ++i)
    if (post[i] > 0.01 * peak && std::fabs(prior[i] / post[i] - 1.0) > worst) worst = std::fabs(prior[i] / post[i] - 1.0), at = th[i];
  char buf[200];
  std::snprintf(buf, sizeof buf, "post = prior with remnants: worst %.2e at %.0f deg (0 deg: %.5f, %.5f mb/sr)", worst,
                at, post[0], prior[0]);
  Check(worst < 2e-3, buf);
  std::printf("%s\n", buf);
}

}  // namespace

int main() {
  TestNarrowBinLimit();
  TestLevelFunction();
  TestPostPrior();
  if (failures) {
    std::printf("transfer_continuum_transfer: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_continuum_transfer: all checks passed\n");
  return 0;
}
