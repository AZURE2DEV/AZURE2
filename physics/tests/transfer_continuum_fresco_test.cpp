// transfer_continuum_fresco (docs/upgrade/transfer.md section 15): the bins and the transfer to them against FRESCO
// 3.5, from tests/data/fresco/dn_bin_*_fr (make_fresco_benchmarks.py; FRESCO is not run here).  15N(d,n)[p + 15N] at
// E_lab = 10 MeV, finite range, post, no remnant, to three bins of p + 15N (kp 5, no spin-orbit), each on [0, RSP]
// with RSP = 40 fm and FRESCO's k rule:
//   nonresonant: l = 1, V = 50 MeV, 1.055-2.945 MeV, weight 1 (ISC = 2);
//   narrow:      l = 2, V = 53 MeV, 1.445-2.355 MeV around a resonance at 1.9015 MeV, Gamma = 95 keV (ISC = 4);
//   wide:        l = 0, V = 47 MeV, 1.0-2.16 MeV inside a resonance of Gamma ~ 3.3 MeV (ISC = 4).
//  1. The bin (FRESCO's listing): energies (6 digits), the phase shift at every 20th node (2e-3 deg: FRESCO prints 3
//     decimals and matches at RSP and RSP - h, we at RSP and RSP - 1 fm), Norm, RMS (4 digits), D0 (l = 0; 2 decimals) and u(r) every 50th point but its last (4 digits).
//  2. dsigma/dOmega from 0 to 90 deg where it is above 1 % of its maximum: ours against FRESCO's refined run (steps
//     halved) to 1 % at NK = 1600 and at NK = 3200 (the same bins: FRESCO's k rule converges as 1/NK, and so does
//     ours with it); FRESCO's rmatch (50 against 60 fm) and rnl (12 against 10 fm) variants against its reference
//     to 2e-3, and its step refinement printed (it is above 2e-3 at two minima: section 15).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/transfer/Transfer.h"
#include "transfer_fresco_common.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;
using namespace fresco_reference;

namespace {

struct BinCase {
  const char *name;
  int l, two_j, two_JF;
  double depth, centre, width;
  BinWeight weight;
};

std::vector<double> Xs(const std::vector<Line> &run) {
  std::vector<double> xs;
  for (const Line &l : Select(run, "xs")) xs.push_back(D(l[2]));
  return xs;
}

/// The worst relative difference of a against b at 0-90 deg where b is above 1 % of its maximum there.
double Worst(const std::vector<double> &th, const std::vector<double> &a, const std::vector<double> &b, double &at) {
  double peak = 0.0, worst = 0.0;
  for (std::size_t i = 0; i < th.size(); ++i)
    if (th[i] <= 90.0) peak = std::max(peak, b[i]);
  for (std::size_t i = 0; i < th.size() && i < a.size() && i < b.size(); ++i)
    if (th[i] <= 90.0 && b[i] > 0.01 * peak && std::fabs(a[i] / b[i] - 1.0) > worst)
      worst = std::fabs(a[i] / b[i] - 1.0), at = th[i];
  return worst;
}

void TestBin(const BinCase &c) {
  const std::string name = std::string("dn_bin_") + c.name + "_fr";
  const auto runs = ReadCase(name);
  if (runs.empty()) return;
  const std::vector<Line> &ref = runs.at("reference");
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 15, 0, 1.25);
  p.real_volume = {c.depth, 1.25, 0.65};
  ContinuumRequest q;
  q.pair = {1, 7, 15.0 / 16.0};
  q.l = c.l, q.two_j = c.two_j, q.two_s = 1, q.step_fm = 0.025, q.r_max_fm = 40.0;  // FRESCO's hcm and RSP
  const ContinuumBin bin = MakeBin(p, q, c.centre, c.width, 1600, c.weight, BinRule::Fresco);
  // 1. the listing
  const Line e = Select(ref, "bin").at(0);
  double worst_phase = 0.0, worst_wf = 0.0, umax = 0.0;
  for (const Line &l : Select(ref, "binphase"))
    worst_phase = std::max(worst_phase, std::fabs(bin.phase_shift_rad[std::stoi(l[1]) - 1] * 180.0 / std::acos(-1.0) - D(l[4])));
  for (double u : bin.u) umax = std::max(umax, std::fabs(u));
  for (const Line &l : Select(ref, "binwf"))  // FRESCO leaves its last node unnormalised
    if (std::stoi(l[1]) < static_cast<int>(bin.u.size()))
      worst_wf = std::max(worst_wf, std::fabs(bin.u[std::stoi(l[1]) - 1] - D(l[2])));
  const BinListing list = ListBin(bin);
  const Line n = Select(ref, "binnorm").at(0);
  char buf[400];
  std::snprintf(buf, sizeof buf,
                "%s: bin %.6f-%.6f MeV (FRESCO %s-%s), phases worst %.1e deg, Norm %.4f (%s), RMS %.4f (%s), D0 %.2f "
                "(%s), u worst %.1e of max",
                name.c_str(), bin.energy_min_MeV, bin.energy_max_MeV, e[1].c_str(), e[2].c_str(), worst_phase, list.norm,
                n[1].c_str(), list.rms_fm, n[2].c_str(), list.d0_MeV_fm32, n[3].c_str(), worst_wf / umax);
  bool ok = std::fabs(bin.energy_min_MeV - D(e[1])) < 6e-7 && std::fabs(bin.energy_max_MeV - D(e[2])) < 6e-7 &&
            worst_phase < 2e-3 && std::fabs(list.norm - D(n[1])) < 6e-5 && std::fabs(list.rms_fm - D(n[2])) < 6e-4 &&
            worst_wf < 1e-3 * umax && std::stoi(e[3]) == 1600;
  if (c.l == 0) ok = ok && std::fabs(list.d0_MeV_fm32 - D(n[3])) < 0.051;
  Check(ok, buf);
  std::printf("%s\n", buf);

  // 2. the transfer
  TransferReaction r;
  r.mass_a_u = 2.0, r.mass_A_u = 15.0, r.mass_s_u = 1.0, r.mass_F_u = 16.0, r.two_s_x = 1;
  r.projectile.radial = RadialOverlapOf(SolveBoundState(NP(), Bound(0, 1, 1.0, 1.0, 0, 1, 1, 2.2246, 0)));
  r.projectile.two_J_core = 1, r.projectile.two_J_composite = 2;
  r.residual.kind = OverlapKind::Bin;
  r.residual.radial = RadialOverlapOf(bin);
  r.residual.two_J_core = 1, r.residual.two_J_composite = c.two_JF;
  const double Ecm = 10.0 * 15.0 / 17.0;
  ScatteringRequest in = Scattering(1, 7, 2.0, 15.0, 0, Ecm, 0.0, 30);
  ScatteringRequest out = Scattering(0, 8, 1.0, 16.0, 0, Ecm, -2.2246 - c.centre, 30);
  in.matching_radius_fm = out.matching_radius_fm = 50.0;
  CheckKinematics(name, ref, "d", in);
  CheckKinematics(name, ref, "n", out);
  const DistortedWaves wi = SolveDistortedWaves(D15N(false), in), wf = SolveDistortedWaves(N16O(), out);
  TransferGrid g;
  g.r_max_fm = 45.0;  // the bin is 0 beyond 40 fm, V_np u_d beyond 10 fm
  std::vector<double> th;
  for (const Line &l : Select(ref, "xs")) th.push_back(std::max(D(l[1]), 1e-3));
  const std::vector<double> ours = FiniteRangeCrossSection(r, wi, wf, g, th);
  r.residual.radial = RadialOverlapOf(MakeBin(p, q, c.centre, c.width, 3200, c.weight, BinRule::Fresco));
  const std::vector<double> ours_nk = FiniteRangeCrossSection(r, wi, wf, g, th);
  const std::vector<double> x_ref = Xs(ref), x_fine = Xs(runs.at("refined")), x_rmatch = Xs(runs.at("variant1")),
                            x_rnl = Xs(runs.at("variant2")), x_nk = Xs(runs.at("variant3"));
  double at1 = 0, at2 = 0, at3 = 0, at4 = 0, at5 = 0;
  const double d_fine = Worst(th, ours, x_fine, at1), d_nk = Worst(th, ours_nk, x_nk, at2);
  const double f_step = Worst(th, x_fine, x_ref, at3), f_rmatch = Worst(th, x_rmatch, x_ref, at4),
               f_rnl = Worst(th, x_rnl, x_ref, at5);
  std::snprintf(buf, sizeof buf,
                "%s: ours - FRESCO %.2e at %.0f deg (NK 1600), %.2e at %.0f deg (NK 3200); FRESCO's own: steps %.1e "
                "(%.0f deg), rmatch %.1e, rnl %.1e; 0 deg: %.4f mb/sr (FRESCO %.4f)",
                name.c_str(), d_fine, at1, d_nk, at2, f_step, at3, f_rmatch, f_rnl, ours[0], x_fine[0]);
  Check(d_fine < 0.01 && d_nk < 0.01 && f_rmatch < 2e-3 && f_rnl < 2e-3, buf);
  std::printf("%s\n", buf);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s <tests/data/fresco>\n", argv[0]);
    return 2;
  }
  data_dir = argv[1];
  const BinCase cases[3] = {{"nonresonant", 1, 1, 2, 50.0, 2.0, 1.89, BinWeight::Flat},
                            {"narrow", 2, 5, 4, 53.0, 1.9, 0.91, BinWeight::TMatrix},
                            {"wide", 0, 1, 2, 47.0, 1.58, 1.16, BinWeight::TMatrix}};
  for (const BinCase &c : cases) TestBin(c);
  if (failures) {
    std::printf("transfer_continuum_fresco: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_continuum_fresco: all checks passed\n");
  return 0;
}
