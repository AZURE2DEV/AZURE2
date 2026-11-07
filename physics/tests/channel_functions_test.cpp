// channel_functions (docs/upgrade/coulomb.md section 5): P, S, phi and S' of an open channel, and eta, rho.
//
//  1. Neutral closed forms (eta = 0): l = 0, P = rho, S = 0, phi = rho; l = 1, P = rho^3/(1 + rho^2),
//     S = -1/(1 + rho^2), phi = rho - atan(rho); S' of l = 1 from d rho^2/dE = rho^2/E.  An error in a formula of
//     Eq. (4) (a wrong power, sign, or F and G swapped) fails here at the 1e-2 level or more.
//  2. Charged channels: P, S, exp(-i phi) from the library's Coulomb functions against the same formulas evaluated
//     on mpmath's F, F', G, G' (data/coulomb.csv), for eta > 0: to 1e-8 (coulsm is good to 1e-9).
//  3. S' (Richardson, h = 1e-3 E) against a five-point difference of S computed with coulcc (another algorithm,
//     S good to 1e-12) at h = 1e-4 E: to 1e-6 relative.
//  4. Kinematics: non-relativistic eta = Z1 Z2 alpha / beta, rho = k a from E = (hbar k)^2 / 2 mu, written
//     independently; relativistic k^2 = (s - (m1 + m2)^2)(s - (m1 - m2)^2) / 4s with s = (m1 + m2 + T)^2, and the
//     non-relativistic limit T/m -> 0.
//  5. E <= 0: invalid_argument and NaN.
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/Constants.h"

using namespace azure2::physics;

namespace {

int failures = 0;

void Check(bool ok, const char *what, double got, double bound) {
  if (!ok) {
    std::printf("FAIL: %s: %.3e (bound %.1e)\n", what, got, bound);
    ++failures;
  }
}

struct Row {
  int l;
  double eta, rho, F, Fp, G, Gp;
};

std::vector<Row> Read(const std::string &path) {
  std::vector<Row> rows;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    Row r{};
    std::string t[6];
    ss >> r.l >> t[0] >> t[1] >> t[2] >> t[3] >> t[4] >> t[5];
    double *v[6] = {&r.eta, &r.rho, &r.F, &r.Fp, &r.G, &r.Gp};
    for (int i = 0; i < 6; ++i) *v[i] = std::strtod(t[i].c_str(), nullptr);
    rows.push_back(r);
  }
  return rows;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s coulomb.csv\n", argv[0]);
    return 2;
  }

  // 1. Neutral closed forms: n + 12C, a = 4 fm, E_c = 0.01-10 MeV (rho = 0.09-2.8).
  const ChannelKinematics n12c{0, 6, 1.008664916 * 12.0 / (1.008664916 + 12.0), 0.0, 0.0, false};
  double e1 = 0.0;
  for (double E : {0.01, 0.1, 0.5, 1.0, 3.0, 10.0}) {
    const double r = rho(n12c, 4.0, E);
    for (coulomb::Method m : {coulomb::Method::coulsm, coulomb::Method::coulcc}) {
      const OpenChannelValues v0 = channel_functions(n12c, 0, 4.0, E, m), v1 = channel_functions(n12c, 1, 4.0, E, m);
      const double phi0 = -std::arg(v0.exp_minus_i_phi), phi1 = -std::arg(v1.exp_minus_i_phi);
      const double d[] = {std::fabs(v0.P / r - 1.0), std::fabs(v0.S), std::fabs(phi0 - r), std::fabs(v0.dS_dE),
                          std::fabs(v1.P / (r * r * r / (1.0 + r * r)) - 1.0), std::fabs(v1.S + 1.0 / (1.0 + r * r)),
                          std::fabs(phi1 - (r - std::atan(r))),
                          std::fabs(v1.dS_dE / ((r * r / E) / ((1.0 + r * r) * (1.0 + r * r))) - 1.0)};
      for (double x : d) e1 = std::fmax(e1, std::isnan(x) ? INFINITY : x);
    }
  }
  std::printf("neutral closed forms: %.2e\n", e1);
  Check(e1 <= 1e-7, "neutral closed forms", e1, 1e-7);  // S' carries the 1e-8 of the difference

  // 2. Charged channels against mpmath's Coulomb functions.
  const std::vector<Row> rows = Read(argv[1]);
  if (rows.size() != 1296) {
    std::printf("FAIL: %s incomplete\n", argv[1]);
    return 1;
  }
  double eP = 0.0, eS = 0.0, ePhi = 0.0;
  for (const Row &r : rows) {
    if (r.eta == 0.0) continue;
    const coulomb::CoulombWaves ref{r.F, r.Fp, r.G, r.Gp, 0.0, coulomb::Status::ok};
    const coulomb::CoulombWaves w = coulomb::coulomb(r.l, r.eta, r.rho);
    const double Pref = penetrability(ref, r.rho), Sref = shift(ref, r.rho);
    eP = std::fmax(eP, std::fabs(penetrability(w, r.rho) / Pref - 1.0));
    eS = std::fmax(eS, std::fabs(shift(w, r.rho) - Sref) / (1.0 + std::fabs(Sref)));
    ePhi = std::fmax(ePhi, std::abs(exp_minus_i_phi(w) - exp_minus_i_phi(ref)));
  }
  std::printf("against mpmath (eta > 0): P %.2e relative, S %.2e, exp(-i phi) %.2e\n", eP, eS, ePhi);
  Check(eP <= 1e-8, "P against mpmath", eP, 1e-8);
  Check(eS <= 1e-8, "S against mpmath", eS, 1e-8);
  Check(ePhi <= 1e-8, "exp(-i phi) against mpmath", ePhi, 1e-8);

  // 3. S' against a five-point difference of coulcc's S: p + 12C (a = 4.5 fm) and alpha + 12C (a = 5.5 fm).
  const ChannelKinematics p12c{1, 6, 1.007276 * 12.0 / (1.007276 + 12.0), 0.0, 0.0, false};
  const ChannelKinematics a12c{2, 6, 4.002603 * 12.0 / (4.002603 + 12.0), 0.0, 0.0, false};
  double e3 = 0.0;
  for (const ChannelKinematics *k : {&p12c, &a12c})
    for (int l = 0; l <= 3; ++l)
      for (double E : {0.3, 1.0, 2.5, 6.0}) {
        const double a = (k == &p12c) ? 4.5 : 5.5, h = 1e-4 * E;
        auto S = [&](double e) {
          return shift(coulomb::coulomb(l, eta(*k, e), rho(*k, a, e), coulomb::Method::coulcc), rho(*k, a, e));
        };
        const double ref = (S(E - 2 * h) - 8 * S(E - h) + 8 * S(E + h) - S(E + 2 * h)) / (12 * h);
        const double got = channel_functions(*k, l, a, E).dS_dE;
        e3 = std::fmax(e3, std::fabs(got - ref) / std::fabs(ref));
      }
  std::printf("S' against coulcc's five-point difference: %.2e relative\n", e3);
  Check(e3 <= 1e-6, "S'", e3, 1e-6);

  // 4. Kinematics.
  const double mu = p12c.mu_u, uc = kAtomicMassUnit_MeV, hc = kHbarC_MeV_fm;
  double e4 = 0.0;
  for (double E : {1e-3, 0.1, 1.0, 10.0}) {
    const double beta = std::sqrt(2.0 * E / (mu * uc)), k = std::sqrt(2.0 * mu * uc * E) / hc;
    e4 = std::fmax(e4, std::fabs(eta(p12c, E) / (6.0 * kFineStructure / beta) - 1.0));
    e4 = std::fmax(e4, std::fabs(rho(p12c, 4.5, E) / (k * 4.5) - 1.0));
    e4 = std::fmax(e4, std::fabs(eta_bound(p12c, E) / eta(p12c, E) - 1.0));
    e4 = std::fmax(e4, std::fabs(rho_bound(p12c, 4.5, E) / rho(p12c, 4.5, E) - 1.0));
  }
  const double m1 = 1.007276 * uc, m2 = 12.0 * uc;
  const ChannelKinematics rel{1, 6, mu, m1, m2, true};
  double e4r = 0.0, e4lim = 0.0;
  for (double T : {1e-3, 1.0, 50.0, 500.0}) {
    const double s = (m1 + m2 + T) * (m1 + m2 + T);
    const double k2 = (s - (m1 + m2) * (m1 + m2)) * (s - (m1 - m2) * (m1 - m2)) / (4.0 * s) / (hc * hc);
    e4r = std::fmax(e4r, std::fabs(rho(rel, 1.0, T) / std::sqrt(k2) - 1.0));
    if (T <= 1e-3) e4lim = std::fmax(std::fabs(eta(rel, T) / eta({1, 6, m1 * m2 / (m1 + m2) / uc, 0, 0, false}, T) - 1.0),
                                     std::fabs(rho(rel, 1.0, T) / rho({1, 6, m1 * m2 / (m1 + m2) / uc, 0, 0, false}, 1.0, T) - 1.0));
  }
  std::printf("kinematics: non-relativistic %.2e, relativistic k %.2e, T/m -> 0 limit %.2e\n", e4, e4r, e4lim);
  Check(e4 <= 1e-14, "non-relativistic eta, rho", e4, 1e-14);
  Check(e4r <= 1e-9, "relativistic k", e4r, 1e-9);  // the closed form cancels at small T
  Check(e4lim <= 1e-6, "relativistic limit", e4lim, 1e-6);

  // 5. Outside the domain.
  for (double E : {0.0, -1.0, (double)NAN}) {
    const OpenChannelValues v = channel_functions(p12c, 0, 4.5, E);
    if (v.status != coulomb::Status::invalid_argument || !std::isnan(v.P) || !std::isnan(v.S) || !std::isnan(v.dS_dE)) {
      std::printf("FAIL: E = %g accepted\n", E);
      ++failures;
    }
  }

  if (failures) std::printf("%d failure(s)\n", failures);
  else std::printf("channel_functions: all checks passed\n");
  return failures ? 1 : 0;
}
