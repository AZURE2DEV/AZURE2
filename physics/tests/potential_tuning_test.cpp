// potential_tuning (docs/upgrade/potential_background.md section 13.1): what a pair's potential gives by itself at low
// energy, against tests/data/potential_tuning.csv (make_potential_tuning_reference.py: scipy's DOP853, mpmath's Coulomb
// functions and digamma, numpy's fit on another window with another degree; nothing of the engine).
//
//  1. h(eta) = Re psi(i eta) - ln eta against mpmath, on both sides of the switch to the asymptotic series (eta = 12).
//  2. The potential's own phase shift (3He + 4He, the Gaussian VD of Tursunov et al.; n + 4He, a made-up Gaussian; l = 0
//     and 1) against scipy, to 1e-6 rad.
//  3. a and r of the effective-range expansion (Coulomb-modified for 3He + 4He l = 0; k^(2l+1) cot delta for n + 4He
//     l = 0, 1): to 3e-5 relative for a, 1e-3 for r (two windows and two fits: what an extrapolation to k = 0 allows),
//     and the stated spread of the engine's fit covers the difference to within a factor 10.
//  4. The number of bound states by the zero-energy solution and S_l(0) at r_m against Sturm's count of scipy's
//     solution at E = -0.01 MeV out to 400 fm, for deep, shallow and unbound waves.
//  5. The depth that gives a0 = 40 fm in VD's branch (two bound s states) against scipy's brentq, to 1e-6 relative; the
//     a0 there; the branch kept.  A charged l = 1 wave is refused with its reason.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/HybridPotential.h"
#include "azure2/physics/PotentialTuning.h"

using namespace azure2::physics;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) ++failures;
}

struct System {
  double m1, m2;
  int z1, z2;
  double rc, r0;
  double Mu() const { return m1 * m2 / (m1 + m2) * kAtomicMassUnit_MeV; }
  HybridPotential Potential(double depth) const {
    HybridPotential p = GaussianPotential(depth, r0);
    p.coulomb_radius_fm = rc;
    return p;
  }
};

const System kHe3{3.0160293, 4.0026032, 2, 2, 3.095, 1.0 / std::sqrt(0.15747)};
const System kNAlpha{1.00866491595, 4.0026032, 0, 2, 0.0, 2.0};

struct Row {
  std::string kind, system;
  int l;
  double x, value;
};

std::vector<Row> ReadRows(const std::string &path) {
  std::vector<Row> rows;
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);  // header
  while (std::getline(in, line)) {
    std::stringstream s(line);
    Row r;
    std::string l, x, v;
    std::getline(s, r.kind, ',');
    std::getline(s, r.system, ',');
    std::getline(s, l, ',');
    std::getline(s, x, ',');
    std::getline(s, v, ',');
    r.l = std::stoi(l);
    r.x = std::stod(x);
    r.value = std::stod(v);
    rows.push_back(r);
  }
  return rows;
}

const System &SystemOf(const std::string &name) { return name == "he3" ? kHe3 : kNAlpha; }

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: potential_tuning_test potential_tuning.csv\n");
    return 2;
  }
  const std::vector<Row> rows = ReadRows(argv[1]);
  Check(rows.size() > 20, "reference rows read (" + std::to_string(rows.size()) + ")");
  std::map<std::string, double> depthOf = {{"he3,0", 67.5}, {"he3,1", 83.589554}, {"nalpha,0", 40.0}, {"nalpha,1", 40.0}};
  for (const Row &r : rows) {
    const System &s = SystemOf(r.system);
    char what[256];
    if (r.kind == "h") {
      const double h = CoulombH(r.x);
      std::snprintf(what, sizeof what, "h(%g) = %.15g against mpmath %.15g", r.x, h, r.value);
      Check(std::fabs(h - r.value) <= 1e-12 * std::fabs(r.value), what);
    } else if (r.kind == "delta") {
      const double v0 = depthOf[r.system + "," + std::to_string(r.l)];
      const double d = PotentialPhaseShift(s.Potential(v0), s.z1, s.z2, s.Mu(), r.l, 0.0, r.x);
      std::snprintf(what, sizeof what, "%s l = %d, E = %g MeV: delta = %.10g rad against scipy %.10g", r.system.c_str(), r.l,
                    r.x, d, r.value);
      Check(std::fabs(d - r.value) <= 1e-6 * std::max(1.0, std::fabs(r.value)) &&
                std::fabs(d - r.value) <= 1e-4 * std::fabs(r.value) + 1e-9,
            what);
    } else if (r.kind == "a" || r.kind == "r") {
      const EffectiveRange e = PotentialEffectiveRange(s.Potential(r.x), s.z1, s.z2, s.Mu(), r.l, 0.0);
      const double got = r.kind == "a" ? e.a : e.r, spread = r.kind == "a" ? e.aSpread : e.rSpread;
      const double tol = r.kind == "a" ? 3e-5 : 1e-3;
      std::snprintf(what, sizeof what, "%s l = %d, V0 = %g MeV: %s = %.8g (spread %.2g; window %.4g-%.4g MeV) against %.8g",
                    r.system.c_str(), r.l, r.x, r.kind.c_str(), got, spread, e.eLow, e.eHigh, r.value);
      Check(e.ok && std::fabs(got - r.value) <= tol * std::fabs(r.value), what);
      std::snprintf(what, sizeof what, "  the stated spread %.2g is within a factor 10 of the difference %.2g or both below "
                    "1e-5 relative", spread, std::fabs(got - r.value));
      Check(std::fabs(got - r.value) <= 10.0 * spread + 1e-5 * std::fabs(r.value), what);
    } else if (r.kind == "bound") {
      const BoundStateCount b = PotentialBoundStates(s.Potential(r.x), s.z1, s.z2, s.Mu(), r.l, 0.0);
      std::snprintf(what, sizeof what, "%s l = %d, V0 = %g MeV: %d bound state(s) (%d nodes, r u'/u %.4g vs S(0) %.4g) against %d",
                    r.system.c_str(), r.l, r.x, b.bound, b.nodes, b.logDerivative, b.exterior, (int)r.value);
      Check(b.bound == (int)r.value, what);
    } else if (r.kind == "depth") {
      const DepthForLength d = DepthForScatteringLength(s.Potential(67.5), s.z1, s.z2, s.Mu(), r.l, 0.0, r.x);
      std::snprintf(what, sizeof what, "%s: the depth for a0 = %g fm is %.10g MeV (a0 %.8g, %d bound, %d nodes) against "
                    "brentq %.10g", r.system.c_str(), r.x, d.depth_MeV, d.a, d.boundStates, d.nodes, r.value);
      Check(d.found && std::fabs(d.depth_MeV - r.value) <= 1e-6 * r.value && std::fabs(d.a - r.x) <= 1e-6 * r.x &&
                d.boundStates == PotentialBoundStates(s.Potential(67.5), s.z1, s.z2, s.Mu(), 0, 0.0).bound,
            what);
    }
  }
  const EffectiveRange p = PotentialEffectiveRange(kHe3.Potential(83.589554), 2, 2, kHe3.Mu(), 1, 0.0);
  Check(!p.ok && p.why.find("l = 0 only") != std::string::npos, "charged l = 1 refused: " + p.why);
  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
