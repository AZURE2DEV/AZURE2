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
//  6. The depth and width for a chosen (a, r): from VD's s wave (67.5 MeV, 2.52 fm) the engine finds back the (70 MeV,
//     2.4 fm) whose a and r scipy computed (40.705 fm, 0.998 fm), from n + 4He (40 MeV, 2 fm) the (35 MeV, 2.2 fm), to
//     3e-4 (the conditioning of (V0, r0) in (a, r) times the 1e-5 between the two expansions' r); the
//     branch kept; an r out of reach (5 fm) reported as not found with the range of r it found, no pair returned; a
//     Woods-Saxon term refused.
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
  // 6. depth and width for (a, r)
  for (const std::string sys : {"he3", "nalpha"}) {
    std::map<std::string, double> v;
    for (const Row &r : rows)
      if (r.system == sys && r.kind.rfind("dw_", 0) == 0) v[r.kind] = r.value;
    if (v.size() != 5) {
      Check(false, sys + ": the dw_ rows of the reference");
      continue;
    }
    const System &s = SystemOf(sys);
    const HybridPotential start = s.Potential(sys == "he3" ? 67.5 : 40.0);
    const DepthWidthForRange d = DepthWidthForEffectiveRange(start, s.z1, s.z2, s.Mu(), 0, 0.0, v["dw_a"], v["dw_r"]);
    char what[400];
    std::snprintf(what, sizeof what, "%s: (a, r) = (%.6g, %.6g) found at V0 = %.8g MeV, r0 = %.8g fm (scipy's %.8g, %.8g); "
                  "a, r there %.8g, %.8g; %d bound (scipy %d) %s", sys.c_str(), v["dw_a"], v["dw_r"], d.depth_MeV,
                  d.width_fm, v["dw_depth"], v["dw_width"], d.a, d.r, d.boundStates, (int)v["dw_bound"], d.why.c_str());
    // the solve reaches the asked a and r to 1e-6 (its own check); the pair agrees with scipy's to 3e-4: (V0, r0) is
    // ill-conditioned in (a, r), and the two expansions' r differ by about 1e-5 (section 13.1)
    Check(d.found && std::fabs(d.a / v["dw_a"] - 1) < 1e-6 && std::fabs(d.r / v["dw_r"] - 1) < 1e-5 &&
              std::fabs(d.depth_MeV / v["dw_depth"] - 1) < 3e-4 && std::fabs(d.width_fm / v["dw_width"] - 1) < 3e-4 &&
              d.boundStates == (int)v["dw_bound"],
          what);
  }
  {
    const DepthWidthForRange d = DepthWidthForEffectiveRange(kHe3.Potential(67.5), 2, 2, kHe3.Mu(), 0, 0.0, 30.95, 5.0);
    Check(!d.found && d.why.find("r runs from") != std::string::npos && d.depth_MeV == 0.0,
          "r = 5 fm with a0 = 30.95 fm out of reach: not found, no pair returned: " + d.why);
    HybridPotential ws = WoodsSaxonPotential(60.0, 2.5, 0.6);
    const DepthWidthForRange e = DepthWidthForEffectiveRange(ws, 0, 2, kNAlpha.Mu(), 0, 0.0, 3.0, 1.5);
    Check(!e.found && e.why.find("Gaussian") != std::string::npos, "a Woods-Saxon term refused: " + e.why);
  }
  const EffectiveRange p = PotentialEffectiveRange(kHe3.Potential(83.589554), 2, 2, kHe3.Mu(), 1, 0.0);
  Check(!p.ok && p.why.find("l = 0 only") != std::string::npos, "charged l = 1 refused: " + p.why);
  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
