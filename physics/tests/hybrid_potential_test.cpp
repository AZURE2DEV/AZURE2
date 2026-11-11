// hybrid_potential (docs/upgrade/spin_orbit_potential.md sections 1.1, 3.2): the potential as a value and the
// inward radial integration.
//
//  1. Values: two potentials built alike compare equal; each field, changed alone, makes them unequal and changes
//     PotentialKey; no potential has the key {0}.  The factories refuse a non-positive parameter.
//  2. V_N = 0 (depth 0): the solution of the radial equation with only the Coulomb (uniform sphere inside R_C) and
//     centrifugal terms, read outside R_C, is the Coulomb H+ = G + iF.  P = rho/|H|^2 and S = rho Re(H' H*)/|H|^2
//     at a = 3.6 fm of p + 12C, l = 0-2, E = 0.5-3 MeV against the library's Coulomb functions (coulsm, 1e-9): to
//     2e-7, the accuracy of the O(step^4) integration (largest 9.3e-8, at l = 2, 3 MeV; 6e-9 in P and 2e-8 in S at
//     l = 1, 0.5 MeV, work/so-note/radial.py part B).  The Coulomb term with e^2 = 1.4399764 MeV fm instead of
//     alpha hbar c fails at 2.3e-5.
//  3. Below, at and above threshold (data/hybrid_bound.csv, make_hybrid_bound_reference.py: scipy DOP853 inward from
//     15 fm with mpmath's Whittaker, Bessel K and Coulomb functions, not the engine's): S(E) of p + 12C with the
//     settings of tests/hybrid_potential (Woods-Saxon 20 MeV, 3.6 fm, 0.6 fm; a = 3.6 fm) at E = -2 ... +0.5 MeV,
//     l = 0-2, to 1e-6 absolute (HybridBoundShift below and at threshold, IntegrateInward above); the same with
//     V_N = 0 (rows "Sc": the Whittaker and Coulomb S, so V_N -> 0 reproduces them); ln W~(a), the ANC <-> gamma
//     factor, to 1e-6; int_a^inf W~^2 dr, the bound-state integral of external capture, to 1e-6 relative.
//  4. S continuous through threshold: |S(E) - S(0)| <= 2e-3 for |E| <= 1 keV (dS/dE ~ 0.5 MeV^-1).
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/HybridPotential.h"
#include "azure2/physics/coulomb/Coulomb.h"

using namespace azure2::physics;

namespace {

int failures = 0;

void Check(bool ok, const char *what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what);
    ++failures;
  }
}

void Values() {
  const HybridPotential a = WoodsSaxonPotential(20.0, 3.6, 0.6), b = WoodsSaxonPotential(20.0, 3.6, 0.6);
  Check(a == b && PotentialKey(a) == PotentialKey(b), "equal potentials compare equal");
  Check(PotentialKey(std::nullopt) == std::vector<double>{0.0}, "no potential: key {0}");
  Check(PotentialKey(a) != PotentialKey(std::nullopt), "a potential's key differs from none");
  auto differs = [&](HybridPotential c, const char *what) {
    Check(!(c == a) && PotentialKey(c) != PotentialKey(a), what);
  };
  HybridPotential c = a;
  c.coulomb_radius_fm = 3.6;
  differs(c, "R_C enters the comparison");
  c = a;
  c.central_depth_MeV = 21.0;
  differs(c, "V0 enters the comparison");
  c = a;
  c.central.radius_fm = 3.7;
  differs(c, "R enters the comparison");
  c = a;
  c.central.diffuseness_fm = 0.65;
  differs(c, "a enters the comparison");
  c = a;
  c.matching_radius_fm = 16.0;
  differs(c, "r_m enters the comparison");
  c = a;
  c.step_fm = 0.03;
  differs(c, "the step enters the comparison");
  const HybridPotential g = GaussianPotential(20.0, 3.6);
  differs(g, "the shape enters the comparison");
  HybridPotential g2 = g;
  g2.central.width_fm = 3.7;
  Check(!(g2 == g) && PotentialKey(g2) != PotentialKey(g), "r0 enters the comparison");
  Check(a.coulomb_radius_fm == 3.5 && a.matching_radius_fm == 15.0 && a.step_fm == 0.0275,
        "defaults R_C = 3.5 fm, r_m = 15 fm, step 0.0275 fm");
  int refused = 0;
  for (auto make : {+[] { WoodsSaxonPotential(0.0, 3.6, 0.6); }, +[] { WoodsSaxonPotential(20.0, -1.0, 0.6); },
                    +[] { WoodsSaxonPotential(20.0, 3.6, 0.0); }, +[] { GaussianPotential(-1.0, 3.0); },
                    +[] { GaussianPotential(20.0, 0.0); }}) {
    try {
      make();
    } catch (const std::invalid_argument &) {
      refused++;
    }
  }
  Check(refused == 5, "the factories refuse non-positive parameters");
}

// p + 12C as pair 1 of tests/hybrid_potential (and make_hybrid_bound_reference.py): masses 1.00728 and 12 u
const double kMu_u = 1.00728 * 12.0 / 13.00728;
const double kMu_MeV = kMu_u * kAtomicMassUnit_MeV;

double CoulombOnly(double tolerance) {
  HybridPotential p = WoodsSaxonPotential(1.0, 3.6, 0.6);
  p.central_depth_MeV = 0.0;
  const double a = 3.6;
  double worst = 0.0;
  for (int l = 0; l <= 2; l++)
    for (double E : {0.5, 1.0, 2.0, 3.0}) {
      const double eta = 6 * kFineStructure * std::sqrt(kMu_MeV / 2.0) / std::sqrt(E);
      const double k = std::sqrt(2.0 * kMu_MeV * E) / kHbarC_MeV_fm;  // fm^-1
      auto outgoing = [&](double r) {
        const coulomb::CoulombWaves w = coulomb::coulomb(l, eta, k * r);
        return std::complex<double>(w.G, w.F);
      };
      const RadialValue h = IntegrateInward(p, 1, 6, kMu_MeV, l, 0.0, a, E, outgoing);
      const double rho = k * a;
      const std::complex<double> dH = h.du_dr / k;
      const double P = rho / std::norm(h.u), S = rho * std::real(dH * std::conj(h.u)) / std::norm(h.u);
      const coulomb::CoulombWaves w = coulomb::coulomb(l, eta, rho);
      const double P0 = rho / (w.F * w.F + w.G * w.G), S0 = rho * (w.F * w.Fp + w.G * w.Gp) / (w.F * w.F + w.G * w.G);
      const double dP = std::fabs(P / P0 - 1.0), dS = std::fabs(S - S0);
      worst = std::max(worst, std::max(dP, dS));
      if (dP > tolerance || dS > tolerance)
        std::printf("FAIL: V_N = 0, l = %d, E = %.1f MeV: |dP/P| %.2e, |dS| %.2e (bound %.0e)\n", l, E, dP, dS,
                    tolerance);
      if (dP > tolerance || dS > tolerance) ++failures;
    }
  return worst;
}

struct RefRow {
  std::string kind;
  int l;
  double E, value;
};

std::vector<RefRow> ReadReference(const char *path) {
  std::vector<RefRow> rows;
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);  // header
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    RefRow r;
    std::string f;
    std::getline(ss, r.kind, ',');
    std::getline(ss, f, ',');
    r.l = std::stoi(f);
    std::getline(ss, f, ',');
    r.E = std::stod(f);
    std::getline(ss, f, ',');
    r.value = std::stod(f);
    rows.push_back(r);
  }
  return rows;
}

// S at E (MeV) of p + 12C with the potential p: HybridBoundShift below and at threshold, rho H'/H above.
double ShiftAt(const HybridPotential &p, int l, double E) {
  if (E <= 0.0) return HybridBoundShift(p, 1, 6, kMu_MeV, l, 0.0, 3.6, -E);
  const double eta = 6 * kFineStructure * std::sqrt(kMu_MeV / 2.0) / std::sqrt(E);
  const double k = std::sqrt(2.0 * kMu_MeV * E) / kHbarC_MeV_fm;
  auto outgoing = [&](double r) {
    const coulomb::CoulombWaves w = coulomb::coulomb(l, eta, k * r);
    return std::complex<double>(w.G, w.F);
  };
  const RadialValue h = IntegrateInward(p, 1, 6, kMu_MeV, l, 0.0, 3.6, E, outgoing);
  return 3.6 * std::real(h.du_dr * std::conj(h.u)) / std::norm(h.u);
}

void AgainstReference(const char *path) {
  const std::vector<RefRow> rows = ReadReference(path);
  if (rows.empty()) {
    std::printf("FAIL: no reference rows in %s\n", path);
    ++failures;
    return;
  }
  const HybridPotential ws = WoodsSaxonPotential(20.0, 3.6, 0.6);
  HybridPotential none = ws;
  none.central_depth_MeV = 0.0;
  double worstS = 0.0, worstSc = 0.0, worstW = 0.0, worstN = 0.0;
  for (const RefRow &r : rows) {
    if (r.kind == "S" || r.kind == "Sc") {
      const double s = ShiftAt(r.kind == "S" ? ws : none, r.l, r.E);
      const double d = std::fabs(s - r.value);
      (r.kind == "S" ? worstS : worstSc) = std::max(r.kind == "S" ? worstS : worstSc, d);
      if (d > 1e-6) {
        std::printf("FAIL: %s l = %d E = %g: %.10f against %.10f\n", r.kind.c_str(), r.l, r.E, s, r.value);
        ++failures;
      }
    } else if (r.kind == "logW") {
      const DecayingValue v = DecayingSolution(ws, 1, 6, kMu_MeV, r.l, 0.0, 3.6, -r.E);
      const double d = std::fabs(v.log_W - r.value);
      worstW = std::max(worstW, d);
      if (d > 1e-6 || v.sign != 1.0) {
        std::printf("FAIL: ln W~(a) l = %d E = %g: %.10f against %.10f\n", r.l, r.E, v.log_W, r.value);
        ++failures;
      }
    } else if (r.kind == "N") {
      // int_a^inf W~^2 dr: Simpson on [a, 15] (W~ from DecayingSolution at each node) and on [15, 60] (W itself, the
      // same function beyond r_m); the tail beyond 60 fm is below 1e-20 of the integral for these bindings.
      auto W2 = [&](double x) {
        const DecayingValue v = DecayingSolution(ws, 1, 6, kMu_MeV, r.l, 0.0, x, -r.E);
        return std::exp(2.0 * v.log_W);
      };
      auto simpson = [&](double a, double b, int n) {
        const double h = (b - a) / n;
        double s = W2(a) + W2(b);
        for (int i = 1; i < n; i++) s += (i % 2 ? 4.0 : 2.0) * W2(a + i * h);
        return s * h / 3.0;
      };
      const double N = simpson(3.6, 15.0, 228) + simpson(15.0, 60.0, 450);
      const double d = std::fabs(N / r.value - 1.0);
      worstN = std::max(worstN, d);
      if (d > 1e-6) {
        std::printf("FAIL: N l = %d E = %g: %.10e against %.10e\n", r.l, r.E, N, r.value);
        ++failures;
      }
    }
  }
  std::printf("against the DOP853 reference: |dS| %.1e, |dS| (V_N = 0) %.1e, |d ln W~(a)| %.1e, |dN/N| %.1e\n", worstS,
              worstSc, worstW, worstN);
}

// S continuous through threshold: |S(E) - S(0)| for |E| <= 1 keV, both sides (the closed band eta > 100 is the
// engine's ChannelFunc; here E = -1e-3 .. 1e-3 MeV on the physics functions).
void Continuity() {
  const HybridPotential ws = WoodsSaxonPotential(20.0, 3.6, 0.6);
  for (int l = 0; l <= 2; l++) {
    const double s0 = ShiftAt(ws, l, 0.0);
    double worst = 0.0;
    for (double E : {-1e-3, -1e-4, -1e-6, 1e-4, 1e-3}) worst = std::max(worst, std::fabs(ShiftAt(ws, l, E) - s0));
    if (worst > 2e-3) {
      std::printf("FAIL: l = %d, S not continuous at threshold: %.2e\n", l, worst);
      ++failures;
    }
  }
}

}  // namespace

int main(int argc, char **argv) {
  Values();
  const double worst = CoulombOnly(2e-7);
  std::printf("V_N = 0 against the Coulomb functions: worst %.2e\n", worst);
  if (argc < 2) {
    std::printf("FAIL: usage: physics_hybrid_potential_test data/hybrid_bound.csv\n");
    return 1;
  }
  AgainstReference(argv[1]);
  Continuity();
  if (failures) {
    std::printf("%d failure(s)\n", failures);
    return 1;
  }
  std::printf("hybrid_potential: all checks passed\n");
  return 0;
}
