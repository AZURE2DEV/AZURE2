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
#include <cmath>
#include <complex>
#include <cstdio>
#include <stdexcept>
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

// p + 12C (atomic masses, u)
const double kMu_u = 1.00782503223 * 12.0 / (1.00782503223 + 12.0);
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
      const RadialValue h = IntegrateInward(p, 1, 6, kMu_MeV, l, a, E, outgoing);
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

}  // namespace

int main() {
  Values();
  const double worst = CoulombOnly(2e-7);
  std::printf("V_N = 0 against the Coulomb functions: worst %.2e\n", worst);
  if (failures) {
    std::printf("%d failure(s)\n", failures);
    return 1;
  }
  std::printf("hybrid_potential: all checks passed\n");
  return 0;
}
