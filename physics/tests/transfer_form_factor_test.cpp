// transfer_form_factor (azure2/physics/transfer/FormFactor.h; docs/upgrade/thm_form_factor.md): the full form factor
// of the THM transfer vertex, against answers that do not come from its own algebra.
//
//  1. Green's identity: W = Int u V j_l(pr) r dr, computed directly from the level function and the well, equals
//     (hbar^2/2mu) u(a) [M_l(pa; B) + Delta(eps)] for every case below, at eps on and off the level energy: 1e-8
//     relative.  It holds only if u solves the radial equation of the well with a u'(a)/u(a) = B.
//  2. The level function's boundary value and node count are those asked for, and Int_0^a u^2 = 1.
//  3. R = W/surface for the four 17O(n,alpha) levels of Oliva's thesis against an independent NumPy implementation
//     (Numerov on 1500 steps, trapezoid overlaps; ~/THM 17O_n_a/oliva2022_repro/scripts/formfactor.py of 2026-11-04,
//     whose numbers are copied below): 1e-4 relative (the two quadratures differ at that level).
//  4. SphericalJ against closed forms (j_0 = sin x/x, j_1, j_2) at small and large x: 1e-13.
//  5. WildermuthNodes: p + 18O (8 protons): l = 0 -> 1s (1), l = 2 -> 0d (0), l = 1 -> 1p (1), l = 3 -> 0f (0);
//     n + 16O (8 neutrons) the same; n + 4He (2 neutrons): l = 1 -> 0p (0), l = 0 -> 1s (1).
//  6. Chebyshev: a series fitted to cos(3x) on [0, 2] from 24 nodes reproduces it to 1e-12 between the nodes.
#include <cmath>
#include <cstdio>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/transfer/FormFactor.h"

using namespace azure2::physics;
namespace ff = azure2::physics::transfer::form_factor;

namespace {

int failures = 0;

void Check(bool ok, const char *what, double got, double want) {
  if (!ok) {
    std::printf("FAIL %s: got %.12g, want %.12g\n", what, got, want);
    failures++;
  }
}

double Rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }

}  // namespace

int main() {
  // n + 17O as in formfactor.py: masses (u), amu (MeV), radius 4.64 fm, r0 1.25 fm, a0 0.65 fm, no Coulomb.
  const double mN = 1.00866492, m17 = 16.99474312, amu = 931.494102, bD = 2.224566;
  ff::Well w;
  w.mu = mN * m17 / (mN + m17) * amu;
  w.a = 4.64;
  w.R = 1.25 * std::cbrt(17.0);
  w.diff = 0.65;
  w.rc = w.R;
  w.zze2 = 0.0;
  w.steps = 1000;
  const double h2m = w.H2M();
  struct Case {
    const char *name;
    int l, nodes;
    double e, B, R;  // E (MeV above threshold), a u'/u, W/surface from formfactor.py
  };
  const Case cases[] = {
      {"1- n1", 1, 1, -0.0054, -1.0049370181563422, 1.1969798495346522},
      {"1- n0", 1, 0, -0.0054, -1.0049370181563422, 1.3384048253377216},
      {"5- n0", 3, 0, 0.09, -2.982238199989906, 1.0949427822511701},
      {"2+ n0", 2, 0, 0.172, -1.9409505703868728, 1.1569430158631027},
      {"2+ n1", 2, 1, 0.172, -1.9409505703868728, 1.1034710457212777},
      {"3- n1", 1, 1, 0.24, -0.8094598593311532, 1.2120215258735039},
      {"3- n0", 1, 0, 0.24, -0.8094598593311532, 1.3653388687635477},
  };
  for (const Case &c : cases) {
    const ff::LevelFunction f = ff::SolveLevel(w, c.l, c.e, c.B, c.nodes);
    if (!f.ok) {
      std::printf("FAIL %s: no level function (%s)\n", c.name, f.why.c_str());
      failures++;
      continue;
    }
    // 2. boundary, nodes, norm
    double L = 0.0;
    int n = -1;
    ff::SolveAtDepth(w, c.l, c.e, f.depth, L, n);
    Check(Rel(L, c.B) < 1e-6, "boundary value", L, c.B);
    Check(n == c.nodes, "node count", n, c.nodes);
    double norm = 0.0;
    const double hstep = w.a / w.steps;
    for (int i = 0; i <= w.steps; i++)
      norm += (i == 0 || i == w.steps ? 1.0 : (i % 2 ? 4.0 : 2.0)) * hstep / 3.0 * f.u[i] * f.u[i];
    Check(std::fabs(norm - 1.0) < 1e-10, "norm", norm, 1.0);
    // 1. Green's identity at the THM momentum and at two others
    for (double eps : {c.e + bD, 0.5, 4.0}) {
      if (!(eps > 0.0)) continue;
      const double p = std::sqrt(eps / h2m);
      const double W = ff::PotentialOverlap(w, f, p);
      const double green = h2m * f.ua * (ff::SurfaceM(c.l, p * w.a, f.B) + ff::VolumeTerm(w, f, eps));
      Check(Rel(green, W) < 1e-8, "Green identity", green, W);
    }
    // 3. against formfactor.py
    const double eps = c.e + bD, p = std::sqrt(eps / h2m);
    const double R = ff::PotentialOverlap(w, f, p) / (h2m * f.ua * ff::SurfaceM(c.l, p * w.a, f.B));
    Check(Rel(R, c.R) < 1e-4, c.name, R, c.R);
  }
  // 4. spherical Bessel functions
  for (double x : {1e-3, 0.3, 2.0, 7.5, 40.0}) {
    const double s = std::sin(x), co = std::cos(x);
    const double j0 = s / x, j1 = s / (x * x) - co / x, j2 = (3.0 / (x * x) - 1.0) * s / x - 3.0 * co / (x * x);
    Check(std::fabs(ff::SphericalJ(0, x) - j0) < 1e-13, "j0", ff::SphericalJ(0, x), j0);
    Check(std::fabs(ff::SphericalJ(1, x) - j1) < 1e-13 * std::max(1.0, 1.0 / x), "j1", ff::SphericalJ(1, x), j1);
    if (x > 0.1) Check(std::fabs(ff::SphericalJ(2, x) - j2) < 1e-12, "j2", ff::SphericalJ(2, x), j2);
  }
  // 5. Wildermuth
  Check(ff::WildermuthNodes(0, 8) == 1, "W p+18O l0", ff::WildermuthNodes(0, 8), 1);
  Check(ff::WildermuthNodes(2, 8) == 0, "W p+18O l2", ff::WildermuthNodes(2, 8), 0);
  Check(ff::WildermuthNodes(1, 8) == 1, "W p+18O l1", ff::WildermuthNodes(1, 8), 1);
  Check(ff::WildermuthNodes(3, 8) == 0, "W p+18O l3", ff::WildermuthNodes(3, 8), 0);
  Check(ff::WildermuthNodes(1, 2) == 0, "W n+4He l1", ff::WildermuthNodes(1, 2), 0);
  Check(ff::WildermuthNodes(0, 2) == 1, "W n+4He l0", ff::WildermuthNodes(0, 2), 1);
  // 6. Chebyshev
  {
    const int n = 24;
    std::vector<double> x = ff::Chebyshev::Nodes(0.0, 2.0, n), v(n);
    for (int k = 0; k < n; k++) v[k] = std::cos(3.0 * x[k]);
    ff::Chebyshev c;
    c.Fit(0.0, 2.0, v);
    for (double t = 0.0; t <= 2.0; t += 0.0137) Check(std::fabs(c(t) - std::cos(3.0 * t)) < 1e-12, "chebyshev", c(t), std::cos(3.0 * t));
  }
  if (failures) {
    std::printf("%d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_form_factor: all checks passed\n");
  return 0;
}
