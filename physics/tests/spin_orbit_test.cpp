// spin_orbit (docs/upgrade/spin_orbit_potential.md sections 2.1-2.2 and 4.2, tests 1-3): the spin-orbit term of the
// hybrid potential, each part against an answer that does not use the engine's formula for it.
//
//  1. Angular factor: <2 l.i> = TwoLDotI(l, 2j, 2i) equals 2 <(i l) j m = j| l.i |(i l) j j> built in an m-scheme from
//     l_z i_z + (l_+ i_- + l_- i_+)/2 and Clebsch-Gordan coefficients, for 2i = 1..4, l = 0..4, every j: catches the
//     factor 2 (l.i instead of 2 l.i) and the sign.  TwoLDotI refuses a j that l and i cannot make.
//  2. The radial term: SpinOrbitRadialTerm(r) = V_so 2.000 fm^2 (1/r) df/dr with df/dr the central difference of the
//     Fermi function (h = 1e-5 fm), to 1e-8 relative: catches the factor 2 of (hbar/m_pi c)^2 and the sign of df/dr;
//     V_so > 0 lowers j = l + i (V(r) of l = 1, j = 3/2 below j = 1/2 at R_so).
//  3. Closed form (radial.py C of the note): a neutral square well (30 MeV, R = 4.5 fm) with the Thomas term of its
//     sharp edge, a delta shell -V_so 2.000 fm^2 <2 l.i> delta(r - R)/R (a jump of u'), solved with Riccati-Bessel
//     functions, against IntegrateInward on Woods-Saxon edges of diffuseness d = 0.01 and 0.005 fm (step d/20)
//     extrapolated linearly to d = 0 (Richardson: the kink of u makes the error O(d)); n + 12C masses, a = 3 fm,
//     E = 2 MeV, l = 1, 2, both j: S, P, phi to 1e-3 (relative to max(1, |S|, |P|)).
//  4. Against data/spin_orbit.csv (make_spin_orbit_reference.py: scipy DOP853 + mpmath Coulomb and Whittaker, not the
//     engine's routines): p + 12C, Woods-Saxon 20/3.6/0.6 with V_so = 6 MeV on the proton, l = 0-2, both j, E = 0.5-3
//     MeV (S, P, phi by IntegrateInward) and E = -1.94351, -0.5 MeV (S by HybridBoundShift, ln W~(a) by
//     DecayingSolution), to 1e-6 (P relative, the others absolute).
//  5. <2 l.i> = 0 (s wave) with a spin-orbit term gives the bits of the central potential (IntegrateInward,
//     DecayingSolution), and RadialPotential is bit-identical on the grid: the term is added last, never as + 0.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/AngularCoupling.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/HybridPotential.h"
#include "azure2/physics/coulomb/Coulomb.h"

using namespace azure2::physics;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    ++failures;
  }
}

// ---------------------------------------------------------------- 1. angular factor
double MSchemeTwoLDotI(int l, int two_i, int two_j) {
  // |(i l) j j> = sum (i mi l ml | j j) |l ml>|i mi> (spin first, AZURE2's order; the diagonal value does not
  // depend on it).  Twice the m-scheme expectation value of l.i.
  auto amp = [&](int two_ml, int two_mi) { return clebsch_gordan(two_i, two_mi, 2 * l, two_ml, two_j, two_j); };
  auto ladder = [](int two_j, int two_m, int step) {  // sqrt(j(j+1) - m(m + step)) for J_+ (step +2) or J_- (-2)
    const double j = two_j / 2.0, m = two_m / 2.0;
    return std::sqrt(j * (j + 1) - m * (m + step / 2.0));
  };
  double sum = 0.0;
  for (int two_ml = -2 * l; two_ml <= 2 * l; two_ml += 2)
    for (int two_mi = -two_i; two_mi <= two_i; two_mi += 2) {
      const double c = amp(two_ml, two_mi);
      if (c == 0.0) continue;
      sum += c * c * (two_ml / 2.0) * (two_mi / 2.0);  // l_z i_z
      // (l_+ i_- + l_- i_+)/2 takes (ml, mi) to (ml +- 1, mi -+ 1)
      if (two_ml + 2 <= 2 * l && two_mi - 2 >= -two_i)
        sum += 0.5 * amp(two_ml + 2, two_mi - 2) * ladder(2 * l, two_ml, 2) * ladder(two_i, two_mi, -2) * c;
      if (two_ml - 2 >= -2 * l && two_mi + 2 <= two_i)
        sum += 0.5 * amp(two_ml - 2, two_mi + 2) * ladder(2 * l, two_ml, -2) * ladder(two_i, two_mi, 2) * c;
    }
  return 2.0 * sum;
}

void Angular() {
  double worst = 0.0;
  int cases = 0;
  for (int two_i = 1; two_i <= 4; two_i++)
    for (int l = 0; l <= 4; l++)
      for (int two_j = std::abs(2 * l - two_i); two_j <= 2 * l + two_i; two_j += 2) {
        const double d = std::fabs(TwoLDotI(l, two_j, two_i) - MSchemeTwoLDotI(l, two_i, two_j));
        worst = std::max(worst, d);
        cases++;
        if (d > 1e-12) std::printf("FAIL: <2 l.i> l = %d 2i = %d 2j = %d: %g against the m-scheme %g\n", l, two_i, two_j,
                                   TwoLDotI(l, two_j, two_i), MSchemeTwoLDotI(l, two_i, two_j));
      }
  Check(worst <= 1e-12, "TwoLDotI equals twice the m-scheme <l.i>");
  Check(TwoLDotI(1, 3, 1) == 1.0 && TwoLDotI(1, 1, 1) == -2.0 && TwoLDotI(2, 5, 1) == 2.0, "spin 1/2: l and -(l + 1)");
  int refused = 0;
  for (auto t : {std::vector<int>{1, 5, 1}, {1, 2, 1}, {0, 3, 1}, {2, 1, 1}}) {
    try {
      TwoLDotI(t[0], t[1], t[2]);
    } catch (const std::invalid_argument &) {
      refused++;
    }
  }
  Check(refused == 4, "TwoLDotI refuses j outside |l - i| .. l + i and a 2j of the wrong parity");
  std::printf("1. angular: %d cases, worst %.1e\n", cases, worst);
}

// ---------------------------------------------------------------- 2. the radial term
void RadialTerm() {
  const SpinOrbitTerm t = WoodsSaxonSpinOrbit(6.0, 3.6, 0.6, SpinOrbitParticle::Light);
  auto fermi = [](double r) { return 1.0 / (1.0 + std::exp((r - 3.6) / 0.6)); };
  double worst = 0.0;
  for (double r : {0.5, 2.0, 3.0, 3.6, 4.2, 6.0, 10.0}) {
    const double h = 1e-5;
    const double expected = 6.0 * 2.000 * (fermi(r + h) - fermi(r - h)) / (2.0 * h) / r;
    worst = std::max(worst, std::fabs(SpinOrbitRadialTerm(t, r) / expected - 1.0));
  }
  Check(worst < 1e-8, "v_so(r) = V_so 2.000 fm^2 (1/r) df/dr (central difference)");
  HybridPotential p = WoodsSaxonPotential(20.0, 3.6, 0.6);
  p.spin_orbit = t;
  const double mu = 1.00728 * 12.0 / 13.00728 * kAtomicMassUnit_MeV;
  Check(RadialPotential(p, 1, 6, mu, 1, TwoLDotI(1, 3, 1), 3.6) < RadialPotential(p, 1, 6, mu, 1, TwoLDotI(1, 1, 1), 3.6),
        "V_so > 0 lowers j = l + 1/2");
  std::printf("2. radial term: worst %.1e\n", worst);
}

// ---------------------------------------------------------------- 3. delta-shell closed form
struct Riccati {
  double F, dF, G, dG;  // F = x j_l(x), G = -x y_l(x) and their x-derivatives
};
Riccati RiccatiAt(int l, double x) {
  // u_0 = sin x, cos x; u_l = (2l - 1)/x u_{l-1} - u_{l-2}; u_l' = u_{l-1} - l u_l / x (both kinds)
  double f0 = std::sin(x), g0 = std::cos(x), f1 = std::sin(x) / x - std::cos(x), g1 = std::cos(x) / x + std::sin(x);
  if (l == 0) return {f0, std::cos(x), g0, -std::sin(x)};
  double fp = f0, gp = g0, f = f1, g = g1;
  for (int k = 2; k <= l; k++) {
    const double fn = (2 * k - 1) / x * f - fp, gn = (2 * k - 1) / x * g - gp;
    fp = f, gp = g, f = fn, g = gn;
  }
  return {f, fp - l * f / x, g, gp - l * g / x};
}

struct Values {
  double S, P, phi;
};

Values FromSolution(std::complex<double> u, std::complex<double> du_drho, double rho) {
  const std::complex<double> L = rho * du_drho / u;
  return {L.real(), L.imag(), std::atan2(u.imag(), u.real())};
}

Values DeltaShell(double mu, double V0, double R, double Vso, double ldots, int l, double a, double E) {
  const double k = std::sqrt(2.0 * mu * E) / kHbarC_MeV_fm;
  const double q = std::sqrt(k * k + 2.0 * mu * V0 / (kHbarC_MeV_fm * kHbarC_MeV_fm));
  const double jump = 2.0 * mu / (kHbarC_MeV_fm * kHbarC_MeV_fm) * (-Vso * 2.000 * ldots / R);  // u'(R+) - u'(R-) = jump u(R)
  const Riccati out = RiccatiAt(l, k * R), in = RiccatiAt(l, q * R);
  const std::complex<double> H(out.G, out.F), dH(out.dG, out.dF);
  // alpha F(qR) + beta G(qR) = H(kR);  q (alpha F'(qR) + beta G'(qR)) = k H'(kR) - jump H(kR)
  const std::complex<double> b1 = H, b2 = k * dH - jump * H;
  const double det = in.F * q * in.dG - in.G * q * in.dF;
  const std::complex<double> alpha = (b1 * q * in.dG - in.G * b2) / det, beta = (in.F * b2 - q * in.dF * b1) / det;
  const Riccati at = RiccatiAt(l, q * a);
  const std::complex<double> u = alpha * at.F + beta * at.G, du = q * (alpha * at.dF + beta * at.dG);
  return FromSolution(u, du / k, k * a);
}

Values Engine(const HybridPotential &p, double mu, int l, double ldots, double a, double E) {
  const double k = std::sqrt(2.0 * mu * E) / kHbarC_MeV_fm;
  auto outgoing = [&](double r) {
    const Riccati w = RiccatiAt(l, k * r);
    return std::complex<double>(w.G, w.F);
  };
  const RadialValue h = IntegrateInward(p, 0, 6, mu, l, ldots, a, E, outgoing);
  return FromSolution(h.u, h.du_dr / k, k * a);
}

void DeltaShellClosedForm() {
  const double mu = 1.00866 * 12.0 / 13.00866 * kAtomicMassUnit_MeV;
  double worst = 0.0;
  for (int l : {1, 2})
    for (int two_j : {2 * l - 1, 2 * l + 1}) {
      const double ldots = TwoLDotI(l, two_j, 1);
      const Values x = DeltaShell(mu, 30.0, 4.5, 7.0, ldots, l, 3.0, 2.0);
      Values v[2];
      int n = 0;
      for (double d : {0.01, 0.005}) {
        HybridPotential p = WoodsSaxonPotential(30.0, 4.5, d);
        p.spin_orbit = WoodsSaxonSpinOrbit(7.0, 4.5, d, SpinOrbitParticle::Light);
        p.step_fm = d / 20.0;
        v[n++] = Engine(p, mu, l, ldots, 3.0, 2.0);
      }
      const Values r{2 * v[1].S - v[0].S, 2 * v[1].P - v[0].P, 2 * v[1].phi - v[0].phi};
      const double scale = std::max(1.0, std::max(std::fabs(x.S), std::fabs(x.P)));
      const double err = std::max({std::fabs(r.S - x.S), std::fabs(r.P - x.P), std::fabs(r.phi - x.phi)}) / scale;
      worst = std::max(worst, err);
      std::printf("3. l = %d 2j = %d: closed form S %.6f P %.6f phi %.6f; Richardson %.6f %.6f %.6f; %.1e\n", l, two_j,
                  x.S, x.P, x.phi, r.S, r.P, r.phi, err);
    }
  Check(worst < 1e-3, "the delta-shell closed form against the solver on a sharp Woods-Saxon edge, Richardson in d");
}

// ---------------------------------------------------------------- 4. against the independent integration
const double kMu_MeV = 1.00728 * 12.0 / 13.00728 * kAtomicMassUnit_MeV;

void AgainstReference(const char *path) {
  HybridPotential p = WoodsSaxonPotential(20.0, 3.6, 0.6);
  p.spin_orbit = WoodsSaxonSpinOrbit(6.0, 3.6, 0.6, SpinOrbitParticle::Light);
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);
  int rows = 0;
  double worst = 0.0;
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string kind, f;
    std::getline(ss, kind, ',');
    std::getline(ss, f, ',');
    const int l = std::stoi(f);
    std::getline(ss, f, ',');
    const int two_j = std::stoi(f);
    std::getline(ss, f, ',');
    const double E = std::stod(f);
    std::getline(ss, f, ',');
    const double ref = std::stod(f);
    const double ldots = TwoLDotI(l, two_j, 1);
    double value = 0.0, err = 0.0;
    if (E > 0.0) {
      const double eta = 6 * kFineStructure * std::sqrt(kMu_MeV / 2.0) / std::sqrt(E);
      const double k = std::sqrt(2.0 * kMu_MeV * E) / kHbarC_MeV_fm;
      auto outgoing = [&](double r) {
        const coulomb::CoulombWaves w = coulomb::coulomb(l, eta, k * r);
        return std::complex<double>(w.G, w.F);
      };
      const RadialValue h = IntegrateInward(p, 1, 6, kMu_MeV, l, ldots, 3.6, E, outgoing);
      const Values v = FromSolution(h.u, h.du_dr / k, k * 3.6);
      value = kind == "S" ? v.S : kind == "P" ? v.P : v.phi;
      err = kind == "P" ? std::fabs(value / ref - 1.0) : std::fabs(value - ref);
    } else {
      value = kind == "S" ? HybridBoundShift(p, 1, 6, kMu_MeV, l, ldots, 3.6, -E)
                          : DecayingSolution(p, 1, 6, kMu_MeV, l, ldots, 3.6, -E).log_W;
      err = std::fabs(value - ref);
    }
    rows++;
    worst = std::max(worst, err);
    if (err > 1e-6) std::printf("FAIL: %s l = %d 2j = %d E = %g: %.10g against %.10g\n", kind.c_str(), l, two_j, E, value, ref);
    if (err > 1e-6) ++failures;
  }
  Check(rows == 80, "80 reference rows read");
  std::printf("4. against make_spin_orbit_reference.py: %d rows, worst %.1e\n", rows, worst);
}

// ---------------------------------------------------------------- 5. no term, the same bits
void ZeroFactorSameBits() {
  const HybridPotential central = WoodsSaxonPotential(20.0, 3.6, 0.6);
  HybridPotential so = central;
  so.spin_orbit = WoodsSaxonSpinOrbit(6.0, 3.6, 0.6, SpinOrbitParticle::Light);
  const double ldots = TwoLDotI(0, 1, 1);  // 0 for the s wave
  Check(ldots == 0.0, "<2 l.i> = 0 for l = 0");
  bool same = true;
  for (double r = 0.0275; r < 16.0; r += 0.0275)
    same = same && RadialPotential(central, 1, 6, kMu_MeV, 0, 0.0, r) == RadialPotential(so, 1, 6, kMu_MeV, 0, ldots, r);
  Check(same, "RadialPotential with <2 l.i> = 0 is bit-identical to the central one");
  const double eta = 6 * kFineStructure * std::sqrt(kMu_MeV / 2.0) / std::sqrt(1.0);
  const double k = std::sqrt(2.0 * kMu_MeV) / kHbarC_MeV_fm;
  auto outgoing = [&](double r) {
    const coulomb::CoulombWaves w = coulomb::coulomb(0, eta, k * r);
    return std::complex<double>(w.G, w.F);
  };
  const RadialValue a = IntegrateInward(central, 1, 6, kMu_MeV, 0, 0.0, 3.6, 1.0, outgoing);
  const RadialValue b = IntegrateInward(so, 1, 6, kMu_MeV, 0, ldots, 3.6, 1.0, outgoing);
  Check(std::memcmp(&a.u, &b.u, sizeof a.u) == 0 && std::memcmp(&a.du_dr, &b.du_dr, sizeof a.du_dr) == 0,
        "IntegrateInward: s wave with a spin-orbit term = central (bits)");
  const DecayingValue c = DecayingSolution(central, 1, 6, kMu_MeV, 0, 0.0, 3.6, 0.5);
  const DecayingValue d = DecayingSolution(so, 1, 6, kMu_MeV, 0, ldots, 3.6, 0.5);
  Check(c.log_W == d.log_W && c.rho_dlog == d.rho_dlog, "DecayingSolution: s wave with a spin-orbit term = central (bits)");
  Check(central != so && PotentialKey(central) != PotentialKey(so) && PotentialKey(central).size() == 9,
        "the term enters the comparison and the key; a central potential keeps its 9 numbers");
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: spin_orbit_test data/spin_orbit.csv\n");
    return 2;
  }
  Angular();
  RadialTerm();
  DeltaShellClosedForm();
  AgainstReference(argv[1]);
  ZeroFactorSameBits();
  std::printf(failures ? "FAILED (%d)\n" : "PASSED\n", failures);
  return failures ? 1 : 0;
}
