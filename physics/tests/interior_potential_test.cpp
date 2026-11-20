// The interior R-function of a hybrid potential (docs/upgrade/interior_potential.md): RegularSolution and
// InteriorChannelFunctions of physics/HybridPotential.h.
//
//  1. The regular solution of alpha + 12C in the Woods-Saxon potential of Koonin, Tombrello and Fox (45.891 MeV,
//     3.55 fm, 0.60 fm, R_C = 3.55 fm), l = 0-2, at a = 4.0, 5.35 and 6.5 fm and E = -1 ... 5 MeV, against an
//     independent integration (data/interior_potential.csv, make_interior_potential_reference.py: scipy DOP853):
//     atan f to 2e-6 and int_0^a u^2 dr / u(a)^2 to 2e-5 relative at the default step (0.0275 fm; largest 1.1e-6
//     and 1.4e-5), the O(step^4) error of the integration in a 46 MeV well: with a quarter of the step both are
//     below 1e-7 (checked for two rows, which also holds the order).
//  2. Without levels the collision function exp(-2 i phi') is the scattering of the potential alone: its phase
//     against the phase shift of the regular solution matched to Coulomb functions at 15 fm (the same file), for the
//     three channel radii and two boundary constants, to 4e-6 rad (largest 2.1e-6, on the l = 1 resonance).
//  3. The pole of R_pot at a = 5.35 fm, b = 0, l = 1 (KTF's Fig. 1: 1.926 MeV, theta^2 = 0.94): energy to 2e-6 MeV,
//     reduced width from the norm (note Eq. (6)) to 1e-5.
//  4. One level and R_pot: U from R = R_pot + gamma^2/(E_lambda - E) directly against the level-matrix form with the
//     primed functions, to 1e-12 (an identity).
//  5. Limits: f = infinity returns the outer functions; P' = 0 at a pole of R_pot; P' >= 0.
//  6. The wave on its grid (RegularWaveOnGrid, InteriorOverlap; note section 6): in a flat well without Coulomb, l = 0,
//     u = sin(K r), the overlap of the waves of two energies against the closed form to 1e-5 (1e-7 at a quarter of the
//     step); the overlap of a wave
//     with itself is the norm of RegularSolution, which is the wave's own f and norm bit for bit; two grids are refused.
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <limits>
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

const double kMu_MeV = 4.00260325413 * 12.0 / 16.00260325413 * kAtomicMassUnit_MeV;
const double kPi = std::acos(-1.0);

HybridPotential Ktf() {
  HybridPotential p = WoodsSaxonPotential(45.891, 3.55, 0.60);
  p.coulomb_radius_fm = 3.55;
  return p;
}

// S, P and exp(-i phi) at a (fm) and E > 0 (MeV) of the outgoing solution of the potential.
ChannelTriple Outer(const HybridPotential &p, int l, double a, double E) {
  const double eta = 12 * kFineStructure * std::sqrt(kMu_MeV / 2.0) / std::sqrt(E);
  const double k = std::sqrt(2.0 * kMu_MeV * E) / kHbarC_MeV_fm;
  auto outgoing = [&](double r) {
    const coulomb::CoulombWaves w = coulomb::coulomb(l, eta, k * r);
    return std::complex<double>(w.G, w.F);
  };
  const RadialValue h = IntegrateInward(p, 2, 6, kMu_MeV, l, 0.0, a, E, outgoing);
  ChannelTriple t;
  t.P = k * a / std::norm(h.u);
  t.S = a * std::real(h.du_dr * std::conj(h.u)) / std::norm(h.u);
  t.exp_minus_i_phi = std::conj(h.u) / std::abs(h.u);
  return t;
}

double Wrap(double x) { return x - kPi * std::floor(x / kPi + 0.5); }  // into [-pi/2, pi/2)

struct RefRow {
  std::string kind;
  int l;
  double a, E, value;
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
    r.a = std::stod(f);
    std::getline(ss, f, ',');
    r.E = std::stod(f);
    std::getline(ss, f, ',');
    r.value = std::stod(f);
    rows.push_back(r);
  }
  return rows;
}

void AgainstReference(const char *path) {
  const std::vector<RefRow> rows = ReadReference(path);
  if (rows.empty()) {
    std::printf("FAIL: no reference rows in %s\n", path);
    ++failures;
    return;
  }
  const HybridPotential p = Ktf();
  double worstF = 0.0, worstN = 0.0, worstD = 0.0, poleE = 0.0;
  for (const RefRow &r : rows) {
    if (r.kind == "f" || r.kind == "norm") {
      const RegularValue v = RegularSolution(p, 2, 6, kMu_MeV, r.l, 0.0, r.a, r.E);
      const double d = r.kind == "f" ? std::fabs(std::atan(v.log_derivative) - std::atan(r.value))
                                     : std::fabs(v.norm_fm / r.value - 1.0);
      (r.kind == "f" ? worstF : worstN) = std::max(r.kind == "f" ? worstF : worstN, d);
      if (d > (r.kind == "f" ? 2e-6 : 2e-5)) {
        std::printf("FAIL: %s l = %d a = %g E = %g: %.10e against %.10e\n", r.kind.c_str(), r.l, r.a, r.E,
                    r.kind == "f" ? v.log_derivative : v.norm_fm, r.value);
        ++failures;
      }
    } else if (r.kind == "delta") {
      for (double a : {4.0, 5.35, 6.5})
        for (double b : {0.0, -1.3}) {
          const double f = RegularSolution(p, 2, 6, kMu_MeV, r.l, 0.0, a, r.E).log_derivative;
          const ChannelTriple t = InteriorChannelFunctions(Outer(p, r.l, a, r.E), f, b);
          const double delta = Wrap(std::arg(t.exp_minus_i_phi));
          const double d = std::fabs(Wrap(delta - r.value));
          worstD = std::max(worstD, d);
          if (d > 4e-6) {
            std::printf("FAIL: phase l = %d E = %g a = %g b = %g: %.9f against %.9f rad\n", r.l, r.E, a, b, delta, r.value);
            ++failures;
          }
        }
    } else if (r.kind == "pole") {
      double lo = 1.0, hi = 3.0;  // f falls through 0 at the pole (df/dE < 0)
      for (int i = 0; i < 60; i++) {
        const double mid = 0.5 * (lo + hi);
        (RegularSolution(p, 2, 6, kMu_MeV, r.l, 0.0, r.a, mid).log_derivative > 0.0 ? lo : hi) = mid;
      }
      poleE = 0.5 * (lo + hi);
      if (std::fabs(poleE - r.value) > 2e-6) {
        std::printf("FAIL: pole of R_pot at %.7f MeV against %.7f\n", poleE, r.value);
        ++failures;
      }
    } else if (r.kind == "gamma2") {
      const RegularValue v = RegularSolution(p, 2, 6, kMu_MeV, r.l, 0.0, r.a, poleE);
      const double gamma2 = kHbarC_MeV_fm * kHbarC_MeV_fm / (2.0 * kMu_MeV * r.a * v.norm_fm);
      const double theta2 = gamma2 / (3.0 * kHbarC_MeV_fm * kHbarC_MeV_fm / (2.0 * kMu_MeV * r.a * r.a));
      std::printf("pole of R_pot (a = %.2f fm, b = 0): E = %.5f MeV, gamma^2 = %.5f MeV, theta^2 = %.3f\n", r.a, poleE,
                  gamma2, theta2);
      if (std::fabs(gamma2 / r.value - 1.0) > 1e-5) {
        std::printf("FAIL: gamma^2 of the pole %.8f against %.8f MeV\n", gamma2, r.value);
        ++failures;
      }
    }
  }
  std::printf("against the DOP853 reference: |d atan f| %.1e, |d norm| %.1e relative, |d phase| %.1e rad\n", worstF, worstN,
              worstD);
}

// The error falls as step^4: two rows of the reference with a quarter of the default step.
void StepOrder(const char *path) {
  HybridPotential p = Ktf();
  p.step_fm = 0.0275 / 4.0;
  for (const RefRow &r : ReadReference(path)) {
    const bool rowF = r.kind == "f" && r.l == 1 && r.a == 6.5 && r.E == 2.0;
    const bool rowN = r.kind == "norm" && r.l == 2 && r.a == 4.0 && r.E == 5.0;
    if (!rowF && !rowN) continue;
    const RegularValue v = RegularSolution(p, 2, 6, kMu_MeV, r.l, 0.0, r.a, r.E);
    const double d = rowF ? std::fabs(v.log_derivative - r.value) : std::fabs(v.norm_fm / r.value - 1.0);
    std::printf("a quarter of the step, %s: %.1e\n", rowF ? "|df|" : "|d norm| relative", d);
    Check(d < 1e-7, "a quarter of the step brings the regular solution within 1e-7 of the reference");
  }
}

void OneLevel() {
  const HybridPotential p = Ktf();
  const double a = 5.35, b = 0.0, El = 3.2, gamma2 = 0.15;
  double worst = 0.0;
  for (double E : {1.0, 1.926, 2.5, 3.1, 4.0}) {
    const ChannelTriple o = Outer(p, 1, a, E);
    const double f = RegularSolution(p, 2, 6, kMu_MeV, 1, 0.0, a, E).log_derivative;
    const std::complex<double> L0(o.S - b, o.P);
    const double R = 1.0 / (f - b) + gamma2 / (El - E);
    const std::complex<double> direct = o.exp_minus_i_phi * o.exp_minus_i_phi * (1.0 - std::conj(L0) * R) / (1.0 - L0 * R);
    const ChannelTriple t = InteriorChannelFunctions(o, f, b);
    const std::complex<double> Ainv = (El - E) - gamma2 * std::complex<double>(t.S - b, t.P);
    const std::complex<double> primed =
        t.exp_minus_i_phi * t.exp_minus_i_phi * (1.0 + std::complex<double>(0.0, 2.0) * t.P * gamma2 / Ainv);
    worst = std::max(worst, std::abs(direct - primed));
    Check(t.P >= 0.0, "P' >= 0");
  }
  std::printf("one level and R_pot, direct against primed: worst %.1e\n", worst);
  Check(worst < 1e-12, "U from R_pot + a pole equals the level-matrix form with the primed functions");
}

void Limits() {
  ChannelTriple o;
  o.S = -0.7;
  o.P = 0.3;
  o.exp_minus_i_phi = std::polar(1.0, -0.4);
  const ChannelTriple same = InteriorChannelFunctions(o, std::numeric_limits<double>::infinity(), 0.0);
  Check(same.S == o.S && same.P == o.P && same.exp_minus_i_phi == o.exp_minus_i_phi, "f = infinity: the outer functions");
  const ChannelTriple far = InteriorChannelFunctions(o, 1e12, 0.0);
  Check(std::fabs(far.S - o.S) < 1e-10 && std::fabs(far.P - o.P) < 1e-10 && std::abs(far.exp_minus_i_phi - o.exp_minus_i_phi) < 1e-10,
        "f -> infinity: the outer functions");
  const ChannelTriple pole = InteriorChannelFunctions(o, 0.25, 0.25);
  Check(pole.P == 0.0 && pole.S == 0.25, "at a pole of R_pot (f = b): P' = 0, S' = b");
  const ChannelTriple below = InteriorChannelFunctions(o, 0.2499, 0.25), above = InteriorChannelFunctions(o, 0.2501, 0.25);
  Check(std::abs(below.exp_minus_i_phi + above.exp_minus_i_phi) < 1e-3,
        "exp(-i phi') changes sign through a pole of R_pot (with P'^1/2 the amplitude is continuous)");
}

void Wave() {
  // neutral particles in a Woods-Saxon that is flat to 1e-39 of its depth over [0, 5 fm]: u = sin(K r)
  const HybridPotential p = WoodsSaxonPotential(30.0, 50.0, 0.5);
  const double mu = 3000.0, a = 5.0, k2 = 2.0 * mu / (kHbarC_MeV_fm * kHbarC_MeV_fm);
  const RegularWave u1 = RegularWaveOnGrid(p, 0, 0, mu, 0, 0.0, a, 1.0), u2 = RegularWaveOnGrid(p, 0, 0, mu, 0, 0.0, a, -4.0);
  const double K1 = std::sqrt(k2 * 31.0), K2 = std::sqrt(k2 * 26.0);
  const double exact = (std::sin((K1 - K2) * a) / (2.0 * (K1 - K2)) - std::sin((K1 + K2) * a) / (2.0 * (K1 + K2))) /
                       (std::sin(K1 * a) * std::sin(K2 * a));
  const double got = InteriorOverlap(u1, u2, [](double) { return 1.0; });
  const int n = static_cast<int>(u1.u.size()) - 1;
  const double fine = InteriorOverlap(RegularWaveOnGrid(p, 0, 0, mu, 0, 0.0, a, 1.0, 4 * n),
                                      RegularWaveOnGrid(p, 0, 0, mu, 0, 0.0, a, -4.0, 4 * n), [](double) { return 1.0; });
  std::printf("overlap of two flat-well waves: %.3e off at the default step, %.3e at a quarter of it\n", got / exact - 1.0,
              fine / exact - 1.0);
  Check(std::fabs(got / exact - 1.0) < 1e-5 && std::fabs(fine / exact - 1.0) < 1e-7,
        "the overlap of two waves of a flat well against the closed form (1e-5; 1e-7 at a quarter of the step)");
  const RegularValue v = RegularSolution(p, 0, 0, mu, 0, 0.0, a, 1.0);
  Check(v.log_derivative == u1.log_derivative && v.norm_fm == u1.norm_fm && u1.u.back() == 1.0 && u1.u.front() == 0.0,
        "RegularSolution is the wave's f and norm");
  Check(std::fabs(InteriorOverlap(u1, u1, [](double) { return 1.0; }) / v.norm_fm - 1.0) < 1e-14, "the overlap of a wave with itself is its norm");
  const RegularWave coarse = RegularWaveOnGrid(p, 0, 0, mu, 0, 0.0, a, 1.0, 100);
  Check(coarse.u.size() == 101 && std::fabs(coarse.log_derivative / u1.log_derivative - 1.0) < 1e-4, "a given number of intervals");
  bool threw = false;
  try {
    InteriorOverlap(u1, coarse, [](double) { return 1.0; });
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  Check(threw, "two waves on different grids are refused");
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("FAIL: usage: physics_interior_potential_test data/interior_potential.csv\n");
    return 1;
  }
  AgainstReference(argv[1]);
  StepOrder(argv[1]);
  OneLevel();
  Limits();
  Wave();
  if (failures) {
    std::printf("%d failure(s)\n", failures);
    return 1;
  }
  std::printf("interior_potential: all checks passed\n");
  return 0;
}
