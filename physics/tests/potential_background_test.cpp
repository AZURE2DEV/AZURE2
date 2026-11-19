// potential_background (docs/upgrade/potential_background.md sections 1-3, 7.2): the interior R-function of a
// potential and its use as the background of an R-matrix channel.
//
//  1. Zero potential, neutral: l(E) = x u_l'(x)/u_l(x), x = k a, u_l the Riccati-Bessel function (l = 0: R = tan x / x
//     at B = 0), l = 0-3, a = 3, 5 fm, E = 0.5, 4 MeV, and below threshold (E = -2 MeV, u_0 = sinh); to 1e-7 relative.
//  2. A constant well over the whole interior (Woods-Saxon with R = a + 4 fm, diffuseness 0.05 fm: V(a) = -V0 to
//     1e-30): l = K a u_l'(Ka)/u_l(Ka), K^2 = 2 mu (E + V0)/hbar^2, l = 0-3, V0 = 40 MeV; to 2e-6 relative (the
//     O(step^4) error of Numerov at K h = 0.05).
//  3. Pure point Coulomb interior (Z1 Z2 = 4, 3He + 4He, R_C = 0): l = rho F_l'/F_l against mpmath's coulombf
//     (data/potential_background.csv, make_potential_background_reference.py); to 2e-6.
//  4. Box states of the zero potential, l = 0, B = 0: E_p = hbar^2 ((n + 1/2) pi / a)^2 / 2mu, gamma_p^2 = hbar^2/(mu a^2),
//     n nodes; and for the deep Gaussian of part 5 the residue -1/(dl/dE) against (hbar^2/2 mu a) u(a)^2 / int u^2
//     (Simpson on the Numerov grid), to 1e-5.  R_pot - gamma_p^2/(E_p - E) is regular at E_p (subtraction).
//  5. The identity (section 2.3), the heart: for the Gaussian cluster potential of 3He + 4He (Dubovichenko's depths per
//     l, R_C = 3.095 fm), the collision matrix of the channel with no levels, from R_pot at a_c and the exterior
//     functions of the same potential (IntegrateInward), U0 e^(-2 i omega) = e^(-2 i phi) (l - L*)/(l - L), equals
//     e^(2 i delta) with delta from scipy's DOP853 integration to 15 fm (data/potential_background.csv), for
//     a_c = 3, 4, 5, 6 fm, l = 0-3, E = 0.5-4 MeV, B_bg = 0 and -2; to 2e-6 in |U0 - e^(2i delta)|.  The same through
//     DressWaves (Eq. (9)): e^(-2 i phi~) of the dressed waves.  With the Coulomb exterior instead of the potential's
//     the identity fails (a check that the test can fail): |dU| > 1e-3 for every l <= 2, E >= 1 MeV.
//  6. Below threshold the background acts through S~ = B + (S - B)/(1 - R0 (S - B)) (ChannelBoundShift), continuous
//     with the dressed open-channel S~ through E = 0 (5e-3 over 1 keV).
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

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

const double kM1 = 3.0160293, kM2 = 4.0026032;
const double kMu_MeV = kM1 * kM2 / (kM1 + kM2) * kAtomicMassUnit_MeV;

HybridPotential ZeroPotential() {
  HybridPotential p = WoodsSaxonPotential(1.0, 3.0, 0.6);
  p.central_depth_MeV = 0.0;
  p.coulomb_radius_fm = 0.0;
  return p;
}

// Riccati-Bessel u_l(x) = x j_l(x) and its derivative, l = 0-3, closed forms.
void Riccati(int l, double x, double &u, double &du) {
  const double s = std::sin(x), c = std::cos(x);
  switch (l) {
    case 0: u = s; du = c; break;
    case 1: u = s / x - c; du = c / x - s / (x * x) + s; break;
    case 2: u = (3.0 / (x * x) - 1.0) * s - 3.0 * c / x;
            du = -6.0 / (x * x * x) * s + (3.0 / (x * x) - 1.0) * c + 3.0 * s / x + 3.0 * c / (x * x); break;
    default: u = (15.0 / (x * x * x) - 6.0 / x) * s - (15.0 / (x * x) - 1.0) * c;
             du = (-45.0 / (x * x * x * x) + 6.0 / (x * x)) * s + (15.0 / (x * x * x) - 6.0 / x) * c +
                  30.0 / (x * x * x) * c + (15.0 / (x * x) - 1.0) * s;
  }
}

double FreeParticle() {
  double worst = 0.0;
  const HybridPotential p = ZeroPotential();
  for (int l = 0; l <= 3; ++l)
    for (double a : {3.0, 5.0})
      for (double E : {0.5, 4.0}) {
        const double k = std::sqrt(2.0 * kMu_MeV * E) / kHbarC_MeV_fm, x = k * a;
        double u, du;
        Riccati(l, x, u, du);
        const double ref = x * du / u, got = InteriorLogDerivative(p, 0, 0, kMu_MeV, l, 0.0, a, E);
        const double d = std::fabs(got / ref - 1.0);
        worst = std::max(worst, d);
        Check(d < 1e-7, "free particle l = " + std::to_string(l) + ": " + std::to_string(got) + " vs " + std::to_string(ref));
        if (l == 0) {
          const double R = InteriorRFunction(p, 0, 0, kMu_MeV, 0, 0.0, a, E, 0.0);
          Check(std::fabs(R / (std::tan(x) / x) - 1.0) < 1e-7, "free particle l = 0: R = tan(ka)/(ka)");
        }
      }
  // below threshold: u_0 = sinh(kappa r), l = kappa a coth(kappa a)
  const double kappa = std::sqrt(2.0 * kMu_MeV * 2.0) / kHbarC_MeV_fm, xa = kappa * 4.0;
  const double got = InteriorLogDerivative(p, 0, 0, kMu_MeV, 0, 0.0, 4.0, -2.0);
  const double d = std::fabs(got / (xa / std::tanh(xa)) - 1.0);
  worst = std::max(worst, d);
  Check(d < 1e-7, "free particle below threshold");
  return worst;
}

double ConstantWell() {
  double worst = 0.0;
  const double V0 = 40.0, a = 4.0;
  HybridPotential p = WoodsSaxonPotential(V0, a + 4.0, 0.05);
  p.coulomb_radius_fm = 0.0;
  for (int l = 0; l <= 3; ++l)
    for (double E : {-5.0, 0.5, 4.0}) {
      const double K = std::sqrt(2.0 * kMu_MeV * (E + V0)) / kHbarC_MeV_fm, x = K * a;
      double u, du;
      Riccati(l, x, u, du);
      const double ref = x * du / u, got = InteriorLogDerivative(p, 0, 0, kMu_MeV, l, 0.0, a, E);
      const double d = std::fabs(got / ref - 1.0);
      worst = std::max(worst, d);
      Check(d < 2e-6, "constant well l = " + std::to_string(l) + " E = " + std::to_string(E));
    }
  return worst;
}

struct Ref {
  std::map<std::tuple<int, double, double>, double> ellC;  // (l, a, E)
  std::map<std::pair<int, double>, double> delta;          // (l, E)
};

Ref ReadReference(const char *path) {
  Ref r;
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string kind, f;
    std::getline(ss, kind, ',');
    std::getline(ss, f, ',');
    const int l = std::stoi(f);
    std::getline(ss, f, ',');
    const double a = std::stod(f);
    std::getline(ss, f, ',');
    const double E = std::stod(f);
    std::getline(ss, f, ',');
    const double v = std::stod(f);
    if (kind == "ellC") r.ellC[{l, a, E}] = v;
    if (kind == "delta") r.delta[{l, E}] = v;
  }
  return r;
}

double PointCoulomb(const Ref &ref) {
  double worst = 0.0;
  const HybridPotential p = ZeroPotential();
  for (const auto &e : ref.ellC) {
    const int l = std::get<0>(e.first);
    const double a = std::get<1>(e.first), E = std::get<2>(e.first);
    const double got = InteriorLogDerivative(p, 2, 2, kMu_MeV, l, 0.0, a, E);
    const double d = std::fabs(got - e.second);
    worst = std::max(worst, d);
    Check(d < 2e-6, "point Coulomb l = " + std::to_string(l) + " a = " + std::to_string(a) + " E = " + std::to_string(E));
  }
  return worst;
}

// The Gaussian cluster potential of 3He + 4He with the depth of partial wave l (the reference script's).
HybridPotential Cluster(int l) {
  const double depth[] = {67.5, 83.83, 69.0, 84.8};
  HybridPotential p = GaussianPotential(depth[l], 1.0 / std::sqrt(0.15747));
  p.coulomb_radius_fm = 3.095;
  return p;
}

double BoxStates() {
  double worst = 0.0;
  // zero potential, l = 0, B = 0
  const HybridPotential zero = ZeroPotential();
  const double a = 4.0, h2m = kHbarC_MeV_fm * kHbarC_MeV_fm / (2.0 * kMu_MeV);
  const std::vector<BoxState> s = InteriorBoxStates(zero, 0, 0, kMu_MeV, 0, 0.0, a, 0.0, 0.0, 60.0);
  Check(s.size() >= 3, "zero potential: at least three box states below 60 MeV");
  for (const BoxState &b : s) {
    const double kn = (b.nodes + 0.5) * M_PI / a;
    const double dE = std::fabs(b.E_MeV / (h2m * kn * kn) - 1.0), dg = std::fabs(b.gamma2_MeV / (2.0 * h2m / (a * a)) - 1.0);
    worst = std::max(worst, std::max(dE, dg));
    Check(dE < 1e-7 && dg < 1e-6, "zero potential box state " + std::to_string(b.nodes));
  }
  // deep Gaussian, l = 0: residue by the derivative against the norm integral (Simpson on a fine grid)
  for (double ac : {3.5, 4.5}) {
    const HybridPotential p = Cluster(0);
    const std::vector<BoxState> st = InteriorBoxStates(p, 2, 2, kMu_MeV, 0, 0.0, ac, 0.0, -80.0, 30.0);
    Check(st.size() >= 3, "Gaussian s wave: at least three box states (two forbidden) below 30 MeV");
    for (std::size_t i = 0; i < st.size(); ++i) Check(st[i].nodes == (int)i, "box states ordered by nodes");
    for (const BoxState &b : st) {
      const int n = 2000;
      double norm = 0.0;
      for (int i = 0; i <= n; ++i) {
        const double r = ac * i / n;
        const InteriorValue v = IntegrateOutward(p, 2, 2, kMu_MeV, 0, 0.0, std::max(r, 1e-3), b.E_MeV);
        const double u = r > 0 ? v.u * std::exp(v.log_scale) : 0.0;
        norm += u * u * ((i == 0 || i == n) ? 1.0 : (i % 2 ? 4.0 : 2.0));
      }
      norm *= ac / n / 3.0;
      const InteriorValue va = IntegrateOutward(p, 2, 2, kMu_MeV, 0, 0.0, ac, b.E_MeV);
      const double ua = va.u * std::exp(va.log_scale);
      const double g2 = h2m / ac * ua * ua / norm;
      const double d = std::fabs(g2 / b.gamma2_MeV - 1.0);
      worst = std::max(worst, d);
      Check(d < 1e-4, "Gaussian residue against the norm integral, nodes " + std::to_string(b.nodes));
      // subtraction leaves R regular at E_p
      InteriorBackground bg;
      bg.subtract_nodes = {b.nodes};
      const double Rm = BackgroundR(p, 2, 2, kMu_MeV, 0, 0.0, ac, b.E_MeV - 1e-3, bg);
      const double Rp = BackgroundR(p, 2, 2, kMu_MeV, 0, 0.0, ac, b.E_MeV + 1e-3, bg);
      Check(std::fabs(Rp - Rm) < 1e-2 * (1.0 + std::fabs(Rp)), "subtracted R is regular at E_p");
    }
  }
  return worst;
}

// U0 e^(-2 i omega) of a channel without levels from the interior R-function and exterior functions.
struct Exterior {
  std::complex<double> H, dH;  // H at a and its rho-derivative
  double rho;
};

Exterior ExteriorOf(const HybridPotential &p, bool hybrid, int l, double a, double E) {
  const double k = std::sqrt(2.0 * kMu_MeV * E) / kHbarC_MeV_fm;
  const double eta = 4.0 * kFineStructure * std::sqrt(kMu_MeV / (2.0 * E));
  Exterior x;
  x.rho = k * a;
  if (!hybrid) {
    const coulomb::CoulombWaves w = coulomb::coulomb(l, eta, x.rho);
    x.H = {w.G, w.F};
    x.dH = {w.Gp, w.Fp};
    return x;
  }
  auto outgoing = [&](double r) {
    const coulomb::CoulombWaves w = coulomb::coulomb(l, eta, k * r);
    return std::complex<double>(w.G, w.F);
  };
  const RadialValue h = IntegrateInward(p, 2, 2, kMu_MeV, l, 0.0, a, E, outgoing);
  x.H = h.u;
  x.dH = h.du_dr / k;
  return x;
}

double Identity(const Ref &ref, double &worstDressed, double &smallestCoulombMiss) {
  double worst = 0.0;
  worstDressed = 0.0;
  smallestCoulombMiss = 1e300;
  for (const auto &e : ref.delta) {
    const int l = e.first.first;
    const double E = e.first.second;
    const std::complex<double> truth = std::exp(std::complex<double>(0.0, 2.0 * e.second));
    const HybridPotential p = Cluster(l);
    for (double a : {3.0, 4.0, 5.0, 6.0})
      for (double B : {0.0, -2.0}) {
        const double R0 = InteriorRFunction(p, 2, 2, kMu_MeV, l, 0.0, a, E, B);
        const Exterior x = ExteriorOf(p, true, l, a, E);
        const std::complex<double> L = x.rho * x.dH / x.H, L0 = L - B;
        const std::complex<double> hsp = std::conj(x.H) / std::abs(x.H);  // e^(-i phi)
        const std::complex<double> U = hsp * hsp * (1.0 - R0 * std::conj(L0)) / (1.0 - R0 * L0);
        const double d = std::abs(U - truth);
        worst = std::max(worst, d);
        Check(d < 2e-6, "identity l = " + std::to_string(l) + " a = " + std::to_string(a) + " E = " + std::to_string(E) +
                            " B = " + std::to_string(B) + ": |dU| = " + std::to_string(d));
        // through the dressed waves
        WavesAtRadius w{x.H.imag(), x.dH.imag(), x.H.real(), x.dH.real()};
        const WavesAtRadius dw = DressWaves(w, x.rho, R0, B);
        const std::complex<double> Hd(dw.G, dw.F);
        const std::complex<double> hspd = std::conj(Hd) / std::abs(Hd);
        const double dd = std::abs(hspd * hspd - truth);
        worstDressed = std::max(worstDressed, dd);
        Check(dd < 2e-6, "identity through DressWaves l = " + std::to_string(l));
        // the Wronskian is kept: F~ G~' - G~ F~' = F G' - G F'
        Check(std::fabs((dw.F * dw.Gp - dw.G * dw.Fp) - (w.F * w.Gp - w.G * w.Fp)) < 1e-10 * (1.0 + std::fabs(w.F * w.Gp)),
              "DressWaves keeps the Wronskian");
        // the Coulomb exterior instead: the identity must fail
        if (l <= 2 && E >= 1.0) {
          const Exterior c = ExteriorOf(p, false, l, a, E);
          const std::complex<double> Lc = c.rho * c.dH / c.H - B;
          const std::complex<double> hc = std::conj(c.H) / std::abs(c.H);
          const std::complex<double> Uc = hc * hc * (1.0 - R0 * std::conj(Lc)) / (1.0 - R0 * Lc);
          smallestCoulombMiss = std::min(smallestCoulombMiss, std::abs(Uc - truth));
        }
      }
  }
  Check(smallestCoulombMiss > 1e-3, "with the Coulomb exterior the identity fails (the test can fail)");
  return worst;
}

// 6. Below threshold: the dressed S (ChannelBoundShift, DressShift) continues the open-channel S~ of the dressed waves
// through E = 0: |S~(+1 keV) - S~(-1 keV)| small (dS/dE ~ 1 MeV^-1), for the s wave of the Gaussian with its background.
double Threshold() {
  HybridPotential p = Cluster(0);
  PartialWave w;
  w.l = 0;
  w.background = InteriorBackground{-0.5, {}};
  p.waves.push_back(w);
  const double a = 4.2, E = 1e-3;
  const double R0 = BackgroundR(p, 2, 2, kMu_MeV, 0, 0.0, a, E, *w.background);
  const Exterior x = ExteriorOf(p, true, 0, a, E);
  WavesAtRadius wv{x.H.imag(), x.dH.imag(), x.H.real(), x.dH.real()};
  const WavesAtRadius d = DressWaves(wv, x.rho, R0, -0.5);
  const double Sabove = x.rho * (d.F * d.Fp + d.G * d.Gp) / (d.F * d.F + d.G * d.G);
  const double Sbelow = ChannelBoundShift(p, 2, 2, kMu_MeV, 0, -1, 0.0, a, E);
  const double S0 = ChannelBoundShift(p, 2, 2, kMu_MeV, 0, -1, 0.0, a, 0.0);
  const double dd = std::max(std::fabs(Sabove - Sbelow), std::fabs(Sabove - S0));
  Check(dd < 5e-3, "dressed S continuous through threshold: " + std::to_string(Sabove) + " / " + std::to_string(S0) +
                       " / " + std::to_string(Sbelow));
  // DressShift is DressWaves on a real wave: S~ of the dressed (W, rho W') equals the formula
  const double S = 1.3, Rr = 0.4, B = -0.5;
  const WavesAtRadius r = DressWaves({1.0, S / 2.0, 0.0, 0.0}, 2.0, Rr, B);
  Check(std::fabs(2.0 * r.Fp / r.F - DressShift(S, Rr, B)) < 1e-14, "DressShift = DressWaves on a real wave");
  return dd;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: physics_potential_background_test <data/potential_background.csv>\n");
    return 2;
  }
  const Ref ref = ReadReference(argv[1]);
  Check(ref.ellC.size() == 16 && ref.delta.size() == 16, "reference file read");
  const double free = FreeParticle();
  const double well = ConstantWell();
  const double coul = PointCoulomb(ref);
  const double box = BoxStates();
  double dressed = 0.0, miss = 0.0;
  const double ident = Identity(ref, dressed, miss);
  const double thr = Threshold();
  std::printf("free particle %.2e, constant well %.2e, point Coulomb %.2e, box states %.2e\n", free, well, coul, box);
  std::printf("identity |U0 - exp(2 i delta)|: %.2e (formula 7), %.2e (dressed waves); Coulomb exterior misses by >= %.2e\n",
              ident, dressed, miss);
  std::printf("threshold: |S~(+1 keV) - S~(-1 keV, 0)| = %.2e\n", thr);
  if (failures) {
    std::printf("%d failure(s)\n", failures);
    return 1;
  }
  std::printf("PASS\n");
  return 0;
}
