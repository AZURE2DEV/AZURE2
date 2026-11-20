// potential_capture (docs/upgrade/potential_background.md section 5): the interior direct-capture kernel of a
// background partial wave and the bound-state rule.
//
//  1. Zero potential, neutral pair: the regular solutions are closed forms, phi_E = sin(kr)/sin(ka) (l = 0) and
//     phi_f = i1(kappa r)/i1(kappa a), i1(x) = cosh x - sinh x / x (l = 1, E_f < 0), R_pot = 1/(ka cot ka - B), so
//     K = R_pot int_0^a phi_f r phi_E dr by a 20000-panel Simpson sum of the closed forms; to 1e-6 relative (Numerov's
//     O(h^4)) for a = 4.0 fm (an even number of intervals) and a = 4.87 fm (odd: the 3/8 rule), B = 0 and -1.
//  2. The same with the lowest box state subtracted (subtract 0): its pole term gamma_p^2 phi_p/(E_p - E) in closed
//     form, k_p a = the root of ka cot ka = B, phi_p = sin(k_p r)/sin(k_p a), gamma_p^2 = -1/(dl/dE); to 1e-5.
//  3. The bound-state rule on the Gaussian potential VD of 3He + 4He (p3/2, 83.589554 MeV, alpha = 0.15747 fm^-2,
//     R_C = 3.095 fm): NearestEigenvalue at a_c = 4 fm equals the binding energy found by scipy's DOP853 shooting of
//     tests/potential_capture/make_reference.py, B = 1.4548867432 MeV, to 1e-6 MeV; DepthBindingAt that energy returns
//     the written depth to 1e-6 relative, and at E = -1.5866 MeV a deeper one; the eigenvalue does not depend on a_c
//     (3.5 and 5.5 fm, 1e-6 MeV).
//  The potential's own ANC (InteriorNormIntegral with the exterior norm) is checked through the engine's set-up note
//  against make_reference.py's 4.0931604450 fm^-1/2 in tests/potential_capture/check.sh.
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/HybridPotential.h"

using namespace azure2::physics;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) ++failures;
}

const double kM1 = 3.0160293, kM2 = 4.0026032;
const double kMu_MeV = kM1 * kM2 / (kM1 + kM2) * kAtomicMassUnit_MeV;

HybridPotential ZeroPotential() {
  HybridPotential p = WoodsSaxonPotential(1.0, 3.0, 0.6);
  p.central_depth_MeV = 0.0;
  p.coulomb_radius_fm = 0.0;
  return p;
}

double Simpson(const std::function<double(double)> &f, double a, double b, int n) {
  const double h = (b - a) / n;
  double s = f(a) + f(b);
  for (int i = 1; i < n; ++i) s += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
  return s * h / 3.0;
}

double Kinetic(double E) { return std::sqrt(2.0 * kMu_MeV * std::fabs(E)) / kHbarC_MeV_fm; }

/// The closed-form kernel of the zero potential (part 1, 2): s wave at E > 0 into the l = 1 regular solution at Ef < 0.
double ClosedKernel(double a, double E, double Ef, double B, bool subtractLowest) {
  const double k = Kinetic(E), kap = Kinetic(Ef);
  auto i1 = [](double x) { return x < 1e-3 ? x * x / 3.0 : std::cosh(x) - std::sinh(x) / x; };
  const double Rpot = 1.0 / (k * a / std::tan(k * a) - B);
  std::function<double(double)> g = [&](double r) { return Rpot * std::sin(k * r) / std::sin(k * a); };
  double kp = 0.0, Ep = 0.0, g2 = 0.0;
  if (subtractLowest) {
    // x cot x = B: the lowest root in (0, pi) for B < 1, by bisection
    double lo = 1e-9, hi = M_PI - 1e-9;
    auto f = [&](double x) { return x / std::tan(x) - B; };
    for (int it = 0; it < 200; ++it) {
      const double m = 0.5 * (lo + hi);
      if ((f(m) > 0.0) == (f(lo) > 0.0)) lo = m; else hi = m;
    }
    kp = 0.5 * (lo + hi) / a;
    Ep = std::pow(kHbarC_MeV_fm * kp, 2.0) / (2.0 * kMu_MeV);
    // gamma^2 = (hbar^2/2 mu a) u(a)^2 / int_0^a u^2, u = sin(k_p r)  (Eq. (4) of the note)
    const double norm = Simpson([&](double r) { return std::pow(std::sin(kp * r), 2.0); }, 0.0, a, 20000);
    g2 = std::pow(kHbarC_MeV_fm, 2.0) / (2.0 * kMu_MeV * a) * std::pow(std::sin(kp * a), 2.0) / norm;
    g = [&, Rpot](double r) {
      return Rpot * std::sin(k * r) / std::sin(k * a) - g2 * std::sin(kp * r) / std::sin(kp * a) / (Ep - E);
    };
  }
  return Simpson([&](double r) { return i1(kap * r) / i1(kap * a) * r * g(r); }, 0.0, a, 20000);
}

void ZeroPotentialKernel() {
  const HybridPotential p = ZeroPotential();
  const auto weight = [](double r) { return r; };
  for (bool sub : {false, true})
    for (double a : {4.0, 4.87})
      for (double B : {0.0, -1.0})
        for (double E : {0.3, 2.0}) {
          InteriorBackground bg;
          bg.boundary = B;
          bg.capture = true;
          if (sub) bg.subtract_nodes = {0};
          const double Ef = -1.5;
          const double K = InteriorCaptureIntegral(p, 0, 0.0, E, bg, p, 1, 0.0, Ef, 0, 0, kMu_MeV, a, weight);
          const double ref = ClosedKernel(a, E, Ef, B, sub);
          char text[200];
          std::snprintf(text, sizeof text, "zero potential%s, a = %.2f fm, B = %g, E = %g MeV: K = %.9g, closed form %.9g (rel %.2e)",
                        sub ? " minus the lowest box pole" : "", a, B, E, K, ref, std::fabs(K / ref - 1.0));
          Check(std::fabs(K / ref - 1.0) < (sub ? 1e-5 : 1e-6), text);
        }
}

HybridPotential VD_p32() {
  HybridPotential p = GaussianPotential(83.589554, 1.0 / std::sqrt(0.15747));
  p.coulomb_radius_fm = 3.095;
  p.matching_radius_fm = 15.0;
  p.step_fm = 0.0275;
  return p;
}

void BoundStateRule() {
  const HybridPotential p = VD_p32();
  const double ldoti = 1.0;  // <2 l.i> of p3/2 with i = 1/2 (TwoLDotI(1, 3, 1))
  Check(std::fabs(TwoLDotI(1, 3, 1) - ldoti) < 1e-15, "<2 l.i> = 1 for p3/2");
  const double scipyB = 1.4548867432;  // make_reference.py (DOP853 shooting, mpmath Whittaker)
  for (double a : {3.5, 4.0, 5.5}) {
    const PotentialEigenvalue ev = NearestEigenvalue(p, 2, 2, kMu_MeV, 1, 0.0, a, -1.5866, 1.0);
    char text[200];
    std::snprintf(text, sizeof text, "VD p3/2 eigenvalue at a_c = %.1f fm: %.10f MeV, scipy -%.10f (diff %.1e MeV)", a,
                  ev.E_MeV, scipyB, ev.E_MeV + scipyB);
    Check(ev.found && std::fabs(ev.E_MeV + scipyB) < 1e-6, text);
  }
  const DepthAtEnergy d = DepthBindingAt(p, 2, 2, kMu_MeV, 1, 0.0, 4.0, -scipyB);
  char text[200];
  std::snprintf(text, sizeof text, "the depth binding at -%.10f MeV: %.9f MeV (written 83.589554)", scipyB, d.depth_MeV);
  Check(d.found && std::fabs(d.depth_MeV / 83.589554 - 1.0) < 1e-6, text);
  const DepthAtEnergy deeper = DepthBindingAt(p, 2, 2, kMu_MeV, 1, 0.0, 4.0, -1.5866);
  std::snprintf(text, sizeof text, "the depth binding at -1.5866 MeV: %.6f MeV (deeper)", deeper.depth_MeV);
  Check(deeper.found && deeper.depth_MeV > 83.589554 && deeper.depth_MeV < 1.05 * 83.589554, text);
  // the matching function vanishes at the eigenvalue and changes sign through it
  const double f1 = BoundMatching(p, 2, 2, kMu_MeV, 1, 0.0, 4.0, -scipyB - 0.01);
  const double f2 = BoundMatching(p, 2, 2, kMu_MeV, 1, 0.0, 4.0, -scipyB + 0.01);
  Check((f1 < 0.0) != (f2 < 0.0), "the matching function changes sign through the eigenvalue");
}

// The poles of R0 (review 2026-11-20, findings 5 and 6), on VD's s wave at a = 4.87 fm, B = 0 (a box pole at 7.59 MeV):
//  - exactly at the bisected pole (to the last bit) R_pot, the dressed waves and the closed-channel factor are finite,
//    and R_pot keeps the sign and size of the pole's neighbourhood;
//  - with the pole subtracted, R0 at E_p and within the 1 keV window is finite and on the line through its values at
//    the window's edges, and continuous with R0 outside (1e-6 MeV and 2e-3 MeV away: the same to 1e-6 of R0);
//  - the kernel of the interior capture: the subtracted part is gamma_p^2 I_p/(E_p - E) with a constant I_p;
//  - the node count includes the last interval: the states found are numbered 0, 1, 2, ... without a gap.
void Poles() {
  HybridPotential p = GaussianPotential(67.5, 1.0 / std::sqrt(0.15747));
  p.coulomb_radius_fm = 3.095;
  p.matching_radius_fm = 15.0;
  p.step_fm = 0.0275;
  const double a = 4.87;
  const std::vector<BoxState> st = InteriorBoxStates(p, 2, 2, kMu_MeV, 0, 0.0, a, 0.0, -60.0, 9.0);
  bool numbered = !st.empty();
  for (size_t i = 0; i < st.size(); ++i) numbered = numbered && st[i].nodes == (int)i;
  Check(numbered, "box states of VD's s wave at 4.87 fm numbered 0, 1, ... (" + std::to_string(st.size()) + " states)");
  if (st.empty()) return;
  const BoxState &top = st.back();
  char text[240];
  std::snprintf(text, sizeof text, "the highest box pole below 9 MeV at E_p = %.15g MeV (expected near 7.5915)", top.E_MeV);
  Check(std::fabs(top.E_MeV - 7.5915) < 1e-3, text);
  const double Rp = InteriorRFunction(p, 2, 2, kMu_MeV, 0, 0.0, a, top.E_MeV, 0.0);
  const double Rn = InteriorRFunction(p, 2, 2, kMu_MeV, 0, 0.0, a, std::nextafter(top.E_MeV, 0.0), 0.0);
  std::snprintf(text, sizeof text, "R_pot at the bisected pole finite: %.3e (at the energy below: %.3e)", Rp, Rn);
  Check(std::isfinite(Rp) && std::isfinite(Rn) && std::fabs(Rp) > 1e8, text);
  WavesAtRadius w;
  w.F = 0.3;
  w.Fp = 0.9;
  w.G = 2.0;
  w.Gp = -1.1;
  const WavesAtRadius d = DressWaves(w, 1.5, Rp, 0.0);
  Check(std::isfinite(d.F) && std::isfinite(d.G) && std::isfinite(d.Fp) && std::isfinite(d.Gp),
        "the dressed waves at the bisected pole are finite");
  // subtracted: the window
  InteriorBackground bg;
  bg.subtract_nodes = {top.nodes};
  const double r0 = BackgroundR(p, 2, 2, kMu_MeV, 0, 0.0, a, top.E_MeV, bg);
  const double rl = BackgroundR(p, 2, 2, kMu_MeV, 0, 0.0, a, top.E_MeV - kSubtractedPoleWindow_MeV, bg);
  const double rh = BackgroundR(p, 2, 2, kMu_MeV, 0, 0.0, a, top.E_MeV + kSubtractedPoleWindow_MeV, bg);
  const double out = BackgroundR(p, 2, 2, kMu_MeV, 0, 0.0, a, top.E_MeV + 2e-3, bg);
  const double in = BackgroundR(p, 2, 2, kMu_MeV, 0, 0.0, a, top.E_MeV + kSubtractedPoleWindow_MeV - 1e-6, bg);
  std::snprintf(text, sizeof text, "subtracted pole: R0(E_p) = %.10g = the window's midpoint %.10g; R0 at 2 keV %.10g, "
                "1 eV inside the edge %.10g against the edge %.10g", r0, 0.5 * (rl + rh), out, in, rh);
  Check(std::isfinite(r0) && std::fabs(r0 - 0.5 * (rl + rh)) < 1e-12 * (1.0 + std::fabs(r0)) &&
            std::fabs(in - rh) < 1e-6 * (1.0 + std::fabs(rh)) && std::fabs(out - rh) < 1e-2 * (1.0 + std::fabs(rh)),
        text);
  // the capture kernel with the pole subtracted: K_unsubtracted - K_subtracted = gamma_p^2 I_p / (E_p - E), I_p the
  // overlap of the box state with the bound state, a constant: (K_un - K_sub)(E) (E_p - E) the same at three energies
  // (the subtraction acts on the capture as on R0; a subtracted pole's capture belongs to the level that replaces it)
  {
    HybridPotential pf = GaussianPotential(83.589554, 1.0 / std::sqrt(0.15747));
    pf.coulomb_radius_fm = 3.095;
    pf.step_fm = 0.0275;
    const auto w1 = [](double r) { return r; };
    InteriorBackground un;
    double c[3];
    const double Es[3] = {top.E_MeV - 0.7, top.E_MeV - 0.05, top.E_MeV + 0.4};
    for (int k = 0; k < 3; ++k) {
      const double ku = InteriorCaptureIntegral(p, 0, 0.0, Es[k], un, pf, 1, 0.0, -1.4548867, 2, 2, kMu_MeV, a, w1);
      const double ks = InteriorCaptureIntegral(p, 0, 0.0, Es[k], bg, pf, 1, 0.0, -1.4548867, 2, 2, kMu_MeV, a, w1);
      c[k] = (ku - ks) * (top.E_MeV - Es[k]);
    }
    std::snprintf(text, sizeof text, "subtracted capture kernel: (K - K_sub)(E_p - E) = %.10g, %.10g, %.10g at three energies",
                  c[0], c[1], c[2]);
    Check(std::fabs(c[1] / c[0] - 1.0) < 1e-6 && std::fabs(c[2] / c[0] - 1.0) < 1e-6, text);
  }
  // the closed-channel factor 1 - R0 (S - B) at a pole below threshold: finite
  const std::vector<BoxState> below = InteriorBoxStates(p, 2, 2, kMu_MeV, 0, 0.0, a, 0.0, -60.0, -0.01);
  if (!below.empty()) {
    const double f = BackgroundWFactor([&] {
      HybridPotential q = p;
      PartialWave pw;
      pw.l = 0;
      pw.background = InteriorBackground();
      q.waves = {pw};
      return q;
    }(), 2, 2, kMu_MeV, 0, -1, 0.0, a, below.back().E_MeV);
    std::snprintf(text, sizeof text, "W~_d/W~ at the bisected pole %.6g MeV below threshold finite: %.3e", below.back().E_MeV, f);
    Check(std::isfinite(f), text);
  }
}

}  // namespace

int main() {
  ZeroPotentialKernel();
  BoundStateRule();
  Poles();
  if (failures) {
    std::printf("potential_capture: %d FAILED\n", failures);
    return 1;
  }
  std::printf("potential_capture: all passed\n");
  return 0;
}
