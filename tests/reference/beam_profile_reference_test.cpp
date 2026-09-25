/*!
 * Reference calculation for the beam-profile experimental effect.
 *
 * The physics regression suite (tests/beam_profile_kernel) pins the chi-squared
 * of a real project against a recorded number: it catches a CHANGE, but it
 * cannot say the number was ever right, because the reference came from the
 * same code.  This test is the independent check.  It drives
 * EPoint::IntegrateTargetEffect directly, with sub-point cross sections chosen
 * by hand, and compares against
 *
 *   - closed-form values, where the fold has one (a constant cross section must
 *     come back exactly, whatever the profile, window or detailed-balance
 *     weight is doing), and
 *   - a second implementation of the documented integral, written here from the
 *     formulas in docs/source/user_guide/experimental_effects.rst and
 *     integrated by dense Simpson quadrature rather than by the engine's
 *     2-point Gauss-Legendre rule over sub-point intervals.
 *
 * What the engine computes, for a point of energy E0 whose detector resolution
 * window is [a, b] and whose sub-points carry sigma(E_i):
 *
 *     Y = Int K(E) sigma(E) dE / Int K(E) dE,
 *     K(E) = G(E) * W(E) * D(E),
 *     G(E) = sum_k w_k / (omega_k sqrt(2 pi)) * exp(-z^2/2) * (1 + erf(alpha_k z / sqrt2)),
 *            z = (E - xi_k) / omega_k,
 *     W(E) = 1/2 [ erf((b - E)/(s sqrt2)) - erf((a - E)/(s sqrt2)) ],  1 if no window,
 *     D(E) = g(E) / g(E0),  g(E) = E / E_gamma(E)^2,  1 unless dbFlag is set,
 *     E_gamma(E) = M (1 - sqrt(1 - 2 (E + Q) / M)).
 *
 * Run:  tests/reference/beam_profile_reference_test          (ctest: beam_profile_reference)
 */
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Config.h"
#include "EData.h"
#include "EPoint.h"
#include "TargetEffect.h"

// Each consumer of the engine defines this itself -- src/AZURE2.cpp for the CLI,
// api/src/pybind_api.cpp for the Python API -- so a program that links
// AZURE2Core has to as well.  CoulFunc reads it at construction.
Config *g_config = nullptr;

namespace {

// Tolerance for anything compared against a closed form evaluated with the
// machine's own pi.  Constants.h defines `pi = 3.141592650`, which is the true
// value truncated at ten digits, a relative error of 1.14e-9; the profile
// weight carries 1/(omega sqrt(2 pi)) and so is off by half of that, 5.7e-10.
// It cancels in the fold itself, where numerator and denominator share it,
// which is why the identities below hold to 1e-12.  Give the constant its full
// precision and this tolerance can go to 1e-12 as well.
const double kPiTolerance = 5.0e-9;

int failures = 0;
int checks = 0;

void check(const std::string &name, bool ok, const std::string &detail = "") {
  checks++;
  std::cout << "  " << (ok ? "ok  " : "FAIL") << "  " << name;
  if (!ok && !detail.empty()) std::cout << "  -- " << detail;
  std::cout << std::endl;
  if (!ok) failures++;
}

void checkClose(const std::string &name, double got, double want, double rtol) {
  double scale = std::max(std::fabs(want), 1.0e-300);
  double rel = std::fabs(got - want) / scale;
  char buf[256];
  snprintf(buf, sizeof(buf), "got %.12g, want %.12g, relative %.3g > %.3g", got, want, rel, rtol);
  check(name, rel <= rtol, buf);
}

// ---------------------------------------------------------------- the reference

struct Component {
  double xi, omega, alpha, weight;
};

struct Kernel {
  std::vector<Component> profile;
  double tpcSigma = 0.0;       // 0 = no resolution window
  double truncation = 0.0;     // 0 = untruncated
  bool photodissociation = false;
  double qGamma = 0.0, mass = 0.0, e0 = 0.0;
  double a = 0.0, b = 0.0;     // window, already shifted
  bool window = false;

  double G(double e) const {
    double sum = 0.0;
    for (const Component &c : profile) {
      double z = (e - c.xi) / c.omega;
      if (truncation > 0.0) {
        double delta = c.alpha / std::sqrt(1.0 + c.alpha * c.alpha);
        double mean = c.xi + c.omega * delta * std::sqrt(2.0 / M_PI);
        double sd = c.omega * std::sqrt(1.0 - 2.0 * delta * delta / M_PI);
        if (std::fabs(e - mean) > truncation * sd) continue;
      }
      sum += c.weight * std::exp(-0.5 * z * z) / (c.omega * std::sqrt(2.0 * M_PI)) *
             (1.0 + std::erf(c.alpha * z / std::sqrt(2.0)));
    }
    return sum;
  }

  double W(double e) const {
    if (!window) return 1.0;
    if (!(tpcSigma > 0.0)) return (e >= a && e <= b) ? 1.0 : 0.0;
    return 0.5 * (std::erf((b - e) / (tpcSigma * std::sqrt(2.0))) -
                  std::erf((a - e) / (tpcSigma * std::sqrt(2.0))));
  }

  double gOf(double e) const {
    double sum = e + qGamma;
    double arg = 1.0 - 2.0 * sum / mass;
    double egamma = (arg > 0.0) ? mass * (1.0 - std::sqrt(arg)) : sum;
    return (egamma > 0.0) ? e / (egamma * egamma) : 0.0;
  }

  double D(double e) const {
    if (!photodissociation) return 1.0;
    double norm = gOf(e0);
    return (norm > 0.0) ? gOf(e) / norm : 1.0;
  }

  double operator()(double e) const { return G(e) * W(e) * D(e); }
};

// Dense Simpson fold of K*sigma / K over [lo, hi]; sigma is evaluated exactly,
// not interpolated, so this shares no machinery with the engine.
template <typename Sigma>
double referenceFold(const Kernel &k, Sigma sigma, double lo, double hi, int intervals = 200000) {
  if (intervals % 2) intervals++;
  double h = (hi - lo) / intervals;
  double num = 0.0, den = 0.0;
  for (int i = 0; i <= intervals; i++) {
    double e = lo + i * h;
    double w = (i == 0 || i == intervals) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    double kk = k(e);
    num += w * kk * sigma(e);
    den += w * kk;
  }
  return num / den;
}

// ---------------------------------------------------------------- the engine

// One targetInt line carrying only the beam-profile block; the leading legacy
// fields are inert, exactly as pyazr's AzrModel.add_target_effect writes them.
std::string targetIntLine(const std::string &profileBlock, int numPoints = 400) {
  std::ostringstream o;
  o << "1  \"1\"  " << numPoints << "  0  0  0  0  \"\"  0  0  0  0  \"\"  0  0  0.04  10  50  " << profileBlock;
  return o.str();
}

TargetEffect makeEffect(const Config &configure, const std::string &line) {
  std::istringstream in(line);
  return TargetEffect(in, configure);
}

// Fold sigma(E) through the engine: build a point with sub-points on a uniform
// grid over [lo, hi], hand it the effect, and call IntegrateTargetEffect.
template <typename Sigma>
double engineFold(const Config &configure, EData &data, int effectNum, Sigma sigma,
                  double e0, double lo, double hi, int nSub,
                  double shift = 0.0, const double *window = nullptr,
                  double qGamma = 0.0, double mass = 0.0) {
  EPoint point(0.0, e0, 1, 2, false, false, false, 0.0, 0, 0);
  point.SetParentData(&data);
  point.SetTargetEffectNum(effectNum);
  if (window) point.SetBinWindow(window[0], window[1]);
  if (mass > 0.0) point.SetPhotoKinematics(qGamma, mass);
  // A shift moves the point and its sub-points; the window stays where it was
  // and the kernel is expected to shift it by the same amount.  Keeping the
  // lab and c.m. energies equal makes that amount exactly `shift`.
  point.SetLabEnergy(e0 + shift);
  point.SetCMEnergy(e0 + shift);
  for (int i = 0; i < nSub; i++) {
    double e = hi - (hi - lo) * i / (nSub - 1);   // engine convention: high to low
    EPoint sub(0.0, e, 1, 2, false, false, false, 0.0, 0, 0);
    sub.SetFitCrossSection(sigma(e));
    point.AddSubPoint(sub);
  }
  point.IntegrateTargetEffect(configure);
  return point.GetFitCrossSection();
}

}  // namespace

int main() {
  std::ostringstream sink;
  Config configure(sink);
  g_config = &configure;

  std::cout << "beam-profile reference calculation" << std::endl;

  // ------------------------------------------------ 1. parsing and the profile
  const std::string oneComponent = "beamprofile  1  3.5  0.30  -1.6  1  0.0733  0  0";
  TargetEffect single = makeEffect(configure, targetIntLine(oneComponent));
  check("a beamprofile block is recognised", single.IsBeamProfile());
  check("and counts as a sub-point effect", single.IsSubPointEffect());
  check("one component is read", single.GetBeamProfile().size() == 1);
  checkClose("its location", single.GetBeamProfile()[0].xi, 3.5, 1e-12);
  checkClose("its scale", single.GetBeamProfile()[0].omega, 0.30, 1e-12);
  checkClose("its skewness", single.GetBeamProfile()[0].alpha, -1.6, 1e-12);
  checkClose("the detector resolution", single.GetBeamTpcSigma(), 0.0733, 1e-12);
  check("the detailed-balance flag is off here", !single.IsBeamPhotodissociation());

  Kernel ref;
  ref.profile = {{3.5, 0.30, -1.6, 1.0}};
  ref.tpcSigma = 0.0733;
  bool weightOk = true;
  double worstWeight = 0.0;
  for (double e = 2.0; e <= 5.0; e += 0.05) {
    double got = single.BeamProfileWeight(e), want = ref.G(e);
    double rel = std::fabs(got - want) / std::max(std::fabs(want), 1e-300);
    worstWeight = std::max(worstWeight, rel);
    if (rel > kPiTolerance) weightOk = false;
  }
  {
    char buf[128];
    snprintf(buf, sizeof(buf), "worst relative deviation %.3g over 2-5 MeV", worstWeight);
    check("BeamProfileWeight matches the closed-form skewed Gaussian", weightOk, buf);
  }

  // two components, unequal weights -- the sum, not the last one
  TargetEffect twin = makeEffect(configure, targetIntLine(
      "beamprofile  2  3.4  0.25  -1.2  0.3  3.8  0.32  -2.0  0.7  0.0733  0  0"));
  Kernel refTwin;
  refTwin.profile = {{3.4, 0.25, -1.2, 0.3}, {3.8, 0.32, -2.0, 0.7}};
  checkClose("a two-component profile sums its components",
             twin.BeamProfileWeight(3.6), refTwin.G(3.6), kPiTolerance);

  // truncation
  TargetEffect cut = makeEffect(configure, targetIntLine(
      "beamprofile  1  3.5  0.30  -1.6  1  0.0733  2  0"));
  Kernel refCut = ref;
  refCut.truncation = 2.0;
  checkClose("truncation keeps the centre", cut.BeamProfileWeight(3.4), refCut.G(3.4), kPiTolerance);
  check("truncation zeroes the far tail", cut.BeamProfileWeight(2.2) == 0.0);
  check("and the untruncated profile does not", single.BeamProfileWeight(2.2) > 0.0);

  double lo = 0.0, hi = 0.0;
  single.BeamProfileSupport(lo, hi);
  check("the support brackets the profile", lo < 3.5 && hi > 3.5 &&
        single.BeamProfileWeight(lo) < 1e-4 * single.BeamProfileWeight(3.4));

  // lab -> c.m. conversion happens once even when several segments share the effect
  TargetEffect shared = makeEffect(configure, targetIntLine(oneComponent));
  shared.ConvertBeamProfileToCM(0.75);
  shared.ConvertBeamProfileToCM(0.75);
  checkClose("the frame conversion is applied once, not once per segment",
             shared.GetBeamProfile()[0].xi, 3.5 * 0.75, 1e-12);

  // ------------------------------------------------ 2. the fold
  EData data;
  data.AddTargetEffect(single);
  const int effect = data.NumTargetEffects();
  const double e0 = 3.5, gridLo = 2.6, gridHi = 4.4;
  const double window[2] = {3.40, 3.60};
  const double qGamma = 7.16191, mass = 14887.0;   // 12C + alpha, MeV

  Kernel kNoWindow = ref;
  kNoWindow.window = false;

  // (a) the exact identity: a constant folds to itself, whatever the kernel is.
  auto constant = [](double) { return 4.2e-6; };
  checkClose("a constant cross section folds to itself (no window)",
             engineFold(configure, data, effect, constant, e0, gridLo, gridHi, 400),
             4.2e-6, 1e-12);
  checkClose("a constant cross section folds to itself (with the window)",
             engineFold(configure, data, effect, constant, e0, gridLo, gridHi, 400, 0.0, window),
             4.2e-6, 1e-12);
  checkClose("a constant cross section folds to itself (window and detailed balance)",
             engineFold(configure, data, effect, constant, e0, gridLo, gridHi, 400, 0.0, window, qGamma, mass),
             4.2e-6, 1e-12);

  // (b) a linear cross section, against the dense reference
  auto linear = [](double e) { return 1.0e-6 * (0.5 + 0.8 * e); };
  Kernel kWin = ref;
  kWin.window = true;
  kWin.a = window[0];
  kWin.b = window[1];
  checkClose("a linear cross section matches the reference fold",
             engineFold(configure, data, effect, linear, e0, gridLo, gridHi, 800, 0.0, window),
             referenceFold(kWin, linear, gridLo, gridHi), 1e-6);

  // (c) curvature: a narrow Gaussian bump inside the window, where the 2-point
  //     rule and the linear interpolation between sub-points both have to work
  auto bump = [](double e) { return 1.0e-6 * (1.0 + 30.0 * std::exp(-0.5 * std::pow((e - 3.52) / 0.05, 2.0))); };
  double coarse = engineFold(configure, data, effect, bump, e0, gridLo, gridHi, 200, 0.0, window);
  double fine = engineFold(configure, data, effect, bump, e0, gridLo, gridHi, 3000, 0.0, window);
  double truth = referenceFold(kWin, bump, gridLo, gridHi);
  checkClose("a curved cross section matches the reference fold (3000 sub-points)", fine, truth, 2e-4);
  check("and the quadrature converges toward it as sub-points are added",
        std::fabs(fine - truth) < std::fabs(coarse - truth),
        "refining the grid did not move the result toward the reference");

  // (d) the window follows the segment's energy shift
  const double shift = 0.12;
  Kernel kShift = kWin;
  kShift.a = window[0] + shift;
  kShift.b = window[1] + shift;
  checkClose("an energy shift moves the resolution window with the point",
             engineFold(configure, data, effect, linear, e0, gridLo + shift, gridHi + shift, 800, shift, window),
             referenceFold(kShift, linear, gridLo + shift, gridHi + shift), 1e-6);

  // (e) the detailed-balance weight, which is what makes a capture cross section
  //     average the way the inverse photodissociation measurement averaged it
  EData dataDb;
  dataDb.AddTargetEffect(makeEffect(configure, targetIntLine(
      "beamprofile  1  3.5  0.30  -1.6  1  0.0733  0  1")));
  const int effectDb = dataDb.NumTargetEffects();
  Kernel kDb = kWin;
  kDb.photodissociation = true;
  kDb.qGamma = qGamma;
  kDb.mass = mass;
  kDb.e0 = e0;
  double withDb = engineFold(configure, dataDb, effectDb, linear, e0, gridLo, gridHi, 800, 0.0, window, qGamma, mass);
  checkClose("the detailed-balance weight matches the reference",
             withDb, referenceFold(kDb, linear, gridLo, gridHi), 1e-6);
  check("and it is not a no-op",
        std::fabs(withDb / referenceFold(kWin, linear, gridLo, gridHi) - 1.0) > 1e-6);

  std::cout << (failures ? "FAILED " : "passed ") << checks - failures << "/" << checks
            << " checks" << std::endl;
  return failures ? 1 : 0;
}
