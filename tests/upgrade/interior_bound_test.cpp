// A bound level in a channel whose potential also acts inside the channel radius (docs/upgrade/interior_potential.md
// section 6).
//
// On tests/upgrade/interior_bound (alpha + 12C, a Woods-Saxon potential with `interior_boundary 0`, a = 5.35 fm), for
// l = 0-3 and binding energies 0.05-30 MeV:
//  (a) the width factor of ParticleChannelValues is mu a/(hbar^2 W^2) (f - b)^2/(f - S)^2, with f from the regular
//      solution, S from the decaying one and W from WhitFunc;
//  (b) the norm identity: dS'/dE, which ParticleChannelValues differentiates numerically from the primed shift, equals
//      (2 mu a/hbar^2) (N_f + ((S - b)/(f - b))^2 int_0^a u^2/u(a)^2) / D^2, N_f = int_a^inf W^2/W(a)^2 (NFIntegral):
//      the norm an external-capture final state is given (ExternalCaptureFinalStates.cpp) is the one the level
//      conversion uses;
//  (c) at a bound state of the potential itself (f = S, found by bisection) the limit of a large reduced width,
//      C^2 = 2 P / (dS'/dE), is the single-particle ANC^2 of the potential, 1/(W(a)^2 (int_0^a u^2/u(a)^2 + N_f)),
//      as the mean of its values 1 keV above and below the state (both factors have a pole on it).
//
//   interior_bound_test <interior_bound.azr>      (ctest: upgrade_interior_bound)
#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>

#include "AChannel.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Constants.h"
#include "ECAmplitudeCache.h"
#include "EngineChannelFunctions.h"
#include "JGroup.h"
#include "NFIntegral.h"
#include "PPair.h"
#include "SetUpForTests.h"
#include "WhitFunc.h"
#include "azure2/physics/HybridPotential.h"

namespace phys = azure2::physics;

namespace {
int failures = 0;
void check(const std::string &what, bool ok) {
  std::printf("  %-4s %s\n", ok ? "ok" : "FAIL", what.c_str());
  if (!ok) failures++;
}
}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) return 2;
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = argv[1];
  if (cfg.ReadConfigFile() < 0) throw std::runtime_error(std::string("cannot read ") + argv[1]);
  cfg.ResolveMatrixRoute();
  InitializeECAmplitudeCache(cfg);
  CompoundSetUp setUp = test_setup::WithoutData(cfg);
  if (setUp.refusal != CompoundSetUp::Refusal::kNone) throw std::runtime_error("the set-up failed");
  const CompoundSystem &cn = *setUp.compound;
  const PPair *pair = nullptr;
  double a = 0.0;
  for (int j = 1; j <= cn.NumJGroups() && !pair; j++)
    for (int ch = 1; ch <= cn.GetJGroup(j)->NumChannels() && !pair; ch++)
      if (cn.GetJGroup(j)->GetChannel(ch)->GetRadType() == 'P') {
        pair = cn.GetPair(cn.GetJGroup(j)->GetChannel(ch)->GetPairNum());
        a = cn.GetJGroup(j)->GetChannel(ch)->GetRadius();
      }
  if (!pair) throw std::runtime_error("no particle channel");
  check("the pair's potential acts inside the radius", pair->GetHybridPotential(1) && pair->GetHybridPotential(1)->interior_boundary);
  const double mu = pair->GetRedMass() * uconv, toNorm = 2.0 * mu * a / (hbarc * hbarc);

  double worstP = 0.0, worstN = 0.0, worstC = 0.0;
  int nPoints = 0, nStates = 0;
  for (int l = 0; l <= 3; l++) {
    const phys::HybridPotential &pot = *pair->GetHybridPotential(l);
    const double b = *pot.interior_boundary;
    const int twoJ = pair->SpinOrbitTwoJ(l);
    const double twoLDotI = pair->SpinOrbitFactor(l, twoJ);
    WhitFunc whit(pair);
    NFIntegral nf(pair, a);
    auto fMinusS = [&](double B) {
      return phys::RegularSolution(pot, pair->GetZ(1), pair->GetZ(2), mu, l, twoLDotI, a, -B).log_derivative -
             phys::HybridBoundShift(pot, pair->GetZ(1), pair->GetZ(2), mu, l, twoLDotI, a, B);
    };
    // The right side of (b) over (2 mu a / hbar^2) / D^2, and the limit of (c).
    auto expected = [&](double B, double &dInverse, double &W) {
      const phys::RegularValue u = phys::RegularSolution(pot, pair->GetZ(1), pair->GetZ(2), mu, l, twoLDotI, a, -B);
      const double S = phys::HybridBoundShift(pot, pair->GetZ(1), pair->GetZ(2), mu, l, twoLDotI, a, B);
      dInverse = (u.log_derivative - b) / (u.log_derivative - S);
      W = whit(l, twoJ, a, B);
      return std::pow((S - b) / (u.log_derivative - b), 2.0) * u.norm_fm + nf(l, twoJ, pair->GetSepE() + pair->GetExE() - B);
    };
    auto pointCheck = [&](double B) {
      double dInverse, W;
      const double sum = expected(B, dInverse, W);
      const phys::ChannelValues v = ParticleChannelValues(*pair, l, twoJ, a, -B);
      const double wantP = mu * a / (hbarc * hbarc) / (W * W) * dInverse * dInverse;
      worstP = std::fmax(worstP, std::fabs(v.P / wantP - 1.0));
      const double wantN = toNorm * sum * dInverse * dInverse;
      const double dN = std::fabs(v.dS_dE / wantN - 1.0);
      if (dN > 1e-4) std::printf("    l=%d B=%.4f MeV: dS'/dE %.8g, identity %.8g (1/D %.4g)\n", l, B, v.dS_dE, wantN, dInverse);
      worstN = std::fmax(worstN, dN);
      nPoints++;
      return 2.0 * v.P / v.dS_dE;  // C^2 in the limit of a large reduced width
    };
    double prevB = 0.05, prev = fMinusS(prevB);
    for (double B = 0.3; B < 30.0; B += 0.25) {
      const double g = fMinusS(B);
      if (std::isfinite(g) && std::isfinite(prev) && g * prev < 0.0) {
        double lo = prevB, hi = B, glo = prev;
        for (int it = 0; it < 60; it++) {
          const double mid = 0.5 * (lo + hi), gm = fMinusS(mid);
          if (gm * glo < 0.0) hi = mid;
          else { lo = mid; glo = gm; }
        }
        const double Bb = 0.5 * (lo + hi);
        if (std::fabs(fMinusS(Bb)) < 1e-6) {  // a zero of f - S (a bound state), not a pole of f
          double dInverse, W;
          const double sum = expected(Bb, dInverse, W);  // at f = S: int_0^a u^2/u(a)^2 + N_f
          const double single = 1.0 / (W * W * sum);
          // The mean of the two sides: the limit varies linearly across the state, and both factors have a pole on it.
          const double c2 = 0.5 * (pointCheck(Bb - 1e-3) + pointCheck(Bb + 1e-3));
          worstC = std::fmax(worstC, std::fabs(c2 / single - 1.0));
          std::printf("    l=%d: bound state of the potential at %.6f MeV below threshold, single-particle ANC^2 %.6g fm^-1\n", l, Bb, single);
          nStates++;
        }
      }
      pointCheck(B);
      prevB = B;
      prev = g;
    }
  }
  char line[240];
  std::snprintf(line, sizeof line, "(a) the width factor carries (f - b)^2/(f - S)^2: %d points, largest relative difference %.2g", nPoints, worstP);
  check(line, worstP < 1e-10);
  std::snprintf(line, sizeof line, "(b) dS'/dE equals the exterior and interior norm integrals over D^2: largest relative difference %.2g", worstN);
  check(line, worstN < 1e-4);
  std::snprintf(line, sizeof line, "(c) %d bound states of the potential: the large-width limit is the single-particle ANC^2 to %.2g", nStates, worstC);
  check(line, nStates > 0 && worstC < 1e-4);
  std::printf("%s\n", failures ? "FAILED" : "all checks passed");
  return failures ? 1 : 0;
}
