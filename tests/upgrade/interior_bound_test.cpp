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
//  (d) the set-up warning (WarnLevelDoublesPotential): InteriorPotentialStates finds the bound states of (c) and the
//      l = 1 resonance of the potential near 2.40 MeV (about 410 keV wide), and of the project's three levels the one
//      at 9.585 MeV (E_c = 2.42 MeV) is reported, once.
//  (e) without levels the channel's phase is the potential's: the engine's exp(-i phi') at a = 5.35 fm against the
//      phase shift of the regular solution integrated to 14 fm, where the nuclear potential has died out, and matched
//      to the waves there (no R_pot, no boundary constant, another radius), l = 0-3 at 0.3-6 MeV.
//
//   interior_bound_test <interior_bound.azr>      (ctest: upgrade_interior_bound)
#include <cmath>
#include <cstdio>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Minuit2/MnUserParameters.h>

#include "AChannel.h"
#include "AZUREParams.h"
#include "CompoundSystemChecks.h"
#include "CoulFunc.h"
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
  int nPoints = 0, nStates = 0, statesMissed = 0;
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
          bool found = false;
          for (const PotentialState &s : InteriorPotentialStates(*pair, l, twoJ, a, -Bb - 1.0, -Bb + 1.0))
            if (std::fabs(s.Ec + Bb) < 1e-6 && s.width == 0.0) found = true;
          if (!found) statesMissed++;
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
  // (e)
  {
    const double kPi = std::acos(-1.0);
    double worst = 0.0;
    int n = 0;
    CoulFunc coul(pair);
    for (int l = 0; l <= 3; l++) {
      const phys::HybridPotential &pot = *pair->GetHybridPotential(l);
      const int twoJ = pair->SpinOrbitTwoJ(l);
      for (double E = 0.3; E < 6.01; E += 0.3) {
        const double rFar = 14.0, rho = std::sqrt(2.0 * mu * E) / hbarc * rFar;
        const double f = phys::RegularSolution(pot, pair->GetZ(1), pair->GetZ(2), mu, l, pair->SpinOrbitFactor(l, twoJ), rFar, E).log_derivative;
        const CoulWaves w = coul(l, twoJ, rFar, E);
        const double direct = std::atan2(rho * w.dF - f * w.F, f * w.G - rho * w.dG);
        const double engine = std::arg(coul.HardSpherePhase(l, twoJ, a, E));
        double d = std::fmod(engine - direct, kPi);  // phases are equal modulo pi
        if (d > 0.5 * kPi) d -= kPi;
        if (d < -0.5 * kPi) d += kPi;
        if (std::fabs(d) > 1e-4) std::printf("    l=%d E=%.2f MeV: engine %.6f, direct %.6f rad\n", l, E, engine, direct);
        worst = std::fmax(worst, std::fabs(d));
        n++;
      }
    }
    std::snprintf(line, sizeof line, "(e) the phase without levels is the potential's phase shift: %d points, largest difference %.2g rad", n, worst);
    check(line, worst < 1e-4);
  }
  // (d)
  const std::vector<PotentialState> resonances = InteriorPotentialStates(*pair, 1, pair->SpinOrbitTwoJ(1), a, 0.05, 10.0);
  std::snprintf(line, sizeof line, "(d) the states of the potential: the %d bound states found again, one l = 1 resonance below 10 MeV (%.4f MeV, %.1f keV)",
                nStates, resonances.empty() ? 0.0 : resonances[0].Ec, resonances.empty() ? 0.0 : resonances[0].width * 1e3);
  check(line, statesMissed == 0 && resonances.size() == 1 && std::fabs(resonances[0].Ec - 2.40) < 0.02 &&
                  std::fabs(resonances[0].width - 0.41) < 0.02);
  AZUREParams params;
  cn.FillMnParams(params, setUp.state, setUp.inputs, &cfg);
  cn.FillCompoundFromParams(params.GetMinuitParams().Params(), setUp.state);
  cn.CalcShiftFunctions(cfg, setUp.state);
  const std::optional<ObservedLevelTable> table = cn.ObservedOfFitValues(cfg, setUp.state, setUp.inputs);
  check("(d) the observed values of the project's levels", table.has_value());
  if (table) {
    sink.str("");
    const int warned = WarnLevelDoublesPotential(cn, setUp.inputs, *table, cfg);
    const std::string text = sink.str();
    check("(d) one level is reported, the one at 9.585 MeV",
          warned == 1 && text.find("E=9.5850") != std::string::npos && text.find("resonance of the potential") != std::string::npos);
    if (warned != 1) std::printf("%s", text.c_str());
  }
  std::printf("%s\n", failures ? "FAILED" : "all checks passed");
  return failures ? 1 : 0;
}
