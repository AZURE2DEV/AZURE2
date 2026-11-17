// The observed-level record carries the channel's P, mu and a (docs/upgrade/data_model.md section 10).
//
// On T1, T2, T3 (tests/upgrade) and tests/13N, in Brune's, the formal and Park's parametrization, the record of
// CompoundSystem::ObservedOfFitValues at the project's values, every particle channel of every level:
//  (a) P equals the second route a reader had before, ParticleChannelValues(pair, l, SpinOrbitTwoJ(J), E - S - e) at
//      the record's E (the same function on the same arguments: identical, the largest difference is printed);
//  (b) mu_u and a_fm are the pair's GetRedMass() and the channel's GetRadius();
//  (c) P is the factor W was formed with: W = s 2 gamma^2 P / N, N = 1 + sum_c gamma_c^2 dS_c/dE over the level's
//      particle channels (Brune's and the formal basis, whose record holds Brune's values), 1 under Park; a record
//      whose P was not the one W used fails here;
//  (d) photon and beta channels carry mu_u = a_fm = 0.
//
//   observed_level_record_test <T1.azr> <T2.azr> <T3.azr> <13N.azr>      (ctest: upgrade_observed_level_record)
#include <cmath>
#include <cstdio>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

#include <Minuit2/MnUserParameters.h>

#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "EngineChannelFunctions.h"
#include "JGroup.h"
#include "PPair.h"
#include "SetUpForTests.h"

namespace {

int failures = 0;

void check(const std::string &what, bool ok) {
  std::printf("  %-4s %s\n", ok ? "ok" : "FAIL", what.c_str());
  if (!ok) failures++;
}

double RelDiff(double got, double want) {
  if (got == want) return 0.0;
  return std::fabs(got - want) / std::fmax(std::fabs(want), 1e-300);
}

enum class Basis { Brune, Park, Formal };

void CheckProject(const std::string &name, const std::string &azr, Basis basis, const char *basisName) {
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = azr;
  if (basis == Basis::Park) cfg.paramMask |= (Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  if (basis == Basis::Formal) cfg.paramMask &= ~(Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  if (cfg.ReadConfigFile() < 0) throw std::runtime_error("cannot read " + azr);
  cfg.ResolveMatrixRoute();
  InitializeECAmplitudeCache(cfg);
  CompoundSetUp setUp = test_setup::WithoutData(cfg);
  if (setUp.refusal != CompoundSetUp::Refusal::kNone) throw std::runtime_error("the set-up failed for " + azr);
  const CompoundSystem &cn = *setUp.compound;
  // The fit values of the project's start, as a session holds them (AZUREAPI::UpdateParameters).
  AZUREParams params;
  cn.FillMnParams(params, setUp.state, setUp.inputs, &cfg);
  cn.FillCompoundFromParams(params.GetMinuitParams().Params(), setUp.state);
  cn.CalcShiftFunctions(cfg, setUp.state);
  const std::optional<ObservedLevelTable> table = cn.ObservedOfFitValues(cfg, setUp.state, setUp.inputs);
  check(name + ", " + basisName + ": the transformation gives a record", table.has_value());
  if (!table) return;
  const bool park = basis == Basis::Park;
  long n = 0, identical = 0, open = 0, closed = 0, others = 0;
  double worstP = 0.0, worstW = 0.0;
  bool pairOk = true, othersOk = true;
  for (int j = 1; j <= cn.NumJGroups(); j++) {
    const JGroup *group = cn.GetJGroup(j);
    for (int la = 1; la <= group->NumLevels(); la++) {
      const ObservedLevel &level = table->Level(j, la);
      // Brune's N of the level, from the record's amplitudes and the second route's dS/dE
      double N = 1.0;
      for (int ch = 1; ch <= group->NumChannels(); ch++) {
        const AChannel *c = group->GetChannel(ch);
        if (c->GetRadType() != 'P') continue;
        const PPair *pair = cn.GetPair(c->GetPairNum());
        const double g = level.gamma[ch - 1];
        N += ParticleChannelValues(*pair, c->GetL(), pair->SpinOrbitTwoJ(group->GetJ()), c->GetRadius(), level.E - pair->GetSepE() - pair->GetExE()).dS_dE * (g * g);
      }
      if (park) N = 1.0;
      for (int ch = 1; ch <= group->NumChannels(); ch++) {
        const AChannel *c = group->GetChannel(ch);
        if (c->GetRadType() != 'P') {
          others++;
          othersOk = othersOk && level.mu_u[ch - 1] == 0.0 && level.a_fm[ch - 1] == 0.0;
          continue;
        }
        const PPair *pair = cn.GetPair(c->GetPairNum());
        const double Ec = level.E - pair->GetSepE() - pair->GetExE();
        const double P = ParticleChannelValues(*pair, c->GetL(), pair->SpinOrbitTwoJ(group->GetJ()), c->GetRadius(), Ec).P;
        n++;
        (Ec > 0.0 ? open : closed)++;
        if (level.P[ch - 1] == P) identical++;
        worstP = std::fmax(worstP, RelDiff(level.P[ch - 1], P));
        pairOk = pairOk && level.mu_u[ch - 1] == pair->GetRedMass() && level.a_fm[ch - 1] == c->GetRadius();
        const double g = level.gamma[ch - 1];
        const double W = (g < 0.0 ? -1.0 : 1.0) * 2.0 * g * g * level.P[ch - 1] / N;
        worstW = std::fmax(worstW, RelDiff(level.W[ch - 1], W));
      }
    }
  }
  char line[300];
  std::snprintf(line, sizeof(line), "%s, %s: (a) P of %ld particle channels (%ld open, %ld closed) = "
                "ParticleChannelValues at the record's E: %ld identical, largest relative difference %.2g",
                name.c_str(), basisName, n, open, closed, identical, worstP);
  check(line, n > 0 && worstP <= 1e-15);
  check(name + ", " + basisName + ": (b) mu_u and a_fm are the pair's GetRedMass() and the channel's GetRadius()", n > 0 && pairOk);
  std::snprintf(line, sizeof(line), "%s, %s: (c) W = s 2 gamma^2 P / N with the record's P, largest relative "
                "difference %.2g (<= 1e-13)", name.c_str(), basisName, worstW);
  check(line, worstW <= 1e-13);
  std::snprintf(line, sizeof(line), "%s, %s: (d) %ld photon or beta channels carry mu_u = a_fm = 0", name.c_str(),
                basisName, others);
  check(line, othersOk);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::fprintf(stderr, "usage: %s T1.azr T2.azr T3.azr 13N.azr\n", argv[0]);
    return 2;
  }
  const char *names[4] = {"T1", "T2", "T3", "13N"};
  const struct {
    Basis basis;
    const char *name;
  } bases[3] = {{Basis::Brune, "Brune"}, {Basis::Formal, "formal"}, {Basis::Park, "Park"}};
  for (int t = 0; t < 4; t++)
    for (const auto &b : bases) {
      try {
        CheckProject(names[t], argv[t + 1], b.basis, b.name);
      } catch (const std::exception &e) {
        std::printf("  FAIL %s, %s threw: %s\n", names[t], b.name, e.what());
        failures++;
      }
    }
  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
