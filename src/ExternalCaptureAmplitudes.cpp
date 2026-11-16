#include "ExternalCaptureAmplitudes.h"

#include <cmath>
#include <iostream>
#include <memory>

#include "AChannel.h"
#include "GSLException.h"
#include "ALevel.h"
#include "CompoundSystem.h"
#include "EngineChannelFunctions.h"
#include "SystemOfCompound.h"
#include "Config.h"
#include "ECIntegral.h"
#include "JGroup.h"
#include "PPair.h"
#include "azure2/physics/ExternalCapture.h"
#include "azure2/physics/LevelParameters.h"

namespace phys = azure2::physics;

double SqrtNFFactor(const CompoundSystem &compound, const JGroup *group, const ALevel *level, double n0,
                    const std::vector<double> &gamma) {
  std::vector<phys::BoundChannel> channels;
  for (int ch = 1; ch <= group->NumChannels() && ch <= level->NumNFIntegrals(); ch++) {
    const PPair *pair = compound.GetPair(group->GetChannel(ch)->GetPairNum());
    channels.push_back({group->GetChannel(ch)->GetRadius(), pair->GetRedMass(), gamma[ch - 1], level->GetNFIntegral(ch)});
  }
  return 1.0 / std::sqrt(phys::bound_state_normalization(n0, channels));
}

std::pair<int, int> ExternalCaptureFinalState(const CompoundSystem &compound, const AChannel *channel) {
  const char type = channel->GetRadType();
  if (type != 'E' && !(type == 'M' && channel->GetL() == 1)) return {0, 0};
  for (int j = 1; j <= compound.NumJGroups(); j++)
    for (int la = 1; la <= compound.GetJGroup(j)->NumLevels(); la++) {
      const ALevel *f = compound.GetJGroup(j)->GetLevel(la);
      if (!f->IsECLevel() || channel->GetPairNum() != f->GetECPairNum()) continue;
      const int L = channel->GetL();
      const unsigned char mask = f->GetECMultMask();
      if ((type == 'E' && L == 1 && (mask & isE1)) || (type == 'M' && L == 1 && (mask & isM1)) ||
          (type == 'E' && L == 2 && (mask & isE2)))
        return {j, la};
      return {0, 0};
    }
  return {0, 0};
}

complex ExternalReducedWidthAt(const CompoundSystem &compound, const JGroup *group, const AChannel *channel,
                               std::pair<int, int> finalState,
                               const ExternalCaptureInputs &values, const Config &configure,
                               ExternalCaptureSlopes *slopes) {
  if (slopes) {
    slopes->level.assign(group->NumChannels(), complex(0.0, 0.0));
    slopes->finalScale.clear();
  }
  if (finalState.first == 0) return complex(0.0, 0.0);
  const char type = channel->GetRadType();
  const int L = channel->GetL();
  const JGroup *finalGroup = compound.GetJGroup(finalState.first);
  const ALevel *finalLevel = finalGroup->GetLevel(finalState.second);
  std::vector<phys::CaptureChannel> in;
  std::vector<const AChannel *> inChannels;
  std::vector<int> inNumber;
  for (int ch = 1; ch <= group->NumChannels(); ch++) {
    const AChannel *c = group->GetChannel(ch);
    if (c->GetRadType() != 'P') continue;
    in.push_back({c->GetPairNum(), c->GetL(), c->GetS(), values.levelGamma[ch - 1]});
    inChannels.push_back(c);
    inNumber.push_back(ch);
  }
  std::vector<phys::FinalChannel> out;
  std::vector<const AChannel *> outChannels;
  std::vector<int> outNumber;
  for (int chp = 1; chp <= finalGroup->NumChannels(); chp++) {
    const AChannel *c = finalGroup->GetChannel(chp);
    if (c->GetRadType() != 'P') continue;
    // s_fc = gamma~_fc N_f^(-1/2) K_c (Eq. 2), in the engine's order.
    out.push_back({c->GetPairNum(), c->GetL(), c->GetS(),
                   values.finalGamma[chp - 1] * values.finalSqrtNF * finalLevel->GetECConversionFactor(chp)});
    outChannels.push_back(c);
    outNumber.push_back(chp);
  }
  // The integrals of the pairs that couple, kept for the slopes (Eq. 4 is bilinear in gamma_lambda c and s_fc).
  std::vector<complex> integrals(in.size() * out.size());
  std::vector<unsigned char> known(in.size() * out.size(), 0);
  const phys::RadialIntegral integral = [&](std::size_t c, std::size_t cf) {
    if (known[c * out.size() + cf]) return integrals[c * out.size() + cf];
    ECIntegral theECIntegral(compound.GetPair(outChannels[cf]->GetPairNum()), configure, outChannels[cf]->GetRadius());
    const complex v = theECIntegral(inChannels[c]->GetL(), outChannels[cf]->GetL(), inChannels[c]->GetS(),
                                    outChannels[cf]->GetS(), group->GetJ(), finalGroup->GetJ(), L, type,
                                    values.levelEnergy, values.finalEnergy, true);
    integrals[c * out.size() + cf] = v;
    known[c * out.size() + cf] = 1;
    return v;
  };
  const phys::Multipole X = type == 'M' ? phys::Multipole::M : phys::Multipole::E;
  const complex value = phys::external_reduced_width(X, L, finalGroup->GetJ(), in, out, integral);
  if (!slopes) return value;
  slopes->finalScale.assign(finalGroup->NumChannels(), complex(0.0, 0.0));
  for (std::size_t c = 0; c < in.size(); c++) {
    std::vector<phys::CaptureChannel> unit = in;
    for (std::size_t k = 0; k < unit.size(); k++) unit[k].gamma = k == c ? 1.0 : 0.0;
    slopes->level[inNumber[c] - 1] = phys::external_reduced_width(X, L, finalGroup->GetJ(), unit, out, integral);
  }
  for (std::size_t cf = 0; cf < out.size(); cf++) {
    std::vector<phys::FinalChannel> unit = out;
    for (std::size_t k = 0; k < unit.size(); k++) unit[k].scale = k == cf ? 1.0 : 0.0;
    slopes->finalScale[outNumber[cf] - 1] = phys::external_reduced_width(X, L, finalGroup->GetJ(), in, unit, integral);
  }
  return value;
}

complex ExternalReducedWidthOfInputs(const CompoundSystem &compound, const LevelState &state, const LevelInputs &inputs, int j,
                                     int la, const AChannel *channel, const Config &configure) {
  const std::pair<int, int> finalState = ExternalCaptureFinalState(compound, channel);
  if (finalState.first == 0) return complex(0.0, 0.0);
  const JGroup *group = compound.GetJGroup(j);
  ExternalCaptureInputs values;
  values.levelEnergy = inputs.E(j, la);
  values.finalEnergy = inputs.E(finalState.first, finalState.second);
  for (int ch = 1; ch <= group->NumChannels(); ch++) values.levelGamma.push_back(inputs.Gamma(j, la, ch));
  for (int ch = 1; ch <= compound.GetJGroup(finalState.first)->NumChannels(); ch++)
    values.finalGamma.push_back(inputs.Gamma(finalState.first, finalState.second, ch));
  values.finalSqrtNF = state.SqrtNF(finalState.first, finalState.second);
  return ExternalReducedWidthAt(compound, group, channel, finalState, values, configure);
}

complex ExternalReducedWidthOfTable(const CompoundSystem &compound, const ObservedLevelTable &table, int j, int la,
                                    const AChannel *channel, const Config &configure) {
  const std::pair<int, int> finalState = ExternalCaptureFinalState(compound, channel);
  if (finalState.first == 0) return complex(0.0, 0.0);
  const ObservedLevel &level = table.Level(j, la);
  const ObservedLevel &finalLevel = table.Level(finalState.first, finalState.second);
  ExternalCaptureInputs values;
  values.levelEnergy = level.E;
  values.finalEnergy = finalLevel.E;
  values.levelGamma = level.gamma;
  values.finalGamma = finalLevel.gamma;
  values.finalSqrtNF = finalLevel.sqrtNF;
  return ExternalReducedWidthAt(compound, compound.GetJGroup(j), channel, finalState, values, configure);
}

void FinalStateOfAnotherPair(const CompoundSystem &compound, const Config &configure, int otherPairKey, int level, int pairNum) {
  const PPair *pair = compound.GetPair(pairNum);
  configure.outStream << "ERROR: external capture into pairs " << otherPairKey
                      << " and " << pair->GetPairKey() << ": both have the final state J^pi = " << pair->GetJ(2)
                      << (pair->GetPi(2) > 0 ? "+" : "-") << " at E_x = " << pair->GetHeavyExE() << " MeV (level " << level
                      << " of that J^pi group).  A photon pair has its own final state: remove the duplicated pair, or"
                      << " give it the energy of another level." << std::endl;
}
