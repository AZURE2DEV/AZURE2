// The index map of the fit's parameter vector (azure2/inference/ParameterIndexMap.h).
#include "azure2/inference/ParameterIndexMap.h"

namespace azure2 {
namespace inference {

void ParamIndexMap::SizeLevelBlock(const std::vector<int> &numLevels, const std::vector<int> &numChannels) {
  eFlat_.assign(numLevels.size(), {});
  gammaFlat_.assign(numLevels.size(), {});
  for (std::size_t g = 0; g < numLevels.size(); g++) {
    eFlat_[g].assign(numLevels[g], -1);
    gammaFlat_[g].assign(numLevels[g], std::vector<int>(numChannels[g], -1));
  }
}

int ParamIndexMap::AppendLevelEnergy(int jGroup, int level) {
  const int i = Append(ParamDesc{ParamKind::LevelEnergy, jGroup, level, -1, -1, false});
  eFlat_[jGroup - 1][level - 1] = i;
  return i;
}

int ParamIndexMap::AppendAmplitude(int jGroup, int level, int channel, bool radiative) {
  const int i = Append(ParamDesc{radiative ? ParamKind::RadiativeDecay : ParamKind::Gamma, jGroup, level, channel, -1, false});
  gammaFlat_[jGroup - 1][level - 1][channel - 1] = i;
  return i;
}

int ParamIndexMap::AppendIsospin(ParamKind kind, int jGroup, int level, int channel, int nucleus) {
  ParamDesc d{kind, jGroup, level, channel, -1, false};
  d.nucleus = nucleus;
  return Append(d);
}

int ParamIndexMap::AppendAbundance(int ordinal, int sample, int nucleus) {
  ParamDesc d{ParamKind::Abundance, -1, -1, -1, -1, false};
  d.sample = sample;
  d.nucleus = nucleus;
  return abundanceIndex_[ordinal] = Append(d);
}

int ParamIndexMap::AppendNorm(int segment) {
  return normIndex_[segment] = Append(ParamDesc{ParamKind::Norm, -1, -1, -1, segment, false});
}

int ParamIndexMap::AppendEnergyShift(int segment) {
  return shiftIndex_[segment] = Append(ParamDesc{ParamKind::EnergyShift, -1, -1, -1, segment, false});
}

int ParamIndexMap::AppendEnergyShiftSqrt(int segment) {
  return sqrtShiftIndex_[segment] = Append(ParamDesc{ParamKind::EnergyShiftSqrt, -1, -1, -1, segment, false});
}

int ParamIndexMap::AppendEnergyShiftFn(int segment, int j) {
  ParamDesc d{ParamKind::EnergyShiftFn, -1, -1, -1, segment, false};
  d.fnTerm = j;
  return fnShiftIndex_[{segment, j}] = Append(d);
}

int ParamIndexMap::AppendScreening(int k) {
  ParamDesc d{ParamKind::Screening, -1, -1, -1, -1, false};
  d.screening = k;
  return screeningIndex_[k] = Append(d);
}

int ParamIndexMap::AppendDensity(int k) {
  ParamDesc d{ParamKind::Density, -1, -1, -1, -1, false};
  d.density = k;
  return densityIndex_[k] = Append(d);
}

int ParamIndexMap::AppendDoppler(int k) {
  ParamDesc d{ParamKind::Doppler, -1, -1, -1, -1, false};
  d.doppler = k;
  return dopplerIndex_[k] = Append(d);
}

int ParamIndexMap::AppendThickness(int k, int segment) {
  ParamDesc d{ParamKind::Thickness, -1, -1, -1, segment, false};
  d.thickness = k;
  return thicknessIndex_[k] = Append(d);
}

int ParamIndexMap::AppendThmCoherent(int k) { return Append(ParamDesc{ParamKind::ThmCoherent, -1, -1, k, -1, false}); }

int ParamIndexMap::AppendApparatus(int k) {
  ParamDesc d{ParamKind::Apparatus, -1, -1, -1, -1, false};
  d.apparatus = k;
  return apparatusIndex_[k] = Append(d);
}

void ParamIndexMap::Pack(const std::vector<bool> &fixed) {
  const int nFull = (int)desc_.size();
  fullToPacked_.assign(nFull, -1);
  packedToFull_.clear();
  for (int f = 0; f < nFull; f++) {
    const bool isFixed = (f < (int)fixed.size()) ? fixed[f] : false;
    desc_[f].fixed = isFixed;
    if (!isFixed) {
      fullToPacked_[f] = (int)packedToFull_.size();
      packedToFull_.push_back(f);
    }
  }
  numPacked_ = (int)packedToFull_.size();
}

}  // namespace inference
}  // namespace azure2
