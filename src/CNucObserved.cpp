// CNuc: the observed quantities of the current fit values (a table), and one fit value set with what derives
// from it, for the model data (ModelData.cpp), which compare observed widths and energies with literature values
// inside a fit.
#include <exception>
#include <vector>

#include "ALevel.h"
#include "CNuc.h"
#include "CNucExternalCapture.h"
#include "Config.h"
#include "JGroup.h"

ObservedLevelTable CNuc::EmptyObservedTable() const {
  std::vector<int> numLevels, numChannels;
  for (int j = 1; j <= this->NumJGroups(); j++) {
    numLevels.push_back(this->GetJGroup(j)->NumLevels());
    numChannels.push_back(this->GetJGroup(j)->NumChannels());
  }
  return ObservedLevelTable(numLevels, numChannels);
}

std::optional<ObservedLevelTable> CNuc::ObservedOfFitValues(const Config &configure, LevelState &state,
                                                          const LevelInputs &inputs) const {
  // In Park's basis TransformToBrune recomputes the shift functions, which set the bound-state normalizations
  // N_f^-1/2 with Park's J; a calculation reads the ones its own evaluation set, so they are put back.
  std::vector<double> sqrtNF;
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++)
      sqrtNF.push_back(state.SqrtNF(j, la));
  std::optional<ObservedLevelTable> table = this->EmptyObservedTable();
  try {
    this->TransformToBrune(configure, Messages::kQuiet, state, inputs, *table);
    this->ObservedWidths(configure, state, *table);
  } catch (const std::exception &) {
    table.reset();
  }
  std::size_t k = 0;
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->GetJGroup(j)->NumLevels(); la++)
      state.SetSqrtNF(j, la, sqrtNF[k++]);
  return table;
}

void CNuc::SetFitLevelParameter(LevelState &state, const LevelParameterLayout::Slot &slot, double value) const {
  const int j = slot.group, la = slot.level;
  const JGroup *group = this->GetJGroup(j);
  const ALevel *level = group->GetLevel(la);
  if (slot.kind == LevelParameterLayout::Slot::Kind::kEnergy) {
    state.SetE(j, la, value);
  } else {
    state.SetGamma(j, la, slot.channel, value);
    state.SetECGamma(j, la, slot.channel, value);
  }
  // What FillCompoundFromParams derives from the level's values.
  std::vector<double> g;
  for (int c = 1; c <= group->NumChannels(); c++) g.push_back(state.Gamma(j, la, c));
  state.SetSqrtNF(j, la, SqrtNFFactor(*this, group, level, 1.0, g));
  if (ecFinalBlocks_) this->ECFinalStatesFromFormal(state);
  if (captureMaskSetup_) this->FormalCaptureMasksFromParams(state);
}
