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

std::optional<ObservedLevelTable> CNuc::ObservedOfFitValues(const Config &configure) {
  // In Park's basis TransformToBrune recomputes the shift functions, which set the bound-state normalizations
  // N_f^-1/2 with Park's J; a calculation reads the ones its own evaluation set, so they are put back.
  std::vector<double> sqrtNF;
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++)
      sqrtNF.push_back(levelState_.SqrtNF(j, la));
  std::optional<ObservedLevelTable> table = this->EmptyObservedTable();
  try {
    this->TransformToBrune(configure, Messages::kQuiet, *table);
    this->ObservedWidths(configure, *table);
  } catch (const std::exception &) {
    table.reset();
  }
  std::size_t k = 0;
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++)
      levelState_.SetSqrtNF(j, la, sqrtNF[k++]);
  return table;
}

void CNuc::SetFitLevelParameter(const LevelParameterLayout::Slot &slot, double value) {
  const int j = slot.group, la = slot.level;
  JGroup *group = this->Group(j);
  ALevel *level = group->GetLevel(la);
  if (slot.kind == LevelParameterLayout::Slot::Kind::kEnergy) {
    levelState_.SetE(j, la, value);
  } else {
    levelState_.SetGamma(j, la, slot.channel, value);
    levelState_.SetECGamma(j, la, slot.channel, value);
  }
  // What FillCompoundFromParams derives from the level's values.
  std::vector<double> g;
  for (int c = 1; c <= group->NumChannels(); c++) g.push_back(levelState_.Gamma(j, la, c));
  levelState_.SetSqrtNF(j, la, SqrtNFFactor(*this, group, level, 1.0, g));
  if (ecFinalBlocks_) this->ECFinalStatesFromFormal();
  if (captureMaskSetup_) this->FormalCaptureMasksFromParams();
}
