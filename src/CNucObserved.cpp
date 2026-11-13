// CNuc: the observed quantities of the current fit values, and one fit value set with what derives from it, for the
// model data (ModelData.cpp), which compare observed widths and energies with literature values inside a fit.
#include <exception>
#include <vector>

#include "ALevel.h"
#include "CNuc.h"
#include "CNucExternalCapture.h"
#include "Config.h"
#include "JGroup.h"

bool CNuc::ObservedOfFitValues(const Config &configure) {
  // TransformToBrune sets the bound-state normalization N_f^-1/2 of the formal basis from the Brune amplitudes; a
  // calculation reads the one FillCompoundFromParams set, so it is put back.  The other slots written (transform,
  // observed, external amplitude, W) are read only by the output and by ObservedWidths itself.
  std::vector<double> sqrtNF;
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++)
      sqrtNF.push_back(this->Group(j)->GetLevel(la)->GetSqrtNFFactor());
  bool ok = true;
  try {
    this->TransformToBrune(configure, false);
    this->ObservedWidths(configure);
  } catch (const std::exception &) {
    ok = false;
  }
  std::size_t k = 0;
  for (int j = 1; j <= this->NumJGroups(); j++)
    for (int la = 1; la <= this->Group(j)->NumLevels(); la++)
      this->Group(j)->GetLevel(la)->SetSqrtNFFactor(sqrtNF[k++]);
  return ok;
}

void CNuc::SetFitLevelParameter(int j, int la, int ch, double value) {
  JGroup *group = this->Group(j);
  ALevel *level = group->GetLevel(la);
  if (ch == 0) {
    level->SetFitE(value);
  } else {
    level->SetFitGamma(ch, value);
    level->SetECGamma(ch, value);
  }
  // What FillCompoundFromParams derives from the level's values.
  std::vector<double> g;
  for (int c = 1; c <= group->NumChannels(); c++) g.push_back(level->GetFitGamma(c));
  level->SetSqrtNFFactor(SqrtNFFactor(*this, group, level, 1.0, g));
  if (ecFinalBlocks_) this->ECFinalStatesFromFormal();
  if (captureMaskSetup_) this->FormalCaptureMasksFromParams();
}
