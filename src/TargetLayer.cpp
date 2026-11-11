// The set-up refusals of a yield on a layer (TargetLayer.h).
#include "TargetLayer.h"

#include <cmath>
#include <iostream>

#include "Config.h"
#include "Segment.h"
#include "TargetEffect.h"
#include "Unit.h"
#include "WorkingCopy.h"
#include "azure2/experiment/LayerChecks.h"
#include "azure2/experiment/LayerExit.h"

int RefuseLayerAtBeamEnergy(const TargetEffect &effect, const Unit &point, const std::string &dataFile, const Config &configure) {
  const azure2::experiment::LayerPoint where{effect.GetTargetName(), dataFile, point.GetLabEnergy()};
  const double eps = effect.StoppingAt(point.GetLabEnergy(), configure, effect.GetStoppingScaleNominal());
  configure.outStream << "ERROR: " << azure2::experiment::LayerStoppingUndefined(where, point.GetLabEnergy(), eps).ToText()
                      << std::endl;
  return -1;
}

namespace {
bool CheckPoint(const WorkingCopy &data, const Unit &point, const TargetEffect &te, const std::string &dataFile,
                std::ostream &out) {
  const int n = point.NumSubPoints();
  const double f = point.GetCMFactor();  // c.m./lab: the sub-points carry c.m. energies and c.m. stopping
  const azure2::experiment::LayerPoint where{te.GetTargetName(), dataFile, point.GetLabEnergy()};
  for (int j = 1; j <= n; j++) {
    const double s = point.GetSubPoint(data, j)->GetStoppingPower();
    if (!(s > 0.0) || !std::isfinite(s)) {
      out << "ERROR: "
          << azure2::experiment::LayerStoppingUndefined(where, point.GetSubPoint(data, j)->GetCMEnergy() / f, s / f).ToText()
          << std::endl;
      return false;
    }
  }
  const bool noSpread = !te.IsGaussian() && !te.IsStraggling();
  const double top = (noSpread && n > 0) ? point.GetSubPoint(data, 1)->GetCMEnergy() : point.GetCMEnergy();
  double eOut = 0.0;
  if (azure2::experiment::LayerExitEnergy(
          n, [&](int j) { return point.GetSubPoint(data, j)->GetCMEnergy(); },
          [&](int j) { return point.GetSubPoint(data, j)->GetStoppingPower(); }, top, te.GetLayoutDensity(), eOut))
    return true;
  out << "ERROR: "
      << azure2::experiment::LayerExitBelowModel(where, te.GetLayoutDensity(), TargetEffect::minIntegrationEnergy / f).ToText()
      << std::endl;
  return false;
}
}  // namespace

bool CheckLayerYields(WorkingCopy &data, const Config &configure) {
  for (const Segment &segment : data.Data().GetSegments()) {
    if (!segment.IsTargetEffect()) continue;
    for (int t = 0; t < segment.NumTerms(); t++) {
      const TargetEffect *effect = data.Data().GetTargetEffect(segment.TermTargetEffectNum(t));
      if (!effect || !effect->IsLayerYield()) continue;
      for (const Unit &point : segment.TermPoints(data, t))
        if (point.IsTargetEffect() && point.NumSubPoints() > 0 &&
            !CheckPoint(data, point, *effect, segment.GetDataFile(), configure.outStream))
          return false;
    }
  }
  return true;
}
