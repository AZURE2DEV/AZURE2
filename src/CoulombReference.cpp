// The Coulomb leaf of a ratio to the Rutherford or Mott cross section, and the S-factor conversion of a point
// (CoulombReference.h).
#include "CoulombReference.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "CNuc.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "Segment.h"
#include "Kinematics.h"
#include "PPair.h"
#include "PointFold.h"
#include "TargetEffect.h"
#include "azure2/physics/AstrophysicalFactor.h"
#include "azure2/physics/CoulombScattering.h"

namespace azure2 {

double CoulombReferenceValue(const WorkingCopy &data, Unit *point, experiment::LeafPart part, CNuc *compound,
                             const Config &configure) {
  if (!experiment::IsCoulombReference(part))
    throw std::logic_error("CoulombReferenceValue: the leaf is not a Coulomb reference");
  PPair *pair = compound->GetPair(compound->GetPairNumFromKey(point->GetEntranceKey()));
  const physics::ChannelKinematics kin = Kinematics::Channel(pair).Physics();
  const int twoS = static_cast<int>(std::lround(2.0 * pair->GetJ(1)));
  const bool mott = part == experiment::LeafPart::Mott;
  auto at = [&](const Unit *p) {
    return mott ? physics::mott_cross_section(kin, twoS, p->GetCMEnergy(), p->GetCMAngle()).value
                : physics::rutherford_cross_section(kin, p->GetCMEnergy(), p->GetCMAngle()).value;
  };
  if (!Segment::HasSubPointIntegration(data, point)) return at(point);
  std::vector<double> sub(static_cast<std::size_t>(point->NumSubPoints()));
  for (int j = 1; j <= point->NumSubPoints(); j++) sub[j - 1] = at(point->GetSubPoint(data, j));
  double folded = FoldValues(*point, data, configure, sub);
  // Unit::Calculate blends a folded value with the bare one at a window edge outside a thin-target integration (a
  // thin target blends inside its kernel, which FoldValues applies).
  const TargetEffect *effect = data.Data().GetTargetEffect(point->GetTargetEffectNum());
  const bool thinTarget = effect->IsTargetIntegration() && !effect->IsThickTarget();
  const double blend = point->GetTargetBlendWeight();
  if (!thinTarget && blend < 1.0) folded = blend * folded + (1.0 - blend) * at(point);
  return folded;
}

double SFactorConversionValue(const Unit *point, CNuc *compound, double shape_perMeV) {
  PPair *pair = compound->GetPair(compound->GetPairNumFromKey(point->GetEntranceKey()));
  return physics::modified_sfactor_conversion(Kinematics::Channel(pair).Physics(), point->GetCMEnergy(), shape_perMeV);
}

}  // namespace azure2
