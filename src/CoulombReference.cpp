// The Coulomb leaf of a ratio to the Rutherford or Mott cross section (CoulombReference.h).
#include "CoulombReference.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "EPoint.h"
#include "ESegment.h"
#include "Kinematics.h"
#include "PPair.h"
#include "TargetEffect.h"
#include "azure2/physics/CoulombScattering.h"

namespace azure2 {

double CoulombReferenceValue(EPoint *point, experiment::LeafPart part, CNuc *compound, const Config &configure) {
  if (!experiment::IsCoulombReference(part))
    throw std::logic_error("CoulombReferenceValue: the leaf is not a Coulomb reference");
  PPair *pair = compound->GetPair(compound->GetPairNumFromKey(point->GetEntranceKey()));
  const physics::ChannelKinematics kin = Kinematics::Channel(pair).Physics();
  const int twoS = static_cast<int>(std::lround(2.0 * pair->GetJ(1)));
  const bool mott = part == experiment::LeafPart::Mott;
  auto at = [&](const EPoint *p) {
    return mott ? physics::mott_cross_section(kin, twoS, p->GetCMEnergy(), p->GetCMAngle()).value
                : physics::rutherford_cross_section(kin, p->GetCMEnergy(), p->GetCMAngle()).value;
  };
  if (!ESegment::HasSubPointIntegration(point)) return at(point);
  std::vector<double> sub(static_cast<std::size_t>(point->NumSubPoints()));
  for (int j = 1; j <= point->NumSubPoints(); j++) sub[j - 1] = at(point->GetSubPoint(j));
  double folded = point->IntegrateSubPointValues(configure, sub);
  // EPoint::Calculate blends a folded value with the bare one at a window edge outside a thin-target integration (a
  // thin target blends inside its kernel, which IntegrateSubPointValues applies).
  const TargetEffect *effect = point->GetParentData()->GetTargetEffect(point->GetTargetEffectNum());
  const bool thinTarget = effect->IsTargetIntegration() && !effect->IsThickTarget();
  const double blend = point->GetTargetBlendWeight();
  if (!thinTarget && blend < 1.0) folded = blend * folded + (1.0 - blend) * at(point);
  return folded;
}

}  // namespace azure2
