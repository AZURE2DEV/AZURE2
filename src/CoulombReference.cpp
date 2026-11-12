// The Coulomb leaf of a ratio to the Rutherford or Mott cross section, and the S-factor conversion of a point
// (CoulombReference.h).
#include "CoulombReference.h"

#include <cmath>
#include <functional>
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
  // at the point's angle and each slot's energy: an elastic slot has no angle of its own (Site::AngleFree)
  const double angle = point->GetCMAngle();
  auto at = [&](const Site *p) {
    return mott ? physics::mott_cross_section(kin, twoS, p->GetCMEnergy(), angle).value
                : physics::rutherford_cross_section(kin, p->GetCMEnergy(), angle).value;
  };
  return FoldReferenceAtPoint(data, point, configure, at);
}

double FoldReferenceAtPoint(const WorkingCopy &data, Unit *point, const Config &configure,
                            const std::function<double(const Site *)> &at) {
  if (!Segment::HasSubPointIntegration(data, point)) return at(point);
  std::vector<double> sub(static_cast<std::size_t>(point->NumSubPoints()));
  for (int j = 1; j <= point->NumSubPoints(); j++) sub[j - 1] = at(point->GetSubPoint(data, j));
  return FoldValues(*point, data, configure, sub);
}

double SFactorConversionValue(const Site *point, CNuc *compound, double shape_perMeV) {
  PPair *pair = compound->GetPair(compound->GetPairNumFromKey(point->GetEntranceKey()));
  return physics::modified_sfactor_conversion(Kinematics::Channel(pair).Physics(), point->GetCMEnergy(), shape_perMeV);
}

}  // namespace azure2
