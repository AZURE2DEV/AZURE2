#ifndef AZURE2_COULOMBREFERENCE_H
#define AZURE2_COULOMBREFERENCE_H

#include "azure2/experiment/Prediction.h"

class CNuc;
class Config;
class EData;
class EPoint;

namespace azure2 {

/*!
 * The value of a Coulomb leaf of a ratio to the Rutherford or Mott cross section (experiment::LeafPart::Rutherford,
 * Mott; data_kinds.md section 4.3) at a point, b/sr, c.m.: physics/CoulombScattering.h at the point's c.m. energy and
 * angle, with the engine's eta and k of the point's entrance pair (Kinematics::Channel: its masses, its charges, the
 * run's kinematics) and, for Mott, the pair's particle spin.  Under a sub-point effect (EPoint::IsTargetEffect with a
 * sub-point kernel) the reference is folded as the point's observable is: evaluated at the point's own sub-points and
 * integrated by the point's own kernel (EPoint::IntegrateSubPointValues, the one quadrature), with the edge blend
 * w fold + (1 - w) bare of EPoint::Calculate where that applies, so the model is the ratio of the two folds,
 * <dsigma/dOmega> / <sigma_C>, not the fold of the ratio.
 *
 * It reads the point's present energy and angle: an energy shift that re-energised the point (ESegment::
 * UpdatePointEnergiesWithShift) moves it as it moves the numerator; no level parameter does.  `data` is the data set
 * that holds the point (its target effect).
 */
double CoulombReferenceValue(const EData &data, EPoint *point, experiment::LeafPart part, CNuc *compound,
                             const Config &configure);

/*!
 * The factor E exp(2 pi eta(E)) (MeV) that takes the cross section of a point to its S factor (experiment::FactorRef::
 * SFactor; physics/AstrophysicalFactor.h): E the point's present c.m. energy (after an energy shift that re-energised
 * it), eta of its entrance pair with the pair's masses and charges.  Under a fold it is taken at the point's energy:
 * the model is E exp(2 pi eta(E)) <sigma>, the S factor of the folded cross section at the point.
 */
double SFactorConversionValue(const EPoint *point, CNuc *compound);

}  // namespace azure2

#endif
