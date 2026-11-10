#include <exception>
#include <functional>
#include <cmath>
#include <limits>
#include "InputFiles.h"
#include "CNuc.h"
#include "Config.h"
#include "DataLine.h"
#include "DataCovariance.h"
#include "EData.h"
#include "ESegment.h"
#include "ModelPass.h"
#include "azure2/experiment/Expression.h"
#include "ExperimentAdapter.h"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include "ThmLineshape.h"
#include "ThmDistortion.h"
#include "ThmAngular.h"
#include "ThmDwVertex.h"

/*!
 * The observable flags and the model-side options of a data set or an extrapolation (ReactionSpec).  A THM segment
 * (modified R-matrix / half-off-shell; docs/THM_IMPLEMENTATION.md) computes the same observables with the THM
 * amplitude.  The vector analyzing power is differential in the c.m. frame, so it takes the angular machinery of a
 * differential cross section although it is a dimensionless ratio.
 */

void ESegment::SetReaction(const azure2::experiment::ReactionSpec &r) {
  using azure2::experiment::Observable;
  entrancekey_ = r.entranceKey;
  exitkey_ = r.exitKey;
  isTHM_ = r.thm;
  isAnalyzingPower_ = (r.observable == Observable::AnalyzingPower);
  iscmdifferential_ = (r.observable == Observable::DifferentialCM || r.observable == Observable::AnalyzingPower);
  isdifferential_ = iscmdifferential_ || r.observable == Observable::DifferentialLab;
  isphase_ = (r.observable == Observable::PhaseShift);
  if (isphase_) {
    j_ = r.phaseJ;
    l_ = r.phaseL;
    phaseS_ = r.phaseS;
  } else {
    j_ = 0.0;
    l_ = 0;
  }
  isTotalCapture_ = (r.observable == Observable::TotalCapture) ? 1 : 0;
  // E1 / E2 data are angle-integrated capture compared only with that multipole's part of the cross section.
  if (r.observable == Observable::CrossSectionE1)
    crossSectionComponent_ = 1;
  else if (r.observable == Observable::CrossSectionE2)
    crossSectionComponent_ = 2;
  else
    crossSectionComponent_ = 0;
  isAngDist_ = false;
  maxAngDistOrder_ = 0;
  // Legendre coefficients of the angular distribution as data: the points are angular-distribution points, and
  // every order the compound allows is computed.
  isCoefficients_ = (r.observable == Observable::AngularCoefficients);
  if (isCoefficients_) maxAngDistOrder_ = 20;
  isAdvanced_ = (r.combination != azure2::experiment::ReactionSpec::Combination::None);
  operationType_ = (r.combination == azure2::experiment::ReactionSpec::Combination::Ratio) ? 1 : 0;
  components_ = r.components;
  naturalToken_ = r.sample;
  isTotal_ = r.total;
  transmission_ = r.transmission;
  yield_ = r.yieldThickness;
  thickness_ = GetNominalThickness();
}

/*!
 * A data segment: points read from a file (Fill) and compared with the model.
 */

ESegment::ESegment(const azure2::experiment::DatasetSpec &spec) {
  SetReaction(spec.reaction);
  min_e_ = spec.minE;
  max_e_ = spec.maxE;
  min_a_ = spec.minA;
  max_a_ = spec.maxA;
  coefficientsAbsolute_ = spec.absolute;
  datafile_ = spec.dataFile;
  dataNorm_ = dataNormNominal_ = spec.norm;
  dataNormError_ = spec.normError;
  energyShift_ = spec.energyShift;
  // The shift prior is centred on no shift (the spec's value is only the start), unless the spec gives a centre
  // (applied by EData::Fill once the groups are resolved, as is the normalization's).
  energyShiftNominal_ = 0.0;
  lastEnergyShift_ = std::numeric_limits<double>::quiet_NaN();  // the points start unshifted: the first fill always applies the shift
  energyShiftError_ = spec.energyShiftError;
  varyEnergyShift_ = spec.varyEnergyShift;
  energyShiftSqrt_ = energyShiftSqrtNominal_ = spec.energyShiftSqrt;
  lastEnergyShiftSqrt_ = std::numeric_limits<double>::quiet_NaN();
  energyShiftSqrtError_ = spec.energyShiftSqrtError;
  varyEnergyShiftSqrt_ = spec.varyEnergyShiftSqrt;
  shiftGroup_ = spec.shiftGroup;
  shiftLeaderPos_ = -1;
  hasShiftGroupFollowers_ = false;
  shiftFnGroup_ = spec.shiftFunctionGroup;
  normGroup_ = spec.normGroup;
  covarianceFile_ = spec.covarianceFile;
  normLeaderPos_ = -1;
  hasNormGroupFollowers_ = false;
  hasSqrtShiftBlock_ = spec.hasSqrtShift;
  shiftFnDeclared_ = spec.shiftFunction;
  shiftFnNominal_ = spec.shiftFunctionValues;
  shiftFnError_ = spec.shiftFunctionErrors;
  shiftFnVary_ = spec.shiftFunctionVary;
  hasShiftFn_ = false;
  shiftFnConfig_ = nullptr;
  shiftFnSourcePos_ = -1;
  shiftFnParamStart_ = -1;
  varyNorm_ = spec.varyNorm;
  thicknessFit_ = spec.thicknessFit;
  thicknessVary_ = spec.thicknessVary ? 1 : 0;
  thicknessSigma_ = spec.thicknessSigma;
  thicknessGroup_ = spec.thicknessGroup;
  isUPOS_ = spec.upos;
  secondaryDecayL_ = spec.secondaryDecayL;
  Ic_ = spec.secondaryFinalSpin;
  delta_ = spec.secondaryMixingRatio;
  targetEffectNum_ = 0;
  isTargetEffect_ = false;

  // Initialize advanced segment composition
  segmentOperationType_ = SUM;  // Default to SUM operation
  componentScaling_ = 1.0;
  componentCalculationMutex_ = std::make_shared<std::mutex>();
}

/*!
 * An extrapolation segment: no data, the points are made on a grid (EData::MakePoints).
 */

ESegment::ESegment(const azure2::experiment::ExtrapolationSpec &spec) {
  SetReaction(spec.reaction);
  min_e_ = spec.minE;
  max_e_ = spec.maxE;
  min_a_ = spec.minA;
  max_a_ = spec.maxA;
  if (spec.reaction.observable == azure2::experiment::Observable::AngularDistribution) {
    isAngDist_ = true;
    maxAngDistOrder_ = spec.maxAngDistOrder;
  }
  datafile_ = "";
  dataNorm_ = dataNormNominal_ = 1.;
  dataNormError_ = 0.;
  energyShift_ = energyShiftNominal_ = 0.0;
  energyShiftError_ = 0.0;
  lastEnergyShift_ = std::numeric_limits<double>::quiet_NaN();  // the points start unshifted: the first fill always applies the shift
  varyEnergyShift_ = false;
  energyShiftSqrt_ = energyShiftSqrtNominal_ = 0.0;
  lastEnergyShiftSqrt_ = std::numeric_limits<double>::quiet_NaN();
  energyShiftSqrtError_ = 0.0;
  shiftGroup_ = "";
  shiftLeaderPos_ = -1;
  hasShiftGroupFollowers_ = false;
  shiftFnGroup_ = "";
  normGroup_ = "";
  normLeaderPos_ = -1;
  hasNormGroupFollowers_ = false;
  hasSqrtShiftBlock_ = false;
  hasShiftFn_ = false;
  shiftFnConfig_ = nullptr;
  shiftFnSourcePos_ = -1;
  shiftFnParamStart_ = -1;
  varyEnergyShiftSqrt_ = false;
  varyNorm_ = false;
  // a test segment has no thickness parameter: it uses its own sample thickness (SetReaction)

  isUPOS_ = false;
  secondaryDecayL_ = 0;
  Ic_ = 0.0;
  delta_ = 0.0;
  targetEffectNum_ = 0;
  isTargetEffect_ = false;

  segmentOperationType_ = (operationType_ == 1) ? RATIO : SUM;
  componentScaling_ = 1.0;
  componentCalculationMutex_ = std::make_shared<std::mutex>();
}

/*!
 * Returns true if a point is with the specified angle and energy ranges of a segment, otherwise
 * returns false.
 */

bool ESegment::IsInSegment(EPoint point) {
  bool b = false;
  if (point.GetLabEnergy() >= this->GetMinEnergy() &&
      point.GetLabEnergy() <= this->GetMaxEnergy()) {
    if (this->IsDifferential()) {
      if (point.GetLabAngle() >= this->GetMinAngle() &&
          point.GetLabAngle() <= this->GetMaxAngle()) b = true;
    } else
      b = true;
  }
  return b;
}

/*!
 * Returns true if the segment is differential cross section, otherwise returns false.
 */

bool ESegment::IsDifferential() const {
  return isdifferential_;
}

/*!
 * Returns true if the segment is differential C.M. cross section, otherwise returns false.
 */

bool ESegment::IsCMDifferential() const {
  return iscmdifferential_;
}

/*!
 * Returns true if the segment is phase shift, otherwise returns false.
 */

bool ESegment::IsPhase() const {
  return isphase_;
}

/*!
 * Returns true if this is a THM (modified R-matrix, half-off-shell) segment,
 * whose entrance vertex uses the transfer form factor instead of the
 * penetrability and whose cross section is the HOES formula (arbitrary units).
 */

bool ESegment::IsTHM() const {
  return isTHM_;
}

/*!
 * A THM segment's overall scale is arbitrary (Typel & Baur, Ann. Phys. 305
 * (2003) 228, Sec. 6.3), so when its norm is declared free it is profiled out
 * analytically rather than carried as a Minuit parameter -- otherwise it would
 * trade against the direct-data normalizations and the joint fit could collapse
 * the direct norm.
 */

bool ESegment::IsProfiledNorm() const {
  return isTHM_ && IsVaryNorm();
}

/*!
 * With the norm n multiplying the data (residual m - n d, error n e), the
 * chi^2 = sum ((m - n d)/(n e))^2 is minimized at n* = S_mm / S_md, giving
 * chi^2_min = S_dd - S_md^2 / S_mm, where S_ab = sum a_i b_i / e_i^2.  The
 * points' model values must be in `results` (azure2::EvaluateModelPoint); n* is
 * returned in `profiledNorm`.
 */

double ESegment::ProfileNormChiSquared(const azure2::Results &results, double &profiledNorm) const {
  double Smm = 0.0, Smd = 0.0, Sdd = 0.0;
  for (const EPoint &pt : points_) {
    double e = pt.GetCMCrossSectionError();
    if (e == 0.0) continue;
    double w = 1.0 / (e * e);
    double m = pt.GetModelValue(results);
    double d = pt.GetCMCrossSection();
    Smm += m * m * w;
    Smd += m * d * w;
    Sdd += d * d * w;
  }
  if (Smd > 0.0 && Smm > 0.0) {
    profiledNorm = Smm / Smd;
    double chi = Sdd - Smd * Smd / Smm;
    return chi > 0.0 ? chi : 0.0;
  }
  // Degenerate (no positive model-data overlap): leave the scale at unity and
  // report the unscaled chi^2 so the fit still sees a finite, honest value.
  profiledNorm = 1.0;
  double chi = 0.0;
  for (const EPoint &pt : points_) chi += azure2::PointResidual(pt, results, 1.0).Chi2();
  return chi;
}

/*!
 * Returns which cross section component the segment is compared against:
 * 0 for the full cross section, 1 for the E1 component only, and 2 for the
 * E2 component only.
 */

int ESegment::GetCrossSectionComponent() const {
  return crossSectionComponent_;
}

/*!
 * Returns the number of total capture segments to be summed.  Should be zero if the segment is
 * not total capture, otherwise the parameter should be the number of segments in the sum
 * (inclusive of the current segment).
 */

int ESegment::IsTotalCapture() const {
  return isTotalCapture_;
}

/*!
 * Returns true if the segment is angular distribution.
 */

bool ESegment::IsAngularDist() const {
  return isAngDist_;
}

/*!
 * Returns true if this segment is an Unobserved Primary, Observed Secondary (UPOS) segment.
 */

bool ESegment::IsUPOS() const {
  return isUPOS_;
}

/*!
 * Returns the angular momentum of the secondary decay for a UPOS segment.
 */

int ESegment::GetSecondaryDecayL() const {
  return secondaryDecayL_;
}

/*!
 * Returns the spin of the final state for a UPOS segment.
 */

double ESegment::GetIc() const {
  return Ic_;
}

/*!
 * Returns the multipole mixing ratio for a UPOS segment.
 */

double ESegment::GetDelta() const {
  return delta_;
}

/*!
 * Returns true if the segment has a corresponding TargetEffect object, otherwise
 * returns false.
 */

bool ESegment::IsTargetEffect() const {
  return isTargetEffect_;
}

/*!
 * Returns true if the normalization parameter for the segment is to be fit, otherwise
 * returns false.
 */

bool ESegment::IsVaryNorm() const {
  return varyNorm_;
}

/*!
 * Returns the number of data point objects in the segment.
 */

int ESegment::NumPoints() const {
  return points_.size();
}

/*!
 * Returns the entrance particle pair key of the segment.
 */

int ESegment::GetEntranceKey() const {
  return entrancekey_;
}

/*!
 * Returns the exit particle pair key of the segment.
 */

int ESegment::GetExitKey() const {
  return exitkey_;
}

/*!
 * Fills the segment with points from the specified data file according to the
 * maximum and minimum energy and angle ranges.  If the data
 * in the lab frame, conversions are performed to the center of mass frame.
 */

int ESegment::Fill(CNuc *theCNuc, EData *theData, const Config &configure) {
  std::string infile = this->GetDataFile();
  std::unique_ptr<std::istream> inStream = InputFiles::Open(infile, configure.cacheInputFiles);
  std::istream &in = *inStream;
  if (!in) return -1;
  int lineNumber = 0;
  int record = 0;  // data rows read so far: the rows of a data covariance file
  std::shared_ptr<DataCovariance::File> covFile;
  std::string canonicalData;
  if (!covarianceFile_.empty()) {
    covFile = theData->CovarianceFileFor(covarianceFile_, configure);
    if (!covFile) return -1;
    canonicalData = DataCovariance::CanonicalPath(infile);
    if (covFile->FileIndex(canonicalData) < 0) {
      configure.outStream << "ERROR: Data file " << infile << " of segment #" << this->GetSegmentKey()
                          << " is not among the \"file\" lines of its covariance file " << covarianceFile_ << "."
                          << std::endl;
      theData->SetCovarianceError();
      return -1;
    }
  }
  while (true) {
    DataLine line(in);
    lineNumber += line.linesConsumed();
    if (line.atEnd()) break;
    if (!line.valid()) {
      // Not a comment and not a data row.  Say exactly where, then let the
      // caller drop the segment with its usual "Could Not Fill" warning --
      // silently skipping a mangled row would hide a corrupt file.
      configure.outStream << "ERROR: Cannot parse line " << lineNumber << " of " << infile << ": \"" << line.raw()
                          << "\"" << std::endl
                          << "       A data line needs at least four numeric columns (energy, angle, value,"
                          << " uncertainty). Lines beginning with '#' and blank lines are ignored." << std::endl;
      return -1;
    }
    record++;
    {
      EPoint NewEPoint(line, this);
      NewEPoint.SetFileRecord(record);
      if (this->IsCoefficients()) {
        // the angle column is the Legendre order k
        const double k = line.angle();
        if (k < 0.0 || k > 20.0 || std::fabs(k - std::round(k)) > 1.0e-9) {
          configure.outStream << "ERROR: line " << lineNumber << " of " << infile << ": an angular-coefficient"
                              << " segment's second column is the Legendre order k (an integer, 0-20), got " << k
                              << "." << std::endl;
          return -1;
        }
        NewEPoint.SetCoefficientOrder((int)std::round(k));
        NewEPoint.SetCoefficientAbsolute(this->IsCoefficientAbsolute());
      }
      if (this->IsInSegment(NewEPoint) && covFile) {
        // The point's uncertainty becomes sqrt(C_ii) (the data value standing in for the
        // model in relvector terms); chi2 itself does not depend on it (DataCovariance.h).
        const int row = covFile->RowOf(canonicalData, record);
        const double c = covFile->Variance(row, line.crossSection());
        if (!(c > 0.0)) {
          configure.outStream << "ERROR: Row " << row + 1 << " of covariance file " << covarianceFile_ << " (data row "
                              << record << " of " << infile << ") has a variance " << c << " <= 0." << std::endl;
          theData->SetCovarianceError();
          return -1;
        }
        const double sigma = std::sqrt(c);
        const double rel = std::fabs(line.error() - sigma) / sigma;
        if (rel > 0.01) {
          if (covFile->numSigmaMismatch++ == 0)
            covFile->firstSigmaMismatch = infile + " data row " + std::to_string(record) + ": error column " +
                                          std::to_string(line.error()) + ", sqrt(C_ii) " + std::to_string(sigma);
          covFile->worstSigmaMismatch = std::max(covFile->worstSigmaMismatch, rel);
        }
        NewEPoint.SetDataError(sigma);
      }
      if (this->IsInSegment(NewEPoint)) {
        this->AddPoint(NewEPoint);
        PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
        PPair *exitPair = theCNuc->GetPairFromKey(this->GetExitKey());  // nullptr for total capture (-1)
        EPoint *converted = this->GetPoint(this->NumPoints());
        converted->SetParentData(theData);
        // Ratios are dimensionless: no cross-section conversion.  IsAdvanced() and GetLegacyOperationType() since
        // components aren't added yet during Fill().
        const bool isRatio = this->IsAdvanced() && this->GetLegacyOperationType() == 1;
        if (!this->ConvertPointToCM(converted, entrancePair, exitPair, configure, !isRatio, BeyondRange::Stop)) {
          // Beyond the kinematic maximum lab angle: not a NaN point, not a dropped segment (audit 2026-10-25, 31).
          configure.outStream << "ERROR: " << "segment #" << this->GetSegmentKey() << ": lab angle " << converted->GetLabAngle() << " deg at E_lab = " << converted->GetLabEnergy()
                              << " MeV: " << converted->GetLabAngleWhy() << "." << std::endl;
          theData->SetKinematicsError();
          return -1;
        }
      }
    }
  }
  return 0;
}

double ESegment::CMEnergyFactor(PPair *entrancePair, PPair *exitPair) {
  return (entrancePair->GetPType() == 20 && exitPair) ? azure2::experiment::DecayCMFactor(PartitionOf(exitPair))
                                                      : azure2::experiment::BeamCMFactor(PartitionOf(entrancePair));
}

void ESegment::ConvertPointEnergyToCM(EPoint *point, PPair *entrancePair, PPair *exitPair) {
  if (entrancePair->GetPType() == 20 && exitPair)
    point->ConvertDecayEnergy(exitPair);
  else
    point->ConvertLabEnergy(entrancePair);
}

bool ESegment::ConvertPointToCM(EPoint *point, PPair *entrancePair, PPair *exitPair, const Config &configure,
                                bool convertValue, BeyondRange beyond) const {
  ConvertPointEnergyToCM(point, entrancePair, exitPair);
  // Lab-frame differential data only: not c.m. data (isDiff 4 / 5), phase shifts, angular distributions.
  const bool labDifferential = this->IsDifferential() && !this->IsPhase() && !this->IsAngularDist() && !this->IsCMDifferential();
  if (!exitPair || !labDifferential) return true;
  // Convert: from the state of the point as loaded (EPoint::KeepLoadFrame), so that where the angle cannot be
  // converted the angle and the factor are the load-time ones, whatever this copy converted before (H4).
  if (beyond == BeyondRange::Convert) point->RestoreLoadFrame();
  bool inRange = true;
  if (exitPair->GetPType() == 0) {
    inRange = (this->GetEntranceKey() == this->GetExitKey()) ? point->ConvertLabAngle(entrancePair)
                                                             : point->ConvertLabAngle(entrancePair, exitPair, configure);
    if (!inRange && beyond == BeyondRange::Stop) return false;
    if (convertValue) point->ConvertCrossSection(entrancePair, exitPair);
  }
  if (exitPair->GetPType() == 10) {
    point->ConvertLabAngleGammas(entrancePair);
    if (convertValue) point->ConvertCrossSectionGammas(entrancePair);
  }
  if (beyond == BeyondRange::Stop) point->KeepLoadFrame();
  return inRange;
}

/*!
 * Returns the orbital angular momentum value for the segment.  Applies only if the segment is phase shift.
 */

int ESegment::GetL() const {
  return l_;
}

/*!
 * Returns the position of the corresponding TargetEffect object in the vector
 *  of the parent EData object.
 */

int ESegment::GetTargetEffectNum() const {
  return targetEffectNum_;
}

/*!
 * Returns the segment key for the current ESegment object.  The segment key is the order of
 * the segment specified in the input file, INCLUDING non-active segments.
 */

int ESegment::GetSegmentKey() const {
  return segmentKey_;
}

/*!
 * Returns the maximum polynomial order if segment is angular distribution.
 */

int ESegment::GetMaxAngDistOrder() const {
  return maxAngDistOrder_;
}

/*!
 * Returns the minimum energy of the segment (lab frame).
 */

double ESegment::GetMinEnergy() const {
  return min_e_;
}

/*!
 * Returns the maximum energy of the segment (lab frame).
 */

double ESegment::GetMaxEnergy() const {
  return max_e_;
}

/*!
 * Returns the minimum angle of the segment (lab frame).
 */

double ESegment::GetMinAngle() const {
  return min_a_;
}

/*!
 * Returns the maximum angle of the segment (lab frame).
 */

double ESegment::GetMaxAngle() const {
  return max_a_;
}

/*!
 * Returns the chi-squared value due to the normalization for each segment after the fitting process.
 */

// double ESegment::GetNormChiSquared() const {
//   return norm_chi_squared_;
// }

/*!
 * Returns the total spin of the segment.  Only applies if the segment is phase shift.
 */

double ESegment::GetJ() const {
  return j_;
}

/*!
 * Returns the normalization parameter for the data segment.
 */

double ESegment::GetNorm() const {
  return dataNorm_;
}

/*!
 * Returns the nominal normalization parameter for the data segment.
 */

double ESegment::GetNominalNorm() const {
  return dataNormNominal_;
}

/*!
 * Returns the normalization error for the data segment.
 */

double ESegment::GetNormError() const {
  return dataNormError_;
}

/*!
 * Returns the name of the data file from which to read.
 */

std::string ESegment::GetDataFile() const {
  return datafile_;
}

/*!
 * Returns the energy shift value for the data segment.
 */

double ESegment::GetEnergyShift() const {
  return energyShift_;
}

/*!
 * Returns the last applied energy shift value for the data segment.
 */
double ESegment::GetLastEnergyShift() const {
  return lastEnergyShift_;
}

/*!
 * Returns the nominal energy shift value for the data segment.
 */

double ESegment::GetNominalEnergyShift() const {
  return energyShiftNominal_;
}

/*!
 * Returns the energy shift error for the data segment.
 */

double ESegment::GetEnergyShiftError() const {
  return energyShiftError_;
}

/*!
 * Returns true if the energy shift is varied in the fit.
 */

bool ESegment::IsVaryEnergyShift() const {
  return varyEnergyShift_;
}

/*!
 * Returns true if this is an advanced segment (sum/ratio).
 */

bool ESegment::IsAdvanced() const {
  return isAdvanced_;
}

/*!
 * Returns the operation type (0 for sum, 1 for ratio).
 */

int ESegment::GetLegacyOperationType() const {
  return operationType_;
}

/*!
 * Sets the energy shift value for the data segment.
 */

void ESegment::SetEnergyShift(double energyShift) {
  energyShift_ = energyShift;
}

/*!
 * Sets the last applied energy shift value for the data segment.
 */
void ESegment::SetLastEnergyShift(double lastEnergyShift) {
  lastEnergyShift_ = lastEnergyShift;
}

/*!
 * Returns the coefficient of the sqrt(E) energy-shift term currently applied.
 */

double ESegment::GetEnergyShiftSqrt() const {
  return energyShiftSqrt_;
}

/*!
 * Returns the sqrt(E) coefficient last applied to the points.
 */

double ESegment::GetLastEnergyShiftSqrt() const {
  return lastEnergyShiftSqrt_;
}

/*!
 * Returns the sqrt(E) coefficient as declared in the input file.
 */

double ESegment::GetNominalEnergyShiftSqrt() const {
  return energyShiftSqrtNominal_;
}

/*!
 * Returns the uncertainty of the sqrt(E) coefficient.
 */

double ESegment::GetEnergyShiftSqrtError() const {
  return energyShiftSqrtError_;
}

/*!
 * Returns true if the sqrt(E) coefficient is a fit parameter.
 */

std::string ESegment::GetShiftGroup() const { return shiftGroup_; }

int ESegment::GetShiftLeaderPos() const { return shiftLeaderPos_; }

bool ESegment::IsShiftGroupFollower() const { return shiftLeaderPos_ >= 0; }

bool ESegment::HasShiftGroupFollowers() const { return hasShiftGroupFollowers_; }

std::string ESegment::GetShiftFnGroup() const { return shiftFnGroup_; }

void ESegment::SetShiftGroupFollower(int pos, double shift) {
  shiftLeaderPos_ = pos;
  varyEnergyShift_ = false;  // the group's parameter belongs to its first member
  energyShift_ = shift;
}

void ESegment::SetHasShiftGroupFollowers(bool has) { hasShiftGroupFollowers_ = has; }

std::string ESegment::GetNormGroup() const { return normGroup_; }

int ESegment::GetNormLeaderPos() const { return normLeaderPos_; }

bool ESegment::IsNormGroupFollower() const { return normLeaderPos_ >= 0; }

bool ESegment::HasNormGroupFollowers() const { return hasNormGroupFollowers_; }

void ESegment::SetNormGroupFollower(int pos, double norm, double nominal) {
  normLeaderPos_ = pos;
  // The group's norm, its nominal value and its penalty belong to the first member.
  // varyNorm_ is kept: a follower with a varied norm keeps its (fixed) parameter, so
  // the parameter layout, and every index after it, does not depend on the groups.
  dataNorm_ = norm;
  dataNormNominal_ = nominal;
  dataNormError_ = 0.0;
}

void ESegment::SetHasNormGroupFollowers(bool has) { hasNormGroupFollowers_ = has; }

bool ESegment::IsVaryEnergyShiftSqrt() const {
  return varyEnergyShiftSqrt_;
}

/*!
 * The complete lab-energy shift of a point whose original lab energy is
 * labEnergy: the constant term plus b*sqrt(E/MeV).  The sqrt term is
 * evaluated at the point's own (unshifted) energy, so the mapping is the
 * same function E -> E + a + b*sqrt(E) for every point and every call.
 */

double ESegment::TotalEnergyShift(double labEnergy) const {
  const double shift = (energyShiftSqrt_ == 0.0) ? energyShift_ : energyShift_ + energyShiftSqrt_ * sqrt(fabs(labEnergy));
  if (!hasShiftFn_) return shift;
  return shift + shiftFn_.Evaluate(*shiftFnConfig_, labEnergy);
}

double ESegment::ShiftFnTerm(double labEnergy) const {
  return hasShiftFn_ ? shiftFn_.Evaluate(*shiftFnConfig_, labEnergy) : 0.0;
}

bool ESegment::OwnsShiftFn() const { return !shiftFnDeclared_.empty(); }
std::string ESegment::GetDeclaredShiftFn() const { return shiftFnDeclared_; }
const std::vector<double> &ESegment::GetNominalShiftFnValues() const { return shiftFnNominal_; }
const std::vector<double> &ESegment::GetShiftFnErrors() const { return shiftFnError_; }
bool ESegment::IsVaryShiftFn(int j) const { return j >= 0 && j < (int)shiftFnVary_.size() && shiftFnVary_[j]; }
bool ESegment::HasSqrtShiftBlock() const { return hasSqrtShiftBlock_; }
bool ESegment::HasShiftFn() const { return hasShiftFn_; }
std::string ESegment::GetShiftFnExpression() const { return hasShiftFn_ ? shiftFn_.GetEquation() : std::string(); }
int ESegment::NumShiftFnParams() const { return (int)shiftFnValues_.size(); }
int ESegment::GetShiftFnSourcePos() const { return shiftFnSourcePos_; }
int ESegment::GetShiftFnParamStart() const { return shiftFnParamStart_; }
void ESegment::SetShiftFnParamStart(int start) { shiftFnParamStart_ = start; }
const std::vector<double> &ESegment::GetShiftFnValues() const { return shiftFnValues_; }
const std::vector<double> &ESegment::GetLastShiftFnValues() const { return lastShiftFnValues_; }
void ESegment::SetLastShiftFnValues(const std::vector<double> &values) { lastShiftFnValues_ = values; }

void ESegment::SetShiftFnValue(int j, double value) {
  if (!hasShiftFn_ || j < 0 || j >= (int)shiftFnValues_.size()) return;
  shiftFnValues_[j] = value;
  shiftFn_.SetParameter(j, value, *shiftFnConfig_);
}

void ESegment::InstallShiftFn(const ESegment &owner, int ownerPos, const Config &configure) {
  const int n = (int)owner.shiftFnNominal_.size();
  shiftFn_ = Equation(owner.shiftFnDeclared_, n, configure);
  shiftFnConfig_ = &configure;
  hasShiftFn_ = true;
  shiftFnSourcePos_ = ownerPos;
  shiftFnValues_.assign(n, 0.0);
  for (int j = 0; j < n; j++) SetShiftFnValue(j, owner.shiftFnNominal_[j]);
  // Never "already applied": the first FillEnergyShiftsFromParams re-energises the points.
  lastShiftFnValues_.assign(n, std::numeric_limits<double>::quiet_NaN());
}

void ESegment::CopyShiftFnFrom(const ESegment &other) {
  if (!other.hasShiftFn_) return;
  if (!hasShiftFn_ || shiftFn_.GetEquation() != other.shiftFn_.GetEquation()) {
    shiftFn_ = other.shiftFn_;
    shiftFnConfig_ = other.shiftFnConfig_;
    hasShiftFn_ = true;
    shiftFnSourcePos_ = other.shiftFnSourcePos_;
    shiftFnValues_.assign(other.shiftFnValues_.size(), 0.0);
  }
  for (int j = 0; j < (int)other.shiftFnValues_.size(); j++)
    if (shiftFnValues_[j] != other.shiftFnValues_[j]) SetShiftFnValue(j, other.shiftFnValues_[j]);
}

// Every energy-shift term is penalised from ZERO (no shift): (dE/sigma)^2, (b/sigma)^2, (a_j/sigma)^2.  The value
// on the segment line is only the starting value.  Before 2026-10-25 the penalty was measured from that value, so
// a refit seeded from an earlier fit's shifts carried those shifts as its prior centre (DeBoer: "It should measure
// the penalty from 0").  The normalization penalty is unchanged (measured from the nominal norm).  On thm an
// explicit "segment_N_energy_shift prior_centre c" row of <parameterSettings> moves the constant term's centre
// to c (GetNominalEnergyShift); without one it is 0, as on dev.
void ESegment::AddShiftPenalty(double &chi) const {
  // A prior the prior covariance replaces is not added (objective.md D3).
  if (varyEnergyShift_ && energyShiftError_ != 0.0 && !shiftPriorReplaced_)
    chi += pow((energyShift_ - energyShiftNominal_) / energyShiftError_, 2.0);  // nominal: 0 or the prior_centre row
  if (varyEnergyShiftSqrt_ && energyShiftSqrtError_ != 0.0 && !shiftSqrtPriorReplaced_) chi += pow(energyShiftSqrt_ / energyShiftSqrtError_, 2.0);
  // The shift function's parameters are penalised once, on the segment that owns them.
  if (OwnsShiftFn() && hasShiftFn_)
    for (int j = 0; j < (int)shiftFnValues_.size(); j++)
      if (shiftFnVary_[j] && shiftFnError_[j] != 0.0 && !IsShiftFnPriorReplaced(j)) chi += pow(shiftFnValues_[j] / shiftFnError_[j], 2.0);
}

void ESegment::ShiftPenaltyParts(double &constant, double &sqrtTerm, double &function) const {
  constant = sqrtTerm = function = 0.0;
  if (varyEnergyShift_ && energyShiftError_ != 0.0 && !shiftPriorReplaced_)
    constant = pow((energyShift_ - energyShiftNominal_) / energyShiftError_, 2.0);
  if (varyEnergyShiftSqrt_ && energyShiftSqrtError_ != 0.0 && !shiftSqrtPriorReplaced_)
    sqrtTerm = pow(energyShiftSqrt_ / energyShiftSqrtError_, 2.0);
  if (OwnsShiftFn() && hasShiftFn_)
    for (int j = 0; j < (int)shiftFnValues_.size(); j++)
      if (shiftFnVary_[j] && shiftFnError_[j] != 0.0 && !IsShiftFnPriorReplaced(j))
        function += pow(shiftFnValues_[j] / shiftFnError_[j], 2.0);
}

void ESegment::SetEnergyShiftSqrt(double energyShiftSqrt) {
  energyShiftSqrt_ = energyShiftSqrt;
}

void ESegment::SetLastEnergyShiftSqrt(double lastEnergyShiftSqrt) {
  lastEnergyShiftSqrt_ = lastEnergyShiftSqrt;
}

namespace {

/*!
 * Applies an energy shift to a point energy.
 *
 * The only thing guarded here is that the result stays strictly positive: the
 * lab->CM conversions and the Coulomb functions are undefined at or below zero.
 * The floor is a fixed fraction of the point's own energy, so it can only ever
 * engage for a shift large enough to push the point through threshold -- for
 * every physically meaningful shift the mapping is the exact identity plus the
 * shift, and therefore continuous in it.
 *
 * This used to be an absolute 0.01 MeV floor, which silently snapped *every*
 * data point below 10 keV onto the same energy the moment a segment's shift
 * became non-zero.  Since the update is skipped entirely while the shift is
 * still zero, that turned chi^2 into a step function of the shift for any data
 * set reaching below 10 keV (3H+d, 3He+d, ...): an infinitesimal shift moved
 * chi^2 by thousands, the finite-difference gradient was meaningless, and the
 * fit either froze at its starting point or ran away.
 */
double ShiftedEnergy(double originalEnergy, double shift) {
  static const double kMinEnergyFraction = 1.e-6;
  double shifted = originalEnergy + shift;
  double floorEnergy = originalEnergy * kMinEnergyFraction;
  return (shifted < floorEnergy) ? floorEnergy : shifted;
}

}  // namespace

/*!
 * Updates the energies of all points in this segment based on the current energy shift.
 */

void ESegment::UpdatePointEnergiesWithShift(CNuc *theCNuc, const Config *configure) {
  // No lock needed since each thread works on independent EData clones
  // after fixing EData::Clone to properly remap component segment pointers

  for (int i = 0; i < NumPoints(); i++) UpdatePointEnergyWithShift(i, theCNuc, configure);
}

std::vector<double> ESegment::ShiftState() const {
  std::vector<double> state{this->GetEnergyShift(), this->GetEnergyShiftSqrt()};
  state.insert(state.end(), this->GetShiftFnValues().begin(), this->GetShiftFnValues().end());
  return state;
}

void ESegment::BeginShiftRefill(CNuc *theCNuc, const Config &configure) {
  // Tables of another system (a session that changed a radius or the Coulomb route) are dropped.
  if (!shiftTables_.empty() && !shiftTables_.front().table->SameSystem(theCNuc, configure)) shiftTables_.clear();
  shiftRows_.assign(this->NumPoints(), {});
}

void ESegment::CommitShiftTable(CNuc *theCNuc, const Config &configure) {
  if ((int)shiftRows_.size() != this->NumPoints()) return;  // no refill since BeginShiftRefill
  std::vector<std::pair<ChannelTable::Key, ChannelTable::Row>> all;
  for (auto &v : shiftRows_) all.insert(all.end(), v.begin(), v.end());
  shiftRows_.clear();
  if (all.empty()) return;
  const std::vector<double> state = ShiftState();
  for (auto it = shiftTables_.begin(); it != shiftTables_.end(); ++it)
    if (it->state == state) {
      shiftTables_.erase(it);
      break;
    }
  shiftTables_.insert(shiftTables_.begin(), ShiftTable{state, ChannelTable::FromRows(theCNuc, configure, std::move(all))});
  if (shiftTables_.size() > kShiftTables) shiftTables_.resize(kShiftTables);
}

ChannelRowSource ESegment::ShiftRowSource(int i) {
  ChannelRowSource rows;
  for (const ShiftTable &t : shiftTables_) rows.tables.push_back(t.table.get());
  if (i >= 0 && i < (int)shiftRows_.size()) rows.used = &shiftRows_[i];
  return rows;
}

/*!
 * The earlier point of this segment that point i (0-based) is mapped onto
 * (EData::MapData), or NULL -- also for a point mapped onto a point of another
 * segment, which this segment's shift step re-energises on its own.
 */
EPoint *ESegment::ShiftHost(int i) {
  if (i < 0 || i >= NumPoints()) return nullptr;
  EPoint *point = GetPoint(i + 1);
  if (!point->IsMapped()) return nullptr;
  EnergyMap map = point->GetMap();
  if (map.point < 1 || map.point > i) return nullptr;
  EPoint *candidate = GetPoint(map.point);
  if (candidate->NumSubPoints() != point->NumSubPoints()) return nullptr;
  for (int m = 0, n = candidate->PlanMapped().size(); m < n; m++)
    if (candidate->PlanMapped()[m] == point) return candidate;
  return nullptr;
}

/*!
 * Re-energises point i and then every point mapped onto it, on the calling
 * thread.  A mapped point copies its energy-dependent values from its host
 * (UpdatePointEnergyWithShift), so the host must be finished first: the
 * parallel re-energise loops (EData::FillEnergyShiftsFromParams, the shift
 * columns of AZURECalc::ResidualJacobian) hand out these groups, never single
 * mapped points.  A point with a host returns at once.
 */
void ESegment::UpdatePointGroupEnergyWithShift(int i, CNuc *theCNuc, const Config *configure) {
  if (ShiftHost(i)) return;
  UpdatePointEnergyWithShift(i, theCNuc, configure);
  EPoint *point = GetPoint(i + 1);
  // A point mapped from another segment (EData::MapData) is not this
  // segment's to move: it is re-energised with its own segment, on its own.
  std::less<const EPoint *> before;
  const EPoint *first = &points_.front(), *last = &points_.back();
  for (int m = 0, n = point->PlanMapped().size(); m < n; m++) {
    const EPoint *mapped = point->PlanMapped()[m];
    if (before(mapped, first) || before(last, mapped)) continue;
    const int j = (int)(mapped - first);
    if (j > i && j < NumPoints() && ShiftHost(j) == point) UpdatePointEnergyWithShift(j, theCNuc, configure);
  }
}

/*!
 * Applies the segment's energy shift to point i (0-based) and its sub-points
 * and recalculates their energy-dependent values.  Writes only that point, so
 * different points may be updated concurrently -- except that a mapped point
 * reads its host (ShiftHost), so it must run after its host and not
 * concurrently with it: UpdatePointGroupEnergyWithShift.
 */
void ESegment::UpdatePointEnergyWithShift(int i, CNuc *theCNuc, const Config *configure) {
  // An ordinary segment's energies are positive and must stay so: points at
  // E <= 0 are left alone and a shift is floored above zero (ShiftedEnergy).
  // A THM (HOES) segment has points and folding sub-points below the entrance
  // threshold, and all of them shift: leaving the sub-points at E <= 0 behind
  // (as the guard did) tore the fold of every point within its reach of
  // threshold, and a point below threshold never moved at all.
  const bool thm = this->IsTHM();
  {
    EPoint *point = GetPoint(i + 1);
    // A point mapped onto an earlier point of this segment (EData::MapData)
    // is at the same energy as its host, which has already been moved (index
    // order in UpdatePointEnergiesWithShift; host first within a group in
    // UpdatePointGroupEnergyWithShift): its energy-dependent values (and its
    // sub-points') are copied from the host instead of recomputed.  Only the
    // angle-dependent Coulomb amplitude is its own.  CopyEDependentValuesFrom
    // refuses (and the values are recomputed) if the host is not at this
    // energy.
    EPoint *host = ShiftHost(i);
    if (point && (thm || point->GetOriginalEnergy() > 0)) {
      double originalEnergy = point->GetOriginalEnergy();
      double shiftedEnergy = thm ? originalEnergy + TotalEnergyShift(originalEnergy)
                                 : ShiftedEnergy(originalEnergy, TotalEnergyShift(originalEnergy));

      // Set the shifted energy
      point->SetLabEnergy(shiftedEnergy);

      // Recalculate energy dependent values
      if (theCNuc && configure) {
        PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
        PPair *exitPair = theCNuc->GetPairFromKey(this->GetExitKey());  // nullptr for total capture (-1)

        // The same conversion as ESegment::Fill (ratios are dimensionless: no cross-section conversion)
        const bool isRatio = this->IsAdvanced() && this->GetLegacyOperationType() == 1;
        if (!this->ConvertPointToCM(point, entrancePair, exitPair, *configure, !isRatio, BeyondRange::Convert)) {
          // An energy shift moved the point past the kinematic maximum lab angle mid-fit: the c.m. angle keeps
          // its value as loaded and the run says so (it cannot stop here; audit 2026-10-25, item 31).
          configure->outStream << "WARNING: segment #" << this->GetSegmentKey() << ": the energy shift moved the point at lab angle "
                               << point->GetLabAngle() << " deg to E_lab = " << point->GetLabEnergy()
                               << " MeV, beyond its kinematic maximum " << point->GetLabAngleMax()
                               << " deg; its c.m. angle is kept at its value as loaded." << std::endl;
        }

        if (this->IsDifferential() || this->IsCMDifferential()) {
          point->ClearLegendrePolynomials();
          // Get TargetEffect with Q-coefficients if applicable
          TargetEffect *effect = NULL;
          EData *parentData = point->GetParentData();
          if (parentData && this->IsTargetEffect()) {
            TargetEffect *te = parentData->GetTargetEffect(this->GetTargetEffectNum());
            if (te && te->IsQCoefficients()) {
              effect = te;
            }
          }
          point->CalcLegendreP(configure->maxLOrder, theCNuc, effect);
        }
        // The sub-points are moved and recalculated below, one by one.
        const ChannelRowSource rows = ShiftRowSource(i);
        if (host && point->CopyEDependentValuesFrom(*host))
          point->CalcCoulombAmplitude(theCNuc, false);
        else
          point->RecalcEDependentValues(theCNuc, *configure, false, &rows);

        // IMPORTANT: Also apply energy shift to all subpoints (used in convolution/target integration)
        for (int j = 1; j <= point->NumSubPoints(); j++) {
          EPoint *subPoint = point->GetSubPoint(j);
          if (subPoint && (thm || subPoint->GetOriginalEnergy() > 0)) {
            double subOriginalEnergy = subPoint->GetOriginalEnergy();
            // Sub-points hold c.m. energies and are the depth / resolution samples
            // of their parent point: they move RIGIDLY with it, by the parent's
            // own lab shift converted to the c.m. frame.  A shift is a calibration
            // of the measured beam energy; the target thickness, the beam profile
            // and the depth grid are physical and keep their width.  Before
            // 2026-10-26 each sub-point was moved by the shift at its own energy,
            // which for a sqrt(E) shift or a shift function stretched the window by
            // (1 + ds/dE) -- a yield factor the fit could trade against b -- and put
            // the resonance-anchored lattice (EData::ReplaceSubPointGrids, which
            // offsets by the parent's shift) off a narrow resonance by s' x depth
            // (audit item 11, tests/shift_sqrt_subpoints).  A constant shift was
            // never affected.
            // A beta-delayed segment (entrance ptype 20) carries the light decay
            // particle's lab energy; c.m. = lab*(M1+M2)/M2 of the EXIT pair
            // (EPoint::ConvertDecayEnergy), not the entrance-pair ratio.
            const double cmFactor = CMEnergyFactor(entrancePair, exitPair);
            // A THM segment's shift is not floored (its energies may lie below threshold), as for the point.
            double energyShiftCM = cmFactor * TotalEnergyShift(originalEnergy);
            double subShiftedEnergy = thm ? subOriginalEnergy + energyShiftCM
                                          : ShiftedEnergy(subOriginalEnergy, energyShiftCM);

            // Set the shifted energy for subpoint
            // Must set both CM and Lab energy since CalcLegendreP uses GetLabEnergy() for Q-coefficients
            subPoint->SetCMEnergy(subShiftedEnergy);
            subPoint->SetLabEnergy(subShiftedEnergy);

            // Recalculate energy dependent values for subpoint
            if (entrancePair->GetPType() == 20) {
              // The c.m. energy was set above. ConvertDecayEnergy would read the lab
              // slot (now holding the c.m. value) and multiply by (M1+M2)/M2 again,
              // pushing every convolution sub-point of a beta-delayed spectrum up by
              // 4/3 and corrupting the folded yield; only the excitation energy needs updating.
              subPoint->ConvertExcitationEnergy(exitPair);
            } else if (this->IsCMDifferential()) {
              // subPoint->ConvertLabEnergy(entrancePair);
            } else if (!this->IsCMDifferential()) {
              // subPoint->ConvertLabEnergy(entrancePair);
            }

            // A sub-point is built at its parent's c.m. angle (EData::SetUpTargetEffects, EPoint(angle, energy,
            // parent) puts that angle in both slots) and carries no data of its own, so the lab-to-c.m. angle and
            // cross-section conversions of the parent (above) must not be applied to it.  Before 2026-10-26 they
            // were, which took the parent's c.m. angle for a lab angle and moved every sub-point by the full
            // lab-c.m. difference (2-5 deg for p + light nuclei, 18 deg for alpha + 12C): a lab-frame differential
            // segment with target/resolution integration and any active shift (nonzero, varied, shiftgroup,
            // sqrtshift, shiftfn) fitted at the wrong angles, and its shift derivatives were wrong with it
            // (16O(p,p) with a 0.1 keV convolution: chi2 13,484 -> 23,097 at a 1e-7 MeV shift).  The parent's
            // c.m. angle at the shifted energy is the right value; tests/shift_subpoint_angle.
            if (this->IsDifferential() && !this->IsCMDifferential()) {
              subPoint->SetLabAngle(point->GetCMAngle());
              subPoint->SetCMAngle(point->GetCMAngle());
            }

            if (this->IsDifferential() || this->IsCMDifferential()) {
              subPoint->ClearLegendrePolynomials();
              // Get TargetEffect with Q-coefficients if applicable
              // Note: Use point->GetParentData() since subpoints don't have parentData_ set
              TargetEffect *subEffect = NULL;
              EData *parentData = point->GetParentData();
              if (parentData && this->IsTargetEffect()) {
                TargetEffect *te = parentData->GetTargetEffect(this->GetTargetEffectNum());
                if (te && te->IsQCoefficients()) {
                  subEffect = te;
                }
              }
              subPoint->CalcLegendreP(configure->maxLOrder, theCNuc, subEffect);
            }
            if (host && subPoint->CopyEDependentValuesFrom(*host->GetSubPoint(j)))
              subPoint->CalcCoulombAmplitude(theCNuc);
            else
              subPoint->RecalcEDependentValues(theCNuc, *configure, true, &rows);
          }
        }
      }
    }
  }
}

/*!
 * Adds a data point to the segment.
 */

void ESegment::AddPoint(EPoint point) {
  // The observable is a property of the segment; stamp it on the point so the
  // calculation does not have to look back up.
  point.SetIsAnalyzingPower(this->IsAnalyzingPower());
  point.SetIsTotal(this->IsTotalTerm());
  points_.push_back(point);
}

/*!
 * Stamps IsTotalTerm on every point and sub-point: a component made from a "total" base segment
 * (EData::CreateComponentSegment) is a total only if its exit pair is its entrance pair.
 */
void ESegment::ApplyTotalToPoints() {
  const bool total = this->IsTotalTerm();
  for (EPoint &point : points_) {
    point.SetIsTotal(total);
    for (int i = 1; i <= point.NumSubPoints(); i++) point.GetSubPoint(i)->SetIsTotal(total);
  }
}

/*!
 * Sets the position of the corresponding TargetEffect object in the vector
 * of the parent EData object.
 */

void ESegment::SetTargetEffectNum(int targetEffectNum) {
  targetEffectNum_ = targetEffectNum;
  isTargetEffect_ = true;
}

/*!
 * Sets the segment key for the current ESegment object.
 */

void ESegment::SetSegmentKey(int segmentKey) {
  segmentKey_ = segmentKey;
}

/*!
 * Sets the normalization parameter for the segment.
 */

void ESegment::SetNorm(double norm) {
  dataNorm_ = norm;
}

/*!
 * Sets the exit pair key to the given value.
 */

void ESegment::SetExitKey(int key) {
  exitkey_ = key;
  for (int i = 1; i <= this->NumPoints(); i++)
    this->GetPoint(i)->SetExitKey(key);
}

/*!
 * Sets the entrance key for the segment and all its points
 */
void ESegment::SetEntranceKey(int key) {
  entrancekey_ = key;
  for (int i = 1; i <= this->NumPoints(); i++)
    this->GetPoint(i)->SetEntranceKey(key);
}

/*!
 * Sets the number of total capture segments to be summed.  Should be zero if the segment is
 * not total capture, otherwise the parameter should be the number of segments in the sum
 * (inclusive of the current segment).
 */

void ESegment::SetIsTotalCapture(int num) {
  isTotalCapture_ = num;
}

/*!
 * Set the flag determining if the normalization is varied.
 */

void ESegment::SetVaryNorm(bool varyNorm) {
  varyNorm_ = varyNorm;
}

/*!
 * Sets the minimum angle for the segment.
 */
void ESegment::SetMinAngle(double minAngle) {
  min_a_ = minAngle;
}

/*!
 * Sets the maximum angle for the segment.
 */
void ESegment::SetMaxAngle(double maxAngle) {
  max_a_ = maxAngle;
}

/*!
 * Returns a pointer to the data point object specified by a position in the EPoint vector.
 */

EPoint *ESegment::GetPoint(int pointNum) {
  EPoint *b = &points_[pointNum - 1];
  return b;
}

/*!
 * Returns a reference to the vector of EPoint objects.
 */

std::vector<EPoint> &ESegment::GetPoints() {
  return points_;
}

// Advanced segment composition methods

/*!
 * Adds a component segment defined by entrance and exit keys
 */
void ESegment::AddComponent(int entranceKey, int exitKey) {
  // This method is now deprecated - use AddComponentSegment instead
  // Components should be full ESegment copies, not just entrance/exit keys
}

/*!
 * Adds a component segment (full ESegment copy)
 */
void ESegment::AddComponentSegment(ESegment *componentSegment) {
  if (componentSegment) {
    componentSegments_.push_back(componentSegment);
  }
}

/*!
 * Sets the operation type for combining components (SUM or RATIO)
 */
void ESegment::SetOperationType(OperationType operation) {
  segmentOperationType_ = operation;
}

/*!
 * Gets the operation type for combining components
 */
OperationType ESegment::GetOperationType() const {
  return segmentOperationType_;
}

/*!
 * Returns true if this segment has additional components
 */
bool ESegment::HasComponents() const {
  return !componentSegments_.empty();
}

/*!
 * Returns a reference to the component segments vector
 */
const std::vector<ESegment *> &ESegment::GetComponentSegments() const {
  return componentSegments_;
}

/*!
 * Clears all components from the segment
 */
void ESegment::ClearComponents() {
  // Clean up component segments (they are owned by EData, so just clear pointers)
  componentSegments_.clear();
}

/*!
 * Sets the per-component scaling factor applied to this segment when it is
 * combined as a component of a parent advanced (SUM) segment.
 */
void ESegment::SetComponentScaling(double scaling) {
  componentScaling_ = scaling;
}

/*!
 * Returns the per-component scaling factor. Defaults to 1.0 for segments that
 * are not used as a scaled component.
 */
double ESegment::GetComponentScaling() const {
  return componentScaling_;
}

/*!
 * The frame factor of a SUM term (ESegment.h).  Each point carries the
 * lab-to-c.m. factor of its own conversion (EPoint::ConvertCrossSection,
 * ConvertCrossSectionGammas; 1 when none was done), refreshed with the energy
 * shift (UpdatePointEnergyWithShift).  Before 2026-10-26 the terms were added
 * as they were: 13C(a,n1+n2) at 6 MeV, -2.9 % at 0 deg to +5.6 % at 150 deg
 * (audit item 7, tests/sum_lab_frame).
 */
double ESegment::TermFrameFactor(const EPoint *basePoint, const EPoint *termPoint) {
  if (!basePoint || !termPoint) return 1.0;
  return azure2::experiment::TermFrameFactor(basePoint->GetCrossSectionKinFactor(), termPoint->GetCrossSectionKinFactor());
}

/*!
 * The lab-to-c.m. factor of one RATIO term at a point (ESegment.h): the term's own
 * pairs, the photon Jacobian for a gamma exit and the particle one otherwise, and
 * 1 for a ratio that is not a lab-frame differential one (c.m., angle-integrated,
 * phase shift) or a term without a proper exit pair (a total-capture term).
 */
double ESegment::RatioTermFrameFactor(EPoint *point, const ESegment *term, CNuc *cnuc) {
  if (!point || !term || !cnuc) return 1.0;
  if (!term->IsDifferential() || term->IsCMDifferential() || term->IsPhase() || term->IsAngularDist()) return 1.0;
  if (!cnuc->IsPairKey(term->GetEntranceKey()) || !cnuc->IsPairKey(term->GetExitKey())) return 1.0;
  PPair *entrancePair = cnuc->GetPair(cnuc->GetPairNumFromKey(term->GetEntranceKey()));
  PPair *exitPair = cnuc->GetPair(cnuc->GetPairNumFromKey(term->GetExitKey()));
  if (entrancePair->GetPType() == 20) return 1.0;
  if (exitPair->GetPType() == 10) return point->CalculateCrossSectionGammaConversionFactor(entrancePair);
  if (exitPair->GetPType() != 0) return 1.0;
  return point->CalculateCrossSectionConversionFactor(entrancePair, exitPair);
}

/*!
 * Transmission of a sample n atoms/barn thick, T = exp(-n s), and dT/ds (azure2::experiment::node::Transmission).
 */
double ESegment::Transmission(double n, double s, double *dTds, double *dTdn) {
  return azure2::experiment::node::Transmission(n, s, dTds, dTdn);
}

/*!
 * SAMMY's self-shielded yield Y0 = (s / st)(1 - exp(-n st)) and its partial derivatives
 * (azure2::experiment::node::YieldY0).
 */
double ESegment::YieldY0(double n, double s, double st, double *dYds, double *dYdst, double *dYdn) {
  return azure2::experiment::node::YieldY0(n, s, st, dYds, dYdst, dYdn);
}

bool ESegment::HasSubPointIntegration(EPoint *point) {
  if (!point || !point->IsTargetEffect() || point->NumSubPoints() <= 0 || !point->GetParentData()) return false;
  TargetEffect *te = point->GetParentData()->GetTargetEffect(point->GetTargetEffectNum());
  return te && te->IsSubPointEffect();
}

bool ESegment::TransformTerms(int pointIndex, CNuc *cnuc, std::vector<TransformTerm> &terms) {
  terms.clear();
  EPoint *basePoint = GetPoint(pointIndex + 1);
  if (!basePoint) return false;
  const bool natural = this->IsNatural();
  auto nucleusOf = [&](int entranceKey) { return cnuc->GetPair(cnuc->GetPairNumFromKey(entranceKey))->GetNucleus(); };
  const int baseNucleus = nucleusOf(this->GetEntranceKey());
  terms.push_back({basePoint, this, natural ? cnuc->GetAbundance(naturalSample_, baseNucleus) : 1.0, 1.0, baseNucleus,
                   false, crossSectionComponent_});
  for (ESegment *comp : componentSegments_) {
    if (!comp) continue;
    EPoint *cp = (pointIndex < comp->NumPoints()) ? comp->GetPoint(pointIndex + 1) : nullptr;
    if (!cp) return false;
    const int nucleus = nucleusOf(comp->GetEntranceKey());
    if (comp->IsYieldTotalTerm()) {
      terms.push_back({cp, comp, 1.0, 1.0, nucleus, true, 0});
      continue;
    }
    const double scale = comp->GetComponentScaling() * TermFrameFactor(basePoint, cp);
    const double w = natural ? scale * cnuc->GetAbundance(naturalSample_, nucleus) : scale;
    terms.push_back({cp, comp, w, scale, nucleus, false, comp->GetCrossSectionComponent()});
  }
  return true;
}

void ESegment::SetThmLineshape(std::shared_ptr<const ThmLineshape> l) {
  thmLineshape_ = l;
  for (EPoint &point : points_) point.SetThmLineshape(l.get());
}

void ESegment::SetThmDistortion(std::shared_ptr<const ThmDistortion> d) {
  thmDistortion_ = d;
  for (EPoint &point : points_) point.SetThmDistortion(d.get());
}

void ESegment::SetThmDwVertex(std::shared_ptr<const ThmDwVertex> v) {
  thmDwVertex_ = v;
  for (EPoint &point : points_) point.SetThmDwVertex(v.get());
}

void ESegment::SetThmCoherent(std::shared_ptr<const ThmCoherentBackground> c) {
  thmCoherent_ = c;
  for (EPoint &point : points_) point.SetThmCoherent(c.get());
}

void ESegment::SetThmAngleWindow(std::shared_ptr<const ThmAngleWindow> w) {
  thmAngle_ = w;
  for (EPoint &point : points_) point.SetThmAngleWindow(w.get());
}

void ESegment::SetThmSpectatorWindow(std::shared_ptr<const ThmSpectatorWindow> w) {
  thmWindow_ = w;
  for (EPoint &point : points_) point.SetThmSpectatorWindow(w.get());
}
