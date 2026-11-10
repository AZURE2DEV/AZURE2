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
  reaction_ = r;
  naturalToken_ = r.sample;
  isTotal_ = r.total;
  transmission_ = r.transmission;
  yield_ = r.yieldThickness;
  thickness_ = GetNominalThickness();
}

/*!
 * A data segment: points read from a file (Fill) and compared with the model.
 */

ESegment::ESegment(const azure2::experiment::DatasetSpec &spec) : parameterNames_(spec.names) {
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
  energyShiftSqrt_ = spec.energyShiftSqrt;
  lastEnergyShiftSqrt_ = std::numeric_limits<double>::quiet_NaN();
  energyShiftSqrtError_ = spec.energyShiftSqrtError;
  varyEnergyShiftSqrt_ = spec.varyEnergyShiftSqrt;
  shiftGroup_ = spec.shiftGroup;
  shiftLeaderPos_ = -1;
  shiftFnGroup_ = spec.shiftFunctionGroup;
  normGroup_ = spec.normGroup;
  covarianceFile_ = spec.covarianceFile;
  normLeaderPos_ = -1;
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
}

/*!
 * An extrapolation segment: no data, the points are made on a grid (EData::MakePoints).
 */

ESegment::ESegment(const azure2::experiment::ExtrapolationSpec &spec) : parameterNames_(spec.names) {
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
  energyShiftSqrt_ = 0.0;
  lastEnergyShiftSqrt_ = std::numeric_limits<double>::quiet_NaN();
  energyShiftSqrtError_ = 0.0;
  shiftGroup_ = "";
  shiftLeaderPos_ = -1;
  shiftFnGroup_ = "";
  normGroup_ = "";
  normLeaderPos_ = -1;
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
        // Ratios are dimensionless: no cross-section conversion (decided by the spec: components are not added yet
        // during Fill()).
        if (!this->ConvertPointToCM(converted, entrancePair, exitPair, configure, !this->IsDimensionlessRatio(),
                                    BeyondRange::Stop)) {
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

bool ESegment::HasSubPointIntegration(EPoint *point) {
  if (!point || !point->IsTargetEffect() || point->NumSubPoints() <= 0 || !point->GetParentData()) return false;
  TargetEffect *te = point->GetParentData()->GetTargetEffect(point->GetTargetEffectNum());
  return te && te->IsSubPointEffect();
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
