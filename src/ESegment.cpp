#include <exception>
#include <functional>
#include <cmath>
#include <limits>
#include "InputFiles.h"
#include "CNuc.h"
#include "Config.h"
#include "Kinematics.h"
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
  columns_ = spec.columns;
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
}

/*!
 * Returns true if a point is with the specified angle and energy ranges of a segment, otherwise
 * returns false.
 */

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

double ESegment::ProfileNormChiSquared(const EData &data, const azure2::Results &results, double &profiledNorm) const {
  double Smm = 0.0, Smd = 0.0, Sdd = 0.0;
  const UnitRange<const EPoint> points = GetPoints(data);
  for (int i = 0; i < points.size(); i++) {
    const EPoint &pt = points[i];
    double e = CMError(i);
    if (e == 0.0) continue;
    double w = 1.0 / (e * e);
    double m = pt.GetModelValue(results);
    double d = CMValue(i);
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
  for (int i = 0; i < points.size(); i++) chi += azure2::PointResidual(data, *this, i, results, 1.0).Chi2();
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
  return own_.count;
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
  // A row inside the window takes sqrt(C_ii) of a covariance file as its uncertainty (the data value standing in for
  // the model in relvector terms); chi2 itself does not depend on it (DataCovariance.h).
  auto covarianceError = [&](azure2::experiment::DataRow &row) -> bool {
    if (!covFile) return true;
    const int fileRow = covFile->RowOf(canonicalData, row.record);
    const double c = covFile->Variance(fileRow, row.value);
    if (!(c > 0.0)) {
      configure.outStream << "ERROR: Row " << fileRow + 1 << " of covariance file " << covarianceFile_ << " (data row "
                          << row.record << " of " << infile << ") has a variance " << c << " <= 0." << std::endl;
      theData->SetCovarianceError();
      return false;
    }
    const double sigma = std::sqrt(c);
    const double rel = std::fabs(row.error - sigma) / sigma;
    if (rel > 0.01) {
      if (covFile->numSigmaMismatch++ == 0)
        covFile->firstSigmaMismatch = infile + " data row " + std::to_string(row.record) + ": error column " +
                                      std::to_string(row.error) + ", sqrt(C_ii) " + std::to_string(sigma);
      covFile->worstSigmaMismatch = std::max(covFile->worstSigmaMismatch, rel);
    }
    row.error = sigma;
    return true;
  };
  azure2::experiment::RowWindow window;
  window.minEnergy = min_e_;
  window.maxEnergy = max_e_;
  window.angles = this->IsDifferential();
  window.minAngle = min_a_;
  window.maxAngle = max_a_;
  auto rows = std::make_shared<azure2::experiment::MeasurementRows>();
  rows->angleIsLegendreOrder = this->IsCoefficients();
  rows->absolute = this->IsCoefficientAbsolute();
  azure2::experiment::MeasurementProvenance provenance;
  std::string error;
  if (!azure2::experiment::ReadMeasurementRows(in, infile, window, columns_, this->IsCoefficients(), covarianceError,
                                               *rows, provenance, error)) {
    // Not a data row, or a bad order: say exactly where, then let the caller drop the segment with its usual "Could
    // Not Fill" warning.  A covariance refusal has been reported already.
    if (!error.empty()) configure.outStream << error << std::endl;
    return -1;
  }
  this->SetMeasurementRows(*theData, rows, provenance, false);
  PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
  PPair *exitPair = theCNuc->GetPairFromKey(this->GetExitKey());  // nullptr for total capture (-1)
  for (int i = 0; i < this->NumPoints(); i++) {
    // Ratios are dimensionless: no cross-section conversion (decided by the spec: components are not added yet).
    if (!this->ConvertRowToCM(*theData, 0, i, entrancePair, exitPair, !this->IsDimensionlessRatio(), BeyondRange::Stop)) {
      // Beyond the kinematic maximum lab angle: not a NaN point, not a dropped segment (audit 2026-10-25, 31).
      const EPoint *converted = this->GetPoint(*theData, i + 1);
      configure.outStream << "ERROR: " << "segment #" << this->GetSegmentKey() << ": lab angle " << converted->GetLabAngle() << " deg at E_lab = " << converted->GetLabEnergy()
                          << " MeV: " << converted->GetLabAngleWhy() << "." << std::endl;
      theData->SetKinematicsError();
      return -1;
    }
  }
  this->KeepLoadFrame(0);
  return 0;
}

void ESegment::SetMeasurementRows(EData &data, std::shared_ptr<const azure2::experiment::MeasurementRows> rows,
                                  const azure2::experiment::MeasurementProvenance &provenance, bool extrapolation) {
  auto m = std::make_shared<azure2::experiment::Measurement>();
  m->key = this->GetSegmentKey();
  m->extrapolation = extrapolation;
  m->rows = std::move(rows);
  m->provenance = provenance;
  measurement_ = m;
  const azure2::experiment::MeasurementRows &r = *measurement_->rows;
  own_.settings = MakeOwnSettings();
  // the points: a new range at the end of the data set's units (an earlier range is left for the compaction)
  own_.first = data.NumUnits();
  own_.count = r.Size();
  own_.frame.rows.clear();
  own_.frame.rows.reserve(r.Size());
  for (int i = 0; i < r.Size(); i++) {
    data.AppendUnit(EPoint(r, i, this, 0));
    own_.frame.rows.push_back(azure2::experiment::UnconvertedRow(r.energy[i], r.angle[i], r.value[i], r.error[i],
                                                            r.binLow[i], r.binHigh[i]));
  }
}

void ESegment::SetGridRows(EData &data, const std::vector<double> &energies, const std::vector<double> &angles) {
  auto rows = std::make_shared<azure2::experiment::MeasurementRows>();
  for (double angle : angles)
    for (double energy : energies)
      rows->Add(energy, angle, 0.0, 0.1, std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), 0);
  azure2::experiment::MeasurementProvenance provenance;
  provenance.rowsRead = provenance.rowsKept = rows->Size();
  this->SetMeasurementRows(data, rows, provenance, true);
}

void ESegment::SetPrediction(std::shared_ptr<const azure2::experiment::Prediction> p) {
  auto m = std::make_shared<azure2::experiment::Measurement>(*measurement_);
  m->prediction = std::move(p);
  measurement_ = m;
}

double ESegment::CMEnergyFactor(PPair *entrancePair, PPair *exitPair) {
  return (entrancePair->GetPType() == 20 && exitPair) ? azure2::experiment::DecayCMFactor(PartitionOf(exitPair))
                                                      : azure2::experiment::BeamCMFactor(PartitionOf(entrancePair));
}

namespace {

// How the rows of a term with these pairs are converted (experiment::FrameRule).
azure2::experiment::FrameRule RuleOf(const ESegment &segment, int t, PPair *entrancePair, PPair *exitPair,
                                     bool convertValue) {
  using azure2::experiment::ExitKind;
  azure2::experiment::FrameRule r;
  r.entrance = PartitionOf(entrancePair);
  if (exitPair) {
    r.exit = PartitionOf(exitPair);
    r.exitKind = exitPair->GetPType() == 0 ? ExitKind::Particle : exitPair->GetPType() == 10 ? ExitKind::Photon
                                                                                              : ExitKind::Other;
  }
  r.elastic = segment.TermEntranceKey(t) == segment.TermExitKey(t);
  r.betaDelayed = entrancePair->GetPType() == 20 && exitPair;
  r.relativistic = Kinematics::Relativistic();
  // Lab-frame differential data only: not c.m. data (isDiff 4 / 5), phase shifts, angular distributions.
  r.labDifferential = segment.IsDifferential() && !segment.IsPhase() && !segment.IsAngularDist() &&
                      !segment.IsCMDifferential();
  r.convertValue = convertValue;
  return r;
}

}  // namespace

bool ESegment::ConvertRowToCM(EData &data, int t, int i, PPair *entrancePair, PPair *exitPair, bool convertValue,
                              BeyondRange beyond) {
  TermUnits &u = Units(t);
  EPoint *point = &data.Unit(u.first + i);
  azure2::experiment::FrameRow &row = u.frame.rows[i];
  // After a shift: from the row as converted at load, so that where the angle cannot be converted the angle and the
  // factor are the load-time ones, whatever this copy converted before (data_model.md section 1.3, H4).
  if (beyond == BeyondRange::Convert && u.loadFrame) row = u.loadFrame->rows[i];
  const azure2::experiment::MeasurementRows &rows = *measurement_->rows;
  azure2::experiment::AngleLimit lim;
  const bool ok = azure2::experiment::ConvertRow(RuleOf(*this, t, entrancePair, exitPair, convertValue),
                                                 point->GetLabEnergy(), point->GetLabAngle(), rows.value[i],
                                                 rows.error[i], rows.binLow[i], rows.binHigh[i],
                                                 beyond == BeyondRange::Stop, row, lim);
  point->NoteAngleLimit(lim);
  point->SetFrame(row);
  return ok;
}

void ESegment::ConvertRowEnergyToCM(EData &data, int t, int i, PPair *entrancePair, PPair *exitPair) {
  TermUnits &u = Units(t);
  EPoint &point = data.Unit(u.first + i);
  azure2::experiment::FrameRow &row = u.frame.rows[i];
  const azure2::experiment::MeasurementRows &rows = *measurement_->rows;
  azure2::experiment::ConvertRowEnergy(RuleOf(*this, t, entrancePair, exitPair, false), point.GetLabEnergy(),
                                       rows.binLow[i], rows.binHigh[i], row);
  point.SetFrame(row);
}

int ESegment::AddLeaf(EData &data, int entranceKey, int exitKey, std::optional<double> fixedAngle, double scaling,
                      int key) {
  Leaf leaf;
  leaf.key = key;
  leaf.entranceKey = entranceKey;
  leaf.exitKey = exitKey;
  leaf.fixedAngle = fixedAngle;
  leaf.scaling = scaling;
  // the leaf's points: copies of the segment's, a new range at the end of the data set's units
  leaf.units.first = data.NumUnits();
  leaf.units.count = own_.count;
  for (int i = 0; i < own_.count; i++) {
    if (data.Unit(own_.first + i).NumSubPoints())
      throw std::logic_error("ESegment::AddLeaf: a leaf is added after the sub-points of its segment");
    data.AppendUnit(EPoint(data.Unit(own_.first + i)));
  }
  leaf.units.frame = own_.frame;
  leaf.units.loadFrame = own_.loadFrame;
  for (int i = 0; i < leaf.units.count; i++) {
    EPoint &point = data.Unit(leaf.units.first + i);
    if (fixedAngle) {
      point.SetCMAngle(*fixedAngle);
      point.SetLabAngle(*fixedAngle);
      leaf.units.frame.rows[i].thetaCM = *fixedAngle;
    }
  }
  leaves_.push_back(std::move(leaf));
  const int t = (int)leaves_.size();
  SetLeafSettings(data, t);
  return t;
}

void ESegment::MakeYieldTotalLeaf(EData &data, int t) {
  leaves_[t - 1].yieldTotal = true;
  SetLeafSettings(data, t);
}

void ESegment::SetLeafSettings(EData &data, int t) {
  // The segment's settings as they stand when the leaf is made, with the leaf's pairs and role.
  auto settings = std::make_shared<TermSettings>(*own_.settings);
  settings->entranceKey = leaves_[t - 1].entranceKey;
  settings->exitKey = leaves_[t - 1].exitKey;
  settings->total = IsTermTotal(t);
  SetTermSettings(data, t, std::move(settings));
}

bool ESegment::IsTermTotal(int t) const {
  return (isTotal_ || IsYieldTotalLeaf(t)) && TermExitKey(t) == TermEntranceKey(t);
}

void ESegment::SetTermTargetEffectNum(int t, int n) {
  if (t == 0)
    SetTargetEffectNum(n);
  else
    leaves_[t - 1].targetEffectNum = n;
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

std::shared_ptr<TermSettings> ESegment::MakeOwnSettings() const {
  auto t = std::make_shared<TermSettings>();
  t->entranceKey = entrancekey_;
  t->exitKey = exitkey_;
  t->differential = IsDifferential();
  t->phase = IsPhase();
  t->thm = IsTHM();
  t->angularDist = IsAngularDist() || IsCoefficients();
  t->analyzingPower = IsAnalyzingPower();
  t->total = IsTotalTerm();
  t->upos = IsUPOS();
  t->coefficientAbsolute = IsCoefficientAbsolute();
  t->maxAngDistOrder = GetMaxAngDistOrder();
  t->j = GetJ();
  t->l = GetL();
  t->phaseS = GetPhaseS();
  t->secondaryDecayL = GetSecondaryDecayL();
  t->ic = GetIc();
  t->delta = GetDelta();
  t->crossSectionComponent = GetCrossSectionComponent();
  t->thmWeight = thmWeight_.get();
  t->thmLineshape = thmLineshape_.get();
  t->thmWindow = thmWindow_.get();
  t->thmDistortion = thmDistortion_.get();
  t->thmDw = thmDwVertex_.get();
  t->thmAngle = thmAngle_.get();
  t->thmCoherent = thmCoherent_.get();
  return t;
}

void ESegment::SetTermSettings(EData &data, int t, std::shared_ptr<const TermSettings> settings) {
  TermUnits &u = Units(t);
  u.settings = std::move(settings);
  for (EPoint &point : TermPoints(data, t)) {
    point.UseSettings(u.settings.get());
    for (EPoint &sub : point.GetSubPoints(data)) sub.UseSettings(u.settings.get());
  }
}

void ESegment::RefreshOwnSettings(EData &data) {
  SetTermSettings(data, 0, MakeOwnSettings());
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

void ESegment::SetExitKey(EData &data, int key) {
  exitkey_ = key;
  RefreshOwnSettings(data);
}

/*!
 * Sets the entrance key for the segment and all its points
 */
void ESegment::SetEntranceKey(EData &data, int key) {
  entrancekey_ = key;
  RefreshOwnSettings(data);
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

UnitRange<EPoint> ESegment::TermPoints(EData &data, int t) const {
  return data.Units(Units(t).first, Units(t).count);
}

UnitRange<const EPoint> ESegment::TermPoints(const EData &data, int t) const {
  return data.Units(Units(t).first, Units(t).count);
}

EPoint *ESegment::GetTermPoint(EData &data, int t, int i) const {
  return &data.Unit(Units(t).first + i - 1);
}

// Advanced segment composition methods

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
 * The frame factor of a SUM term (ESegment.h).  Each point carries the
 * lab-to-c.m. factor of its own conversion (EPoint::ConvertCrossSection,
 * ConvertCrossSectionGammas; 1 when none was done), refreshed with the energy
 * shift (UpdatePointEnergyWithShift).  Before 2026-10-26 the terms were added
 * as they were: 13C(a,n1+n2) at 6 MeV, -2.9 % at 0 deg to +5.6 % at 150 deg
 * (audit item 7, tests/sum_lab_frame).
 */
double ESegment::TermFrameFactor(int t, int i) const {
  return azure2::experiment::TermFrameFactor(Jacobian(0, i), Jacobian(t, i));
}

/*!
 * The lab-to-c.m. factor of one RATIO term at a point (ESegment.h): the term's own
 * pairs, the photon Jacobian for a gamma exit and the particle one otherwise, and
 * 1 for a ratio that is not a lab-frame differential one (c.m., angle-integrated,
 * phase shift) or a term without a proper exit pair (a total-capture term).
 */
double ESegment::RatioTermFrameFactor(int t, EPoint *point, CNuc *cnuc) const {
  if (!point || !cnuc) return 1.0;
  if (!IsDifferential() || IsCMDifferential() || IsPhase() || IsAngularDist()) return 1.0;
  if (!cnuc->IsPairKey(TermEntranceKey(t)) || !cnuc->IsPairKey(TermExitKey(t))) return 1.0;
  PPair *entrancePair = cnuc->GetPair(cnuc->GetPairNumFromKey(TermEntranceKey(t)));
  PPair *exitPair = cnuc->GetPair(cnuc->GetPairNumFromKey(TermExitKey(t)));
  if (entrancePair->GetPType() == 20) return 1.0;
  if (exitPair->GetPType() == 10) return point->CalculateCrossSectionGammaConversionFactor(entrancePair);
  if (exitPair->GetPType() != 0) return 1.0;
  return point->CalculateCrossSectionConversionFactor(entrancePair, exitPair);
}

bool ESegment::HasSubPointIntegration(const EData &data, const EPoint *point) {
  if (!point || !point->IsTargetEffect() || point->NumSubPoints() <= 0) return false;
  const TargetEffect *te = data.GetTargetEffect(point->GetTargetEffectNum());
  return te && te->IsSubPointEffect();
}

void ESegment::SetThmLineshape(EData &data, std::shared_ptr<const ThmLineshape> l) {
  thmLineshape_ = l;
  RefreshOwnSettings(data);
}

void ESegment::SetThmDistortion(EData &data, std::shared_ptr<const ThmDistortion> d) {
  thmDistortion_ = d;
  RefreshOwnSettings(data);
}

void ESegment::SetThmDwVertex(EData &data, std::shared_ptr<const ThmDwVertex> v) {
  thmDwVertex_ = v;
  RefreshOwnSettings(data);
}

void ESegment::SetThmCoherent(EData &data, std::shared_ptr<const ThmCoherentBackground> c) {
  thmCoherent_ = c;
  RefreshOwnSettings(data);
}

void ESegment::SetThmAngleWindow(EData &data, std::shared_ptr<const ThmAngleWindow> w) {
  thmAngle_ = w;
  RefreshOwnSettings(data);
}

void ESegment::SetThmSpectatorWindow(EData &data, std::shared_ptr<const ThmSpectatorWindow> w) {
  thmWindow_ = w;
  RefreshOwnSettings(data);
}
