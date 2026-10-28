#ifndef ESEGMENT_H
#define ESEGMENT_H

#include "EPoint.h"
#include "Equation.h"
#include <vector>
#include <mutex>
#include <memory>

class EData;
class ExtrapLine;
class SegLine;
class CNuc;
class Config;

enum OperationType {
  SUM = 0,
  RATIO = 1
};

/// An AZURE data segment

/*!
 * An AZURE data segment is specified by an entrance and exit particle pair key, as well as a range of
 * energy and angle value and a data file name. The segment also contains flags specifing they type of data point it contains. The ESegment object is the container object for a vector of EData objects.
 */

class ESegment {
 public:
  /// Build a data segment from a <segmentsData> line.
  ESegment(SegLine);
  /// Build an extrapolation segment from a <segmentsTest> line -- points on a grid, no data file.
  ESegment(ExtrapLine);
  /// Does this point fall inside the segment's energy and angle ranges?
  bool IsInSegment(EPoint);
  /// Differential cross section in the lab frame?
  bool IsDifferential() const;
  /// Differential cross section in the centre-of-mass frame?
  bool IsCMDifferential() const;
  /// Phase shift?
  bool IsPhase() const;
  /// Does the segment carry target effects?
  bool IsTargetEffect() const;
  /// Is the normalization a free fit parameter?
  bool IsVaryNorm() const;
  /// Angular distribution (Legendre coefficients)?
  bool IsAngularDist() const;
  /*!
   * Vector analyzing power segment (isDiff = 7). The comparison quantity is
   * A_y rather than a cross section, so the values in the data file are
   * dimensionless and bounded by one.
   */
  bool IsAnalyzingPower() const { return isAnalyzingPower_; };
  /// Unobserved-primary, observed-secondary reaction?
  bool IsUPOS() const;
  /// Number of segments summed for total capture; 0 if this is not a total-capture segment.
  int IsTotalCapture() const;
  /// Returns the cross section component this segment is compared against (0 = full, 1 = E1 only, 2 = E2 only).
  int GetCrossSectionComponent() const;
  /// Number of points in the segment.
  int NumPoints() const;
  /// Entrance particle pair key.
  int GetEntranceKey() const;
  /// Exit particle pair key.
  int GetExitKey() const;
  /// Read the data file and keep the points inside the energy and angle ranges, converting from the lab frame where needed.
  int Fill(CNuc *, EData *, const Config &);
  /// Orbital angular momentum. Phase-shift segments only.
  int GetL() const;
  /// 1-based position of the TargetEffect in the parent EData.
  int GetTargetEffectNum() const;
  /// Segment key -- its position in the input file, counting inactive segments too.
  int GetSegmentKey() const;
  /// Highest Legendre order. Angular-distribution segments only.
  int GetMaxAngDistOrder() const;
  /// Angular momentum of the secondary decay. UPOS segments only.
  int GetSecondaryDecayL() const;
  /// Lowest energy of the segment, lab frame.
  double GetMinEnergy() const;
  /// Highest energy of the segment, lab frame.
  double GetMaxEnergy() const;
  /// Smallest angle of the segment, lab frame.
  double GetMinAngle() const;
  /// Largest angle of the segment, lab frame.
  double GetMaxAngle() const;
  /// Chi-squared contributed by this segment.
  double GetSegmentChiSquared() const;
  /// Energy step between generated points. Extrapolation segments only.
  double GetEStep() const;
  /// Angle step between generated points. Extrapolation segments only.
  double GetAStep() const;
  /// Total spin. Phase-shift segments only.
  double GetJ() const;
  /// Channel spin of a phase-shift segment ("spin <s>"), or -1 when the line gives none; see EData::Initialize.
  double GetPhaseS() const { return phaseS_; }
  /// Current normalization applied to the data.
  double GetNorm() const;
  /// Normalization as declared in the input file, the value the fit penalty pulls towards.
  double GetNominalNorm() const;
  /// Fractional systematic uncertainty on the normalization.
  double GetNormError() const;
  /// Spin of the final state. UPOS segments only.
  double GetIc() const;
  /// Multipole mixing ratio. UPOS segments only.
  double GetDelta() const;
  /// Current energy shift applied to the segment's points.
  double GetEnergyShift() const;
  /// Shift last applied, so a change can be undone before applying the new one.
  double GetLastEnergyShift() const;
  /// Energy shift as declared in the input file.
  double GetNominalEnergyShift() const;
  /// Uncertainty on the energy shift.
  double GetEnergyShiftError() const;
  /// Is the energy shift a free fit parameter?
  bool IsVaryEnergyShift() const;
  /// Coefficient b of the sqrt(E) energy-shift term, E' = E + a + b*sqrt(E/MeV) (MeV^1/2).
  double GetEnergyShiftSqrt() const;
  /// sqrt(E) coefficient last applied to the points.
  double GetLastEnergyShiftSqrt() const;
  /// sqrt(E) coefficient as declared in the input file (the penalty's reference).
  double GetNominalEnergyShiftSqrt() const;
  /// Uncertainty on the sqrt(E) coefficient (MeV^1/2).
  double GetEnergyShiftSqrtError() const;
  /// Is the sqrt(E) coefficient a free fit parameter?
  bool IsVaryEnergyShiftSqrt() const;
  /// Shift-group label ("" if none): segments with one label share one constant shift.
  std::string GetShiftGroup() const;
  /// Position (0-based, in EData's segment vector) of the group member whose shift
  /// parameter this segment uses; -1 for a segment that uses its own.
  int GetShiftLeaderPos() const;
  /// Does this segment take its constant shift from another segment of its group?
  bool IsShiftGroupFollower() const;
  /// Does another segment take its constant shift from this one?
  bool HasShiftGroupFollowers() const;
  /// Make this segment a follower of the group member at position pos: its own
  /// shift stops being a fit parameter and starts at the given value.
  void SetShiftGroupFollower(int pos, double shift);
  void SetHasShiftGroupFollowers(bool);
  /// Shift-function-group label ("" if none): segments with one label apply the shift
  /// function of the first of them but keep their own constant shifts.
  std::string GetShiftFnGroup() const;
  /// Normalization-group label ("" if none): segments with one label share one norm.
  std::string GetNormGroup() const;
  /// Position (0-based, in EData's segment vector) of the group member whose norm this
  /// segment uses; -1 for a segment that uses its own.
  int GetNormLeaderPos() const;
  /// Does this segment take its normalization from another segment of its group?
  bool IsNormGroupFollower() const;
  /// Does another segment take its normalization from this one?
  bool HasNormGroupFollowers() const;
  /// Make this segment a follower of the group member at position pos: it takes that
  /// member's norm and nominal norm, and carries no norm penalty of its own.
  void SetNormGroupFollower(int pos, double norm, double nominal);
  void SetHasNormGroupFollowers(bool);
  /// Data-covariance file of the segment ("" if none): the segments naming one file
  /// share one covariance matrix (EData::ResolveCovarianceGroups).
  const std::string &GetCovarianceFile() const { return covarianceFile_; }
  /// Index of the segment's covariance group in EData (-1 if it has none).
  int GetCovarianceGroup() const { return covGroup_; }
  void SetCovarianceGroup(int g) { covGroup_ = g; }
  /// (R^-1 r)_i of point \p pointIdx (0-based) at the last EData::CovarianceChi2Correction:
  /// the point's share of d(chi2)/d(r_i) / 2 with the correlation matrix R of its group.
  double GetCovarianceSeed(int pointIdx) const {
    return pointIdx >= 0 && pointIdx < (int)covSeed_.size() ? covSeed_[pointIdx] : 0.0;
  }
  void SetCovarianceSeed(int pointIdx, double z) {
    if ((int)covSeed_.size() <= pointIdx) covSeed_.resize(pointIdx + 1, 0.0);
    covSeed_[pointIdx] = z;
  }
  /// The whole lab-energy shift at a point's original lab energy: a + b*sqrt(E).
  double TotalEnergyShift(double labEnergy) const;
  /// Add this segment's energy-shift penalties, ((x - nominal)/sigma)^2 for each varied
  /// shift term with sigma != 0, to chi.  The one definition every chi^2 path uses; the
  /// index form (gradient, curvature, least-squares priors) is ForEachShiftPenalty in
  /// AZURECalc.cpp and must list the same terms.
  void AddShiftPenalty(double &chi) const;
  /// Shift function ("shiftfn" block): does this segment's own line declare one?  Such a
  /// segment owns the function's fit parameters (and their penalties).
  bool OwnsShiftFn() const;
  /// Expression, starting values (the penalties' nominal values), uncertainties and vary
  /// flags declared on this segment's line (owners only).
  std::string GetDeclaredShiftFn() const;
  const std::vector<double> &GetNominalShiftFnValues() const;
  const std::vector<double> &GetShiftFnErrors() const;
  bool IsVaryShiftFn(int j) const;
  /// Does this line carry a "sqrtshift" block?
  bool HasSqrtShiftBlock() const;
  /// Is a shift function applied to this segment's points (its own, or its group leader's)?
  bool HasShiftFn() const;
  /// Expression of the applied shift function ("" if none).
  std::string GetShiftFnExpression() const;
  /// Number of parameters of the applied shift function.
  int NumShiftFnParams() const;
  /// Position (0-based, in EData's segment vector) of the segment that owns the applied
  /// shift function's parameters (this segment or its group leader); -1 if none.
  int GetShiftFnSourcePos() const;
  /// Owners: index of the function's first parameter within EData's shift-function block.
  int GetShiftFnParamStart() const;
  void SetShiftFnParamStart(int);
  /// Current and last-applied parameter values of the applied shift function.
  const std::vector<double> &GetShiftFnValues() const;
  const std::vector<double> &GetLastShiftFnValues() const;
  void SetShiftFnValue(int j, double value);
  void SetLastShiftFnValues(const std::vector<double> &);
  /// Apply the shift function owned by the segment at position ownerPos (this one, or its
  /// group leader) to this segment; its parameters start at the owner's declared values.
  void InstallShiftFn(const ESegment &owner, int ownerPos, const Config &);
  /// Copy the applied shift function (expression and current values) of another segment,
  /// for the components of a composite or total-capture segment.
  void CopyShiftFnFrom(const ESegment &);
  /// f(x) alone, at original lab energy x (MeV); 0 without a shift function.
  double ShiftFnTerm(double labEnergy) const;
  /// Is this a composite segment, combining others by sum or ratio?
  bool IsAdvanced() const;
  /// Combination type as an integer: 0 sum, 1 ratio.
  int GetLegacyOperationType() const;
  /// Semicolon-separated component list, as written in the input file.
  std::string GetComponentsList() const;
  /// Path of the data file this segment reads.
  std::string GetDataFile() const;
  /// Append a point.
  void AddPoint(EPoint);
  /// Record this segment's chi-squared.
  void SetSegmentChiSquared(double);
  /// Point the segment at a TargetEffect in the parent EData.
  void SetTargetEffectNum(int);
  void SetSegmentKey(int);
  /// Set the normalization applied to the data.
  void SetNorm(double);
  /// Set the energy shift; UpdatePointEnergiesWithShift applies it to the points.
  void SetEnergyShift(double);
  void SetLastEnergyShift(double);
  /// Set the sqrt(E) coefficient; UpdatePointEnergiesWithShift applies it to the points.
  void SetEnergyShiftSqrt(double);
  void SetLastEnergyShiftSqrt(double);
  /// Re-apply the current energy shift to every point, undoing the previous one.
  void UpdatePointEnergiesWithShift(CNuc *theCNuc = NULL, const Config *configure = NULL);
  /// The same for one point (0-based index) and its sub-points.
  void UpdatePointEnergyWithShift(int i, CNuc *theCNuc, const Config *configure);
  /// Point i (0-based) and, after it, every point mapped onto it; nothing if point i is itself
  /// mapped onto an earlier point (its host does it).  The unit for parallel re-energising.
  void UpdatePointGroupEnergyWithShift(int i, CNuc *theCNuc, const Config *configure);
  /// The earlier point of this segment that point i (0-based) is mapped onto, or NULL.
  EPoint *ShiftHost(int i);
  /// Change the exit pair key.
  void SetExitKey(int);
  /// Change the entrance pair key, for this segment and every point in it.
  void SetEntranceKey(int);
  /// Number of segments to sum for total capture; 0 if not total capture.
  void SetIsTotalCapture(int);
  /// Free or fix the normalization in the fit.
  void SetVaryNorm(bool);
  void SetMinAngle(double);
  void SetMaxAngle(double);
  /// Point \p i, 1-based.
  EPoint *GetPoint(int);
  /// The points, by reference.
  std::vector<EPoint> &GetPoints();

  // Advanced segment composition methods
  /// Add a component named by its entrance and exit keys.
  void AddComponent(int entranceKey, int exitKey);
  /// Add a component as a complete segment.
  void AddComponentSegment(ESegment *componentSegment);
  /// How the components combine: SUM or RATIO.
  void SetOperationType(OperationType operation);
  /// How the components combine.
  OperationType GetOperationType() const;
  /// Does this segment combine other segments?
  bool HasComponents() const;
  /// A natural-target segment: its own term and every component is weighted
  /// by the isotopic abundance of the term's entrance pair's compound nucleus
  /// in the segment's sample (the trailing "natural <1|sample>" token of the
  /// line, resolved to a 1-based sample of the <abundances> block by
  /// EData::Fill; 0 = not a natural segment).
  bool IsNatural() const { return naturalSample_ > 0; }
  /// "total 1" on the line: every term of the segment whose exit pair is its entrance pair is the total cross
  /// section of that pair, from the optical theorem (GenMatrixFunc::CalculateBareCrossSection); a term to another
  /// exit pair (a SUM component adding capture, say) stays what it is.  Neutral entrance pairs only.
  bool IsTotal() const { return isTotal_; }
  /// An angular-coefficient data segment (isDiff 9): rows E_lab, k, value, error; the value
  /// is the Legendre coefficient a_k = A_k / A_0 of dsigma/dOmega_cm = A_0 [1 + sum_k a_k P_k]
  /// (or the absolute A_k with "absolute 1").  The points are angular-distribution points
  /// whose fit value is their own order's coefficient.
  bool IsCoefficients() const { return isCoefficients_; }
  bool IsCoefficientAbsolute() const { return coefficientsAbsolute_; }
  /// IsTotal for this term: the segment is "total" and its exit pair is its entrance pair.
  bool IsTotalTerm() const { return isTotal_ && this->GetExitKey() == this->GetEntranceKey(); }
  void ApplyTotalToPoints();
  int GetNaturalSample() const { return naturalSample_; }
  void SetNaturalSample(int sample) { naturalSample_ = sample; }
  /// The token as read ("0", "1" or a sample name), before resolution.
  const std::string &GetNaturalToken() const { return naturalToken_; }
  /// Components are never weighted themselves (the base segment does it).
  void SetNatural(bool natural) { if (!natural) naturalSample_ = 0; }
  /// The component segments, by reference.
  const std::vector<ESegment *> &GetComponentSegments() const;
  /// Drop every component.
  void ClearComponents();

  // Per-component scaling factor applied when this segment is combined as a
  // component of a parent advanced segment. Defaults to 1.0 (no scaling).
  /// Scale factor applied to this segment when it is summed into a parent composite segment.
  void SetComponentScaling(double scaling);
  /// Scale factor as a component; 1.0 unless it is one.
  double GetComponentScaling() const;
  /// The factor that takes a SUM term's c.m. cross section into the frame of
  /// the base point: the base point's lab-to-c.m. conversion factor over the
  /// term's.  A lab-frame differential segment is compared in the base pair's
  /// c.m. frame (sigma_cm = J_base sigma_lab, ESegment::Fill), while a
  /// component to another exit pair is evaluated in its own c.m. frame with its
  /// own Jacobian J_term, so the lab-frame sum  sum_i s_i sigma_cm,i / J_i
  /// reads  sum_i s_i (J_base / J_i) sigma_cm,i  there.  1 when the frames
  /// agree: c.m. or angle-integrated segments, a component of the same pair,
  /// gamma terms of one entrance pair.
  static double TermFrameFactor(const EPoint *basePoint, const EPoint *termPoint);

  /// The lab-to-c.m. factor of one RATIO term (sigma_cm = J sigma_lab) at a point,
  /// with the TERM's own entrance and exit pair and the Jacobian of its own exit
  /// type (photon or particle).  1 unless the ratio is a lab-frame differential
  /// one (the term's points were then never converted, so nothing is stored on
  /// them: the RATIO data value is dimensionless and ESegment::Fill skips the
  /// cross-section conversion).  Before 2026-10-27 the denominator's factor was
  /// computed with the BASE segment's pairs (audit item 15).
  static double RatioTermFrameFactor(EPoint *point, const ESegment *term, CNuc *cnuc);

  // Calculate theoretical cross section including components
  /// Theoretical value at a point, combining the components if this is a composite segment.
  /*!
   * `precalculated`: every non-mapped point of this segment and of its components
   * has already been Calculate()d for the current parameters (AZURECalc does
   * that in one parallel pass), so only the combination is done here.
   */
  double CalculateTheoreticalCrossSection(int pointIndex, CNuc *cnuc, const Config &configure, EData *edata,
                                          bool precalculated = false);

 private:
  bool isdifferential_;
  bool iscmdifferential_;
  bool isphase_;
  double phaseS_ = -1.0;
  /// Which cross section component to compare against: 0 = full, 1 = E1 only, 2 = E2 only.
  int crossSectionComponent_;
  bool isTargetEffect_;
  bool varyNorm_;
  bool isAngDist_;
  bool isAnalyzingPower_;
  bool isUPOS_;
  int secondaryDecayL_;
  double Ic_;
  double delta_;
  int isTotalCapture_;
  int entrancekey_;
  int exitkey_;
  int l_;
  int targetEffectNum_;
  int segmentKey_;
  int maxAngDistOrder_;
  double min_e_;
  double max_e_;
  double min_a_;
  double max_a_;
  double e_step_;
  double a_step_;
  double segment_chi_squared_;
  double j_;
  double dataNorm_;
  double dataNormNominal_;
  double dataNormError_;
  double energyShift_;
  double lastEnergyShift_;  // To track last applied energy shift
  double energyShiftNominal_;
  double energyShiftError_;
  bool varyEnergyShift_;
  double energyShiftSqrt_;
  double lastEnergyShiftSqrt_;  // To track last applied sqrt(E) coefficient
  double energyShiftSqrtNominal_;
  double energyShiftSqrtError_;
  bool varyEnergyShiftSqrt_;
  std::string shiftGroup_;
  int shiftLeaderPos_;
  bool hasShiftGroupFollowers_;
  std::string shiftFnGroup_;
  std::string normGroup_;
  int normLeaderPos_;
  bool hasNormGroupFollowers_;
  std::string covarianceFile_;
  int covGroup_ = -1;
  std::vector<double> covSeed_;
  bool hasSqrtShiftBlock_;
  // Shift function declared on this line (owner) ...
  std::string shiftFnDeclared_;
  std::vector<double> shiftFnNominal_;
  std::vector<double> shiftFnError_;
  std::vector<bool> shiftFnVary_;
  // ... and the one applied to the points (own or the group leader's).
  bool hasShiftFn_;
  Equation shiftFn_;
  std::vector<double> shiftFnValues_;
  std::vector<double> lastShiftFnValues_;
  const Config *shiftFnConfig_;
  int shiftFnSourcePos_;
  int shiftFnParamStart_;
  bool isAdvanced_;
  int operationType_;
  std::string naturalToken_ = "0";
  int naturalSample_ = 0;
  bool isTotal_ = false;  // IsTotal
  bool isCoefficients_ = false;
  bool coefficientsAbsolute_ = false;
  std::string componentsList_;
  std::string datafile_;
  std::vector<EPoint> points_;

  // Advanced segment composition
  std::vector<ESegment *> componentSegments_;
  OperationType segmentOperationType_;
  double componentScaling_;
  mutable std::shared_ptr<std::mutex> componentCalculationMutex_;
};

#endif
