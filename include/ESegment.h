#ifndef ESEGMENT_H
#define ESEGMENT_H

#include "EPoint.h"
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
  /// Re-apply the current energy shift to every point, undoing the previous one.
  void UpdatePointEnergiesWithShift(CNuc *theCNuc = NULL, const Config *configure = NULL);
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

  // Calculate theoretical cross section including components
  /// Theoretical value at a point, combining the components if this is a composite segment.
  double CalculateTheoreticalCrossSection(int pointIndex, CNuc *cnuc, const Config &configure, EData *edata);

 private:
  bool isdifferential_;
  bool iscmdifferential_;
  bool isphase_;
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
  bool isAdvanced_;
  int operationType_;
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
