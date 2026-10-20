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
struct ThmWeightTable;
class ThmDistortion;
class ThmDwVertex;
struct ThmAngleWindow;
struct ThmCoherentBackground;
struct ThmLineshape;
struct ThmSpectatorWindow;

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
  /// Returns true if this is a THM (modified R-matrix, half-off-shell) segment.
  bool IsTHM() const;
  /// True when the normalization is profiled (eliminated analytically) rather
  /// than fitted: a THM segment carries an arbitrary overall scale (Typel &
  /// Baur 2003), so a free THM norm is set to its closed-form optimum each
  /// evaluation instead of competing with the direct-data norms as a Minuit
  /// parameter.
  bool IsProfiledNorm() const;
  /// Set the normalization to the value that minimizes this segment's chi^2 at
  /// the current model (points must already hold their fit cross sections) and
  /// return that minimum chi^2.  Used for profiled (THM) norms.
  double ProfileNormChiSquared();
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
  /// Centre of the normalization prior: the norm field of the input file, or
  /// its explicit prior_centre in <parameterSettings> when it has one.
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
  /// Centre of the energy-shift prior: the shift field of the input file, or
  /// its explicit prior_centre in <parameterSettings>.
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
  /// Move the centre of the normalization prior away from the norm field (an
  /// explicit "segment_N_norm prior_centre c" row of <parameterSettings>).
  /// The current normalization (the start value) is left alone.
  void SetNominalNorm(double norm) { dataNormNominal_ = norm; }
  /// Same for the energy-shift prior ("segment_N_energy_shift prior_centre c").
  void SetNominalEnergyShift(double shift) { energyShiftNominal_ = shift; }
  /// Attach the <thm> weight table of this segment; points created afterwards
  /// carry it.  The segment keeps the table alive.
  void SetThmWeight(std::shared_ptr<const ThmWeightTable> w) { thmWeight_ = w; }
  const ThmWeightTable *GetThmWeight() const { return thmWeight_.get(); }
  /// Attach the Coulomb line shape of this segment's THM experiment
  /// (lineshape=on) to the segment and to every point it has (and sub-points
  /// created afterwards).  The segment keeps it alive.
  void SetThmLineshape(std::shared_ptr<const ThmLineshape> l);
  const ThmLineshape *GetThmLineshape() const { return thmLineshape_.get(); }
  /// Attach the spectator-momentum window of this segment's THM experiment
  /// (ps=...), as SetThmLineshape.
  void SetThmSpectatorWindow(std::shared_ptr<const ThmSpectatorWindow> w);
  const ThmSpectatorWindow *GetThmSpectatorWindow() const { return thmWindow_.get(); }
  /// Attach the distortion factor R(E) of this segment's THM experiment
  /// (distortion=...), as SetThmLineshape.
  void SetThmDistortion(std::shared_ptr<const ThmDistortion> d);
  const ThmDistortion *GetThmDistortion() const { return thmDistortion_.get(); }
  void SetThmDwVertex(std::shared_ptr<const ThmDwVertex> v);
  const ThmDwVertex *GetThmDwVertex() const { return thmDwVertex_.get(); }
  /// Attach the angular window of this segment's THM experiment
  /// (theta=thmin-thmax), as SetThmLineshape.
  void SetThmAngleWindow(std::shared_ptr<const ThmAngleWindow> w);
  const ThmAngleWindow *GetThmAngleWindow() const { return thmAngle_.get(); }
  /// THM coherent background (cbackground=), as SetThmLineshape.
  void SetThmCoherent(std::shared_ptr<const ThmCoherentBackground> c);
  const ThmCoherentBackground *GetThmCoherent() const { return thmCoherent_.get(); }
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
  bool isTHM_;
  std::shared_ptr<const ThmWeightTable> thmWeight_;
  std::shared_ptr<const ThmLineshape> thmLineshape_;
  std::shared_ptr<const ThmSpectatorWindow> thmWindow_;
  std::shared_ptr<const ThmDistortion> thmDistortion_;
  std::shared_ptr<const ThmDwVertex> thmDwVertex_;
  std::shared_ptr<const ThmAngleWindow> thmAngle_;
  std::shared_ptr<const ThmCoherentBackground> thmCoherent_;
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
  std::string normGroup_;
  int normLeaderPos_;
  bool hasNormGroupFollowers_;
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
