#ifndef SEGMENT_H
#define SEGMENT_H

#include "Unit.h"
#include "Equation.h"
#include "SegmentState.h"
#include <vector>
#include <mutex>
#include <memory>
#include <optional>
#include "azure2/experiment/Frames.h"
#include "azure2/experiment/Measurement.h"
#include "azure2/experiment/SetupSpecs.h"

namespace azure2 {
namespace experiment {
struct Prediction;
}  // namespace experiment
}  // namespace azure2

class WorkingCopy;
class CNuc;
class Config;
struct ThmWeightTable;
class ThmDistortion;
class ThmDwVertex;
struct ThmAngleWindow;
struct ThmCoherentBackground;
struct ThmLineshape;
struct ThmSpectatorWindow;

class Segment;
typedef std::vector<Segment>::const_iterator SegmentIterator;

enum OperationType {
  SUM = 0,
  RATIO = 1
};

/// A data set of the project: its measurement, its settings and the terms of its prediction (its points are units of
/// the data set that holds it, WorkingCopy::Units).

class Segment {
 public:
  /// Build a data segment from its spec (experiment/'s SetupSpecs.h).
  explicit Segment(const azure2::experiment::DatasetSpec &);
  /// Build an extrapolation segment from its spec -- points on a grid, no data file.
  explicit Segment(const azure2::experiment::ExtrapolationSpec &);
  /// Differential cross section in the lab frame?
  bool IsDifferential() const;
  /// Differential cross section in the centre-of-mass frame?
  bool IsCMDifferential() const;
  /// Phase shift?
  bool IsPhase() const;
  /// Returns true if this is a THM (modified R-matrix, half-off-shell) segment.
  bool IsTHM() const;
  /*!
   * The lowest c.m. energy (MeV) of the pair `entrance` that a fold of this data set reaches: the integration floor
   * TargetEffect::minIntegrationEnergy, or for half-off-shell data, which reach below the entrance threshold, the
   * compound system's E_x = minIntegrationEnergy (THM.md section 1): minIntegrationEnergy - (S + E_x of the pair).
   */
  double IntegrationFloor(const PPair *entrance) const;
  /// Whether an energy shift of this data set keeps its energies positive (Segment::UpdatePointEnergyWithShift): every
  /// data set but half-off-shell data, whose energies lie below the entrance threshold as well and all move.
  bool ShiftIsFloored() const { return !IsTHM(); }
  /// True when the normalization is profiled (eliminated analytically) rather
  /// than fitted: a THM segment carries an arbitrary overall scale (Typel &
  /// Baur 2003), so a free THM norm is set to its closed-form optimum each
  /// evaluation instead of competing with the direct-data norms as a Minuit
  /// parameter.
  bool IsProfiledNorm() const;
  /// The normalization that minimizes this segment's chi^2 at the model values in `results` (profiledNorm) and that
  /// minimum chi^2 (returned).  Used for profiled (THM) norms; the segment's own norm is not changed.
  double ProfileNormChiSquared(const WorkingCopy &data, const azure2::Results &results, double &profiledNorm) const;
  /// The segment's position (0-based) among the data set's segments: its state in a working copy (WorkingCopy::State).
  int Pos() const { return pos_; }
  void SetPos(int pos) { pos_ = pos; }
  /// Does the segment carry target effects?
  bool IsTargetEffect() const;
  /// Is the normalization a free fit parameter?
  bool IsVaryNorm() const;
  /*!
   * The diagonal priors of this segment's parameters that the prior covariance (<config> priorCovariance=) replaces:
   * its multivariate Gaussian takes their place (docs/upgrade/objective.md D3), set once at set-up
   * (ReplaceDiagonalPriors, Objective.h).
   */
  /// The three parts of AddShiftPenalty separately: the constant shift, the sqrt(E) coefficient, the shift function.
  void ShiftPenaltyParts(const SegmentState &state, double &constant, double &sqrtTerm, double &function) const;
  void ReplaceNormPrior() { normPriorReplaced_ = true; }
  bool IsNormPriorReplaced() const { return normPriorReplaced_; }
  void ReplaceShiftPrior() { shiftPriorReplaced_ = true; }
  bool IsShiftPriorReplaced() const { return shiftPriorReplaced_; }
  void ReplaceShiftSqrtPrior() { shiftSqrtPriorReplaced_ = true; }
  bool IsShiftSqrtPriorReplaced() const { return shiftSqrtPriorReplaced_; }
  void ReplaceShiftFnPrior(int j) {
    if (j >= 0) {
      if ((int)shiftFnPriorReplaced_.size() <= j) shiftFnPriorReplaced_.resize(j + 1, false);
      shiftFnPriorReplaced_[j] = true;
    }
  }
  bool IsShiftFnPriorReplaced(int j) const { return j >= 0 && j < (int)shiftFnPriorReplaced_.size() && shiftFnPriorReplaced_[j]; }
  void ReplaceThicknessPrior() { thicknessPriorReplaced_ = true; }
  bool IsThicknessPriorReplaced() const { return thicknessPriorReplaced_; }
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
  /// Read the data file and keep the rows inside the energy and angle ranges (the segment's measurement), one point
  /// per row of the working copy `data` (set-up), converted from the lab frame where needed.
  int Fill(CNuc *, WorkingCopy *, const Config &);
  /// What the segment measured: its rows (shared, never written after set-up), their provenance, its prediction.
  const azure2::experiment::Measurement &GetMeasurement() const { return *measurement_; }
  bool HasMeasurement() const { return measurement_ != nullptr; }
  /// Make `rows` the segment's measurement and one point per row (a data file's rows, or the grid of an
  /// extrapolation), in the frame as read; ConvertRowToCM converts them.
  void SetMeasurementRows(WorkingCopy &data, std::shared_ptr<const azure2::experiment::MeasurementRows> rows,
                          const azure2::experiment::MeasurementProvenance &provenance, bool extrapolation);
  /// The grid of an extrapolation as its measurement: one row per (angle, energy), energies varying fastest, value
  /// 0 and uncertainty 0.1 (lab MeV, deg).
  void SetGridRows(WorkingCopy &data, const std::vector<double> &energies, const std::vector<double> &angles);
  /// Row i (0-based) of the measurement as read: the value, its uncertainty (after a covariance file), its data-file
  /// record.
  double LabValue(int i) const { return measurement_->rows->value[i]; }
  double LabError(int i) const { return measurement_->rows->error[i]; }
  int FileRecord(int i) const { return measurement_->rows->record[i]; }
  /// dOmega_lab/dOmega_cm of row i (0-based) at the present energy shift of the segment's `state` in a working copy
  /// (1 where the value is not converted: c.m. data, angle-integrated data, a ratio); the Jacobian of term t (the
  /// leaf's own exit pair) for a SUM or RATIO frame factor.
  static double Jacobian(const SegmentState &state, int t, int i) { return state.terms[t].frame.rows[i].jacobian; }
  /// The value and uncertainty of row i in the c.m. frame at the present shift, J y and J dy with J = Jacobian(state,
  /// 0, i): computed for display (the written c.m. data columns, the API's c.m. data), never stored.
  double CMValue(const SegmentState &state, int i) const { return LabValue(i) * Jacobian(state, 0, i); }
  double CMError(const SegmentState &state, int i) const { return LabError(i) * Jacobian(state, 0, i); }
  /// What ConvertRowToCM does with a lab angle beyond the kinematic maximum: stop (at load: the value is not
  /// converted), or convert the value at the load-time angle (an energy shift moved the point there mid-fit).
  enum class BeyondRange { Stop, Convert };
  /// The lab -> c.m. conversion of row i (0-based) of term t at the point's present lab energy
  /// (docs/upgrade/data_structure.md section 1.2, experiment::ConvertRow): the one sequence of the load (Stop), of an
  /// energy shift (Convert) and of the leaves, written into the term's frame and point.  With convertValue (false for
  /// a ratio: dimensionless) the value and its uncertainty carry dOmega_lab/dOmega_cm.  False where the lab angle is
  /// beyond the kinematic maximum, with the maximum and the reason in `limit`; the caller reports it.
  bool ConvertRowToCM(WorkingCopy &data, int t, int i, PPair *entrancePair, PPair *exitPair, bool convertValue,
                      BeyondRange beyond, azure2::experiment::AngleLimit &limit) const;
  /// How the rows of term t with these pairs are converted (experiment::FrameRule, the rule of ConvertRowToCM; the
  /// value converted).  experiment::AngleIsLab of it: the term's angle is given in the lab.
  azure2::experiment::FrameRule TermFrameRule(int t, PPair *entrancePair, PPair *exitPair) const;
  /// E_c.m., E_x and the c.m. bin edges of row i of term t alone (the leaf grids are laid out before the leaves'
  /// angles are converted).
  void ConvertRowEnergyToCM(WorkingCopy &data, int t, int i, PPair *entrancePair, PPair *exitPair) const;
  /// The frame of term t as it stands now in `data` becomes the frame an energy shift starts from
  /// (BeyondRange::Convert); set-up.
  void KeepLoadFrame(const WorkingCopy &data, int t);

  // ---- The terms of the segment's prediction (docs/upgrade/data_structure.md section 1.5).  Term 0 is the
  // segment's own; term k >= 1 is its k-th leaf: a component of a SUM, RATIO or total-capture segment, or the total
  // cross section of a yield segment.  A leaf is evaluated at the rows of the segment's measurement with its own
  // pairs; its points are calculation units that hold no data, and it takes every setting (shift, flags, THM) from
  // the segment.
  int NumTerms() const { return 1 + (int)leaves_.size(); }
  /// The points of term t: units [TermFirstUnit, + NumTermPoints) of `data`, a working copy of the data set that
  /// holds the segment (WorkingCopy::Units; every function here that takes `data` means one).  One point per row of the
  /// measurement in every term.
  UnitRange<Unit> TermPoints(WorkingCopy &data, int t) const;
  UnitRange<const Unit> TermPoints(const WorkingCopy &data, int t) const;
  int NumTermPoints(int) const { return NumPoints(); }
  int TermFirstUnit(const WorkingCopy &data, int t) const;
  /// Point i (1-based) of term t.
  Unit *GetTermPoint(WorkingCopy &data, int t, int i) const;
  /// What term t's points compute (shared by them and their sub-points; immutable once set up).
  const std::shared_ptr<const TermSettings> &TermSettingsOf(int t) const { return Fixed(t).settings; }
  int TermEntranceKey(int t) const { return t ? leaves_[t - 1].entranceKey : entrancekey_; }
  int TermExitKey(int t) const { return t ? leaves_[t - 1].exitKey : exitkey_; }
  /// The key of the term in caches and the intEC file: the segment key, or the leaf's (negative, from the data set's
  /// counter; never a data set of the file).
  int TermKey(int t) const { return t ? leaves_[t - 1].key : segmentKey_; }
  int TermTargetEffectNum(int t) const { return t ? leaves_[t - 1].targetEffectNum : targetEffectNum_; }
  void SetTermTargetEffectNum(int t, int n);
  /// The total cross section of the entrance pair (optical theorem): the segment is "total" (or the leaf is a yield
  /// segment's total term) and the term's exit pair is its entrance pair.
  bool IsTermTotal(int t) const;
  /// The total-cross-section leaf of a yield segment: it enters Y0 as s_t and is not part of the segment's sum.
  bool IsYieldTotalLeaf(int t) const { return t > 0 && leaves_[t - 1].yieldTotal; }
  /// The scale factor of leaf t in the segment's sum (1 for term 0).
  double LeafScaling(int t) const { return t ? leaves_[t - 1].scaling : 1.0; }
  /// Add a leaf with these pairs, a fixed lab angle (deg; a RATIO denominator) or none, a scale factor and its key:
  /// its points are the segment's (at the rows, with the frame at load) with the leaf's pairs.  Returns its term.
  int AddLeaf(WorkingCopy &data, int entranceKey, int exitKey, std::optional<double> fixedAngle, double scaling, int key);
  /// Make leaf t the total-cross-section term of a yield segment.
  void MakeYieldTotalLeaf(WorkingCopy &data, int t);
  /// E_c.m. / E_lab of the segment's points, non-relativistic: (M1 + M2) / M2 of the decay pair for a beta-delayed
  /// point, M2 / (M1 + M2) of the entrance pair otherwise (the rigid shift of sub-points, the target-effect grids).
  static double CMEnergyFactor(PPair *entrancePair, PPair *exitPair);
  /// Orbital angular momentum. Phase-shift segments only.
  int GetL() const;
  /// 1-based position of the TargetEffect in the parent WorkingCopy.
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
  /// Total spin. Phase-shift segments only.
  double GetJ() const;
  /// Channel spin of a phase-shift segment ("spin <s>"), or -1 when the line gives none; see WorkingCopy::Initialize.
  double GetPhaseS() const { return phaseS_; }
  /// The normalization, energy shifts and sample thickness a working copy starts from (WorkingCopy::AddSegment): the
  /// spec's values (a group follower: its leader's).
  double StartNorm() const { return startNorm_; }
  double StartEnergyShift() const { return startEnergyShift_; }
  double StartEnergyShiftSqrt() const { return startEnergyShiftSqrt_; }
  double StartThickness() const { return startThickness_; }
  /// Centre of the normalization prior: the norm field of the input file, or
  /// its explicit prior_centre in <parameterSettings> when it has one.
  double GetNominalNorm() const;
  /// Fractional systematic uncertainty on the normalization.
  double GetNormError() const;
  /// Spin of the final state. UPOS segments only.
  double GetIc() const;
  /// Multipole mixing ratio. UPOS segments only.
  double GetDelta() const;
  /// Centre of the energy-shift prior: 0 (dev 0dd28e5; the shift field of the input file is only the
  /// start), or its explicit prior_centre in <parameterSettings>.
  double GetNominalEnergyShift() const;
  /// Uncertainty on the energy shift.
  double GetEnergyShiftError() const;
  /// Is the energy shift a free fit parameter?
  bool IsVaryEnergyShift() const;
  /// Uncertainty on the sqrt(E) coefficient (MeV^1/2).
  double GetEnergyShiftSqrtError() const;
  /// The centre of the sqrt(E) term's prior (MeV^1/2): as written, 0 unless written.
  double GetEnergyShiftSqrtCentre() const { return energyShiftSqrtCentre_; }
  /// Is the sqrt(E) coefficient a free fit parameter?
  bool IsVaryEnergyShiftSqrt() const;
  /// Shift-group label ("" if none): segments with one label share one constant shift.
  std::string GetShiftGroup() const;
  /// Position (0-based, in WorkingCopy's segment vector) of the group member whose shift
  /// parameter this segment uses; -1 for a segment that uses its own.
  int GetShiftLeaderPos() const;
  /// Does this segment take its constant shift from another segment of its group?
  bool IsShiftGroupFollower() const;
  /// Make this segment a follower of the group member at position pos: its own
  /// shift stops being a fit parameter and starts at the given value.
  void SetShiftGroupFollower(int pos, double shift);
  void SetStartNorm(double norm) { startNorm_ = norm; }
  /// Shift-function-group label ("" if none): segments with one label apply the shift
  /// function of the first of them but keep their own constant shifts.
  std::string GetShiftFnGroup() const;
  /// Normalization-group label ("" if none): segments with one label share one norm.
  std::string GetNormGroup() const;
  /// Position (0-based, in WorkingCopy's segment vector) of the group member whose norm this
  /// segment uses; -1 for a segment that uses its own.
  int GetNormLeaderPos() const;
  /// Does this segment take its normalization from another segment of its group?
  bool IsNormGroupFollower() const;
  /// Make this segment a follower of the group member at position pos: it takes that
  /// member's norm and nominal norm, and carries no norm penalty of its own.
  void SetNormGroupFollower(int pos, double norm, double nominal);
  /// Data-covariance file of the segment ("" if none): the segments naming one file
  /// share one covariance matrix (WorkingCopy::ResolveCovarianceGroups).
  const std::string &GetCovarianceFile() const { return covarianceFile_; }
  /// Index of the segment's covariance group in WorkingCopy (-1 if it has none).
  int GetCovarianceGroup() const { return covGroup_; }
  void SetCovarianceGroup(int g) { covGroup_ = g; }
  /// The whole lab-energy shift of the segment's `state` at a point's original lab energy: a + b*sqrt(E) + f(E).
  double TotalEnergyShift(const SegmentState &state, double labEnergy) const;
  /// Add this segment's energy-shift penalties, ((x - nominal)/sigma)^2 for each varied
  /// shift term with sigma != 0, to chi.  The one definition every chi^2 path uses; the
  /// index form (gradient, curvature, least-squares priors) is ForEachShiftPenalty in
  /// AZURECalc.cpp and must list the same terms.
  void AddShiftPenalty(const SegmentState &state, double &chi) const;
  /// Shift function ("shiftfn" block): does this segment's own line declare one?  Such a
  /// segment owns the function's fit parameters (and their penalties).
  bool OwnsShiftFn() const;
  /// Expression, starting values (the penalties' nominal values), uncertainties and vary
  /// flags declared on this segment's line (owners only).
  std::string GetDeclaredShiftFn() const;
  const std::vector<double> &GetNominalShiftFnValues() const;
  const std::vector<double> &GetShiftFnErrors() const;
  /// The centres of the shift function's priors, one per parameter: as written, 0 unless written.
  const std::vector<double> &GetShiftFnCentres() const;
  bool IsVaryShiftFn(int j) const;
  /// Does this line carry a "sqrtshift" block?
  bool HasSqrtShiftBlock() const;
  /// Is a shift function applied to this segment's points (its own, or its group leader's)?
  bool HasShiftFn() const;
  /// Expression of the applied shift function ("" if none).
  std::string GetShiftFnExpression() const;
  /// Number of parameters of the applied shift function.
  int NumShiftFnParams() const { return numShiftFnParams_; }
  /// Position (0-based, in WorkingCopy's segment vector) of the segment that owns the applied
  /// shift function's parameters (this segment or its group leader); -1 if none.
  int GetShiftFnSourcePos() const;
  /// Owners: index of the function's first parameter within WorkingCopy's shift-function block.
  int GetShiftFnParamStart() const;
  void SetShiftFnParamStart(int);
  /// Set parameter j of the applied shift function in the segment's `state`.
  void SetShiftFnValue(SegmentState &state, int j, double value) const;
  /// Apply the shift function owned by the segment at position ownerPos (this one, or its
  /// group leader) to this segment and to its `state`; its parameters start at the owner's declared values.
  void InstallShiftFn(SegmentState &state, const Segment &owner, int ownerPos, const Config &);
  /// f(x) alone, at original lab energy x (MeV), of the segment's `state`; 0 without a shift function.
  double ShiftFnTerm(const SegmentState &state, double labEnergy) const;
  /// Is this a composite segment, combining others by sum or ratio?
  bool IsAdvanced() const;
  /// Combination type as an integer: 0 sum, 1 ratio.
  int GetLegacyOperationType() const;
  /// The terms of a composite segment after the leading one.
  const std::vector<azure2::experiment::ComponentSpec> &GetComponents() const { return reaction_.components; }
  /// The reaction, observable and model options as the spec gave them (before the set-up resolved pairs and terms).
  const azure2::experiment::ReactionSpec &GetReaction() const { return reaction_; }
  /// Is the value a ratio (a RATIO of two reactions, or a ratio to the Rutherford or Mott cross section)?  It is the
  /// same in the lab and the c.m. frame: the data are not converted at load and the lab output does not divide it.
  bool IsDimensionlessRatio() const {
    return (IsAdvanced() && GetLegacyOperationType() == 1) ||
           reaction_.coulombReference != azure2::experiment::CoulombReference::None || !reaction_.reference.Empty();
  }
  /// What the segment predicts at its points (src/Predictions.cpp, WorkingCopy::BuildPredictions): built once at set-up,
  /// immutable, shared by the working copies; null before it is built.
  const azure2::experiment::Prediction *GetPrediction() const {
    return measurement_ ? measurement_->prediction.get() : nullptr;
  }
  /// The model value is an S factor, c E sigma exp(2 pi eta) with c the scale's constant (ReactionSpec::scale; 1 for a
  /// data set, 1 / (unit factor) for an extrapolation in another unit than MeV b).
  bool ModelIsSFactor() const { return reaction_.scale.sFactor; }
  /// g (1/MeV) of a modified S factor exp(2 pi eta + g E); 0 for the S factor and every other quantity.
  double ModelSFactorShape() const { return reaction_.scale.sFactorShape_perMeV; }
  double ModelScaleFactor() const { return reaction_.scale.factor; }
  /// Give the measurement its prediction (a new measurement that shares the rows: the points read them).
  void SetPrediction(std::shared_ptr<const azure2::experiment::Prediction> p);
  /// Path of the data file this segment reads.
  std::string GetDataFile() const;
  /// Point the segment at a TargetEffect in the parent WorkingCopy.
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
  void SetThmLineshape(WorkingCopy &data, std::shared_ptr<const ThmLineshape> l);
  const ThmLineshape *GetThmLineshape() const { return thmLineshape_.get(); }
  /// Attach the spectator-momentum window of this segment's THM experiment
  /// (ps=...), as SetThmLineshape.
  void SetThmSpectatorWindow(WorkingCopy &data, std::shared_ptr<const ThmSpectatorWindow> w);
  const ThmSpectatorWindow *GetThmSpectatorWindow() const { return thmWindow_.get(); }
  /// Attach the distortion factor R(E) of this segment's THM experiment
  /// (distortion=...), as SetThmLineshape.
  void SetThmDistortion(WorkingCopy &data, std::shared_ptr<const ThmDistortion> d);
  const ThmDistortion *GetThmDistortion() const { return thmDistortion_.get(); }
  void SetThmDwVertex(WorkingCopy &data, std::shared_ptr<const ThmDwVertex> v);
  const ThmDwVertex *GetThmDwVertex() const { return thmDwVertex_.get(); }
  /// Attach the angular window of this segment's THM experiment
  /// (theta=thmin-thmax), as SetThmLineshape.
  void SetThmAngleWindow(WorkingCopy &data, std::shared_ptr<const ThmAngleWindow> w);
  const ThmAngleWindow *GetThmAngleWindow() const { return thmAngle_.get(); }
  /// THM coherent background (cbackground=), as SetThmLineshape.
  void SetThmCoherent(WorkingCopy &data, std::shared_ptr<const ThmCoherentBackground> c);
  const ThmCoherentBackground *GetThmCoherent() const { return thmCoherent_.get(); }
  /// Re-apply the current energy shift of the segment's state in `data` to every point of every term, undoing the
  /// previous one.
  void UpdatePointEnergiesWithShift(WorkingCopy &data, CNuc *theCNuc, const Config *configure) const;
  /// The same for point i (0-based) of term t: the unit of parallel re-energising (it writes that point and the slots
  /// it writes only).  A folded point on the shared mesh moves alone: its nodes stay (ShiftedFold.h).
  void UpdatePointEnergyWithShift(WorkingCopy &data, int t, int i, CNuc *theCNuc, const Config *configure) const;
  /// Change the exit pair key.
  void SetExitKey(WorkingCopy &data, int);
  /// Change the entrance pair key, for this segment and every point in it.
  void SetEntranceKey(WorkingCopy &data, int);
  /// Number of segments to sum for total capture; 0 if not total capture.
  void SetIsTotalCapture(int);
  void SetMinAngle(double);
  void SetMaxAngle(double);
  /// Point \p i (1-based) of the segment's own term, and its points.
  Unit *GetPoint(WorkingCopy &data, int i) const { return GetTermPoint(data, 0, i); }
  UnitRange<Unit> GetPoints(WorkingCopy &data) const { return TermPoints(data, 0); }
  UnitRange<const Unit> GetPoints(const WorkingCopy &data) const { return TermPoints(data, 0); }

  /// How the components combine: SUM or RATIO.
  void SetOperationType(OperationType operation);
  /// How the components combine.
  OperationType GetOperationType() const;
  /// Does this segment combine other terms (it has leaves)?
  bool HasComponents() const { return !leaves_.empty(); }
  /// A natural-target segment: its own term and every component is weighted
  /// by the isotopic abundance of the term's entrance pair's compound nucleus
  /// in the segment's sample (the trailing "natural <1|sample>" token of the
  /// line, resolved to a 1-based sample of the <abundances> block by
  /// WorkingCopy::Fill; 0 = not a natural segment).
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
  /// The Legendre coefficients of an angular-distribution point (an extrapolation of quantity legendre) as the
  /// output and the API give them: the engine's A_L / A_0 (L = 0 ... K) mapped to the segment's expansion
  /// (experiment/LegendreExpansion.h: times A_0 for an absolute one, times the factor of L) and its unit.
  std::vector<double> LegendreCoefficients(const Unit &point, const azure2::Results &results) const;
  /// "transmission <n>": the segment is the transmission T = exp(-n s) of a sample n atoms/barn thick, s the
  /// segment's (assembled) cross section in barns -- the total cross section of the entrance pair.  "yield <n>":
  /// the self-shielded capture or reaction yield Y0 = (s / s_t)(1 - exp(-n s_t)) (SAMMY's Y0, manual III.D.1), s the
  /// segment's cross section and s_t the total cross section of the entrance pair, from a hidden "total 1"
  /// component (IsYieldTotalTerm).  With a sub-point target effect (resolution, Gaussian spread) the transform is
  /// applied at every sub-point before the kernel, as SAMMY does: broadening a cross section and then
  /// exponentiating fills in the dips of a thick sample.  azure2::EvaluateModelPoint; the points
  /// keep the plain cross sections, so they can still share calculations with other segments.
  bool IsTransmission() const { return transmission_ > 0.0; }
  bool IsYield() const { return yield_ > 0.0; }
  bool IsTransformed() const { return transmission_ > 0.0 || yield_ > 0.0; }
  /// n in atoms/barn of a transmission or yield segment, 0 otherwise: the value the transform uses -- the line's n,
  /// a thickness group's leader's (WorkingCopy::ResolveThicknessGroups), or the fit parameter segment_<key>_thickness
  /// (WorkingCopy::FillThicknessFromParams): the value in the segment's `state`.
  double GetSampleThickness(const SegmentState &state) const { return IsTransformed() ? state.thickness : 0.0; }
  /// The line's own n (the prior centre of the fit form), 0 for a segment without the transform.
  double GetNominalThickness() const { return transmission_ > 0.0 ? transmission_ : yield_; }
  /// The fit form "transmission n vary sigma" / "yield n vary sigma": the line has a parameter
  /// segment_<key>_thickness (fixed for vary 0, and for a thickness-group follower), sigma its prior width (0 = none).
  bool HasThicknessFit() const { return thicknessFit_; }
  bool IsThicknessVaried() const { return thicknessFit_ && thicknessVary_ == 1 && thicknessLeaderPos_ < 0; }
  double GetThicknessSigma() const { return thicknessSigma_; }
  /// The centre of its prior (atoms/b): as written, else the start value (docs/upgrade/parameter_kinds.md section 4).
  double GetThicknessCentre() const { return thicknessCentre_; }
  /// "thicknessgroup <label>" ("" = none); the follower's leader position in WorkingCopy's segment list (-1: none).
  const std::string &GetThicknessGroup() const { return thicknessGroup_; }
  /// The names of the segment's fit parameters, as its spec gives them (WorkingCopy::FillMnParams).
  const azure2::experiment::ParameterNames &GetParameterNames() const { return parameterNames_; }
  int GetThicknessLeaderPos() const { return thicknessLeaderPos_; }
  bool IsThicknessGroupFollower() const { return thicknessLeaderPos_ >= 0; }
  void SetThicknessGroupFollower(int pos, double n) {
    thicknessLeaderPos_ = pos;
    startThickness_ = n;
  }
  /// The prior ((n - n0)/sigma)^2 of a free thickness parameter that owns its value (no follower), 0 otherwise (a
  /// fixed parameter pays no prior, docs/upgrade/objective.md D5).
  double ThicknessPenalty(const SegmentState &state) const {
    if (!IsThicknessVaried() || thicknessPriorReplaced_ || !(thicknessSigma_ > 0.0)) return 0.0;
    const double d = (state.thickness - thicknessCentre_) / thicknessSigma_;
    return d * d;
  }
  /// The value the fit compares with the data is formed by the segment, not read off its points: components, a
  /// sample's abundance weights or a transmission / yield transform.
  bool NeedsAssembly() const { return HasComponents() || IsNatural() || IsTransformed(); }
  /// True when the point, a point of `data`, is integrated over sub-points (Unit::IntegrateTargetEffect).
  static bool HasSubPointIntegration(const WorkingCopy &data, const Unit *point);
  /// IsTotal for the segment's own term: the segment is "total" and its exit pair is its entrance pair.
  bool IsTotalTerm() const { return IsTermTotal(0); }
  int GetNaturalSample() const { return naturalSample_; }
  void SetNaturalSample(int sample) { naturalSample_ = sample; }
  /// The token as read ("0", "1" or a sample name), before resolution.
  const std::string &GetNaturalToken() const { return naturalToken_; }
  /// The factor that takes a SUM term's c.m. cross section into the frame of the base point at row i: the base's
  /// lab-to-c.m. factor over the term's (experiment::TermFrameFactor).  A lab-frame differential segment is compared
  /// in the base pair's c.m. frame (sigma_cm = J_base sigma_lab), while a component to another exit pair is evaluated
  /// in its own c.m. frame with its own Jacobian J_term, so the lab-frame sum  sum_i s_i sigma_cm,i / J_i  reads
  /// sum_i s_i (J_base / J_i) sigma_cm,i  there.  1 when the frames agree.
  static double TermFrameFactor(const SegmentState &state, int t, int i);

  /// The lab-to-c.m. factor of one RATIO term (sigma_cm = J sigma_lab) at a point,
  /// with the TERM's own entrance and exit pair and the Jacobian of its own exit
  /// type (photon or particle).  1 unless the ratio is a lab-frame differential
  /// one (the term's points were then never converted, so nothing is stored on
  /// them: the RATIO data value is dimensionless and Segment::Fill skips the
  /// cross-section conversion).  Before 2026-10-27 the denominator's factor was
  /// computed with the BASE segment's pairs (audit item 15).
  double RatioTermFrameFactor(int t, Unit *point, CNuc *cnuc) const;

  /// Free-shift refills (docs/upgrade/channel_table.md section 2.1, interim until the objective has the analytic
  /// shift derivative): BeginShiftRefill, single-threaded before the parallel re-energise loop, keeps the tables of
  /// the segment's last kShiftTables shift states for UpdatePointEnergyWithShift to read (exact keys) and has every
  /// point record the rows it used; CommitShiftTable, single-threaded after the loop, turns those rows into the table
  /// of the present state (the shift values, bitwise), replacing the oldest; nothing writes a stored table.
  void BeginShiftRefill(WorkingCopy &data, int t, CNuc *, const Config &) const;
  void CommitShiftTable(WorkingCopy &data, int t, CNuc *, const Config &) const;
  static constexpr std::size_t kShiftTables = 3;

 private:
  /// The observable flags and the model-side options both constructors take from the reaction spec.
  void SetReaction(const azure2::experiment::ReactionSpec &);
  /// What a term holds that does not change after the set-up: its settings and its frame at load.  Its points, its
  /// present frame and its shift tables are a working copy's (TermState).
  struct TermFixed {
    std::shared_ptr<const TermSettings> settings;  ///< what the term's points compute (TermSettings.h)
    std::shared_ptr<const azure2::experiment::MeasurementFrame> loadFrame;  ///< KeepLoadFrame
  };
  /// A leaf: its pairs, key, fixed angle, scale, role and target effect.
  struct Leaf {
    int key = 0, entranceKey = 0, exitKey = 0;
    std::optional<double> fixedAngle;  ///< deg lab
    double scaling = 1.0;
    bool yieldTotal = false;
    int targetEffectNum = 0;
    TermFixed fixed;
  };
  TermFixed own_;
  std::vector<Leaf> leaves_;
  int pos_ = -1;  // Pos
  TermFixed &Fixed(int t) { return t ? leaves_[t - 1].fixed : own_; }
  /// The settings of the segment's own term from its members; given to its points (SetTermSettings).
  std::shared_ptr<TermSettings> MakeOwnSettings() const;
  void SetTermSettings(WorkingCopy &data, int t, std::shared_ptr<const TermSettings> settings);
  /// After a set-up change of the segment's pairs or THM attachments.
  void RefreshOwnSettings(WorkingCopy &data);
  /// Leaf t's settings: the segment's with the leaf's pairs and role.
  void SetLeafSettings(WorkingCopy &data, int t);
  const TermFixed &Fixed(int t) const { return t ? leaves_[t - 1].fixed : own_; }
  static ChannelRowSource ShiftRowSource(TermState &term, int i);
  static std::vector<double> ShiftState(const SegmentState &state);
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
  double j_;
  double startNorm_;
  double dataNormNominal_;
  double dataNormError_;
  double startEnergyShift_;
  double energyShiftNominal_;
  double energyShiftError_;
  bool varyEnergyShift_;
  double startEnergyShiftSqrt_;
  double energyShiftSqrtError_;
  double energyShiftSqrtCentre_ = 0.0;  // the prior's centre (MeV^1/2), 0 unless written
  bool varyEnergyShiftSqrt_;
  std::string shiftGroup_;
  int shiftLeaderPos_;
  std::string shiftFnGroup_;
  std::string normGroup_;
  int normLeaderPos_;
  std::string covarianceFile_;
  int covGroup_ = -1;
  bool hasSqrtShiftBlock_;
  // Shift function declared on this line (owner) ...
  std::string shiftFnDeclared_;
  std::vector<double> shiftFnNominal_;
  std::vector<double> shiftFnError_;
  std::vector<double> shiftFnCentre_;  // the priors' centres, 0 unless written
  std::vector<bool> shiftFnVary_;
  bool normPriorReplaced_ = false, shiftPriorReplaced_ = false, shiftSqrtPriorReplaced_ = false;
  bool thicknessPriorReplaced_ = false;
  std::vector<bool> shiftFnPriorReplaced_;
  // ... and the one applied to the points (own or the group leader's).
  bool hasShiftFn_;
  int numShiftFnParams_ = 0;
  std::string shiftFnExpression_;
  const Config *shiftFnConfig_;
  int shiftFnSourcePos_;
  int shiftFnParamStart_;
  bool isAdvanced_;
  int operationType_;
  std::string naturalToken_ = "0";
  int naturalSample_ = 0;
  bool isTotal_ = false;  // IsTotal
  double transmission_ = 0.0;  // IsTransmission
  double yield_ = 0.0;         // IsYield
  double startThickness_ = 0.0;  // StartThickness
  bool thicknessFit_ = false;  // HasThicknessFit
  int thicknessVary_ = 0;
  double thicknessSigma_ = 0.0;
  double thicknessCentre_ = 0.0;  // the prior's centre (atoms/b): as written, else the start (GetNominalThickness)
  std::string thicknessGroup_;
  azure2::experiment::ParameterNames parameterNames_;
  int thicknessLeaderPos_ = -1;  // IsThicknessGroupFollower
  bool isCoefficients_ = false;
  bool coefficientsAbsolute_ = false;
  azure2::experiment::ReactionSpec reaction_;  ///< GetReaction
  std::shared_ptr<const azure2::experiment::Measurement> measurement_;  ///< GetMeasurement (shared by every copy)
  std::string datafile_;
  azure2::experiment::ColumnConversion columns_;  ///< the rows of datafile_ into the engine's units and frame (format 3)
  OperationType segmentOperationType_;
};

/*!
 * A term of a segment by position: the segment and the term (0 its own, k its k-th leaf).  What the hand-over and
 * the grid re-placement queue; it answers what they ask of a term, the term's pairs, key, target effect and points.
 */
struct TermRef {
  const Segment *segment = nullptr;
  int term = 0;
  bool operator==(const TermRef &o) const { return segment == o.segment && term == o.term; }
  bool operator!=(const TermRef &o) const { return !(*this == o); }
  int NumPoints() const { return segment->NumTermPoints(term); }
  /// Point i, 1-based, in the units of `data`, the data set that holds the segment.
  Unit *GetPoint(WorkingCopy &data, int i) const { return segment->GetTermPoint(data, term, i); }
  UnitRange<Unit> GetPoints(WorkingCopy &data) const { return segment->TermPoints(data, term); }
  int GetEntranceKey() const { return segment->TermEntranceKey(term); }
  int GetExitKey() const { return segment->TermExitKey(term); }
  int GetSegmentKey() const { return segment->TermKey(term); }
  int GetTargetEffectNum() const { return segment->TermTargetEffectNum(term); }
};

#endif
