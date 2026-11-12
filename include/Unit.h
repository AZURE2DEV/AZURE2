#ifndef UNIT_H
#define UNIT_H

#include "Site.h"

class Unit;
typedef Unit *UnitIterator;

/// An AZURE data point
/*!
 * A data point: a calculation site (Site) with its measurement row, its target-effect fold over its slots and its
 * result and model records.  Its entrance and exit pair, energy, angle, conversions and the flags that set the type
 * of data (angle integrated or differential) are those of its site.
 */
class Unit : public Site {
 public:
  /// Compute everything that depends on energy but not on the fit parameters. Run once before fitting.
  /// `data`: the data set that holds the point (its target effect may carry Q coefficients), null outside one.
  void Initialize(WorkingCopy *data, CNuc *, const Config &);
  /// Take the energies, c.m. angle and c.m. bin edges of the point's row from its segment's frame
  /// (Segment::ConvertToCM, experiment::ConvertRow).
  void SetFrame(const azure2::experiment::FrameRow &row);
  /// The point of row `row` (0-based) of a measurement, taking its conventions from the segment.  The point reads
  /// what was measured there (energy, value, uncertainty, bins, file record, Legendre order) from the rows, which the
  /// segment's measurement keeps alive and never changes.
  Unit(const azure2::experiment::MeasurementRows &rows, int row, const Segment *segment, int term);
  /// A point at a given c.m. energy and angle that is no row of a measurement (a sub-point of a target effect), with
  /// the conventions of term `term` of the segment (Segment::TermEntranceKey ...).
  Unit(double angle, double energy, const Segment *segment, int term);
  /// Build with every convention given explicitly, rather than inherited from a segment.
  Unit(double, double, int, int, bool, bool, bool, double, int, int);
  /// Does the point carry target effects?
  bool IsTargetEffect() const;
  /// How the adaptive sub-point grid of a point was laid out (WorkingCopy::SetUpTargetEffects), so WorkingCopy::
  /// ReplaceSubPointGrids can lay it out again for other level parameters; held by the data set (DataSet::PointSetUp).
  /// Energies are c.m., unshifted.
  struct SubGridSpec {
    bool valid = false;
    double startEnergy = 0.0, endEnergy = 0.0;
    int maxPoints = 0, entranceKey = 0;
    double baseEnergyStep = 0.0, resonanceWidthMultiplier = 0.0, pointsPerWidth = 0.0;
    bool targetIntegration = false;
    double stoppingConversion = 1.0;  ///< c.m./lab factor of the stopping power
    /// On the shared mesh of its nucleus (docs/upgrade/energy_grid.md section 2.2): the mesh (DataSet::Meshes) and
    /// the support's index in it, -1 for a point with a grid of its own; E_x - E_cm of its pair (MeV).
    int mesh = -1, meshSupport = -1;
    double frameShift = 0.0;
  };
  /// What the set-up worked out for this point that its evaluations read again (its sub-point grid's layout, its
  /// user-defined resolution kernel): record SetUpIndex of the data set's DataSet::PointSetUps, -1 for none.
  int SetUpIndex() const { return set_up_; }
  void SetSetUpIndex(int i) { set_up_ = i; }
  /// The point's own lab -> c.m. factor and entrance-pair charges, with which every target-effect kernel
  /// converts the line's lab energies (beam sigma, beam profile, Bohr straggling).  Set by
  /// WorkingCopy::SetUpTargetEffects for each point and component point; a line shared by several pairs
  /// used to be converted once with the first pair's factor (2026-10-27).
  void SetTargetKinematics(double cmFactor, int z1, int z2) { cm_factor_ = cmFactor; bohr_z1_ = z1; bohr_z2_ = z2; }
  double GetCMFactor() const { return cm_factor_; }
  int GetBeamCharge() const { return bohr_z1_; }
  int GetTargetCharge() const { return bohr_z2_; }
  /// Number of sub-points, the samples a target-effect integral is evaluated on.
  int NumSubPoints() const { return sub_count_; }
  /// 1-based position of the TargetEffect in the parent WorkingCopy.
  int GetTargetEffectNum() const;
  /// Weight in [0,1] with which the target effect applies to this point;
  /// 1 everywhere unless the effect declares energy ranges.
  double GetTargetBlendWeight() const { return target_blend_weight_; };
  void SetTargetBlendWeight(double w) { target_blend_weight_ = w; };
  /// Thin-target integration restricted to energy windows (2026-11-05): outside the windows (blend weight 0) the
  /// point is evaluated at the single depth-averaged energy E_bar = E - D of its own target instead of being a
  /// plain point at its data (beam) energy, which was off by the centroid shift ~Delta/2.  Such a point carries
  /// one sub-point at E_bar (fixed density) or three spanning the range of E_bar a fitted density can reach,
  /// and SetTargetFallbackOnly(true) marks it.  D (c.m. MeV, TargetEffect::MeanEnergyLoss) is set on every thin-target
  /// point: edge blending and the automatic decision fall back to sigma(E_bar) as well.
  void SetTargetFallbackOnly(bool f) { target_fallback_only_ = f; }
  bool IsTargetFallbackOnly() const { return target_fallback_only_; }
  /// Inside a TargetFallbackScope.
  bool IsTargetFallbackNow() const { return target_fallback_now_; }
  double GetTargetFallbackDepth() const { return target_fallback_depth_; }
  void SetTargetFallbackDepth(double d) { target_fallback_depth_ = d; }
  /// While alive, the target integration of the point returns sigma(E_bar) from its sub-points (an automatic
  /// tolerance decided the integration is not needed at this point).
  class TargetFallbackScope {
   public:
    explicit TargetFallbackScope(Unit *p) : p_(p), prev_(p->target_fallback_now_) { p_->target_fallback_now_ = true; }
    ~TargetFallbackScope() { p_->target_fallback_now_ = prev_; }
    TargetFallbackScope(const TargetFallbackScope &) = delete;
    TargetFallbackScope &operator=(const TargetFallbackScope &) = delete;
   private:
    Unit *p_;
    bool prev_;
  };
  /// Where the point's model record lies in Results::model (WorkingCopy::LayoutResults; kNoRecord for a component point, a
  /// sub-point or a point outside a data set), and how many angular-distribution coefficients it holds.  The model
  /// value is what the data set predicts at the point (ModelPass.h); the leaf values above are what the point itself
  /// calculates.
  void SetModelRecord(std::size_t slot, int nAng) { model_slot_ = slot, num_model_ang_ = (std::uint16_t)nAng; }
  double GetModelValue(const azure2::Results &r) const { return r.model[ModelRecord(r)]; }
  void SetModelValue(azure2::Results &r, double v) const { r.model[ModelRecord(r)] = v; }
  /// The point's covariance seed (M^-1 r)_i (Results::ObjectiveValues::covarianceSeed), at its model record.
  double GetCovarianceSeed(const azure2::Results &r) const { return r.objective.covarianceSeed.at(ModelRecord(r)); }
  void SetCovarianceSeed(azure2::Results &r, double z) const { r.objective.covarianceSeed.at(ModelRecord(r)) = z; }
  /// The coefficients a_k of an angular-distribution point as the model record holds them (count, then values).
  int GetNumModelAngularDists(const azure2::Results &r) const {
    if (!num_model_ang_) return 0;
    const double n = r.model[ModelRecord(r) + 1];  // NaN (poisoned, not written): the capacity, as the leaf's count
    return std::isnan(n) ? (int)num_model_ang_ : (int)n;
  }
  double GetModelAngularDist(const azure2::Results &r, int k) const { return r.model[ModelRecord(r) + 2 + k]; }
  /// Copy the leaf record's angular-distribution count and coefficients into the model record (EvaluateModelPoint).
  void CopyAngularDistsToModel(azure2::Results &r) const {
    if (!num_model_ang_) return;
    const double *from = Rec(r) + azure2::Results::AngularOffset();
    std::copy(from, from + 1 + std::min<int>(num_model_ang_, num_ang_slots_), r.model.begin() + ModelRecord(r) + 1);
  }
  /// Calculate the point (its sub-points and its fold) for the current parameters, outside the energy pass
  /// (WorkingCopy::CalculatePointsInParallel), writing the values into `results`, the buffer of the evaluation laid out
  /// by `data`, the data set that holds the point (its target effects give the fold).  One solve at the point's energy,
  /// or one per sub-point (azure2::SolveAtOneEnergy), then FoldSubPoints.  Null `data` for a point outside any data
  /// set (no target effect: a reaction-rate point).
  void Calculate(WorkingCopy *data, CNuc *, const Config &configure, azure2::Results &results);
  /// Before its sub-points are calculated: with external capture, their energy-dependent values at their energies.
  void PrepareSubPoints(WorkingCopy *data, CNuc *, const Config &configure);
  /// The fold of this point from its sub-points' values already in `results`: the bare probe, the automatic tolerance
  /// and edge blending, the thin-target fallback, the integrations.
  void FoldSubPoints(WorkingCopy *data, CNuc *, const Config &configure, azure2::Results &results);
  /// The fold of a thin-target integration, with its sigma(E_bar) fallback (see the .cpp).
  void FoldThinTarget(WorkingCopy &data, azure2::Results &results, const Config &configure, const TargetEffect *effect);
  /// The task of the data set's calculation plan that calculates this point (CalculationPlan.h; -1: none, as for a
  /// sub-point, a probe copy or a point outside a data set).  Set when the plan is built (WorkingCopy::BuildPlan).
  int PlanTask() const { return plan_task_; }
  void SetPlanTask(int task) { plan_task_ = task; }
  /// Take a copy of the point out of the plan: it calculates itself alone (a probe).
  void DetachFromPlan() { plan_task_ = -1; }
  void SetTargetEffectNum(int);
  /// The point's own site (its energy, angle and tables), without its slots.
  Site &OwnSite() { return *this; }
  const Site &OwnSite() const { return *this; }
  /// The Site functions of the same names for the point's own site and then for each slot of its fold, in order.
  void CalcLegendreP(WorkingCopy *data, int maxL, CNuc *, const TargetEffect *);
  void CalcEDependentValues(WorkingCopy *data, CNuc *, const Config &, const ChannelRowSource *rows = nullptr);
  void ChannelKeys(WorkingCopy *data, CNuc *, std::vector<ChannelTable::Key> &);
  void CalcCoulombAmplitude(WorkingCopy *data, CNuc *);
  void CalculateECAmplitudes(WorkingCopy *data, CNuc *, const Config &, const ChannelTable *table = nullptr);
  /// Site::RecalcEDependentValues for the point and, when the point's own values were not current, each slot
  /// emptied and recomputed at its present energy (rows and Coulomb amplitude).
  void RecalcEDependentValues(WorkingCopy *data, CNuc *, const Config &, const ChannelRowSource *rows = nullptr);
  /// Make `sub` a sub-point of this point (WorkingCopy::AddSubPoint places it): it shares the term's settings and the
  /// point's screening.
  void MakeSubPoint(Site &sub) const;
  /// The point's slots, a range of the slots of its working copy (WorkingCopy::SlotAt); set by the working copy.
  int SubPointFirst() const { return sub_first_; }
  void SetSubPointRange(int first, int count) { sub_first_ = first, sub_count_ = count; }
  /// The point's fold of its observable slots (FoldPoint, include/PointFold.h, on each slot the data set compares).
  void IntegrateTargetEffectForObservable(const WorkingCopy &data, const Config &, azure2::Results &results);
  /// Integrate the E1 and E2 components of the cross section over the same
  /// sub-points as the total, for a segment compared against one component.
  void IntegrateTargetEffectComponents(const WorkingCopy &data, const Config &, azure2::Results &results);
  /// Per-point energy window (columns 5-6 of the data file) of a beam-profile
  /// effect: the detector-reconstructed energy slice the point was built from.
  bool HasBinWindow() const { return bin_low_cm_ == bin_low_cm_ && bin_high_cm_ == bin_high_cm_; };
  /// Set the per-point energy window directly, in the c.m. frame.  A row of a measurement takes it from columns 5-6
  /// of the data file through its frame; this is for a point built outside that path, as the beam-profile
  /// reference test does.
  void SetBinWindow(double low, double high) {
    bin_low_cm_ = low;
    bin_high_cm_ = high;
  };
  double GetBinLowCM() const { return bin_low_cm_; };
  double GetBinHighCM() const { return bin_high_cm_; };
  /// Kinematics of the inverse photodissociation reaction for the
  /// detailed-balance weight: Q_gamma = S_entrance - Ex_final (MeV) and the
  /// compound-nucleus mass (MeV).
  void SetPhotoKinematics(double qGamma, double compoundMassMeV) { photo_qgamma_ = qGamma; photo_mcn_ = compoundMassMeV; };
  double GetPhotoQ() const { return photo_qgamma_; }
  double GetPhotoCompoundMass() const { return photo_mcn_; }
  /// Slot i (1-based) of the point in `data`, the working copy that holds it (WorkingCopy::SlotAt); null beyond the
  /// last.
  Site *GetSubPoint(WorkingCopy &data, int i);
  const Site *GetSubPoint(const WorkingCopy &data, int i) const;
  UnitRange<Site> GetSubPoints(WorkingCopy &data);
 private:
  double cm_factor_ = 1.0;   ///< c.m./lab energy factor of the point's entrance pair (SetTargetKinematics)
  int bohr_z1_ = 0;          ///< beam charge of the entrance pair, for Bohr straggling
  int bohr_z2_ = 0;          ///< target charge of the entrance pair
  int set_up_ = -1;  ///< SetUpIndex
  int targetEffectNum_;
  double target_blend_weight_ = 1.0;
  bool target_fallback_only_ = false;   ///< see SetTargetFallbackOnly
  bool target_fallback_now_ = false;    ///< see TargetFallbackScope
  double target_fallback_depth_ = 0.0;  ///< D = E - E_bar, c.m. MeV (see SetTargetFallbackOnly)
  std::shared_ptr<const TermSettings> own_settings_;
  // Beam-profile effect: per-point energy window (c.m.; NaN when the data file has no such columns) and the
  // inverse reaction kinematics for the detailed-balance weight.
  double bin_low_cm_ = std::numeric_limits<double>::quiet_NaN();
  double bin_high_cm_ = std::numeric_limits<double>::quiet_NaN();
  double photo_qgamma_ = 0.0;
  double photo_mcn_ = 0.0;
  int plan_task_ = -1;  // see PlanTask
  int sub_first_ = 0, sub_count_ = 0;  ///< the slots: [first, first + count) of WorkingCopy::SlotAt
  std::uint16_t num_model_ang_ = 0;  // see SetModelRecord
  // Where the model record lies (see SetModelRecord).
  std::size_t model_slot_ = kNoRecord;
  std::size_t ModelRecord(const azure2::Results &r) const {
    if (model_slot_ == kNoRecord || model_slot_ + azure2::Results::ModelRecordSize(num_model_ang_) > r.model.size())
      NoModelRecord(r);
    return model_slot_;
  }
  [[noreturn]] void NoModelRecord(const azure2::Results &r) const;
};

#endif