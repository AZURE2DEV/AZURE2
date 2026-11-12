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
    /// An energy shift on the shared mesh (include/ShiftedFold.h): whether the data set's shift is active at the
    /// set-up; the point's c.m. energy the support [endEnergy, startEnergy] was laid out at; the lowest c.m. energy a
    /// node may have and whether the support's low edge is hard there; the reach R_i and the c.m. shift at the start
    /// values Delta_i(nu_0) (MeV).
    bool shifted = false;
    double centreEnergy = 0.0, floorEnergy = 0.0;
    bool lowHard = false;
    double shiftReach = 0.0, startShift = 0.0;
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
  /// The coefficients a_k = B_k / (f A_0) of the point's angular distribution, k = 0 .. max_ang_dist_order, from the
  /// Legendre coefficients its record holds, into the model record (EvaluateModelPoint).  A record without them (a
  /// failed point) gives none.
  void CopyAngularDistsToModel(azure2::Results &r) const;
  /// Calculate the point (its sub-points and its fold) for the current parameters, outside the energy pass
  /// (WorkingCopy::CalculatePointsInParallel), writing the values into `results`, the buffer of the evaluation laid out
  /// by `data`, the data set that holds the point (its target effects give the fold).  One solve at the point's energy,
  /// or one per sub-point (azure2::SolveAtOneEnergy), then FoldSubPoints.  Null `data` for a point outside any data
  /// set (no target effect: a reaction-rate point).
  void Calculate(WorkingCopy *data, CNuc *, const Config &configure, azure2::Results &results);
  /// Before its sub-points are calculated: with external capture, their energy-dependent values at their energies.
  void PrepareSubPoints(WorkingCopy *data, CNuc *, const Config &configure);
  /// The task of the data set's calculation plan that calculates this point (CalculationPlan.h; -1: none, as for a
  /// sub-point or a point outside a data set).  Set when the plan is built (WorkingCopy::BuildPlan).
  int PlanTask() const { return plan_task_; }
  void SetPlanTask(int task) { plan_task_ = task; }
  void SetTargetEffectNum(int);
  /// The point's own site (its energy, angle and tables), without its slots.
  Site &OwnSite() { return *this; }
  const Site &OwnSite() const { return *this; }
  /// The Site functions of the same names for the point's own site and then for each slot of its fold, in order.
  void CalcLegendreP(WorkingCopy *data, int maxL, CNuc *, const TargetEffect *);
  /// A row set of the point's own over its slots (RowSet.h, not cached), where it reads angle-free slots; after its
  /// own rows and whenever its slots move (a grid of its own moved by an energy shift).  WorkingCopy::BuildRowSets
  /// replaces it by a shared one.  `labFrame`: its term's conversion where its angle is given in the lab, else null.
  void CalcNodeRows(const WorkingCopy &data, int maxL, CNuc *, const TargetEffect *,
                    const azure2::experiment::FrameRule *labFrame);
  /// The rows of the point's angular identity over the slots [first, first + count) of `data` (RowSet.h; its own
  /// Legendre rows, or per slot where they depend on the energy, f_C per slot for elastic scattering; with `labFrame`,
  /// the conversion of a lab angle, the c.m. angle and the Jacobian of every slot).  Throws where the lab angle is
  /// beyond the kinematic maximum at a slot.
  azure2::RowSet RowSetOver(const WorkingCopy &data, int first, int count, int maxL, CNuc *, const TargetEffect *,
                            const azure2::experiment::FrameRule *labFrame) const;
  /// The factor of a fold read through a lab-frame row set (RowSet.h): the point's own dOmega_lab/dOmega_cm as its
  /// frame holds it (SetFrame; 1 for a ratio), so that the fold, a lab value, is the c.m. value the point's comparison
  /// divides by the same Jacobian.  1 for every other point.
  double LabFoldFactor() const { return (RowSetOf() && RowSetOf()->labFrame) ? frame_jacobian_ : 1.0; }
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
  /// Whether the point's fold is a ratio of two folds with the same kernel, (K.(w r))/(K.w): a kind with two values per
  /// slot (SiteObservableRow::values, A_y over its cross section) or a normalised Legendre coefficient (a_k over A_0).
  bool FoldsRatio() const;
  /// The weight w and the ratio's element r of the point's slot `slot` (1-based) for FoldsRatio: the slot's values
  /// (SlotValue, SlotAnalyzingPower; for a coefficient A_0 and a_k).
  void RatioPartsAt(const WorkingCopy &data, int slot, const azure2::Results &results, double &weight,
                    double &element) const;
  /// Where the cotangent of the weight w goes (PointTBar's xsComponent, ayCrossSection): the cross section of an A_y
  /// slot (0, true), A_0 of a coefficient (-1, false); the element's goes to the observable (0, false).
  void RatioWeightRoute(int &component, bool &ayCrossSection) const;
  /// The fold of this point from its slots' values already in `results` (FoldPoint, include/PointFold.h; a ratio,
  /// FoldsRatio, is the ratio of two folds).
  void FoldSubPoints(const WorkingCopy &data, const Config &, azure2::Results &results);
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
  std::shared_ptr<const TermSettings> own_settings_;
  // Beam-profile effect: per-point energy window (c.m.; NaN when the data file has no such columns) and the
  // inverse reaction kinematics for the detailed-balance weight.
  double bin_low_cm_ = std::numeric_limits<double>::quiet_NaN();
  double bin_high_cm_ = std::numeric_limits<double>::quiet_NaN();
  double photo_qgamma_ = 0.0;
  double photo_mcn_ = 0.0;
  double frame_jacobian_ = 1.0;  ///< dOmega_lab/dOmega_cm of the point's frame row (SetFrame; LabFoldFactor)
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