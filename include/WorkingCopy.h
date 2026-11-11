#ifndef WORKINGCOPY_H
#define WORKINGCOPY_H

#include "CalculationPlan.h"
#include "DataSet.h"
#include "SegmentState.h"
#include "ChannelTable.h"
#include "Segment.h"
#include "TargetEffect.h"
#include "NodeTableCache.h"
#include "NodeTable.h"
#include "ThmExperiment.h"
#include "ThmGroup.h"
#include "EvaluationCount.h"
#include "Results.h"
#include "ParameterKind.h"
#include <memory>
#include <atomic>
#include <ctime>
#include <deque>
#include <map>
#include <optional>
#include "DataCovariance.h"
#include <memory>
#include <fstream>
#include <ios>

// THM types WorkingCopy holds by pointer or takes by reference (ThmNorm.cpp).
struct ThmLineshape;
struct ThmSpectatorWindow;
struct ThmAngleWindow;
struct ThmVertexReport;
struct ThmLineshapeReport;
struct ThmDistortionReport;
class ThmDistortion;
class ThmDwVertex;

class CNuc;
class AZUREParams;
class PriorCovariance;
struct BandData;
namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT

/*!
 * A working copy of the data of a project (docs/upgrade/data_structure.md section 1.4): the data set it evaluates
 * (DataSet, immutable and shared by every copy through Data()), and what an evaluation changes -- the state of every
 * segment and term (SegmentState: norms, shifts, frames, shift tables), the present values of the target effects'
 * fit parameters, the calculation units (one flat vector), the result layout.  A copy (Clone) is a plain copy that
 * shares the data set.  The set-up builds the data set through this object (Draft) and freezes it at the end of
 * Initialize.
 */

class WorkingCopy {
 public:
  using LeafRef = DataSet::LeafRef;
  WorkingCopy();
  WorkingCopy(const WorkingCopy &);
  WorkingCopy &operator=(const WorkingCopy &) = delete;
  /// The data set this copy evaluates: never written after the set-up (DataSet.h).
  const DataSet &Data() const { return *set_; }
  /// The shared pointer to it (two copies of one working copy share one data set).
  const std::shared_ptr<const DataSet> &DataSetPointer() const { return set_; }
  /// The data set while the set-up builds it; throws once it is frozen.
  DataSet &Draft();
  /// The end of the set-up: from now on the data set is only read (Initialize calls it last).
  void Freeze() { draft_.reset(); }
  bool IsFrozen() const { return !draft_; }
  /*!
   * A change of a setting after the set-up (a prior covariance the caller gives, a parameter layout laid out again):
   * this copy gets a new data set, a copy of the old one with `change` applied, and then reads only that.  Copies
   * made before keep the old one.  During the set-up the change is made to the draft.
   */
  template <class F>
  void ReviseDataSet(F change) {
    if (draft_) {
      change(*draft_);
      return;
    }
    auto next = std::make_shared<DataSet>(*set_);
    change(*next);
    set_ = std::move(next);
  }
  /// The state of segment `segment` (one of Data()'s) in this copy.
  SegmentState &State(const Segment &segment) { return segmentStates_[segment.Pos()]; }
  const SegmentState &State(const Segment &segment) const { return segmentStates_[segment.Pos()]; }
  /// The present values of target effect n's (1-based) fit parameters in this copy.
  TargetEffectValues &EffectValues(int n) { return effectValues_[n - 1]; }
  const TargetEffectValues &EffectValues(int n) const { return effectValues_[n - 1]; }
  /// The present value of THM coherent-background parameter k (0-based).
  double ThmCoherentValue(int k) const { return thmCoherentValues_[k]; }
  /// Build the segments from the data-set specs (each filled from its data file) and their effects.  Returns -1 on
  /// failure (with the message).
  int Fill(const Config &, CNuc *, const azure2::experiment::DataSpecs &);
  /// Build points from the extrapolation grids instead of data files, for a run with no data.  Returns -1 on failure.
  int MakePoints(const Config &, CNuc *, const azure2::experiment::ExtrapolationSpecs &);
  /// TEMPORARY (step 3d): the same from the old sections of the project (src/ProjectSetup.cpp).  A format-3 project
  /// that the set-up refuses returns kFillRefused: its data files cannot be read (the files and the refusal line
  /// printed) or a set-up check refused it (its ERROR line printed, the run's last); the caller prints nothing more.
  int Fill(const Config &, CNuc *);
  static constexpr int kFillRefused = -2;
  int MakePoints(const Config &, CNuc *);
  /*!
   * The prior centres of the data sets' normalizations and energy shifts (DatasetSpec::normCentre, shiftCentre):
   * the prior of a varied normalization or shift is centred there instead of on the start value (normalization)
   * or 0 (shift).  Applied after the groups are resolved.  A centre for a data set without a segment (not active,
   * or skipped at load) is reported and ignored; specs.priorCentreError (a row the old-line reader refused) is
   * reported and returns -1.
   */
  int ApplyPriorCentres(const Config &, const azure2::experiment::DataSpecs &specs);
  /// Push each screening effect's current U_e onto the points and sub-points of its segments.
  void ApplyScreening();
  /// Set the screening effects' U_e from a parameter vector, then ApplyScreening (called by
  /// FillEnergyShiftsFromParams, so every evaluation sees the vector's values).
  void FillScreeningFromParams(const vector_r &);
  /// Sum of the screening priors ((U_e - nominal)/sigma)^2.
  double ScreeningPenalty() const;
  /// Set the fitted densities from a parameter vector and recompute the target thickness of every point of
  /// their lines (called by FillEnergyShiftsFromParams, so every evaluation sees the vector's values).
  void FillDensitiesFromParams(const vector_r &, const Config &);
  /// Sum of the density priors ((N - N0)/sigma)^2.
  double DensityPenalty() const;
  /// The apparatus parameters from the parameter vector (ApparatusParameters.cpp): each row writes the values of the
  /// effects of its parameters, then what follows them (the stopping scale: the sub-points' stopping, the thickness).
  void FillApparatusFromParams(const vector_r &, const Config &);
  /// The sum of the priors ((x - x0)/sigma)^2 of the apparatus parameters of row `row` (-1: of every row).
  double ApparatusPenalty(int row = -1) const;
  /// What follows a change of the stopping scale of the effects marked in `changed` (0-based): the stopping of their
  /// sub-points and the thickness of their points (StoppingScale.cpp, the row's follow-up).
  void FollowStoppingScales(const std::vector<char> &changed, const Config &);
  /// Set the fitted temperatures from a parameter vector (clamped to [T0/4, layout temperature], with a warning
  /// printed once), called by FillEnergyShiftsFromParams so every evaluation sees the vector's values.
  void FillDopplerFromParams(const vector_r &, const Config *configure);
  /// Sum of the Doppler priors ((T - T0)/sigma)^2.
  double DopplerPenalty() const;
  /// Set the fitted thicknesses from a parameter vector (a value below 0 is taken as 0); thickness-group followers
  /// take their leader's.  Called by FillEnergyShiftsFromParams.
  void FillThicknessFromParams(const vector_r &);
  /// Sum of the thickness priors ((n - n0)/sigma)^2.
  double ThicknessPenalty() const;
  /// Link the thicknessgroup members to the group's first segment in file order; -1 (with an ERROR) when a member has
  /// no transmission / yield block.  Also lays out the thickness parameters (NumThicknessParams).
  int ResolveThicknessGroups(const Config &);
  /// Build the TargetEffect objects from their specs and give each linked segment its copy.
  int SetUpTargetEffects(const Config &, CNuc *,
                         const azure2::experiment::SpecList<azure2::experiment::TargetEffectSpec> &);
  /// Is this a fit? AZURECalc clones the compound nucleus and data per thread only when it is.
  bool IsFit() const;
  /// Is this a Minos error-analysis call? Suppresses the transformation and file output.
  bool IsErrorAnalysis() const;
  void SetFit(bool);
  void SetErrorAnalysis(bool);
  /// Count one more objective evaluation. Returns the count before this one (0 for the first), taken in one atomic
  /// operation, so concurrent callers (Minuit2's OpenMP numerical gradient) each get a distinct number.
  int Iterate();
  /// Initialize every point: Unit::Initialize over the whole data set.
  int Initialize(CNuc *, const Config &);
  /// Append a segment.
  void AddSegment(Segment);
  /// Print the data as read, or as generated for a run without data.
  void PrintData(const Config &);
  /// Unit::CalcLegendreP for every point.
  void CalcLegendreP(int, CNuc *);
  void PrintLegendreP(const Config &);
  /// Unit::CalcEDependentValues for every point -- penetrabilities, shifts and Coulomb phases at each energy.
  int CalcEDependentValues(CNuc *, const Config &);
  void PrintEDependentValues(const Config &, CNuc *);
  /// Unit::CalcCoulombAmplitude for every point.
  void CalcCoulombAmplitude(CNuc *);
  void PrintCoulombAmplitude(const Config &, CNuc *);
  /// Write AZUREOut_*, chiSquared.out and the rest of the run's output files.  The values written are those of
  /// `results`, the buffer of the last evaluation of this data set: its model values and what its objective stored
  /// (Results::objective).  The THM profiles are taken again from the model values first (a snapshot written during
  /// a fit has not been through the objective).
  void WriteOutputFiles(const Config &, azure2::Results &results, bool = false, const BandData * = nullptr);
  /// The normalization the data of segment i (1-based) are compared with in `results`: the profiled n* of a THM
  /// segment (stored by the objective or the writer), else the segment's own (set by the hand-over).
  double NormUsed(int i, const azure2::Results &results) const {
    const Segment &segment = *Data().GetSegment(i);
    return segment.IsProfiledNorm() ? results.objective.norm[i] : State(segment).norm;
  }
  /// The prior parameter covariance of the run (<config> priorCovariance=, PriorCovariance.h; DataSet::
  /// GetPriorCovariance), null without one.  Set by the CLI (AZUREMain) and pyazr (AZUREAPI) once the parameters are
  /// known: a revision of the data set (ReviseDataSet).
  void SetPriorCovariance(std::shared_ptr<const PriorCovariance> prior);
  /// External-capture amplitudes for every point that has an EC component.
  int CalculateECAmplitudes(CNuc *, const Config &);
  /// How many external-capture amplitudes this model expects in an intEC file.
  /// Mirrors CalculateECAmplitudes exactly, so a mismatch means the cached file belongs to a different model.
  long long CountECAmplitudes(CNuc *, const Config &);
  /// Hash of everything the external-capture integrals depend on: the energies they are evaluated at,
  /// the pairs and final states that enter them, and the Coulomb-function options.  Stored beside an
  /// intEC file (see ECSignaturePath) so that a file built for another grid is not reused.
  std::string ECSignature(CNuc *, const Config &);
  /// The sidecar that records the signature of an intEC file: the file name with ".sig" appended.
  static std::string ECSignaturePath(const std::string &integralsFile);
  /// The leaves after the segments: their frames with their own pairs, energy-dependent values, Legendre
  /// polynomials, Coulomb amplitudes and external-capture amplitudes (appended to the intEC file).
  int InitializeLeaves(CNuc *, const Config &);
  /// Validates a component's pair keys (both present, one compound nucleus) and flags its entrance pair; false with the error printed otherwise.
  bool CheckComponentPairs(CNuc *, const Segment &baseSegment, int entranceKey, int exitKey, const Config &);
  /// A leaf of `base` (Segment::AddLeaf) with these pairs, a fixed lab angle (deg; a RATIO denominator) or none and a
  /// scale, its key from the data set's counter; returns its term.
  int AddLeaf(Segment &base, int entranceKey, int exitKey, std::optional<double> fixedAngle, double scaling);
  /// A transmission or yield segment after its components exist: checks, plain components, the yield's total term.
  bool FinishTransformedSegment(CNuc *, Segment *, int key, bool isTest, const Config &);
  /// Every term of a transmission or yield segment on the base point's sub-point grid (after the target effects).
  bool CheckTransformedSegments(const Config &);
  /// The expression tree of every segment (src/Predictions.cpp), once its components, pairs and sample are resolved;
  /// false, with a message, if one cannot be built.
  bool BuildPredictions(CNuc *, const Config &);
  /// The calculation plan of the data set, one task per point, built and checked (CalculationPlan.h).
  void BuildPlan();
  /// The node table of the stage (NodeTable.h), built at the end of Initialize and shared by the copies; null before.
  const azure2::NodeTable *Nodes() const { return nodes_.table.get(); }
  /// AZURE_NODE_CHECK=1 (the point pass, before its tasks): the node table and slots against the units of the plan
  /// tasks `work`, bit for bit; prints the counts to stderr and returns the number that differ (src/NodeTable.cpp).
  long CheckNodeTable(CNuc *cnuc, const std::vector<int> &work);
  /// The slots of every plan task in this copy (NodeTable.h): the stage table's until the hand-over moves a unit off it.
  const azure2::NodeSlots *Slots() const { return nodes_.slots.get(); }
  /// The copy's overflow table (NodeTable.h): the nodes of its units that the stage table does not hold.
  const azure2::NodeOverflow &Overflow() const { return nodes_.overflow; }
  /// The end of the hand-over (src/NodeTable.cpp): the slots of the plan tasks `changedTasks` (their units moved: a
  /// re-placed grid, an energy shift) and of every task already on the overflow table looked up again, the overflow
  /// table placed again from nothing for them.  With nothing changed and the rows not released, nothing is done.
  void PlaceOverflowNodes(CNuc *cnuc, const std::vector<int> &changedTasks, bool rowsReleased);
  /// Let go of the overflow table's rows before the units are re-energised (they are holders of the units' rows).
  void ReleaseOverflowRows();
  /// The calculation plan: which task calculates which point (CalculationPlan.h).  Positions only: a shift or a
  /// re-placed sub-point grid leaves it as it is (the energy tasks it gives are rebuilt from the units).
  const azure2::CalculationPlan &Plan() const { return plan_; }
  /// A point of the plan, by its unit.
  Unit *PlanPoint(azure2::PointRef r) { return &UnitAt(r); }
  /// The units of every point of every term in data order: the segments' own points, then the leaves' (LeafRefs).
  std::vector<azure2::PointRef> PlanUnits() const;

  /// The calculation units of the working copy (data_structure.md section 1.4), src/WorkingCopyUnits.cpp: every point
  /// of every term and every sub-point in one vector; a term's points and a point's sub-points are ranges of it.
  Unit &UnitAt(int i) { return units_[i]; }
  const Unit &UnitAt(int i) const { return units_[i]; }
  int NumUnits() const { return (int)units_.size(); }
  /// The index of a unit of this working copy (UnitAt(UnitIndex(u)) is u).
  int UnitIndex(const Unit *u) const { return (int)(u - units_.data()); }
  UnitRange<Unit> Units(int first, int count) { return UnitRange<Unit>(units_.data() + first, count); }
  UnitRange<const Unit> Units(int first, int count) const {
    return UnitRange<const Unit>(units_.data() + first, count);
  }
  /// Append a unit (a point of a term, at set-up); returns its index.
  int AppendUnit(Unit unit) {
    units_.push_back(std::move(unit));
    return (int)units_.size() - 1;
  }
  /// Give `point` (a unit of this data set) the sub-point `sub` after those already given; held until
  /// CommitSubPoints, so that no unit moves while the set-up holds pointers to points.
  void AddSubPoint(Unit &point, Unit sub);
  /// Place the sub-points given since the last commit (each point's replace any it had) and lay the units out again
  /// (CompactUnits).  Single-threaded.
  void CommitSubPoints();
  /// Per work item of a grid re-placement: (a point's unit, its new sub-points).  PlaceSubPointGrids puts a grid with
  /// as many nodes as the old one in its place, any other in a new layout of the units.  Single-threaded.
  using SubPointGrids = std::vector<std::vector<std::pair<int, std::vector<Unit>>>>;
  void PlaceSubPointGrids(SubPointGrids &grids);
  /// Lay the units out again (see Units), dropping every unit no term and no point refers to.
  void CompactUnits();
  /// Append a target-effect definition.
  void AddTargetEffect(TargetEffect);
  /// Seed the Minuit parameter array with the normalizations and energy shifts.
  /// Link the members of each "shiftgroup" to the group's first segment (called by Fill).
  /// Also installs and checks the shift functions ("shiftfn"); returns -1 on an invalid one.
  int ResolveShiftGroups(const Config &);
  /// Install and check the shift functions (called by ResolveShiftGroups).
  int ResolveShiftFunctions(const Config &);
  /// Set each shift-group follower's parameter to its leader's value (for the parameter files).
  void SyncShiftGroupParams(ROOT::Minuit2::MnUserParameters &) const;
  /// Link the members of each "normgroup" to the group's first segment (called by Fill).
  void ResolveNormGroups(const Config &);
  /// Set each norm-group follower's (fixed) norm parameter to its leader's value.
  void SyncNormGroupParams(ROOT::Minuit2::MnUserParameters &) const;
  /// The data covariance file \p path, read on first use (segment Fill); null (with an
  /// ERROR printed and SetCovarianceError) if it cannot be read.
  std::shared_ptr<DataCovariance::File> CovarianceFileFor(const std::string &path, const Config &);
  void SetCovarianceError() { covarianceError_ = true; }
  /// A data point beyond the kinematic maximum lab angle (Segment::Fill): the run stops instead of dropping the segment.
  void SetKinematicsError() { kinematicsError_ = true; }
  /// Gather the points of the segments carrying "covariance <file>" into one group per
  /// file and factor each group's matrix (called by Fill).  Returns -1 on an error.
  int ResolveCovarianceGroups(const Config &);
  /// From the points' model values in `results`: the chi2 of every covariance group, r^T M^-1 r,
  /// minus the sum of r_i^2 the per-point loops counted for those points (so adding the
  /// return value to their total gives the covariance chi2).  Sets each covariance
  /// segment's chi2 in `results` (and in \p segChis, indexed by segment number, when given) to the sum over its
  /// points of r_i (M^-1 r)_i, and the points' seeds (M^-1 r)_i.
  double CovarianceChi2Correction(azure2::Results &results, std::vector<double> *segChis = nullptr);
  /// Whiten the data rows of a residual vector (and of a row-major Jacobian with \p nCols
  /// columns, when given): r <- W r within every covariance group (W^T W = M^-1).  Rows are the
  /// points in segment order, as ComputeResidualJacobian lays them out.
  void WhitenResidualRows(std::vector<double> &residuals, std::vector<double> *jacobian, int nCols);
  /// covariance.out: each group's chi2 and the diagonal-only sum, at the current fit values.
  void WriteCovarianceSummary(const Config &, const azure2::Results &results);
  /// Forget the frozen model of the relvector terms; the next evaluation freezes it again.
  void RefreezeCovariance();
  /// Calculate the covariance segments' points at the current parameters and set their
  /// seeds, for a gradient pass that needs (M^-1 r) before it reaches the points.
  void PrepareCovarianceSeeds(CNuc *, const Config &, azure2::Results &results);
  void FillMnParams(AZUREParams &);
  /// Write the normalizations back from a Minuit parameter vector.
  void FillNormsFromParams(const vector_r &);
  /// Write the energy shifts back from a Minuit parameter vector.  The last step of the hand-over of every
  /// evaluation, so it ends with BeginEvaluation(results): `results` is the buffer of the evaluation that follows,
  /// owned by its caller.
  void FillEnergyShiftsFromParams(const vector_r &, azure2::Results &results, CNuc *theCNuc = nullptr,
                                  const Config *configure = nullptr);
  /// The start of an evaluation into the caller's buffer `results`, after the hand-over: the records are laid out
  /// again if a point or sub-point is not where the layout puts it (a new point, a re-placed sub-point grid, a copied
  /// segment), each moving its values in `results` to its new record; `results` is sized for the layout (a new buffer
  /// starts at 0); the model values and the objective's per-segment values are sized too; then, under
  /// AZURE_POISON_RESULTS=1, all of them are quiet NaN, so a value read before this evaluation writes it shows as
  /// NaN.  O(points), no allocation when nothing changed.
  void BeginEvaluation(azure2::Results &results);
  /// Size a buffer laid out as the data set's records are (a new buffer: every value 0) for its leaf records, model
  /// values and the objective's per-segment and per-experiment values.
  void SizeResults(azure2::Results &results) const {
    results.slot.resize(resultsSize_, 0.0);
    results.model.resize(modelSize_, 0.0);
    results.objective.norm.resize(Data().NumSegments() + 1, 0.0);
    results.objective.segmentChi2.resize(Data().NumSegments() + 1, 0.0);
    results.objective.thmProfiles.resize(Data().NumThmGroups());
    results.objective.covarianceSeed.resize(Data().HasCovariance() ? modelSize_ : 0, 0.0);
  }
  /// Lay the adaptive sub-point grids out again for the levels theCNuc holds now
  /// (its fit values, as observed parameters) and the segments' energy shifts, and
  /// rebuild the sub-points of every point whose grid changed.  Sub-points of a
  /// segment in \p shiftActive are left for the energy-shift step to energise and
  /// their segments are appended to \p reshift; the others are energised here.
  /// Returns the number of points whose grid changed, and appends their plan tasks to \p replacedTasks.  Called by
  /// FillEnergyShiftsFromParams.
  int ReplaceSubPointGrids(CNuc *theCNuc, const Config &configure, const std::vector<TermRef> &shiftActive,
                           std::vector<TermRef> &reshift, std::vector<int> &replacedTasks);
  /// Drop the last segment.
  void DeleteLastSegment();
  WorkingCopy *Clone() const;

  using ThmGroup = ::ThmGroup;                  // ThmGroup.h
  using ThmCoherentParam = ::ThmCoherentParam;  // ThmGroup.h
  /// Reads the coherent-background block of a full parameter vector (if it
  /// has one) into the parameters and, if given, the compound the HOES model
  /// reads them from.  Called by FillEnergyShiftsFromParams.
  void FillThmCoherentFromParams(const vector_r &p, CNuc *theCNuc);
  /*!
   * Profiles group g from the model values of its segments' points: stores in `results` (Results::objective) the
   * profile, every segment's norm (the shared n*) and its chi^2 (that of its own points), and returns the group's
   * chi^2.
   */
  double ProfileThmGroup(int g, azure2::Results &results);
  /// Standardized residuals (f - d)/e of segment i under its group's profile in `results`, one per point.
  void ThmGroupResiduals(int g, int i, const azure2::Results &results, std::vector<double> &out);
  /// The background (model units) that the output adds to the model value of a point of segment i at c.m. energy
  /// E under the profile in `results`; 0 outside a group with one.
  double ThmBackgroundAt(int i, double energy, const azure2::Results &results);
  /// The lowest and highest c.m. energy (MeV) of the points of a group's segments: the grids of the set-up
  /// (ThmNorm.cpp) and the ranges of the reports (ThmReports.cpp).
  static void ThmGroupEnergyRange(WorkingCopy &data, const ThmGroup &group, double &eLo, double &eHi);
  /// Profile of every experiment, trivial ones included, for the output and pyazr.
  std::vector<ThmExperimentReport> ThmExperimentReports(const azure2::Results &results);
  /*!
   * The Coulomb line shape of THM experiment `name` (lineshape=on) at the
   * c.m. energies `energies` (MeV, x + A) and the current fit parameters of
   * `compound`: E_sF, eta_0, and per exit pair of its segments zeta, the
   * eta_sb estimate and, per level of the J groups coupling entrance and
   * exit, the pole and |N_C|^2.  False (and `why`) if there is no such
   * experiment or it has no line shape.
   */
  bool ThmLineshapeTable(const std::string &name, const std::vector<double> &energies, CNuc *compound,
                         const Config &configure, ThmLineshapeReport &out, std::string &why);
  /*!
   * The entrance vertex of THM experiment `name` at the c.m. energies
   * `energies` (MeV, x + A) and the current parameters: the nodes and weights
   * of its spectator-momentum window (one node for ps=delta), rho = p_xA a/hbar c
   * per node, and per entrance channel and level the window average of
   * |M_l|^2 and its quasi-free value (p_s = 0), with the boundary the vertex
   * uses for that level.  False (and `why`) if there is no such experiment.
   */
  bool ThmVertexTable(const std::string &name, const std::vector<double> &energies, CNuc *compound,
                      const Config &configure, ThmVertexReport &out, std::string &why);
  /*!
   * The distortion factor of THM experiment `name` (distortion=coulomb or
   * optical) at the c.m. energies `energies` (MeV, x + A): E_sF, k_sF,
   * eta_sF, the spectator angle, q, |M|^2, |M_PW|^2, the l summed, R
   * evaluated directly and R as the model uses it (interpolated on the grid).
   * For distortion=table only the model's w(E).  False (and `why`) if there is
   * no such experiment or it has no distortion.
   */
  bool ThmDistortionTable(const std::string &name, const std::vector<double> &energies, ThmDistortionReport &out,
                          std::string &why);
  /*!
   * Calculates every point of \p segments (segments of this data set) and of their leaves (src/EnergyPass.cpp):
   * the energy tasks of their plan tasks (CalculationPlan::EnergyTasks: one level-matrix solve per compound energy,
   * points and sub-points alike) in one OpenMP pass, largest first, then the folds of the points with sub-points in
   * a second.  A plan task whose calculation throws gets NaN on its points (Unit::MarkFailed).  Returns the number of
   * such plan tasks, or -1 when nothing was calculated: no OpenMP, already inside a parallel region, or
   * AZURE_SERIAL_CHI2 set.  After a pass that returns 0 the caller reads the points with
   * azure2::EvaluateModelPoint(..., LeafValues::Ready).
   */
  int CalculatePointsInParallel(const std::vector<const Segment *> &segments, CNuc *theCNuc, const Config &configure,
                                azure2::Results &results);
  /// The energy tasks of the point pass over \p segments (the ones CalculatePointsInParallel runs): the analytic
  /// gradient back-propagates through the same tasks (src/GradEnergyTasks.cpp).
  const azure2::EnergyTaskList &EnergyTasksOf(const std::vector<const Segment *> &segments, const Config &configure);

 private:
  void LayOutUnits(std::vector<std::vector<Unit> *> *fresh);  // CompactUnits, PlaceSubPointGrids
  /// The plan tasks of the points of \p segments and their leaves, in data order (src/EnergyPass.cpp).
  std::vector<int> PassWork(const std::vector<const Segment *> &segments) const;
  /// The data set (shared, immutable) and, during the set-up only, the same object writable (Draft).
  std::shared_ptr<const DataSet> set_;
  std::shared_ptr<DataSet> draft_;
  std::vector<SegmentState> segmentStates_;     // State: one per segment of the data set, by Pos
  std::vector<TargetEffectValues> effectValues_;  // EffectValues: one per target effect
  std::vector<double> thmCoherentValues_;       // ThmCoherentValue: one per coherent-background parameter
  std::vector<Unit> units_;  // Unit
  std::vector<std::vector<Unit>> pendingSubPoints_;  // AddSubPoint: per point's unit, until the commit
  azure2::CalculationPlan plan_;  // built by BuildPlan; positions only, so a copy copies it as it is
  std::map<std::string, std::shared_ptr<DataCovariance::File>> covFiles_;  // only while filling
  bool covarianceError_ = false;              // a covariance file could not be read (CovarianceFileFor)
  bool kinematicsError_ = false;              // a lab angle beyond the kinematic maximum (audit 2026-10-25, item 31)
  EvaluationCount iterations_;  // objective evaluations (EvaluationCount.h)
  /// Group the densityfit / densitygroup lines into parameters (SetUpTargetEffects); -1 on an error.
  int SetUpDensityParams(const Config &);
  /// Groups the effects' declarations into apparatus parameters by their key (ApparatusParameters.cpp).
  void SetUpApparatusParams(const Config &);
  /// The set-up's call (src/TargetSetup.cpp, the node-layout branch's file) under its former name.
  void SetUpStoppingScaleParams(const Config &configure) { SetUpApparatusParams(configure); }
  /// The set-up refusals of a library stopping (evaluation, table): its data, and its range at every sub-point
  /// (StoppingChecks.h); the source of each target's stopping in the log, once per target.  false after a refusal.
  bool CheckTargetStoppings(const Config &);
  bool CheckSubPointStoppings(const Config &);
  bool isFit_;
  bool isErrorAnalysis_;
  std::streampos ecReadPos_;  // File offset where component-segment EC integrals begin in the intEC file
  // Whether the intEC file is being read back rather than written.  Decided once, in
  // CalculateECAmplitudes, and followed by the component-segment pass so the two cannot disagree.
  bool ecUsePrevious_;
  std::string ecSignature_;   // Signature of the calculation being set up (see ECSignature)
  std::string ecOutputFile_;  // The intEC file being written, when not reading one back
  /// Builds thmGroups_ from the <thm> experiments once the segments are read; -1 after an ERROR.
  int BuildThmGroups(const Config &, CNuc *, int numLines);
  /// Writes output/thm_experiments.out.
  void WriteThmExperiments(const Config &, const azure2::Results &results);
  /// The parts of WriteOutputFiles (src/DataOutputFiles.cpp), in the order it calls them.  The AZUREOut_* files with
  /// their .band siblings (when `band` holds gradients) and the sample decompositions; kinflag 0 writes the c.m.
  /// frame, 1 the lab frame (kinflag.dat).
  void WriteSegmentFiles(const Config &, const azure2::Results &results, int kinflag, const BandData *band);
  /// chiSquared.out: one row per segment, in the order of the segment files, then the totals.
  void WriteChiSquaredFile(const Config &, const azure2::Results &results);
  /// normalizations.out, when a segment's normalization is varied.
  void WriteNormalizationsFile(const Config &, const azure2::Results &results);
  /// shifts.out, when a segment's energy shift is varied, follows a group or has a shift function.
  void WriteShiftsFile(const Config &, const azure2::Results &results);
  /// Slots of a buffer laid out for this data set: one record per point and sub-point of its segments and component
  /// segments (Results.h, LayoutResults; the data set holds the layout, never the values), and its model values
  /// (Results::model).
  std::size_t resultsSize_ = 0, modelSize_ = 0;
  /// One record per point and sub-point, in the order of the segments' points then the leaves', each point before its
  /// sub-points: each gets its offset (Unit::SetResultRecord); and one model record per point of segments_
  /// (Unit::SetModelRecord).  With `results`, the buffer of the layout the points
  /// had, a point's record keeps the values it held there; a point new to the layout starts at 0, as a newly
  /// constructed point did.
  void LayoutResults(azure2::Results *results);
  /// Calls f(point, angular slots, natural slots) for every point and sub-point in layout order.
  template <class F>
  void ForEachResultUnit(F f);
  void FillEnergyShiftsFromParamsImpl(const vector_r &, CNuc *theCNuc, const Config *configure);
  /// The external-capture amplitudes of the leaves (InitializeLeaves): computed and appended to `out`, or read from
  /// `in` where the data segments' amplitudes end (usePrevious); -1 when stopped or a GSL error ends the loop.
  int LeafECAmplitudes(CNuc *, const Config &, bool usePrevious, std::ifstream &in, std::ofstream &out);
  /// Redraws the progress bar of the external-capture loops (CalculateECAmplitudes, InitializeLeaves)
  /// from inside their OpenMP loops, under the named critical section ec_progress.
  static void ECProgressRedraw(const Config &configure, const char *segmentKeyOut, const std::atomic<int> &pointIndex,
                               int numPoints, time_t &startTime);
  // A copy of the compound nucleus that holds the observed level parameters the
  // sub-point grids are laid out from (ReplaceSubPointGrids).  Never shared: a
  // copy or Clone of the WorkingCopy starts without one.
  struct GridNucleus {
    std::unique_ptr<CNuc> nuc;
    GridNucleus();
    GridNucleus(const GridNucleus &);
    GridNucleus &operator=(const GridNucleus &);
    ~GridNucleus();
  };
  GridNucleus gridNuc_;
  // The tables of re-placed sub-point nodes, shared by this working copy and its copies, and written by them during
  // an evaluation: on purpose, each entry under its shard's mutex, never a table in place (NodeTableCache.h).
  std::shared_ptr<NodeTableCache> nodeCache_ = std::make_shared<NodeTableCache>();
  /// The stage's channel table (docs/upgrade/channel_table.md): built by CalcEDependentValues and
  /// InitializeLeaves before their point loops, read-only afterwards, replaced, never written.
  std::shared_ptr<const ChannelTable> channelTable_;
  azure2::NodeState nodes_;  // Nodes(), Slots(), Overflow(): built at the end of Initialize, placed by the hand-over
  std::shared_ptr<const ChannelTable> BuildChannelTable(CNuc *, const Config &, const std::vector<Unit *> &points);
};

#include "SubPointAccess.h"

#endif
