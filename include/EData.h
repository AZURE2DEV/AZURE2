#ifndef EDATA_H
#define EDATA_H

#include "ESegment.h"
#include "TargetEffect.h"
#include "EDataIterator.h"
#include "ThmExperiment.h"
#include <memory>
#include <deque>
#include <map>
#include "DataCovariance.h"
#include <memory>
#include <ios>

// THM types EData holds by pointer or takes by reference (EDataThm.cpp).
struct ThmLineshape;
struct ThmSpectatorWindow;
struct ThmAngleWindow;
struct ThmVertexReport;
struct ThmLineshapeReport;
struct ThmDistortionReport;
class ThmDistortion;
class ThmDwVertex;

class CNuc;
class PriorCovariance;
struct BandData;
namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT

/// An AZURE data object

/*!
 * The EData object is the top level data object in AZURE.  It is the container object for a vector of ESegment objects.
 */

class EData {
 public:
  EData();
  /// Number of data segments.
  int NumSegments() const;
  /// Build the segments from the segment and data files. Returns -1 on failure.
  int Fill(const Config &, CNuc *);
  /// Build points from <segmentsTest> grids instead of data files, for a run with no data. Returns -1 on failure.
  int MakePoints(const Config &, CNuc *);
  /// Fit iterations used so far.
  int Iterations() const;
  /// Number of target-effect definitions attached to the data.
  int NumTargetEffects() const;
  /// Index at which the normalization parameters start in the Minuit vector.
  int GetNormParamOffset() const;
  /// Index at which the energy-shift parameters start in the Minuit vector.
  int GetEnergyShiftParamOffset() const;
  /*!
   * Explicit prior centres of normalizations and energy shifts: rows
   *   segment_N_norm prior_centre c
   *   segment_N_energy_shift prior_centre c
   * of <parameterSettings> (N the segment key, the 1-based <segmentsData> line
   * counting inactive lines).  Such a row moves only the centre of that prior
   * (for a norm also the scale of its percentage width); the norm/shift field
   * of the segment line stays the start value.  Without a row the field is
   * both, as it always was.  Three-token rows are skipped by the older
   * <parameterSettings> readers (engine and GUI), which then fall back to the
   * field.  Returns -1 on a malformed row.
   */
  int ReadPriorCentres(const Config &);
  /// Index at which the sqrt(E) energy-shift coefficients start in the Minuit vector.
  int GetEnergyShiftSqrtParamOffset() const;
  /// Index at which the shift-function parameters start in the Minuit vector (after the sqrt(E)
  /// coefficients): the parameters of each segment that declares a "shiftfn" block, in segment order.
  int GetEnergyShiftFnParamOffset() const;
  /*!
   * Electron screening.  The <targetInt> lines with a `screening` block, in line order, each with one
   * parameter `screening_<n>` (n = the effect number: its place among the active <targetInt> lines),
   * U_e in eV, after the shift-function parameters.  None without a screening block.  A line listing
   * several segments has one copy of its effect per segment (ReadTargetEffectsFile); the parameter is
   * the line's, shared by every copy.
   */
  int NumScreeningParams() const { return (int)screeningEffects_.size(); }
  /// Effect number of screening parameter k (1-based): its place among the active <targetInt> lines.
  int GetScreeningEffect(int k) const { return screeningEffects_[k - 1]; }
  /// The target effect of screening parameter k (the first copy of its line; every copy holds the same U_e).
  TargetEffect *GetScreeningTargetEffect(int k) { return this->GetTargetEffect(screeningCopies_[k - 1].front()); }
  const TargetEffect *GetScreeningTargetEffect(int k) const {
    return const_cast<EData *>(this)->GetTargetEffect(screeningCopies_[k - 1].front());
  }
  int GetScreeningParamOffset() const { return screeningParamOffset_; }
  /// Push each screening effect's current U_e onto the points and sub-points of its segments.
  void ApplyScreening();
  /// Set the screening effects' U_e from a parameter vector, then ApplyScreening (called by
  /// FillEnergyShiftsFromParams, so every evaluation sees the vector's values).
  void FillScreeningFromParams(const vector_r &);
  /// Sum of the screening priors ((U_e - nominal)/sigma)^2.
  double ScreeningPenalty() const;
  /*!
   * Fitted target densities (TargetEffect::HasDensityFit).  One parameter `density_<n>` (atoms/cm^2) per
   * ungrouped `densityfit` line and per `densitygroup`, n = the effect number of the line (of the group's first
   * line: its place among the active <targetInt> lines), after the screening parameters.  None without a
   * densityfit block.  As for screening, a line listing several segments has one copy of its effect per segment;
   * the parameter is the line's, shared by every copy.
   */
  int NumDensityParams() const { return (int)densityEffects_.size(); }
  /// Effect number of density parameter k (1-based): the line, or the first line of its group.
  int GetDensityEffect(int k) const { return densityEffects_[k - 1]; }
  /// The target effect of density parameter k (the first copy of its line; every copy holds the same density).
  TargetEffect *GetDensityTargetEffect(int k) { return this->GetTargetEffect(densityLeadCopies_[k - 1]); }
  const TargetEffect *GetDensityTargetEffect(int k) const {
    return const_cast<EData *>(this)->GetTargetEffect(densityLeadCopies_[k - 1]);
  }
  int GetDensityParamOffset() const { return densityParamOffset_; }
  /// The density parameter (1-based) effect \p n follows, 0 if its density is not fitted.
  int DensityParamOfEffect(int n) const {
    return (n >= 1 && n <= (int)densityParamOfEffect_.size()) ? densityParamOfEffect_[n - 1] : 0;
  }
  /// Set the fitted densities from a parameter vector and recompute the target thickness of every point of
  /// their lines (called by FillEnergyShiftsFromParams, so every evaluation sees the vector's values).
  void FillDensitiesFromParams(const vector_r &, const Config &);
  /// Sum of the density priors ((N - N0)/sigma)^2.
  double DensityPenalty() const;
  /*!
   * Fitted Doppler temperatures (TargetEffect::HasDopplerFit, `doppler T vary sigma [Tmax]`).  One parameter
   * `doppler_<n>` (K) per such line, n = the line number (its place among the active <targetInt> lines, as
   * screening_<n>), after the density parameters.  None without the fit form.
   */
  int NumDopplerParams() const { return (int)dopplerEffects_.size(); }
  /// Line number of Doppler parameter k (1-based): the n of doppler_<n>.  THM gives every segment a line lists its
  /// own copy of the effect; the parameter moves all of them (dopplerCopies_), as screening_<n> does.
  int GetDopplerEffect(int k) const { return dopplerEffects_[k - 1]; }
  /// The first copy of Doppler parameter k's line (its temperature, sigma and range are the line's).
  TargetEffect *GetDopplerTargetEffect(int k) { return this->GetTargetEffect(dopplerCopies_[k - 1].front()); }
  const TargetEffect *GetDopplerTargetEffect(int k) const {
    return const_cast<EData *>(this)->GetTargetEffect(dopplerCopies_[k - 1].front());
  }
  int GetDopplerParamOffset() const { return dopplerParamOffset_; }
  /// The Doppler parameter (1-based) of effect \p n, 0 if its temperature is not a fit parameter.
  int DopplerParamOfEffect(int n) const;
  /// Set the fitted temperatures from a parameter vector (clamped to [T0/4, layout temperature], with a warning
  /// printed once), called by FillEnergyShiftsFromParams so every evaluation sees the vector's values.
  void FillDopplerFromParams(const vector_r &, const Config *configure);
  /// Sum of the Doppler priors ((T - T0)/sigma)^2.
  double DopplerPenalty() const;
  /*!
   * Fitted sample thicknesses (ESegment::HasThicknessFit, `transmission n vary sigma` / `yield n vary sigma`).  One
   * parameter `segment_<key>_thickness` (atoms/barn) per data segment whose line has the fit form, in segment order,
   * after the Doppler temperatures; a thickness-group follower's is fixed at its leader's value.  None without the
   * fit form.
   */
  int NumThicknessParams() const { return (int)thicknessSegments_.size(); }
  /// Position (1-based, as GetSegment) of the segment owning thickness parameter k (1-based).
  int GetThicknessSegment(int k) const { return thicknessSegments_[k - 1]; }
  int GetThicknessParamOffset() const { return thicknessParamOffset_; }
  /// The thickness parameter (1-based) whose value segment \p pos (1-based) uses: its own, its group leader's for a
  /// follower; 0 when that thickness is no parameter.
  int ThicknessParamOfSegment(int pos) const;
  /// Set the fitted thicknesses from a parameter vector (a value below 0 is taken as 0); thickness-group followers
  /// take their leader's.  Called by FillEnergyShiftsFromParams.
  void FillThicknessFromParams(const vector_r &);
  /// Sum of the thickness priors ((n - n0)/sigma)^2.
  double ThicknessPenalty() const;
  /// Link the thicknessgroup members to the group's first segment in file order; -1 (with an ERROR) when a member has
  /// no transmission / yield block.  Also lays out the thickness parameters (NumThicknessParams).
  int ResolveThicknessGroups(const Config &);
  /// Read the target-effects input file and build the TargetEffect objects.
  int ReadTargetEffectsFile(const Config &, CNuc *);
  /// Is this a fit? AZURECalc clones the compound nucleus and data per thread only when it is.
  bool IsFit() const;
  /// Is this a Minos error-analysis call? Suppresses the transformation and file output.
  bool IsErrorAnalysis() const;
  /// Does a segment with this key exist?
  bool IsSegmentKey(int);
  void SetFit(bool);
  void SetErrorAnalysis(bool);
  /// Count one more fit iteration.
  void Iterate();
  /// Reset the iteration counter to zero.
  void ResetIterations();
  /// Initialize every point: EPoint::Initialize over the whole data set.
  int Initialize(CNuc *, const Config &);
  /// Append a segment.
  void AddSegment(ESegment);
  /// Print the data as read, or as generated for a run without data.
  void PrintData(const Config &);
  /// EPoint::CalcLegendreP for every point.
  void CalcLegendreP(int, CNuc *);
  void PrintLegendreP(const Config &);
  /// EPoint::CalcEDependentValues for every point -- penetrabilities, shifts and Coulomb phases at each energy.
  int CalcEDependentValues(CNuc *, const Config &);
  void PrintEDependentValues(const Config &, CNuc *);
  /// EPoint::CalcCoulombAmplitude for every point.
  void CalcCoulombAmplitude(CNuc *);
  void PrintCoulombAmplitude(const Config &, CNuc *);
  /// Write AZUREOut_*, chiSquared.out and the rest of the run's output files.
  void WriteOutputFiles(const Config &, bool = false, const BandData * = nullptr);
  /// Park formalism: the J > 0 penalty of the last objective evaluation, for
  /// chiSquared.out (set by AZURECalc::operator(); 0 otherwise).
  void SetParkPenalty(double penalty) { parkPenalty_ = penalty; }
  double GetParkPenalty() const { return parkPenalty_; }
  /// The abundance penalty (priors + sum-to-one) of the last evaluation, reported in chiSquared.out.
  void SetAbundancePenalty(double penalty) { abundancePenalty_ = penalty; }
  double GetAbundancePenalty() const { return abundancePenalty_; }
  /// chi2 of the measured scattering lengths (<scatteringLengths>) at the last evaluation.
  void SetScatteringLengthChi2(double chi2) { scatteringLengthChi2_ = chi2; }
  double GetScatteringLengthChi2() const { return scatteringLengthChi2_; }
  /// The prior parameter covariance of the run (<config> priorCovariance=, PriorCovariance.h), shared by every copy
  /// (Clone); null without one.  Set by the CLI (AZUREMain) and pyazr (AZUREAPI) once the parameters are known.
  void SetPriorCovariance(std::shared_ptr<const PriorCovariance> prior) { priorCov_ = std::move(prior); }
  const PriorCovariance *GetPriorCovariance() const { return priorCov_.get(); }
  /// Its chi2 term at the last evaluation, reported in chiSquared.out (Total-Prior-Covariance-Chi-Squared).
  void SetPriorCovarianceChi2(double chi2) { priorCovChi2_ = chi2; }
  double GetPriorCovarianceChi2() const { return priorCovChi2_; }
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
  /// Set up the entrance/exit pairs of the component segments in the compound nucleus.
  int InitializeComponentSegments(CNuc *, const Config &);
  /// Validates a component's pair keys (both present, one compound nucleus) and flags its entrance pair; false with the error printed otherwise.
  bool CheckComponentPairs(CNuc *, const ESegment &baseSegment, int entranceKey, int exitKey, const Config &);
  ESegment *CreateComponentSegment(const ESegment &baseSegment, int entranceKey, int exitKey);
  /// A transmission or yield segment after its components exist: checks, plain components, the yield's total term.
  bool FinishTransformedSegment(CNuc *, ESegment *, int key, bool isTest, const Config &);
  /// Every term of a transmission or yield segment on the base point's sub-point grid (after the target effects).
  bool CheckTransformedSegments(const Config &);
  ESegment *CreateComponentSegment(const ESegment &baseSegment, int entranceKey, int exitKey, double fixedAngle);
  /// Map points at equal energies onto one another, so a shared energy is calculated once.
  /// The compound nucleus tells which segments have a particle exit; null is conservative.
  void MapData(CNuc *theCNuc = nullptr);
  /// The mapping inside segments with components (MapData calls it last).
  void MapComponentSegments(CNuc *theCNuc);
  /// Append a target-effect definition.
  void AddTargetEffect(TargetEffect);
  void SetNormParamOffset(int);
  void SetEnergyShiftParamOffset(int);
  void SetEnergyShiftSqrtParamOffset(int);
  void SetEnergyShiftFnParamOffset(int);
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
  /// A data point beyond the kinematic maximum lab angle (ESegment::Fill): the run stops instead of dropping the segment.
  void SetKinematicsError() { kinematicsError_ = true; }
  /// Gather the points of the segments carrying "covariance <file>" into one group per
  /// file and factor each group's matrix (called by Fill).  Returns -1 on an error.
  int ResolveCovarianceGroups(const Config &);
  bool HasCovariance() const { return !covGroups_.empty(); }
  int NumCovarianceGroups() const { return (int)covGroups_.size(); }
  DataCovariance::Group &GetCovarianceGroup(int g) const { return *covGroups_[g]; }
  /// From the points' current fit values: the chi2 of every covariance group, r^T M^-1 r,
  /// minus the sum of r_i^2 the per-point loops counted for those points (so adding the
  /// return value to their total gives the covariance chi2).  Sets each covariance
  /// segment's chi2 to the sum over its points of r_i (M^-1 r)_i (also in \p segChis,
  /// indexed by segment number, when given) and the points' seeds (M^-1 r)_i.
  double CovarianceChi2Correction(std::vector<double> *segChis = nullptr);
  /// Whiten the data rows of a residual vector (and of a row-major Jacobian with \p nCols
  /// columns, when given): r <- W r within every covariance group (W^T W = M^-1).  Rows are the
  /// points in segment order, as ComputeResidualJacobian lays them out.
  void WhitenResidualRows(std::vector<double> &residuals, std::vector<double> *jacobian, int nCols);
  /// covariance.out: each group's chi2 and the diagonal-only sum, at the current fit values.
  void WriteCovarianceSummary(const Config &);
  /// Forget the frozen model of the relvector terms; the next evaluation freezes it again.
  void RefreezeCovariance();
  /// True when any data-covariance group has relvector terms (a model frozen at the first evaluation).
  bool HasRelativeCovariance() const;
  /// Calculate the covariance segments' points at the current parameters and set their
  /// seeds, for a gradient pass that needs (M^-1 r) before it reaches the points.
  void PrepareCovarianceSeeds(CNuc *, const Config &);
  /// The segment (1-based) whose norm segment s uses: its norm-group leader, or s itself.
  int NormOwner(int s) const;
  void FillMnParams(ROOT::Minuit2::MnUserParameters &);
  /// Write the normalizations back from a Minuit parameter vector.
  void FillNormsFromParams(const vector_r &);
  /// Write the energy shifts back from a Minuit parameter vector.
  void FillEnergyShiftsFromParams(const vector_r &, EData *data = nullptr, CNuc *theCNuc = nullptr, const Config *configure = nullptr);
  /// Lay the adaptive sub-point grids out again for the levels theCNuc holds now
  /// (its fit values, as observed parameters) and the segments' energy shifts, and
  /// rebuild the sub-points of every point whose grid changed.  Sub-points of a
  /// segment in \p shiftActive are left for the energy-shift step to energise and
  /// their segments are appended to \p reshift; the others are energised here.
  /// Returns the number of points whose grid changed.  Called by
  /// FillEnergyShiftsFromParams.
  int ReplaceSubPointGrids(CNuc *theCNuc, const Config &configure, const std::vector<ESegment *> &shiftActive,
                           std::vector<ESegment *> &reshift);
  /// Drop the last segment.
  void DeleteLastSegment();
  /// Segment \p i, 1-based.
  ESegment *GetSegment(int);
  /// The segment with this key, or null if there is none.
  ESegment *GetSegmentFromKey(int);
  EData *Clone() const;

  /*!
   * A THM experiment (<thm> experiment[<name>] ..., ThmExperiment.h) as the
   * data see it: the segments that share one profiled norm and background.
   */
  struct ThmGroup {
    std::string name;
    int terms = 0;              ///< background terms, 0..3
    std::vector<int> segments;  ///< 1-based indices into the segments, ascending
    /// One segment and no background: exactly the per-segment profile of
    /// ESegment::ProfileNormChiSquared, which then handles it.
    bool trivial = false;
    ThmProfile profile;         ///< the last profile (ProfileThmGroup)
    /// Coulomb line shape (lineshape=on), or null; shared with its segments.
    std::shared_ptr<const ThmLineshape> lineshape;
    /// Spectator-momentum window (ps=hulthen|gauss|table), or null; shared
    /// with its segments.
    std::shared_ptr<const ThmSpectatorWindow> window;
    /// Distortion factor R(E) (distortion=...), or null; shared with its segments.
    std::shared_ptr<const ThmDistortion> distortion;
    /// Distorted-wave entrance vertex (vertexModel=dw), or null; shared with
    /// its segments.  It replaces R(E): distortion is then null.
    std::shared_ptr<const ThmDwVertex> dwVertex;
    /// Angular window of the exit pair (theta=thmin-thmax), or null; shared
    /// with its segments.
    std::shared_ptr<const ThmAngleWindow> angle;
    /// Entrance pair key of its segments (0 if they differ).
    int pairKey = 0;
    /// Coherent background (cbackground=), or null; shared with its segments.
    std::shared_ptr<const ThmCoherentBackground> coherent;
  };
  /*!
   * The fit parameters of the THM coherent backgrounds (cbackground=), the
   * last block of the parameter vector after the energy shifts (none without
   * the key): Re c0, Im c0 [, Re c1, Im c1] of each combination, in the order
   * of the experiments and their terms.  `value` is the current value.
   */
  struct ThmCoherentParam {
    std::string name;
    std::string experiment;
    double value = 0.0;
    bool fixed = false;
  };
  int NumThmCoherentParams() const { return (int)thmCoherentParams_.size(); }
  const ThmCoherentParam &GetThmCoherentParam(int k) const { return thmCoherentParams_[k]; }
  int GetThmCoherentParamOffset() const { return thmCoherentParamOffset_; }
  /// Reads the coherent-background block of a full parameter vector (if it
  /// has one) into the parameters and, if given, the compound the HOES model
  /// reads them from.  Called by FillEnergyShiftsFromParams.
  void FillThmCoherentFromParams(const vector_r &p, CNuc *theCNuc);
  int NumThmGroups() const { return (int)thmGroups_.size(); }
  const ThmGroup &GetThmGroup(int g) const { return thmGroups_[g]; }
  /// The (non-trivial) THM experiment that profiles segment i jointly with
  /// others or with a background, or -1: then the segment is handled as
  /// before (its own profile, or its fixed norm).  Only a segment whose norm
  /// is profiled (ESegment::IsProfiledNorm) belongs.
  int ThmGroupOf(int i);
  /// Is i the last segment of group g, where its profile is taken once the
  /// models of all its segments are in?
  bool IsLastOfThmGroup(int g, int i);
  /*!
   * Profiles group g from the fit cross sections of its segments' points:
   * stores the profile, sets every segment's norm to the shared n* and its
   * chi^2 to that of its own points, and returns the group's chi^2.
   */
  double ProfileThmGroup(int g);
  /// Standardized residuals (f - d)/e of segment i under its group's last profile, one per point.
  void ThmGroupResiduals(int g, int i, std::vector<double> &out);
  /// The background (model units) that the output adds to the fit cross
  /// section of a point of segment i at c.m. energy E; 0 outside a group with one.
  double ThmBackgroundAt(int i, double energy);
  /// Profile of every experiment, trivial ones included, for the output and pyazr.
  std::vector<ThmExperimentReport> ThmExperimentReports();
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
  TargetEffect *GetTargetEffect(int);
  EDataIterator begin();
  EDataIterator end();
  std::vector<ESegment> &GetSegments();
  /*!
   * Calculates every point of \p segments and of their component segments in
   * one OpenMP pass, largest items first.  A host point also fills the points
   * mapped onto it, so mapped points are not items of their own.  A point whose
   * calculation throws gets a fit cross section of 0.  Returns the number of
   * such points, or -1 when nothing was calculated: no OpenMP, already inside a
   * parallel region, point mapping switched off on this thread, or
   * AZURE_SERIAL_CHI2 set.  The caller then reads the points with
   * ESegment::CalculateTheoreticalCrossSection(..., precalculated = true).
   */
  static int CalculatePointsInParallel(const std::vector<ESegment *> &segments, CNuc *theCNuc,
                                       const Config &configure);

 private:
  std::vector<TargetEffect> targetEffects_;
  std::vector<ESegment> segments_;
  std::deque<ESegment> componentSegments_;  // Separate storage for component segments (deque avoids pointer invalidation)
  bool componentError_ = false;
  std::map<std::string, std::shared_ptr<DataCovariance::File>> covFiles_;  // only while filling
  std::vector<std::shared_ptr<DataCovariance::Group>> covGroups_;         // shared by every copy
  bool covarianceError_ = false;              // a component failed CheckComponentPairs
  bool kinematicsError_ = false;              // a lab angle beyond the kinematic maximum (audit 2026-10-25, item 31)
  int iterations_;
  int normParamOffset_;
  int energyShiftParamOffset_;
  int energyShiftSqrtParamOffset_;
  int energyShiftFnParamOffset_ = 0;
  int screeningParamOffset_ = -1;  // -1 until FillMnParams has laid the vector out
  std::vector<int> screeningEffects_;  // line numbers (effect numbers) with a screening block, in order
  std::vector<std::vector<int>> screeningCopies_;  // per screening parameter: the target effects of its line
  int densityParamOffset_ = -1;        // -1 until FillMnParams has laid the vector out
  std::vector<int> densityEffects_;    // line number (effect number) holding each density parameter, in order
  std::vector<int> densityLeadCopies_;     // per density parameter: the first target effect (copy) of that line
  std::vector<int> densityParamOfEffect_;  // per target effect (copy, 0-based): its density parameter, 0 = none
  int dopplerParamOffset_ = -1;        // -1 until FillMnParams has laid the vector out
  std::vector<int> dopplerEffects_;    // line number holding each Doppler parameter, in order (as screening_)
  std::vector<std::vector<int>> dopplerCopies_;  // per Doppler parameter: the target effects (copies) of its line
  int thicknessParamOffset_ = -1;      // -1 until FillMnParams has laid the vector out
  std::vector<int> thicknessSegments_;  // positions (1-based) of the segments with a thickness parameter, in order
  /// Group the densityfit / densitygroup lines into parameters (ReadTargetEffectsFile); -1 on an error.
  int SetUpDensityParams(const Config &);
  bool isFit_;
  bool isErrorAnalysis_;
  std::streampos ecReadPos_;  // File offset where component-segment EC integrals begin in the intEC file
  // Whether the intEC file is being read back rather than written.  Decided once, in
  // CalculateECAmplitudes, and followed by the component-segment pass so the two cannot disagree.
  bool ecUsePrevious_;
  std::string ecSignature_;   // Signature of the calculation being set up (see ECSignature)
  std::string ecOutputFile_;  // The intEC file being written, when not reading one back
  std::vector<ThmGroup> thmGroups_;  // THM experiments (BuildThmGroups), empty without them
  std::vector<ThmCoherentParam> thmCoherentParams_;  // cbackground= parameters, empty without the key
  int thmCoherentParamOffset_ = -1;                 // their first index in the parameter vector (FillMnParams)
  /// Builds thmGroups_ from the <thm> experiments once the segments are read; -1 after an ERROR.
  int BuildThmGroups(const Config &, CNuc *, int numLines);
  /// Writes output/thm_experiments.out.
  void WriteThmExperiments(const Config &);
  double parkPenalty_ = 0.0;
  double abundancePenalty_ = 0.0;
  double scatteringLengthChi2_ = 0.0;
  std::shared_ptr<const PriorCovariance> priorCov_;  // shared by every copy (read-only once set up)
  double priorCovChi2_ = 0.0;
  std::vector<std::string> nucleusLabels_;   // one per compound nucleus, for the natural decomposition file
  // A copy of the compound nucleus that holds the observed level parameters the
  // sub-point grids are laid out from (ReplaceSubPointGrids).  Never shared: a
  // copy or Clone of the EData starts without one.
  struct GridNucleus {
    std::unique_ptr<CNuc> nuc;
    GridNucleus();
    GridNucleus(const GridNucleus &);
    GridNucleus &operator=(const GridNucleus &);
    ~GridNucleus();
  };
  GridNucleus gridNuc_;
};

#endif
