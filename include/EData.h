#ifndef EDATA_H
#define EDATA_H

#include "ESegment.h"
#include "TargetEffect.h"
#include "EDataIterator.h"
#include <deque>
#include <map>
#include "DataCovariance.h"
#include <memory>
#include <ios>

class CNuc;
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
  /// Index at which the sqrt(E) energy-shift coefficients start in the Minuit vector.
  int GetEnergyShiftSqrtParamOffset() const;
  /// Index at which the shift-function parameters start in the Minuit vector (after the sqrt(E)
  /// coefficients): the parameters of each segment that declares a "shiftfn" block, in segment order.
  int GetEnergyShiftFnParamOffset() const;
  /*!
   * Electron screening.  The target effects with a `screening` block, in effect order, each with one
   * parameter `screening_<n>` (n = the effect number: its place among the active <targetInt> lines),
   * U_e in eV, after the shift-function parameters.  None without a screening block.
   */
  int NumScreeningParams() const { return (int)screeningEffects_.size(); }
  /// Effect number of screening parameter k (1-based).
  int GetScreeningEffect(int k) const { return screeningEffects_[k - 1]; }
  int GetScreeningParamOffset() const { return screeningParamOffset_; }
  /// Push each screening effect's current U_e onto the points and sub-points of its segments.
  void ApplyScreening();
  /// Set the screening effects' U_e from a parameter vector, then ApplyScreening (called by
  /// FillEnergyShiftsFromParams, so every evaluation sees the vector's values).
  void FillScreeningFromParams(const vector_r &);
  /// Sum of the screening priors ((U_e - nominal)/sigma)^2.
  double ScreeningPenalty() const;
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
  TargetEffect *GetTargetEffect(int);
  EDataIterator begin();
  EDataIterator end();
  std::vector<ESegment> &GetSegments();

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
  std::vector<int> screeningEffects_;  // effect numbers with a screening block, in order
  bool isFit_;
  bool isErrorAnalysis_;
  std::streampos ecReadPos_;  // File offset where component-segment EC integrals begin in the intEC file
  // Whether the intEC file is being read back rather than written.  Decided once, in
  // CalculateECAmplitudes, and followed by the component-segment pass so the two cannot disagree.
  bool ecUsePrevious_;
  std::string ecSignature_;   // Signature of the calculation being set up (see ECSignature)
  std::string ecOutputFile_;  // The intEC file being written, when not reading one back
  double parkPenalty_ = 0.0;
  double abundancePenalty_ = 0.0;
  double scatteringLengthChi2_ = 0.0;
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
