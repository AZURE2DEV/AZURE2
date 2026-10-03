#ifndef EDATA_H
#define EDATA_H

#include "ESegment.h"
#include "TargetEffect.h"
#include "EDataIterator.h"
#include "ThmExperiment.h"
#include "ThmLineshape.h"
#include "ThmDistortion.h"
#include "ThmDwVertex.h"
#include "ThmAngular.h"
#include <memory>
#include <deque>
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
  ESegment *CreateComponentSegment(const ESegment &baseSegment, int entranceKey, int exitKey);
  ESegment *CreateComponentSegment(const ESegment &baseSegment, int entranceKey, int exitKey, double fixedAngle);
  /// Map points at equal energies onto one another, so a shared energy is calculated once.
  void MapData();
  /// Append a target-effect definition.
  void AddTargetEffect(TargetEffect);
  void SetNormParamOffset(int);
  void SetEnergyShiftParamOffset(int);
  /// Seed the Minuit parameter array with the normalizations and energy shifts.
  void FillMnParams(ROOT::Minuit2::MnUserParameters &);
  /// Write the normalizations back from a Minuit parameter vector.
  void FillNormsFromParams(const vector_r &);
  /// Write the energy shifts back from a Minuit parameter vector.
  void FillEnergyShiftsFromParams(const vector_r &, EData *data = nullptr, CNuc *theCNuc = nullptr, const Config *configure = nullptr);
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

 private:
  std::vector<TargetEffect> targetEffects_;
  std::vector<ESegment> segments_;
  std::deque<ESegment> componentSegments_;  // Separate storage for component segments (deque avoids pointer invalidation)
  int iterations_;
  int normParamOffset_;
  int energyShiftParamOffset_;
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
};

#endif
