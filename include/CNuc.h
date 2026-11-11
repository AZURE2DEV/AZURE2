#ifndef CNUC_H
#define CNUC_H

#include <string>
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "JGroup.h"
#include "PPair.h"
#include "PathwayTable.h"
#include "ModelData.h"
#include "azure2/experiment/SetupSpecs.h"

namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT
class Config;
class AZUREParams;
namespace azure2 {
namespace io {
class LevelRow;
}
}  // namespace azure2

/// An AZURE compound nucleus

/*!
 * The compound nucleus is the fundamental concept of R-Matrix theory.  As such, the CNuc object
 * in AZURE is the top level container object for all structure and reaction objects.  Specifically,
 * the CNuc object is the container object for vectors of PPair and JGroup objects, within which all other
 * nuclear data objects are contained.
 */

/// One compound nucleus of a project: the (Z, A) its pairs form.
struct CompoundNucleus {
  int Z = 0;
  int A = 0;
};

/// The isotopic abundance of one compound nucleus's target in one sample.
struct SampleAbundance {
  double abundance = 1.0;   ///< current value (a fit parameter when vary is set)
  double nominal = 1.0;     ///< value read from the file, centre of the prior
  double sigma = 0.0;       ///< prior width; 0 = no prior
  bool vary = false;
  bool priorReplaced = false;  ///< the prior covariance replaces this prior (docs/upgrade/objective.md D3)
  std::string isotope;  ///< format 3: the isotope as written, which names the parameter ("" = the old name)
};

/// A sample: one set of isotopic (or stoichiometric) abundances, one entry per
/// compound nucleus.  The <abundances> block holds one unnamed sample (the
/// lines without a "<name>:" prefix) and any number of named ones; a segment
/// chooses its sample with "natural 1" (the unnamed, or the only, sample) or
/// "natural <name>".  Each sample has its own sum-to-one penalty width: the
/// block's "<name>: sum <s>" line, else the config's abundanceSumSigma.
struct AbundanceSample {
  std::string name;                      ///< "" for the unnamed sample
  std::vector<SampleAbundance> entries;   ///< index nucleus-1
  double sumSigma = -1.0;                ///< < 0: use the config's abundanceSumSigma
};

class CNuc {
 public:
  // -- structure ------------------------------------------------------------
  /// Number of compound nuclei the pairs form.  One for every project written
  /// before several nuclei could share a file.
  int NumNuclei() const;
  /// Nucleus \p i, 1-based.
  CompoundNucleus *GetNucleus(int);
  const CompoundNucleus *GetNucleus(int) const;
  /// "Z=10 A=23" style label of nucleus \p i.
  std::string NucleusLabel(int) const;
  /// Did the project carry an <abundances> block?  Only then are the
  /// abundances Minuit parameters (one per nucleus and sample, after the
  /// level block), so a project without the block keeps its parameter layout.
  bool HasAbundanceBlock() const { return hasAbundanceBlock_; }
  /// The model data (ModelData.h): the rows of the objective's term T9, read by Fill.
  const ModelDataSet &GetModelData() const { return *modelData_; }
  /// The lines of the optional <scatteringLengths> block as read (ScatteringLength::ReadConstraints), for the report.
  const std::vector<ScatteringLength::Constraint> &ScatteringLengthConstraints() const {
    return modelData_->scatteringLengthLines;
  }
  /// Number of samples the <abundances> block defines (0 without the block).
  int NumSamples() const { return (int)samples_.size(); }
  /// Sample \p s, 1-based.
  AbundanceSample *GetSample(int s) { return &samples_[s - 1]; }
  const AbundanceSample *GetSample(int s) const { return &samples_[s - 1]; }
  /// 1-based index of the sample with this name, 0 if there is none.
  int FindSample(const std::string &name) const;
  /// The sample "natural 1" means: the unnamed one if present, else the only one; 0 if ambiguous or none.
  int DefaultSample() const;
  /// Number of abundance parameters FillMnParams emits: samples x nuclei with the block, else 0.
  int NumAbundanceParams() const { return hasAbundanceBlock_ ? (int)(samples_.size() * nuclei_.size()) : 0; }
  /// 1-based ordinal of the abundance parameter of (sample, nucleus), both 1-based: samples are blocks of nuclei.
  int AbundanceOrdinal(int sample, int nucleus) const { return (sample - 1) * (int)nuclei_.size() + nucleus; }
  /// The (sample, nucleus) of a 1-based abundance ordinal.
  void AbundanceOf(int ordinal, int &sample, int &nucleus) const;
  /// Minuit name: "abundance_<Z>_<A>" for the unnamed sample, "abundance_<sample>_<Z>_<A>" otherwise.
  std::string AbundanceParamName(int sample, int nucleus) const;
  /// Current abundance of nucleus \p nucleus in sample \p sample (1 when there is no block).
  double GetAbundance(int sample, int nucleus) const;
  /// Sum of the current abundances of one sample.
  double AbundanceSum(int sample) const;
  /// The sum-to-one width of a sample: its own "sum" line, else \p globalSumSigma.
  double SampleSumSigma(int sample, double globalSumSigma) const;
  /// The abundance penalty added to chi-squared, summed over samples: every
  /// entry with a prior contributes ((w - w0)/sigma)^2, and each sample whose
  /// width s is > 0 contributes ((sum w) - 1)^2 / s^2.
  double AbundancePenalty(double globalSumSigma) const;
  double AbundanceSumPenalty(double globalSumSigma) const;  // sum-to-one terms only (MCMC: priors come from BuildAutoPriors)
  /// Does a pair with this key exist? Keys are what the input files write; the
  /// PPair vector position may differ, so use GetPairNumFromKey to convert.
  bool IsPairKey(int);
  int NumPairs() const;
  int NumJGroups() const;
  /// 1-based position of the pair, or 0 if absent.
  int IsPair(PPair);
  /// 1-based position of the \f$J^\pi\f$ group, or 0 if absent.
  int IsJGroup(JGroup);
  /// Convert an input-file pair key to its 1-based position in the PPair vector.
  int GetPairNumFromKey(int);
  /// Pair \p i, 1-based.
  PPair *GetPair(int);
  const PPair *GetPair(int) const;
  /// The pair with this key, or nullptr when no pair has it (a total-capture exit key -1, a typo): the
  /// GetPair(GetPairNumFromKey(key)) idiom read pairs_[-1] for those (audit 2026-10-25, item 32).
  PPair *GetPairFromKey(int key);
  /// Keys of the capture (gamma, pair type 10) pairs of the compound nucleus of the pair with key
  /// \p entranceKey, in pair order: what the exit key -1 (total capture) stands for.  Empty when the
  /// entrance key is unknown or the nucleus has no capture pair.
  std::vector<int> CapturePairKeys(int entranceKey);
  /// Why a reaction rate cannot be calculated for these keys ("" when it can).  The exit key -1 is the
  /// total radiative capture: the sum over CapturePairKeys(entranceKey).
  std::string RatePairsProblem(int entranceKey, int exitKey);
  /// \f$J^\pi\f$ group \p i, 1-based.
  JGroup *GetJGroup(int);
  const JGroup *GetJGroup(int) const;
  void AddPair(PPair);
  void AddJGroup(JGroup);

  // -- building -------------------------------------------------------------
  /// Build the whole object from the nuclear and external-capture input files.
  /// Returns -1 if a file cannot be read. \p radii overrides channel radii, pair number (1-based) -> fm, as the
  /// API's SetRadius does; it holds every override of the session, so a rebuild keeps them all.
  int Fill(const Config &, const std::map<int, double> &radii = std::map<int, double>());
  /// Adds a level-less J^pi group, with its particle channels, for every J^pi the particle pairs allow up to
  /// \p maxL that no active level provides. Such a group contributes hard-sphere scattering only.
  int AddHardSphereGroups(const Config &, int maxL);
  /// Read the external-capture file into the ECLevel vector and check each
  /// final state against the nuclear file.
  int ParseExternalCapture(const Config &, std::map<int, int> &, const std::vector<azure2::io::LevelRow> &rows);
  /// Boundary conditions, the physical-to-formal transformation, the reaction
  /// pathways and the angular-distribution coefficients -- everything derived
  /// that a calculation needs. Run after Fill, and again after anything that
  /// changes the Coulomb functions.
  void Initialize(const Config &);

  // -- derived quantities ---------------------------------------------------
  /// Boundary condition per channel, evaluated at the energy of the group's
  /// first level.
  void CalcBoundaryConditions(const Config &);
  /// Build the pathway table (Pathways()); CNucPathways.cpp.
  void BuildPathways(const Config &);
  /// Recompute the shift functions at the current level energies. Needed every
  /// iteration under the Brune parameterization, where they move with the fit.
  /// In the Park formalism also dS/dE and J, and, unless `secondDerivatives`
  /// is false (a chi-squared alone needs none), d^2S/dE^2 for the gradient.
  void CalcShiftFunctions(const Config &, bool secondDerivatives = true);
  /// Forget every level's cached shift functions (after a change of the
  /// nuclear potential, which moves them at fixed energy).
  void InvalidateShiftFunctions();
  /// Rescale the level amplitudes of a Minuit parameter set between Brune's
  /// and Park's normalizations (gamma_Park = gamma_Brune sqrt J), in place.
  /// `from`/`to` are AZUREParams::Basis values; false if the conversion is
  /// not level-wise (either side standard) or a level has J <= 0.
  bool ConvertAmplitudeBasis(ROOT::Minuit2::MnUserParameters &, int from, int to, const Config &);
  /// Park formalism: the penalty added to chi-squared for levels whose
  /// overlap J = 1 - sum gamma^2 dS/dE is not positive, sum (J/kParkNormScale)^2
  /// over J < 0.  Zero for every valid parameter set, and in Brune mode.
  double ParkNormPenalty();
  /// Scale of the Park J > 0 penalty: a level at J = -kParkNormScale costs one
  /// unit of chi-squared, so the wall is steep on the scale of J but smooth.
  static constexpr double kParkNormScale = 1.0e-3;
  /// External reduced width amplitude for one channel.

  // -- parameter transformations --------------------------------------------
  /// Observed widths and energies to formal R-matrix parameters. Returns false
  /// if a level could not be transformed.
  bool TransformIn(const Config &);
  /// Why the last TransformIn returned false (empty if it did not).
  const std::string &TransformError() const { return transformError_; }
  /// Has TransformIn run on the current level scheme?  Until it has, the
  /// level gammas are the input values (partial widths in eV, ANCs, gamma
  /// widths in eV), not reduced width amplitudes -- consumers that need to
  /// estimate a resonance width before initialisation completes must check.
  bool IsTransformedIn() const { return transformedIn_; };
  /// Formal basis with external capture to a final state that shares its J^pi block with other levels: the final
  /// state's capture amplitude is its Brune amplitude, recomputed with the fit parameters (cnuc_split.md section 8).
  bool HasFormalECFinalBlocks() const { return ecFinalBlocks_ != nullptr; }
  /// The --ignore-externals mask of a channel capture in the formal basis (external_capture.md section 3): the
  /// projector Pi = T^-1 P T over the R-matrix levels `levels` of group `group`, for the photon channel `channel`,
  /// at the current fit parameters.  The masked entrance and exit vectors are Pi gamma.
  struct CaptureMask {
    int group = 0, channel = 0;
    std::vector<int> levels;
    std::vector<double> projector;  ///< row-major, levels x levels; empty: Pi = P (no masked level, or all decoupled)
  };
  /// The formal-basis mask of (group, photon channel), or nullptr: the mask is then diagonal, a level whose input
  /// width in the channel is zero (ALevel::ZeroInputWidth) is left out (Brune, Park, a block of one level).
  const CaptureMask *FormalCaptureMask(int group, int channel) const;
  /// Formal basis, --ignore-externals, channel capture through a block of several levels: the masks are recomputed
  /// with the fit parameters.  True when the adjoint does not differentiate that dependence: a slot has a projector,
  /// or a masked level has a free particle amplitude (freeGamma(group, level, channel), 1-based).
  bool FormalCaptureMasksOutsideAdjoint(const std::function<bool(int, int, int)> &freeGamma);
  /// Formal parameters back to observed energies and partial widths, for
  /// reporting. Under Brune the two are the same and this is a near no-op.
  void TransformOut(const Config &);
  /// The transformed physical parameters as a flat vector.
  vector_r GetTransformParams(const Config &configure);
  /// The observed widths of a parameter vector in the levels' observed slots, as TransformOut computes them but
  /// without its messages (for WarnLargeRadiativeWidths, CNucChecks.h); the fit state is set back to the vector.
  /// False, with nothing set, if the transformation fails.
  bool ObservedWidthsAt(const Config &, const vector_r &);
  /// The observed quantities of the current fit values (TransformOut without its messages) in the levels' observed and
  /// transform slots, as the model data read them (ModelData.cpp) inside an evaluation: nothing a calculation reads
  /// is changed (the bound-state normalizations TransformToBrune sets are put back).  False if the transformation
  /// fails (the slots are then not meaningful).  CNucObserved.cpp.
  bool ObservedOfFitValues(const Config &);
  /// Set one fit value, as FillCompoundFromParams sets it with everything derived from it: the energy of level la of
  /// J-group j (channel 0) or its amplitude in channel ch (1-based).  For the numerical derivatives of the model data;
  /// setting the value back restores the state bit for bit.  CNucObserved.cpp.
  void SetFitLevelParameter(int j, int la, int ch, double value);

  // -- fitting --------------------------------------------------------------
  /// Seed the parameter vector from the current level parameters, each with its record (what it is).
  void FillMnParams(AZUREParams &, const Config *config = nullptr);
  /// Write a Minuit parameter vector back into the levels.
  void FillCompoundFromParams(const vector_r &);
  /// As above, for a vector already in physical rather than formal parameters.
  void FillCompoundFromParamsPhysical(const vector_r &);

  /// The reaction pathways of the run (PathwayTable.h), built by Initialize and shared by every copy (Clone).
  /// Throws std::logic_error before Initialize.
  const azure2::PathwayTable &Pathways() const;
  /// Ask Initialize to build the capture analyzing-power terms of (entrance pair, exit pair), both 1-based, in the
  /// pathway table.  Set-up state, like PPair::SetEntrance: it must be called before Initialize.
  void RequestCaptureAy(int entrancePair, int exitPair);

  /// Deep copy. AZURECalc gives each thread its own, since a calculation
  /// writes the fitted parameters into the object it runs on.
  CNuc *Clone() const;

  /// Current values of the THM coherent backgrounds (`cbackground=`,
  /// ThmCoherentBackground in ThmExperiment.h), in parameter order: the data
  /// own the parameters (WorkingCopy::FillThmCoherentFromParams), the HOES model
  /// (THMMatrixFunc) reads them here, on the compound it is handed.  Empty
  /// without the key.
  const std::vector<double> &ThmCoherentValues() const { return thmCoherentValues_; }
  void SetThmCoherentValues(const std::vector<double> &v) { thmCoherentValues_ = v; }

 private:
  struct ECFinalBlocks;
  /// Formal basis with external capture: the blocks whose final states need their Brune amplitude (Initialize).
  void SetUpECFinalBlocks(const Config &);
  /// Their final states' external-capture amplitudes and N_f^1/2 from the current fit values (FillCompoundFromParams).
  void ECFinalStatesFromFormal();
  std::shared_ptr<const ECFinalBlocks> ecFinalBlocks_;
  struct CaptureMaskSetup;
  /// Formal basis with --ignore-externals: the (block, photon channel) pairs of channel capture (Initialize, after
  /// the pathways are sorted).
  void SetUpFormalCaptureMasks(const Config &);
  /// Their projectors from the current fit values (FillCompoundFromParams).
  void FormalCaptureMasksFromParams();
  bool MaskedLevelCoupled(JGroup *, int) const;
  std::shared_ptr<const CaptureMaskSetup> captureMaskSetup_;
  std::vector<CaptureMask> captureMasks_;
  /// TransformOut, first part: the Brune (Park) E and gamma of the fit parameters in the levels' transform slots
  /// (formal basis: physics/ brune_from_formal); report prints the warnings of levels that do not convert.
  void TransformToBrune(const Config &, bool report);
  /// TransformOut, second part: the observed quantities from the transform slots (physics/ observed_levels).
  void ObservedWidths(const Config &);
  std::string transformError_;
  /// 1-based index of the nucleus with this (Z, A), appending it if new.
  int FindOrAddNucleus(int z, int a);
  /// TEMPORARY (step 3d, src/ProjectSetup.cpp): the old <abundances> section through SetAbundances.
  int ReadAbundances(const Config &);
  /// The abundances of the samples from their spec into samples_; -1 (with the reader's message) when the reader
  /// stopped.
  int SetAbundances(const Config &, const azure2::experiment::AbundanceSpecs &);
  /// The pathway table (Initialize) and the capture-A_y requests it is built with.
  std::shared_ptr<const azure2::PathwayTable> pathways_;
  std::vector<std::pair<int, int> > captureAyRequests_;
  std::vector<PPair> pairs_;
  std::vector<JGroup> jgroups_;
  std::vector<CompoundNucleus> nuclei_;
  std::vector<AbundanceSample> samples_;
  bool hasAbundanceBlock_ = false;
  /// Setup state, shared by the copies (Clone).
  std::shared_ptr<const ModelDataSet> modelData_ = std::make_shared<ModelDataSet>();
  bool transformedIn_ = false;
  std::vector<double> thmCoherentValues_;
};

extern double DoubleFactorial(int);

#endif
