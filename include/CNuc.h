#ifndef CNUC_H
#define CNUC_H

#include <string>
#include <functional>
#include <exception>
#include <memory>
#include <vector>
#include <map>
#include <optional>
#include "JGroup.h"
#include "PPair.h"
#include "PathwayTable.h"
#include "ModelData.h"
#include "GammaDecayMap.h"
#include "ObservedLevelTable.h"
#include "LevelInputs.h"
#include "LevelState.h"
#include "LevelParameterLayout.h"
#include "azure2/experiment/SetupSpecs.h"

namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT
class Config;
class AZUREParams;
struct CNucSystem;
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
 *
 * It holds the set-up only: the structure, the pathways and the set-up records, written by Fill, AddHardSphereGroups,
 * ParseExternalCapture and Initialize and read-only after them.  The values of an evaluation are the caller's: a
 * LevelState (and the LevelInputs of a transformation) per owner, passed to every fill, transformation and reader.
 * One object is shared, constant, by every evaluation, thread and walker of a run or a session; it is not copied.
 */

/// One compound nucleus of a project: the (Z, A) its pairs form.
struct CompoundNucleus {
  int Z = 0;
  int A = 0;
};

class CNuc {
 public:
  CNuc() = default;
  CNuc(const CNuc &) = delete;
  CNuc &operator=(const CNuc &) = delete;
  // -- structure ------------------------------------------------------------
  /// Number of compound nuclei the pairs form.  One for every project written
  /// before several nuclei could share a file.
  int NumNuclei() const;
  /// Nucleus \p i, 1-based.
  const CompoundNucleus *GetNucleus(int) const;
  /// "Z=10 A=23" style label of nucleus \p i.
  std::string NucleusLabel(int) const;
  /// The model data (ModelData.h): the rows of the objective's term T9, read by Fill.
  const ModelDataSet &GetModelData() const { return *modelData_; }
  /// The gamma_decay records (GammaDecayMap.h), read by Fill; empty without <gammaDecays>.
  const GammaDecaySet &GetGammaDecays() const { return *gammaDecays_; }
  /// Does a pair with this key exist? Keys are what the input files write; the
  /// PPair vector position may differ, so use GetPairNumFromKey to convert.
  bool IsPairKey(int) const;
  int NumPairs() const;
  int NumJGroups() const;
  /// 1-based position of the \f$J^\pi\f$ group, or 0 if absent.
  int IsJGroup(JGroup) const;
  /// Convert an input-file pair key to its 1-based position in the PPair vector.
  int GetPairNumFromKey(int) const;
  /// Pair \p i, 1-based.  The structure (pairs, J^pi groups, channels) is read-only after Fill and Initialize: only
  /// CNuc's own set-up writes it (Pair, Group, below).
  const PPair *GetPair(int) const;
  /// The pair with this key, or nullptr when no pair has it (a total-capture exit key -1, a typo): the
  /// GetPair(GetPairNumFromKey(key)) idiom read pairs_[-1] for those (audit 2026-10-25, item 32).
  const PPair *GetPairFromKey(int key) const;
  /// Keys of the capture (gamma, pair type 10) pairs of the compound nucleus of the pair with key
  /// \p entranceKey, in pair order: what the exit key -1 (total capture) stands for.  Empty when the
  /// entrance key is unknown or the nucleus has no capture pair.
  std::vector<int> CapturePairKeys(int entranceKey) const;
  /// Why a reaction rate cannot be calculated for these keys ("" when it can).  The exit key -1 is the
  /// total radiative capture: the sum over CapturePairKeys(entranceKey).
  std::string RatePairsProblem(int entranceKey, int exitKey) const;
  /// \f$J^\pi\f$ group \p i, 1-based (read-only, as GetPair).
  const JGroup *GetJGroup(int) const;
  /// Where each level's energy and amplitudes sit in the parameter vector (LevelParameterLayout.h); built by Fill.
  /// Every walk of the level block reads it.
  const LevelParameterLayout &ParameterLayout() const { return *parameterLayout_; }
  /// Take the fit energies and amplitudes of `from` into `to` (states of this structure), nothing derived from them
  /// (WorkingCopy::ReplaceSubPointGrids' state).
  void CopyFitEnergiesAndAmplitudes(const LevelState &from, LevelState &to) const;

  // -- building -------------------------------------------------------------
  /// Build the whole object from the nuclear and external-capture input files.
  /// Returns -1 if a file cannot be read. \p radii overrides channel radii, pair number (1-based) -> fm, as the
  /// API's SetRadius does; it holds every override of the session, so a rebuild keeps them all.  `state` and `inputs`
  /// are the caller's (the run's, the session's): shaped for the levels, the inputs made from the <levels> rows and
  /// the gamma_decay records, the state given the input energies and zero-width flags.  This object keeps neither.
  int Fill(const Config &, LevelState &state, LevelInputs &inputs,
           const std::map<int, double> &radii = std::map<int, double>());
  /// Adds a level-less J^pi group, with its particle channels, for every J^pi the particle pairs allow up to
  /// \p maxL that no active level provides. Such a group contributes hard-sphere scattering only.
  int AddHardSphereGroups(const Config &, int maxL);
  /// Read the external-capture file into the ECLevel vector and check each
  /// final state against the nuclear file.
  int ParseExternalCapture(const Config &, std::map<int, int> &, const std::vector<azure2::io::LevelRow> &rows);
  /// Boundary conditions, the physical-to-formal transformation, the reaction
  /// pathways and the angular-distribution coefficients -- everything derived
  /// that a calculation needs. Run after Fill, and again after anything that
  /// changes the Coulomb functions.  The transformation writes the caller's `state` and `inputs` (those Fill made).
  void Initialize(const Config &, LevelState &state, LevelInputs &inputs);

  // -- derived quantities ---------------------------------------------------
  /// Boundary condition per channel, evaluated at the energy of the group's
  /// first level.
  void CalcBoundaryConditions(const Config &, const LevelInputs &inputs);
  /// Build the pathway table (Pathways()); CNucPathways.cpp.
  void BuildPathways(const Config &);
  /// Which energy derivatives of the shift functions CalcShiftFunctions takes in Park's basis: a chi-squared alone
  /// needs no second derivatives, the gradient does.
  enum class ShiftOrders { kWithSecondDerivatives, kFirstDerivativesOnly };
  /// Recompute the shift functions of `state` at its level energies. Needed every
  /// iteration under the Brune parameterization, where they move with the fit.
  /// In the Park formalism also dS/dE and J, and d^2S/dE^2 as `orders` asks.
  void CalcShiftFunctions(const Config &, LevelState &state,
                          ShiftOrders orders = ShiftOrders::kWithSecondDerivatives) const;
  /// Rescale the level amplitudes of a Minuit parameter set between Brune's
  /// and Park's normalizations (gamma_Park = gamma_Brune sqrt J), in place.
  /// `from`/`to` are AZUREParams::Basis values; false if the conversion is
  /// not level-wise (either side standard) or a level has J <= 0.
  bool ConvertAmplitudeBasis(ROOT::Minuit2::MnUserParameters &, int from, int to, const Config &) const;
  /// External reduced width amplitude for one channel.

  // -- parameter transformations --------------------------------------------
  /// Observed widths and energies of `inputs` to the fit's basis (Brune, Park or formal), in `inputs`, with the
  /// levels' N_f^1/2 and Park's J in `state`.  Returns false if a level could not be transformed;
  /// inputs.TransformError() says why.
  bool TransformIn(const Config &, LevelState &state, LevelInputs &inputs) const;
  /// Formal basis with external capture to a final state that shares its J^pi block with other levels: the final
  /// state's capture amplitude is its Brune amplitude, recomputed with the fit parameters (cnuc_split.md section 8).
  bool HasFormalECFinalBlocks() const { return ecFinalBlocks_ != nullptr; }
  /// The --ignore-externals mask of a channel capture in the formal basis (external_capture.md section 3): a slot per
  /// (group, photon channel) with the R-matrix levels of the block the projector Pi = T^-1 P T acts on.  The projector
  /// at the current fit parameters is the state's (LevelState::CaptureProjector(slot)); the masked entrance and exit
  /// vectors are Pi gamma.  Set up by Initialize.
  struct CaptureMaskSlot {
    std::size_t block = 0;  ///< the block of the structure System
    int group = 0, channel = 0;
    std::vector<int> levels;
  };
  /// The formal-basis mask slot of (group, photon channel), or -1: the mask is then diagonal, a level whose input
  /// width in the channel is zero (LevelState::ZeroInputWidth) is left out (Brune, Park, a block of one level).
  int FormalCaptureMaskSlot(int group, int channel) const;
  const CaptureMaskSlot &FormalCaptureMaskOf(int slot) const;
  /// Formal basis, --ignore-externals, channel capture through a block of several levels: the masks are recomputed
  /// with the fit parameters.  True when the adjoint does not differentiate that dependence: a slot has a projector,
  /// or a masked level has a free particle amplitude (freeGamma(group, level, channel), 1-based).
  bool FormalCaptureMasksOutsideAdjoint(const LevelState &state,
                                        const std::function<bool(int, int, int)> &freeGamma) const;
  /// The Brune (Park) and observed values of the levels at the fit values of `state`, for the reports; prints the
  /// warnings of levels that do not convert.  Changes nothing of `state` but Park's shift functions (TransformToBrune).
  ObservedLevelTable TransformOut(const Config &, LevelState &state, const LevelInputs &inputs) const;
  /// The observed parameters of `table` as a flat vector, as a <levels> line holds them (a gamma_decay record's
  /// numbers at its slots); the level block only.
  vector_r GetTransformParams(const LevelState &state, const ObservedLevelTable &table) const;
  /// The observed values of a parameter vector, as TransformOut computes them but without its messages (for
  /// WarnLargeRadiativeWidths, CNucChecks.h); `state` is left at the vector.  Empty if the transformation fails.
  std::optional<ObservedLevelTable> ObservedWidthsAt(const Config &, const vector_r &, LevelState &state,
                                                     const LevelInputs &inputs) const;
  /// The observed values of the fit values of `state` (TransformOut without its messages), as the model data read
  /// them (ModelData.cpp) inside an evaluation: nothing a calculation reads is changed.  Empty if the transformation
  /// fails.  CNucObserved.cpp.
  std::optional<ObservedLevelTable> ObservedOfFitValues(const Config &, LevelState &state,
                                                        const LevelInputs &inputs) const;
  /// Set the fit value of one slot of the level block (ParameterLayout), as FillCompoundFromParams sets it with
  /// everything derived from it: a level's energy or its amplitude in one channel.  For the numerical derivatives of
  /// the model data; setting the value back restores the state bit for bit.  CNucObserved.cpp.
  void SetFitLevelParameter(LevelState &state, const LevelParameterLayout::Slot &slot, double value) const;

  // -- fitting --------------------------------------------------------------
  /// Seed the parameter vector from the levels' input values `inputs`, each with its record (what it is); a width
  /// whose input was zero (state.ZeroInputWidth) is fixed.
  void FillMnParams(AZUREParams &, const LevelState &state, const LevelInputs &inputs,
                    const Config *config = nullptr) const;
  /// Write a parameter vector into `state`, with everything derived from it but the shift functions.
  void FillCompoundFromParams(const vector_r &, LevelState &state) const;
  /// The vector of observed parameters (GetTransformParams' form: a gamma_decay record's numbers at its slots) as the
  /// levels' input values in `inputs`, with N_f^1/2 in `state`; a record's numbers become its input numbers and its
  /// partial widths the level's input widths (GammaDecaySet::WithPartialWidths; gamma_parametrisation.md section 2.8).
  void FillCompoundFromParamsPhysical(const vector_r &, LevelState &state, LevelInputs &inputs) const;

  /// The physics/ System of the structure (CNucSystem.h, CNucSystemContent::Structure: R-matrix blocks, B_c = 0, no
  /// level values) and the same with every J-group and level (CNucSystemContent::AllGroups).  Built once, when Fill
  /// has made the structure; the structure does not change after Fill but for B_c,
  /// which neither holds.  Throws, as BuildSystemFromCNuc did, when the System cannot represent the structure, and
  /// std::logic_error before Fill.
  const CNucSystem &SystemStructure() const;
  const CNucSystem &SystemAllGroups() const;
  /// The handle of SystemStructure(), for a set-up object that keeps it (GammaDecaySet, the external-capture set-up).
  std::shared_ptr<const CNucSystem> SystemStructureHandle() const;
  /// The reaction pathways of the run (PathwayTable.h), built by Initialize.
  /// Throws std::logic_error before Initialize.
  const azure2::PathwayTable &Pathways() const;
  /// Ask Initialize to build the capture analyzing-power terms of (entrance pair, exit pair), both 1-based, in the
  /// pathway table.  Set-up state, like RequestEntrance: it must be called before Initialize.
  void RequestCaptureAy(int entrancePair, int exitPair);
  /// Ask Initialize to build the pathways of entrance pair `pair` (1-based).  Set-up state: before Initialize.
  void RequestEntrance(int pair);

 private:
  /// The set-up's writable views of the structure (Fill, AddHardSphereGroups, ParseExternalCapture, Initialize).
  PPair *Pair(int);
  JGroup *Group(int);
  /// 1-based position of the pair, or 0 if absent.
  int IsPair(PPair);
  void AddPair(PPair);
  void AddJGroup(JGroup);
  struct ECFinalBlocks;
  /// Formal basis with external capture: the blocks whose final states need their Brune amplitude (Initialize).
  void SetUpECFinalBlocks(const Config &, const LevelInputs &inputs);
  /// Their final states' external-capture amplitudes and N_f^1/2 from the current fit values (FillCompoundFromParams).
  void ECFinalStatesFromFormal(LevelState &state) const;
  std::shared_ptr<const ECFinalBlocks> ecFinalBlocks_;
  struct CaptureMaskSetup;
  /// Formal basis with --ignore-externals: the (block, photon channel) pairs of channel capture (Initialize, after
  /// the pathways are sorted).
  void SetUpFormalCaptureMasks(const Config &, LevelState &state);
  /// Their projectors from the current fit values (FillCompoundFromParams).
  void FormalCaptureMasksFromParams(LevelState &state) const;
  bool MaskedLevelCoupled(const LevelState &state, int j, int la) const;
  std::shared_ptr<const CaptureMaskSetup> captureMaskSetup_;
  /// Whether TransformToBrune prints the warnings of levels that do not convert.
  enum class Messages { kPrint, kQuiet };
  /// A table with an entry for every level and channel, all values zero.
  ObservedLevelTable EmptyObservedTable() const;
  /// TransformOut, first part: the Brune (Park) E and gamma of the fit parameters, with N_f^-1/2, in `table` (formal
  /// basis: physics/ brune_from_formal).
  void TransformToBrune(const Config &, Messages, LevelState &state, const LevelInputs &inputs,
                        ObservedLevelTable &table) const;
  /// TransformOut, second part: the observed quantities from the Brune (Park) values of `table` (physics/
  /// observed_levels).
  void ObservedWidths(const Config &, const LevelState &state, ObservedLevelTable &table) const;
  /// 1-based index of the nucleus with this (Z, A), appending it if new.
  int FindOrAddNucleus(int z, int a);
  /// The pathway table (Initialize) and the capture-A_y requests it is built with.
  std::shared_ptr<const azure2::PathwayTable> pathways_;
  std::vector<std::pair<int, int> > captureAyRequests_;
  std::vector<int> entranceRequests_;  ///< RequestEntrance, 1-based pair numbers
  /// Fill's last step: SystemStructure and SystemAllGroups, or the exception their build threw.
  void BuildSystems(const Config &);
  std::shared_ptr<const CNucSystem> systemStructure_, systemAllGroups_;
  std::exception_ptr systemError_;
  std::vector<PPair> pairs_;
  std::vector<JGroup> jgroups_;
  std::vector<CompoundNucleus> nuclei_;
  /// Set-up state.
  std::shared_ptr<const ModelDataSet> modelData_ = std::make_shared<ModelDataSet>();
  std::shared_ptr<const GammaDecaySet> gammaDecays_ = std::make_shared<GammaDecaySet>();
  /// The level block of the parameter vector (ParameterLayout); built by Fill.
  std::shared_ptr<const LevelParameterLayout> parameterLayout_ = std::make_shared<LevelParameterLayout>();
};

extern double DoubleFactorial(int);

#endif
