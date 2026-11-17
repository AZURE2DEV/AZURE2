// AZUREAPI, the session's model: made and released, set up from the project (Initialize, SetUpModel), set up
// again when a radius, a potential or a mode changes (SetRadius, Rebuild, RebuildImpl, ReplaceModel), its prior
// covariance, and its parameter vectors with their bounds and index lists (UpdateParameters, SetParamsAllRWA,
// TransformAllRWAParameters, ReadFitSettings, GetParameterBounds, Get*Indices).  One unit because these are what
// replaces or re-reads compound_, data_, the level state and the vectors every other call takes as given.  The
// other jobs of the class are files of their own: AZUREAPIReports.cpp (what a calculation reports),
// AZUREAPIChi2.cpp (the objective's value), AZUREAPIGradient.cpp (its derivatives), AZUREAPIObjective.cpp (the
// objective with its breakdown, as rows, as a log-probability), AZUREAPIFit.cpp (the program's fitters),
// AZUREAPIModelData.cpp, AZUREAPIChannelFunctions.cpp, AZUREAPIDiagnostics.cpp.
#include <algorithm>
#include "Objective.h"
#include "ParameterLimitsManager.h"
#include "AZUREAPI.h"
#include "CompoundSystemSetUp.h"
#include "AdaptiveExtrapolation.h"
#include "ChannelFunc.h"
#include "ModelPass.h"
#include "ScatteringLength.h"
#include "AZUREParams.h"

#include "GSLException.h"
#include <stdexcept>
#include "ECAmplitudeCache.h"

#include "Config.h"
#include "CompoundSystem.h"
#include "LevelReport.h"
#include "CaptureAyRequests.h"
#include "ProjectSetup.h"
#include "PPair.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "Unit.h"
#include "TargetEffect.h"
#include "AMatrixFunc.h"
#include "AZUREGrad.h"
#include "CovarianceBand.h"
#include "PriorCovarianceFile.h"
#include "AZURECalc.h"
#include "CoulFunc.h"
#include "ECIntegral.h"
#include "Profile.h"
#include "InputFiles.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <new>
#include <cmath>

AZUREAPI::AZUREAPI(Config &configure) :
  configure_(configure), data_(nullptr), state_(std::make_unique<LevelState>()),
  inputs_(std::make_unique<LevelInputs>()) {
  // A session rebuilds its model from the files (SetRadius, potential changes, mode switches): read each
  // input file once and print each setup Note once, not on every rebuild (Config::cacheInputFiles, Note).
  configure_.cacheInputFiles = true;
  if (!configure_.noteLog) configure_.noteLog = std::make_shared<Config::NoteLog>();
}

AZUREAPI::~AZUREAPI() {
  // A session used to be a process, so the OS reclaimed these; in-process the
  // compound nucleus and its data are ~10 MB per model and must be returned.
  compound_.reset();
  if (data_ != nullptr) delete data_;
}

int AZUREAPI::Initialize() {
  // Initialize caches for performance
  InitializeECAmplitudeCache(configure());

  configure_.paramMask |= Config::USE_EXTERNAL_CAPTURE;

  std::string file;
  if (configure().paramMask & Config::CALCULATE_WITH_DATA)
    file = configure().outputdir + "intEC.dat";
  else
    file = configure().outputdir + "intEC.extrap";

  std::ifstream in(configure().OpenPath(file).c_str());
  if (!in)
    configure_.paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  else
    configure_.paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  configure_.integralsfile = file;

  // Initialize EC Integral caching system

  std::string cacheFile;
  if (configure().paramMask & Config::CALCULATE_WITH_DATA) {
    cacheFile = configure().outputdir + "intEC_cache.dat";
  } else {
    cacheFile = configure().outputdir + "intEC_cache.extrap";
  }

  // Initialize() is called again on every mode switch: a new model, as RebuildImpl makes, with the session's radius
  // overrides.
  if (SetUpModel(ModelSetUpOccasion::kFirst, radiusOverrides_) != CompoundSetUp::Refusal::kNone) return -1;

  UpdateParameters();
  if (!SetUpPriorCovariance()) return -1;
  if (data()->Initialize(compound(), configure()) == -1) return -1;
  AZUREProfile::Report(configure().profile, std::cerr, "api set-up");
  // The session's buffer for the new layout, at 0 until an evaluation writes it (a read before any evaluation, such
  // as GetThmExperiments, sees 0 as it did when the data set held the buffer).
  results_ = azure2::Results();
  data()->SizeResults(results_);

  return 0;
}

bool AZUREAPI::SetUpPriorCovariance() {
  std::shared_ptr<const PriorCovariance> prior;
  if (data() && (configure().paramMask & Config::CALCULATE_WITH_DATA) &&
      !(configure().paramMask & Config::CALCULATE_REACTION_RATE) &&
      LoadPriorCovariance(configure(), names_, fixed_, AZUREParams::BasisForConfig(configure()), prior) < 0)
    return false;
  if (data()) data()->SetPriorCovariance(prior);
  // The prior covariance replaces the diagonal prior of each of its parameters (docs/upgrade/objective.md D3).
  if (prior && prior->Size() > 0) {
    const std::vector<std::string> replaced = ReplaceDiagonalPriors(compound(), data(), *prior);
    std::string list;
    for (const std::string &name : replaced) list += (list.empty() ? "" : ", ") + name;
    if (!list.empty())
      configure().Note("Note: prior covariance: its multivariate Gaussian replaces the prior of " + list +
                       " of its own (normalization error, shift error or sigma); a nuisance prior on a listed "
                       "parameter is not applied either.");
  }
  return true;
}

vector_r AZUREAPI::PriorCovarianceResidualsRWA(const vector_r &rwaParams) const {
  vector_r r;
  const PriorCovariance *pc = data() ? data()->Data().GetPriorCovariance() : nullptr;
  if (!pc) return r;
  vector_r full = all_rwa_;
  for (int i = 0, k = 0; i < (int)all_rwa_.size(); ++i)
    if (!fixed_[i] && k < (int)rwaParams.size()) full[i] = rwaParams[k++];
  pc->Residuals(full, r);
  return r;
}

std::vector<std::string> AZUREAPI::PriorCovarianceNames() const {
  const PriorCovariance *pc = data() ? data()->Data().GetPriorCovariance() : nullptr;
  return pc ? pc->Names() : std::vector<std::string>();
}

bool AZUREAPI::UpdateParameters() {
  all_.clear();
  all_rwa_.clear();
  names_.clear();
  records_.clear();
  fixed_.clear();
  values_.clear();
  values_rwa_.clear();

  AZUREParams params;
  compound()->FillMnParams(params, state(), inputs(), &configure());
  data()->FillMnParams(params);

  compound()->FillCompoundFromParams(params.GetMinuitParams().Params(), state());

  compound()->CalcShiftFunctions(configure(), state());
  const ObservedLevelTable observed = compound()->TransformOut(configure(), state(), inputs());

  records_ = params.Records();
  for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
    // if( !params.GetMinuitParams().Parameter(i).IsFixed( ) ){
    names_.push_back(params.GetMinuitParams().Parameter(i).GetName());
    //}
    all_rwa_.push_back(params.GetMinuitParams().Parameter(i).Value());
    fixed_.push_back(params.GetMinuitParams().Parameter(i).IsFixed());
    if (!fixed_.back()) {
      values_rwa_.push_back(params.GetMinuitParams().Parameter(i).Value());
    }
  }

  all_ = compound()->GetTransformParams(state(), observed);
  data()->Data().Abundances().AppendTo(all_rwa_, all_);

  // Norms and shifts are missing from all_ and need to be added using the all_rwa_ values
  for (int i = all_.size(); i < all_rwa_.size(); ++i) {
    all_.push_back(all_rwa_[i]);
  }

  for (int i = 0; i < all_.size(); ++i) {
    // std::cout << "all_[" << i << "] = " << all_[i] << std::endl;
    // std::cout << "fixed_[" << i << "] = " << fixed_[i] << std::endl;
    // std::cout << "names_[" << i << "] = " << names_[i] << std::endl;
    if (!fixed_[i]) {
      values_.push_back(all_[i]);
    }
  }

  return true;
}

/*!
 * Replaces every parameter value, fixed ones included, with p (RWA basis, the
 * params_all_rwa() layout) and refreshes the physical copies, without
 * calculating anything.
 *
 * The free-vector entry points (CalculateChi2RWA, UpdateSegmentsRWA,
 * TransformAllRWAParameters, ...) take the FIXED parameters from all_rwa_, i.e.
 * from the .azr.  The CLI's external parameter file sets fixed level
 * parameters too (AZUREParams::ReadUserParameters), so a session given only the
 * free part of a param.sav computes a different model whenever the file's
 * fixed values differ from the .azr's: 11B+alpha 10-12-26, 10 fixed widths,
 * chi2 9925.03 against the CLI's 9868.04.  This is how a session takes them.
 */
bool AZUREAPI::SetParamsAllRWA(const vector_r &p) {
  if (p.size() != all_rwa_.size()) return false;
  all_rwa_ = p;

  // Physical values: the R-matrix block through the same transform as
  // UpdateParameters; norms and shifts are the same in both bases.
  compound()->FillCompoundFromParams(all_rwa_, state());
  data()->FillAbundancesFromVector(all_rwa_);
  compound()->CalcShiftFunctions(configure(), state());
  all_ = compound()->GetTransformParams(state(), compound()->TransformOut(configure(), state(), inputs()));
  data()->Data().Abundances().AppendTo(all_rwa_, all_);
  for (size_t i = all_.size(); i < all_rwa_.size(); ++i) all_.push_back(all_rwa_[i]);

  values_.clear();
  values_rwa_.clear();
  for (size_t i = 0; i < all_rwa_.size(); ++i) {
    if (!fixed_[i]) {
      values_.push_back(all_[i]);
      values_rwa_.push_back(all_rwa_[i]);
    }
  }
  return true;
}

// Transform RWA parameters to physical values
vector_r AZUREAPI::TransformAllRWAParameters(const vector_r &p, bool includeFixed) {
  // p holds every parameter, fixed ones included; a short p updates a prefix.
  vector_r params = all_rwa_;
  for (int i = 0; i < p.size() && i < params.size(); ++i) {
    params[i] = p[i];
  }

  const CompoundSystem *localCompound = NULL;
  WorkingCopy *localData = NULL;
  localCompound = compound();
  localData = data();

  localCompound->FillCompoundFromParams(params, state());
  localData->FillAbundancesFromVector(params);

  vector_r transformedParams = localCompound->GetTransformParams(state(), localCompound->TransformOut(configure(), state(), inputs()));
  data()->Data().Abundances().AppendTo(params, transformedParams);

  if (includeFixed) return transformedParams;
  // Get only non fixed parameters
  vector_r transformed;
  for (int i = 0; i < transformedParams.size(); ++i) {
    if (!fixed_[i]) {
      transformed.push_back(transformedParams[i]);
    }
  }

  return transformed;
}

bool AZUREAPI::CalculateExternalCapture() {
  configure_.paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  data()->CalculateECAmplitudes(compound(), configure());
  configure_.paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  return true;
}

// Set AZURE2 to calculate data points
void AZUREAPI::SetData() {
  configure_.paramMask |= Config::CALCULATE_WITH_DATA;
}

// Set AZURE2 to calculate extrapolations
void AZUREAPI::SetExtrap() {
  configure_.paramMask &= ~Config::CALCULATE_WITH_DATA;
}


// Set the channel radius of one particle pair (1-based) and rebuild.
//
// The radius enters the penetrabilities, shift functions, boundary conditions,
// Wigner limits, the ANC <-> reduced-width conversion and the lower limit of
// every external-capture integral, so nothing can be reused: the compound
// nucleus and data objects are rebuilt from the .azr with the override, the
// EC-integral cache is bypassed (USE_PREVIOUS_INTEGRALS cleared) so the
// integrals are recomputed on the new radius rather than read back from
// output/intEC*, and the parameter bookkeeping -- values, names, fixed flags
// and the transformed physical parameters -- is refilled to match.
bool AZUREAPI::SetRadius(int idx, double r) {
  // The new override joins the earlier ones; it is kept only if the rebuild worked.
  std::map<int, double> radii = radiusOverrides_;
  radii[idx] = r;
  if (!RebuildImpl(radii)) return false;
  radiusOverrides_ = radii;
  return true;
}

bool AZUREAPI::Rebuild() {
  return RebuildImpl(radiusOverrides_);
}

// Shared body of SetRadius and Rebuild.  Every radius override of the session
// goes to CompoundSystem::Fill, so a rebuild for one reason (a potential change, another
// pair's radius) does not undo an earlier SetRadius.
// The one place compound_ and data_ are replaced (Initialize, RebuildImpl): with them goes everything derived from
// them and built lazily (the nuisance priors' manager and its parameter list), so it is built again for the new model.
void AZUREAPI::ReplaceModel() {
  compound_.reset();
  delete data_;
  data_ = new WorkingCopy();
  *state_ = LevelState();  // in place: a manager made later refers to the session's levels
  *inputs_ = LevelInputs();
  nuisance_.reset();
  nuisanceParams_.reset();
  nuisanceRead_ = false;
}

CompoundSetUp::Refusal AZUREAPI::SetUpModel(ModelSetUpOccasion occasion, const std::map<int, double> &radii) {
  // What a first set-up and a rebuild differ in; everything else is one set-up.
  //  - The project's check files (compound nucleus, data): the first set-up writes them; a rebuild makes the same
  //    project's model again and writes none.
  const bool writesCheckFiles = occasion == ModelSetUpOccasion::kFirst;
  //  - The external-capture integrals: a rebuild computes them again (its radii or potentials change them) instead of
  //    reading intEC* back: USE_PREVIOUS_INTEGRALS is cleared once the data are set up, and RebuildImpl sets it
  //    again after the data's initialization or a refusal of the last part.
  const bool recomputesIntegrals = occasion == ModelSetUpOccasion::kRebuild;
  //  - A refused last part: the first set-up says it in the program's words (its message, a blank line, "Calculation
  //    was aborted."); a rebuild, asked for by a caller who handles its false, prints the message only.
  const bool saysAbort = occasion == ModelSetUpOccasion::kFirst;

  ReplaceModel();
  // The session is silent on a refusal (it sets Config::USE_API; the caller raises), but for the rate's pairs.
  const DataSetUpStep dataSetUp = [this, writesCheckFiles, recomputesIntegrals](
                                      const CompoundSystem &structure,
                                      const LevelInputs &levelInputs) -> std::optional<CompoundSetUpRequests> {
    if (writesCheckFiles && ((configure().screenCheckMask | configure().fileCheckMask) & Config::CHECK_COMPOUND_NUCLEUS))
      WriteCompoundNucleusCheck(structure, levelInputs, configure());
    CompoundSetUpRequests requests;
    if (!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
      // The data object from the segments and the data files: the structure gives the pair keys, the centre-of-mass
      // conversions and the S-factor conversions.
      if (configure().paramMask & Config::CALCULATE_WITH_DATA) {
        if (data()->Fill(configure(), &structure, levelInputs) != 0) {
          return std::nullopt;
        } else if (data()->Data().NumSegments() == 0 && structure.GetModelData().NumResiduals() == 0) {
          return std::nullopt;
        }
        RequestDataEntrances(structure, configure(), requests);
      } else {
        if (data()->MakePoints(configure(), &structure, levelInputs) != 0 || SessionRefusesAdaptiveGrids(configure())) {
          return std::nullopt;
        } else if (data()->Data().NumSegments() == 0) {
          return std::nullopt;
        }
        RequestExtrapolationEntrances(structure, configure(), requests);
      }
      if (writesCheckFiles && ((configure().fileCheckMask | configure().screenCheckMask) & Config::CHECK_DATA))
        data()->PrintData(configure());
    } else {
      // The program's check of the rate's pairs (SetUpProgramCompound): exit -1 is the total radiative capture.
      const std::string rateProblem =
          structure.RatePairsProblem(configure().rateParams.entrancePair, configure().rateParams.exitPair);
      if (!rateProblem.empty()) {
        configure().outStream << rateProblem << std::endl;
        return std::nullopt;
      }
      requests.AddEntrance(structure.GetPairNumFromKey(configure().rateParams.entrancePair));
    }
    if (recomputesIntegrals) configure_.paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
    // The capture-A_y terms the data sets need, built with the pathway table
    RequestCaptureAyTerms(structure, *data(), requests);
    return requests;
  };
  CompoundSetUp setUp = SetUpCompoundSystem(configure(), dataSetUp, radii);
  if (setUp.refusal == CompoundSetUp::Refusal::kFinish) {
    if (saysAbort)
      PrintSetUpRefusal(setUp, configure().outStream);
    else
      configure().outStream << setUp.finishError << std::endl;
  }
  if (setUp.refusal != CompoundSetUp::Refusal::kNone) return setUp.refusal;
  compound_ = std::move(setUp.compound);
  *state_ = std::move(setUp.state);  // in place: a manager made later refers to the session's levels
  *inputs_ = std::move(setUp.inputs);
  return CompoundSetUp::Refusal::kNone;
}

bool AZUREAPI::RebuildImpl(const std::map<int, double> &radii) {
  const CompoundSetUp::Refusal refusal = SetUpModel(ModelSetUpOccasion::kRebuild, radii);
  if (refusal == CompoundSetUp::Refusal::kFinish) configure_.paramMask |= Config::USE_PREVIOUS_INTEGRALS;
  if (refusal != CompoundSetUp::Refusal::kNone) return false;
  if (data()->Initialize(compound(), configure()) == -1) {
    configure_.paramMask |= Config::USE_PREVIOUS_INTEGRALS;
    return false;
  }
  results_ = azure2::Results();  // as after the first Initialize
  data()->SizeResults(results_);
  configure_.paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  // Without this the client keeps the pre-change parameter vector, names,
  // fixed flags and Wigner limits -- silently mismatched with the new model.
  UpdateParameters();
  if (!SetUpPriorCovariance()) return false;

  return true;
}

void AZUREAPI::ReadFitSettings() {
  if (nuisanceRead_) return;
  nuisanceRead_ = true;
  // A scratch parameter set: the limits manager applies the settings to it, not to the session's parameters.
  nuisanceParams_.reset(new AZUREParams);
  compound()->FillMnParams(*nuisanceParams_, state(), inputs(), &configure());
  data()->FillMnParams(*nuisanceParams_);
  nuisance_.reset(new ParameterLimitsManager(&configure(), compound(), state(), inputs(), data(), nuisanceParams_.get()));
  if (nuisance_->ReadParameterSettings()) nuisance_->ApplyAllParameterSettings(*nuisanceParams_);
}

std::pair<double, double> AZUREAPI::GetParameterBounds(int index) {
  ReadFitSettings();
  const ROOT::Minuit2::MinuitParameter &limits =
      static_cast<const AZUREParams &>(*nuisanceParams_).GetMinuitParams().Parameter(index);
  return {limits.HasLowerLimit() ? limits.LowerLimit() : -std::numeric_limits<double>::infinity(),
          limits.HasUpperLimit() ? limits.UpperLimit() : std::numeric_limits<double>::infinity()};
}

// Function to get indeces of non-fixed normalization parameters
vector_r AZUREAPI::GetNormalizationIndices() {
  int k = 0;
  vector_r indices;
  int totalParams = all_rwa_.size();
  for (int i = 0; i < totalParams; ++i) {
    if (!fixed_[i]) {
      if (records_[i].kind == ParamKind::Norm) {
        indices.push_back(k);
      }
      ++k;
    }
  }

  return indices;
}

// Function to get indeces of non-fixed energy shift parameters
vector_r AZUREAPI::GetEnergyShiftIndices() {
  int k = 0;
  vector_r indices;
  int totalParams = all_rwa_.size();
  for (int i = 0; i < totalParams; ++i) {
    if (!fixed_[i]) {
      // Every energy-shift term: the constant, the sqrt(E) coefficient and the shift-function parameters
      const ParamKind kind = records_[i].kind;
      if (kind == ParamKind::EnergyShift || kind == ParamKind::EnergyShiftSqrt || kind == ParamKind::EnergyShiftFn) {
        indices.push_back(k);
      }
      ++k;
    }
  }

  return indices;
}
