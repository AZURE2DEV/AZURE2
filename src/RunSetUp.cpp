// The set-up stages of a program run (ProgramRun.h): the compound nucleus and the data (the sampler's set-up too),
// the parameter vector, and for a data run the objective.
#include "ProgramRun.h"

#include <iostream>
#include <memory>
#include <optional>
#include <string>

#include "AdaptiveExtrapolation.h"
#include "CompoundSystem.h"
#include "CompoundSystemSetUp.h"
#include "CompoundSystemChecks.h"
#include "LevelReport.h"
#include "CaptureAyRequests.h"
#include "ProjectSetup.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "GSLException.h"
#include "Objective.h"
#include "PriorCovarianceFile.h"
#include "Profile.h"
#include "WorkingCopy.h"

CompoundSetUp SetUpProgramCompound(const Config &configure, WorkingCopy &data, DataCheckAt dataCheck) {
  configure.outStream << "Filling Compound Nucleus..." << std::endl;
  const DataSetUpStep dataSetUp = [&configure, &data, dataCheck](
                                      const CompoundSystem &structure,
                                      const LevelInputs &inputs) -> std::optional<CompoundSetUpRequests> {
    if ((configure.screenCheckMask | configure.fileCheckMask) & Config::CHECK_COMPOUND_NUCLEUS)
      WriteCompoundNucleusCheck(structure, inputs, configure);
    CompoundSetUpRequests requests;
    if (!(configure.paramMask & Config::CALCULATE_REACTION_RATE)) {
      // The data object from the segments and the data files: the structure gives the pair keys, the centre-of-mass
      // conversions and the S-factor conversions.
      configure.outStream << "Filling Data Structures..." << std::endl;
      if (configure.paramMask & Config::CALCULATE_WITH_DATA) {
        if (const int filled = data.Fill(configure, &structure, inputs); filled != 0) {
          if (filled != WorkingCopy::kFillRefused) configure.outStream << "Could not fill data object from file." << std::endl;
          return std::nullopt;
        } else if (data.Data().NumSegments() == 0 && structure.GetModelData().NumResiduals() == 0) {  // model data alone run
          configure.outStream << "There is no data provided." << std::endl;
          return std::nullopt;
        }
        RequestDataEntrances(structure, configure, requests);
      } else {
        if (const int made = data.MakePoints(configure, &structure, inputs); made != 0) {
          if (made != WorkingCopy::kFillRefused) configure.outStream << "Could not fill data object from file." << std::endl;
          return std::nullopt;
        } else if (data.Data().NumSegments() == 0) {
          configure.outStream << "Extrapolation segments produce no data." << std::endl;
          return std::nullopt;
        }
        RequestExtrapolationEntrances(structure, configure, requests);
      }
      // The capture-A_y terms the data sets need are built with the pathway table (the compound's last part).
      RequestCaptureAyTerms(structure, data, requests);
      // The data check lists the points that are calculated: on a grid refined later, it is written after that.
      if (((configure.fileCheckMask | configure.screenCheckMask) & Config::CHECK_DATA) && dataCheck == DataCheckAt::kDataSetUp)
        data.PrintData(configure);
    } else {
      const std::string rateProblem =
          structure.RatePairsProblem(configure.rateParams.entrancePair, configure.rateParams.exitPair);
      if (!rateProblem.empty()) {
        configure.outStream << rateProblem << std::endl;
        return std::nullopt;
      }
      requests.AddEntrance(structure.GetPairNumFromKey(configure.rateParams.entrancePair));
    }
    return requests;
  };
  CompoundSetUp setUp = SetUpCompoundSystem(configure, dataSetUp);
  PrintSetUpRefusal(setUp, configure.outStream);
  return setUp;
}

int SetUpRun(ProgramRun &run) {
  // An extrapolation run whose project has an adaptive grid refines it once the parameters are known (below).
  run.refinesAdaptiveGrids =
      !(run.configure.paramMask & (Config::CALCULATE_WITH_DATA | Config::CALCULATE_REACTION_RATE)) &&
      HasAdaptiveGrids(run.configure);

  // The compound system from the <levels> rows, the run's data set up in the middle, then its last part.
  CompoundSetUp setUp = SetUpProgramCompound(
      run.configure, *run.data, run.refinesAdaptiveGrids ? DataCheckAt::kAfterGridRefinement : DataCheckAt::kDataSetUp);
  if (setUp.refusal != CompoundSetUp::Refusal::kNone) return -1;
  run.compound = std::move(setUp.compound);
  run.state = std::move(setUp.state);
  run.inputs = std::move(setUp.inputs);

  // Create new parameters for minuit, fill them from compound nucleus object and data file.
  run.compound->FillMnParams(run.params, run.state, run.inputs, &run.configure);
  run.data->FillMnParams(run.params);
  if (run.configure.paramMask & Config::USE_PREVIOUS_PARAMETERS) {
    run.configure.outStream << "Reading User Parameter File..." << std::endl;
    if (!run.params.ReadUserParameters(run.configure) || !run.params.ReconcileBasis(run.compound.get(), run.configure)) return -1;
  }
  // An adaptive extrapolation grid is a function of the parameters (docs/upgrade/project_format3.md section 3.8): it is
  // refined now, at the run's parameters, and the extrapolations are made again on the refined grids before anything
  // is computed on them (AdaptiveExtrapolation.h).
  if (run.refinesAdaptiveGrids) {
    std::unique_ptr<WorkingCopy> refined = MakeRefinedExtrapolations(run.configure, run.compound.get(), run.state, run.inputs, run.params.GetMinuitParams().Params());
    if (!refined) return -1;
    delete run.data;
    run.data = refined.release();
    if ((run.configure.fileCheckMask | run.configure.screenCheckMask) & Config::CHECK_DATA) run.data->PrintData(run.configure);
  }
  // Prior parameter covariance (<config> priorCovariance=): matched onto the parameter vector as it stands now (a
  // parameter file can fix or free level parameters), before anything is written, so a refused file leaves no
  // output.  Only a run with data has a chi-squared to add it to.
  if ((run.configure.paramMask & Config::CALCULATE_WITH_DATA) &&
      !(run.configure.paramMask & Config::CALCULATE_REACTION_RATE) && !run.configure.priorCovariance.empty()) {
    const ROOT::Minuit2::MnUserParameters &mp = run.params.GetMinuitParams();
    std::vector<std::string> names;
    std::vector<bool> fixedFlags;
    for (unsigned int i = 0; i < mp.Params().size(); i++) {
      names.push_back(mp.GetName(i));
      fixedFlags.push_back(mp.Parameter(i).IsFixed());
    }
    if (LoadPriorCovariance(run.configure, names, fixedFlags, AZUREParams::BasisForConfig(run.configure), run.priorCov) < 0)
      return -1;
  }
  run.data->SetPriorCovariance(run.priorCov);
  if (!(run.configure.paramMask & Config::USE_PREVIOUS_PARAMETERS)) {
    run.configure.outStream << "Creating New param.par File..." << std::endl;
    run.params.WriteUserParameters(run.configure, false);
  }
  return 0;
}

int SetUpObjective(ProgramRun &run) {
  // Initialize shared EC amplitude cache for interpolation
  if (run.configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
    InitializeECAmplitudeCache(run.configure);
    run.configure.outStream << "Initialized and cleared EC amplitude cache..." << std::endl;
  }

  // Initialize data object
  if (run.data->Initialize(run.compound.get(), run.configure) == -1) {
    run.configure.outStream << std::endl
                          << "Calculation was aborted." << std::endl;
    return -1;
  }

  // Create ParameterLimitsManager and apply settings
  run.configure.outStream << "Applying parameter settings with Parameter Manager..." << std::endl;
  run.limits.reset(new ParameterLimitsManager(&run.configure, run.compound.get(), run.state, run.inputs, run.data, &run.params));
  run.limits->ReadParameterSettings();
  run.limits->ApplyAllParameterSettings(run.params);

  // Warn once if a level has a radiative width comparable to its particle width,
  // where the R-Matrix treatment of capture breaks down.
  if (const std::optional<ObservedLevelTable> observed =
          run.compound->ObservedWidthsAt(run.configure, run.params.GetMinuitParams().Params(), run.state, run.inputs))
  {
    WarnLargeRadiativeWidths(*run.compound.get(), run.inputs, *observed, run.configure);
    WarnLevelDoublesPotential(*run.compound.get(), run.inputs, *observed, run.configure);
  }

  // AZURE_PROFILE: the set-up's stages (the channel table, the points' energy-dependent values) on their own, not
  // in the first evaluation's report.
  AZUREProfile::Report(run.configure.profile, std::cerr, "set-up");
  // The objective: the chi^2 of the data, the priors and the constraints (Objective.h), on the run's data.
  run.objective.reset(new AZURECalc(run.data, run.results, run.compound.get(), run.state, run.inputs, *run.fitterStates, run.configure, run.limits.get()));
  // The prior covariance replaces the diagonal prior of each of its parameters (docs/upgrade/objective.md D3).
  if (run.priorCov && run.priorCov->Size() > 0) {
    const std::vector<std::string> replaced = ReplaceDiagonalPriors(run.compound.get(), run.data, *run.priorCov);
    std::string list;
    for (const std::string &name : replaced) list += (list.empty() ? "" : ", ") + name;
    if (!list.empty())
      run.configure.Note("Note: prior covariance: its multivariate Gaussian replaces the prior of " + list +
                         " of its own (normalization error, shift error or sigma); a nuisance prior on a listed "
                         "parameter is not applied either.");
  }
  run.objective->SetErrorDef(1.0);
  // The fixed flags the minimizers and the band honour: the .azr's, then an
  // external param.par's, then the <parameterSettings> (MIGRAD reads them off
  // run.params; see AZURECalc::FixedMask).
  run.fitFixedMask.assign(run.params.GetMinuitParams().Params().size(), false);
  for (size_t fi = 0; fi < run.fitFixedMask.size(); fi++)
    run.fitFixedMask[fi] = run.params.GetMinuitParams().Parameter(fi).IsFixed();
  run.objective->SetFixedMask(run.fitFixedMask);
  return 0;
}
