// The set-up stages of a program run (ProgramRun.h): the compound nucleus, the data, the parameter vector, and for a
// data run the objective.
#include "ProgramRun.h"

#include <iostream>
#include <memory>
#include <string>

#include "AdaptiveExtrapolation.h"
#include "CNuc.h"
#include "CNucChecks.h"
#include "CNucReport.h"
#include "CaptureAyRequests.h"
#include "ProjectSetup.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "GSLException.h"
#include "Objective.h"
#include "PriorCovariance.h"
#include "Profile.h"
#include "WorkingCopy.h"

int SetUpRun(ProgramRun &run) {
  // Fill compound nucleus from nucfile
  run.configure.outStream << "Filling Compound Nucleus..." << std::endl;
  if (run.compound->Fill(run.configure) == -1) {
    run.configure.outStream << "Could not fill compound nucleus from file."
                          << std::endl;
    return -1;
  } else if (run.compound->NumPairs() == 0 || run.compound->NumJGroups() == 0) {
    run.configure.outStream << "No nuclear data exists. Calculation not possible." << std::endl;
    return -1;
  }
  if ((run.configure.screenCheckMask | run.configure.fileCheckMask) &
      Config::CHECK_COMPOUND_NUCLEUS) WriteCompoundNucleusCheck(*run.compound, run.configure);

  // An extrapolation run whose project has an adaptive grid refines it once the parameters are known (below).
  run.refinesAdaptiveGrids =
      !(run.configure.paramMask & (Config::CALCULATE_WITH_DATA | Config::CALCULATE_REACTION_RATE)) &&
      HasAdaptiveGrids(run.configure);

  if (!(run.configure.paramMask & Config::CALCULATE_REACTION_RATE)) {
    // Fill the data object from the segments and data file
    //   Compound object is passed to the function for pair key verification and
    //   center of mass conversions, s-factor conversions, etc.
    run.configure.outStream << "Filling Data Structures..." << std::endl;
    if (run.configure.paramMask & Config::CALCULATE_WITH_DATA) {
      if (const int filled = run.data->Fill(run.configure, run.compound); filled != 0) {
        if (filled != WorkingCopy::kFillRefused) run.configure.outStream << "Could not fill data object from file." << std::endl;
        return -1;
      } else if (run.data->Data().NumSegments() == 0 && run.compound->GetModelData().NumResiduals() == 0) {  // model data alone run
        run.configure.outStream << "There is no data provided." << std::endl;
        return -1;
      }
      RequestDataEntrances(*run.compound, run.configure);
    } else {
      if (const int made = run.data->MakePoints(run.configure, run.compound); made != 0) {
        if (made != WorkingCopy::kFillRefused) run.configure.outStream << "Could not fill data object from file." << std::endl;
        return -1;
      } else if (run.data->Data().NumSegments() == 0) {
        run.configure.outStream << "Extrapolation segments produce no data." << std::endl;
        return -1;
      }
      RequestExtrapolationEntrances(*run.compound, run.configure);
    }
    // The capture-A_y terms the data sets need are built with the pathway table (CNuc::Initialize below).
    RequestCaptureAyTerms(*run.compound, *run.data);
    // The data check lists the points that are calculated: an extrapolation run with an adaptive grid writes it after
    // the refinement below, on the refined grid.
    if (((run.configure.fileCheckMask | run.configure.screenCheckMask) & Config::CHECK_DATA) && !run.refinesAdaptiveGrids)
      run.data->PrintData(run.configure);

    // Cache will be populated dynamically as integrals are calculated
  } else {
    const std::string rateProblem =
        run.compound->RatePairsProblem(run.configure.rateParams.entrancePair, run.configure.rateParams.exitPair);
    if (!rateProblem.empty()) {
      run.configure.outStream << rateProblem << std::endl;
      return -1;
    } else {
      run.compound->GetPair(run.compound->GetPairNumFromKey(run.configure.rateParams.entrancePair))->SetEntrance();
      run.compound->RequestEntrance(run.compound->GetPairNumFromKey(run.configure.rateParams.entrancePair));
    }
  }

  // Initialize compound nucleus object
  try {
    run.compound->Initialize(run.configure);
  } catch (GSLException e) {
    run.configure.outStream << e.what() << std::endl;
    run.configure.outStream << std::endl
                          << "Calculation was aborted." << std::endl;
    return -1;
  }

  // Create new parameters for minuit, fill them from compound nucleus object and data file.
  run.compound->FillMnParams(run.params, &run.configure);
  run.data->FillMnParams(run.params);
  if (run.configure.paramMask & Config::USE_PREVIOUS_PARAMETERS) {
    run.configure.outStream << "Reading User Parameter File..." << std::endl;
    if (!run.params.ReadUserParameters(run.configure) || !run.params.ReconcileBasis(run.compound, run.configure)) return -1;
  }
  // An adaptive extrapolation grid is a function of the parameters (docs/upgrade/project_format3.md section 3.8): it is
  // refined now, at the run's parameters, and the extrapolations are made again on the refined grids before anything
  // is computed on them (AdaptiveExtrapolation.h).
  if (run.refinesAdaptiveGrids) {
    std::unique_ptr<WorkingCopy> refined = MakeRefinedExtrapolations(run.configure, run.compound, run.params.GetMinuitParams().Params());
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
    if (PriorCovariance::Load(run.configure, names, fixedFlags, AZUREParams::BasisForConfig(run.configure), run.priorCov) < 0)
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
  if (run.data->Initialize(run.compound, run.configure) == -1) {
    run.configure.outStream << std::endl
                          << "Calculation was aborted." << std::endl;
    return -1;
  }

  // Create ParameterLimitsManager and apply settings
  run.configure.outStream << "Applying parameter settings with Parameter Manager..." << std::endl;
  run.limits.reset(new ParameterLimitsManager(&run.configure, run.compound, run.data, &run.params));
  run.limits->ReadParameterSettings();
  run.limits->ApplyAllParameterSettings(run.params);

  // Warn once if a level has a radiative width comparable to its particle width,
  // where the R-Matrix treatment of capture breaks down.
  if (run.compound->ObservedWidthsAt(run.configure, run.params.GetMinuitParams().Params()))
    WarnLargeRadiativeWidths(*run.compound, run.configure);

  // AZURE_PROFILE: the set-up's stages (the channel table, the points' energy-dependent values) on their own, not
  // in the first evaluation's report.
  AZUREProfile::Report(run.configure.profile, std::cerr, "set-up");
  // The objective: the chi^2 of the data, the priors and the constraints (Objective.h), on the run's data.
  run.objective.reset(new AZURECalc(run.data, run.results, run.compound, run.configure, run.limits.get()));
  // The prior covariance replaces the diagonal prior of each of its parameters (docs/upgrade/objective.md D3).
  if (run.priorCov && run.priorCov->Size() > 0) {
    const std::vector<std::string> replaced = ReplaceDiagonalPriors(run.compound, run.data, *run.priorCov);
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
