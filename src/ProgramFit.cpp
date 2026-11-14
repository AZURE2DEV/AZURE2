// The program's fit stage (ProgramRun.h): the fitter of the program's switches (Fitters.h), its progress on standard
// output, HESSE for the band, MINOS, and the files of a fit: fitStatus.out, covariance_matrix.out, param.sav.
#include "ProgramRun.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <Minuit2/FunctionMinimum.h>

#include "CompoundSystem.h"
#include "Config.h"
#include "FitErrors.h"
#include "FitEvaluations.h"
#include "Fitters.h"
#include "ParameterLabel.h"
#include "WorkingCopy.h"

ProgramProgress::ProgramProgress(AZURECalc &objective, std::ostream &out) : objective_(objective), out_(out) {
  runControl_.progress = [this](const FitProgress &record) { Print(record); };
  runControl_.checkpoint = [this](const std::vector<double> &p) { objective_.WriteIterationOutput(p); };
  observer_.evaluated = [this](long evaluation, double value, const std::vector<double> &p, bool errorAnalysis) {
    FitProgress record;
    record.evaluation = evaluation;
    record.objective = value;
    record.errorAnalysis = errorAnalysis;
    record.x = &p;
    HandOverEvaluation(runControl_, record);
  };
  objective_.Observe(&observer_);
}

ProgramProgress::~ProgramProgress() { objective_.Observe(nullptr); }

FitControl ProgramProgress::Control(const Config &configure) const {
  FitControl control = runControl_;
  control.maxCalls = configure.maxCalls;
  control.callLimitSource = " (--max-calls)";
  control.refreezePasses = configure.refreezeCovariance;
  control.log = &configure.outStream;
  return control;
}

void ProgramProgress::Print(const FitProgress &record) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (record.kind == FitProgress::Kind::kIteration) {
    out_ << "\r\t" << record.label << ": " << std::setw(4) << record.iteration << "  Chi-Squared: " << record.objective
         << "        ";
    out_.flush();
    return;
  }
  // MINOS's evaluations (RunMinos under a FitEvaluations of stage "minos"): every 10th, counted from MINOS's start.
  if (std::string(record.stage) == "minos") {
    if (record.calls % 10 == 0)
      out_ << "\r\tMINOS evaluation: " << std::setw(6) << record.calls << " Chi-Squared: " << record.objective;
    out_.flush();
    return;
  }
  // An evaluation: every 10th of the run outside MINOS.
  if (record.errorAnalysis || record.evaluation == 0) return;
  if (record.evaluation % 10 == 0)
    out_ << "\r\tIteration: " << std::setw(6) << record.evaluation << " Chi-Squared: " << record.objective;
  out_.flush();
}

void ReportParametersAtLimits(const ProgramRun &run, const AZUREParams &params) {
  // A parameter pinned against one of its limits is the most common reason a
  // fit cannot produce a covariance, and the one the user can act on. Name any
  // such parameter physically -- the Minuit name alone does not say which
  // level or channel it belongs to.
  const ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  int nAtLimit = 0;
  for (unsigned int i = 0; i < mp.Params().size(); i++) {
    const ROOT::Minuit2::MinuitParameter &par = mp.Parameter(i);
    if (par.IsFixed() || !par.HasLimits()) continue;
    const double value = mp.Value(i);
    const double scale = std::max(std::fabs(value), 1.e-30);
    bool atLimit = false;
    if (par.HasLowerLimit() && std::fabs(value - par.LowerLimit()) <= 1.e-6 * scale) atLimit = true;
    if (par.HasUpperLimit() && std::fabs(value - par.UpperLimit()) <= 1.e-6 * scale) atLimit = true;
    if (!atLimit) continue;
    if (!nAtLimit) run.configure.outStream << "  Parameter(s) sitting at a limit:" << std::endl;
    std::string text = mp.GetName(i);
    const std::string label = AZURELabel::Parameter(run.compound.get(), run.inputs, run.data, (int)i);
    if (!label.empty()) text += "  =  " + label;
    run.configure.outStream << "    " << text << " = " << value << std::endl;
    nAtLimit++;
  }
  if (!nAtLimit)
    run.configure.outStream << "  No parameter is at a limit, so the flat direction is a "
                               "degeneracy between correlated parameters."
                            << std::endl;
}

void CaptureMinuitCovariance(ProgramRun &run, const std::vector<double> &covData, const AZUREParams &atLimits) {
  if (!(run.configure.paramMask & Config::CALCULATE_COVARIANCE_BAND)) return;
  const int nMn = run.params.GetMinuitParams().Params().size();
  std::vector<bool> fixedMask(nMn);
  for (int fi = 0; fi < nMn; fi++) fixedMask[fi] = run.params.GetMinuitParams().Parameter(fi).IsFixed();
  ParamIndexMap pmap = BuildParamIndexMap(run.compound.get(), run.data, fixedMask);
  run.bandCov = BuildBandCovarianceFromMinuit(covData, pmap);
  // Saved later (after the optional reduced-chi^2 scaling) in the band block.
  if (!run.bandCov.empty())
    run.haveBandCov = true;
  else {
    const int np = pmap.NumPacked();
    if ((int)covData.size() == np * (np + 1) / 2 && BandPackedColumns(pmap).empty()) {
      // No level parameter is free (or no parameter at all): reported where the band is skipped.
    } else if (covData.empty()) {
      run.configure.outStream << "Note: no covariance available for the uncertainty band -- the fit "
                                 "did not produce one. This usually means the minimum is not well "
                                 "determined (the fit did not fully converge, a parameter is sitting "
                                 "at a limit, or two parameters are degenerate)."
                              << std::endl;
      ReportParametersAtLimits(run, atLimits);
    } else
      run.configure.outStream << "Note: covariance/parameter size mismatch for the uncertainty band ("
                              << covData.size() << " covariance entries vs " << np * (np + 1) / 2
                              << " expected for " << np << " free parameters); band skipped." << std::endl;
  }
}

namespace {

// fitStatus.out: how the fit ended, for scripts and for anyone who reads only the output files.  Every fitter, the
// same keys: the fitter asked and the one that ran (rows of Fitters()), the outcome in the shared words (FitOutcome),
// the fitter's own status word, EDM and its goal ("none" for a fitter without them), calls, the call limit.
void WriteFitStatus(const Config &configure, const FitterRow &fitter, const FitResult &result,
                    const FitControl &control) {
  std::ofstream out(configure.OpenPath(configure.outputdir + "fitStatus.out").c_str());
  if (!out) return;
  const FitterRow &ran = FitterThatRan(fitter, result);
  const unsigned int noCallLimit = 1000000000u;
  std::string callLimitText = "none";
  if (control.maxCalls > 0)
    callLimitText = std::to_string(std::min(control.maxCalls, noCallLimit)) + control.callLimitSource +
                    (ran.callLimit ? "" : " ignored");
  out << "minimizer      " << ran.name << std::endl
      << "fitter         " << fitter.name << std::endl
      << "status         " << FitOutcome(result) << std::endl
      << "fitter_status  " << result.status << std::endl
      << "valid          " << (result.converged ? 1 : 0) << std::endl;
  if (result.edmGoal > 0.0)
    out << "edm            " << result.edm << std::endl << "edm_goal       " << result.edmGoal << std::endl;
  else
    out << "edm            none" << std::endl << "edm_goal       none" << std::endl;
  out << "calls          " << result.calls << std::endl << "call_limit     " << callLimitText << std::endl;
}

// covariance_matrix.out of a MINOS run: the parameters, Minuit's covariance and correlation matrices.
void WriteCovarianceMatrixFile(const Config &configure, AZUREParams &params, const ROOT::Minuit2::FunctionMinimum &min) {
  char filename[256];
  snprintf(filename, sizeof(filename), "%scovariance_matrix.out", configure.outputdir.c_str());
  std::ofstream out;
  out.open(configure.OpenPath(filename));
  if (!out) {
    std::cout << "Could not save " << configure.outputdir.c_str() << "covariance_matrix.out file." << std::endl;
    return;
  }
  // Write header information (params hold the fit's end point: AdoptMigradEndPoint)
  out << "Parameter List" << std::endl << std::endl;
  for (int i = 0; i < (int)params.GetMinuitParams().Params().size(); i++) {
    out << std::setw(20) << params.GetMinuitParams().GetName(i) << " (" << std::setw(3) << i << ")" << std::scientific
        << std::setw(20) << params.GetMinuitParams().Value(i);
    if (params.GetMinuitParams().Parameter(i).IsFixed())
      out << " Fixed" << std::endl;
    else
      out << " Fitted by Minuit" << std::endl;
  }
  std::cout << min.UserCovariance();
  std::vector<double> CovarianceData = min.UserCovariance().Data();
  std::vector<int> parameterTable;
  parameterTable.reserve(params.GetMinuitParams().Params().size());
  // Header for Covariance Matrix
  out << std::endl << "Covariance Matrix" << std::endl << std::endl;
  out << std::setw(5) << " ";
  for (int i = 0; i < (int)params.GetMinuitParams().Params().size(); i++) {
    if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
      out << std::setw(15) << i;
      parameterTable.push_back(i);
    }
  }
  out << std::endl;
  int size = parameterTable.size();
  std::cout << "covariance length " << CovarianceData.size() << std::endl;
  // Covariance matrix (Minuit's packed lower triangle)
  for (int i = 0; i < size; i++) {
    out << std::setw(5) << parameterTable[i];
    for (int j = 0; j < size; j++) {
      const int k = i <= j ? ((j + 1) * j) / 2 + i : ((i + 1) * i) / 2 + j;
      std::cout << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k << std::endl;
      out << std::setw(15) << CovarianceData[k];
    }
    out << std::endl;
  }
  // Header for Correlation Matrix
  out << std::endl << "Correlation Matrix" << std::endl << std::endl;
  out << std::setw(5) << " ";
  for (int i = 0; i < (int)params.GetMinuitParams().Params().size(); i++)
    if (!params.GetMinuitParams().Parameter(i).IsFixed()) out << std::setw(15) << i;
  out << std::endl;
  // Correlation matrix
  for (int i = 0; i < size; i++) {
    out << std::setw(5) << parameterTable[i];
    for (int j = 0; j < size; j++) {
      const int k = i <= j ? ((j + 1) * j) / 2 + i : ((i + 1) * i) / 2 + j;
      int jdiag = ((j + 2) * (j + 1)) / 2 - 1;
      int idiag = ((i + 2) * (i + 1)) / 2 - 1;
      std::cout << k << "  " << CovarianceData[k] << std::endl;
      out << std::setw(15) << std::fixed << CovarianceData[k] / sqrt(fabs(CovarianceData[jdiag] * CovarianceData[idiag]));
    }
    out << std::endl;
  }
  out.flush();
  out.close();
}

// The error analysis of a MIGRAD row: HESSE for the band (falling back to MIGRAD's matrix), then MINOS with its files.
void MigradErrorAnalysis(ProgramRun &run, FitResult &result, const AZUREParams &start, const FitControl &control) {
  const Config &configure = run.configure;
  ROOT::Minuit2::FunctionMinimum &min = *result.minimum;
  // MIGRAD's error matrix is approximate; refine it with HESSE for the band,
  // falling back to the MIGRAD matrix if HESSE gives no valid covariance.
  // With no free parameter there is nothing to refine, and MnHesse on such a state ends in a segmentation
  // fault inside Minuit2 (found by tests/fit_modes; the code is the same at the branch point with dev).
  // The parameters named as sitting at a limit are those of `start` (the run's parameters before the fit).
  if ((configure.paramMask & Config::CALCULATE_COVARIANCE_BAND) && !HasFreeParameterForHesse(result)) {
    CaptureMinuitCovariance(run, min.UserCovariance().Data(), start);
  } else if (configure.paramMask & Config::CALCULATE_COVARIANCE_BAND) {
    std::vector<double> migradCov = min.UserCovariance().Data();
    configure.outStream << std::endl << "Computing accurate covariance (HESSE) for the uncertainty band..." << std::endl;
    if (RunHesse(*run.objective, result)) {
      CaptureMinuitCovariance(run, min.UserCovariance().Data(), start);
    } else {
      configure.outStream << "HESSE did not return a valid covariance; using the MIGRAD "
                             "error matrix (uncertainties may be approximate)."
                          << std::endl;
      ReportParametersAtLimits(run, start);
      CaptureMinuitCovariance(run, migradCov, start);
    }
  }
  if (configure.paramMask & Config::PERFORM_ERROR_ANALYSIS) {
    configure.outStream << std::endl
                        << "Performing parameter error analysis with Up=" << configure.chiVariance << "." << std::endl;
    {
      FitEvaluations evaluations(*run.objective, control, result.fitter.c_str(), "minos");
      RunMinos(*run.objective, run.params, configure.chiVariance, result, configure.outStream);
    }
    // The fit's end point with Minuit's errors, for the files and the program's parameters.
    AdoptMigradEndPoint(run.params, result);
    WriteCovarianceMatrixFile(configure, run.params, min);
    run.params.WriteParameterErrors(result.minos, configure);
  }
}

}  // namespace

void RunFitStage(ProgramRun &run, const ProgramProgress &progress) {
  const Config &configure = run.configure;
  if (!(configure.paramMask & Config::PERFORM_FIT)) {
    if (configure.paramMask & Config::USE_AMATRIX)
      configure.outStream << "Performing A-Matrix Calculation..." << std::endl;
    else
      configure.outStream << "Performing R-Matrix Calculation..." << std::endl;
    return;
  }
  const std::string fitStatusFile = configure.outputdir + "fitStatus.out";
  std::remove(configure.OpenPath(fitStatusFile).c_str());  // a stale file from an earlier fit must not describe this one
  if (configure.paramMask & Config::USE_AMATRIX)
    configure.outStream << "Performing A-Matrix Fit..." << std::endl;
  else
    configure.outStream << "Performing R-Matrix Fit..." << std::endl;
  run.data->SetFit(true);

  const FitterRow &fitter = FitterOfConfig(configure);
  const FitControl control = progress.Control(configure);
  const AZUREParams start = run.params;
  FitResult result = RunFitter(fitter, *run.objective, run.params, control);
  WriteFitStatus(configure, fitter, result, control);

  if (result.fellBackToMigrad) {
    // A least-squares row without its analytic Jacobian ran MIGRAD: its matrix for the band, no HESSE, no MINOS.
    CaptureMinuitCovariance(run, result.minimum->UserCovariance().Data(), run.params);
  } else if (fitter.minos) {
    MigradErrorAnalysis(run, result, start, control);
    AdoptMigradEndPoint(run.params, result);
  } else if (fitter.residualRows) {
    if ((configure.paramMask & Config::CALCULATE_COVARIANCE_BAND) && !result.leastSquaresBand.empty()) {
      run.bandCov = result.leastSquaresBand;
      run.haveBandCov = true;  // saved later, after optional reduced-chi^2 scaling
    }
  } else if (configure.paramMask & Config::PERFORM_ERROR_ANALYSIS) {
    configure.outStream << std::endl << "Error analysis is not available with NLopt minimizer." << std::endl;
    configure.outStream << "Please use Minuit2 for error analysis." << std::endl;
  }
  run.data->SyncShiftGroupParams(run.params.GetMinuitParams());
  run.data->SyncNormGroupParams(run.params.GetMinuitParams());
  run.params.WriteUserParameters(configure, true);
}
