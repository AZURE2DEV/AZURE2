#include <cmath>
#include "AZURECalc.h"
#include "AZUREMain.h"
#include "ScatteringLength.h"
#include "ParameterLabel.h"
#include "AZUREParams.h"
#include "Config.h"
#include "CovarianceBand.h"
#include "PriorCovariance.h"
#include "ReactionRate.h"
#include "ParameterLimitsManager.h"
#include "ECAmplitudeCache.h"
#include "CoulFuncCache.h"
#include "Minuit2/MnPrint.h"
#include "GSLException.h"
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnHesse.h>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <cstdio>
#include <iomanip>

#ifdef USE_NLOPT
#include <nlopt.h>

// Structure to pass data to the nlopt objective function
struct NLoptData {
  const AZURECalc *theFunc;
  std::vector<double> *full_params;         // Full parameter vector (includes fixed params)
  std::vector<unsigned> *variable_indices;  // Indices of variable parameters
};

// NLopt objective function wrapper
// This wrapper receives only the VARIABLE parameters from nlopt,
// reconstructs the full parameter vector with fixed parameters,
// and passes it to AZURECalc
double nlopt_objective(unsigned n, const double *x, double *grad, void *f_data) {
  NLoptData *data = static_cast<NLoptData *>(f_data);

  // Reconstruct full parameter vector by inserting variable params at correct positions
  for (unsigned j = 0; j < n; j++) {
    unsigned i = (*data->variable_indices)[j];
    (*data->full_params)[i] = x[j];
  }

  // Calculate chi-squared with full parameter vector
  double chi_squared = (*data->theFunc)(*data->full_params);

  // Gradient computation not implemented (using derivative-free algorithms)
  if (grad) {
    // NLopt requires gradient to be set to NULL or computed
    // For derivative-free algorithms, this won't be called
  }

  return chi_squared;
}
#endif

int AZUREMain::operator()() {
  // Clear Coulomb function cache at the start of each new calculation
  // This ensures the cache doesn't retain stale data from previous calculations
  // or nuclear potential changes from the GUI
  if (g_coulFuncCache) {
    g_coulFuncCache->Clear();
    configure().outStream << "Cleared Coulomb function cache." << std::endl;
  }

  // Fill compound nucleus from nucfile
  configure().outStream << "Filling Compound Nucleus..." << std::endl;
  if (compound()->Fill(configure()) == -1) {
    configure().outStream << "Could not fill compound nucleus from file."
                          << std::endl;
    return -1;
  } else if (compound()->NumPairs() == 0 || compound()->NumJGroups() == 0) {
    configure().outStream << "No nuclear data exists. Calculation not possible." << std::endl;
    return -1;
  }
  if ((configure().screenCheckMask | configure().fileCheckMask) &
      Config::CHECK_COMPOUND_NUCLEUS) compound()->PrintNuc(configure());

  if (!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    // Fill the data object from the segments and data file
    //   Compound object is passed to the function for pair key verification and
    //   center of mass conversions, s-factor conversions, etc.
    configure().outStream << "Filling Data Structures..." << std::endl;
    if (configure().paramMask & Config::CALCULATE_WITH_DATA) {
      if (data()->Fill(configure(), compound()) == -1) {
        configure().outStream << "Could not fill data object from file." << std::endl;
        return -1;
      } else if (data()->NumSegments() == 0) {
        configure().outStream << "There is no data provided." << std::endl;
        return -1;
      }
    } else {
      if (data()->MakePoints(configure(), compound()) == -1) {
        configure().outStream << "Could not fill data object from file." << std::endl;
        return -1;
      } else if (data()->NumSegments() == 0) {
        configure().outStream << "Extrapolation segments produce no data." << std::endl;
        return -1;
      }
    }
    if ((configure().fileCheckMask | configure().screenCheckMask) & Config::CHECK_DATA)
      data()->PrintData(configure());

    // Cache will be populated dynamically as integrals are calculated
  } else {
    const std::string rateProblem =
        compound()->RatePairsProblem(configure().rateParams.entrancePair, configure().rateParams.exitPair);
    if (!rateProblem.empty()) {
      configure().outStream << rateProblem << std::endl;
      return -1;
    } else {
      compound()->GetPair(compound()->GetPairNumFromKey(configure().rateParams.entrancePair))->SetEntrance();
    }
  }

  // Initialize compound nucleus object
  try {
    compound()->Initialize(configure());
  } catch (GSLException e) {
    configure().outStream << e.what() << std::endl;
    configure().outStream << std::endl
                          << "Calculation was aborted." << std::endl;
    return -1;
  }

  // Create new parameters for minuit, fill them from compound nucleus object and data file.
  AZUREParams params;
  compound()->FillMnParams(params.GetMinuitParams(), &configure());
  data()->FillMnParams(params.GetMinuitParams());
  if (configure().paramMask & Config::USE_PREVIOUS_PARAMETERS) {
    configure().outStream << "Reading User Parameter File..." << std::endl;
    params.ReadUserParameters(configure());
    if (!params.ReconcileBasis(compound(), configure())) return -1;
  }
  // Prior parameter covariance (<config> priorCovariance=): matched onto the parameter vector as it stands now (a
  // parameter file can fix or free level parameters), before anything is written, so a refused file leaves no
  // output.  Only a run with data has a chi-squared to add it to.
  std::shared_ptr<const PriorCovariance> priorCov;
  if ((configure().paramMask & Config::CALCULATE_WITH_DATA) &&
      !(configure().paramMask & Config::CALCULATE_REACTION_RATE) && !configure().priorCovariance.empty()) {
    const ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
    std::vector<std::string> names;
    std::vector<bool> fixedFlags;
    for (unsigned int i = 0; i < mp.Params().size(); i++) {
      names.push_back(mp.GetName(i));
      fixedFlags.push_back(mp.Parameter(i).IsFixed());
    }
    if (PriorCovariance::Load(configure(), names, fixedFlags, AZUREParams::BasisForConfig(configure()), priorCov) < 0)
      return -1;
  }
  data()->SetPriorCovariance(priorCov);
  if (!(configure().paramMask & Config::USE_PREVIOUS_PARAMETERS)) {
    configure().outStream << "Creating New param.par File..." << std::endl;
    params.WriteUserParameters(configure(), false);
  }

  if (!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    // Initialize shared EC amplitude cache for interpolation
    if (configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
      InitializeECAmplitudeCache();
      // Clear cache to ensure fresh start for this calculation
      if (g_ecAmplitudeCache) {
        g_ecAmplitudeCache->Clear();
      }
      configure().outStream << "Initialized and cleared EC amplitude cache..." << std::endl;
    }

    // Initialize data object
    if (data()->Initialize(compound(), configure()) == -1) {
      configure().outStream << std::endl
                            << "Calculation was aborted." << std::endl;
      return -1;
    }

    // Create ParameterLimitsManager and apply settings
    configure().outStream << "Applying parameter settings with Parameter Manager..." << std::endl;
    ParameterLimitsManager limitsManager(&configure(), compound(), data(), &params);
    limitsManager.ReadParameterSettings();
    limitsManager.ApplyAllParameterSettings(params.GetMinuitParams());

    // Warn once if a level has a radiative width comparable to its particle width,
    // where the R-Matrix treatment of capture breaks down.
    compound()->CheckRadiativeWidths(configure(), params.GetMinuitParams().Params());

    // Declare a new instance of FCNBase
    AZURECalc theFunc(data(), compound(), configure(), &limitsManager);
    // A prior-covariance parameter that also carries a diagonal prior of its own gets both (on purpose or not).
    if (priorCov && priorCov->Size() > 0) {
      const std::vector<int> own = theFunc.IndicesWithOwnPrior();
      std::string list;
      for (int f : priorCov->FullIndex())
        if (std::binary_search(own.begin(), own.end(), f))
          list += (list.empty() ? "" : ", ") + params.GetMinuitParams().GetName(f);
      if (!list.empty())
        configure().Note("Note: prior covariance: " + list + " also carr" +
                         (list.find(',') == std::string::npos ? "ies" : "y") +
                         " a prior of its own (normalization error, shift error, sigma or nuisance prior); both are "
                         "applied, which counts the same information twice if the covariance already holds it.");
    }
    // theFunc.InitializePools();
    theFunc.SetErrorDef(1.0);

    // Parameter covariance for the cross-section band, saved for reuse by later
    // extrapolation runs.  Any fit (MIGRAD or LM) yields it; MINOS is not needed.
    BandCovariance bandCov;
    bool haveBandCov = false;
    std::string covPath = configure().outputdir + "covariance.dat";

    // A parameter pinned against one of its limits is the most common reason a
    // fit cannot produce a covariance, and the one the user can act on. Name any
    // such parameter physically -- the Minuit name alone does not say which
    // level or channel it belongs to.
    auto reportParametersAtLimits = [&]() {
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
        if (!nAtLimit)
          configure().outStream << "  Parameter(s) sitting at a limit:" << std::endl;
        std::string text = mp.GetName(i);
        const std::string label = AZURELabel::Parameter(compound(), data(), (int)i);
        if (!label.empty()) text += "  =  " + label;
        configure().outStream << "    " << text << " = " << value << std::endl;
        nAtLimit++;
      }
      if (!nAtLimit)
        configure().outStream << "  No parameter is at a limit, so the flat direction is a "
                                 "degeneracy between correlated parameters."
                              << std::endl;
    };

    // Build the band covariance from a Minuit covariance (MIGRAD paths).
    auto captureMinuitCov = [&](const std::vector<double> &covData) {
      if (!(configure().paramMask & Config::CALCULATE_COVARIANCE_BAND)) return;
      const int nMn = params.GetMinuitParams().Params().size();
      std::vector<bool> fixedMask(nMn);
      for (int fi = 0; fi < nMn; fi++) fixedMask[fi] = params.GetMinuitParams().Parameter(fi).IsFixed();
      ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), fixedMask);
      bandCov = BuildBandCovarianceFromMinuit(covData, pmap);
      // Saved later (after the optional reduced-chi^2 scaling) in the band block.
      if (!bandCov.empty())
        haveBandCov = true;
      else {
        const int np = pmap.NumPacked();
        if (covData.empty()) {
          configure().outStream << "Note: no covariance available for the uncertainty band -- the fit "
                                   "did not produce one. This usually means the minimum is not well "
                                   "determined (the fit did not fully converge, a parameter is sitting "
                                   "at a limit, or two parameters are degenerate)."
                                << std::endl;
          reportParametersAtLimits();
        } else
          configure().outStream << "Note: covariance/parameter size mismatch for the uncertainty band ("
                                << covData.size() << " covariance entries vs " << np * (np + 1) / 2
                                << " expected for " << np << " free parameters); band skipped." << std::endl;
      }
    };


    // Call limit of a MIGRAD or NLopt fit: --max-calls=N, or none, so that MIGRAD's EDM criterion (NLopt: its
    // tolerances) ends the fit.  Before 2026-10-31 it was a fixed 50,000: the 19F (15N+a) refit, 145 free
    // parameters, stopped there with no message, and its chiSquared.out and param.sav looked like a converged fit.
    // Minuit2 scales the limit by 1.3 in an int, so "none" is 10^9 calls.
    const unsigned int noCallLimit = 1000000000u;
    const unsigned int migradCalls = configure().maxCalls > 0 ? std::min(configure().maxCalls, noCallLimit) : noCallLimit;
    const std::string callLimitText =
        configure().maxCalls > 0 ? std::to_string(migradCalls) + " (--max-calls)" : std::string("none");
    const std::string fitStatusFile = configure().outputdir + "fitStatus.out";
    // fitStatus.out: how the minimizer ended, for scripts and for anyone who reads only the output files.
    auto writeFitStatus = [&](const std::string &minimizer, const std::string &status, bool valid, double edm,
                              double edmGoal, long calls) {
      std::ofstream out(fitStatusFile.c_str());
      if (!out) return;
      out << "minimizer   " << minimizer << std::endl
          << "status      " << status << std::endl
          << "valid       " << (valid ? 1 : 0) << std::endl;
      if (edmGoal > 0.0) out << "edm         " << edm << std::endl << "edm_goal    " << edmGoal << std::endl;
      out << "calls       " << calls << std::endl
          << "call_limit  " << callLimitText << std::endl;
    };
    // Say how MIGRAD ended; a fit that stopped at the call limit or with an invalid minimum is flagged.
    auto reportMigrad = [&](const ROOT::Minuit2::FunctionMinimum &m) {
      const double goal = 0.002 * 0.1 * m.Up();  // MnMigrad's default tolerance 0.1 (VariableMetricBuilder)
      std::string status;
      std::ostringstream why;
      why << std::scientific << std::setprecision(3);
      if (m.IsValid()) {
        status = "converged";
        configure().outStream << std::endl << "MIGRAD converged: valid minimum, EDM " << std::scientific << std::setprecision(3)
                              << m.Edm() << " (goal " << goal << ")" << std::defaultfloat << std::setprecision(6)
                              << ", " << m.NFcn() << " function calls, call limit " << callLimitText << "."
                              << std::endl;
      } else {
        if (m.HasReachedCallLimit()) {
          status = "call_limit";
          why << "stopped at the call limit (" << m.NFcn() << " function calls, limit " << callLimitText << ")";
        } else {
          status = "invalid";
          why << "invalid minimum";
          if (m.IsAboveMaxEdm()) why << ", EDM above its goal";
          if (!m.HasValidCovariance()) why << ", no valid covariance";
          if (m.HasMadePosDefCovar()) why << ", covariance forced positive definite";
          if (m.HesseFailed()) why << ", HESSE failed";
          why << "; " << m.NFcn() << " function calls";
        }
        why << "; EDM " << m.Edm() << " (goal " << goal << ")";
        configure().outStream << std::endl << "WARNING: MIGRAD did NOT converge: " << why.str()
                              << ".  The parameters and outputs written are its last point, not a minimum; rerun"
                              << " the fit from param.sav" << (status == "call_limit" ? " or raise --max-calls" : "")
                              << "." << std::endl;
      }
      writeFitStatus("MIGRAD", status, m.IsValid(), m.Edm(), goal, m.NFcn());
    };

    if (configure().paramMask & Config::PERFORM_FIT) {
      std::remove(fitStatusFile.c_str());  // a stale file from an earlier fit must not describe this one
      // Call minimizer for function minimization, write minimized parameters to params
      if (configure().paramMask & Config::USE_AMATRIX)
        configure().outStream << "Performing A-Matrix Fit..." << std::endl;
      else
        configure().outStream << "Performing R-Matrix Fit..." << std::endl;
      data()->SetFit(true);

#ifdef USE_NLOPT
      if (configure().paramMask & Config::USE_NLOPT_MINIMIZER) {
        // Use NLopt for minimization
        configure().outStream << "Using NLopt minimizer..." << std::endl;

        // Build full parameter vector and identify variable parameters
        unsigned n_total = params.GetMinuitParams().Params().size();
        std::vector<double> full_params(n_total);
        std::vector<unsigned> variable_indices;

        // Extract all parameters and identify which are variable
        for (unsigned i = 0; i < n_total; i++) {
          full_params[i] = params.GetMinuitParams().Value(i);
          if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
            variable_indices.push_back(i);
          }
        }

        // Only pass variable parameters to nlopt
        unsigned n_variable = variable_indices.size();
        std::vector<double> x(n_variable);
        std::vector<double> lb(n_variable);
        std::vector<double> ub(n_variable);

        // Extract initial values and bounds for variable parameters only
        for (unsigned j = 0; j < n_variable; j++) {
          unsigned i = variable_indices[j];
          x[j] = params.GetMinuitParams().Value(i);

          if (params.GetMinuitParams().Parameter(i).HasLowerLimit())
            lb[j] = params.GetMinuitParams().Parameter(i).LowerLimit();
          else
            lb[j] = -1e10;  // Large negative number for unbounded

          if (params.GetMinuitParams().Parameter(i).HasUpperLimit())
            ub[j] = params.GetMinuitParams().Parameter(i).UpperLimit();
          else
            ub[j] = 1e10;  // Large positive number for unbounded
        }

        configure().outStream << "Optimizing " << n_variable << " parameters (out of "
                              << n_total << " total)..." << std::endl;

        // Map algorithm index to nlopt algorithm enum
        nlopt_algorithm algorithms[] = {
            NLOPT_LN_SBPLX,      // 0: SBPLX (Subplex)
            NLOPT_LN_COBYLA,     // 1: COBYLA
            NLOPT_LN_BOBYQA,     // 2: BOBYQA
            NLOPT_LN_NEWUOA,     // 3: NEWUOA
            NLOPT_LN_PRAXIS,     // 4: PRAXIS
            NLOPT_LN_NELDERMEAD  // 5: Nelder-Mead
        };
        const char *algorithm_names[] = {
            "SBPLX", "COBYLA", "BOBYQA", "NEWUOA", "PRAXIS", "Nelder-Mead"};

        int alg_index = configure().nloptAlgorithm;

        configure().outStream << "Using NLopt " << algorithm_names[alg_index]
                              << " algorithm..." << std::endl;

        // Create NLopt optimizer with selected algorithm
        nlopt_opt opt = nlopt_create(algorithms[alg_index], n_variable);

        // Set up objective function with wrapper data
        NLoptData nlopt_data;
        nlopt_data.theFunc = &theFunc;
        nlopt_data.full_params = &full_params;
        nlopt_data.variable_indices = &variable_indices;
        nlopt_set_min_objective(opt, nlopt_objective, &nlopt_data);

        // Set bounds
        nlopt_set_lower_bounds(opt, lb.data());
        nlopt_set_upper_bounds(opt, ub.data());

        // Set stopping criteria
        nlopt_set_maxeval(opt, configure().maxCalls > 0 ? static_cast<int>(migradCalls) : 0);  // 0: no limit
        nlopt_set_ftol_rel(opt, 1e-8);
        nlopt_set_xtol_rel(opt, 1e-8);

        // Perform optimization
        double minf;
        nlopt_result result = nlopt_optimize(opt, x.data(), &minf);

        const long nloptCalls = nlopt_get_numevals(opt);
        if (result == NLOPT_MAXEVAL_REACHED) {
          configure().outStream << std::endl
                                << "WARNING: NLopt did NOT converge: stopped at the call limit (" << nloptCalls
                                << " function calls, limit " << callLimitText << "); final chi-squared " << minf
                                << ".  Rerun the fit from param.sav or raise --max-calls." << std::endl;
          writeFitStatus("NLopt", "call_limit", false, 0.0, 0.0, nloptCalls);
        } else if (result < 0) {
          configure().outStream << std::endl;
          configure().outStream << "NLopt optimization failed: "
                                << nlopt_result_to_string(result) << std::endl;
          writeFitStatus("NLopt", std::string("failed_") + nlopt_result_to_string(result), false, 0.0, 0.0, nloptCalls);
        } else {
          writeFitStatus("NLopt", std::string("converged_") + nlopt_result_to_string(result), true, 0.0, 0.0,
                         nloptCalls);
          configure().outStream << std::endl;
          configure().outStream << "NLopt optimization succeeded: "
                                << nlopt_result_to_string(result) << std::endl;
          configure().outStream << "Final chi-squared: " << minf << std::endl;
        }

        // Update only the variable parameters with optimized values
        // Fixed parameters remain unchanged in full_params
        for (unsigned j = 0; j < n_variable; j++) {
          unsigned i = variable_indices[j];
          params.GetMinuitParams().SetValue(i, x[j]);
        }

        nlopt_destroy(opt);

        // Note: Error analysis with NLopt is not directly available
        // Skip Minos error analysis for NLopt
        if (configure().paramMask & Config::PERFORM_ERROR_ANALYSIS) {
          configure().outStream << std::endl
                                << "Error analysis is not available with NLopt minimizer." << std::endl;
          configure().outStream << "Please use Minuit2 for error analysis." << std::endl;
        }

        data()->SyncShiftGroupParams(params.GetMinuitParams());
        data()->SyncNormGroupParams(params.GetMinuitParams());
        params.WriteUserParameters(configure(), true);
      } else {
#endif
        if (configure().paramMask & Config::USE_GSL_LM_MINIMIZER) {
          // GSL trust-region least-squares (geodesic-accelerated); same analytic
          // Jacobian as the LM path, but factorizes J directly (QR).
          configure().outStream << "Using GSL trust-region least-squares minimizer with analytic Jacobian." << std::endl;
          bool wantBand = configure().paramMask & Config::CALCULATE_COVARIANCE_BAND;
          double gslChi2 = theFunc.RunGSLNonlinear(params, 200, wantBand ? &bandCov : nullptr);
          if (gslChi2 < 0.0) {
            // The analytic Jacobian is not available for this model -> MIGRAD.
            configure().outStream << "Analytic Jacobian unsupported for this model; "
                                     "falling back to MIGRAD."
                                  << std::endl;
            const ROOT::Minuit2::FCNBase &fb = theFunc;
            ROOT::Minuit2::MnMigrad migrad(fb, params.GetMinuitParams());
            ROOT::Minuit2::FunctionMinimum min = migrad(migradCalls);
            reportMigrad(min);
            params.GetMinuitParams() = min.UserParameters();
            captureMinuitCov(min.UserCovariance().Data());
          } else if (wantBand && !bandCov.empty()) {
            haveBandCov = true;  // saved later, after optional reduced-chi^2 scaling
          }
          data()->SyncShiftGroupParams(params.GetMinuitParams());
          data()->SyncNormGroupParams(params.GetMinuitParams());
          params.WriteUserParameters(configure(), true);
        } else if (configure().paramMask & Config::USE_LM_MINIMIZER) {
          // Levenberg-Marquardt / Gauss-Newton using the analytic residual
          // Jacobian (selected as an alternative to MIGRAD).
          configure().outStream << "Using Levenberg-Marquardt minimizer with analytic Jacobian." << std::endl;
          // LM yields the covariance (J^T J + priors)^{-1} directly for the band.
          bool wantBand = configure().paramMask & Config::CALCULATE_COVARIANCE_BAND;
          double lmChi2 = theFunc.RunLevenbergMarquardt(params, 200, wantBand ? &bandCov : nullptr);
          if (lmChi2 < 0.0) {
            // The analytic Jacobian is not available for this model -> MIGRAD.
            configure().outStream << "Analytic Jacobian unsupported for this model; "
                                     "falling back to MIGRAD."
                                  << std::endl;
            const ROOT::Minuit2::FCNBase &fb = theFunc;
            ROOT::Minuit2::MnMigrad migrad(fb, params.GetMinuitParams());
            ROOT::Minuit2::FunctionMinimum min = migrad(migradCalls);
            reportMigrad(min);
            params.GetMinuitParams() = min.UserParameters();
            captureMinuitCov(min.UserCovariance().Data());
          } else if (wantBand && !bandCov.empty()) {
            haveBandCov = true;  // saved later, after optional reduced-chi^2 scaling
          }
          data()->SyncShiftGroupParams(params.GetMinuitParams());
          data()->SyncNormGroupParams(params.GetMinuitParams());
          params.WriteUserParameters(configure(), true);
        } else {
          // Use Minuit2 MIGRAD for minimization (default).  When the analytic-
          // gradient runtime option is on, pass theFunc as an FCNGradientBase so
          // Migrad uses AZURECalc::Gradient; otherwise pass it as a plain FCNBase
          // so Minuit falls back to its own numerical gradient.
          auto runMigrad = [&](const ROOT::Minuit2::MnUserParameters &start) {
            if (configure().paramMask & Config::USE_ANALYTIC_GRADIENT) {
              configure().outStream << "Using analytic gradient for Minuit2." << std::endl;
              std::vector<double> p0 = start.Params();
              if (theFunc.EnableG2(p0))
                configure().outStream << "Seeding MIGRAD with the analytic Gauss-Newton G2." << std::endl;
              else
                configure().outStream << "Analytic G2 unavailable for this model; MIGRAD seeds numerically." << std::endl;
              ROOT::Minuit2::MnMigrad migrad(theFunc, start);
              return migrad(migradCalls);
            } else {
              const ROOT::Minuit2::FCNBase &fb = theFunc;
              ROOT::Minuit2::MnMigrad migrad(fb, start);
              return migrad(migradCalls);
            }
          };
          ROOT::Minuit2::FunctionMinimum min = runMigrad(params.GetMinuitParams());
          reportMigrad(min);
          // relvector data-covariance terms use the model frozen at the FIRST evaluation of the run,
          // i.e. at the starting parameters, so MIGRAD minimized a chi^2 whose covariance belongs to
          // the start.  Re-freeze it at the fitted parameters: the outputs (chiSquared.out, ...) then
          // match a calculation from the param.sav written below; with --refreeze-covariance[=N]
          // MIGRAD is rerun with the re-frozen covariance until the objective stops changing.
          // (12C+a_onefile 10-29-26: fits that started 15k-90k above their minimum reported a chi^2
          // 50-600 away from a mode-1 calculation of their own param.sav.)
          if (data()->HasRelativeCovariance()) {
            for (int pass = 0;; pass++) {
              const double frozen = min.Fval();
              data()->RefreezeCovariance();
              const double refrozen = theFunc(min.UserParameters().Params());
              configure().outStream << "Data covariance: relvector model re-frozen at the fitted parameters; objective "
                                    << frozen << " (model frozen at the start of this pass) -> " << refrozen << "."
                                    << std::endl;
              if (std::fabs(refrozen - frozen) < 0.5) break;
              if (pass >= configure().refreezeCovariance) {
                configure().outStream << "WARNING: the fitted parameters minimize the chi^2 with the start's covariance "
                                         "model, not with their own"
                                      << (configure().refreezeCovariance > 0 ? " (pass limit reached)" : "")
                                      << "; rerun from param.sav or use --refreeze-covariance." << std::endl;
                break;
              }
              configure().outStream << "Refitting with the re-frozen covariance (pass " << pass + 1 << " of at most "
                                    << configure().refreezeCovariance << ")..." << std::endl;
              min = runMigrad(min.UserParameters());
              reportMigrad(min);
            }
          }
          // MIGRAD's error matrix is approximate; refine it with HESSE for the band,
          // falling back to the MIGRAD matrix if HESSE gives no valid covariance.
          if (configure().paramMask & Config::CALCULATE_COVARIANCE_BAND) {
            std::vector<double> migradCov = min.UserCovariance().Data();
            configure().outStream << std::endl
                                  << "Computing accurate covariance (HESSE) for the uncertainty band..."
                                  << std::endl;
            ROOT::Minuit2::MnHesse hesse;
            try {
              hesse(theFunc, min);
            } catch (...) {}
            if (min.HasValidCovariance()) {
              captureMinuitCov(min.UserCovariance().Data());
            } else {
              configure().outStream << "HESSE did not return a valid covariance; using the MIGRAD "
                                       "error matrix (uncertainties may be approximate)."
                                    << std::endl;
              reportParametersAtLimits();
              captureMinuitCov(migradCov);
            }
          }
          if (configure().paramMask & Config::PERFORM_ERROR_ANALYSIS) {
            configure().outStream << std::endl
                                  << "Performing parameter error analysis with Up=" << configure().chiVariance << "." << std::endl;
            data()->SetErrorAnalysis(true);
            theFunc.SetErrorDef(configure().chiVariance);
            ROOT::Minuit2::MnMinos minos(theFunc, min);
            std::vector<std::pair<double, double>> errors;
            for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
              configure().outStream << "\tParameter " << i + 1 << "..." << std::endl;
              if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
                std::pair<double, double> error = minos(i);
                errors.push_back(error);
              } else
                errors.push_back(std::pair<double, double>(0., 0.));
            }

            // New output of the covariance matrix
            char filename[256];
            snprintf(filename, sizeof(filename), "%scovariance_matrix.out", configure().outputdir.c_str());
            std::ofstream out;
            out.open(filename);
            if (out) {
              // Write header information
              params.GetMinuitParams() = min.UserParameters();
              out << "Parameter List" << std::endl
                  << std::endl;
              for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
                out << std::setw(20) << params.GetMinuitParams().GetName(i)
                    << " (" << std::setw(3) << i << ")"
                    << std::scientific << std::setw(20) << params.GetMinuitParams().Value(i);
                if (params.GetMinuitParams().Parameter(i).IsFixed()) {
                  out << " Fixed" << std::endl;
                } else {
                  out << " Fitted by Minuit" << std::endl;
                }
              }

              std::cout << min.UserCovariance();

              std::vector<double> CovarianceData;
              CovarianceData = min.UserCovariance().Data();
              std::vector<int> parameterTable;
              parameterTable.reserve(params.GetMinuitParams().Params().size());

              // Header for Covariance Matrix
              out << std::endl
                  << "Covariance Matrix" << std::endl
                  << std::endl;
              out << std::setw(5) << " ";
              for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
                if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
                  out << std::setw(15) << i;
                  parameterTable.push_back(i);
                }
              }
              out << std::endl;
              int size = parameterTable.size();

              std::cout << "covariance length " << CovarianceData.size() << std::endl;

              // Covariance matrix
              for (int i = 0; i < size; i++) {
                out << std::setw(5) << parameterTable[i];
                for (int j = 0; j < size; j++) {
                  int k;
                  if (i <= j) {
                    k = ((j + 1) * j) / 2 + i;
                  }
                  if (i > j) {
                    k = ((i + 1) * i) / 2 + j;
                  }
                  std::cout << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k << std::endl;
                  out << std::setw(15) << CovarianceData[k];
                }
                out << std::endl;
              }

              // Header for Correlation Matrix
              out << std::endl
                  << "Correlation Matrix" << std::endl
                  << std::endl;
              out << std::setw(5) << " ";
              for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
                if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
                  out << std::setw(15) << i;
                }
              }
              out << std::endl;

              // Correlation matrix
              for (int i = 0; i < size; i++) {
                out << std::setw(5) << parameterTable[i];
                for (int j = 0; j < size; j++) {
                  int k;
                  if (i <= j) {
                    k = ((j + 1) * j) / 2 + i;
                  }
                  if (i > j) {
                    k = ((i + 1) * i) / 2 + j;
                  }
                  int jdiag = ((j + 2) * (j + 1)) / 2 - 1;
                  int idiag = ((i + 2) * (i + 1)) / 2 - 1;
                  std::cout << k << "  " << CovarianceData[k] << std::endl;
                  out << std::setw(15) << std::fixed << CovarianceData[k] / sqrt(fabs(CovarianceData[jdiag] * CovarianceData[idiag]));
                }
                out << std::endl;
              }

              out.flush();
              out.close();
            } else
              std::cout << "Could not save " << configure().outputdir.c_str()
                        << "covariance_matrix.out file." << std::endl;

            // ECS this line added to set the azure variables to the minimised fit parameters
            params.GetMinuitParams() = min.UserParameters();
            params.WriteParameterErrors(errors, configure());
          }
          params.GetMinuitParams() = min.UserParameters();
          data()->SyncShiftGroupParams(params.GetMinuitParams());
          data()->SyncNormGroupParams(params.GetMinuitParams());
          params.WriteUserParameters(configure(), true);
        }  // end MIGRAD (else of Levenberg-Marquardt selector)
#ifdef USE_NLOPT
      }
#endif
    } else {
      if (configure().paramMask & Config::USE_AMATRIX)
        configure().outStream << "Performing A-Matrix Calculation..." << std::endl;
      else
        configure().outStream << "Performing R-Matrix Calculation..." << std::endl;
    }
    data()->SetFit(false);
    data()->SetErrorAnalysis(false);
    double chiSquared = theFunc(params.GetMinuitParams().Params());
    if (configure().paramMask & Config::CALCULATE_WITH_DATA) {
      configure().outStream << std::endl
                            << std::endl;
      for (ESegmentIterator segment = data()->GetSegments().begin();
           segment < data()->GetSegments().end(); segment++) {
        configure().outStream << "Segment #"
                              << segment->GetSegmentKey()
                              << " Chi-Squared/N: "
                              << segment->GetSegmentChiSquared() / segment->NumPoints()
                              << std::endl;
        // With component segments, no need to skip - total capture is handled within the segment
      }
      configure().outStream << "Total Chi-Squared: "
                            << chiSquared << std::endl
                            << std::endl;
    }

    // Write Output Files
    configure().outStream << "Writing output files..." << std::endl;
    // Band covariance: this run's fit if available, else one saved by a prior fit.
    BandData bandData;
    const BandData *bandPtr = nullptr;
    if (configure().paramMask & Config::CALCULATE_COVARIANCE_BAND) {
      // Optionally inflate to reduced chi^2 = 1 (Birge ratio), then persist.
      if (haveBandCov) {
        if ((configure().paramMask & Config::SCALE_COVARIANCE_BY_CHI2) &&
            (configure().paramMask & Config::CALCULATE_WITH_DATA)) {
          double chi2data = 0.;
          int nData = 0;
          for (int s = 1; s <= data()->NumSegments(); s++) {
            ESegment *seg = data()->GetSegment(s);
            if (!seg) continue;
            chi2data += seg->GetSegmentChiSquared();
            nData += seg->NumPoints();
          }
          // Degrees of freedom use the total number of free fit parameters
          // (including norms/shifts), even though the covariance matrix itself
          // holds only the R-matrix sub-block (bandCov.size()).
          int nFree = 0;
          const int nMnAll = params.GetMinuitParams().Params().size();
          for (int i = 0; i < nMnAll; i++)
            if (!params.GetMinuitParams().Parameter(i).IsFixed()) nFree++;
          int ndf = nData - nFree;
          double redChi2 = (ndf > 0) ? chi2data / ndf : 1.0;
          double factor = (redChi2 > 1.0) ? redChi2 : 1.0;  // inflate only
          configure().outStream << "Calculating covariance scaling: chi-Squared = " << chi2data
                                << ", N = " << nData << ", free params = " << nFree
                                << ", nu = " << ndf << ", reduced chi-squared = " << redChi2 << "." << std::endl;
          if (factor > 1.0) {
            const int nMat = bandCov.size();
            for (int a = 0; a < nMat; a++)
              for (int b = 0; b < nMat; b++) bandCov.M[a][b] *= factor;
            bandCov.scale = factor;
            configure().outStream << "Covariance scaled by " << factor << "." << std::endl;
          } else {
            configure().outStream << "Reduced chi-squared <= 1; covariance left unscaled."
                                  << std::endl;
          }
        }
        // Name the columns ("#columns" header, which LoadBandCovariance and numpy.loadtxt skip): the free level
        // energies and reduced-width amplitudes in parameter order, as param.sav holds them (pyazr.prior_from_fit).
        {
          const ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
          std::vector<std::string> names;
          for (unsigned int i = 0; i < mp.Params().size(); i++) {
            const std::string &nm = mp.GetName(i);
            if (!mp.Parameter(i).IsFixed() && (nm.rfind("energy_", 0) == 0 || nm.rfind("width_", 0) == 0))
              names.push_back(nm);
          }
          if ((int)names.size() == bandCov.size()) bandCov.names = names;
        }
        if (SaveBandCovariance(covPath, bandCov))
          configure().outStream << "Saved parameter covariance for cross-section bands." << std::endl;
      }
      BandCovariance cov;
      bool fitRun = configure().paramMask & Config::PERFORM_FIT;  // fit or MINOS
      if (haveBandCov)
        cov = bandCov;
      else if (fitRun) {
        // The fit ran but produced no usable covariance (already reported above);
        // don't fall back to a stale saved file.
        configure().outStream << "Uncertainty band skipped: this fit produced no usable covariance."
                              << std::endl;
      } else if (!LoadBandCovariance(covPath, cov)) {
        // Extrapolation / no-fit run: the band needs a covariance from a prior fit.
        configure().outStream << "No saved parameter covariance found (" << covPath
                              << "); writing the normal output without an uncertainty band. "
                                 "To get a band here, first run a fit (or MINOS) with the "
                                 "\"Uncertainty band\" option enabled, then rerun this extrapolation."
                              << std::endl;
      }
      if (!cov.empty()) {
        if (BuildBandData(compound(), data(), configure(), cov, bandData)) {
          bandPtr = &bandData;
          configure().outStream << "Writing cross-section uncertainty bands..." << std::endl;
        } else {
          configure().outStream << "Cross-section bands unavailable: the analytic model "
                                   "sensitivities are not supported for this configuration."
                                << std::endl;
        }
        // BuildBandData computes the per-point sensitivities on the real
        // compound/data objects (its gradient map is keyed by the real EPoint*),
        // and that forward+adjoint solve overwrites each point's stored best-fit
        // cross section.  Re-run the best-fit forward pass so the plotted and
        // written segments are the best fit, not the last gradient probe; the
        // band itself lives in bandData and is unaffected.
        theFunc(params.GetMinuitParams().Params());
      }
    }
    data()->WriteOutputFiles(configure(), false, bandPtr);
    // s-wave scattering lengths of every particle pair at the parameters just written, and the
    // measured ones of the <scatteringLengths> block: the table goes to scatteringLengths.out only,
    // the run log (GUI execution screen) names the file.
    try {
      std::string report = ScatteringLength::Report(compound(), configure(), compound()->ScatteringLengthConstraints());
      std::ofstream slOut((configure().outputdir + "scatteringLengths.out").c_str());
      slOut << report;
      configure().outStream << "Scattering lengths written to " << configure().outputdir << "scatteringLengths.out"
                            << std::endl;
    } catch (const std::exception &e) {
      configure().outStream << "Scattering lengths not calculated: " << e.what() << std::endl;
    }
  } else {
    // Calculate Reaction Rate
    //  This uses the adaptive integration routines of GSL.  As the energy stepsize is
    //  unknown until integration, AZURE is called not in segment control mode but as
    //  an energy dependent function.  As every data point must be reinitialized (new energy
    //  dependent terms calculated) the routine is slow.  This should be a tradeoff
    //  for good accuracy.
    configure().outStream << "Performing reaction rate calculation..." << std::endl;
    ReactionRate reactionRate(compound(), params.GetMinuitParams().Params(), configure(),
                              configure().rateParams.entrancePair, configure().rateParams.exitPair);
    // Warn once if a level has a radiative width comparable to its particle width,
    // where the R-Matrix treatment of capture breaks down.
    compound()->CheckRadiativeWidths(configure(), params.GetMinuitParams().Params());
    if (configure().paramMask & Config::USE_BRUNE_FORMALISM)
      compound()->CalcShiftFunctions(configure());
    if (configure().rateParams.useFile)
      reactionRate.CalculateFileRates();
    else
      reactionRate.CalculateRates();
    reactionRate.WriteRates();
  }

  configure().outStream << "Performing final parameter transformation..." << std::endl;
  try {
    compound()->TransformOut(configure());
  } catch (GSLException e) {
    configure().outStream << e.what() << std::endl;
    configure().outStream << "Problem with output transformation.  Aborting."
                          << std::endl;
    return -1;
  }
  compound()->PrintTransformParams(configure());

  // Cleanup shared EC amplitude cache
  if (configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
    if (g_ecAmplitudeCache) {
      // g_ecAmplitudeCache->PrintStats();
    }
    CleanupECAmplitudeCache();
  }

  return 0;
}