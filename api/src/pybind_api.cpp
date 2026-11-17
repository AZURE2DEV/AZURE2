// pybind11 bindings for the AZURE2 R-matrix engine.
//
// This replaces the old TCP/IP socket server ("--use-api <port>").  A Python
// session is a real AZUREAPI object living in the interpreter process: no
// subprocess, no sockets, no port bookkeeping.
//
// The bound class is ``_azure2.Session``.  One Session owns one Config and one
// AZUREAPI; constructing it reads the .azr file and runs the (expensive)
// initialization, so an invalid model fails at construction time.
//
// Several Sessions may be live at once, also on several Python threads at the
// same time: every value a calculation keeps between calls belongs to the
// session's Config (its scratch, its external-capture amplitudes, its notices)
// or to its AZUREAPI, and the engine reads nothing of another session
// (docs/upgrade/evaluation.md; tests/pyazr/session_isolation_test.py).  The
// engine calls release the GIL.  One Session is not reentrant: drive each from
// one thread at a time.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>

#include <cstdio>
#include <fstream>
#include <sstream>
#include "azure2/io/ProjectFormat.h"
#include "azure2/io/Document.h"
#include "AZUREAPI.h"
#include "CompoundSystem.h"
#include "InferenceTypes.h"
#include "CommandLineValues.h"
#include "Config.h"
#include "GSLException.h"
#include "ECAmplitudeCache.h"
#include <gsl/gsl_errno.h>
#include "pybind_notices.h"
#include "ProjectOptions.h"
#include "ReactionRate.h"
#include "CompoundSystemChecks.h"
#include <optional>
#include "azure2/io/ConfigSection.h"
#include "pybind_engine_call.h"
#include "pybind_arrays.h"

namespace py = pybind11;

// Raised when the .azr cannot be loaded or the model will not initialize.
class AZURE2Error : public std::runtime_error {
 public:
  explicit AZURE2Error(const std::string &what) :
    std::runtime_error(what) {}
};

// A free-parameter vector the engine cannot take (wrong length, NaN, infinity): pyazr.ParameterVectorError, which is
// an AZURE2Error too.
class ParameterVectorError : public AZURE2Error {
 public:
  explicit ParameterVectorError(const std::string &what) : AZURE2Error(what) {}
};

// pybind_records.cpp: the parameter and pair records with their named fields
py::list ParameterRecordDicts(AZUREAPI &api);
py::list PairRecordDicts(AZUREAPI &api);
// pybind_evaluation.cpp: one evaluation read as the writer reads it, the objective's rows and log-probability
py::list PointOutputFields();
py::dict EvaluationDict(AZUREAPI &api, const vector_r &free);
py::dict ObjectiveRowsDict(AZUREAPI &api, const vector_r &free, bool withJacobian,
                           const std::vector<int> &skip);
double LogProbability(AZUREAPI &api, const vector_r &free);
py::tuple SamplerLogProbability(AZUREAPI &api, const vector_r &x0, const vector_r &x);
// pybind_fit.cpp: the engine's fitters on the session (fit_with_engine), StopRequest, fitters()
py::dict FitWithEngineDict(AZUREAPI &api, const std::string &method, bool minos, const std::vector<double> &x0,
                           unsigned int maxCalls, const py::object &progress, const py::object &stop);
void BindFit(py::module_ &m);
py::dict SampleWithEngineDict(AZUREAPI &api, const py::kwargs &options);  // pybind_sampler.cpp: the engine's sampler
std::string ReadProjectText(const std::string &path);  // pybind_project.cpp: the bytes of a project file
void BindDataRows(py::module_ &m);
std::vector<std::string> OwnerFieldsOf(ParamKind kind);
// pybind_potentials.cpp: the view of a session's hybrid potentials
pybind11::tuple SessionPotentialView(const Config &config, int pair);
// pybind_thm_reports.cpp: the THM reports of a session as dicts
py::list ThmExperimentsList(AZUREAPI &api);
py::dict ThmLineshapeDict(const ThmLineshapeReport &r);
py::dict ThmVertexDict(const ThmVertexReport &r);
py::dict ThmDistortionDict(const ThmDistortionReport &r);

namespace {

// What a session is given besides the project file: the options that are not project state (pybind_describe.cpp
// says why for each), and the options of the calculation its call named (keyword -> value; the rows of
// io::ConfigKeys(), docs/upgrade/project_format3.md section 3.12).  An option left out is not given: the project's
// <config> key, else the default, decides (ApplyProjectOptions, when the session reads the project).
struct RuntimeOptions {
  bool data_mode = true;            // CALCULATE_WITH_DATA
  std::string coulomb_table;        // Config::coulombTablePath ("" = the default file)
  std::string base_directory;       // Config::baseDirectory: the project's relative paths ("" = the project file's directory)
  std::map<std::string, bool> options;  // the options of the calculation the call gave
};

void apply_options(Config &config, const RuntimeOptions &opt) {
  if (opt.data_mode)
    config.paramMask |= Config::CALCULATE_WITH_DATA;
  else
    config.paramMask &= ~Config::CALCULATE_WITH_DATA;
  config.baseDirectory = opt.base_directory;
  // A relative table file is the session's directory's, as it was when pyazr entered that directory per call.
  config.coulombTablePath = opt.coulomb_table.empty() ? opt.coulomb_table : config.OpenPath(opt.coulomb_table);
  for (const auto &given : opt.options)
    if (!azure2::io::ConfigOptionNamed(given.first)) throw py::value_error(NotAnOptionText(given.first));
  // In the table's order (use_brune before use_park: use_park=True with use_brune=False gives both on, as before)
  for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys()) {
    const auto given = k.option ? opt.options.find(k.name) : opt.options.end();
    if (given != opt.options.end()) GiveOption(config, k.name, given->second, OptionSource::kKeyword);
  }
}

}  // namespace

void AddModelData(py::dict &terms, AZUREAPI &api, const vector_r &free);  // pybind_model_data.cpp
class Session {
 public:
  // rateEntrance != 0: a reaction-rate session (the rate's pair keys; exit -1 the total radiative capture), set up as
  // the program's rate command sets up its run (RunPreparation.cpp, kRate): the levels, no data.
  Session(const std::string &configfile, const RuntimeOptions &opt, int rateEntrance = 0, int rateExit = 0) :
    config_(nullptr),
    api_(nullptr) {
    // AZURE_REFINE_QUADRATURE: a bad value raises (ValueError), before any allocation
    const int refine = azure2::RefineQuadratureFromEnvironment();
    config_ = new Config(noticeStream_);  // standard output, but the notices of HeldNotices() (pybind_notices.h)
    config_->refineQuadrature = refine;
    config_->configfile = configfile;
    config_->paramMask |= Config::USE_API;  // silence the "Calculating ..." chatter
    apply_options(*config_, opt);
    if (rateEntrance != 0) {
      config_->paramMask &= ~Config::CALCULATE_WITH_DATA;
      config_->paramMask |= Config::CALCULATE_REACTION_RATE;
      config_->rateParams.entrancePair = rateEntrance;
      config_->rateParams.exitPair = rateExit;
    }


    // The binding is the API: mirror the checks main() performs for a headless run.
    int configStatus = config_->ReadConfigFile();
    if (configStatus == Config::kConfigUnreadable) {
      fail("Could not open " + configfile + " or no <config> block found.");
    } else if (configStatus == Config::kConfigNotFormat3) {
      fail(azure2::io::NotFormat3Message(ReadProjectText(configfile), configfile));
    } else if (configStatus == Config::kConfigRefused && !config_->refusedOptions.empty()) {
      fail(config_->refusedOptions + configfile + ": the options of this call contradict each other or the project's <config> keys.");
    } else if (configStatus == Config::kConfigRefused) {
      fail(configfile + " (format 3) is refused (see the ERROR lines above; AZURE2 --validate lists them).");
    } else if (configStatus < 0) {
      fail("Malformed <thm> section in " + configfile + " (see the ERROR above).");
    }
    // R-matrix route + ignore_externals (the default): A-matrix with a warning, as the CLI and the GUI.  Before
    // 2026-11-01 pyazr ran the R-matrix route and silently kept the internal-capture term ignore_externals drops
    // (13N capture 3e-5 off the A-matrix route).
    config_->ResolveMatrixRoute();

    InitializeECAmplitudeCache(*config_);

    api_ = new AZUREAPI(*config_);
    if (api_->Initialize() != 0) {
      std::string msg = "AZURE2 could not initialize " + configfile +
          "; check the output above.";
      fail(msg);
    }
  }

  ~Session() {
    if (api_ != nullptr) delete api_;
    if (config_ != nullptr) delete config_;
  }

  Session(const Session &) = delete;
  Session &operator=(const Session &) = delete;

  [[noreturn]] void fail(const std::string &msg) {
    if (api_ != nullptr) {
      delete api_;
      api_ = nullptr;
    }
    if (config_ != nullptr) {
      delete config_;
      config_ = nullptr;
    }
    throw AZURE2Error(msg);
  }

  // -- output files and external capture --------------------------------------

  bool write_output_files() {
    return api_->WriteOutputFiles();
  }

  bool calculate_external_capture() {
    return api_->CalculateExternalCapture();
  }

  // -- options ---------------------------------------------------------------

  // The options of the calculation in effect: {name: (value, source)}, source "default", "project" or "keyword".
  py::dict option_sources() const {
    py::dict out;
    for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys())
      if (k.option) out[k.name] = py::make_tuple(OptionIsOn(*config_, k.name), OptionSourceName(OptionSourceOf(*config_, k.name)));
    return out;
  }

  // -- hybrid nuclear potential, per particle pair ---------------------------
  //
  // The session's Config holds one setting per pair (the project's potential records) and a default for the pairs
  // that name none; CompoundSystem::Fill gives each pair its potential when the session is made.
  // (enabled, type, V0, R, a, r0, has_own_setting, coulomb_radius, spin_orbit, matching_radius, step) for the given pair, resolved against the default.
  py::tuple get_potential(int pair) { return SessionPotentialView(*config_, pair); }

  // -- data bookkeeping ------------------------------------------------------

  int update_data() {
    return api_->UpdateData();
  }

  // -- parameter metadata ----------------------------------------------------

  /// One dict per parameter, in params_names() order: its record's named fields (pybind_records.cpp).
  py::list parameter_records() const { return ParameterRecordDicts(*api_); }
  /// One dict per particle pair, in pair-number order (pybind_records.cpp).
  py::list pair_records() const { return PairRecordDicts(*api_); }

  // -- parameter vectors -----------------------------------------------------

  py::array_t<double> params_values() const { return to_array(api_->params_values()); }
  py::array_t<double> params_values_rwa() const { return to_array(api_->params_values_rwa()); }
  py::array_t<double> params_all_rwa() const { return to_array(api_->params_all_rwa()); }
  std::string params_names(int i) const { return api_->params_names(i); }
  std::vector<bool> params_fixed() const { return api_->params_fixed(); }

  // -- data arrays -----------------------------------------------------------

  py::array_t<double> data_energies(int i) const { return to_array(api_->data_energies(i)); }
  py::array_t<double> data_angles(int i) const { return to_array(api_->data_angles(i)); }
  py::array_t<double> data_segments(int i) const { return to_array(api_->data_segments(i)); }
  py::array_t<double> data_segments_errors(int i) const { return to_array(api_->data_segments_errors(i)); }
  py::array_t<double> data_excitation_energies(int i) const { return to_array(api_->data_excitation_energies(i)); }
  py::array_t<double> data_conv(int i) const { return to_array(api_->data_conv(i)); }
  py::array_t<double> norms() {
    api_->UpdateNorms();
    return to_array(api_->norms());
  }

  // -- calculated arrays -----------------------------------------------------

  py::array_t<double> calculated_segments(int i) const { return to_array(api_->calculated_segments(i)); }
  py::array_t<double> calculated_segments_e1(int i) const { return to_array(api_->calculated_segments_e1(i)); }
  py::array_t<double> calculated_segments_e2(int i) const { return to_array(api_->calculated_segments_e2(i)); }
  py::array_t<double> calculated_energies(int i) const { return to_array(api_->calculated_energies(i)); }
  py::array_t<double> calculated_angles(int i) const { return to_array(api_->calculated_angles(i)); }
  py::array_t<double> calculated_angular_dists(int i) const { return to_array(api_->calculated_angular_dists(i)); }
  py::array_t<double> calculated_excitation_energies(int i) const { return to_array(api_->calculated_excitation_energies(i)); }
  py::array_t<double> calculated_conv(int i) const { return to_array(api_->calculated_conv(i)); }

  // -- argument checks ---------------------------------------------------------
  // The engine copies p[k] for every free parameter without looking at p's
  // length.  A short vector therefore read past its end and left the norms and
  // shifts at whatever the engine last held: cross sections and S-factors still
  // right (they use no norm), chi2 garbage, no error (12C+alpha, 2026-10-08).
  // A NaN or an infinity reached the Coulomb functions and failed there, naming no parameter (or gave a NaN chi2).
  void expect_free(const vector_r &v, const char *what) const {
    std::vector<bool> fixed = api_->params_fixed();
    size_t nfree = 0;
    for (bool f : fixed) if (!f) nfree++;
    if (v.size() != nfree)
      throw ParameterVectorError(std::string(what) + ": expected the free parameter vector (" +
                                 std::to_string(nfree) + " values, in the order of the fit vector), got " +
                                 std::to_string(v.size()) + ".  " + kVectorHint);
    std::vector<int> slot;  // free index -> index among all parameters
    for (size_t k = 0; k < fixed.size(); ++k) if (!fixed[k]) slot.push_back((int)k);
    expect_finite(v, slot, what, "free parameter vector");
  }
  void expect_all(const vector_r &v, const char *what) const {
    size_t nall = api_->params_fixed().size();
    if (v.size() != nall)
      throw ParameterVectorError(std::string(what) + ": expected every parameter (" + std::to_string(nall) +
                                 " values, fixed ones included), got " + std::to_string(v.size()) + ".  " +
                                 "Session._engine.params_all_rwa() is such a vector.");
    std::vector<int> slot(nall);
    for (size_t k = 0; k < nall; ++k) slot[k] = (int)k;
    expect_finite(v, slot, what, "parameter vector");
  }
  static constexpr const char *kVectorHint =
      "Session.parameters.x0 is such a vector (the start values); Session.parameters.free names its entries in order.";
  // The first entry that is not finite, by its index and the parameter's name.
  void expect_finite(const vector_r &v, const std::vector<int> &slot, const char *what, const char *vector) const {
    size_t bad = 0, count = 0;
    for (size_t i = 0; i < v.size(); ++i)
      if (!std::isfinite(v[i]) && count++ == 0) bad = i;
    if (count == 0) return;
    std::ostringstream text;
    text << what << ": entry " << bad << " of the " << vector << " (" << api_->params_names(slot[bad]) << ") is "
         << v[bad] << (count > 1 ? " (" + std::to_string(count) + " entries are not finite)" : std::string())
         << "; every entry must be a finite number.  " << kVectorHint;
    throw ParameterVectorError(text.str());
  }

  // -- segment updates (return the number of segments) ------------------------

  int update_segments(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "update_segments");
    return api_->UpdateSegments(v);
  }
  int update_segments_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "update_segments_rwa");
    return api_->UpdateSegmentsRWA(v);
  }
  void set_params_all_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_all(v, "set_params_all_rwa");
    api_->SetParamsAllRWA(v);
  }
  int update_segments_all_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_all(v, "update_segments_all_rwa");
    return api_->UpdateSegmentsAllRWA(v);
  }

  // -- transforms ------------------------------------------------------------

  py::array_t<double> transform_all_rwa(const py::array_t<double, py::array::forcecast> &p, bool include_fixed) {
    vector_r v = to_vector(p), out;
    expect_all(v, "transform_all_rwa");
    {
      EngineCallScope engine;
      out = api_->TransformAllRWAParameters(v, include_fixed);
    }
    return to_array(out);
  }

  // -- chi-squared and derivatives -------------------------------------------

  double calculate_chi2_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "calculate_chi2_rwa");
    return api_->CalculateChi2RWA(v);
  }
  double calculate_chi2_physical(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "calculate_chi2_physical");
    return api_->CalculateChi2Physical(v);
  }
  py::array_t<double> scattering_lengths_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "scattering_lengths_rwa");
    return to_array(api_->CalculateScatteringLengthsRWA(v));
  }
  py::array_t<double> park_norms(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "park_norms");
    return to_array(api_->ParkNorms(v));
  }
  py::array_t<double> observed_levels_all_rwa(const py::array_t<double, py::array::forcecast> &p) { return to_array(api_->ObservedLevelRowsAllRWA(to_vector(p))); }
  py::array_t<double> prior_covariance_residuals(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "prior_covariance_residuals");
    return to_array(api_->PriorCovarianceResidualsRWA(v));
  }
  std::vector<std::string> prior_covariance_names() {
    return api_->PriorCovarianceNames();
  }
  int refreeze_covariance() {
    return api_->RefreezeCovariance();
  }

  // -- the reaction rate (a session made with rate_entrance) ----------------------------------------------------
  // N_A<sigma v> at every parameter value of `p` (all parameters, the params_all_rwa() layout), as the program's rate
  // command computes it (AZUREMain.cpp RunReactionRate): on the temperatures min + i step (i = 0 .. (max - min) / step,
  // sorted) when `temperature_file` is empty, else on the file's temperatures in its order.  (T9, rate) pairs.
  std::vector<std::pair<double, double>> reaction_rate(const py::array_t<double, py::array::forcecast> &p,
                                                       double minTemp, double maxTemp, double tempStep,
                                                       const std::string &temperatureFile) {
    if (!(config_->paramMask & Config::CALCULATE_REACTION_RATE))
      throw py::value_error("reaction_rate: the session was not made for a rate (Session(..., rate_entrance=))");
    vector_r v = to_vector(p);
    expect_all(v, "reaction_rate");
    std::vector<std::pair<double, double>> out;
    {
      EngineCallScope engine;
      Config &configure = *config_;
      configure.rateParams.useFile = !temperatureFile.empty();
      configure.rateParams.temperatureFile = temperatureFile;
      configure.rateParams.minTemp = minTemp;
      configure.rateParams.maxTemp = maxTemp;
      configure.rateParams.tempStep = tempStep;
      ReactionRate rate(api_->compound(), api_->state(), api_->inputs(), v, configure, configure.rateParams.entrancePair,
                        configure.rateParams.exitPair);
      if (configure.rateParams.useFile)
        rate.CalculateFileRates();
      else
        rate.CalculateRates();
      for (const RateData &r : rate.Rates()) out.emplace_back(r.temperature, r.rate);
    }
    return out;
  }
  py::array_t<double> calculate_chi2_grad_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p), out;
    expect_free(v, "calculate_chi2_grad_rwa");
    {
      EngineCallScope engine;
      out = api_->CalculateChi2GradRWA(v);
    }
    return to_array(out);
  }
  // The fit's objective (docs/upgrade/objective.md): its value, its gradient, and its terms (names, values, and the
  // normalization and shift priors per segment, index 0 unused).
  double calculate_objective_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "calculate_objective_rwa");
    return api_->CalculateObjectiveRWA(v);
  }
  py::array_t<double> calculate_objective_grad_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p), out;
    expect_free(v, "calculate_objective_grad_rwa");
    {
      EngineCallScope engine;
      out = api_->CalculateObjectiveGradRWA(v);
    }
    return to_array(out);
  }
  py::dict objective_terms_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "objective_terms_rwa");
    std::vector<double> terms;
    ObjectiveContext::SegmentPriors seg;
    const double total = api_->CalculateObjectiveRWA(v, &terms, &seg);
    py::dict d;
    d["total"] = total;
    d["names"] = api_->ObjectiveTermNames();
    d["values"] = to_array(terms);
    d["norm"] = to_array(seg.norm);
    d["shift"] = to_array(seg.shift);
    d["shift_sqrt"] = to_array(seg.shiftSqrt);
    d["shift_fn"] = to_array(seg.shiftFn);
    AddModelData(d, *api_, v);  // the <modelData> data, one record each
    return d;
  }
  // pyazr's Evaluation and Objective (pybind_evaluation.cpp; docs/upgrade/pyazr3.md section 7).
  py::dict evaluation_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "evaluation_rwa");
    return EvaluationDict(*api_, v);
  }
  py::dict objective_rows_rwa(const py::array_t<double, py::array::forcecast> &p, bool jacobian,
                              const std::vector<int> &skip) {
    vector_r v = to_vector(p);
    expect_free(v, "objective_rows_rwa");
    for (int c : skip)
      if (c < 0 || c >= (int)v.size())
        throw AZURE2Error("objective_rows_rwa: skip column " + std::to_string(c) + " is outside the " +
                          std::to_string(v.size()) + " free parameters.");
    return ObjectiveRowsDict(*api_, v, jacobian, skip);
  }
  std::vector<std::string> objective_term_names() const { return api_->ObjectiveTermNames(); }
  py::dict fit_with_engine(const std::string &method, bool minos, const std::vector<double> &x0, unsigned int maxCalls,
                           const py::object &progress, const py::object &stop) {
    return FitWithEngineDict(*api_, method, minos, x0, maxCalls, progress, stop);
  }
  py::dict sample_with_engine(const py::kwargs &options) { return SampleWithEngineDict(*api_, options); }
  double log_prob_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p);
    expect_free(v, "log_prob_rwa");
    return LogProbability(*api_, v);
  }
  py::tuple sampler_log_probability(const py::array_t<double, py::array::forcecast> &x0,
                                    const py::array_t<double, py::array::forcecast> &p) {
    vector_r s = to_vector(x0), v = to_vector(p);
    expect_free(s, "sampler_log_probability");
    expect_free(v, "sampler_log_probability");
    return SamplerLogProbability(*api_, s, v);
  }
  py::array_t<double> calculate_residual_jacobian_rwa(const py::array_t<double, py::array::forcecast> &p,
                                                      const std::vector<int> &skip) {
    vector_r v = to_vector(p), out;
    expect_free(v, "calculate_residual_jacobian_rwa");
    for (int c : skip)
      if (c < 0 || c >= (int)v.size())
        throw AZURE2Error("calculate_residual_jacobian_rwa: skip column " + std::to_string(c) + " is outside the " +
                          std::to_string(v.size()) + " free parameters.");
    {
      EngineCallScope engine;
      out = api_->CalculateResidualJacobianRWA(v, skip);
    }
    return to_array(out);
  }
  py::array_t<double> calculate_residuals_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    {
      EngineCallScope engine;
      out = api_->CalculateResidualsRWA(v);
    }
    return to_array(out);
  }
  // The THM experiments as the last chi-squared evaluation profiled them.
  py::list thm_experiments() { return ThmExperimentsList(*api_); }
  // The Coulomb line shape of a THM experiment at the current parameters.
  py::dict thm_lineshape(const std::string &name, py::array_t<double, py::array::forcecast> e) {
    ThmLineshapeReport r;
    std::string why;
    if (!api_->GetThmLineshape(name, to_vector(e), r, why)) throw AZURE2Error(why);
    return ThmLineshapeDict(r);
  }
  // The window-averaged THM entrance vertex of an experiment at the current parameters.
  py::dict thm_vertex(const std::string &name, py::array_t<double, py::array::forcecast> e) {
    ThmVertexReport r;
    std::string why;
    const vector_r energies = to_vector(e);
    bool ok;
    {
      // The DW vertex and R(E) can take seconds: without the GIL.
      EngineCallScope engine;
      ok = api_->GetThmVertex(name, energies, r, why);
    }
    if (!ok) throw AZURE2Error(why);
    return ThmVertexDict(r);
  }
  // The distortion factor R(E) of a THM experiment (distortion=...).
  py::dict thm_distortion(const std::string &name, py::array_t<double, py::array::forcecast> e) {
    ThmDistortionReport r;
    std::string why;
    const vector_r energies = to_vector(e);
    bool ok;
    {
      EngineCallScope engine;
      ok = api_->GetThmDistortion(name, energies, r, why);
    }
    if (!ok) throw AZURE2Error(why);
    return ThmDistortionDict(r);
  }
  py::array_t<double> calculate_model_gradients_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p), out;
    expect_free(v, "calculate_model_gradients_rwa");
    {
      EngineCallScope engine;
      out = api_->CalculateModelGradientsRWA(v);
    }
    return to_array(out);
  }

  // -- external region and caches --------------------------------------------

  py::array_t<double> coulomb_functions(const py::array_t<double, py::array::forcecast> &request) {
    vector_r v = to_vector(request), out;
    {
      EngineCallScope engine;
      out = api_->GetCoulombFunctions(v);
    }
    return to_array(out);
  }
  py::array_t<double> ec_integrals(const py::array_t<double, py::array::forcecast> &request) {
    vector_r v = to_vector(request), out;
    {
      EngineCallScope engine;
      out = api_->GetECIntegrals(v);
    }
    return to_array(out);
  }
  py::array_t<double> cache_stats() const {
    vector_r out;
    {
      EngineCallScope engine;
      out = api_->GetCacheStats();
    }
    return to_array(out);
  }

 private:
  Config *config_;
  AZUREAPI *api_;
  azure2_binding::NoticeBuffer notices_{std::cout};
  std::ostream noticeStream_{&notices_};

 public:
  /// The notices held back since the last call (pybind_notices.h): one dict each with the engine's line ("text"), the
  /// file it concerns ("file", a path) and the consequence ("consequence").
  py::list take_notices() {
    py::list out;
    for (size_t i : notices_.Take()) {
      const azure2_binding::HeldNotice &n = azure2_binding::HeldNotices()[i];
      py::dict d;
      d["text"] = n.text;
      d["file"] = NoticeFile(n.file);
      d["consequence"] = n.consequence;
      out.append(d);
    }
    return out;
  }

 private:
  // The path a notice's file key names.
  std::string NoticeFile(const std::string &key) const {
    if (config_ == nullptr) return "";
    if (key == "ec_amplitudes")  // WorkingCopy::CalculateECAmplitudes
      return config_->OpenPath(config_->outputdir +
                               ((config_->paramMask & Config::CALCULATE_WITH_DATA) ? "intEC.dat" : "intEC.extrap"));
    // DataOutputFiles.cpp
    if (key == "normalizations") return config_->OpenPath(config_->outputdir + "normalizations.out");
    if (key == "shifts") return config_->OpenPath(config_->outputdir + "shifts.out");
    return "";
  }
};

void BindPhysics(py::module_ &m);  // pybind_physics.cpp: the free functions of physics/
void BindThm(py::module_ &m);      // pybind_thm.cpp: <thm> blocks as the engine reads them
void BindSchema(py::module_ &m);   // pybind_schema.cpp: the table of project format 3
void BindConvert(py::module_ &m);  // pybind_convert.cpp: format-3 validation and the converter old -> format 3
void BindProject(py::module_ &m);  // pybind_project.cpp: the project file as a document
void BindDescription(py::module_ &m);  // pybind_describe.cpp: what pyazr.describe() and pyazr.kinds() read
void BindRecordFields(py::module_ &m);  // pybind_records.cpp: the fields of the parameter and pair records
void BindPrecision(py::module_ &m);  // pybind_precision.cpp: the x87 control words (a diagnostic, a test hook)

namespace {

// pyazr's error classes (pyazr/errors.py) for the engine's refusals: AZURE2Error is pyazr.AZURE2Error (also named
// pyazr.EngineError), ParameterVectorError pyazr.ParameterVectorError, so that the class a C++ refusal raises and the
// class a Python check raises are one.  The module is part of the package pyazr (its __name__ "pyazr._azure2"); the
// classes are looked up there when one is raised.
py::object ErrorClass(const char *name) {
  return py::module_::import("pyazr.errors").attr(name);
}

void BindErrors(py::module_ &m) {
  m.attr("AZURE2Error") = ErrorClass("AZURE2Error");
  m.attr("ParameterVectorError") = ErrorClass("ParameterVectorError");
  py::register_local_exception_translator([](std::exception_ptr p) {
    try {
      if (p) std::rethrow_exception(p);
    } catch (const ParameterVectorError &e) {
      PyErr_SetString(ErrorClass("ParameterVectorError").ptr(), e.what());
    } catch (const AZURE2Error &e) {
      PyErr_SetString(ErrorClass("AZURE2Error").ptr(), e.what());
    }
  });
}

}  // namespace

PYBIND11_MODULE(_azure2, m) {
  m.doc() = "In-process bindings for the AZURE2 R-matrix engine.";

  // AZURE2.cpp installs this for the executable; the module has to do it for
  // itself, and it matters more here.  GSL's default handler calls abort(), so
  // an integral that will not converge -- which the hybrid nuclear potential
  // can produce -- would take the whole interpreter down instead of raising.
  gsl_set_error_handler(&GSLException::GSLErrorHandler);

  BindErrors(m);
  BindPhysics(m);  // pybind_physics.cpp: the table of physics functions (physics_functions(), physics_evaluate)
  BindThm(m);  // pybind_thm.cpp: the engine's reading of a <thm> block
  BindFit(m);  // pybind_fit.cpp: the engine's fitters (fitters(), StopRequest)
  BindDataRows(m);  // pybind_data_rows.cpp: a data set's file as the engine reads it (data_rows)
  BindPrecision(m);  // pybind_precision.cpp: _x87_control_words, _set_x87_control_word

  m.def(
      "parameter_kinds",
      []() {
        py::list out;
        for (const ParamKindRow &row : ParamKindTable()) {
          py::dict d;
          d["key"] = row.key;
          d["unit"] = row.unit;
          d["owner"] = OwnerFieldsOf(row.kind);
          d["own_prior"] = row.ownPrior;
          d["segment_field"] = row.segmentField;
          d["vary_from_project"] = row.varyFromProject;
          d["doc"] = row.doc;
          out.append(d);
        }
        return out;
      },
      "The kinds of fit parameter (ParameterKind.h), one dict per kind in their order in the vector: key (the "
      "Parameter.kind of pyazr), unit, owner (the record fields that say what it belongs to, "
      "parameter_record_fields()), own_prior, segment_field, vary_from_project, doc.");
  BindRecordFields(m);  // pybind_records.cpp: parameter_record_fields, pair_record_fields
  m.def("point_output_fields", &PointOutputFields,
        "The numbers a point carries beside its model value (pybind_evaluation.cpp's table): name, unit ('model': the "
        "data set's engine unit), when (the data sets it applies to), doc, per_point ('number' or 'list').");
  BindSchema(m);  // pybind_schema.cpp: the format-3 table (schema(), as AZURE2 --schema)
  BindConvert(m);  // pybind_convert.cpp: validate_format3_text, convert_format3_text
  BindProject(m);  // pybind_project.cpp: the project file as a lossless document (_azure2.Document)
  // The project file's format (io::FormatOf), on text.
  m.def(
      "project_format", [](const std::string &text) { return azure2::io::FormatOf(azure2::io::Document::Parse(text)); },
      py::arg("text"), "The format of a project's text: 2 or 3 with format=2 or format=3 in <config>, else 1.");
  py::register_exception<GSLException>(m, "GSLError", PyExc_ArithmeticError);

  py::class_<RuntimeOptions>(m, "RuntimeOptions")
      .def(py::init<>())
      .def_readwrite("data_mode", &RuntimeOptions::data_mode)
      .def_readwrite("coulomb_table", &RuntimeOptions::coulomb_table)
      .def_readwrite("base_directory", &RuntimeOptions::base_directory)
      .def_readwrite("options", &RuntimeOptions::options,
                     "The options of the calculation the call gave, {keyword: bool}; one left out is not given.")
;

  py::class_<Session>(m, "Session")
      .def(py::init<const std::string &, const RuntimeOptions &, int, int>(), PrecisionGuard(),
           py::arg("configfile"), py::arg("options") = RuntimeOptions(), py::arg("rate_entrance") = 0,
           py::arg("rate_exit") = 0,
           "Load ``configfile`` and build the model.  Raises AZURE2Error if "
           "the file is missing or the model will not initialize.  rate_entrance (not 0) and rate_exit: a session "
           "for the reaction rate of those pairs (reaction_rate), set up as the program's rate command.")
      .def("write_output_files", &Session::write_output_files,
           py::call_guard<EngineCallScope>(),
           "Write AZUREOut_*, chiSquared.out and the rest into the output "
           "directory, from the cross sections currently on the points.")
      .def("calculate_external_capture", &Session::calculate_external_capture,
           py::call_guard<EngineCallScope>(),
           "Recompute every external-capture integral from scratch.")
      .def("reaction_rate", &Session::reaction_rate, PrecisionGuard(), py::arg("params"), py::arg("min_temp"),
           py::arg("max_temp"), py::arg("temp_step"), py::arg("temperature_file"),
           "N_A<sigma v> (cm^3 mol^-1 s^-1) of a rate session (Session(..., rate_entrance, rate_exit)) at every "
           "parameter value (params_all_rwa() layout), as the program's rate command: [(T9, rate)] on min_temp + i "
           "temp_step (GK) up to max_temp, sorted, or on the temperatures of temperature_file in its order when that is not empty.")
      .def("option_sources", &Session::option_sources,
           "The options of the calculation in effect: {name: (value, source)}, source 'default', 'project' (the "
           "<config> key) or 'keyword' (given to the session).")
      .def("get_potential", &Session::get_potential, py::arg("pair"),
           "(enabled, type, V0, R, a, r0, has_own_setting, coulomb_radius, spin_orbit [] or [V_so, R_so, a_so, on], matching_radius, step) for a pair.")
      .def("update_data", &Session::update_data, py::call_guard<EngineCallScope>(),
           "Refill the data arrays; returns the number of segments.")
      .def("parameter_records", &Session::parameter_records,
           "Every parameter as a dict of named fields (parameter_record_fields()), in params_names() order.")
      .def("pair_records", &Session::pair_records,
           "Every particle pair as a dict of named fields (pair_record_fields()), in pair-number order.")
      .def("params_values", &Session::params_values)
      .def("params_values_rwa", &Session::params_values_rwa)
      .def("params_all_rwa", &Session::params_all_rwa)
      .def("params_names", &Session::params_names, py::arg("i"))
      .def("params_fixed", &Session::params_fixed)
      .def("data_energies", &Session::data_energies, py::arg("segment"))
      .def("data_angles", &Session::data_angles, py::arg("segment"))
      .def("data_segments", &Session::data_segments, py::arg("segment"))
      .def("data_segments_errors", &Session::data_segments_errors, py::arg("segment"))
      .def("data_excitation_energies", &Session::data_excitation_energies, py::arg("segment"))
      .def("data_conv", &Session::data_conv, py::arg("segment"))
      .def("norms", &Session::norms)
      .def("calculated_segments", &Session::calculated_segments, py::arg("segment"))
      .def("calculated_segments_e1", &Session::calculated_segments_e1, py::arg("segment"))
      .def("calculated_segments_e2", &Session::calculated_segments_e2, py::arg("segment"))
      .def("calculated_energies", &Session::calculated_energies, py::arg("segment"))
      .def("calculated_angles", &Session::calculated_angles, py::arg("segment"))
      .def("calculated_angular_dists", &Session::calculated_angular_dists, py::arg("segment"))
      .def("calculated_excitation_energies", &Session::calculated_excitation_energies, py::arg("segment"))
      .def("calculated_conv", &Session::calculated_conv, py::arg("segment"))
      .def("update_segments", &Session::update_segments,
           py::call_guard<EngineCallScope>(), py::arg("params"))
      .def("update_segments_rwa", &Session::update_segments_rwa,
           py::call_guard<EngineCallScope>(), py::arg("params"))
      .def("set_params_all_rwa", &Session::set_params_all_rwa,
           py::call_guard<EngineCallScope>(), py::arg("params"),
           "Set every parameter, fixed ones included (params_all_rwa layout), "
           "without calculating.  Fixed values are what the free-vector calls use.")
      .def("update_segments_all_rwa", &Session::update_segments_all_rwa,
           py::call_guard<EngineCallScope>(), py::arg("params"))
      .def("transform_all_rwa", &Session::transform_all_rwa, PrecisionGuard(), py::arg("params"),
           py::arg("include_fixed") = false)
      .def("calculate_chi2_rwa", &Session::calculate_chi2_rwa,
           py::call_guard<EngineCallScope>(), py::arg("params"))
      .def("scattering_lengths_rwa", &Session::scattering_lengths_rwa, PrecisionGuard(), py::arg("params"),
           "s-wave scattering lengths of every particle pair at a free RWA vector: 6 numbers per channel "
           "(pair key, channel spin, parity, Re a0 fm, Im a0 fm, has-levels flag).")
      .def("park_norms", &Session::park_norms, PrecisionGuard(), py::arg("params"), "Park overlap J of every R-matrix level (1 in Brune mode).")
      .def("observed_levels_all_rwa", &Session::observed_levels_all_rwa, PrecisionGuard(), py::arg("params"))
      .def("prior_covariance_residuals", &Session::prior_covariance_residuals, PrecisionGuard(), py::arg("params"),
           "Prior parameter covariance (<config> priorCovariance=): the whitened residuals L^-1 (p - mu) at the free "
           "RWA parameters (empty without one); their squared sum is the prior's chi2 term.")
      .def("prior_covariance_names", &Session::prior_covariance_names,
           "The prior covariance's parameter names, in its row order (fixed ones marginalized); empty without one.")
      .def("refreeze_covariance", &Session::refreeze_covariance,
           "Data covariances: forget the model frozen for the relvector terms; the next chi2 evaluation "
           "freezes it at its parameters.  Returns the number of covariance groups with relvector terms.")
      .def("calculate_chi2_physical", &Session::calculate_chi2_physical,
           py::call_guard<EngineCallScope>(), py::arg("params"))
      .def("calculate_objective_rwa", &Session::calculate_objective_rwa, PrecisionGuard(), py::arg("params"),
           "The fit's objective at the free RWA parameters: data chi^2 plus every prior and constraint.")
      .def("calculate_objective_grad_rwa", &Session::calculate_objective_grad_rwa, PrecisionGuard(), py::arg("params"),
           "[objective, gradient...] at the free RWA parameters.")
      .def("objective_terms_rwa", &Session::objective_terms_rwa, PrecisionGuard(), py::arg("params"),
           "The objective's terms: total, names, values, and per-segment norm/shift priors (index 0 unused).")
      .def("calculate_chi2_grad_rwa", &Session::calculate_chi2_grad_rwa, PrecisionGuard(), py::arg("params"))
      .def("evaluation_rwa", &Session::evaluation_rwa, PrecisionGuard(), py::arg("params"),
           "One evaluation at the free vector: the objective, its terms, the data chi^2 and every data set's arrays, "
           "read as the output writer reads them (pyazr.evaluate.Evaluation).")
      .def("objective_rows_rwa", &Session::objective_rows_rwa, PrecisionGuard(), py::arg("params"), py::arg("jacobian"),
           py::arg("skip") = std::vector<int>(),
           "The objective as least-squares rows (sum r^2 = the objective) and, with jacobian, their Jacobian; the "
           "free columns listed in skip are not finite-differenced (left zero).")
      .def("fit_with_engine", &Session::fit_with_engine, PrecisionGuard(), py::arg("method"), py::arg("minos"),
           py::arg("x0"), py::arg("max_calls"), py::arg("progress"), py::arg("stop"),
           "The program's fit on the session (pyazr/session_fit.py): a row of fitters(), without the GIL.")
      .def("sample_with_engine", &Session::sample_with_engine, PrecisionGuard(),
           "The program's ensemble sampler on the session (pybind_sampler.cpp: its keywords), without the GIL.")
      .def("objective_term_names", &Session::objective_term_names,
           "The names of the objective's terms (ObjectiveTerm::Name), in the engine's order.")
      .def("log_prob_rwa", &Session::log_prob_rwa, PrecisionGuard(), py::arg("params"),
           "-1/2 x the objective, -infinity outside the ensemble sampler's support.")
      .def("sampler_log_probability", &Session::sampler_log_probability, PrecisionGuard(), py::arg("x0"),
           py::arg("params"),
           "The engine's ensemble sampler set up at x0 (fit basis): (ln P at params, its constant).  A test's tie.")
      .def("calculate_residual_jacobian_rwa", &Session::calculate_residual_jacobian_rwa, PrecisionGuard(),
           py::arg("params"), py::arg("skip") = std::vector<int>(),
           "Residuals and Jacobian; the free columns listed in skip are not finite-differenced (left zero).")
      .def("calculate_residuals_rwa", &Session::calculate_residuals_rwa, PrecisionGuard(), py::arg("params"),
           "Standardized residuals (forward pass; THM norms profiled).")
      .def("thm_experiments", &Session::thm_experiments,
           "THM experiments (<thm> experiment[...]) as the last chi-squared evaluation "
           "profiled them: list of dicts (name, segments, background, points, chi2, status, "
           "value = [norm, b0, b1, b2], covariance = 4x4 row-major).")
      .def("thm_lineshape", &Session::thm_lineshape, PrecisionGuard(), py::arg("name"), py::arg("energies"),
           "Coulomb line shape of THM experiment `name` (lineshape=on) at c.m. energies, with the "
           "level poles of the last evaluation's parameters.")
      .def("thm_vertex", &Session::thm_vertex, PrecisionGuard(), py::arg("name"), py::arg("energies"),
           "THM entrance vertex of experiment `name` at c.m. energies: the spectator-momentum "
           "window at every energy (nodes p_s, weights, T_s, theta_cm, rho: lists over E) and, per "
           "entrance channel and level, <|M_l|^2> over the window and |M_l|^2 at p_s = 0.")
      .def("thm_distortion", &Session::thm_distortion, PrecisionGuard(), py::arg("name"), py::arg("energies"),
           "Distortion factor R(E) of THM experiment `name` (distortion=coulomb|optical|table) at "
           "c.m. energies: E_sF, eta_sF, the spectator angle, |M|^2, |M_PW|^2, R directly and as the model uses it (R_model).")
      .def("calculate_model_gradients_rwa", &Session::calculate_model_gradients_rwa, PrecisionGuard(),
           py::arg("params"))
      .def("coulomb_functions", &Session::coulomb_functions, PrecisionGuard(), py::arg("request"))
      .def("ec_integrals", &Session::ec_integrals, PrecisionGuard(), py::arg("request"))
      .def("cache_stats", &Session::cache_stats)
      .def("take_notices", &Session::take_notices,
           "The engine's notices held back since the last call (pybind_notices.h), cleared: dicts with text, file, "
           "consequence.  pyazr warns with them (pyazr.EngineWarning).");
  BindDescription(m);  // pybind_describe.cpp: build_info, coulomb_table_info, session_options
}
