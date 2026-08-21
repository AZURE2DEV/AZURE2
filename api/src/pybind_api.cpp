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
// Several Sessions may be live at once, but a Session is not reentrant and the
// engine keeps process-wide state (see ConfigScope below), so drive them from
// one thread.  For parallelism give each *process* its own Session.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <cstring>
#include <iostream>
#include <string>

#include "AZUREAPI.h"
#include "Config.h"
#include "GSLException.h"
#include "CoulFuncCache.h"
#include "ECAmplitudeCache.h"
#include "NuclearPotentialManager.h"
#include <gsl/gsl_errno.h>

namespace py = pybind11;

// Defined by AZURE2.cpp for the CLI/GUI executable; here the module owns it,
// because the binding links only the core library, not the executable.
Config* g_config = nullptr;

// Raised when the .azr cannot be loaded or the model will not initialize.
class AZURE2Error : public std::runtime_error {
 public:
  explicit AZURE2Error(const std::string& what) : std::runtime_error(what) {}
};

namespace {

// CoulFunc reads g_config at construction to decide whether the hybrid
// potential model is on, and construction happens deep inside a calculation.
// With several live sessions a single global is not enough: this publishes the
// calling session's Config for the duration of the call and puts back whatever
// was there before, so sessions cannot silently reconfigure each other.
struct ConfigScope {
  explicit ConfigScope(Config* c) : previous_(g_config) { g_config = c; }
  ~ConfigScope() { g_config = previous_; }
  ConfigScope(const ConfigScope&) = delete;
  ConfigScope& operator=(const ConfigScope&) = delete;
 private:
  Config* previous_;
};

// Convert a std::vector<double> to a new float64 numpy array (a copy, so the
// array owns its memory and outlives the C++ object).
py::array_t<double> to_array(const vector_r& v) {
  return py::array_t<double>(static_cast<py::ssize_t>(v.size()), v.data());
}

// The mirror image: a flat copy of any float64-coercible array.  The engine
// takes vector_r&, so it needs its own buffer either way.
vector_r to_vector(const py::array_t<double, py::array::forcecast>& a) {
  vector_r v(a.size());
  if (a.size() > 0) std::memcpy(v.data(), a.data(), v.size() * sizeof(double));
  return v;
}

// The runtime options the CLI used to carry as command-line flags.  Defaults
// mirror Config::Reset() plus the CLI's defaults (--no-transform off, etc.).
struct RuntimeOptions {
  bool data_mode = true;          // CALCULATE_WITH_DATA
  bool use_brune = true;          // USE_BRUNE_FORMALISM
  bool ignore_externals = true;   // IGNORE_ZERO_WIDTHS
  bool transform = true;          // TRANSFORM_PARAMETERS
  bool use_long_wavelength = true;// USE_LONGWAVELENGTH_APPROX
  bool use_gsl_coul = false;      // USE_GSL_COULOMB_FUNC
  bool use_rmc = false;           // USE_RMC_FORMALISM
};

void apply_options(Config& config, const RuntimeOptions& opt) {
  if (opt.data_mode) config.paramMask |= Config::CALCULATE_WITH_DATA;
  else config.paramMask &= ~Config::CALCULATE_WITH_DATA;
  if (opt.use_brune) config.paramMask |= Config::USE_BRUNE_FORMALISM;
  else config.paramMask &= ~Config::USE_BRUNE_FORMALISM;
  if (opt.ignore_externals) config.paramMask |= Config::IGNORE_ZERO_WIDTHS;
  else config.paramMask &= ~Config::IGNORE_ZERO_WIDTHS;
  if (opt.transform) config.paramMask |= Config::TRANSFORM_PARAMETERS;
  else config.paramMask &= ~Config::TRANSFORM_PARAMETERS;
  if (opt.use_long_wavelength) config.paramMask |= Config::USE_LONGWAVELENGTH_APPROX;
  else config.paramMask &= ~Config::USE_LONGWAVELENGTH_APPROX;
  if (opt.use_gsl_coul) config.paramMask |= Config::USE_GSL_COULOMB_FUNC;
  else config.paramMask &= ~Config::USE_GSL_COULOMB_FUNC;
  if (opt.use_rmc) config.paramMask |= Config::USE_RMC_FORMALISM;
  else config.paramMask &= ~Config::USE_RMC_FORMALISM;
}

}  // namespace

class Session {
 public:
  Session(const std::string& configfile, const RuntimeOptions& opt)
      : config_(nullptr), api_(nullptr) {
    config_ = new Config(std::cout);
    config_->configfile = configfile;
    config_->paramMask |= Config::USE_API;   // silence the "Calculating ..." chatter
    apply_options(*config_, opt);

    ConfigScope guard(config_);

    // The API stub used to guard the CLI's reads against --use-api; the binding
    // is the API, so mirror the checks main() performed for a headless run.
    if (config_->ReadConfigFile() == -1) {
      std::string msg = "Could not open " + configfile + " or no <config> block found.";
      fail(msg);
    }

    InitializeCoulFuncCache();
    InitializeECAmplitudeCache();

    api_ = new AZUREAPI(*config_);
    if (api_->Initialize() != 0) {
      std::string msg = "AZURE2 could not initialize " + configfile +
                        "; check the output above.";
      fail(msg);
    }
  }

  ~Session() {
    if (api_ != nullptr) delete api_;
    // Only disown the global if it is still ours: a session destroyed while
    // another is live must not pull the config out from under it.
    if (g_config == config_) g_config = nullptr;
    if (config_ != nullptr) delete config_;
  }

  Session(const Session&) = delete;
  Session& operator=(const Session&) = delete;

  [[noreturn]] void fail(const std::string& msg) {
    if (api_ != nullptr) { delete api_; api_ = nullptr; }
    if (g_config == config_) g_config = nullptr;
    if (config_ != nullptr) { delete config_; config_ = nullptr; }
    throw AZURE2Error(msg);
  }

  // Every method below that runs the engine takes a ConfigScope; the plain
  // getters do not, because they only hand back vectors already computed.

  // -- lifecycle -----------------------------------------------------------

  bool rebuild() {
    ConfigScope guard(config_);
    return api_->Rebuild();
  }

  bool initialize() {
    ConfigScope guard(config_);
    return api_->Initialize() == 0;
  }
  bool calculate_external_capture() {
    ConfigScope guard(config_);
    return api_->CalculateExternalCapture();
  }

  // -- modes / rebuilds ------------------------------------------------------

  void set_data() { api_->SetData(); }
  void set_extrap() { api_->SetExtrap(); }
  bool set_radius(int idx, double radius) {
    ConfigScope guard(config_);
    return api_->SetRadius(idx, radius);
  }

  // -- hybrid nuclear potential, per particle pair ---------------------------
  //
  // A potential belongs to a pair: it bends that channel's radial wave
  // functions and no other.  NuclearPotentialManager holds one setting per
  // pair plus a default for the pairs that name none, and CoulFunc reads its
  // own pair's setting at construction -- so a change only reaches the model
  // once it is rebuilt, which is what reinitializing does.

  // pair = 0 addresses the default that every unnamed pair falls back to.
  void set_potential(int pair, const std::string& type, bool enabled,
                     double V0, double R, double a, double r0) {
    NuclearPotentialSetting s;
    s.enabled = enabled;
    s.type = type;
    s.V0 = V0;
    s.R = R;
    s.a = a;
    s.r0 = r0;
    NuclearPotentialManager& mgr = NuclearPotentialManager::instance();
    try {
      if (pair) mgr.setSetting(pair, s);
      else mgr.setDefaultSetting(s);
    } catch (const std::exception& e) {
      throw AZURE2Error(e.what());
    }
    // The master switch follows the pairs: with no pair enabled the hybrid
    // model is off, and CoulFunc short-circuits on it.
    config_->useHybridMethod = mgr.isAnyEnabled();
    if (config_->useHybridMethod) config_->paramMask |= Config::USE_HYBRID_COULOMB;
    else config_->paramMask &= ~Config::USE_HYBRID_COULOMB;
  }

  void clear_potential(int pair) {
    NuclearPotentialManager& mgr = NuclearPotentialManager::instance();
    if (pair) mgr.clearPairSetting(pair);
    else mgr.resetToDefault();
    config_->useHybridMethod = mgr.isAnyEnabled();
    if (config_->useHybridMethod) config_->paramMask |= Config::USE_HYBRID_COULOMB;
    else config_->paramMask &= ~Config::USE_HYBRID_COULOMB;
  }

  // (enabled, type, V0, R, a, r0, has_own_setting) for the given pair,
  // resolved against the default.
  py::tuple get_potential(int pair) const {
    const NuclearPotentialManager& mgr = NuclearPotentialManager::instance();
    NuclearPotentialSetting s = pair ? mgr.getSetting(pair) : mgr.getDefaultSetting();
    bool own = pair ? mgr.hasPairSetting(pair) : true;
    return py::make_tuple(s.enabled, s.type, s.V0, s.R, s.a, s.r0, own);
  }

  std::vector<int> configured_potential_pairs() const {
    return NuclearPotentialManager::instance().configuredPairs();
  }

  // -- data bookkeeping ------------------------------------------------------

  int update_data() {
    ConfigScope guard(config_);
    return api_->UpdateData();
  }
  void update_norms() { api_->UpdateNorms(); }
  bool update_parameters() { return api_->UpdateParameters(); }

  // -- parameter metadata ----------------------------------------------------

  py::array_t<double> parameter_info() const { return to_array(api_->GetParameterInfo()); }
  py::array_t<double> pairs_info() const { return to_array(api_->GetPairsInfo()); }
  py::array_t<double> normalization_indices() const { return to_array(api_->GetNormalizationIndices()); }
  py::array_t<double> energy_shift_indices() const { return to_array(api_->GetEnergyShiftIndices()); }

  // -- parameter vectors -----------------------------------------------------

  py::array_t<double> params_values() const { return to_array(api_->params_values()); }
  py::array_t<double> params_values_rwa() const { return to_array(api_->params_values_rwa()); }
  py::array_t<double> params_all() const { return to_array(api_->params_all()); }
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
  py::array_t<double> norms() { api_->UpdateNorms(); return to_array(api_->norms()); }
  py::array_t<double> norms_errors() { api_->UpdateNorms(); return to_array(api_->norms_errors()); }

  // -- calculated arrays -----------------------------------------------------

  py::array_t<double> calculated_segments(int i) const { return to_array(api_->calculated_segments(i)); }
  py::array_t<double> calculated_segments_e1(int i) const { return to_array(api_->calculated_segments_e1(i)); }
  py::array_t<double> calculated_segments_e2(int i) const { return to_array(api_->calculated_segments_e2(i)); }
  py::array_t<double> calculated_energies(int i) const { return to_array(api_->calculated_energies(i)); }
  py::array_t<double> calculated_angles(int i) const { return to_array(api_->calculated_angles(i)); }
  py::array_t<double> calculated_angular_dists(int i) const { return to_array(api_->calculated_angular_dists(i)); }
  py::array_t<double> calculated_excitation_energies(int i) const { return to_array(api_->calculated_excitation_energies(i)); }
  py::array_t<double> calculated_conv(int i) const { return to_array(api_->calculated_conv(i)); }

  // -- segment updates (return the number of segments) ------------------------

  int update_segments(py::array_t<double, py::array::forcecast> p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->UpdateSegments(v);
  }
  int update_segments_rwa(py::array_t<double, py::array::forcecast> p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->UpdateSegmentsRWA(v);
  }
  int update_segments_all_rwa(py::array_t<double, py::array::forcecast> p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->UpdateSegmentsAllRWA(v);
  }

  // -- transforms ------------------------------------------------------------

  py::array_t<double> transform_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    { py::gil_scoped_release release; ConfigScope guard(config_);
      out = api_->TransformRWAParameters(v); }
    return to_array(out);
  }
  py::array_t<double> transform_all_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    { py::gil_scoped_release release; ConfigScope guard(config_);
      out = api_->TransformAllRWAParameters(v); }
    return to_array(out);
  }

  // -- chi-squared and derivatives -------------------------------------------

  double calculate_chi2_rwa(py::array_t<double, py::array::forcecast> p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->CalculateChi2RWA(v);
  }
  double calculate_chi2_physical(py::array_t<double, py::array::forcecast> p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->CalculateChi2Physical(v);
  }
  py::array_t<double> calculate_chi2_grad_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    { py::gil_scoped_release release; ConfigScope guard(config_);
      out = api_->CalculateChi2GradRWA(v); }
    return to_array(out);
  }
  py::array_t<double> calculate_residual_jacobian_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    { py::gil_scoped_release release; ConfigScope guard(config_);
      out = api_->CalculateResidualJacobianRWA(v); }
    return to_array(out);
  }
  py::array_t<double> calculate_model_gradients_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    { py::gil_scoped_release release; ConfigScope guard(config_);
      out = api_->CalculateModelGradientsRWA(v); }
    return to_array(out);
  }

  // -- external region and caches --------------------------------------------

  py::array_t<double> coulomb_functions(py::array_t<double, py::array::forcecast> request) {
    vector_r v = to_vector(request), out;
    { py::gil_scoped_release release; ConfigScope guard(config_);
      out = api_->GetCoulombFunctions(v); }
    return to_array(out);
  }
  py::array_t<double> ec_integrals(py::array_t<double, py::array::forcecast> request) {
    vector_r v = to_vector(request), out;
    { py::gil_scoped_release release; ConfigScope guard(config_);
      out = api_->GetECIntegrals(v); }
    return to_array(out);
  }
  py::array_t<double> cache_stats() const {
    vector_r out;
    { py::gil_scoped_release release; out = api_->GetCacheStats(); }
    return to_array(out);
  }

 private:
  Config* config_;
  AZUREAPI* api_;
};

PYBIND11_MODULE(_azure2, m) {
  m.doc() = "In-process bindings for the AZURE2 R-matrix engine.";

  // AZURE2.cpp installs this for the executable; the module has to do it for
  // itself, and it matters more here.  GSL's default handler calls abort(), so
  // an integral that will not converge -- which the hybrid nuclear potential
  // can produce -- would take the whole interpreter down instead of raising.
  gsl_set_error_handler(&GSLException::GSLErrorHandler);

  py::register_exception<AZURE2Error>(m, "AZURE2Error", PyExc_RuntimeError);
  py::register_exception<GSLException>(m, "GSLError", PyExc_ArithmeticError);

  py::class_<RuntimeOptions>(m, "RuntimeOptions")
      .def(py::init<>())
      .def_readwrite("data_mode", &RuntimeOptions::data_mode)
      .def_readwrite("use_brune", &RuntimeOptions::use_brune)
      .def_readwrite("ignore_externals", &RuntimeOptions::ignore_externals)
      .def_readwrite("transform", &RuntimeOptions::transform)
      .def_readwrite("use_long_wavelength", &RuntimeOptions::use_long_wavelength)
      .def_readwrite("use_gsl_coul", &RuntimeOptions::use_gsl_coul)
      .def_readwrite("use_rmc", &RuntimeOptions::use_rmc);

  py::class_<Session>(m, "Session")
      .def(py::init<const std::string&, const RuntimeOptions&>(),
           py::arg("configfile"), py::arg("options") = RuntimeOptions(),
           "Load ``configfile`` and build the model.  Raises AZURE2Error if "
           "the file is missing or the model will not initialize.")
      .def("rebuild", &Session::rebuild,
           py::call_guard<py::gil_scoped_release>(),
           "Rebuild the model from the .azr, recomputing the external-capture "
           "integrals instead of reading back intEC.dat.  Needed after "
           "anything that changes the Coulomb functions.")
      .def("initialize", &Session::initialize,
           py::call_guard<py::gil_scoped_release>(),
           "Re-run the (expensive) model initialization; needed after a mode "
           "switch.")
      .def("calculate_external_capture", &Session::calculate_external_capture,
           py::call_guard<py::gil_scoped_release>(),
           "Recompute every external-capture integral from scratch.")
      .def("set_data", &Session::set_data, "Switch to the <segmentsData> mode.")
      .def("set_extrap", &Session::set_extrap, "Switch to the <segmentsTest> mode.")
      .def("set_radius", &Session::set_radius, py::arg("pair"), py::arg("radius"),
           py::call_guard<py::gil_scoped_release>(),
           "Rebuild the model with one pair's channel radius changed (fm).")
      .def("set_potential", &Session::set_potential, py::arg("pair"),
           py::arg("type"), py::arg("enabled"), py::arg("V0"), py::arg("R"),
           py::arg("a"), py::arg("r0"),
           "Set the hybrid nuclear potential for one particle pair; pair=0 "
           "sets the default the unnamed pairs fall back to.  Call "
           "initialize() afterwards for it to reach the model.")
      .def("clear_potential", &Session::clear_potential, py::arg("pair"),
           "Drop a pair's own potential so it follows the default again; "
           "pair=0 resets the default and every pair.")
      .def("get_potential", &Session::get_potential, py::arg("pair"),
           "(enabled, type, V0, R, a, r0, has_own_setting) for a pair.")
      .def("configured_potential_pairs", &Session::configured_potential_pairs,
           "The pair keys carrying a potential of their own.")
      .def("update_data", &Session::update_data, py::call_guard<py::gil_scoped_release>(),
           "Refill the data arrays; returns the number of segments.")
      .def("update_norms", &Session::update_norms, "Refill the normalization arrays.")
      .def("update_parameters", &Session::update_parameters, "Refill the parameter bookkeeping.")
      .def("parameter_info", &Session::parameter_info)
      .def("pairs_info", &Session::pairs_info)
      .def("normalization_indices", &Session::normalization_indices)
      .def("energy_shift_indices", &Session::energy_shift_indices)
      .def("params_values", &Session::params_values)
      .def("params_values_rwa", &Session::params_values_rwa)
      .def("params_all", &Session::params_all)
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
      .def("norms_errors", &Session::norms_errors)
      .def("calculated_segments", &Session::calculated_segments, py::arg("segment"))
      .def("calculated_segments_e1", &Session::calculated_segments_e1, py::arg("segment"))
      .def("calculated_segments_e2", &Session::calculated_segments_e2, py::arg("segment"))
      .def("calculated_energies", &Session::calculated_energies, py::arg("segment"))
      .def("calculated_angles", &Session::calculated_angles, py::arg("segment"))
      .def("calculated_angular_dists", &Session::calculated_angular_dists, py::arg("segment"))
      .def("calculated_excitation_energies", &Session::calculated_excitation_energies, py::arg("segment"))
      .def("calculated_conv", &Session::calculated_conv, py::arg("segment"))
      .def("update_segments", &Session::update_segments,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("update_segments_rwa", &Session::update_segments_rwa,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("update_segments_all_rwa", &Session::update_segments_all_rwa,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("transform_rwa", &Session::transform_rwa, py::arg("params"))
      .def("transform_all_rwa", &Session::transform_all_rwa, py::arg("params"))
      .def("calculate_chi2_rwa", &Session::calculate_chi2_rwa,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("calculate_chi2_physical", &Session::calculate_chi2_physical,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("calculate_chi2_grad_rwa", &Session::calculate_chi2_grad_rwa,
           py::arg("params"))
      .def("calculate_residual_jacobian_rwa", &Session::calculate_residual_jacobian_rwa,
           py::arg("params"))
      .def("calculate_model_gradients_rwa", &Session::calculate_model_gradients_rwa,
           py::arg("params"))
      .def("coulomb_functions", &Session::coulomb_functions, py::arg("request"))
      .def("ec_integrals", &Session::ec_integrals, py::arg("request"))
      .def("cache_stats", &Session::cache_stats);
}
