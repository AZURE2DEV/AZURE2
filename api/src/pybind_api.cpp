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

#include <cstdio>
#include <fstream>
#include <sstream>
#include "azure2/io/ProjectFormat.h"
#include "azure2/io/Document.h"
#include "AZUREAPI.h"
#include "CNuc.h"
#include "ParameterKind.h"
#include "Config.h"
#include "GSLException.h"
#include "ECAmplitudeCache.h"
#include <gsl/gsl_errno.h>

namespace py = pybind11;

// Defined by AZURE2.cpp for the CLI/GUI executable; here the module owns it,
// because the binding links only the core library, not the executable.
Config *g_config = nullptr;

// Raised when the .azr cannot be loaded or the model will not initialize.
class AZURE2Error : public std::runtime_error {
 public:
  explicit AZURE2Error(const std::string &what) :
    std::runtime_error(what) {}
};

// pybind_potentials.cpp: a session's hybrid potentials
azure2::experiment::HybridPotentialSetting SessionPotentialSetting(const std::string &type, bool enabled, double V0, double R, double a, double r0,
                                                                   double coulombRadius_fm, double matchingRadius_fm, double step_fm, const std::vector<double> &spinOrbit, CNuc *compound, int pair);
pybind11::tuple SessionPotentialView(const Config &config, int pair);
pybind11::tuple SessionPotentialRecords(const Config &config);

namespace {

/// The bytes of a project file ("" when it cannot be read), for the refusal of one that is not format 3.
std::string ReadProjectText(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

// CoulFunc reads g_config at construction to decide whether the hybrid
// potential model is on, and construction happens deep inside a calculation.
// With several live sessions a single global is not enough: this publishes the
// calling session's Config for the duration of the call and puts back whatever
// was there before, so sessions cannot silently reconfigure each other.  It is
// not a lock: calls run with the GIL released, so two Python threads driving
// sessions at once would race on g_config -- hence "one thread" above.
struct ConfigScope {
  // Besides the Config, every call swaps in the session's own external-capture
  // amplitude cache.  `g_ecAmplitudeCache` is process-global and its keys
  // (k-group, EC m-group, entrance, exit, segment key) say nothing about which
  // session they belong to, so with one shared object a second session on the
  // same model -- an extrapolation grid beside a data session, say -- first
  // deleted the first session's entries (InitializeECAmplitudeCache) and then
  // filled the cache under colliding keys, and the first session interpolated
  // its capture amplitudes from the wrong grid: residuals off by up to 28 %,
  // silently, for as long as the process lived (8Be+alpha, 2026-10-07).
  // `slot` is the session's cache pointer; whatever the call leaves in the
  // global (Initialize re-creates it) is captured back into the slot on exit.
  // The session's hybrid nuclear potentials need no swap: they are values in
  // the session's Config, and each pair of its CNuc holds its own.
  ConfigScope(Config *c, ECAmplitudeCache **slot) :
    previous_(g_config), previousCache_(g_ecAmplitudeCache), slot_(slot) {
    g_config = c;
    g_ecAmplitudeCache = *slot_;
  }
  ~ConfigScope() {
    *slot_ = g_ecAmplitudeCache;
    g_ecAmplitudeCache = previousCache_;
    g_config = previous_;
  }
  ConfigScope(const ConfigScope &) = delete;
  ConfigScope &operator=(const ConfigScope &) = delete;

 private:
  Config *previous_;
  ECAmplitudeCache *previousCache_;
  ECAmplitudeCache **slot_;
};

// Convert a std::vector<double> to a new float64 numpy array (a copy, so the
// array owns its memory and outlives the C++ object).
py::array_t<double> to_array(const vector_r &v) {
  return py::array_t<double>(static_cast<py::ssize_t>(v.size()), v.data());
}

// The mirror image: a flat copy of any float64-coercible array.  The engine
// takes vector_r&, so it needs its own buffer either way.  Bindings that
// release the GIL for the whole call (py::call_guard) take their array by
// const reference: a by-value py::array is reference-counted inside the
// GIL-free call (copy in, destroy out), and a temporary converted from a
// Python list then crashed (segfault on calculate_chi2_rwa([...]), October
// 2026).  The others convert their argument with the GIL held and release it
// only around the engine call (an inner block: gil_scoped_release, then
// ConfigScope, so that the scope ends before the GIL is taken back); they may
// take the array by value.
vector_r to_vector(const py::array_t<double, py::array::forcecast> &a) {
  vector_r v(a.size());
  if (a.size() > 0) std::memcpy(v.data(), a.data(), v.size() * sizeof(double));
  return v;
}

// The runtime options the CLI used to carry as command-line flags.  Defaults
// mirror Config::Reset() plus the CLI's defaults (--no-transform off, etc.).
struct RuntimeOptions {
  bool data_mode = true;            // CALCULATE_WITH_DATA
  bool use_brune = true;            // USE_BRUNE_FORMALISM
  bool ignore_externals = true;     // IGNORE_ZERO_WIDTHS
  bool transform = true;            // TRANSFORM_PARAMETERS
  bool use_long_wavelength = true;  // USE_LONGWAVELENGTH_APPROX
  bool use_gsl_coul = false;        // USE_GSL_COULOMB_FUNC
  bool use_rmc = false;             // USE_RMC_FORMALISM
  bool use_park = false;            // USE_PARK_FORMALISM (implies Brune)
  bool use_coulomb_table = true;    // Config::useCoulombTable
  bool replace_grids = true;        // Config::replaceSubPointGrids
  bool use_ec_table = true;         // Config::useECTable
  bool ec_convection_current = true;  // Config::ecConvectionCurrent (exact-form external capture only)
  std::string coulomb_table;        // Config::coulombTablePath ("" = the default file)
};

void apply_options(Config &config, const RuntimeOptions &opt) {
  if (opt.data_mode)
    config.paramMask |= Config::CALCULATE_WITH_DATA;
  else
    config.paramMask &= ~Config::CALCULATE_WITH_DATA;
  if (opt.use_brune)
    config.paramMask |= Config::USE_BRUNE_FORMALISM;
  else
    config.paramMask &= ~Config::USE_BRUNE_FORMALISM;
  if (opt.ignore_externals)
    config.paramMask |= Config::IGNORE_ZERO_WIDTHS;
  else
    config.paramMask &= ~Config::IGNORE_ZERO_WIDTHS;
  if (opt.transform)
    config.paramMask |= Config::TRANSFORM_PARAMETERS;
  else
    config.paramMask &= ~Config::TRANSFORM_PARAMETERS;
  if (opt.use_long_wavelength)
    config.paramMask |= Config::USE_LONGWAVELENGTH_APPROX;
  else
    config.paramMask &= ~Config::USE_LONGWAVELENGTH_APPROX;
  if (opt.use_gsl_coul)
    config.paramMask |= Config::USE_GSL_COULOMB_FUNC;
  else
    config.paramMask &= ~Config::USE_GSL_COULOMB_FUNC;
  if (opt.use_rmc)
    config.paramMask |= Config::USE_RMC_FORMALISM;
  else
    config.paramMask &= ~Config::USE_RMC_FORMALISM;
  config.useCoulombTable = opt.use_coulomb_table;
  config.coulombTablePath = opt.coulomb_table;
  config.replaceSubPointGrids = opt.replace_grids;
  config.useECTable = opt.use_ec_table;
  config.ecConvectionCurrent = opt.ec_convection_current;
  // Park shares Brune's per-level shift functions, so it switches Brune on.
  if (opt.use_park)
    config.paramMask |= (Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  else
    config.paramMask &= ~Config::USE_PARK_FORMALISM;
}

}  // namespace

void AddModelData(py::dict &terms, const AZUREAPI &api, const vector_r &free);  // pybind_model_data.cpp
class Session {
 public:
  Session(const std::string &configfile, const RuntimeOptions &opt) :
    config_(nullptr),
    api_(nullptr) {
    config_ = new Config(std::cout);
    config_->configfile = configfile;
    config_->paramMask |= Config::USE_API;  // silence the "Calculating ..." chatter
    apply_options(*config_, opt);

    ConfigScope guard(config_, &ecCache_);

    // The binding is the API: mirror the checks main() performs for a headless run.
    int configStatus = config_->ReadConfigFile();
    if (configStatus == Config::kConfigUnreadable) {
      fail("Could not open " + configfile + " or no <config> block found.");
    } else if (configStatus == Config::kConfigNotFormat3) {
      fail(azure2::io::NotFormat3Message(ReadProjectText(configfile), configfile));
    } else if (configStatus == Config::kConfigRefused) {
      fail(configfile + " (format 3) is refused (see the ERROR lines above; AZURE2 --validate lists them).");
    } else if (configStatus < 0) {
      fail("Malformed <thm> section in " + configfile + " (see the ERROR above).");
    }
    // R-matrix route + ignore_externals (the default): A-matrix with a warning, as the CLI and the GUI.  Before
    // 2026-11-01 pyazr ran the R-matrix route and silently kept the internal-capture term ignore_externals drops
    // (13N capture 3e-5 off the A-matrix route).
    config_->ResolveMatrixRoute();

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
    if (ecCache_ != nullptr) {
      if (g_ecAmplitudeCache == ecCache_) g_ecAmplitudeCache = nullptr;
      delete ecCache_;
    }
  }

  Session(const Session &) = delete;
  Session &operator=(const Session &) = delete;

  [[noreturn]] void fail(const std::string &msg) {
    if (api_ != nullptr) {
      delete api_;
      api_ = nullptr;
    }
    if (g_config == config_) g_config = nullptr;
    if (config_ != nullptr) {
      delete config_;
      config_ = nullptr;
    }
    throw AZURE2Error(msg);
  }

  // Every method below that runs the engine takes a ConfigScope; the plain
  // getters do not, because they only hand back vectors already computed.

  // -- lifecycle -----------------------------------------------------------

  bool rebuild() {
    ConfigScope guard(config_, &ecCache_);
    return api_->Rebuild();
  }

  bool write_output_files() {
    ConfigScope guard(config_, &ecCache_);
    return api_->WriteOutputFiles();
  }

  bool initialize() {
    ConfigScope guard(config_, &ecCache_);
    return api_->Initialize() == 0;
  }
  bool calculate_external_capture() {
    ConfigScope guard(config_, &ecCache_);
    return api_->CalculateExternalCapture();
  }

  // -- modes / rebuilds ------------------------------------------------------

  void set_data() { api_->SetData(); }
  void set_extrap() { api_->SetExtrap(); }
  bool set_radius(int idx, double radius) {
    ConfigScope guard(config_, &ecCache_);
    return api_->SetRadius(idx, radius);
  }

  // -- hybrid nuclear potential, per particle pair ---------------------------
  //
  // A potential belongs to a pair: it bends that channel's radial wave
  // functions and no other.  The session's Config holds one setting per pair
  // plus a default for the pairs that name none (HybridPotentialSettings), and
  // CNuc::Fill gives each pair its potential -- so a change only reaches the
  // model once it is rebuilt, which pyazr does after every change.

  // pair = 0 addresses the default that every unnamed pair falls back to.
  void set_potential(int pair, const std::string &type, bool enabled, double V0, double R, double a, double r0,
                     double coulomb_radius, double matching_radius, double step, const std::vector<double> &spin_orbit) {
    try {
      azure2::experiment::HybridPotentialSetting s = SessionPotentialSetting(type, enabled, V0, R, a, r0, coulomb_radius, matching_radius, step, spin_orbit, api_->compound(), pair);
      (pair ? config_->potentials.pairs[pair] : config_->potentials.defaults) = s;
    } catch (const std::exception &e) {
      throw AZURE2Error(std::string(pair ? "Failed to set potential for pair: " : "Failed to set default potential: ") +
                        e.what());
    }
    PotentialsChanged();
  }

  void clear_potential(int pair) {
    if (pair)
      config_->potentials.pairs.erase(pair);
    else
      config_->potentials = azure2::experiment::HybridPotentialSettings();
    PotentialsChanged();
  }

  // (enabled, type, V0, R, a, r0, has_own_setting, coulomb_radius, spin_orbit, matching_radius, step) for the given pair, resolved against the default.
  py::tuple get_potential(int pair) { return SessionPotentialView(*config_, pair); }
  // (changed, record lines, notes): the session's potentials as format-3 records (pybind_potentials.cpp)
  py::tuple potential_records() { return SessionPotentialRecords(*config_); }

  std::vector<int> configured_potential_pairs() { return config_->potentials.ConfiguredPairs(); }

  // -- data bookkeeping ------------------------------------------------------

  int update_data() {
    ConfigScope guard(config_, &ecCache_);
    return api_->UpdateData();
  }
  void update_norms() { api_->UpdateNorms(); }
  bool update_parameters() { return api_->UpdateParameters(); }

  // -- parameter metadata ----------------------------------------------------

  py::array_t<double> parameter_info() const { return to_array(api_->GetParameterInfo()); }
  /// One dict per parameter, in params_names() order: its kind (ParamKindTable key) and what it belongs to, from the
  /// record written where it was created (ParameterKind.h).  Fields that do not apply are absent.
  py::list parameter_records() const {
    py::list out;
    CNuc *nuc = api_->compound();
    for (const ParameterRecord &r : api_->GetParameterRecords()) {
      py::dict d;
      d["kind"] = r.Row().key;
      d["code"] = (int)r.kind;
      d["ordinal"] = r.ordinal;
      switch (r.kind) {
        case ParamKind::LevelEnergy:
        case ParamKind::Gamma:
          d["jgroup"] = r.jGroup;
          d["level"] = r.level;
          d["level_number"] = r.levelNumber;
          if (r.kind == ParamKind::Gamma) d["channel"] = r.channel;
          break;
        case ParamKind::Norm:
        case ParamKind::EnergyShift:
        case ParamKind::EnergyShiftSqrt:
        case ParamKind::EnergyShiftFn:
        case ParamKind::Thickness:
          d["segment"] = r.segment;
          d["segment_key"] = r.segmentKey;
          if (r.kind == ParamKind::EnergyShiftFn) d["fn_term"] = r.fnTerm;
          break;
        case ParamKind::Screening:
        case ParamKind::Density: case ParamKind::Doppler: case ParamKind::Apparatus:
          d["effect"] = r.effect;
          if (r.kind == ParamKind::Apparatus) d["apparatus"] = r.apparatus, d["unit"] = r.apparatusUnit, d["owner"] = r.owner;
          break;
        case ParamKind::Abundance:
          d["sample"] = r.sample;
          d["sample_name"] = r.sampleName;
          d["nucleus"] = r.nucleus;
          d["nucleus_Z"] = r.nucleusZ;
          d["nucleus_A"] = r.nucleusA;
          break;
        case ParamKind::ThmCoherent: {
          d["thm_experiment"] = r.thmExperiment;
          d["jgroup"] = r.jGroup;
          d["channel"] = r.channel;
          d["exit_channel"] = r.exitChannel;
          d["thm_term"] = r.thmTerm;
          // The combination as cbackground= writes it: J^pi, exit pair key, (s,l) -> (s',l').
          JGroup *jg = nuc->GetJGroup(r.jGroup);
          AChannel *in = jg->GetChannel(r.channel);
          AChannel *ex = jg->GetChannel(r.exitChannel);
          d["thm_jpi"] = ThmSpinText(jg->GetJ()) + (jg->GetPi() > 0 ? "+" : "-");
          d["thm_exit_key"] = nuc->GetPair(ex->GetPairNum())->GetPairKey();
          d["thm_channels"] = ThmSpinText(in->GetS()) + "," + std::to_string(in->GetL()) + "," +
                              ThmSpinText(ex->GetS()) + "," + std::to_string(ex->GetL());
          break;
        }
      }
      out.append(d);
    }
    return out;
  }
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
  py::array_t<double> norms() {
    api_->UpdateNorms();
    return to_array(api_->norms());
  }
  py::array_t<double> norms_errors() {
    api_->UpdateNorms();
    return to_array(api_->norms_errors());
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
  void expect_free(const vector_r &v, const char *what) const {
    std::vector<bool> fixed = api_->params_fixed();
    size_t nfree = 0;
    for (bool f : fixed) if (!f) nfree++;
    if (v.size() != nfree)
      throw AZURE2Error(std::string(what) + ": expected the free parameter vector (" +
                        std::to_string(nfree) + " values, params_rwa / params order), got " +
                        std::to_string(v.size()) + ".");
  }
  void expect_all(const vector_r &v, const char *what) const {
    size_t nall = api_->params_fixed().size();
    if (v.size() != nall)
      throw AZURE2Error(std::string(what) + ": expected every parameter (" + std::to_string(nall) +
                        " values, fixed ones included), got " + std::to_string(v.size()) + ".");
  }

  // -- segment updates (return the number of segments) ------------------------

  int update_segments(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_free(v, "update_segments");
    return api_->UpdateSegments(v);
  }
  int update_segments_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_free(v, "update_segments_rwa");
    return api_->UpdateSegmentsRWA(v);
  }
  void set_params_all_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_all(v, "set_params_all_rwa");
    api_->SetParamsAllRWA(v);
  }
  int update_segments_all_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_all(v, "update_segments_all_rwa");
    return api_->UpdateSegmentsAllRWA(v);
  }

  // -- transforms ------------------------------------------------------------

  py::array_t<double> transform_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->TransformRWAParameters(v);
    }
    return to_array(out);
  }
  py::array_t<double> transform_all_rwa(const py::array_t<double, py::array::forcecast> &p, bool include_fixed) {
    vector_r v = to_vector(p), out;
    expect_all(v, "transform_all_rwa");
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->TransformAllRWAParameters(v, include_fixed);
    }
    return to_array(out);
  }
  py::array_t<double> transform_all_rwa_full(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p), out;
    expect_all(v, "transform_all_rwa_full");
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->TransformAllRWAParametersFull(v);
    }
    return to_array(out);
  }

  // -- chi-squared and derivatives -------------------------------------------

  double calculate_chi2_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_free(v, "calculate_chi2_rwa");
    return api_->CalculateChi2RWA(v);
  }
  double calculate_chi2_physical(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_free(v, "calculate_chi2_physical");
    return api_->CalculateChi2Physical(v);
  }
  py::array_t<double> scattering_lengths_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_free(v, "scattering_lengths_rwa");
    return to_array(api_->CalculateScatteringLengthsRWA(v));
  }
  py::array_t<double> park_norms(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_free(v, "park_norms");
    return to_array(api_->ParkNorms(v));
  }
  py::array_t<double> prior_covariance_residuals(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_free(v, "prior_covariance_residuals");
    return to_array(api_->PriorCovarianceResidualsRWA(v));
  }
  std::vector<std::string> prior_covariance_names() {
    ConfigScope guard(config_, &ecCache_);
    return api_->PriorCovarianceNames();
  }
  int refreeze_covariance() {
    ConfigScope guard(config_, &ecCache_);
    return api_->RefreezeCovariance();
  }
  py::array_t<double> calculate_chi2_grad_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p), out;
    expect_free(v, "calculate_chi2_grad_rwa");
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->CalculateChi2GradRWA(v);
    }
    return to_array(out);
  }
  // The fit's objective (docs/upgrade/objective.md): its value, its gradient, and its terms (names, values, and the
  // normalization and shift priors per segment, index 0 unused).
  double calculate_objective_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
    vector_r v = to_vector(p);
    expect_free(v, "calculate_objective_rwa");
    return api_->CalculateObjectiveRWA(v);
  }
  py::array_t<double> calculate_objective_grad_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p), out;
    expect_free(v, "calculate_objective_grad_rwa");
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->CalculateObjectiveGradRWA(v);
    }
    return to_array(out);
  }
  py::dict objective_terms_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_, &ecCache_);
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
  py::array_t<double> calculate_residual_jacobian_rwa(const py::array_t<double, py::array::forcecast> &p,
                                                      const std::vector<int> &skip) {
    vector_r v = to_vector(p), out;
    expect_free(v, "calculate_residual_jacobian_rwa");
    for (int c : skip)
      if (c < 0 || c >= (int)v.size())
        throw AZURE2Error("calculate_residual_jacobian_rwa: skip column " + std::to_string(c) + " is outside the " +
                          std::to_string(v.size()) + " free parameters.");
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->CalculateResidualJacobianRWA(v, skip);
    }
    return to_array(out);
  }
  py::array_t<double> calculate_residuals_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->CalculateResidualsRWA(v);
    }
    return to_array(out);
  }
  py::array_t<double> current_norms() {
    ConfigScope guard(config_, &ecCache_);
    return to_array(api_->GetCurrentNorms());
  }
  // The THM experiments as the last chi-squared evaluation profiled them.
  py::list thm_experiments() {
    ConfigScope guard(config_, &ecCache_);
    py::list out;
    for (const ThmExperimentReport &r : api_->GetThmExperiments()) {
      py::dict d;
      d["name"] = r.name;
      d["segments"] = r.segments;
      d["background"] = r.background;
      d["points"] = r.points;
      d["chi2"] = r.chi2;
      d["status"] = r.status;
      d["value"] = std::vector<double>(r.value, r.value + 4);
      d["covariance"] = std::vector<double>(r.covariance, r.covariance + 16);
      if (!r.coherentNames.empty()) {
        d["cbkg_names"] = r.coherentNames;
        d["cbkg_values"] = r.coherentValues;
      }
      out.append(d);
    }
    return out;
  }
  // The Coulomb line shape of a THM experiment at the current parameters.
  py::dict thm_lineshape(const std::string &name, py::array_t<double, py::array::forcecast> e) {
    ConfigScope guard(config_, &ecCache_);
    ThmLineshapeReport r;
    std::string why;
    if (!api_->GetThmLineshape(name, to_vector(e), r, why)) throw AZURE2Error(why);
    py::dict d;
    d["experiment"] = r.experiment;
    d["spectator"] = r.spectator;
    d["Zs"] = r.Zs;
    d["ZF"] = r.ZF;
    d["E_aA"] = r.eAA;
    d["B"] = r.bind;
    d["E"] = to_array(r.energy);
    d["E_sF"] = to_array(r.esf);
    d["eta_0"] = to_array(r.eta0);
    py::list exits;
    for (const ThmLineshapeReport::Exit &x : r.exits) {
      py::dict xd;
      xd["pair"] = x.pairKey;
      xd["Zb"] = x.Zb;
      xd["ZB"] = x.ZB;
      xd["mb"] = x.mb;
      xd["mB"] = x.mB;
      xd["zeta"] = to_array(x.zeta);
      xd["eta_sb"] = to_array(x.etaSb);
      py::list levels;
      for (const ThmLineshapeReport::Level &l : x.levels) {
        py::dict ld;
        ld["jgroup"] = l.jgroup;
        ld["level"] = l.level;
        ld["J"] = l.J;
        ld["pi"] = l.pi;
        ld["E_level"] = l.energy;
        ld["Gamma"] = l.width;
        ld["NC2"] = to_array(l.nc2);
        levels.append(ld);
      }
      xd["levels"] = levels;
      exits.append(xd);
    }
    d["exits"] = exits;
    return d;
  }
  // The window-averaged THM entrance vertex of an experiment at the current parameters.
  py::dict thm_vertex(const std::string &name, py::array_t<double, py::array::forcecast> e) {
    ThmVertexReport r;
    std::string why;
    const vector_r energies = to_vector(e);
    bool ok;
    {
      // The DW vertex and R(E) can take seconds: without the GIL.
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      ok = api_->GetThmVertex(name, energies, r, why);
    }
    if (!ok) throw AZURE2Error(why);
    py::dict d;
    d["experiment"] = r.experiment;
    d["window"] = r.window;
    d["pair"] = r.pairKey;
    d["mu_sx"] = r.muSx;
    d["B"] = r.bind;
    d["radius"] = r.radius;
    // The nodes at each energy (they follow the accepted directions): lists
    // over E of arrays over the nodes; empty lists with the DW vertex.
    auto rows = [&](const std::vector<std::vector<double>> &v) {
      py::list l;
      for (const std::vector<double> &row : v) l.append(to_array(row));
      return l;
    };
    d["p_s"] = rows(r.p);
    d["weights"] = rows(r.weight);
    d["T_s"] = rows(r.es);
    if (!r.theta.empty()) d["theta_cm"] = rows(r.theta);
    d["E"] = to_array(r.energy);
    {
      // Where the window (or the DW grid) does not reach E, the engine uses
      // the nodes of the nearest data point / grid energy (as for a folding
      // sub-point): flagged False here.
      py::array_t<bool> reached(r.reached.size());
      auto m = reached.mutable_unchecked<1>();
      for (size_t i = 0; i < r.reached.size(); i++) m(i) = r.reached[i] != 0;
      d["reached"] = reached;
    }
    py::list rho;
    for (const std::vector<double> &row : r.rho) rho.append(to_array(row));
    d["rho"] = rho;
    d["model"] = r.model;
    if (r.model == "dw") {
      py::list q, w;
      for (const std::vector<double> &row : r.dwQ) q.append(to_array(row));
      for (const std::vector<double> &row : r.dwWeight) w.append(to_array(row));
      d["dw_q"] = q;
      d["dw_weights"] = w;
      d["dw_q_delta"] = to_array(r.dwQDelta);
      d["dw_p_delta"] = to_array(r.dwPDelta);
      if (!r.dwTheta.empty()) {
        py::list th, aq, aw;
        for (const std::vector<double> &row : r.dwTheta) th.append(to_array(row));
        for (const std::vector<double> &row : r.dwAngleQ) aq.append(to_array(row));
        for (const std::vector<double> &row : r.dwAngleWeight) aw.append(to_array(row));
        d["angle_theta_cm"] = th;
        d["angle_q"] = aq;
        d["angle_weights"] = aw;
      }
    }
    py::list channels;
    for (const ThmVertexReport::Channel &c : r.channels) {
      py::dict cd;
      cd["jgroup"] = c.jgroup;
      cd["channel"] = c.channel;
      cd["J"] = c.J;
      cd["pi"] = c.pi;
      cd["l"] = c.l;
      cd["s"] = c.s;
      py::list levels;
      for (const ThmVertexReport::Level &l : c.levels) {
        py::dict ld;
        ld["level"] = l.level;
        ld["boundary"] = l.boundary;
        ld["M2"] = to_array(l.m2);
        ld["M2_qf"] = to_array(l.m2qf);
        if (!l.m2pw.empty()) ld["M2_pw"] = to_array(l.m2pw);
        levels.append(ld);
      }
      cd["levels"] = levels;
      channels.append(cd);
    }
    d["channels"] = channels;
    return d;
  }
  // The distortion factor R(E) of a THM experiment (distortion=...).
  py::dict thm_distortion(const std::string &name, py::array_t<double, py::array::forcecast> e) {
    ThmDistortionReport r;
    std::string why;
    const vector_r energies = to_vector(e);
    bool ok;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      ok = api_->GetThmDistortion(name, energies, r, why);
    }
    if (!ok) throw AZURE2Error(why);
    py::dict d;
    d["experiment"] = r.experiment;
    d["kind"] = r.kind;
    d["description"] = r.description;
    d["E"] = to_array(r.energy);
    d["R_model"] = to_array(r.rModel);
    if (r.kind == "table") return d;
    d["E_ref"] = r.eRef;
    d["E_aA"] = r.eAA;
    d["B"] = r.bind;
    d["k_aA"] = r.kAA;
    d["eta_aA"] = r.etaAA;
    d["kappa"] = r.kappa;
    d["eta_b"] = r.etaB;
    d["beta"] = r.beta;
    d["E_sF"] = to_array(r.esf);
    d["k_sF"] = to_array(r.ksf);
    d["eta_sF"] = to_array(r.etasf);
    d["theta_cm"] = to_array(r.thetaCm);
    d["x"] = to_array(r.x);
    d["q"] = to_array(r.q);
    d["M2"] = to_array(r.m2);
    d["M2_PW"] = to_array(r.mpw2);
    d["R"] = to_array(r.r);
    d["lmax"] = r.lmax;
    return d;
  }
  py::array_t<double> calculate_model_gradients_rwa(const py::array_t<double, py::array::forcecast> &p) {
    vector_r v = to_vector(p), out;
    expect_free(v, "calculate_model_gradients_rwa");
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->CalculateModelGradientsRWA(v);
    }
    return to_array(out);
  }

  // -- external region and caches --------------------------------------------

  py::array_t<double> coulomb_functions(const py::array_t<double, py::array::forcecast> &request) {
    vector_r v = to_vector(request), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->GetCoulombFunctions(v);
    }
    return to_array(out);
  }
  py::array_t<double> ec_integrals(const py::array_t<double, py::array::forcecast> &request) {
    vector_r v = to_vector(request), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_, &ecCache_);
      out = api_->GetECIntegrals(v);
    }
    return to_array(out);
  }
  py::array_t<double> cache_stats() const {
    vector_r out;
    {
      py::gil_scoped_release release;
      out = api_->GetCacheStats();
    }
    return to_array(out);
  }

 private:
  /// After an edit of config_->potentials: the switch follows the pairs (on when any pair has the model), and the
  /// shift functions cached on the levels, computed with the old potential, are dropped
  /// (CNuc::CalcShiftFunctions skips a level that has not moved).
  void PotentialsChanged() {
    ConfigScope guard(config_, &ecCache_);
    config_->useHybridMethod = config_->potentials.AnyEnabled();
    if (config_->useHybridMethod)
      config_->paramMask |= Config::USE_HYBRID_COULOMB;
    else
      config_->paramMask &= ~Config::USE_HYBRID_COULOMB;
    if (api_ && api_->compound()) api_->compound()->InvalidateShiftFunctions();
  }

  Config *config_;
  ECAmplitudeCache *ecCache_ = nullptr;  // this session's external-capture amplitude cache
  AZUREAPI *api_;
};

void BindPhysics(py::module_ &m);  // pybind_physics.cpp: the free functions of physics/
void BindThm(py::module_ &m);      // pybind_thm.cpp: <thm> blocks as the engine reads them
void BindSchema(py::module_ &m);   // pybind_schema.cpp: the table of project format 3
void BindConvert(py::module_ &m);  // pybind_convert.cpp: format-3 validation and the converter old -> format 3
void BindProject(py::module_ &m);  // pybind_project.cpp: the project file as a document

PYBIND11_MODULE(_azure2, m) {
  m.doc() = "In-process bindings for the AZURE2 R-matrix engine.";

  // AZURE2.cpp installs this for the executable; the module has to do it for
  // itself, and it matters more here.  GSL's default handler calls abort(), so
  // an integral that will not converge -- which the hybrid nuclear potential
  // can produce -- would take the whole interpreter down instead of raising.
  gsl_set_error_handler(&GSLException::GSLErrorHandler);

  py::register_exception<AZURE2Error>(m, "AZURE2Error", PyExc_RuntimeError);
  BindPhysics(m);
  BindThm(m);  // pybind_thm.cpp: the engine's reading of a <thm> block

  m.def(
      "parameter_kinds",
      []() {
        py::list out;
        for (const ParamKindRow &row : ParamKindTable()) {
          py::dict d;
          d["code"] = (int)row.kind;
          d["key"] = row.key;
          d["unit"] = row.unit;
          d["owner"] = row.owner;
          d["own_prior"] = row.ownPrior;
          d["segment_field"] = row.segmentField;
          d["vary_from_project"] = row.varyFromProject;
          d["doc"] = row.doc;
          out.append(d);
        }
        return out;
      },
      "The kinds of fit parameter (ParameterKind.h), one dict per kind in the order of their codes: code, key (the "
      "Parameter.kind of pyazr), unit, owner fields, own_prior, segment_field, vary_from_project, doc.");
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
      .def_readwrite("use_brune", &RuntimeOptions::use_brune)
      .def_readwrite("ignore_externals", &RuntimeOptions::ignore_externals)
      .def_readwrite("transform", &RuntimeOptions::transform)
      .def_readwrite("use_long_wavelength", &RuntimeOptions::use_long_wavelength)
      .def_readwrite("use_gsl_coul", &RuntimeOptions::use_gsl_coul)
      .def_readwrite("use_coulomb_table", &RuntimeOptions::use_coulomb_table)
      .def_readwrite("coulomb_table", &RuntimeOptions::coulomb_table)
      .def_readwrite("replace_grids", &RuntimeOptions::replace_grids)
      .def_readwrite("use_ec_table", &RuntimeOptions::use_ec_table)
      .def_readwrite("ec_convection_current", &RuntimeOptions::ec_convection_current)
      .def_readwrite("use_rmc", &RuntimeOptions::use_rmc)
      .def_readwrite("use_park", &RuntimeOptions::use_park);

  py::class_<Session>(m, "Session")
      .def(py::init<const std::string &, const RuntimeOptions &>(),
           py::arg("configfile"), py::arg("options") = RuntimeOptions(),
           "Load ``configfile`` and build the model.  Raises AZURE2Error if "
           "the file is missing or the model will not initialize.")
      .def("write_output_files", &Session::write_output_files,
           py::call_guard<py::gil_scoped_release>(),
           "Write AZUREOut_*, chiSquared.out and the rest into the output "
           "directory, from the cross sections currently on the points.")
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
      .def("set_data", &Session::set_data, "Switch to the data run (the <data> records).")
      .def("set_extrap", &Session::set_extrap, "Switch to the extrapolation run (the <extrapolations> records).")
      .def("set_radius", &Session::set_radius, py::arg("pair"), py::arg("radius"),
           py::call_guard<py::gil_scoped_release>(),
           "Rebuild the model with one pair's channel radius changed (fm).")
      .def("set_potential", &Session::set_potential, py::arg("pair"),
           py::arg("type"), py::arg("enabled"), py::arg("V0"), py::arg("R"),
           py::arg("a"), py::arg("r0"), py::arg("coulomb_radius"), py::arg("matching_radius"), py::arg("step"), py::arg("spin_orbit"),
           "Set the hybrid nuclear potential for one particle pair (lengths fm); pair=0 "
           "sets the default the unnamed pairs fall back to.  spin_orbit: [] or [V_so MeV, R_so fm, a_so fm, "
           "on (0 light, 1 heavy)].  Call initialize() afterwards for it to reach the model.")
      .def("clear_potential", &Session::clear_potential, py::arg("pair"),
           "Drop a pair's own potential so it follows the default again; "
           "pair=0 resets the default and every pair.")
      .def("get_potential", &Session::get_potential, py::arg("pair"),
           "(enabled, type, V0, R, a, r0, has_own_setting, coulomb_radius, spin_orbit [] or [V_so, R_so, a_so, on], matching_radius, step) for a pair.")
      .def("potential_records", &Session::potential_records,
           "(changed, lines, notes): the potentials as format-3 records (a format-3 project only).")
      .def("configured_potential_pairs", &Session::configured_potential_pairs,
           "The pair keys carrying a potential of their own.")
      .def("update_data", &Session::update_data, py::call_guard<py::gil_scoped_release>(),
           "Refill the data arrays; returns the number of segments.")
      .def("update_norms", &Session::update_norms, "Refill the normalization arrays.")
      .def("update_parameters", &Session::update_parameters, "Refill the parameter bookkeeping.")
      .def("parameter_info", &Session::parameter_info)
      .def("parameter_records", &Session::parameter_records,
           "What each parameter is: one dict per parameter (kind, owner fields), in params_names() order.")
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
      .def("set_params_all_rwa", &Session::set_params_all_rwa,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"),
           "Set every parameter, fixed ones included (params_all_rwa layout), "
           "without calculating.  Fixed values are what the free-vector calls use.")
      .def("update_segments_all_rwa", &Session::update_segments_all_rwa,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("transform_rwa", &Session::transform_rwa, py::arg("params"))
      .def("transform_all_rwa", &Session::transform_all_rwa, py::arg("params"),
           py::arg("include_fixed") = false)
      .def("transform_all_rwa_full", &Session::transform_all_rwa_full, py::arg("params"),
           "As transform_all_rwa, returning every R-matrix parameter's physical value, fixed ones included.")
      .def("calculate_chi2_rwa", &Session::calculate_chi2_rwa,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("scattering_lengths_rwa", &Session::scattering_lengths_rwa, py::arg("params"),
           "s-wave scattering lengths of every particle pair at a free RWA vector: 6 numbers per channel "
           "(pair key, channel spin, parity, Re a0 fm, Im a0 fm, has-levels flag).")
      .def("park_norms", &Session::park_norms, py::arg("params"),
           "Park overlap J of every R-matrix level (1 in Brune mode).")
      .def("prior_covariance_residuals", &Session::prior_covariance_residuals, py::arg("params"),
           "Prior parameter covariance (<config> priorCovariance=): the whitened residuals L^-1 (p - mu) at the free "
           "RWA parameters (empty without one); their squared sum is the prior's chi2 term.")
      .def("prior_covariance_names", &Session::prior_covariance_names,
           "The prior covariance's parameter names, in its row order (fixed ones marginalized); empty without one.")
      .def("refreeze_covariance", &Session::refreeze_covariance,
           "Data covariances: forget the model frozen for the relvector terms; the next chi2 evaluation "
           "freezes it at its parameters.  Returns the number of covariance groups with relvector terms.")
      .def("calculate_chi2_physical", &Session::calculate_chi2_physical,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("calculate_objective_rwa", &Session::calculate_objective_rwa, py::arg("params"),
           "The fit's objective at the free RWA parameters: data chi^2 plus every prior and constraint.")
      .def("calculate_objective_grad_rwa", &Session::calculate_objective_grad_rwa, py::arg("params"),
           "[objective, gradient...] at the free RWA parameters.")
      .def("objective_terms_rwa", &Session::objective_terms_rwa, py::arg("params"),
           "The objective's terms: total, names, values, and per-segment norm/shift priors (index 0 unused).")
      .def("calculate_chi2_grad_rwa", &Session::calculate_chi2_grad_rwa,
           py::arg("params"))
      .def("calculate_residual_jacobian_rwa", &Session::calculate_residual_jacobian_rwa,
           py::arg("params"), py::arg("skip") = std::vector<int>(),
           "Residuals and Jacobian; the free columns listed in skip are not finite-differenced (left zero).")
      .def("calculate_residuals_rwa", &Session::calculate_residuals_rwa, py::arg("params"),
           "Standardized residuals (forward pass; THM norms profiled).")
      .def("current_norms", &Session::current_norms,
           "Norm each segment carries now (profiled THM norms at their optimum).")
      .def("thm_experiments", &Session::thm_experiments,
           "THM experiments (<thm> experiment[...]) as the last chi-squared evaluation "
           "profiled them: list of dicts (name, segments, background, points, chi2, status, "
           "value = [norm, b0, b1, b2], covariance = 4x4 row-major).")
      .def("thm_lineshape", &Session::thm_lineshape, py::arg("name"), py::arg("energies"),
           "Coulomb line shape of THM experiment `name` (lineshape=on) at c.m. energies, with the "
           "level poles of the last evaluation's parameters.")
      .def("thm_vertex", &Session::thm_vertex, py::arg("name"), py::arg("energies"),
           "THM entrance vertex of experiment `name` at c.m. energies: the spectator-momentum "
           "window at every energy (nodes p_s, weights, T_s, theta_cm, rho: lists over E) and, per "
           "entrance channel and level, <|M_l|^2> over the window and |M_l|^2 at p_s = 0.")
      .def("thm_distortion", &Session::thm_distortion, py::arg("name"), py::arg("energies"),
           "Distortion factor R(E) of THM experiment `name` (distortion=coulomb|optical|table) at "
           "c.m. energies: E_sF, eta_sF, the spectator angle, |M|^2, |M_PW|^2, R directly and as the "
           "model uses it (R_model).")
      .def("calculate_model_gradients_rwa", &Session::calculate_model_gradients_rwa,
           py::arg("params"))
      .def("coulomb_functions", &Session::coulomb_functions, py::arg("request"))
      .def("ec_integrals", &Session::ec_integrals, py::arg("request"))
      .def("cache_stats", &Session::cache_stats);
}
