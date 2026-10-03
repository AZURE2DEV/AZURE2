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
Config *g_config = nullptr;

// Raised when the .azr cannot be loaded or the model will not initialize.
class AZURE2Error : public std::runtime_error {
 public:
  explicit AZURE2Error(const std::string &what) :
    std::runtime_error(what) {}
};

namespace {

// CoulFunc reads g_config at construction to decide whether the hybrid
// potential model is on, and construction happens deep inside a calculation.
// With several live sessions a single global is not enough: this publishes the
// calling session's Config for the duration of the call and puts back whatever
// was there before, so sessions cannot silently reconfigure each other.
struct ConfigScope {
  explicit ConfigScope(Config *c) :
    previous_(g_config) { g_config = c; }
  ~ConfigScope() { g_config = previous_; }
  ConfigScope(const ConfigScope &) = delete;
  ConfigScope &operator=(const ConfigScope &) = delete;

 private:
  Config *previous_;
};

// Convert a std::vector<double> to a new float64 numpy array (a copy, so the
// array owns its memory and outlives the C++ object).
py::array_t<double> to_array(const vector_r &v) {
  return py::array_t<double>(static_cast<py::ssize_t>(v.size()), v.data());
}

// The mirror image: a flat copy of any float64-coercible array.  The engine
// takes vector_r&, so it needs its own buffer either way.  Bindings that
// release the GIL take their array by const reference: a by-value py::array
// is reference-counted inside the GIL-free call (copy in, destroy out), and a
// temporary converted from a Python list then crashed (segfault on
// calculate_chi2_rwa([...]), October 2026).
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
}

}  // namespace

class Session {
 public:
  Session(const std::string &configfile, const RuntimeOptions &opt) :
    config_(nullptr),
    api_(nullptr) {
    config_ = new Config(std::cout);
    config_->configfile = configfile;
    config_->paramMask |= Config::USE_API;  // silence the "Calculating ..." chatter
    apply_options(*config_, opt);

    ConfigScope guard(config_);

    // The API stub used to guard the CLI's reads against --use-api; the binding
    // is the API, so mirror the checks main() performed for a headless run.
    int configStatus = config_->ReadConfigFile();
    if (configStatus == -1) {
      std::string msg = "Could not open " + configfile + " or no <config> block found.";
      fail(msg);
    } else if (configStatus < 0) {
      fail("Malformed <potential> or <thm> block in " + configfile + " (see the ERROR above).");
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
    ConfigScope guard(config_);
    return api_->Rebuild();
  }

  bool write_output_files() {
    ConfigScope guard(config_);
    return api_->WriteOutputFiles();
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
  void set_potential(int pair, const std::string &type, bool enabled,
                     double V0, double R, double a, double r0) {
    NuclearPotentialSetting s;
    s.enabled = enabled;
    s.type = type;
    s.V0 = V0;
    s.R = R;
    s.a = a;
    s.r0 = r0;
    NuclearPotentialManager &mgr = NuclearPotentialManager::instance();
    try {
      if (pair)
        mgr.setSetting(pair, s);
      else
        mgr.setDefaultSetting(s);
    } catch (const std::exception &e) {
      throw AZURE2Error(e.what());
    }
    // The master switch follows the pairs: with no pair enabled the hybrid
    // model is off, and CoulFunc short-circuits on it.
    config_->useHybridMethod = mgr.isAnyEnabled();
    if (config_->useHybridMethod)
      config_->paramMask |= Config::USE_HYBRID_COULOMB;
    else
      config_->paramMask &= ~Config::USE_HYBRID_COULOMB;
  }

  void clear_potential(int pair) {
    NuclearPotentialManager &mgr = NuclearPotentialManager::instance();
    if (pair)
      mgr.clearPairSetting(pair);
    else
      mgr.resetToDefault();
    config_->useHybridMethod = mgr.isAnyEnabled();
    if (config_->useHybridMethod)
      config_->paramMask |= Config::USE_HYBRID_COULOMB;
    else
      config_->paramMask &= ~Config::USE_HYBRID_COULOMB;
  }

  // (enabled, type, V0, R, a, r0, has_own_setting) for the given pair,
  // resolved against the default.
  py::tuple get_potential(int pair) const {
    const NuclearPotentialManager &mgr = NuclearPotentialManager::instance();
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

  // -- segment updates (return the number of segments) ------------------------

  int update_segments(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->UpdateSegments(v);
  }
  int update_segments_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->UpdateSegmentsRWA(v);
  }
  int update_segments_all_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->UpdateSegmentsAllRWA(v);
  }

  // -- transforms ------------------------------------------------------------

  py::array_t<double> transform_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_);
      out = api_->TransformRWAParameters(v);
    }
    return to_array(out);
  }
  py::array_t<double> transform_all_rwa(py::array_t<double, py::array::forcecast> p, bool include_fixed) {
    vector_r v = to_vector(p), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_);
      out = api_->TransformAllRWAParameters(v, include_fixed);
    }
    return to_array(out);
  }

  // -- chi-squared and derivatives -------------------------------------------

  double calculate_chi2_rwa(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->CalculateChi2RWA(v);
  }
  double calculate_chi2_physical(const py::array_t<double, py::array::forcecast> &p) {
    ConfigScope guard(config_);
    vector_r v = to_vector(p);
    return api_->CalculateChi2Physical(v);
  }
  py::array_t<double> calculate_chi2_grad_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_);
      out = api_->CalculateChi2GradRWA(v);
    }
    return to_array(out);
  }
  py::array_t<double> calculate_residual_jacobian_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_);
      out = api_->CalculateResidualJacobianRWA(v);
    }
    return to_array(out);
  }
  py::array_t<double> calculate_residuals_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_);
      out = api_->CalculateResidualsRWA(v);
    }
    return to_array(out);
  }
  py::array_t<double> current_norms() {
    ConfigScope guard(config_);
    return to_array(api_->GetCurrentNorms());
  }
  // The THM experiments as the last chi-squared evaluation profiled them.
  py::list thm_experiments() {
    ConfigScope guard(config_);
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
    ConfigScope guard(config_);
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
    ConfigScope guard(config_);
    ThmVertexReport r;
    std::string why;
    if (!api_->GetThmVertex(name, to_vector(e), r, why)) throw AZURE2Error(why);
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
    ConfigScope guard(config_);
    ThmDistortionReport r;
    std::string why;
    if (!api_->GetThmDistortion(name, to_vector(e), r, why)) throw AZURE2Error(why);
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
  py::array_t<double> calculate_model_gradients_rwa(py::array_t<double, py::array::forcecast> p) {
    vector_r v = to_vector(p), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_);
      out = api_->CalculateModelGradientsRWA(v);
    }
    return to_array(out);
  }

  // -- external region and caches --------------------------------------------

  py::array_t<double> coulomb_functions(py::array_t<double, py::array::forcecast> request) {
    vector_r v = to_vector(request), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_);
      out = api_->GetCoulombFunctions(v);
    }
    return to_array(out);
  }
  py::array_t<double> ec_integrals(py::array_t<double, py::array::forcecast> request) {
    vector_r v = to_vector(request), out;
    {
      py::gil_scoped_release release;
      ConfigScope guard(config_);
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
  Config *config_;
  AZUREAPI *api_;
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
      .def("transform_all_rwa", &Session::transform_all_rwa, py::arg("params"),
           py::arg("include_fixed") = false)
      .def("calculate_chi2_rwa", &Session::calculate_chi2_rwa,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("calculate_chi2_physical", &Session::calculate_chi2_physical,
           py::call_guard<py::gil_scoped_release>(), py::arg("params"))
      .def("calculate_chi2_grad_rwa", &Session::calculate_chi2_grad_rwa,
           py::arg("params"))
      .def("calculate_residual_jacobian_rwa", &Session::calculate_residual_jacobian_rwa,
           py::arg("params"))
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
