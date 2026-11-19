// azure2::front's physics (Engine.h, "Physics"; docs/upgrade/physics_tab.md section 7): the rows of the physics and
// stopping tables (api/physics_table, the rows pyazr.physics and pyazr.stopping read) and one evaluated element-wise;
// a PhysicsSession, the engine's session of a saved project (AZUREAPI, set up as pyazr's Session and front::ReadParameters set it up), whose calls
// hand on what the session computes: the channel functions (AZUREAPI::GetCoulombFunctions, pyazr's
// Session.channel_functions), a J^pi block (AZUREAPI::JpiBlockAllRWA, pyazr's Session.jpi_block), the observed
// levels (AZUREAPI::ObservedLevelRowsAllRWA) and the THM tables (AZUREAPI::GetThm*, pyazr's session_thm.py).  No
// formula here: the arrays are copied into series.
#include "azure2/front/Engine.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>

#include <gsl/gsl_errno.h>

#include "AChannel.h"
#include "ALevel.h"
#include "AZUREAPI.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "GSLException.h"
#include "JGroup.h"
#include "LevelState.h"
#include "PPair.h"
#include "ThmExperiment.h"
#include "azure2/physics/LevelParameters.h"
#include "azure2/physics_table/PhysicsTable.h"
#include "azure2/physics_table/StoppingTable.h"
#include "azure2/stopping/Stopping.h"
#include "azure2/support/FloatingPointPrecision.h"

namespace azure2 {
namespace front {

namespace table = azure2::physics_table;

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

std::string Text(double v) {
  char b[64];
  std::snprintf(b, sizeof(b), "%.10g", v);
  return b;
}

/// "1/2", "3": a spin in words.
std::string Spin(double s) {
  const long twice = std::lround(2.0 * s);
  return twice % 2 ? std::to_string(twice) + "/2" : std::to_string(twice / 2);
}

PhysicsArgument ArgumentOf(const table::ArgSpec &s) {
  PhysicsArgument a;
  a.name = s.name, a.unit = s.unit, a.doc = s.doc;
  a.kind = table::KindName(s.kind);
  a.domain = table::DomainText(s);
  a.lower = s.lo, a.upper = s.hi, a.lowerOpen = s.lo_open, a.upperOpen = s.hi_open;
  a.choices = s.choices;
  a.hasDefault = s.has_default, a.defaultValue = s.default_value;
  a.hasExample = s.has_example, a.example = s.example_value;
  a.hasSweep = s.has_sweep, a.sweepLower = s.sweep_lo, a.sweepUpper = s.sweep_hi;
  a.sweepLog = s.sweep_spacing == table::Spacing::log;
  return a;
}

PhysicsSeries Series(const std::string &name, const std::string &unit, const std::string &label,
                     std::vector<std::pair<std::string, double>> keys, std::size_t n) {
  PhysicsSeries s;
  s.name = name, s.unit = unit, s.label = label, s.keys = std::move(keys);
  s.values.assign(n, kNaN);
  return s;
}

const char *KindWord(unsigned char kind) {
  using azure2::physics::ObservedKind;
  switch (static_cast<ObservedKind>(kind)) {
    case ObservedKind::Width: return "width";
    case ObservedKind::ANC: return "anc";
    case ObservedKind::MagneticMoment: return "magnetic_moment";
    case ObservedKind::QuadrupoleMoment: return "quadrupole_moment";
    case ObservedKind::Beta: return "beta";
  }
  return "width";
}

std::string PairKind(int ptype) { return ptype == 0 ? "particle" : ptype == 10 ? "photon" : "beta"; }

}  // namespace

// ---------------------------------------------------------------- the physics table

std::vector<PhysicsFunction> PhysicsFunctions() {
  std::vector<PhysicsFunction> out;
  for (const table::PhysicsRow &r : table::PhysicsTable()) {
    PhysicsFunction f;
    f.name = r.name, f.summary = r.summary, f.doc = r.doc, f.library = r.library, f.axis = r.axis, f.group = r.group;
    for (const table::ArgSpec &s : r.args) f.arguments.push_back(ArgumentOf(s));
    for (const table::OutSpec &s : r.outs) f.returns.push_back({s.name, s.unit, s.doc});
    out.push_back(f);
  }
  return out;
}

bool EvaluatePhysics(const std::string &name, const std::vector<std::vector<double>> &arguments, PhysicsValues &out,
                     std::string &why) {
  out = PhysicsValues();
  const table::PhysicsRow *row = table::FindPhysicsRow(name);
  if (!row) {
    why = "no physics function named '" + name + "' (PhysicsFunctions() lists them)";
    return false;
  }
  const std::size_t na = row->args.size(), no = row->outs.size();
  if (arguments.size() != na) {
    why = name + ": " + std::to_string(na) + " argument columns expected, " + std::to_string(arguments.size()) +
          " given";
    return false;
  }
  std::size_t n = na ? 1 : 0;
  for (const std::vector<double> &c : arguments)
    if (c.size() != 1) {
      if (n != 1 && c.size() != n) {
        why = name + ": the argument columns must have one length (or length 1, repeated)";
        return false;
      }
      n = c.size();
    }
  for (const std::vector<double> &c : arguments)
    if (c.empty()) n = 0;
  out.values.assign(no, std::vector<double>(n, kNaN));
  out.status.assign(n, 0);
  const azure2::support::ExtendedPrecisionScope fpPrecision;
  std::vector<double> a(na), o(no);
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t k = 0; k < na; k++) a[k] = arguments[k].size() == 1 ? arguments[k][0] : arguments[k][i];
    out.status[i] = table::EvaluateElement(*row, a.data(), o.data());
    for (std::size_t j = 0; j < no; j++) out.values[j][i] = out.status[i] == table::kOk ? o[j] : kNaN;
  }
  return true;
}

std::string PhysicsStatusText(const std::string &name, int status) {
  const table::PhysicsRow *row = table::FindPhysicsRow(name);
  if (status == table::kOk) return "computed";
  if (status == table::kFailed) return "not converged or failed";
  if (status == table::kJointly) return "the arguments together are outside the domain";
  if (row && status >= 1 && status <= (int)row->args.size()) {
    const table::ArgSpec &s = row->args[status - 1];
    return s.name + " outside its domain (" + table::DomainText(s) + ")";
  }
  return "status " + std::to_string(status) + " (the function's doc says what it means)";
}

// ---------------------------------------------------------------- the stopping table

std::vector<PhysicsFunction> StoppingFunctions() {
  std::vector<PhysicsFunction> out;
  for (const table::StoppingRow &r : table::StoppingTable()) {
    PhysicsFunction f;
    f.name = r.name, f.summary = r.summary, f.doc = r.doc, f.library = r.library, f.axis = r.axis;
    f.group = "stopping";
    f.arguments.push_back(ArgumentOf(table::StoppingFormula()));
    for (const table::ArgSpec &s : table::StoppingNumericArgs(r)) f.arguments.push_back(ArgumentOf(s));
    for (const table::OutSpec &s : r.outs) f.returns.push_back({s.name, s.unit, s.doc});
    out.push_back(f);
  }
  return out;
}

std::string StoppingCitation() { return azure2::stopping::kSrimCitation; }

bool EvaluateStopping(const std::string &name, const std::string &formula,
                      const std::vector<std::vector<double>> &arguments, PhysicsValues &out, std::string &why) {
  out = PhysicsValues();
  const table::StoppingRow *row = table::FindStoppingRow(name);
  if (!row) {
    why = "no stopping function named '" + name + "' (StoppingFunctions() lists them)";
    return false;
  }
  const std::size_t na = table::StoppingNumericArgs(*row).size(), no = row->outs.size();
  if (arguments.size() != na) {
    why = name + ": " + std::to_string(na) + " numeric argument columns expected (after the formula), " +
          std::to_string(arguments.size()) + " given";
    return false;
  }
  std::size_t n = 1;
  for (const std::vector<double> &c : arguments)
    if (c.size() != 1) {
      if (n != 1 && c.size() != n) {
        why = name + ": the argument columns must have one length (or length 1, repeated)";
        return false;
      }
      n = c.size();
    }
  for (const std::vector<double> &c : arguments)
    if (c.empty()) n = 0;
  const azure2::support::ExtendedPrecisionScope fpPrecision;
  try {
    table::StoppingMaterial(formula, "", 1.0);  // a formula the library refuses is refused here, before any element
  } catch (const std::exception &e) {
    why = name + ": the formula '" + formula + "' is refused: " + e.what();
    return false;
  }
  std::vector<std::vector<double>> full(na);
  std::vector<const double *> columns(na);
  for (std::size_t k = 0; k < na; k++) {
    full[k] = arguments[k].size() == 1 ? std::vector<double>(n, arguments[k][0]) : arguments[k];
    columns[k] = full[k].data();
  }
  std::vector<double> values(no * n);
  out.status.assign(n, 0);
  table::EvaluateStoppingRow(*row, formula, columns, n, values.data(), out.status.data());
  out.values.assign(no, std::vector<double>(n));
  for (std::size_t j = 0; j < no; j++)
    for (std::size_t i = 0; i < n; i++) out.values[j][i] = values[j * n + i];
  return true;
}

double PhysicsSeries::Key(const std::string &key, double absent) const {
  for (const auto &k : keys)
    if (k.first == key) return k.second;
  return absent;
}

const PhysicsSeries *PhysicsCurves::Find(const std::string &name,
                                         const std::vector<std::pair<std::string, double>> &match) const {
  for (const PhysicsSeries &s : series) {
    if (s.name != name) continue;
    bool all = true;
    for (const auto &m : match) all = all && s.Key(m.first, kNaN) == m.second;
    if (all) return &s;
  }
  return nullptr;
}

// ---------------------------------------------------------------- the session

struct PhysicsSession::Impl {
  std::ostringstream notes;  // the session's Config writes here (declared before it: it outlives the Config)
  std::unique_ptr<Config> config;
  std::unique_ptr<AZUREAPI> api;
  vector_r x;  // every parameter, fixed ones included (params_all_rwa layout)
};

PhysicsSession::PhysicsSession() : impl_(new Impl) {}
PhysicsSession::~PhysicsSession() = default;
PhysicsSession::PhysicsSession(PhysicsSession &&) noexcept = default;
PhysicsSession &PhysicsSession::operator=(PhysicsSession &&) noexcept = default;

bool PhysicsSession::IsOpen() const { return impl_ && impl_->api; }

bool PhysicsSession::Open(const std::string &projectFile, const std::string &parameterFile, bool extrapolations,
                          PhysicsSession &out, std::string &why) {
  out = PhysicsSession();
  Impl &s = *out.impl_;
  const azure2::support::ExtendedPrecisionScope fpPrecision;
  gsl_set_error_handler(&GSLException::GSLErrorHandler);
  // A session as pyazr's Session and front::ReadParameters make it: the project read as the program reads it, the
  // matrix route resolved, the model and the data (or the extrapolations) set up.
  s.config.reset(new Config(s.notes));
  Config &config = *s.config;
  config.configfile = projectFile;
  config.paramMask |= Config::USE_API;
  if (extrapolations)
    config.paramMask &= ~Config::CALCULATE_WITH_DATA;
  else
    config.paramMask |= Config::CALCULATE_WITH_DATA;
  const int status = config.ReadConfigFile();
  if (status < 0) {
    why = projectFile + " cannot be set up (status " + std::to_string(status) + "): " + s.notes.str();
    out = PhysicsSession();
    return false;
  }
  config.ResolveMatrixRoute();
  try {
    InitializeECAmplitudeCache(config);
    s.api.reset(new AZUREAPI(config));
    if (s.api->Initialize() != 0) {
      why = "the engine could not set up " + projectFile + ": " + s.notes.str();
      out = PhysicsSession();
      return false;
    }
    s.x = s.api->params_all_rwa();
    if (!parameterFile.empty()) {
      // The program's --parameter-file (src/ProgramMain.cpp): the vector as the set-up makes it, then the file read
      // into it by AZUREParams::ReadUserParameters and its basis reconciled.
      AZUREParams params;
      s.api->compound()->FillMnParams(params, s.api->state(), s.api->inputs(), &config);
      s.api->data()->FillMnParams(params);
      config.paramfile = parameterFile;
      if (!params.ReadUserParameters(config) || !params.ReconcileBasis(s.api->compound(), config)) {
        why = parameterFile + " cannot be applied to " + projectFile + ": " + s.notes.str();
        out = PhysicsSession();
        return false;
      }
      const std::vector<double> values = params.GetMinuitParams().Params();
      if (values.size() != s.x.size() || !s.api->SetParamsAllRWA(vector_r(values.begin(), values.end()))) {
        why = parameterFile + ": the parameter vector does not match the session's (" + std::to_string(values.size()) +
              " against " + std::to_string(s.x.size()) + " values)";
        out = PhysicsSession();
        return false;
      }
      s.x = s.api->params_all_rwa();
    }
  } catch (const std::exception &e) {
    why = std::string("the engine could not set up the project: ") + e.what();
    out = PhysicsSession();
    return false;
  }
  return true;
}

std::string PhysicsSession::Parametrization() const {
  if (!IsOpen()) return "";
  const Config &c = impl_->api->configure();
  if (c.paramMask & Config::USE_PARK_FORMALISM) return "park";
  if (c.paramMask & Config::USE_BRUNE_FORMALISM) return "brune";
  return "standard";
}

std::vector<PhysicsPair> PhysicsSession::Pairs() const {
  std::vector<PhysicsPair> out;
  if (!IsOpen()) return out;
  const CompoundSystem *c = impl_->api->compound();
  for (int i = 1; i <= c->NumPairs(); i++) {
    const PPair *p = c->GetPair(i);
    PhysicsPair q;
    q.key = p->GetPairKey();
    q.kind = PairKind(p->GetPType());
    q.z1 = p->GetZ(1), q.z2 = p->GetZ(2);
    q.m1 = p->GetM(1), q.m2 = p->GetM(2);
    q.j1 = p->GetJ(1), q.j2 = p->GetJ(2);
    q.separation = p->GetSepE(), q.excitation = p->GetExE();
    q.radius = p->GetPairRadius();
    q.identical = p->IsIdentical();
    out.push_back(q);
  }
  return out;
}

std::vector<PhysicsGroup> PhysicsSession::Groups() const {
  std::vector<PhysicsGroup> out;
  if (!IsOpen()) return out;
  const CompoundSystem *c = impl_->api->compound();
  for (int j = 1; j <= c->NumJGroups(); j++) {
    const JGroup *g = c->GetJGroup(j);
    PhysicsGroup q;
    q.jgroup = j, q.J = g->GetJ(), q.parity = g->GetPi();
    for (int la = 1; la <= g->NumLevels(); la++) q.levels += g->GetLevel(la)->IsInRMatrix() ? 1 : 0;
    for (int ch = 1; ch <= g->NumChannels(); ch++) {
      const AChannel *a = g->GetChannel(ch);
      const PPair *p = c->GetPair(a->GetPairNum());
      PhysicsChannel k;
      k.channel = ch, k.pair = p->GetPairKey(), k.radiation = std::string(1, a->GetRadType());
      k.l = a->GetL(), k.s = a->GetS();
      k.radius = p->GetPType() == 0 ? a->GetRadius() : 0.0;
      q.channels.push_back(k);
    }
    out.push_back(q);
  }
  return out;
}

std::vector<PhysicsPole> PhysicsSession::Poles() {
  std::vector<PhysicsPole> out;
  if (!IsOpen()) return out;
  const azure2::support::ExtendedPrecisionScope fpPrecision;
  AZUREAPI &api = *impl_->api;
  const vector_r rows = api.ObservedLevelRowsAllRWA(impl_->x);  // fills the state at x
  const CompoundSystem *c = api.compound();
  for (std::size_t r = 0; r + 10 <= rows.size(); r += 10) {
    const int j = (int)rows[r], la = (int)rows[r + 1], ch = (int)rows[r + 2];
    const JGroup *g = c->GetJGroup(j);
    if (out.empty() || out.back().jgroup != j || out.back().level != la) {
      PhysicsPole p;
      p.jgroup = j, p.level = la, p.J = g->GetJ(), p.parity = g->GetPi();
      p.energy = rows[r + 3];
      p.formalEnergy = api.state().E(j, la);
      out.push_back(p);
    }
    const AChannel *a = g->GetChannel(ch);
    const PPair *pair = c->GetPair(a->GetPairNum());
    PhysicsPole::Width w;
    w.channel = ch, w.pair = pair->GetPairKey(), w.radiation = std::string(1, a->GetRadType());
    w.l = a->GetL(), w.s = a->GetS();
    w.kind = KindWord((unsigned char)rows[r + 4]);
    w.width = w.kind == std::string("width") ? rows[r + 5] : 0.0;
    w.value = rows[r + 6];
    w.penetrability = rows[r + 7];
    out.back().totalWidth += std::fabs(w.width);
    out.back().widths.push_back(w);
  }
  return out;
}

bool PhysicsSession::ChannelFunctions(int pair, int l, double radius, int twoJ, const std::vector<double> &eCm,
                                      PhysicsCurves &out, std::string &why) {
  out = PhysicsCurves();
  if (!IsOpen()) {
    why = "no session";
    return false;
  }
  AZUREAPI &api = *impl_->api;
  const CompoundSystem *c = api.compound();
  if (!c->IsPairKey(pair) || c->GetPair(c->GetPairNumFromKey(pair))->GetPType() != 0) {
    why = "pair " + std::to_string(pair) + " is no particle pair of the project (Pairs() lists them)";
    return false;
  }
  if (l < 0) {
    why = "l must be a whole number >= 0";
    return false;
  }
  vector_r request = {(double)pair, (double)l, (double)twoJ, radius, (double)eCm.size()};
  request.insert(request.end(), eCm.begin(), eCm.end());
  vector_r response;
  try {
    const azure2::support::ExtendedPrecisionScope fpPrecision;
    response = api.GetCoulombFunctions(request);
  } catch (const std::exception &e) {
    why = e.what();
    return false;
  }
  if (response.size() != 1 + 7 * eCm.size()) {
    why = "the engine gave no channel functions for pair " + std::to_string(pair);
    return false;
  }
  out.axis = "e_cm", out.axisUnit = "MeV", out.axisValues = eCm;
  const std::string label = "pair " + std::to_string(pair) + ", l " + std::to_string(l);
  const std::vector<std::pair<std::string, double>> keys = {{"pair", pair}, {"l", l}};
  const char *names[7] = {"F", "Fp", "G", "Gp", "P", "S", "delta_hs"};
  const char *units[7] = {"", "", "", "", "", "", "rad"};
  for (int k = 0; k < 7; k++) out.series.push_back(Series(names[k], units[k], label, keys, eCm.size()));
  PhysicsSeries open = Series("open", "", label, keys, eCm.size());
  for (std::size_t i = 0; i < eCm.size(); i++) {
    for (int k = 0; k < 7; k++) out.series[k].values[i] = response[1 + 7 * i + k];
    // as pyazr's channel_functions: E > 0 and the engine's functions defined there
    open.values[i] = eCm[i] > 0 && (response[1 + 7 * i] != 0 || response[1 + 7 * i + 2] != 0) ? 1.0 : 0.0;
  }
  out.series.push_back(open);
  return true;
}

bool PhysicsSession::Block(int jgroup, const std::vector<double> &excitation, PhysicsCurves &out, std::string &why) {
  out = PhysicsCurves();
  if (!IsOpen()) {
    why = "no session";
    return false;
  }
  JpiBlockReport r;
  try {
    const azure2::support::ExtendedPrecisionScope fpPrecision;
    r = impl_->api->JpiBlockAllRWA(impl_->x, jgroup, vector_r(excitation.begin(), excitation.end()));
  } catch (const std::exception &e) {
    why = e.what();
    return false;
  }
  const std::size_t ne = r.ex.size(), nc = r.channels.size(), n = r.levels.size();
  out.axis = "E_x", out.axisUnit = "MeV", out.axisValues = r.ex;
  out.facts.push_back({"J", Spin(r.J) + (r.pi < 0 ? "-" : "+")});
  out.facts.push_back({"parametrization", r.park ? "park" : r.brune ? "brune" : "standard"});
  for (std::size_t i = 0; i < n; i++)
    out.facts.push_back({"E_level " + std::to_string(r.levels[i]), Text(r.levelEnergy[i]) + " MeV"});
  const double jg = jgroup;
  auto channelLabel = [&](std::size_t c) {
    const JpiBlockReport::Channel &k = r.channels[c];
    return "channel " + std::to_string(k.channel) + " (pair " + std::to_string(k.pairKey) + ", l " +
           std::to_string(k.l) + ", s " + Spin(k.s) + ")";
  };
  for (std::size_t c = 0; c < nc; c++) {
    const JpiBlockReport::Channel &k = r.channels[c];
    const std::vector<std::pair<std::string, double>> keys = {
        {"jgroup", jg}, {"channel", (double)k.channel}, {"pair", (double)k.pairKey}, {"l", (double)k.l}, {"s", k.s}};
    PhysicsSeries ec = Series("E_c", "MeV", channelLabel(c), keys, ne), P = Series("P", "", channelLabel(c), keys, ne),
                  S = Series("S", "", channelLabel(c), keys, ne), Lre = Series("L_re", "", channelLabel(c), keys, ne),
                  Lim = Series("L_im", "", channelLabel(c), keys, ne);
    for (std::size_t e = 0; e < ne; e++) {
      ec.values[e] = r.ec[e * nc + c];
      if (k.radType != 'P') continue;
      P.values[e] = r.P[e * nc + c], S.values[e] = r.S[e * nc + c];
      Lre.values[e] = r.lo[e * nc + c].real(), Lim.values[e] = r.lo[e * nc + c].imag();
    }
    for (PhysicsSeries *s : {&ec, &P, &S, &Lre, &Lim}) out.series.push_back(*s);
  }
  for (std::size_t c = 0; c < nc; c++)
    for (std::size_t cp = c; cp < nc; cp++) {
      if (r.channels[c].radType != 'P' || r.channels[cp].radType != 'P') continue;
      const std::vector<std::pair<std::string, double>> keys = {
          {"jgroup", jg}, {"channel", (double)r.channels[c].channel}, {"channel_out", (double)r.channels[cp].channel}};
      const std::string label = channelLabel(c) + " -> " + channelLabel(cp);
      PhysicsSeries Ure = Series("U_re", "", label, keys, ne), Uim = Series("U_im", "", label, keys, ne),
                    Uabs = Series("U_abs", "", label, keys, ne), Tre = Series("T_re", "", label, keys, ne),
                    Tim = Series("T_im", "", label, keys, ne), delta = Series("delta", "deg", label, keys, ne),
                    eta = Series("eta", "", label, keys, ne);
      for (std::size_t e = 0; e < ne; e++) {
        const std::complex<double> U = r.U[e * nc * nc + c * nc + cp], T = r.T[e * nc * nc + c * nc + cp];
        if (std::isnan(U.real())) continue;
        Ure.values[e] = U.real(), Uim.values[e] = U.imag(), Uabs.values[e] = std::abs(U);
        Tre.values[e] = T.real(), Tim.values[e] = T.imag();
        if (c == cp) delta.values[e] = r.delta[e * nc + c], eta.values[e] = r.eta[e * nc + c];
      }
      for (PhysicsSeries *s : {&Ure, &Uim, &Uabs, &Tre, &Tim}) out.series.push_back(*s);
      if (c == cp) {
        out.series.push_back(delta);
        out.series.push_back(eta);
      }
    }
  for (std::size_t l = 0; l < n; l++)
    for (std::size_t m = l; m < n; m++) {
      const std::vector<std::pair<std::string, double>> keys = {
          {"jgroup", jg}, {"level", (double)r.levels[l]}, {"level_out", (double)r.levels[m]}};
      PhysicsSeries a = Series("A_abs", "MeV^-1",
                               "level " + std::to_string(r.levels[l]) + ", level " + std::to_string(r.levels[m]), keys,
                               ne);
      for (std::size_t e = 0; e < ne; e++) a.values[e] = std::abs(r.A[e * n * n + l * n + m]);
      out.series.push_back(a);
    }
  const std::vector<std::pair<std::string, double>> keys = {{"jgroup", jg}};
  PhysicsSeries dre = Series("det_re", "MeV^n", "det A^-1", keys, ne), dim = Series("det_im", "MeV^n", "det A^-1", keys, ne),
                dlog = Series("log10_abs_det", "", "log10 |det A^-1|", keys, ne);
  for (std::size_t e = 0; e < ne; e++) {
    dre.values[e] = r.detInverse[e].real(), dim.values[e] = r.detInverse[e].imag();
    dlog.values[e] = std::log10(std::abs(r.detInverse[e]));
    if (r.singular[e]) out.notes.push_back("A^-1 singular at E_x = " + Text(r.ex[e]) + " MeV (A, U, T set to 0)");
  }
  out.series.push_back(dre);
  out.series.push_back(dim);
  out.series.push_back(dlog);
  for (std::size_t c = 0; c < nc; c++)
    if (r.channels[c].radType != 'P')
      out.notes.push_back(channelLabel(c) + " is a photon or beta channel: it enters A (R-matrix capture only) but has "
                          "no U or T element here");
  return true;
}

std::vector<std::string> PhysicsSession::ThmExperiments() const {
  std::vector<std::string> out;
  if (!IsOpen()) return out;
  for (const ThmExperiment &e : impl_->api->configure().thm.experiments) out.push_back(e.name);
  return out;
}

namespace {

/// The session's chi-squared at its vector, before the THM tables that read the last evaluation (pyazr's
/// session_thm.py does the same: calculate_chi2_rwa(x) first).
void EvaluateAtVector(AZUREAPI &api) { api.CalculateChi2RWA(api.params_values_rwa()); }

}  // namespace

bool PhysicsSession::ThmVertex(const std::string &experiment, const std::vector<double> &energies, PhysicsCurves &out,
                               std::string &why) {
  out = PhysicsCurves();
  if (!IsOpen()) {
    why = "no session";
    return false;
  }
  ThmVertexReport r;
  try {
    const azure2::support::ExtendedPrecisionScope fpPrecision;
    EvaluateAtVector(*impl_->api);
    if (!impl_->api->GetThmVertex(experiment, vector_r(energies.begin(), energies.end()), r, why)) return false;
  } catch (const std::exception &e) {
    why = e.what();
    return false;
  }
  const std::size_t ne = r.energy.size();
  out.axis = "E", out.axisUnit = "MeV", out.axisValues = r.energy;
  out.facts = {{"experiment", r.experiment}, {"window", r.window}, {"model", r.model},
               {"pair", std::to_string(r.pairKey)}, {"B", Text(r.bind) + " MeV"}, {"radius", Text(r.radius) + " fm"}};
  PhysicsSeries reached = Series("reached", "", "reached", {}, ne);
  for (std::size_t e = 0; e < ne && e < r.reached.size(); e++) reached.values[e] = r.reached[e] ? 1.0 : 0.0;
  out.series.push_back(reached);
  for (const ThmVertexReport::Channel &c : r.channels)
    for (const ThmVertexReport::Level &l : c.levels) {
      const std::vector<std::pair<std::string, double>> keys = {
          {"jgroup", (double)c.jgroup}, {"channel", (double)c.channel}, {"J", c.J}, {"pi", (double)c.pi},
          {"l", (double)c.l}, {"s", c.s}, {"level", (double)l.level}};
      const std::string label = "J " + Spin(c.J) + (c.pi < 0 ? "-" : "+") + ", channel " + std::to_string(c.channel) +
                                ", l " + std::to_string(c.l) + ", level " + std::to_string(l.level);
      const std::pair<const char *, const std::vector<double> *> fields[3] = {
          {"M2", &l.m2}, {"M2_qf", &l.m2qf}, {"M2_pw", &l.m2pw}};
      for (const auto &f : fields) {
        if (f.second->empty()) continue;
        PhysicsSeries s = Series(f.first, "", label, keys, ne);
        for (std::size_t e = 0; e < ne && e < f.second->size(); e++) s.values[e] = (*f.second)[e];
        out.series.push_back(s);
      }
    }
  for (std::size_t e = 0; e < ne && e < r.reached.size(); e++)
    if (!r.reached[e]) {
      out.notes.push_back("the window does not reach some energies: the vertex there is the nearest data point's");
      break;
    }
  return true;
}

bool PhysicsSession::ThmLineshape(const std::string &experiment, const std::vector<double> &energies,
                                  PhysicsCurves &out, std::string &why) {
  out = PhysicsCurves();
  if (!IsOpen()) {
    why = "no session";
    return false;
  }
  ThmLineshapeReport r;
  try {
    const azure2::support::ExtendedPrecisionScope fpPrecision;
    EvaluateAtVector(*impl_->api);
    if (!impl_->api->GetThmLineshape(experiment, vector_r(energies.begin(), energies.end()), r, why)) return false;
  } catch (const std::exception &e) {
    why = e.what();
    return false;
  }
  const std::size_t ne = r.energy.size();
  out.axis = "E", out.axisUnit = "MeV", out.axisValues = r.energy;
  out.facts = {{"experiment", r.experiment}, {"spectator", r.spectator}, {"E_aA", Text(r.eAA) + " MeV"},
               {"B", Text(r.bind) + " MeV"}, {"Zs", std::to_string(r.Zs)}, {"ZF", std::to_string(r.ZF)}};
  PhysicsSeries esf = Series("E_sF", "MeV", "E_sF", {}, ne), eta0 = Series("eta_0", "", "eta_0", {}, ne);
  for (std::size_t e = 0; e < ne; e++) {
    if (e < r.esf.size()) esf.values[e] = r.esf[e];
    if (e < r.eta0.size()) eta0.values[e] = r.eta0[e];
  }
  out.series.push_back(esf);
  out.series.push_back(eta0);
  for (const ThmLineshapeReport::Exit &x : r.exits) {
    const std::string pl = "exit pair " + std::to_string(x.pairKey);
    PhysicsSeries zeta = Series("zeta", "", pl, {{"pair", (double)x.pairKey}}, ne),
                  etaSb = Series("eta_sb", "", pl, {{"pair", (double)x.pairKey}}, ne);
    for (std::size_t e = 0; e < ne; e++) {
      if (e < x.zeta.size()) zeta.values[e] = x.zeta[e];
      if (e < x.etaSb.size()) etaSb.values[e] = x.etaSb[e];
    }
    out.series.push_back(zeta);
    out.series.push_back(etaSb);
    for (const ThmLineshapeReport::Level &l : x.levels) {
      PhysicsSeries s = Series("NC2", "",
                               pl + ", J " + Spin(l.J) + (l.pi < 0 ? "-" : "+") + ", level " + std::to_string(l.level) +
                                   " (E " + Text(l.energy) + " MeV, Gamma " + Text(l.width) + " MeV)",
                               {{"pair", (double)x.pairKey}, {"jgroup", (double)l.jgroup}, {"level", (double)l.level},
                                {"J", l.J}, {"pi", (double)l.pi}},
                               ne);
      for (std::size_t e = 0; e < ne && e < l.nc2.size(); e++) s.values[e] = l.nc2[e];
      out.series.push_back(s);
    }
  }
  return true;
}

bool PhysicsSession::ThmDistortion(const std::string &experiment, const std::vector<double> &energies,
                                   PhysicsCurves &out, std::string &why) {
  out = PhysicsCurves();
  if (!IsOpen()) {
    why = "no session";
    return false;
  }
  ThmDistortionReport r;
  try {
    const azure2::support::ExtendedPrecisionScope fpPrecision;
    if (!impl_->api->GetThmDistortion(experiment, vector_r(energies.begin(), energies.end()), r, why)) return false;
  } catch (const std::exception &e) {
    why = e.what();
    return false;
  }
  const std::size_t ne = r.energy.size();
  out.axis = "E", out.axisUnit = "MeV", out.axisValues = r.energy;
  out.facts = {{"experiment", r.experiment}, {"kind", r.kind},
               {"description", r.description}, {"E_ref", Text(r.eRef) + " MeV"},
               {"E_aA", Text(r.eAA) + " MeV"}, {"B", Text(r.bind) + " MeV"},
               {"k_aA", Text(r.kAA) + " fm^-1"}, {"eta_aA", Text(r.etaAA)},
               {"kappa", Text(r.kappa) + " fm^-1"}, {"eta_b", Text(r.etaB)}, {"beta", Text(r.beta)}};
  const std::tuple<const char *, const char *, const std::vector<double> *> fields[] = {
      {"R", "", &r.r},          {"R_model", "", &r.rModel}, {"M2", "fm^6", &r.m2},     {"M2_PW", "fm^6", &r.mpw2},
      {"E_sF", "MeV", &r.esf},  {"k_sF", "fm^-1", &r.ksf},  {"eta_sF", "", &r.etasf},  {"theta_cm", "deg", &r.thetaCm},
      {"x", "", &r.x},          {"q", "fm^-1", &r.q}};
  for (const auto &f : fields) {
    PhysicsSeries s = Series(std::get<0>(f), std::get<1>(f), std::get<0>(f), {}, ne);
    for (std::size_t e = 0; e < ne && e < std::get<2>(f)->size(); e++) s.values[e] = (*std::get<2>(f))[e];
    out.series.push_back(s);
  }
  PhysicsSeries lmax = Series("lmax", "", "lmax", {}, ne);
  for (std::size_t e = 0; e < ne && e < r.lmax.size(); e++) lmax.values[e] = r.lmax[e];
  out.series.push_back(lmax);
  return true;
}

}  // namespace front
}  // namespace azure2
