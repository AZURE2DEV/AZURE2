// azure2::front's fit parameters (Engine.h; docs/upgrade/qt_gui.md section 7): the parameter vector of a session of the
// project (AZUREAPI, as pyazr's Session sets it up: names, values, fixed flags, limits, records), a parameter file as
// the engine reads it, and the rows of <parameterSettings> and of <mcmc> <parameters> in the grammar of their readers
// (ParameterLimitsManager::ReadParameterSettings, src/SamplerRun.cpp).
#include "azure2/front/Engine.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>

#include <gsl/gsl_errno.h>

#include "AZUREAPI.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "GSLException.h"
#include "azure2/support/FloatingPointPrecision.h"

namespace azure2 {
namespace front {

namespace {

std::vector<std::string> Tokens(const std::string &line) {
  std::istringstream s(line);
  std::vector<std::string> out;
  std::string t;
  while (s >> t) out.push_back(t);
  return out;
}

std::string Trim(const std::string &s) {
  const std::size_t a = s.find_first_not_of(" \t\r");
  if (a == std::string::npos) return std::string();
  return s.substr(a, s.find_last_not_of(" \t\r") - a + 1);
}

std::string Number(double v) {
  char b[64];
  std::snprintf(b, sizeof(b), "%.10g", v);
  return b;
}

double ToDouble(const std::string &s) {
  char *end = nullptr;
  const double v = std::strtod(s.c_str(), &end);
  return end && *end == '\0' ? v : std::numeric_limits<double>::quiet_NaN();
}

/// The first line of the file's section `name` (1-based), 0 when it has none.
int SectionLine(const Project &project, const std::string &name) {
  std::istringstream text(project.Text());
  std::string line;
  int n = 0;
  while (std::getline(text, line)) {
    n++;
    if (Trim(line) == "<" + name + ">") return n;
  }
  return 0;
}

/// A row of <parameterSettings> read as ParameterLimitsManager reads it (9, 12 or 14 tokens), or false.
bool ReadSettingRow(const std::string &raw, ParameterSetting &row) {
  const std::string line = Trim(raw);
  if (line.empty() || line[0] == '#' || line[0] == '<') return false;
  const std::vector<std::string> p = Tokens(line);
  std::size_t names = 0;
  if (p.size() == 9) names = 1;
  else if (p.size() == 12) names = 4;
  else if (p.size() == 14) names = 6;
  else return false;
  row = ParameterSetting();
  for (std::size_t i = 0; i < names; i++) row.name += (i ? " " : "") + p[i];
  row.nominal = ToDouble(p[names]);
  row.lower = ToDouble(p[names + 1]);
  row.upper = ToDouble(p[names + 2]);
  row.error = ToDouble(p[names + 3]);
  row.fitError = ToDouble(p[names + 4]);
  row.nuisance = p[names + 5] == "1";
  row.category = p[names + 6];
  row.index = std::atoi(p[names + 7].c_str());
  return true;
}

std::string SettingText(const ParameterSetting &r) {
  return r.name + " " + Number(r.nominal) + " " + Number(r.lower) + " " + Number(r.upper) + " " + Number(r.error) + " " +
         Number(r.fitError) + " " + (r.nuisance ? "1" : "0") + " " + (r.category.empty() ? "level" : r.category) + " " +
         std::to_string(r.index);
}

}  // namespace

bool ReadParameters(const std::string &projectFile, std::vector<Parameter> &out, std::string &why) {
  return ReadParameters(projectFile, std::string(), out, why);
}

bool ReadParameters(const std::string &projectFile, const std::string &scratchDirectory, std::vector<Parameter> &out,
                    std::string &why) {
  out.clear();
  const azure2::support::ExtendedPrecisionScope fpPrecision;
  gsl_set_error_handler(&GSLException::GSLErrorHandler);
  std::ostringstream notes;
  // A session as pyazr's (api/src/pybind_api.cpp, Session): the project read as the program reads it, the matrix
  // route resolved, the model and the data set up.
  Config config(notes);
  config.configfile = projectFile;
  config.paramMask |= Config::USE_API;
  const int status = config.ReadConfigFile();
  if (status < 0) {
    why = projectFile + " cannot be set up (status " + std::to_string(status) + "): " + notes.str();
    return false;
  }
  if (!scratchDirectory.empty()) {  // what the set-up writes (an intEC.dat) goes there, not into the outputs
    config.outputdir = scratchDirectory + (scratchDirectory.back() == '/' ? "" : "/");
    config.checkdir = config.outputdir;
  }
  config.ResolveMatrixRoute();
  try {
    InitializeECAmplitudeCache(config);
    AZUREAPI api(config);
    if (api.Initialize() != 0) {
      why = "the engine could not set up " + projectFile + ": " + notes.str();
      return false;
    }
    const vector_r values = api.params_all_rwa();
    const std::vector<bool> fixed = api.params_fixed();
    const std::vector<ParameterRecord> &records = api.GetParameterRecords();
    for (std::size_t i = 0; i < values.size(); i++) {
      Parameter p;
      p.name = api.params_names((int)i);
      p.value = values[i];
      p.fixed = i < fixed.size() && fixed[i];
      if (i < records.size()) {
        p.kind = records[i].Row().key;
        p.unit = records[i].Row().unit;
        p.levelNumber = records[i].levelNumber;
        p.channel = records[i].channel;
      }
      const std::pair<double, double> bounds = api.GetParameterBounds((int)i);
      p.lower = bounds.first;
      p.upper = bounds.second;
      out.push_back(p);
    }
  } catch (const std::exception &e) {
    why = std::string("the engine could not set up the project: ") + e.what();
    return false;
  }
  return true;
}

bool ReadParameterFile(const std::string &path, std::vector<ParameterValue> &out, std::string &why) {
  out.clear();
  std::ifstream in(path);
  if (!in) {
    why = "cannot read " + path;
    return false;
  }
  std::string line;
  while (std::getline(in, line)) {
    const std::vector<std::string> t = Tokens(line);
    if (t.size() < 2 || t[0][0] == '#') continue;
    ParameterValue v;
    v.name = t[0];
    v.value = ToDouble(t[1]);
    v.error = t.size() > 2 ? ToDouble(t[2]) : std::numeric_limits<double>::quiet_NaN();
    out.push_back(v);
  }
  return true;
}

std::vector<ParameterSetting> ParameterSettings(const Project &project) {
  std::vector<ParameterSetting> out;
  const int first = SectionLine(project, "parameterSettings");
  const std::vector<std::string> body = project.SectionBody("parameterSettings");
  for (std::size_t i = 0; i < body.size(); i++) {
    ParameterSetting row;
    if (!ReadSettingRow(body[i], row)) continue;
    row.line = first ? first + 1 + (int)i : 0;
    out.push_back(row);
  }
  return out;
}

bool SetParameterSetting(Project &project, const ParameterSetting &row, std::string &why, bool remove) {
  const std::vector<std::string> names = Tokens(row.name);
  if (names.empty() || row.name.find('#') != std::string::npos) {
    why = "a parameter setting needs a name without '#'";
    return false;
  }
  // The reader takes 1, 4 ("Level N Energy (MeV)") or 6 ("Level N Channel M Width (eV)") name tokens.
  if (names.size() != 1 && names.size() != 4 && names.size() != 6) {
    why = "\"" + row.name + "\": a parameter setting is named by one token (as param.sav) or the Level forms";
    return false;
  }
  for (double v : {row.nominal, row.lower, row.upper, row.error, row.fitError})
    if (!std::isfinite(v)) {
      why = "\"" + row.name + "\": every number of a parameter setting must be finite (0 0 for no limits)";
      return false;
    }
  std::vector<std::string> body = project.SectionBody("parameterSettings");
  bool found = false;
  for (std::size_t i = 0; i < body.size(); i++) {
    ParameterSetting old;
    if (!ReadSettingRow(body[i], old) || old.name != Trim(row.name)) continue;
    found = true;
    if (remove) body.erase(body.begin() + (long)i);
    else body[i] = SettingText(row);
    break;
  }
  if (!found && !remove) body.push_back(SettingText(row));
  if (found || !remove) project.SetSectionBody("parameterSettings", body);
  return true;
}

std::vector<SamplerPrior> SamplerPriors(const Project &project) {
  std::vector<SamplerPrior> out;
  bool in = false;
  for (const std::string &raw : project.SectionBody("mcmc")) {
    const std::string line = Trim(raw);
    if (line == "<parameters>") in = true;
    else if (line == "</parameters>") in = false;
    else if (in && !line.empty()) {
      const std::vector<std::string> p = Tokens(line);
      if (p.size() < 5) continue;
      SamplerPrior r;
      for (std::size_t i = 0; i + 4 < p.size(); i++) r.name += (i ? " " : "") + p[i];
      r.value = std::atof(p[p.size() - 4].c_str());
      r.mean = std::atof(p[p.size() - 3].c_str());
      r.sigma = std::atof(p[p.size() - 2].c_str());
      r.use = std::atoi(p[p.size() - 1].c_str()) == 1;
      out.push_back(r);
    }
  }
  return out;
}

void SetSamplerPriors(Project &project, const std::vector<SamplerPrior> &rows) {
  std::vector<std::string> body, block;
  block.push_back("<parameters>");
  for (const SamplerPrior &r : rows)
    block.push_back(r.name + " " + Number(r.value) + " " + Number(r.mean) + " " + Number(r.sigma) + " " +
                    (r.use ? "1" : "0"));
  block.push_back("</parameters>");
  bool in = false, placed = false;
  for (const std::string &raw : project.SectionBody("mcmc")) {
    const std::string line = Trim(raw);
    if (line == "<parameters>") {
      in = true;
      if (!placed) body.insert(body.end(), block.begin(), block.end());
      placed = true;
      continue;
    }
    if (line == "</parameters>") {
      in = false;
      continue;
    }
    if (!in) body.push_back(raw);
  }
  if (!placed && !rows.empty()) body.insert(body.end(), block.begin(), block.end());
  if (rows.empty() && !project.HasSection("mcmc")) return;
  project.SetSectionBody("mcmc", body);
}

bool HasCommand(const std::string &command) {
  for (const Command &c : Commands())
    if (c.name == command && c.built) return true;
  return false;
}

}  // namespace front
}  // namespace azure2
