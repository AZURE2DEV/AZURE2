// A session rebuilt at another channel radius (AZUREAPI::SetRadius, the rebuild) sets its model up as a session opened
// at that radius does (AZUREAPI::SetUpModel, docs/upgrade/cnuc_split.md section 12.3):
//   - a session for the reaction rate (tests/reaction_rate/models/p12C.azr, entrance 1, exit 2), rebuilt with pair 1
//     at 4.0 fm: accepted, and its rate on a temperature grid equals the rate of a fresh session of the project with
//     4.0 fm written for pair 1 (and differs from the rate at the project's 3.6 fm);
//   - a project of model data alone (tests/upgrade/T2/T2_true.azr without its data sets, with a <modelData> block),
//     rebuilt with pair 1 at 5.0 fm: accepted, and its model values equal a fresh session's at 5.0 fm (and differ from
//     those at the project's 4.85 fm).
//
//   session_rebuild_test <tests dir>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "AZUREAPI.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "ReactionRate.h"

namespace fs = std::filesystem;

namespace {
int failures = 0;
void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), ok || detail.empty() ? "" : "  -- ",
              ok ? "" : detail.c_str());
  if (!ok) failures++;
}

std::string ReadText(const fs::path &file) {
  std::ifstream in(file);
  std::stringstream text;
  text << in.rdbuf();
  return text.str();
}

// `text` with the channel radius (column 28 of a <levels> row) of every row of pair key `pairKey` (column 6) set to
// `radius` (written as given).
std::string WithRadius(const std::string &text, int pairKey, const std::string &radius) {
  std::istringstream in(text);
  std::ostringstream out;
  std::string line;
  bool inLevels = false;
  while (std::getline(in, line)) {
    if (line == "<levels>") inLevels = true;
    if (line == "</levels>") inLevels = false;
    std::istringstream fields(line);
    std::vector<std::string> f;
    for (std::string x; fields >> x;) f.push_back(x);
    if (inLevels && f.size() >= 28 && std::atoi(f[5].c_str()) == pairKey) {
      f[27] = radius;
      line.clear();
      for (size_t i = 0; i < f.size(); i++) line += (i ? " " : "") + f[i];
    }
    out << line << '\n';
  }
  return out.str();
}

// A project directory under `root`: `name`/t.azr with `text`, and its output/ and checks/ directories.
fs::path MakeProject(const fs::path &root, const std::string &name, const std::string &text) {
  const fs::path dir = root / name;
  fs::create_directories(dir / "output");
  fs::create_directories(dir / "checks");
  std::ofstream(dir / "t.azr") << text;
  return dir;
}

// A session as the binding opens one (api/src/pybind_api.cpp, Session): a rate session when rateEntrance != 0.
struct Session {
  std::ostringstream sink;
  Config cfg{sink};
  std::unique_ptr<AZUREAPI> api;
  bool opened = false;
  Session(const fs::path &dir, int rateEntrance, int rateExit) {
    fs::current_path(dir);
    cfg.configfile = (dir / "t.azr").string();
    cfg.paramMask |= Config::USE_API | Config::CALCULATE_WITH_DATA;
    if (rateEntrance != 0) {
      cfg.paramMask &= ~Config::CALCULATE_WITH_DATA;
      cfg.paramMask |= Config::CALCULATE_REACTION_RATE;
      cfg.rateParams.entrancePair = rateEntrance;
      cfg.rateParams.exitPair = rateExit;
    }
    if (cfg.ReadConfigFile() != 0) return;
    cfg.ResolveMatrixRoute();
    InitializeECAmplitudeCache(cfg);
    api.reset(new AZUREAPI(cfg));
    opened = api->Initialize() == 0;
  }
  // N_A<sigma v> on T9 = 0.1, 0.4, 0.7, 1.0 at the session's parameters, as the binding's reaction_rate computes it.
  std::vector<double> Rate() {
    cfg.rateParams.useFile = false;
    cfg.rateParams.minTemp = 0.1;
    cfg.rateParams.maxTemp = 1.0;
    cfg.rateParams.tempStep = 0.3;
    ReactionRate rate(api->compound(), api->state(), api->inputs(), api->params_all_rwa(), cfg,
                      cfg.rateParams.entrancePair, cfg.rateParams.exitPair);
    rate.CalculateRates();
    std::vector<double> out;
    for (const RateData &r : rate.Rates()) out.push_back(r.rate);
    return out;
  }
  // The model value of every <modelData> datum (and its imaginary part) at the session's parameters.
  std::vector<double> ModelValues() {
    std::vector<double> out;
    for (const AZUREAPI::ModelDatumValue &v : api->ModelDataRWA(api->params_values_rwa())) {
      out.push_back(v.model);
      out.push_back(v.modelImag);
    }
    return out;
  }
};

std::string List(const std::vector<double> &v) {
  std::string s;
  char buffer[32];
  for (double x : v) {
    std::snprintf(buffer, sizeof buffer, "%s%.12g", s.empty() ? "" : " ", x);
    s += buffer;
  }
  return s;
}

// The largest relative distance between `a` and `b` (infinite when their lengths differ or either is empty).
double Distance(const std::vector<double> &a, const std::vector<double> &b) {
  if (a.size() != b.size() || a.empty()) return INFINITY;
  double d = 0.0;
  for (size_t i = 0; i < a.size(); i++) {
    const double scale = std::max(std::fabs(a[i]), std::fabs(b[i]));
    if (scale > 0.0) d = std::max(d, std::fabs(a[i] - b[i]) / scale);
  }
  return d;
}
}  // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: session_rebuild_test <tests dir>\n");
    return 2;
  }
  const fs::path tests = fs::absolute(argv[1]);
  // A new directory under the system's temporary directory (std::filesystem: mkdtemp is POSIX only).
  fs::path root;
  {
    std::error_code ec;
    const fs::path base = fs::temp_directory_path(ec);
    std::random_device device;
    for (int attempt = 0; !ec && attempt < 16 && root.empty(); ++attempt) {
      const fs::path candidate = base / ("session_rebuild." + std::to_string(device()));
      std::error_code made;
      if (fs::create_directory(candidate, made) && !made) root = candidate;
    }
  }
  if (root.empty()) {
    std::fprintf(stderr, "session_rebuild_test: no temporary directory\n");
    return 2;
  }

  std::printf("a rate session rebuilt at another radius (p12C, entrance 1, exit 2; pair 1 at 4.0 fm, the file's 3.6)\n");
  {
    const std::string text = ReadText(tests / "reaction_rate/models/p12C.azr");
    const fs::path asGiven = MakeProject(root, "rate", text);
    const fs::path atNew = MakeProject(root, "rate_4fm", WithRadius(text, 1, "4.0"));
    std::vector<double> before, rebuilt, fresh;
    {
      Session s(asGiven, 1, 2);
      check("opened", s.opened, s.sink.str().substr(0, 300));
      if (s.opened) {
        before = s.Rate();
        const bool accepted = s.api->SetRadius(1, 4.0);
        check("rebuilt at 4.0 fm: accepted", accepted, s.sink.str().substr(0, 300));
        if (accepted) rebuilt = s.Rate();
      }
    }
    {
      Session s(atNew, 1, 2);
      check("a fresh session at 4.0 fm opened", s.opened, s.sink.str().substr(0, 300));
      if (s.opened) fresh = s.Rate();
    }
    std::printf("        N_A<sigma v> at T9 = 0.1 0.4 0.7 1.0: 3.6 fm %s; rebuilt at 4.0 fm %s; fresh at 4.0 fm %s\n",
                List(before).c_str(), List(rebuilt).c_str(), List(fresh).c_str());
    check("the rate after the rebuild equals a fresh session's at 4.0 fm", Distance(rebuilt, fresh) == 0.0);
    check("and the radius is in it: the rate at 3.6 fm differs", Distance(before, fresh) > 1e-6);
  }

  std::printf("a project of model data alone rebuilt at another radius (T2 without data sets; pair 1 at 5.0 fm, the "
              "file's 4.85)\n");
  {
    // T2_true.azr without its data sets and what refers to them (the THM experiment, the targets, the prior
    // covariance), with model data: a resonance strength, a total width, a partial width and a scattering length.
    std::istringstream in(ReadText(tests / "upgrade/T2/T2_true.azr"));
    std::ostringstream out;
    bool skip = false;
    for (std::string line; std::getline(in, line);) {
      if (line.find("priorCovariance") != std::string::npos) continue;
      if (line == "</data>" || line == "</thm>" || line == "</targets>") skip = false;
      if (!skip) out << line << '\n';
      if (line == "<data>" || line == "<thm>" || line == "<targets>") skip = true;
    }
    out << "<modelData>\n"
           "datum wg_L5\n    quantity resonance_strength\n    level 5\n    J 3-\n    reaction 15N(p,a)12C\n"
           "    value 2.7318028 keV\n    error 0.1 keV\n"
           "datum G_L5\n    quantity total_width\n    level 5\n    J 3-\n    value 15.622 keV\n    error 1 keV\n"
           "datum Gp_L5\n    quantity partial_width\n    level 5\n    J 3-\n    pair p+15N\n    value 1802 eV\n"
           "    error 100 eV\n"
           "datum a_p15N\n    quantity scattering_length\n    pair p+15N\n    s 1\n    value 0 fm\n    error 100 fm\n"
           "</modelData>\n";
    const fs::path asGiven = MakeProject(root, "model", out.str());
    const fs::path atNew = MakeProject(root, "model_5fm", WithRadius(out.str(), 1, "5.0"));
    std::vector<double> before, rebuilt, fresh;
    {
      Session s(asGiven, 0, 0);
      check("opened", s.opened, s.sink.str().substr(0, 300));
      if (s.opened) {
        before = s.ModelValues();
        const bool accepted = s.api->SetRadius(1, 5.0);
        check("rebuilt at 5.0 fm: accepted", accepted, s.sink.str().substr(0, 300));
        if (accepted) rebuilt = s.ModelValues();
      }
    }
    {
      Session s(atNew, 0, 0);
      check("a fresh session at 5.0 fm opened", s.opened, s.sink.str().substr(0, 300));
      if (s.opened) fresh = s.ModelValues();
    }
    std::printf("        model values (each with its imaginary part): 4.85 fm %s; rebuilt at 5.0 fm %s; fresh at 5.0 fm "
                "%s\n", List(before).c_str(), List(rebuilt).c_str(), List(fresh).c_str());
    check("the model values after the rebuild equal a fresh session's at 5.0 fm", Distance(rebuilt, fresh) == 0.0);
    check("and the radius is in them: the values at 4.85 fm differ", Distance(before, fresh) > 1e-6);
  }

  fs::current_path(tests);
  fs::remove_all(root);
  std::printf("%s (%d failures)\n", failures ? "FAILED" : "passed", failures);
  return failures ? 1 : 0;
}
