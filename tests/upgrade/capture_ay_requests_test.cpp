// The capture analyzing-power terms are built at set-up for the reactions the data sets name (cnuc_split.md section
// 11.5, decision 2; include/CaptureAyRequests.h).  Every project of the tree given is loaded as the engine loads it
// (CNuc::Fill, EData::Fill or MakePoints, RequestCaptureAyTerms, CNuc::Initialize), in a scratch copy, and:
//  (a) every term of an analyzing-power data set with a photon exit has its exit's terms in the pathway table;
//  (b) no beta entrance has capture-A_y terms;
//  (c) the table is the same with the requests in one order and in the reverse order (every entrance and every photon
//      pair requested, both ways, compared row by row);
//  (d) for the projects of (a): loaded without the requests, the capture-A_y points are refused (NaN, the message of
//      GenMatrixFunc::CaptureAyNotBuiltMessage reported once), and with them the points are finite.
// Arguments: the tests/ directory and a scratch directory.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "AZUREParams.h"
#include "CNuc.h"
#include "CaptureAyRequests.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "Segment.h"
#include "ModelPass.h"
#include "PathwayTable.h"
#include "Results.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) failures++;
}

struct Loaded {
  std::unique_ptr<std::ostringstream> sink;
  std::unique_ptr<Config> cfg;
  std::unique_ptr<CNuc> compound;
  std::unique_ptr<EData> data;
  bool ok = false;
};

enum class Requests { FromData, None, AllForward, AllReverse };

// Load as AZUREMain does, in the current directory.
Loaded Load(const std::string &file, Requests requests, bool withData) {
  Loaded l;
  l.sink.reset(new std::ostringstream);
  l.cfg.reset(new Config(*l.sink));
  l.cfg->configfile = file;
  g_config = l.cfg.get();
  l.compound.reset(new CNuc);
  l.data.reset(new EData);
  if (l.cfg->ReadConfigFile() < 0) return l;
  if (withData)
    l.cfg->paramMask |= Config::CALCULATE_WITH_DATA;
  else
    l.cfg->paramMask &= ~Config::CALCULATE_WITH_DATA;
  l.cfg->fileCheckMask = l.cfg->screenCheckMask = 0;
  l.cfg->ResolveMatrixRoute();
  InitializeECAmplitudeCache();
  if (l.compound->Fill(*l.cfg) == -1) return l;
  const int filled = (l.cfg->paramMask & Config::CALCULATE_WITH_DATA) ? l.data->Fill(*l.cfg, l.compound.get())
                                                                       : l.data->MakePoints(*l.cfg, l.compound.get());
  if (filled == -1) return l;
  CNuc &c = *l.compound;
  if (requests == Requests::FromData) RequestCaptureAyTerms(c, *l.data);
  if (requests == Requests::AllForward || requests == Requests::AllReverse) {
    std::vector<std::pair<int, int> > all;
    for (int a = 1; a <= c.NumPairs(); a++)
      for (int p = 1; p <= c.NumPairs(); p++)
        if (c.GetPair(a)->GetPType() != 10 && c.GetPair(p)->GetPType() == 10) all.push_back({a, p});
    if (requests == Requests::AllReverse) std::reverse(all.begin(), all.end());
    for (const auto &r : all) c.RequestCaptureAy(r.first, r.second);
  }
  try {
    c.Initialize(*l.cfg);
  } catch (const std::exception &) {
    return l;
  }
  l.ok = true;
  return l;
}

// The capture-A_y part of the table, one line per exit and per term, doubles as hexadecimal floats.
std::vector<std::string> DumpCaptureAy(CNuc &c) {
  std::vector<std::string> out;
  char buf[256];
  for (int a = 1; a <= c.NumPairs(); a++) {
    const azure2::EntrancePathways &P = c.Pathways().Entrance(a);
    for (const azure2::ExitRow &x : P.exits) {
      std::snprintf(buf, sizeof buf, "%d -> %d: built %d, first %d, %d terms", a, x.pair, (int)x.has_capture_ay,
                    x.first_capture_ay, x.num_capture_ay);
      out.push_back(buf);
    }
    for (const azure2::CaptureAyTerm &t : P.capture_ay) {
      std::snprintf(buf, sizeof buf, "  %d %d/%d/%d %d/%d/%d %a %a", t.kOrder, t.kGroup1, t.path1, (int)t.isEC1, t.kGroup2,
                    t.path2, (int)t.isEC2, t.ak, t.bk);
      out.push_back(buf);
    }
  }
  return out;
}

// The capture-A_y terms of the data sets: (number of such terms, number of them whose exit has its terms built).
std::pair<int, int> CaptureAyTerms(Loaded &l) {
  int terms = 0, built = 0;
  CNuc &c = *l.compound;
  for (int i = 1; i <= l.data->Data().NumSegments(); i++) {
    const ESegment &s = *l.data->Data().GetSegment(i);
    if (!s.IsAnalyzingPower()) continue;
    for (int t = 0; t < s.NumTerms(); t++) {
      if (!c.IsPairKey(s.TermEntranceKey(t)) || !c.IsPairKey(s.TermExitKey(t))) continue;
      const int a = c.GetPairNumFromKey(s.TermEntranceKey(t)), p = c.GetPairNumFromKey(s.TermExitKey(t));
      if (c.GetPair(p)->GetPType() != 10) continue;
      terms++;
      const azure2::EntrancePathways &P = c.Pathways().Entrance(a);
      const int e = P.ExitOf(p);
      if (e && P.exits[e - 1].has_capture_ay) built++;
    }
  }
  return {terms, built};
}

// Evaluate every point once at the start parameters; the capture-A_y points' values.
std::vector<double> CaptureAyValues(Loaded &l) {
  AZUREParams params;
  l.compound->FillMnParams(params, l.cfg.get());
  l.data->FillMnParams(params);
  if (l.data->Initialize(l.compound.get(), *l.cfg) == -1) return {};
  azure2::Results results;
  l.data->SizeResults(results);
  const vector_r p = params.GetMinuitParams().Params();
  l.compound->FillCompoundFromParams(p);
  l.data->FillNormsFromParams(p);
  l.data->FillEnergyShiftsFromParams(p, results, l.compound.get(), l.cfg.get());
  if (l.cfg->paramMask & Config::USE_BRUNE_FORMALISM) l.compound->CalcShiftFunctions(*l.cfg, false);
  azure2::EvaluateModel(*l.data, l.compound.get(), *l.cfg, results, azure2::LeafValues::Calculate);
  std::vector<double> values;
  CNuc &c = *l.compound;
  for (int i = 1; i <= l.data->Data().NumSegments(); i++) {
    const ESegment &s = *l.data->Data().GetSegment(i);
    if (!s.IsAnalyzingPower() || !c.IsPairKey(s.GetExitKey())) continue;
    if (c.GetPair(c.GetPairNumFromKey(s.GetExitKey()))->GetPType() != 10) continue;
    for (int k = 1; k <= s.NumPoints(); k++) values.push_back(s.GetPoint(*l.data, k)->GetFitCrossSection(results));
  }
  return values;
}

}  // namespace

int main(int argc, char **argv) {
  namespace fs = std::filesystem;
  if (argc < 3) {
    std::fprintf(stderr, "usage: capture_ay_requests_test <tests dir> <scratch dir>\n");
    return 2;
  }
  const fs::path tests = fs::absolute(argv[1]), scratch = fs::absolute(argv[2]);
  std::vector<fs::path> files;
  for (const auto &e : fs::recursive_directory_iterator(tests))
    if (e.is_regular_file() && e.path().extension() == ".azr" && e.path().string().find("/output/") == std::string::npos)
      files.push_back(e.path());
  std::sort(files.begin(), files.end());
  fs::remove_all(scratch);
  int projects = 0, captureProjects = 0, captureTerms = 0, betaEntrances = 0;
  std::vector<std::pair<fs::path, bool> > withCaptureAy;  // (file, data mode)
  for (const fs::path &file : files) {
    // A scratch copy of the project's folder (the data are read relative to it; nothing is written into tests/).
    const fs::path dir = scratch / fs::relative(file.parent_path(), tests);
    fs::create_directories(dir);
    fs::copy(file.parent_path(), dir, fs::copy_options::recursive | fs::copy_options::overwrite_existing);
    fs::create_directories(dir / "output");
    fs::create_directories(dir / "checks");
    fs::current_path(dir);
    for (int withData = 1; withData >= 0; withData--) {
      const std::string name = file.filename().string(),
                        tag = fs::relative(file, tests).string() + (withData ? " (data)" : " (extrapolation)");
      Loaded l = Load(name, Requests::FromData, withData);
      if (!l.ok) continue;
      projects++;
      // (a)
      const std::pair<int, int> n = CaptureAyTerms(l);
      if (n.first > 0) {
        captureProjects++;
        captureTerms += n.first;
        withCaptureAy.push_back({file, withData == 1});
        Check(n.second == n.first, tag + ": " + std::to_string(n.second) + " of " + std::to_string(n.first) +
                                       " capture-A_y terms of the data sets have their terms built");
      }
      // (b) and (c)
      Loaded forward = Load(name, Requests::AllForward, withData), reverse = Load(name, Requests::AllReverse, withData);
      if (forward.ok && reverse.ok) {
        const std::vector<std::string> f = DumpCaptureAy(*forward.compound), r = DumpCaptureAy(*reverse.compound);
        if (f != r) Check(false, tag + ": the table depends on the order of the requests");
        CNuc &c = *forward.compound;
        for (int a = 1; a <= c.NumPairs(); a++) {
          if (c.GetPair(a)->GetPType() != 20) continue;
          betaEntrances++;
          if (!c.Pathways().Entrance(a).capture_ay.empty()) Check(false, tag + ": a beta entrance has capture-A_y terms");
        }
      }
    }
  }
  std::printf("%d project loads (data and extrapolation mode); %d with capture-A_y data sets (%d terms):\n", projects,
              captureProjects, captureTerms);
  for (const auto &f : withCaptureAy)
    std::printf("    %s (%s)\n", fs::relative(f.first, tests).c_str(), f.second ? "data" : "extrapolation");
  Check(captureProjects > 0, "some project has capture-A_y data");
  Check(betaEntrances > 0, std::to_string(betaEntrances) + " beta entrances seen, none with capture-A_y terms");
  // (d)
  bool refusalMessage = false;
  for (const auto &item : withCaptureAy) {
    const fs::path &file = item.first;
    const bool withData = item.second;
    const fs::path dir = scratch / fs::relative(file.parent_path(), tests);
    fs::current_path(dir);
    const std::string tag = fs::relative(file, tests).string() + (withData ? " (data)" : " (extrapolation)");
    Loaded without = Load(file.filename().string(), Requests::None, withData);
    const std::vector<double> refused = CaptureAyValues(without);
    bool allNaN = !refused.empty();
    for (double v : refused) allNaN = allNaN && std::isnan(v);
    // EPoint::MarkFailed reports the first failure of the process only: the message is looked for once, below.
    refusalMessage = refusalMessage || without.sink->str().find("has no terms") != std::string::npos;
    Check(allNaN, tag + ": without the requests the " + std::to_string(refused.size()) +
                      " capture-A_y points are refused (NaN)");
    Loaded with = Load(file.filename().string(), Requests::FromData, withData);
    const std::vector<double> values = CaptureAyValues(with);
    bool finite = !values.empty();
    for (double v : values) finite = finite && std::isfinite(v);
    Check(finite && with.sink->str().find("has no terms") == std::string::npos,
          tag + ": with them the " + std::to_string(values.size()) + " capture-A_y points are finite");
  }
  Check(refusalMessage, "the refusal names the missing terms (\"has no terms\")");
  std::printf("%s: %d failed\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
