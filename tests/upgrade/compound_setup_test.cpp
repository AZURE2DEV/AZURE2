// The one set-up of the compound system and the data (include/CompoundSystemSetUp.h,
// docs/upgrade/cnuc_split.md section 12.3, the one set-up as built):
//   - in the middle (the data's set-up), what the last part writes refuses to be read, each accessor with a
//     std::logic_error naming itself: the pathway table, the formal external-capture blocks, the capture masks, B_c,
//     a final state's energy; after the set-up all of them read;
//   - the requests keep their first order and drop repeats; the entrances asked for are the pathway table's;
//   - the refusals come back as the step that refused: the <levels> (a final state missing), the data, the last part
//     (no standard R-matrix for 7Li_p_ay in the formal basis, with its message), and no compound with them;
//   - the structure alone (InspectCompoundStructure) is never finished and carries the file's input values.
//
//   compound_setup_test <tests dir>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>

#include "SetUpForTests.h"
#include "ALevel.h"
#include "AChannel.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "JGroup.h"
#include "PathwayTable.h"

namespace {
int failures = 0;
void check(const std::string &name, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s%s\n", ok ? "ok  " : "FAIL", name.c_str(), ok || detail.empty() ? "" : "  -- ",
              ok ? "" : detail.c_str());
  if (!ok) failures++;
}

// The message of the std::logic_error `read` throws, or "" when it does not throw one.
std::string LogicErrorOf(const std::function<void()> &read) {
  try {
    read();
  } catch (const std::logic_error &e) {
    return e.what();
  } catch (...) {
    return "(another exception)";
  }
  return "";
}

// The first external-capture final state of `c`: (group, level), or (0, 0).
std::pair<int, int> FirstFinalState(const CompoundSystem &c) {
  for (int j = 1; j <= c.NumJGroups(); j++)
    for (int la = 1; la <= c.GetJGroup(j)->NumLevels(); la++)
      if (c.GetJGroup(j)->GetLevel(la)->IsECLevel()) return {j, la};
  return {0, 0};
}

// Every accessor of what the last part writes, read from `c`: the error message of each ("" when it read).
std::vector<std::pair<std::string, std::string> > ReadLastPart(const CompoundSystem &c) {
  const std::pair<int, int> fs = FirstFinalState(c);
  return {
      {"CompoundSystem::Pathways", LogicErrorOf([&] { c.Pathways(); })},
      {"CompoundSystem::HasFormalECFinalBlocks", LogicErrorOf([&] { c.HasFormalECFinalBlocks(); })},
      {"CompoundSystem::FormalCaptureMaskSlot", LogicErrorOf([&] { c.FormalCaptureMaskSlot(1, 1); })},
      {"AChannel::GetBoundaryCondition", LogicErrorOf([&] { c.GetJGroup(1)->GetChannel(1)->GetBoundaryCondition(); })},
      {"ALevel::GetECEnergy", LogicErrorOf([&] { c.GetJGroup(fs.first)->GetLevel(fs.second)->GetECEnergy(); })},
  };
}

struct Project {
  std::ostringstream sink;
  Config cfg{sink};
  bool read = false;
  // `azr` relative to `dir`, the directory the project's data paths start from.
  Project(const std::filesystem::path &dir, const std::string &azr, bool formal) {
    std::filesystem::current_path(dir);
    cfg.configfile = azr;
    cfg.paramMask |= Config::USE_EXTERNAL_CAPTURE;
    if (formal) cfg.paramMask &= ~(Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
    read = cfg.ReadConfigFile() >= 0;
    cfg.fileCheckMask = cfg.screenCheckMask = 0;
    cfg.ResolveMatrixRoute();
    InitializeECAmplitudeCache(cfg);
  }
};
}  // namespace

int main(int argc, char **argv) {
  namespace fs = std::filesystem;
  if (argc != 2) {
    std::fprintf(stderr, "usage: compound_setup_test <tests dir>\n");
    return 2;
  }
  const fs::path tests = fs::absolute(argv[1]);

  std::printf("the requests\n");
  {
    CompoundSetUpRequests r;
    for (int p : {3, 1, 3, 2, 1}) r.AddEntrance(p);
    for (std::pair<int, int> q : {std::pair<int, int>{1, 4}, {2, 4}, {1, 4}}) r.AddCaptureAy(q.first, q.second);
    check("entrances in first order, each once", r.entrancePairs == std::vector<int>{3, 1, 2});
    check("capture-A_y pairs in first order, each once",
          r.captureAyPairs == std::vector<std::pair<int, int> >{{1, 4}, {2, 4}});
  }

  std::printf("in the middle and after (16O(p,g) with external capture, validate_run_refusals/models/capture.f3, formal basis)\n");
  {
    Project p(tests / "validate_run_refusals/models", "capture.f3", true);
    check("reads the project", p.read, p.sink.str().substr(0, 300));
    WorkingCopy data;
    std::vector<std::pair<std::string, std::string> > middle;
    int entrance = 0;
    CompoundSetUp setUp = SetUpCompoundSystem(p.cfg, [&](const CompoundSystem &structure, const LevelInputs &inputs)
                                                         -> std::optional<CompoundSetUpRequests> {
      check("the structure has a final state", FirstFinalState(structure).first > 0);
      middle = ReadLastPart(structure);
      if (data.Fill(p.cfg, &structure, inputs) != 0) return std::nullopt;
      CompoundSetUpRequests requests;
      RequestDataEntrances(structure, p.cfg, requests);
      entrance = requests.entrancePairs.empty() ? 0 : requests.entrancePairs.front();
      return requests;
    });
    for (const auto &m : middle)
      check("in the middle, " + m.first + " refuses, naming itself", m.second.find(m.first) != std::string::npos,
            m.second.empty() ? "it read" : m.second);
    check("set up", setUp.refusal == CompoundSetUp::Refusal::kNone && setUp.compound != nullptr,
          p.sink.str().substr(0, 300));
    if (setUp.compound) {
      for (const auto &m : ReadLastPart(*setUp.compound))
        check("after the set-up, " + m.first + " reads", m.second.empty(), m.second);
      check("the requested entrance has its pathways", entrance > 0 && setUp.compound->GetPair(entrance)->IsEntrance());
      check("the inputs are transformed", setUp.inputs.TransformedIn());
    }
  }

  std::printf("the refusals\n");
  {
    Project p(tests / "ec_final_state_listed", "models/missing_level.azr", false);
    bool middleRan = false;
    const CompoundSetUp setUp = SetUpCompoundSystem(p.cfg, [&](const CompoundSystem &, const LevelInputs &) {
      middleRan = true;
      return std::optional<CompoundSetUpRequests>(CompoundSetUpRequests());
    });
    check("a final state that is not a level: the <levels>, before the data",
          setUp.refusal == CompoundSetUp::Refusal::kUnreadable && !middleRan && !setUp.compound);
  }
  {
    Project p(tests / "upgrade/T1", "T1.azr", false);
    const CompoundSetUp setUp =
        SetUpCompoundSystem(p.cfg, [](const CompoundSystem &, const LevelInputs &) { return std::optional<CompoundSetUpRequests>(); });
    check("the data refused: the data, no compound", setUp.refusal == CompoundSetUp::Refusal::kData && !setUp.compound);
  }
  {
    Project p(tests / "7Li_p_ay", "7Li_p_ay.azr", true);
    const CompoundSetUp setUp = test_setup::WithoutData(p.cfg);
    check("no standard R-matrix (7Li_p_ay, formal): the last part, with its message, no compound",
          setUp.refusal == CompoundSetUp::Refusal::kFinish && !setUp.compound &&
              setUp.finishError.find("Exception thrown") != std::string::npos,
          setUp.finishError.substr(0, 200));
  }

  std::printf("the structure alone\n");
  {
    Project p(tests / "upgrade/T1", "T1.azr", false);
    std::string pathways;
    bool transformed = true;
    const bool read = InspectCompoundStructure(p.cfg, [&](const CompoundSystem &structure, const LevelState &,
                                                          const LevelInputs &inputs) {
      pathways = LogicErrorOf([&] { structure.Pathways(); });
      transformed = inputs.TransformedIn();
    });
    check("read, never finished, the file's input values",
          read && pathways.find("CompoundSystem::Pathways") != std::string::npos && !transformed, pathways);
  }

  std::printf("%s (%d failures)\n", failures ? "FAILED" : "passed", failures);
  return failures ? 1 : 0;
}
