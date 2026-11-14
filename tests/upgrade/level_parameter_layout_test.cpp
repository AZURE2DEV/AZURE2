// The level block of the parameter vector (include/LevelParameterLayout.h).
//  (a) A layout made by hand (groups of 2 levels x 3 channels, none, 1 level x 1 channel): every index, name and
//      level number against the order of the format (per level its energy, then one amplitude per channel).
//  (b) Every project given, loaded with CompoundSystem::Fill: the format's order written out in the test (groups, levels, per
//      level its energy, then its channels in order) against every slot of the layout, its index and name, and the
//      names CompoundSystem::FillMnParams gives the level block (a gamma_decay record's slot excepted: the record names its own
//      numbers).
// Arguments: .azr files.  Reads them in place and writes nothing.
#include <cstdio>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <Minuit2/MnUserParameters.h>

#include "SetUpForTests.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "JGroup.h"
#include "LevelParameterLayout.h"

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

using Slot = LevelParameterLayout::Slot;

void ByHand() {
  const LevelParameterLayout layout({2, 0, 1}, {3, 4, 1});
  Check(layout.Size() == 2 * 4 + 0 + 1 * 2, "size by hand");
  const std::vector<std::string> names = {"energy_1", "width_1_1", "width_1_2", "width_1_3", "energy_2", "width_2_1",
                                          "width_2_2", "width_2_3", "energy_3", "width_3_1"};
  for (int i = 0; i < layout.Size() && i < (int)names.size(); i++)
    Check(LevelParameterLayout::Name(layout.Slots()[i]) == names[i], "name of slot " + std::to_string(i));
  Check(layout.Of(1, 1).energy == 0 && layout.Of(1, 2).energy == 4 && layout.Of(3, 1).energy == 8, "energy indices");
  Check(layout.Of(1, 2).amplitude == std::vector<int>({5, 6, 7}), "amplitudes of level 2 of group 1");
  Check(layout.IndexOf(3, 1, 1) == 9 && layout.IndexOf(3, 1, 0) == 8, "IndexOf");
  Check(layout.Of(3, 1).levelNumber == 3 && layout.Slots()[9].group == 3 && layout.Slots()[9].channel == 1,
        "level number and slot of the third level");
  Check(layout.Levels().size() == 3, "three levels");
}

void OnProject(const std::string &file) {
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = file;
  if (cfg.ReadConfigFile() < 0) {
    Check(false, file + ": cannot read");
    return;
  }
  cfg.fileCheckMask = cfg.screenCheckMask = 0;
  cfg.ResolveMatrixRoute();
  InitializeECAmplitudeCache(cfg);
  // The structure finished, with the levels as its set-up reads them from the file (the values FillMnParams seeds).
  const CompoundSetUp setUp = test_setup::WithoutData(cfg);
  LevelState state;
  LevelInputs inputs;
  if (setUp.refusal != CompoundSetUp::Refusal::kNone || !test_setup::LevelsBeforeTheLastPart(cfg, state, inputs)) {
    Check(false, file + ": set-up failed");
    return;
  }
  const CompoundSystem &compound = *setUp.compound;
  const LevelParameterLayout &layout = compound.ParameterLayout();
  AZUREParams params;
  compound.FillMnParams(params, state, inputs, &cfg);
  const ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  // The format's order written out here, independently of the layout: the slots and the names FillMnParams gives.
  int index = 0, levelNumber = 1;
  auto slotIs = [&](int j, int la, int ch, const std::string &name) {
    const std::string where = file + ": slot " + std::to_string(index);
    if (index >= layout.Size() || index >= (int)p.Params().size()) {
      Check(false, where + " missing");
      return;
    }
    const Slot &slot = layout.Slots()[index];
    Check(slot.group == j && slot.level == la && slot.channel == ch && slot.levelNumber == levelNumber,
          where + ": group, level, channel");
    Check(layout.IndexOf(j, la, ch) == index, where + ": IndexOf");
    Check(LevelParameterLayout::Name(slot) == name, where + ": name " + LevelParameterLayout::Name(slot));
    if (!compound.GetGammaDecays().OwnsSlot(j, la, ch)) Check(p.Name(index) == name, where + ": FillMnParams' name");
  };
  for (int j = 1; j <= compound.NumJGroups(); j++) {
    const JGroup *group = compound.GetJGroup(j);
    for (int la = 1; la <= group->NumLevels(); la++, levelNumber++) {
      slotIs(j, la, 0, "energy_" + std::to_string(levelNumber));
      index++;
      for (int ch = 1; ch <= group->NumChannels(); ch++, index++)
        slotIs(j, la, ch, "width_" + std::to_string(levelNumber) + "_" + std::to_string(ch));
    }
  }
  Check(layout.Size() == index, file + ": size");
  std::printf("  %s: %d slots\n", file.c_str(), layout.Size());
}

}  // namespace

int main(int argc, char **argv) {
  ByHand();
  for (int a = 1; a < argc; a++) OnProject(argv[a]);
  std::printf("%s (%d failures)\n", failures ? "FAILED" : "passed", failures);
  return failures ? 1 : 0;
}
