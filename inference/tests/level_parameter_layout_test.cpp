// inference_level_parameter_layout: the level block of the parameter vector (azure2/inference/LevelParameterLayout.h)
// for three J groups of unequal size, every index, name and level number written out by hand.  The order is the
// parameter file's: group after group, level after level, the energy and then one amplitude per channel; an
// amplitude placed before its energy, or levels numbered per group, fails here.
#include "azure2/inference/LevelParameterLayout.h"

#include <string>
#include <vector>

#include "check.h"

using azure2::inference::LevelParameterLayout;

int main() {
  // Group 1: 2 levels of 3 channels; group 2: no level (2 channels); group 3: 1 level of 1 channel.
  const LevelParameterLayout layout({2, 0, 1}, {3, 2, 1});
  check("size: 2 x (1 + 3) + 0 + 1 x (1 + 1) = 10", layout.Size() == 10 && layout.Slots().size() == 10);
  check("three levels", layout.Levels().size() == 3);

  const std::vector<std::string> names = {"energy_1", "width_1_1", "width_1_2", "width_1_3", "energy_2",
                                          "width_2_1", "width_2_2", "width_2_3", "energy_3", "width_3_1"};
  const int group[] = {1, 1, 1, 1, 1, 1, 1, 1, 3, 3}, level[] = {1, 1, 1, 1, 2, 2, 2, 2, 1, 1};
  const int channel[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1}, number[] = {1, 1, 1, 1, 2, 2, 2, 2, 3, 3};
  bool slots = true, named = true;
  for (int i = 0; i < layout.Size(); i++) {
    const LevelParameterLayout::Slot &s = layout.Slots()[i];
    slots = slots && s.group == group[i] && s.level == level[i] && s.channel == channel[i] && s.levelNumber == number[i] &&
            (s.kind == LevelParameterLayout::Slot::Kind::kEnergy) == (channel[i] == 0);
    named = named && LevelParameterLayout::Name(s) == names[i];
  }
  check("every slot: group, level, channel, level number, kind", slots);
  check("every name: energy_N, width_N_ch with N over all groups", named);

  check("IndexOf: energies", layout.IndexOf(1, 1, 0) == 0 && layout.IndexOf(1, 2, 0) == 4 && layout.IndexOf(3, 1, 0) == 8);
  check("IndexOf: amplitudes", layout.IndexOf(1, 1, 1) == 1 && layout.IndexOf(1, 1, 3) == 3 && layout.IndexOf(1, 2, 2) == 6 &&
                                   layout.IndexOf(3, 1, 1) == 9);
  // The gamma_decay records rely on it: a level's amplitude in channel ch sits at its energy's index + ch.
  bool contiguous = true;
  for (const LevelParameterLayout::Level &l : layout.Levels())
    for (std::size_t ch = 1; ch <= l.amplitude.size(); ch++) contiguous = contiguous && l.amplitude[ch - 1] == l.energy + (int)ch;
  check("amplitude of channel ch at energy + ch", contiguous);
  const LevelParameterLayout::Level &last = layout.Of(3, 1);
  check("Of(3, 1) past an empty group", last.group == 3 && last.level == 1 && last.levelNumber == 3 && last.energy == 8 &&
                                            last.amplitude.size() == 1);
  check("an empty layout", LevelParameterLayout().Size() == 0 && LevelParameterLayout({}, {}).Levels().empty());
  return result();
}
