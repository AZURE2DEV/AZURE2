#ifndef AZURE2_INFERENCE_LEVELPARAMETERLAYOUT_H
#define AZURE2_INFERENCE_LEVELPARAMETERLAYOUT_H

#include <string>
#include <vector>

namespace azure2 {
namespace inference {

/*!
 * Where each level's parameters sit in the fit's parameter vector: the level block, which opens the vector (the
 * abundances and the data's parameters follow it).  Every walk of the level block reads it (CompoundSystem::FillMnParams, the
 * fills from a fit or an observed vector, GetTransformParams, ConvertAmplitudeBasis, BuildParamIndexMap, the labels
 * of ParameterLabel.cpp, the gamma_decay records' slots, the model data's derivative columns), so that the order is
 * written once.
 *
 * The order is the format's (param.par, parameters.out, the parameter names and records): J-group after J-group,
 * level after level; per level its energy, then one amplitude per channel of its group, in channel order.  So a
 * level's amplitude in channel ch sits at EnergyIndex + ch (the gamma_decay records' numbers rely on it).
 *
 * Built once from the structure (CompoundSystem::Fill, when the groups, levels and channels are final); a value, shared by the
 * copies of the compound.  Plain data: built from the groups' sizes alone (inference/, docs/upgrade/inference.md).
 */
class LevelParameterLayout {
 public:
  /// One entry of the level block.
  struct Slot {
    enum class Kind { kEnergy, kAmplitude };
    Kind kind = Kind::kEnergy;
    int group = 0;        ///< J-group, 1-based
    int level = 0;        ///< level of the group, 1-based
    int channel = 0;      ///< channel of the group, 1-based; 0 for the energy
    int levelNumber = 0;  ///< the level's number over all groups, 1-based (the N of energy_N and width_N_M)
  };
  /// The slots of one level.
  struct Level {
    int group = 0, level = 0, levelNumber = 0;
    int energy = 0;                ///< index of the energy in the vector, 0-based
    std::vector<int> amplitude;    ///< index of the amplitude of channel ch at [ch - 1], 0-based
  };

  LevelParameterLayout() = default;
  /// The layout of numLevels[j - 1] levels of numChannels[j - 1] channels in J-group j.
  LevelParameterLayout(const std::vector<int> &numLevels, const std::vector<int> &numChannels);

  /// Number of entries of the level block.
  int Size() const { return (int)slots_.size(); }
  /// Every entry, in the vector's order.
  const std::vector<Slot> &Slots() const { return slots_; }
  /// Every level, in the vector's order.
  const std::vector<Level> &Levels() const { return levels_; }
  /// Level la of J-group j, both 1-based.
  const Level &Of(int j, int la) const { return levels_[first_[j - 1] + la - 1]; }
  /// The index of a slot (j, la, ch) in the vector: the energy for ch = 0, else the amplitude of channel ch.
  int IndexOf(int j, int la, int ch) const { return ch == 0 ? Of(j, la).energy : Of(j, la).amplitude[ch - 1]; }
  /// The parameter's name in the vector: energy_N for an energy, width_N_ch for an amplitude (N: levelNumber).
  static std::string Name(const Slot &slot);

 private:
  std::vector<Slot> slots_;
  std::vector<Level> levels_;
  std::vector<int> first_;  ///< index in levels_ of the first level of J-group j at [j - 1]
};

}  // namespace inference
}  // namespace azure2

#endif
