#ifndef AZURE2_IO_LEVELGROUPING_H
#define AZURE2_IO_LEVELGROUPING_H

#include <map>
#include <set>
#include <vector>

namespace azure2 {
namespace io {

/*!
 * Which level of a J^pi group a <levels> line belongs to: the one rule that the engine (CNuc::Fill) and pyazr
 * (AzrModel) both use (docs/upgrade/project.md, P2).  A line joins the level of its group that already holds a line
 * with the same level ID within 1 keV; failing that, a level at the identical energy (reported: two IDs read as one
 * level); otherwise it opens a new level (reported when another level lies within 1 keV: kept as two levels).
 * Lines are fed in file order, active lines only.
 */
class LevelGrouper {
 public:
  enum class Kind {
    Joined,       ///< the level of its ID
    NewLevel,     ///< a new level
    SameEnergy,   ///< merged into a level of another ID at the identical energy (warn)
    NewNearby     ///< a new level within 1 keV of a level of another ID (warn)
  };
  struct Slot {
    int level = 0;   ///< 1-based level number within the group
    Kind kind = Kind::NewLevel;
    int other = 0;   ///< SameEnergy / NewNearby: the other level (1-based)
    int otherID = 0; ///< ... and the first (smallest) level ID it holds
  };
  /// The level of a line of group `group` (any key the caller uses for its J^pi groups) at `energy` (MeV) with
  /// level ID `id`.
  Slot Assign(int group, double energy, int id);
  /// The energy of level `level` (1-based) of `group`, as first set.
  double Energy(int group, int level) const { return groups_.at(group).at(level - 1).energy; }

 private:
  struct Level {
    double energy;
    std::set<int> ids;
  };
  std::map<int, std::vector<Level> > groups_;
};

}  // namespace io
}  // namespace azure2

#endif
