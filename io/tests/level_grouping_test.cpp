// io::LevelGrouper (LevelGrouping.h; docs/upgrade/project.md P2), the level rule of CNuc::Fill and pyazr's Project:
// a line joins the level of its ID within 1 keV; otherwise one at the identical energy (SameEnergy); otherwise a new
// level, reported (NewNearby) when another level of the group is within 1 keV; groups are independent.
#include <cstdio>
#include <string>

#include "azure2/io/LevelGrouping.h"

using azure2::io::LevelGrouper;

namespace {
int failures = 0, checks = 0;
void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}
}  // namespace

int main() {
  LevelGrouper g;
  check("first line: new level 1", g.Assign(1, 5.0, 10).level == 1);
  check("same ID, same energy: joined", g.Assign(1, 5.0, 10).kind == LevelGrouper::Kind::Joined);
  check("same ID within 1 keV: joined", g.Assign(1, 5.0005, 10).level == 1);
  LevelGrouper::Slot s = g.Assign(1, 5.0, 11);
  check("other ID at the identical energy: merged, reported", s.level == 1 && s.kind == LevelGrouper::Kind::SameEnergy &&
                                                                  s.otherID == 10);
  check("the merged ID then joins", g.Assign(1, 5.0, 11).kind == LevelGrouper::Kind::Joined);
  s = g.Assign(1, 5.0007, 12);
  check("other ID within 1 keV: a new level, reported", s.level == 2 && s.kind == LevelGrouper::Kind::NewNearby &&
                                                            s.other == 1);
  check("same ID beyond 1 keV: a new level", g.Assign(1, 6.0, 10).level == 3);
  check("another group is independent", g.Assign(2, 5.0, 10).level == 1);
  check("a later line of an earlier level finds it", g.Assign(1, 5.0, 10).level == 1);
  check("the level energy is the first line's", g.Energy(1, 2) == 5.0007);
  std::printf("%s: %d checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
