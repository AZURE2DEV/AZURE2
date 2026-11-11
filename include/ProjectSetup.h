#ifndef PROJECTSETUP_H
#define PROJECTSETUP_H

#include <set>
#include <utility>
#include <vector>

#include "azure2/io/LevelRows.h"

class Config;

/*!
 * What the set-up reads of a project before the compound nucleus is filled (src/ProjectSetup.cpp): the <levels>
 * rows, read once (a format-3 project: by Config::ReadConfigFile, with the rest of the project), and the exit pairs of
 * its active data sets and extrapolations, from their specs (and the reaction-rate pair): the pairs whose external
 * capture a run can use (CNuc::ParseExternalCapture refuses a missing final state of one of them).
 */
/// false when the project has no <levels> section
bool ProjectLevelSection(const Config &, azure2::io::LevelSection &out);
std::set<int> ProjectExitPairKeys(const Config &);
/// The (entrance, exit) pair keys of every term (the leading one and the components of a SUM or RATIO) of the active
/// data sets (`data`) or extrapolations, in order; false when the project
/// has none of those sections (an old project) or its specs could not be read.
bool ProjectReactionPairs(const Config &, bool data, std::vector<std::pair<int, int>> &pairs);
/// An active data set with a THM observable.
bool ProjectHasThmData(const Config &);

#endif
