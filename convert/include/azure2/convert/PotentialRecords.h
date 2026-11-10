#ifndef AZURE2_CONVERT_POTENTIALRECORDS_H
#define AZURE2_CONVERT_POTENTIALRECORDS_H

#include <map>
#include <set>
#include <string>
#include <vector>

#include "azure2/experiment/HybridPotentialSettings.h"
#include "azure2/experiment/LevelTables.h"
#include "azure2/io/Record.h"

namespace azure2 {
namespace convert {

/// The `potential` records of a project's hybrid potentials, and what was left out.
struct PotentialRecordsResult {
  io::RecordSet records;
  /// One line per setting that writes no record: a section that switches the model off, a default that applies to no
  /// particle pair, a section for a key that is no pair of <levels>, a photon or beta pair the default would reach.
  std::vector<std::string> notes;
};

/*!
 * The `potential` records (format 3, <potentials>) of hybrid-potential settings (the old <potential> block as
 * ReadPotentialBlock reads it, or a pyazr session's): one record per particle pair whose setting -- its own, else the
 * default -- has the model on, in pair-key order, written by the converter to format 3 and by pyazr's writer.  The
 * default has no record of its own: it is expanded into the pairs it applies to.  A photon or beta pair gets none
 * (the engine computes no channel function of a potential for it).  Every number is written in fm or MeV as the
 * shortest decimal that reads back to the same double, coulomb_radius, matching_radius and step included, so the
 * record holds exactly the potential the settings make.  The pair is named by its label; `pair_key` is added where
 * another pair of `pairs` has the same particles and states.  A record's name is the one `names` gives its pair key
 * (a project's own record), else the label without '+' ("p12C"), made unique against `taken` (and each other) with a
 * suffix _2, _3, ...
 */
PotentialRecordsResult PotentialRecords(const experiment::HybridPotentialSettings &settings,
                                        const std::vector<experiment::PairInfo> &pairs, const std::set<std::string> &taken,
                                        const std::map<int, std::string> &names = {});

}  // namespace convert
}  // namespace azure2

#endif
