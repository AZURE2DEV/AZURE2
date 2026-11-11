#ifndef AZURE2_EXPERIMENT_LEVELEDITS_H
#define AZURE2_EXPERIMENT_LEVELEDITS_H

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/LevelTables.h"

namespace azure2 {
namespace io {
class Document;
}
namespace experiment {

/*!
 * The <levels> section of a project as a person reads it -- levels named by J^pi and energy, channels by pair and
 * (l, s) -- and its edits as line replacements of the document (docs/upgrade/pyazr_on_records.md, "Edits of the
 * nuclear part").  The rows are read by io::LevelRow and grouped by LevelTablesOfRows (the engine's CNuc::Fill
 * numbering); a row is written by io::EditLevelRow, which keeps every byte of the row but the tokens it replaces.
 * Rows an edit does not touch keep their bytes.  Numbers are written as the shortest decimal that reads back as the
 * same double (ShortestNumberText).  Library azure2_experiment_reader (it reads the document).
 */

/// One row of <levels> (one channel of a level) and where it stands.
struct LevelRowLine {
  int line = 0;  ///< 0-based index of the row's line in the document
  LevelRowSlot slot;  ///< the engine's (group, level, channel); zeros for an inactive row
  LevelChannelInfo channel;  ///< pair key, l, 2s
  std::string widthText;     ///< the width column as written
  double width = 0.0;        ///< eV (open channel, photon), fm^-1/2 (closed channel: an ANC), MeV^1/2 with widthIsAmplitude
  bool widthIsAmplitude = false;  ///< gammaIsRWA = 1
  bool fixed = false;             ///< channelFix
};

/// One level: the rows of one level of a J^pi group (active rows), or consecutive inactive rows of one level ID.
struct LevelEntry {
  /// "<J^pi> <E_x as written>", e.g. "3/2- 3.5032"; " #2", " #3" ... appended to the later levels sharing a label
  std::string label;
  int twoJ = 0, parity = 0;
  std::string energyText;  ///< E_x (MeV) as written in the level's first row
  double energy = 0.0;     ///< MeV
  bool energyFixed = false;
  bool active = false;
  std::vector<int> ids;      ///< the level IDs of its rows (column 9)
  int group = 0, level = 0;  ///< the engine's indices (1-based); 0 for an inactive level
  int nucleusZ = 0, nucleusA = 0;
  std::vector<int> rows;  ///< indices into LevelRowTable::rows, in file order
};

/// The section as read.
struct LevelRowTable {
  std::vector<LevelRowLine> rows;
  std::vector<LevelEntry> levels;  ///< in the order of their first rows
  LevelTables tables;              ///< LevelTablesOfRows of the rows (pairs, the engine's groups)
  std::size_t open = 0, close = 0;  ///< 0-based lines of the <levels> and </levels> tags
};

/// Reads the <levels> section of `doc`: false with `why` when there is none, it is not closed, or a row does not read.
bool ReadLevelRowTable(const io::Document &doc, LevelRowTable &out, std::string &why);

/// "3/2-", "2+": J^pi as the label writes it.
std::string JpiText(int twoJ, int parity);

/// A level named by its label, by J^pi and E_x (within 1 keV, the engine's rule for one level, io::LevelGrouper):
/// "3/2- 3.50", or by the engine's indices (group > 0).
struct LevelAddress {
  std::string label;
  int group = 0, level = 0;
};
/// The index in table.levels of the level `address` names; -1 with `why` when none or more than one is named.
int FindLevel(const LevelRowTable &table, const LevelAddress &address, std::string &why);

/// Replace the document's lines [first, last) (0-based) by `lines`.
struct LineReplacement {
  int first = 0, last = 0;
  std::vector<std::string> lines;
};

/// The outcome of an edit: the replacements to apply in the order given (descending lines), or `error`.
struct LevelEdit {
  std::vector<LineReplacement> replacements;
  std::vector<std::pair<std::size_t, std::string>> unplaced;  ///< LevelValueEdits: (index of a value with no row, why)
  std::string error;
};

struct LevelChange {
  std::optional<double> energy;  ///< E_x (MeV), on every row of the level
  std::optional<bool> energyFixed;
  std::optional<bool> active;
};
struct ChannelChange {
  std::optional<double> width;  ///< the width column in the row's convention (eV, fm^-1/2, or MeV^1/2 with gammaIsRWA)
  std::optional<bool> fixed;
};
/// The pair columns, on every row of the pair (the separation energy also in the entrance-pair copy of the rows whose
/// entrance pair it is).
struct PairChange {
  std::optional<double> radius;  ///< fm
  std::optional<double> lightMass, heavyMass;  ///< u
  std::optional<double> separationEnergy;      ///< MeV
};
struct NewChannel {
  int pairKey = 0, l = 0, twoS = 0;
  double width = 0.0;
  bool fixed = false;
};
/// A level to add: its rows are those of its J^pi group (the group's channels, those not given with width 0, fixed),
/// or for a J^pi the project does not have, one row per channel cloned from a row of its pair (no channel given: the
/// group's channels, all width 0, fixed).  It gets the next free level ID and goes at the end of the section.
struct NewLevel {
  int twoJ = 0, parity = 1;
  double energy = 0.0;
  bool energyFixed = true;
  std::vector<NewChannel> channels;
};

LevelEdit SetLevel(const io::Document &doc, const LevelAddress &address, const LevelChange &change);
LevelEdit SetChannel(const io::Document &doc, const LevelAddress &address, const LevelChannelInfo &channel,
                     const ChannelChange &change);
LevelEdit SetPair(const io::Document &doc, int pairKey, const PairChange &change);
LevelEdit AddLevel(const io::Document &doc, const NewLevel &level);
LevelEdit RemoveLevel(const io::Document &doc, const LevelAddress &address);

/// A fitted value of a level parameter, named by the engine's indices of the fit parameter: the level's energy, or
/// the width of its channel (pair key, l, 2s).
struct LevelValue {
  enum class Kind { Energy, Width } kind = Kind::Energy;
  int group = 0, level = 0;
  LevelChannelInfo channel;  ///< Width only
  double value = 0.0;        ///< MeV for an energy; a width in the row's convention
};
/// The edits that write `values` into the rows they come from (pyazr's save_fit; the counterpart of
/// ParameterValueEdit for level parameters).  A value with no row (a channel the level does not list) is reported in
/// `unplaced` and not written.
LevelEdit LevelValueEdits(const io::Document &doc, const std::vector<LevelValue> &values);

}  // namespace experiment
}  // namespace azure2

#endif
