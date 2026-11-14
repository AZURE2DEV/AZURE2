#ifndef AZURE2_EXPERIMENT_LEVELEDITS_H
#define AZURE2_EXPERIMENT_LEVELEDITS_H

#include <map>
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
  bool pairRow = false;           ///< a pair row (AddPair): the row of a pair no level names, not a level's
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
  /// gammaIsRWA (column 33): the width column is a reduced-width amplitude (MeV^1/2) instead of a partial width or
  /// ANC; a particle channel only.  The value is not converted: give the width in the new convention with it.
  std::optional<bool> widthIsAmplitude;
};
/// The pair columns, on every row of the pair (the separation energy also in the entrance-pair copy of the rows whose
/// entrance pair it is).
struct PairChange {
  std::optional<double> radius;  ///< fm
  std::optional<double> lightMass, heavyMass;  ///< u
  std::optional<double> separationEnergy;      ///< MeV
  /// The excitation energies (MeV) of the light (e1, the optional last column) and the heavy partner (e2); a light
  /// excitation of 0 is not written on a row that does not write the column.
  std::optional<double> lightExcitation, heavyExcitation;
  /// The THM binding energy B (MeV, >= 0; 0 = none) of the composite projectile whose transferred particle is the
  /// pair's light partner (column 32); a particle pair only.  A row that does not write the THM columns gets them
  /// only for a value other than 0.
  std::optional<double> bindingEnergy;
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

/// A pair to add: the pair columns of a <levels> row.  A photon pair's light partner is the photon (J 1, parity +1,
/// Z 0, mass 0, g 0, whatever is given); a beta pair's light partner is the electron or positron as given.
struct NewPair {
  int key = 0;   ///< 0: the next free key (one more than the largest)
  int kind = 0;  ///< pType: 0 particle pair, 10 photon (capture to a level of the compound nucleus), 20 beta
  int twoJLight = 1, twoJHeavy = 0;  ///< 2j of the light and heavy partner (photon pair: heavy = the final level)
  int parityLight = 1, parityHeavy = 1;
  int zLight = 0, zHeavy = 0;
  double massLight = 0.0, massHeavy = 0.0;  ///< u
  double gLight = 0.0, gHeavy = 0.0;        ///< g-factors (magnetic moments in nuclear magnetons over j)
  double lightExcitation = 0.0, heavyExcitation = 0.0;  ///< MeV (photon pair: E_x of the final level)
  double separationEnergy = 0.0;  ///< MeV
  double radius = 0.0;            ///< fm, the channel radius (a particle pair; 0 for photon and beta)
  int externalCapture = 0;        ///< ecMultMask, the external-capture multipole mask (column 31)
  double bindingEnergy = 0.0;     ///< MeV, the THM binding energy (a particle pair; 0 = none)
};
/*!
 * The <levels> file has no place for a pair but its rows, so a pair no level names yet is held by one row of its own:
 * inactive (isActive 0) with level ID 0, J 0+, E_x 0, l 0, s 0, width 0 (a "pair row").  ReadLevelRowTable lists its
 * pair and no level of it; the engine reads it as any inactive row (its pair is part of the system only with
 * includeAllJPi).  AddLevel removes the pair row of a pair its new rows name.  The pair row goes at the end of the
 * section; `key` gets the key used.
 */
LevelEdit AddPair(const io::Document &doc, const NewPair &pair, int &key);
/*!
 * Removes pair `pairKey`: its pair row, and with `removeChannels` the channel rows of the levels that name it.  Refused
 * while a level names it (without `removeChannels`), naming the levels and channels, and when a level would be left
 * without a channel (remove the level first).  Data sets and records that name the pair are the validator's to refuse.
 */
LevelEdit RemovePair(const io::Document &doc, int pairKey, bool removeChannels);

/// The limits of AllowedChannels (the Qt LevelsTab's spin boxes and their defaults).
struct ChannelLimits {
  int lMax = 2;                ///< the largest orbital angular momentum of a particle channel
  int multipolarityMax = 2;    ///< the largest multipolarity L of a photon channel
  int multipolaritiesMax = 2;  ///< the most photon channels per photon pair (the lowest multipolarities first)
};
/*!
 * The channels a level of 2J `twoJ` and parity `parity` has with the pairs of the section (keys in `pairKeys`, or all
 * when empty), in pair order then (s, l) ascending: the rule of the Qt GUI's LevelsTab::calculateChannels.
 * - particle pair (j1 pi1, j2 pi2): every channel spin s = |j1 - j2| ... j1 + j2 and every whole l = |s - J| ... s + J
 *   with pi1 pi2 (-1)^l = pi and l <= lMax;
 * - photon pair (final level J_f): s = J_f and the multipolarities L = 1 ... multipolarityMax with |L - J_f| <= J <=
 *   L + J_f, at most multipolaritiesMax of them (E or M follows from the parities: the engine's, not a column);
 * - beta pair (parent j2 pi2): pi = pi2 and J = j2: s = 0 with l = 1 (Gamow-Teller) and l = 0 (Fermi); |J - j2| = 1:
 *   s = 1, l = 1 (Gamow-Teller).
 * Refused (`error`) for 2J < 0, a parity other than +1 or -1, negative limits, or an unknown key.
 */
struct AllowedChannelList {
  std::vector<LevelChannelInfo> channels;
  std::string error;
};
AllowedChannelList AllowedChannels(const io::Document &doc, int twoJ, int parity, const ChannelLimits &limits,
                                   const std::vector<int> &pairKeys);

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

/// The <levels> rows of several projects in one section (pyazr's merge_projects: the levels of several compound
/// nuclei in one project).
struct MergedLevelRows {
  std::vector<std::string> rows;             ///< every row of every member, member by member, in file order
  std::vector<std::map<int, int>> pairKeys;  ///< per member: the pair keys its rows name -> the keys in `rows`
  std::string error;                         ///< the member whose section does not read, and why; "" when merged
};
/*!
 * Member k's rows with each pair key (the ir column) moved past the keys of the members before it -- by the sum of
 * their largest keys -- and each level ID past the largest ID before it; every other byte of a row is kept (inactive
 * rows too).  The aa column is kept as written: the engine does not read it.
 */
MergedLevelRows MergeLevelRows(const std::vector<const io::Document *> &members);

}  // namespace experiment
}  // namespace azure2

#endif
