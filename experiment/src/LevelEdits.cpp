// The <levels> section as levels and channels, and its edits as line replacements (LevelEdits.h).
#include "azure2/experiment/LevelEdits.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>

#include "azure2/experiment/ProjectReader.h"
#include "azure2/experiment/ValueText.h"
#include "azure2/io/Document.h"
#include "azure2/io/LevelRowText.h"
#include "azure2/io/LevelRows.h"

namespace azure2 {
namespace experiment {

namespace {

// The table with what the edits need of each row: the row as read and its text.
struct Rows {
  LevelRowTable table;
  std::vector<io::LevelRow> row;
  std::vector<io::LevelRowText> text;
};

bool BlankLine(const std::string &s) { return s.find_first_not_of(" \t\r") == std::string::npos; }

int Column(const char *name) { return io::LevelColumnIndex(name); }

std::string Flag(bool on) { return on ? "1" : "0"; }

// The compound nucleus (Z, A) of a row's pair (PPair's constructor; LevelTablesOfRows groups by it).
std::pair<int, int> CompoundOf(const io::LevelRow &r) {
  if (r.pType() == 20) return {r.z2() - r.z1(), (int)std::floor(r.m2() + 0.5)};
  if (r.pType() == 10) return {r.z2(), (int)std::floor(r.m2() + 0.5)};
  return {r.z1() + r.z2(), (int)std::floor(r.m1() + r.m2() + 0.5)};
}

LevelEdit Refused(const std::string &why) {
  LevelEdit e;
  e.error = why;
  return e;
}

// One replacement per edited row, descending lines.
LevelEdit RowEdits(const Rows &t, const std::map<int, std::map<int, std::string>> &byRow) {
  LevelEdit e;
  for (auto it = byRow.rbegin(); it != byRow.rend(); ++it) {
    const io::LevelRowText &text = t.text[it->first];
    const std::string line = io::EditLevelRow(text, it->second);
    const int at = t.table.rows[it->first].line;
    if (line != text.line) e.replacements.push_back({at, at + 1, {line}});
  }
  std::sort(e.replacements.begin(), e.replacements.end(),
            [](const LineReplacement &a, const LineReplacement &b) { return a.first > b.first; });
  return e;
}

bool Finite(const std::optional<double> &x) { return !x || std::isfinite(*x); }

// The index of the first row of a pair (any row, active or not), or -1.
int PairRow(const Rows &t, int key) {
  for (std::size_t r = 0; r < t.row.size(); r++)
    if (t.row[r].ir() == key) return (int)r;
  return -1;
}

std::string ChannelText(const LevelChannelInfo &c) {
  std::ostringstream s;
  s << "pair " << c.pairKey << ", l = " << c.l << ", s = " << (c.twoS % 2 ? std::to_string(c.twoS) + "/2" : std::to_string(c.twoS / 2));
  return s.str();
}

// J as written in a label: "3/2" or "1.5" or "2".
bool ParseTwoJ(const std::string &j, int &twoJ) {
  char *end = nullptr;
  const std::size_t slash = j.find('/');
  if (slash != std::string::npos) {
    if (slash == 0 || j.substr(slash + 1) != "2") return false;
    const long v = std::strtol(j.substr(0, slash).c_str(), &end, 10);
    if (*end || v < 0) return false;
    twoJ = (int)v;
    return true;
  }
  const double x = std::strtod(j.c_str(), &end);
  if (j.empty() || *end || x < 0 || std::fabs(2.0 * x - std::lround(2.0 * x)) > 1e-9) return false;
  twoJ = (int)std::lround(2.0 * x);
  return true;
}

}  // namespace

std::string JpiText(int twoJ, int parity) {
  const std::string j = twoJ % 2 ? std::to_string(twoJ) + "/2" : std::to_string(twoJ / 2);
  return j + (parity < 0 ? "-" : "+");
}

namespace {

bool ReadRows(const io::Document &doc, Rows &t, std::string &why) {
  t = Rows();
  LevelRowTable &out = t.table;
  const io::Document::Section *s = doc.Find("levels");
  if (!s || s->close == io::Document::Section::npos) {
    why = "the project has no closed <levels> section";
    return false;
  }
  out.open = s->open;
  out.close = s->close;
  for (std::size_t i = s->open + 1; i < s->close; i++) {
    const std::string &line = doc.Line(i);
    if (BlankLine(line)) continue;  // ReadLevelSection skips blank lines (blanks and tabs)
    std::istringstream in(line);
    io::LevelRow row(in);
    io::LevelRowText text;
    if ((in.rdstate() & (std::ios::failbit | std::ios::badbit)) || !io::SplitLevelRow(line, text, why)) {
      why = "<levels> line " + std::to_string(i + 1) + " is not a row: '" + line + "'";
      return false;
    }
    t.row.push_back(row);
    t.text.push_back(text);
    LevelRowLine r;
    r.line = (int)i;
    r.channel = LevelChannelInfo{row.ir(), row.l(), (int)std::lround(2.0 * row.s())};
    r.widthText = text.Token(Column("gamma"));
    r.width = row.gamma();
    r.widthIsAmplitude = row.gammaIsRWA() == 1;
    r.fixed = row.channelFix() != 0;
    out.rows.push_back(r);
  }
  out.tables = LevelTablesOfRows(t.row);
  for (std::size_t k = 0; k < out.rows.size(); k++) out.rows[k].slot = out.tables.rowSlots[k];
  // Levels: an active row joins the level of its (group, level); an inactive row the level of the row before it when
  // that one is inactive with the same level ID and J^pi.
  std::map<std::pair<int, int>, int> ofSlot;
  int lastInactive = -1;
  for (std::size_t k = 0; k < out.rows.size(); k++) {
    const LevelRowSlot &slot = out.rows[k].slot;
    const io::LevelRow &row = t.row[k];
    const bool active = slot.group > 0;
    int index = -1;
    if (active) {
      const auto it = ofSlot.find({slot.group, slot.level});
      if (it != ofSlot.end()) index = it->second;
      lastInactive = -1;
    } else if (lastInactive >= 0) {
      const io::LevelRow &first = t.row[out.levels[lastInactive].rows.front()];
      if (first.levelID() == row.levelID() && first.levelJ() == row.levelJ() && first.levelPi() == row.levelPi())
        index = lastInactive;
    }
    if (index < 0) {
      LevelEntry l;
      l.twoJ = (int)std::lround(2.0 * row.levelJ());
      l.parity = row.levelPi();
      l.energyText = t.text[k].Token(Column("levelE"));
      l.energy = row.levelE();
      l.energyFixed = row.levelFix() != 0;
      l.active = active;
      l.group = slot.group;
      l.level = slot.level;
      const std::pair<int, int> nucleus = CompoundOf(row);
      l.nucleusZ = nucleus.first;
      l.nucleusA = nucleus.second;
      index = (int)out.levels.size();
      out.levels.push_back(l);
      if (active) ofSlot[{slot.group, slot.level}] = index;
    }
    if (!active) lastInactive = index;
    LevelEntry &l = out.levels[index];
    l.rows.push_back((int)k);
    if (std::find(l.ids.begin(), l.ids.end(), row.levelID()) == l.ids.end()) l.ids.push_back(row.levelID());
  }
  std::map<std::string, int> seen;
  for (LevelEntry &l : out.levels) {
    const std::string label = JpiText(l.twoJ, l.parity) + " " + l.energyText;
    const int n = ++seen[label];
    l.label = n == 1 ? label : label + " #" + std::to_string(n);
  }
  return true;
}

}  // namespace

bool ReadLevelRowTable(const io::Document &doc, LevelRowTable &out, std::string &why) {
  Rows t;
  const bool ok = ReadRows(doc, t, why);
  out = t.table;
  return ok;
}

int FindLevel(const LevelRowTable &t, const LevelAddress &a, std::string &why) {
  if (a.group > 0) {
    for (std::size_t k = 0; k < t.levels.size(); k++)
      if (t.levels[k].group == a.group && t.levels[k].level == a.level) return (int)k;
    why = "no level has the engine's indices (group " + std::to_string(a.group) + ", level " + std::to_string(a.level) + ")";
    return -1;
  }
  for (std::size_t k = 0; k < t.levels.size(); k++)
    if (t.levels[k].label == a.label) return (int)k;
  // "<J><+|-> <E_x>": J as 3/2 or 1.5
  std::istringstream in(a.label);
  std::string jpi, energyText, rest;
  in >> jpi >> energyText;
  char *end = nullptr;
  const double energy = std::strtod(energyText.c_str(), &end);
  int twoJ = -1, parity = 0;
  if (!jpi.empty() && (jpi.back() == '+' || jpi.back() == '-')) {
    parity = jpi.back() == '+' ? 1 : -1;
    if (!ParseTwoJ(jpi.substr(0, jpi.size() - 1), twoJ)) twoJ = -1;
  }
  std::string labels;
  for (const LevelEntry &l : t.levels) labels += (labels.empty() ? "" : "; ") + l.label;
  if (twoJ < 0 || energyText.empty() || *end || (in >> rest)) {
    why = "no level is named '" + a.label + "' (a level is named by its label or by J^pi and E_x: '3/2- 3.50'; the levels: " + labels + ")";
    return -1;
  }
  std::vector<int> found;
  for (std::size_t k = 0; k < t.levels.size(); k++)
    if (t.levels[k].twoJ == twoJ && t.levels[k].parity == parity && std::fabs(t.levels[k].energy - energy) <= 1e-3)
      found.push_back((int)k);
  if (found.size() == 1) return found[0];
  if (found.empty()) {
    why = "no " + JpiText(twoJ, parity) + " level within 1 keV of " + energyText + " MeV (the levels: " + labels + ")";
  } else {
    why = "'" + a.label + "' names " + std::to_string(found.size()) + " levels:";
    for (int k : found) why += " '" + t.levels[k].label + "'";
    why += "; give the label";
  }
  return -1;
}

LevelEdit SetLevel(const io::Document &doc, const LevelAddress &address, const LevelChange &change) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(why);
  const int k = FindLevel(t.table, address, why);
  if (k < 0) return Refused(why);
  if (!Finite(change.energy)) return Refused("the energy of a level is a finite number (MeV)");
  std::map<int, std::map<int, std::string>> byRow;
  for (int r : t.table.levels[k].rows) {
    std::map<int, std::string> &tokens = byRow[r];
    if (change.energy) tokens[Column("levelE")] = ShortestNumberText(*change.energy);
    if (change.energyFixed) tokens[Column("levelFix")] = Flag(*change.energyFixed);
    if (change.active) tokens[Column("isActive")] = Flag(*change.active);
  }
  return RowEdits(t, byRow);
}

LevelEdit SetChannel(const io::Document &doc, const LevelAddress &address, const LevelChannelInfo &channel,
                     const ChannelChange &change) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(why);
  const int k = FindLevel(t.table, address, why);
  if (k < 0) return Refused(why);
  if (!Finite(change.width)) return Refused("a width is a finite number");
  const LevelEntry &level = t.table.levels[k];
  std::string listed;
  for (int r : level.rows) {
    const LevelChannelInfo &c = t.table.rows[r].channel;
    listed += (listed.empty() ? "" : "; ") + ChannelText(c);
    if (c.pairKey != channel.pairKey || c.l != channel.l || c.twoS != channel.twoS) continue;
    std::map<int, std::map<int, std::string>> byRow;
    if (change.width) byRow[r][Column("gamma")] = ShortestNumberText(*change.width);
    if (change.fixed) byRow[r][Column("channelFix")] = Flag(*change.fixed);
    return RowEdits(t, byRow);
  }
  return Refused("level '" + level.label + "' has no row for " + ChannelText(channel) + " (its channels: " + listed + ")");
}

LevelEdit SetPair(const io::Document &doc, int pairKey, const PairChange &change) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(why);
  const int first = PairRow(t, pairKey);
  if (first < 0) {
    std::string keys;
    for (const PairInfo &p : t.table.tables.pairs)
      keys += (keys.empty() ? "" : ", ") + std::to_string(p.key) + " (" + PairLabel(p) + ")";
    return Refused("no pair has the key " + std::to_string(pairKey) + " (the pairs: " + keys + ")");
  }
  for (const std::optional<double> &x : {change.radius, change.lightMass, change.heavyMass})
    if (x && !(std::isfinite(*x) && *x > 0.0)) return Refused("a channel radius and a mass are positive numbers");
  if (!Finite(change.separationEnergy)) return Refused("a separation energy is a finite number (MeV)");
  if (change.radius && t.row[first].pType() != 0)
    return Refused("pair " + std::to_string(pairKey) + " is a photon or beta pair: it has no channel radius");
  std::map<int, std::map<int, std::string>> byRow;
  for (std::size_t r = 0; r < t.row.size(); r++) {
    const io::LevelRow &row = t.row[r];
    if (row.ir() == pairKey) {
      std::map<int, std::string> &tokens = byRow[(int)r];
      if (change.radius) tokens[Column("chRad")] = ShortestNumberText(*change.radius);
      if (change.lightMass) tokens[Column("m1")] = ShortestNumberText(*change.lightMass);
      if (change.heavyMass) tokens[Column("m2")] = ShortestNumberText(*change.heavyMass);
      if (change.separationEnergy) tokens[Column("sepE")] = ShortestNumberText(*change.separationEnergy);
    }
    if (change.separationEnergy && row.aa() == pairKey)
      byRow[(int)r][Column("entranceSepE")] = ShortestNumberText(*change.separationEnergy);
  }
  return RowEdits(t, byRow);
}

LevelEdit AddLevel(const io::Document &doc, const NewLevel &n) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(why);
  if (n.twoJ < 0 || (n.parity != 1 && n.parity != -1) || !std::isfinite(n.energy))
    return Refused("a level has 2J >= 0, parity +1 or -1 and a finite energy (MeV)");
  // The compound nucleus of the channels' pairs (no channel given: of the J^pi group that exists), and the J^pi group
  // of that nucleus.
  std::set<std::pair<int, int>> nuclei;
  for (const NewChannel &c : n.channels) {
    const int p = PairRow(t, c.pairKey);
    if (p < 0) return Refused("no pair has the key " + std::to_string(c.pairKey) + " (a new level names existing pairs)");
    if (!std::isfinite(c.width)) return Refused("a width is a finite number");
    nuclei.insert(CompoundOf(t.row[p]));
  }
  if (n.channels.empty())
    for (const LevelEntry &l : t.table.levels)
      if (l.active && l.twoJ == n.twoJ && l.parity == n.parity) nuclei.insert({l.nucleusZ, l.nucleusA});
  if (nuclei.empty())
    return Refused("no " + JpiText(n.twoJ, n.parity) + " level exists, so the new level needs its channels");
  if (nuclei.size() != 1)
    return Refused(n.channels.empty() ? "the " + JpiText(n.twoJ, n.parity) + " levels of several compound nuclei exist: "
                                        "name a channel of the new level"
                                      : "the channels of a level belong to pairs of one compound nucleus");
  int maxID = 0;
  for (const io::LevelRow &r : t.row) maxID = std::max(maxID, r.levelID());
  const std::string id = std::to_string(maxID + 1);
  const LevelEntry *group = nullptr;
  for (const LevelEntry &l : t.table.levels)
    if (l.active && l.twoJ == n.twoJ && l.parity == n.parity && std::make_pair(l.nucleusZ, l.nucleusA) == *nuclei.begin()) {
      group = &l;
      break;
    }
  // The template rows (index into the rows) and the channel given for each: one per channel of the J^pi group (its
  // first row of the channel), or one per given channel (the pair's first row).
  std::vector<std::pair<int, const NewChannel *>> rows;
  if (group) {
    std::map<int, int> ofChannel;  // the group's channel index -> its first row
    for (std::size_t r = 0; r < t.table.rows.size(); r++) {
      const LevelRowSlot &slot = t.table.rows[r].slot;
      if (slot.group == group->group && !ofChannel.count(slot.channel)) ofChannel[slot.channel] = (int)r;
    }
    std::vector<bool> used(n.channels.size(), false);
    std::string listed;
    for (const auto &kv : ofChannel) {
      const LevelChannelInfo &c = t.table.rows[kv.second].channel;
      listed += (listed.empty() ? "" : "; ") + ChannelText(c);
      const NewChannel *given = nullptr;
      for (std::size_t i = 0; i < n.channels.size(); i++)
        if (n.channels[i].pairKey == c.pairKey && n.channels[i].l == c.l && n.channels[i].twoS == c.twoS) {
          given = &n.channels[i];
          used[i] = true;
        }
      rows.push_back({kv.second, given});
    }
    for (std::size_t i = 0; i < n.channels.size(); i++)
      if (!used[i])
        return Refused("the " + JpiText(n.twoJ, n.parity) + " levels share one channel set (" + listed + "); " +
                       ChannelText(LevelChannelInfo{n.channels[i].pairKey, n.channels[i].l, n.channels[i].twoS}) +
                       " is not in it");
  } else {
    for (const NewChannel &c : n.channels) rows.push_back({PairRow(t, c.pairKey), &c});
  }
  std::vector<std::string> lines;
  for (const auto &rc : rows) {
    std::map<int, std::string> tokens;
    tokens[Column("levelJ")] = ShortestNumberText(0.5 * n.twoJ);
    tokens[Column("levelPi")] = n.parity > 0 ? "1" : "-1";
    tokens[Column("levelE")] = ShortestNumberText(n.energy);
    tokens[Column("levelFix")] = Flag(n.energyFixed);
    tokens[Column("levelID")] = id;
    tokens[Column("isActive")] = "1";
    tokens[Column("gamma")] = ShortestNumberText(rc.second ? rc.second->width : 0.0);
    tokens[Column("channelFix")] = Flag(rc.second ? rc.second->fixed : true);
    if (!group) {
      tokens[Column("l2")] = std::to_string(2 * rc.second->l);
      tokens[Column("s2")] = std::to_string(rc.second->twoS);
    }
    lines.push_back(io::EditLevelRow(t.text[rc.first], tokens));
  }
  lines.push_back("");
  const int close = (int)t.table.close;
  if (close - 1 > (int)t.table.open && !BlankLine(doc.Line(close - 1))) lines.insert(lines.begin(), "");
  LevelEdit e;
  e.replacements.push_back({close, close, lines});
  return e;
}

LevelEdit RemoveLevel(const io::Document &doc, const LevelAddress &address) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(why);
  const int k = FindLevel(t.table, address, why);
  if (k < 0) return Refused(why);
  std::vector<int> lines;
  for (int r : t.table.levels[k].rows) lines.push_back(t.table.rows[r].line);
  // Contiguous rows with a blank line after them and a blank line (or the tag) before: the blank line goes too.
  const int open = (int)t.table.open, close = (int)t.table.close;
  const bool contiguous = lines.back() - lines.front() + 1 == (int)lines.size();
  if (contiguous && lines.back() + 1 < close && BlankLine(doc.Line(lines.back() + 1)) &&
      (lines.front() - 1 == open || BlankLine(doc.Line(lines.front() - 1))))
    lines.push_back(lines.back() + 1);
  LevelEdit e;
  for (auto it = lines.rbegin(); it != lines.rend(); ++it) e.replacements.push_back({*it, *it + 1, {}});
  return e;
}

LevelEdit LevelValueEdits(const io::Document &doc, const std::vector<LevelValue> &values) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(why);
  std::map<int, std::map<int, std::string>> byRow;
  std::vector<std::pair<std::size_t, std::string>> unplaced;
  for (std::size_t i = 0; i < values.size(); i++) {
    const LevelValue &v = values[i];
    const int k = FindLevel(t.table, LevelAddress{"", v.group, v.level}, why);
    if (k < 0 || !std::isfinite(v.value)) {
      unplaced.push_back({i, k < 0 ? why : "a value that is not finite"});
      continue;
    }
    const std::string text = ShortestNumberText(v.value);
    if (v.kind == LevelValue::Kind::Energy) {
      for (int r : t.table.levels[k].rows) byRow[r][Column("levelE")] = text;
      continue;
    }
    bool placed = false;
    for (int r : t.table.levels[k].rows) {
      const LevelChannelInfo &c = t.table.rows[r].channel;
      if (c.pairKey == v.channel.pairKey && c.l == v.channel.l && c.twoS == v.channel.twoS) {
        byRow[r][Column("gamma")] = text;
        placed = true;
        break;
      }
    }
    if (!placed) unplaced.push_back({i, "level '" + t.table.levels[k].label + "' has no row for " + ChannelText(v.channel)});
  }
  LevelEdit e = RowEdits(t, byRow);
  e.unplaced = unplaced;
  return e;
}

}  // namespace experiment
}  // namespace azure2
