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

LevelEdit Refused(const EditRefusal &why) {
  LevelEdit e;
  e.refusal = why;
  e.error = why.Text();
  return e;
}

// The section does not read (ReadRows's `why`).
EditRefusal Unreadable(const std::string &why) {
  return {why, "a closed <levels> section whose lines are rows of the <levels> columns or blank",
          "correct the section in the project's text before editing it"};
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
    r.pairRow = row.isActive() != 1 && row.levelID() == 0;
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
    if (out.rows[k].pairRow) {  // holds a pair, not a level
      lastInactive = -1;
      continue;
    }
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

std::string EditRefusal::Text() const {
  std::string s = what;
  if (!allowed.empty()) s += (s.empty() ? "" : ": ") + allowed;
  if (!fix.empty()) s += (s.empty() ? "" : "; ") + fix;
  return s;
}

int FindLevel(const LevelRowTable &t, const LevelAddress &a, EditRefusal &why) {
  std::string labels;
  for (const LevelEntry &l : t.levels) labels += (labels.empty() ? "" : "; ") + l.label;
  const std::string listed = "the levels: " + (labels.empty() ? std::string("none") : labels);
  if (a.group > 0) {
    for (std::size_t k = 0; k < t.levels.size(); k++)
      if (t.levels[k].group == a.group && t.levels[k].level == a.level) return (int)k;
    why = {"no level has the engine's indices (group " + std::to_string(a.group) + ", level " +
               std::to_string(a.level) + ")",
           "the indices of an active level (" + listed + ")", "name the level by its label (Project.levels())"};
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
  if (twoJ < 0 || energyText.empty() || *end || (in >> rest)) {
    why = {"no level is named '" + a.label + "'",
           "a level is named by its label or by J^pi and E_x, '3/2- 3.50' (" + listed + ")",
           "give one of the labels (Project.levels())"};
    return -1;
  }
  std::vector<int> found;
  for (std::size_t k = 0; k < t.levels.size(); k++)
    if (t.levels[k].twoJ == twoJ && t.levels[k].parity == parity && std::fabs(t.levels[k].energy - energy) <= 1e-3)
      found.push_back((int)k);
  if (found.size() == 1) return found[0];
  if (found.empty()) {
    why = {"no " + JpiText(twoJ, parity) + " level within 1 keV of " + energyText + " MeV",
           "J^pi and E_x of a level within 1 keV (" + listed + ")", "give one of the labels (Project.levels())"};
  } else {
    std::string named;
    for (int k : found) named += (named.empty() ? "'" : ", '") + t.levels[k].label + "'";
    why = {"'" + a.label + "' names " + std::to_string(found.size()) + " levels (" + named + ")",
           "a name of one level", "give the label of the one meant"};
  }
  return -1;
}

namespace {

// Appends the removals of the row lines `lines` (0-based, ascending) to `e`: contiguous rows with a blank line after
// them and a blank line (or the tag) before take the blank line with them, so that an added block (AddLevel, AddPair)
// removed again leaves the section as it was.
void RemoveRowLines(const io::Document &doc, const Rows &t, std::vector<int> lines, LevelEdit &e) {
  const int open = (int)t.table.open, close = (int)t.table.close;
  const bool contiguous = lines.back() - lines.front() + 1 == (int)lines.size();
  if (contiguous && lines.back() + 1 < close && BlankLine(doc.Line(lines.back() + 1)) &&
      (lines.front() - 1 == open || BlankLine(doc.Line(lines.front() - 1))))
    lines.push_back(lines.back() + 1);
  for (auto it = lines.rbegin(); it != lines.rend(); ++it) e.replacements.push_back({*it, *it + 1, {}});
}

std::string PairsListed(const Rows &t) {
  std::string keys;
  for (const PairInfo &p : t.table.tables.pairs)
    keys += (keys.empty() ? "" : ", ") + std::to_string(p.key) + " (" + PairLabel(p) + ")";
  return keys.empty() ? "none" : keys;
}

bool HalfInteger(int twoJ) { return twoJ >= 0; }

}  // namespace

LevelEdit AddPair(const io::Document &doc, const NewPair &n, int &key) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(Unreadable(why));
  if (n.kind != 0 && n.kind != 10 && n.kind != 20)
    return Refused({"the pair kind " + std::to_string(n.kind), "particle (0), photon (10) or beta (20)",
                    "give one of them"});
  const bool particle = n.kind == 0, photon = n.kind == 10;
  for (double x : {n.massLight, n.massHeavy, n.gLight, n.gHeavy, n.lightExcitation, n.heavyExcitation,
                   n.separationEnergy, n.radius, n.bindingEnergy})
    if (!std::isfinite(x))
      return Refused({"a number of the pair is not finite", "finite masses, g-factors and energies",
                      "give a number for each"});
  if (!HalfInteger(n.twoJLight) || !HalfInteger(n.twoJHeavy))
    return Refused({"the spins 2j = " + std::to_string(n.twoJLight) + ", " + std::to_string(n.twoJHeavy),
                    "a spin is a multiple of 1/2, >= 0", "give J as 0, 1/2, 1, ..."});
  if (std::abs(n.parityLight) != 1 || std::abs(n.parityHeavy) != 1)
    return Refused({"the parities " + std::to_string(n.parityLight) + ", " + std::to_string(n.parityHeavy),
                    "a parity is +1 or -1", "give +1 or -1"});
  if (!(n.massHeavy > 0.0) || (!photon && !(n.massLight > 0.0)))
    return Refused({"the masses " + ShortestNumberText(n.massLight) + ", " + ShortestNumberText(n.massHeavy) + " u",
                    "the masses of a pair's partners are > 0 (u; a photon's is 0, given by the pair kind)",
                    "give the atomic masses in u"});
  if (particle && !(n.radius > 0.0))
    return Refused({"the channel radius " + ShortestNumberText(n.radius) + " fm",
                    "a particle pair's channel radius is a number > 0 (fm)", "give radius= in fm"});
  if (!particle && n.radius != 0.0)
    return Refused({"the channel radius " + ShortestNumberText(n.radius) + " fm of a photon or beta pair",
                    "a photon or beta pair has no channel radius", "give 0 or none"});
  if (n.lightExcitation < 0.0 || n.heavyExcitation < 0.0)
    return Refused({"the excitation energies " + ShortestNumberText(n.lightExcitation) + ", " +
                        ShortestNumberText(n.heavyExcitation) + " MeV",
                    "an excitation energy is >= 0 (MeV, 0 for the ground state)", "give it in MeV"});
  if (!particle && n.lightExcitation != 0.0)
    return Refused({"a light excitation of a photon or beta pair",
                    "only a particle pair's light partner can be excited",
                    "give light_excitation=0"});
  if (n.bindingEnergy < 0.0 || (!particle && n.bindingEnergy != 0.0))
    return Refused({"the THM binding energy " + ShortestNumberText(n.bindingEnergy) + " MeV",
                    "a number >= 0 (MeV; 0 = none) on a particle pair", "give binding_energy >= 0 or none"});
  if (n.externalCapture < 0)
    return Refused({"the external-capture mask " + std::to_string(n.externalCapture), "a whole number >= 0",
                    "give 0 for none"});
  int maxKey = 0;
  for (const io::LevelRow &r : t.row) maxKey = std::max(maxKey, r.ir());
  key = n.key > 0 ? n.key : maxKey + 1;
  if (n.key < 0)
    return Refused({"the pair key " + std::to_string(n.key), "a whole number >= 1 (0 or none: the next free key)",
                    "give another key or none"});
  if (PairRow(t, key) >= 0)
    return Refused({"pair key " + std::to_string(key) + " is taken (the pairs: " + PairsListed(t) + ")",
                    "a key no pair has", "give another key or none for the next free one (" +
                    std::to_string(maxKey + 1) + ")"});
  // The unused columns aa and entranceSepE as the first row writes them (the engine reads neither).
  const std::string aa = t.row.empty() ? "1" : t.text[0].Token(Column("aa"));
  const std::string entranceSepE = t.row.empty() ? ShortestNumberText(n.separationEnergy)
                                                 : t.text[0].Token(Column("entranceSepE"));
  auto spin = [](int twoJ) { return ShortestNumberText(0.5 * twoJ); };
  const std::vector<std::string> columns = {
      "0", "1", "0", "1", aa, std::to_string(key), "0", "0", "0", "0", "1", "0",
      photon ? "1" : spin(n.twoJLight), photon ? "1" : std::to_string(n.parityLight), spin(n.twoJHeavy),
      std::to_string(n.parityHeavy), ShortestNumberText(n.heavyExcitation), photon ? "0" : ShortestNumberText(n.massLight),
      ShortestNumberText(n.massHeavy), photon ? "0" : std::to_string(n.zLight), std::to_string(n.zHeavy), entranceSepE,
      ShortestNumberText(n.separationEnergy), "0", "0", "0.0", std::to_string(n.kind), ShortestNumberText(n.radius),
      photon ? "0" : ShortestNumberText(n.gLight), ShortestNumberText(n.gHeavy), std::to_string(n.externalCapture)};
  std::string line;
  for (const std::string &c : columns) line += (line.empty() ? "" : " ") + c;
  if (n.bindingEnergy != 0.0) line += "  " + ShortestNumberText(n.bindingEnergy) + " 0";
  if (n.lightExcitation != 0.0) line += " " + ShortestNumberText(n.lightExcitation);
  std::vector<std::string> lines{line, ""};
  const int close = (int)t.table.close;
  if (close - 1 > (int)t.table.open && !BlankLine(doc.Line(close - 1))) lines.insert(lines.begin(), "");
  LevelEdit e;
  e.replacements.push_back({close, close, lines});
  return e;
}

LevelEdit RemovePair(const io::Document &doc, int pairKey, bool removeChannels) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(Unreadable(why));
  if (PairRow(t, pairKey) < 0)
    return Refused({"no pair has the key " + std::to_string(pairKey),
                    "the key of a pair (the pairs: " + PairsListed(t) + ")",
                    "give one of the keys (Project.pairs())"});
  std::vector<int> lines;
  std::string users, bare;
  for (const LevelEntry &l : t.table.levels) {
    std::string named;
    int others = 0;
    for (int r : l.rows) {
      if (t.table.rows[r].channel.pairKey != pairKey) {
        others++;
        continue;
      }
      const LevelChannelInfo &c = t.table.rows[r].channel;
      named += (named.empty() ? "" : "; ") + std::string("l = ") + std::to_string(c.l) + ", s = " +
               (c.twoS % 2 ? std::to_string(c.twoS) + "/2" : std::to_string(c.twoS / 2));
      lines.push_back(t.table.rows[r].line);
    }
    if (named.empty()) continue;
    users += (users.empty() ? "" : ", ") + std::string("'") + l.label + "' (" + named + ")";
    if (others == 0) bare += (bare.empty() ? "" : ", ") + std::string("'") + l.label + "'";
  }
  const std::string label = PairLabel(*PairByKey(pairKey, t.table.tables.pairs));
  if (!users.empty() && !removeChannels)
    return Refused({"the pair (" + label + ") is a channel of the levels " + users,
                    "a pair is removed when no level names it",
                    "remove those levels or channels first, or remove the channels with the pair "
                    "(remove_channels=True)"});
  if (!bare.empty())
    return Refused({"the pair (" + label + ") is the only channel of the levels " + bare,
                    "a level keeps at least one channel", "remove those levels first (remove_level)"});
  for (std::size_t r = 0; r < t.row.size(); r++)
    if (t.table.rows[r].pairRow && t.row[r].ir() == pairKey) lines.push_back(t.table.rows[r].line);
  std::sort(lines.begin(), lines.end());
  LevelEdit e;
  // Each run of consecutive lines is removed as a block (with its blank line, as RemoveLevel does).
  for (std::size_t a = 0; a < lines.size();) {
    std::size_t b = a + 1;
    while (b < lines.size() && lines[b] == lines[b - 1] + 1) b++;
    RemoveRowLines(doc, t, std::vector<int>(lines.begin() + a, lines.begin() + b), e);
    a = b;
  }
  std::sort(e.replacements.begin(), e.replacements.end(),
            [](const LineReplacement &x, const LineReplacement &y) { return x.first > y.first; });
  return e;
}

AllowedChannelList AllowedChannels(const io::Document &doc, int twoJ, int parity, const ChannelLimits &limits,
                                   const std::vector<int> &pairKeys) {
  AllowedChannelList out;
  Rows t;
  auto refuse = [&out](const EditRefusal &why) {
    out.refusal = why;
    out.error = why.Text();
    return out;
  };
  std::string why;
  if (!ReadRows(doc, t, why)) return refuse(Unreadable(why));
  if (twoJ < 0 || std::abs(parity) != 1)
    return refuse({"2J = " + std::to_string(twoJ) + ", parity " + std::to_string(parity),
                   "a level has 2J >= 0 and parity +1 or -1", "give J as 0, 1/2, 1, ... and parity +1 or -1"});
  if (limits.lMax < 0 || limits.multipolarityMax < 0 || limits.multipolaritiesMax < 0)
    return refuse({"a negative limit", "l_max, multipolarity_max and multipolarities_max are whole numbers >= 0",
                   "give each >= 0 (their defaults are 2)"});
  for (int k : pairKeys)
    if (!PairByKey(k, t.table.tables.pairs))
      return refuse({"no pair has the key " + std::to_string(k),
                     "the key of a pair (the pairs: " + PairsListed(t) + ")",
                     "give keys of the pairs, or none for all"});
  for (const PairInfo &p : t.table.tables.pairs) {
    if (!pairKeys.empty() && std::find(pairKeys.begin(), pairKeys.end(), p.key) == pairKeys.end()) continue;
    if (p.kind == 0) {
      for (int twoS = std::abs(p.twoJHeavy - p.twoJLight); twoS <= p.twoJHeavy + p.twoJLight; twoS += 2)
        for (int twoL = std::abs(twoS - twoJ); twoL <= twoS + twoJ; twoL += 2) {
          const int l = twoL / 2;
          if (twoL % 2 == 0 && p.lightParity * p.heavyParity * (l % 2 ? -1 : 1) == parity && l <= limits.lMax)
            out.channels.push_back({p.key, l, twoS});
        }
    } else if (p.kind == 20) {
      if (p.heavyParity == parity && p.twoJHeavy == twoJ) {
        out.channels.push_back({p.key, 1, 0});
        out.channels.push_back({p.key, 0, 0});
      } else if (p.heavyParity == parity && std::abs(p.twoJHeavy - twoJ) == 2) {
        out.channels.push_back({p.key, 1, 2});
      }
    } else {
      int count = 0;
      for (int L = 1; L <= limits.multipolarityMax && count < limits.multipolaritiesMax; L++)
        if (std::abs(2 * L - p.twoJHeavy) <= twoJ && twoJ <= 2 * L + p.twoJHeavy) {
          out.channels.push_back({p.key, L, p.twoJHeavy});
          count++;
        }
    }
  }
  return out;
}

LevelEdit SetLevel(const io::Document &doc, const LevelAddress &address, const LevelChange &change) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(Unreadable(why));
  EditRefusal unknown;
  const int k = FindLevel(t.table, address, unknown);
  if (k < 0) return Refused(unknown);
  if (!Finite(change.energy))
    return Refused({"the energy " + ShortestNumberText(*change.energy) + " MeV", "a finite number (MeV)",
                    "give E_x in MeV"});
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
  if (!ReadRows(doc, t, why)) return Refused(Unreadable(why));
  EditRefusal unknown;
  const int k = FindLevel(t.table, address, unknown);
  if (k < 0) return Refused(unknown);
  if (!Finite(change.width))
    return Refused({"the width " + ShortestNumberText(*change.width), "a finite number in the row's convention",
                    "give the width in eV, fm^-1/2 (an ANC) or MeV^1/2 (an amplitude row)"});
  const LevelEntry &level = t.table.levels[k];
  std::string listed;
  for (int r : level.rows) {
    const LevelChannelInfo &c = t.table.rows[r].channel;
    listed += (listed.empty() ? "" : "; ") + ChannelText(c);
    if (c.pairKey != channel.pairKey || c.l != channel.l || c.twoS != channel.twoS) continue;
    if (change.widthIsAmplitude && t.row[r].pType() != 0)
      return Refused({"width_is_amplitude on a photon or beta channel",
                      "the flag is for a particle channel (a photon or beta channel has no reduced-width amplitude "
                      "input)",
                      "set width_is_amplitude on a channel of a particle pair"});
    std::map<int, std::map<int, std::string>> byRow;
    if (change.width) byRow[r][Column("gamma")] = ShortestNumberText(*change.width);
    if (change.fixed) byRow[r][Column("channelFix")] = Flag(*change.fixed);
    if (change.widthIsAmplitude && (*change.widthIsAmplitude || t.text[r].hasThm))
      byRow[r][Column("gammaIsRWA")] = Flag(*change.widthIsAmplitude);
    return RowEdits(t, byRow);
  }
  return Refused({"level '" + level.label + "' has no row for " + ChannelText(channel),
                  "a channel of the level (its channels: " + listed + ")",
                  "give one of its channels, or add the channel to the levels of its J^pi"});
}

LevelEdit SetPair(const io::Document &doc, int pairKey, const PairChange &change) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(Unreadable(why));
  const int first = PairRow(t, pairKey);
  if (first < 0)
    return Refused({"no pair has the key " + std::to_string(pairKey),
                    "the key of a pair (the pairs: " + PairsListed(t) + ")",
                    "give one of the keys (Project.pairs())"});
  struct Positive {
    const std::optional<double> &value;
    const char *what, *unit;
  };
  for (const Positive &x : {Positive{change.radius, "channel radius", "fm"}, Positive{change.lightMass, "light mass", "u"},
                            Positive{change.heavyMass, "heavy mass", "u"}})
    if (x.value && !(std::isfinite(*x.value) && *x.value > 0.0))
      return Refused({"the " + std::string(x.what) + " " + ShortestNumberText(*x.value) + " " + x.unit +
                          " is out of range",
                      "it is a number > 0 (" + std::string(x.unit) + ")", "give a positive " + std::string(x.what)});
  if (!Finite(change.separationEnergy))
    return Refused({"the separation energy " + ShortestNumberText(*change.separationEnergy) + " MeV is out of range",
                    "it is a finite number (MeV)", "give the separation energy in MeV"});
  for (const Positive &x : {Positive{change.lightExcitation, "light excitation energy", "MeV"},
                            Positive{change.heavyExcitation, "heavy excitation energy", "MeV"}})
    if (x.value && !(std::isfinite(*x.value) && *x.value >= 0.0))
      return Refused({"the " + std::string(x.what) + " " + ShortestNumberText(*x.value) + " MeV is out of range",
                      "it is a number >= 0 (MeV, 0 for the ground state)", "give it in MeV"});
  if (change.radius && t.row[first].pType() != 0)
    return Refused({"a channel radius for a photon or beta pair", "a channel radius belongs to a particle pair",
                    "set radius on a particle pair"});
  if (change.bindingEnergy && t.row[first].pType() != 0)
    return Refused({"a THM binding energy for a photon or beta pair",
                    "a photon or beta pair cannot be the entrance pair of a THM experiment, so it has no THM binding "
                    "energy",
                    "set binding_energy on a particle pair"});
  if (change.bindingEnergy && !(std::isfinite(*change.bindingEnergy) && *change.bindingEnergy >= 0.0))
    return Refused({"the THM binding energy " + ShortestNumberText(*change.bindingEnergy) + " MeV is out of range",
                    "it is a number >= 0 (MeV; 0 = none), B of the Trojan horse a = x + s",
                    "give B in MeV, or 0 for none"});
  std::map<int, std::map<int, std::string>> byRow;
  for (std::size_t r = 0; r < t.row.size(); r++) {
    const io::LevelRow &row = t.row[r];
    if (row.ir() == pairKey) {
      std::map<int, std::string> &tokens = byRow[(int)r];
      if (change.radius) tokens[Column("chRad")] = ShortestNumberText(*change.radius);
      if (change.lightMass) tokens[Column("m1")] = ShortestNumberText(*change.lightMass);
      if (change.heavyMass) tokens[Column("m2")] = ShortestNumberText(*change.heavyMass);
      if (change.separationEnergy) tokens[Column("sepE")] = ShortestNumberText(*change.separationEnergy);
      if (change.heavyExcitation) tokens[Column("e2")] = ShortestNumberText(*change.heavyExcitation);
      if (change.lightExcitation && (*change.lightExcitation != 0.0 || t.text[r].hasE1))
        tokens[Column("e1")] = ShortestNumberText(*change.lightExcitation);
      if (change.bindingEnergy && (*change.bindingEnergy != 0.0 || t.text[r].hasThm))
        tokens[Column("bindingE")] = ShortestNumberText(*change.bindingEnergy);
    }
    if (change.separationEnergy && row.aa() == pairKey)
      byRow[(int)r][Column("entranceSepE")] = ShortestNumberText(*change.separationEnergy);
  }
  return RowEdits(t, byRow);
}

LevelEdit AddLevel(const io::Document &doc, const NewLevel &n) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(Unreadable(why));
  if (n.twoJ < 0 || (n.parity != 1 && n.parity != -1) || !std::isfinite(n.energy))
    return Refused({"2J = " + std::to_string(n.twoJ) + ", parity " + std::to_string(n.parity) + ", E_x " +
                        ShortestNumberText(n.energy) + " MeV",
                    "a level has 2J >= 0, parity +1 or -1 and a finite energy (MeV)",
                    "give J as 0, 1/2, 1, ..., parity +1 or -1 and E_x in MeV"});
  // The compound nucleus of the channels' pairs (no channel given: of the J^pi group that exists), and the J^pi group
  // of that nucleus.
  std::set<std::pair<int, int>> nuclei;
  for (const NewChannel &c : n.channels) {
    const int p = PairRow(t, c.pairKey);
    if (p < 0)
      return Refused({"no pair has the key " + std::to_string(c.pairKey),
                      "a new level names existing pairs (the pairs: " + PairsListed(t) + ")",
                      "give one of the keys, or add the pair first (add_pair)"});
    if (!std::isfinite(c.width))
      return Refused({"the width " + ShortestNumberText(c.width) + " of the channel of pair " +
                          std::to_string(c.pairKey),
                      "a finite number in the row's convention", "give the width in eV, fm^-1/2 or MeV^1/2"});
    nuclei.insert(CompoundOf(t.row[p]));
  }
  if (n.channels.empty())
    for (const LevelEntry &l : t.table.levels)
      if (l.active && l.twoJ == n.twoJ && l.parity == n.parity) nuclei.insert({l.nucleusZ, l.nucleusA});
  if (nuclei.empty())
    return Refused({"no channels for a " + JpiText(n.twoJ, n.parity) + " level, and no " + JpiText(n.twoJ, n.parity) +
                        " level exists",
                    "a level of a new J^pi needs its channels",
                    "give its channels, or channels=\"allowed\" (Project.allowed_channels)"});
  if (nuclei.size() != 1)
    return Refused(n.channels.empty()
                       ? EditRefusal{"the " + JpiText(n.twoJ, n.parity) + " levels of several compound nuclei exist",
                                     "the levels of one compound nucleus", "name a channel of the new level"}
                       : EditRefusal{"the channels name pairs of several compound nuclei",
                                     "the channels of a level belong to pairs of one compound nucleus",
                                     "give channels of one compound nucleus"});
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
        return Refused({ChannelText(LevelChannelInfo{n.channels[i].pairKey, n.channels[i].l, n.channels[i].twoS}) +
                            " is not a channel of the " + JpiText(n.twoJ, n.parity) + " levels",
                        "the " + JpiText(n.twoJ, n.parity) + " levels share one channel set (" + listed + ")",
                        "give channels of the set; the others are added with width 0, fixed"});
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
  // The pair rows of the pairs the new rows name: the pair is now held by a level.
  std::set<int> named;
  for (const auto &rc : rows) named.insert(t.row[rc.first].ir());
  for (std::size_t r = t.row.size(); r-- > 0;)
    if (t.table.rows[r].pairRow && named.count(t.row[r].ir())) RemoveRowLines(doc, t, {t.table.rows[r].line}, e);
  return e;
}

LevelEdit RemoveLevel(const io::Document &doc, const LevelAddress &address) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(Unreadable(why));
  EditRefusal unknown;
  const int k = FindLevel(t.table, address, unknown);
  if (k < 0) return Refused(unknown);
  std::vector<int> lines;
  for (int r : t.table.levels[k].rows) lines.push_back(t.table.rows[r].line);
  LevelEdit e;
  RemoveRowLines(doc, t, lines, e);
  return e;
}

LevelEdit LevelValueEdits(const io::Document &doc, const std::vector<LevelValue> &values) {
  Rows t;
  std::string why;
  if (!ReadRows(doc, t, why)) return Refused(Unreadable(why));
  std::map<int, std::map<int, std::string>> byRow;
  std::vector<std::pair<std::size_t, std::string>> unplaced;
  for (std::size_t i = 0; i < values.size(); i++) {
    const LevelValue &v = values[i];
    EditRefusal unknown;
    const int k = FindLevel(t.table, LevelAddress{"", v.group, v.level}, unknown);
    if (k < 0 || !std::isfinite(v.value)) {
      unplaced.push_back({i, k < 0 ? unknown.what : "a value that is not finite"});
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

MergedLevelRows MergeLevelRows(const std::vector<const io::Document *> &members) {
  MergedLevelRows out;
  int keyOffset = 0, idOffset = 0;
  for (std::size_t k = 0; k < members.size(); k++) {
    Rows t;
    std::string why;
    if (!ReadRows(*members[k], t, why)) {
      out.error = "project " + std::to_string(k + 1) + ": " + why;
      return out;
    }
    std::map<int, int> keys;
    int maxKey = 0, maxID = 0;
    for (const io::LevelRow &r : t.row) {
      keys[r.ir()] = r.ir() + keyOffset;
      maxKey = std::max(maxKey, r.ir());
      maxID = std::max(maxID, r.levelID());
    }
    for (std::size_t r = 0; r < t.row.size(); r++)
      out.rows.push_back(io::EditLevelRow(t.text[r], {{Column("ir"), std::to_string(keys[t.row[r].ir()])},
                                                      {Column("levelID"), std::to_string(t.row[r].levelID() + idOffset)}}));
    out.pairKeys.push_back(keys);
    keyOffset += maxKey;
    idOffset += maxID;
  }
  return out;
}

}  // namespace experiment
}  // namespace azure2
