// The facade's part for the Potential panel of the Pairs tab and the Levels tab's check of the partial waves
// (docs/upgrade/qt_gui.md section 19).  Every number is the engine's: a pair's potential as the engine's reader builds
// the record (experiment::ReadFormat3Text, BuildPotentials), the pair as the engine makes it from its <levels> row
// (PPair(io::LevelRow)), the wave's arguments as the session's calls take them (PotentialWaveOf, api/include/
// PotentialWave.h), and physics/'s functions of the potential (InteriorBoxStates, SubtractedBoxStates,
// PotentialEffectiveRange, DepthForScatteringLength, DepthBindingAt, NearestEigenvalue, PotentialBoundStateANC).  The
// partial waves' rule is the engine's set-up (CompoundSystem::Fill and AddHardSphereGroups, src/CompoundSystemBuild.cpp).
#include "azure2/front/Engine.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <tuple>

#include "Constants.h"
#include "PPair.h"
#include "PotentialWave.h"
#include "azure2/experiment/Build.h"
#include "azure2/experiment/ChannelRadii.h"
#include "azure2/experiment/LevelTables.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/Document.h"
#include "azure2/io/LevelRows.h"
#include "azure2/physics/HybridPotential.h"
#include "azure2/physics/PotentialTuning.h"
#include "azure2/physics/Widths.h"

namespace azure2 {
namespace front {

namespace {

namespace ex = azure2::experiment;
namespace ph = azure2::physics;

const char *const kSection = "potentials";

/// A <levels> row as the engine's reader takes it (the tokens joined by blanks); false when it does not read.
bool ReadLevelRow(const LevelRow &row, std::unique_ptr<io::LevelRow> &out) {
  std::string text;
  for (const std::string &t : row.tokens) {
    if (t.empty()) continue;
    if (!text.empty()) text += ' ';
    text += t;
  }
  std::istringstream in(text);
  out = std::make_unique<io::LevelRow>(in);
  return !in.fail();
}

/// The compound nucleus of a row as the engine and the Pairs tab infer it from the pair's columns (PPair): Z1 + Z2,
/// round(M1 + M2) of a particle pair, Z2, round(M2) of a photon pair, Z2 - Z1, round(M2) of a beta pair.
std::pair<int, int> NucleusOf(const io::LevelRow &r) {
  if (r.pType() == 10) return {r.z2(), (int)std::lround(r.m2())};
  if (r.pType() == 20) return {r.z2() - r.z1(), (int)std::lround(r.m2())};
  return {r.z1() + r.z2(), (int)std::lround(r.m1() + r.m2())};
}

std::vector<io::LevelRow> ReadLevelRows(const std::vector<LevelRow> &rows) {
  std::vector<io::LevelRow> out;
  for (const LevelRow &r : rows) {
    std::unique_ptr<io::LevelRow> row;
    if (ReadLevelRow(r, row)) out.push_back(*row);
  }
  return out;
}

int TwoJOf(const std::string &spin) {
  if (spin.empty()) return -1;
  const std::size_t slash = spin.find('/');
  if (slash != std::string::npos) return (int)std::lround(2.0 * std::stod(spin.substr(0, slash)) / std::stod(spin.substr(slash + 1)));
  return (int)std::lround(2.0 * std::stod(spin));
}

/// The project with only what a pair's potential needs: its <config>, <levels> and <channelRadii>, and in
/// <potentials> the one record `record`, on; every other section of records emptied (no data file is read).
Project PotentialCopy(const Project &project, const Record &record) {
  Project copy = Project::FromText(project.Text(), project.Path());
  for (const std::string &section : copy.Sections())
    if (section != "config" && section != "levels" && section != "channelRadii" && !copy.Records(section).empty())
      copy.SetSectionBody(section, {});
  copy.SetSectionBody(kSection, {});
  std::string why;
  copy.SetRecord(kSection, record, why);
  return copy;
}

ex::Format3Project Built(const Project &copy) {
  ex::Format3ReadOptions read;
  read.fileName = copy.Path();
  read.projectDir = copy.Directory();
  read.readDataFiles = false;
  read.dataFiles = ex::DataFiles::NotUsed;
  read.setUpChecks = false;
  return ex::ReadFormat3Text(copy.Text(), read);
}

std::string ErrorsOf(const ex::Format3Project &built) {
  std::string out;
  for (const ex::Diagnostic &d : built.diagnostics)
    if (d.severity == ex::Severity::Error)  // without the copy's file and line: the record and field name it
      out += (out.empty() ? "" : "; ") + d.code + " (" + d.record + (d.field.empty() ? "" : "." + d.field) + "): " + d.message;
  return out;
}

/// A wave of a pair's potential set up as the engine sets it up: the pair from its <levels> row with the record's
/// potential, and the wave's arguments (PotentialWaveOf).
struct WaveSetUp {
  std::unique_ptr<PPair> pair;
  std::unique_ptr<io::LevelRow> row;
  PotentialWaveArgs args;
  ex::Format3Project built;
  double radius = 0.0;
  std::string radiusNote;
  int twoJ = -1;  ///< the 2J asked (the lowest of the wave's when it names none and the potential depends on J)
};

bool SetUpWave(const Project &project, const Record &record, int pairKey, int l, int twoJ, WaveSetUp &w, std::string &why) {
  w.built = Built(PotentialCopy(project, record));
  const std::string errors = ErrorsOf(w.built);
  if (!errors.empty()) {
    why = errors;
    return false;
  }
  const auto setting = w.built.built.potentials.pairs.find(pairKey);
  if (setting == w.built.built.potentials.pairs.end()) {
    why = "the potential " + record.name + " is not the potential of pair " + std::to_string(pairKey);
    return false;
  }
  if (!w.built.levelSection) {
    why = "the project has no <levels>";
    return false;
  }
  for (const io::LevelRow &r : w.built.levelSection->rows)
    if (r.ir() == pairKey && r.pType() == 0) {
      w.row = std::make_unique<io::LevelRow>(r);
      break;
    }
  if (!w.row) {
    why = "pair " + std::to_string(pairKey) + " is no particle pair of <levels>";
    return false;
  }
  try {
    w.pair = std::make_unique<PPair>(*w.row, PairFunctionSettings());
    w.pair->SetHybridPotential(setting->second.Potential());
    // a wave without J of a potential that depends on J (a `wave` names a J): its lowest J, |l - s| of the smallest s
    if (twoJ < 0 && ph::DependsOnJ(*w.pair->GetHybridPotential())) {
      const int a = (int)std::lround(2 * w.row->j1()), b = (int)std::lround(2 * w.row->j2());
      for (int twoS = std::abs(a - b); twoS <= a + b; twoS += 2)
        if (twoJ < 0 || std::abs(2 * l - twoS) < twoJ) twoJ = std::abs(2 * l - twoS);
    }
    w.twoJ = twoJ;
    w.args = PotentialWaveOf(*w.pair, l, twoJ, "the potential panel");
  } catch (const std::exception &e) {
    why = e.what();
    return false;
  }
  // the channel radius: the pair's, or a <channelRadii> record's for (l, s, J); channel spins s that differ: the first
  w.radius = w.pair->GetPairRadius();
  const int twoJ1 = (int)std::lround(2 * w.row->j1()), twoJ2 = (int)std::lround(2 * w.row->j2());
  twoJ = w.twoJ;
  bool first = true;
  for (int twoS = std::abs(twoJ1 - twoJ2); twoS <= twoJ1 + twoJ2; twoS += 2) {
    if (twoJ >= 0 && (twoJ < std::abs(2 * l - twoS) || twoJ > 2 * l + twoS)) continue;
    bool ambiguous = false;
    const ex::ChannelRadiusEntry *e = ex::ChannelRadiusFor(w.built.built.channelRadii, pairKey, l, twoS, twoJ, &ambiguous);
    const double radius = e ? e->radius_fm : w.pair->GetPairRadius();
    if (first) {
      w.radius = radius;
      if (e) w.radiusNote = "the channel's own radius (channel_radius " + e->record + ")";
      first = false;
    } else if (radius != w.radius) {
      w.radiusNote += (w.radiusNote.empty() ? "" : "; ") + std::string("its channel spins have different radii: the first is shown");
    }
  }
  return true;
}

}  // namespace

bool ParsePotentialWave(const std::string &text, PotentialWave &out, std::string &why) {
  out = PotentialWave();
  ValueParts parts;
  if (!ParseFieldValue("potential", "wave", text, "", parts, why)) return false;
  for (const NamedValue &p : parts.parts) {
    if (p.name == "l") out.l = std::stoi(p.value.text);
    else if (p.name == "J") {
      out.J = p.value.text;
      out.twoJ = TwoJOf(out.J);
    } else if (p.name == "depth") {
      out.depth = p.value.number;
      out.depthUnit = p.value.unit;
    } else if (p.name == "width") {
      out.width = p.value.number;
      out.widthUnit = p.value.unit;
    } else if (p.name == "background") out.background = true;
    else if (p.name == "boundary") out.boundary = p.value.text.empty() ? p.value.number : p.value.text;
    else if (p.name == "subtract") out.subtract = p.value.text;
    else if (p.name == "capture") out.capture = true;
    else if (p.name == "only") out.boundStateOnly = true;
  }
  return true;
}

std::string PotentialWaveText(const PotentialWave &wave) {
  // the parts in the alternative's order, each filled as its type wants it (the enum words from the schema's part)
  const std::vector<KindInfo> kinds = RecordKinds();
  const FieldInfo *field = nullptr;
  for (const KindInfo &k : kinds)
    if (k.kind == "potential")
      for (const FieldInfo &f : k.fields)
        if (f.key == "wave") field = &f;
  if (!field || field->value.alternatives.empty()) return std::string();
  ValueParts parts;
  parts.type = field->value.type;
  parts.alternative = field->value.alternatives.front().keyword;
  for (const PartInfo &part : field->value.alternatives.front().parts) {
    ValueParts v;
    v.type = part.value.type;
    const std::string word = part.value.values.empty() ? std::string() : part.value.values.front();
    bool given = true;
    if (part.name == "l") v.text = std::to_string(wave.l);
    else if (part.name == "J") given = !wave.J.empty(), v.text = wave.J;
    else if (part.name == "depth") given = !wave.depth.empty(), v.number = wave.depth, v.unit = wave.depthUnit;
    else if (part.name == "width") given = !wave.width.empty(), v.number = wave.width, v.unit = wave.widthUnit;
    else if (part.name == "background") given = wave.background, v.text = word;
    else if (part.name == "boundary") given = wave.background && !wave.boundary.empty(), v.text = v.number = wave.boundary;
    else if (part.name == "subtract") given = wave.background && !wave.subtract.empty(), v.text = wave.subtract;
    else if (part.name == "capture") given = wave.background && wave.capture, v.text = word;
    else if (part.name == "only") given = wave.boundStateOnly, v.text = word;
    else given = false;
    if (given) parts.parts.push_back({part.name, v});
  }
  return ValueText(field->value, parts);
}

PairPotential PotentialOfPair(const Project &project, int pairKey) {
  PairPotential out;
  for (bool on : {true, false}) {
    for (const Record &r : on ? project.Records(kSection) : project.RecordsOff(kSection)) {
      if (r.kind != "potential") continue;
      const ex::Format3Project built = Built(PotentialCopy(project, r));
      const auto it = built.built.potentialRecords.find(pairKey);
      if (it == built.built.potentialRecords.end() || it->second.first != r.name) continue;
      out.found = true;
      out.on = on;
      out.record = r;
      return out;
    }
  }
  return out;
}

bool PotentialWaveOf(const Project &project, const Record &record, int pairKey, int l, int twoJ, double eMin, double eMax,
                     PotentialWaveReport &out, std::string &why) {
  out = PotentialWaveReport();
  WaveSetUp w;
  if (!SetUpWave(project, record, pairKey, l, twoJ, w, why)) return false;
  const PotentialWaveArgs &a = w.args;
  const int z1 = w.pair->GetZ(1), z2 = w.pair->GetZ(2);
  out.l = l;
  out.twoJ = w.twoJ;
  out.parity = w.row->pi1() * w.row->pi2() * (l % 2 == 0 ? 1 : -1);
  out.radius = w.radius;
  out.radiusNote = w.radiusNote;
  out.threshold = w.row->sepE() + w.row->e2() + w.row->e1();
  out.depth = a.potential.central_depth_MeV;
  out.width = a.potential.central.shape == ph::FormShape::Gaussian ? a.potential.central.width_fm : 0.0;
  const ph::InteriorBackground *bg = ph::BackgroundOf(*w.pair->GetHybridPotential(), l, a.twoJ);
  out.background = bg != nullptr;
  out.boundary = bg ? bg->boundary : 0.0;
  try {
    const std::vector<ph::BoxState> states =
        ph::InteriorBoxStates(a.potential, z1, z2, a.mu, l, a.ldoti, out.radius, out.boundary, eMin, eMax);
    const double gammaW2 = ph::wigner_limit(w.pair->GetRedMass(), out.radius);
    for (const ph::BoxState &b : states) {
      PotentialPole p;
      p.nodes = b.nodes;
      p.energy = b.E_MeV;
      p.excitation = b.E_MeV + out.threshold;
      p.gamma2 = b.gamma2_MeV;
      p.wignerRatio = b.gamma2_MeV / gammaW2;
      p.subtracted = bg && std::find(bg->subtract_nodes.begin(), bg->subtract_nodes.end(), b.nodes) != bg->subtract_nodes.end();
      out.poles.push_back(p);
    }
  } catch (const std::exception &e) {
    why = e.what();
    return false;
  }
  // the set-up's checks of a background (CheckPotentialWaves, src/CompoundSystemChecks.cpp), with its words
  if (bg && out.radius > a.potential.matching_radius_fm)
    out.radiusError = "an interior background at the channel radius " + std::to_string(out.radius) +
                      " fm, beyond the matching radius " + std::to_string(a.potential.matching_radius_fm) +
                      " fm, where the potential does not act: it would be ignored.  Raise matching_radius or lower the "
                      "channel radius.";
  if (bg && !bg->subtract_nodes.empty()) {
    try {
      ph::SubtractedBoxStates(a.potential, z1, z2, a.mu, l, a.ldoti, out.radius, *bg);
    } catch (const std::exception &e) {
      out.subtractError = "subtract names a box state the potential does not have there (" + std::string(e.what()) +
                          ").  Name the node number of a state the box has (potential_background.md section 3).";
    }
  }
  const ph::BoundStateCount count = ph::PotentialBoundStates(a.potential, z1, z2, a.mu, l, a.ldoti);
  out.boundStates = count.bound;
  out.zeroEnergyNodes = count.nodes;
  const ph::EffectiveRange e = ph::PotentialEffectiveRange(a.potential, z1, z2, a.mu, l, a.ldoti);
  out.effectiveRange = e.ok;
  out.effectiveRangeWhy = e.why;
  out.a = e.a;
  out.r = e.r;
  out.aSpread = e.aSpread;
  out.rSpread = e.rSpread;
  out.eLow = e.eLow;
  out.eHigh = e.eHigh;
  // the project's levels of this J^pi with a channel of the pair and l
  for (const LevelRow &row : project.Levels()) {
    std::unique_ptr<io::LevelRow> r;
    if (!ReadLevelRow(row, r) || r->ir() != pairKey || r->pType() != 0 || r->l() != l || r->levelPi() != out.parity) continue;
    if (twoJ >= 0 && (int)std::lround(2 * r->levelJ()) != twoJ) continue;
    const ChannelQuantities q = ChannelAtLevel(row, out.radius);
    PotentialLevel lv;
    lv.J = r->levelJ();
    lv.parity = r->levelPi();
    lv.excitation = r->levelE();
    lv.energy = r->levelE() - out.threshold;
    lv.active = r->isActive() == 1;
    lv.s = r->s();
    lv.kind = q.kind;
    lv.value = r->gamma();
    out.levels.push_back(lv);
  }
  return true;
}

bool PotentialDepthForScatteringLength(const Project &project, const Record &record, int pairKey, int l, int twoJ,
                                       double a, PotentialDepth &out, std::string &why) {
  out = PotentialDepth();
  WaveSetUp w;
  if (!SetUpWave(project, record, pairKey, l, twoJ, w, why)) return false;
  const ph::DepthForLength d = ph::DepthForScatteringLength(w.args.potential, w.pair->GetZ(1), w.pair->GetZ(2), w.args.mu,
                                                            l, w.args.ldoti, a);
  out.found = d.found;
  out.why = d.why;
  out.depth = d.depth_MeV;
  out.a = d.a;
  out.r = d.r;
  out.boundStates = d.boundStates;
  out.nodes = d.nodes;
  return true;
}

bool PotentialDepthForBinding(const Project &project, const Record &record, int pairKey, int l, int twoJ, double energy,
                              PotentialDepth &out, std::string &why) {
  out = PotentialDepth();
  if (!(energy < 0.0)) {
    why = "the energy must be below the pair's threshold (E_c < 0)";
    return false;
  }
  WaveSetUp w;
  if (!SetUpWave(project, record, pairKey, l, twoJ, w, why)) return false;
  const PotentialWaveArgs &a = w.args;
  const int z1 = w.pair->GetZ(1), z2 = w.pair->GetZ(2);
  // the radius of the bound-state rule (GetPotentialBoundState): the matching radius when the channel's is beyond it
  const double a_fm = std::min(w.radius, a.potential.matching_radius_fm);
  const ph::DepthAtEnergy d = ph::DepthBindingAt(a.potential, z1, z2, a.mu, l, a.ldoti, a_fm, energy);
  if (!d.found) {
    out.why = "no depth within 30 % of " + std::to_string(a.potential.central_depth_MeV) + " MeV binds the wave at " +
              std::to_string(energy) + " MeV";
    return true;
  }
  ph::HybridPotential q = a.potential;
  q.central_depth_MeV = d.depth_MeV;
  const ph::PotentialEigenvalue ev = ph::NearestEigenvalue(q, z1, z2, a.mu, l, a.ldoti, a_fm, energy, std::min(0.5, -energy));
  out.found = ev.found;
  out.depth = d.depth_MeV;
  if (!ev.found) {
    out.why = "the depth " + std::to_string(d.depth_MeV) + " MeV was found but its eigenvalue not (window 0.5 MeV)";
    return true;
  }
  out.energy = ev.E_MeV;
  out.anc = ph::PotentialBoundStateANC(q, z1, z2, a.mu, l, a.ldoti, a_fm, ev.E_MeV);
  ph::IntegrateOutward(q, z1, z2, a.mu, l, a.ldoti, q.matching_radius_fm, ev.E_MeV, &out.nodes);
  out.boundStates = ph::PotentialBoundStates(q, z1, z2, a.mu, l, a.ldoti).bound;
  return true;
}

namespace {

/// The largest l the hard-sphere groups take (CompoundSystem::Fill): maxOrbitalMomentum when given, else the largest l
/// of a particle row of <levels> (0 when none).
int KeyMaxL(const std::vector<io::LevelRow> &rows, int maxOrbitalMomentum) {
  if (maxOrbitalMomentum >= 0) return maxOrbitalMomentum;
  int maxL = 0;
  for (const io::LevelRow &r : rows)
    if (r.pType() == 0 && r.l() > maxL) maxL = r.l();
  return maxL;
}

}  // namespace

std::vector<PartialWaveState> PartialWaveStates(const std::vector<LevelRow> &rows, bool includeAllJPi,
                                                int maxOrbitalMomentum, int pairKey, int l, int twoJ) {
  std::vector<PartialWaveState> out;
  const std::vector<io::LevelRow> read = ReadLevelRows(rows);
  const io::LevelRow *pair = nullptr;
  for (const io::LevelRow &r : read)
    if (r.ir() == pairKey && r.pType() == 0) pair = &r;
  if (!pair) return out;
  const int parity = pair->pi1() * pair->pi2() * (l % 2 == 0 ? 1 : -1);
  std::set<int> twoJs;
  if (twoJ >= 0) {
    twoJs.insert(twoJ);
  } else {
    const int a = (int)std::lround(2 * pair->j1()), b = (int)std::lround(2 * pair->j2());
    for (int twoS = std::abs(a - b); twoS <= a + b; twoS += 2)
      for (int j = std::abs(2 * l - twoS); j <= 2 * l + twoS; j += 2) twoJs.insert(j);
  }
  const bool keyCovers = l <= KeyMaxL(read, maxOrbitalMomentum);
  for (int j : twoJs) {
    PartialWaveState s;
    s.twoJ = j;
    s.parity = parity;
    for (const io::LevelRow &r : read)
      s.level = s.level || (r.isActive() == 1 && r.ir() == pairKey && r.l() == l && r.levelPi() == parity &&
                            (int)std::lround(2 * r.levelJ()) == j);
    // the key adds a group only for a J^pi without one (an active level of the J^pi without the channel keeps it out)
    bool group = false;
    for (const io::LevelRow &r : read)
      group = group || (r.isActive() == 1 && r.levelPi() == parity && (int)std::lround(2 * r.levelJ()) == j &&
                        NucleusOf(r) == NucleusOf(*pair));
    s.key = keyCovers && !group;
    s.included = s.level || (includeAllJPi && s.key);
    out.push_back(s);
  }
  return out;
}

std::vector<MissingPartialWave> PartialWavesWithoutLevel(const std::vector<LevelRow> &rows, int lMax,
                                                         int maxOrbitalMomentum, int &keyMaxL) {
  std::vector<MissingPartialWave> out;
  const std::vector<io::LevelRow> read = ReadLevelRows(rows);
  keyMaxL = read.empty() ? -1 : KeyMaxL(read, maxOrbitalMomentum);
  // one row per particle pair, in key order
  std::map<int, const io::LevelRow *> pairs;
  for (const io::LevelRow &r : read)
    if (r.pType() == 0 && !pairs.count(r.ir())) pairs[r.ir()] = &r;
  // J^pi -> its lowest l and a pair, per compound nucleus (Z, A of the pair's two particles)
  struct Key {
    int Z, A, twoJ, parity;
    bool operator<(const Key &o) const {
      return std::tie(Z, A, twoJ, parity) < std::tie(o.Z, o.A, o.twoJ, o.parity);
    }
  };
  std::set<int> identicalPairs;
  {
    std::string text = "<levels>\n";
    for (const LevelRow &r : rows) {
      std::string line;
      for (const std::string &t : r.tokens)
        if (!t.empty()) line += (line.empty() ? "" : " ") + t;
      text += line + "\n";
    }
    text += "</levels>\n";
    bool ok = false;
    const ex::LevelTables tables = ex::LevelTablesOfRows(ex::LevelRowsOf(io::Document::Parse(text), ok));
    for (const ex::PairInfo &p : tables.pairs)
      if (p.identical) identicalPairs.insert(p.key);
  }
  std::map<Key, MissingPartialWave> formed;
  for (const auto &kv : pairs) {
    const io::LevelRow &p = *kv.second;
    const int Z = NucleusOf(p).first, A = NucleusOf(p).second;
    const int a = (int)std::lround(2 * p.j1()), b = (int)std::lround(2 * p.j2());
    // identical particles: only l + s even (the engine's set-up; the pair's flag as the engine's tables give it)
    const bool identical = identicalPairs.count(kv.first) > 0;
    for (int twoS = std::abs(a - b); twoS <= a + b; twoS += 2)
      for (int l = 0; l <= lMax; l++) {
        if (identical && (l + twoS / 2) % 2 != 0) continue;
        const int parity = p.pi1() * p.pi2() * (l % 2 == 0 ? 1 : -1);
        for (int j = std::abs(2 * l - twoS); j <= 2 * l + twoS; j += 2) {
          const Key k{Z, A, j, parity};
          auto it = formed.find(k);
          if (it == formed.end() || l < it->second.l)
            formed[k] = {j, parity, Z, A, kv.first, l, false};
        }
      }
  }
  for (auto &kv : formed) {
    MissingPartialWave &m = kv.second;
    bool carried = false;
    // an active level of the J^pi in the nucleus (each row of a level carries the level's J^pi and its pair's columns)
    for (const io::LevelRow &r : read)
      carried = carried || (r.isActive() == 1 && (int)std::lround(2 * r.levelJ()) == m.twoJ && r.levelPi() == m.parity &&
                            NucleusOf(r) == std::make_pair(m.nucleusZ, m.nucleusA));
    if (carried) continue;
    m.key = keyMaxL >= 0 && m.l <= keyMaxL;
    out.push_back(m);
  }
  return out;
}

}  // namespace front
}  // namespace azure2
