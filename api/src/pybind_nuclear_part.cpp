// The nuclear part of a format-3 project for pyazr's Project (docs/upgrade/pyazr_on_records.md, "Edits of the nuclear
// part"): the <levels> rows as levels and pairs (experiment/LevelEdits.h), the <config> entries (io/ConfigSection.h)
// and the <thm> lines (ThmBlockEdits.h), read and edited by the engine's one grammar of each.  An edit returns the
// replacements of the document's lines, (first, last, [lines]) with 0-based [first, last), to be applied in the order
// given; pyazr applies them to its Document and validates the project as it would stand.  A refusal is a ValueError
// with the reason.  Texts are bytes.
#include <algorithm>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ThmBlockEdits.h"
#include "ThmExperiment.h"
#include "ThmKeys.h"
#include "azure2/experiment/LevelEdits.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/ConfigSection.h"
#include "azure2/io/Document.h"

namespace py = pybind11;

namespace {

namespace io = azure2::io;
namespace ex = azure2::experiment;

py::list Replacements(const std::vector<ex::LineReplacement> &edits) {
  py::list out;
  for (const ex::LineReplacement &r : edits) {
    py::list lines;
    for (const std::string &l : r.lines) lines.append(py::bytes(l));
    out.append(py::make_tuple(r.first, r.last, lines));
  }
  return out;
}

// A refused edit as a ValueError whose attributes what, allowed and fix are the refusal's parts (None when empty);
// pyazr's FormatError takes them from there.
[[noreturn]] void ThrowRefusal(const ex::EditRefusal &r) {
  auto part = [](const std::string &s) { return s.empty() ? py::object(py::none()) : py::object(py::str(s)); };
  py::object error = py::reinterpret_borrow<py::object>(PyExc_ValueError)(r.Text());
  error.attr("what") = part(r.what);
  error.attr("allowed") = part(r.allowed);
  error.attr("fix") = part(r.fix);
  PyErr_SetObject(PyExc_ValueError, error.ptr());
  throw py::error_already_set();
}

py::list Done(const ex::LevelEdit &e) {
  if (!e.error.empty()) ThrowRefusal(e.refusal);
  return Replacements(e.replacements);
}

// A level named by its label (str) or by the engine's indices ((group, level)).
ex::LevelAddress Address(const py::object &level) {
  ex::LevelAddress a;
  if (py::isinstance<py::str>(level) || py::isinstance<py::bytes>(level)) {
    a.label = level.cast<std::string>();
  } else {
    const auto gl = level.cast<std::pair<int, int>>();
    a.group = gl.first;
    a.level = gl.second;
    if (a.group < 1) throw py::value_error("a level's engine indices start at 1");
  }
  return a;
}

py::dict LevelTable(const io::Document &doc) {
  ex::LevelRowTable t;
  std::string why;
  if (!ex::ReadLevelRowTable(doc, t, why)) throw py::value_error(why);
  py::list levels, pairs;
  for (const ex::LevelEntry &l : t.levels) {
    py::dict d;
    d["label"] = l.label;
    d["two_j"] = l.twoJ;
    d["parity"] = l.parity;
    d["jpi"] = ex::JpiText(l.twoJ, l.parity);
    d["energy_text"] = l.energyText;
    d["energy"] = l.energy;
    d["fixed"] = l.energyFixed;
    d["active"] = l.active;
    d["ids"] = l.ids;
    d["group"] = l.group;
    d["level"] = l.level;
    d["nucleus"] = py::make_tuple(l.nucleusZ, l.nucleusA);
    py::list channels;
    for (int r : l.rows) {
      const ex::LevelRowLine &row = t.rows[r];
      py::dict c;
      c["pair"] = row.channel.pairKey;
      c["l"] = row.channel.l;
      c["two_s"] = row.channel.twoS;
      c["width_text"] = row.widthText;
      c["width"] = row.width;
      c["fixed"] = row.fixed;
      c["amplitude"] = row.widthIsAmplitude;
      c["channel"] = row.slot.channel;
      c["line"] = row.line + 1;
      channels.append(c);
    }
    d["channels"] = channels;
    levels.append(d);
  }
  for (const ex::PairInfo &p : t.tables.pairs) {
    py::dict d;
    d["key"] = p.key;
    d["label"] = ex::PairLabel(p);
    d["kind"] = p.kind;
    d["light"] = py::make_tuple(p.lightZ, p.lightA);
    d["heavy"] = py::make_tuple(p.heavyZ, p.heavyA);
    d["light_mass"] = p.lightMass;
    d["heavy_mass"] = p.heavyMass;
    d["light_excitation"] = p.lightExcitation;
    d["heavy_excitation"] = p.heavyExcitation;
    d["separation_energy"] = p.separationEnergy;
    d["threshold"] = p.threshold;
    d["radius"] = p.radius;
    d["two_j_light"] = p.twoJLight;
    d["two_j_heavy"] = p.twoJHeavy;
    d["heavy_parity"] = p.heavyParity;
    d["light_parity"] = p.lightParity;
    d["binding_energy"] = p.bindingEnergy;
    pairs.append(d);
  }
  py::dict out;
  out["levels"] = levels;
  out["pairs"] = pairs;
  return out;
}

// The <config> lines of a document and the 0-based line of its tag; ValueError when there is no section.
std::vector<std::string> ConfigLines(const io::Document &doc, std::size_t &open) {
  std::vector<std::string> lines = io::ConfigLinesOf(doc, open);
  if (lines.empty()) throw py::value_error("the project has no <config> section");
  return lines;
}

// The lines between <thm> and </thm> (0-based line of the tag in `open`); empty with open = npos without a section.
std::vector<std::string> ThmBody(const io::Document &doc, std::size_t &open) {
  const io::Document::Section *s = doc.Find("thm");
  open = std::string::npos;
  if (!s) return {};
  if (s->close == io::Document::Section::npos) throw py::value_error("the <thm> section is not closed");
  open = s->open;
  return doc.Body(*s);
}

py::list ThmReplacements(const io::Document &doc, std::size_t open, const std::vector<std::string> &body,
                         const std::vector<ThmLineEdit> &edits) {
  if (open == std::string::npos) {  // a new section after </config> (or at the start)
    const io::Document::Section *c = doc.Find("config");
    const int at = c && c->close != io::Document::Section::npos ? (int)c->close + 1 : 0;
    std::vector<std::string> lines{"<thm>"};
    for (const std::string &l : ApplyThmLineEdits(body, edits)) lines.push_back(l);
    lines.push_back("</thm>");
    return Replacements({ex::LineReplacement{at, at, lines}});
  }
  // The body replaced as a whole only where it changed: one replacement per edit, in the order given.
  std::vector<ex::LineReplacement> out;
  for (const ThmLineEdit &e : edits) {
    const int line = (int)open + 1 + e.line;
    if (e.kind == ThmLineEdit::Kind::Replace)
      out.push_back({line, line + 1, {e.text}});
    else if (e.kind == ThmLineEdit::Kind::Remove)
      out.push_back({line, line + 1, {}});
    else
      out.push_back({line, line, {e.text}});
  }
  return Replacements(out);
}

// A key of <thm> as a field of a record kind of the schema (the shape pyazr.schema() gives a field).
py::dict ThmKeyField(const ThmKeyRow &r) {
  py::dict f;
  f["type"] = r.type;
  f["doc"] = r.doc;
  f["required"] = false;
  if (!r.title.empty() || !r.unit.empty()) f["title"] = r.title;
  if (!r.unit.empty()) f["unit"] = r.unit;
  if (!r.values.empty()) f["values"] = r.values;
  if (!r.defaultValue.empty()) f["default"] = r.defaultValue;
  if (!r.grammar.empty()) f["grammar"] = r.grammar;
  if (!r.presentWith.empty()) {
    py::dict when;
    when["if"] = "field_present";
    when["field"] = r.presentWith;
    when["text"] = r.presentWith + " is given";
    f["applies_when"] = when;
  }
  return f;
}

template <class Key>
py::dict ThmKeyKindDict(const ThmKeyKind &kind, const std::vector<Key> &keys) {
  py::dict fields;
  for (const Key &k : keys) fields[py::str(k.row.key)] = ThmKeyField(k.row);
  py::dict d;
  d["title"] = kind.title;
  d["doc"] = kind.doc;
  d["fields"] = fields;
  return d;
}

}  // namespace

void BindNuclearPart(py::module_ &m) {
  m.def("jpi_text", &ex::JpiText, py::arg("two_j"), py::arg("parity"),
        "J^pi as the engine writes it (experiment::JpiText): two_j = 2J, parity +1 or -1; \"1/2+\", \"2-\".");
  m.def("level_table", &LevelTable, py::arg("doc"),
        "The <levels> section as levels and pairs (experiment::ReadLevelRowTable): {'levels': [dicts with label, "
        "two_j, parity, jpi, energy_text, energy (MeV), fixed, active, ids, group and level (the engine's indices, 0 "
        "for an inactive level), nucleus (Z, A), channels: [dicts with pair, l, two_s, width_text, width, fixed, "
        "amplitude, channel (the group's channel index), line (1-based)]] (a pair row, experiment::AddPair, is no level), 'pairs': [dicts with key, label, kind, "
        "light and heavy (Z, A), light_mass, heavy_mass (u), light_excitation, heavy_excitation, separation_energy, "
        "threshold (MeV), radius (fm), two_j_light, two_j_heavy, heavy_parity, light_parity, binding_energy (MeV)]}.  ValueError when a row does not read.");
  m.def(
      "reaction_label",
      [](const io::Document &doc, int entrance, int exit) {
        ex::LevelRowTable t;
        std::string why;
        if (!ex::ReadLevelRowTable(doc, t, why)) throw py::value_error(why);
        const ex::PairInfo *in = ex::PairByKey(entrance, t.tables.pairs), *out = ex::PairByKey(exit, t.tables.pairs);
        if (!in || (exit != -1 && !out))
          throw py::value_error("no pair has the key " + std::to_string(in ? exit : entrance) + " (the pairs: " +
                                [&] {
                                  std::string l;
                                  for (const std::string &x : ex::PairListing(t.tables.pairs)) l += (l.empty() ? "" : "; ") + x;
                                  return l;
                                }() + ")");
        return ex::ReactionLabel(*in, exit == -1 ? nullptr : out);
      },
      py::arg("doc"), py::arg("entrance"), py::arg("exit"),
      "The reaction label from pair `entrance` to pair `exit` (keys of <levels>; exit -1: the total capture), as the "
      "converter writes it (experiment::ReactionLabel): '12C(p,g)13N'.");
  m.def(
      "find_level",
      [](const io::Document &doc, const py::object &level) {
        ex::LevelRowTable t;
        std::string why;
        if (!ex::ReadLevelRowTable(doc, t, why)) throw py::value_error(why);
        ex::EditRefusal unknown;
        const int k = ex::FindLevel(t, Address(level), unknown);
        if (k < 0) ThrowRefusal(unknown);
        return k;
      },
      py::arg("doc"), py::arg("level"),
      "The index in level_table()['levels'] of the level named by its label, by 'J^pi E_x' (within 1 keV) or by the "
      "engine's (group, level); ValueError when none or more than one.");
  m.def(
      "set_level_edits",
      [](const io::Document &doc, const py::object &level, std::optional<double> energy, std::optional<bool> fixed,
         std::optional<bool> active) {
        ex::LevelChange c;
        c.energy = energy;
        c.energyFixed = fixed;
        c.active = active;
        return Done(ex::SetLevel(doc, Address(level), c));
      },
      py::arg("doc"), py::arg("level"), py::kw_only(), py::arg("energy") = py::none(), py::arg("fixed") = py::none(),
      py::arg("active") = py::none(),
      "The replacements that set a level's E_x (MeV), its fixed flag or its active flag on every row of the level; "
      "the level by its label (str) or the engine's (group, level).");
  m.def(
      "set_channel_edits",
      [](const io::Document &doc, const py::object &level, int pair, int l, int twoS, std::optional<double> width,
         std::optional<bool> fixed, std::optional<bool> widthIsAmplitude) {
        ex::ChannelChange c;
        c.width = width;
        c.fixed = fixed;
        c.widthIsAmplitude = widthIsAmplitude;
        return Done(ex::SetChannel(doc, Address(level), ex::LevelChannelInfo{pair, l, twoS}, c));
      },
      py::arg("doc"), py::arg("level"), py::arg("pair"), py::arg("l"), py::arg("two_s"), py::kw_only(),
      py::arg("width") = py::none(), py::arg("fixed") = py::none(), py::arg("width_is_amplitude") = py::none(),
      "The replacements that set the width column (the row's convention: eV, fm^-1/2, or MeV^1/2 with gammaIsRWA), "
      "the fixed flag or the amplitude flag (gammaIsRWA, a particle channel only) of one channel (pair key, l, 2s) of "
      "a level.");
  m.def(
      "set_pair_edits",
      [](const io::Document &doc, int pair, std::optional<double> radius, std::optional<double> lightMass,
         std::optional<double> heavyMass, std::optional<double> separationEnergy, std::optional<double> lightExcitation,
         std::optional<double> heavyExcitation, std::optional<double> bindingEnergy) {
        ex::PairChange c;
        c.radius = radius;
        c.lightMass = lightMass;
        c.heavyMass = heavyMass;
        c.separationEnergy = separationEnergy;
        c.lightExcitation = lightExcitation;
        c.heavyExcitation = heavyExcitation;
        c.bindingEnergy = bindingEnergy;
        return Done(ex::SetPair(doc, pair, c));
      },
      py::arg("doc"), py::arg("pair"), py::kw_only(), py::arg("radius") = py::none(), py::arg("light_mass") = py::none(),
      py::arg("heavy_mass") = py::none(), py::arg("separation_energy") = py::none(),
      py::arg("light_excitation") = py::none(), py::arg("heavy_excitation") = py::none(),
      py::arg("binding_energy") = py::none(),
      "The replacements that set a pair's channel radius (fm), masses (u), separation energy, excitation energies "
      "of its light and heavy partner or THM binding energy (MeV, a particle pair only) on every row of the pair (the "
      "separation energy also in the entrance copy of the rows whose entrance pair it is).");
  m.def(
      "add_level_edits",
      [](const io::Document &doc, int twoJ, int parity, double energy, bool fixed,
         const std::vector<std::tuple<int, int, int, double, bool>> &channels) {
        ex::NewLevel n;
        n.twoJ = twoJ;
        n.parity = parity;
        n.energy = energy;
        n.energyFixed = fixed;
        for (const auto &c : channels)
          n.channels.push_back(ex::NewChannel{std::get<0>(c), std::get<1>(c), std::get<2>(c), std::get<3>(c), std::get<4>(c)});
        return Done(ex::AddLevel(doc, n));
      },
      py::arg("doc"), py::arg("two_j"), py::arg("parity"), py::arg("energy"), py::kw_only(), py::arg("fixed"),
      py::arg("channels"),
      "The replacements that add a level (2J, parity, E_x in MeV, fixed) with channels [(pair, l, 2s, width, fixed)]: "
      "the rows of its J^pi group (channels not given: width 0, fixed), or one row per channel cloned from a row of "
      "its pair; the next free level ID; at the end of <levels>.");
  m.def(
      "add_pair_edits",
      [](const io::Document &doc, int key, int kind, int twoJLight, int parityLight, int zLight, double massLight,
         double gLight, int twoJHeavy, int parityHeavy, int zHeavy, double massHeavy, double gHeavy,
         double lightExcitation, double heavyExcitation, double separationEnergy, double radius, int externalCapture,
         double bindingEnergy) {
        ex::NewPair n;
        n.key = key;
        n.kind = kind;
        n.twoJLight = twoJLight;
        n.parityLight = parityLight;
        n.zLight = zLight;
        n.massLight = massLight;
        n.gLight = gLight;
        n.twoJHeavy = twoJHeavy;
        n.parityHeavy = parityHeavy;
        n.zHeavy = zHeavy;
        n.massHeavy = massHeavy;
        n.gHeavy = gHeavy;
        n.lightExcitation = lightExcitation;
        n.heavyExcitation = heavyExcitation;
        n.separationEnergy = separationEnergy;
        n.radius = radius;
        n.externalCapture = externalCapture;
        n.bindingEnergy = bindingEnergy;
        int used = 0;
        const ex::LevelEdit e = ex::AddPair(doc, n, used);
        return py::make_tuple(Done(e), used);
      },
      py::arg("doc"), py::kw_only(), py::arg("key"), py::arg("kind"), py::arg("two_j_light"), py::arg("parity_light"),
      py::arg("z_light"), py::arg("mass_light"), py::arg("g_light"), py::arg("two_j_heavy"), py::arg("parity_heavy"),
      py::arg("z_heavy"), py::arg("mass_heavy"), py::arg("g_heavy"), py::arg("light_excitation"),
      py::arg("heavy_excitation"), py::arg("separation_energy"), py::arg("radius"), py::arg("external_capture"),
      py::arg("binding_energy"),
      "(replacements, key): the pair row of a new pair (experiment::AddPair: inactive, level ID 0) at the end of "
      "<levels>; key 0 takes the next free key.  Masses in u, energies in MeV, the radius in fm; kind 0 particle, 10 "
      "photon, 20 beta.");
  m.def(
      "remove_pair_edits",
      [](const io::Document &doc, int pair, bool removeChannels) { return Done(ex::RemovePair(doc, pair, removeChannels)); },
      py::arg("doc"), py::arg("pair"), py::kw_only(), py::arg("remove_channels"),
      "The replacements that remove pair `pair` (experiment::RemovePair): its pair row and, with remove_channels, the "
      "channel rows of the levels that name it; ValueError while a level names it without remove_channels, or when a "
      "level would be left without a channel.");
  m.def(
      "allowed_channels",
      [](const io::Document &doc, int twoJ, int parity, int lMax, int multipolarityMax, int multipolaritiesMax,
         const std::vector<int> &pairs) {
        ex::ChannelLimits limits;
        limits.lMax = lMax;
        limits.multipolarityMax = multipolarityMax;
        limits.multipolaritiesMax = multipolaritiesMax;
        const ex::AllowedChannelList a = ex::AllowedChannels(doc, twoJ, parity, limits, pairs);
        if (!a.error.empty()) ThrowRefusal(a.refusal);
        py::list out;
        for (const ex::LevelChannelInfo &c : a.channels) out.append(py::make_tuple(c.pairKey, c.l, c.twoS));
        return out;
      },
      py::arg("doc"), py::arg("two_j"), py::arg("parity"), py::kw_only(), py::arg("l_max"),
      py::arg("multipolarity_max"), py::arg("multipolarities_max"), py::arg("pairs"),
      "[(pair key, l, 2s)]: the channels a level of 2J and parity has with the pairs of <levels> (pairs: keys, or [] "
      "for all), by the rule of experiment::AllowedChannels (the Qt LevelsTab's).");
  m.def(
      "remove_level_edits", [](const io::Document &doc, const py::object &level) { return Done(ex::RemoveLevel(doc, Address(level))); },
      py::arg("doc"), py::arg("level"), "The replacements that remove every row of a level.");
  m.def(
      "level_value_edits",
      [](const io::Document &doc, const std::vector<std::tuple<std::string, int, int, int, int, int, double>> &values) {
        std::vector<ex::LevelValue> v;
        for (const auto &x : values) {
          ex::LevelValue lv;
          const std::string &kind = std::get<0>(x);
          if (kind != "energy" && kind != "width") throw py::value_error("a level value is an 'energy' or a 'width'");
          lv.kind = kind == "energy" ? ex::LevelValue::Kind::Energy : ex::LevelValue::Kind::Width;
          lv.group = std::get<1>(x);
          lv.level = std::get<2>(x);
          lv.channel = ex::LevelChannelInfo{std::get<3>(x), std::get<4>(x), std::get<5>(x)};
          lv.value = std::get<6>(x);
          v.push_back(lv);
        }
        const ex::LevelEdit e = ex::LevelValueEdits(doc, v);
        return py::make_tuple(Done(e), e.unplaced);
      },
      py::arg("doc"), py::arg("values"),
      "The replacements that write fitted level values [(kind 'energy' | 'width', group, level, pair, l, 2s, value)] "
      "into the rows they come from (experiment::LevelValueEdits), and the values with no row: (replacements, "
      "[(index into values, why)]).");
  m.def(
      "merge_level_rows",
      [](const std::vector<const io::Document *> &members) {
        const ex::MergedLevelRows r = ex::MergeLevelRows(members);
        if (!r.error.empty()) throw py::value_error(r.error);
        return py::make_tuple(r.rows, r.pairKeys);
      },
      py::arg("members"),
      "The <levels> rows of several projects (documents) in one section (experiment::MergeLevelRows): (rows, per "
      "member {its pair key: the merged key}); each member's pair keys moved past the largest keys of the members "
      "before it, its level IDs past the largest ID before it.  ValueError when a member's section does not read.");

  m.def(
      "config_entries",
      [](const io::Document &doc) {
        std::size_t open = 0;
        const io::ConfigSection c = io::ReadConfigSection(ConfigLines(doc, open));
        py::list entries;
        for (const io::ConfigEntry &e : c.entries) {
          py::dict d;
          d["name"] = e.name;
          d["value"] = py::bytes(e.value);
          d["line"] = (int)open + e.line + 1;
          d["positional"] = e.positional;
          entries.append(d);
        }
        return py::make_tuple(entries, c.closed);
      },
      py::arg("doc"),
      "The <config> entries as the engine reads them (io::ReadConfigSection): ([dicts with name, value, line "
      "(1-based), positional], the section is complete).");
  m.def(
      "config_keys",
      []() {
        static const char *const kinds[] = {"bool", "int", "number", "path", "version"};
        py::list out;
        for (const io::ConfigKey &k : io::ConfigKeys()) out.append(py::make_tuple(k.name, kinds[k.kind], k.doc));
        return py::make_tuple(io::ConfigPositionalNames(), out);
      },
      "(the positional entries' names, [(key, kind, doc)]) of <config>.");
  m.def(
      "config_table",
      []() {
        static const char *const kinds[] = {"bool", "int", "number", "path", "version"};
        py::list out;
        for (const io::ConfigPositional &p : io::ConfigPositionals()) {
          py::dict d;
          d["name"] = p.name;
          d["positional"] = true;
          d["kind"] = p.token ? "choice" : "line";
          d["values"] = p.values;
          d["title"] = p.title;
          d["group"] = p.group;
          d["doc"] = p.doc;
          d["rare"] = p.rare;
          out.append(d);
        }
        for (const io::ConfigKey &k : io::ConfigKeys()) {
          py::dict d;
          d["name"] = k.name;
          d["positional"] = false;
          d["kind"] = kinds[k.kind];
          d["minimum"] = k.minimum;
          d["title"] = k.title;
          d["group"] = k.group;
          d["doc"] = k.doc;
          d["absent"] = k.absent;
          d["rare"] = k.rare;
          d["option"] = k.option;
          if (k.option) {
            d["default"] = k.defaultOn;
            d["flag"] = k.flag;
            d["negation"] = k.negation;
            d["effect"] = k.effect;
          }
          out.append(d);
        }
        return out;
      },
      "Every entry of <config> as the engine's tables hold it (io::ConfigPositionals(), io::ConfigKeys()), in file "
      "order: name, positional, kind (choice | line for a positional entry; bool | int | number | path | version for "
      "a key), values (a choice's, blank-separated), minimum (an int key's), title, group, doc, absent (what a missing "
      "key means), rare, option (and default, flag, negation, effect of an option).");
  m.def(
      "config_edits",
      [](const io::Document &doc, const std::string &name, const std::string &value, bool remove) {
        std::size_t open = 0;
        const std::vector<std::string> lines = ConfigLines(doc, open);
        std::vector<io::ConfigLineEdit> edits;
        std::string why;
        if (!io::ConfigEdits(lines, name, value, remove ? io::ConfigChange::Remove : io::ConfigChange::Set, edits, why))
          throw py::value_error(why);
        std::vector<ex::LineReplacement> out;
        for (const io::ConfigLineEdit &e : edits) {
          const int line = (int)open + e.line;
          if (e.insert)
            out.push_back({line, line, {e.text}});
          else
            out.push_back({line, line + 1, e.remove ? std::vector<std::string>() : std::vector<std::string>{e.text}});
        }
        std::sort(out.begin(), out.end(), [](const ex::LineReplacement &a, const ex::LineReplacement &b) { return a.first > b.first; });
        return Replacements(out);
      },
      py::arg("doc"), py::arg("name"), py::arg("value"), py::kw_only(), py::arg("remove") = false,
      "The replacements that set a <config> entry (io::ConfigEdits) or, with remove, remove a key.");

  m.def(
      "thm_keys",
      []() {
        py::dict out;
        out[py::str(ThmOptionKind().name)] = ThmKeyKindDict(ThmOptionKind(), ThmOptionKeys());
        out[py::str(ThmExperimentKind().name)] = ThmKeyKindDict(ThmExperimentKind(), ThmExperimentKeys());
        return out;
      },
      "The keys of <thm> as two record kinds of a schema, '<thm>' (the global options) and 'experiment' (the keys of "
      "an experiment[<name>] line): {kind: {title, doc, fields: {key: {type, doc, required, title, unit, values, "
      "default, grammar, applies_when}}}}, from the tables the engine's parser reads (ThmOptionKeys, "
      "ThmExperimentKeys).");
  m.def(
      "thm_entries",
      [](const io::Document &doc) {
        std::size_t open = 0;
        const std::vector<std::string> body = ThmBody(doc, open);
        py::list out;
        for (const ThmEntry &e : ThmBlockEntries(body))
          out.append(py::make_tuple(e.experiment, e.key, e.value, open == std::string::npos ? 0 : (int)open + 2 + e.line));
        return py::make_tuple(open != std::string::npos, out);
      },
      py::arg("doc"),
      "The <thm> section's key=value entries (ThmBlockEntries): (the section exists, [(experiment or '', key, value, "
      "line (1-based))]).");
  m.def(
      "thm_option_edits",
      [](const io::Document &doc, const std::string &key, const std::string &value, bool remove) {
        std::size_t open = 0;
        const std::vector<std::string> body = ThmBody(doc, open);
        std::vector<ThmLineEdit> edits;
        std::string why;
        if (!ThmOptionEdits(body, key, value, remove ? ThmChange::Remove : ThmChange::Set, edits, why)) throw py::value_error(why);
        return ThmReplacements(doc, open, body, edits);
      },
      py::arg("doc"), py::arg("key"), py::arg("value"), py::kw_only(), py::arg("remove") = false,
      "The replacements that set (or remove) a global option of <thm> (ThmOptionEdits); a project without <thm> gets "
      "the section after </config>.");
  m.def(
      "thm_experiment_edits",
      [](const io::Document &doc, const std::string &name, const std::string &key, const std::string &value, bool remove) {
        std::size_t open = 0;
        const std::vector<std::string> body = ThmBody(doc, open);
        std::vector<ThmLineEdit> edits;
        std::string why;
        if (!ThmExperimentEdits(body, name, key, value, remove ? ThmChange::Remove : ThmChange::Set, edits, why))
          throw py::value_error(why);
        return ThmReplacements(doc, open, body, edits);
      },
      py::arg("doc"), py::arg("name"), py::arg("key"), py::arg("value"), py::kw_only(), py::arg("remove") = false,
      "The replacements that set (or remove) a key of experiment[name] (ThmExperimentEdits); with remove and key '' "
      "the experiment's lines.");
  m.def(
      "check_format3_thm",
      [](const std::string &text) {
        ex::Format3ReadOptions read;
        read.readDataFiles = false;
        read.setUpChecks = false;
        const ex::Format3Project p = ex::ReadFormat3Text(text, read);
        const io::Document doc = io::Document::Parse(text);
        const io::Document::Section *s = doc.Find("thm");
        if (!p.Ok()) return std::string();  // the project's own errors say what is wrong
        if (!s) return p.built.thm.experiments.empty() ? std::string() : "the project has no <thm> section, which its THM data name";  // nothing to check, or the project's own errors say what is wrong
        std::vector<std::string> section{"<thm>"};
        for (const std::string &l : doc.Body(*s)) section.push_back(l);
        section.push_back("</thm>");
        ThmOptions thm;
        std::vector<std::string> idle, warnings;
        return ReadFormat3ThmSection(section, p.built.thm, "", false, thm, idle, warnings);
      },
      py::arg("text"),
      "The <thm> section of a format-3 project text as the engine reads it (ReadFormat3ThmSection: members and "
      "weights from the records; files not read): \"\" or the engine's message.");
}
