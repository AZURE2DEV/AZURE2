// The engine's reading of the <thm> block, for Python without a session (written for pyazr's old file class, deleted
// in pyazr 3; parse_thm_cbackground and thm_nuclides still have callers): the <thm> lines are read and checked here,
// so pyazr and AZURE2 take exactly the same blocks (docs/upgrade/project.md section 5).  No files are read except by
// read_thm_weight_table.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "ThmExperiment.h"
#include "ThmOptical.h"
#include "ThmOptions.h"

namespace py = pybind11;

namespace {

py::tuple NuclideTuple(const ThmNuclide &n) { return py::make_tuple(n.name, n.Z, n.A, n.mass); }

py::dict CoherentTermDict(const ThmExperiment::CoherentTerm &t) {
  py::dict d;
  d["J"] = t.J;
  d["parity"] = t.parity;
  d["exit"] = t.exitKey;
  d["channels"] = t.hasChannels ? py::object(py::make_tuple(t.s, t.l, t.sp, t.lp)) : py::object(py::none());
  d["form"] = t.form == 2 ? "linear" : "const";
  if (t.hasValues) {
    py::list values;
    for (int k = 0; k < 2 * t.form; k++) values.append(py::make_tuple(t.value[k], t.fixed[k]));
    d["values"] = values;
  } else {
    d["values"] = py::none();
  }
  return d;
}

// One experiment as read_thm_block returns it: the keys given (in order) and, per key, the value -- the text as written
// for the keys whose value is a word or a specification, the number for Ebeam / psNodes / spectatorAngleNodes, the
// sorted segment list, (text, Z, A, mass) for the three nuclei, on/off as a bool, and cbackground canonical.
py::dict ExperimentDict(const ThmExperiment &x) {
  py::dict d;
  py::list keys;
  for (const std::string &key : x.keys) keys.append(key);
  d["keys"] = keys;
  for (size_t i = 0; i < x.keys.size(); i++) {
    const std::string &key = x.keys[i], &text = x.values[i];
    if (key == "segments")
      d["segments"] = x.segments;
    else if (key == "background")
      d["background"] = ThmExperiment::BackgroundName(x.backgroundTerms);
    else if (key == "beam")
      d["beam"] = NuclideTuple(x.beam);
    else if (key == "target")
      d["target"] = NuclideTuple(x.target);
    else if (key == "spectator")
      d["spectator"] = NuclideTuple(x.spectator);
    else if (key == "Ebeam")
      d["Ebeam"] = x.beamEnergy;
    else if (key == "lineshape")
      d["lineshape"] = x.lineshape;
    else if (key == "psNodes")
      d["psNodes"] = x.psNodes;
    else if (key == "spectatorAngleNodes")
      d["spectatorAngleNodes"] = x.angleNodes;
    else if (key == "cbackground")
      d["cbackground"] = FormatThmCoherentBackground(x.cbackground);
    else
      d[py::str(key)] = text;
  }
  return d;
}

// The block as read_thm_block returns it: the five global options (later lines win), the per-pair spectator energies,
// the weight files, and the experiments by name in order of first appearance.  ValueError with AZURE2's message.
py::dict ReadThmBlock(const std::vector<std::string> &lines) {
  std::vector<ThmOptionLine> options;
  std::vector<ThmExperiment> experiments;
  const std::string why = ParseThmBlockLines(lines, options, experiments);
  if (!why.empty()) throw py::value_error(why);
  py::dict s, byPair, weight, weightTest, named;
  s["entranceL"] = "incoherent";
  s["vertex"] = "constant";
  s["kinematics"] = "lacognata";
  s["coulombIntegral"] = false;
  s["spectatorEnergy"] = 0.0;
  for (const ThmOptionLine &o : options) {
    switch (o.kind) {
      case ThmOptionLine::VERTEX: s["vertex"] = o.word; break;
      case ThmOptionLine::KINEMATICS: s["kinematics"] = o.word; break;
      case ThmOptionLine::ENTRANCE_L: s["entranceL"] = o.word; break;
      case ThmOptionLine::COULOMB_INTEGRAL: s["coulombIntegral"] = o.flag; break;
      case ThmOptionLine::SPECTATOR_ENERGY:
        if (o.perPair)
          byPair[py::int_(o.pairKey)] = o.energy;
        else
          s["spectatorEnergy"] = o.energy;
        break;
      case ThmOptionLine::WEIGHT: (o.test ? weightTest : weight)[py::int_(o.segment)] = o.file; break;
    }
  }
  for (const ThmExperiment &x : experiments) named[py::str(x.name)] = ExperimentDict(x);
  s["spectatorByPair"] = byPair;
  s["weight"] = weight;
  s["weightTest"] = weightTest;
  s["experiments"] = named;
  return s;
}

// One global line key=value: (kind, word, flag, energy, pair key or None, test, segment, file), or None if AZURE2 does not
// take it.
py::object ThmOption(const std::string &key, const std::string &value) {
  ThmOptionLine o;
  if (!ParseThmOptionLine(key, value, o)) return py::none();
  static const char *kinds[] = {"vertex", "kinematics", "entranceL", "coulombIntegral", "spectatorEnergy", "weight"};
  py::dict d;
  d["kind"] = kinds[o.kind];
  d["word"] = o.word;
  d["flag"] = o.flag;
  d["energy"] = o.energy;
  d["pair"] = o.perPair ? py::object(py::int_(o.pairKey)) : py::object(py::none());
  d["test"] = o.test;
  d["segment"] = o.segment;
  d["file"] = o.file;
  return d;
}

}  // namespace

void BindThm(py::module_ &m) {
  m.def("check_thm_block", &CheckThmBlockLines, py::arg("lines"), py::arg("global_rules") = false,
        "The lines of a <thm> block as AZURE2 reads them: \"\" or the message it reports.");
  m.def("read_thm_block", &ReadThmBlock, py::arg("lines"),
        "The lines of a <thm> block read as AZURE2 reads them (options and experiments, see pybind_thm.cpp); "
        "ValueError with AZURE2's message.");
  m.def("thm_option", &ThmOption, py::arg("key"), py::arg("value"),
        "One global <thm> key=value as AZURE2 reads it (a dict), or None if it is refused.");
  m.def(
      "parse_thm_cbackground",
      [](const std::string &value) {
        std::vector<ThmExperiment::CoherentTerm> terms;
        const std::string why = ParseThmCoherentBackground(value, terms);
        if (!why.empty()) throw py::value_error(why);
        py::list out;
        for (const auto &t : terms) out.append(CoherentTermDict(t));
        return py::make_tuple(out, FormatThmCoherentBackground(terms));
      },
      py::arg("value"), "A cbackground= value: (terms as dicts, the canonical value); ValueError if refused.");
  m.def("thm_spin_text", &ThmSpinText, py::arg("j"), "A spin as AZURE2 writes it in parameter names: 1/2, 2.");
  m.def(
      "thm_nuclides",
      []() {
        py::list out;
        for (const ThmNuclide &n : ThmNuclide::Table()) out.append(NuclideTuple(n));
        return out;
      },
      "The built-in nuclei of experiment lines: (name, Z, A, nuclear mass in u).");
  m.def(
      "thm_global_opticals",
      []() {
        py::list out;
        for (const ThmGlobalOptical &g : ThmGlobalOpticals()) {
          py::dict d;
          d["name"] = g.name;
          d["projectiles"] = g.projectiles;
          d["reference"] = g.reference;
          d["a_range"] = py::make_tuple(g.aMin, g.aMax);
          d["e_range"] = py::make_tuple(g.eMin, g.eMax);
          out.append(d);
        }
        return out;
      },
      "The built-in global optical potentials: name, projectiles, reference, target A and lab-energy (MeV) ranges.");
  m.def(
      "thm_global_optical_for",
      [](const std::string &name, int Z, int A) {
        const int k = ThmGlobalOpticalIndex(name);
        return k >= 0 && ThmGlobalOpticalFor(k, Z, A);
      },
      py::arg("name"), py::arg("Z"), py::arg("A"), "Whether the global optical potential is made for (Z, A).");
  m.def(
      "read_thm_weight_table",
      [](const std::string &path) {
        ThmWeightTable table;
        return table.Read(path);
      },
      py::arg("path"), "\"\" if the weight file is usable as AZURE2 reads it, else why not.");
}
