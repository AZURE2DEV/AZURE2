// The THM tables and grammars of the engine for Python without a session: a cbackground= value as the engine reads it
// (pyazr.snapshot), the built-in nuclei and global optical potentials.  The <thm> section of a project is checked by
// check_format3_thm (pybind_nuclear_part.cpp, ReadFormat3ThmSection).  No files are read.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "ThmExperiment.h"
#include "azure2/physics/transfer/GlobalOptical.h"

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

}  // namespace

void BindThm(py::module_ &m) {
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
        for (const azure2::physics::transfer::GlobalOpticalSet &g : azure2::physics::transfer::GlobalOpticalSets()) {
          py::dict d;
          d["name"] = g.name;
          d["projectiles"] = g.projectiles;
          d["reference"] = g.reference;
          d["a_range"] = py::make_tuple(g.a_min, g.a_max);
          d["e_range"] = py::make_tuple(g.e_min_MeV, g.e_max_MeV);
          out.append(d);
        }
        return out;
      },
      "The built-in global optical potentials: name, projectiles, reference, target A and lab-energy (MeV) ranges.");
}
