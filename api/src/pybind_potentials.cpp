// The hybrid potentials of a pyazr session (Session, pybind_api.cpp): the setting set_nuclear_potential makes, the view
// nuclear_potential reads, and the session's potentials as `potential` records of format 3 (pyazr's save_fit of a
// format-3 project), written by the converter's writer (convert/PotentialRecords.h), so both write the same records.
#include <stdexcept>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Config.h"
#include "azure2/convert/PotentialRecords.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/RecordReader.h"

namespace py = pybind11;
using azure2::experiment::HybridPotentialSetting;

HybridPotentialSetting SessionPotentialSetting(const std::string &type, bool enabled, double V0, double R, double a,
                                               double r0, double coulombRadius_fm) {
  HybridPotentialSetting s;
  s.enabled = enabled;
  if (type == "WoodsSaxon")
    s.shape = azure2::physics::FormShape::WoodsSaxon;
  else if (type == "Gaussian")
    s.shape = azure2::physics::FormShape::Gaussian;
  else
    throw std::invalid_argument("unknown potential type '" + type + "'; expected 'WoodsSaxon' or 'Gaussian'");
  s.V0 = V0;
  s.R = R;
  s.a = a;
  s.r0 = r0;
  s.coulomb_radius_fm = coulombRadius_fm;
  s.Potential();  // throws std::invalid_argument when a parameter is not positive
  return s;
}

py::tuple SessionPotentialView(const Config &config, int pair) {
  const azure2::experiment::HybridPotentialSettings &p = config.potentials;
  const HybridPotentialSetting &s = pair ? p.For(pair) : p.defaults;
  const bool own = pair ? p.pairs.count(pair) > 0 : true;
  return py::make_tuple(s.enabled, s.ShapeName(), s.V0, s.R, s.a, s.r0, own, s.coulomb_radius_fm);
}

py::tuple SessionPotentialRecords(const Config &config) {
  if (!config.format3)
    throw std::invalid_argument("the session's project is not in format 3: its potentials are the <potential> block");
  const azure2::experiment::Format3Project &f = *config.format3;
  // changed: a particle pair's potential differs from the one the project's records give
  bool changed = false;
  for (const azure2::experiment::PairInfo &pair : f.levels.pairs)
    if (pair.kind == 0 && config.potentials.PotentialFor(pair.key) != f.built.potentials.PotentialFor(pair.key)) changed = true;
  // the other records' names are taken; a pair keeps the name of its record
  std::set<std::string> taken;
  for (const azure2::experiment::ValidatedRecord &r : f.validated.records)
    if (r.kind != "potential") taken.insert(r.name);
  std::map<int, std::string> names;
  for (const auto &r : f.built.potentialRecords) names[r.first] = r.second.first;
  const azure2::convert::PotentialRecordsResult out =
      azure2::convert::PotentialRecords(config.potentials, f.levels.pairs, taken, names);
  return py::make_tuple(changed, azure2::io::WriteRecords(out.records), out.notes);
}
