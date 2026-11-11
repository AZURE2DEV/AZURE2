// The hybrid potentials of a pyazr session (Session, pybind_api.cpp): the setting set_nuclear_potential makes, the view
// nuclear_potential reads, and the session's potentials as `potential` records of format 3 (pyazr's save_fit of a
// format-3 project), written by the converter's writer (convert/PotentialRecords.h), so both write the same records.
#include <stdexcept>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "CNuc.h"
#include "Config.h"
#include "PPair.h"
#include "azure2/experiment/SpinOrbitChecks.h"
#include "azure2/convert/PotentialRecords.h"
#include "azure2/experiment/ProjectReader.h"
#include "azure2/io/RecordReader.h"

namespace py = pybind11;
using azure2::experiment::HybridPotentialSetting;

// spinOrbit: {} (none) or {V_so (MeV), R_so (fm), a_so (fm), on (0 light, 1 heavy)}; refused on the default (pair 0,
// which would give the term to every pair) and where SpinOrbitPairRefusal refuses it for the pair of `compound`.
HybridPotentialSetting SessionPotentialSetting(const std::string &type, bool enabled, double V0, double R, double a,
                                               double r0, double coulombRadius_fm, const std::vector<double> &spinOrbit,
                                               CNuc *compound, int pair) {
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
  if (!spinOrbit.empty()) {
    if (spinOrbit.size() != 4) throw std::invalid_argument("spin_orbit needs [V_so, R_so, a_so, on]");
    if (!pair) throw std::invalid_argument("a spin-orbit term belongs to one pair: name the pair (the default would give it to every pair)");
    const azure2::physics::SpinOrbitParticle on =
        spinOrbit[3] != 0.0 ? azure2::physics::SpinOrbitParticle::Heavy : azure2::physics::SpinOrbitParticle::Light;
    s.spin_orbit = azure2::physics::WoodsSaxonSpinOrbit(spinOrbit[0], spinOrbit[1], spinOrbit[2], on);
    if (compound && compound->IsPairKey(pair)) {
      const PPair *p = compound->GetPair(compound->GetPairNumFromKey(pair));
      const azure2::experiment::SpinOrbitRefusal r = azure2::experiment::SpinOrbitPairRefusal(
          (int)std::lround(2.0 * p->GetJ(1)), (int)std::lround(2.0 * p->GetJ(2)), on, std::to_string(pair));
      if (r) throw std::invalid_argument("ERROR (" + r.code + "): " + r.message);
    }
  }
  s.Potential();  // throws std::invalid_argument when a parameter is not positive
  return s;
}

py::tuple SessionPotentialView(const Config &config, int pair) {
  const azure2::experiment::HybridPotentialSettings &p = config.potentials;
  const HybridPotentialSetting &s = pair ? p.For(pair) : p.defaults;
  const bool own = pair ? p.pairs.count(pair) > 0 : true;
  std::vector<double> spinOrbit;
  if (s.spin_orbit)
    spinOrbit = {s.spin_orbit->depth_MeV, s.spin_orbit->form.radius_fm, s.spin_orbit->form.diffuseness_fm,
                 s.spin_orbit->acts_on == azure2::physics::SpinOrbitParticle::Heavy ? 1.0 : 0.0};
  return py::make_tuple(s.enabled, s.ShapeName(), s.V0, s.R, s.a, s.r0, own, s.coulomb_radius_fm, spinOrbit);
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
