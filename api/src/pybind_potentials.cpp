// The hybrid potentials of a pyazr session (Session, pybind_api.cpp): the view pyazr's nuclear_potential reads, each
// pair's setting as the session's Config resolved it from the project's potential records.
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Config.h"

namespace py = pybind11;
using azure2::experiment::HybridPotentialSetting;

py::tuple SessionPotentialView(const Config &config, int pair, int l) {
  const azure2::experiment::HybridPotentialSettings &p = config.potentials;
  // The record of this l when the pair has one (a record with an l list), else the pair's record or the default.
  const auto ofL = (pair && l >= 0) ? p.byL.find({pair, l}) : p.byL.end();
  const HybridPotentialSetting &s = ofL != p.byL.end() ? ofL->second : pair ? p.For(pair) : p.defaults;
  const bool own = ofL != p.byL.end() || (pair ? p.pairs.count(pair) > 0 : true);
  std::vector<double> spinOrbit;
  if (s.spin_orbit)
    spinOrbit = {s.spin_orbit->depth_MeV, s.spin_orbit->form.radius_fm, s.spin_orbit->form.diffuseness_fm,
                 s.spin_orbit->acts_on == azure2::physics::SpinOrbitParticle::Heavy ? 1.0 : 0.0};
  return py::make_tuple(s.enabled, s.ShapeName(), s.V0, s.R, s.a, s.r0, own, s.coulomb_radius_fm, spinOrbit,
                        s.matching_radius_fm, s.step_fm,
                        s.interior_boundary ? py::object(py::float_(*s.interior_boundary)) : py::object(py::none()));
}
