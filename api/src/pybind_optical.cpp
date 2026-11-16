// Elastic scattering on an optical potential for pyazr (pyazr.physics.optical_scattering): the distorted waves of
// physics/transfer (docs/upgrade/transfer.md sections 2 and 4) from plain numbers, without a session or a project.
// No physics in this file: it fills an OpticalPotential and a ScatteringRequest from the caller's numbers, calls
// SolveDistortedWaves, ElasticOverRutherford / ElasticCrossSection and ReactionCrossSection, and hands back arrays.
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/physics/transfer/DistortedWaves.h"
#include "azure2/physics/transfer/OpticalPotential.h"

namespace py = pybind11;
namespace transfer = azure2::physics::transfer;

namespace {

using RealArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

/// One term of the potential: (depth in MeV, radius parameter in fm, diffuseness in fm).
transfer::WoodsSaxonTerm TermOf(const py::handle &value, const std::string &name) {
  const auto numbers = py::cast<std::vector<double>>(value);
  if (numbers.size() != 3)
    throw std::invalid_argument("optical_scattering: the term '" + name + "' of the potential must be (depth_MeV, "
                                "r0_fm, diffuseness_fm), three numbers");
  for (double x : numbers)
    if (!std::isfinite(x))
      throw std::invalid_argument("optical_scattering: a number of the term '" + name + "' is not finite");
  return transfer::WoodsSaxonTerm{numbers[0], numbers[1], numbers[2]};
}

/// The potential of OpticalPotential.h's Eq. (2.1) from a dict of its members; an unknown key is refused.
transfer::OpticalPotential PotentialOf(const py::dict &d) {
  transfer::OpticalPotential p;
  for (const auto &item : d) {
    const std::string key = py::cast<std::string>(item.first);
    if (key == "radius_convention") {
      const std::string c = py::cast<std::string>(item.second);
      if (c == "absolute")
        p.convention = transfer::RadiusConvention::Absolute;
      else if (c == "target_mass")
        p.convention = transfer::RadiusConvention::TargetMass;
      else if (c == "both_masses")
        p.convention = transfer::RadiusConvention::BothMasses;
      else
        throw std::invalid_argument("optical_scattering: radius_convention '" + c + "' is not known; allowed: "
                                    "'absolute' (R = r0), 'target_mass' (R = r0 A_t^1/3), 'both_masses' "
                                    "(R = r0 (A_p^1/3 + A_t^1/3))");
    } else if (key == "target_mass_number") {
      p.target_mass_number = py::cast<double>(item.second);
    } else if (key == "projectile_mass_number") {
      p.projectile_mass_number = py::cast<double>(item.second);
    } else if (key == "coulomb_r0_fm") {
      p.coulomb_r0_fm = py::cast<double>(item.second);
    } else if (key == "real_volume") {
      p.real_volume = TermOf(item.second, key);
    } else if (key == "imaginary_volume") {
      p.imaginary_volume = TermOf(item.second, key);
    } else if (key == "imaginary_surface") {
      p.imaginary_surface = TermOf(item.second, key);
    } else if (key == "real_spin_orbit") {
      p.real_spin_orbit = TermOf(item.second, key);
    } else if (key == "imaginary_spin_orbit") {
      p.imaginary_spin_orbit = TermOf(item.second, key);
    } else {
      throw std::invalid_argument("optical_scattering: the potential has no member '" + key + "'; allowed: "
                                  "radius_convention, target_mass_number, projectile_mass_number, coulomb_r0_fm, "
                                  "real_volume, imaginary_volume, imaginary_surface, real_spin_orbit, "
                                  "imaginary_spin_orbit");
    }
  }
  return p;
}

py::dict OpticalScattering(int z1, int z2, double m1_u, double m2_u, double e_cm, const py::dict &potential,
                           const RealArray &angles_deg, int two_s, int l_max, double step_fm,
                           double matching_radius_fm) {
  if (!(m1_u > 0.0) || !(m2_u > 0.0) || !std::isfinite(m1_u) || !std::isfinite(m2_u))
    throw std::invalid_argument("optical_scattering: m1_u and m2_u must be finite masses > 0 (u)");
  if (angles_deg.ndim() != 1)
    throw std::invalid_argument("optical_scattering: angles_deg must be one-dimensional");
  for (py::ssize_t i = 0; i < angles_deg.shape(0); i++) {
    const double theta = angles_deg.data()[i];
    if (!(theta > 0.0 && theta <= 180.0))
      throw std::invalid_argument("optical_scattering: angle " + std::to_string(i) + " = " + std::to_string(theta) +
                                  " deg is outside 0 < theta <= 180 (c.m.)");
  }
  const transfer::OpticalPotential p = PotentialOf(potential);
  transfer::ScatteringRequest q;
  q.pair.z1 = z1;
  q.pair.z2 = z2;
  q.pair.mu_u = m1_u * m2_u / (m1_u + m2_u);  // as the rows of pybind_physics.cpp make a pair
  q.two_s = two_s;
  q.energy_MeV = e_cm;
  q.step_fm = step_fm;
  q.matching_radius_fm = matching_radius_fm;
  q.l_max = l_max;
  q.keep_waves = false;
  const transfer::DistortedWaves w = transfer::SolveDistortedWaves(p, q);
  const py::ssize_t n = static_cast<py::ssize_t>(w.waves.size());
  py::array_t<int> l(n), twoJ(n);
  py::array_t<std::complex<double>> S(n);
  for (py::ssize_t i = 0; i < n; i++) {
    l.mutable_data()[i] = w.waves[i].l;
    twoJ.mutable_data()[i] = w.waves[i].two_j;
    S.mutable_data()[i] = w.waves[i].S;
  }
  const bool charged = z1 * z2 != 0;
  py::array_t<double> elastic(angles_deg.shape(0));
  for (py::ssize_t i = 0; i < angles_deg.shape(0); i++)
    elastic.mutable_data()[i] = charged ? transfer::ElasticOverRutherford(w, angles_deg.data()[i])
                                        : transfer::ElasticCrossSection(w, angles_deg.data()[i]);
  py::array_t<double> sigma(static_cast<py::ssize_t>(w.sigma.size()));
  for (std::size_t i = 0; i < w.sigma.size(); i++) sigma.mutable_data()[i] = w.sigma[i];
  py::dict out;
  out["k"] = w.k_fm;
  out["eta"] = w.eta;
  out["coulomb_phase"] = sigma;
  out["l"] = l;
  out["two_j"] = twoJ;
  out["S"] = S;
  out[charged ? "ratio_to_rutherford" : "dsigma_domega_mb_sr"] = elastic;
  out["reaction_mb"] = transfer::ReactionCrossSection(w);
  out["spin_orbit"] = w.spin_orbit;
  out["description"] = transfer::Describe(p);
  return out;
}

}  // namespace

void BindOptical(py::module_ &m) {
  m.def("optical_scattering", &OpticalScattering, py::arg("z1"), py::arg("z2"), py::arg("m1_u"), py::arg("m2_u"),
        py::arg("e_cm"), py::arg("potential"), py::arg("angles_deg"), py::arg("two_s") = 0, py::arg("l_max") = 20,
        py::arg("step_fm") = 0.02, py::arg("matching_radius_fm") = 20.0,
        "Elastic scattering of a projectile of spin two_s / 2 (0, 1/2 or 1) on a target whose spin is ignored, in "
        "the optical potential of physics/transfer/OpticalPotential.h Eq. (2.1) (a dict of its members; a term is "
        "(depth_MeV, r0_fm, diffuseness_fm)), at the c.m. energy e_cm (MeV), masses in u: a dict of k (fm^-1), eta, "
        "coulomb_phase (sigma_l, rad, l = 0 ... l_max), l, two_j and S (the nuclear S-matrix element of each "
        "partial wave), ratio_to_rutherford at angles_deg (c.m.; a neutral pair: dsigma_domega_mb_sr), reaction_mb, "
        "spin_orbit and description (the potential in words).  pyazr.physics.optical_scattering is the documented "
        "call.");
}
