// The stopping/ library in the _azure2 module (submodule _azure2.stopping; pyazr.stopping wraps it): SRIM-2013's
// stopping powers ported from LunaSRIM, energy loss, range and straggling of a layer, provenance
// (docs/upgrade/stopping.md).  A thin binding: arrays in, arrays out, the C++ functions element-wise; no second
// implementation.  Units as stopping/'s Material.h: energies MeV (lab, of the projectile), stopping cross sections
// eV/(1e15 atoms/cm^2) per atom of the material, areal densities atoms/cm^2 of the material, straggling widths MeV.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/EnergyLoss.h"
#include "azure2/stopping/Stopping.h"
#include "azure2/stopping/Straggling.h"

namespace py = pybind11;
namespace st = azure2::stopping;

namespace {

using Array = py::array_t<double, py::array::forcecast>;

/// The material of a formula, its phase ("" = DefaultPhase, "solid", "gas") and SRIM's compound correction.
st::Material MakeMaterial(const std::string &formula, const std::string &phase, double compoundCorrection) {
  st::Material m = st::MaterialFromFormula(formula);
  if (!phase.empty()) m.phase = st::PhaseFromName(phase);
  m.compound_correction = compoundCorrection;
  return m;
}

/// The projectile of charge z1 and atomic mass mass_u (u).
st::Projectile MakeProjectile(int z1, double massU) {
  if (!(massU > 0.0)) throw std::invalid_argument("stopping: the projectile's mass must be > 0 u");
  return st::Ion(z1, massU, "Z1=" + std::to_string(z1));
}

std::vector<double> Values(const Array &a) {
  const auto v = a.unchecked<1>();
  std::vector<double> out((size_t)v.shape(0));
  for (py::ssize_t i = 0; i < v.shape(0); i++) out[(size_t)i] = v(i);
  return out;
}

void SameLength(const std::vector<double> &a, const std::vector<double> &b, const char *what) {
  if (a.size() != b.size())
    throw std::invalid_argument(std::string("stopping: ") + what + " must have the same length (got " +
                                std::to_string(a.size()) + " and " + std::to_string(b.size()) + ")");
}

py::array_t<double> ToArray(const std::vector<double> &v) {
  py::array_t<double> out((py::ssize_t)v.size());
  auto o = out.mutable_unchecked<1>();
  for (size_t i = 0; i < v.size(); i++) o((py::ssize_t)i) = v[i];
  return out;
}

py::dict MaterialInfo(const std::string &formula, const std::string &phase, double compoundCorrection) {
  const st::Material m = MakeMaterial(formula, phase, compoundCorrection);
  py::list constituents;
  for (const st::Constituent &c : m.constituents) constituents.append(py::make_tuple(c.z, c.atom_fraction));
  py::dict d;
  d["name"] = m.name;
  d["constituents"] = constituents;
  d["phase"] = st::PhaseName(m.phase);
  d["compound_correction"] = m.compound_correction;
  d["mean_atomic_weight_u"] = st::MeanAtomicWeight(m);
  return d;
}

py::dict Stopping(const std::string &formula, int z1, double massU, const Array &energies, const std::string &phase,
                  double compoundCorrection) {
  const st::Material m = MakeMaterial(formula, phase, compoundCorrection);
  const st::Projectile p = MakeProjectile(z1, massU);
  const std::vector<double> e = Values(energies);
  std::vector<double> el(e.size()), nu(e.size()), tot(e.size()), unc(e.size());
  for (size_t i = 0; i < e.size(); i++) {
    const st::StoppingValue v = st::Stopping(m, p, e[i]);
    el[i] = v.electronic;
    nu[i] = v.nuclear;
    tot[i] = v.total;
    unc[i] = v.rel_uncertainty;
  }
  py::dict d;
  d["electronic"] = ToArray(el);
  d["nuclear"] = ToArray(nu);
  d["total"] = ToArray(tot);
  d["rel_uncertainty"] = ToArray(unc);
  return d;
}

py::array_t<double> ArealDensityBetween(const std::string &formula, int z1, double massU, const Array &eIn, const Array &eOut,
                                        const std::string &phase, double compoundCorrection) {
  const st::Material m = MakeMaterial(formula, phase, compoundCorrection);
  const st::Projectile p = MakeProjectile(z1, massU);
  const std::vector<double> a = Values(eIn), b = Values(eOut);
  SameLength(a, b, "energy_in and energy_out");
  std::vector<double> out(a.size());
  for (size_t i = 0; i < a.size(); i++) out[i] = st::ArealDensityBetween(m, p, a[i], b[i]);
  return ToArray(out);
}

py::array_t<double> EnergyAfter(const std::string &formula, int z1, double massU, const Array &eIn, const Array &areal,
                                const std::string &phase, double compoundCorrection) {
  const st::Material m = MakeMaterial(formula, phase, compoundCorrection);
  const st::Projectile p = MakeProjectile(z1, massU);
  const std::vector<double> a = Values(eIn), n = Values(areal);
  SameLength(a, n, "energy_in and areal_density");
  std::vector<double> out(a.size());
  for (size_t i = 0; i < a.size(); i++) out[i] = st::EnergyAfter(m, p, a[i], n[i]);
  return ToArray(out);
}

py::array_t<double> Range(const std::string &formula, int z1, double massU, const Array &energies, const std::string &phase,
                          double compoundCorrection) {
  const st::Material m = MakeMaterial(formula, phase, compoundCorrection);
  const st::Projectile p = MakeProjectile(z1, massU);
  const std::vector<double> e = Values(energies);
  std::vector<double> out(e.size());
  for (size_t i = 0; i < e.size(); i++) out[i] = st::Range(m, p, e[i]);
  return ToArray(out);
}

py::dict Straggling(const std::string &formula, int z1, double massU, const Array &eIn, const Array &areal,
                    const Array &beamSigma, const std::string &model, const std::string &phase, double compoundCorrection) {
  const st::Material m = MakeMaterial(formula, phase, compoundCorrection);
  const st::Projectile p = MakeProjectile(z1, massU);
  const st::StragglingModel sm = st::StragglingModelFromName(model);
  const std::vector<double> a = Values(eIn), n = Values(areal), s0 = Values(beamSigma);
  SameLength(a, n, "energy_in and areal_density");
  SameLength(a, s0, "energy_in and beam_sigma");
  std::vector<double> out(a.size()), strag(a.size()), beam(a.size()), total(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    const st::StragglingWidths w = st::PropagateStraggling(sm, m, p, a[i], n[i], s0[i]);
    out[i] = w.energy_out_MeV;
    strag[i] = w.straggling;
    beam[i] = w.beam;
    total[i] = w.total;
  }
  py::dict d;
  d["energy_out"] = ToArray(out);
  d["straggling"] = ToArray(strag);
  d["beam"] = ToArray(beam);
  d["total"] = ToArray(total);
  return d;
}

std::string Provenance(const std::string &formula, int z1, double massU, const std::string &phase, double compoundCorrection) {
  return st::StoppingProvenance(MakeMaterial(formula, phase, compoundCorrection), MakeProjectile(z1, massU));
}

py::tuple LightProjectile(const std::string &name) {
  st::Projectile p;
  if (name == "p") p = st::Proton();
  else if (name == "d") p = st::Deuteron();
  else if (name == "t") p = st::Triton();
  else if (name == "3He") p = st::Helium3();
  else if (name == "a") p = st::Alpha();
  else throw std::invalid_argument("stopping: projectile '" + name + "' is not one of p, d, t, 3He, a (give z1 and mass_u)");
  return py::make_tuple(p.z, p.mass_u);
}

}  // namespace

void BindStopping(py::module_ &parent) {
  py::module_ m = parent.def_submodule("stopping", "stopping/: SRIM-2013 stopping, energy loss and straggling (cite SRIM)");
  m.attr("citation") = st::kSrimCitation;
  m.attr("implementation") = st::kSrimImplementation;
  m.def("material", &MaterialInfo, py::arg("formula"), py::arg("phase") = "", py::arg("compound_correction") = 1.0);
  m.def("light_projectile", &LightProjectile, py::arg("name"));
  m.def("stopping", &Stopping, py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("energy"), py::arg("phase") = "",
        py::arg("compound_correction") = 1.0);
  m.def("areal_density_between", &ArealDensityBetween, py::arg("formula"), py::arg("z1"), py::arg("mass_u"),
        py::arg("energy_in"), py::arg("energy_out"), py::arg("phase") = "", py::arg("compound_correction") = 1.0);
  m.def("energy_after", &EnergyAfter, py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("energy_in"),
        py::arg("areal_density"), py::arg("phase") = "", py::arg("compound_correction") = 1.0);
  m.def("range", &Range, py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("energy"), py::arg("phase") = "",
        py::arg("compound_correction") = 1.0);
  m.def("straggling", &Straggling, py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("energy_in"),
        py::arg("areal_density"), py::arg("beam_sigma"), py::arg("model") = "trim", py::arg("phase") = "",
        py::arg("compound_correction") = 1.0);
  m.def("provenance", &Provenance, py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("phase") = "",
        py::arg("compound_correction") = 1.0);
}
