// The stopping/ library in the _azure2 module (submodule _azure2.stopping; pyazr.stopping wraps it): SRIM-2013's
// stopping powers ported from LunaSRIM, energy loss, range and straggling of a layer, provenance
// (docs/upgrade/stopping.md).  A thin binding: arrays in, arrays out, the C++ functions element-wise; no second
// implementation.  The same functions are also rows of a table of the physics table's shape (pybind_rows.h: arguments
// with unit, domain, example and sweep; returned fields with units; and SRIM's citation), evaluated element-wise
// with a status per element (functions(), evaluate(): pyazr.kinds()["stopping"], the web UI's Physics tab).
// Units as stopping/'s Material.h: energies MeV (lab, of the projectile), stopping cross sections
// eV/(1e15 atoms/cm^2) per atom of the material, areal densities atoms/cm^2 of the material, straggling widths MeV.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "pybind_rows.h"

#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/EnergyLoss.h"
#include "azure2/stopping/Stopping.h"
#include "azure2/stopping/Straggling.h"
#include "pybind_engine_call.h"

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

// ---------------------------------------------------------------- the table of stopping functions

namespace rows = azure2_rows;
using rows::ArgSpec;
using rows::OutSpec;
using rows::Spacing;

/// What a row's adapter is given besides the element's numbers: the material and the projectile of the element.
struct StoppingCall {
  const st::Material &material;
  const st::Projectile &projectile;
};

/// The adapter of a stopping row: the element's numeric arguments after the material and projectile ones
/// (a[0 .. n)), its fields o[0 .. n_outs); returns a status (pybind_rows.h).
using StoppingAdapter = int (*)(const StoppingCall &c, const double *a, double *o);

/// A row: the material (formula as text, phase, compound correction) and projectile (z1, mass_u) arguments come
/// first in every row (kCommon), then the row's own; returns; the adapter.
struct StoppingRow {
  std::string name, summary, doc, library, axis;
  std::vector<ArgSpec> args;  ///< the row's own numeric arguments, after the common ones
  std::vector<OutSpec> outs;
  StoppingAdapter call;
};

const std::vector<std::string> kPhases = {"default", "solid", "gas"};  // "default": SRIM's phase of the material
constexpr std::size_t kCommon = 4;                                     // phase, compound_correction, z1, mass_u
const char *const kStoppingUnit = "eV/(1e15 atoms/cm^2)";
const char *const kArealUnit = "atoms/cm^2";

ArgSpec Formula() {
  return rows::Text("formula", "chemical formula of the target material (counts may be decimal: Ta2O5, TiN, CH4)",
                    "Ta2O5");
}
std::vector<ArgSpec> Common() {
  return {rows::Choice("phase", "phase of the material (default: SRIM's: gas for a gaseous element, solid otherwise)",
                       kPhases, 0),
          rows::Real("compound_correction", "", "SRIM's factor of Bragg's rule (1: Bragg's rule)").gt(0).otherwise(1),
          rows::Integer("z1", "charge number of the projectile").ge(1).le(92).example(1),
          rows::Real("mass_u", "u", "atomic mass of the projectile").gt(0).example(st::Proton().mass_u)};
}
ArgSpec Energy(const std::string &name, const std::string &doc) {
  return rows::Real(name, "MeV", doc).gt(0).example(1).sweep(0.001, 10, Spacing::log);
}
ArgSpec Areal() {
  return rows::Real("areal_density", kArealUnit, "areal density of the layer (atoms of the material)").ge(0)
      .example(1e18).sweep(0, 1e19, Spacing::linear);
}

std::vector<StoppingRow> BuildStoppingTable() {
  std::vector<StoppingRow> t;
  t.push_back({"stopping_power", "Stopping cross section of a material for an ion (eV/(1e15 atoms/cm^2)).",
               "Per atom of the material (all atoms counted), electronic (SRIM-2013) and nuclear (the universal "
               "ZBL potential), and SRIM's published accuracy for the ion class as a relative 1-sigma uncertainty "
               "(docs/upgrade/stopping.md).", "stopping/Stopping.h: Stopping", "energy",
               {Energy("energy", "lab kinetic energy of the projectile")},
               {{"electronic", kStoppingUnit, "electronic stopping"}, {"nuclear", kStoppingUnit, "nuclear stopping"},
                {"total", kStoppingUnit, "electronic + nuclear"},
                {"rel_uncertainty", "", "relative 1-sigma uncertainty"}},
               [](const StoppingCall &c, const double *a, double *o) {
                 const st::StoppingValue v = st::Stopping(c.material, c.projectile, a[0]);
                 o[0] = v.electronic, o[1] = v.nuclear, o[2] = v.total, o[3] = v.rel_uncertainty;
                 return rows::kOk;
               }});
  t.push_back({"range", "Range R(E) = N(E, 0) of an ion in a material (atoms/cm^2).",
               "The areal density over which the ion stops: the integral of dE / epsilon from 0 to E.",
               "stopping/EnergyLoss.h: Range", "energy", {Energy("energy", "lab kinetic energy of the projectile")},
               {{"range", kArealUnit, "range"}},
               [](const StoppingCall &c, const double *a, double *o) {
                 o[0] = st::Range(c.material, c.projectile, a[0]);
                 return rows::kOk;
               }});
  t.push_back({"energy_after", "Energy of an ion after a layer of given areal density (MeV).",
               "E_out with N = int_{E_out}^{E_in} dE / epsilon.  Status -2 where the layer is at least the range.",
               "stopping/EnergyLoss.h: EnergyAfter", "areal_density",
               {Energy("energy_in", "lab kinetic energy before the layer"), Areal()},
               {{"energy_out", "MeV", "energy after the layer"}},
               [](const StoppingCall &c, const double *a, double *o) {
                 o[0] = st::EnergyAfter(c.material, c.projectile, a[0], a[1]);
                 return rows::kOk;
               }});
  t.push_back({"areal_density_between", "Areal density between two energies of an ion (atoms/cm^2).",
               "N(E_in, E_out) = int_{E_out}^{E_in} dE / epsilon.  Status -2 for E_out > E_in.",
               "stopping/EnergyLoss.h: ArealDensityBetween", "energy_out",
               {Energy("energy_in", "lab kinetic energy before the layer"),
                rows::Real("energy_out", "MeV", "lab kinetic energy after the layer").ge(0).example(0.5)
                    .sweep(0, 1, Spacing::linear)},
               {{"areal_density", kArealUnit, "areal density"}},
               [](const StoppingCall &c, const double *a, double *o) {
                 o[0] = st::ArealDensityBetween(c.material, c.projectile, a[0], a[1]);
                 return rows::kOk;
               }});
  t.push_back({"energy_straggling", "Energy after a layer and its spread: the layer's straggling and the beam's (MeV).",
               "TRIM's or Bohr's straggling with Tschalaer's propagation; the beam's entrance spread beam_sigma "
               "carried through the layer; total = both in quadrature (1 sigma).  Status -2 where the layer is at "
               "least the range.",
               "stopping/Straggling.h: PropagateStraggling", "areal_density",
               {Energy("energy_in", "lab kinetic energy before the layer"), Areal(),
                rows::Real("beam_sigma", "MeV", "the beam's energy spread at the entrance (1 sigma)").ge(0)
                    .otherwise(0),
                rows::Choice("model", "straggling model", {"trim", "bohr"}, 0)},
               {{"energy_out", "MeV", "energy after the layer"}, {"straggling", "MeV", "the layer's straggling"},
                {"beam", "MeV", "the entrance spread at the exit"}, {"total", "MeV", "both in quadrature"}},
               [](const StoppingCall &c, const double *a, double *o) {
                 const st::StragglingModel m = st::StragglingModelFromName(a[3] == 0 ? "trim" : "bohr");
                 const st::StragglingWidths w = st::PropagateStraggling(m, c.material, c.projectile, a[0], a[1], a[2]);
                 o[0] = w.energy_out_MeV, o[1] = w.straggling, o[2] = w.beam, o[3] = w.total;
                 return rows::kOk;
               }});
  return t;
}

const std::vector<StoppingRow> &StoppingTable() {
  static const std::vector<StoppingRow> table = BuildStoppingTable();  // immutable after the first call
  return table;
}

const StoppingRow &FindStoppingRow(const std::string &name) {
  for (const StoppingRow &r : StoppingTable())
    if (r.name == name) return r;
  throw std::invalid_argument("stopping: no function named '" + name + "' (pyazr.kinds()['stopping'] lists them)");
}

/// Every numeric argument of a row in order: the common ones, then the row's own.
std::vector<ArgSpec> NumericArgs(const StoppingRow &r) {
  std::vector<ArgSpec> all = Common();
  all.insert(all.end(), r.args.begin(), r.args.end());
  return all;
}

py::list DescribeStopping() {
  py::list out;
  for (const StoppingRow &r : StoppingTable()) {
    py::list args, outs;
    args.append(rows::DescribeArgument(Formula()));
    for (const ArgSpec &a : NumericArgs(r)) args.append(rows::DescribeArgument(a));
    for (const OutSpec &o : r.outs) outs.append(rows::DescribeReturn(o));
    py::dict d;
    d["name"] = r.name;
    d["summary"] = r.summary;
    d["doc"] = r.doc;
    d["library"] = r.library;
    d["axis"] = r.axis;
    d["group"] = "stopping";
    d["citation"] = st::kSrimCitation;
    d["implementation"] = st::kSrimImplementation;
    d["arguments"] = args;
    d["returns"] = outs;
    out.append(d);
  }
  return out;
}

using Column = py::array_t<double, py::array::c_style | py::array::forcecast>;

/// Row `name` element-wise for the material `formula` on equal-length 1-d arrays, one per numeric argument in the
/// row's order (the common ones first): (values (n_returns, n), status (n,)), NaN where the status is not 0.  The
/// material and projectile are made once per distinct (phase, compound correction, z1, mass) in a row of elements.
py::tuple EvaluateStopping(const std::string &name, const std::string &formula, const std::vector<Column> &columns) {
  const StoppingRow &row = FindStoppingRow(name);
  const std::vector<ArgSpec> spec = NumericArgs(row);
  const std::size_t na = spec.size(), no = row.outs.size();
  if (columns.size() != na)
    throw std::invalid_argument("stopping " + name + ": " + std::to_string(na) + " argument arrays expected, " +
                                std::to_string(columns.size()) + " given");
  const py::ssize_t n = columns[0].size();
  std::vector<const double *> in(na);
  for (std::size_t k = 0; k < na; k++) {
    if (columns[k].ndim() != 1 || columns[k].size() != n)
      throw std::invalid_argument("stopping " + name + ": the argument arrays must be 1-d of one length");
    in[k] = columns[k].data();
  }
  st::MaterialFromFormula(formula);  // a formula the library refuses is refused here, before any element
  py::array_t<double> values({static_cast<py::ssize_t>(no), n});
  py::array_t<int> status(n);
  double *v = values.mutable_data();
  int *stp = status.mutable_data();
  {
    EngineCallScope engine;
    std::vector<double> a(na), o(no);
    std::vector<double> made(kCommon, rows::kNaN);  // the common arguments of the material and projectile made
    st::Material material;
    st::Projectile projectile;
    for (py::ssize_t i = 0; i < n; i++) {
      for (std::size_t k = 0; k < na; k++) a[k] = in[k][i];
      int s = rows::CheckDomain(spec, a.data());
      if (s == rows::kOk) {
        try {
          if (!std::equal(made.begin(), made.end(), a.begin())) {
            const int phase = static_cast<int>(a[0]);
            material = MakeMaterial(formula, phase == 0 ? "" : kPhases[static_cast<std::size_t>(phase)], a[1]);
            projectile = MakeProjectile(static_cast<int>(a[2]), a[3]);
            std::copy(a.begin(), a.begin() + kCommon, made.begin());
          }
          s = row.call(StoppingCall{material, projectile}, a.data() + kCommon, o.data());
        } catch (const std::domain_error &) {
          s = rows::kJointly;  // the layer at least the range, an energy order the integral refuses
        } catch (const std::invalid_argument &) {
          s = rows::kJointly;
        } catch (const std::exception &) {
          s = rows::kFailed;
        }
      }
      for (std::size_t j = 0; j < no; j++) v[j * n + i] = s == rows::kOk ? o[j] : rows::kNaN;
      stp[i] = s;
    }
  }
  return py::make_tuple(values, status);
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
  m.def("material", &MaterialInfo, PrecisionGuard(), py::arg("formula"), py::arg("phase") = "", py::arg("compound_correction") = 1.0);
  m.def("light_projectile", &LightProjectile, py::arg("name"));
  m.def("stopping", &Stopping, PrecisionGuard(), py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("energy"), py::arg("phase") = "",
        py::arg("compound_correction") = 1.0);
  m.def("areal_density_between", &ArealDensityBetween, PrecisionGuard(), py::arg("formula"), py::arg("z1"), py::arg("mass_u"),
        py::arg("energy_in"), py::arg("energy_out"), py::arg("phase") = "", py::arg("compound_correction") = 1.0);
  m.def("energy_after", &EnergyAfter, PrecisionGuard(), py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("energy_in"),
        py::arg("areal_density"), py::arg("phase") = "", py::arg("compound_correction") = 1.0);
  m.def("range", &Range, PrecisionGuard(), py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("energy"), py::arg("phase") = "",
        py::arg("compound_correction") = 1.0);
  m.def("straggling", &Straggling, PrecisionGuard(), py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("energy_in"),
        py::arg("areal_density"), py::arg("beam_sigma"), py::arg("model") = "trim", py::arg("phase") = "",
        py::arg("compound_correction") = 1.0);
  m.def("functions", &DescribeStopping,
        "The table of stopping functions, one dict per row of the physics table's shape (pybind_rows.h): name, "
        "summary, doc, library, axis, group, citation (SRIM's, owed by every use), implementation, arguments "
        "(formula, a text, then the numeric ones: phase, compound_correction, z1, mass_u and the row's own) and "
        "returns; plain JSON data.");
  m.def("evaluate", &EvaluateStopping, PrecisionGuard(), py::arg("name"), py::arg("formula"), py::arg("arguments"),
        "Row `name` for the material `formula` element-wise on equal-length 1-d float arrays, one per numeric "
        "argument in the row's order (a choice as its index): (values, status), values of shape (n_returns, n), NaN "
        "where status != 0; status 0 computed, k > 0 numeric argument k (1-based) outside its domain, -2 the "
        "arguments together outside the domain (a layer at least the range), -1 failed.");
  m.def("provenance", &Provenance, PrecisionGuard(), py::arg("formula"), py::arg("z1"), py::arg("mass_u"), py::arg("phase") = "",
        py::arg("compound_correction") = 1.0);
}
