// The rows of the stopping table (include/azure2/physics_table/StoppingTable.h): one row per function of stopping/,
// read by the binding (api/src/pybind_stopping.cpp: pyazr.stopping, pyazr.kinds()["stopping"]) and by the front ends'
// facade (api/front: StoppingFunctions, EvaluateStopping).  No physics in this file: an adapter calls the library.
#include "azure2/physics_table/StoppingTable.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/EnergyLoss.h"
#include "azure2/stopping/Stopping.h"
#include "azure2/stopping/Straggling.h"
#include "azure2/stopping/Tables.h"

namespace azure2 {
namespace physics_table {

namespace {

namespace st = azure2::stopping;

const char *const kStoppingUnit = "eV/(1e15 atoms/cm^2)";
const char *const kArealUnit = "atoms/cm^2";

ArgSpec FormulaArgument() {
  return Text("formula", "chemical formula of the target material (counts may be decimal: Ta2O5, TiN, CH4)",
                    "Ta2O5");
}
std::vector<ArgSpec> CommonArguments() {
  return {Choice("phase", "phase of the material (default: SRIM's: gas for a gaseous element, solid otherwise)",
                       StoppingPhases(), 0),
          Real("compound_correction", "", "SRIM's factor of Bragg's rule (1: Bragg's rule)").gt(0).otherwise(1),
          Integer("z1", "charge number of the projectile").ge(1).le(92).example(1),
          Real("mass_u", "u", "atomic mass of the projectile").gt(0).example(st::Proton().mass_u)};
}
ArgSpec Energy(const std::string &name, const std::string &doc) {
  return Real(name, "MeV", doc).gt(0).example(1).sweep(0.001, 10, Spacing::log);
}
ArgSpec Areal() {
  return Real("areal_density", kArealUnit, "areal density of the layer (atoms of the material)").ge(0)
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
                 return kOk;
               }});
  t.push_back({"range", "Range R(E) = N(E, 0) of an ion in a material (atoms/cm^2).",
               "The areal density over which the ion stops: the integral of dE / epsilon from 0 to E.",
               "stopping/EnergyLoss.h: Range", "energy", {Energy("energy", "lab kinetic energy of the projectile")},
               {{"range", kArealUnit, "range"}},
               [](const StoppingCall &c, const double *a, double *o) {
                 o[0] = st::Range(c.material, c.projectile, a[0]);
                 return kOk;
               }});
  t.push_back({"energy_after", "Energy of an ion after a layer of given areal density (MeV).",
               "E_out with N = int_{E_out}^{E_in} dE / epsilon.  Status -2 where the layer is at least the range.",
               "stopping/EnergyLoss.h: EnergyAfter", "areal_density",
               {Energy("energy_in", "lab kinetic energy before the layer"), Areal()},
               {{"energy_out", "MeV", "energy after the layer"}},
               [](const StoppingCall &c, const double *a, double *o) {
                 o[0] = st::EnergyAfter(c.material, c.projectile, a[0], a[1]);
                 return kOk;
               }});
  t.push_back({"areal_density_between", "Areal density between two energies of an ion (atoms/cm^2).",
               "N(E_in, E_out) = int_{E_out}^{E_in} dE / epsilon.  Status -2 for E_out > E_in.",
               "stopping/EnergyLoss.h: ArealDensityBetween", "energy_out",
               {Energy("energy_in", "lab kinetic energy before the layer"),
                Real("energy_out", "MeV", "lab kinetic energy after the layer").ge(0).example(0.5)
                    .sweep(0, 1, Spacing::linear)},
               {{"areal_density", kArealUnit, "areal density"}},
               [](const StoppingCall &c, const double *a, double *o) {
                 o[0] = st::ArealDensityBetween(c.material, c.projectile, a[0], a[1]);
                 return kOk;
               }});
  t.push_back({"energy_straggling", "Energy after a layer and its spread: the layer's straggling and the beam's (MeV).",
               "TRIM's or Bohr's straggling with Tschalaer's propagation; the beam's entrance spread beam_sigma "
               "carried through the layer; total = both in quadrature (1 sigma).  Status -2 where the layer is at "
               "least the range.",
               "stopping/Straggling.h: PropagateStraggling", "areal_density",
               {Energy("energy_in", "lab kinetic energy before the layer"), Areal(),
                Real("beam_sigma", "MeV", "the beam's energy spread at the entrance (1 sigma)").ge(0)
                    .otherwise(0),
                Choice("model", "straggling model", {"trim", "bohr"}, 0)},
               {{"energy_out", "MeV", "energy after the layer"}, {"straggling", "MeV", "the layer's straggling"},
                {"beam", "MeV", "the entrance spread at the exit"}, {"total", "MeV", "both in quadrature"}},
               [](const StoppingCall &c, const double *a, double *o) {
                 const st::StragglingModel m = st::StragglingModelFromName(a[3] == 0 ? "trim" : "bohr");
                 const st::StragglingWidths w = st::PropagateStraggling(m, c.material, c.projectile, a[0], a[1], a[2]);
                 o[0] = w.energy_out_MeV, o[1] = w.straggling, o[2] = w.beam, o[3] = w.total;
                 return kOk;
               }});
  return t;
}

}  // namespace

const std::vector<std::string> &StoppingPhases() {
  static const std::vector<std::string> phases = {"default", "solid", "gas"};  // "default": SRIM's phase
  return phases;
}

ArgSpec StoppingFormula() { return FormulaArgument(); }
std::vector<ArgSpec> StoppingCommon() { return CommonArguments(); }

const std::vector<StoppingRow> &StoppingTable() {
  static const std::vector<StoppingRow> table = BuildStoppingTable();  // immutable after the first call
  return table;
}

const StoppingRow *FindStoppingRow(const std::string &name) {
  for (const StoppingRow &r : StoppingTable())
    if (r.name == name) return &r;
  return nullptr;
}

std::vector<ArgSpec> StoppingNumericArgs(const StoppingRow &row) {
  std::vector<ArgSpec> all = StoppingCommon();
  all.insert(all.end(), row.args.begin(), row.args.end());
  return all;
}

st::Material StoppingMaterial(const std::string &formula, const std::string &phase, double compoundCorrection) {
  st::Material m = st::MaterialFromFormula(formula);
  if (!phase.empty()) m.phase = st::PhaseFromName(phase);
  m.compound_correction = compoundCorrection;
  return m;
}

st::Projectile StoppingProjectile(int z1, double massU) {
  if (!(massU > 0.0)) throw std::invalid_argument("stopping: the projectile's mass must be > 0 u");
  return st::Ion(z1, massU, "Z1=" + std::to_string(z1));
}

void EvaluateStoppingRow(const StoppingRow &row, const std::string &formula, const std::vector<const double *> &columns,
                         std::size_t n, double *v, int *stp) {
  const std::vector<ArgSpec> spec = StoppingNumericArgs(row);
  const std::size_t na = spec.size(), no = row.outs.size();
  std::vector<double> a(na), o(no);
  std::vector<double> made(kStoppingCommon, kNaN);  // the common arguments of the material and projectile made
  st::Material material;
  st::Projectile projectile;
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t k = 0; k < na; k++) a[k] = columns[k][i];
    int s = CheckDomain(spec, a.data());
    if (s == kOk) {
      try {
        if (!std::equal(made.begin(), made.end(), a.begin())) {
          const int phase = static_cast<int>(a[0]);
          material = StoppingMaterial(formula, phase == 0 ? "" : StoppingPhases()[static_cast<std::size_t>(phase)], a[1]);
          projectile = StoppingProjectile(static_cast<int>(a[2]), a[3]);
          std::copy(a.begin(), a.begin() + kStoppingCommon, made.begin());
        }
        s = row.call(StoppingCall{material, projectile}, a.data() + kStoppingCommon, o.data());
      } catch (const std::domain_error &) {
        s = kJointly;  // the layer at least the range, an energy order the integral refuses
      } catch (const std::invalid_argument &) {
        s = kJointly;
      } catch (const std::exception &) {
        s = kFailed;
      }
    }
    for (std::size_t j = 0; j < no; j++) v[j * n + i] = s == kOk ? o[j] : kNaN;
    stp[i] = s;
  }
}

}  // namespace physics_table
}  // namespace azure2
