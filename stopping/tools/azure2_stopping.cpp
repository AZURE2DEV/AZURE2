// azure2-stopping: the stopping module from the command line (docs/upgrade/stopping.md section 11).
//
//   azure2-stopping stopping   PROJECTILE MATERIAL --from E1 --to E2 [--points n]
//   azure2-stopping range      PROJECTILE MATERIAL --from E1 --to E2 [--points n]
//   azure2-stopping layer      PROJECTILE MATERIAL --energy E --areal-density N [--beam-sigma S] [--straggling trim|bohr]
//   azure2-stopping provenance PROJECTILE MATERIAL
//
//   PROJECTILE: --projectile p|d|t|3He|a, or --ion-z Z --ion-mass M (u) for any ion Z = 1..92
//   MATERIAL:   --material F (a chemical formula, Ta2O5) [--phase solid|gas] [--compound-correction C]
//
// The phase defaults to SRIM's (gas for a single natural-gas element, solid otherwise); a gaseous compound needs
// --phase gas.  Energies in MeV (lab), stopping in eV/(1e15 atoms/cm^2), areal densities in atoms/cm^2, widths in MeV
// (1 sigma).  Energy grids are logarithmic.  Output: comment lines (starting with '#') naming SRIM-2013 with the
// citation it asks for, the material and phase, then tab-separated columns with a header line naming each column and
// its unit.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <map>
#include <stdexcept>
#include <string>

#include "azure2/stopping/EnergyLoss.h"
#include "azure2/stopping/Stopping.h"
#include "azure2/stopping/Straggling.h"

using namespace azure2::stopping;

namespace {

int Usage() {
  std::fprintf(stderr,
               "usage: azure2-stopping stopping|range PROJECTILE MATERIAL --from E1 --to E2 [--points n]\n"
               "       azure2-stopping layer PROJECTILE MATERIAL --energy E --areal-density N [--beam-sigma S] "
               "[--straggling trim|bohr]\n"
               "       azure2-stopping provenance PROJECTILE MATERIAL\n"
               "PROJECTILE: --projectile p|d|t|3He|a, or --ion-z Z --ion-mass M\n"
               "MATERIAL:   --material F [--phase solid|gas] [--compound-correction C]\n"
               "units: MeV, eV/(1e15 atoms/cm^2), atoms/cm^2\n");
  return 2;
}

double Number(const std::map<std::string, std::string> &opt, const std::string &key) {
  const auto it = opt.find(key);
  if (it == opt.end()) throw std::invalid_argument("missing " + key);
  char *end = nullptr;
  const double v = std::strtod(it->second.c_str(), &end);
  if (*end != '\0') throw std::invalid_argument(key + " '" + it->second + "' is not a number");
  return v;
}

Projectile ProjectileFrom(const std::map<std::string, std::string> &opt) {
  if (opt.count("--ion-z") || opt.count("--ion-mass")) {
    const double z = Number(opt, "--ion-z");
    if (z != std::floor(z)) throw std::invalid_argument("--ion-z must be an integer");
    return Ion(static_cast<int>(z), Number(opt, "--ion-mass"), "Z=" + opt.at("--ion-z"));
  }
  const auto it = opt.find("--projectile");
  if (it == opt.end()) throw std::invalid_argument("missing --projectile (or --ion-z and --ion-mass)");
  const std::string &name = it->second;
  if (name == "p") return Proton();
  if (name == "d") return Deuteron();
  if (name == "t") return Triton();
  if (name == "3He") return Helium3();
  if (name == "a") return Alpha();
  throw std::invalid_argument("projectile '" + name + "': expected p, d, t, 3He or a (other ions: --ion-z, --ion-mass)");
}

Material MaterialFrom(const std::map<std::string, std::string> &opt) {
  const auto it = opt.find("--material");
  if (it == opt.end()) throw std::invalid_argument("missing --material");
  Material m = MaterialFromFormula(it->second);
  if (opt.count("--phase")) m.phase = PhaseFromName(opt.at("--phase"));
  if (opt.count("--compound-correction")) m.compound_correction = Number(opt, "--compound-correction");
  return m;
}

void Header(const Projectile &p, const Material &m) {
  std::printf("# %s\n# cite SRIM: %s\n# %s (Z1 = %d, %.8g u) in %s, phase %s, compound correction %g\n",
              kSrimImplementation, kSrimCitation, p.label.c_str(), p.z, p.mass_u, m.name.c_str(),
              PhaseName(m.phase).c_str(), m.compound_correction);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) return Usage();
  const std::string command = argv[1];
  if ((argc - 2) % 2 != 0) return Usage();
  std::map<std::string, std::string> opt;
  for (int i = 2; i + 1 < argc; i += 2) opt[argv[i]] = argv[i + 1];
  try {
    const Projectile p = ProjectileFrom(opt);
    const Material m = MaterialFrom(opt);
    if (command == "provenance") {
      std::printf("%s", StoppingProvenance(m, p).c_str());
      return 0;
    }
    if (command == "stopping" || command == "range") {
      const double e1 = Number(opt, "--from"), e2 = Number(opt, "--to");
      const int n = opt.count("--points") ? static_cast<int>(Number(opt, "--points")) : 50;
      if (!(e1 > 0.0) || !(e2 > e1) || n < 2) throw std::invalid_argument("need 0 < --from < --to and --points >= 2");
      Header(p, m);
      if (command == "stopping") {
        std::printf("E_MeV\teps_electronic_eV_per_1e15\teps_nuclear_eV_per_1e15\teps_total_eV_per_1e15\trel_uncertainty\n");
      } else {
        std::printf("E_MeV\trange_atoms_per_cm2\n");
      }
      for (int i = 0; i < n; ++i) {
        const double e = e1 * std::pow(e2 / e1, static_cast<double>(i) / (n - 1));
        if (command == "stopping") {
          const StoppingValue v = Stopping(m, p, e);
          std::printf("%.10g\t%.10g\t%.10g\t%.10g\t%.3g\n", e, v.electronic, v.nuclear, v.total, v.rel_uncertainty);
        } else {
          std::printf("%.6g\t%.6g\n", e, Range(m, p, e));
        }
      }
      return 0;
    }
    if (command == "layer") {
      const double e = Number(opt, "--energy"), n = Number(opt, "--areal-density");
      const double beam = opt.count("--beam-sigma") ? Number(opt, "--beam-sigma") : 0.0;
      const StragglingModel model =
          opt.count("--straggling") ? StragglingModelFromName(opt.at("--straggling")) : StragglingModel::kTrim;
      const StragglingWidths w = PropagateStraggling(model, m, p, e, n, beam);
      Header(p, m);
      std::printf("# straggling model: %s\n", StragglingModelName(model).c_str());
      std::printf("E_in_MeV\tareal_density_atoms_per_cm2\tE_out_MeV\tdelta_E_MeV\tsigma_straggling_MeV\tsigma_beam_MeV\tsigma_total_MeV\n");
      std::printf("%.8g\t%.6g\t%.8g\t%.6g\t%.4g\t%.4g\t%.4g\n", e, n, w.energy_out_MeV, e - w.energy_out_MeV, w.straggling,
                  w.beam, w.total);
      return 0;
    }
    return Usage();
  } catch (const std::exception &ex) {
    std::fprintf(stderr, "azure2-stopping: %s\n", ex.what());
    return 1;
  }
}
