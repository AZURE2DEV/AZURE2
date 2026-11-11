// azure2-stopping: the stopping module from the command line (docs/upgrade/stopping.md section 11).
//
//   azure2-stopping stopping   --projectile P --material F [--phase solid|liquid|gas] --from E1 --to E2 [--points n]
//   azure2-stopping range      --projectile P --material F --from E1 --to E2 [--points n]
//   azure2-stopping layer      --projectile P --material F --energy E --areal-density N [--beam-sigma S]
//   azure2-stopping provenance --projectile P --material F
//
// P is p, d, t, 3He or a; F a chemical formula (Ta2O5).  Energies in MeV (lab), stopping in eV/(1e15 atoms/cm^2),
// areal densities in atoms/cm^2, widths in MeV (1 sigma).  Energy grids are logarithmic.  Output: tab-separated
// columns with a header line naming each column and its unit.
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
               "usage: azure2-stopping stopping|range --projectile P --material F [--phase solid|liquid|gas] "
               "--from E1 --to E2 [--points n]\n"
               "       azure2-stopping layer --projectile P --material F --energy E --areal-density N [--beam-sigma S]\n"
               "       azure2-stopping provenance --projectile P --material F\n"
               "P: p d t 3He a; F: a formula (Ta2O5); MeV, eV/(1e15 atoms/cm^2), atoms/cm^2\n");
  return 2;
}

Projectile ProjectileNamed(const std::string &name) {
  if (name == "p") return Proton();
  if (name == "d") return Deuteron();
  if (name == "t") return Triton();
  if (name == "3He") return Helium3();
  if (name == "a") return Alpha();
  throw std::invalid_argument("projectile '" + name + "': expected p, d, t, 3He or a");
}

double Number(const std::map<std::string, std::string> &opt, const std::string &key) {
  const auto it = opt.find(key);
  if (it == opt.end()) throw std::invalid_argument("missing " + key);
  char *end = nullptr;
  const double v = std::strtod(it->second.c_str(), &end);
  if (*end != '\0') throw std::invalid_argument(key + " '" + it->second + "' is not a number");
  return v;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) return Usage();
  const std::string command = argv[1];
  std::map<std::string, std::string> opt;
  for (int i = 2; i + 1 < argc; i += 2) opt[argv[i]] = argv[i + 1];
  if ((argc - 2) % 2 != 0) return Usage();
  try {
    if (!opt.count("--projectile") || !opt.count("--material")) return Usage();
    const Projectile p = ProjectileNamed(opt["--projectile"]);
    const Phase phase = opt.count("--phase") ? PhaseFromName(opt["--phase"]) : Phase::kSolid;
    const Material m = MaterialFromFormula(opt["--material"], phase, 0.0);
    if (command == "provenance") {
      std::printf("%s", StoppingProvenance(m, p).c_str());
      return 0;
    }
    if (command == "stopping" || command == "range") {
      const double e1 = Number(opt, "--from"), e2 = Number(opt, "--to");
      const int n = opt.count("--points") ? static_cast<int>(Number(opt, "--points")) : 50;
      if (!(e1 > 0.0) || !(e2 > e1) || n < 2) throw std::invalid_argument("need 0 < --from < --to and --points >= 2");
      if (command == "stopping") {
        std::printf("E_MeV\teps_electronic_eV_per_1e15\teps_nuclear_eV_per_1e15\teps_total_eV_per_1e15\trel_uncertainty\n");
      } else {
        std::printf("E_MeV\trange_atoms_per_cm2\n");
      }
      for (int i = 0; i < n; ++i) {
        const double e = e1 * std::pow(e2 / e1, static_cast<double>(i) / (n - 1));
        if (command == "stopping") {
          const StoppingValue v = Stopping(m, p, e);
          std::printf("%.6g\t%.6g\t%.6g\t%.6g\t%.3g\n", e, v.electronic, v.nuclear, v.total, v.rel_uncertainty);
        } else {
          std::printf("%.6g\t%.6g\n", e, Range(m, p, e));
        }
      }
      return 0;
    }
    if (command == "layer") {
      const double e = Number(opt, "--energy"), n = Number(opt, "--areal-density");
      const double beam = opt.count("--beam-sigma") ? Number(opt, "--beam-sigma") : 0.0;
      const StragglingWidths w = PropagateStraggling(StragglingModel{}, m, p, e, n, beam);
      std::printf("E_in_MeV\tareal_density_atoms_per_cm2\tE_out_MeV\tdelta_E_MeV\tsigma_bohr_MeV\tsigma_beam_MeV\tsigma_total_MeV\n");
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
