// The engine's stopping per active nucleus (TargetStopping, SetupSpecs.h StoppingSpec): epsilon_active = 1e-21 x the
// material's stopping per formula unit / the active atoms per formula unit, checked on Ta2O5 with 18O active and on
// TiN with 14N active against the same number put together by hand from the library's per-atom values of each element
// in the solid (Bragg's rule: per formula unit sum_i n_i eps_i), and TiN against the value azure2-stopping prints.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "TargetStopping.h"
#include "azure2/stopping/Stopping.h"

namespace {
int failures = 0;
void Check(const std::string &what, double got, double want, double tol) {
  const double rel = std::fabs(got / want - 1.0);
  const bool ok = rel <= tol;
  std::printf("  %s  %s: %.10e against %.10e (%.1e, tolerance %.0e)\n", ok ? "ok  " : "FAIL", what.c_str(), got, want, rel, tol);
  if (!ok) failures++;
}

// The library's stopping per atom of element z alone, in the solid (eV/(1e15 atoms/cm^2)).
double ElementInSolid(int z, const azure2::stopping::Projectile &p, double e) {
  azure2::stopping::Material m = azure2::stopping::MaterialFromFractions("e", {{z, 1.0}});
  m.phase = azure2::stopping::Phase::kSolid;
  return azure2::stopping::Stopping(m, p, e).total;
}

azure2::experiment::StoppingSpec Spec(const std::string &material, std::vector<std::pair<int, double>> composition, int activeZ) {
  azure2::experiment::StoppingSpec s;
  s.source = azure2::experiment::StoppingSpec::Source::Evaluation;
  s.material = material;
  s.composition = composition;
  s.activeZ = activeZ;
  s.projectileZ = 1;
  s.projectileA = 1;
  s.projectileMass = 1.007825;
  return s;
}
}  // namespace

int main() {
  const azure2::stopping::Projectile p = azure2::stopping::Ion(1, 1.007825, "1H");
  for (double e : {0.1, 0.5, 2.0}) {
    // Ta2O5, 18O active: per formula unit 2 eps_Ta + 5 eps_O, five O atoms
    const TargetStopping ta2o5(Spec("Ta2O5", {{73, 2.0}, {8, 5.0}}, 8));
    const double hand = 1e-21 * (2.0 * ElementInSolid(73, p, e) + 5.0 * ElementInSolid(8, p, e)) / 5.0;
    Check("Ta2O5, 18O active, E = " + std::to_string(e) + " MeV", ta2o5.LibraryPerActiveNucleus(e), hand, 1e-14);
    // TiN, 14N active: per formula unit eps_Ti + eps_N, one N atom
    const TargetStopping tin(Spec("TiN", {{22, 1.0}, {7, 1.0}}, 7));
    const double handTiN = 1e-21 * (ElementInSolid(22, p, e) + ElementInSolid(7, p, e)) / 1.0;
    Check("TiN, 14N active, E = " + std::to_string(e) + " MeV", tin.LibraryPerActiveNucleus(e), handTiN, 1e-14);
  }
  // the value the command-line tool prints for p (1.007825 u) in TiN at 0.5 MeV: 12.00133921 eV/(1e15 atoms/cm^2) per
  // atom of TiN, i.e. 2 x 12.00133921e-21 MeV cm^2 per N atom
  const TargetStopping tin(Spec("TiN", {{22, 1.0}, {7, 1.0}}, 7));
  Check("TiN at 0.5 MeV against azure2-stopping's 12.00133921 per atom", tin.LibraryPerActiveNucleus(0.5),
        2.0 * 12.00133921e-21, 1e-9);
  if (failures) {
    std::printf("stopping_convention: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("stopping_convention: all checks passed\n");
  return 0;
}
