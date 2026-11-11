// The comparison of this module's stopping with ERYA's (docs/upgrade/stopping.md section 9): not a test, a table.
// Run from external/ (where ERYA finds erya/data/SRIM2013.xml), e.g. through tools/compare_erya.sh.
//
// ERYA (external/erya, SRIMUtilities) has proton coefficients only (SRIM-2013 eight-parameter fits, energy argument
// the proton energy in keV, result MeV cm^2 per atom).  For alphas it has no data; the column "erya" for alphas is
// 4 x ERYA's proton value at the same velocity (the bare-charge scaling the engine would have used with ERYA's
// equation), and is labelled so.  ERYA's fits are of the electronic stopping;
// this module's electronic stopping is compared with them, its nuclear stopping printed beside.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "SRIMUtilities.h"
#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/Stopping.h"

using namespace azure2::stopping;

int main() {
  SRIMUtilities srim;
  if (!srim.isDataLoaded()) {
    std::fprintf(stderr, "SRIM2013.xml not found: run from external/\n");
    return 2;
  }
  const std::vector<std::string> materials = {"H", "He", "C", "N", "O", "Ne", "Al", "Ti", "Ta", "Au", "Ta2O5", "TiN"};
  const std::vector<double> energies = {0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0};
  std::printf("projectile\tmaterial\tE_MeV\tthis_electronic\tthis_nuclear\terya\trel_diff_electronic\terya_note\n");
  for (const Projectile &p : {Proton(), Alpha()}) {
    for (const std::string &formula : materials) {
      const Material m = MaterialFromFormula(formula, Phase::kSolid, 0.0);
      std::vector<CompoundElement> elements;
      for (const auto &c : m.constituents) elements.emplace_back(c.z, c.atom_fraction);
      for (double e : energies) {
        const StoppingValue v = Stopping(m, p, e);
        // ERYA: proton energy in keV at the same velocity; MeV cm^2 per atom -> eV/(1e15 atoms/cm^2)
        const double e_p_keV = e * Proton().mass_u / p.mass_u * 1e3;
        double erya = srim.calculateCompoundStoppingPower(e_p_keV, elements, 0) / kStoppingToMeVcm2;
        if (p.z == 2) erya *= 4.0;
        std::printf("%s\t%s\t%g\t%.5g\t%.4g\t%.5g\t%+.4f\t%s\n", p.label.c_str(), formula.c_str(), e, v.electronic, v.nuclear,
                    erya, v.electronic / erya - 1.0, p.z == 2 ? "4 x ERYA(p, same velocity)" : "ERYA(p)");
      }
    }
  }
  return 0;
}
