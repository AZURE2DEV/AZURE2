// SRIM-2013's electronic stopping (docs/upgrade/stopping.md section 4) against LunaSRIM's own values
// (srim_reference.csv: p, d, 3He, alpha, 7Li, 12C, 14N, 16O in 25 elements, gases in both phases, every branch with
// each boundary and 1e-3 on each side), the phase factor rules, continuity at the branch boundaries, velocity
// scaling and the stated uncertainties.
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/stopping/Electronic.h"
#include "check.h"

using namespace azure2::stopping;

int main(int argc, char **argv) {
  const std::string ref = argc > 1 ? argv[1] : "srim_reference.csv";
  check::Largest largest;
  int n = 0;
  for (const auto &r : check::Rows(ref, "se")) {
    const int z1 = std::stoi(r[0]), z2 = std::stoi(r[2]);
    const double m1 = std::stod(r[1]), e_keV = std::stod(r[4]), want = std::stod(r[5]);
    const Phase phase = r[3] == "1" ? Phase::kGas : Phase::kSolid;
    const double got = ElectronicStopping(Ion(z1, m1, "ion"), z2, e_keV * 1e-3, phase, 1.0).value;
    check::RelQuiet("S_e Z1=" + r[0] + " M1=" + r[1] + " Z2=" + r[2] + " phase " + r[3] + " E=" + r[4] + " keV", got, want, 1e-10);
    largest.Track("Z1=" + r[0] + " M1=" + r[1], got, want);
    ++n;
  }
  check::True("electronic stopping: " + std::to_string(n) + " values against LunaSRIM to 1e-10", n > 5000);
  largest.Summary("electronic stopping against LunaSRIM, per projectile");

  // The phase factor: 1 for noble gases and an element in its natural phase, 0.92 for a natural gas as a solid,
  // 1/0.92 for a natural solid as a gas.
  check::Rel("phase factor: Ne as a solid is 1 (noble gas)", SrimPhaseFactor(10, Phase::kSolid), 1.0, 0.0);
  check::Rel("phase factor: O as a gas is 1", SrimPhaseFactor(8, Phase::kGas), 1.0, 0.0);
  check::Rel("phase factor: O as a solid is 0.92", SrimPhaseFactor(8, Phase::kSolid), 0.92, 0.0);
  check::Rel("phase factor: C as a gas is 1/0.92", SrimPhaseFactor(6, Phase::kGas), 1.0 / 0.92, 0.0);
  check::Rel("fade-out: a factor 1 stays 1", SrimFadeOut(1.0, 100.0, 1), 1.0, 0.0);
  check::Rel("fade-out: 1 - 0.08 / (1 + exp(1.48 (sqrt(1/25) - 7))) at 1 keV/u (H)", SrimFadeOut(0.92, 1.0, 1),
             1 - 0.08 / (1 + std::exp(1.48 * (0.2 - 7))), 1e-15);
  check::Rel("fade-out: nearly gone at 10 MeV/u (H)", SrimFadeOut(0.92, 1e4, 1), 1.0, 1e-9);
  check::Rel("fade-out: half at sqrt(e/25) = 4 for heavy ions (400 keV/u)", SrimFadeOut(0.9, 400.0, 6), 0.95, 1e-14);

  // Continuous at every branch boundary (the formulas join with kinks, not jumps).
  const SrimData &d = SrimTables();
  const int failures_before = check::failures;
  for (int z1 : {1, 2, 3, 6}) {
    for (int z2 : {1, 6, 73}) {
      for (double b : {2.0, 1000.0, 1200.0, 2000.0, 2400.0, 30000.0}) {
        const double lo = SrimElementStopping(d, z1, z2, b * (1 - 1e-12), Phase::kSolid);
        const double hi = SrimElementStopping(d, z1, z2, b * (1 + 1e-12), Phase::kSolid);
        check::RelQuiet("continuous at " + std::to_string(b) + " keV/u, Z1=" + std::to_string(z1) + " Z2=" + std::to_string(z2),
                        hi, lo, 1e-9);
      }
    }
  }
  check::True("continuous at 2, 1000, 1200, 2000, 2400, 30000 keV/u (Z1 = 1, 2, 3, 6 in H, C, Ta)", check::failures == failures_before);
  // Velocity-proportional below 2 keV/u; Bethe alone above 1.2 x the end of the fit.
  check::Rel("p in Al: S(0.5 keV/u) = S(2 keV/u) / 2", SrimElementStopping(d, 1, 13, 0.5, Phase::kSolid),
             0.5 * SrimElementStopping(d, 1, 13, 2.0, Phase::kSolid), 1e-14);
  check::Rel("p in Al at 2 MeV/u is SRIM's Bethe", SrimElementStopping(d, 1, 13, 2000.0, Phase::kSolid),
             SrimBethe(ElementData(13), 1, 2000.0, Phase::kSolid), 0.0);

  // Same velocity, same charge, same stopping: d at E (m_d/m_p) equals p at E.
  for (double e : {0.001, 0.05, 1.1, 6.0}) {
    check::Rel("d at the velocity of a " + std::to_string(e) + " MeV proton", ElectronicStopping(Deuteron(), 26, e * Deuteron().mass_u / Proton().mass_u, Phase::kSolid, 1.0).value,
               ElectronicStopping(Proton(), 26, e, Phase::kSolid, 1.0).value, 1e-12);
  }
  // Uncertainty: SRIM's published accuracy per ion class.
  check::Rel("uncertainty of H ions 4.0 %", ElectronicStopping(Proton(), 13, 0.5, Phase::kSolid, 1.0).rel_uncertainty, 0.040, 0.0);
  check::Rel("uncertainty of He ions 3.9 %", ElectronicStopping(Alpha(), 13, 0.5, Phase::kSolid, 1.0).rel_uncertainty, 0.039, 0.0);
  check::Rel("uncertainty of Li ions 4.8 %", ElectronicRelUncertainty(3), 0.048, 0.0);
  check::Rel("uncertainty of Be-U ions 5.8 %", ElectronicRelUncertainty(6), 0.058, 0.0);
  check::Throws<std::domain_error>("E = 0 refused", [] { ElectronicStopping(Proton(), 13, 0.0, Phase::kSolid, 1.0); });
  check::Throws<std::domain_error>("E above 1e7 keV/u refused", [] { ElectronicStopping(Proton(), 13, 2e4, Phase::kSolid, 1.0); });
  check::Throws<std::invalid_argument>("a compound correction of 0 refused", [] { ElectronicStopping(Proton(), 13, 1.0, Phase::kSolid, 0.0); });
  check::Throws<std::out_of_range>("Z2 = 93 refused", [] { ElectronicStopping(Proton(), 93, 1.0, Phase::kSolid, 1.0); });
  return check::Done();
}
