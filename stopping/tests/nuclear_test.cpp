// SRIM-2013's nuclear stopping (docs/upgrade/stopping.md section 5) against LunaSRIM's own values (SNUC03.dat
// pair-specific fit and the high-energy limit, on each side of the change of form at reduced energy 30), the
// conversion to eV/(1e15 atoms/cm^2) against the definition s_n = eps S_n / (pi a^2 gamma E), and the jump at 30.
#include <cmath>
#include <string>

#include "azure2/physics/Constants.h"
#include "azure2/stopping/Nuclear.h"
#include "check.h"

using namespace azure2::stopping;
namespace phys = azure2::physics;

int main(int argc, char **argv) {
  const std::string ref = argc > 1 ? argv[1] : "srim_reference.csv";
  check::Largest largest;
  int n = 0;
  for (const auto &r : check::Rows(ref, "sn")) {
    const int z1 = std::stoi(r[0]), z2 = std::stoi(r[2]);
    const double m1 = std::stod(r[1]), m2 = std::stod(r[3]), e_keV = std::stod(r[4]), want = std::stod(r[5]);
    const double got = NuclearStopping(Ion(z1, m1, "ion"), z2, m2, e_keV * 1e-3);
    check::RelQuiet("S_n Z1=" + r[0] + " Z2=" + r[2] + " E=" + r[4] + " keV", got, want, 1e-10);
    largest.Track("Z1=" + r[0] + " M1=" + r[1], got, want);
    ++n;
  }
  check::True("nuclear stopping: " + std::to_string(n) + " values against LunaSRIM to 1e-10", n > 1000);
  largest.Summary("nuclear stopping against LunaSRIM, per projectile");

  // The definition: eps = a M2 E / (Z1 Z2 e^2 (M1 + M2)), S_n = pi a^2 gamma E s_n / eps, a = 0.8854 a0 /
  // (Z1^0.23 + Z2^0.23), gamma = 4 M1 M2 / (M1 + M2)^2.  SRIM's 32.53 and 8.462 are these constants rounded.
  const double a0_cm = phys::kHbarC_MeV_fm / (phys::kFineStructure * phys::kElectronMass_MeV) * 1e-13;
  const double e2_eV_cm = phys::kFineStructure * phys::kHbarC_MeV_fm * 1e-13 * 1e6;
  const SrimData &d = SrimTables();
  struct Case { Projectile p; int z2; double e_MeV; };
  for (const Case &c : {Case{Proton(), 13, 0.010}, Case{Alpha(), 79, 0.1}, Case{Ion(6, 12.0, "12C"), 8, 2.0}}) {
    const double z1 = c.p.z, m1 = c.p.mass_u, m2 = ElementData(c.z2).atomic_weight_u;
    const double a = 0.8854 * a0_cm / (std::pow(z1, 0.23) + std::pow(c.z2, 0.23));
    const double e_eV = c.e_MeV * 1e6;
    const double eps = a * m2 * e_eV / (z1 * c.z2 * e2_eV_cm * (m1 + m2));
    const double gamma = 4 * m1 * m2 / ((m1 + m2) * (m1 + m2));
    const double sn = SrimReducedNuclearStopping(NuclearCoefficients(d, c.p.z, c.z2), eps);
    const double sn_definition = 3.14159265358979 * a * a * gamma * e_eV * sn / eps * 1e15;
    check::Rel("reduced energy of " + c.p.label + " on Z2=" + std::to_string(c.z2) + " from its definition",
               ZblReducedEnergy(c.p.z, m1, c.z2, m2, c.e_MeV * 1e3), eps, 2e-3);
    check::Rel("S_n of " + c.p.label + " on Z2=" + std::to_string(c.z2) + " from the definition of s_n",
               NuclearStopping(c.p, c.z2, m2, c.e_MeV), sn_definition, 2e-3);
  }
  // The change of form at reduced energy 30: NuclearFormChange_MeV is where eps = 30, and S_n jumps there (by +18 %
  // for H on H: the fit is above the limit ln(eps)/(2 eps)); the integrals split there (Stopping.h).
  const double m_h = ElementData(1).atomic_weight_u;
  const double e30 = NuclearFormChange_MeV(Proton(), 1, m_h);
  check::Rel("eps = 30 at NuclearFormChange_MeV", ZblReducedEnergy(1, Proton().mass_u, 1, m_h, e30 * 1e3), 30.0, 1e-14);
  check::Rel("S_n of p on H jumps by 17.9 % at eps = 30",
             NuclearStopping(Proton(), 1, m_h, e30 * (1 - 1e-9)) / NuclearStopping(Proton(), 1, m_h, e30 * (1 + 1e-9)) - 1,
             0.1786, 1e-3);
  check::Throws<std::domain_error>("S_n at E = 0 refused", [] { NuclearStopping(Proton(), 13, 26.98, 0.0); });
  return check::Done();
}
