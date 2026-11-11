// ZBL nuclear stopping (docs/upgrade/stopping.md section 5): the analytic fit against the universal potential's own
// nuclear stopping (classical scattering integral, make_reference.py), against an independent evaluation of the fit,
// and the conversion to eV/(1e15 atoms/cm^2) against the definition s_n = epsilon S_n / (pi a^2 gamma E).
#include <cmath>
#include <string>

#include "azure2/physics/Constants.h"
#include "azure2/stopping/Nuclear.h"
#include "azure2/stopping/Tables.h"
#include "check.h"

using namespace azure2::stopping;
namespace phys = azure2::physics;

int main(int argc, char **argv) {
  const std::string ref = argc > 1 ? argv[1] : "reference.csv";
  for (const auto &r : check::Reference(ref, "zbl_fit")) {
    check::Rel("ZBL fit s_n(" + std::to_string(r[0]) + ") against an independent evaluation", ZblReducedStopping(r[0]), r[1], 1e-12);
  }
  // The fit represents the universal potential to 1.1 % up to epsilon = 3; at epsilon = 30, where it joins the
  // high-energy form, 4.3 % below the integral (a property of the published fit, stated in stopping.md section 5).
  for (const auto &r : check::Reference(ref, "zbl_universal")) {
    const double tol = r[0] <= 3.0 ? 0.012 : 0.05;
    check::Rel("ZBL fit s_n(" + std::to_string(r[0]) + ") against the universal potential's scattering integral",
               ZblReducedStopping(r[0]), r[1], tol);
  }
  // The definition: epsilon = a M2 E_lab / (Z1 Z2 e^2 (M1 + M2)), S_n = pi a^2 gamma E_lab s_n / epsilon,
  // a = 0.8854 a0 / (Z1^0.23 + Z2^0.23), gamma = 4 M1 M2 / (M1 + M2)^2.  The fit's 32.53 and 8.462 are these
  // constants rounded to four digits.
  const double a0_cm = phys::kHbarC_MeV_fm / (phys::kFineStructure * phys::kElectronMass_MeV) * 1e-13;
  const double e2_eV_cm = phys::kFineStructure * phys::kHbarC_MeV_fm * 1e-13 * 1e6;
  struct Case { Projectile p; int z2; double m2; double e_MeV; };
  for (const Case &c : {Case{Proton(), 13, Element(13).atomic_weight_u, 0.010},
                        Case{Alpha(), 79, Element(79).atomic_weight_u, 0.1},
                        Case{Ion(6, 12.0, "12C"), 8, Element(8).atomic_weight_u, 2.0}}) {
    const double z1 = c.p.z, m1 = c.p.mass_u;
    const double a = 0.8854 * a0_cm / (std::pow(z1, 0.23) + std::pow(c.z2, 0.23));
    const double e_eV = c.e_MeV * 1e6;
    const double eps = a * c.m2 * e_eV / (z1 * c.z2 * e2_eV_cm * (m1 + c.m2));
    const double gamma = 4 * m1 * c.m2 / ((m1 + c.m2) * (m1 + c.m2));
    const double sn_definition = 3.14159265358979 * a * a * gamma * e_eV * ZblReducedStopping(eps) / eps * 1e15;
    check::Rel("reduced energy of " + c.p.label + " on Z2=" + std::to_string(c.z2) + " from its definition",
               ZblReducedEnergy(c.p.z, m1, c.z2, c.m2, c.e_MeV * 1e3), eps, 2e-3);
    check::Rel("S_n of " + c.p.label + " on Z2=" + std::to_string(c.z2) + " from the definition of s_n",
               NuclearStopping(c.p, c.z2, c.m2, c.e_MeV), sn_definition, 2e-3);
  }
  check::Throws<std::domain_error>("S_n at E = 0 refused", [] { NuclearStopping(Proton(), 13, 26.98, 0.0); });
  return check::Done();
}
