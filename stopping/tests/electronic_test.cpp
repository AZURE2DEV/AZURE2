// Electronic stopping (docs/upgrade/stopping.md section 4): the fit forms and the Bethe terms against an independent
// evaluation (make_reference.py; Bloch's correction against the digamma function), the velocity-proportional limit,
// continuity at the two joins, velocity scaling for d, t, 3He and effective-charge scaling for heavier ions.
#include <cmath>
#include <string>

#include "azure2/physics/Constants.h"
#include "azure2/stopping/Electronic.h"
#include "azure2/stopping/Tables.h"
#include "check.h"

using namespace azure2::stopping;

int main(int argc, char **argv) {
  const std::string ref = argc > 1 ? argv[1] : "reference.csv";
  for (const auto &r : check::Reference(ref, "fit_p")) {
    check::Rel("AZ_p fit Z2=" + std::to_string(int(r[0])) + " T=" + std::to_string(r[1]) + " keV",
               ProtonFitStopping(Coefficients(int(r[0]), 1), r[1]), r[2], 1e-12);
  }
  for (const auto &r : check::Reference(ref, "fit_he")) {
    check::Rel("Z_He fit Z2=" + std::to_string(int(r[0])) + " T=" + std::to_string(r[1]) + " MeV",
               HeliumFitStopping(Coefficients(int(r[0]), 2), r[1]), r[2], 1e-12);
  }
  for (const auto &r : check::Reference(ref, "bloch")) {
    check::Rel("Bloch L2(y=" + std::to_string(r[0]) + ") = psi(1) - Re psi(1+iy)", BlochCorrection(r[0]), r[1], 1e-9);
  }
  for (const auto &r : check::Reference(ref, "shell")) {
    check::Rel("Barkas-Berger C(I=" + std::to_string(r[0]) + ", eta=" + std::to_string(r[1]) + ")",
               ShellCorrectionBarkasBerger(r[0], r[1]), r[2], 1e-12);
  }
  for (const auto &r : check::Reference(ref, "bethe")) {
    check::Rel("Bethe z1=" + std::to_string(int(r[0])) + " Z2=" + std::to_string(int(r[1])) + " beta2=" + std::to_string(r[3]),
               BetheStopping(int(r[0]), int(r[1]), r[2], r[3]), r[4], 1e-9);
  }

  // Velocity-proportional below 10 keV (p) and 1 keV (He): S(T) = S(T0) sqrt(T/T0).
  check::Rel("AZ_p: S(2.5 keV) = S(10 keV)/2", ProtonFitStopping(Coefficients(13, 1), 2.5),
             0.5 * ProtonFitStopping(Coefficients(13, 1), 10.0), 1e-14);
  check::Rel("Z_He: S(0.25 keV) = S(1 keV)/2", HeliumFitStopping(Coefficients(79, 2), 2.5e-4),
             0.5 * HeliumFitStopping(Coefficients(79, 2), 1e-3), 1e-14);

  // The joined evaluation is continuous at the end of the fit (1 MeV/u) and at the start of Bethe (3 MeV/u).
  for (const Projectile &p : {Proton(), Alpha()}) {
    for (int z2 : {1, 6, 13, 79}) {
      const double t1 = Coefficients(z2, p.z).fit_max_MeV_u * p.mass_u;
      const double t2 = kBetheFrom_MeV_u * p.mass_u;
      const std::string tag = p.label + " in Z2=" + std::to_string(z2);
      check::Rel("continuous at the end of the fit, " + tag, ElectronicStopping(p, z2, t1 * (1 + 1e-9)).value,
                 ElectronicStopping(p, z2, t1 * (1 - 1e-9)).value, 1e-7);
      check::Rel("continuous at the start of Bethe, " + tag, ElectronicStopping(p, z2, t2 * (1 + 1e-9)).value,
                 ElectronicStopping(p, z2, t2 * (1 - 1e-9)).value, 1e-7);
    }
  }
  // Below the end of the fit the evaluation is the fit, above 3 MeV/u it is Bethe.
  const Projectile p = Proton();
  check::Rel("p in Al at 0.5 MeV is the fit", ElectronicStopping(p, 13, 0.5).value,
             ProtonFitStopping(Coefficients(13, 1), 0.5 / p.mass_u * 1e3), 1e-12);
  check::Rel("p in Al at 8 MeV is Bethe", ElectronicStopping(p, 13, 8.0).value,
             BetheStopping(1, 13, Element(13).mean_excitation_eV, Beta2(p, 8.0)), 1e-12);
  check::Rel("uncertainty of the fit region is the row's", ElectronicStopping(p, 13, 0.5).rel_uncertainty,
             Coefficients(13, 1).rel_uncertainty, 1e-15);
  check::Rel("uncertainty of the Bethe region", ElectronicStopping(p, 13, 8.0).rel_uncertainty, kBetheRelUncertainty, 1e-15);

  // Same velocity, same charge: d at E (m_d/m_p) and t at E (m_t/m_p) equal p at E; 3He equals alpha.
  for (double e : {0.05, 0.8, 2.0, 6.0}) {
    const double sp = ElectronicStopping(Proton(), 26, e).value;
    check::Rel("d at the same velocity as p (" + std::to_string(e) + " MeV p)",
               ElectronicStopping(Deuteron(), 26, e * Deuteron().mass_u / Proton().mass_u).value, sp, 1e-12);
    check::Rel("t at the same velocity as p (" + std::to_string(e) + " MeV p)",
               ElectronicStopping(Triton(), 26, e * Triton().mass_u / Proton().mass_u).value, sp, 1e-12);
    const double sa = ElectronicStopping(Alpha(), 26, 4 * e).value;
    check::Rel("3He at the same velocity as alpha (" + std::to_string(4 * e) + " MeV alpha)",
               ElectronicStopping(Helium3(), 26, 4 * e * Helium3().mass_u / Alpha().mass_u).value, sa, 1e-12);
  }
  // Heavier ions: (q Z1)^2 times the proton value at the same velocity.
  const Projectile c12 = Ion(6, 12.0, "12C");
  const double e12 = 12.0;  // MeV: 1 MeV/u
  const double q = EffectiveChargeFraction(6, Beta2(c12, e12));
  const double vr = std::sqrt(Beta2(c12, e12)) / azure2::physics::kFineStructure / std::pow(6.0, 2.0 / 3.0);
  check::Rel("Pierce-Blann q/Z1 = 1 - exp(-0.95 v/(v0 Z1^2/3))", q, 1.0 - std::exp(-0.95 * vr), 1e-14);
  check::Rel("12C at 1 MeV/u: (6q)^2 x p at the same velocity", ElectronicStopping(c12, 13, e12).value,
             36 * q * q * ElectronicStopping(Proton(), 13, e12 * Proton().mass_u / 12.0).value, 1e-12);
  check::Rel("12C uncertainty (v > v0 Z1^2/3)", ElectronicStopping(c12, 13, e12).rel_uncertainty, kHeavyIonRelUncertainty, 1e-15);
  check::Throws<std::domain_error>("E = 0 refused", [] { ElectronicStopping(Proton(), 13, 0.0); });
  return check::Done();
}
