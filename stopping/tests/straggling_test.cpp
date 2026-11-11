// Straggling (docs/upgrade/stopping.md section 10): Bohr's value, TRIM's model (mean target charge, relativistic
// factor, low-energy limit), the thin-layer and constant-stopping limits of Tschalaer's propagation, the propagation
// against a closed form for epsilon = c E^p, the beam spread, and both models after a layer against LunaSRIM.
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/stopping/Stopping.h"
#include "azure2/stopping/Straggling.h"
#include "check.h"

using namespace azure2::stopping;

namespace {
/// The straggling after a layer against LunaSRIM's step integration of dOmega^2/dx = W - 2 S' Omega^2 (its S' is a
/// numerical difference of an interpolated table); the largest difference found is 3e-6 (stopping.md section 9).
constexpr double kLayerStragglingTolerance = 1e-5;
}  // namespace

int main(int argc, char **argv) {
  const std::string ref = argc > 1 ? argv[1] : "srim_reference.csv";
  // Bohr: Omega^2 [keV^2] = 0.26 Z1^2 Z2 N [1e18 atoms/cm^2] (the rounded form quoted e.g. in the SIMNRA guide).
  check::Rel("Bohr: 0.26 keV^2 per 1e18 atoms/cm^2 for Z1 = Z2 = 1", BohrSpecificVariance(1, 1.0) * 1e18 * 1e6, 0.26, 3e-3);
  check::Rel("Bohr scales as Z1^2 Z2", BohrSpecificVariance(2, 13.0), 4 * 13 * BohrSpecificVariance(1, 1.0), 1e-14);
  const Material m = MaterialFromFormula("Ta2O5");
  check::Rel("Bohr of a compound uses the mean Z2 = (2x73 + 5x8)/7", MeanTargetCharge(m), (2 * 73 + 5 * 8) / 7.0, 1e-15);
  check::Rel("Bohr model of Ta2O5 for alphas", SpecificStragglingVariance(StragglingModel::kBohr, m, Alpha(), 2.0),
             BohrSpecificVariance(2, MeanTargetCharge(m)), 1e-14);

  // TRIM: the integer mean Z (round half to even after each element), Bohr with 177^2/12 and the relativistic factor,
  // the low-energy limit k E.
  check::True("TRIM mean Z of Ta2O5: round(2/7 73) = 21, round(21 + 5/7 8) = 27", TrimMeanTargetCharge(m) == 27);
  check::True("TRIM mean Z of TiN: 11, then round(11 + 3.5) = 14 (half to even)", TrimMeanTargetCharge(MaterialFromFormula("TiN")) == 14);
  const Material al = MaterialFromFormula("Al");
  const double e = 2.0;  // MeV
  const double q = 1 - 1 / std::pow(1 + e * 1e6 / (Proton().mass_u * 1000 * 931494), 2);
  const double w = 177.0 * std::sqrt(13.0) * (1 / q - 0.5) / (1 / q - 1);
  check::Rel("TRIM for p in Al at 2 MeV: (177 Z2^1/2 relativistic factor)^2 / 12, eV^2 A^2 -> MeV^2 cm^2",
             TrimSpecificVariance(al, Proton(), e), w * w / 12 * 1e-28, 1e-14);
  const double kPi = 3.14159265358979323846;
  check::Rel("TRIM / (Bohr x relativistic factor^2) = 177^2/12 / (4 pi e^4), e^2 = 14.39964 eV A (p in Al, 2 MeV)",
             TrimSpecificVariance(al, Proton(), e) / (BohrSpecificVariance(1, 13.0) * std::pow((1 / q - 0.5) / (1 / q - 1), 2)),
             177.0 * 177.0 / 12 / (4 * kPi * 14.39964 * 14.39964), 1e-6);
  const double k = std::sqrt(3.2 * 12 / (Proton().mass_u * 1000) * std::pow(14.0, 8.0 / 3.0));
  check::Rel("TRIM below the Bohr value: w = k E (p in Al at 5 eV)", TrimSpecificVariance(al, Proton(), 5e-6),
             std::pow(k * 5.0, 2) / 12 * 1e-28, 1e-14);
  check::True("model names: trim, bohr", StragglingModelFromName("trim") == StragglingModel::kTrim &&
                                             StragglingModelFromName("bohr") == StragglingModel::kBohr);
  check::Throws<std::invalid_argument>("an unknown model (chu) refused", [] { StragglingModelFromName("chu"); });

  // Constant stopping: no propagation factor; sigma^2 = k N, the beam spread unchanged and added in quadrature.
  const StoppingFunction flat = [](double) { return 40.0; };
  const double kb = BohrSpecificVariance(1, 13.0);
  const auto wf = PropagateStraggling(flat, [kb](double) { return kb; }, 2.0, 5e18, 0.003);
  check::Rel("constant S: straggling^2 = dOmega^2/dN N", wf.straggling * wf.straggling, kb * 5e18, 1e-9);
  check::Rel("constant S: beam spread unchanged", wf.beam, 0.003, 1e-14);
  check::Rel("constant S: total in quadrature", wf.total, std::hypot(std::sqrt(kb * 5e18), 0.003), 1e-9);

  // Thin layer on a real material: sigma^2 = dOmega^2/dN N.
  const auto thin = PropagateStraggling(StragglingModel::kTrim, m, Proton(), 1.5, 1e16, 0.0);
  check::Rel("thin layer (1e16 atoms/cm^2): TRIM", thin.straggling * thin.straggling,
             SpecificStragglingVariance(StragglingModel::kTrim, m, Proton(), 1.5) * 1e16, 1e-3);

  // epsilon = c E^p with S = c 1e-21 E^p (MeV cm^2), dOmega^2/dN = k:
  //   sigma^2 = S(E_out)^2 k (E_in^(1-3p) - E_out^(1-3p)) / ((1 - 3p) (c 1e-21)^3),  beam x (E_out/E_in)^p.
  for (double p : {0.5, -0.7}) {
    const double c = 30.0;
    const StoppingFunction eps = [c, p](double x) { return c * std::pow(x, p); };
    const double ein = 3.0, n = 0.5 * Range(eps, ein);
    const auto r = PropagateStraggling(eps, [kb](double) { return kb; }, ein, n, 0.01);
    const double eout = EnergyAfter(eps, ein, n);
    const double s = c * 1e-21;
    const double sigma2 = std::pow(s * std::pow(eout, p), 2) * kb *
                          (std::pow(ein, 1 - 3 * p) - std::pow(eout, 1 - 3 * p)) / ((1 - 3 * p) * s * s * s);
    check::Rel("Tschalaer against the closed form, p = " + std::to_string(p), r.straggling * r.straggling, sigma2, 1e-8);
    check::Rel("beam spread carried as S(E_out)/S(E_in), p = " + std::to_string(p), r.beam, 0.01 * std::pow(eout / ein, p), 1e-12);
  }

  // Both models after a layer against LunaSRIM.
  for (const auto &r : check::Rows(ref, "layer")) {
    Material layer = MaterialFromFormula(r[0]);
    layer.phase = r[1] == "1" ? Phase::kGas : Phase::kSolid;
    const Projectile p = Ion(std::stoi(r[2]), std::stod(r[3]), "ion");
    const double e_in = std::stod(r[4]) * 1e-3, n = std::stod(r[5]);
    const std::string tag = "Z1=" + r[2] + " (" + r[4] + " keV) in " + r[5] + " atoms/cm^2 " + r[0] + " against LunaSRIM";
    check::Rel("TRIM straggling of " + tag, PropagateStraggling(StragglingModel::kTrim, layer, p, e_in, n, 0.0).straggling,
               std::stod(r[7]) * 1e-3, kLayerStragglingTolerance);
    check::Rel("Bohr straggling of " + tag, PropagateStraggling(StragglingModel::kBohr, layer, p, e_in, n, 0.0).straggling,
               std::stod(r[8]) * 1e-3, kLayerStragglingTolerance);
  }
  return check::Done();
}
