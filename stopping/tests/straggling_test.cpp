// Straggling (docs/upgrade/stopping.md section 10): Bohr's value, the thin-layer and constant-stopping limits of
// Tschalaer's propagation, the propagation against a closed form for epsilon = c E^p, the beam spread, Chu's table
// and the charge-state term of Yang, O'Connor and Wang with given constants.
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/stopping/Straggling.h"
#include "check.h"

using namespace azure2::stopping;

int main(int, char **) {
  // Bohr: Omega^2 [keV^2] = 0.26 Z1^2 Z2 N [1e18 atoms/cm^2] (the rounded form quoted e.g. in the SIMNRA guide).
  check::Rel("Bohr: 0.26 keV^2 per 1e18 atoms/cm^2 for Z1 = Z2 = 1", BohrSpecificVariance(1, 1.0) * 1e18 * 1e6, 0.26, 3e-3);
  check::Rel("Bohr scales as Z1^2 Z2", BohrSpecificVariance(2, 13.0), 4 * 13 * BohrSpecificVariance(1, 1.0), 1e-14);
  const Material m = MaterialFromFormula("Ta2O5", Phase::kSolid, 8.2);
  check::Rel("Bohr of a compound uses the mean Z2 = (2x73 + 5x8)/7", MeanTargetCharge(m), (2 * 73 + 5 * 8) / 7.0, 1e-15);
  const StragglingModel bohr;
  check::Rel("Bohr model of Ta2O5 for alphas", SpecificStragglingVariance(bohr, m, Alpha(), 2.0),
             BohrSpecificVariance(2, MeanTargetCharge(m)), 1e-14);

  // Constant stopping: no propagation factor; sigma^2 = k N, the beam spread unchanged and added in quadrature.
  const StoppingFunction flat = [](double) { return 40.0; };
  const double k = BohrSpecificVariance(1, 13.0);
  const auto w = PropagateStraggling(flat, [k](double) { return k; }, 2.0, 5e18, 0.003);
  check::Rel("constant S: straggling^2 = dOmega^2/dN N", w.straggling * w.straggling, k * 5e18, 1e-9);
  check::Rel("constant S: beam spread unchanged", w.beam, 0.003, 1e-14);
  check::Rel("constant S: total in quadrature", w.total, std::hypot(std::sqrt(k * 5e18), 0.003), 1e-9);

  // Thin layer on a real material: Bohr's sigma^2 = dOmega^2/dN N.
  const auto thin = PropagateStraggling(bohr, m, Proton(), 1.5, 1e16, 0.0);
  check::Rel("thin layer (1e16 atoms/cm^2): Bohr", thin.straggling * thin.straggling,
             SpecificStragglingVariance(bohr, m, Proton(), 1.5) * 1e16, 1e-3);

  // epsilon = c E^p with S = c 1e-21 E^p (MeV cm^2), dOmega^2/dN = k:
  //   sigma^2 = S(E_out)^2 k (E_in^(1-3p) - E_out^(1-3p)) / ((1 - 3p) (c 1e-21)^3),  beam x (E_out/E_in)^p.
  for (double p : {0.5, -0.7}) {
    const double c = 30.0;
    const StoppingFunction eps = [c, p](double e) { return c * std::pow(e, p); };
    const double ein = 3.0, n = 0.5 * Range(eps, ein);
    const auto r = PropagateStraggling(eps, [k](double) { return k; }, ein, n, 0.01);
    const double eout = EnergyAfter(eps, ein, n);
    const double s = c * 1e-21;
    const double sigma2 = std::pow(s * std::pow(eout, p), 2) * k *
                          (std::pow(ein, 1 - 3 * p) - std::pow(eout, 1 - 3 * p)) / ((1 - 3 * p) * s * s * s);
    check::Rel("Tschalaer against the closed form, p = " + std::to_string(p), r.straggling * r.straggling, sigma2, 1e-8);
    check::Rel("beam spread carried as S(E_out)/S(E_in), p = " + std::to_string(p), r.beam, 0.01 * std::pow(eout / ein, p), 1e-12);
  }

  // Chu's factor from a table, and the Yang charge-state term with given constants.
  ChuTable chu;
  chu.energy_keV_u = {100.0, 500.0, 1000.0};
  chu.z2 = {8, 73};
  chu.h = {{0.6, 0.8, 0.9}, {0.3, 0.5, 0.7}};
  chu.source = "test table";
  check::Rel("Chu factor interpolated linearly", ChuFactor(chu, 73, 300.0), 0.4, 1e-14);
  check::Rel("Chu factor 1 above the table (Bohr's limit)", ChuFactor(chu, 8, 2000.0), 1.0, 0.0);
  check::Throws<std::out_of_range>("Chu factor of a Z2 the table lacks refused", [&] { ChuFactor(chu, 13, 300.0); });
  const StragglingModel chu_model{StragglingKind::kChu, &chu, nullptr};
  const double e = 0.3 * Proton().mass_u;  // 300 keV/u
  check::Rel("Chu model of Ta2O5 for protons", SpecificStragglingVariance(chu_model, m, Proton(), e),
             (2 * 0.4 * BohrSpecificVariance(1, 73) + 5 * 0.7 * BohrSpecificVariance(1, 8)) / 7, 1e-12);
  YangCoefficients yang;
  yang.b = {0.2, 0.3, 0.5, 2.0};
  yang.c = {1.0, 0.1, 0.4, 3.0};
  yang.source = "test constants";
  const double g = 0.5 * (1 - std::exp(-2.0 * 0.3));
  check::Rel("Yang term, hydrogen ions", YangChargeStateTerm(yang, 1, 13, 0.3), 0.2 * g / ((0.3 - 0.3) * (0.3 - 0.3) + g * g), 1e-14);
  const double eps_r = 0.5 / (std::pow(2.0, 1.5) * std::sqrt(13.0));
  const double gh = 0.4 * (1 - std::exp(-3.0 * eps_r));
  check::Rel("Yang term, helium", YangChargeStateTerm(yang, 2, 13, 0.5),
             std::pow(2.0, 4.0 / 3.0) / std::cbrt(13.0) * gh / ((eps_r - 0.1) * (eps_r - 0.1) + gh * gh), 1e-14);
  const StragglingModel yang_model{StragglingKind::kChuYang, &chu, &yang};
  const Material o = MaterialFromFormula("O", Phase::kSolid, 1.0);
  check::Rel("Chu-Yang model for protons: H + Yang term (gamma = 1)", SpecificStragglingVariance(yang_model, o, Proton(), e),
             (ChuFactor(chu, 8, 300.0) + YangChargeStateTerm(yang, 1, 8, 0.3)) * BohrSpecificVariance(1, 8), 1e-12);
  check::Throws<std::invalid_argument>("Chu model without a table refused", [&] {
    SpecificStragglingVariance(StragglingModel{StragglingKind::kChu, nullptr, nullptr}, m, Proton(), 1.0);
  });
  check::Throws<std::invalid_argument>("Yang model without constants refused", [&] {
    SpecificStragglingVariance(StragglingModel{StragglingKind::kChuYang, &chu, nullptr}, m, Proton(), 1.0);
  });
  return check::Done();
}
