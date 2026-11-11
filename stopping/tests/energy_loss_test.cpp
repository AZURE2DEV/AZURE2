// What a target needs (docs/upgrade/stopping.md section 8): N(E_in, E_out), E_out(E_in, N), the range and their
// derivatives with respect to a stopping scale, against closed forms for epsilon = c E^p, and on a real material
// against each other and central differences.
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/stopping/EnergyLoss.h"
#include "azure2/stopping/Stopping.h"
#include "check.h"

using namespace azure2::stopping;

int main(int, char **) {
  // epsilon = c E^p (eV/(1e15 atoms/cm^2), E in MeV): N = 1e21 (E_in^(1-p) - E_out^(1-p)) / (c (1 - p)).
  struct Power { double c, p; };
  for (const Power &w : {Power{20.0, 0.5}, Power{60.0, -0.8}, Power{35.0, 0.0}, Power{10.0, 0.3}}) {
    const StoppingFunction eps = [w](double e) { return w.c * std::pow(e, w.p); };
    const auto n_exact = [w](double ein, double eout) {
      return 1e21 * (std::pow(ein, 1 - w.p) - std::pow(eout, 1 - w.p)) / (w.c * (1 - w.p));
    };
    const std::string tag = "eps = " + std::to_string(w.c) + " E^" + std::to_string(w.p);
    check::Rel("N(2, 0.3) closed form, " + tag, ArealDensityBetween(eps, 2.0, 0.3), n_exact(2.0, 0.3), 1e-9);
    check::Rel("range(2) closed form, " + tag, Range(eps, 2.0), n_exact(2.0, 0.0), 1e-9);
    const double n = 0.6 * n_exact(2.0, 0.0);
    const double e_out = std::pow(std::pow(2.0, 1 - w.p) - n * w.c * (1 - w.p) / 1e21, 1 / (1 - w.p));
    check::Rel("E_out(2, 0.6 R) closed form, " + tag, EnergyAfter(eps, 2.0, n), e_out, 1e-9);
    // scale derivatives against central differences of the scaled stopping
    const double h = 1e-5;
    const StoppingFunction up = [&](double e) { return (1 + h) * eps(e); };
    const StoppingFunction down = [&](double e) { return (1 - h) * eps(e); };
    check::Rel("dN/ds = -N, " + tag, ArealDensityScaleDerivative(eps, 2.0, 0.3),
               (ArealDensityBetween(up, 2.0, 0.3) - ArealDensityBetween(down, 2.0, 0.3)) / (2 * h), 1e-6);
    check::Rel("dR/ds = -R, " + tag, RangeScaleDerivative(eps, 2.0), (Range(up, 2.0) - Range(down, 2.0)) / (2 * h), 1e-6);
    check::Rel("dE_out/ds = -N eps(E_out), " + tag, EnergyAfterScaleDerivative(eps, 2.0, n),
               (EnergyAfter(up, 2.0, n) - EnergyAfter(down, 2.0, n)) / (2 * h), 1e-6);
  }
  const StoppingFunction flat = [](double) { return 50.0; };
  check::Rel("thin layer: Delta E = eps N 1e-21", 1.0 - EnergyAfter(flat, 1.0, 1e18), 50.0 * 1e18 * 1e-21, 1e-12);
  check::Rel("N = 0 leaves the energy", EnergyAfter(flat, 1.0, 0.0), 1.0, 0.0);
  check::Throws<std::domain_error>("a layer thicker than the range refused", [&] { EnergyAfter(flat, 1.0, 3e19); });
  check::Throws<std::invalid_argument>("E_out > E_in refused", [&] { ArealDensityBetween(flat, 1.0, 2.0); });
  check::Throws<std::domain_error>("a non-positive stopping refused",
                                   [] { Range([](double e) { return e - 0.5; }, 1.0); });

  // A real material (with its kinks at 10 keV/u, 1 and 3 MeV/u): E_out and N are inverse to each other, and the
  // stopping-scale derivative of E_out matches a central difference.
  const Material m = MaterialFromFormula("Ta2O5", Phase::kSolid, 8.2);
  for (const Projectile &p : {Proton(), Alpha()}) {
    const double n = 0.5 * Range(m, p, 4.0);
    const double e_out = EnergyAfter(m, p, 4.0, n);
    check::Rel(p.label + " in Ta2O5: N(4 MeV, E_out(4 MeV, N)) = N", ArealDensityBetween(m, p, 4.0, e_out), n, 1e-9);
    const StoppingFunction eps = MaterialStopping(m, p);
    const double h = 1e-5;
    const StoppingFunction up = [&](double e) { return (1 + h) * eps(e); };
    const StoppingFunction down = [&](double e) { return (1 - h) * eps(e); };
    check::Rel(p.label + " in Ta2O5: dE_out/ds", EnergyAfterScaleDerivative(eps, 4.0, n),
               (EnergyAfter(up, 4.0, n) - EnergyAfter(down, 4.0, n)) / (2 * h), 1e-5);
  }
  return check::Done();
}
