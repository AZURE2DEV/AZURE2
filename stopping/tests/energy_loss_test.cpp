// What a target needs (docs/upgrade/stopping.md section 8): N(E_in, E_out), E_out(E_in, N), the range and their
// derivatives with respect to a stopping scale, against closed forms for epsilon = c E^p, and on a real material
// against each other and central differences; the integrals across SRIM's branch boundaries (break points); the
// energy after a layer against LunaSRIM's step integration (lunasrim.straggling.energy_loss).
#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/stopping/EnergyLoss.h"
#include "azure2/stopping/Nuclear.h"
#include "azure2/stopping/Stopping.h"
#include "check.h"

using namespace azure2::stopping;

namespace {
/// The energy loss against LunaSRIM's integration: its log-log interpolation of S on 200 points per decade and its
/// numerical steps limit it to a few 1e-6 (largest found 4.6e-6; stopping.md section 9).
constexpr double kLayerTolerance = 1e-5;
}  // namespace

int main(int argc, char **argv) {
  const std::string ref = argc > 1 ? argv[1] : "srim_reference.csv";
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

  // A real material (with its kinks at 2 keV/u, 1 and 1.2 MeV/u and the nuclear jumps): E_out and N are inverse to
  // each other, and the stopping-scale derivative of E_out matches a central difference.
  const Material m = MaterialFromFormula("Ta2O5");
  for (const Projectile &p : {Proton(), Alpha()}) {
    const double n = 0.5 * Range(m, p, 4.0);
    const double e_out = EnergyAfter(m, p, 4.0, n);
    check::Rel(p.label + " in Ta2O5: N(4 MeV, E_out(4 MeV, N)) = N", ArealDensityBetween(m, p, 4.0, e_out), n, 1e-9);
    const StoppingFunction eps = MaterialStopping(m, p);
    const auto breaks = StoppingBreakpoints(m, p);
    const double h = 1e-5;
    const StoppingFunction up = [&](double e) { return (1 + h) * eps(e); };
    const StoppingFunction down = [&](double e) { return (1 - h) * eps(e); };
    check::Rel(p.label + " in Ta2O5: dE_out/ds", EnergyAfterScaleDerivative(eps, 4.0, n, breaks),
               (EnergyAfter(up, 4.0, n, breaks) - EnergyAfter(down, 4.0, n, breaks)) / (2 * h), 1e-5);
  }

  // Break points: every branch boundary of the stopping is one, and the integral across each equals the sum of the
  // integrals on its two sides (each to 1e-10).  The nuclear stopping jumps at reduced energy 30 (Nuclear.h; the total
  // stopping by up to 1.2 %, alpha in H), the electronic stopping has kinks: the adaptive quadrature also converges
  // across them without the break points, to the same value (measured: within 5e-12).
  const Material hydrogen = MaterialFromFormula("H");
  check::Rel("alpha in H: the total stopping jumps by -1.2 % at the nuclear form change",
             Stopping(hydrogen, Alpha(), NuclearFormChange_MeV(Alpha(), 1, ElementData(1).atomic_weight_u) * (1 + 1e-12)).total /
                     Stopping(hydrogen, Alpha(), NuclearFormChange_MeV(Alpha(), 1, ElementData(1).atomic_weight_u) * (1 - 1e-12)).total - 1,
             -0.0123, 0.02);
  for (const Projectile &p : {Proton(), Alpha(), Ion(6, 12.0, "12C")}) {
    const auto breaks = StoppingBreakpoints(m, p);
    const double per_u = p.mass_u * 1e-3;
    bool all = true;
    for (double b_keV_u : {2.0, p.z == 1 ? 1000.0 : 2000.0, p.z == 1 ? 1200.0 : 2400.0}) {
      bool found = false;
      for (double b : breaks) found = found || std::fabs(b - b_keV_u * per_u) <= 1e-12 * b;
      all = all && found;
    }
    for (int z : {73, 8}) {
      bool found = false;
      const double b30 = NuclearFormChange_MeV(p, z, ElementData(z).atomic_weight_u);
      for (double b : breaks) found = found || b == b30;
      all = all && found;
    }
    check::True(p.label + " in Ta2O5: break points at 2 keV/u, the end of the fit, 1.2 x it, and eps = 30 of Ta and O", all);
    const StoppingFunction eps = MaterialStopping(m, p);
    for (double b : breaks) {
      const double whole = ArealDensityBetween(eps, 1.5 * b, b / 1.5, breaks);
      const double sides = ArealDensityBetween(eps, 1.5 * b, b) + ArealDensityBetween(eps, b, b / 1.5);
      check::Rel(p.label + " in Ta2O5: N across the break at " + std::to_string(b * 1e3) + " keV = sum of its two sides",
                 whole, sides, 1e-9);
      check::Rel(p.label + " in Ta2O5: N across the break at " + std::to_string(b * 1e3) + " keV without break points",
                 ArealDensityBetween(eps, 1.5 * b, b / 1.5), whole, 1e-9);
    }
  }

  // The energy after a layer against LunaSRIM (RK4 in depth on a log-log interpolated table, 2000 steps; this module:
  // adaptive quadrature to 1e-10 and a root).  Compared: the energy loss E_in - E_out.
  for (const auto &r : check::Rows(ref, "layer")) {
    Material layer = MaterialFromFormula(r[0]);
    layer.phase = r[1] == "1" ? Phase::kGas : Phase::kSolid;
    const Projectile p = Ion(std::stoi(r[2]), std::stod(r[3]), "ion");
    const double e_in = std::stod(r[4]) * 1e-3, n = std::stod(r[5]), e_out = std::stod(r[6]) * 1e-3;
    check::Rel("energy loss of Z1=" + r[2] + " (" + r[4] + " keV) in " + r[5] + " atoms/cm^2 " + r[0] + " against LunaSRIM",
               e_in - EnergyAfter(layer, p, e_in, n), e_in - e_out, kLayerTolerance);
  }
  return check::Done();
}
