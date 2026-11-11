#ifndef AZURE2_STOPPING_ENERGYLOSS_H
#define AZURE2_STOPPING_ENERGYLOSS_H

#include <functional>
#include <vector>

#include "azure2/stopping/Material.h"

namespace azure2 {
namespace stopping {

/*!
 * What a target needs of a stopping (docs/upgrade/stopping.md section 8).  Energies in MeV (lab), epsilon in
 * eV/(1e15 atoms/cm^2), areal densities in atoms/cm^2.  Every function exists twice: on a Material and projectile
 * (the evaluation of Stopping(), total), and on any function epsilon(E) (so that a closed form tests it, and a
 * caller can pass a scaled or user stopping).
 *
 *   N(E_in, E_out) = int_{E_out}^{E_in} dE / epsilon(E)            the areal density between two energies
 *   E_out(E_in, N): N(E_in, E_out) = N                              the energy after an areal density
 *   R(E) = N(E, 0)                                                   the range (path length, areal density)
 *
 * The integral is computed in u = sqrt(E) (dE = 2u du), which removes the 1/sqrt(E) end point of a
 * velocity-proportional stopping, with adaptive Gauss-Kronrod 7-15 quadrature to the relative tolerance
 * kQuadratureRelTol; E_out is the root of N(E_in, E) - N by safeguarded Newton steps (dN/dE_out = -1/epsilon).
 * The integration interval is split at the energies \p breaks_MeV where epsilon changes form (SRIM's nuclear stopping
 * jumps at reduced energy 30, its electronic stopping has kinks: Stopping.h StoppingBreakpoints); the functions on a
 * Material pass them.
 *
 * The stopping-scale nuisance of ARCHITECTURE.md section 7 multiplies epsilon by s.  At s = 1:
 *   dN/ds = -N,   dR/ds = -R,   dE_out/ds = -N epsilon(E_out)   (so d(Delta E)/ds = N epsilon(E_out)).
 */
constexpr double kQuadratureRelTol = 1.0e-10;

/// epsilon(E): eV/(1e15 atoms/cm^2) at lab energy E (MeV).
using StoppingFunction = std::function<double(double)>;

/// N(E_in, E_out) for any epsilon (E_in >= E_out >= 0).  Throws std::invalid_argument for E_out > E_in or E_out < 0,
/// std::domain_error if epsilon is not positive inside the interval.
double ArealDensityBetween(const StoppingFunction &epsilon, double energy_in_MeV, double energy_out_MeV,
                           const std::vector<double> &breaks_MeV = {});
/// E_out after areal density \p areal_density (>= 0).  Throws std::domain_error if the projectile stops (N >= R(E_in)).
double EnergyAfter(const StoppingFunction &epsilon, double energy_in_MeV, double areal_density,
                   const std::vector<double> &breaks_MeV = {});
/// R(E) = N(E, 0).
double Range(const StoppingFunction &epsilon, double energy_MeV, const std::vector<double> &breaks_MeV = {});

/// The same on a material: epsilon = Stopping(material, projectile, E).total.
StoppingFunction MaterialStopping(const Material &material, const Projectile &projectile);
double ArealDensityBetween(const Material &material, const Projectile &projectile, double energy_in_MeV,
                           double energy_out_MeV);
double EnergyAfter(const Material &material, const Projectile &projectile, double energy_in_MeV, double areal_density);
double Range(const Material &material, const Projectile &projectile, double energy_MeV);

/// The derivatives with respect to the stopping scale s (epsilon -> s epsilon) at s = 1.
double ArealDensityScaleDerivative(const StoppingFunction &epsilon, double energy_in_MeV, double energy_out_MeV,
                                   const std::vector<double> &breaks_MeV = {});
double EnergyAfterScaleDerivative(const StoppingFunction &epsilon, double energy_in_MeV, double areal_density,
                                  const std::vector<double> &breaks_MeV = {});
double RangeScaleDerivative(const StoppingFunction &epsilon, double energy_MeV, const std::vector<double> &breaks_MeV = {});

/// int_{a}^{b} f(E) dE for E in [a, b] with 0 <= a <= b, in u = sqrt(E), adaptive Gauss-Kronrod 7-15 on each piece
/// between the break points inside (a, b), to the relative tolerance \p rel_tol of the whole (absolute floor 1e-300).
/// The quadrature every function above uses; exposed for the straggling propagation and its tests.
double IntegrateInSqrtE(const std::function<double(double)> &f, double a_MeV, double b_MeV, double rel_tol,
                        const std::vector<double> &breaks_MeV = {});

}  // namespace stopping
}  // namespace azure2

#endif
