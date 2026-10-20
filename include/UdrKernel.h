#ifndef UDRKERNEL_H
#define UDRKERNEL_H

#include <cmath>
#include <vector>

/*!
 * A user-defined numerical resolution function (SAMMY's "UDR") evaluated
 * for one data point: the probability density R(tau) that a neutron (or any
 * time-of-flight projectile) is detected with a time-of-flight delay tau
 * (microseconds) relative to its true flight time, after convolution with
 * the burst and channel widths and, if requested, re-centred on its
 * centroid.  R is piecewise linear between the tabulated tau values and
 * normalised to unit integral.
 *
 * A positive delay makes the projectile look slower than it is, so the
 * true energies that contribute to a point at nominal energy E0 satisfy
 * t(E') = t(E0) - tau: a delay maps to a HIGHER true energy.
 */
struct UdrKernel {
  std::vector<double> tau;  ///< delay, microseconds, ascending
  std::vector<double> r;    ///< probability density, per microsecond
  double tauMin = 0.0;      ///< support: first tau with non-negligible density
  double tauMax = 0.0;      ///< support: last tau with non-negligible density
};

namespace UdrKinematics {
const double speedOfLightMPerMicroS = 299.792458;  ///< m / microsecond

/// Flight time (microseconds) of a projectile of rest energy mass (MeV) and
/// laboratory kinetic energy e (MeV) over flightPath metres, relativistic.
inline double TimeOfFlight(double e, double mass, double flightPath) {
  double gamma = 1.0 + e / mass;
  double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
  return flightPath / (beta * speedOfLightMPerMicroS);
}

/// Inverse of TimeOfFlight: laboratory kinetic energy (MeV) for a flight
/// time t (microseconds); 0 when t is not a physical flight time.
inline double EnergyFromTime(double t, double mass, double flightPath) {
  if (t <= 0.0) return 0.0;
  double beta = flightPath / (t * speedOfLightMPerMicroS);
  if (beta >= 1.0) return 0.0;
  double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  return mass * (gamma - 1.0);
}
}  // namespace UdrKinematics

#endif
