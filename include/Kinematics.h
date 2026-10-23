#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <cmath>

#include "Config.h"
#include "Constants.h"
#include "PPair.h"

/*!
 * Two-body kinematics of a particle pair: the channel wave number k and Sommerfeld
 * parameter eta at a centre-of-mass kinetic energy T, and the lab <-> c.m. energy maps.
 *
 * AZURE2 works per channel with T_c = E_x - S_c (the excitation energy less the
 * channel's separation energy).  Non-relativistically T_c is m2 E_lab / m and
 * k^2 = 2 mu T_c.  With relativistic kinematics ("relativistic=1" in the <config>
 * block) T_c is the invariant kinetic energy sqrt(s) - m1 - m2 and, with m = m1 + m2,
 *
 *   s = (m + T)^2 = m^2 + 2 m2 E_lab                      (target, particle 2, at rest)
 *   k^2 = (s - m^2)(s - (m1 - m2)^2) / (4 s)
 *       = T (T + 2m) (T + 2 m1) (T + 2 m2) / (4 (m + T)^2)    (negative below threshold)
 *   eta = Z1 Z2 alpha / beta = Z1 Z2 alpha (s - m1^2 - m2^2) / (2 k sqrt(s))
 *       = Z1 Z2 alpha (2 m1 m2 + 2 m T + T^2) / (2 k (m + T))
 *
 * (Haeck, Hale, Paris and Thompson, LA-UR-19-31084 (2019), the ENDF KRL=1 kinematics
 * of the LANL code EDA).  Both tend to the non-relativistic forms as T/m -> 0, and the
 * threshold quantities (k^2/T and eta k at T = 0) are exactly the non-relativistic ones.
 * The R-matrix denominators stay E_lambda - E_x, so level energies remain excitation
 * energies; EData::WriteRelativisticLevels writes them in the E(s) convention of EDA.
 *
 * Masses are in MeV/c^2 (PPair masses in u times uconv), energies in MeV, k in fm^-1.
 */
namespace Kinematics {

/// Is relativistic kinematics on for this run?
inline bool Relativistic() { return g_config && g_config->relativistic; }

/// k^2 (fm^-2) at c.m. kinetic energy T (MeV); negative below threshold.
inline double K2(double m1, double m2, double T) {
  const double m = m1 + m2;
  return T * (T + 2.0 * m) * (T + 2.0 * m1) * (T + 2.0 * m2) / (4.0 * (m + T) * (m + T)) / (hbarc * hbarc);
}

/// |k| (fm^-1): the wave number above threshold, the decay constant kappa below.
inline double K(double m1, double m2, double T) { return std::sqrt(std::fabs(K2(m1, m2, T))); }

/// eta at T (above threshold), or eta_b = Z1 Z2 alpha mu_eff / kappa below it (T < 0).
inline double Eta(double z1z2, double m1, double m2, double T) {
  const double m = m1 + m2;
  const double p = K(m1, m2, T) * hbarc;  // MeV
  return z1z2 * fstruc * (2.0 * m1 * m2 + 2.0 * m * T + T * T) / (2.0 * p * (m + T));
}

/// c.m. kinetic energy from the lab kinetic energy of particle 1 on particle 2 at rest.
inline double CMFromLab(double m1, double m2, double eLab) {
  const double m = m1 + m2;
  return 2.0 * m2 * eLab / (std::sqrt(m * m + 2.0 * m2 * eLab) + m);
}

/// Lab kinetic energy of particle 1 (particle 2 at rest) for c.m. kinetic energy T.
inline double LabFromCM(double m1, double m2, double T) { return T * (2.0 * (m1 + m2) + T) / (2.0 * m2); }

/// c.m. kinetic energy from the kinetic energy of particle 1 of a two-body decay at rest.
inline double CMFromDecay(double m1, double m2, double e1) {
  const double p2 = e1 * (e1 + 2.0 * m1);
  return e1 + p2 / (std::sqrt(m2 * m2 + p2) + m2);
}

/*!
 * The kinematic quantities of one particle pair.  The non-relativistic branches are
 * AZURE2's expressions as they were before relativistic kinematics existed, operation
 * for operation, so that a run without "relativistic=1" is bit-identical.
 */
struct Channel {
  int z1 = 0, z2 = 0;
  double redmass = 0.0;  ///< u
  double m1 = 0.0, m2 = 0.0;  ///< MeV/c^2
  bool rel = false;

  Channel() = default;
  explicit Channel(const PPair *pair, bool relativistic = Relativistic())
      : z1(pair->GetZ(1)), z2(pair->GetZ(2)), redmass(pair->GetRedMass()),
        m1(pair->GetM(1) * uconv), m2(pair->GetM(2) * uconv), rel(relativistic) {}
  Channel(int z1_, int z2_, double redmass_, double m1u, double m2u, bool relativistic)
      : z1(z1_), z2(z2_), redmass(redmass_), m1(m1u * uconv), m2(m2u * uconv), rel(relativistic) {}

  /// eta at c.m. kinetic energy E > 0.
  double Eta(double energy) const {
    if (!rel) return sqrt(uconv / 2.) * fstruc * z1 * z2 * sqrt(redmass / energy);
    return Kinematics::Eta((double)z1 * z2, m1, m2, energy);
  }
  /// rho = k a at c.m. kinetic energy E > 0.
  double Rho(double radius, double energy) const {
    if (!rel) return sqrt(2. * uconv) / hbarc * radius * sqrt(redmass * energy);
    return Kinematics::K(m1, m2, energy) * radius;
  }
  /// eta_b at binding energy B > 0 (the channel's T = -B).
  double EtaBound(double binding) const {
    if (!rel) return sqrt(uconv / 2.) * fstruc * z1 * z2 * sqrt(redmass / binding);
    return Kinematics::Eta((double)z1 * z2, m1, m2, -binding);
  }
  /// kappa a at binding energy B > 0.
  double RhoBound(double radius, double binding) const {
    if (!rel) return sqrt(2. * uconv) / hbarc * radius * sqrt(redmass * binding);
    return Kinematics::K(m1, m2, -binding) * radius;
  }
  /// (c / v)^(1/2), v the relative velocity at c.m. kinetic energy E > 0: the unit-flux factor of the
  /// external-capture overlap (non-relativistically (mu c^2 / 2E)^(1/4)).
  double FluxFactor(double energy) const {
    if (!rel) return pow(redmass * uconv / 2. / energy, 0.25);
    const double m = m1 + m2, p = Kinematics::K(m1, m2, energy) * hbarc;
    const double beta = 2.0 * p * (m + energy) / (2.0 * m1 * m2 + 2.0 * m * energy + energy * energy);
    return pow(1.0 / beta, 0.5);
  }
  /// The geometric factor pi / k^2 (fm^2) at c.m. kinetic energy E > 0.
  double PiOverK2(double energy) const {
    if (!rel) return pi * pow(hbarc, 2.) / (2 * redmass * uconv * energy);
    return pi / Kinematics::K2(m1, m2, energy);
  }
};

}  // namespace Kinematics

#endif
