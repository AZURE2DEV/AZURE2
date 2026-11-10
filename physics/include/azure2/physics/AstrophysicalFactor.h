#ifndef AZURE2_PHYSICS_ASTROPHYSICALFACTOR_H
#define AZURE2_PHYSICS_ASTROPHYSICALFACTOR_H

#include "azure2/physics/ChannelFunctions.h"

namespace azure2 {
namespace physics {

/*!
 * The astrophysical S factor of a charged entrance pair (docs/upgrade/data_kinds.md section 4.1; LEXFOR pp. A.6-A.7):
 *
 *   S(E) = E sigma(E) exp(2 pi eta(E)),   eta = Z1 Z2 alpha sqrt(mu c^2 / (2 E)),
 *
 * E the c.m. kinetic energy (MeV), sigma in b, S in MeV b; mu the reduced mass of the atomic masses of the pair (u, as
 * the engine's <levels>), alpha and the atomic mass unit from Constants.h.  eta is the non-relativistic Sommerfeld
 * parameter of ChannelFunctions.h whatever kin.relativistic says: the S factor is a convention for the data, which
 * must not change with the kinematics of a run.  sfactor_conversion is the factor E exp(2 pi eta) (MeV); it is
 * defined for E > 0 and Z1 Z2 > 0 (NaN otherwise).  Mass numbers in place of the masses change 2 pi eta and so S by
 * exp(2 pi Delta eta) - 1: +2.9 % for 3He + alpha at 93 keV, +6.7 % for p + 12C at 100 keV.
 */
double sfactor_conversion(const ChannelKinematics &kin, double E_MeV);

/// S (MeV b) from sigma (b) at E (MeV, c.m.).
inline double sfactor_from_cross_section(const ChannelKinematics &kin, double E_MeV, double sigma_b) {
  return sigma_b * sfactor_conversion(kin, E_MeV);
}
/// sigma (b) from S (MeV b) at E (MeV, c.m.).
inline double cross_section_from_sfactor(const ChannelKinematics &kin, double E_MeV, double S_MeVb) {
  return S_MeVb / sfactor_conversion(kin, E_MeV);
}

}  // namespace physics
}  // namespace azure2

#endif
