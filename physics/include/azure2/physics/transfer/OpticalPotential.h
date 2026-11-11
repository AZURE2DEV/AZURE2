#ifndef AZURE2_PHYSICS_TRANSFER_OPTICALPOTENTIAL_H
#define AZURE2_PHYSICS_TRANSFER_OPTICALPOTENTIAL_H

#include <complex>
#include <string>

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The optical potential of a pair of nuclei (docs/upgrade/transfer.md section 2), as a value: explicit numbers, no
 * global parameter set.  At the distance r (fm) between the two centres,
 *
 *   U(r) = V_C(r) - V f(r; R_V, a_V) - i W f(r; R_W, a_W) - i 4 W_D f(r; R_D, a_D) [1 - f(r; R_D, a_D)]
 *          + [V_so + i W_so] (hbar/m_pi c)^2 (1/r) d f(r; R_so, a_so)/dr <2 l.s>,                            (2.1)
 *
 *   f(r; R, a) = 1 / (1 + exp((r - R)/a)),   (hbar/m_pi c)^2 = 2.000 fm^2,
 *   <2 l.s> = j(j+1) - l(l+1) - s(s+1),
 *
 * all in MeV; V, W, W_D > 0 attract and absorb.  The surface term is the derivative form, 4 f (1 - f) = -4 a df/dr
 * (= 1 at r = R_D), FRESCO's TYPE=2 SHAPE=0.  The spin-orbit term is FRESCO's TYPE=3 SHAPE=0 (Thompson, Comput.
 * Phys. Rep. 7 (1988) 167; the FRESCO namelist manual, "&POT"): the form factor times 2 l.s, s the spin of the
 * projectile; for s = 1/2, 2 l.s = l.sigma, so a strength quoted for l.s without the factor 2 is twice V_so here.
 * Its radial factor is physics/HybridPotential.h's SpinOrbitRadialTerm, the same convention.  V_C is the potential
 * of a uniformly charged sphere of radius R_C, Z1 Z2 e^2 (3 - r^2/R_C^2)/(2 R_C) inside and Z1 Z2 e^2 / r outside
 * (R_C = 0: a point charge), with e^2 = alpha hbar c of physics/Constants.h.
 *
 * Radii.  Every radius parameter r0 (fm) of one potential becomes a radius R (fm) by the potential's one radius
 * convention: R = r0 (Absolute), R = r0 A_t^{1/3} (TargetMass) or R = r0 (A_p^{1/3} + A_t^{1/3}) (BothMasses), with
 * the mass numbers A_p, A_t the potential states (FRESCO's TYPE=0 p1 = A_t, p2 = A_p; A_p = 0 is TargetMass).
 *
 * A term with depth 0 contributes nothing whatever its geometry; a term with a non-zero depth needs R > 0, a > 0.
 */

/// How a radius parameter r0 (fm) becomes a radius R (fm).
enum class RadiusConvention {
  Absolute,    ///< R = r0
  TargetMass,  ///< R = r0 A_t^{1/3}
  BothMasses   ///< R = r0 (A_p^{1/3} + A_t^{1/3})
};

/// One Woods-Saxon term: its depth (MeV) and the geometry of its form factor (r0 converted by the convention).
struct WoodsSaxonTerm {
  double depth_MeV = 0.0;       ///< V, W, W_D, V_so or W_so (MeV), sign as in Eq. (2.1)
  double r0_fm = 0.0;           ///< radius parameter (fm)
  double diffuseness_fm = 0.0;  ///< a (fm)
};
bool operator==(const WoodsSaxonTerm &a, const WoodsSaxonTerm &b);

/// The potential of Eq. (2.1).
struct OpticalPotential {
  RadiusConvention convention = RadiusConvention::Absolute;
  double target_mass_number = 0.0;      ///< A_t of the radius convention
  double projectile_mass_number = 0.0;  ///< A_p of the radius convention (BothMasses only)
  double coulomb_r0_fm = 0.0;           ///< R_C parameter (fm); 0: point charge
  WoodsSaxonTerm real_volume;           ///< V
  WoodsSaxonTerm imaginary_volume;      ///< W
  WoodsSaxonTerm imaginary_surface;     ///< W_D
  WoodsSaxonTerm real_spin_orbit;       ///< V_so
  WoodsSaxonTerm imaginary_spin_orbit;  ///< W_so
};
bool operator==(const OpticalPotential &a, const OpticalPotential &b);

/// R (fm) of a radius parameter r0 (fm) under the potential's convention.  Throws std::invalid_argument for a
/// mass number the convention needs that is not positive.
double RadiusOf(const OpticalPotential &p, double r0_fm);

/// The values of the terms of Eq. (2.1) at one radius, each in MeV.
struct PotentialTerms {
  double coulomb_MeV = 0.0;                    ///< V_C (for the charge product the caller gave)
  std::complex<double> volume_MeV;             ///< -V f - i W f
  std::complex<double> surface_MeV;            ///< -i 4 W_D f (1 - f)
  std::complex<double> spin_orbit_form_MeV;    ///< (V_so + i W_so) 2 fm^2 (1/r) df/dr, before <2 l.s>
};

/// The terms at r > 0 (fm) for the charge product z1 z2.  Throws std::invalid_argument for a term with a non-zero
/// depth and a non-positive radius or diffuseness.
PotentialTerms TermsAt(const OpticalPotential &p, int z1z2, double r_fm);

/// U(r) of Eq. (2.1) (MeV) at r > 0 (fm) for the charge product z1 z2 and <2 l.s> = two_l_dot_s.
std::complex<double> PotentialAt(const OpticalPotential &p, int z1z2, double two_l_dot_s, double r_fm);

/// The potential in words and the radii it uses (one line per term), for a listing next to FRESCO's.
std::string Describe(const OpticalPotential &p);

/// A potential with no term (U = V_C = 0): the plane-wave limit.
OpticalPotential NoPotential();

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
