#ifndef THMFUNC_H
#define THMFUNC_H

#include "Constants.h"

class PPair;

/// Trojan Horse Method (modified R-matrix) transfer form factor.

/*!
 * In the modified R-matrix formalism (Mukhamedzhanov et al.), the entrance
 * penetrability of conventional R-matrix is replaced at the channel surface by
 * the half-off-energy-shell (HOES) transfer form factor
 *
 *   M_l(E) = (b - 1) j_l(rho) - rho dj_l/drho ,   rho = p r / (hbar c)
 *
 * with the half-off-shell momentum p = sqrt(2 mu (E + B)), where B is the
 * binding energy of the transferred particle in the Trojan-Horse nucleus, mu
 * and r the entrance channel reduced mass and radius, b the entrance channel
 * boundary constant, and j_l the spherical Bessel function.
 *
 * Only entrance-pair channels use M_l; every other (exit) channel keeps the
 * ordinary R-matrix sqrt(penetrability). This mirrors mrmpy's
 * ``MRMModel._form_factor`` (see docs/THM_IMPLEMENTATION.md).
 */

/// Spherical Bessel function j_l(x) (thin wrapper for consistency/testing).
double ThmSphericalBessel(int l, double x);

/// Half-off-shell momentum rho = sqrt(2 mu (E + B)) r / (hbar c) (dimensionless).
/*!
 * \param mu     entrance channel reduced mass in MeV/c^2 (PPair::GetRedMass()*uconv)
 * \param E      entrance channel c.m. energy (MeV)
 * \param B      THM binding energy of the transferred particle (MeV)
 * \param radius entrance channel radius (fm)
 */
double ThmRho(double mu, double E, double B, double radius);

/// Boundary-independent pieces of M_l: jl = j_l(rho) and rhoDjl = rho dj_l/drho.
/*!
 * These depend only on the point energy, binding energy, reduced mass and
 * channel radius, so they can be cached per point.  The boundary-dependent
 * assembly M_l = (b - 1) jl - rhoDjl happens at the entrance vertex, where b
 * is the per-level shift function S_c(E_lambda) under the Brune formalism
 * (mrmpy vertex_boundary="per_level") or the fixed channel boundary constant
 * otherwise.  Parameters as in ThmRho, plus the orbital angular momentum l.
 */
void ThmBesselParts(int l, double mu, double E, double B, double radius,
                    double &jl, double &rhoDjl);

/// THM transfer form factor M_l(E).
/*!
 * \param l      orbital angular momentum of the entrance channel
 * \param b      entrance channel boundary constant
 * \param mu     entrance channel reduced mass in MeV/c^2
 * \param E      entrance channel c.m. energy (MeV)
 * \param B      THM binding energy of the transferred particle (MeV)
 * \param radius entrance channel radius (fm)
 *
 * The derivative dj_l/drho is taken by a forward finite difference with step
 * 1e-6, matching mrmpy for bit-level cross-validation.
 */
double ThmFormFactor(int l, double b, double mu, double E, double B,
                     double radius);

/// External Coulomb term of the THM entrance vertex.
/*!
 *   C_l(E) = 2 eta k \int_a^inf dr O_l(k r) / O_l(k a) j_l(p r)
 *
 * (Tribble et al. 2014 eq. 2.79, last line; Mukhamedzhanov et al. 2017 eq. 27;
 * Typel & Baur 2003 eq. A.4), with O_l = G_l + i F_l the outgoing Coulomb wave
 * of the entrance pair at the on-shell energy E > 0 and, for E < 0, the
 * decaying Whittaker function in its place.  2 eta k = 2 Z1 Z2 alpha mu c^2/(hbar c)
 * does not depend on E.  Dimensionless; 0 for a neutral pair.
 *
 * \param pair   entrance particle pair (charges, reduced mass, channel radius a)
 * \param l      orbital angular momentum
 * \param E      entrance c.m. energy (MeV), may be negative
 * \param p      half-off-shell momentum (fm^-1)
 * \param useGSL Coulomb functions from GSL instead of the built-in routine
 */
complex ThmCoulombTerm(PPair *pair, int l, double E, double p, bool useGSL);

#endif
