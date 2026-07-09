#ifndef THMFUNC_H
#define THMFUNC_H

///Trojan Horse Method (modified R-matrix) transfer form factor.

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

///Spherical Bessel function j_l(x) (thin wrapper for consistency/testing).
double ThmSphericalBessel(int l, double x);

///Half-off-shell momentum rho = sqrt(2 mu (E + B)) r / (hbar c) (dimensionless).
/*!
 * \param mu     entrance channel reduced mass in MeV/c^2 (PPair::GetRedMass()*uconv)
 * \param E      entrance channel c.m. energy (MeV)
 * \param B      THM binding energy of the transferred particle (MeV)
 * \param radius entrance channel radius (fm)
 */
double ThmRho(double mu, double E, double B, double radius);

///THM transfer form factor M_l(E).
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

#endif
