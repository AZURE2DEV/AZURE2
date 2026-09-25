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
 * decaying Whittaker function W_{-eta, l+1/2}(2 kappa r) in its place.
 * 2 eta k = 2 Z1 Z2 alpha mu c^2/(hbar c) (fm^-1) does not depend on E.
 * Dimensionless; exactly 0 for a neutral pair.
 *
 * Sign and normalisation.  The term is ADDED to the surface part
 * (B - 1) j_l(rho) - rho j_l'(rho) (EPoint::GetThmFormFactor), with a plus sign
 * and no further factor.  All three references agree: Tribble 2.79 and
 * Mukhamedzhanov 27 write j_l [(B - 1) - D] + 2 Z1 Z2 e^2 mu \int O/O(a) j_l
 * (hbar = c = 1, so 2 Z1 Z2 e^2 mu = 2 eta k), and Typel & Baur A.4, rewritten
 * with z_l(x) = x j_l(x), B = a O'/O and D = a j_l'/j_l, reads
 *   J_l = k p O_l(ka)/(k^2 - p^2) [ (B - 1) j_l - rho j_l' + C_l ],
 * the same bracket (their eta > 0 for repulsion, eq. A.1).
 *
 * E = 0 exactly returns the common limit of both sides, O(r)/O(a) ->
 * sqrt(r/a) K_{2l+1}(sqrt(8 eta k r)) / K_{2l+1}(sqrt(8 eta k a)) (real).  For
 * E <= 0 the result is real.  O_l is the point-Coulomb wave: a hybrid nuclear
 * potential, if enabled for the pair, is not applied here.
 *
 * Numerics (details in src/ThmFunc.cpp): Gauss-Legendre panels on the real
 * axis; for E > 0 the conditionally convergent tail beyond the asymptotic
 * region is integrated exactly along complex rays (Cauchy/Jordan), and O_l on
 * [a, R] is obtained by integrating the Coulomb equation inwards from its
 * asymptotic expansion at R (deep below the barrier, where G_l(ka) may
 * overflow, G_l is integrated inwards from inside the barrier instead); for
 * E < 0 WhitFunc is used.  No library Coulomb function is called, so useGSL
 * has no effect (GSL's G_l/F_l are off by up to ~2% for l >= 2, eta ~ 2,
 * rho ~ 1-4, which would otherwise leak into C_l).  Agrees with an independent
 * mpmath evaluation to ~1e-13 relative (tests/reference/thm_coulomb_term_test);
 * ~0.05-1 ms per call, and repeats of the same (pair, l, E, p) -- one per J
 * group and channel -- come from a small per-thread memo.
 *
 * \param pair   entrance particle pair (charges, reduced mass, channel radius a)
 * \param l      orbital angular momentum
 * \param E      entrance c.m. energy (MeV), may be negative
 * \param p      half-off-shell momentum (fm^-1)
 * \param useGSL kept for the interface; has no effect (see above)
 */
complex ThmCoulombTerm(PPair *pair, int l, double E, double p, bool useGSL);

#endif
