#ifndef AZURE2_PHYSICS_OBSERVABLECOUPLING_H
#define AZURE2_PHYSICS_OBSERVABLECOUPLING_H

namespace azure2 {
namespace physics {

/*!
 * Coupling coefficients of the observables: the Legendre coefficients of angular distributions and of capture
 * analyzing powers, per pair of pathways (docs/upgrade/cnuc_split.md section 9).  All angular momenta are passed
 * doubled (two_j = 2j) except orbital angular momenta and multipolarities (l, L) and the Legendre order (k), which
 * are integers and passed as they are.  The coefficients use AngularCoupling's 3j, 6j and 9j; x^ = (2x + 1)^(1/2).
 * The channel spin couples the light particle first, (I_1 I_2 | s), as everywhere in AZURE2.
 */

/// Blatt-Biedenharn Zbar(l1 j1 l2 j2; s k) = l1^ l2^ j1^ j2^ (l1 0 l2 0 | k 0) W(l1 j1 l2 j2; s k)   (9.1).
double z_bar(int l1, int two_j1, int l2, int two_j2, int two_s, int k);

/// The photon form, Zbar(L1 j1 L2 j2; jf k) = L1^ L2^ j1^ j2^ (L1 1 L2 -1 | k 0) W(L1 j1 L2 j2; jf k)   (9.2).
double z_bar_photon(int L1, int two_j1, int L2, int two_j2, int two_jf, int k);

/*!
 * Interference coefficient of two pathways (l1, j1 -> l1') and (l2, j2 -> l2') of a particle exit, entrance channel
 * spin s and exit s': (-1)^(s' - s) / 4 Zbar(l1 j1 l2 j2; s k) Zbar(l1' j1 l2' j2; s' k).
 */
double particle_interference(int l1, int l1p, int two_j1, int l2, int l2p, int two_j2, int two_s, int two_sp, int k);

/*!
 * Interference coefficient of two pathways (l1, j1 -> photon L1) and (l2, j2 -> L2) to a final state of spin jf,
 * entrance channel spin s: (-1)^(1 + s - jf) / 4 Zbar(l1 j1 l2 j2; s k) Zbar_photon(L1 j1 L2 j2; jf k).
 */
double photon_interference(int l1, int L1, int two_j1, int l2, int L2, int two_j2, int two_s, int two_jf, int k);

/*!
 * Secondary-particle coefficient of an unobserved primary (Brune and deBoer, PRC 102, 024628 (2020), Eq. 22) for
 * the exit channel spins s' (first amplitude) and s2' (second), in the light-first channel-spin order, which
 * multiplies each amplitude by (-1)^(I_b + I_B - s') and the product by (-1)^(s' - s2'); I_b, I_B are the spins
 * two_jb (light) and two_jB (heavy) of the exit pair.  l1' is the exit orbital angular momentum of the first pathway.
 */
double upos_interference(int l1, int l1p, int two_j1, int l2, int two_j2, int two_s, int two_sp, int two_sp2,
                         int two_jb, int two_jB, int k);

/// One capture pathway t = {p L b l s} of Seyler and Weller: photon mode p (1 electric, 0 magnetic), multipolarity
/// L, compound spin b, entrance orbital angular momentum l, channel spin s.
struct CapturePathway {
  int l;
  int L;
  int two_b;
  int two_s;
  int p;
};

/// The pair (a_k, b_k) of one pair of capture pathways at Legendre order k.
struct CaptureAyCoefficients {
  double a = 0.0;
  double b = 0.0;
};

/*!
 * a_k and b_k of R. G. Seyler and H. R. Weller, Phys. Rev. C 20 (1979) 453, Eqs. (20) and (21), for pathways t1, t2,
 * projectile spin x, target spin a, residual spin c; b_k carries the (-1)^(s1 - s2) of the light-first channel spin.
 * a_k needs s1 = s2, b_k does not; b_0 = 0.  Both are 0 when the parity bracket of Eq. (15) vanishes, when any of
 * (l1 0 l2 0 | k 0), (L1 1 L2 -1 | k 0), W(L1 b1 L2 b2; c k) or (b_k only) W(x s1 x s2; a 1) is below 1e-12 in
 * magnitude, or (b_k only) when the 9j is below 1e-14.  The 9j is the bare symbol, which Seyler and Weller call
 * X(l s b; l' s' b'; k 1 k), not the [prod (2j+1)]^(1/2)-normalized Fano X, as their worked example on p. 458 shows.
 */
CaptureAyCoefficients capture_ay_coefficients(const CapturePathway &t1, const CapturePathway &t2, int k, int two_x,
                                              int two_a, int two_c);

}  // namespace physics
}  // namespace azure2

#endif
