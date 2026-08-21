#ifndef CAPTUREAYTERM_H
#define CAPTUREAYTERM_H

///One \f$(t,t')\f$ pathway pair of a capture angular distribution.

/*!
 * The Legendre coefficients of a particle-capture-\f$\gamma\f$ angular
 * distribution, R. G. Seyler and H. R. Weller, Phys. Rev. C \b 20 (1979) 453,
 * Eqs. (20) and (21) in the channel-spin representation:
 *
 * \f[ \sigma(\theta,\phi) = N \sum_k \left[ a_k P_k(\cos\theta)
 *                                 + b_k P_k^1(\cos\theta)\, p_y \right] \f]
 * so that
 * \f[ A_y(\theta) = \frac{\sum_k b_k P_k^1(\cos\theta)}
 *                        {\sum_k a_k P_k(\cos\theta)} . \f]
 *
 * Both \f$a_k\f$ and \f$b_k\f$ are sums over pairs of reaction pathways
 * \f$t = \{p L b l s\}\f$, weighted by \f$\mathrm{Re}\,R R'^*\f$ and
 * \f$\mathrm{Re}\,(i R R'^*)\f$ respectively, with \f$R\f$ the collision-matrix
 * element AZURE2 already forms as its T-matrix element. One CaptureAyTerm holds
 * the geometry of one such pair: which two pathways, at which Legendre order,
 * and the two coefficients.
 *
 * Unlike the Interference objects that carry the unpolarized angular
 * distribution, a term may pair pathways from *different* KGroups: \f$a_k\f$
 * requires \f$s = s'\f$, but \f$b_k\f$ does not, and that channel-spin mixing
 * is exactly the information a polarization measurement adds. Hence each side
 * carries its own KGroup index.
 */

struct CaptureAyTerm {
  int kOrder;      //!< Legendre order \f$k\f$
  int kGroup1;     //!< KGroup of the first pathway (1-based)
  int path1;       //!< MGroup (isEC1 false) or ECMGroup (true) index, 1-based
  bool isEC1;      //!< first pathway is external capture
  int kGroup2;     //!< KGroup of the second pathway (1-based)
  int path2;       //!< MGroup or ECMGroup index, 1-based
  bool isEC2;      //!< second pathway is external capture
  double ak;       //!< Eq. (20) coefficient, multiplying \f$\mathrm{Re}\,T T'^*\f$
  double bk;       //!< Eq. (21) coefficient, multiplying \f$\mathrm{Re}\,(i T T'^*)\f$
};

#endif
