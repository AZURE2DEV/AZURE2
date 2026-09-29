#ifndef THM_LINESHAPE_H
#define THM_LINESHAPE_H

#include "Constants.h"

#include <atomic>
#include <string>
#include <vector>

class CNuc;
class JGroup;
class ALevel;
class Config;

/*!
 * The Coulomb line-shape factor N_C of a THM experiment (<thm> experiment[...]
 * lineshape=on; docs/source/theory/thm_implementation.rst, "Coulomb line
 * shape").
 *
 * In A + a(x+s) -> s + F*(x+A) -> s + b + B the spectator s leaves in the
 * Coulomb field of the resonance F*, and after F* decays in that of b and B.
 * Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020) 233 (arXiv:2007.13331), eqs.
 * (55)-(57), and Mukhamedzhanov, EPJA 58 (2022) 71 (arXiv:2203.04094), eqs.
 * (32)-(34): the resonance factor 1/(E_0 - E - i Gamma/2) of the THM amplitude
 * becomes (E_0 - E - i Gamma/2)^(-1 - i zeta) times level-independent Coulomb
 * factors, with
 *
 *   zeta = eta_sb + eta_sB - eta_0,   eta_ij = Z_i Z_j alpha mu_ij / k_ij,
 *
 * eta_0 = Z_s Z_F alpha mu_sF / k_sF the s-F* Sommerfeld parameter of the
 * intermediate state.  AZURE2 takes the paper's own tractable limit (2020 p. 15
 * case 2, 2022 eq. 34): B much heavier than s and b, so k_sB ~ k_sF, and
 * |eta_sb| << 1, dropped:
 *
 *   zeta(E) = Z_s alpha (Z_B mu_sB - Z_F mu_sF) / k_sF(E)      (< 0 when Z_B < Z_F)
 *
 * with E_sF = E_aA - B_xs - E (energy conservation in the three-body final
 * state; E_aA the beam-target c.m. energy, B_xs the binding of x in a, E the
 * x + A c.m. energy from the x + A threshold).  eta_0 is taken at the point's
 * k_sF(E) rather than at the pole momentum k_0: the two agree at E = E_lambda
 * (narrow resonance), and this keeps zeta common to all levels at one energy,
 * so the arbitrary energy unit of the complex power is a common phase.
 *
 * The factor applied to level lambda (pole E_lambda - i Gamma_lambda/2) is
 *
 *   N_C = exp(pi zeta/2) (E_lambda - E - i Gamma_lambda/2)^(-i zeta)      (MeV)
 *   |N_C|^2 = exp[2 zeta arctan(2 (E_lambda - E)/Gamma_lambda)]    (2020 eq. 62)
 *
 * normalized to 1 at the pole; the level-independent factors of eq. (56) and
 * the sinh/hypergeometric prefactor of eq. (62) are smooth and left to the
 * arbitrary THM normalization.
 */
struct ThmLineshape {
  std::string experiment;
  std::string spectator;  ///< name, for the output
  int Zs = 0, ZF = 0;
  double ms = 0.0, mF = 0.0;  ///< u (nuclear masses; F = x + A)
  double eAA = 0.0;           ///< beam-target c.m. energy (MeV)
  double bind = 0.0;          ///< B_xs (MeV)
  /// Exit pairs of the experiment's segments: b the lighter nucleus, B the
  /// heavier; E_bB = E + q (q = entrance minus exit threshold, MeV).
  struct Exit {
    int pairKey = 0;
    int Zb = 0, ZB = 0;
    double mb = 0.0, mB = 0.0;  ///< u
    double q = 0.0;
  };
  std::vector<Exit> exits;
  /// Set once an energy with E_sF <= 0 has been reported.
  mutable std::atomic<bool> warned{false};

  /// E_sF(E) = E_aA - B_xs - E (MeV).
  double EsF(double energy) const { return eAA - bind - energy; }
  /// Reduced mass of s + F (MeV).
  double MuSF() const { return ms * mF / (ms + mF) * uconv; }
  /// k_sF (fm^-1) at E_sF(E), E_sF clamped at 1 keV.
  double KsF(double energy) const;
  /// eta_0 = Z_s Z_F alpha mu_sF / k_sF.
  double Eta0(double energy) const;
  /// eta_sB = Z_s Z_B alpha mu_sB / k_sF (paper eq. 34: k_sB ~ k_sF).
  double EtaSB(double energy, int ZB, double mB) const;
  /// zeta = eta_sB - eta_0 of the exit pair whose heavy nucleus is (Z_B, m_B).
  double Zeta(double energy, int ZB, double mB) const { return EtaSB(energy, ZB, mB) - Eta0(energy); }
  /*!
   * Size of the neglected eta_sb (validity of the case-2 limit): Z_s Z_b alpha
   * c <1/|v_s - v_b|> in the F* rest frame, the average over an isotropic
   * direction of b, <1/|v_s - v_b|> = 1/max(v_s, v_b).  eBB is the b + B
   * relative energy (MeV); masses in u.
   */
  double EtaSbEstimate(double energy, double eBB, int Zb, double mb, double mB) const;
};

/// N_C for zeta, x = E_lambda - E and Gamma (MeV).  Gamma <= 0: the limit
/// Gamma -> 0+ (a step exp(+-pi zeta)); the phase is -zeta ln|x - i Gamma/2|.
complex ThmLineshapeFactor(double zeta, double x, double Gamma);
/// |N_C|^2 = exp[2 zeta arctan(2 x / Gamma)], the same conventions.
double ThmLineshapeFactorSq(double zeta, double x, double Gamma);

/*!
 * Pole of a level as the line-shape factor uses it: the level energy
 * (compound excitation, MeV) and total width (MeV) of the fit parameters.
 * Under Brune (required by lineshape=on) GetFitE is the observed energy and
 * Gamma = sum_c 2 gamma_c^2 P_c(E_lambda) / (1 + sum_c gamma_c^2 dS_c/dE)
 * over the open particle channels (the normalization over all particle
 * channels), plus the radiative widths, as CNuc::TransformOut writes them to
 * parameters.out (internal part; a ground-state moment and closed channels
 * carry no width).  Memoized per thread on the level and its parameters.
 */
double ThmLevelWidth(CNuc *compound, JGroup *jgroup, ALevel *level, const Config &configure);

/// The line shape of one experiment on an energy grid (EData::ThmLineshapeTable).
struct ThmLineshapeReport {
  std::string experiment, spectator;
  int Zs = 0, ZF = 0;
  double eAA = 0.0, bind = 0.0;
  std::vector<double> energy;  ///< E, c.m. of x + A (MeV)
  std::vector<double> esf;     ///< E_sF(E) (MeV)
  std::vector<double> eta0;    ///< eta_0(E)
  struct Level {
    int jgroup = 0, level = 0;  ///< 1-based, as parameters.out
    double J = 0.0;
    int pi = 1;
    double energy = 0.0;  ///< E_lambda, c.m. of x + A (MeV)
    double width = 0.0;   ///< Gamma_lambda (MeV)
    std::vector<double> nc2;  ///< |N_C|^2 on the grid
  };
  struct Exit {
    int pairKey = 0;
    int Zb = 0, ZB = 0;          ///< light (b) and heavy (B) nucleus of the pair
    double mb = 0.0, mB = 0.0;   ///< u
    std::vector<double> zeta, etaSb;
    std::vector<Level> levels;
  };
  std::vector<Exit> exits;
};

#endif
