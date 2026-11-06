#ifndef THM_LINESHAPE_H
#define THM_LINESHAPE_H

#include "Constants.h"
#include "ThmDistortion.h"

#include <atomic>
#include <functional>
#include <memory>
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
 * carry no width).  Under Park (--use-park, which implies Brune) the
 * amplitudes are the observed ones, gamma_Park^2 = gamma_Brune^2 / (1 + sum),
 * and Gamma = sum_c 2 gamma_c^2 P_c: the same width.  Memoized per thread on
 * the level, its parameters and the parametrization.
 */
double ThmLevelWidth(CNuc *compound, JGroup *jgroup, ALevel *level, const Config &configure);

struct ThmExperiment;

/*!
 * The spectator-momentum window of a THM experiment with the plane-wave
 * vertex (<thm> experiment[...] ps=..., with or without spectatorAngles=;
 * docs/source/theory/thm_implementation.rst, "Spectator-momentum window" and
 * "Experimental acceptance").
 *
 * The off-shell x-A momentum of the entrance vertex depends on the
 * spectator momentum p_s = p_sx (the s-x relative momentum in a, = the
 * spectator momentum in the rest frame of a):
 *
 *   p_xA^2 / 2 mu_xA = E + B_xs + p_s^2 / 2 mu_sx
 *
 * (Mukhamedzhanov et al., PRC 96 (2017) 024623, eq. 31; Typel & Baur, Ann.
 * Phys. 305 (2003) 228, eq. 11).  At fixed E the three-body phase space is
 * d Omega_sF (|k_sF| is fixed by energy conservation), and the TH double
 * differential cross section per d Omega_sF dE is |phi_a(p_s)|^2 times the
 * HOES cross section at p_xA(p_s) (2017 eq. 34); the direction of k_sF fixes
 * p_s = |k_sF - beta k_aA|.  The data are the yield of an energy bin over the
 * Monte Carlo integral of KF |phi|^2 over the same events (KF times the
 * Jacobian of the measured variables is the phase-space density, a function
 * of E alone), so the HOES "datum" at E is the mean of the HOES cross
 * section over the accepted directions with the event weight
 *
 *   A(theta) |phi(q)|^2 d cos(theta_cm),   d cos(theta_cm) = q dq / (beta k_sF k_aA),
 *
 * A the acceptance (spectatorAngles=; 1 without), q restricted to the ps
 * window [pmin, pmax] and to what the kinematics reach at E,
 * |k_sF - beta k_aA| <= q <= k_sF + beta k_aA.  Not |phi|^2 p^2 dp: that is
 * the measure of events integrated over E as well (d^3 p_s = d^3 k_sF).
 * The nodes are Gauss-Legendre in cos(theta_cm) on every accepted interval
 * (ThmDistortion::AngleNodes, the same directions R(E) and the DW vertex
 * average over), so they depend on E: every point and folding sub-point
 * takes its own (EPoint::ThmPsTable), node k adding T_k = q_k^2 / 2 mu_sx to
 * E + B in the vertex (and in the Coulomb term C_l), as spectatorEnergy does.
 */
struct ThmSpectatorWindow {
  std::string experiment;
  std::string description;  ///< for the output: distribution, cut, acceptance, nodes
  double muSx = 0.0;        ///< reduced mass of s + x (MeV)
  double pMin = 0.0, pMax = 0.0;  ///< |p_s| cut (MeV/c)
  /// |phi(p)|^2 of the spectator momentum distribution (p in MeV/c, any scale).
  std::function<double(double)> phi2;
  /// Kinematics and accepted directions (ThmDistortion::Setup without waves).
  std::shared_ptr<const ThmDistortion> acc;
  /// Energies of the data points (sorted): a folding sub-point beyond the
  /// reach of the window takes the nodes of the nearest of them.
  std::vector<double> dataE;
  struct Node {
    double theta = 0.0;   ///< c.m. angle of the spectator to the beam (deg)
    double p = 0.0;       ///< |p_s| = q (MeV/c)
    double weight = 0.0;  ///< normalized (sum 1 over the nodes of one energy)
    double es = 0.0;      ///< T = p^2 / 2 mu_sx (MeV)
  };
  /// The accepted nodes at E (weight > 0); false if no direction is accepted.
  bool NodesAt(double energy, std::vector<Node> &out) const;
  /// NodesAt, or (no direction accepted at E) the nodes of the nearest data
  /// energy; *moved is set then.
  void Nodes(double energy, std::vector<Node> &out, bool *moved = nullptr) const;
  /// <T_s> = sum_k w_k T_k at E (MeV), 0 if nothing is accepted there.
  double MeanEs(double energy) const;
  /// The |p_s| the kinematics reach at E (MeV/c): |k_sF - beta k_aA| to
  /// k_sF + beta k_aA; false if E_sF <= 0.
  bool Reach(double energy, double &qLo, double &qHi) const;
  /// The event weight per unit |p_s| at E, unnormalized: A(theta(q)) |phi(q)|^2 q
  /// (d cos theta_cm = q dq / beta k_sF k_aA), 0 outside the acceptance or the
  /// reach (q in MeV/c).  What the nodes integrate; for plots.
  double Density(double energy, double q) const;
};

/*!
 * The window of experiment x (not PS_DELTA) for the s + x reduced mass muSx
 * (MeV), the kinematics k (ThmDistortion::Kinematics, as BuildThmGroups fills
 * them) and the data energies (c.m. of x + A, MeV), every one of which must
 * have an accepted direction.  |phi|^2: Hulthen, (1/(a^2 + q^2) -
 * 1/(b^2 + q^2))^2 with q = p/hbar c (Tribble et al., RPP 77 (2014) 106901,
 * eq. 4.4); gauss, exp(-4 ln 2 p^2/FWHM^2); a table gives |phi(p)|^2 itself
 * (linear between rows; its range is the cut).  "" or what is wrong.
 */
std::string BuildThmSpectatorWindow(const ThmExperiment &x, double muSx, const ThmDistortion::Kinematics &k,
                                    const std::vector<double> &dataE, ThmSpectatorWindow &out);

/// |phi(p)|^2 (p in MeV/c) of experiment x's ps distribution, as
/// BuildThmSpectatorWindow uses it; null for ps=delta.
std::function<double(double)> ThmPsDistribution(const ThmExperiment &x);

#endif
