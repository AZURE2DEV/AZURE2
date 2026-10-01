#ifndef THM_DW_VERTEX_H
#define THM_DW_VERTEX_H

#include "Constants.h"
#include "ThmDistortion.h"

#include <atomic>
#include <functional>
#include <string>
#include <vector>

struct ThmExperiment;

/*!
 * The distorted-wave entrance vertex of a THM experiment (`vertexModel=dw` on
 * its experiment line; docs/source/theory/thm_implementation.rst,
 * "Distorted-wave entrance vertex").
 *
 * Prior-form DWBA of a + A -> s + F*, the x-A interior by the R matrix
 * (Mukhamedzhanov, PRC 84 (2011) 044616; Mukhamedzhanov, Kadyrov & Pang, EPJA
 * 56 (2020) 233, eqs. 28-32): the vertex of an entrance channel (l) is the
 * surface term at r_xA = a of the source
 *
 *   S(r) = Int d^3u phi(u) chi^(-)*_{k_sF}(alpha r + u) chi^(+)_{k_aA}(r + beta u),
 *
 * r = r_xA, u = r_sx, alpha = m_A/m_F, beta = m_s/m_a:
 *
 *   V_lm(B) = (B - 1) S_lm(a) - a S_lm'(a),   S_lm = Int dOmega Y*_lm S,
 *
 * normalized by 4 pi phi~(q) (the plane-wave source is phi~(q) e^{i p.r}, so
 * without distortion V_lm/4pi phi~ = i^l Y*_lm(p^) M_l(p), the plane-wave
 * vertex).  The angle-integrated observable is incoherent in l and m, so
 * |M_l|^2 becomes (4 pi/(2l+1)) sum_m |V_lm/4 pi phi~|^2 = c^+ G c with
 * c = (B - 1, -1) and the 2x2 Gram matrix G of (s_m, d_m) = (S_lm(a),
 * a S_lm'(a))/4 pi phi~.  G = L^+ L splits the vertex into two incoherent
 * components M^(k)(B) = a_k (B - 1) - d_k (Factor).
 *
 * The surface term only: the external prior term (the three-body remnant
 * outside a, whose plane-wave limit is C_l of coulombIntegral=1) is not
 * computed, so vertexModel=dw is the counterpart of coulombIntegral=0.
 *
 * Numerics: reduced amplitudes h^l_{Ls La}(a) (independent of the directions
 * of k_aA and k_sF) by a 2D quadrature in (u, cos theta_u), the waves of
 * ThmDistortion (Numerov + COUL; plane waves as Riccati-Bessel functions) read
 * by six-point Lagrange interpolation; G at every node on a 20 keV grid,
 * interpolated by cubic Lagrange.
 */
class ThmDwVertex {
 public:
  std::string experiment;
  std::string description;  ///< settings, for the output
  ThmDistortion dist;       ///< settings, kinematics, channels and bound state (ThmDistortion::Setup)
  double radius = 0.0;      ///< channel radius a of the entrance pair (fm)
  std::vector<int> lvals;   ///< entrance orbital momenta, ascending
  double alpha = 0.0, beta = 0.0;
  // Spectator-momentum window (ps=): the nodes are put on the part of
  // [pMin, pMax] (MeV/c) that the kinematics reach at each energy.
  bool window = false;
  double pMin = 0.0, pMax = 0.0;
  int psNodes = 1;
  std::function<double(double)> psWeight;
  /// Spectator-direction window (spectatorAngles=, ThmDistortion::AngleNodes):
  /// the nodes are the accepted directions at each energy, weighted by
  /// d cos(theta_cm) x acceptance x |phi~(q)|^2; a ps window only cuts q.
  bool angles = false;
  // Quadrature and partial waves (for the output).
  int laMax = 0, lsMax = 0, uNodes = 0, cNodes = 0;
  double uMax = 0.0;
  // The grid.
  double gridLo = 0.0, gridStep = 0.02;
  int nE = 0;
  int nNodes = 1;  ///< nodes per energy (1 without a window)
  /// Per grid energy e and node k: weight w[e nNodes + k] (sum 1 over k), q in
  /// MeV/c, and per l index li the Gram entries g11, g22, Re g12, Im g12 at
  /// G[((e nNodes + k) nl + li) 4 + c].
  std::vector<double> w, q, G;
  /// spectatorAngles: nAng directions per grid energy e (nNodes is then 1,
  /// with the weight-averaged G, since the model is linear in G): normalized
  /// weight aw, q (MeV/c) and the c.m. angle of the spectator to the beam
  /// (deg) ath at [e nAng + k], for the report.
  int nAng = 0;
  std::vector<double> aw, aq, ath;
  /// The same at the spectatorAngle direction (the vertex without a window):
  /// Gd[(e nl + li) 4 + c], qd, pd (fm^-1), thd (deg, between k_sF and k_aA).
  std::vector<double> Gd, qd, pd, thd;
  std::vector<char> valid;  ///< grid energies where the window is reachable
  double buildSeconds = 0.0;
  mutable std::atomic<bool> warned{false};

  /// The vertex at one energy (interpolated): per node its weight and per l
  /// index the two components M^(k)(B) = a[k] (B - 1) - d[k].
  struct At {
    int nodes = 0;
    int nl = 0;
    std::vector<double> weight;
    std::vector<complex> a, d;  ///< [(k nl + li) 2 + comp]
    bool outside = false;
  };

  /*!
   * Builds the tables for E in [eLo, eHi] (MeV, the grid) for the entrance
   * orbital momenta lvals and channel radius a; every energy in `points` (the
   * data) must have a reachable window.  "" or what is wrong.
   */
  std::string Build(const ThmExperiment &x, const ThmDistortion::Kinematics &k, double a,
                    const std::vector<int> &lvals, double eLo, double eHi, const std::vector<double> &points);
  /// The vertex at E, from the grid (end values beyond it, outside set).
  void Evaluate(double energy, At &out) const;
  /// The Gram entries at E: of every node (window) and at the spectatorAngle
  /// direction (gd, nl x 4), with the node weights and q (MeV/c).
  void Interpolate(double energy, std::vector<double> &weight, std::vector<double> &qk, std::vector<double> &g,
                   std::vector<double> &gd, double &qDelta, double &pDelta, bool *outside = nullptr) const;
  /// spectatorAngles: the directions at the grid energy nearest to E (weight,
  /// q in MeV/c, c.m. angle to the beam in deg); empty without the window.
  void AngleNodesAt(double energy, std::vector<double> &weight, std::vector<double> &qk,
                    std::vector<double> &theta) const;
  /// Index of l in lvals, or -1.
  int LIndex(int l) const;
  /// G = L^+ L for g = (g11, g22, Re g12, Im g12): the components
  /// M^(k)(B) = a[k] (B - 1) - d[k], k = 0, 1.
  static void Factor(const double g[4], complex a[2], complex d[2]);
  /// c^+ G c, c = (B - 1, -1): the vertex squared for boundary B.
  static double Vertex2(const double g[4], complex B);
  /// |M^(0)(B)|^2 + |M^(1)(B)|^2 from Factor: what the HOES model uses (equal
  /// to Vertex2 up to the clamp of a slightly negative pivot, which cubic
  /// interpolation of a rank-one G can produce at the 1e-9 level).
  static double Vertex2Factored(const double g[4], complex B);
};

/// Clebsch-Gordan <j1 m1 j2 m2|J M> for integer arguments (long double Racah sum).
double ThmCG(int j1, int m1, int j2, int m2, int J, int M);

#endif
