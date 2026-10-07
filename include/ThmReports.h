#ifndef THM_REPORTS_H
#define THM_REPORTS_H

#include <memory>
#include <string>
#include <vector>

struct ThmSpectatorWindow;

/*!
 * What EData's THM tables (ThmVertexTable, ThmLineshapeTable,
 * ThmDistortionTable) report, for thm_experiments.out, pyazr (thm_vertex,
 * thm_lineshape, thm_distortion) and the GUI diagnostics: plain data, apart
 * from the physics objects that compute them.
 */

/// The window-averaged entrance vertex of one experiment (EData::ThmVertexTable).
struct ThmVertexReport {
  std::string experiment;
  std::string window;  ///< description ("delta" without a window)
  int pairKey = 0;
  double muSx = 0.0;   ///< MeV (0 without a window)
  double bind = 0.0;   ///< B_xs as the vertex uses it (the entrance pair's), MeV
  double radius = 0.0; ///< channel radius a (fm)
  std::vector<double> energy;         ///< E (MeV)
  /// The nodes at each E, [E][node]: |p_s| (MeV/c), normalized weight, T_s
  /// (MeV) and the c.m. angle of the spectator to the beam (deg).  A delta:
  /// one node, p = 0, T = spectatorEnergy (theta empty).  With the DW vertex
  /// the plane-wave nodes are not used (empty).
  std::vector<std::vector<double>> p, weight, es, theta;
  std::vector<std::vector<double>> rho;  ///< [E][node] p_xA a / hbar c
  /// 1 where the vertex at E is E's own, 0 where the window (or the DW grid)
  /// does not reach E and the nodes and vertex are those of the nearest data
  /// point or grid energy, as for a folding sub-point (ThmSpectatorWindow::Nodes,
  /// ThmDwVertex::Interpolate).  Always 1 for a delta experiment.
  std::vector<char> reached;
  /// The plane-wave window itself (null without one), e.g. for its Density.
  std::shared_ptr<const ThmSpectatorWindow> windowObject;
  /// "pw" (the plane-wave vertex) or "dw" (vertexModel=dw, ThmDwVertex.h).
  /// For dw, m2 is the DW vertex averaged over the nodes the DW vertex uses
  /// (the reachable part of the window at each E), m2qf the DW vertex at the
  /// spectatorAngle direction, and the nodes are per energy:
  std::string model = "pw";
  std::vector<std::vector<double>> dwQ, dwWeight;  ///< [E][node] q (MeV/c), weight
  std::vector<double> dwQDelta, dwPDelta;          ///< spectatorAngle node: q (MeV/c), p (fm^-1)
  /// spectatorAngles=: the accepted directions at the grid energy nearest to
  /// each E, [E][node]: the c.m. angle of the spectator to the beam (deg), q
  /// (MeV/c) and the normalized weight (d cos theta_cm x acceptance x
  /// |phi~(q)|^2); dwQ/dwWeight are then the one node of the averaged vertex.
  std::vector<std::vector<double>> dwTheta, dwAngleQ, dwAngleWeight;
  struct Level {
    int level = 0;          ///< 1-based in the J group
    double boundary = 0.0;  ///< Re of the vertex boundary B used for this level
    std::vector<double> m2;   ///< <|M_l|^2> over the window, on the grid
    std::vector<double> m2qf; ///< |M_l|^2 at p_s = 0 (quasi-free), on the grid
    /// vertexModel=dw: the plane-wave |M_l|^2 at the spectatorAngle node's
    /// p = |k_aA - alpha k_sF| (the DW vertex without distortion); empty for pw.
    std::vector<double> m2pw;
  };
  struct Channel {
    int jgroup = 0, channel = 0;  ///< 1-based
    double J = 0.0;
    int pi = 1;
    int l = 0;
    double s = 0.0;
    std::vector<Level> levels;
  };
  std::vector<Channel> channels;
};

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

/// What the output file and pyazr report for one experiment's distortion.
struct ThmDistortionReport {
  std::string experiment, kind, description;
  double eRef = 0.0, eAA = 0.0, bind = 0.0, kAA = 0.0, etaAA = 0.0, kappa = 0.0, etaB = 0.0, beta = 0.0;
  std::vector<double> energy, esf, ksf, etasf, thetaCm, x, q, m2, mpw2, r, rModel;
  std::vector<int> lmax;
};

#endif
