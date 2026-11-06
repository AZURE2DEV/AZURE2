#ifndef THMDIAGNOSTICS_H
#define THMDIAGNOSTICS_H

#include <QList>
#include <QString>
#include <QThread>
#include <QVector>
#include <limits>

/*!
 * What the Diagnostics page of the THM workspace shows for one THM data
 * segment, computed by the engine (AZUREAPI, the library pyazr drives) for the
 * parameters of the project as it stands, on a grid over the segment's data:
 *
 *   vertex     |M_l(E)|^2 of the entrance vertex at the quasi-free point, per
 *              entrance J^pi group and orbital l, with its nodes;
 *   HOES       the HOES cross section of the segment's entrance -> exit
 *              channel (no resolution, weight or line shape) next to the
 *              on-shell angle-integrated cross section of the same channel;
 *   line shape with lineshape=on: zeta(E) and |N_C|^2 per level
 *              (AZUREAPI::GetThmLineshape);
 *   weight     with weight[k]: w(E) (ThmWeightTable, the engine's reader);
 *   window     with a spectator-momentum window (ps=): the event weight
 *              w(p) = |phi(p)|^2 p^2 over [pmin, pmax] with the engine's
 *              Gauss-Legendre nodes, and in the vertex panel the window
 *              average <|M_l|^2> (AZUREAPI::GetThmVertex, EData::ThmVertexTable);
 *   distortion with distortion=: R(E) as the model uses it and, for coulomb
 *              and optical, |M|^2 and |M_PW|^2 normalized at E_ref
 *              (AZUREAPI::GetThmDistortion, EData::ThmDistortionTable);
 *   angular    with an exit-angle window (theta=): dsigma/dOmega(theta) at one
 *              energy over 0-180 degrees, and its average over the window.
 *              AZUREAPI has no call for the distribution at arbitrary angles,
 *              so the engine runs a copy of the project whose data are one
 *              energy in several single-angle experiments (theta=t-t, the
 *              segment's experiment otherwise: reaction, line shape, window;
 *              no resolution, weight, distortion or background -- the last
 *              three only scale it at one energy).
 *
 * Nothing here is written into the project: the engine runs on a copy in a
 * temporary directory.  ComputeThmDiagnostics does the work and may take a
 * few seconds; ThmDiagnosticsThread runs it off the GUI thread.
 */
struct ThmDiagnosticsRequest {
  QString projectText;  ///< the project as a save would write it
  QString projectDir;   ///< where its relative file names resolve (the engine reads data from the cwd)
  unsigned int paramMask = 0;  ///< Config::paramMask of the GUI (formalism, Coulomb functions, ...)
  int segment = 0;      ///< <segmentsData> line (1-based, counting inactive lines)
  int points = 201;     ///< grid points over the data range
  /// Angular distribution (an experiment with theta=): the c.m. energy (MeV;
  /// NaN: the middle of the data) and the number of angles over 0-180.
  double angularEnergy = std::numeric_limits<double>::quiet_NaN();
  int angularPoints = 19;
};

struct ThmDiagnosticsCurve {
  QString label;
  QVector<double> y;
  QVector<double> yWindow;  ///< vertex: <|M_l|^2> over the ps window (the engine's), or empty
  QVector<double> nodes;  ///< zeros of a real vertex (MeV)
  double boundary = 0.0;  ///< vertex: the real boundary B_c used (constant / perlevel)
  double pole = 0.0, width = 0.0;  ///< line shape: E_lambda (x + A c.m.) and Gamma_lambda (MeV)
};

struct ThmDiagnosticsResult {
  QString error;  ///< "" on success, else what went wrong (with the engine's messages)
  int segment = 0, entranceKey = 0, exitKey = 0;
  QString experiment;  ///< the segment's THM experiment, or ""
  double eLo = 0.0, eHi = 0.0;  ///< c.m. range of the segment's points (MeV)
  QVector<double> energy;       ///< the grid (c.m. of the entrance pair, MeV)

  // Entrance vertex.
  QString vertexMode;           ///< constant | perlevel | onshell
  bool vertexComplex = false;   ///< onshell boundary or Coulomb integral: no real nodes
  double binding = 0.0;         ///< B + spectator energy used in p (MeV)
  struct VertexGroup {
    QString jpi;                ///< e.g. "2+"
    QList<ThmDiagnosticsCurve> curves;  ///< one per entrance l; label "l = 1"
  };
  QList<VertexGroup> vertex;

  // HOES and on-shell (extrapolation run on the grid).
  QVector<double> hoesEnergy, hoes;        ///< HOES model, scaled by hoesScale
  QVector<double> onShellEnergy, onShell;  ///< on-shell sigma (b)
  double hoesScale = 1.0;                  ///< arbitrary THM scale matched to the on-shell curve

  // Line shape (lineshape=on).
  bool lineshape = false;
  QVector<double> zeta;
  QList<ThmDiagnosticsCurve> nc2;  ///< |N_C|^2 of the levels with a pole near the data
  int nc2Hidden = 0;               ///< levels further away (|N_C|^2 flat there)

  // Weight table (weight[k]).
  QString weightFile;
  QVector<double> weight;

  // Spectator-momentum window (ps=), from AZUREAPI::GetThmVertex.
  bool window = false;
  QString windowText;           ///< the engine's description
  double muSx = 0.0, meanTs = 0.0;  ///< MeV (<T_s> at windowE)
  double windowE = 0.0;             ///< the energy of the nodes below (the middle of the grid, MeV)
  QVector<double> nodeP, nodeWeight, nodeTs;  ///< nodes at windowE (MeV/c), normalized weights, T_s (MeV)
  /// The event weight per unit p_s at windowE, A |phi(p)|^2 p over the
  /// accepted p_s (ThmSpectatorWindow::Density), normalized to unit area.
  QVector<double> windowP, windowW;
  QVector<double> nodeW;             ///< the same at the nodes

  // Distortion factor (distortion=), from AZUREAPI::GetThmDistortion.
  bool distortion = false;
  QString distortionKind;         ///< coulomb | optical | table
  QString distortionText;         ///< the engine's description
  QString distortionError;        ///< the engine's reason if it could not be tabulated on the grid
  double distortionRef = 0.0;     ///< E_ref (MeV; coulomb, optical)
  QVector<double> distortionR;    ///< R(E) as the model uses it (interpolated; a table: w(E))
  QVector<double> distortionDirect;  ///< R(E) evaluated directly (coulomb, optical)
  QVector<double> distortionM2, distortionPW2;  ///< |M|^2 and |M_PW|^2 over their values at E_ref
  bool distortionRatioPW = true;  ///< dwpw (R = |M|^2/|M_PW|^2 ...) or dw

  // Angular distribution (theta= window), from single-angle engine runs.
  bool angular = false;
  QString angularError;           ///< the engine's reason if it could not be computed
  double thetaMin = 0.0, thetaMax = 180.0;  ///< the experiment's window (degrees)
  double angularEnergy = 0.0;     ///< c.m. energy of the distribution (MeV)
  QVector<double> angle, dsdo;    ///< theta (degrees) and dsigma/dOmega (HOES, per sr)
  double windowMean = 0.0;        ///< dsigma/dOmega averaged over the window, as the model of a point
};

ThmDiagnosticsResult ComputeThmDiagnostics(const ThmDiagnosticsRequest &request);

class ThmDiagnosticsThread : public QThread {
  Q_OBJECT
 public:
  explicit ThmDiagnosticsThread(const ThmDiagnosticsRequest &request, QObject *parent = 0) :
    QThread(parent), request_(request) {}
  const ThmDiagnosticsResult &result() const { return result_; }

 protected:
  void run() override { result_ = ComputeThmDiagnostics(request_); }

 private:
  ThmDiagnosticsRequest request_;
  ThmDiagnosticsResult result_;
};

#endif
