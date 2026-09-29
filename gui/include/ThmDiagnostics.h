#ifndef THMDIAGNOSTICS_H
#define THMDIAGNOSTICS_H

#include <QList>
#include <QString>
#include <QThread>
#include <QVector>

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
 *              average <|M_l|^2> (AZUREAPI::GetThmVertex, EData::ThmVertexTable).
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
  double muSx = 0.0, meanTs = 0.0;  ///< MeV
  QVector<double> nodeP, nodeWeight, nodeTs;  ///< nodes (MeV/c), normalized weights, T_s (MeV)
  QVector<double> windowP, windowW;  ///< w(p) on [pmin, pmax], normalized to unit area (per MeV/c)
  QVector<double> nodeW;             ///< the same w(p) at the nodes
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
