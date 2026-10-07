#include "ThmDiagnostics.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QRegularExpression>
#include <QStringList>
#include <QTemporaryDir>
#include <QTextStream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <sstream>

#include "ALevel.h"
#include "AChannel.h"
#include "AZUREAPI.h"
#include "Constants.h"
#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "ESegment.h"
#include "EPoint.h"
#include "JGroup.h"
#include "NuclearPotentialManager.h"
#include "PPair.h"
#include "ThmExperiment.h"
#include "ThmFunc.h"
#include "ThmLineshape.h"
#include "ThmReports.h"
#include "ThmVertexBoundary.h"

namespace {

/// Replaces the body of <tag> ... </tag> (lines of their own); false if absent.
bool replaceBlock(QString &text, const QString &tag, const QString &body) {
  const QString open = "<" + tag + ">\n", close = "</" + tag + ">";
  int a = text.startsWith(open) ? 0 : text.indexOf("\n" + open);
  if (a < 0) return false;
  if (a > 0) a += 1;
  a += open.size();
  const int b = text.indexOf("\n" + close, a - 1);
  if (b < 0) return false;
  text.replace(a, b + 1 - a, body);
  return true;
}

/// The lines of the <thm> block (without the tags), or none.
QStringList thmLines(const QString &text, int *from = nullptr, int *to = nullptr) {
  QRegularExpression rx("(^|\\n)[ \\t]*<thm>[^\\n]*\\n");
  QRegularExpressionMatch m = rx.match(text);
  if (!m.hasMatch()) return QStringList();
  const int a = m.capturedEnd();
  const int b = text.indexOf(QRegularExpression("\\n[ \\t]*</thm>"), a - 1);
  if (b < 0) return QStringList();
  if (from) *from = a;
  if (to) *to = b + 1;
  return text.mid(a, b + 1 - a).split('\n', Qt::SkipEmptyParts);
}

void setThmLines(QString &text, const QStringList &lines) {
  int a = -1, b = -1;
  thmLines(text, &a, &b);
  if (a < 0) return;
  text.replace(a, b - a, lines.isEmpty() ? QString() : lines.join('\n') + '\n');
}

/// weight[k]= / weightTest[k]= paths, and ps=table:<file> and
/// distortion=table:<file> of experiment lines, made absolute (the copy of the
/// project is not in the project's directory, and the engine resolves them
/// against it).
QStringList absoluteWeights(const QStringList &lines, const QString &dir) {
  QStringList out;
  QRegularExpression rx("^(\\s*weight(?:Test)?\\s*\\[\\s*\\d+\\s*\\]\\s*=\\s*)([^#\\s][^#]*?)(\\s*(#.*)?)$");
  QRegularExpression table("(^|[ \\t])(?:ps|distortion|spectatorAngles)=(?:cm:)?table:([^ \\t#]+)");
  for (const QString &line : lines) {
    QRegularExpressionMatch m = rx.match(line);
    if (m.hasMatch() && QFileInfo(m.captured(2)).isRelative()) {
      out << m.captured(1) + QDir(dir).absoluteFilePath(m.captured(2)) + m.captured(3);
      continue;
    }
    const int hash = line.indexOf('#');
    const QString code = hash < 0 ? line : line.left(hash);
    QString changed = line;
    if (code.trimmed().startsWith("experiment[")) {
      // From the last match back, so that the earlier positions stay valid.
      QList<QRegularExpressionMatch> matches;
      for (QRegularExpressionMatchIterator it = table.globalMatch(code); it.hasNext();) matches.prepend(it.next());
      for (const QRegularExpressionMatch &t : matches)
        if (QFileInfo(t.captured(2)).isRelative())
          changed.replace(t.capturedStart(2), t.capturedLength(2), QDir(dir).absoluteFilePath(t.captured(2)));
    }
    out << changed;
  }
  return out;
}

QString jpiText(double J, int pi) {
  const int twice = (int)std::lround(2.0 * J);
  const QString j = twice % 2 ? QString("%1/2").arg(twice) : QString::number(twice / 2);
  return j + (pi < 0 ? "-" : "+");
}

QString lastLines(const std::string &log, int n = 12) {
  QStringList lines = QString::fromStdString(log).split('\n', Qt::SkipEmptyParts);
  while (lines.size() > n) lines.removeFirst();
  return lines.join('\n');
}

/// A Config for the copy at `file`, with the GUI's formalism flags.
struct Engine {
  std::ostringstream log;
  Config config;
  std::unique_ptr<AZUREAPI> api;
  Engine() : config(log) {}
  QString start(const QString &file, const QString &outDir, unsigned int mask, bool withData) {
    // The copy's <potential> block is the GUI's own NuclearPotentialManager
    // written out (AZURESetup::writeProject).  Read back, Config::
    // ReadPotentialBlock would reset and refill that process-wide singleton
    // from this worker thread while the GUI thread may read it (the
    // workspace's other pages stay live during a computation).  So the block
    // is taken out, its useAdaptiveGrid kept, and the engine uses the
    // manager as it is -- the same settings, never written here.
    bool adaptiveGrid = true;
    {
      QFile f(file);
      if (!f.open(QIODevice::ReadOnly)) return QObject::tr("Cannot read %1.").arg(file);
      QString text = QString::fromUtf8(f.readAll());
      f.close();
      QRegularExpression block("(^|\\n)[ \\t]*<potential>[^\\n]*\\n(.*?\\n)?[ \\t]*</potential>[^\\n]*(\\n|$)",
                               QRegularExpression::DotMatchesEverythingOption);
      QRegularExpressionMatch m = block.match(text);
      if (m.hasMatch()) {
        QRegularExpressionMatch grid = QRegularExpression("(^|\\n)[ \\t]*useAdaptiveGrid[ \\t]*=[ \\t]*(\\d+)").match(m.captured(0));
        if (grid.hasMatch()) adaptiveGrid = grid.captured(2).toInt() == 1;
        text.replace(m.capturedStart(), m.capturedLength(), m.captured(1));
        if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) return QObject::tr("Cannot write %1.").arg(file);
        f.write(text.toUtf8());
        f.close();
      }
    }
    config.configfile = QDir::toNativeSeparators(file).toStdString();
    config.paramMask = mask;
    config.paramMask &= ~(Config::PERFORM_FIT | Config::PERFORM_ERROR_ANALYSIS | Config::CALCULATE_REACTION_RATE |
                          Config::CALCULATE_COVARIANCE_BAND | Config::USE_PREVIOUS_PARAMETERS |
                          Config::USE_PREVIOUS_INTEGRALS | Config::CALCULATE_WITH_DATA);
    if (withData) config.paramMask |= Config::CALCULATE_WITH_DATA;
    config.paramMask |= Config::USE_API;
    const int status = config.ReadConfigFile();
    if (status == -1) return QObject::tr("AZURE2 could not read the copy of the project.");
    if (status < 0) return QObject::tr("AZURE2 refuses the project:\n%1").arg(lastLines(log.str()));
    // What ReadPotentialBlock sets from the block.
    config.useAdaptiveGrid = adaptiveGrid;
    config.useHybridMethod = NuclearPotentialManager::instance().isAnyEnabled();
    // Nothing of this run goes to the project's output or checks directories.
    config.outputdir = QDir::toNativeSeparators(outDir).toStdString() + "/";
    config.checkdir = config.outputdir;
    config.screenCheckMask = 0;
    config.fileCheckMask = 0;
    api.reset(new AZUREAPI(config));
    if (api->Initialize() != 0) return QObject::tr("AZURE2 could not set up the project:\n%1").arg(lastLines(log.str()));
    return QString();
  }
};

/// The <segmentsData> line `key` (1-based, every line counted) of a project text.
QString segmentsDataLine(const QString &text, int key) {
  const int a = text.indexOf("<segmentsData>\n");
  const int b = text.indexOf("</segmentsData>", a);
  if (a < 0 || b < 0) return QString();
  const QStringList lines = text.mid(a + 15, b - a - 15).split('\n', Qt::SkipEmptyParts);
  int k = 0;
  for (const QString &line : lines)
    if (!line.trimmed().isEmpty() && ++k == key) return line;
  return QString();
}

/*!
 * dsigma/dOmega(theta) at one lab energy of the segment's entrance pair: a
 * copy of the project whose data are that energy (two points: an experiment
 * needs more points than profiled parameters) in one THM segment per angle,
 * each in its own experiment with theta=t-t and the segment's experiment's
 * other keys (reaction, line shape, spectator window); the last segment has
 * the experiment's own window.  No folding, weight, distortion or background:
 * at one energy they only scale the curve.  Fills r.angle, r.dsdo and
 * r.windowMean; returns "" or the reason.
 */
QString angularDistribution(const ThmDiagnosticsRequest &request, const ThmExperiment &experiment, double eLab,
                            const QString &dir, ThmDiagnosticsResult &r) {
  if (experiment.vertexDW)
    return QObject::tr("the angular distribution is the fixed-angle observable (theta=), which is not available "
                       "with vertexModel=dw.");
  QString text = request.projectText;
  const QStringList fields = segmentsDataLine(text, request.segment).split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
  if (fields.size() < 8) return QObject::tr("segment %1 not found in <segmentsData>.").arg(request.segment);
  // The experiment's keys but those the copy sets or leaves out.
  const QStringList dropped = {"segments", "background", "theta", "distortion", "opticalAA", "opticalSF",
                               "spectatorAngle", "distortionRef", "distortionRatio", "boundState",
                               "spectatorAngles", "spectatorAngleNodes"};
  QStringList keep, thm;
  for (const QString &line : absoluteWeights(thmLines(text), request.projectDir)) {
    const QString code = line.left(line.indexOf('#')).trimmed();
    if (code.startsWith("weight")) continue;
    if (!code.startsWith("experiment[")) {
      thm << line;
      continue;
    }
    if (code.mid(11, code.indexOf(']') - 11) != QString::fromStdString(experiment.name)) continue;
    for (const QString &token : code.mid(code.indexOf(']') + 1).split(QRegularExpression("[ \t]+"), Qt::SkipEmptyParts))
      if (!dropped.contains(token.left(token.indexOf('=')))) keep << token;
  }
  const int n = std::max(2, request.angularPoints);
  const int segments = n + 1;
  QString data, dataLines;
  QTextStream ds(&data), ls(&dataLines);
  ds.setRealNumberPrecision(17);
  ls.setRealNumberPrecision(17);
  const QString file = QDir(dir).filePath("angular.dat");
  if (file.contains(QRegularExpression("\\s"))) return QObject::tr("The temporary directory has a space in its path.");
  ds << eLab << " 90 1 0.1\n" << eLab << " 90 1 0.1\n";
  ds.flush();
  const double width = std::max(1.0e-6, 1.0e-9 * std::fabs(eLab));
  for (int k = 1; k <= segments; k++) {
    ls << "1 " << fields[1] << " " << fields[2] << " " << eLab - width << " " << eLab + width << " 0 180 " << fields[7]
       << " 1 1 0 0 0 0 " << file << "\n";
    const QString window = k <= n ? QString("%1-%1").arg(180.0 * (k - 1) / (n - 1), 0, 'g', 17)
                                  : QString("%1-%2").arg(experiment.thetaMin, 0, 'g', 17).arg(experiment.thetaMax, 0, 'g', 17);
    thm << QString("experiment[a%1] segments=%1 theta=%2 %3").arg(k).arg(window, keep.join(' ')).trimmed();
  }
  ls.flush();
  {
    QFile f(file);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) return QObject::tr("Cannot write %1.").arg(file);
    f.write(data.toUtf8());
  }
  setThmLines(text, thm);
  if (!replaceBlock(text, "segmentsData", dataLines) || !replaceBlock(text, "segmentsTest", QString()) ||
      !replaceBlock(text, "targetInt", QString()))
    return QObject::tr("The project has no <segmentsData>, <segmentsTest> or <targetInt> block.");
  const QString projectFile = QDir(dir).filePath("diagnostics_angular.azr");
  {
    QFile f(projectFile);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) return QObject::tr("Cannot write %1.").arg(projectFile);
    f.write(text.toUtf8());
  }
  QDir(dir).mkpath("angular_run");
  Engine engine;
  const QString why = engine.start(projectFile, QDir(dir).filePath("angular_run"), request.paramMask, true);
  if (!why.isEmpty()) return why;
  vector_r p = engine.api->params_values_rwa();
  if (engine.api->UpdateSegmentsRWA(p) != segments)
    return QObject::tr("AZURE2 could not evaluate the angles:\n%1").arg(lastLines(engine.log.str()));
  for (int k = 0; k < segments; k++) {
    const vector_r v = engine.api->calculated_segments(k), e = engine.api->calculated_energies(k);
    if (v.empty() || e.empty()) return QObject::tr("AZURE2 could not evaluate the angles:\n%1").arg(lastLines(engine.log.str()));
    if (k == 0) r.angularEnergy = e[0];
    if (k < n) {
      r.angle << 180.0 * k / (n - 1);
      r.dsdo << v[0];
    } else {
      r.windowMean = v[0];
    }
  }
  return QString();
}

}  // namespace

ThmDiagnosticsResult ComputeThmDiagnostics(const ThmDiagnosticsRequest &request) {
  ThmDiagnosticsResult r;
  r.segment = request.segment;
  QTemporaryDir tmp;
  if (!tmp.isValid()) {
    r.error = QObject::tr("No temporary directory for the engine run.");
    return r;
  }
  // The engine reads data files relative to the working directory, as a run
  // from the main window does (the project's directory).
  // The main window keeps the working directory there (AZURESetup::readFile
  // and writeFile), so this is normally a no-op: the process-wide directory
  // is only touched if it is somewhere else.
  const QString oldCwd = QDir::currentPath();
  const bool moved = QDir(oldCwd) != QDir(request.projectDir);
  if (moved) QDir::setCurrent(request.projectDir);
  struct RestoreCwd {
    QString dir;
    bool moved;
    ~RestoreCwd() {
      if (moved) QDir::setCurrent(dir);
    }
  } restore{oldCwd, moved};

  // 1. The project as it is, with data: the segment's points, the compound
  //    nucleus at the current parameters, the line shape and the weight.
  QString text = request.projectText;
  setThmLines(text, absoluteWeights(thmLines(text), request.projectDir));
  const QString dataFile = tmp.filePath("diagnostics.azr");
  {
    QFile f(dataFile);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
      r.error = QObject::tr("Cannot write %1.").arg(dataFile);
      return r;
    }
    f.write(text.toUtf8());
  }
  QDir(tmp.path()).mkpath("data_run");
  Engine data;
  QString why = data.start(dataFile, tmp.filePath("data_run"), request.paramMask, true);
  if (!why.isEmpty()) {
    r.error = why;
    return r;
  }
  AZUREAPI &api = *data.api;
  CNuc *compound = api.compound();
  ESegment *segment = nullptr;
  for (ESegment &s : api.data()->GetSegments())
    if (s.GetSegmentKey() == request.segment) segment = &s;
  if (!segment || !segment->IsTHM() || segment->GetPoints().empty()) {
    r.error = QObject::tr("Data segment %1 is not an active THM segment with points.").arg(request.segment);
    return r;
  }
  r.entranceKey = segment->GetEntranceKey();
  r.exitKey = segment->GetExitKey();
  r.eLo = 1e300;
  r.eHi = -1e300;
  for (EPoint &p : segment->GetPoints()) {
    r.eLo = std::min(r.eLo, p.GetCMEnergy());
    r.eHi = std::max(r.eHi, p.GetCMEnergy());
  }
  double lo = r.eLo, hi = r.eHi;
  if (!(hi > lo)) {
    lo -= 0.05;
    hi += 0.05;
  }
  const int n = std::max(2, request.points);
  std::vector<double> grid(n);
  for (int i = 0; i < n; i++) grid[i] = lo + (hi - lo) * i / (n - 1);
  r.energy = QVector<double>(grid.begin(), grid.end());
  const Config &cfg = data.config;
  const bool useGSL = !!(cfg.paramMask & Config::USE_GSL_COULOMB_FUNC);

  // Entrance vertex M_l = (B_c - 1) j_l(pa) - pa j_l'(pa) [+ C_l] at the
  // quasi-free point (EPoint::CalcEDependentValues), with the boundary
  // THMMatrixFunc::CalculateTHMCrossSection chooses.
  const int aa = compound->GetPairNumFromKey(r.entranceKey);
  PPair *pair = compound->GetPair(aa);
  const double mu = pair->GetRedMass() * uconv;
  const double bind = pair->GetBindingEnergy() + cfg.thm.SpectatorEnergy(r.entranceKey);
  const double radius = pair->GetChRad();
  const bool coulomb = cfg.thm.coulombIntegral && pair->GetZ(1) * pair->GetZ(2) != 0;
  r.binding = bind;
  r.vertexMode = cfg.thm.vertex == Config::ThmOptions::CONSTANT   ? "constant"
                 : cfg.thm.vertex == Config::ThmOptions::ON_SHELL ? "onshell"
                                                                  : "perlevel";
  const bool onShell = cfg.thm.vertex == Config::ThmOptions::ON_SHELL;
  r.vertexComplex = onShell || (coulomb && hi > 0.0);

  // The segment's experiment and its spectator-momentum window: the nodes and
  // <|M_l|^2> over them from the engine (EData::ThmVertexTable).
  const ThmExperiment *experiment = nullptr;
  for (const ThmExperiment &x : cfg.thm.experiments)
    if (std::find(x.segments.begin(), x.segments.end(), request.segment) != x.segments.end()) experiment = &x;
  ThmVertexReport windowReport;
  if (experiment && experiment->psKind != ThmExperiment::PS_DELTA) {
    std::string whyNot;
    if (!api.GetThmVertex(experiment->name, grid, windowReport, whyNot)) {
      r.error = QObject::tr("Spectator-momentum window: %1").arg(QString::fromStdString(whyNot));
      return r;
    }
    // The plane-wave window (the DW vertex averages over the directions itself).
    const ThmSpectatorWindow *w = windowReport.windowObject.get();
    if (w) {
      r.window = true;
      r.windowText = QString::fromStdString(windowReport.window);
      r.muSx = windowReport.muSx;
      // The nodes follow the accepted directions, which depend on E: those at
      // the middle of the grid, and the event weight per unit p_s there,
      // A |phi(p)|^2 p (d cos theta_cm = p dp/(beta k_sF k_aA)), unit area.
      const int mid = n / 2;
      r.windowE = grid[mid];
      r.nodeP = QVector<double>(windowReport.p[mid].begin(), windowReport.p[mid].end());
      r.nodeWeight = QVector<double>(windowReport.weight[mid].begin(), windowReport.weight[mid].end());
      r.nodeTs = QVector<double>(windowReport.es[mid].begin(), windowReport.es[mid].end());
      for (int k = 0; k < r.nodeTs.size(); k++) r.meanTs += r.nodeWeight[k] * r.nodeTs[k];
      double qLo = 0.0, qHi = 0.0;
      if (w->Reach(r.windowE, qLo, qHi)) {
        const double p0 = std::max(qLo, w->pMin), p1 = std::min(qHi, w->pMax);
        if (p1 > p0) {
          const int m = 401;
          double area = 0.0;
          for (int i = 0; i < m; i++) {
            const double p = p0 + (p1 - p0) * i / (m - 1);
            r.windowP << p;
            r.windowW << w->Density(r.windowE, p);
            if (i > 0) area += 0.5 * (r.windowW[i] + r.windowW[i - 1]) * (p1 - p0) / (m - 1);
          }
          if (area > 0.0) {
            for (double &v : r.windowW) v /= area;
            for (double p : r.nodeP) r.nodeW << w->Density(r.windowE, p) / area;
          } else {
            r.windowP.clear();
            r.windowW.clear();
          }
        }
      }
    }
  }

  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    ThmDiagnosticsResult::VertexGroup group;
    group.jpi = jpiText(jg->GetJ(), jg->GetPi());
    QList<int> ls;
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      AChannel *c = jg->GetChannel(ch);
      if (c->GetPairNum() != aa || pair->GetPType() != 0 || ls.contains(c->GetL())) continue;
      // The boundary the model uses, for the lowest level of the J group.
      const ThmVertexBoundary boundaryOf(cfg, pair, jg, ch);
      const int lowestIndex = boundaryOf.LowestLevel();
      if (!lowestIndex) break;  // no level in the R matrix
      const int l = c->GetL();
      ls << l;
      const double b = onShell ? 0.0 : boundaryOf.Level(jg->GetLevel(lowestIndex));
      ThmDiagnosticsCurve curve;
      curve.label = QString("l = %1").arg(l);
      curve.boundary = b;
      std::vector<double> re(n);
      for (int i = 0; i < n; i++) {
        const double e = grid[i];
        if (!(e + bind > 0.0)) {
          curve.y << 0.0;
          re[i] = 0.0;
          continue;
        }
        double jl, rhoDjl;
        ThmBesselParts(l, mu, e, bind, radius, jl, rhoDjl);
        const complex boundary = onShell ? boundaryOf.OnShellAt(e) : complex(b, 0.0);
        complex m = (boundary - 1.0) * jl - rhoDjl;
        if (coulomb) m += ThmCoulombTerm(pair, l, e, ThmRho(mu, e, bind, 1.0), useGSL);
        curve.y << std::norm(m);
        re[i] = m.real();
      }
      if (r.window) {
        // The engine's window average for this channel and the level whose boundary is used.
        for (const ThmVertexReport::Channel &rc : windowReport.channels)
          if (rc.jgroup == j && rc.channel == ch)
            for (const ThmVertexReport::Level &lv : rc.levels)
              if (lv.level == lowestIndex) curve.yWindow = QVector<double>(lv.m2.begin(), lv.m2.end());
      }
      if (!r.vertexComplex) {
        // Nodes: sign changes on the grid, refined by bisection of the same M_l.
        auto value = [&](double e) {
          double jl, rhoDjl;
          ThmBesselParts(l, mu, e, bind, radius, jl, rhoDjl);
          double m = (b - 1.0) * jl - rhoDjl;
          if (coulomb) m += ThmCoulombTerm(pair, l, e, ThmRho(mu, e, bind, 1.0), useGSL).real();
          return m;
        };
        for (int i = 1; i < n; i++) {
          if (!(grid[i - 1] + bind > 0.0) || re[i - 1] == 0.0 || (re[i - 1] > 0.0) == (re[i] > 0.0)) continue;
          double a = grid[i - 1], z = grid[i], fa = re[i - 1];
          for (int it = 0; it < 60 && z - a > 1e-12; it++) {
            const double mid = 0.5 * (a + z), fm = value(mid);
            if ((fm > 0.0) == (fa > 0.0)) {
              a = mid;
              fa = fm;
            } else {
              z = mid;
            }
          }
          curve.nodes << 0.5 * (a + z);
        }
      }
      group.curves << curve;
    }
    if (!group.curves.isEmpty()) r.vertex << group;
  }

  // Line shape of the segment's experiment, if on.
  if (experiment) {
    const ThmExperiment &x = *experiment;
    {
      r.experiment = QString::fromStdString(x.name);
      if (x.lineshape) {
        ThmLineshapeReport report;
        std::string whyNot;
        if (!api.GetThmLineshape(x.name, grid, report, whyNot)) {
          r.error = QObject::tr("Line shape: %1").arg(QString::fromStdString(whyNot));
          return r;
        }
        for (const ThmLineshapeReport::Exit &e : report.exits) {
          if (e.pairKey != r.exitKey) continue;
          r.lineshape = true;
          r.zeta = QVector<double>(e.zeta.begin(), e.zeta.end());
          // The levels with a pole inside the data range: elsewhere |N_C|^2
          // is nearly flat, exp(+-pi zeta), and only crowds the plot.
          QList<const ThmLineshapeReport::Level *> near;
          for (const ThmLineshapeReport::Level &lv : e.levels)
            if (lv.energy >= lo && lv.energy <= hi)
              near << &lv;
            else
              r.nc2Hidden++;
          std::sort(near.begin(), near.end(), [](const ThmLineshapeReport::Level *a, const ThmLineshapeReport::Level *b) {
            return a->energy < b->energy;
          });
          for (int k = 0; k < near.size(); k++) {
            if (k >= 4) {
              r.nc2Hidden++;
              continue;
            }
            ThmDiagnosticsCurve c;
            c.label = QString("%1, %2 MeV").arg(jpiText(near[k]->J, near[k]->pi)).arg(near[k]->energy, 0, 'f', 2);
            c.y = QVector<double>(near[k]->nc2.begin(), near[k]->nc2.end());
            c.pole = near[k]->energy;
            c.width = near[k]->width;
            r.nc2 << c;
          }
        }
      }
    }
  }

  // Distortion factor of the segment's experiment (EData::ThmDistortionTable).
  if (experiment && experiment->distortion != ThmExperiment::DIST_NONE) {
    r.distortion = true;
    r.distortionRatioPW = experiment->distortionRatioPW;
    ThmDistortionReport report;
    std::string whyNot;
    if (!api.GetThmDistortion(experiment->name, grid, report, whyNot)) {
      r.distortionError = QString::fromStdString(whyNot);
    } else {
      r.distortionKind = QString::fromStdString(report.kind);
      r.distortionText = QString::fromStdString(report.description);
      r.distortionR = QVector<double>(report.rModel.begin(), report.rModel.end());
      if (report.kind != "table") {
        r.distortionRef = report.eRef;
        r.distortionDirect = QVector<double>(report.r.begin(), report.r.end());
        ThmDistortionReport ref;
        if (api.GetThmDistortion(experiment->name, std::vector<double>(1, report.eRef), ref, whyNot) &&
            !ref.m2.empty() && ref.m2[0] > 0.0 && ref.mpw2[0] > 0.0)
          for (size_t i = 0; i < report.m2.size(); i++) {
            r.distortionM2 << report.m2[i] / ref.m2[0];
            r.distortionPW2 << report.mpw2[i] / ref.mpw2[0];
          }
      }
    }
  }

  // Weight table of the segment.
  auto w = cfg.thm.weightBySegment.find(request.segment);
  if (w != cfg.thm.weightBySegment.end() && w->second) {
    r.weightFile = QString::fromStdString(w->second->name);
    for (double e : grid) r.weight << (*w->second)(e);
  }

  // 2. HOES and on-shell cross sections on the grid: the same project with
  //    two extrapolation segments of the segment's channel and nothing else
  //    (no resolution, no weight, no experiment, hence no line shape).
  text = request.projectText;
  QStringList kept;
  for (const QString &line : thmLines(text)) {
    const QString code = line.left(line.indexOf('#')).trimmed();
    if (code.startsWith("weight") || code.startsWith("experiment[")) continue;
    kept << line;
  }
  setThmLines(text, kept);
  const double toLab = (pair->GetM(1) + pair->GetM(2)) / pair->GetM(2);
  const double step = (hi - lo) / (n - 1) * toLab;
  QString tests;
  QTextStream ts(&tests);
  ts.setRealNumberPrecision(12);
  ts << "1 " << r.entranceKey << " " << r.exitKey << " " << lo * toLab << " " << hi * toLab << " " << step
     << " 0 0 0 10\n";
  const bool withOnShell = r.entranceKey != r.exitKey;
  if (withOnShell)
    ts << "1 " << r.entranceKey << " " << r.exitKey << " " << lo * toLab << " " << hi * toLab << " " << step
       << " 0 0 0 0\n";
  ts.flush();
  if (!replaceBlock(text, "segmentsTest", tests) || !replaceBlock(text, "targetInt", QString())) {
    r.error = QObject::tr("The project has no <segmentsTest> or <targetInt> block.");
    return r;
  }
  const QString extrapFile = tmp.filePath("diagnostics_extrap.azr");
  {
    QFile f(extrapFile);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
      r.error = QObject::tr("Cannot write %1.").arg(extrapFile);
      return r;
    }
    f.write(text.toUtf8());
  }
  QDir(tmp.path()).mkpath("extrap_run");
  Engine extrap;
  why = extrap.start(extrapFile, tmp.filePath("extrap_run"), request.paramMask, false);
  if (!why.isEmpty()) {
    r.error = why;
    return r;
  }
  vector_r p = extrap.api->params_values_rwa();
  if (extrap.api->UpdateSegmentsRWA(p) < (withOnShell ? 2 : 1)) {
    r.error = QObject::tr("AZURE2 could not evaluate the extrapolation:\n%1").arg(lastLines(extrap.log.str()));
    return r;
  }
  const vector_r he = extrap.api->calculated_energies(0), hx = extrap.api->calculated_segments(0);
  r.hoesEnergy = QVector<double>(he.begin(), he.end());
  r.hoes = QVector<double>(hx.begin(), hx.end());
  if (withOnShell) {
    const vector_r oe = extrap.api->calculated_energies(1), ox = extrap.api->calculated_segments(1);
    r.onShellEnergy = QVector<double>(oe.begin(), oe.end());
    r.onShell = QVector<double>(ox.begin(), ox.end());
    // The THM scale is arbitrary: match it to the on-shell curve (geometric mean of the ratio).
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < r.hoes.size() && i < r.onShell.size(); i++)
      if (r.hoes[i] > 0.0 && r.onShell[i] > 0.0) {
        sum += std::log(r.onShell[i] / r.hoes[i]);
        count++;
      }
    r.hoesScale = count ? std::exp(sum / count) : 1.0;
    for (double &v : r.hoes) v *= r.hoesScale;
  }

  // 3. Angular distribution of an experiment with a theta window.
  if (experiment && experiment->hasTheta) {
    r.angular = true;
    r.thetaMin = experiment->thetaMin;
    r.thetaMax = experiment->thetaMax;
    // By default the middle of the data, rounded (to the power of ten below a
    // twentieth of their range, the step of the energy box).
    const double step = std::pow(10.0, std::floor(std::log10(std::max(r.eHi - r.eLo, 1e-3) / 20.0)));
    r.angularEnergy = std::isfinite(request.angularEnergy) ? request.angularEnergy
                                                           : std::round(0.5 * (r.eLo + r.eHi) / step) * step;
    r.angularError = angularDistribution(request, *experiment, r.angularEnergy * toLab, tmp.path(), r);
  }
  return r;
}
