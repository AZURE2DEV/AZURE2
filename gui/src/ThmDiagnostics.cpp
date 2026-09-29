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
#include <memory>
#include <sstream>

#include "ALevel.h"
#include "AChannel.h"
#include "AZUREAPI.h"
#include "CNuc.h"
#include "ChannelFunc.h"
#include "Config.h"
#include "EData.h"
#include "ESegment.h"
#include "EPoint.h"
#include "JGroup.h"
#include "PPair.h"
#include "ThmExperiment.h"
#include "ThmFunc.h"
#include "ThmLineshape.h"

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

/// weight[k]= / weightTest[k]= paths made absolute (the copy of the project
/// is not in the project's directory, and the engine resolves them against it).
QStringList absoluteWeights(const QStringList &lines, const QString &dir) {
  QStringList out;
  QRegularExpression rx("^(\\s*weight(?:Test)?\\s*\\[\\s*\\d+\\s*\\]\\s*=\\s*)([^#\\s][^#]*?)(\\s*(#.*)?)$");
  for (const QString &line : lines) {
    QRegularExpressionMatch m = rx.match(line);
    if (m.hasMatch() && QFileInfo(m.captured(2)).isRelative())
      out << m.captured(1) + QDir(dir).absoluteFilePath(m.captured(2)) + m.captured(3);
    else
      out << line;
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
  const QString oldCwd = QDir::currentPath();
  QDir::setCurrent(request.projectDir);
  struct RestoreCwd {
    QString dir;
    ~RestoreCwd() { QDir::setCurrent(dir); }
  } restore{oldCwd};

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
  const bool brune = !!(cfg.paramMask & Config::USE_BRUNE_FORMALISM);

  // Entrance vertex M_l = (B_c - 1) j_l(pa) - pa j_l'(pa) [+ C_l] at the
  // quasi-free point (EPoint::CalcEDependentValues), with the boundary
  // THMMatrixFunc::CalculateTHMCrossSection chooses.
  const int aa = compound->GetPairNumFromKey(r.entranceKey);
  PPair *pair = compound->GetPair(aa);
  const double threshold = pair->GetSepE() + pair->GetExE();
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
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    ALevel *lowest = nullptr;
    for (int la = 1; la <= jg->NumLevels(); la++) {
      ALevel *level = jg->GetLevel(la);
      if (level->IsInRMatrix() && (!lowest || level->GetFitE() < lowest->GetFitE())) lowest = level;
    }
    if (!lowest) continue;
    ThmDiagnosticsResult::VertexGroup group;
    group.jpi = jpiText(jg->GetJ(), jg->GetPi());
    QList<int> ls;
    for (int ch = 1; ch <= jg->NumChannels(); ch++) {
      AChannel *c = jg->GetChannel(ch);
      if (c->GetPairNum() != aa || pair->GetPType() != 0 || ls.contains(c->GetL())) continue;
      const int l = c->GetL();
      ls << l;
      double b = 0.0;
      if (cfg.thm.vertex == Config::ThmOptions::CONSTANT)
        b = ChannelFunc(pair, useGSL).Shift(l, lowest->GetFitE() - threshold);
      else if (!onShell)
        b = brune ? lowest->GetShiftFunction(ch) : c->GetBoundaryCondition();
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
        complex boundary(b, 0.0);
        if (onShell) {
          ChannelFunc f(pair, useGSL);
          boundary = complex(f.Shift(l, e), f.Penetrability(l, e));
        }
        complex m = (boundary - 1.0) * jl - rhoDjl;
        if (coulomb) m += ThmCoulombTerm(pair, l, e, ThmRho(mu, e, bind, 1.0), useGSL);
        curve.y << std::norm(m);
        re[i] = m.real();
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
  for (const ThmExperiment &x : cfg.thm.experiments)
    if (std::find(x.segments.begin(), x.segments.end(), request.segment) != x.segments.end()) {
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
  return r;
}
