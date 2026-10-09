// Headless test of the MCMC tab's parameter list (Load Parameters).
//
// The tab shows, for every varying parameter, the prior the sampler will use.
// The automatic priors of the segment normalizations, energy shifts and
// sqrt(E) shift coefficients come from the errors quoted in the segment
// lines; AZURECalcMCMC::BuildAutoPriors() classifies the parameters by the
// form of their names (segment_<k>_norm, _energy_shift, _energy_shift_sqrt,
// cbkg_...).  The tab classified by substring, so a
// segment_<k>_energy_shift_sqrt row was taken for an energy shift and shown
// with the shift's prior (the engine used the right one).
//
// Checked on tests/energy_shift_sqrt with segment 1's norm (5 %), shift
// (0.004 +- 0.002 MeV) and sqrt(E) coefficient (0.003 +- 0.001) free: each
// gets its own kind and prior; the name classification itself on every form.
// With Use Reduced Widths the tab reads param.par: by name, as AZURE2 does
// (a positional read took the #parametrization line for a parameter).
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QDir>
#include <QFile>
#include <QMetaObject>
#include <QProcess>
#include <QString>
#include <QTemporaryDir>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "AZURESetup.h"
#include "Config.h"
#include "MCMCTab.h"
struct SegPairs {int firstPair; int secondPair;};

// Defined by AZURE2.cpp, which belongs to the executable rather than the GUI
// library, so this test supplies its own. They are never called from here.
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

#ifndef AZURE2_SOURCE_DIR
#error "AZURE2_SOURCE_DIR must be defined"
#endif

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

// The project text with CR removed (a checkout may have CRLF line ends).
static QString slurp(const QString& path) {
  QFile f(path);
  if(!f.open(QIODevice::ReadOnly)) return QString();
  return QString::fromUtf8(f.readAll()).remove('\r');
}
static void spit(const QString& path, const QString& text) {
  QFile f(path);
  if(f.open(QIODevice::WriteOnly | QIODevice::Truncate)) f.write(text.toUtf8());
}

int main(int argc, char** argv) {
  QTemporaryDir settingsDir;
  qputenv("XDG_CONFIG_HOME", settingsDir.path().toUtf8());
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("AZURE2-tests");
  QCoreApplication::setApplicationName("mcmc_priors_test");

  // 1. The kind of a parameter from its name, as the engine takes it.
  {
    int key = -1;
    ok("segment_3_norm: norm of segment 3", MCMCTab::categoryOfName("segment_3_norm", false, &key) == "norm" && key == 3);
    ok("segment_12_energy_shift: shift", MCMCTab::categoryOfName("segment_12_energy_shift", false, &key) == "shift" &&
                                             key == 12);
    ok("segment_2_energy_shift_sqrt: the sqrt(E) coefficient",
       MCMCTab::categoryOfName("segment_2_energy_shift_sqrt", false, &key) == "shift_sqrt" && key == 2);
    ok("cbkg_ names are coherent background, whatever they contain",
       MCMCTab::categoryOfName("cbkg_norm_1/2+_2_1/2,0,1/2,0_re0", false) == "cbkg" &&
           MCMCTab::categoryOfName("cbkg_shift_1/2+_2_1/2,0,1/2,0_im1", true) == "cbkg");
    ok("levels", MCMCTab::categoryOfName("energy_1", false) == "level" &&
                     MCMCTab::categoryOfName("width_2_1", true) == "level_rwa");
  }

  // 2. The tab's list for a project with the three segment parameters free.
  QTemporaryDir work;
  const QString src = QString(AZURE2_SOURCE_DIR) + "/tests/energy_shift_sqrt";
  QString text = slurp(src + "/energy_shift_sqrt.azr");
  const QString seg1 = "1 1 1 0 100 40 40 4 1 0 0 0.004 0 0 data/xs_40.dat 0 0 sqrtshift 0.003 0 0";
  ok("found tests/energy_shift_sqrt and its segment 1", text.contains(seg1));
  text.replace(seg1, "1 1 1 0 100 40 40 4 1 1 5 0.004 0.002 1 data/xs_40.dat 0 0 sqrtshift 0.003 0.001 1");
  QDir(work.path()).mkpath("data");
  QDir(work.path()).mkpath("output");
  QDir(work.path()).mkpath("checks");
  for(const QString& f : QDir(src + "/data").entryList(QDir::Files))
    QFile::copy(src + "/data/" + f, work.filePath("data/" + f));
  const QString path = work.filePath("sqrt.azr");
  spit(path, text);

  AZURESetup w;
  w.open(path);
  ok("project opens", QString::fromStdString(w.GetConfig().configfile).endsWith("sqrt.azr"));
  w.GetConfig().paramMask |= Config::CALCULATE_WITH_DATA;
  MCMCTab* tab = w.findChild<MCMCTab*>();
  ok("MCMC tab found", tab != nullptr);
  if(tab) {
    QMetaObject::invokeMethod(tab, "loadFromPhysical");
    int norms = 0, shifts = 0, sqrts = 0;
    bool normPrior = false, shiftPrior = false, sqrtPrior = false;
    for(const MCMCParameter& p : tab->parameters()) {
      if(p.category == "norm") {
        norms++;
        normPrior = p.autoPrior && p.useGaussianPrior && std::fabs(p.priorStd - 0.05 * p.priorMean) < 1e-12 * p.priorMean;
      } else if(p.category == "shift") {
        shifts++;
        shiftPrior = p.autoPrior && p.useGaussianPrior && p.priorMean == 0.004 && p.priorStd == 0.002;
      } else if(p.category == "shift_sqrt") {
        sqrts++;
        sqrtPrior = p.autoPrior && p.useGaussianPrior && p.priorMean == 0.003 && p.priorStd == 0.001;
      }
    }
    ok("one norm, one shift, one sqrt(E) coefficient", norms == 1 && shifts == 1 && sqrts == 1,
       QString("%1 %2 %3 of %4").arg(norms).arg(shifts).arg(sqrts).arg(tab->parameters().size()));
    ok("norm: its 5 % prior", normPrior);
    ok("energy shift: 0.004 +- 0.002 MeV", shiftPrior);
    ok("sqrt(E) coefficient: its own prior, 0.003 +- 0.001", sqrtPrior);
  }

  // 3. Use Reduced Widths: param.par read by name, as AZURE2 reads it.  The
  //    file starts with its #parametrization line (since October 2026); a
  //    positional read took that line for the first parameter and shifted
  //    every value by one row.  The engine's own param.par, the same without
  //    the line (a file of an older AZURE2) and with its rows reversed must
  //    give the same values.
  if(tab) {
    QProcess engine;
    engine.setWorkingDirectory(work.path());
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    if(!env.contains("OMP_NUM_THREADS")) env.insert("OMP_NUM_THREADS", "1");
    engine.setProcessEnvironment(env);
    engine.start(AZURE2_BINARY, QStringList() << "--no-gui" << "--no-readline" << "sqrt.azr");
    engine.waitForStarted(30000);
    engine.write("1\n\n\n7\n");
    engine.closeWriteChannel();
    engine.waitForFinished(300000);
    const QString par = slurp(work.filePath("output/param.par"));
    ok("the engine wrote a tagged param.par", par.startsWith(QString("#parametrization").rightJustified(20)), par.left(80));
    QStringList rows = par.split('\n', Qt::SkipEmptyParts);
    const QString tag = rows.takeFirst();
    QStringList reversed = rows;
    std::reverse(reversed.begin(), reversed.end());
    auto loaded = [&](const QString& fileText) {
      spit(work.filePath("output/param.par"), fileText);
      QMetaObject::invokeMethod(tab, "loadFromReduced");
      QVector<double> v;
      for(const MCMCParameter& p : tab->parameters()) v << p.value;
      return v;
    };
    const QVector<double> asWritten = loaded(par);
    const QVector<double> untagged = loaded(rows.join('\n') + '\n');
    const QVector<double> backwards = loaded(tag + '\n' + reversed.join('\n') + '\n');
    ok("RWA load: some free parameters", !untagged.isEmpty());
    ok("RWA load: the tagged file gives the untagged file's values", asWritten == untagged);
    ok("RWA load: rows in another order, the same values", backwards == untagged);
  }

  std::cout << (fails ? "FAILED" : "PASSED") << " (" << fails << " failure(s))" << std::endl;
  return fails ? 1 : 0;
}
