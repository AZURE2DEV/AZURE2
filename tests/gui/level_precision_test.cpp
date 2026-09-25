// Headless check that saving a project through the GUI leaves its numbers alone.
//
// The GUI rebuilds the .azr from its tabs on every save.  Its writers streamed
// doubles into a QTextStream, whose default is six significant digits, so a
// save rounded the model: tests/7Li_p_a's 55469980 eV width came back as
// 5.547e+07, its -130789.99999999999 as -130790, and the saved project was a
// different evaluation from the one that was opened.  The writers now go
// through roundTripNumber (gui/include/RoundTripNumber.h).
//
// Checked here, on tests/7Li_p_a:
//   - roundTripNumber itself: exact read-back, and six-digit output kept for
//     values that already fit;
//   - open + save leaves every number of every <levels> line exactly equal
//     (both files parsed to doubles; the GUI writes the levels sorted and
//     renumbers the level-id field, so lines are compared as a set with that
//     field left out, and the optional fields 32-33 it adds are read as 0);
//   - the engine's total chi-squared of the saved file equals the original's.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QDir>
#include <QFile>
#include <QProcess>
#include <QRegExp>
#include <QString>
#include <QStringList>
#include <QTemporaryDir>
#include <QTextStream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include "AZURESetup.h"
#include "Config.h"
#include "RoundTripNumber.h"
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
#ifndef AZURE2_BINARY
#error "AZURE2_BINARY must be defined"
#endif

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

static QString slurp(const QString& path) {
  QFile f(path);
  if(!f.open(QIODevice::ReadOnly)) return QString();
  return QString::fromUtf8(f.readAll());
}

// Every non-blank line of <levels>, as doubles, level-id field (9th) removed.
typedef std::vector<double> Row;
static std::vector<Row> levelRows(const QString& text, bool& parsed) {
  std::vector<Row> rows;
  parsed = true;
  const int a = text.indexOf("<levels>"), b = text.indexOf("</levels>");
  if(a < 0 || b < a) { parsed = false; return rows; }
  const QStringList lines = text.mid(a + 8, b - a - 8).split('\n');
  for(const QString& line : lines) {
    const QStringList tok = line.simplified().split(' ', Qt::SkipEmptyParts);
    if(tok.isEmpty()) continue;
    Row row;
    for(int i = 0; i < tok.size(); i++) {
      bool good = false;
      const double v = tok[i].toDouble(&good);
      if(!good) parsed = false;
      if(i != 8) row.push_back(v);
    }
    // Fields 32 and 33 (THM binding energy, RWA flag) are optional, 0 when
    // absent; the GUI always writes them.
    while(row.size() < 32) row.push_back(0.);
    rows.push_back(row);
  }
  std::sort(rows.begin(), rows.end());
  return rows;
}

static void copyDir(const QString& from, const QString& to) {
  QDir().mkpath(to);
  const QStringList files = QDir(from).entryList(QDir::Files);
  for(const QString& f : files) QFile::copy(from + "/" + f, to + "/" + f);
}

// Runs "Calculate Segments From Data" on dir/name and returns the stdout total
// (the one printed with 12 significant digits), or an empty string.
static QString engineChi2(const QString& dir, const QString& name) {
  QDir(dir).mkpath("output");
  QDir(dir).mkpath("checks");
  QProcess p;
  p.setWorkingDirectory(dir);
  p.setProcessChannelMode(QProcess::MergedChannels);
  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  if(!env.contains("OMP_NUM_THREADS")) env.insert("OMP_NUM_THREADS", "1");
  p.setProcessEnvironment(env);
  p.start(AZURE2_BINARY, QStringList() << "--no-gui" << "--no-readline" << name);
  if(!p.waitForStarted(30000)) return QString();
  p.write("1\n\n\n7\n");
  p.closeWriteChannel();
  if(!p.waitForFinished(600000)) { p.kill(); return QString(); }
  const QString out = QString::fromUtf8(p.readAll());
  QRegExp rx("Total Chi-Squared: ([-+0-9.eE]+)");
  return rx.indexIn(out) >= 0 ? rx.cap(1) : QString();
}

int main(int argc, char** argv) {
  QTemporaryDir settingsDir;
  qputenv("XDG_CONFIG_HOME", settingsDir.path().toUtf8());
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("AZURE2-tests");
  QCoreApplication::setApplicationName("level_precision_test");

  // 1. The helper.
  {
    const double values[] = {55469980., -130789.99999999999, 27.494, 0.1, 1. / 3., 5e6,
                             1.00727646688, -1.2e-300, 6.02214076e23, 0.};
    bool allExact = true;
    QString bad;
    for(double v : values)
      if(roundTripNumber(v).toDouble() != v) { allExact = false; bad += roundTripNumber(v) + " "; }
    ok("roundTripNumber reads back exactly", allExact, bad);
    ok("six-digit values keep their six-digit text",
       roundTripNumber(20.2) == "20.2" && roundTripNumber(5e6) == "5e+06" &&
           roundTripNumber(100000.) == "100000" && roundTripNumber(1.00728) == "1.00728",
       roundTripNumber(20.2) + " " + roundTripNumber(5e6) + " " + roundTripNumber(100000.));
    ok("more digits only where needed", roundTripNumber(55469980.) == "55469980" ||
                                            roundTripNumber(55469980.) == "5.546998e+07",
       roundTripNumber(55469980.));
  }

  // 2. Open + save tests/7Li_p_a, compare the <levels> numbers.
  const QString srcDir = QString(AZURE2_SOURCE_DIR) + "/tests/7Li_p_a";
  const QString original = slurp(srcDir + "/7Li_p_a.azr");
  ok("found tests/7Li_p_a", !original.isEmpty());

  QTemporaryDir work;
  const QString origDir = work.filePath("orig"), guiDir = work.filePath("gui");
  copyDir(srcDir + "/data", origDir + "/data");
  copyDir(srcDir + "/data", guiDir + "/data");
  QFile::copy(srcDir + "/7Li_p_a.azr", origDir + "/run.azr");
  QFile::copy(srcDir + "/7Li_p_a.azr", guiDir + "/run.azr");

  {
    AZURESetup w;
    w.open(guiDir + "/run.azr");
    w.saveProject();
  }
  const QString saved = slurp(guiDir + "/run.azr");
  // For a failure: LEVEL_PRECISION_KEEP=path keeps a copy of the saved file.
  if(qEnvironmentVariableIsSet("LEVEL_PRECISION_KEEP"))
    QFile::copy(guiDir + "/run.azr", qEnvironmentVariable("LEVEL_PRECISION_KEEP"));
  ok("GUI wrote the project", saved.contains("</levels>") && saved != original);

  bool p1 = false, p2 = false;
  const std::vector<Row> a = levelRows(original, p1), b = levelRows(saved, p2);
  ok("both <levels> sections parse as numbers", p1 && p2);
  ok("same number of level lines", a.size() == b.size() && !a.empty(),
     QString("%1 vs %2").arg(a.size()).arg(b.size()));
  int differing = 0;
  QString firstDiff;
  for(size_t i = 0; i < std::min(a.size(), b.size()); i++) {
    if(a[i] == b[i]) continue;
    differing++;
    if(firstDiff.isEmpty()) {
      for(size_t k = 0; k < std::min(a[i].size(), b[i].size()); k++)
        if(a[i][k] != b[i][k]) {
          firstDiff = QString("field %1: %2 -> %3").arg(k < 8 ? k + 1 : k + 2)
                          .arg(roundTripNumber(a[i][k])).arg(roundTripNumber(b[i][k]));
          break;
        }
      if(firstDiff.isEmpty()) firstDiff = "different field count";
    }
  }
  ok("every level-line number is exactly the one opened", differing == 0,
     QString("%1 lines differ, first %2").arg(differing).arg(firstDiff));
  ok("including the ones six digits would round", saved.contains(QRegExp("\\s(55469980|5\\.546998e\\+07)\\s")));

  // 3. The engine sees the same project.
  const QString c0 = engineChi2(origDir, "run.azr");
  const QString c1 = engineChi2(guiDir, "run.azr");
  std::cout << "  chi2 original " << c0.toStdString() << "   GUI-saved " << c1.toStdString() << std::endl;
  ok("engine produced a chi-squared for both", !c0.isEmpty() && !c1.isEmpty());
  ok("engine chi-squared of the GUI-saved file equals the original's", !c0.isEmpty() && c0 == c1);

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
