// The GUI warns when a project is opened whose last run used the GSL Coulomb functions (2026-11-06).
// Those functions were a speed compromise; they are wrong by up to several percent at some (eta, rho) and give
// a NaN shift just above charged-particle thresholds.  About 1,100 archive projects have the option set in
// <lastRun> (bit Config::USE_GSL_COULOMB_FUNC), and the GUI silently used it again on loading them.
//
//   - readProjectFile() (headless) shows nothing;
//   - open() shows a non-modal QMessageBox "gslCoulombWarning" naming the option; a plain project shows none;
//   - "Use the default functions" clears the option, and the project is written without it;
//   - "Keep GSL" keeps it.
#include <QApplication>
#include <QDir>
#include <QFile>
#include <QMessageBox>
#include <QPushButton>
#include <QSettings>
#include <QTemporaryDir>
#include <QTextStream>
#include <iostream>
#include "AZURESetup.h"
#include "Config.h"

struct SegPairs {int firstPair; int secondPair;};

// Defined by AZURE2.cpp, which belongs to the executable rather than the GUI library, so this test supplies its own.
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

namespace {
int failures = 0;
void ok(const std::string &what, bool cond) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what << std::endl;
  if (!cond) failures++;
}
QString readText(const QString &path) {
  QFile f(path);
  if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) return QString();
  return QTextStream(&f).readAll();
}
bool writeText(const QString &path, const QString &text) {
  QFile f(path);
  if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
  QTextStream(&f) << text;
  return true;
}
// The run mask: the first line after <lastRun>.
unsigned lastRunMask(const QString &text) {
  const int i = text.indexOf("<lastRun>");
  if (i < 0) return 0;
  return text.mid(i + 9).trimmed().section('\n', 0, 0).trimmed().toUInt();
}
QMessageBox *warning(AZURESetup &setup) {
  QMessageBox *box = setup.findChild<QMessageBox *>("gslCoulombWarning");
  return (box && box->isVisible()) ? box : nullptr;
}
}  // namespace

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  QTemporaryDir tmp;
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, tmp.path() + "/settings");
  const QString dir = tmp.path();
  QDir().mkpath(dir + "/output");
  QDir().mkpath(dir + "/checks");

  const QString source = readText(QString(AZURE2_TEST_SOURCE_DIR) + "/tests/13N/13N.azr");
  const unsigned plainMask = lastRunMask(source);
  ok("fixture 13N.azr read, its last run without GSL Coulomb functions",
     source.contains("<lastRun>") && plainMask != 0 && !(plainMask & Config::USE_GSL_COULOMB_FUNC));
  QString gslText = source;
  const unsigned gslMask = plainMask | Config::USE_GSL_COULOMB_FUNC;
  gslText.replace(QString("<lastRun>\n%1\n").arg(plainMask), QString("<lastRun>\n%1\n").arg(gslMask));
  const QString gsl = dir + "/gsl.azr", plain = dir + "/plain.azr", switched = dir + "/switched.azr", kept = dir + "/kept.azr";
  ok("fixtures written", lastRunMask(gslText) == gslMask && writeText(gsl, gslText) && writeText(plain, source));

  {  // headless reading
    AZURESetup setup;
    ok("readProjectFile reads the GSL project and its option",
       setup.readProjectFile(gsl) && (setup.GetConfig().paramMask & Config::USE_GSL_COULOMB_FUNC));
    app.processEvents();
    ok("readProjectFile shows no dialog", warning(setup) == nullptr);
  }
  {  // open(): the warning, "Use the default functions"
    AZURESetup setup;
    setup.open(gsl);
    app.processEvents();
    QMessageBox *box = warning(setup);
    ok("open() on the GSL project shows the warning", box != nullptr);
    ok("the warning is non-modal", box && !box->isModal());
    ok("it names the option and the alternative",
       box && box->text().contains("GSL Coulomb functions") && box->text().contains("Runtime Options") &&
           box->text().contains("--gsl-coul"));
    QPushButton *useDefault = box ? box->findChild<QPushButton *>("gslCoulombUseDefault") : nullptr;
    ok("it offers \"Use the default functions\"", useDefault != nullptr);
    if (useDefault) useDefault->click();
    app.processEvents();
    ok("clicking it clears the option", !(setup.GetConfig().paramMask & Config::USE_GSL_COULOMB_FUNC));
    ok("and closes the warning", warning(setup) == nullptr);
    ok("the project is written without the option",
       setup.writeProjectFile(switched) && !(lastRunMask(readText(switched)) & Config::USE_GSL_COULOMB_FUNC));
    setup.open(switched);
    app.processEvents();
    ok("reopening it shows no warning", warning(setup) == nullptr);
  }
  {  // "Keep GSL"
    AZURESetup setup;
    setup.open(gsl);
    app.processEvents();
    QMessageBox *box = warning(setup);
    QPushButton *keep = box ? box->findChild<QPushButton *>("gslCoulombKeep") : nullptr;
    ok("it offers \"Keep GSL\"", keep != nullptr);
    if (keep) keep->click();
    app.processEvents();
    ok("Keep GSL keeps the option", (setup.GetConfig().paramMask & Config::USE_GSL_COULOMB_FUNC) != 0);
    ok("and it is written back", setup.writeProjectFile(kept) && (lastRunMask(readText(kept)) & Config::USE_GSL_COULOMB_FUNC));
  }
  {  // a plain project
    AZURESetup setup;
    setup.open(plain);
    app.processEvents();
    ok("open() on the plain project shows no warning", warning(setup) == nullptr);
  }

  if (failures) {
    std::cout << failures << " failure(s)" << std::endl;
    return 1;
  }
  std::cout << "all passed" << std::endl;
  return 0;
}
