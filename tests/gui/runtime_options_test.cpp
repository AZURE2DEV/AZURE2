// Relativistic kinematics is a Runtime Option (Configure > Runtime Options), not a Particle Pairs control:
//   1. the Runtime Options dialog holds the "Relativistic kinematics" check box, in its own box that says the
//      setting is saved in the project (unlike the other runtime options, the command line and pyazr read it);
//   2. the Particle Pairs tab no longer has it;
//   3. a project with relativistic=1 and relativisticReferencePair=2 reads into the Config and is written back
//      with both keys; switching it off in the Config (what the dialog does) drops both keys from the file;
//      a plain project read afterwards is non-relativistic.
//   4. the prior parameter covariance file (2026-11-05): the dialog's "Prior Parameter Covariance" box holds the
//      file field and its Choose... button; priorCovariance=<path> reads into the Config (comment stripped, a blank
//      inside the path kept), is written back and read by the engine's own <config> reader; cleared, the key is
//      dropped; a plain project has none.
// Before 2026-11-01 the check box sat below the pairs table and the file was written from it, so (1) and (2) and
// the switch-off in (3) fail on older builds.
//
// With AZURE2_OPTIONS_DIALOG_PNG=<file.png> (under xvfb-run -a) the dialog is also saved as an image.
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QCheckBox>
#include <QDir>
#include <QFile>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSettings>
#include <QTemporaryDir>
#include <QTextStream>
#include <iostream>
#include "AZURESetup.h"
#include "Config.h"
#include "EditOptionsDialog.h"
#include "PairsTab.h"
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
bool hasKey(const QString &text, const QString &key) {
  for (const QString &line : text.split('\n'))
    if (line.section('#', 0, 0).trimmed().startsWith(key + "=")) return true;
  return false;
}
}  // namespace

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  QTemporaryDir tmp;
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, tmp.path() + "/settings");
  QSettings::setPath(QSettings::NativeFormat, QSettings::UserScope, tmp.path() + "/settings");
  const QString dir = tmp.path();
  QDir().mkpath(dir + "/output");
  QDir().mkpath(dir + "/checks");

  // 1. The dialog.
  {
    EditOptionsDialog dialog;
    QCheckBox *check = dialog.findChild<QCheckBox *>("relativisticCheck");
    ok("Runtime Options dialog has the Relativistic kinematics check box", check != nullptr);
    QGroupBox *box = dialog.findChild<QGroupBox *>("kinematicsBox");
    ok("it sits in its own Kinematics box", check && box && box->isAncestorOf(check));
    QLabel *note = dialog.findChild<QLabel *>("kinematicsNote");
    ok("the box says it is saved in the project", box && note && box->isAncestorOf(note) && note->text().contains("Saved in the project"));
    const QString png = qEnvironmentVariable("AZURE2_OPTIONS_DIALOG_PNG");
    if (!png.isEmpty()) {
      dialog.adjustSize();
      dialog.show();
      app.processEvents();
      ok("dialog image saved", dialog.grab().save(png));
    }
  }

  // 2. The pairs tab.
  {
    PairsTab pairsTab;
    ok("Particle Pairs tab no longer has it", pairsTab.findChild<QCheckBox *>("relativisticCheck") == nullptr);
  }

  // 3. Through the project file.
  const QString source = readText(QString(AZURE2_TEST_SOURCE_DIR) + "/tests/13N/13N.azr");
  ok("fixture 13N.azr read", source.contains("</config>"));
  QString relativistic = source;
  relativistic.replace("</config>", "relativistic=1\nrelativisticReferencePair=2\n</config>");
  const QString in = dir + "/rel.azr", plain = dir + "/plain.azr", a = dir + "/a.azr", b = dir + "/b.azr", c = dir + "/c.azr";
  ok("fixtures written", writeText(in, relativistic) && writeText(plain, source));
  {
    AZURESetup setup;
    ok("relativistic project read", setup.readProjectFile(in));
    ok("relativistic=1 reaches the Config", setup.GetConfig().relativistic);
    ok("relativisticReferencePair=2 reaches the Config", setup.GetConfig().relativisticReferencePair == 2);
    ok("written", setup.writeProjectFile(a));
    const QString out = readText(a);
    ok("relativistic=1 written back", hasKey(out, "relativistic"));
    ok("relativisticReferencePair=2 written back", out.contains("relativisticReferencePair=2"));
    setup.GetConfig().relativistic = false;  // what unticking it in Runtime Options does
    ok("written after switching it off", setup.writeProjectFile(b));
    const QString off = readText(b);
    ok("switched off: no relativistic key", !hasKey(off, "relativistic"));
    ok("switched off: no relativisticReferencePair key", !hasKey(off, "relativisticReferencePair"));
    ok("plain project read in the same window", setup.readProjectFile(plain));
    ok("plain project is non-relativistic", !setup.GetConfig().relativistic);
    ok("plain project written", setup.writeProjectFile(c));
    ok("plain project written without the key", !hasKey(readText(c), "relativistic"));
  }

  // 4. The prior parameter covariance file (priorCovariance=<path> in <config>, 2026-11-05).
  {
    EditOptionsDialog dialog;
    QLineEdit *edit = dialog.findChild<QLineEdit *>("priorCovarianceEdit");
    QPushButton *choose = dialog.findChild<QPushButton *>("priorCovarianceChoose");
    QGroupBox *box = dialog.findChild<QGroupBox *>("priorCovarianceBox");
    ok("Runtime Options dialog has the prior covariance file field and its Choose... button",
       edit && choose && box && box->isAncestorOf(edit) && box->isAncestorOf(choose) && choose->text() == "Choose...");
  }
  {
    QString withPrior = source;
    withPrior.replace("</config>", "priorCovariance=priors/fit A.txt #a comment\n</config>");
    const QString pin = dir + "/prior.azr", pa = dir + "/pa.azr", pb = dir + "/pb.azr", pc = dir + "/pc.azr";
    ok("prior fixture written", writeText(pin, withPrior));
    AZURESetup setup;
    ok("project with priorCovariance read", setup.readProjectFile(pin));
    ok("priorCovariance reaches the Config (comment stripped, inner blank kept)",
       setup.GetConfig().priorCovariance == "priors/fit A.txt");
    ok("written", setup.writeProjectFile(pa));
    const QString out = readText(pa);
    ok("priorCovariance written back", out.contains("priorCovariance=priors/fit A.txt #"));
    Config check(std::cout);
    check.configfile = pa.toStdString();
    ok("the engine's reader reads the written key", check.ReadConfigFile() == 0 && check.priorCovariance == "priors/fit A.txt");
    setup.GetConfig().priorCovariance.clear();  // what clearing the field in Runtime Options does
    ok("written after clearing it", setup.writeProjectFile(pb));
    ok("cleared: no priorCovariance key", !hasKey(readText(pb), "priorCovariance"));
    ok("plain project read in the same window", setup.readProjectFile(plain));
    ok("plain project has no prior covariance", setup.GetConfig().priorCovariance.empty());
    ok("plain project written", setup.writeProjectFile(pc));
    ok("plain project written without the key", !hasKey(readText(pc), "priorCovariance"));
  }

  if (failures) {
    std::cout << failures << " failure(s)" << std::endl;
    return 1;
  }
  std::cout << "all passed" << std::endl;
  return 0;
}
