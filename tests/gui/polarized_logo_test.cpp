// The logo follows the project (2026-11-06): the plain AZURE2 logo by default, the one with the polarization
// vector over the A while any data or test segment is a polarization observable (analyzing power, data type 7).
//   - a project without one (tests/shift_fn): the logo button shows the plain logo (tests/13N has six
//     analyzing-power segments, 12C(p,p) A_y, so it shows the polarized one);
//   - a project with one (tests/13N_capture_ay): the polarized logo, also in the About box;
//   - removing its analyzing-power segment switches back; reading the plain project again too.
#include <QApplication>
#include <QDialog>
#include <QIcon>
#include <QImage>
#include <QLabel>
#include <QPushButton>
#include <QSettings>
#include <QTemporaryDir>
#include <QTimer>
#include <QToolButton>
#include <iostream>
#include "AZURESetup.h"
#include "Config.h"
#include "SegmentsDataModel.h"

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
QImage image(const QIcon &icon) { return icon.pixmap(64, 64).toImage().convertToFormat(QImage::Format_ARGB32); }
// Which logo the button shows: "plain", "polarized" or "other".
std::string shown(AZURESetup &setup) {
  QToolButton *b = setup.findChild<QToolButton *>("logoButton");
  if (!b) return "no logo button";
  const QImage im = image(b->icon());
  const QImage plain = image(QIcon(":/azure-icon.png")), pol = image(QIcon(":/azure-icon-polarized.png"));
  if (!pol.isNull() && im == pol) return "polarized";
  if (!plain.isNull() && im == plain) return "plain";
  return "other";
}
// The About box's text, opened from the logo button (closed again by a timer).
QString aboutText(AZURESetup &setup) {
  QString text;
  QTimer::singleShot(200, [&text]() {
    for (QWidget *w : QApplication::topLevelWidgets())
      if (QDialog *d = qobject_cast<QDialog *>(w))
        if (d->isVisible() && d->windowTitle().contains("About")) {
          for (QLabel *l : d->findChildren<QLabel *>()) text += l->text();
          d->accept();
        }
  });
  if (QToolButton *b = setup.findChild<QToolButton *>("logoButton")) b->click();
  return text;
}
}  // namespace

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  QTemporaryDir tmp;
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, tmp.path() + "/settings");
  const QString src = QString(AZURE2_TEST_SOURCE_DIR) + "/tests/";

  AZURESetup setup;
  ok("a new window shows the plain logo", shown(setup) == "plain");
  ok("tests/shift_fn (no polarization segment) read", setup.readProjectFile(src + "shift_fn/shift_fn.azr"));
  ok("  plain logo", shown(setup) == "plain");
  ok("tests/13N_capture_ay (an analyzing-power segment) read", setup.readProjectFile(src + "13N_capture_ay/13N_capture_ay.azr"));
  ok("  polarized logo on the logo button", shown(setup) == "polarized");
  ok("  polarized logo in the About box", aboutText(setup).contains("azure-icon-polarized.png"));
  SegmentsDataModel *data = setup.getSegmentsTab()->getSegmentsDataModel();
  int ay = -1;
  for (int i = 0; i < data->getLines().size(); i++)
    if (data->getLines().at(i).dataType == 7) ay = i;
  ok("  its analyzing-power segment found", ay >= 0);
  if (ay >= 0) data->removeRows(ay, 1);
  ok("removing that segment brings back the plain logo", shown(setup) == "plain");
  ok("  and the plain logo in the About box", aboutText(setup).contains("azure-icon.png") &&
                                              !aboutText(setup).contains("azure-icon-polarized.png"));
  ok("the analyzing-power project again: polarized", setup.readProjectFile(src + "13N_capture_ay/13N_capture_ay.azr") && shown(setup) == "polarized");
  ok("tests/shift_fn again: plain", setup.readProjectFile(src + "shift_fn/shift_fn.azr") && shown(setup) == "plain");
  ok("tests/13N (six 12C(p,p) A_y segments): polarized", setup.readProjectFile(src + "13N/13N.azr") && shown(setup) == "polarized");

  if (failures) {
    std::cout << failures << " failure(s)" << std::endl;
    return 1;
  }
  std::cout << "all passed" << std::endl;
  return 0;
}
