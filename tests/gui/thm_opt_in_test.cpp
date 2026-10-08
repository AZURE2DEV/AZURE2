// Headless test of the THM opt-in (Configure > Runtime Options > Use Trojan
// Horse Method), which follows the nuclear-potential switch:
//  1. a new project and a classic one (tests/15N_p_a) open with THM off: no
//     THM Workspace entry, no THM tick in the segment dialogs, no THM
//     Background tab; the Runtime Options box is unticked;
//  2. a project opens with THM on when it has THM content: a <thm> block
//     only, a THM segment, a binding energy only, an amplitude width only;
//  3. on and off through the Runtime Options dialog: switching off a project
//     with THM content asks first (Cancel keeps it on) and hides the controls
//     without touching the content; switching on shows them again;
//  4. every save in between is byte for byte the project's own text, classic
//     or THM, and toggling changes no byte;
//  5. Park's parametrization and THM together: ticking Park on a THM project
//     keeps THM on (Brune implied), both stay ticked, saved and reopened;
//     unticking Park restores the configuration and the file.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QDialog>
#include <QFile>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QPushButton>
#include <QTabWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <functional>
#include <iostream>

#include "AZURESetup.h"
#include "ChannelsModel.h"
#include "Config.h"
#include "EditOptionsDialog.h"
#include "FittingTab.h"
#include "PairsModel.h"
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

static QAction* menuAction(AZURESetup& w, const QString& text) {
  for(QAction* a : w.menuBar()->actions())
    if(a->menu())
      for(QAction* b : a->menu()->actions())
        if(b->text().remove('&').contains(text)) return b;
  return nullptr;
}

// Runs `act` with `onModal` called on the first modal dialog it opens.
static void withModal(const std::function<void()>& act, const std::function<void(QWidget*)>& onModal) {
  bool seen = false;
  QTimer t;
  t.setInterval(20);
  QObject::connect(&t, &QTimer::timeout, [&]() {
    QWidget* m = QApplication::activeModalWidget();
    if(!m || seen) return;
    seen = true;
    onModal(m);
  });
  t.start();
  act();
  t.stop();
}

// Is the THM tick of the Add Data Segment dialog shown?  (The dialog is closed.)
static bool segmentTickShown(AZURESetup& w) {
  bool shown = false, found = false;
  withModal([&]() { w.getSegmentsTab()->addSegDataLine(); },
            [&](QWidget* m) {
              for(QCheckBox* c : m->findChildren<QCheckBox*>())
                if(c->text().contains("THM")) {
                  found = true;
                  shown = c->isVisible();
                }
              if(QDialog* d = qobject_cast<QDialog*>(m)) d->reject();
            });
  return found && shown;
}

static bool thmBackgroundTab(AZURESetup& w) {
  for(QTabWidget* t : w.getFittingTab()->findChildren<QTabWidget*>())
    for(int i = 0; i < t->count(); i++)
      if(t->tabText(i).contains("THM")) return true;
  return false;
}

// Everything a THM user sees, off or on together.
static void expectThm(AZURESetup& w, bool on, const QString& where, bool background = false) {
  const QString p = where + ": ";
  ok(qPrintable(p + (on ? "THM on" : "THM off")), w.thmEnabled() == on);
  QAction* ws = menuAction(w, "THM Workspace");
  ok(qPrintable(p + "THM Workspace entry " + (on ? "shown" : "hidden")), ws && ws->isVisible() == on);
  ok(qPrintable(p + "segment dialog THM tick " + (on ? "shown" : "hidden")), segmentTickShown(w) == on);
  if(background || !on)
    ok(qPrintable(p + "THM Background tab " + (on ? "shown" : "hidden")), thmBackgroundTab(w) == on);
}

// Runtime Options: sets the THM box to `on` and accepts; `answer` is the
// button pressed on the warning, if one comes (0: none expected).
static bool runtimeOptions(AZURESetup& w, bool on, QMessageBox::StandardButton answer, bool* boxWas = nullptr) {
  QAction* a = menuAction(w, "Runtime Options");
  if(!a) return false;
  bool warned = false;
  QTimer t;
  t.setInterval(20);
  QObject::connect(&t, &QTimer::timeout, [&]() {
    QWidget* m = QApplication::activeModalWidget();
    if(EditOptionsDialog* d = qobject_cast<EditOptionsDialog*>(m)) {
      if(boxWas) *boxWas = d->useThmCheck->isChecked();
      d->useThmCheck->setChecked(on);
      d->accept();
    } else if(QMessageBox* b = qobject_cast<QMessageBox*>(m)) {
      warned = true;
      b->button(answer ? answer : QMessageBox::Cancel)->click();
    }
  });
  t.start();
  a->trigger();
  t.stop();
  return warned;
}

static QString save(AZURESetup& w) {
  w.saveProject();
  return slurp(QString::fromStdString(w.GetConfig().configfile));
}

int main(int argc, char** argv) {
  // Keep the recent-file list the GUI writes out of the user's settings.
  QTemporaryDir settingsDir;
  qputenv("XDG_CONFIG_HOME", settingsDir.path().toUtf8());
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("AZURE2-tests");
  QCoreApplication::setApplicationName("thm_opt_in_test");

  QTemporaryDir work;
  const QString src = QString(AZURE2_SOURCE_DIR) + "/tests/";
  AZURESetup w;

  // 1. Off for a new project and for a classic one.
  expectThm(w, false, "new project");
  bool boxWas = true;
  runtimeOptions(w, false, QMessageBox::NoButton, &boxWas);
  ok("new project: Runtime Options box unticked", !boxWas);

  // The GUI's own text of a classic project, then its round trips.
  const QString classicPath = work.filePath("classic.azr");
  spit(classicPath, slurp(src + "15N_p_a/15N_p_a.azr"));
  w.open(classicPath);
  const QString classic = save(w);
  ok("classic: saved", classic.contains("</targetInt>"));
  w.open(classicPath);
  expectThm(w, false, "classic");
  ok("classic: open + save byte for byte", save(w) == classic);
  ok("classic: no warning when switched on", !runtimeOptions(w, true, QMessageBox::Ok));
  expectThm(w, true, "classic, switched on");
  ok("classic, switched on: save byte for byte", save(w) == classic);
  ok("classic: no warning when switched off (no THM content)", !runtimeOptions(w, false, QMessageBox::Ok));
  expectThm(w, false, "classic, switched off");
  ok("classic, switched off: save byte for byte", save(w) == classic);
  w.open(classicPath);
  ok("classic: reopened with THM off", !w.thmEnabled());

  // 2. On when the project has THM content, whichever it is.
  {
    // A <thm> block only.
    const QString path = work.filePath("block.azr");
    spit(path, classic + "<thm>\nkinematics=kf3body\n</thm>\n");
    w.open(path);
    const QString saved = save(w);
    ok("<thm> block: THM on", w.thmEnabled());
    ok("<thm> block: content named", w.thmContent().contains("<thm>"), w.thmContent());
    w.open(path);
    ok("<thm> block: open + save byte for byte", save(w) == saved);
  }
  {
    // A binding energy only, set on the classic project and saved.
    w.open(classicPath);
    PairsModel* pairs = w.getPairsTab()->getPairsModel();
    pairs->setData(pairs->index(0, 15), 2.2246, Qt::EditRole);
    const QString path = work.filePath("binding.azr");
    w.saveProject();  // to classic.azr
    spit(path, slurp(classicPath));
    spit(classicPath, classic);
    w.open(path);
    ok("binding energy: THM on", w.thmEnabled() && w.thmContent() == "THM binding energies", w.thmContent());
  }
  {
    // A width entered as an amplitude only.
    w.open(classicPath);
    ChannelsModel* channels = w.getLevelsTab()->getChannelsModel();
    int row = -1;
    for(int i = 0; i < channels->getChannels().size() && row < 0; i++)
      if(channels->getChannels().at(i).radType == QChar('P')) row = i;
    ok("amplitude: a particle channel", row >= 0);
    if(row >= 0) channels->setData(channels->index(row, 7), 1, Qt::EditRole);
    const QString path = work.filePath("amplitude.azr");
    w.saveProject();
    spit(path, slurp(classicPath));
    spit(classicPath, classic);
    w.open(path);
    ok("amplitude width: THM on", w.thmEnabled() && w.thmContent().contains("amplitudes"), w.thmContent());
  }

  // 3. A THM project (tests/7Li_p_a: a THM segment and a binding energy) with a
  //    coherent background, so that the Fitting tab has its THM tab.
  const QString thmPath = work.filePath("thm.azr");
  QString thmIn = slurp(src + "7Li_p_a/7Li_p_a.azr");
  if(!thmIn.endsWith('\n')) thmIn += '\n';
  spit(thmPath, thmIn + "<thm>\nexperiment[E1] segments=1 cbackground=2+:4=0.50,-0.25\n</thm>\n");
  w.open(thmPath);
  const QString thm = save(w);
  w.open(thmPath);
  expectThm(w, true, "THM project", true);
  const QString content = w.thmContent();
  ok("THM project: content named", content.startsWith("a <thm> block, ") && content.contains("THM segment") &&
                                       content.contains("binding energies"), content);
  ok("THM project: open + save byte for byte", save(w) == thm);
  runtimeOptions(w, false, QMessageBox::NoButton, &boxWas);
  ok("THM project: Runtime Options box ticked", boxWas);
  ok("THM project: switching off warns", runtimeOptions(w, false, QMessageBox::Cancel));
  expectThm(w, true, "THM project, Cancel", true);
  ok("THM project: switching off, Ok", runtimeOptions(w, false, QMessageBox::Ok));
  expectThm(w, false, "THM project, switched off");
  ok("THM project, switched off: the content is kept", w.thmContent() == content);
  ok("THM project, switched off: save byte for byte", save(w) == thm);
  ok("THM project, switched off: still a THM segment",
     w.getSegmentsTab()->getSegmentsDataModel()->getLines().value(0).isTHM == 1);
  ok("THM project: no warning when switched on", !runtimeOptions(w, true, QMessageBox::Ok));
  expectThm(w, true, "THM project, switched on again", true);
  ok("THM project, on again: save byte for byte", save(w) == thm);
  w.open(thmPath);
  ok("THM project: reopened with THM on", w.thmEnabled());

  // 5. Park's parametrization and THM together (THM under --use-park).
  {
    auto options = [&](const std::function<void(EditOptionsDialog*)>& set) {
      QAction* a = menuAction(w, "Runtime Options");
      if(!a) return;
      QTimer t;
      t.setInterval(20);
      QObject::connect(&t, &QTimer::timeout, [&]() {
        if(EditOptionsDialog* d = qobject_cast<EditOptionsDialog*>(QApplication::activeModalWidget())) {
          set(d);
          d->accept();
        }
      });
      t.start();
      a->trigger();
      t.stop();
    };
    const unsigned int park = Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM;
    const unsigned int before = w.GetConfig().paramMask;
    bool thmBox = false, bruneBox = false, bruneEnabled = true;
    options([&](EditOptionsDialog* d) {
      d->useParkCheck->setChecked(true);
      thmBox = d->useThmCheck->isChecked();
      bruneBox = d->useBruneCheck->isChecked();
      bruneEnabled = d->useBruneCheck->isEnabled();
    });
    ok("Park ticked on a THM project: the THM box stays ticked, Brune implied", thmBox && bruneBox && !bruneEnabled);
    ok("Park + THM: Park and Brune in the configuration", (w.GetConfig().paramMask & park) == park);
    expectThm(w, true, "Park + THM", true);
    bool parkShown = false;
    thmBox = false;
    options([&](EditOptionsDialog* d) {
      parkShown = d->useParkCheck->isChecked();
      thmBox = d->useThmCheck->isChecked();
    });
    ok("Park + THM: both boxes ticked when the dialog opens again", parkShown && thmBox);
    const QString withPark = save(w);
    w.open(thmPath);
    ok("Park + THM: reopened with both", w.thmEnabled() && (w.GetConfig().paramMask & park) == park);
    ok("Park + THM: open + save byte for byte", save(w) == withPark);
    options([&](EditOptionsDialog* d) { d->useParkCheck->setChecked(false); });
    ok("Park unticked: the configuration as before", w.GetConfig().paramMask == before);
    expectThm(w, true, "Park unticked", true);
    ok("Park unticked: save byte for byte as before Park", save(w) == thm);
  }

  // A segment line's sqrtshift block (dev 1c3e7e3) survives open + save.
  {
    const QString path = work.filePath("sqrtshift.azr");
    spit(path, slurp(src + "energy_shift_sqrt/energy_shift_sqrt.azr"));
    w.open(path);
    const QString saved = save(w);
    ok("sqrtshift: both blocks kept", saved.count("sqrtshift") == 2 && saved.contains("sqrtshift 0.003 0 0") &&
                                          saved.contains("sqrtshift -0.002 0 0"));
    w.open(path);
    ok("sqrtshift: open + save byte for byte", save(w) == saved);
  }

  // A THM project followed by a classic one: THM goes off.
  w.open(classicPath);
  expectThm(w, false, "classic after THM");
  ok("classic after THM: save byte for byte", save(w) == classic);

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
