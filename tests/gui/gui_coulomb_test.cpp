// The Coulomb table of this computer in the window (docs/upgrade/qt_gui.md section 19.12), without a display and
// without this machine's own table: HOME is an empty folder of the build tree and AZURE2_COULOMB_TABLE is unset, as on
// CI (the file runs look for is then <HOME>/.azure2/coulomb_table.bin, the engine's own rule).
//   gui_coulomb_test <source dir> <screenshot dir>
//  (a) no table: the start page and the Overview each show the quiet line with a Build... button (not red); the Run
//      menu has Build the Coulomb Table... and Check the Coulomb Table... in standard view; Tools' entry opens the same
//      dialog; the line's button opens the dialog (not modal);
//  (b) a file there that is no table: the line says it is not used, with the loader's reason;
//  (c) no build starts while a calculation of the window runs;
//  (d) a build (a table of 2000 rows, never the full 88 MB one) advances its progress bar, is cancelled: no file at
//      the path nor at <path>.new, the line still shown;
//  (e) a small table (l <= 2, eta 0.1-3, rho 0.1-10, 30 nodes per decade, as ctest coulomb_table_small) built through
//      the facade from the dialog: the file is written, progress reached 100 %, the engine's loader takes it although
//      it found no file before (CoulombTable::ForgetMissing), both lines are gone and the next calculation of T1 in
//      this process takes values from it (it took none before); the threads of a build are the command line's
//      default and never above 4.
// PASS/FAIL per check; exit 1 on a failure.
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <QAction>
#include <QApplication>
#include <QDir>
#include <QElapsedTimer>
#include <QLabel>
#include <QMenu>
#include <QMenuBar>
#include <QProgressBar>
#include <QPushButton>
#include <QSettings>

#include "AZURESetup.h"
#include "CoulombTableUi.h"
#include "OverviewPage.h"
#include "PythonTools.h"
#include "StartPage.h"

namespace fs = std::filesystem;
namespace front = azure2::front;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
  if (!ok) failures++;
}

void Settle(int ms = 200) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < ms) QApplication::processEvents(QEventLoop::AllEvents, 20);
}

/// Process events until `done` or `ms` passed; true when done.
template <typename F>
bool WaitFor(F done, int ms) {
  QElapsedTimer t;
  t.start();
  while (!done() && t.elapsed() < ms) QApplication::processEvents(QEventLoop::AllEvents, 20);
  return done();
}

/// The notice on the start page and the one on the Overview.
CoulombTableNotice *NoticeIn(QWidget *page) {
  const auto all = page->findChildren<CoulombTableNotice *>();
  return all.isEmpty() ? nullptr : all.first();
}

bool Shown(QWidget *w) { return w && !w->isHidden(); }

/// The colour a label is drawn in is not a red (the line is a hint, not an error).
bool NotRed(QLabel *l) {
  const QColor c = l->palette().color(l->foregroundRole());
  return !(c.red() > c.green() + 60 && c.red() > c.blue() + 60);
}

QAction *MenuAction(QMainWindow *w, const QString &menuTitle, const QString &name) {
  for (QAction *m : w->menuBar()->actions())
    if (m->menu() && m->text() == menuTitle)
      for (QAction *a : m->menu()->actions())
        if (a->objectName() == name) return a;
  return nullptr;
}

struct Quiet : front::RunObserver {
  void Message(const std::string &, bool) override {}
};

/// A calculation of the project through the facade (the window's run); what the engine wrote to stderr, where a run
/// that took values from the Coulomb table says "Coulomb table: N lookups, M outside the table".
std::string CalculateStderr(const fs::path &project) {
  for (const char *d : {"output", "checks"}) fs::create_directories(project.parent_path() / d);
  front::RunRequest req;
  req.command = "calculate";
  req.project = project.string();
  req.threads = 2;
  Quiet quiet;
  std::ostringstream err;
  std::streambuf *old = std::cerr.rdbuf(err.rdbuf());
  front::Run(req, quiet, nullptr);
  std::cerr.rdbuf(old);
  return err.str();
}

bool UsedTable(const std::string &log) { return log.find(" lookups, ") != std::string::npos; }

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: gui_coulomb_test <source dir> <screenshot dir>" << std::endl;
    return 2;
  }
  // Before anything reads the environment (the engine's record is made once, at its first use): an empty HOME of our
  // own and no AZURE2_COULOMB_TABLE.  This machine's ~/.azure2 is never looked at.
  const fs::path work = fs::absolute("gui_coulomb_work");
  fs::remove_all(work);
  const fs::path home = work / "home";
  fs::create_directories(home);
  qputenv("HOME", QByteArray::fromStdString(home.string()));
  qputenv("USERPROFILE", QByteArray::fromStdString(home.string()));
  qunsetenv("AZURE2_COULOMB_TABLE");

  QApplication app(argc, argv);
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, true);
  const fs::path source = argv[1];
  const QString shots = QString::fromLocal8Bit(argv[2]);
  QDir().mkpath(shots);
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, QString::fromStdString((work / "settings").string()));

  const fs::path table = home / ".azure2" / "coulomb_table.bin";
  {
    const front::CoulombTableState s = front::CoulombTableStatus();
    Check(fs::path(s.path) == table, "the table runs look for is <HOME>/.azure2/coulomb_table.bin: " + s.path);
    Check(!s.exists && !s.usable, "no table in the empty HOME");
  }

  AZURESetup w;
  w.resize(1280, 800);
  w.show();
  Settle();
  CoulombTableControl *control = w.coulombTable();
  Check(control != nullptr, "the window has its Coulomb table control");

  // (a) the line on the start page and on the Overview, the menu entries
  CoulombTableNotice *start = NoticeIn(w.startPage());
  Check(Shown(start), "(a) the start page shows the Coulomb table line");
  if (start) {
    Check(start->label()->text().contains("Coulomb table not built") && start->label()->text().contains("88 MB"),
          "(a) its words: " + start->label()->text().toStdString());
    Check(start->label()->property("role").toString() == "hint" && NotRed(start->label()),
          "(a) the line is a hint, not drawn in red");
    Check(start->buildButton()->isEnabled() && start->buildButton()->text().startsWith("Build"),
          "(a) its Build... button is there and enabled");
  }
  w.grab().save(shots + "/coulomb_start_page.png");
  const fs::path t1 = work / "T1";
  fs::copy(source / "tests/upgrade/T1", t1, fs::copy_options::recursive | fs::copy_options::skip_symlinks);
  QString why;
  Check(w.readProjectFile(QString::fromStdString((t1 / "T1.azr").string()), why), "(a) T1 opens: " + why.toStdString());
  Settle(400);
  OverviewPage *overview = nullptr;
  for (int i = 0; i < w.tabs()->count() && !overview; i++) overview = qobject_cast<OverviewPage *>(w.tabs()->widget(i));
  CoulombTableNotice *onOverview = overview ? NoticeIn(overview) : nullptr;
  Check(Shown(onOverview), "(a) the Overview shows the Coulomb table line");
  const std::string before = CalculateStderr(t1 / "T1.azr");
  Check(!UsedTable(before), "(a) a calculation of T1 runs without a table");
  w.grab().save(shots + "/coulomb_overview.png");

  Check(!w.expertMode(), "(a) standard view");
  QAction *build = MenuAction(&w, QObject::tr("&Run"), "buildCoulombTable");
  QAction *checkEntry = MenuAction(&w, QObject::tr("&Run"), "checkCoulombTable");
  Check(build && build->isVisible() && QString(build->text()).remove('&').contains("Build the Coulomb Table"),
        "(a) Run > Build the Coulomb Table... in standard view");
  Check(checkEntry && checkEntry->isVisible() && QString(checkEntry->text()).remove('&').contains("Check the Coulomb Table"),
        "(a) Run > Check the Coulomb Table... in standard view");
  if (checkEntry) checkEntry->trigger();
  Settle();
  CoulombTableDialog *dialog = control->dialog();
  Check(dialog && dialog->isVisible() && !dialog->isModal(), "(a) Check opens the table's dialog, not modal");
  if (dialog) {
    Check(dialog->stateLabel()->text().contains(QDir::toNativeSeparators(QString::fromStdString(table.string()))) &&
              dialog->stateLabel()->text().contains("None there yet"),
          "(a) the dialog names the path and says there is none: " + dialog->stateLabel()->text().toStdString());
    dialog->grab().save(shots + "/coulomb_dialog.png");
    dialog->hide();
  }
  emit w.tools()->coulombTableRequested();
  Settle();
  Check(dialog && dialog->isVisible() && control->dialog() == dialog, "(a) Tools' entry opens the same dialog");
  if (dialog) dialog->hide();
  if (onOverview) onOverview->buildButton()->click();
  Settle();
  Check(dialog && dialog->isVisible(), "(a) the line's Build... opens the dialog");
  if (!dialog) {
    std::cout << "FAILED" << std::endl;
    return 1;
  }
  const int threads = dialog->buildThreads();
  Check(threads >= 1 && threads <= 4 && threads == std::min(front::DefaultCoulombTableBuild().threads, 4),
        "(e) a build's threads: the command line's default, at most 4 (" + std::to_string(threads) + ")");

  // (b) a file that is no table
  fs::create_directories(table.parent_path());
  { std::ofstream(table) << "not a table\n"; }
  control->refresh();
  Settle();
  Check(Shown(start) && start->label()->text().contains("not used") && start->label()->text().contains("not a Coulomb table"),
        "(b) a file that is no table: the line says why: " + (start ? start->label()->text().toStdString() : std::string()));
  fs::remove(table);
  control->refresh();
  Settle();
  Check(Shown(start) && start->label()->text().contains("not built"), "(b) removed: the line says not built again");

  // (c) no build while a calculation runs
  control->setCalculationRunning(true);
  Check(!dialog->buildButton()->isEnabled() && !start->buildButton()->isEnabled(),
        "(c) while a calculation runs, Build is disabled in the dialog and on the line");
  why.clear();
  const bool refused = !dialog->startBuild(&why);
  Check(refused && why.contains("calculation"), "(c) startBuild refuses: " + why.toStdString());
  control->setCalculationRunning(false);
  Check(dialog->buildButton()->isEnabled(), "(c) enabled again after the calculation");

  // (d) a build advances, is cancelled, leaves nothing
  front::CoulombTableBuild big = front::DefaultCoulombTableBuild();
  big.lMax = 3;
  big.perDecade = 100;  // 4 x 501 rows of 501 nodes: minutes to the end, a second to its first per cent
  big.threads = 2;
  dialog->setSpec(big);
  Check(dialog->startBuild(&why), "(d) a build starts");
  Check(dialog->building() && !start->buildButton()->isEnabled(), "(d) while it runs the line's Build... is disabled");
  const bool advanced = WaitFor([&]() { return dialog->progressBar()->value() >= 1; }, 60000);
  Check(advanced && dialog->progressBar()->isVisible(),
        "(d) the progress bar advances (" + std::to_string(dialog->progressBar()->value()) + " %)");
  dialog->grab().save(shots + "/coulomb_building.png");
  bool cancelOk = true;
  QString cancelError;
  QObject::connect(dialog, &CoulombTableDialog::buildFinished, [&](bool ok, const QString &e) {
    cancelOk = ok;
    cancelError = e;
  });
  dialog->cancelButton()->click();
  Check(WaitFor([&]() { return !dialog->building(); }, 60000), "(d) Cancel ends the build");
  Check(!cancelOk && cancelError == "cancelled", "(d) it ends as cancelled: " + cancelError.toStdString());
  Check(!fs::exists(table) && !fs::exists(table.string() + ".new"), "(d) no file at the path nor at <path>.new");
  Check(Shown(start) && Shown(onOverview), "(d) the lines are still shown");

  // (e) a small table built from the dialog: written, taken by the loader, the lines gone
  front::CoulombTableBuild small = front::DefaultCoulombTableBuild();
  small.lMax = 2;
  small.etaMin = 0.1;
  small.etaMax = 3;
  small.rhoMin = 0.1;
  small.rhoMax = 10;
  small.perDecade = 30;
  small.threads = 2;
  dialog->setSpec(small);
  int lastValue = -1;
  bool okEnd = false;
  QObject::disconnect(dialog, &CoulombTableDialog::buildFinished, nullptr, nullptr);
  QObject::connect(dialog, &CoulombTableDialog::buildFinished, [&](bool ok, const QString &) {
    okEnd = ok;
    lastValue = dialog->progressBar()->value();
  });
  dialog->buildButton()->click();
  Check(WaitFor([&]() { return lastValue >= 0; }, 120000), "(e) the small build ends");
  Check(okEnd && lastValue == 100, "(e) it succeeds, progress at 100 %");
  Check(fs::is_regular_file(table) && fs::file_size(table) > 1000 && !fs::exists(table.string() + ".new"),
        "(e) the file is written at the path (no .new left)");
  const front::CoulombTableState after = front::CoulombTableStatus();
  Check(after.usable && after.lMax == 2,
        "(e) the engine's loader takes it although it found none before (l up to " + std::to_string(after.lMax) + ")");
  Settle();
  Check(!Shown(start) && !Shown(onOverview), "(e) both lines are gone");
  const std::string afterRun = CalculateStderr(t1 / "T1.azr");
  Check(UsedTable(afterRun), "(e) the next calculation of T1 takes values from the new table (same process)");
  Check(dialog->stateLabel()->text().contains("Found and used"), "(e) the dialog says it is found and used");
  dialog->grab().save(shots + "/coulomb_built.png");
  w.grab().save(shots + "/coulomb_overview_after.png");

  fs::remove_all(work);
  std::cout << (failures ? "FAILED" : "all passed") << std::endl;
  return failures ? 1 : 0;
}
