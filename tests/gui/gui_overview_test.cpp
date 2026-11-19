// The Overview page without a display (docs/upgrade/qt_gui.md section 15):
//   gui_overview_test <source dir> <AZURE2 binary> <screenshot dir>
//  (a) T1 and T2 opened: the Overview is the first tab and is shown; its cards' numbers are the facade's
//      (SummarizeProject: pairs, levels per J^pi, data sets, basis); without outputs "No run yet" with its two buttons;
//  (b) a calculation run through the window: the page shows the run's total chi^2 and one row per data set with the
//      chi^2 of the run (RunOutcome), the free parameters of the worker (ReadRunParameters), chi^2/(N - p);
//  (c) a click on a data set's row opens the Plot tab with that data set; a step's link goes to its tab;
//  (d) an edit is not stale until saved; saved (the Save action), the stale mark appears;
//  (e) screenshots light and dark at 1400 x 900; at that size T2's page needs no vertical scrolling;
//  (f) changes on disk (section 16): the window's own save is no change; the project edited outside without unsaved
//      edits is re-read (tab kept, the new value in the tab); with unsaved edits the bar, both sides kept, the
//      differences, Keep mine; an outside run's outputs read by the Overview; a data file rewritten: stale; following
//      off.
// PASS/FAIL per check; exit 1 on a failure.
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <functional>
#include <iostream>
#include <thread>

#include <QAction>
#include <QDialog>
#include <QPlainTextEdit>
#include <QRegularExpression>
#include <QStatusBar>
#include <QToolBar>
#include <QApplication>
#include <QDir>
#include <QElapsedTimer>
#include <QLabel>
#include <QPushButton>
#include <QScrollArea>
#include <QScrollBar>
#include <QSettings>
#include <QSpinBox>
#include <QTabWidget>
#include <QTableWidget>

#include "AZURESetup.h"
#include "AzureTheme.h"
#include "FittingTab.h"
#include "OverviewCharts.h"
#include "OverviewPage.h"
#include "PlotTab.h"
#include "ProjectHolder.h"
#include "ProjectWatcher.h"
#include "RecordTabs.h"
#include "RecordsPanel.h"
#include "RunTab.h"

namespace fs = std::filesystem;
namespace front = azure2::front;

namespace {

/// std::system for a command line that may start with a quoted program and name /dev/null: cmd.exe drops the first
/// and the last quote of a line that starts with one (an outer pair keeps the inner ones) and calls the sink NUL.
int Shell(std::string command) {
#ifdef _WIN32
  for (std::size_t at = command.find("/dev/null"); at != std::string::npos; at = command.find("/dev/null"))
    command.replace(at, 9, "NUL");
  if (!command.empty() && command[0] == '"') command = "\"" + command + "\"";
#endif
  return std::system(command.c_str());
}

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
  if (!ok) failures++;
}

fs::path CopyProject(const fs::path &projectFile, const fs::path &to) {
  fs::remove_all(to);
  fs::create_directories(to);
  fs::copy(projectFile.parent_path(), to, fs::copy_options::recursive | fs::copy_options::skip_symlinks);
  for (const char *d : {"output", "checks"}) {
    fs::remove_all(to / d);
    fs::create_directories(to / d);
  }
  return to / projectFile.filename();
}

void Settle(int rounds = 5) {
  for (int i = 0; i < rounds; i++) QApplication::processEvents();
}

/// Turn the event loop until `done` or `ms` have passed.
bool WaitFor(const std::function<bool()> &done, int ms) {
  QElapsedTimer timer;
  timer.start();
  while (!done() && timer.elapsed() < ms) QApplication::processEvents(QEventLoop::AllEvents, 50);
  return done();
}

RunOutcome RunAndWait(AZURESetup &w, const QString &command) {
  RunOutcome outcome;
  bool ended = false;
  QMetaObject::Connection c = QObject::connect(&w, &AZURESetup::runFinished, [&](const RunOutcome &o) {
    outcome = o;
    ended = true;
  });
  QString why;
  if (!w.runCommand(command, why)) {
    QObject::disconnect(c);
    std::cout << "  not run: " << why.toStdString() << std::endl;
    outcome.exitCode = -100;
    return outcome;
  }
  WaitFor([&ended]() { return ended; }, 600000);
  QObject::disconnect(c);
  if (!ended) outcome.exitCode = -101;
  return outcome;
}

QString LabelText(QWidget *w, const char *name) {
  QLabel *l = w->findChild<QLabel *>(name);
  return l ? l->text() : QString();
}

/// The cards against the facade's summary of the open project.
void CardsAgainstFacade(AZURESetup &w, OverviewPage *page, const std::string &name) {
  const front::ProjectSummary s = front::SummarizeProject(w.holder()->project());
  const front::ProjectSummary &shown = page->projectSummary();
  Check(w.tabs()->widget(0) == page && w.tabs()->currentWidget() == page, name + ": the Overview is the first tab, shown on opening");
  Check(shown.pairs.size() == s.pairs.size() && shown.levels == s.levels && shown.dataSets == s.dataSets &&
            shown.activeDataSets == s.activeDataSets && shown.basis == s.basis,
        name + ": the page's summary is the facade's (" + std::to_string(s.pairs.size()) + " pairs, " +
            std::to_string(s.levels) + " levels, " + std::to_string(s.dataSets) + " data sets, " + s.basis + ")");
  BarChart *levels = page->findChild<BarChart *>("overviewLevelBars");
  const QString head = LabelText(page, "overviewProjectHead");
  bool pairsShown = true;
  for (const front::PairChoice &p : s.pairs) pairsShown = pairsShown && head.contains(QString::fromStdString(p.label).toHtmlEscaped());
  Check(levels && levels->bars() == (int)s.levelGroups.size() && pairsShown &&
            LabelText(page, "overviewLevelsTitle").contains(QString::number(s.levels)),
        name + ": the pairs, " + std::to_string(s.levelGroups.size()) + " J^pi bars and the level count shown");
  const QString text = LabelText(page, "overviewProjectText");
  Check(text.contains(QString("<b>%1</b> of %2 active").arg(s.activeDataSets).arg(s.dataSets)) &&
            text.contains(QString::fromStdString(s.basis)),
        name + ": the data sets and the basis in the project card");
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: gui_overview_test <source dir> <AZURE2 binary> <screenshot dir>" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, true);
  const fs::path source = argv[1];
  const QString shots = QString::fromLocal8Bit(argv[3]);
  const fs::path work = fs::absolute("gui_overview_work");
  fs::remove_all(work);
  fs::create_directories(work);
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, QString::fromStdString((work / "settings").string()));
  QCoreApplication::setOrganizationName("jina-overview-test");
  QCoreApplication::setApplicationName("azure2-overview-test");
  QDir().mkpath(shots);
  qputenv("AZURE2_THEME", "");
  applyAzureTheme(app);
  setAzureTheme(ThemeChoice::Light);

  AZURESetup w;
  w.resize(1400, 900);
  w.show();
  Settle();
  OverviewPage *page = w.findChild<OverviewPage *>();
  Check(page != nullptr, "the window has an Overview page");
  if (!page) return 1;

  // (a) T1 without outputs
  const fs::path t1 = CopyProject(source / "tests/upgrade/T1/T1.azr", work / "T1");
  QString why;
  Check(w.readProjectFile(QString::fromStdString(t1.string()), why), "T1 opens (" + why.toStdString() + ")");
  Settle();
  CardsAgainstFacade(w, page, "T1");
  bool noRunStep = false;
  for (const OverviewPage::Step &s : page->steps()) noRunStep = noRunStep || s.text.contains("no run yet");
  Check(page->noRunShown() && page->findChild<QPushButton *>("overviewNoRun_calculate") &&
            page->findChild<QPushButton *>("overviewNoRun_fit") && noRunStep,
        "T1 without outputs: 'No run yet' with Calculate and Fit, and the step says so");
  Check(WaitFor([page]() { return !page->parametersBusy() && page->haveParameters(); }, 120000),
        "T1: the parameters' worker ends");
  std::vector<front::Parameter> all;
  std::string unused;
  front::ReadParameters(t1.string(), all, unused);
  int free = 0;
  for (const front::Parameter &p : all) free += !p.fixed;
  Check(page->runParameters().free == free && free > 0 &&
            LabelText(page, "overviewProjectText").contains(QString("<b>%1</b> of %2").arg(free).arg(all.size())),
        "T1: " + std::to_string(free) + " free parameters of " + std::to_string(all.size()) + " in the project card");
  w.grab().save(shots + "/overview_norun_light.png");

  // (b) a calculation through the window
  const RunOutcome outcome = RunAndWait(w, "calculate");
  Check(outcome.exitCode == 0 && outcome.haveChiSquared, "T1: calculated through the window");
  w.tabs()->setCurrentWidget(page);
  Settle();
  Check(WaitFor([page]() { return !page->parametersBusy(); }, 120000), "T1: the parameters' worker ends after the run");
  Settle();
  const front::RunSummary &run = page->runSummary();
  QTableWidget *table = page->dataSetTable();
  bool rows = !page->noRunShown() && run.command == "calculate" && table->rowCount() == outcome.dataSetChiSquared.size();
  double sum = 0;
  for (int i = 0; rows && i < outcome.dataSetChiSquared.size(); i++) {
    const front::DataSetResult &d = run.dataSets[i];
    rows = d.key == outcome.dataSetChiSquared[i].first &&
           std::fabs(d.chiSquared - outcome.dataSetChiSquared[i].second) <= 1e-5 * std::max(1.0, d.chiSquared);
    sum += d.chiSquared;
  }
  Check(rows, "T1: " + std::to_string(table->rowCount()) + " rows, each the run's chi^2 of its data set");
  bool named = true;
  for (int i = 0; i < table->rowCount(); i++) named = named && !table->item(i, 0)->text().startsWith("data set");
  Check(named, "T1: every row named by its record");
  Check(page->tileText("chi2") == QString::number(run.totalChiSquared, 'g', 8) && std::fabs(run.totalChiSquared - sum) < 1e-5 * sum &&
            page->tileText("points") == QString::number(run.points),
        "T1: the chi^2 tile " + page->tileText("chi2").toStdString() + " = the sum of the rows, N " +
            page->tileText("points").toStdString());
  const double dof = run.totalChiSquared / (run.points - free);
  Check(page->tileText("free") == QString::number(free) && page->tileText("dof") == QString::number(dof, 'g', 4) &&
            !page->staleShown(),
        "T1: free parameters " + page->tileText("free").toStdString() + ", chi^2/(N - p) " +
            page->tileText("dof").toStdString() + ", not stale");
  w.grab().save(shots + "/overview_T1_light.png");

  // (c) a row opens the Plot tab with its data set; a step's link goes to its tab
  table->sortItems(4, Qt::DescendingOrder);
  const QString first = table->item(0, 0)->data(Qt::UserRole + 1).toString();
  emit table->cellClicked(0, 0);
  Settle();
  Check(w.tabs()->currentWidget() == w.getPlotTab() && w.getPlotTab()->availableRecords().contains(first),
        "a click on the worst row (" + first.toStdString() + ") opens the Plot tab");
  page->followLink("tab:fitting");
  Check(w.tabs()->currentWidget() == (QWidget *)w.getFittingTab(), "a step's link to Fitting Settings");
  w.tabs()->setCurrentWidget(page);
  Settle();

  // (d) the stale mark
  std::this_thread::sleep_for(std::chrono::milliseconds(30));
  front::Record r = w.holder()->records("data")[0];
  for (front::Field &f : r.fields)
    if (f.key == "norm") f.value = "1.01 fixed";
  Check(w.holder()->setRecord("data", r, why), "T1: a norm edited");
  WaitFor([]() { return false; }, 400);  // the page's quiet time
  Check(!page->staleShown(), "an edit not saved: not stale");
  w.action("save")->trigger();
  Settle();
  bool staleStep = false;
  for (const OverviewPage::Step &s : page->steps()) staleStep = staleStep || s.text.contains("older than the project");
  Check(page->staleShown() && staleStep, "saved with the Save action: the stale mark and its step appear");
  w.grab().save(shots + "/overview_T1_stale_light.png");

  // (f) changes on disk (section 16; added by the director's message of 2026-11-20)
  {
    ProjectWatcher *watcher = w.findChild<ProjectWatcher *>("projectWatcher");
    QToolBar *bar = w.findChild<QToolBar *>("diskChangeBar");
    Check(watcher && bar && watcher->enabled() && w.findChild<QAction *>("followDisk") &&
              w.findChild<QAction *>("followDisk")->isChecked(),
          "the watcher, its bar, and View > Follow Changes on Disk on by default");
    if (!watcher || !bar) return 1;
    int changes = 0;
    QObject::connect(watcher, &ProjectWatcher::projectChanged, [&changes]() { changes++; });
    const fs::path file = t1;
    const std::string dataName = w.holder()->records("data")[0].name;
    auto normInWindow = [&w]() { return w.holder()->records("data")[0].Value("norm"); };
    auto outsideEdit = [&file](const std::string &norm) {  // another program (pyazr, a script) edits and saves
      front::Project other;
      std::string why;
      front::Project::Open(file.string(), other, why);
      front::Record r = other.Records("data")[0];
      for (front::Field &f : r.fields)
        if (f.key == "norm") f.value = norm;
      other.SetRecord("data", r, why);
      other.Save("", why);
    };
    auto diskNorm = [&file]() {
      front::Project other;
      std::string why;
      front::Project::Open(file.string(), other, why);
      return other.Records("data")[0].Value("norm");
    };
    // the window's own save is no change
    w.action("save")->trigger();
    WaitFor([]() { return false; }, 1500);
    Check(changes == 0 && bar->isHidden() && !w.holder()->project().Modified(), "the window's own Save: no reload, no bar");

    // no unsaved edit: re-read at once, the tab kept, the new value shown without a click
    w.tabs()->setCurrentWidget(w.getSegmentsTab());
    outsideEdit("1.02 fixed");
    const bool reloaded = WaitFor([&]() { return normInWindow() == "1.02 fixed"; }, 10000);
    bool shown = false;
    if (RecordsPanel *panel = w.getSegmentsTab()->records("data"))
      for (QAbstractItemModel *m : panel->findChildren<QAbstractItemModel *>())
        for (int r = 0; r < m->rowCount() && !shown; r++)
          for (int c = 0; c < m->columnCount() && !shown; c++) shown = m->index(r, c).data().toString().startsWith("1.02");
    Check(reloaded && shown && w.tabs()->currentWidget() == (QWidget *)w.getSegmentsTab() &&
              w.statusBar()->currentMessage().contains("reloaded") && bar->isHidden(),
          "changed outside, no unsaved edit: reloaded, the Segments tab kept and showing 1.02 (" +
              w.statusBar()->currentMessage().toStdString() + ")");

    // unsaved edits: the bar, and neither side lost
    front::Record mine = w.holder()->records("data")[0];
    for (front::Field &f : mine.fields)
      if (f.key == "norm") f.value = "1.03 fixed";
    w.holder()->setRecord("data", mine, why);
    outsideEdit("1.04 fixed");
    Check(WaitFor([bar]() { return !bar->isHidden(); }, 10000) && normInWindow() == "1.03 fixed" && diskNorm() == "1.04 fixed",
          "changed outside with unsaved edits: the bar; the window keeps 1.03, the disk 1.04");
    w.tabs()->setCurrentWidget(page);
    Settle();
    w.grab().save(shots + "/overview_disk_bar_light.png");
    w.findChild<QAction *>("diskDifferences")->trigger();
    Settle();
    QDialog *diff = w.findChild<QDialog *>("diskDifferencesDialog");
    QPlainTextEdit *diffText = diff ? diff->findChild<QPlainTextEdit *>() : nullptr;
    const QString d = diffText ? diffText->toPlainText() : QString();
    Check(d.contains(QRegularExpression("\\n- +norm +1\\.04 fixed")) && d.contains(QRegularExpression("\\n\\+ +norm +1\\.03 fixed")),
          "Show differences: the disk's line and the window's line");
    if (diff) {
      diff->grab().save(shots + "/overview_disk_differences_light.png");
      diff->close();
    }
    w.findChild<QAction *>("diskKeep")->trigger();
    WaitFor([]() { return false; }, 1500);
    Check(bar->isHidden() && normInWindow() == "1.03 fixed" && diskNorm() == "1.04 fixed" && w.holder()->project().Modified(),
          "Keep mine: the bar goes, both sides as they were");
    changes = 0;
    w.action("save")->trigger();
    WaitFor([]() { return false; }, 1500);
    Check(changes == 0 && diskNorm() == "1.03 fixed" && bar->isHidden(), "saved: the window's version on disk, no reload");

    // outputs written by an outside run: the Overview reads them without a click
    w.tabs()->setCurrentWidget(page);
    outsideEdit("1.2 fixed");
    WaitFor([&]() { return normInWindow() == "1.2 fixed"; }, 10000);
    const QString before = page->tileText("chi2");
    const std::string line = "\"" + std::string(argv[2]) + "\" calculate --threads=2 \"" + file.string() + "\" > \"" +
                             (work / "outside.txt").string() + "\" 2>&1";
    Check(Shell(line) == 0, "an outside run of the program");
    const bool read = WaitFor([&]() { return page->tileText("chi2") != before && !page->staleShown(); }, 15000);
    Check(read, "the outside run's chi^2 shown without a click: " + before.toStdString() + " -> " +
                    page->tileText("chi2").toStdString());
    // a data file written again: the stale mark names the data set
    const fs::path dataFile = file.parent_path() / w.holder()->records("data")[0].Value("file");
    std::string content;
    {
      std::ifstream in(dataFile, std::ios::binary);
      content.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(30));
    {
      std::ofstream out(dataFile, std::ios::binary);
      out << content;
    }
    const bool dataStale = WaitFor([&]() { return page->staleShown(); }, 10000);
    Check(dataStale && page->findChild<QLabel *>("overviewStale")->text().contains(QString::fromStdString(dataName)),
          "a data file rewritten: the stale mark names " + dataName);
    // following off: an outside edit is not read
    w.findChild<QAction *>("followDisk")->setChecked(false);
    outsideEdit("1.3 fixed");
    WaitFor([]() { return false; }, 1500);
    Check(normInWindow() == "1.2 fixed" && QSettings().value("view/followDisk").toBool() == false,
          "Follow Changes on Disk off: not re-read, the choice kept");
    w.findChild<QAction *>("followDisk")->setChecked(true);
  }

  // a fit through the window (T1_start, at most 80 calls): its status in the engine's words and the convergence chart
  {
    const fs::path start = CopyProject(source / "tests/upgrade/T1/T1_start.azr", work / "T1_start");
    Check(w.readProjectFile(QString::fromStdString(start.string()), why), "T1_start opens");
    QSpinBox *calls = nullptr;
    for (QSpinBox *s : w.getFittingTab()->findChildren<QSpinBox *>())
      if (s->toolTip().contains("--max-calls")) calls = s;
    if (calls) calls->setValue(80);
    const RunOutcome fit = RunAndWait(w, "fit");
    w.tabs()->setCurrentWidget(page);
    WaitFor([page]() { return !page->parametersBusy(); }, 120000);
    Settle();
    ChiSquaredTrace *trace = page->findChild<ChiSquaredTrace *>("overviewTrace");
    std::string fileStatus;
    std::ifstream in(start.parent_path() / "output" / "fitStatus.out");
    for (std::string line; std::getline(in, line);)
      if (line.rfind("status ", 0) == 0) fileStatus = line.substr(line.find_first_not_of(' ', 6));
    Check(calls && fit.exitCode == 0 && page->runSummary().command == "fit" && page->runSummary().status == fileStatus &&
              LabelText(page, "overviewRunStatus").contains(QString::fromStdString(fileStatus)) && trace &&
              trace->points() > 0 && trace->isVisible(),
          "a fit through the window: command fit, status '" + fileStatus + "' shown, " +
              std::to_string(trace ? trace->points() : 0) + " points of the convergence chart");
    w.grab().save(shots + "/overview_T1_fit_light.png");
    if (calls) calls->setValue(0);
  }

  // (a) and (e): T2, calculated, no vertical scrolling at 1400 x 900; light and dark
  const fs::path t2 = CopyProject(source / "tests/upgrade/T2/T2.azr", work / "T2");
  Check(w.readProjectFile(QString::fromStdString(t2.string()), why), "T2 opens (" + why.toStdString() + ")");
  Settle();
  CardsAgainstFacade(w, page, "T2");
  const RunOutcome t2run = RunAndWait(w, "calculate");
  Check(t2run.exitCode == 0, "T2: calculated through the window");
  w.tabs()->setCurrentWidget(page);
  WaitFor([page]() { return !page->parametersBusy(); }, 300000);
  Settle(10);
  Check(page->dataSetTable()->rowCount() == t2run.dataSetChiSquared.size(), "T2: a row per data set of the run");
  QScrollArea *scroll = page->findChild<QScrollArea *>();
  Check(scroll && scroll->verticalScrollBar()->maximum() == 0 && scroll->horizontalScrollBar()->maximum() == 0,
        "T2 at 1400 x 900: no scrolling (vertical range " +
            std::to_string(scroll ? scroll->verticalScrollBar()->maximum() : -1) + ")");
  w.grab().save(shots + "/overview_T2_light.png");
  setAzureTheme(ThemeChoice::Dark);
  Settle(10);
  w.grab().save(shots + "/overview_T2_dark.png");
  w.resize(1000, 800);
  Settle(10);
  w.grab().save(shots + "/overview_T2_narrow_dark.png");
  setAzureTheme(ThemeChoice::Light);
  w.resize(1400, 900);
  Settle();

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
