// The window's shell without a display (docs/upgrade/qt_gui.md section 11):
//   gui_shell_test <source dir> <AZURE2 binary> <screenshot dir>
//  (a) the look: the light and the dark palette, the View menu's choice kept in the settings;
//  (b) the start page while no project is open (New, Open, the example projects found), the tool bar's actions and
//      their state, the status bar (file, modified mark, validity), the Diagnostics dock hidden on a valid project and
//      opened by an error, a double click on a diagnostic going to its tab, the recent projects, a dropped .azr;
//  (c) the Calculate tab's editors made from the facade's option rows (units, ranges, an option enabled only with the
//      flag it needs), a value set in an editor reaching the command line; calculate from the tool bar: the summary
//      (a row per data set with its name, the total = the engine's), sortable, Open Plot; a fit's chi^2 curve;
//  (d) the Plot tab: check-box lists with All and None, wheel zoom, pan, the cursor's read-out, CSV export;
//  (e) screenshots of the start page and the tabs in the light and the dark look.
// PASS/FAIL per check; exit 1 on a failure.
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>

#include <QAction>
#include <QApplication>
#include <QFileInfo>
#include <QCheckBox>
#include <QDir>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QDropEvent>
#include <QElapsedTimer>
#include <QFile>
#include <QLabel>
#include <QListWidget>
#include <QMimeData>
#include <QMouseEvent>
#include <QProgressBar>
#include <QPushButton>
#include <QSettings>
#include <QSpinBox>
#include <QStackedWidget>
#include <QStatusBar>
#include <QTabWidget>
#include <QTableWidget>
#include <QToolBar>
#include <QToolButton>

#include "AZURESetup.h"
#include "AzureTheme.h"
#include "DiagnosticsPanel.h"
#include "FittingTab.h"
#include "MCMCTab.h"
#include "PlotTab.h"
#include "PlotWidget.h"
#include "RecordTabs.h"
#include "RunTab.h"
#include "StartPage.h"

namespace fs = std::filesystem;

namespace {

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

/// Trigger `start`, then turn the event loop until the window's run ends.
RunOutcome RunAndWait(AZURESetup &w, const std::function<void()> &start) {
  RunOutcome outcome;
  bool ended = false;
  QMetaObject::Connection c = QObject::connect(&w, &AZURESetup::runFinished, [&](const RunOutcome &o) {
    outcome = o;
    ended = true;
  });
  start();
  if (!w.running()) {
    QObject::disconnect(c);
    outcome.exitCode = -100;
    return outcome;
  }
  QElapsedTimer timer;
  timer.start();
  while (!ended && timer.elapsed() < 600000) QApplication::processEvents(QEventLoop::AllEvents, 50);
  QObject::disconnect(c);
  if (!ended) outcome.exitCode = -101;
  return outcome;
}

void Settle() {
  for (int i = 0; i < 5; i++) QApplication::processEvents();
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: gui_shell_test <source dir> <AZURE2 binary> <screenshot dir>" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, true);
  const fs::path source = argv[1];
  const QString shots = QString::fromLocal8Bit(argv[3]);
  const fs::path work = fs::absolute("gui_shell_work");
  fs::remove_all(work);
  fs::create_directories(work);
  // the settings of this test only (recent files, window state, the look), in its work directory
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, QString::fromStdString((work / "settings").string()));
  QCoreApplication::setOrganizationName("jina-shell-test");
  QCoreApplication::setApplicationName("azure2-shell-test");
  QDir().mkpath(shots);

  // (a) the look
  qputenv("AZURE2_THEME", "");
  applyAzureTheme(app);
  setAzureTheme(ThemeChoice::Dark);
  const bool dark = azureThemeIsDark() && app.palette().color(QPalette::Window).lightness() < 100 &&
                    QSettings().value("view/theme").toString() == "dark";
  setAzureTheme(ThemeChoice::Light);
  const bool light = !azureThemeIsDark() && app.palette().color(QPalette::Window).lightness() > 200;
  Check(dark && light, "the look: the dark palette (window lightness < 100) and the light one (> 200), kept in the settings");
  Check(azurePalette(true).color(QPalette::Text).lightness() > 180 && azurePalette(false).color(QPalette::Text).lightness() < 80,
        "the look: readable text in both palettes");
  Check(app.styleSheet().contains("QHeaderView::section") && app.styleSheet().contains("primary=\"true\"") &&
            !app.styleSheet().contains("font-size"),
        "the look: one style sheet (headers, primary buttons), no font size in it");

  AZURESetup w;
  w.resize(1400, 900);
  w.show();
  Settle();

  // (b) the start page, the tool bar, the status bar
  QStackedWidget *central = qobject_cast<QStackedWidget *>(w.centralWidget());
  Check(!w.hasProject() && central && central->currentWidget() == w.startPage(),
        "the start page is shown while no project is open");
  QStringList toolNames;
  for (QAction *a : w.mainToolBar()->actions())
    if (!a->isSeparator()) toolNames << a->objectName().mid(7);
  Check(toolNames == QStringList({"new", "open", "save", "validate", "calculate", "fit", "stop", "plot"}),
        "the tool bar: " + toolNames.join(" ").toStdString());
  bool icons = true;
  for (const QString &n : toolNames) {
    const bool ok = !w.action(n)->icon().isNull() && !w.action(n)->icon().pixmap(QSize(22, 22)).isNull() &&
                    !w.action(n)->toolTip().isEmpty();
    if (!ok) std::cout << "  action " << n.toStdString() << ": icon null " << w.action(n)->icon().isNull() << ", tip '"
                       << w.action(n)->toolTip().toStdString() << "'" << std::endl;
    icons = icons && ok;
  }
  Check(icons && w.action("calculate")->shortcut() == QKeySequence(Qt::Key_F5) &&
            w.action("fit")->shortcut() == QKeySequence(Qt::Key_F6) && w.action("validate")->shortcut() == QKeySequence(Qt::Key_F7),
        "the tool bar's actions have icons, tool tips and shortcuts (F5 calculate, F6 fit, F7 validate)");
  Check(w.action("new")->isEnabled() && w.action("open")->isEnabled() && !w.action("save")->isEnabled() &&
            !w.action("calculate")->isEnabled() && !w.action("stop")->isEnabled(),
        "without a project: New and Open enabled; Save, Calculate and Stop not");
  Check(w.statusFile()->text().contains("No project"), "the status bar says no project is open");
  const QStringList examples = StartPage::findExamples();
  bool t1Listed = false;
  for (const QString &e : examples) t1Listed = t1Listed || e.endsWith("tests/upgrade/T1/T1.azr");
  // the list's project lines (its group headings, qt_gui.md section 19.9, carry no file)
  int projectLines = 0;
  for (int i = 0; i < w.startPage()->exampleList()->count(); i++)
    projectLines += !w.startPage()->exampleList()->item(i)->data(Qt::UserRole).toString().isEmpty();
  Check(examples.size() >= 4 && t1Listed && projectLines == examples.size(),
        "the start page lists " + std::to_string(examples.size()) + " example projects, T1 among them");
  QImage startLight = w.grab().toImage();
  startLight.save(shots + "/shell_start_light.png");

  // a project opened: the tabs, the status bar, the actions
  const fs::path t1 = CopyProject(source / "tests/upgrade/T1/T1.azr", work / "T1");
  QString why;
  Check(w.readProjectFile(QString::fromStdString(t1.string()), why), "T1 opened: " + why.toStdString());
  Settle();
  Check(w.hasProject() && central->currentWidget() == w.tabs() && w.statusFile()->text().contains("T1.azr") &&
            !w.statusFile()->text().contains("modified") && w.action("save")->isEnabled() &&
            w.action("calculate")->isEnabled(),
        "T1 open: the tabs, the file in the status bar, Save and Calculate enabled");
  const bool t1Clean = w.getDiagnostics()->errorCount() == 0 && w.getDiagnostics()->warningCount() == 0;
  Check(t1Clean && !w.diagnosticsDockWidget()->isVisible() && w.statusValidity()->isVisible() &&
            w.statusValidity()->text().startsWith("Valid"),
        "T1 valid: the Diagnostics dock out of the way, \"" + w.statusValidity()->text().toStdString() + "\" in the status bar");
  // a broken unit: the modified mark, the dock opens, its double click goes to the record's tab
  std::vector<azure2::front::Record> data = w.holder()->records("data");
  azure2::front::Record broken = data.at(0);
  for (azure2::front::Field &f : broken.fields)
    if (f.key == "unit") f.value = "barn_per_furlong";
  Check(w.holder()->setRecord("data", broken, why), "T1: a data set's unit broken through the project");
  Settle();
  Check(w.isWindowModified() && w.statusFile()->text().contains("modified"), "the modified mark in the title and the status bar");
  w.action("validate")->trigger();
  Settle();
  Check(w.getDiagnostics()->errorCount() > 0 && w.diagnosticsDockWidget()->isVisible() &&
            w.statusValidity()->text().contains("error"),
        "an error: the dock opens, \"" + w.statusValidity()->text().toStdString() + "\" in the status bar");
  w.grab().save(shots + "/shell_diagnostics_light.png");
  w.tabs()->setCurrentIndex(0);
  int errorRow = -1;
  for (int i = 0; i < w.getDiagnostics()->rowCount() && errorRow < 0; i++) errorRow = i;
  w.getDiagnostics()->activateRow(errorRow);
  Settle();
  Check(w.tabs()->currentWidget() == w.getSegmentsTab(), "a double click on the diagnostic shows the Segments tab");
  w.statusValidity()->click();
  Check(w.diagnosticsDockWidget()->isVisible(), "the validity in the status bar opens the Diagnostics");
  // back to the valid record, saved: the dock closes again; the project closed: the start page and its recent list
  Check(w.holder()->setRecord("data", data.at(0), why) && w.writeProjectFile(QString::fromStdString(t1.string()), why),
        "T1 restored and saved");
  w.action("validate")->trigger();
  Settle();
  Check(!w.diagnosticsDockWidget()->isVisible() && !w.isWindowModified(), "valid again: the dock hides, no modified mark");
  Check(w.closeProject() && central->currentWidget() == w.startPage() && !w.action("save")->isEnabled(),
        "File > Close: the start page again");
  bool recent = false;
  for (int i = 0; i < w.startPage()->recentList()->count(); i++)
    recent = recent || w.startPage()->recentList()->item(i)->data(Qt::UserRole).toString() == QString::fromStdString(t1.string());
  Check(recent, "the start page lists T1 among the recent projects");
  // a dropped .azr is opened
  Check(AZURESetup::acceptsDrop({QUrl::fromLocalFile("/x/a.azr")}) && !AZURESetup::acceptsDrop({QUrl::fromLocalFile("/x/a.txt")}),
        "a dropped .azr is taken, a .txt is not");
  {
    Check(w.acceptDrops(), "the window takes drops");
    Check(w.openDropped({QUrl::fromLocalFile("/x/notes.txt"), QUrl::fromLocalFile(QString::fromStdString(t1.string()))}),
          "a drop of a .txt and T1.azr: T1 is taken");
    Settle();
    Check(w.hasProject() && QFileInfo(w.holder()->path()) == QFileInfo(QString::fromStdString(t1.string())), "T1 dropped on the window: opened");
  }

  // (c) the Calculate tab's editors from the facade's rows
  RunTab *run = w.getRunTab();
  run->setCommand("rate");
  QDoubleSpinBox *tmin = qobject_cast<QDoubleSpinBox *>(run->editor("--temperature-min="));
  Check(tmin && tmin->suffix() == " GK" && tmin->value() == tmin->minimum() && !run->options().join(" ").contains("temperature"),
        "rate: --temperature-min a spin box in GK, not set until given");
  Check(run->setOption("--temperature-min=", "0.01") && run->options().contains("--temperature-min=0.01"),
        "rate: 0.01 GK in the editor gives --temperature-min=0.01 (" + run->options().join(" ").toStdString() + ")");
  run->setOption("--temperature-min=", "");
  QLabel *hint = run->findChild<QLabel *>("hint_temperature-min");
  Check(hint && hint->text().contains("temperature") && hint->text().contains("--temperature-file"),
        "rate: the option's help in one line beside it, with its other way");
  run->setCommand("fit");
  QCheckBox *band = qobject_cast<QCheckBox *>(run->editor("--covariance-band"));
  QWidget *scale = run->editor("--scale-covariance");
  const bool scaleOff = scale && !scale->isEnabled();
  if (band) band->setChecked(true);
  Check(band && scaleOff && scale->isEnabled(), "fit: --scale-covariance enabled only with --covariance-band");
  band->setChecked(false);
  run->setCommand("minos");
  Check(run->setOption("--delta-chi2=", "0") && !run->options().join(" ").contains("delta") &&
            run->setOption("--delta-chi2=", "2.71") && run->options().contains("--delta-chi2=2.71"),
        "minos: --delta-chi2 above 0 (0 is not set), 2.71 reaches the command line");
  run->setOption("--delta-chi2=", "");
  Check(RunTab::lineSeverity("Warning: the band was not written") == 1 && RunTab::lineSeverity("ERROR: refused") == 2 &&
            RunTab::lineSeverity("Segment #1 Chi-Squared/N: 17.4") == 0,
        "the log: warnings and errors highlighted, other lines not");

  // calculate from the tool bar: the summary
  if (QSpinBox *threads = run->findChild<QSpinBox *>("threads")) threads->setValue(2);
  const RunOutcome calc = RunAndWait(w, [&]() { w.action("calculate")->trigger(); });
  Settle();
  QTableWidget *chi = run->chiTable();
  const QStringList names = w.getPlotTab()->dataRecords();
  bool named = chi->rowCount() == calc.dataSetChiSquared.size() && chi->rowCount() == names.size();
  for (int i = 0; named && i < chi->rowCount(); i++) named = chi->item(i, 1)->text() == names[i];
  Check(calc.exitCode == 0 && named && chi->rowCount() >= 7,
        "calculate from the tool bar: a row per data set with its record's name (" + std::to_string(chi->rowCount()) + ")");
  Check(calc.haveChiSquared && run->totalLabel()->text().contains(QString::number(calc.totalChiSquared, 'g', 12)),
        "the total chi-squared of the engine in the summary: " + run->totalLabel()->text().toStdString());
  chi->sortByColumn(3, Qt::DescendingOrder);
  bool sorted = true;
  for (int i = 1; i < chi->rowCount(); i++)
    sorted = sorted && chi->item(i - 1, 3)->data(Qt::UserRole).toDouble() >= chi->item(i, 3)->data(Qt::UserRole).toDouble();
  Check(sorted, "the summary sorts by chi-squared");
  chi->sortByColumn(0, Qt::AscendingOrder);
  Check(!w.statusJob()->isVisible() && w.statusBar()->currentMessage().contains("total chi-squared"),
        "the status bar after the run: \"" + w.statusBar()->currentMessage().toStdString() + "\"");
  w.tabs()->setCurrentWidget(run);
  Settle();
  w.grab().save(shots + "/shell_calculate_light.png");
  run->plotButton()->click();
  Check(w.tabs()->currentWidget() == w.getPlotTab(), "Open Plot shows the Plot tab");

  // (d) the Plot tab
  PlotTab *plot = w.getPlotTab();
  PlotWidget *pw = plot->plot();
  plot->checkAll(false, true);
  Check(plot->checkedRecords().size() == names.size() && pw->series().size() == names.size(),
        "Plot: All checks every data set (" + std::to_string(pw->series().size()) + " series)");
  plot->checkAll(false, false);
  Check(pw->series().isEmpty() && plot->checkedRecords().isEmpty(), "Plot: None clears the plot");
  plot->selectRecords({"cap_g0", "el_lab090"});
  plot->setCurveStyle("cap_g0", QColor("#123456"), 2);
  bool styled = false;
  for (const PlotSeries &ps : pw->series()) styled = styled || (ps.label == "cap_g0" && ps.color == QColor("#123456") && ps.marker == 2);
  plot->setCurveStyle("cap_g0", QColor(), -1);
  bool restored = false;
  for (const PlotSeries &ps : pw->series()) restored = restored || (ps.label == "cap_g0" && ps.color != QColor("#123456"));
  Check(styled && restored, "Plot: a curve's own colour and marker, and back to the palette's");
  Settle();
  const QRectF whole = pw->view();
  const QPointF centre = pw->rect().center();
  pw->wheelZoom(centre, 1);
  const QRectF in = pw->view();
  Check(in.width() < whole.width() * 0.85 && in.width() > whole.width() * 0.75 && whole.contains(in.center()),
        "Plot: a wheel step zooms in about the cursor (width " + std::to_string(in.width() / whole.width()) + " of the whole)");
  pw->wheelZoom(centre, 1, Qt::ControlModifier);
  Check(pw->view().height() == in.height() && pw->view().width() < in.width(), "Plot: Ctrl + wheel zooms x only");
  const QRectF before = pw->view();
  pw->panBy(QPoint(40, 0));
  Check(pw->view().left() < before.left() && std::fabs(pw->view().width() - before.width()) < 1e-9 * before.width(),
        "Plot: a drag to the right pans the view to lower x");
  QMouseEvent move(QEvent::MouseMove, centre, pw->mapToGlobal(centre), Qt::NoButton, Qt::NoButton, Qt::NoModifier);
  QApplication::sendEvent(pw, &move);
  QPointF at;
  Check(pw->dataAt(centre, at) && pw->cursorText().startsWith("x = ") && plot->cursorLabel()->text().startsWith("x = "),
        "Plot: the cursor's coordinates read out (" + pw->cursorText().toStdString() + ")");
  pw->resetZoom();
  const QString csvFile = QString::fromStdString((work / "plot.csv").string());
  Check(plot->exportPlot(csvFile, why), "Plot: CSV export " + why.toStdString());
  {
    QFile f(csvFile);
    Check(f.open(QIODevice::ReadOnly | QIODevice::Text), "Plot: the CSV file reads back");
    int dataRows = 0;
    bool inData = false;
    for (const QByteArray &line : f.readAll().split('\n')) {
      if (line.startsWith("# ")) inData = line.endsWith(": data");
      else if (inData && !line.isEmpty() && !line.startsWith("x,")) dataRows++;
    }
    int points = 0;
    for (const PlotSeries &s : pw->series()) points += s.dataX.size();
    Check(points > 0 && dataRows == points, "Plot: the CSV holds the " + std::to_string(points) + " data points drawn (" +
                                                std::to_string(dataRows) + " rows)");
  }
  plot->checkAll(false, true);
  Settle();
  w.grab().save(shots + "/shell_plot_light.png");

  // a fit with a call limit: the chi^2 curve; the job in the status bar while it runs
  {
    const fs::path start = CopyProject(source / "tests/upgrade/T1/T1_start.azr", work / "fit_T1");
    w.readProjectFile(QString::fromStdString(start.string()), why);
    w.getFittingTab()->findChild<QSpinBox *>()->setValue(40);
    bool jobShown = false;
    RunOutcome fit;
    {
      bool ended = false;
      QMetaObject::Connection c = QObject::connect(&w, &AZURESetup::runFinished, [&](const RunOutcome &o) {
        fit = o;
        ended = true;
      });
      w.action("fit")->trigger();
      QElapsedTimer timer;
      timer.start();
      while (!ended && timer.elapsed() < 600000) {
        QApplication::processEvents(QEventLoop::AllEvents, 50);
        jobShown = jobShown || (w.statusJob()->isVisible() && w.action("stop")->isEnabled() && !w.action("fit")->isEnabled());
      }
      QObject::disconnect(c);
    }
    Check(fit.exitCode == 0 && jobShown, "fit from the tool bar: the job and Stop in the status bar while it ran");
    Check(run->trace()->points() >= 2, "fit: the best chi-squared drawn as a curve (" + std::to_string(run->trace()->points()) +
                                           " points)");
    w.tabs()->setCurrentWidget(run);
    Settle();
    w.grab().save(shots + "/shell_fit_light.png");
  }

  // the sampler's results: a chain written here (3 steps, 2 walkers, 2 parameters; two of the four moves stay put)
  if (MCMCTab *mcmc = w.getMCMCTab()) {
    {
      std::ofstream f(fs::path(w.holder()->directory().toStdString()) / "output" / "samples.mcmc");
      f << "step,walker,logprob,loglikelihood,logprior,param0,param1\n"
        << "0,0,-5,-5,0,1,10\n0,1,-6,-6,0,3,12\n1,0,-5,-5,0,1,10\n1,1,-4,-4,0,5,14\n2,0,-3,-3,0,2,11\n2,1,-4,-4,0,5,14\n";
    }
    mcmc->setBurnIn(0);
    const bool read = mcmc->readResults(why);
    QTableWidget *res = mcmc->resultsTable();
    const double mean0 = res->rowCount() == 2 ? res->item(0, 1)->data(Qt::UserRole).toDouble() : NAN;
    const double median1 = res->rowCount() == 2 ? res->item(1, 4)->data(Qt::UserRole).toDouble() : NAN;
    Check(read && res->rowCount() == 2 && std::fabs(mean0 - 17.0 / 6) < 1e-12 && median1 == 11.5 &&
              std::fabs(mcmc->acceptanceFraction() - 0.5) < 1e-12,
          "sampler results: mean 17/6, median 11.5, acceptance 0.5 from a written chain (" + why.toStdString() + ")");
    Check(mcmc->showResultPlot("trace:1") && mcmc->resultsPlot()->series().size() == 2 &&
              mcmc->showResultPlot("histogram:0") && mcmc->showResultPlot("logprob") &&
              mcmc->resultsPlot()->series().at(0).modelY.at(2) == -3.5,
          "sampler results: trace (a curve per walker), histogram, the walkers' mean log posterior");
    mcmc->setBurnIn(1);
    Check(std::fabs(mcmc->resultsTable()->item(0, 1)->data(Qt::UserRole).toDouble() - 13.0 / 4) < 1e-12,
          "sampler results: a burn-in of one step leaves the first out (mean 13/4)");
    w.tabs()->setCurrentWidget(mcmc);
    if (QTabWidget *pages = mcmc->findChild<QTabWidget *>("mcmcPages")) pages->setCurrentIndex(2);
    mcmc->showResultPlot("trace:0");
    Settle();
    w.grab().save(shots + "/shell_mcmc_results_light.png");
  }

  // (e) the dark look
  setAzureTheme(ThemeChoice::Dark);
  Settle();
  for (const auto &tab : {std::make_pair(static_cast<QWidget *>(run), QString("calculate")),
                          std::make_pair(static_cast<QWidget *>(plot), QString("plot")),
                          std::make_pair(static_cast<QWidget *>(w.getSegmentsTab()), QString("segments")),
                          std::make_pair(static_cast<QWidget *>(w.getFittingTab()), QString("fitting"))}) {
    w.tabs()->setCurrentWidget(tab.first);
    Settle();
    w.grab().save(shots + "/shell_" + tab.second + "_dark.png");
  }
  const QImage plotDark = pw->grab().toImage();
  Check(QColor(plotDark.pixel(5, 5)).lightness() < 100, "the plot drawn in the dark palette on screen");
  w.closeProject();
  Settle();
  w.grab().save(shots + "/shell_start_dark.png");
  setAzureTheme(ThemeChoice::Light);
  Settle();
  w.tabs()->setCurrentWidget(w.getSegmentsTab());
  w.readProjectFile(QString::fromStdString(t1.string()), why);
  Settle();
  w.grab().save(shots + "/shell_segments_light.png");
  w.tabs()->setCurrentWidget(w.getFittingTab());
  Settle();
  w.grab().save(shots + "/shell_fitting_light.png");

  std::cout << (failures ? "FAILED: " : "OK: ") << failures << " failure(s)" << std::endl;
  return failures ? 1 : 0;
}
