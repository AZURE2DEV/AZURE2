// The Qt setup program without a display (docs/upgrade/qt_gui.md section 5):
//   gui_window_test <source dir> <AZURE2 binary> <screenshot dir>
//  (a) T1, T2, T3 and 13N opened in the window and saved byte for byte; every record of every section written back
//      through its form (Apply) without moving a byte;
//  (b) a channel's width edited in the Levels tab: one row of <levels> changes, the engine reads the new width;
//      a field broken in a form: the engine's code in the diagnostics panel; an old project refused on opening;
//  (c) calculate run from the Calculate tab in the worker thread: the program's output files byte for byte, the
//      window responsive meanwhile (its event loop turns); a fit with a call limit; a fit stopped from the window;
//  (d) every tab grabbed into the screenshot directory.
// PASS/FAIL per check; exit 1 on a failure.
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include <QApplication>
#include <QComboBox>
#include <QDir>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>
#include <QTableWidget>
#include <QTimer>

#include "AZURESetup.h"
#include "DiagnosticsPanel.h"
#include "FittingTab.h"
#include "LevelsTab.h"
#include "RecordTabs.h"
#include "RecordsPanel.h"
#include "RunTab.h"

namespace fs = std::filesystem;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
  if (!ok) failures++;
}

std::string Bytes(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
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

/// Run the window's run of its Calculate tab and wait for its end, turning the event loop; how many times the loop
/// turned while the run was in progress says the window stayed responsive.
RunOutcome RunAndWait(AZURESetup &w, int &turns, std::function<void()> during = nullptr) {
  RunOutcome outcome;
  bool ended = false;
  QMetaObject::Connection c = QObject::connect(&w, &AZURESetup::runFinished, [&](const RunOutcome &o) {
    outcome = o;
    ended = true;
  });
  QString why;
  if (!w.startRun(why)) {
    std::cout << "  the run did not start: " << why.toStdString() << std::endl;
    QObject::disconnect(c);
    outcome.exitCode = -100;
    return outcome;
  }
  turns = 0;
  QElapsedTimer timer;
  timer.start();
  while (!ended && timer.elapsed() < 600000) {
    QApplication::processEvents(QEventLoop::AllEvents, 50);
    turns++;
    if (during) during();
  }
  QObject::disconnect(c);
  if (!ended) outcome.exitCode = -101;
  return outcome;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: gui_window_test <source dir> <AZURE2 binary> <screenshot dir>" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("jina-test");
  QCoreApplication::setApplicationName("azure2-gui-test");
  const fs::path source = argv[1];
  const std::string binary = argv[2];
  const QString shots = QString::fromLocal8Bit(argv[3]);
  const fs::path work = fs::absolute("gui_test_work");
  fs::create_directories(work);

  AZURESetup w;
  w.resize(1400, 900);
  w.show();

  // (a) round trips through the window
  const std::vector<fs::path> projects = {source / "tests/upgrade/T1/T1.azr", source / "tests/upgrade/T2/T2.azr",
                                          source / "tests/upgrade/T3/T3.azr", source / "tests/13N/13N.azr"};
  for (const fs::path &f : projects) {
    const std::string name = f.stem().string();
    const fs::path copy = CopyProject(f, work / ("rt_" + name));
    QString why;
    Check(w.readProjectFile(QString::fromStdString(copy.string()), why), name + ": opened (" + why.toStdString() + ")");
    Check(w.getDiagnostics()->errorCount() == 0, name + ": the engine's diagnostics: no error");
    const fs::path out = work / ("rt_" + name) / (name + "_saved.azr");
    Check(w.writeProjectFile(QString::fromStdString(out.string()), why) && Bytes(out) == Bytes(f),
          name + ": saved byte for byte (" + why.toStdString() + ")");
    // every record through its form
    int applied = 0;
    bool allApplied = true;
    for (SectionsTab *tab : {w.getSegmentsTab(), w.getEffectsTab()})
      for (const char *section : {"data", "extrapolations", "modelData", "targets"}) {
        RecordsPanel *panel = tab->records(section);
        if (!panel) continue;
        for (const azure2::front::Record &r : w.holder()->records(section)) {
          allApplied = allApplied && panel->selectRecord(QString::fromStdString(r.name)) && panel->applyForm(why);
          applied++;
        }
      }
    const std::string text = w.holder()->project().Text(), original = Bytes(f);
    if (text != original) {
      std::size_t at = 0;
      while (at < text.size() && at < original.size() && text[at] == original[at]) at++;
      const std::size_t from = original.rfind('\n', at) + 1;
      std::cout << "  first difference, the file:   " << original.substr(from, original.find('\n', at) - from) << "\n"
                << "                  the window: " << text.substr(from, text.find('\n', at) - from) << std::endl;
    }
    Check(allApplied && text == original,
          name + ": " + std::to_string(applied) + " records written back through their forms, not a byte moved");
  }

  // (b) edits
  {
    const fs::path copy = CopyProject(projects[0], work / "edit_T1");
    QString why;
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    const std::vector<azure2::front::LevelRow> before = w.holder()->project().Levels();
    ChannelsModel *channels = w.getLevelsTab()->getChannelsModel();
    channels->setData(channels->index(0, 6), 1.91, Qt::EditRole);
    Check(w.writeProjectFile(QString::fromStdString(copy.string()), why), "T1: saved after a width was edited");
    azure2::front::Project back;
    std::string reason;
    azure2::front::Project::Open(copy.string(), back, reason);
    const std::vector<azure2::front::LevelRow> after = back.Levels();
    int changed = 0;
    for (std::size_t i = 0; i < before.size() && i < after.size(); i++) changed += before[i].tokens != after[i].tokens;
    bool found = false;
    for (const azure2::front::LevelRow &r : after) found = found || (r.tokens.size() > 11 && r.tokens[11] == "1.91");
    Check(found && changed == 1 && before.size() == after.size(),
          "T1: the edited width read back by the engine (" + std::to_string(changed) + " row(s) differ)");
    Check(back.Validate().ok, "T1: the edited project is valid");
    // a field broken in the form: the engine names it
    RecordsPanel *data = w.getSegmentsTab()->records("data");
    data->selectRecord("el_lab090");
    if (QComboBox *q = qobject_cast<QComboBox *>(data->fieldWidget("unit"))) q->setCurrentText("barn_per_furlong");
    data->applyForm(why);
    w.getDiagnostics()->validate();
    Check(w.getDiagnostics()->errorCount() > 0 && !w.getDiagnostics()->codes().isEmpty(),
          "T1: a wrong unit: the engine's diagnostic " + w.getDiagnostics()->codes().join(",").toStdString());
    Check(!w.readProjectFile(QString::fromStdString((source / "tests/convert/fixtures/tests/13N/13N.azr").string()), why) &&
              why.contains("--convert"),
          "an old project refused on opening, with the conversion");
  }

  // (c) runs
  int turns = 0;
  {
    const fs::path copy = CopyProject(projects[0], work / "run_T1");
    const fs::path cli = CopyProject(projects[0], work / "cli_T1");
    QString why;
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    w.getRunTab()->setCommand("calculate");
    const RunOutcome o = RunAndWait(w, turns);
    Check(o.exitCode == 0 && o.haveChiSquared, "T1: calculate from the window ends with exit 0 and a chi^2");
    Check(turns > 3, "T1: the window's event loop turned " + std::to_string(turns) + " times during the run");
    const std::string command = "\"" + binary + "\" calculate --threads=2 \"" + cli.string() + "\" > /dev/null 2>&1";
    Check(std::system(command.c_str()) == 0, "T1: the program's calculate");
    int same = 0, total = 0;
    for (const auto &e : fs::directory_iterator(cli.parent_path() / "output")) {
      total++;
      same += Bytes(e.path()) == Bytes(copy.parent_path() / "output" / e.path().filename());
    }
    Check(total > 0 && same == total, "T1: the window's output files = the program's: " + std::to_string(same) + " of " +
                                          std::to_string(total) + " identical");
    Check(w.getRunTab()->chiTable()->rowCount() > 1 && w.getRunTab()->log().contains("Total Chi-Squared"),
          "T1: the chi^2 table and the program's lines in the Calculate tab");
  }
  {
    const fs::path copy = CopyProject(source / "tests/upgrade/T1/T1_start.azr", work / "fit_T1");
    QString why;
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    w.getRunTab()->setCommand("fit");
    w.getFittingTab()->findChild<QSpinBox *>()->setValue(40);  // the call limit
    const RunOutcome o = RunAndWait(w, turns);
    Check(o.exitCode == 0 && fs::exists(copy.parent_path() / "output" / "fitStatus.out"),
          "T1_start: a fit with a call limit of 40 from the window, fitStatus.out written");
    w.getFittingTab()->findChild<QSpinBox *>()->setValue(0);
    bool asked = false;
    QLabel *progress = w.getRunTab()->findChild<QLabel *>("runProgress");
    const RunOutcome s = RunAndWait(w, turns, [&]() {
      // stop once the fit has begun (its progress line shows)
      if (!asked && w.running() && progress && progress->text().startsWith("fit")) {
        asked = true;
        w.stopRun();
      }
    });
    const std::string status = Bytes(copy.parent_path() / "output" / "fitStatus.out");
    Check(s.exitCode == 0 && s.stopped && status.find("stopped") != std::string::npos,
          "T1_start: a fit stopped from the window ends at its best point (fitStatus.out: stopped)");
  }

  // (d) screenshots of every tab, T2 (data sets, targets, THM) open
  {
    const fs::path copy = CopyProject(projects[1], work / "shots_T2");
    QString why;
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    w.getSegmentsTab()->records("data")->selectRecord("s02_diff090");
    w.getEffectsTab()->records("targets")->selectRecord("tgt_s05_implanted");
    w.saveTabScreenshots(shots);
    const int n = (int)QDir(shots).entryList({"tab*.png"}).size();
    Check(n == w.tabs()->count(), "screenshots of " + std::to_string(n) + " tabs in " + shots.toStdString());
  }

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
