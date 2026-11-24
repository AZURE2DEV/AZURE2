// The window for an occasional user (docs/upgrade/qt_gui.md section 19), without a display:
//   gui_casual_test <source dir> <screenshot dir>
//  (a) every shipped example (examples/*/ *.azr) and the test projects of the start page (tests/upgrade/T1-T3,
//      tests/13N) open through the window's own open path with no modal dialog at all; the pairs read from <levels>
//      keep their keys (c12c12_tumino2018 names one pair under two keys: each key keeps its pair, and the rows the
//      Levels tab would write back carry the masses and charges of the file);
//  (b) the notes bar: 13N (valid, 50 warnings) opens with one note of the engine's check and the Diagnostics dock
//      closed; Details shows every note in a window that is not modal; Dismiss hides the bar; a valid project opens
//      with no bar; opening another project clears the notes;
//  (d) the converter's repair of an old line: the old 13N with one data line short of its shift's vary flag; the
//      window shows the line, its number and the two readings; the first taken and converted again gives the project
//      the unmodified old 13N gives, byte for byte; the old file is unchanged, the repaired copy beside it, a note says
//      so;
// PASS/FAIL per check; exit 1 on a failure.
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>

#include <QApplication>
#include <QDialog>
#include <QDir>
#include <QElapsedTimer>
#include <QFileInfo>
#include <QMessageBox>
#include <QSettings>
#include <QTimer>

#include "AZURESetup.h"
#include "NotesBar.h"
#include "ConvertRepairDialog.h"
#include <QFile>
#include <QLineEdit>
#include <QPushButton>
#include <fstream>
#include <sstream>
#include <QAction>
#include <QDockWidget>
#include "LevelsTab.h"
#include "PairsModel.h"
#include "PairsTab.h"

namespace fs = std::filesystem;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  std::cout << (ok ? "PASS " : "FAIL ") << what << std::endl;
  if (!ok) failures++;
}

void Settle(int ms = 300) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < ms) QApplication::processEvents(QEventLoop::AllEvents, 20);
}

/// A copy of a project's folder (the window may write nothing into the source tree).
fs::path CopyProject(const fs::path &projectFile, const fs::path &to) {
  fs::remove_all(to);
  fs::create_directories(to);
  fs::copy(projectFile.parent_path(), to, fs::copy_options::recursive | fs::copy_options::skip_symlinks);
  return to / projectFile.filename();
}

/// Counts and closes every modal dialog that appears while it lives (a dialog the window would show a user).
struct ModalWatch {
  QTimer timer;
  QStringList seen;
  ModalWatch() {
    timer.setInterval(20);
    QObject::connect(&timer, &QTimer::timeout, [this]() {
      QWidget *m = QApplication::activeModalWidget();
      if (!m) return;
      QString text = m->windowTitle();
      if (QMessageBox *b = qobject_cast<QMessageBox *>(m)) text += ": " + b->text();
      seen << text;
      if (QDialog *d = qobject_cast<QDialog *>(m)) d->reject();
      else m->close();
    });
    timer.start();
  }
};

double Num(const std::vector<std::string> &t, int c) {
  return c >= 0 && c < (int)t.size() && !t[c].empty() ? std::stod(t[c]) : 0.0;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: gui_casual_test <source dir> <screenshot dir>" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, true);
  const fs::path source = argv[1];
  const QString shots = QString::fromLocal8Bit(argv[2]);
  QDir().mkpath(shots);
  const fs::path work = fs::absolute("gui_casual_work");
  fs::remove_all(work);
  fs::create_directories(work);
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, QString::fromStdString((work / "settings").string()));

  // (a) every example opens with no modal dialog; the pairs keep their keys
  std::vector<fs::path> projects;
  for (const auto &dir : fs::directory_iterator(source / "examples"))
    if (dir.is_directory())
      for (const auto &f : fs::directory_iterator(dir.path()))
        if (f.path().extension() == ".azr") projects.push_back(f.path());
  for (const char *t : {"tests/upgrade/T1/T1.azr", "tests/upgrade/T2/T2.azr", "tests/upgrade/T3/T3.azr", "tests/13N/13N.azr"})
    projects.push_back(source / t);
  std::sort(projects.begin(), projects.end());
  Check(projects.size() >= 14, "the examples and test projects found (" + std::to_string(projects.size()) + ")");
  const int kIr = azure2::front::LevelColumnIndex("ir"), kM1 = azure2::front::LevelColumnIndex("m1"),
            kM2 = azure2::front::LevelColumnIndex("m2"), kZ1 = azure2::front::LevelColumnIndex("z1"),
            kZ2 = azure2::front::LevelColumnIndex("z2");
  for (const fs::path &p : projects) {
    const std::string name = p.stem().string();
    const fs::path copy = CopyProject(p, work / ("open_" + name));
    AZURESetup w;
    w.resize(1280, 720);
    w.show();
    Settle(100);
    ModalWatch watch;
    w.open(QString::fromStdString(copy.string()));
    Settle(1500);
    Check(w.hasProject(), name + ": opened");
    Check(watch.seen.isEmpty(), name + ": no modal dialog on opening (" + std::to_string(watch.seen.size()) + ": " +
                                    watch.seen.join(" / ").left(200).toStdString() + ")");
    // the pairs at their keys: as many pairs as the largest key, and every row the Levels tab writes back carries the
    // masses and charges of the file's row
    std::vector<azure2::front::LevelRow> file, back;
    for (const auto &r : w.holder()->project().Levels())
      if (!r.tokens.empty()) file.push_back(r);  // the blank lines between the levels are no rows
    for (const auto &r : w.getLevelsTab()->levelRows())
      if (!r.tokens.empty()) back.push_back(r);
    int maxKey = 0;
    for (const auto &r : file)
      if (!r.tokens.empty()) maxKey = std::max(maxKey, (int)Num(r.tokens, kIr));
    Check(w.getPairsTab()->getPairsModel()->numPairs() == maxKey,
          name + ": " + std::to_string(w.getPairsTab()->getPairsModel()->numPairs()) + " pairs for the largest key " +
              std::to_string(maxKey));
    bool same = back.size() == file.size();
    for (std::size_t i = 0; same && i < file.size(); i++)
      for (int c : {kIr, kM1, kM2, kZ1, kZ2}) same = same && Num(file[i].tokens, c) == Num(back[i].tokens, c);
    Check(same, name + ": the rows written back keep each channel's pair (key, masses, charges)");
  }

  // (b) the notes bar
  {
    const fs::path n13 = CopyProject(source / "tests/13N/13N.azr", work / "notes_13N");
    const fs::path t1 = CopyProject(source / "tests/upgrade/T1/T1.azr", work / "notes_T1");
    AZURESetup w;
    w.resize(1280, 720);
    w.show();
    Settle(100);
    ModalWatch watch;
    w.open(QString::fromStdString(n13.string()));
    Settle(800);
    NotesBar *bar = w.notesBar();
    Check(bar && bar->isVisible() && bar->notes().size() == 1 && bar->notes().first().contains("warning"),
          "13N: one note of the engine's check (" + (bar ? bar->notes().join(" / ").toStdString() : std::string()) + ")");
    Check(!w.diagnosticsDockWidget()->isVisible(), "13N: warnings only, the Diagnostics dock stays closed");
    w.grab().save(shots + "/notes_13N_1280x720.png");
    bar->detailsAction()->trigger();
    Settle(200);
    Check(bar->detailsWindow() && bar->detailsWindow()->isVisible() && !bar->detailsWindow()->isModal() &&
              !QApplication::activeModalWidget(),
          "Details: a window that does not block the main one");
    bar->dismissAction()->trigger();
    Settle(100);
    Check(!bar->isVisible() && w.hasProject(), "Dismiss hides the bar");
    w.open(QString::fromStdString(t1.string()));
    Settle(800);
    Check(!bar->isVisible() && bar->notes().isEmpty(), "T1 (valid): no notes bar; the notes of 13N are gone");
    Check(watch.seen.isEmpty(), "no modal dialog in (b)");
  }

  // (d) the converter's repair of a line
  {
    const fs::path dir = work / "repair_13N";
    fs::remove_all(dir);
    fs::create_directories(dir);
    fs::copy(source / "tests/13N/data", dir / "data", fs::copy_options::recursive);
    std::ifstream in(source / "tests/convert/fixtures/tests/13N/13N.azr", std::ios::binary);
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string original = buffer.str();
    // the data line of meyer_84: "... 0 0.005 0 data/meyer_84.dat" -> without the shift's vary flag
    const std::string good = "0.005          0               data/meyer_84.dat", bad = "0.005          data/meyer_84.dat";
    std::string broken = original;
    const std::size_t at = broken.find(good);
    Check(at != std::string::npos, "(d) the meyer_84 line found in the old 13N");
    if (at != std::string::npos) broken.replace(at, good.size(), bad);
    std::ofstream(dir / "13N.azr", std::ios::binary) << broken;
    std::ofstream(dir / "ref.azr", std::ios::binary) << original;
    std::string refReport;
    const bool refOk = azure2::front::ConvertProject((dir / "ref.azr").string(), (dir / "ref.f3.azr").string(), refReport);
    Check(refOk, "(d) the unmodified old 13N converts");

    AZURESetup w;
    w.resize(1280, 720);
    w.show();
    Settle(100);
    int dialogs = 0, readings = 0, lineNumber = 0;
    QString shown, taken;
    QTimer driver;
    driver.setInterval(50);
    QObject::connect(&driver, &QTimer::timeout, [&]() {
      ConvertRepairDialog *d = qobject_cast<ConvertRepairDialog *>(QApplication::activeModalWidget());
      if (!d) return;
      dialogs++;
      readings = d->readingButtons().size();
      shown = d->lineText();
      d->grab().save(shots + "/convert_repair_dialog.png");
      if (!d->readingButtons().isEmpty()) d->readingButtons().first()->click();
      taken = d->lineText();
      d->retryButton()->click();
    });
    driver.start();
    QString report;
    const bool ok = w.convertOldProject(QString::fromStdString((dir / "13N.azr").string()),
                                        QString::fromStdString((dir / "13N.f3.azr").string()), report);
    driver.stop();
    (void)lineNumber;
    Check(dialogs == 1 && readings == 2 && shown.contains("0.005 data/meyer_84.dat"),
          "(d) one repair window: the line as written, two readings (" + std::to_string(dialogs) + " window(s), " +
              std::to_string(readings) + " readings)");
    Check(taken.contains("0.005 0 data/meyer_84.dat"), "(d) the first reading (the shift fixed) put into the line");
    Check(ok && w.hasProject() && w.holder()->path().endsWith("13N.f3.azr"), "(d) converted again and opened: " + report.toStdString());
    std::ifstream a(dir / "13N.f3.azr", std::ios::binary), b(dir / "ref.f3.azr", std::ios::binary);
    std::stringstream sa, sb;
    sa << a.rdbuf();
    sb << b.rdbuf();
    Check(!sa.str().empty() && sa.str() == sb.str(), "(d) the converted project is byte for byte that of the unmodified old 13N");
    std::ifstream o(dir / "13N.azr", std::ios::binary);
    std::stringstream so;
    so << o.rdbuf();
    Check(so.str() == broken, "(d) the old file is unchanged");
    Check(fs::exists(dir / "13N.repaired.azr"), "(d) the repaired copy beside it");
    Check(w.notesBar()->notes().join(" ").contains("repaired"), "(d) a note says which line was repaired");
  }

  std::cout << (failures ? "FAILED: " : "ALL PASSED") << (failures ? std::to_string(failures) : "") << std::endl;
  return failures ? 1 : 0;
}
