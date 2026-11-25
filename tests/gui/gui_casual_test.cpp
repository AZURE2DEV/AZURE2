// The window for an occasional user (docs/upgrade/qt_gui.md section 19), without a display:
//   gui_casual_test <source dir> <screenshot dir>
//  (a) every shipped example (examples/*/ *.azr) and the test projects of the start page (tests/upgrade/T1-T3,
//      tests/13N) open through the window's own open path with no modal dialog at all; the pairs read from <levels>
//      keep their keys (c12c12_tumino2018 names one pair under two keys: each key keeps its pair, and the rows the
//      Levels tab would write back carry the masses and charges of the file);
//  (b) the notes bar: 13N (valid, 50 warnings) opens with one note of the engine's check and the Diagnostics dock
//      closed; Details shows every note in a window that is not modal; Dismiss hides the bar; a valid project opens
//      with no bar; opening another project clears the notes;
//  (c) the record tabs at 1280x720 (the PI: "the bottom part covers everything"): 13N's data-set table (16 records,
//      the notes bar shown) shows at least eight rows on open; the Form button closes the form and the table takes
//      its height; the choice is kept for the next window; the Levels tab's level list shows at least ten rows and its
//      empty Gamma Decays / Isospin panels are not shown in standard view;
//  (h) an old project of each era converted from the window under a decimal-comma locale: the bytes of the
//      conversion under "C";
//  (i) a converted project whose data files are elsewhere: one note with Locate Data Files..., located by name;
//  (d) the converter's repair of an old line: the old 13N with one data line short of its shift's vary flag; the
//      window shows the line, its number and the two readings; the first taken and converted again gives the project
//      the unmodified old 13N gives, byte for byte; the old file is unchanged, the repaired copy beside it, a note says
//      so;
//  (e) the expert switch: standard view by default; T1's Nuclear Potential and MCMC tabs (no records, no <mcmc>) and
//      Model Data sub-tab hidden, shown in expert view; a project with a potential record shows its tab in standard
//      view; the Calculate tab's --capture-amplitudes hidden in standard view, shown in expert view and when set;
//      the switch changes no byte of the project (saved before and after, compared) and is kept in the settings;
//  (f) the second review (qt_gui.md section 19.8): natural_neon's target of kind none shows no field marked required
//      that does not apply, its data-set table fits 1280 px; Tumino's channel table shows its Value column and its
//      extrapolation grids rounded for display (the record kept); the THM experiment's expert keys hidden unless given;
//      a phase wrapped at 180 drawn continuous; the multipolarity limits shown when a photon pair is added and the
//      Parameter Settings page when <parameterSettings> gets a line (standard view); the start page's test projects
//      in expert view only, a README title; New Project's first step; the notes bar's warning; Physics > Functions'
//      title line without a gap; the Overview's Calculate and Fit once before a run.
// PASS/FAIL per check; exit 1 on a failure.
#include <cstdlib>
#include <clocale>
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
#include "RecordExtras.h"
#include "RecordTabs.h"
#include "RecordsPanel.h"
#include "RunTab.h"
#include <QSplitter>
#include <QTabWidget>
#include <QTableView>
#include <QHeaderView>
#include <QToolButton>
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
#include "ChannelsModel.h"
#include "FittingTab.h"
#include "RecordTableModel.h"
#include "OverviewPage.h"
#include "PhysicsTab.h"
#include "PlotWidget.h"
#include "StartPage.h"
#include "OverviewCharts.h"
#include <QRegularExpression>
#include <QTextDocument>
#include <QGroupBox>
#include <QSet>
#include "ThmTab.h"
#include "ValueEditor.h"
#include <QListWidget>
#include <QScrollBar>
#include <QSpinBox>
#include <cmath>

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

/// The visible labels and placeholders under `root` that show a file's tag or a grammar ("<thm>", "<A><Symbol>",
/// "<value> <unit> ..."): what a user of the standard view should never read (qt_gui.md section 19.9).
QStringList AngleTexts(QWidget *root) {
  QStringList out;
  const QRegularExpression tag("<[A-Za-z][A-Za-z0-9_ ]*>");
  for (QLabel *l : root->findChildren<QLabel *>()) {
    if (!l->isVisibleTo(root) || l->text().isEmpty()) continue;
    QString plain = l->text();
    if (Qt::mightBeRichText(plain)) {
      QTextDocument d;
      d.setHtml(plain);
      plain = d.toPlainText();
    }
    if (tag.match(plain).hasMatch()) out << plain.left(60);
  }
  for (QLineEdit *e : root->findChildren<QLineEdit *>())
    if (e->isVisibleTo(root) && tag.match(e->placeholderText()).hasMatch()) out << "placeholder " + e->placeholderText().left(60);
  for (QGroupBox *b : root->findChildren<QGroupBox *>())
    if (b->isVisibleTo(root) && tag.match(b->title()).hasMatch()) out << "box " + b->title().left(60);
  return out;
}

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

  // (0) a decimal-comma locale of the user does not reach the readers: with LC_NUMERIC of such a language strtod
  // stops at the point and every "5.441" of a project is refused.  The window resets it; the rest of this test then
  // runs as the program would for such a user.  (Skipped where the machine has no such locale.)
  {
    const char *comma = nullptr;
    for (const char *name : {"de_DE.UTF-8", "de_DE.utf8", "it_IT.UTF-8", "pl_PL.UTF-8", "fr_FR.UTF-8", "de_DE", "German_Germany.1252", "de-DE"})
      if (std::setlocale(LC_NUMERIC, name) && std::string(std::localeconv()->decimal_point) == ",") {
        comma = name;
        break;
      }
    if (!comma) {
      std::setlocale(LC_NUMERIC, "C");
      std::cout << "skip (0): no decimal-comma locale on this machine" << std::endl;
    } else {
      AZURESetup probe;
      Check(std::string(std::localeconv()->decimal_point) == "." && std::strtod("5.441", nullptr) == 5.441,
            std::string("(0) with the locale ") + comma + " the window reads 5.441 as 5.441");
    }
  }

  // (a) every example opens with no modal dialog; the pairs keep their keys
  std::vector<fs::path> projects;
  for (const auto &dir : fs::directory_iterator(source / "examples"))
    if (dir.is_directory())
      for (const auto &f : fs::directory_iterator(dir.path()))
        if (f.path().extension() == ".azr") projects.push_back(f.path());
  for (const char *t : {"tests/upgrade/T1/T1.azr", "tests/upgrade/T2/T2.azr", "tests/upgrade/T3/T3.azr", "tests/13N/13N.azr"})
    projects.push_back(source / t);
  // every project of the examples' manifest (qt_gui.md section 19.9), read where the start page reads it
  std::vector<fs::path> manifest;
  const QList<StartPage::Example> listed = StartPage::readManifest(QString::fromStdString(source.string()));
  int entries = 0;
  QSet<QString> reactions;
  bool unique = true;
  for (const StartPage::Example &e : listed) {
    if (e.project.isEmpty()) continue;
    manifest.push_back(fs::path(e.project.toStdString()));
    if (!e.variant) {
      entries++;
      unique = unique && !reactions.contains(e.reaction);
      reactions.insert(e.reaction);
    }
    if (std::find(projects.begin(), projects.end(), manifest.back()) == projects.end()) projects.push_back(manifest.back());
  }
  Check(entries >= 10 && unique, "(a) the manifest lists " + std::to_string(entries) + " reactions, each once, and " +
                                     std::to_string(manifest.size() - entries) + " variants");
  for (const fs::path &m : manifest) Check(fs::exists(m), "(a) manifest: " + m.string() + " exists");
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
    if (std::find(manifest.begin(), manifest.end(), p) != manifest.end()) {
      // a project the start page offers: no note at all (no warning of the check), and it validates clean
      const azure2::front::Validation v = w.holder()->project().Validate();
      int problems = 0;
      for (const azure2::front::Diagnostic &d : v.diagnostics) problems += d.severity == "error" || d.severity == "warning";
      Check(w.notesBar()->notes().isEmpty() && v.ok && problems == 0,
            name + " (manifest): no note on opening, validates with no error or warning (" +
                w.notesBar()->notes().join(" / ").toStdString() + ")");
    }
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

  // (c) the record tabs at 1280x720
  {
    QSettings().setValue("view/expertOptions", false);
    const fs::path n13 = CopyProject(source / "tests/13N/13N.azr", work / "layout_13N");
    {
      AZURESetup w;
      w.resize(1280, 720);
      w.show();
      Settle(100);
      w.open(QString::fromStdString(n13.string()));
      Settle(500);
      w.tabs()->setCurrentWidget(w.getSegmentsTab());
      Settle(400);
      RecordsPanel *data = w.getSegmentsTab()->records("data");
      const int rows = data ? data->visibleRows() : 0;
      w.grab().save(shots + "/layout_13N_data_1280x720.png");
      Check(data && data->recordCount() == 16 && rows >= 8,
            "(c) 13N at 1280x720: the data-set table shows " + std::to_string(rows) + " rows on open (at least 8)");
      Check(data && data->formShown(), "(c) the form below the table is shown");
      data->setFormShown(false);
      Settle(200);
      const int closedRows = data->visibleRows();
      Check(!data->formShown() && closedRows > rows, "(c) Form closed: the table takes the height (" + std::to_string(closedRows) + " rows)");
      w.tabs()->setCurrentWidget(w.getLevelsTab());
      Settle(400);
      QTableView *levels = w.getLevelsTab()->findChild<QTableView *>("levelsView");
      const int levelRows = levels ? levels->viewport()->height() / std::max(1, levels->verticalHeader()->defaultSectionSize()) : 0;
      Check(levelRows >= 10, "(c) the Levels tab's list shows " + std::to_string(levelRows) + " rows (at least 10)");
      RecordExtras *extras = w.getLevelsTab()->extraPanels();
      Check(extras && !extras->tabs()->isVisible(), "(c) the empty Gamma Decays / Isospin panels are not shown in standard view");
      w.grab().save(shots + "/layout_13N_levels_1280x720.png");
    }
    {
      // the next window: the closed form is remembered
      AZURESetup w;
      w.resize(1280, 720);
      w.show();
      Settle(100);
      w.open(QString::fromStdString(n13.string()));
      Settle(500);
      w.tabs()->setCurrentWidget(w.getSegmentsTab());
      Settle(400);
      RecordsPanel *data = w.getSegmentsTab()->records("data");
      Check(data && !data->formShown(), "(c) the next window keeps the form closed");
      data->setFormShown(true);
      Settle(200);
      Check(data->formShown() && data->visibleRows() >= 5, "(c) the form opened again, the table keeps at least five rows");
    }
  }

  // (e) the expert switch
  {
    QSettings().remove("view/expertOptions");
    const fs::path t1 = CopyProject(source / "tests/upgrade/T1/T1.azr", work / "expert_T1");
    const fs::path ex = CopyProject(source / "experiment/tests/data/section_examples.azr", work / "expert_sections");
    auto bytes = [](const fs::path &p) {
      std::ifstream in(p, std::ios::binary);
      std::stringstream b;
      b << in.rdbuf();
      return b.str();
    };
    const std::string before = bytes(t1);
    AZURESetup w;
    w.resize(1280, 720);
    w.show();
    Settle(100);
    Check(!w.expertMode(), "(e) standard view by default");
    w.open(QString::fromStdString(t1.string()));
    Settle(500);
    QTabWidget *tabs = w.tabs();
    auto tabShown = [tabs](const QString &title) {
      for (int i = 0; i < tabs->count(); i++)
        if (tabs->tabText(i).remove('&') == title) return tabs->isTabVisible(i);
      return false;
    };
    QTabWidget *segmentTabs = w.getSegmentsTab()->findChild<QTabWidget *>();
    auto subShown = [segmentTabs](const QString &title) {
      for (int i = 0; segmentTabs && i < segmentTabs->count(); i++)
        if (segmentTabs->tabText(i) == title) return segmentTabs->isTabVisible(i);
      return false;
    };
    w.getRunTab()->setCommand("calculate");
    Check(!tabShown("Nuclear Potential") && !tabShown("MCMC") && !subShown("Model Data") && tabShown("Segments"),
          "(e) T1, standard view: Nuclear Potential, MCMC and Model Data hidden");
    Check(!w.getRunTab()->optionShown("--capture-amplitudes=") && w.getRunTab()->optionShown("--parameter-file="),
          "(e) Calculate, standard view: --capture-amplitudes hidden, --parameter-file shown");
    w.grab().save(shots + "/expert_off_T1.png");
    w.setExpertMode(true);
    Settle(200);
    Check(tabShown("Nuclear Potential") && tabShown("MCMC") && subShown("Model Data") &&
              w.getRunTab()->optionShown("--capture-amplitudes="),
          "(e) expert view: the hidden tabs, sub-tab and option shown");
    Check(QSettings().value("view/expertOptions").toBool(), "(e) the switch kept in the settings");
    w.grab().save(shots + "/expert_on_T1.png");
    w.setExpertMode(false);
    Settle(200);
    w.getRunTab()->setOption("--capture-amplitudes=", "intEC.dat");
    w.getRunTab()->setCommand("fit");
    Settle(100);
    Check(w.getRunTab()->optionShown("--capture-amplitudes="), "(e) standard view: an expert option that is set is shown");
    w.getRunTab()->setOption("--capture-amplitudes=", "");
    Check(!w.isWindowModified() && !w.holder()->project().Modified(), "(e) the switch marks nothing modified");
    QString why;
    Check(w.writeProjectFile(QString::fromStdString(t1.string()), why) && bytes(t1) == before,
          "(e) saved after switching both ways: the file is byte for byte as before " + why.toStdString());
    w.open(QString::fromStdString(ex.string()));
    Settle(500);
    Check(w.getSegmentsTab() && tabShown("Nuclear Potential") && subShown("Model Data"),
          "(e) a project with potential and model-data records: their tab and sub-tab shown in standard view");
    QSettings().setValue("view/expertOptions", false);
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
    // the line as in the file, its runs of blanks kept (the window showed it collapsed before qt_gui.md 19.8)
    Check(dialogs == 1 && readings == 2 && shown.contains(QString::fromStdString(bad)),
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

    // "Convert again" with the line unedited: the copy holds the line exactly as the file (tabs and runs of blanks)
    const fs::path dir2 = work / "repair_13N_unedited";
    fs::remove_all(dir2);
    fs::create_directories(dir2);
    fs::copy(source / "tests/13N/data", dir2 / "data", fs::copy_options::recursive);
    std::string tabbed = original;
    const std::size_t at2 = tabbed.find(good);
    const std::string badTab = "0.005  \t  data/meyer_84.dat";
    if (at2 != std::string::npos) tabbed.replace(at2, good.size(), badTab);
    std::ofstream(dir2 / "13N.azr", std::ios::binary) << tabbed;
    int rounds = 0;
    QString shownTab;
    QTimer driver2;
    driver2.setInterval(50);
    QObject::connect(&driver2, &QTimer::timeout, [&]() {
      ConvertRepairDialog *d = qobject_cast<ConvertRepairDialog *>(QApplication::activeModalWidget());
      if (!d) return;
      if (++rounds == 1) {
        shownTab = d->lineText();
        d->retryButton()->click();  // unedited
      } else {
        d->reject();
      }
    });
    driver2.start();
    QString report2;
    const bool ok2 = w.convertOldProject(QString::fromStdString((dir2 / "13N.azr").string()),
                                         QString::fromStdString((dir2 / "13N.f3.azr").string()), report2);
    driver2.stop();
    std::ifstream c(dir2 / "13N.repaired.azr", std::ios::binary);
    std::stringstream sc;
    sc << c.rdbuf();
    Check(shownTab.contains(QString::fromStdString(badTab)), "(d) the line shown with its tab and blanks as in the file");
    Check(!ok2 && rounds == 2 && sc.str() == tabbed,
          "(d) Convert again unedited: the copy is byte for byte the old file (" + std::to_string(rounds) + " rounds)");
  }

  // (f) the second review (qt_gui.md section 19.8)
  {
    QSettings().setValue("view/expertOptions", false);
    // a target of kind none (natural_neon): its active nucleus, areal density, stopping and straggling do not apply,
    // so no field is marked required and no red line is shown
    const fs::path neon = CopyProject(source / "examples/natural_neon/natural_neon.azr", work / "review_neon");
    AZURESetup w;
    w.resize(1280, 720);
    w.show();
    Settle(100);
    ModalWatch watch;
    w.open(QString::fromStdString(neon.string()));
    Settle(800);
    w.tabs()->setCurrentWidget(w.getEffectsTab());
    Settle(300);
    RecordsPanel *targets = w.getEffectsTab()->records("targets");
    if (targets) targets->selectRecord("natural");
    Settle(300);
    QStringList red;
    for (const char *key : {"active", "areal_density", "stopping", "straggling"})
      if (FieldEditor *e = targets ? targets->fieldEditor(key) : nullptr)
        if (!e->message().isEmpty()) red << QString(key) + ": " + e->message().left(40);
    Check(targets && !targets->fieldApplies("active") && !targets->fieldApplies("areal_density") && red.isEmpty(),
          "(f) natural_neon, target of kind none: active nucleus and areal density do not apply, no red line (" +
              red.join(" / ").toStdString() + ")");
    Check(targets && !targets->fieldRowShown("active") && !targets->fieldRowShown("areal_density"),
          "(f) natural_neon: the rows that do not apply are not shown while empty");
    int isotopeDocs = 0;  // the isotope's doc once (its box), not on every entry
    if (targets)
      for (QLabel *l : targets->findChildren<QLabel *>())
        isotopeDocs += l->isVisibleTo(targets) && l->text().contains("abundance weight of a target nucleus");
    Check(isotopeDocs <= 1, "(f) natural_neon: the isotope's help line shown " + std::to_string(isotopeDocs) + " time(s) (at most once)");
    const QStringList angles = targets ? AngleTexts(targets) : QStringList{"no panel"};
    Check(angles.isEmpty(), "(f) natural_neon target form, standard view: no tag or grammar shown (" + angles.join(" | ").toStdString() + ")");
    w.grab().save(shots + "/review_neon_target.png");
    // the data-set table at 1280 px: no sideways scrolling
    w.tabs()->setCurrentWidget(w.getSegmentsTab());
    Settle(400);
    QTableView *dataTable = w.getSegmentsTab()->records("data") ? w.getSegmentsTab()->records("data")->table() : nullptr;
    Check(dataTable && dataTable->horizontalScrollBar()->maximum() == 0,
          "(f) natural_neon at 1280x720: the data-set table fits without sideways scrolling (" +
              std::to_string(dataTable ? dataTable->horizontalScrollBar()->maximum() : -1) + ")");
    w.grab().save(shots + "/review_neon_data.png");
    // the Overview: Calculate and Fit once on the page before a run
    OverviewPage *overview = w.findChild<OverviewPage *>();
    QPushButton *stepCalc = overview ? overview->findChild<QPushButton *>("overviewAction_runcalculate") : nullptr;
    Check(stepCalc && stepCalc->isHidden() && overview->findChild<QPushButton *>("overviewNoRun_calculate"),
          "(f) Overview before a run: Calculate and Fit once (in The last run), not twice");
    Check(watch.seen.isEmpty(), "(f) natural_neon: no modal dialog");
  }
  {
    // Tumino: the channel table's Value column readable; the extrapolation grids' display rounded, the file kept
    const fs::path tumino = CopyProject(source / "examples/c12c12_tumino2018/c12c12_tumino2018.azr", work / "review_tumino");
    AZURESetup w;
    w.resize(1280, 720);
    w.show();
    Settle(100);
    w.open(QString::fromStdString(tumino.string()));
    Settle(800);
    w.tabs()->setCurrentWidget(w.getLevelsTab());
    Settle(300);
    QTableView *channels = w.getLevelsTab()->findChild<QTableView *>("channelsView");
    const int valueWidth = channels ? channels->columnWidth(ChannelsModel::ValueColumn) : 0;
    Check(valueWidth >= 70, "(f) Tumino at 1280x720: the channel table's Value column is " + std::to_string(valueWidth) +
                                " px wide (at least 70)");
    RecordsPanel *extrap = w.getSegmentsTab()->records("extrapolations");
    QString shown, stored;
    if (extrap)
      for (int c = 0; c < extrap->model()->columnCount(); c++) {
        const QString t = extrap->model()->data(extrap->model()->index(0, c)).toString();
        if (t.startsWith("list ")) shown = t;
      }
    for (const azure2::front::Record &r : w.holder()->project().Records("extrapolations"))
      if (r.name == "x1") stored = QString::fromStdString(r.Value("energies"));
    Check(shown.contains(" 5.5 MeV") && !shown.contains("5.4999"),
          "(f) Tumino: the grid is shown \"" + shown.toStdString() + "\"");
    // the form of x1: the grid's line shows the rounded numbers, the field still gives the file's text
    if (extrap) {
      extrap->selectRecord("x1");
      Settle(200);
    }
    FieldEditor *energies = extrap ? extrap->fieldEditor("energies") : nullptr;
    bool roundedLine = false;
    if (energies)
      for (QLineEdit *l : energies->findChildren<QLineEdit *>())
        roundedLine = roundedLine || (l->isVisibleTo(energies) && l->text().contains(" 5.5") && !l->text().contains("5.4999"));
    Check(energies && roundedLine && energies->text() == stored,
          "(f) Tumino x1 form: the grid shown rounded, the field's value the file's text until edited");
    Check(stored.endsWith("5.49999999999997 MeV lab"),
          "(f) Tumino: the record keeps the grid as written (" + stored.right(30).toStdString() + ")");
    w.grab().save(shots + "/review_tumino_levels.png");
    // the THM experiment: the keys of an occasional user shown, the others only when given (or in expert view)
    ThmTab *thm = w.getThmTab();
    if (thm) thm->selectExperiment("E1");
    Check(thm && thm->rowShown("beam", true) && thm->rowShown("Ebeam", true) && thm->rowShown("ps", true) &&
              !thm->rowShown("psNodes", true) && !thm->rowShown("opticalAA", true),
          "(f) THM, standard view: beam, Ebeam, ps shown; psNodes, opticalAA (not given) hidden");
    w.setExpertMode(true);
    Settle(100);
    Check(thm && thm->rowShown("psNodes", true) && thm->rowShown("opticalAA", true), "(f) THM, expert view: every key");
    w.setExpertMode(false);
    Settle(100);
    const QStringList thmAngles = thm ? AngleTexts(thm) : QStringList{"no tab"};
    Check(thmAngles.isEmpty(), "(f) THM tab, standard view: no tag or grammar shown (" + thmAngles.join(" | ").toStdString() + ")");
    w.setExpertMode(true);
    Settle(100);
    w.setExpertMode(false);
    Settle(100);
    // Physics > R-matrix: the E_x fields show their whole value
    w.tabs()->setCurrentIndex(w.tabs()->count() - 1);
    if (PhysicsTab *physics = w.findChild<PhysicsTab *>()) physics->showPage(PhysicsTab::RMatrix);
    // The page fills its fields when its session has read the project: on a slow machine that takes far more than a
    // fixed pause, so wait for the two fields (up to 3 minutes) instead.
    for (int waited = 0; waited < 180000; waited += 500) {
      Settle(500);
      int filled = 0;
      for (QLineEdit *e : w.findChildren<QLineEdit *>())
        filled += e->isVisible() && (e->objectName() == "sweepFrom" || e->objectName() == "sweepTo") && !e->text().isEmpty();
      if (filled == 2 && waited >= 3500) break;
    }
    int cut = 0, fields = 0;
    for (QLineEdit *e : w.findChildren<QLineEdit *>()) {
      if (!e->isVisible() || (e->objectName() != "sweepFrom" && e->objectName() != "sweepTo") || e->text().isEmpty()) continue;
      fields++;
      cut += e->fontMetrics().horizontalAdvance(e->text()) + 6 > e->width();
    }
    Check(fields == 2 && cut == 0, "(f) Tumino, Physics > R-matrix: the E_x fields show their values whole (" +
                                       std::to_string(fields) + " fields, " + std::to_string(cut) + " cut)");
    w.grab().save(shots + "/review_tumino_rmatrix.png");
  }
  {
    // the legend where it hides no point, else below the axes; pole labels that never overlap
    PlotWidget plot;
    plot.resize(700, 420);
    plot.show();
    QList<PlotSeries> full, corner;
    for (int k = 0; k < 6; k++) {
      PlotSeries a;
      a.label = QString("channel %1 (pair 6, l 4, s 0) -> channel %1 (pair 6, l 4, s 0)").arg(k + 1);
      for (int i = 0; i <= 100; i++) {
        a.modelX << 14.0 + 0.035 * i;
        a.modelY << (k % 2 ? 1.0 : -1.0) * (i % 2 ? 1.0 : 0.2);  // fills the panel, every corner
      }
      full << a;
    }
    plot.setSeries(full);
    plot.setMarkers({14.7694, 15.1080, 15.1445, 15.1640, 15.2316, 15.3549, 15.3866, 15.8883, 16.2731, 16.3362},
                    {"14.769448 MeV", "15.108045 MeV", "15.144539 MeV", "15.164 MeV", "15.2316 MeV", "15.3549 MeV",
                     "15.3866 MeV", "15.888346 MeV", "16.273096 MeV", "16.336232 MeV"});
    plot.grab();
    Check(plot.legendBelowAxes(), "(f) plot: curves in every corner, the legend below the axes");
    const QList<QRect> labels = plot.markerLabelRects();
    bool apart = !labels.isEmpty();
    for (int i = 0; i < labels.size(); i++)
      for (int j = i + 1; j < labels.size(); j++) apart = apart && !labels[i].intersects(labels[j]);
    Check(apart, "(f) plot: " + std::to_string(labels.size()) + " pole labels drawn, none over another");
    plot.grab().save(shots + "/review_plot_legend_below.png");
    PlotSeries low;
    low.label = "low";
    for (int i = 0; i <= 50; i++) low.modelX << i, low.modelY << (i < 10 ? 1.0 : 0.0);
    corner << low;
    plot.setSeries(corner);
    plot.setMarkers({});
    plot.grab();
    Check(!plot.legendRect().isEmpty() && !plot.legendBelowAxes(), "(f) plot: a free corner, the legend inside the panel");
    // the chi^2/N bars: one outlier cuts the scale at ten times the median, the others keep their lengths
    bool capped = false;
    const double scale = BarChart::ScaleTop({1.0, 1.2, 0.8, 0.9, 500.0}, &capped);
    bool plain = true;
    const double top = BarChart::ScaleTop({1.0, 1.2, 0.8, 2.0}, &plain);
    Check(capped && std::fabs(scale - 10.0) < 1e-12 && !plain && std::fabs(top - 2.0) < 1e-12,
          "(f) chi^2/N bars: an outlier caps the scale at ten times the median (" + std::to_string(scale) + "), else the largest");
  }
  {
    // a phase near 0 that the convention [0, 180) wraps to 180: drawn continuous
    const QVector<double> wrapped = {0.5, 179.8, 0.2, 179.9, 1.0, 60.0, 120.0, 170.0, 10.0};
    const QVector<double> u = physicsqt::UnwrappedPhaseDegrees(wrapped);
    bool smooth = u.size() == wrapped.size();
    for (int i = 1; smooth && i < u.size(); i++) smooth = std::fabs(u[i] - u[i - 1]) <= 90;
    Check(smooth && std::fabs(u[1] + 0.2) < 1e-9 && std::fabs(u[8] - 190.0) < 1e-9,
          "(f) the phase drawn continuous (-0.2 for 179.8, the last step through 180 to 190)");
  }
  {
    // the multipolarity limits: shown in standard view as soon as a photon pair is added
    const fs::path o18 = CopyProject(source / "examples/o18_lacognata2008/o18_lacognata2008.azr", work / "review_o18");
    AZURESetup w;
    w.resize(1280, 720);
    w.show();
    Settle(100);
    w.open(QString::fromStdString(o18.string()));
    Settle(600);
    QSpinBox *mult = w.getLevelsTab()->findChild<QSpinBox *>("maxMultSpin");
    bool photon = false;
    for (const PairsData &p : w.getPairsTab()->getPairsModel()->getPairs()) photon = photon || p.pairType == 10;
    Check(mult && !photon && mult->isHidden(), "(f) o18: no photon pair, the multipolarity limits hidden in standard view");
    PairsData g = w.getPairsTab()->getPairsModel()->getPairs().first();
    g.pairType = 10;
    g.lightJ = 1, g.lightPi = 1, g.lightZ = 0, g.lightM = 0, g.lightG = 0;
    w.getPairsTab()->addPair(g, w.getPairsTab()->getPairsModel()->numPairs(), false);
    Settle(200);
    Check(mult && !mult->isHidden(), "(f) a photon pair added in standard view: the multipolarity limits shown at once");
    // a row of <parameterSettings> written while its page is hidden: the page is shown
    QTabWidget *pages = w.getFittingTab()->findChild<QTabWidget *>("fittingPages");
    auto pageShown = [pages](const QString &prefix) {
      for (int i = 0; pages && i < pages->count(); i++)
        if (pages->tabText(i).startsWith(prefix)) return pages->isTabVisible(i);
      return false;
    };
    const bool before = pageShown("Parameter Settings");
    w.holder()->setSectionBody("parameterSettings", {"# a comment", "R1.E fixed"});
    Settle(100);
    Check(!before && pageShown("Parameter Settings"),
          "(f) Fitting Settings: the text page hidden while empty, shown once <parameterSettings> holds a line");
  }
  {
    // the start page: the manifest's groups as headings, one line per reaction with its source, the variants
    // indented under it, every project line openable (qt_gui.md section 19.9)
    AZURESetup w;
    w.resize(1280, 720);
    w.show();
    Settle(200);
    QListWidget *list = w.startPage()->exampleList();
    int headings = 0, lines = 0, variants = 0;
    bool tumino = false, o18variant = false;
    for (int i = 0; i < list->count(); i++) {
      const QListWidgetItem *it = list->item(i);
      const QString file = it->data(Qt::UserRole).toString();
      if (file.isEmpty()) {
        headings += it->flags() == Qt::NoItemFlags;
        continue;
      }
      lines++;
      const bool variant = it->text().trimmed().startsWith(QChar(0x21B3));
      variants += variant;
      tumino = tumino || (it->text().contains("Tumino et al. 2018") && file.endsWith("c12c12_tumino2018.azr"));
      o18variant = o18variant || (variant && file.endsWith("o18_lacognata2008.azr"));
    }
    Check(headings == 4 && lines == (int)StartPage::findExamples().size() && variants >= 5 && tumino && o18variant,
          "(f) start page: " + std::to_string(headings) + " groups, " + std::to_string(lines) + " project lines (" +
              std::to_string(variants) + " variants under their reaction), names from the manifest");
    w.grab().save(shots + "/review_start.png");
    // New Project: the next step is a particle pair, not a run
    QMetaObject::invokeMethod(&w, "newProject");
    Settle(500);
    OverviewPage *overview = w.findChild<OverviewPage *>();
    QString first;
    if (overview && !overview->steps().isEmpty()) first = overview->steps().first().text;
    bool runStep = false;
    for (const OverviewPage::Step &st : overview ? overview->steps() : QList<OverviewPage::Step>()) runStep = runStep || st.text.contains("no run yet");
    Check(first.contains("particle pair") && !runStep, "(f) New Project: the first step is \"" + first.toStdString() + "\"");
    w.grab().save(shots + "/review_new_overview.png");
    // the notes bar: a warning goes first and changes the icon
    w.notesBar()->clearNotes();
    w.notesBar()->addNote("an information");
    w.notesBar()->addNote("a warning", QString(), false, true);
    Check(w.notesBar()->hasWarning() && w.notesBar()->notes().first() == "a warning", "(f) notes bar: a warning first");
    w.notesBar()->clearNotes();
    Check(!w.notesBar()->hasWarning() && w.notesBar()->notes().isEmpty(), "(f) notes bar: cleared, no warning left");
    // Physics > Functions: no blank band under the function's title
    w.tabs()->setCurrentIndex(w.tabs()->count() - 1);
    Settle(300);
    QLabel *doc = w.findChild<QLabel *>("functionsDoc");
    Check(doc && doc->height() <= 3 * doc->fontMetrics().height(),
          "(f) Physics > Functions: the title line is " + std::to_string(doc ? doc->height() : -1) + " px high (no gap)");
  }

  // (g) a project the start page offers, saved after every tab, sub-tab and record form was shown, is byte for byte the
  //     save made right after opening it: showing a value (rounded numbers, qt_gui.md section 19.9) never edits it
  {
    for (const StartPage::Example &e : StartPage::readManifest(QString::fromStdString(source.string()))) {
      if (e.project.isEmpty()) continue;
      const fs::path p = e.project.toStdString();
      const std::string name = p.stem().string();
      const fs::path copy = CopyProject(p, work / ("visit_" + name));
      AZURESetup w;
      w.resize(1280, 720);
      w.show();
      Settle(100);
      ModalWatch watch;
      QString why;
      w.open(QString::fromStdString(copy.string()));
      Settle(600);
      const QString a = QString::fromStdString((copy.parent_path() / "a.azr").string());
      const QString b = QString::fromStdString((copy.parent_path() / "b.azr").string());
      const bool savedA = w.writeProjectFile(a, why);
      int forms = 0;
      for (int i = 0; i < w.tabs()->count(); i++) {
        if (!w.tabs()->isTabVisible(i)) continue;
        w.tabs()->setCurrentIndex(i);
        Settle(80);
        for (QTabWidget *inner : w.tabs()->widget(i)->findChildren<QTabWidget *>())
          for (int j = 0; j < inner->count(); j++) {
            inner->setCurrentIndex(j);
            Settle(40);
          }
        for (RecordsPanel *r : w.tabs()->widget(i)->findChildren<RecordsPanel *>())
          for (const QString &record : r->recordNames()) {
            forms += r->selectRecord(record);
            Settle(20);
          }
      }
      const bool savedB = w.writeProjectFile(b, why);
      std::ifstream fa(a.toStdString(), std::ios::binary), fb(b.toStdString(), std::ios::binary);
      std::stringstream sa, sb;
      sa << fa.rdbuf();
      sb << fb.rdbuf();
      Check(savedA && savedB && !sa.str().empty() && sa.str() == sb.str() && watch.seen.isEmpty(),
            "(g) " + name + ": saved after showing every tab and " + std::to_string(forms) +
                " record forms, byte for byte the save made on opening " + why.toStdString());
    }
  }

  // (h) an old project converted from the window under a decimal-comma locale of the user: the bytes the converter
  //     writes under "C" (the window resets LC_NUMERIC when it is made; the converter reads with strtod and streams).
  //     The forms of 2025-08 and 2026-05 (decimals in norms, shifts, angles, scalings) and of 2025-01 (coefficients).
  {
    const char *comma = nullptr;
    for (const char *name : {"de_DE.UTF-8", "de_DE.utf8", "it_IT.UTF-8", "pl_PL.UTF-8", "fr_FR.UTF-8", "de_DE"})
      if (std::setlocale(LC_NUMERIC, name) && std::string(std::localeconv()->decimal_point) == ",") {
        comma = name;
        break;
      }
    std::setlocale(LC_NUMERIC, "C");
    if (!comma) {
      std::cout << "skip (h): no decimal-comma locale on this machine" << std::endl;
    } else {
      for (const char *form : {"energy_shift_2025_08", "components_scaling_2026_05", "target_2025_01_coefficients"}) {
        const fs::path dir = work / (std::string("locale_") + form);
        fs::remove_all(dir);
        fs::create_directories(dir);
        fs::copy(source / "tests/convert/fixtures/history" / form, dir, fs::copy_options::recursive);
        const fs::path old = dir / (std::string(form) + ".azr");
        std::string report;
        const bool refOk = azure2::front::ConvertProject(old.string(), (dir / "ref.f3.azr").string(), report);
        std::setlocale(LC_NUMERIC, comma);  // as QApplication sets it from the user's environment
        bool ok = false;
        {
          AZURESetup w;  // resets LC_NUMERIC, as the program's window does
          QString why;
          ok = w.convertOldProject(QString::fromStdString(old.string()), QString::fromStdString((dir / "gui.f3.azr").string()), why);
        }
        std::setlocale(LC_NUMERIC, "C");
        std::ifstream a(dir / "ref.f3.azr", std::ios::binary), b(dir / "gui.f3.azr", std::ios::binary);
        std::stringstream sa, sb;
        sa << a.rdbuf();
        sb << b.rdbuf();
        Check(refOk && ok && !sa.str().empty() && sa.str() == sb.str(),
              std::string("(h) ") + form + " converted from the window with LC_NUMERIC " + comma +
                  ": the bytes of the conversion under C");
      }
    }
  }

  // (i) a converted project whose data files are not where it says (an old project from another computer, a moved
  //     folder): one note naming them with "Locate Data Files...", no Diagnostics dock, no box; the folder chosen, each
  //     file found by its name, its data set pointed to it, the project valid
  {
    const fs::path dir = work / "locate_2020";
    fs::remove_all(dir);
    fs::create_directories(dir);
    fs::copy(source / "tests/convert/fixtures/history/gui_2020", dir, fs::copy_options::recursive);
    std::string report;
    const bool converted = azure2::front::ConvertProject((dir / "gui_2020.azr").string(), (dir / "gui_2020.f3.azr").string(), report);
    fs::create_directories(dir / "elsewhere");
    fs::rename(dir / "data", dir / "elsewhere" / "measured");
    fs::create_directories(dir / "output");
    fs::create_directories(dir / "checks");
    AZURESetup w;
    w.resize(1280, 720);
    w.show();
    Settle(100);
    ModalWatch watch;
    w.open(QString::fromStdString((dir / "gui_2020.f3.azr").string()));
    Settle(800);
    const QStringList missing = w.missingDataFiles();
    Check(converted && missing.size() == 2, "(i) the converted 2020 project, its data folder moved: 2 data files missing (" +
                                                std::to_string(missing.size()) + ")");
    Check(w.notesBar()->isVisible() && w.notesBar()->notes().join(" ").contains("2 data file(s) not found") &&
              w.notesBar()->locateAction()->isVisible(),
          "(i) one note names them, with Locate Data Files...");
    Check(!w.diagnosticsDockWidget()->isVisible() && watch.seen.isEmpty(), "(i) no Diagnostics dock and no box for missing data files");
    w.grab().save(shots + "/locate_data_files_note.png");
    QStringList still;
    const int found = w.locateDataFiles(QString::fromStdString((dir / "elsewhere").string()), &still);
    Check(found == 2 && still.isEmpty() && w.missingDataFiles().isEmpty() && w.holder()->project().Modified(),
          "(i) located under the chosen folder: both data sets pointed to their files, the project modified");
    const azure2::front::Validation v = w.holder()->project().Validate();
    bool dataError = false;
    for (const azure2::front::Diagnostic &d : v.diagnostics) dataError = dataError || d.code == "E_DATA_FILE";
    Check(!dataError, "(i) after locating: no E_DATA_FILE");
  }

  std::cout << (failures ? "FAILED: " : "ALL PASSED") << (failures ? std::to_string(failures) : "") << std::endl;
  return failures ? 1 : 0;
}
