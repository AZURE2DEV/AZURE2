// The Qt setup program without a display (docs/upgrade/qt_gui.md section 10):
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
#include <deque>
#include <functional>
#include <map>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include <QAbstractButton>
#include <QApplication>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QMessageBox>
#include <QCheckBox>
#include <QFile>
#include <QFileInfo>
#include <QImage>
#include <QComboBox>
#include <QDir>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QSpinBox>
#include <QTabWidget>
#include <QTableWidget>
#include <QTimer>

#include "AZURESetup.h"
#include "ChannelsModel.h"
#include "ConfigDialogs.h"
#include "LevelsModel.h"
#include "DiagnosticsPanel.h"
#include "FittingTab.h"
#include "LevelsTab.h"
#include "MCMCTab.h"
#include "PairsModel.h"
#include "PairsTab.h"
#include "PlotTab.h"
#include "PlotWidget.h"
#include "RecordTabs.h"
#include "RecordsPanel.h"
#include "PythonTools.h"
#include "RunTab.h"
#include "ThmTab.h"

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
RunOutcome RunAndWait(AZURESetup &w, int &turns, std::function<void()> during = nullptr,
                      std::function<bool(QString &)> start = nullptr) {
  RunOutcome outcome;
  bool ended = false;
  QMetaObject::Connection c = QObject::connect(&w, &AZURESetup::runFinished, [&](const RunOutcome &o) {
    outcome = o;
    ended = true;
  });
  QString why;
  if (!(start ? start(why) : w.startRun(why))) {
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

/// A fresh folder `to` with output/ and checks/ and the files `names` of `from` (as tests/cli_commands' fixture).
fs::path Fixture(const fs::path &from, const std::vector<std::string> &names, const fs::path &to) {
  fs::remove_all(to);
  fs::create_directories(to / "output");
  fs::create_directories(to / "checks");
  for (const std::string &n : names) fs::copy(from / n, to / n, fs::copy_options::recursive);
  return to / names[0];
}

/// Every file of a's output directory equals b's, byte for byte; `report` names the ones that differ.
bool SameOutputs(const fs::path &a, const fs::path &b, std::string &report) {
  int same = 0, total = 0;
  report.clear();
  for (const auto &e : fs::directory_iterator(a / "output")) {
    total++;
    if (Bytes(e.path()) == Bytes(b / "output" / e.path().filename())) same++;
    else report += " " + e.path().filename().string();
  }
  report = std::to_string(same) + " of " + std::to_string(total) + " identical" + (report.empty() ? "" : ", differ:" + report);
  return total > 0 && same == total;
}

/// Answers the modal dialogs a menu item opens, one handler per dialog in order (a file chooser given its file, a
/// message box its button), polled from the event loop.
struct DialogDriver {
  std::deque<std::function<bool(QWidget *)>> answers;
  QTimer timer;
  DialogDriver() {
    timer.setInterval(100);
    QObject::connect(&timer, &QTimer::timeout, [this]() {
      QWidget *modal = QApplication::activeModalWidget();
      if (modal && !answers.empty() && answers.front()(modal)) answers.pop_front();
      else if (modal && std::getenv("AZURE2_TEST_DIALOGS"))
        std::cout << "  (a dialog without an answer: " << modal->metaObject()->className() << " "
                  << modal->windowTitle().toStdString() << ")" << std::endl;
    });
    timer.start();
  }
  void file(const QString &path) {
    answers.push_back([path](QWidget *w) {
      QFileDialog *d = qobject_cast<QFileDialog *>(w);
      if (!d) return false;
      // the file name typed into the dialog's name field and its Open/Save button pressed, as a user does
      QLineEdit *name = d->findChild<QLineEdit *>("fileNameEdit");
      QDialogButtonBox *buttons = d->findChild<QDialogButtonBox *>("buttonBox");
      if (!name || !buttons) return false;
      name->setText(path);
      for (QAbstractButton *b : buttons->buttons())
        if (buttons->buttonRole(b) == QDialogButtonBox::AcceptRole) {
          b->click();
          return true;
        }
      return false;
      return true;
    });
  }
  void button(QMessageBox::StandardButton b) {
    answers.push_back([b](QWidget *w) {
      QMessageBox *m = qobject_cast<QMessageBox *>(w);
      if (!m) return false;
      if (QAbstractButton *x = m->button(b)) x->click();
      else m->accept();
      return true;
    });
  }
};

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: gui_window_test <source dir> <AZURE2 binary> <screenshot dir>" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, true);  // the dialogs the tests answer are Qt's
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
        if (!tab->records(section)) continue;
        for (const azure2::front::Record &r : w.holder()->records(section)) {
          RecordsPanel *panel = tab->panelFor(section, QString::fromStdString(r.name));
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
    data->setFieldText("unit", "barn_per_furlong");
    data->applyForm(why);
    w.getDiagnostics()->validate();
    Check(w.getDiagnostics()->errorCount() > 0 && !w.getDiagnostics()->codes().isEmpty(),
          "T1: a wrong unit: the engine's diagnostic " + w.getDiagnostics()->codes().join(",").toStdString());
    Check(!w.readProjectFile(QString::fromStdString((source / "tests/convert/fixtures/tests/13N/13N.azr").string()), why) &&
              why.contains("--convert"),
          "an old project refused on opening, with the conversion");
  }

  // (b2) the THM tab's form: <thm> edited by the engine's edits of its lines, the keys from its tables
  {
    const fs::path copy = CopyProject(projects[1], work / "thm_T2");
    QString why;
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    const std::string original = w.holder()->project().Text();
    ThmTab *thm = w.getThmTab();
    Check(thm->experiments() == QStringList{"LC09"} && thm->selectExperiment("LC09"), "T2: the THM tab lists experiment[LC09]");
    QWidget *well = thm->experimentEditor("ffWell");
    const bool wellOffBefore = well && !well->isEnabled();
    const bool ff = thm->setExperimentKey("LC09", "vertexModel", "ff", why);
    bool found = false;
    for (const azure2::front::ThmEntry &e : azure2::front::ThmEntries(w.holder()->project()))
      found = found || (e.experiment == "LC09" && e.key == "vertexModel" && e.value == "ff");
    Check(ff && found && wellOffBefore && well->isEnabled(),
          "T2: vertexModel=ff set through the form; its field ffWell enabled only then (" + why.toStdString() + ")");
    Check(thm->setExperimentKey("LC09", "ffWell", "1.3,0.6", why), "T2: ffWell=1.3,0.6 set");
    Check(!thm->setExperimentKey("LC09", "vertexModel", "zz", why) && !why.isEmpty(),
          "T2: vertexModel=zz refused by the engine's reader: " + why.toStdString());
    Check(thm->addExperiment("second", why) && thm->experiments().size() == 2 && thm->removeExperiment("second", why) &&
              thm->experiments().size() == 1,
          "T2: an experiment added and removed");
    Check(thm->setOption("kinematics", "triple", why) && thm->setOption("kinematics", "", why),
          "T2: a global option set and removed");
    thm->setExperimentKey("LC09", "ffWell", "", why);
    thm->setExperimentKey("LC09", "vertexModel", "", why);
    const std::string back = w.holder()->project().Text();
    Check(back == original, "T2: every key removed again: the file as it was");
    if (back != original) std::cout << back.substr(back.find("<thm>"), 300) << std::endl;
  }

  // (b3) editing paths: records added, duplicated, renamed and removed in each list; a pair added and removed; levels
  // added, edited and removed with their channels; saved, reopened and validated
  {
    const fs::path copy = CopyProject(projects[1], work / "edit_T2");
    QString why;
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    struct List {
      RecordsPanel *panel;
      std::string section, kind;
    };
    const std::vector<List> lists = {{w.getSegmentsTab()->records("data"), "data", "dataset"},
                                     {w.getSegmentsTab()->records("extrapolations"), "extrapolations", "extrapolation"},
                                     {w.getEffectsTab()->records("targets"), "targets", "target"}};
    std::map<std::string, int> expected;
    for (const List &l : lists) {
      RecordsPanel *panel = l.panel;
      // the records of the section (the Experimental Effects tab lists each kind of <targets> on its own)
      const int before = panel ? (int)w.holder()->records(QString::fromStdString(l.section)).size() : -1;
      std::string first;
      for (const azure2::front::Record &r : w.holder()->records(QString::fromStdString(l.section)))
        if (first.empty() && r.kind == l.kind) first = r.name;
      QString copyName;
      // a list without a record of the kind (T2 has no extrapolation): only the add and the remove
      const bool duplicated = first.empty() || (panel && panel->selectRecord(QString::fromStdString(first)) &&
                                                panel->duplicateShown(copyName, why));
      const QString renamed = QString::fromStdString(first) + "_dup";
      const bool renamedOk = first.empty() || (duplicated && panel->renameShown(renamed, why) &&
                                               panel->recordNames().contains(renamed) && !panel->recordNames().contains(copyName));
      const int kept = panel ? panel->recordCount() : -1;
      const bool added = panel && panel->addNamedRecord(QString::fromStdString(l.kind), "added_one", why) &&
                         panel->recordCount() == kept + 1;
      const bool removed = panel && panel->selectRecord("added_one") && panel->removeShown(why) &&
                           panel->recordCount() == kept;
      Check(duplicated && renamedOk && added && removed,
            "T2 " + l.section + ": " + (first.empty() ? std::string("(none of kind ") + l.kind + ")" :
                                        first + " duplicated and renamed " + renamed.toStdString()) +
                ", a new " + l.kind + " added and removed (" + why.toStdString() + ")");
      expected[l.section] = before + (first.empty() ? 0 : 1);
    }
    Check(w.writeProjectFile(QString::fromStdString(copy.string()), why), "T2: saved after the record edits");
    azure2::front::Project back;
    std::string reason;
    const bool opened = azure2::front::Project::Open(copy.string(), back, reason);
    bool counts = opened;
    std::string found;
    for (const auto &e : expected) {
      counts = counts && (int)back.Records(e.first).size() == e.second;
      found += " " + e.first + " " + std::to_string(back.Records(e.first).size()) + "/" + std::to_string(e.second);
    }
    const azure2::front::Validation v = back.Validate();
    std::string errors;
    for (const azure2::front::Diagnostic &d : v.diagnostics)
      if (d.severity == "error") errors += " " + d.text;
    Check(opened && counts && v.ok, "T2: reopened with the new records (" + found + " ), valid" + errors);
  }
  {
    const fs::path copy = CopyProject(projects[0], work / "levels_T1");
    QString why;
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    PairsModel *pairs = w.getPairsTab()->getPairsModel();
    const int numPairs = pairs->numPairs();
    PairsData pair = pairs->getPairs().first();
    pair.channelRadius += 0.25;
    w.getPairsTab()->addPair(pair, numPairs, false);
    const bool pairAdded = pairs->numPairs() == numPairs + 1;
    w.getPairsTab()->removePairAt(numPairs);
    Check(pairAdded && pairs->numPairs() == numPairs, "T1: a pair added and removed (" + std::to_string(numPairs) + " pairs)");
    LevelsTab *levels = w.getLevelsTab();
    LevelsModel *lm = levels->getLevelsModel();
    const int numLevels = lm->getLevels().size();
    const int numChannels = levels->getChannelsModel()->getChannels().size();
    levels->addLevel(LevelsData{1, 0, 0.5, 1, 6.0}, false);
    levels->addLevel(LevelsData{1, 0, 2.5, -1, 7.0}, false);
    int newChannels = 0;
    for (const ChannelsData &c : levels->getChannelsModel()->getChannels()) newChannels += c.levelIndex >= numLevels;
    Check(lm->getLevels().size() == numLevels + 2 && newChannels > 0 &&
              levels->getChannelsModel()->getChannels().size() == numChannels + newChannels,
          "T1: two levels added (1/2+ at 6 MeV, 5/2- at 7 MeV) with " + std::to_string(newChannels) + " channels");
    lm->setData(lm->index(numLevels, 4), 6.5, Qt::EditRole);
    levels->removeLevelAt(numLevels + 1);
    int kept = 0;
    for (const ChannelsData &c : levels->getChannelsModel()->getChannels()) kept += c.levelIndex >= numLevels;
    Check(lm->getLevels().size() == numLevels + 1 && kept > 0 && kept < newChannels,
          "T1: the 1/2+ level moved to 6.5 MeV, the 5/2- level removed with its channels (" + std::to_string(kept) +
              " channels left of the new level)");
    Check(w.writeProjectFile(QString::fromStdString(copy.string()), why), "T1: saved after the level edits (" + why.toStdString() + ")");
    azure2::front::Project back;
    std::string reason;
    azure2::front::Project::Open(copy.string(), back, reason);
    bool found = false;
    for (const azure2::front::LevelRow &r : back.Levels())
      found = found || (r.tokens.size() > 2 && r.tokens[0] == "0.5" && r.tokens[1] == "1" && r.tokens[2] == "6.5");
    const azure2::front::Validation v = back.Validate();
    std::string errors;
    for (const azure2::front::Diagnostic &d : v.diagnostics)
      if (d.severity == "error") errors += " " + d.text;
    Check(found && v.ok, "T1: reopened with the level at 6.5 MeV, valid" + errors);
  }
  // the Configure dialogs, the formalism, File > New, Open Recent and Convert (through its dialogs)
  {
    const fs::path copy = CopyProject(projects[0], work / "config_T1");
    QString why;
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    const std::string original = w.holder()->project().Text();
    const QJsonArray table = QJsonDocument::fromJson(QByteArray::fromStdString(azure2::front::ConfigTableJson())).array();
    struct Dialog {
      QString title;
      QStringList groups, exclude;
      std::string shot;
    };
    for (const Dialog &d : {Dialog{"Checks", {"Checks"}, {}, "dialog_checks"},
                            Dialog{"Directories", {"Directories"}, {}, "dialog_directories"},
                            Dialog{"Runtime Options", {}, {"Checks", "Directories", "File"}, "dialog_runtime_options"}}) {
      ConfigDialog dialog(w.holder(), d.title, d.groups, d.exclude, &w);
      // the first entry of the dialog the test can change: a check's word, a directory's text, an option's box
      QString name, from, to;
      for (const QJsonValue &v : table) {
        const QJsonObject row = v.toObject();
        if (!dialog.editor(row["name"].toString()) || row["kind"].toString() == "version") continue;
        name = row["name"].toString();
        break;
      }
      QWidget *e = dialog.editor(name);
      const QString asRead = w.holder()->configValue(name, QString(QChar(1)));  // the line as it stands, or none
      if (QComboBox *c = qobject_cast<QComboBox *>(e)) {
        from = c->currentText();
        c->setCurrentIndex((c->currentIndex() + 1) % c->count());
        to = c->currentText();
      } else if (QCheckBox *b = qobject_cast<QCheckBox *>(e)) {
        from = b->isChecked() ? "1" : "0";
        b->setChecked(!b->isChecked());
        to = b->isChecked() ? "1" : "0";
      } else if (QLineEdit *l = qobject_cast<QLineEdit *>(e)) {
        from = l->text();
        l->setText("output2/");
        to = "output2/";
      }
      dialog.show();
      QApplication::processEvents();
      dialog.grab().save(shots + "/" + QString::fromStdString(d.shot) + ".png");
      const bool written = dialog.write(why) && w.holder()->configValue(name) == to;
      Check(!name.isEmpty() && written && w.holder()->project().Validate().format3,
            "Configure > " + d.title.toStdString() + ": " + name.toStdString() + " " + from.toStdString() + " -> " +
                to.toStdString() + " written by the engine's edit (" + why.toStdString() + ")");
      if (asRead == QString(QChar(1))) w.holder()->setConfig(name, QString(), why, true);
      else w.holder()->setConfig(name, asRead, why);
    }
    QAction *rMatrix = nullptr, *aMatrix = nullptr;
    for (QAction *a : w.findChildren<QAction *>()) {
      if (a->text() == "&R-Matrix") rMatrix = a;
      if (a->text() == "&A-Matrix") aMatrix = a;
    }
    const std::string amatrix = w.holder()->project().Config().front().value;
    if (rMatrix) rMatrix->trigger();
    const std::string r = w.holder()->project().Config().front().value;
    if (aMatrix) aMatrix->trigger();
    Check(rMatrix && aMatrix && amatrix == "true" && r == "false" && w.holder()->project().Config().front().value == "true",
          "Configure > Formalism: R-Matrix writes the first <config> line false, A-Matrix true again");
    // every entry back to its value (an edited line keeps the engine's layout of an edit: "true  #..." for "true #...")
    azure2::front::Project before = azure2::front::Project::FromText(original, copy.string());
    bool same = before.Config().size() == w.holder()->project().Config().size();
    for (std::size_t i = 0; same && i < before.Config().size(); i++)
      same = before.Config()[i].name == w.holder()->project().Config()[i].name &&
             before.Config()[i].value == w.holder()->project().Config()[i].value;
    Check(same, "T1: every Configure edit undone: the <config> entries as read");
    w.writeProjectFile(QString::fromStdString(copy.string()), why);

    // File > New
    QMetaObject::invokeMethod(&w, "newProject");
    Check(w.holder()->path().endsWith("untitled.azr") && w.holder()->records("data").empty() &&
              w.holder()->project().Levels().empty(),
          "File > New: an empty project (untitled.azr)");
    // Open Recent: the item of the T1 copy saved above
    QAction *recent = nullptr;
    for (QAction *a : w.findChildren<QAction *>())
      if (a->data().toString() == QFileInfo(QString::fromStdString(copy.string())).absoluteFilePath()) recent = a;
    if (recent) recent->trigger();
    Check(recent && w.holder()->path() == QFileInfo(QString::fromStdString(copy.string())).absoluteFilePath(),
          "File > Open Recent: the T1 copy opened from its item");
    // Convert through its dialogs: the old file, the new name, the converter's report, then the new file opens
    const fs::path conv = work / "convert_13N";
    fs::remove_all(conv);
    fs::create_directories(conv / "output");
    fs::create_directories(conv / "checks");
    fs::copy(source / "tests/convert/fixtures/tests/13N/13N.azr", conv / "13N.azr");
    fs::copy(source / "tests/13N/data", conv / "data", fs::copy_options::recursive);
    DialogDriver driver;
    driver.file(QString::fromStdString((conv / "13N.azr").string()));
    driver.file(QString::fromStdString((conv / "13N.f3.azr").string()));
    driver.button(QMessageBox::Ok);
    QMetaObject::invokeMethod(&w, "convertProject");
    const azure2::front::Validation v = w.holder()->project().Validate();
    std::string errors;
    for (const azure2::front::Diagnostic &d : v.diagnostics)
      if (d.severity == "error") errors += " " + d.text;
    Check(fs::exists(conv / "13N.f3.azr") && w.holder()->path().endsWith("13N.f3.azr") && v.ok,
          "File > Convert: the old 13N through the dialogs, the converted file opened and valid" + errors);
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

    // the Plot tab after the run: every data set, as many data points as the output files have rows
    PlotTab *plot = w.getPlotTab();
    const QStringList records = plot->availableRecords();
    plot->selectRecords(records);
    int rows = 0;
    for (const auto &e : fs::directory_iterator(copy.parent_path() / "output")) {
      if (e.path().filename().string().rfind("AZUREOut_aa=", 0) != 0 || e.path().extension() != ".out") continue;
      std::ifstream in(e.path());
      std::string line;
      while (std::getline(in, line)) rows += line.find_first_not_of(" \t") != std::string::npos;
    }
    int points = 0;
    for (const PlotSeries &s : plot->plot()->series()) points += s.dataX.size();
    Check(records.size() >= 7 && rows > 0 && points == rows,
          "T1: the plot holds " + std::to_string(points) + " data points of " + std::to_string(records.size()) +
              " data sets; the output files " + std::to_string(rows) + " rows");
    w.tabs()->setCurrentWidget(plot);
    plot->selectRecords({"cap_g0"});
    QApplication::processEvents();
    const QImage image = plot->plot()->grab().toImage();
    int ink = 0;
    for (int y = 0; y < image.height(); y += 2)
      for (int x = 0; x < image.width(); x += 2) ink += image.pixel(x, y) != qRgb(255, 255, 255);
    Check(ink > 500 && plot->plot()->dataPointsDrawn() > 0,
          "T1: the plot of cap_g0 is drawn (" + std::to_string(ink) + " coloured pixels of every fourth)");
    QDir().mkpath(shots);
    image.save(shots + "/plot_T1_cap_g0.png");
    // log axes, a zoom and its reset, residuals; the three export formats
    const QRectF whole = plot->plot()->view();
    plot->plot()->zoomTo(QRectF(whole.left(), whole.top(), whole.width() / 2, whole.height()));
    const bool zoomed = plot->plot()->view() != whole;
    plot->plot()->resetZoom();
    Check(zoomed && plot->plot()->view() == whole, "T1: the plot zooms and resets");
    bool exported = true;
    for (const char *suffix : {"png", "svg", "pdf"}) {
      const QString file = QString::fromStdString((work / "plot_T1.").string()) + suffix;
      exported = exported && plot->exportPlot(file, why) && QFileInfo(file).size() > 1000;
    }
    Check(exported, "T1: the plot exported to PNG, SVG and PDF (" + why.toStdString() + ")");
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

    // the parameter table: a session's parameters, a limit written as a <parameterSettings> row, the fit's values
    // beside the start, and the fitted values as the next start (--parameter-file)
    FittingTab *fitting = w.getFittingTab();
    QTableWidget *table = fitting->parameterTable();
    const bool readParameters = fitting->readParameters(QString::fromStdString(copy.string()), why);
    Check(readParameters && table->rowCount() > 10,
          "T1_start: the parameter table has " + std::to_string(table->rowCount()) + " rows (" + why.toStdString() + ")");
    int energyRow = -1;
    for (int r = 0; r < table->rowCount(); r++)
      if (table->item(r, 0)->text() == "energy_1") energyRow = r;
    Check(energyRow >= 0, "T1_start: energy_1 in the table");
    if (energyRow >= 0) {
      table->item(energyRow, 4)->setText("0.3");
      table->item(energyRow, 5)->setText("0.7");
    }
    const std::vector<azure2::front::ParameterSetting> settings = azure2::front::ParameterSettings(w.holder()->project());
    Check(settings.size() == 1 && settings[0].name == "energy_1" && settings[0].lower == 0.3 && settings[0].upper == 0.7 &&
              w.holder()->project().Validate().ok,
          "T1_start: the limits of energy_1 written as its <parameterSettings> row, the project valid");
    const bool readFitted = fitting->readFitted(fitting->fittedFile(), why);
    Check(readFitted && energyRow >= 0 && !table->item(energyRow, 9)->text().isEmpty(),
          "T1_start: the fitted value of energy_1 beside its start (" + (energyRow >= 0 ? table->item(energyRow, 9)->text().toStdString() : std::string()) + ")");
    if (energyRow >= 0) {  // the row removed again: the next runs are those of the file as it was
      table->item(energyRow, 4)->setText("");
      table->item(energyRow, 5)->setText("");
    }
    Check(azure2::front::ParameterSettings(w.holder()->project()).empty(), "T1_start: the emptied row removed");
    Check(fitting->takeFittedAsStart(fitting->fittedFile(), why) &&
              w.getRunTab()->options().join(" ").contains("--parameter-file="),
          "T1_start: the fitted values as the next start: --parameter-file= in the run's options");
    // the window's calculate from that start = the program's calculate --parameter-file=<the same file>, on a copy
    const fs::path cliStart = work / "cli_start_T1";
    fs::remove_all(cliStart);
    fs::copy(copy.parent_path(), cliStart, fs::copy_options::recursive);
    w.getRunTab()->setCommand("calculate");
    const RunOutcome c = RunAndWait(w, turns);
    const std::string command = "cd \"" + cliStart.string() + "\" && \"" + binary +
                                "\" calculate --threads=2 --parameter-file=T1_start.start.sav T1_start.azr > run.log 2>&1";
    const int cliStatus = std::system(command.c_str());
    std::string report;
    const bool same = SameOutputs(copy.parent_path(), cliStart, report);
    Check(c.exitCode == 0 && cliStatus == 0 && same,
          "T1_start: calculate from the fitted start = AZURE2 calculate --parameter-file=<it>: " + report);
    w.getRunTab()->presetOption("--parameter-file=", "");
  }

  // (c2) the other commands from the window, each against the program's run of the same command line on a copy:
  // every output file byte for byte (tests/cli_commands' projects and options)
  {
    struct Case {
      std::string name, command;
      fs::path from;
      std::vector<std::string> files;
      std::vector<std::pair<std::string, std::string>> options;  // the window's editors: option, value
      std::string cli;                                          // the program's options
    };
    const fs::path tests = source / "tests";
    // the program's runs below take --threads=2: so does the window's
    if (QSpinBox *threads = w.getRunTab()->findChild<QSpinBox *>("threads")) threads->setValue(2);
    const std::vector<Case> cases = {
        {"extrapolate_band", "extrapolate", tests / "upgrade/T1", {"T1.azr", "data"},
         {{"--covariance-band", "1"}}, "--covariance-band"},
        {"fit_band", "fit", tests / "identical_pp_res", {"identical_pp_res.azr", "data"},
         {{"--covariance-band", "1"}, {"--scale-covariance", "1"}}, "--covariance-band --scale-covariance"},
        {"minos", "minos", tests / "identical_pp_res", {"identical_pp_res.azr", "data"}, {{"--delta-chi2=", "1"}},
         "--delta-chi2=1"},
        {"rate", "rate", tests / "reaction_rate/models", {"p12C.azr"},
         {{"--rate-entrance=", "1"}, {"--rate-exit=", "2"}, {"--temperature-min=", "0.01"},
          {"--temperature-max=", "2"}, {"--temperature-step=", "0.1"}},
         "--rate-entrance=1 --rate-exit=2 --temperature-min=0.01 --temperature-max=2 --temperature-step=0.1"},
    };
    for (const Case &c : cases) {
      const fs::path mine = Fixture(c.from, c.files, work / ("window_" + c.name));
      const fs::path theirs = Fixture(c.from, c.files, work / ("cli_" + c.name));
      QString why;
      w.readProjectFile(QString::fromStdString(mine.string()), why);
      w.getRunTab()->setCommand(QString::fromStdString(c.command));
      bool set = true;
      for (const auto &o : c.options)
        set = set && w.getRunTab()->setOption(QString::fromStdString(o.first), QString::fromStdString(o.second));
      const RunOutcome o = RunAndWait(w, turns);
      const std::string command = "cd \"" + theirs.parent_path().string() + "\" && \"" + binary + "\" " + c.command +
                                  " " + c.cli + " --threads=2 " + theirs.filename().string() + " > run.log 2>&1";
      const int status = std::system(command.c_str());
      std::string report;
      const bool same = SameOutputs(mine.parent_path(), theirs.parent_path(), report);
      Check(set && o.exitCode == 0 && status == 0 && same,
            c.name + ": " + c.command + " from the window = AZURE2 " + c.command + " " + c.cli + ": " + report +
                (set ? "" : " (an option's editor was not set)"));
      for (const auto &opt : c.options) w.getRunTab()->setOption(QString::fromStdString(opt.first), opt.second == "1" ? "0" : "");
      if (c.name == "rate") {
        QComboBox *entrance = w.getRunTab()->findChild<QComboBox *>("option_rate-entrance");
        Check(entrance && entrance->count() == 1 + (int)w.getPairsTab()->getPairsModel()->numPairs(),
              "rate: the entrance pair is a drop-down of the project's " +
                  std::to_string(w.getPairsTab()->getPairsModel()->numPairs()) + " pairs");
        w.getRunTab()->setOption("--rate-entrance=", "");
        w.getRunTab()->setOption("--rate-exit=", "");
      }
      if (c.name == "fit_band") {
        const azure2::front::Curves curves = azure2::front::ReadCurves(w.holder()->project());
        int banded = 0;
        for (const azure2::front::Curve &k : curves.curves) banded += k.haveBand;
        Check(banded > 0 && banded == (int)curves.curves.size(),
              "fit_band: every data set's curve has its band (" + std::to_string(banded) + ")");
        w.getPlotTab()->selectRecords(w.getPlotTab()->availableRecords());
        w.tabs()->setCurrentWidget(w.getPlotTab());
        QApplication::processEvents();
        w.getPlotTab()->plot()->grab().save(shots + "/plot_identical_pp_res_band.png");
      }
    }
    // mcmc from its tab, where the engine has the command
    if (w.getMCMCTab()) {
      const fs::path mine = Fixture(tests / "mcmc_physical_norm", {"13N_fast.azr"}, work / "window_mcmc");
      const fs::path theirs = Fixture(tests / "mcmc_physical_norm", {"13N_fast.azr"}, work / "cli_mcmc");
      fs::copy(tests / "13N/data", mine.parent_path() / "data", fs::copy_options::recursive);
      fs::copy(tests / "13N/data", theirs.parent_path() / "data", fs::copy_options::recursive);
      QString why;
      w.readProjectFile(QString::fromStdString(mine.string()), why);
      MCMCTab *mcmc = w.getMCMCTab();
      mcmc->setWalkers(4);
      mcmc->setSteps(100);
      mcmc->setSpread(0.5);
      mcmc->setThreads(1);
      mcmc->setBasis("rwa");
      mcmc->setFreshStart(true);
      const RunOutcome o = RunAndWait(w, turns, nullptr, [&w](QString &why) { return w.startMCMC(why); });
      const std::string command = "cd \"" + theirs.parent_path().string() + "\" && \"" + binary +
                                  "\" mcmc --walkers=4 --steps=100 --parameter-spread=0.5 --energy-spread=1 "
                                  "--sampler-threads=1 --basis=rwa --overwrite-samples --threads=2 13N_fast.azr > run.log 2>&1";
      const int status = std::system(command.c_str());
      std::string report;
      const bool same = SameOutputs(mine.parent_path(), theirs.parent_path(), report);
      Check(o.exitCode == 0 && status == 0 && same, "mcmc: the MCMC tab's run = AZURE2 mcmc: " + report);
      Check(mcmc->lastStep() == 100 && mcmc->acceptance().startsWith("Acceptance fraction:"),
            "mcmc: the progress reached step " + std::to_string(mcmc->lastStep()) + " of 100; " +
                mcmc->acceptance().toStdString());
      mcmc->setFreshStart(false);
      w.getRunTab()->setOption("--overwrite-samples", "0");
    } else {
      Check(!azure2::front::HasCommand("mcmc"), "no MCMC tab: the engine has no mcmc command");
    }
  }

  // (e) the Tools menu: pyazr's commands as child processes; ENDF-6 export of T1 through the code its item calls
  {
    PythonTools *tools = w.tools();
    QString why;
    PythonTools::setInterpreter("/no/such/python");
    QMenu menu;
    tools->fillMenu(&menu, w.holder()->path());
    const QList<QAction *> offItems = menu.actions();
    Check(!offItems.isEmpty() && !offItems.first()->isEnabled() && offItems.first()->toolTip().contains("does not start"),
          "Tools: every item disabled with its reason when the interpreter does not start: " +
              (offItems.isEmpty() ? std::string() : offItems.first()->toolTip().toStdString()));
    const char *python = std::getenv("AZURE2_TEST_PYTHON");
    PythonTools::setInterpreter(python ? python : "python3");
    tools->fillMenu(&menu, w.holder()->path());
    const bool enabled = !menu.actions().isEmpty() && menu.actions().first()->isEnabled();
    Check(enabled, "Tools: the items enabled with " + PythonTools::interpreter().toStdString() + " (" +
                       (menu.actions().isEmpty() ? std::string() : menu.actions().first()->toolTip().toStdString()) + ")");
    const fs::path copy = CopyProject(projects[0], work / "endf_T1");
    const fs::path endf = copy.parent_path() / "T1.endf";
    int exitCode = -100;
    QStringList lines;
    QMetaObject::Connection c1 = QObject::connect(tools, &PythonTools::finished, [&](int code) { exitCode = code; });
    QMetaObject::Connection c2 =
        QObject::connect(tools, &PythonTools::line, [&](const QString &l, bool) { lines << l; });
    const bool started = enabled && tools->start({"endf", "write", QString::fromStdString(copy.string()),
                                                  QString::fromStdString(endf.string())},
                                                 QString::fromStdString(copy.parent_path().string()), why);
    QElapsedTimer timer;
    timer.start();
    while (started && exitCode == -100 && timer.elapsed() < 600000) QApplication::processEvents(QEventLoop::AllEvents, 50);
    QObject::disconnect(c1);
    QObject::disconnect(c2);
    Check(started && exitCode == 0 && fs::exists(endf) && fs::file_size(endf) > 1000 &&
              w.getRunTab()->log().contains("endf write"),
          "Tools: Export ENDF-6 of T1 (python -m pyazr endf write) exit " + std::to_string(exitCode) + ", " +
              (fs::exists(endf) ? std::to_string(fs::file_size(endf)) : std::string("no")) + " bytes, " +
              std::to_string(lines.size()) + " lines in the log pane" + (lines.isEmpty() ? "" : ": " + lines.last().toStdString()));
  }

  // (d) screenshots of every tab, T2 (data sets, targets, THM) open
  {
    const fs::path copy = CopyProject(projects[1], work / "shots_T2");
    QString why;
    w.readProjectFile(QString::fromStdString(copy.string()), why);
    w.getSegmentsTab()->records("data")->selectRecord("s02_diff090");
    w.getEffectsTab()->records("targets")->selectRecord("tgt_s05_implanted");
    // after a calculation (the Calculate tab's results, the Plot tab), the parameter tables filled
    w.getRunTab()->setCommand("calculate");
    const RunOutcome o = RunAndWait(w, turns);
    Check(o.exitCode == 0, "T2: calculate from the window for the screenshots");
    w.getPlotTab()->selectRecords({"s02_diff090", "s03_diff150"});
    w.getFittingTab()->readParameters(QString::fromStdString(copy.string()), why);
    if (w.getMCMCTab()) w.getMCMCTab()->readParameters(QString::fromStdString(copy.string()), why);
    w.getThmTab()->selectExperiment("LC09");
    for (const QString &old : QDir(shots).entryList({"tab*.png"})) QFile::remove(shots + "/" + old);
    w.saveTabScreenshots(shots);
    const int n = (int)QDir(shots).entryList({"tab*.png"}).size();
    Check(n == w.tabs()->count(), "screenshots of " + std::to_string(n) + " tabs in " + shots.toStdString());
    // the sub-pages that a tab grab does not show, the Tools menu, About
    auto shot = [&](QWidget *page, QWidget *tab, const QString &name) {
      if (QTabWidget *pages = qobject_cast<QTabWidget *>(page->parentWidget() ? page->parentWidget()->parentWidget() : nullptr))
        pages->setCurrentWidget(page);
      w.tabs()->setCurrentWidget(tab);
      QApplication::processEvents();
      w.grab().save(shots + "/" + name + ".png");
    };
    if (QTabWidget *pages = w.getThmTab()->findChild<QTabWidget *>("thmPages")) {
      pages->setCurrentIndex(1);
      shot(pages->currentWidget(), w.getThmTab(), "sub_thm_text");
      pages->setCurrentIndex(0);
    }
    if (QTabWidget *pages = w.getFittingTab()->findChild<QTabWidget *>("fittingPages"))
      for (int i = 1; i < pages->count(); i++) {
        pages->setCurrentIndex(i);
        shot(pages->currentWidget(), w.getFittingTab(), QString("sub_fitting_%1").arg(i));
      }
    if (w.getMCMCTab())
      if (QTabWidget *pages = w.getMCMCTab()->findChild<QTabWidget *>("mcmcPages")) {
        pages->setCurrentIndex(1);
        shot(pages->currentWidget(), w.getMCMCTab(), "sub_mcmc_settings");
      }
    QMenu tools;
    w.tools()->fillMenu(&tools, w.holder()->path());
    tools.popup(w.mapToGlobal(QPoint(200, 30)));
    QApplication::processEvents();
    tools.grab().save(shots + "/menu_tools.png");
    tools.close();
    Check(QFileInfo(shots + "/menu_tools.png").size() > 0, "a screenshot of the Tools menu");
  }

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
