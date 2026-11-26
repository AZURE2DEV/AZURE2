// The Levels and Channels and the Particle Pairs tabs on the go (docs/upgrade/qt_gui.md section 11), headless:
//   gui_levels_test <source dir> [<screenshot dir>]
// Each action is driven through the widgets (keys typed into the add bar and the cells' editors, Del, Ctrl+D,
// Ctrl+Z), then the project is saved and its text read back:
//  (a) T1, T3 and 13N opened and saved without an edit: byte for byte; the first level selected, its channels shown;
//  (b) a level added from the add bar; a refused energy; an energy, a J^pi and a width edited in place, a refused
//      value leaving the model as it was; the value column's meaning by kind (Gamma in eV, an ANC in fm^-1/2);
//  (c) the details pane's numbers = the facade's; theta^2 in the channel table;
//  (d) two levels removed with Del after one confirmation; a duplicate (Ctrl+D) with the widths;
//  (e) undo of everything back to the opened text, byte for byte, and redo to the edited one;
//  (f) the search, the J^pi filter with its counts and "only included"; bulk fix of the selection;
//  (g) the missing-J^pi text names what to do and its link fills the add bar;
//  (h) the Particle Pairs tab: a cell edited in place, a pair's removal naming the data sets that use it.
// PASS/FAIL per check; exit 1 on a failure.
#include <cmath>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>

#include <QAbstractButton>
#include <QApplication>
#include <QComboBox>
#include <QDir>
#include <QDialogButtonBox>
#include <QElapsedTimer>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QTableView>
#include <QTest>
#include <QTimer>
#include <QUndoStack>

#include "AZURESetup.h"
#include "ChannelDetails.h"
#include "ChannelsModel.h"
#include "LevelsModel.h"
#include "LevelsTab.h"
#include "LevelsViews.h"
#include "PairsModel.h"
#include "PairsTab.h"
#include "ProjectHolder.h"

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
  return to / projectFile.filename();
}

void Spin(int ms = 50) {
  QElapsedTimer t;
  t.start();
  while (t.elapsed() < ms) QApplication::processEvents(QEventLoop::AllEvents, 10);
}

/// Answers the message boxes in order (the button given), counting them.
struct Answers {
  std::deque<QMessageBox::StandardButton> buttons;
  int answered = 0;
  QTimer timer;
  Answers() {
    timer.setInterval(30);
    QObject::connect(&timer, &QTimer::timeout, [this]() {
      QMessageBox *m = qobject_cast<QMessageBox *>(QApplication::activeModalWidget());
      if (!m) return;
      lastText = m->text() + "\n" + m->informativeText();
      const QMessageBox::StandardButton b = buttons.empty() ? QMessageBox::Cancel : buttons.front();
      if (!buttons.empty()) buttons.pop_front();
      if (QAbstractButton *x = m->button(b)) x->click();
      else m->reject();
      answered++;
    });
    timer.start();
  }
  QString lastText;
};

/// The rows of the project's <levels> as the engine reads them, whose column `name` is `value`.
int RowsWith(const azure2::front::Project &p, const std::string &name, const std::string &value) {
  const int c = azure2::front::LevelColumnIndex(name);
  int n = 0;
  for (const azure2::front::LevelRow &r : p.Levels()) n += (int)r.tokens.size() > c && r.tokens[c] == value;
  return n;
}

/// The model row of the level of J, parity and energy, or -1.
int LevelRow(LevelsTab *t, double j, int pi, double e) {
  const QList<LevelsData> levels = t->getLevelsModel()->getLevels();
  for (int i = 0; i < levels.size(); i++)
    if (levels[i].jValue == j && levels[i].piValue == pi && std::fabs(levels[i].energy - e) < 1e-9) return i;
  return -1;
}

/// The key that opens a cell's editor: F2, and Return on macOS (QAbstractItemView::EditKeyPressed is the platform's).
#ifdef Q_OS_MACOS
const Qt::Key kEditKey = Qt::Key_Return;
#else
const Qt::Key kEditKey = Qt::Key_F2;
#endif

/// The window active again (a message box takes the activation and offscreen does not give it back) and `w` focused.
void Focus(QWidget *w) {
  w->window()->activateWindow();
  Spin(10);
  w->setFocus();
}

/// Open the editor of a cell of a view with the edit key, type `text` and press Enter.
void EditCell(QTableView *view, const QModelIndex &index, const QString &text) {
  Focus(view);
  view->setCurrentIndex(index);
  QTest::keyClick(view, kEditKey);
  Spin(20);
  QWidget *editor = QApplication::focusWidget();
  if (!qobject_cast<QLineEdit *>(editor))
    std::cout << "  (no editor: focus on " << (editor ? editor->metaObject()->className() : "nothing") << ")" << std::endl;
  if (QLineEdit *e = qobject_cast<QLineEdit *>(editor)) {
    e->selectAll();
    QTest::keyClicks(e, text);
    QTest::keyClick(e, Qt::Key_Return);
  }
  Spin(30);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: gui_levels_test <source dir>" << std::endl;
    return 2;
  }
  QApplication app(argc, argv);
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs, true);
  QCoreApplication::setOrganizationName("jina-test");
  QCoreApplication::setApplicationName("azure2-gui-levels-test");
  const fs::path source = argv[1];
  const fs::path work = fs::absolute("gui_levels_work");
  fs::create_directories(work);

  AZURESetup w;
  w.resize(1400, 900);
  w.show();
  w.activateWindow();
  LevelsTab *levels = w.getLevelsTab();
  w.tabs()->setCurrentWidget(levels);
  Spin(100);
  QTableView *levelsView = levels->getLevelsView();
  QTableView *channelsView = levels->getChannelsView();
  LevelsModel *lm = levels->getLevelsModel();
  ChannelsModel *cm = levels->getChannelsModel();
  Answers answers;

  // (a) opened and saved without an edit; the first level selected and its channels shown
  for (const fs::path &f : {source / "tests/upgrade/T1/T1.azr", source / "tests/upgrade/T3/T3.azr", source / "tests/13N/13N.azr"}) {
    const std::string name = f.stem().string();
    const fs::path copy = CopyProject(f, work / ("rt_" + name));
    QString why;
    Check(w.readProjectFile(QString::fromStdString(copy.string()), why), name + ": opened");
    Spin(100);
    const QList<int> selected = levels->selectedLevelRows();
    const int first = levels->getLevelsProxy()->mapToSource(levels->getLevelsProxy()->index(0, 0)).row();
    Check(selected.size() == 1 && selected.first() == first && channelsView->model()->rowCount() > 0,
          name + ": the first level shown is selected and its " + std::to_string(channelsView->model()->rowCount()) +
              " channels are on screen");
    Check(levels->undoStack()->count() == 0, name + ": nothing to undo after opening");
    const fs::path out = copy.parent_path() / (name + "_saved.azr");
    Check(w.writeProjectFile(QString::fromStdString(out.string()), why) && Bytes(out) == Bytes(f),
          name + ": saved without an edit, byte for byte");
  }

  // (b) T1: add, edit in place, refuse
  const fs::path t1 = source / "tests/upgrade/T1/T1.azr";
  const fs::path copy = CopyProject(t1, work / "edit_T1");
  QString why;
  w.readProjectFile(QString::fromStdString(copy.string()), why);
  w.tabs()->setCurrentWidget(levels);  // an opened project shows the Overview first (qt_gui.md section 15)
  Spin(100);
  const std::string original = w.holder()->project().Text();
  const int numLevels = lm->getLevels().size();
  {
    levels->addJSpin()->setValue(2.5);
    if (levels->addParityToggle()->parity() != -1) levels->addParityToggle()->click();
    QLineEdit *e = levels->addEnergyEdit();
    e->setFocus();
    e->clear();
    QTest::keyClicks(e, "7.5");
    QTest::keyClick(e, Qt::Key_Return);
    Spin(50);
    const int row = LevelRow(levels, 2.5, -1, 7.5);
    int channels = 0;
    for (const ChannelsData &c : cm->getChannels()) channels += c.levelIndex == row;
    Check(row >= 0 && lm->getLevels().size() == numLevels + 1 && channels > 0 && levels->currentLevelRow() == row &&
              e->text().isEmpty(),
          "add bar: 5/2- at 7.5 MeV added with " + std::to_string(channels) + " channels, selected, the field cleared");
    QTest::keyClicks(e, "abc");
    QTest::keyClick(e, Qt::Key_Return);
    Spin(20);
    Check(lm->getLevels().size() == numLevels + 1 && levels->statusLabel()->text().contains("not a number"),
          "add bar: 'abc' refused, the reason shown (" + levels->statusLabel()->text().toStdString() + ")");
    e->clear();
    QTest::keyClicks(e, "7.5005");
    QTest::keyClick(e, Qt::Key_Return);
    Spin(20);
    Check(lm->getLevels().size() == numLevels + 1 && levels->statusLabel()->text().contains("1 keV"),
          "add bar: a 5/2- level 0.5 keV from another refused (the engine's 1 keV)");
    e->clear();
  }
  {
    // the energy in place: F2, type, Enter; then an invalid one
    const int row = LevelRow(levels, 2.5, -1, 7.5);
    const QModelIndex at = levels->getLevelsProxy()->mapFromSource(lm->index(row, 4));
    EditCell(levelsView, at, "7.75");
    Check(LevelRow(levels, 2.5, -1, 7.75) == row, "in place: the new level's energy typed as 7.75 MeV");
    EditCell(levelsView, levels->getLevelsProxy()->mapFromSource(lm->index(row, 4)), "seven");
    Check(LevelRow(levels, 2.5, -1, 7.75) == row && levels->statusLabel()->text().contains("not a number"),
          "in place: 'seven' refused, the energy unchanged");
    // J^pi in place: the editor's spin box and parity toggle, Enter
    levelsView->setFocus();
    levelsView->setCurrentIndex(levels->getLevelsProxy()->mapFromSource(lm->index(row, 2)));
    QTest::keyClick(levelsView, kEditKey);
    Spin(20);
    HalfSpinBox *j = levelsView->viewport()->findChild<HalfSpinBox *>();
    ParityToggle *pi = levelsView->viewport()->findChild<ParityToggle *>();
    if (j && pi) {
      j->setValue(1.5);
      pi->click();  // - -> +
      QTest::keyClick(j, Qt::Key_Return);
    }
    Spin(30);
    int channels = 0;
    for (const ChannelsData &c : cm->getChannels()) channels += c.levelIndex == row;
    Check(j && pi && LevelRow(levels, 1.5, 1, 7.75) == row && channels > 0,
          "in place: J^pi 5/2- -> 3/2+ with the J spin box and the parity toggle; " + std::to_string(channels) +
              " channels of 3/2+");
  }
  {
    // a width in place: the 1/2+ level at 2.365 MeV, p + 12C, l = 0
    const int row = LevelRow(levels, 0.5, 1, 2.365);
    levels->selectLevelRows({row});
    Spin(30);
    int ch = -1;
    const QList<ChannelsData> channels = cm->getChannels();
    for (int i = 0; i < channels.size(); i++)
      if (channels[i].levelIndex == row && channels[i].pairIndex == 0 && channels[i].lValue == 0) ch = i;
    const QSortFilterProxyModel *cp = qobject_cast<const QSortFilterProxyModel *>(channelsView->model());
    const QModelIndex at = cp->mapFromSource(cm->index(ch, ChannelsModel::ValueColumn));
    Check(cm->data(cm->index(ch, ChannelsModel::MeaningColumn), Qt::DisplayRole).toString() == QString::fromUtf8("Γ [eV]"),
          "the value of an open particle channel is Γ [eV]");
    EditCell(channelsView, at, "34000");
    Check(cm->getChannels()[ch].reducedWidth == 34000, "in place: the width typed as 34000 eV");
    EditCell(channelsView, cp->mapFromSource(cm->index(ch, ChannelsModel::ValueColumn)), "1e");
    Check(cm->getChannels()[ch].reducedWidth == 34000, "in place: '1e' refused, the width unchanged");
    // (c) the details pane: P of the facade at the level's energy, theta^2 in the table
    channelsView->setCurrentIndex(cp->mapFromSource(cm->index(ch, 3)));
    channelsView->selectRow(cp->mapFromSource(cm->index(ch, 3)).row());
    Spin(30);
    const azure2::front::ChannelQuantities q = levels->channelQuantities(ch);
    const QString details = levels->findChild<ChannelDetails *>()->details->text();
    Check(details.contains(QString::number(q.P, 'g', 6)) && std::fabs(q.P - 0.0112659697575860) < 1e-9,
          "details: the penetrability of the facade, 0.0112660 at 2.365 MeV (" + std::to_string(q.P) + ")");
    Check(cm->data(cm->index(ch, ChannelsModel::RatioColumn), Qt::DisplayRole).toString() ==
              QString::number(q.wignerRatio, 'g', 3),
          "the table's theta^2 = the facade's (" + std::to_string(q.wignerRatio) + ")");
    // a closed channel: the 1/2- level at 0 MeV
    const int bound = LevelRow(levels, 0.5, -1, 0.0);
    levels->selectLevelRows({bound});
    Spin(30);
    QString meaning;
    const QList<ChannelsData> all = cm->getChannels();
    for (int i = 0; i < all.size(); i++)
      if (all[i].levelIndex == bound && all[i].pairIndex == 0) meaning = cm->data(cm->index(i, ChannelsModel::MeaningColumn), Qt::DisplayRole).toString();
    Check(meaning == "C [fm^-1/2]", "the value of a closed channel is an ANC, C [fm^-1/2] (" + meaning.toStdString() + ")");
  }
  {
    // the rows the tab writes (LevelsTab::levelRows, what a save hands to the engine's row writer); the project is
    // saved only after the undo of (e), so that the undo is compared with the text as opened
    azure2::front::Project rows = azure2::front::Project::FromText(original, copy.string());
    std::string reason;
    Check(rows.SetLevels(levels->levelRows(), reason) && RowsWith(rows, "levelE", "7.75") > 0 &&
              RowsWith(rows, "gamma", "34000") == 1 && rows.Validate().ok,
          "T1's rows after the edits: the 3/2+ level at 7.75 MeV and the width 34000 eV, read by the engine, valid");
  }

  // (d) remove two levels with Del after one confirmation; duplicate with Ctrl+D
  {
    const int n = lm->getLevels().size();
    const int a = LevelRow(levels, 1.5, 1, 7.75), b = LevelRow(levels, 0.5, 1, 2.365);
    levels->selectLevelRows({a, b});
    levelsView->setFocus();
    answers.buttons.push_back(QMessageBox::Yes);
    const int asked = answers.answered;
    QTest::keyClick(levelsView, Qt::Key_Delete);
    Spin(100);
    Check(lm->getLevels().size() == n - 2 && answers.answered == asked + 1 && answers.lastText.contains("2.365") &&
              answers.lastText.contains("7.75"),
          "Del: two levels removed after one confirmation naming both");
    const int src = LevelRow(levels, 1.5, -1, 1.0);
    levels->selectLevelRows({src});
    Focus(levelsView);
    QTest::keyClick(levelsView, Qt::Key_D, Qt::ControlModifier);
    Spin(50);
    QWidget *editor = QApplication::focusWidget();
    if (qobject_cast<QLineEdit *>(editor)) QTest::keyClick(editor, Qt::Key_Escape);
    const int copyRow = LevelRow(levels, 1.5, -1, 1.01);
    bool sameWidths = copyRow >= 0;
    for (const ChannelsData &c : cm->getChannels()) {
      if (c.levelIndex != copyRow) continue;
      for (const ChannelsData &o : cm->getChannels())
        if (o.levelIndex == src && o.pairIndex == c.pairIndex && o.lValue == c.lValue && o.sValue == c.sValue)
          sameWidths = sameWidths && o.reducedWidth == c.reducedWidth && o.isFixed == c.isFixed;
    }
    Check(copyRow >= 0 && sameWidths && levels->currentLevelRow() == copyRow,
          "Ctrl+D: 3/2- duplicated 10 keV above (1.01 MeV) with its widths, the copy selected");
  }

  // (e) undo to the opened text, redo to the edited one
  {
    Spin(50);
    QUndoStack *undo = levels->undoStack();
    const int steps = undo->count();
    QStringList before;
    for (const LevelsData &l : lm->getLevels()) before << QString("%1 %2 %3").arg(l.jValue).arg(l.piValue).arg(l.energy);
    Focus(levelsView);
    for (int i = 0; i < steps; i++) QTest::keyClick(levelsView, Qt::Key_Z, Qt::ControlModifier);
    Spin(50);
    const QString back = QString::fromStdString((work / "edit_T1" / "T1_undone.azr").string());
    Check(steps >= 6 && !undo->canUndo() && w.writeProjectFile(back, why) && Bytes(back.toStdString()) == Bytes(t1),
          "Ctrl+Z x " + std::to_string(steps) + " (add, energy, J^pi, width, remove, duplicate): the project as opened, "
          "byte for byte");
    for (int i = 0; i < steps; i++) QTest::keyClick(levelsView, Qt::Key_Z, Qt::ControlModifier | Qt::ShiftModifier);
    Spin(50);
    QStringList after;
    for (const LevelsData &l : lm->getLevels()) after << QString("%1 %2 %3").arg(l.jValue).arg(l.piValue).arg(l.energy);
    const QString again = QString::fromStdString((work / "edit_T1" / "T1_redone.azr").string());
    azure2::front::Project redone;
    std::string reason;
    Check(!undo->canRedo() && after == before && w.writeProjectFile(again, why) &&
              azure2::front::Project::Open(again.toStdString(), redone, reason) && RowsWith(redone, "levelE", "1.01") > 0 &&
              redone.Validate().ok,
          "Ctrl+Shift+Z x " + std::to_string(steps) + ": every edit again; saved with the duplicate at 1.01 MeV, valid");
  }

  // (f) search, J^pi filter, only included, bulk fix
  {
    LevelsFilterProxy *proxy = levels->getLevelsProxy();
    int threeHalves = 0;
    for (const LevelsData &l : lm->getLevels()) threeHalves += l.jValue == 1.5 && l.piValue == -1;
    QTest::keyClicks(levels->searchEdit(), "3/2-");
    Spin(20);
    Check(proxy->rowCount() == threeHalves && threeHalves >= 2, "search '3/2-': the " + std::to_string(threeHalves) + " 3/2- levels");
    levels->searchEdit()->clear();
    QComboBox *jpi = levels->jpiFilterCombo();
    const int k = jpi->findText(QString("%1 (%2)").arg(AddLevelDialog::spinParityLabel(1.5, -1)).arg(threeHalves));
    jpi->setCurrentIndex(k);
    Spin(20);
    Check(k > 0 && proxy->rowCount() == threeHalves, "J^pi filter: '3/2- (" + std::to_string(threeHalves) + ")' shows them");
    jpi->setCurrentIndex(0);
    const int all = proxy->rowCount();
    const int out = LevelRow(levels, 1.5, -1, 1.01);
    lm->setData(lm->index(out, 0), Qt::Unchecked, Qt::CheckStateRole);
    levels->onlyIncludedCheck()->setChecked(true);
    Spin(20);
    Check(proxy->rowCount() == all - 1, "only included: the level left out is hidden");
    levels->onlyIncludedCheck()->setChecked(false);
    QList<int> rows;
    for (int i = 0; i < lm->getLevels().size(); i++) rows << i;
    levels->setFixed(rows, false, true);
    bool fixed = true;
    for (const ChannelsData &c : cm->getChannels()) fixed = fixed && c.isFixed == 1;
    levels->setFixed(rows.mid(0, 1), true, false);
    Check(fixed && lm->getLevels()[0].isFixed == 0, "bulk: every width fixed, the first level's energy freed");
  }

  // (g) the missing J^pi: what to do, and the link fills the add bar
  {
    levels->setIncludeAllJPi(false);
    Spin(20);
    // the label is one line (qt_gui.md section 19.6); what to do is its tool tip
    const QString text = levels->missingLabel()->text();
    Check(levels->missingLabel()->toolTip().contains("switch the option on") && text.contains("href='jpi:"),
          "missing J^pi: says what to do (its tip), each J^pi a link");
    const int at = text.indexOf("href='jpi:");
    const QString link = text.mid(at + 6, text.indexOf('\'', at + 6) - at - 6);
    QMetaObject::invokeMethod(levels->missingLabel(), "linkActivated", Q_ARG(QString, link));
    Spin(20);
    const QStringList p = link.split(':');
    Check(p.size() >= 3 && (int)std::lround(2 * levels->addJSpin()->value()) == p[1].toInt() &&
              levels->addParityToggle()->parity() == p[2].toInt(),
          "missing J^pi: the link " + link.toStdString() + " fills the add bar");
  }

  // (h) the Particle Pairs tab
  {
    PairsTab *pairs = w.getPairsTab();
    PairsModel *pm = pairs->getPairsModel();
    w.tabs()->setCurrentWidget(pairs);
    Spin(30);
    const std::vector<std::string> uses = azure2::front::PairUses(w.holder()->project(), 1);
    std::string list;
    for (const std::string &u : uses) list += " " + u;
    Check(!uses.empty(), "pair 1 is named by the data sets:" + list);
    const QStringList goes = pairs->pairConsequences(0);
    Check(goes.join(" ").contains("channel(s) of") && goes.join(" ").contains("data: "),
          "removing pair 1 names its channels and data sets: " + goes.join(" | ").toStdString());
    // the radius of pair 1 in place: F2, type, Enter; a negative one refused
    const double radius = pm->getPairs()[0].channelRadius;
    EditCell(pairs->getPairsView(), pm->index(0, 12), "4.1");
    Check(pm->getPairs()[0].channelRadius == 4.1, "pairs in place: pair 1's radius typed as 4.1 fm (was " +
                                                      std::to_string(radius) + ")");
    EditCell(pairs->getPairsView(), pm->index(0, 12), "-1");
    Check(pm->getPairs()[0].channelRadius == 4.1, "pairs in place: a radius of -1 fm refused");
    levels->recordUndo();
    Focus(pairs->getPairsView());
    QTest::keyClick(pairs->getPairsView(), Qt::Key_Z, Qt::ControlModifier);
    Spin(30);
    Check(pm->getPairs()[0].channelRadius == radius, "pairs: Ctrl+Z gives the radius back");
  }

  // (i) screenshots (a second argument: the directory) of T1 as opened, a channel of the 1/2+ level selected
  if (argc >= 3) {
    const QString dir = QString::fromLocal8Bit(argv[2]);
    QDir().mkpath(dir);
    const fs::path shot = CopyProject(t1, work / "shot_T1");
    w.readProjectFile(QString::fromStdString(shot.string()), why);
    Spin(100);
    w.tabs()->setCurrentWidget(levels);
    levels->selectLevelRows({LevelRow(levels, 0.5, 1, 2.365)});
    Spin(50);
    channelsView->selectRow(0);
    levels->setIncludeAllJPi(false);
    Spin(100);
    w.grab().save(dir + "/tab1_levels_and_channels.png");
    levels->jpiFilterCombo()->setCurrentIndex(1);
    levels->searchEdit()->setText("");
    Spin(50);
    w.grab().save(dir + "/tab1_levels_filtered.png");
    levels->jpiFilterCombo()->setCurrentIndex(0);
    w.tabs()->setCurrentWidget(w.getPairsTab());
    Spin(100);
    w.grab().save(dir + "/tab0_particle_pairs.png");
    std::cout << "screenshots in " << dir.toStdString() << std::endl;
  }

  std::cout << (failures ? "FAILED: " + std::to_string(failures) : std::string("all passed")) << std::endl;
  return failures ? 1 : 0;
}
