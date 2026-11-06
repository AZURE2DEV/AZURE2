// Tied systematics in the GUI: the "normgroup <label>", "shiftgroup <label>" and "shiftfngroup <label>" blocks of
// a data-segment line (SegLine.h; EData::ResolveNormGroups / ResolveShiftGroups / ResolveShiftFunctions).  Every
// active segment with one label shares the normalization / constant shift / shift function of the group's first
// active segment in file order.
//
//   1. SegmentKeywords: the labels read from a tail (block by block, so a value spelled like a keyword is not taken
//      for one), and a label replaced, removed or added with every other token kept;
//   2. a project whose lines carry the blocks among sqrtshift, covariance, shiftfn and sample blocks: the Segments
//      tab's Groups column and tooltips, the group leaders (inactive lines are not members);
//   3. the data-segment dialog shows each line's labels in its three drop-downs (listing the project's labels),
//      refuses a label with a blank, and an edit that changes, clears or adds a label rewrites that block only;
//      an edit that changes nothing leaves the tail byte for byte; a line added from the dialog gets its blocks;
//   4. the Fitting tab shows a group member's normalization / energy shift as tied ("tied to segment 1 (normgroup
//      D)"), greyed, without a check box and with the leader's value, and follows the Segments tab;
//   5. the project written back by the GUI's writer carries the edited lines.
//
// Runs without a display (QT_QPA_PLATFORM=offscreen).  With AZURE2_SEG_GROUPS_PNG=<stem> it saves the edit dialog,
// the Segments tab and the Fitting tab's Normalization table as <stem>_dialog.png, <stem>_segments.png and
// <stem>_fitting.png (use xvfb-run for real fonts).
#include <QApplication>
#include <QComboBox>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QGroupBox>
#include <QLineEdit>
#include <QTableView>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTextStream>
#include <QTimer>
#include <cstdlib>
#include <functional>
#include <iostream>
#include "AZURESetup.h"
#include "AddSegDataDialog.h"
#include "Config.h"
#include "FittingTab.h"
#include "SegmentKeywords.h"
#include "SegmentsDataModel.h"
#include "SegmentsTab.h"
struct SegPairs {int firstPair; int secondPair;};
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

namespace {
int failures = 0;
void ok(const std::string& what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if (!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if (!cond) failures++;
}

const char* pngStem() { return std::getenv("AZURE2_SEG_GROUPS_PNG"); }

// Runs \p onDialog on the data-segment dialog the next call opens, then accepts it (a safety net rejects any
// dialog still open after 10 s, so the test never hangs).
void whenDialogOpens(const std::function<void(AddSegDataDialog*)>& onDialog) {
  QTimer::singleShot(0, [onDialog]() {
    for (QWidget* w : QApplication::topLevelWidgets())
      if (AddSegDataDialog* d = qobject_cast<AddSegDataDialog*>(w))
        if (d->isVisible()) {
          onDialog(d);
          d->accept();
          return;
        }
  });
  QTimer::singleShot(10000, []() {
    for (QWidget* w : QApplication::topLevelWidgets())
      if (QDialog* d = qobject_cast<QDialog*>(w))
        if (d->isVisible()) d->reject();
  });
}

void editRow(SegmentsTab* tab, int row, const std::function<void(AddSegDataDialog*)>& onDialog) {
  QTableView* view = tab->findChild<QTableView*>("segmentsDataView");
  view->clearSelection();
  view->selectRow(row);
  whenDialogOpens(onDialog);
  tab->editSegDataLine();
}

const FittingParameter* param(FittingTab* tab, const QString& name) {
  for (const FittingParameter& p : tab->getFittingParameters())
    if (p.name == name) return &p;
  return nullptr;
}

int tableRow(QTableWidget* table, const QString& name) {
  for (int r = 0; r < table->rowCount(); r++)
    if (table->item(r, 0) && table->item(r, 0)->text() == name) return r;
  return -1;
}

// The <segmentsData> lines of a project file, whitespace-normalised.
QStringList segmentLines(const QString& path) {
  QFile f(path);
  QStringList out;
  if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) return out;
  bool in = false;
  for (const QString& line : QString::fromUtf8(f.readAll()).split('\n')) {
    if (line.trimmed() == "<segmentsData>") { in = true; continue; }
    if (line.trimmed() == "</segmentsData>") break;
    if (in && !line.trimmed().isEmpty()) out << line.simplified();
  }
  return out;
}
}  // namespace

int main(int argc, char** argv) {
  QApplication app(argc, argv);

  // ---- 1. the tail helpers ---------------------------------------------------------------------------------------
  {
    using namespace SegmentKeywords;
    const QString t = "sqrtshift 0 0 0 normgroup D shiftgroup A";
    ok("labels read from a tail", groupLabel(t, "normgroup") == "D" && groupLabel(t, "shiftgroup") == "A" &&
                                      groupLabel(t, "shiftfngroup").isEmpty());
    ok("a covariance file named 'normgroup' is not a normgroup block",
       groupLabel("covariance normgroup normgroup D", "normgroup") == "D" &&
           groupLabel("covariance normgroup", "normgroup").isEmpty());
    ok("a shiftfn block is skipped whole (its quoted expression may hold blanks)",
       groupLabel("shiftfn \"a0 * x\" 1 0 0.02 1 shiftfngroup F", "shiftfngroup") == "F");
    ok("change a label in place", withGroupLabel(t, "normgroup", "E") == "sqrtshift 0 0 0 normgroup E shiftgroup A",
       withGroupLabel(t, "normgroup", "E"));
    ok("clear a label", withGroupLabel(t, "shiftgroup", "") == "sqrtshift 0 0 0 normgroup D",
       withGroupLabel(t, "shiftgroup", ""));
    ok("add a label at the end", withGroupLabel("sqrtshift 0 0 0", "normgroup", "D") == "sqrtshift 0 0 0 normgroup D");
    ok("add a label before a sample block",
       withGroupLabel("covariance c.txt sample 1", "shiftgroup", "A") == "covariance c.txt shiftgroup A sample 1",
       withGroupLabel("covariance c.txt sample 1", "shiftgroup", "A"));
    ok("add to an empty tail", withGroupLabel("", "normgroup", "D") == "normgroup D");
    const QString fn = "shiftfn \"a0 * (sqrt(x)-sqrt(1.7))\" 1 0 0.02 1 shiftfngroup F";
    ok("a quoted expression with blanks survives an edit",
       withGroupLabel(fn, "normgroup", "D") == fn + " normgroup D", withGroupLabel(fn, "normgroup", "D"));
    ok("labels with a blank or a quote are refused",
       !isValidGroupLabel("a b") && !isValidGroupLabel("a\"") && isValidGroupLabel("Dayras_1976"));
  }

  // ---- 2. a project with groups ----------------------------------------------------------------------------------
  QTemporaryDir tmp;
  const QString source = QString(AZURE2_TEST_SOURCE_DIR) + "/tests/norm_groups/norm_groups.azr";
  const QString project = tmp.filePath("groups.azr");
  const QStringList lines = {
      "1 1 1 0 100 40 40 4 1.2 1 10 0.004 0.01 1 data/a.dat 0 0 sqrtshift 0 0 0 normgroup D shiftgroup A",
      "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/b.dat 0 0 normgroup D shiftgroup A",
      "1 1 1 0 100 60 60 4 1 0 1 0 0 0 data/c.dat 0 0 covariance cov/c.txt shiftfngroup F shiftfn "
      "\"a0*(sqrt(x)-sqrt(1.7))\" 1 0 0.02 1 sample 1",
      "1 1 1 0 100 120 120 4 1 0 0 0 0 0 data/d.dat 0 0 shiftfngroup F",
      "0 1 1 0 100 150 150 4 1 0 0 0 0 0 data/e.dat 0 0 normgroup D"};
  {
    QFile in(source);
    ok("read the norm_groups project", in.open(QIODevice::ReadOnly | QIODevice::Text));
    QString text = QString::fromUtf8(in.readAll());
    const int a = text.indexOf("<segmentsData>\n") + QString("<segmentsData>\n").size();
    const int b = text.indexOf("</segmentsData>");
    text.replace(a, b - a, lines.join("\n") + "\n");
    QFile out(project);
    out.open(QIODevice::WriteOnly | QIODevice::Text);
    out.write(text.toUtf8());
  }
  AZURESetup setup;
  ok("open the project", setup.readProjectFile(project));
  SegmentsTab* tab = setup.getSegmentsTab();
  SegmentsDataModel* model = tab->getSegmentsDataModel();
  ok("five data lines", model->rowCount(QModelIndex()) == 5);
  ok("the model counts the Groups column", model->columnCount(QModelIndex()) == SegmentsDataModel::GroupsColumn + 1);
  auto cell = [&](int row, int role) { return model->data(model->index(row, SegmentsDataModel::GroupsColumn), role).toString(); };
  ok("Groups column of line 1", cell(0, Qt::DisplayRole) == "norm: D, shift: A", cell(0, Qt::DisplayRole));
  ok("Groups column of line 3", cell(2, Qt::DisplayRole) == "shift fn: F", cell(2, Qt::DisplayRole));
  ok("Groups column header", model->headerData(SegmentsDataModel::GroupsColumn, Qt::Horizontal, Qt::DisplayRole) == "Groups");
  ok("line 1 leads, line 2 follows (normgroup and shiftgroup)",
     model->groupLeader(0, "normgroup") == -1 && model->groupLeader(1, "normgroup") == 0 &&
         model->groupLeader(1, "shiftgroup") == 0 && model->groupLeader(3, "shiftfngroup") == 2);
  ok("an inactive line is no member", model->groupLeader(4, "normgroup") == -1);
  ok("tooltip of the leader names its followers", cell(0, Qt::ToolTipRole).contains("used by segment 2"), cell(0, Qt::ToolTipRole));
  ok("tooltip of a follower names its leader", cell(1, Qt::ToolTipRole).contains("of segment 1"), cell(1, Qt::ToolTipRole));
  ok("a follower's own norm is marked unused",
     model->data(model->index(1, 9), Qt::ToolTipRole).toString().contains("Not used"));
  ok("labels offered for normgroup", model->groupLabels("normgroup") == QStringList{"D"});

  // ---- 3. the dialog ---------------------------------------------------------------------------------------------
  bool shown = false;
  editRow(tab, 0, [&](AddSegDataDialog* d) {
    shown = d->normGroupCombo->currentText() == "D" && d->shiftGroupCombo->currentText() == "A" &&
            d->shiftFnGroupCombo->currentText().isEmpty() && d->groupsBox->isVisible();
    ok("the drop-downs list the project's labels",
       d->normGroupCombo->findText("D") > 0 && d->shiftGroupCombo->findText("A") > 0 &&
           d->shiftFnGroupCombo->findText("F") > 0 && d->normGroupCombo->itemText(0).isEmpty());
    ok("shift-function group disabled while an energy-shift group is set", !d->shiftFnGroupCombo->isEnabled());
    int pos = 0;
    QString blank = "a b";
    ok("a label with a blank is refused as typed",
       d->normGroupCombo->validator() &&
           d->normGroupCombo->validator()->validate(blank, pos) == QValidator::Invalid);
    if (const char* stem = pngStem()) {
      d->resize(d->sizeHint());
      d->grab().save(QString("%1_dialog.png").arg(stem));
    }
    d->setGroupLabel("shiftgroup", "");  // clear
    ok("clearing the energy-shift group enables the shift-function group", d->shiftFnGroupCombo->isEnabled());
  });
  ok("edit dialog showed line 1's labels", shown);
  ok("shiftgroup cleared, every other token kept", model->getLines()[0].extraTokens == "sqrtshift 0 0 0 normgroup D",
     model->getLines()[0].extraTokens);

  editRow(tab, 1, [](AddSegDataDialog* d) { d->normGroupCombo->setEditText("E"); });
  ok("normgroup changed in place", model->getLines()[1].extraTokens == "normgroup E shiftgroup A",
     model->getLines()[1].extraTokens);

  editRow(tab, 3, [](AddSegDataDialog* d) { d->setGroupLabel("normgroup", "D"); });
  ok("normgroup added at the end", model->getLines()[3].extraTokens == "shiftfngroup F normgroup D",
     model->getLines()[3].extraTokens);

  const QString before = model->getLines()[2].extraTokens;
  editRow(tab, 2, [](AddSegDataDialog* d) {
    ok("line 3: energy-shift group disabled while a shift-function group is set",
       d->shiftFnGroupCombo->currentText() == "F" && !d->shiftGroupCombo->isEnabled());
  });
  ok("an edit that changes nothing leaves the tail byte for byte", model->getLines()[2].extraTokens == before,
     model->getLines()[2].extraTokens);

  whenDialogOpens([](AddSegDataDialog* d) {
    d->lowEnergyText->setText("0");
    d->highEnergyText->setText("100");
    d->dataFileText->setText("data/f.dat");
    d->setGroupLabel("normgroup", "D");
    d->shiftGroupCombo->setEditText("B");
  });
  tab->addSegDataLine();
  ok("a line added from the dialog", model->rowCount(QModelIndex()) == 6);
  if (model->rowCount(QModelIndex()) == 6)
    ok("... carries its group blocks", model->getLines()[5].extraTokens == "normgroup D shiftgroup B",
       model->getLines()[5].extraTokens);
  ok("Groups column follows the edits", cell(1, Qt::DisplayRole) == "norm: E, shift: A" &&
                                            cell(3, Qt::DisplayRole) == "norm: D, shift fn: F",
     cell(1, Qt::DisplayRole) + " | " + cell(3, Qt::DisplayRole));

  // ---- 4. the Fitting tab ----------------------------------------------------------------------------------------
  FittingTab* fitting = setup.getFittingTab();
  QMetaObject::invokeMethod(fitting, "refreshParameters");
  QTableWidget* norms = fitting->findChild<QTableWidget*>("normParamsTable");
  QTableWidget* shifts = fitting->findChild<QTableWidget*>("shiftParamsTable");
  const FittingParameter* n1 = param(fitting, "segment_1_norm");
  const FittingParameter* n4 = param(fitting, "segment_4_norm");
  const FittingParameter* n2 = param(fitting, "segment_2_norm");
  ok("leader's norm is a parameter of its own", n1 && n1->tiedTo.isEmpty());
  ok("a follower's norm is tied to the leader, with its value",
     n4 && n4->tiedTo == "tied to segment 1 (normgroup D)" && n4->value == 1.2 && n4->tiedLeader == 0,
     n4 ? n4->tiedTo + QString(" %1").arg(n4->value) : QString("missing"));
  ok("a lone label ties nothing (normgroup E has one active line)", n2 && n2->tiedTo.isEmpty());
  ok("the follower's own norm in the Segments tab is untouched", model->getLines()[3].dataNorm == 1.0);
  const int r4 = norms ? tableRow(norms, "segment_4_norm") : -1;
  ok("the follower's row is in the Normalization table", r4 >= 0);
  if (r4 >= 0) {
    QTableWidgetItem* nuis = norms->item(r4, 6);
    ok("its last cell says what it is tied to", nuis && nuis->text() == "tied to segment 1 (normgroup D)");
    ok("no check box: it cannot be varied on its own", nuis && !nuis->data(Qt::CheckStateRole).isValid() &&
                                                          !(nuis->flags() & Qt::ItemIsUserCheckable));
    bool greyed = true;
    for (int c = 0; c < norms->columnCount(); c++)
      if (norms->item(r4, c) && (norms->item(r4, c)->flags() & (Qt::ItemIsEnabled | Qt::ItemIsEditable))) greyed = false;
    ok("greyed and read-only in every column", greyed);
    ok("the value shown is the leader's", norms->item(r4, 1)->text() == "1.2", norms->item(r4, 1)->text());
    ok("tooltip explains it", norms->item(r4, 1)->toolTip().contains("tied to segment 1"));
  }
  if (const char* stem = pngStem()) {
    QTabWidget* tabs = fitting->findChild<QTabWidget*>();
    if (tabs) tabs->setCurrentWidget(norms);
    fitting->resize(1100, 420);
    fitting->grab().save(QString("%1_fitting.png").arg(stem));
    tab->resize(1400, 520);
    tab->grab().save(QString("%1_segments.png").arg(stem));
  }

  // the leader's norm changes in the Segments tab: the tied row follows
  model->setData(model->index(0, 9), 1.3, Qt::EditRole);
  n4 = param(fitting, "segment_4_norm");
  ok("the tied row follows the leader's value", n4 && n4->value == 1.3 && r4 >= 0 && norms->item(r4, 1)->text() == "1.3",
     n4 ? QString::number(n4->value) : QString());
  // the follower leaves the group: a parameter of its own again, with its own value and a check box
  model->setLineExtras(3, "shiftfngroup F");
  n4 = param(fitting, "segment_4_norm");
  const int r4b = norms ? tableRow(norms, "segment_4_norm") : -1;
  ok("a line that leaves its group is independent again", n4 && n4->tiedTo.isEmpty() && n4->value == 1.0 && r4b >= 0 &&
                                                              norms->item(r4b, 6)->data(Qt::CheckStateRole).isValid());
  // energy shifts: line 2 follows line 1 in shift group A once line 1 is in it again
  model->setLineExtras(0, "sqrtshift 0 0 0 normgroup D shiftgroup A");
  const FittingParameter* s2 = param(fitting, "segment_2_energy_shift");
  ok("a shift-group follower's energy shift is tied",
     s2 && s2->tiedTo == "tied to segment 1 (shiftgroup A)" && s2->value == 0.004,
     s2 ? s2->tiedTo + QString(" %1").arg(s2->value) : QString("missing"));
  const int rs2 = shifts ? tableRow(shifts, "segment_2_energy_shift") : -1;
  ok("... and shown so", rs2 >= 0 && shifts->item(rs2, 6)->text() == "tied to segment 1 (shiftgroup A)");
  // deactivating the leader hands the group to the next active line
  model->setData(model->index(0, 0), Qt::Unchecked, Qt::CheckStateRole);
  s2 = param(fitting, "segment_2_energy_shift");
  ok("an inactive leader hands the group on (line 2 leads shift group A)", s2 && s2->tiedTo.isEmpty());
  model->setData(model->index(0, 0), Qt::Checked, Qt::CheckStateRole);

  // ---- 5. written back -------------------------------------------------------------------------------------------
  const QString saved = tmp.filePath("saved.azr");
  ok("project written", setup.writeProjectFile(saved));
  QStringList want = lines;
  for (QString& l : want) l = l.simplified();
  want[0] = "1 1 1 0 100 40 40 4 1.3 1 10 0.004 0.01 1 data/a.dat 0 0 sqrtshift 0 0 0 normgroup D shiftgroup A";
  want[1] = "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/b.dat 0 0 normgroup E shiftgroup A";
  want << "1 1 1 0 100 0 180 0 1 0 0 0 0 0 data/f.dat 0 0 normgroup D shiftgroup B";
  const QStringList got = segmentLines(saved);
  ok("six lines written", got.size() == 6, QString::number(got.size()));
  for (int i = 0; i < std::min(got.size(), want.size()); i++)
    ok(QString("line %1 written as expected").arg(i + 1).toStdString(), got[i] == want[i], got[i]);
  // and read back by the GUI's own reader, unchanged
  AZURESetup again;
  ok("written project reopens", again.readProjectFile(saved));
  QStringList tails;
  for (const SegmentsDataData& l : again.getSegmentsTab()->getSegmentsDataModel()->getLines()) tails << l.extraTokens;
  QStringList tailsBefore;
  for (const SegmentsDataData& l : model->getLines()) tailsBefore << l.extraTokens;
  ok("tails read back identical", tails == tailsBefore, tails.join(" | "));

  // the engine's own group test projects, opened and saved unedited: segment lines unchanged
  for (const QString& name : {QString("norm_groups/norm_groups.azr"), QString("shift_groups/shift_groups.azr")}) {
    const QString file = QString(AZURE2_TEST_SOURCE_DIR) + "/tests/" + name;
    AZURESetup s;
    const QString out = tmp.filePath(QFileInfo(name).fileName());
    ok(QString("%1 opened and saved").arg(name).toStdString(), s.readProjectFile(file) && s.writeProjectFile(out));
    const QStringList a = segmentLines(file), b = segmentLines(out);
    ok(QString("%1: segment lines reproduced").arg(name).toStdString(), !a.isEmpty() && a == b, b.join(" | "));
  }

  if (failures) {
    std::cout << "FAIL: " << failures << " check(s)" << std::endl;
    return 1;
  }
  std::cout << "all checks passed" << std::endl;
  return 0;
}
