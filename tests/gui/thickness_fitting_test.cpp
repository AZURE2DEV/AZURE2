// The fitted sample thickness in the GUI: the transmission / self-shielded-yield block of a segment line,
// "transmission n [vary sigma]" / "yield n [vary sigma]", and "thicknessgroup <label>" (SegLine.h,
// EData::ResolveThicknessGroups).
//
//   1. SegmentKeywords: the block read from a tail (plain and fit form; a group label after the fit form is still
//      found), and the block set, changed or removed with every other token kept;
//   2. the data-segment dialog's Sample row shows a line's block and group; an edit that changes nothing leaves the
//      tail byte for byte, a change rewrites only the block (or the group), a line without one gets one; the
//      Segments table shows the observable and the thickness, the Groups column the thickness group;
//   3. the test-segment dialog's Sample row (kind and n, no fit parameter) keeps a test line's fit form as read;
//   4. the Fitting tab's Sample Thicknesses rows: one segment_<key>_thickness per active line with the fit form
//      (Value n, Error the prior sigma, Nuisance = vary, always a Minuit index), a thicknessgroup follower tied to
//      its leader with the leader's value; edits in the table go back to the line; a .sav value comes back as n;
//   5. the project written by the GUI carries the edited lines.
//
// Runs without a display (QT_QPA_PLATFORM=offscreen).
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QFile>
#include <QLineEdit>
#include <QTableView>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTextStream>
#include <QTimer>
#include <functional>
#include <iostream>
#include "AZURESetup.h"
#include "AddSegDataDialog.h"
#include "AddSegTestDialog.h"
#include "Config.h"
#include "FittingTab.h"
#include "SegmentKeywords.h"
#include "SegmentsDataModel.h"
#include "SegmentsTab.h"
#include "SegmentsTestModel.h"
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

// Runs \p onDialog on the dialog of type D the next call opens, then accepts it (a safety net rejects any dialog
// still open after 10 s).
template <class D>
void whenDialogOpens(const std::function<void(D*)>& onDialog) {
  QTimer::singleShot(0, [onDialog]() {
    for (QWidget* w : QApplication::topLevelWidgets())
      if (D* d = qobject_cast<D*>(w))
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
  whenDialogOpens<AddSegDataDialog>(onDialog);
  tab->editSegDataLine();
}

const FittingParameter* param(FittingTab* tab, const QString& name) {
  for (const FittingParameter& p : tab->getFittingParameters())
    if (p.name == name) return &p;
  return nullptr;
}

int tableRow(QTableWidget* table, const QString& name) {
  for (int r = 0; table && r < table->rowCount(); r++)
    if (table->item(r, 0) && table->item(r, 0)->text() == name) return r;
  return -1;
}

QStringList blockLines(const QString& path, const QString& block) {
  QFile f(path);
  QStringList out;
  if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) return out;
  bool in = false;
  for (const QString& line : QString::fromUtf8(f.readAll()).split('\n')) {
    if (line.trimmed() == "<" + block + ">") { in = true; continue; }
    if (line.trimmed() == "</" + block + ">") break;
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
    SampleBlock b = sampleBlock("total 1 transmission 0.15 1 0.01 thicknessgroup s");
    ok("fit form read (kind, n, vary, sigma as written)",
       b.kind == "transmission" && b.n == "0.15" && b.fit && b.vary && b.sigma == "0.01");
    ok("a group label after the fit form is found", groupLabel("transmission 0.15 1 0.01 thicknessgroup s", "thicknessgroup") == "s");
    b = sampleBlock("yield 0.2 normgroup D");
    ok("plain form read", b.kind == "yield" && b.n == "0.2" && !b.fit);
    ok("no block", sampleBlock("total 1 normgroup D").kind.isEmpty());
    SampleBlock f;
    f.kind = "transmission";
    f.n = "0.2";
    f.fit = true;
    f.vary = true;
    f.sigma = "0.02";
    ok("plain -> fit form in place, the rest kept",
       withSampleBlock("total 1 transmission 0.2 normgroup D", f) == "total 1 transmission 0.2 1 0.02 normgroup D",
       withSampleBlock("total 1 transmission 0.2 normgroup D", f));
    ok("fit form -> plain", withSampleBlock("transmission 0.15 1 0.01 sample 1", SampleBlock{"transmission", "0.15"}) ==
                                "transmission 0.15 sample 1");
    ok("a new block goes before a sample block",
       withSampleBlock("covariance c.txt sample 1", SampleBlock{"yield", "0.1"}) == "covariance c.txt yield 0.1 sample 1");
    ok("removed", withSampleBlock("total 1 transmission 0.15 1 0 thicknessgroup s", SampleBlock()) == "total 1 thicknessgroup s");
    ok("thicknessgroup is a group keyword", groupKeywords().contains("thicknessgroup") && isKeyword("thicknessgroup"));
  }

  // ---- 2. the data-segment dialog --------------------------------------------------------------------------------
  QTemporaryDir tmp;
  const QString source = QString(AZURE2_TEST_SOURCE_DIR) + "/tests/reaction_rate/models/n12C.azr";
  const QString project = tmp.filePath("thickness.azr");
  const QStringList lines = {
      "1 1 1 0.29 0.66 0 180 0 1 0 0 0 0 0 data/d.dat 0 0 total 1 transmission 0.15 1 0.01 thicknessgroup s normgroup D",
      "1 1 2 0.29 0.66 0 180 0 1 0 0 0 0 0 data/g.dat 0 0 sqrtshift 0 0 0 yield 0.12 1 0 thicknessgroup s",
      "1 1 1 0.29 0.66 0 180 0 1 0 0 0 0 0 data/d2.dat 0 0 total 1 transmission 0.2",
      "1 1 2 0.29 0.66 0 180 0 1 0 0 0 0 0 data/g2.dat 0 0 normgroup D"};
  const QStringList tests = {"1 1 1 0.30 0.65 0.001 0 0 0 0 total 1 transmission 0.15 1 0.01"};
  {
    QFile in(source);
    ok("read n12C.azr", in.open(QIODevice::ReadOnly | QIODevice::Text));
    QString text = QString::fromUtf8(in.readAll());
    text.replace("<segmentsData>\n</segmentsData>", "<segmentsData>\n" + lines.join("\n") + "\n</segmentsData>");
    text.replace("<segmentsTest>\n</segmentsTest>", "<segmentsTest>\n" + tests.join("\n") + "\n</segmentsTest>");
    // n12C.azr's <parameterSettings> carry errors for its old segments 3 and 10, which the Fitting tab would lay on
    // these lines: dropped.
    const int ps = text.indexOf("<parameterSettings>\n"), pe = text.indexOf("</parameterSettings>");
    if (ps >= 0 && pe > ps) text.remove(ps + 20, pe - ps - 20);
    QFile out(project);
    out.open(QIODevice::WriteOnly | QIODevice::Text);
    out.write(text.toUtf8());
  }
  AZURESetup setup;
  ok("open the project", setup.readProjectFile(project));
  SegmentsTab* tab = setup.getSegmentsTab();
  SegmentsDataModel* model = tab->getSegmentsDataModel();
  ok("four data lines", model->rowCount(QModelIndex()) == 4);
  auto tail = [&](int row) { return model->getLines()[row].extraTokens; };
  const QString t0 = tail(0), t1 = tail(1);
  auto cell = [&](int row, int col) { return model->data(model->index(row, col), Qt::DisplayRole).toString(); };
  ok("the type column shows the observable and n",
     cell(0, 7).contains("Transmission") && cell(0, 7).contains("n = 0.15 atoms/b, fitted") &&
         cell(1, 7).contains("Self-Shielded Yield") && cell(2, 7).contains("n = 0.2 atoms/b") &&
         !cell(2, 7).contains("fitted"),
     cell(0, 7) + " | " + cell(1, 7) + " | " + cell(2, 7));
  ok("the Groups column shows the thickness group",
     cell(0, SegmentsDataModel::GroupsColumn) == "norm: D, thickness: s" &&
         cell(1, SegmentsDataModel::GroupsColumn) == "thickness: s",
     cell(0, SegmentsDataModel::GroupsColumn));
  ok("line 2 follows line 1 in thickness group s", model->groupLeader(1, "thicknessgroup") == 0);

  bool shown = false;
  editRow(tab, 0, [&](AddSegDataDialog* d) {
    SampleRow* r = d->sampleRow;
    shown = r->kindCombo->currentData().toString() == "transmission" && r->thicknessText->text() == "0.15" &&
            r->fitCheck->isChecked() && r->sigmaText->text() == "0.01" && r->groupCombo->currentText() == "s" &&
            r->groupCombo->findText("s") > 0;
  });
  ok("the Sample row shows line 1's block and group", shown);
  ok("an edit that changes nothing leaves the tail byte for byte", tail(0) == t0, tail(0));

  editRow(tab, 1, [&](AddSegDataDialog* d) {
    ok("line 2: vary 1, prior 0 shown as Fit checked, sigma 0",
       d->sampleRow->fitCheck->isChecked() && d->sampleRow->sigmaText->text() == "0");
    d->sampleRow->thicknessText->setText("0.13");
  });
  ok("n changed in place, every other token kept", tail(1) == "sqrtshift 0 0 0 yield 0.13 1 0 thicknessgroup s", tail(1));

  editRow(tab, 2, [](AddSegDataDialog* d) {
    d->sampleRow->fitCheck->setChecked(true);
    d->sampleRow->sigmaText->setText("0.02");
  });
  ok("a plain block becomes the fit form", tail(2) == "total 1 transmission 0.2 1 0.02", tail(2));

  editRow(tab, 3, [](AddSegDataDialog* d) {
    ok("a line without a block: None, the fields disabled",
       d->sampleRow->kindCombo->currentIndex() == 0 && !d->sampleRow->thicknessText->isEnabled() &&
           !d->sampleRow->groupCombo->isEnabled());
    d->sampleRow->kindCombo->setCurrentIndex(2);  // Self-shielded yield
    d->sampleRow->thicknessText->setText("0.1");
    d->sampleRow->fitCheck->setChecked(true);
    d->setGroupLabel("thicknessgroup", "t");
  });
  ok("a block and a group added after the other blocks", tail(3) == "normgroup D yield 0.1 1 0 thicknessgroup t", tail(3));

  editRow(tab, 3, [](AddSegDataDialog* d) { d->sampleRow->kindCombo->setCurrentIndex(0); });
  ok("None removes the block and its group", tail(3) == "normgroup D", tail(3));
  editRow(tab, 3, [](AddSegDataDialog* d) {
    d->sampleRow->kindCombo->setCurrentIndex(2);
    d->sampleRow->thicknessText->setText("0.1");
    d->sampleRow->fitCheck->setChecked(true);
    d->setGroupLabel("thicknessgroup", "t");
  });

  // ---- 3. the test-segment dialog --------------------------------------------------------------------------------
  SegmentsTestModel* tmodel = tab->getSegmentsTestModel();
  const QString x0 = tmodel->getLines()[0].extraTokens;
  QTableView* tview = tab->findChild<QTableView*>("segmentsTestView");
  ok("the test segments' view", tview != nullptr);
  if (tview) {
    tview->selectRow(0);
    bool testShown = false;
    whenDialogOpens<AddSegTestDialog>([&](AddSegTestDialog* d) {
      testShown = d->sampleRow->kindCombo->currentData().toString() == "transmission" &&
                  d->sampleRow->thicknessText->text() == "0.15" && d->sampleRow->fitCheck == nullptr;
    });
    tab->editSegTestLine();
    ok("the test dialog shows kind and n, no Fit box", testShown);
    ok("unchanged: the test line's tail (fit form included) is kept", tmodel->getLines()[0].extraTokens == x0,
       tmodel->getLines()[0].extraTokens);
    tview->selectRow(0);
    whenDialogOpens<AddSegTestDialog>([](AddSegTestDialog* d) { d->sampleRow->thicknessText->setText("0.25"); });
    tab->editSegTestLine();
    ok("a new n on a test line keeps its fit form", tmodel->getLines()[0].extraTokens == "total 1 transmission 0.25 1 0.01",
       tmodel->getLines()[0].extraTokens);
  }

  // ---- 4. the Fitting tab ----------------------------------------------------------------------------------------
  FittingTab* fitting = setup.getFittingTab();
  QMetaObject::invokeMethod(fitting, "refreshParameters");
  QTableWidget* table = fitting->findChild<QTableWidget*>("thicknessParamsTable");
  const FittingParameter* p1 = param(fitting, "segment_1_thickness");
  const FittingParameter* p2 = param(fitting, "segment_2_thickness");
  const FittingParameter* p3 = param(fitting, "segment_3_thickness");
  const FittingParameter* p4 = param(fitting, "segment_4_thickness");
  ok("four thickness rows (every line has the fit form now)", p1 && p2 && p3 && p4 && table && table->rowCount() == 4);
  if (p1) ok("segment_1_thickness: 0.15, prior 0.01, fitted, a Minuit index",
             p1->value == 0.15 && p1->error == 0.01 && p1->useAsNuisance && p1->minuitIndex >= 0 && p1->tiedTo.isEmpty());
  if (p2)
    ok("segment_2_thickness is tied to segment 1 with its value",
       p2->tiedTo == "tied to segment 1 (thicknessgroup s)" && p2->value == 0.15 && p2->tiedLeader == 0 &&
           p2->minuitIndex >= 0,
       p2->tiedTo + QString(" %1").arg(p2->value));
  if (p3) ok("segment_3_thickness: 0.2, prior 0.02", p3->value == 0.2 && p3->error == 0.02 && p3->tiedTo.isEmpty());
  if (p4) ok("segment_4_thickness leads its own group t", p4->value == 0.1 && p4->tiedTo.isEmpty());
  const int r2 = tableRow(table, "segment_2_thickness");
  if (r2 >= 0) {
    QTableWidgetItem* nuis = table->item(r2, 6);
    ok("the follower's row: no check box, greyed, the leader named",
       nuis && nuis->text() == "tied to segment 1 (thicknessgroup s)" && !(nuis->flags() & Qt::ItemIsUserCheckable) &&
           !(table->item(r2, 1)->flags() & Qt::ItemIsEnabled));
  }
  const int r1 = tableRow(table, "segment_1_thickness");
  ok("the leader's row is in the table", r1 >= 0);
  if (r1 >= 0) {
    table->item(r1, 1)->setText("0.18");
    ok("a value edit reaches line 1's block, the rest of the tail kept",
       tail(0) == "total 1 transmission 0.18 1 0.01 thicknessgroup s normgroup D", tail(0));
    p2 = param(fitting, "segment_2_thickness");
    ok("the follower shows the leader's new value", p2 && p2->value == 0.18 && r2 >= 0 && table->item(r2, 1)->text() == "0.18");
    table->item(r1, 4)->setText("0.03");
    table->item(r1, 6)->setCheckState(Qt::Unchecked);
    ok("Error is the prior sigma, Nuisance the vary flag",
       tail(0) == "total 1 transmission 0.18 0 0.03 thicknessgroup s normgroup D", tail(0));
    p1 = param(fitting, "segment_1_thickness");
    ok("a fixed thickness keeps its Minuit index (a segment parameter)", p1 && p1->minuitIndex >= 0);
    table->item(r1, 6)->setCheckState(Qt::Checked);
  }
  // the Segments tab changes a line: the row follows
  model->setLineExtras(2, "total 1 transmission 0.21 1 0.02");
  p3 = param(fitting, "segment_3_thickness");
  ok("the row follows the line", p3 && p3->value == 0.21);
  model->setLineExtras(2, "total 1 transmission 0.21");
  ok("a plain block has no row", param(fitting, "segment_3_thickness") == nullptr);
  model->setLineExtras(2, "total 1 transmission 0.21 1 0.02");

  // ---- 5. written back -------------------------------------------------------------------------------------------
  const QString saved = tmp.filePath("saved.azr");
  ok("project written", setup.writeProjectFile(saved));
  const QStringList got = blockLines(saved, "segmentsData");
  const QStringList want = {
      "1 1 1 0.29 0.66 0 180 0 1 0 0 0 0 0 data/d.dat 0 0 total 1 transmission 0.18 1 0.03 thicknessgroup s normgroup D",
      "1 1 2 0.29 0.66 0 180 0 1 0 0 0 0 0 data/g.dat 0 0 sqrtshift 0 0 0 yield 0.13 1 0 thicknessgroup s",
      "1 1 1 0.29 0.66 0 180 0 1 0 0 0 0 0 data/d2.dat 0 0 total 1 transmission 0.21 1 0.02",
      "1 1 2 0.29 0.66 0 180 0 1 0 0 0 0 0 data/g2.dat 0 0 normgroup D yield 0.1 1 0 thicknessgroup t"};
  ok("four lines written", got.size() == 4, QString::number(got.size()));
  for (int i = 0; i < std::min(got.size(), want.size()); i++)
    ok(QString("line %1 written as expected").arg(i + 1).toStdString(), got[i] == want[i], got[i]);
  const QStringList gotTests = blockLines(saved, "segmentsTest");
  ok("the test line written with its new n", gotTests.size() == 1 && gotTests[0].endsWith("total 1 transmission 0.25 1 0.01"),
     gotTests.join(" | "));

  std::cout << (failures ? "FAILED" : "PASSED") << std::endl;
  return failures ? 1 : 0;
}
