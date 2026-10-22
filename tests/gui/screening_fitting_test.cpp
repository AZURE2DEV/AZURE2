// The Fitting tab's Screening rows: one screening_<n> (eV) per active Experimental Effects line with
// electron screening, n its place among the active lines (EData::FillMnParams), Value = U_e, Error = the
// prior's sigma, Nuisance = fitted.  Edits in either tab reach the other, a fitted U_e takes a Minuit
// index, and the project file written back carries the values.
//
// Project: tests/electron_screening/grad.azr (two inactive 13N lines, then the screening line on 6,7 with
// U_e 1000 eV, fitted, prior 200 eV).
#include <QApplication>
#include <QFile>
#include <QTemporaryDir>
#include <QTextStream>
#include <QTableWidget>
#include <iostream>
#include "AZURESetup.h"
#include "Config.h"
#include "FittingTab.h"
#include "TargetIntModel.h"
#include "TargetIntTab.h"
struct SegPairs {int firstPair; int secondPair;};
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

namespace {
int failures = 0;
void ok(const std::string &what, bool cond) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what << std::endl;
  if (!cond) failures++;
}
const FittingParameter *find(FittingTab *tab, const QString &name) {
  for (const FittingParameter &p : tab->getFittingParameters())
    if (p.name == name) return &p;
  return nullptr;
}
int countCategory(FittingTab *tab, const QString &category) {
  int n = 0;
  for (const FittingParameter &p : tab->getFittingParameters())
    if (p.category == category) n++;
  return n;
}
}  // namespace

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  const QString project = QString(AZURE2_TEST_SOURCE_DIR) + "/tests/electron_screening/grad.azr";
  QTemporaryDir tmp;

  AZURESetup setup;
  ok("read the screening project", setup.readProjectFile(project));
  FittingTab *fitting = setup.getFittingTab();
  TargetIntModel *effects = setup.getTargetIntTab()->getTargetIntModel();
  QMetaObject::invokeMethod(fitting, "refreshParameters");  // a slot

  const FittingParameter *p = find(fitting, "screening_1");
  ok("one screening row, named after the line's place among the active lines (screening_1)",
     p && countCategory(fitting, "screening") == 1);
  if (p) {
    ok("value = U_e 1000 eV, error = prior 200 eV, fitted", p->value == 1000.0 && p->error == 200.0 && p->useAsNuisance);
    ok("a fitted U_e takes a Minuit index", p->minuitIndex >= 0);
  }

  // The Experimental Effects tab changes U_e and fixes it: the row follows, and loses its index.
  const int row = 2;  // the screening line is the third <targetInt> line
  effects->setData(effects->index(row, 36), 750.0, Qt::EditRole);
  effects->setData(effects->index(row, 37), false, Qt::EditRole);
  p = find(fitting, "screening_1");
  ok("row follows the Experimental Effects tab (U_e 750, fixed, no index)",
     p && p->value == 750.0 && !p->useAsNuisance && p->minuitIndex < 0);

  // Switching screening off removes the row; on again brings it back.
  effects->setData(effects->index(row, 35), false, Qt::EditRole);
  ok("no screening, no row", countCategory(fitting, "screening") == 0);
  effects->setData(effects->index(row, 35), true, Qt::EditRole);
  ok("screening back, row back", find(fitting, "screening_1") != nullptr);

  // Activating an earlier line renumbers it among the active lines.
  effects->setData(effects->index(0, 0), Qt::Checked, Qt::CheckStateRole);
  ok("an earlier line made active renames the row screening_2",
     find(fitting, "screening_2") != nullptr && find(fitting, "screening_1") == nullptr);
  effects->setData(effects->index(0, 0), Qt::Unchecked, Qt::CheckStateRole);
  ok("and back to screening_1", find(fitting, "screening_1") != nullptr);

  // Edits in the Fitting tab's Screening table go back to the Experimental Effects line.
  QTableWidget *table = nullptr;
  int trow = -1;
  for (QTableWidget *t : fitting->findChildren<QTableWidget *>())
    for (int r = 0; r < t->rowCount(); r++)
      if (t->item(r, 0) && t->item(r, 0)->text() == "screening_1") { table = t; trow = r; }
  ok("the screening row is in a table of its own", table && trow >= 0 && table->rowCount() == 1);
  if (table && trow >= 0) {
    table->item(trow, 1)->setText("820");
    table->item(trow, 4)->setText("150");
    table->item(trow, 6)->setCheckState(Qt::Checked);
    TargetIntData d = effects->getLines().at(row);
    ok("Fitting-tab edits reach the line (U_e 820, prior 150, fitted)",
       d.screeningUe == 820.0 && d.screeningSigma == 150.0 && d.screeningVary);
    p = find(fitting, "screening_1");
    ok("and the row takes a Minuit index again", p && p->minuitIndex >= 0);
    table->item(trow, 1)->setText("750");
    table->item(trow, 4)->setText("200");
    table->item(trow, 6)->setCheckState(Qt::Unchecked);
  }

  // The written project carries what the tabs hold.
  const QString out = tmp.filePath("out.azr");
  ok("project written", setup.writeProjectFile(out));
  QFile f(out);
  f.open(QIODevice::ReadOnly | QIODevice::Text);
  const QString text = QTextStream(&f).readAll();
  ok("the line is written with U_e 750, fixed, prior 200", text.contains("screening 750 0 200"));

  std::cout << (failures ? "FAILED" : "PASSED") << std::endl;
  return failures ? 1 : 0;
}
