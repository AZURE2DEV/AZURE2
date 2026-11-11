// The Fitting tab's Target Densities rows: one density_<n> (atoms/cm^2) per active Experimental Effects line that
// holds a densityfit parameter (a densitygroup's first line; its members follow it), n its place among the active
// lines (EData::SetUpDensityParams), Value = the density, Error = the prior's sigma in atoms/cm^2, Nuisance = fitted.
// Edits in either tab reach the other (a density to every line of its group), a percentage prior on the line
// survives a change of the density, a fitted density takes a Minuit index, and the project file carries the values.
//
// Project: tests/target_density/fit_group.azr (lines on segments 3 and 4 in densitygroup t; the first holds
// densityfit 1 0 1.5e18 at 1e18, the second carries its own 2e18 and only the label).
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
  const QString project = QString(AZURE2_TEST_SOURCE_DIR) + "/tests/target_density/fit_group.azr";
  QTemporaryDir tmp;

  AZURESetup setup;
  ok("read the density project", setup.readProjectFile(project));
  FittingTab *fitting = setup.getFittingTab();
  TargetIntModel *effects = setup.getTargetIntTab()->getTargetIntModel();
  QMetaObject::invokeMethod(fitting, "refreshParameters");

  const FittingParameter *p = find(fitting, "density_1");
  ok("one density row for the group (density_1; the member has none)", p && countCategory(fitting, "density") == 1);
  if (p) {
    ok("value 1e18, no prior, fitted", p->value == 1e18 && p->error == 0.0 && p->useAsNuisance);
    ok("a fitted density takes a Minuit index", p->minuitIndex >= 0);
  }

  // The Experimental Effects tab: a 2 % prior and fixed -> the row follows and loses its index.
  effects->setData(effects->index(0, 41), "2%", Qt::EditRole);
  effects->setData(effects->index(0, 40), false, Qt::EditRole);
  p = find(fitting, "density_1");
  ok("row follows the line (prior 2 % = 2e16, fixed, no index)",
     p && p->error == 2e16 && !p->useAsNuisance && p->minuitIndex < 0);
  effects->setData(effects->index(0, 40), true, Qt::EditRole);

  // No densityfit on the leader: no row; back: the row is back.
  effects->setData(effects->index(0, 39), false, Qt::EditRole);
  ok("no densityfit, no row", countCategory(fitting, "density") == 0);
  effects->setData(effects->index(0, 39), true, Qt::EditRole);
  ok("densityfit back, row back", find(fitting, "density_1") != nullptr);

  // Edits in the Fitting tab's table go back to the lines.
  QTableWidget *table = nullptr;
  int trow = -1;
  for (QTableWidget *t : fitting->findChildren<QTableWidget *>())
    for (int r = 0; r < t->rowCount(); r++)
      if (t->item(r, 0) && t->item(r, 0)->text() == "density_1") { table = t; trow = r; }
  ok("the density row is in a table of its own", table && trow >= 0 && table->rowCount() == 1);
  if (table && trow >= 0) {
    table->item(trow, 1)->setText("1.05e18");
    QList<TargetIntData> lines = effects->getLines();
    ok("a new value reaches both lines of the group", lines.at(0).density == 1.05e18 && lines.at(1).density == 1.05e18);
    ok("the percentage prior is kept", lines.at(0).densitySigma == "2%");
    table->item(trow, 4)->setText("3e16");
    table->item(trow, 6)->setCheckState(Qt::Unchecked);
    lines = effects->getLines();
    ok("an Error edit is written in atoms/cm^2, the Nuisance box fixes it",
       lines.at(0).densitySigma == "3e+16" && !lines.at(0).densityVary);
    p = find(fitting, "density_1");
    ok("and the row loses its Minuit index", p && p->minuitIndex < 0);
    table->item(trow, 6)->setCheckState(Qt::Checked);
  }

  // The written project carries what the tabs hold.
  const QString out = tmp.filePath("out.azr");
  ok("project written", setup.writeProjectFile(out));
  QFile f(out);
  f.open(QIODevice::ReadOnly | QIODevice::Text);
  const QString text = QTextStream(&f).readAll();
  ok("leader written with densityfit 1 3e+16 1.5e+18 densitygroup t",
     text.contains("densityfit 1 3e+16 1.5e+18 densitygroup t\n"));
  ok("member written with the label only", text.count("densityfit") == 1 && text.count("densitygroup t") == 2);
  const QString block = text.mid(text.indexOf("<targetInt>"), text.indexOf("</targetInt>") - text.indexOf("<targetInt>"));
  ok("both lines carry the new density", block.count("1.05e+18") == 2);

  std::cout << (failures ? "FAILED" : "PASSED") << std::endl;
  return failures ? 1 : 0;
}
