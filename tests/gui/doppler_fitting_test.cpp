// The Fitting tab's Doppler Temperatures rows: one doppler_<n> (K) per active Experimental Effects line with the fit
// form 'doppler T vary sigma [Tmax]', n its place among the active lines (EData::FillMnParams), Value = T, Error =
// the prior's sigma in K, Nuisance = fitted.  A plain 'doppler T' line has no row.  Edits in either tab reach the
// other, a fitted T takes a Minuit index, and the project file carries the values.
//
// Project: tests/reaction_rate/models/n12C.azr with its two inactive <targetInt> lines and three active ones added
// here: "doppler 300 1 20 900" (effect 1), a plain "doppler 310" (effect 2) and "doppler 320 0 0" (effect 3).
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
  QTemporaryDir tmp;

  // The project: n12C.azr with three active Doppler lines after its two inactive ones.
  QFile src(QString(AZURE2_TEST_SOURCE_DIR) + "/tests/reaction_rate/models/n12C.azr");
  ok("read n12C.azr", src.open(QIODevice::ReadOnly | QIODevice::Text));
  QString text = QTextStream(&src).readAll();
  const QString added = "1 \"1,2\" 400 0 0 0 0 \"\" 0 0 0 0 \"\" 0 doppler 300 1 20 900\n"
                        "1 \"3\" 400 0 0 0 0 \"\" 0 0 0 0 \"\" 0 doppler 310\n"
                        "1 \"4\" 400 0 0 0 0 \"\" 0 0 0 0 \"\" 0 doppler 320 0 0\n";
  text.replace("</targetInt>", added + "</targetInt>");
  const QString project = tmp.filePath("doppler.azr");
  {
    QFile f(project);
    f.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream(&f) << text;
  }

  AZURESetup setup;
  ok("read the Doppler project", setup.readProjectFile(project));
  FittingTab *fitting = setup.getFittingTab();
  TargetIntModel *effects = setup.getTargetIntTab()->getTargetIntModel();
  QMetaObject::invokeMethod(fitting, "refreshParameters");

  // the rows of the model: the two inactive lines, then effects 1-3
  const int row1 = 2;
  ok("five lines, the Doppler ones from row 2", effects->getLines().size() == 5 && effects->getLines().at(row1).dopplerFit);
  const FittingParameter *p = find(fitting, "doppler_1");
  const FittingParameter *p3 = find(fitting, "doppler_3");
  ok("two Doppler rows, doppler_1 and doppler_3 (the plain line has none)",
     p && p3 && !find(fitting, "doppler_2") && countCategory(fitting, "doppler") == 2);
  if (p)
    ok("doppler_1: 300 K, prior 20 K, fitted, with a Minuit index",
       p->value == 300.0 && p->error == 20.0 && p->useAsNuisance && p->minuitIndex >= 0);
  if (p3)
    ok("doppler_3: 320 K, no prior, fixed, no Minuit index",
       p3->value == 320.0 && p3->error == 0.0 && !p3->useAsNuisance && p3->minuitIndex < 0);

  // The Experimental Effects tab: a 30 K prior and fixed -> the row follows and loses its index.
  effects->setData(effects->index(row1, 48), 30.0, Qt::EditRole);
  effects->setData(effects->index(row1, 47), false, Qt::EditRole);
  p = find(fitting, "doppler_1");
  ok("row follows the line (prior 30 K, fixed, no index)", p && p->error == 30.0 && !p->useAsNuisance && p->minuitIndex < 0);
  effects->setData(effects->index(row1, 47), true, Qt::EditRole);
  p = find(fitting, "doppler_1");
  ok("fitted again: an index again", p && p->useAsNuisance && p->minuitIndex >= 0);

  // No fit form on the line: no row; back: the row is back.
  effects->setData(effects->index(row1, 46), false, Qt::EditRole);
  ok("plain line, no row", find(fitting, "doppler_1") == nullptr && countCategory(fitting, "doppler") == 1);
  effects->setData(effects->index(row1, 46), true, Qt::EditRole);
  ok("fit form back, row back", find(fitting, "doppler_1") != nullptr);

  // Edits in the Fitting tab's table go back to the line (a fitted T is written back like a density).
  QTableWidget *table = nullptr;
  int trow = -1;
  for (QTableWidget *t : fitting->findChildren<QTableWidget *>())
    for (int r = 0; r < t->rowCount(); r++)
      if (t->item(r, 0) && t->item(r, 0)->text() == "doppler_1") { table = t; trow = r; }
  ok("the Doppler rows are in a table of their own", table && trow >= 0 && table->rowCount() == 2);
  if (table && trow >= 0) {
    table->item(trow, 1)->setText("450");
    QList<TargetIntData> lines = effects->getLines();
    ok("a new value reaches the line's temperature", lines.at(row1).dopplerTemperature == 450.0);
    table->item(trow, 4)->setText("35");
    table->item(trow, 6)->setCheckState(Qt::Unchecked);
    lines = effects->getLines();
    ok("an Error edit is the prior sigma in K, the Nuisance box fixes it",
       lines.at(row1).dopplerSigma == 35.0 && !lines.at(row1).dopplerVary);
    p = find(fitting, "doppler_1");
    ok("and the row loses its Minuit index", p && p->minuitIndex < 0);
    table->item(trow, 6)->setCheckState(Qt::Checked);
  }

  // The written project carries what the tabs hold.
  const QString out = tmp.filePath("out.azr");
  ok("project written", setup.writeProjectFile(out));
  QFile f(out);
  f.open(QIODevice::ReadOnly | QIODevice::Text);
  const QString written = QTextStream(&f).readAll();
  ok("effect 1 written as 'doppler 450 1 35 900'", written.contains(" doppler 450 1 35 900\n"));
  ok("effect 2 stays plain, effect 3 keeps its fit form",
     written.contains(" doppler 310\n") && written.contains(" doppler 320 0 0\n"));

  std::cout << (failures ? "FAILED" : "PASSED") << std::endl;
  return failures ? 1 : 0;
}
