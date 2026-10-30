// Headless checks on the AME2020 mass table in the GUI (MassTable) and the places that fill from it: the pair
// dialog's isotope fields and "Fill masses and S" button, the S column's tooltip and warning colour, and the
// Pairs tab's "Masses from AME2020..." dialog (MassCheckDialog).
//
// The GUI embeds pyazr/data/mass_1.mas20, the file pyazr.masses reads; the reference numbers below are the ones
// tests/pyazr/masses_test.py checks, so the two front ends are pinned to the same values.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.  With AZURE2_MASS_TABLE_PNG=<stem> it
// saves the pair dialog and the masses dialog as <stem>_pair.png and <stem>_masses.png for a layout review.

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QTableWidget>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "AddPairDialog.h"
#include "Config.h"
#include "MassCheckDialog.h"
#include "MassTable.h"
#include "PairsModel.h"
#include "PairsTab.h"
struct SegPairs {int firstPair; int secondPair;};

Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if (!detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if (!cond) fails++;
}

static PairsData pair(int type, int z1, double m1, int z2, double m2, double S) {
  PairsData p{};
  p.lightJ = type == 10 ? 1.0 : 0.5; p.lightPi = 1; p.lightZ = z1; p.lightM = m1; p.lightG = 0.;
  p.heavyJ = 0.; p.heavyPi = 1; p.heavyZ = z2; p.heavyM = m2; p.heavyG = 0.;
  p.excitationEnergy = 0.; p.seperationEnergy = S; p.channelRadius = type == 10 ? 0. : 4.0;
  p.pairType = type; p.ecMultMask = 0;
  return p;
}

static void addTo(PairsModel& model, const PairsData& p) {
  const int row = model.numPairs();
  model.insertRows(row, 1);
  const QVariant v[PairsData::SIZE] = {p.lightJ, p.lightPi, p.lightZ, p.lightM, p.lightG, p.heavyJ, p.heavyPi, p.heavyZ,
                                       p.heavyM, p.heavyG, p.excitationEnergy, p.seperationEnergy, p.channelRadius,
                                       p.pairType, p.ecMultMask};
  for (int c = 0; c < PairsData::SIZE; c++) model.setData(model.index(row, c), v[c]);
}

static void typeIsotope(QLineEdit* edit, const QString& text) {
  edit->setText(text);
  edit->setModified(true);   // as a keystroke would
  emit edit->editingFinished();
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  const MassTable& mt = MassTable::instance();

  std::cout << "== table" << std::endl;
  ok("the embedded AME2020 table is read in full", mt.isLoaded() && mt.size() == 3558, QString::number(mt.size()));
  struct Ref { int Z, A; double excess, mass; bool estimated; };
  const Ref refs[] = {{0, 1, 8071.31806, 1.00866491590, false}, {1, 1, 7288.971064, 1.007825031898, false},
                      {2, 4, 2424.91587, 4.00260325413, false}, {6, 12, 0.0, 12.0, false},
                      {7, 13, 5345.481, 13.005738609, false}, {8, 16, -4737.00217, 15.99491461926, false},
                      {3, 3, 28667.0, 3.030775, true}};
  bool rowsOk = true;
  for (const Ref& r : refs) {
    const MassEntry* e = mt.lookup(r.Z, r.A);
    rowsOk = rowsOk && e && e->massExcess == r.excess && e->atomicMass == r.mass && e->estimated == r.estimated;
  }
  ok("reference rows (n, 1H, 4He, 12C, 13N, 16O, estimated 3Li) equal pyazr's", rowsOk);
  int Z = -1, A = -1;
  bool names = true;
  const struct { const char* text; int Z, A; } cases[] = {{"p", 1, 1}, {"n", 0, 1}, {"a", 2, 4}, {"alpha", 2, 4},
                                                          {"helion", 2, 3}, {"3He", 2, 3}, {"12C", 6, 12}, {"C12", 6, 12},
                                                          {"c-12", 6, 12}, {"N15", 7, 15}, {"1n", 0, 1}};
  for (const auto& c : cases) names = names && mt.parseNuclide(c.text, &Z, &A) && Z == c.Z && A == c.A;
  ok("isotope names and shorthands parse as in pyazr", names);
  ok("\"h\", \"Xx12\" and \"12\" are refused",
     !mt.parseNuclide("h", &Z, &A) && !mt.parseNuclide("Xx12", &Z, &A) && !mt.parseNuclide("12", &Z, &A));
  double S = 0.;
  const bool s1 = mt.separationEnergy(1, 1, 6, 12, &S) && S == 1.943490064;
  const bool s2 = mt.separationEnergy(2, 4, 6, 12, &S) && S == 7.16191804;
  const bool s3 = mt.separationEnergy(0, 1, 6, 12, &S) && S == 4.94630873;
  ok("S(p+12C), S(a+12C), S(n+12C) equal pyazr's to the last digit", s1 && s2 && s3);

  // the 13N test file's pair: bare-proton mass, S 20 eV high
  const PairMassCheck c13 = mt.check(pair(0, 1, 1.00728, 6, 12.0, 1.94351));
  ok("check: a bare proton is identified as 1H, S +0.020 keV, mass offset 5.45e-4 u",
     c13.fillable && c13.hasS && std::fabs(c13.dSkeV - 0.019936) < 1e-6 && std::fabs(c13.massOffset - 5.45031898e-4) < 1e-12,
     QString("dS %1 keV, offset %2 u").arg(c13.dSkeV, 0, 'f', 6).arg(c13.massOffset));
  const PairMassCheck cg = mt.check(pair(10, 0, 0.0, 7, 13.0057, 0.0));
  ok("check: a gamma pair is fillable through its heavy particle and S is not compared",
     cg.fillable && !cg.hasS && cg.heavy && cg.heavy->A == 13);
  ok("check: a nuclide the table lacks makes the pair unfillable", !mt.check(pair(0, 6, 1.0, 6, 12.0, 1.0)).fillable);

  std::cout << "== pair dialog" << std::endl;
  {
    AddPairDialog d;
    typeIsotope(d.lightIsotopeText, "p");
    typeIsotope(d.heavyIsotopeText, "12C");
    ok("isotope fields fill Z and the atomic mass",
       d.lightZText->text() == "1" && d.lightMText->text() == "1.007825031898" && d.heavyZText->text() == "6" &&
           d.heavyMText->text() == "12",
       d.lightZText->text() + " " + d.lightMText->text() + " " + d.heavyZText->text() + " " + d.heavyMText->text());
    d.fillButton->click();
    ok("the fill button sets S from the table", d.seperationEnergyText->text() == "1.943490064", d.seperationEnergyText->text());
    ok("the check line says consistent", d.massCheckLabel->text().contains("consistent"), d.massCheckLabel->text());
    d.seperationEnergyText->setText("1.9485");
    ok("an S 5 keV off is flagged with its difference", d.massCheckLabel->text().contains("differs by 5.01 keV"),
       d.massCheckLabel->text());
    d.lightMText->setText("1.00728");
    emit d.lightIsotopeText->editingFinished();   // focus passing through, no new entry
    ok("an isotope field left alone does not overwrite a typed mass", d.lightMText->text() == "1.00728");
    ok("a mass off the table is reported, not coloured below 1e-3 u",
       d.massCheckLabel->text().contains("differ from AME2020") &&
           !d.massCheckLabel->text().section("<br/>", 1).contains("color"),
       d.massCheckLabel->text());
    d.lightMText->setText("1.08664");
    ok("a mistyped mass (0.078 u off) is coloured", d.massCheckLabel->text().section("<br/>", 1).contains("color"),
       d.massCheckLabel->text());
    d.lightMText->setText("1.00728");
    if (const char* stem = std::getenv("AZURE2_MASS_TABLE_PNG")) {
      d.show();
      d.grab().save(QString("%1_pair.png").arg(stem));
    }
    typeIsotope(d.heavyIsotopeText, "12Fe");
    ok("an isotope the table lacks is reported and changes nothing",
       d.heavyMText->text() == "12" && d.massCheckLabel->text().contains("no nuclide"), d.massCheckLabel->text());
  }
  {
    AddPairDialog d;
    d.pairTypeCombo->setCurrentIndex(1);   // particle, gamma
    ok("the light isotope field is off for a gamma pair", !d.lightIsotopeText->isEnabled());
    d.heavyZText->setText("7");
    d.heavyMText->setText("13.0057");
    d.fillButton->click();
    ok("a gamma pair gets the compound nucleus's mass, S stays 0",
       d.heavyMText->text() == "13.005738609" && d.seperationEnergyText->text() == "0.0",
       d.heavyMText->text() + " / " + d.seperationEnergyText->text());
  }

  std::cout << "== pairs table and masses dialog" << std::endl;
  PairsTab tab;
  PairsModel* model = tab.getPairsModel();
  addTo(*model, pair(0, 1, 1.00728, 6, 12.0, 1.94351));        // bare proton, S 20 eV high
  addTo(*model, pair(0, 0, 1.0, 6, 12.0, 4.95130873));         // mass numbers, S 5 keV high
  addTo(*model, pair(10, 0, 0.0, 7, 13.0057, 0.0));            // gamma pair, 4-decimal mass
  addTo(*model, pair(0, 2, 4.00260325413, 6, 12.0, 7.16191804));  // consistent already (another nucleus)
  const QModelIndex s1i = model->index(0, 11), s2i = model->index(1, 11), s3i = model->index(2, 11);
  ok("S column: a 0.02 keV pair is not coloured, its tooltip says consistent",
     !model->data(s1i, Qt::ForegroundRole).isValid() && model->data(s1i, Qt::ToolTipRole).toString().contains("Consistent"));
  ok("S column: a 5 keV pair is coloured and its tooltip gives the table's S",
     model->data(s2i, Qt::ForegroundRole).isValid() && model->data(s2i, Qt::ToolTipRole).toString().contains("4.946309"),
     model->data(s2i, Qt::ToolTipRole).toString());
  ok("S column: no tooltip for a gamma pair", !model->data(s3i, Qt::ToolTipRole).isValid());

  {
    MassCheckDialog dialog(model);
    QTableWidget* t = dialog.table();
    if (const char* stem = std::getenv("AZURE2_MASS_TABLE_PNG")) {
      dialog.show();
      dialog.grab().save(QString("%1_masses.png").arg(stem));
    }
    ok("mass cells: a mass number (0.0087 u off) is coloured, a bare proton (5e-4 u) is not",
       t->item(1, 2)->foreground().style() != Qt::NoBrush && t->item(0, 2)->foreground().style() == Qt::NoBrush);
    ok("the dialog lists every pair, the differing ones ticked",
       t->rowCount() == 4 && t->item(0, 0)->checkState() == Qt::Checked && t->item(1, 0)->checkState() == Qt::Checked &&
           t->item(2, 0)->checkState() == Qt::Checked && t->item(3, 0)->checkState() == Qt::Unchecked);
    QList<QPair<int, PairsData>> f = dialog.filledPairs();
    const bool values = f.size() == 3 && f[0].first == 0 && f[0].second.lightM == 1.007825031898 &&
                        f[0].second.seperationEnergy == 1.943490064 && f[1].second.lightM == 1.0086649159 &&
                        f[1].second.seperationEnergy == 4.94630873 && f[2].second.heavyM == 13.005738609 &&
                        f[2].second.seperationEnergy == 0.0;
    ok("filled values: table masses and S, gamma pair mass only", values);
    dialog.separationCheck()->setChecked(false);
    f = dialog.filledPairs();
    ok("without the S box only masses change", f.size() == 3 && f[0].second.seperationEnergy == 1.94351 &&
                                                f[1].second.seperationEnergy == 4.95130873);
    dialog.separationCheck()->setChecked(true);
    t->item(1, 0)->setCheckState(Qt::Unchecked);
    f = dialog.filledPairs();
    ok("an unticked pair is left alone", f.size() == 2 && f[0].first == 0 && f[1].first == 2);
    t->item(1, 0)->setCheckState(Qt::Checked);
    for (const auto& change : dialog.filledPairs()) tab.editPair(change.second, change.first, true);
  }
  const QList<PairsData> after = model->getPairs();
  ok("applied through PairsTab::editPair, the model holds the table's values",
     after[0].lightM == 1.007825031898 && after[0].seperationEnergy == 1.943490064 && after[1].lightM == 1.0086649159 &&
         after[1].seperationEnergy == 4.94630873 && after[2].heavyM == 13.005738609 && after[3].lightM == 4.00260325413);
  bool consistent = true;
  for (const PairsData& p : after) {
    const PairMassCheck c = mt.check(p);
    consistent = consistent && c.massOffset == 0.0 && (!c.hasS || c.dSkeV == 0.0);
  }
  ok("every pair now agrees with the table exactly", consistent);
  ok("and nothing is coloured", !model->data(model->index(1, 11), Qt::ForegroundRole).isValid());

  if (fails) {
    std::cout << fails << " check(s) FAILED" << std::endl;
    return 1;
  }
  std::cout << "all checks passed" << std::endl;
  return 0;
}
