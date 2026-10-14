// Several compound nuclei in one project, GUI side (headless):
//   1. the compound nucleus of a pair (PairsModel::nucleusOf) and its labels;
//   2. the abundance table follows the pairs, and the <abundances> block
//      round-trips through PairsTab::readAbundances / writeAbundances;
//   3. levels belong to a nucleus: J^pi options and channels per nucleus,
//      two 0+ ground states of different nuclei stay distinct on load, and
//      writeNuclearFile writes one block per nucleus;
//   4. the "natural 1" keyword tail: SegmentKeywords, the segment dialogs'
//      check box (hidden with one nucleus, refused for a ratio or a phase
//      shift) and the "nat" flag on the row label;
//   5. the Fitting tab's abundance parameters sit between the levels and
//      the normalizations and count only when fitted.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QTextStream>
#include <QString>
#include <QStringList>
#include <iostream>
#include "PairsTab.h"
#include "PairsModel.h"
#include "AbundancesModel.h"
#include "LevelsTab.h"
#include "LevelsModel.h"
#include "ChannelsModel.h"
#include "SegmentsTab.h"
#include "SegmentsDataModel.h"
#include "SegmentsTestModel.h"
#include "SegmentKeywords.h"
#include "AddSegDataDialog.h"
#include "AddSegTestDialog.h"
#include "AddLevelDialog.h"
#include "FittingTab.h"
#include "Config.h"
struct SegPairs {int firstPair; int secondPair;};

// Defined by AZURE2.cpp, which belongs to the executable rather than the GUI
// library, so this test supplies its own. They are never called from here.
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

static QStringList tokens(const QString& text) {
  QStringList out;
  for (const QString& line : text.split('\n', Qt::SkipEmptyParts)) {
    QString t = line.simplified();
    if (!t.isEmpty() && !t.startsWith("</")) out << t;
  }
  return out;
}

static void addPair(PairsModel& pairs, int row, const QVariant (&v)[PairsData::SIZE]) {
  pairs.insertRows(row, 1);
  for (int c = 0; c < PairsData::SIZE; c++) pairs.setData(pairs.index(row, c), v[c]);
}

static QString jpiList(const QList<SpinParityOption>& options) {
  QStringList out;
  for (const SpinParityOption& o : options) out << AddLevelDialog::spinParityLabel(o.jValue, o.piValue);
  return out.join(" ");
}

// Pairs: p+12C (-> 13N), gamma+13N (-> 13N), p+13C (-> 14N).  Columns as PairsData.
static const QVariant P_12C[PairsData::SIZE] = {0.5, 1, 1, 1.00782, 0., 0.0, 1, 6, 12.0, 0., 0., 1.9435, 5.0, 0, 0};
static const QVariant G_13N[PairsData::SIZE] = {1.0, -1, 0, 0., 0., 0.5, -1, 7, 13.0057, 0., 0., 0., 0., 10, 0};
static const QVariant P_13C[PairsData::SIZE] = {0.5, 1, 1, 1.00782, 0., 0.5, -1, 6, 13.00335, 0., 0., 7.5506, 5.0, 0, 0};

int main(int argc, char** argv) {
  QApplication app(argc, argv);

  // ---- 1. nuclei of the pairs ---------------------------------------------
  PairsTab pairsTab;
  PairsModel& pairs = *pairsTab.getPairsModel();
  addPair(pairs, 0, P_12C);
  {
    ok("one nucleus from p+12C", pairs.nuclei().size() == 1 && pairs.nuclei().first() == Nucleus{7, 13});
    ok("single nucleus: no abundance block", !pairsTab.hasAbundanceBlock());
  }
  addPair(pairs, 1, G_13N);
  ok("the capture pair to 13N joins the same nucleus", pairs.nuclei().size() == 1);
  addPair(pairs, 2, P_13C);
  QList<Nucleus> nuclei = pairs.nuclei();
  ok("p+13C adds 14N", nuclei.size() == 2 && nuclei.at(1) == Nucleus{7, 14});
  ok("html label", PairsModel::nucleusLabel(Nucleus{7, 13}) == "<sup>13</sup>N", PairsModel::nucleusLabel(Nucleus{7, 13}));
  ok("plain label (unicode superscript)", PairsModel::nucleusLabel(Nucleus{7, 13}, false) == QString::fromUtf8("\u00b9\u00b3N"),
     PairsModel::nucleusLabel(Nucleus{7, 13}, false));
  ok("long label names the first particle pair", pairs.nucleusLongLabel(Nucleus{7, 14}).startsWith(QString::fromUtf8("\u00b9\u2074N (")) &&
     pairs.nucleusLongLabel(Nucleus{7, 14}).contains(QString::fromUtf8("\u00b9\u00b3C+p")), pairs.nucleusLongLabel(Nucleus{7, 14}));
  ok("nucleusIndex", pairs.nucleusIndex(Nucleus{7, 14}) == 1 && pairs.nucleusIndex(Nucleus{8, 16}) == -1);

  // ---- 2. abundances --------------------------------------------------------
  AbundancesModel* abundances = pairsTab.getAbundancesModel();
  ok("abundance table has one row per nucleus", abundances->rowCount() == 2, QString::number(abundances->rowCount()));
  ok("two nuclei: the block is written even at the defaults", pairsTab.hasAbundanceBlock() && abundances->isDefault());
  {
    QString in("7 13 0.6 0 0\n7 14 0.4 1 0.05\nenr: 7 13 0.1 0 0\nenr: 7 14 0.9 1 0.02\nenr: sum 0.01\n</abundances>\n");
    QTextStream inStream(&in);
    ok("block read", pairsTab.readAbundances(inStream));
    ok("two samples", abundances->samples() == (QStringList() << QString() << "enr"), abundances->samples().join(","));
    AbundanceEntry e2 = abundances->entry("enr", Nucleus{7, 14});
    ok("the named sample's entry", e2.value == 0.9 && e2.vary && e2.sigma == 0.02 && e2.sample == "enr");
    ok("the named sample's own sum width", abundances->sampleSumSigma("enr") == 0.01 && abundances->sampleSumSigma(QString()) == -1.0);
    ok("parameter names", AbundancesModel::parameterName(e2) == "abundance_enr_7_14" &&
       AbundancesModel::parameterName(abundances->entry(Nucleus{7, 13})) == "abundance_7_13");
    abundances->setCurrentSample("enr");
    ok("the table shows the chosen sample", abundances->data(abundances->index(1, AbundancesModel::ValueCol), Qt::DisplayRole).toDouble() == 0.9);
    abundances->setCurrentSample(QString());
    AbundanceEntry a = abundances->entry(Nucleus{7, 13}), b = abundances->entry(Nucleus{7, 14});
    ok("entries", a.value == 0.6 && !a.vary && a.sigma == 0.0 && b.value == 0.4 && b.vary && b.sigma == 0.05,
       QString("%1 %2 %3 / %4 %5 %6").arg(a.value).arg(a.vary).arg(a.sigma).arg(b.value).arg(b.vary).arg(b.sigma));
    QString out;
    QTextStream outStream(&out);
    pairsTab.writeAbundances(outStream);
    QStringList written = tokens(out);
    for (QString& line : written) line = line.section('#', 0, 0).trimmed();  // the writer comments each line with the nucleus
    ok("block written back, sample by sample", written == (QStringList() << "7 13 0.6 0 0" << "7 14 0.4 1 0.05" << "enr: 7 13 0.1 0 0" << "enr: 7 14 0.9 1 0.02" << "enr: sum 0.01"), tokens(out).join(" | "));
    ok("a bad sample name is refused", !abundances->addSample("two words") && !abundances->addSample("a:b"));
    ok("sum sigma default off", pairsTab.abundanceSumSigma() == 0.0);
    pairsTab.setAbundanceSumSigma(0.01);
    ok("sum sigma set", pairsTab.abundanceSumSigma() == 0.01);
    ok("a negative value is refused", !abundances->setEntry(Nucleus{7, 13}, -1.0, false, 0.0));
  }

  // ---- 3. levels per nucleus ------------------------------------------------
  LevelsTab levelsTab;
  levelsTab.setPairsModel(&pairs);
  QObject::connect(&levelsTab, SIGNAL(readNewPair(PairsData, int, bool)), &pairsTab, SLOT(addPair(PairsData, int, bool)));
  QObject::connect(&levelsTab, SIGNAL(readExistingPair(PairsData, int, bool)), &pairsTab, SLOT(editPair(PairsData, int, bool)));
  {
    QString n13 = jpiList(levelsTab.spinParityOptions(Nucleus{7, 13}));
    QString n14 = jpiList(levelsTab.spinParityOptions(Nucleus{7, 14}));
    const QString plus = QString(QChar(0x207A)), minus = QString(QChar(0x207B));
    ok("13N options are half-integer (p+12C, 1/2+ on 0+)",
       n13.startsWith(QString("1/2%1 1/2%2 3/2%1 3/2%2 5/2%1 5/2%2").arg(plus).arg(minus)), n13);
    ok("14N options are integer (p+13C, 1/2+ on 1/2-)",
       n14.startsWith(QString("0%1 0%2 1%1 1%2 2%1 2%2").arg(plus).arg(minus)), n14);
    LevelsData half = {1, 0, 0.5, 1, 2.365, 7, 13};
    LevelsData one = {1, 0, 1.0, 1, 8.25, 7, 14};
    bool pure = true;
    for (const ChannelsData& c : levelsTab.calculateChannels(half, -1)) pure &= (c.pairIndex == 0 || c.pairIndex == 1);
    for (const ChannelsData& c : levelsTab.calculateChannels(one, -1)) pure &= (c.pairIndex == 2);
    ok("channels of a level come from its own nucleus's pairs", pure);
    ok("13N 1/2+ has its p and gamma channels", levelsTab.calculateChannels(half, -1).size() == 1 + 1,
       QString::number(levelsTab.calculateChannels(half, -1).size()));
    LevelsData anyHalf = {1, 0, 0.5, 1, 2.365};
    ok("a level without a nucleus takes channels from every pair (here only 13N's fit a 1/2+)",
       levelsTab.calculateChannels(anyHalf, -1).size() == 2, QString::number(levelsTab.calculateChannels(anyHalf, -1).size()));
  }
  {
    // 13N 1/2+ and 3/2- (pair 1), 14N 0+ and 1- (pair 3), all at 9.0 MeV, level IDs 1-4.
    const QString block =
        "0.5 1 9.0 0 1 1 1 0 1 1 0 0.1 0.5 1 0 1 0 1.00782 12.0 1 6 1.9435 1.9435 0 0 0.0 0 5.0 0 0 0\n"
        "0 1 9.0 0 1 3 1 2 2 1 0 0.2 0.5 1 0.5 -1 0 1.00782 13.00335 1 6 1.9435 7.5506 0 0 0.0 0 5.0 0 0 0\n"
        "1.5 -1 9.0 0 1 1 1 2 3 1 0 0.3 0.5 1 0 1 0 1.00782 12.0 1 6 1.9435 1.9435 0 0 0.0 0 5.0 0 0 0\n"
        "1 -1 9.0 0 1 3 2 0 4 1 0 0.4 0.5 1 0.5 -1 0 1.00782 13.00335 1 6 1.9435 7.5506 0 0 0.0 0 5.0 0 0 0\n"
        "</levels>\n";
    QString in(block);
    QTextStream inStream(&in);
    ok("two-nucleus levels block read", levelsTab.readNuclearFile(inStream));
    QList<LevelsData> levels = levelsTab.getLevelsModel()->getLevels();
    ok("four levels", levels.size() == 4, QString::number(levels.size()));
    int n13 = 0, n14 = 0;
    for (const LevelsData& l : levels) { if (l.nucleus() == Nucleus{7, 13}) n13++; if (l.nucleus() == Nucleus{7, 14}) n14++; }
    ok("levels carry the nucleus of their channels' pair", n13 == 2 && n14 == 2, QString("%1 %2").arg(n13).arg(n14));
    ok("the nucleus column shows the nucleus",
       levelsTab.getLevelsModel()->data(levelsTab.getLevelsModel()->index(1, LevelsModel::NucleusColumn), Qt::DisplayRole).toString() == "<sup>14</sup>N");
    QString out;
    QTextStream outStream(&out);
    levelsTab.writeNuclearFile(outStream);
    QStringList got = tokens(out);
    ok("written as one block per nucleus (13N levels first, then 14N)", got.size() == 4 &&
       got[0].startsWith("0.5 1 9") && got[1].startsWith("1.5 -1 9") && got[2].startsWith("0 1 9") && got[3].startsWith("1 -1 9"),
       got.join(" | "));
    // the same (J, pi, E) is another level in another nucleus
    LevelsData in14 = {1, 0, 1.0, -1, 9.0, 7, 14};
    LevelsData in13 = {1, 0, 1.0, -1, 9.0, 7, 13};
    LevelsData any = {1, 0, 1.0, -1, 9.0};
    ok("isLevel compares the nucleus (an unassigned one matches)", levelsTab.getLevelsModel()->isLevel(in14) == 3 &&
       levelsTab.getLevelsModel()->isLevel(in13) == -1 && levelsTab.getLevelsModel()->isLevel(any) == 3,
       QString("%1 %2 %3").arg(levelsTab.getLevelsModel()->isLevel(in14)).arg(levelsTab.getLevelsModel()->isLevel(in13))
           .arg(levelsTab.getLevelsModel()->isLevel(any)));
    QStringList missing = levelsTab.missingSpinParities();
    ok("missing J^pi are reported per nucleus", !missing.isEmpty() && missing.first().startsWith("<sup>13</sup>N "), missing.join(", "));
  }
  {
    // the level dialog: the nucleus box is shown with two nuclei and switches the J^pi list
    AddLevelDialog dialog;
    QStringList labels;
    QList<QList<SpinParityOption> > perNucleus;
    for (const Nucleus& n : nuclei) { labels << pairs.nucleusLongLabel(n); perNucleus << levelsTab.spinParityOptions(n); }
    dialog.setSpinParityOptions(levelsTab.spinParityOptions());
    dialog.setNucleusOptions(labels, perNucleus, 0, true);
    QComboBox* nucleusCombo = dialog.findChild<QComboBox*>("nucleusCombo");
    QComboBox* jpiCombo = dialog.findChild<QComboBox*>("spinParityCombo");
    ok("level dialog lists both nuclei", nucleusCombo && nucleusCombo->count() == 2);
    if (nucleusCombo && jpiCombo) {
      int n13 = jpiCombo->count();
      nucleusCombo->setCurrentIndex(1);
      ok("choosing 14N switches the J^pi list", jpiCombo->count() == perNucleus[1].size() && jpiCombo->count() != n13,
         QString("%1 -> %2").arg(n13).arg(jpiCombo->count()));
      ok("nucleusIndex follows the box", dialog.nucleusIndex() == 1);
    }
  }

  // ---- 4. natural segments ------------------------------------------------
  ok("isNatural", SegmentKeywords::isNatural("natural 1") && !SegmentKeywords::isNatural("natural 0") &&
     !SegmentKeywords::isNatural("") && SegmentKeywords::isNatural("sqrtshift 0.1 0.2 1 natural 1"));
  ok("withNatural adds the block last", SegmentKeywords::withNatural("sqrtshift 0.1 0.2 1", true) == "sqrtshift 0.1 0.2 1 natural 1" &&
     SegmentKeywords::withNatural("", true) == "natural 1");
  ok("withNatural(false) removes it and keeps the rest", SegmentKeywords::withNatural("natural 1 sqrtshift 0.1 0.2 1", false) == "sqrtshift 0.1 0.2 1" &&
     SegmentKeywords::withNatural("natural 1", false).isEmpty());
  {
    AddSegDataDialog dialog;
    dialog.setPairsModel(&pairs);
    dialog.setSamples(QStringList() << QString() << "enr");
    ok("data dialog: natural box shown with two nuclei", !dialog.naturalCombo->isHidden() && dialog.naturalCombo->isEnabled());
    ok("data dialog lists No, the default sample and enr", dialog.naturalCombo->count() == 3, QString::number(dialog.naturalCombo->count()));
    dialog.setNaturalToken("enr");
    ok("token round trip", dialog.naturalToken() == "enr", dialog.naturalToken());
    dialog.setNaturalToken("1");
    dialog.advancedModeCheck->setChecked(true);
    dialog.operationCombo->setCurrentIndex(1);  // ratio
    ok("ratio clears and disables it", !dialog.naturalCombo->isEnabled() && dialog.naturalToken().isEmpty());
    dialog.operationCombo->setCurrentIndex(0);  // sum
    ok("sum allows it again", dialog.naturalCombo->isEnabled());
    dialog.setDataTypeCode(2);  // phase shift
    ok("phase shift refuses it", !dialog.naturalCombo->isEnabled());
    dialog.setDataTypeCode(0);
    ok("total cross section allows it", dialog.naturalAllowed());
    PairsModel one;
    addPair(one, 0, P_12C);
    AddSegDataDialog single;
    single.setPairsModel(&one);
    ok("data dialog: hidden with one nucleus", single.naturalCombo->isHidden());
    AddSegTestDialog test;
    test.setPairsModel(&pairs);
    ok("test dialog: shown with two nuclei", !test.naturalCombo->isHidden());
    test.setDataTypeCode(3);  // angular-distribution coefficients
    ok("test dialog: refused for angular-distribution coefficients", !test.naturalCombo->isEnabled());
  }
  SegmentsTab segmentsTab;
  segmentsTab.setPairsModel(&pairs);
  {
    // (the second line differs in its energy range: the tab refuses a duplicate line with a message box)
    QString in("1 1 1 1.0 3.0 0 0 0 1 0 0 0 0 0 data.dat 1 0 1 3 3 -999 natural 1\n"
               "1 1 1 1.5 2.5 0 0 0 1 0 0 0 0 0 plain.dat 0 0\n</segmentsData>\n");
    QTextStream inStream(&in);
    ok("data segments read", segmentsTab.readSegDataFile(inStream));
    SegmentsDataModel* model = segmentsTab.getSegmentsDataModel();
    QString natLabel = model->data(model->index(0, 1), Qt::DisplayRole).toString();
    QString plainLabel = model->data(model->index(1, 1), Qt::DisplayRole).toString();
    ok("natural row is flagged 'nat' in its label", natLabel.contains("<i>nat</i>"), natLabel);
    QString in2("1 2 2 1.0 2.0 0 0 0 1 0 0 0 0 0 enr.dat 0 0 natural enr\n</segmentsData>\n");
    QTextStream inStream2(&in2);
    segmentsTab.readSegDataFile(inStream2);
    QString sampleLabel = model->data(model->index(2, 1), Qt::DisplayRole).toString();
    ok("a named-sample row is flagged nat:<sample>", sampleLabel.contains("<i>nat:enr</i>"), sampleLabel);
    ok("SegmentKeywords sample token", SegmentKeywords::naturalSample("sqrtshift 0.1 0.2 1 natural enr") == "enr" &&
       SegmentKeywords::naturalSample("natural 1") == "1" && SegmentKeywords::naturalSample("natural 0").isEmpty() &&
       SegmentKeywords::withNaturalToken("natural 1", "enr") == "natural enr");
    ok("plain row is not", !plainLabel.contains("nat"), plainLabel);
    QString out;
    QTextStream outStream(&out);
    outStream.setFieldAlignment(QTextStream::AlignLeft);
    segmentsTab.writeSegDataFile(outStream);
    ok("natural written last on its line", tokens(out).first().endsWith(" natural 1"), tokens(out).first());
  }

  // ---- 5. fitting parameters -----------------------------------------------
  {
    FittingTab fitting;
    fitting.setTabReferences(&levelsTab, &segmentsTab);
    fitting.setAbundancesModel(abundances);
    fitting.populateFromCurrentGUIState();
    const QList<FittingParameter>& params = fitting.getFittingParameters();
    int firstAbundance = -1, lastLevel = -1, firstNorm = -1, nAbundance = 0;
    for (int i = 0; i < params.size(); i++) {
      if (params[i].category == "abundance") { nAbundance++; if (firstAbundance < 0) firstAbundance = i; }
      if (params[i].category == "level") lastLevel = i;
      if (params[i].category == "norm" && firstNorm < 0) firstNorm = i;
    }
    ok("one abundance parameter per nucleus and sample", nAbundance == 4, QString::number(nAbundance));
    ok("abundances sit after the levels and before the normalizations",
       firstAbundance > lastLevel && (firstNorm < 0 || firstAbundance < firstNorm),
       QString("%1 %2 %3").arg(lastLevel).arg(firstAbundance).arg(firstNorm));
    if (nAbundance == 4) {
      const FittingParameter& a = params[firstAbundance];
      const FittingParameter& b = params[firstAbundance + 1];
      ok("named as in param.sav, samples as blocks", a.name == "abundance_7_13" && b.name == "abundance_7_14" &&
         params[firstAbundance + 2].name == "abundance_enr_7_13" && params[firstAbundance + 3].name == "abundance_enr_7_14",
         a.name + " " + b.name + " " + params[firstAbundance + 2].name);
      ok("value, sigma and vary from the table", a.value == 0.6 && !a.useAsNuisance && b.value == 0.4 && b.useAsNuisance && b.error == 0.05);
      ok("only the fitted abundance takes a Minuit slot", a.minuitIndex == -1 && b.minuitIndex >= 0,
         QString("%1 %2").arg(a.minuitIndex).arg(b.minuitIndex));
    }
    FittingTab plain;
    PairsTab onePairTab;
    addPair(*onePairTab.getPairsModel(), 0, P_12C);
    plain.setTabReferences(&levelsTab, &segmentsTab);
    plain.setAbundancesModel(onePairTab.getAbundancesModel());
    plain.populateFromCurrentGUIState();
    int n = 0;
    for (const FittingParameter& p : plain.getFittingParameters()) if (p.category == "abundance") n++;
    ok("no abundance parameter without a block", n == 0, QString::number(n));
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
