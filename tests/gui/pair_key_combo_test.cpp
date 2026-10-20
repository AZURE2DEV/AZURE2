// Headless checks on the particle-pair drop-downs that replaced the typed pair
// keys (PairKeyCombo) in the segment dialogs and the reaction-rate box.
//
// The .azr still holds 1-based keys; the box shows the pairs by name and hands
// the key back. Entrance items read "2: 14N+n [0.000 MeV]", exit items read as
// the reaction with the chosen entrance pair, "2: 11B(a,n)14N [0.000 MeV]".
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.
// With AZURE2_PAIR_COMBO_PNG=<stem> it saves the data-segment dialog as
// <stem>_data.png for a layout review.

#include <QApplication>
#include <QComboBox>
#include <QString>
#include <cstdlib>
#include <iostream>
#include "PairKeyCombo.h"
#include "PairsModel.h"
#include "AddSegDataDialog.h"
#include "AddSegTestDialog.h"
#include "Config.h"
struct SegPairs {int firstPair; int secondPair;};

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

static void addPair(PairsModel& pairs, int row, const QVariant (&v)[PairsData::SIZE]) {
  pairs.insertRows(row, 1);
  for (int c = 0; c < PairsData::SIZE; c++) pairs.setData(pairs.index(row, c), v[c]);
}

// Unicode the labels use
static const QString a = QString(QChar(0x03B1)), g = QString(QChar(0x03B3));
static QString sup(const char* digits) {
  static const QString s = QString::fromUtf8("⁰¹²³⁴⁵⁶⁷⁸⁹");
  QString out;
  for (const char* c = digits; *c; c++) out += s.at(*c - '0');
  return out;
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  PairsModel pairs;
  // 11B+a, 14N+n, 14C+p, 15N+gamma (the 11B+a project's pairs)
  const QVariant p1[PairsData::SIZE] = {0.0, 1, 2, 4.0026, 0., 1.5, -1, 5, 11.0093, 0., 0., 10.9912, 5.5, 0, 0};
  const QVariant p2[PairsData::SIZE] = {0.5, 1, 0, 1.00866, 0., 1.0, 1, 7, 14.0031, 0., 0., 10.8333, 4.0, 0, 0};
  const QVariant p3[PairsData::SIZE] = {0.5, 1, 1, 1.00782, 0., 0.0, 1, 6, 14.0032, 0., 0., 10.2075, 4.0, 0, 0};
  const QVariant p4[PairsData::SIZE] = {1.0, -1, 0, 0., 0., 0.5, -1, 7, 15.0001, 0., 0., 0., 0., 10, 0};
  addPair(pairs, 0, p1); addPair(pairs, 1, p2); addPair(pairs, 2, p3); addPair(pairs, 3, p4);

  // 1. plain-text conversion of the model's HTML labels
  ok("plainText: 11B+a", PairsModel::plainText("<center><sup>11</sup>B+&alpha; [0.000 MeV]</center>") == sup("11") + "B+" + a + " [0.000 MeV]",
     PairsModel::plainText("<center><sup>11</sup>B+&alpha; [0.000 MeV]</center>"));
  ok("plainText: beta-", PairsModel::plainText("&beta;<sup>-</sup>") == QString(QChar(0x03B2)) + QChar(0x207B));
  ok("getPairLabel: 14N+n", pairs.getPairLabel(pairs.getPairs().at(1)) == sup("14") + "N+n [0.000 MeV]", pairs.getPairLabel(pairs.getPairs().at(1)));

  // 2. an entrance box lists the pairs by name, in key order, and returns the key
  PairKeyCombo entrance(PairKeyCombo::Entrance);
  entrance.setPairsModel(&pairs);
  ok("entrance: 4 items", entrance.count() == 4, QString::number(entrance.count()));
  ok("entrance: nothing selected until a key is set", entrance.key() == 0);
  ok("entrance: item 1 reads '1: 11B+a [0.000 MeV]'", entrance.itemText(0) == "1: " + sup("11") + "B+" + a + " [0.000 MeV]", entrance.itemText(0));
  ok("entrance: item 4 reads '4: 15N+g [0.000 MeV]'", entrance.itemText(3) == "4: " + sup("15") + "N+" + g + " [0.000 MeV]", entrance.itemText(3));
  entrance.setKey(3);
  ok("entrance: setKey(3) -> key() 3", entrance.key() == 3 && entrance.currentIndex() == 2);

  // 3. an exit box reads as reactions with the entrance pair, and follows it
  PairKeyCombo exit(PairKeyCombo::Exit);
  exit.setPairsModel(&pairs);
  QObject::connect(&entrance, SIGNAL(keyChanged(int)), &exit, SLOT(setEntranceKey(int)));
  entrance.setKey(1);
  exit.setKey(2);
  ok("exit: with entrance 1 item 2 reads '2: 11B(a,n)14N [0.000 MeV]'",
     exit.itemText(1) == "2: " + sup("11") + "B(" + a + ",n)" + sup("14") + "N [0.000 MeV]", exit.itemText(1));
  ok("exit: capture item reads '4: 11B(a,g)15N [0.000 MeV]'",
     exit.itemText(3) == "4: " + sup("11") + "B(" + a + "," + g + ")" + sup("15") + "N [0.000 MeV]", exit.itemText(3));
  entrance.setKey(2);
  ok("exit: relabelled when the entrance changes, selection kept",
     exit.key() == 2 && exit.itemText(1) == "2: " + sup("14") + "N(n,n)" + sup("14") + "N [0.000 MeV]", exit.itemText(1));

  // 4. a key beyond the pairs list (inconsistent file) stays selectable and marked; 0 clears
  exit.setKey(7);
  ok("exit: undefined key 7 kept and marked", exit.key() == 7 && exit.currentText().contains("no such pair"), exit.currentText());
  exit.setKey(0);
  ok("exit: setKey(0) clears", exit.key() == 0 && exit.currentIndex() == -1);

  // 5. a pair added after the box was filled shows up on refresh (what dropping the list does)
  const QVariant p5[PairsData::SIZE] = {0.5, 1, 1, 1.00782, 0., 1.0, -1, 6, 14.0032, 0., 6.0938, 10.2075, 4.0, 0, 0};
  addPair(pairs, 4, p5);
  entrance.refresh();
  ok("entrance: 5 items after a pair is added, key kept", entrance.count() == 5 && entrance.key() == 2);
  ok("entrance: excited residual carries its energy", entrance.itemText(4).endsWith("[6.094 MeV]"), entrance.itemText(4));

  // 6. the dialogs round-trip the keys the Segments tab hands them
  {
    AddSegDataDialog dialog;
    dialog.setPairsModel(&pairs);
    ok("data dialog: fresh dialog starts on pair 1/1", dialog.entrancePairCombo->key() == 1 && dialog.exitPairCombo->key() == 1);
    dialog.entrancePairCombo->setKey(1);
    dialog.exitPairCombo->setKey(2);
    ok("data dialog: keys read back", dialog.entrancePairCombo->key() == 1 && dialog.exitPairCombo->key() == 2);
    ok("data dialog: exit shows the reaction", dialog.exitPairCombo->currentText().contains("B(" + a + ",n)"), dialog.exitPairCombo->currentText());
    dialog.lowEnergyText->setText("0.5"); dialog.highEnergyText->setText("3.0");
    if (const char* stem = std::getenv("AZURE2_PAIR_COMBO_PNG")) {
      dialog.show(); QApplication::processEvents();
      dialog.grab().save(QString("%1_data.png").arg(stem));
    }
  }
  {
    AddSegTestDialog dialog;
    dialog.setPairsModel(&pairs);
    dialog.entrancePairCombo->setKey(2);
    dialog.exitPairCombo->setKey(4);
    ok("test dialog: keys read back", dialog.entrancePairCombo->key() == 2 && dialog.exitPairCombo->key() == 4);
    ok("test dialog: 14N(n,g)15N", dialog.exitPairCombo->currentText() == "4: " + sup("14") + "N(n," + g + ")" + sup("15") + "N [0.000 MeV]", dialog.exitPairCombo->currentText());
  }

  std::cout << (fails ? "FAILED" : "PASSED") << " (" << fails << " failures)" << std::endl;
  return fails ? 1 : 0;
}
