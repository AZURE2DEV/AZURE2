// Headless checks on the observable controls of the segment dialogs.
//
// The observable list is short (Angle Integrated, Differential, Phase Shift,
// Analyzing Power, Angular Coefficients, and for test segments Angular Distribution Coefficients);
// its variants are separate controls that map onto the unchanged isDiff codes:
//   - Differential + Frame Lab / C.M.            -> data 1 / 4,  test 1 / 5
//   - Angle Integrated + gamma exit pair +
//     Multipolarity All / E1 only / E2 only      -> data 0 / 5 / 6 (no test codes)
//   - Angle Integrated + exit item "Total capture (all gamma pairs)"
//                                                -> data 3, test 4 (exit key -1)
// Every isDiff code a dialog can write must come back from an edit unchanged,
// with the same exit key.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.
// With AZURE2_SEG_OBSERVABLE_PNG=<stem> it saves the data dialog as
// <stem>_capture.png (gamma exit, E2 only), <stem>_total.png (total capture)
// and <stem>_differential.png (C.M. differential).
//
// The angle box is titled after the frame of the angles: "C.M. Angle" for a
// C.M. differential cross section and an analyzing power, "Lab Angle" otherwise
// (it used to say "Lab" for every observable).

#include <QApplication>
#include <QComboBox>
#include <QStandardItemModel>
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
static void ok(const QString& what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what.toStdString();
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

static void addPair(PairsModel& pairs, int row, const QVariant (&v)[PairsData::SIZE]) {
  pairs.insertRows(row, 1);
  for (int c = 0; c < PairsData::SIZE; c++) pairs.setData(pairs.index(row, c), v[c]);
}

// 11B+a (particle pair, key 1) and 15N+gamma (pair type 10, key 2)
static const QVariant P_ALPHA[PairsData::SIZE] = {0.0, 1, 2, 4.0026, 0., 1.5, -1, 5, 11.0093, 0., 0., 10.9912, 5.5, 0, 0};
static const QVariant P_GAMMA[PairsData::SIZE] = {1.0, -1, 0, 0., 0., 0.5, -1, 7, 15.0001, 0., 0., 0., 0., 10, 0};

static const QString LAB = "Lab Angle [degrees]", CM = "C.M. Angle [degrees]";

static bool itemEnabled(QComboBox* combo, int code) {
  QStandardItemModel* items = qobject_cast<QStandardItemModel*>(combo->model());
  int i = combo->findData(code);
  return items && i >= 0 && items->item(i)->isEnabled();
}

// What the Segments tab does on an edit: exit key first, then the code; on accept the
// exit key is -1 for total capture, the combo's key otherwise.
template <class Dialog>
static void load(Dialog& d, int exitKey, int code) {
  d.entrancePairCombo->setKey(1);
  d.exitPairCombo->setKey(exitKey);
  d.setDataTypeCode(code);
}
template <class Dialog>
static int writtenExit(Dialog& d, int totalCode) {
  return d.dataTypeCode() == totalCode ? -1 : d.exitPairCombo->key();
}

static void shot(QWidget& w, const char* suffix) {
  if (const char* stem = std::getenv("AZURE2_SEG_OBSERVABLE_PNG")) {
    w.show();
    w.adjustSize();
    QApplication::processEvents();
    w.grab().save(QString("%1_%2.png").arg(stem).arg(suffix));
    w.hide();
  }
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  PairsModel pairs;
  addPair(pairs, 0, P_ALPHA);
  addPair(pairs, 1, P_GAMMA);
  PairsModel noGamma;
  addPair(noGamma, 0, P_ALPHA);

  // ---------------- data dialog ----------------
  {
    AddSegDataDialog d;
    d.setPairsModel(&pairs);
    ok("data: observable list is Angle Integrated, Differential, Phase Shift, Analyzing Power, Angular Coefficients, Scattering Length",
       d.dataTypeCombo->count() == 6 && d.dataTypeCombo->itemData(0).toInt() == 0 && d.dataTypeCombo->findData(1) >= 0 &&
           d.dataTypeCombo->findData(2) >= 0 && d.dataTypeCombo->findData(7) >= 0 && d.dataTypeCombo->findData(9) == 4 &&
           d.dataTypeCombo->findData(AddSegDataDialog::ScatteringLengthCode) == 5,
       QString::number(d.dataTypeCombo->count()));
    ok("data: no total-capture, C.M. or E1/E2 items in the observable list",
       d.dataTypeCombo->findData(3) < 0 && d.dataTypeCombo->findData(4) < 0 && d.dataTypeCombo->findData(5) < 0 && d.dataTypeCombo->findData(6) < 0);
    ok("data: total-capture exit item present with a gamma pair", d.exitPairCombo->findData(PairKeyCombo::TotalCapture) >= 0);
    ok("data: total-capture item reads 'Total capture (all γ pairs)'",
       d.exitPairCombo->itemText(d.exitPairCombo->findData(PairKeyCombo::TotalCapture)) == QString::fromUtf8("Total capture (all γ pairs)"),
       d.exitPairCombo->itemText(d.exitPairCombo->findData(PairKeyCombo::TotalCapture)));
    ok("data: component exit combo has no total-capture item", d.componentExitCombo->findData(PairKeyCombo::TotalCapture) < 0);

    // particle exit: no multipolarity
    d.exitPairCombo->setKey(1);
    d.dataTypeCombo->setCurrentIndex(d.dataTypeCombo->findData(0));
    ok("data: particle exit, Angle Integrated -> multipolarity hidden, isDiff 0",
       d.multipolarityCombo->isHidden() && d.multipolarityLabel->isHidden() && d.dataTypeCode() == 0, QString::number(d.dataTypeCode()));
    ok("data: frame hidden for Angle Integrated", d.frameCombo->isHidden() && d.frameLabel->isHidden());
    ok("data: Angle Integrated -> angle box titled Lab", d.angleBox->title() == LAB, d.angleBox->title());

    // gamma exit: multipolarity appears and maps to 0/5/6
    d.exitPairCombo->setKey(2);
    ok("data: gamma exit, Angle Integrated -> multipolarity shown", !d.multipolarityCombo->isHidden() && !d.multipolarityLabel->isHidden());
    const int mult[3] = {0, 5, 6};
    for (int m : mult) {
      d.multipolarityCombo->setCurrentIndex(d.multipolarityCombo->findData(m));
      ok(QString("data: gamma exit, multipolarity %1 -> isDiff %1").arg(m), d.dataTypeCode() == m, QString::number(d.dataTypeCode()));
    }
    shot(d, "capture");
    // multipolarity hidden (and read as All) for other observables and for a particle exit
    d.dataTypeCombo->setCurrentIndex(d.dataTypeCombo->findData(2));
    ok("data: gamma exit, Phase Shift -> multipolarity hidden, isDiff 2", d.multipolarityCombo->isHidden() && d.dataTypeCode() == 2,
       QString::number(d.dataTypeCode()));
    d.dataTypeCombo->setCurrentIndex(d.dataTypeCombo->findData(0));
    d.exitPairCombo->setKey(1);
    ok("data: E2 chosen then particle exit -> hidden, read as All (isDiff 0)", d.multipolarityCombo->isHidden() && d.dataTypeCode() == 0,
       QString::number(d.dataTypeCode()));
    d.multipolarityCombo->setCurrentIndex(d.multipolarityCombo->findData(0));

    // Differential: frame drop-down
    d.dataTypeCombo->setCurrentIndex(d.dataTypeCombo->findData(1));
    ok("data: Differential -> frame shown, Lab by default, isDiff 1",
       !d.frameCombo->isHidden() && d.dataTypeCode() == 1 && d.lowAngleText->isEnabled(), QString::number(d.dataTypeCode()));
    ok("data: Differential + Lab -> angle box titled Lab", d.angleBox->title() == LAB, d.angleBox->title());
    d.frameCombo->setCurrentIndex(d.frameCombo->findData(4));
    ok("data: Differential + C.M. -> angle box titled C.M.", d.angleBox->title() == CM, d.angleBox->title());
    ok("data: Differential + C.M. -> isDiff 4, angles enabled", d.dataTypeCode() == 4 && d.lowAngleText->isEnabled(), QString::number(d.dataTypeCode()));
    shot(d, "differential");
    d.dataTypeCombo->setCurrentIndex(d.dataTypeCombo->findData(7));
    ok("data: Analyzing Power -> frame hidden, isDiff 7, angles enabled, vary norm disabled",
       d.frameCombo->isHidden() && d.dataTypeCode() == 7 && d.lowAngleText->isEnabled() && !d.varyNormCheck->isEnabled());
    ok("data: Analyzing Power -> angle box titled C.M.", d.angleBox->title() == CM, d.angleBox->title());
    d.dataTypeCombo->setCurrentIndex(d.dataTypeCombo->findData(2));
    ok("data: Phase Shift -> J and l shown", !d.phaseJValueText->isHidden() && !d.phaseLValueText->isHidden() && !d.lowAngleText->isEnabled());
    ok("data: Phase Shift -> angle box titled Lab", d.angleBox->title() == LAB, d.angleBox->title());

    // total capture: restricts the observable list to Angle Integrated
    d.dataTypeCombo->setCurrentIndex(d.dataTypeCombo->findData(1));
    d.exitPairCombo->setKey(PairKeyCombo::TotalCapture);
    ok("data: total capture switches the observable to Angle Integrated, isDiff 3",
       d.observableCode() == 0 && d.dataTypeCode() == 3, QString::number(d.dataTypeCode()));
    ok("data: total capture disables the other observables",
       itemEnabled(d.dataTypeCombo, 0) && !itemEnabled(d.dataTypeCombo, 1) && !itemEnabled(d.dataTypeCombo, 2) && !itemEnabled(d.dataTypeCombo, 7));
    ok("data: total capture hides multipolarity and frame", d.multipolarityCombo->isHidden() && d.frameCombo->isHidden());
    ok("data: total capture writes exit key -1", writtenExit(d, 3) == -1);
    d.exitPairCombo->refresh();  // what dropping the list does
    ok("data: total capture survives a refresh of the exit list", d.exitPairCombo->key() == PairKeyCombo::TotalCapture && d.dataTypeCode() == 3);
    d.entrancePairCombo->setKey(2);
    ok("data: total capture survives an entrance change", d.exitPairCombo->key() == PairKeyCombo::TotalCapture);
    d.entrancePairCombo->setKey(1);  // 11B+alpha: the picture shows 11B(a,g) total capture
    ok("data: total capture -> angle box titled Lab", d.angleBox->title() == LAB, d.angleBox->title());
    shot(d, "total");
    d.exitPairCombo->setKey(2);
    ok("data: leaving total capture re-enables the observables", itemEnabled(d.dataTypeCombo, 1) && itemEnabled(d.dataTypeCombo, 7));
  }
  {
    AddSegDataDialog d;
    d.setPairsModel(&noGamma);
    ok("data: no total-capture exit item without a gamma pair", d.exitPairCombo->findData(PairKeyCombo::TotalCapture) < 0,
       QString::number(d.exitPairCombo->count()));
    d.exitPairCombo->refresh();
    ok("data: ... nor after a refresh", d.exitPairCombo->findData(PairKeyCombo::TotalCapture) < 0);
  }
  // edit round trips: (exit key, isDiff) -> controls -> (exit key, isDiff)
  {
    struct Case { int exit; int code; int observable; const char* what; };
    const Case cases[] = {
      {1, 0, 0, "angle integrated, particle exit"}, {2, 0, 0, "angle integrated, gamma exit"},
      {1, 1, 1, "lab differential"},                {-1, 3, 0, "total capture"},
      {1, 4, 1, "C.M. differential"},               {2, 5, 0, "E1 only"},
      {2, 6, 0, "E2 only"},                         {1, 7, 7, "analyzing power"},
      {1, 2, 2, "phase shift"},                     {1, 9, 9, "a code the dialog does not list (9)"},
    };
    for (const Case& c : cases) {
      AddSegDataDialog d;
      d.setPairsModel(&pairs);
      load(d, c.exit, c.code);
      const QString tag = QString("data edit isDiff %1 (%2)").arg(c.code).arg(c.what);
      ok(tag + ": observable restored", d.observableCode() == c.observable, QString::number(d.observableCode()));
      ok(tag + ": isDiff and exit key written back", d.dataTypeCode() == c.code && writtenExit(d, 3) == c.exit,
         QString("%1 / %2").arg(d.dataTypeCode()).arg(writtenExit(d, 3)));
      if (c.code == 1 || c.code == 4)
        ok(tag + ": frame restored", !d.frameCombo->isHidden() && d.frameCombo->currentData().toInt() == c.code);
      ok(tag + ": angle box title", d.angleBox->title() == ((c.code == 4 || c.code == 7) ? CM : LAB), d.angleBox->title());
      if (c.code == 5 || c.code == 6)
        ok(tag + ": multipolarity restored", !d.multipolarityCombo->isHidden() && d.multipolarityCombo->currentData().toInt() == c.code);
      if (c.code == 3)
        ok(tag + ": total-capture exit item selected", d.exitPairCombo->key() == PairKeyCombo::TotalCapture);
    }
    // a total-capture line in a project without a gamma pair keeps its item
    AddSegDataDialog d;
    d.setPairsModel(&noGamma);
    load(d, -1, 3);
    ok("data edit isDiff 3 without a gamma pair: item kept, isDiff 3, exit -1",
       d.exitPairCombo->key() == PairKeyCombo::TotalCapture && d.dataTypeCode() == 3 && writtenExit(d, 3) == -1);
  }

  // ---------------- test dialog ----------------
  {
    AddSegTestDialog d;
    d.setPairsModel(&pairs);
    ok("test: observable list is Angle Integrated, Differential, Phase Shift, Ang. Dist. Coefficients, Analyzing Power",
       d.dataTypeCombo->count() == 5 && d.dataTypeCombo->itemData(0).toInt() == 0 && d.dataTypeCombo->findData(3) >= 0 &&
           d.dataTypeCombo->findData(7) >= 0 && d.dataTypeCombo->findData(4) < 0 && d.dataTypeCombo->findData(5) < 0);
    ok("test: total-capture exit item present with a gamma pair", d.exitPairCombo->findData(PairKeyCombo::TotalCapture) >= 0);
    d.exitPairCombo->setKey(2);
    ok("test: gamma exit, Angle Integrated -> isDiff 0, frame hidden", d.dataTypeCode() == 0 && d.frameCombo->isHidden());
    d.dataTypeCombo->setCurrentIndex(d.dataTypeCombo->findData(1));
    ok("test: Differential -> frame shown, Lab, isDiff 1", !d.frameCombo->isHidden() && d.dataTypeCode() == 1 && d.angleStepText->isEnabled());
    d.frameCombo->setCurrentIndex(d.frameCombo->findData(5));
    ok("test: Differential + C.M. -> isDiff 5", d.dataTypeCode() == 5, QString::number(d.dataTypeCode()));
    ok("test: Differential + C.M. -> angle box titled C.M.", d.angleBox->title() == CM, d.angleBox->title());
    d.dataTypeCombo->setCurrentIndex(d.dataTypeCombo->findData(3));
    ok("test: Angular Distribution Coefficients -> order shown, isDiff 3", !d.angDistSpin->isHidden() && d.dataTypeCode() == 3);
    d.exitPairCombo->setKey(PairKeyCombo::TotalCapture);
    ok("test: total capture switches to Angle Integrated, isDiff 4", d.observableCode() == 0 && d.dataTypeCode() == 4,
       QString::number(d.dataTypeCode()));
    ok("test: total capture disables the other observables", !itemEnabled(d.dataTypeCombo, 1) && !itemEnabled(d.dataTypeCombo, 3));
    ok("test: total capture writes exit key -1", writtenExit(d, 4) == -1);
  }
  {
    AddSegTestDialog d;
    d.setPairsModel(&noGamma);
    ok("test: no total-capture exit item without a gamma pair", d.exitPairCombo->findData(PairKeyCombo::TotalCapture) < 0);
  }
  {
    struct Case { int exit; int code; int observable; };
    const Case cases[] = {{1, 0, 0}, {2, 0, 0}, {1, 1, 1}, {1, 2, 2}, {1, 3, 3}, {-1, 4, 0}, {1, 5, 1}, {1, 7, 7}};
    for (const Case& c : cases) {
      AddSegTestDialog d;
      d.setPairsModel(&pairs);
      load(d, c.exit, c.code);
      const QString tag = QString("test edit isDiff %1").arg(c.code);
      ok(tag + ": observable restored", d.observableCode() == c.observable, QString::number(d.observableCode()));
      ok(tag + ": isDiff and exit key written back", d.dataTypeCode() == c.code && writtenExit(d, 4) == c.exit,
         QString("%1 / %2").arg(d.dataTypeCode()).arg(writtenExit(d, 4)));
      if (c.code == 1 || c.code == 5)
        ok(tag + ": frame restored", !d.frameCombo->isHidden() && d.frameCombo->currentData().toInt() == c.code);
      ok(tag + ": angle box title", d.angleBox->title() == ((c.code == 5 || c.code == 7) ? CM : LAB), d.angleBox->title());
    }
  }

  std::cout << (fails ? "FAILED" : "PASSED") << " (" << fails << " failures)" << std::endl;
  return fails ? 1 : 0;
}
