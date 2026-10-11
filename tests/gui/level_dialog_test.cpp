// Headless checks on the Add/Edit Level dialog's J^pi list.
//
// The dialog offers only the J^pi values that give a level at least one channel
// under the current pairs and limits (LevelsTab::spinParityOptions), falls back
// to every J^pi up to 10 when no pairs are defined, and keeps an edited level's
// own J^pi selectable even when it gives no channels.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.
// With AZURE2_LEVEL_DIALOG_PNG=<stem> it also saves the dialogs as
// <stem>_add.png, <stem>_edit.png and <stem>_nopairs.png for a layout review.

#include <QApplication>
#include <QComboBox>
#include <QPushButton>
#include <QString>
#include <cstdlib>
#include <iostream>
#include "LevelsTab.h"
#include "AddLevelDialog.h"
#include "PairsModel.h"
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

// Expected labels written with ASCII signs; the dialog shows them as superscripts.
static QString SP(const char* ascii) {
  return QString::fromUtf8(ascii).replace("+", QString(QChar(0x207A))).replace("-", QString(QChar(0x207B)));
}

static QString labels(const QList<SpinParityOption>& options) {
  QStringList out;
  for (const SpinParityOption& o : options) out << AddLevelDialog::spinParityLabel(o.jValue, o.piValue);
  return out.join(" ");
}

static QString comboLabels(AddLevelDialog& dialog) {
  QComboBox* combo = dialog.findChild<QComboBox*>();
  QStringList out;
  for (int i = 0; i < combo->count(); i++) out << combo->itemText(i);
  return out.join("|");
}

static void save(QWidget& w, const char* suffix) {
  const char* stem = std::getenv("AZURE2_LEVEL_DIALOG_PNG");
  if (!stem) return;
  w.show();
  QApplication::processEvents();
  w.grab().save(QString("%1_%2.png").arg(stem).arg(suffix));
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  PairsModel pairs;
  LevelsTab tab;
  tab.setPairsModel(&pairs);

  // 1. No pairs: nothing gives channels, so the dialog falls back to every J^pi up to 10.
  ok("no pairs -> no options from the tab", tab.spinParityOptions().isEmpty());
  {
    AddLevelDialog dialog;
    dialog.setSpinParityOptions(tab.spinParityOptions());
    QComboBox* combo = dialog.findChild<QComboBox*>();
    ok("fallback lists 0+ .. 10- (42 entries)", combo->count() == 42, QString::number(combo->count()));
    save(dialog, "nopairs");
  }

  // 2. n + 14N (1/2+ on 1+), maximum l = 2 (the tab's default).
  pairs.insertRows(0, 1);
  const QVariant n14N[PairsData::SIZE] = {0.5, 1, 0, 1.0087, 0., 1.0, 1, 7, 14.0031, 0., 0., 10.833, 5.0, 0, 0};
  for (int c = 0; c < PairsData::SIZE; c++) pairs.setData(pairs.index(0, c), n14N[c]);
  QList<SpinParityOption> options = tab.spinParityOptions();
  // s = 1/2, 3/2; l = 0 (+), 1 (-), 2 (+)
  ok("n+14N, l<=2 -> 1/2+ 1/2- 3/2+ 3/2- 5/2+ 5/2- 7/2+",
     labels(options) == SP("1/2+ 1/2- 3/2+ 3/2- 5/2+ 5/2- 7/2+"), labels(options));
  bool counted = true;
  for (const SpinParityOption& o : options) {
    LevelsData level = {1, 0, o.jValue, o.piValue, 0.};
    counted &= (o.numChannels == tab.calculateChannels(level, -1).size());
  }
  ok("channel counts match calculateChannels", counted);
  bool integerJ = false;
  for (const SpinParityOption& o : options) integerJ |= (int(2 * o.jValue) % 2 == 0);
  ok("no integer J offered for n+14N", !integerJ);

  // 2b. A capture pair to the 15N 1/2- ground state (E1/M1, maximum multipolarity 2) must not open integer J:
  // the triangle condition alone accepted J = 1 for a 1/2 final state (11B+a project, 2026-10-11).
  pairs.insertRows(1, 1);
  const QVariant g15N[PairsData::SIZE] = {1.0, -1, 0, 0., 0., 0.5, -1, 7, 15.0001, 0., 0., 0., 0., 10, 0};
  for (int c = 0; c < PairsData::SIZE; c++) pairs.setData(pairs.index(1, c), g15N[c]);
  options = tab.spinParityOptions();
  ok("capture to 1/2- adds only half-integer J (3/2+ 3/2- 5/2+ 5/2- 1/2+ 1/2- 7/2+ ...)",
     labels(options) == SP("1/2+ 1/2- 3/2+ 3/2- 5/2+ 5/2- 7/2+"), labels(options));
  {
    LevelsData oneMinus = {1, 0, 1.0, -1, 0.};
    ok("J = 1- gets no gamma channel to 1/2-", tab.calculateChannels(oneMinus, -1).isEmpty());
    LevelsData threeHalvesPlus = {1, 0, 1.5, 1, 0.};
    // 3 particle channels (s=1/2 l=2; s=3/2 l=0,2) + E1 and M2 to 1/2-
    ok("J = 3/2+ gets its E1 and M2 channels", tab.calculateChannels(threeHalvesPlus, -1).size() == 3 + 2,
       QString::number(tab.calculateChannels(threeHalvesPlus, -1).size()));
  }
  pairs.removeRows(1, 1);
  options = tab.spinParityOptions();

  // 3. Add: the list is the tab's, Accept waits for a numerical energy.
  {
    AddLevelDialog dialog;
    dialog.setSpinParityOptions(options);
    QPushButton* accept = nullptr;
    for (QPushButton* b : dialog.findChildren<QPushButton*>()) if (b->text() == "Accept") accept = b;
    ok("Accept disabled with no energy", accept && !accept->isEnabled());
    dialog.energyText->setText("abc");
    ok("Accept disabled for a non-number", accept && !accept->isEnabled());
    dialog.energyText->setText("-0.25");
    ok("Accept enabled for a number (bound levels too)", accept && accept->isEnabled());
    ok("first option selected", dialog.jValue() == 0.5 && dialog.piValue() == 1);
    dialog.energyText->setText("10.4");
    save(dialog, "add");
  }

  // 4. Edit: an existing J^pi in the list is selected as it is.
  {
    AddLevelDialog dialog;
    dialog.setSpinParityOptions(options, 2.5, -1, true);
    ok("edit selects 5/2-", dialog.jValue() == 2.5 && dialog.piValue() == -1);
    ok("no extra entry when the J^pi is in the list", dialog.findChild<QComboBox*>()->count() == options.size());
  }

  // 5. Edit: a J^pi with no channels (9/2-, read from a file) stays selectable, in order.
  {
    AddLevelDialog dialog;
    dialog.setWindowTitle("Edit a Level");
    dialog.setSpinParityOptions(options, 4.5, -1, true);
    dialog.energyText->setText("11.2");
    ok("edit keeps 9/2- selected", dialog.jValue() == 4.5 && dialog.piValue() == -1, comboLabels(dialog));
    ok("9/2- added after 7/2+ and marked", comboLabels(dialog).endsWith(SP("7/2+|9/2-  (no channels)")), comboLabels(dialog));
    save(dialog, "edit");
  }
  {
    AddLevelDialog dialog;
    dialog.setSpinParityOptions(options, 1.0, 1, true);
    ok("1+ inserted between 1/2- and 3/2+", comboLabels(dialog).startsWith(SP("1/2+|1/2-|1+  (no channels)|3/2+")),
       comboLabels(dialog));
  }

  std::cout << (fails ? "FAILED" : "PASSED") << " (" << fails << " failures)" << std::endl;
  return fails ? 1 : 0;
}
