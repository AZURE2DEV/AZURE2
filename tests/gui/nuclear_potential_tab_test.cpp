// Headless checks on the Nuclear Potential tab's per-pair selector.
//
// The tab is constructed before a project exists: a project's particle pairs
// are built while <levels> is read, which happens *after* the <potential>
// block. A selector filled once at construction therefore stays empty forever,
// which is exactly the bug this test exists to stop coming back -- the tab
// shipped once showing nothing but "Default" and no way to reach a channel.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QComboBox>
#include <QCheckBox>
#include <QLineEdit>
#include <QPushButton>
#include <QTextStream>
#include <QString>
#include <iostream>
#include "NuclearPotentialTab.h"
#include "NuclearPotentialManager.h"
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

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  NuclearPotentialManager& mgr = NuclearPotentialManager::instance();
  mgr.resetToDefault();

  PairsModel pairs;
  NuclearPotentialTab tab;
  tab.setPairsModel(&pairs);

  QComboBox* combo = tab.findChild<QComboBox*>();
  ok("selector exists", combo != nullptr);
  if(!combo) return 1;

  std::cout << "\n1. before any pairs exist (this is what the user saw)\n";
  ok("only the default entry", combo->count() == 1,
     QString("count=%1").arg(combo->count()));

  std::cout << "\n2. pairs arrive later, as they do when <levels> is read\n";
  pairs.insertRows(0, 3, QModelIndex());
  ok("selector picked them up", combo->count() == 4,
     QString("count=%1").arg(combo->count()));
  for(int i = 0; i < combo->count(); i++)
    std::cout << "        [" << i << "] " << combo->itemText(i).toStdString()
              << "  data=" << combo->itemData(i).toInt() << std::endl;
  ok("keys are 1-based and positional",
     combo->itemData(0).toInt() == 0 && combo->itemData(1).toInt() == 1 &&
     combo->itemData(3).toInt() == 3);
  ok("selector is enabled", combo->isEnabled());

  std::cout << "\n3. select pair 2 and switch it on\n";
  combo->setCurrentIndex(combo->findData(2));
  QCheckBox* enabled = tab.findChild<QCheckBox*>();
  ok("checkbox exists", enabled != nullptr);
  enabled->setChecked(true);
  ok("pair 2 on in the manager", mgr.isPairEnabled(2));
  ok("pair 2 has its own setting", mgr.hasPairSetting(2));
  ok("pair 1 untouched", !mgr.isPairEnabled(1) && !mgr.hasPairSetting(1));
  ok("label shows its state",
     combo->itemText(combo->findData(2)).contains("[on]"),
     combo->itemText(combo->findData(2)));

  std::cout << "\n4. give it a depth and commit\n";
  QList<QLineEdit*> edits = tab.findChildren<QLineEdit*>();
  edits.at(0)->setText("33.5");                      // Woods-Saxon V0
  QList<QPushButton*> buttons = tab.findChildren<QPushButton*>();
  QPushButton* apply = nullptr;
  for(QPushButton* b : buttons) if(b->text().contains("Apply")) apply = b;
  ok("apply button exists", apply != nullptr);
  apply->click();
  ok("V0 stored for pair 2", qFuzzyCompare(mgr.getSetting(2).V0, 33.5),
     QString::number(mgr.getSetting(2).V0));
  ok("default still untouched", qFuzzyCompare(mgr.getDefaultSetting().V0, 150.0),
     QString::number(mgr.getDefaultSetting().V0));

  std::cout << "\n5. switching pairs keeps each one's own values\n";
  combo->setCurrentIndex(combo->findData(3));
  ok("pair 3 inherits, not its own", !mgr.hasPairSetting(3));
  combo->setCurrentIndex(combo->findData(1));
  ok("just visiting pair 3 did not pin a setting to it", !mgr.hasPairSetting(3));
  ok("just visiting pair 1 did not pin a setting to it", !mgr.hasPairSetting(1));
  combo->setCurrentIndex(combo->findData(2));
  ok("pair 2 still 33.5 after the round trip",
     qFuzzyCompare(mgr.getSetting(2).V0, 33.5), QString::number(mgr.getSetting(2).V0));
  ok("its widget shows 33.5", edits.at(0)->text().toDouble() == 33.5,
     edits.at(0)->text());

  std::cout << "\n6. deleting pair 1 shifts pair 2's setting down to key 1\n";
  tab.onPairRemoved(1);
  pairs.removeRows(0, 1, QModelIndex());
  ok("setting moved to key 1", mgr.hasPairSetting(1) &&
     qFuzzyCompare(mgr.getSetting(1).V0, 33.5));
  ok("nothing left at key 2", !mgr.hasPairSetting(2));
  ok("selector shrank", combo->count() == 3,
     QString("count=%1").arg(combo->count()));

  std::cout << "\n7. the <potential> block round-trips\n";
  QString buf;
  {
    QTextStream out(&buf);
    out << "useHybridPotential=" << (mgr.getDefaultEnabled() ? 1 : 0) << "\n";
    tab.writePotentialSettings(out);
    out << "</potential>\n";
  }
  std::cout << "        --- written ---\n";
  for(const QString& l : buf.split("\n")) if(!l.isEmpty()) std::cout << "        " << l.toStdString() << "\n";
  mgr.resetToDefault();
  Config cfg(std::cout);
  {
    QTextStream in(&buf);
    tab.readPotentialSettings(in, cfg);
  }
  ok("pair 1 came back on", mgr.isPairEnabled(1));
  ok("pair 1 V0 came back", qFuzzyCompare(mgr.getSetting(1).V0, 33.5),
     QString::number(mgr.getSetting(1).V0));

  std::cout << "\n8. a file with no pair= line still reads as it always did\n";
  QString legacy = "useHybridPotential=1\npotentialType=0\nV0=77\nR=3.1\na=0.5\n</potential>\n";
  mgr.resetToDefault();
  {
    QTextStream in(&legacy);
    tab.readPotentialSettings(in, cfg);
  }
  ok("default carries it", qFuzzyCompare(mgr.getDefaultSetting().V0, 77.0) &&
     mgr.getDefaultEnabled(), QString::number(mgr.getDefaultSetting().V0));
  ok("no pair has its own", mgr.configuredPairs().empty());
  ok("every pair inherits it", mgr.isPairEnabled(1) && mgr.isPairEnabled(2));

  std::cout << std::endl;
  if(fails) { std::cout << "FAILED: " << fails << " check(s)" << std::endl; return 1; }
  std::cout << "all GUI checks passed" << std::endl;
  return 0;
}
