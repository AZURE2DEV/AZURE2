// Headless round-trip checks on the Experimental Effects tab's optional
// energy-range tokens.
//
// A targetInt line may end with a quoted lab-energy range list, a blend width
// and an automatic-application tolerance.  The GUI has to read them, carry
// them in its model, and write them back -- but only when they are used, so a
// file that never touched the feature stays byte-identical.  Both directions
// are checked here through the tab's own readFile/writeFile.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QTextStream>
#include <QString>
#include <iostream>
#include "TargetIntTab.h"
#include "TargetIntModel.h"
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

static QString roundTrip(TargetIntTab& tab, const QString& block, bool& readOk) {
  tab.reset();
  QString in(block);
  QTextStream inStream(&in);
  readOk = tab.readFile(inStream);
  QString out;
  QTextStream outStream(&out);
  tab.writeFile(outStream);
  return out;
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  TargetIntTab tab;

  // 1. A line with the new tokens survives a read/write round trip.
  bool readOk = false;
  QString withTokens =
      "1 \"1,2\" 50 1 0.030 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 5 50 "
      "\"1.95-2.55,2.8-3.0\" 0.12 0.002\n</targetInt>\n";
  QString out = roundTrip(tab, withTokens, readOk);
  ok("read a line carrying ranges/blend/tolerance", readOk);
  ok("ranges token written back", out.contains("\"1.95-2.55,2.8-3.0\""), out);
  ok("blend width written back", out.contains(" 0.12 "), out);
  ok("tolerance written back", out.contains(" 0.002"), out);

  QList<TargetIntData> lines = tab.getTargetIntModel()->getLines();
  ok("model row present", lines.size() == 1);
  if (lines.size() == 1) {
    ok("model carries ranges", lines.at(0).applyRanges == "1.95-2.55,2.8-3.0",
       lines.at(0).applyRanges);
    ok("model carries blend width", lines.at(0).transitionWidth == 0.12);
    ok("model carries tolerance", lines.at(0).autoTolerance == 0.002);
  }

  // 2. Reading the GUI's own output again reproduces the same fields.
  bool readOk2 = false;
  QString out2 = roundTrip(tab, out + "</targetInt>\n", readOk2);
  ok("re-read the written line", readOk2);
  ok("second round trip stable", out2 == out, out2);

  // 3. A legacy line without the tokens gains nothing on write: files that
  //    never used the feature stay byte-identical in this respect.
  QString legacy = "1 \"1,2\" 50 1 0.030 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 5 50\n</targetInt>\n";
  QString out3 = roundTrip(tab, legacy, readOk);
  ok("read a legacy line", readOk);
  // Exactly the six quotes of the three legacy string fields; a written
  // ranges token would add two more.
  ok("no ranges token invented", out3.count('"') == 6, out3);
  ok("line ends at points-per-width", out3.trimmed().endsWith("50"), out3);

  // 4. Tolerance-only form: empty ranges token, zero width, tolerance set.
  QString tolOnly = "1 \"1\" 50 1 0.030 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 5 50 \"\" 0 0.005\n</targetInt>\n";
  QString out4 = roundTrip(tab, tolOnly, readOk);
  ok("read a tolerance-only line", readOk);
  lines = tab.getTargetIntModel()->getLines();
  ok("empty ranges stay empty", lines.size() == 1 && lines.at(0).applyRanges.isEmpty());
  ok("tolerance-only round trip keeps tolerance", out4.contains(" 0.005"), out4);

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
