// Headless round-trip checks on the Experimental Effects tab's optional
// trailing tokens.
//
// A targetInt line may end with a quoted lab-energy range list, a blend width
// and an automatic-application tolerance, and after those with a beam-profile
// block introduced by the `beamprofile` keyword.  The GUI has to read them,
// carry them in its model, and write them back -- but only when they are used,
// so a file that never touched either feature stays byte-identical.  Both
// directions are checked here through the tab's own readFile/writeFile.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QTextStream>
#include <QString>
#include <cmath>
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

  // 5. A beam-profile block with no ranges tokens before it: the keyword is
  //    the first non-numeric thing on the line after the grid parameters.
  //    This is the shape pyazr writes, so it is the one that matters most.
  QString beamOnly =
      "1 \"7\" 150 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50 "
      "beamprofile 1 3.7802596656384173 0.3760611 -2.06 1 0.07334525 0 1\n</targetInt>\n";
  QString out5 = roundTrip(tab, beamOnly, readOk);
  ok("read a beam-profile line", readOk);
  lines = tab.getTargetIntModel()->getLines();
  ok("beam-profile row present", lines.size() == 1, out5);
  if (lines.size() == 1) {
    const TargetIntData& d = lines.at(0);
    ok("model flags the beam profile", d.isBeamProfile);
    ok("one component, four numbers", d.beamProfile.size() == 4,
       QString::number(d.beamProfile.size()));
    if (d.beamProfile.size() == 4) {
      // 12 significant digits survive the write; the engine's own parser reads
      // the same text, so anything coarser would move the profile on a save.
      ok("location kept to 1e-9", std::fabs(d.beamProfile.at(0) / 3.7802596656384173 - 1.) < 1e-9,
         QString::number(d.beamProfile.at(0), 'g', 17));
      ok("scale kept", d.beamProfile.at(1) == 0.3760611);
      ok("skewness kept", d.beamProfile.at(2) == -2.06);
      ok("weight kept", d.beamProfile.at(3) == 1.);
    }
    ok("resolution sigma kept", d.beamTpcSigma == 0.07334525);
    ok("truncation kept", d.beamTruncation == 0.);
    ok("detailed-balance flag kept", d.beamPhotodissociation);
    // The grid parameters sit just before the keyword and must not be eaten.
    ok("width multiplier still read", d.resonanceWidthMultiplier == 20.);
    ok("points per width still read", d.pointsPerWidth == 50.);
    ok("no ranges invented", d.applyRanges.isEmpty());
  }
  ok("keyword written back", out5.contains("beamprofile 1 "), out5);
  ok("skewness written back", out5.contains("-2.06"), out5);
  ok("no ranges token invented alongside", out5.count('"') == 6, out5);
  bool readOk5 = false;
  QString out5b = roundTrip(tab, out5 + "</targetInt>\n", readOk5);
  ok("re-read the written beam-profile line", readOk5);
  ok("beam-profile round trip stable", out5b == out5, out5b);

  // 6. Several components, and ranges tokens in front of the keyword: both
  //    optional blocks on one line, in the order the engine parses them.
  QString beamMulti =
      "1 \"3-5\" 200 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 10 25 \"1.9-2.6\" 0.1 0.001 "
      "beamprofile 2 2.5 0.3 -1.5 0.4 2.9 0.31 -1.7 0.6 0.055 2 0\n</targetInt>\n";
  QString out6 = roundTrip(tab, beamMulti, readOk);
  ok("read a two-component line with ranges", readOk);
  lines = tab.getTargetIntModel()->getLines();
  if (lines.size() == 1) {
    const TargetIntData& d = lines.at(0);
    ok("two components, eight numbers", d.beamProfile.size() == 8,
       QString::number(d.beamProfile.size()));
    ok("second component read", d.beamProfile.size() == 8 && d.beamProfile.at(4) == 2.9
                                   && d.beamProfile.at(7) == 0.6);
    ok("truncation read", d.beamTruncation == 2.);
    ok("detailed balance off", !d.beamPhotodissociation);
    ok("ranges still read alongside", d.applyRanges == "1.9-2.6", d.applyRanges);
    ok("blend width still read", d.transitionWidth == 0.1);
  }
  ok("both blocks written back", out6.contains("\"1.9-2.6\"") && out6.contains("beamprofile 2 "), out6);
  bool readOk6 = false;
  QString out6b = roundTrip(tab, out6 + "</targetInt>\n", readOk6);
  ok("two-component round trip stable", readOk6 && out6b == out6, out6b);

  // 7. A block that declares more components than it supplies is dropped
  //    rather than half-applied: the engine would otherwise fold the data
  //    with a profile the tab never showed.
  QString beamShort =
      "1 \"1\" 50 0 0 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 5 50 "
      "beamprofile 2 2.5 0.3 -1.5 0.4\n</targetInt>\n";
  QString out7 = roundTrip(tab, beamShort, readOk);
  lines = tab.getTargetIntModel()->getLines();
  ok("truncated block leaves the profile off",
     lines.size() == 1 && !lines.at(0).isBeamProfile && lines.at(0).beamProfile.isEmpty());
  ok("truncated block is not written back", !out7.contains("beamprofile"), out7);

  // 8. A legacy line still gains no beam-profile token (case 3 checks the
  //    ranges half of the same promise).
  QString out8 = roundTrip(tab, legacy, readOk);
  ok("no beam-profile token invented", !out8.contains("beamprofile"), out8);

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
