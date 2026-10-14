// Headless round trip of <segmentsData> and <segmentsTest> lines through the
// Segments tab's own readers and writers.
//
// A segment line may end with keyword blocks the GUI does not edit --
// "sqrtshift b bError vary" (sqrt(E) energy-shift term) and "mixture w_self
// n {file w}..." (natural-abundance targets), "shiftgroup <label>" and
// "shiftfn \"<f(x)>\" N v s vary ..." (energy-shift function) -- which it has to carry through
// a load/save verbatim.  The sqrtshift block used to be lost on every load:
// the reader split it off, but the row was then rebuilt column by column and
// the field had no column.  Both blocks are now stored on the row directly.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QTextStream>
#include <QString>
#include <QStringList>
#include <iostream>
#include "SegmentsTab.h"
#include "SegmentsDataModel.h"
#include "SegmentsTestModel.h"
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

// Whitespace-normalised tokens of a block, line by line.
static QStringList tokens(const QString& text) {
  QStringList out;
  for (const QString& line : text.split('\n', Qt::SkipEmptyParts)) {
    QString t = line.simplified();
    if (!t.isEmpty()) out << t;
  }
  return out;
}

static void addPair(PairsModel& pairs, int row, const QVariant (&v)[PairsData::SIZE]) {
  pairs.insertRows(row, 1);
  for (int c = 0; c < PairsData::SIZE; c++) pairs.setData(pairs.index(row, c), v[c]);
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  // The tab labels its rows through the pairs model, so it needs one (17O+a's n+20Ne and a+17O pairs).
  PairsModel pairs;
  const QVariant p1[PairsData::SIZE] = {0.0, 1, 2, 4.0026, 0., 2.5, 1, 8, 16.9991, 0., 0., 7.347, 1.0, 0, 0};
  const QVariant p2[PairsData::SIZE] = {0.5, 1, 0, 1.00866, 0., 0.0, 1, 10, 19.9924, 0., 0., 7.347, 1.0, 0, 0};
  const QVariant p3[PairsData::SIZE] = {0.5, 1, 0, 1.00866, 0., 2.0, 1, 10, 19.9924, 1.63365, 0., 7.347, 1.0, 0, 0};
  addPair(pairs, 0, p1); addPair(pairs, 1, p2); addPair(pairs, 2, p3);
  SegmentsTab tab;
  tab.setPairsModel(&pairs);

  // ---- data segments -------------------------------------------------------
  const QString dataBlock =
      "1 1 1 0.5 3 0 180 1 1 0 5 0.004 0.002 0 data/xs_40.dat 0 0 sqrtshift 0.003 0.001 0\n"
      "1 2 2 0 2.5 0 180 0 1 0 0 0 0 0 data/Junghans_ntotal.dat 1 0 2 2 1 -999 2 3 -999 0 mixture 0.9048 2 \"tables/22Ne.dat\" 0.0925 \"tables/21Ne.dat\" 0.0027\n"
      "1 1 1 1 3 0 180 4 1 0 1 0 0 0 data/plain.dat 0 0\n"
      "0 1 2 0 5 0 180 0 1 1 7.2 0 0 0 data/both.dat 0 0 sqrtshift 0.001 0.002 1 mixture 0.5 1 \"t.dat\" 0.5\n"
      "1 1 1 0.5 3 0 180 1 1 0 5 0 0.01 1 data/g.dat 0 0 shiftgroup G shiftfn \"a0*(sqrt(x)-sqrt(1.7))\" 1 0 0.02 1\n"
      "</segmentsData>\n";
  {
    QString in(dataBlock);
    QTextStream inStream(&in);
    ok("data block read", tab.readSegDataFile(inStream));
    QList<SegmentsDataData> lines = tab.getSegmentsDataModel()->getLines();
    ok("five data lines", lines.size() == 5, QString::number(lines.size()));
    if (lines.size() == 5) {
      ok("shiftgroup and shiftfn blocks kept on the row",
         lines[4].extraTokens == "shiftgroup G shiftfn \"a0*(sqrt(x)-sqrt(1.7))\" 1 0 0.02 1", lines[4].extraTokens);
      ok("data file of the shift-function row", lines[4].dataFile == "data/g.dat", lines[4].dataFile);
      ok("sqrtshift block kept on the row", lines[0].extraTokens == "sqrtshift 0.003 0.001 0", lines[0].extraTokens);
      ok("mixture block kept on a composite row",
         lines[1].extraTokens == "mixture 0.9048 2 \"tables/22Ne.dat\" 0.0925 \"tables/21Ne.dat\" 0.0027", lines[1].extraTokens);
      ok("composite components still parsed", lines[1].isAdvanced == 1 && lines[1].componentsList.split(";").size() == 2,
         lines[1].componentsList);
      ok("plain row has no extras", lines[2].extraTokens.isEmpty(), lines[2].extraTokens);
      ok("both blocks kept together", lines[3].extraTokens == "sqrtshift 0.001 0.002 1 mixture 0.5 1 \"t.dat\" 0.5",
         lines[3].extraTokens);
      ok("data file of the mixture row", lines[1].dataFile == "data/Junghans_ntotal.dat", lines[1].dataFile);
    }
    QString out;
    QTextStream outStream(&out);
    outStream.setFieldAlignment(QTextStream::AlignLeft);  // as AZURESetup::writeFile sets it
    tab.writeSegDataFile(outStream);
    QStringList got = tokens(out), want = tokens(dataBlock);
    want.removeLast();  // the closing tag is written by the caller
    ok("data lines written back token for token", got == want, got.join(" | "));
  }

  // ---- test segments -------------------------------------------------------
  const QString testBlock =
      "1 1 1 1 3 0.002 120 120 1 1 mixture 0.5 1 \"table.dat\" 0.5\n"
      "1 2 2 0 2.5 0.01 0 180 0 0 1 0 2 2 1 -999 2 3 -999 mixture 0.9048 1 \"tables/22Ne.dat\" 0.0925\n"
      "1 1 1 1.5 3 0.002 120 120 1 1 0\n"
      "</segmentsTest>\n";
  {
    QString in(testBlock);
    QTextStream inStream(&in);
    ok("test block read", tab.readSegTestFile(inStream));
    QList<SegmentsTestData> lines = tab.getSegmentsTestModel()->getLines();
    ok("three test lines", lines.size() == 3, QString::number(lines.size()));
    if (lines.size() == 3) {
      ok("mixture after a plain test line (no isAdvanced token)", lines[0].extraTokens == "mixture 0.5 1 \"table.dat\" 0.5"
                                                                   && lines[0].isAdvanced == 0, lines[0].extraTokens);
      ok("mixture after a composite test line", lines[1].extraTokens == "mixture 0.9048 1 \"tables/22Ne.dat\" 0.0925"
                                                && lines[1].isAdvanced == 1, lines[1].extraTokens);
      ok("plain test line has no extras", lines[2].extraTokens.isEmpty(), lines[2].extraTokens);
    }
    QString out;
    QTextStream outStream(&out);
    outStream.setFieldAlignment(QTextStream::AlignLeft);
    tab.writeSegTestFile(outStream);
    QStringList got = tokens(out);
    // The writer always emits the isAdvanced token, so the first line gains a "0" before its block.
    QStringList want;
    want << "1 1 1 1 3 0.002 120 120 1 1 0 mixture 0.5 1 \"table.dat\" 0.5"
         << "1 2 2 0 2.5 0.01 0 180 0 0 1 0 2 2 1 -999 2 3 -999 mixture 0.9048 1 \"tables/22Ne.dat\" 0.0925"
         << "1 1 1 1.5 3 0.002 120 120 1 1 0";
    ok("test lines written back with their blocks last", got == want, got.join(" | "));
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
