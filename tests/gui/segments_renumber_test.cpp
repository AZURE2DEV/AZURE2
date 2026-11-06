// Headless check that what names a data segment by its line number follows the
// segment when lines are moved, added or deleted in the Segments tab.
//
// The <thm> block's experiment[...] segments= and weight[k], the Experimental
// Effects segment lists and the segment_N_norm / segment_N_energy_shift rows
// of <parameterSettings> (settings and prior_centre) all count lines of
// <segmentsData>.  Before, a move or a delete left them on the old numbers,
// so a prior centre, a nuisance setting or a THM experiment silently passed
// to another data set.  Checked on examples/f19_pag_thm with a richer <thm>
// block, through AZURESetup's own open and save: each edit is saved, the
// saved file is reopened and saved again byte for byte, and an unedited
// project saves as it was.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QDir>
#include <QFile>
#include <QStatusBar>
#include <QString>
#include <QStringList>
#include <QTemporaryDir>
#include <iostream>
#include "AZURESetup.h"
#include "Config.h"
#include "SegmentsTab.h"
#include "ThmSettings.h"
struct SegPairs {int firstPair; int secondPair;};

// Defined by AZURE2.cpp, which belongs to the executable rather than the GUI
// library, so this test supplies its own. They are never called from here.
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

#ifndef AZURE2_SOURCE_DIR
#error "AZURE2_SOURCE_DIR must be defined"
#endif

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

// LF text (a Windows checkout gives the projects CRLF endings).
static QString slurp(const QString& path) {
  QFile f(path);
  if(!f.open(QIODevice::ReadOnly)) return QString();
  return QString::fromUtf8(f.readAll()).remove('\r');
}
static void spit(const QString& path, const QString& text) {
  QFile f(path);
  if(f.open(QIODevice::WriteOnly | QIODevice::Truncate)) f.write(text.toUtf8());
}

// The lines of section <name>, trimmed, blank ones dropped.
static QStringList section(const QString& text, const QString& name) {
  QStringList rows;
  bool in = false;
  for(const QString& raw : text.split('\n')) {
    const QString line = raw.trimmed();
    if(line == "<" + name + ">") { in = true; continue; }
    if(line == "</" + name + ">") { in = false; continue; }
    if(in && !line.isEmpty()) rows << line;
  }
  return rows;
}
// The data file of each <segmentsData> line, in order.
static QStringList dataFiles(const QString& text) {
  QStringList files;
  for(const QString& line : section(text, "segmentsData"))
    files << line.section(QRegExp("\\s+"), 14, 14);
  return files;
}
// The <parameterSettings> row of `name` with its settings (value .. nuisance), or its prior centre.
static QString settingOf(const QString& text, const QString& name) {
  for(const QString& row : section(text, "parameterSettings"))
    if(row.startsWith(name + " ") && !row.contains("prior_centre")) return row.section(' ', 1, 6);
  return QString();
}
static QString centreOf(const QString& text, const QString& name) {
  for(const QString& row : section(text, "parameterSettings"))
    if(row.startsWith(name + " prior_centre ")) return row.section(' ', 2, 2);
  return QString();
}
// The Experimental Effects segment lists.
static QStringList effectLists(const QString& text) {
  QStringList lists;
  for(const QString& line : section(text, "targetInt")) lists << line.section(QRegExp("\\s+"), 1, 1);
  return lists;
}

// Saves, then reopens and saves the saved file: the two must be the same.
static QString saveAndCheckRoundTrip(AZURESetup& w, const QString& path, const char* what) {
  w.saveProject();
  const QString saved = slurp(path);
  const QString again = path + ".again.azr";
  QFile::remove(again);
  QFile::copy(path, again);
  w.open(again);
  w.saveProject();
  ok(what, slurp(again) == saved);
  w.open(path);
  return saved;
}

int main(int argc, char** argv) {
  QTemporaryDir settingsDir;
  qputenv("XDG_CONFIG_HOME", settingsDir.path().toUtf8());
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("AZURE2-tests");
  QCoreApplication::setApplicationName("segments_renumber_test");

  QTemporaryDir work;
  const QString example = QString(AZURE2_SOURCE_DIR) + "/examples/f19_pag_thm";
  QString text = slurp(example + "/f19_pag_thm.azr");
  ok("found examples/f19_pag_thm", text.contains("experiment[SU2025] segments=1"));
  QDir(work.path()).mkpath("data");
  for(const QString& f : QDir(example + "/data").entryList(QDir::Files))
    QFile::copy(example + "/data/" + f, work.filePath("data/" + f));
  QDir(work.path()).mkpath("output");
  QDir(work.path()).mkpath("checks");
  spit(work.filePath("w1.dat"), "0.1 1\n0.5 1\n");
  spit(work.filePath("w3.dat"), "0.1 2\n0.5 2\n");

  // Segments: 1 thm_su2025 (THM), 2 juna_zhang2022, 3 spyrou2000_S.  More
  // references than the example has: weights on 1 and 3, a second
  // experiment on 2 and 3, an effect on 1 and 3, settings and a distinct
  // prior centre on segment 3.
  text.replace("experiment[SU2025] segments=1 background=linear",
               "# the THM run\nweight[1]=w1.dat\nweight[3]=w3.dat\n"
               "experiment[SU2025] segments=1 background=linear");
  text.replace("Ebeam=55\n</thm>", "Ebeam=55\nexperiment[DIRECT] segments=2,3   # kept with its comment\n</thm>");
  text.replace("1  \"1\"  40", "1  \"1,3\"  40");
  text.replace("segment_3_norm prior_centre 1", "segment_3_norm prior_centre 1.05\n"
               "segment_3_norm 1.092788 0.5 1.5 0.1 0 1 norm 0");
  spit(work.filePath("p0.in"), text);

  AZURESetup w;
  const QString path = work.filePath("project.azr");

  std::cout << "\n1. unedited: open + save twice gives the same file\n";
  QFile::copy(work.filePath("p0.in"), path);
  w.open(path);
  ThmSettings thm;
  ok("the <thm> block parses", w.thmSettings(thm) && thm.weight.size() == 2,
     QString::number(thm.weight.size()));
  const QString first = saveAndCheckRoundTrip(w, path, "an unedited project saves byte for byte");
  ok("first save keeps the references",
     first.contains("experiment[DIRECT] segments=2,3   # kept with its comment") &&
     first.contains("weight[3]=w3.dat") && effectLists(first) == QStringList() << "\"1,3\"" &&
     centreOf(first, "segment_3_norm") == "1.05" && settingOf(first, "segment_3_norm").startsWith("1.092788 0.5 1.5 "),
     settingOf(first, "segment_3_norm"));

  std::cout << "\n2. move segment 3 up (spyrou2000_S becomes segment 2)\n";
  ok("move accepted", w.getSegmentsTab()->moveDataSegment(2, 1));
  const QString moved = saveAndCheckRoundTrip(w, path, "reopened and saved: the same file");
  ok("lines moved", dataFiles(moved) == QStringList() << "data/thm_su2025.dat" << "data/spyrou2000_S.dat"
                                                      << "data/juna_zhang2022.dat",
     dataFiles(moved).join(","));
  ok("weight[3] is weight[2], weight[1] stays",
     moved.contains("weight[2]=w3.dat") && moved.contains("weight[1]=w1.dat") && !moved.contains("weight[3]"));
  ok("experiment lines unchanged (same sets)",
     moved.contains("experiment[DIRECT] segments=2,3   # kept with its comment") &&
     moved.contains("experiment[SU2025] segments=1 background=linear"));
  ok("effect list 1,3 is 1,2", effectLists(moved) == QStringList() << "\"1,2\"", effectLists(moved).join(" "));
  ok("prior centres follow their segments",
     centreOf(moved, "segment_2_norm") == "1.05" && centreOf(moved, "segment_3_norm") == "1",
     centreOf(moved, "segment_2_norm") + " " + centreOf(moved, "segment_3_norm"));
  ok("norm settings follow their segment",
     settingOf(moved, "segment_2_norm").startsWith("1.092788 0.5 1.5 ") &&
     settingOf(moved, "segment_3_norm").startsWith("0.904768 0 0 "),
     settingOf(moved, "segment_2_norm") + " | " + settingOf(moved, "segment_3_norm"));
  ok("comments of the block kept", moved.contains("# the THM run"));

  std::cout << "\n3. add a segment (appended: nothing renumbered)\n";
  SegmentsDataData extra = w.getSegmentsTab()->getSegmentsDataModel()->getLines().at(1);
  extra.dataFile = "data/juna_zhang2022.dat";
  extra.lowEnergy = 0.06;
  w.getSegmentsTab()->addSegDataLine(extra);
  const QString added = saveAndCheckRoundTrip(w, path, "reopened and saved: the same file");
  ok("four segments", dataFiles(added).size() == 4, dataFiles(added).join(","));
  ok("references as before",
     section(added, "thm") == section(moved, "thm") && effectLists(added) == effectLists(moved) &&
     centreOf(added, "segment_2_norm") == "1.05" && settingOf(added, "segment_2_norm").startsWith("1.092788 0.5 1.5 "));

  std::cout << "\n4. move the new segment to the top\n";
  for(int row = 3; row > 0; row--) w.getSegmentsTab()->moveDataSegment(row, row - 1);
  const QString top = saveAndCheckRoundTrip(w, path, "reopened and saved: the same file");
  ok("experiments follow", top.contains("experiment[SU2025] segments=2 ") &&
                           top.contains("experiment[DIRECT] segments=3,4   # kept with its comment"));
  ok("weights follow", top.contains("weight[2]=w1.dat") && top.contains("weight[3]=w3.dat") &&
                       !top.contains("weight[1]"));
  ok("effect list follows", effectLists(top) == QStringList() << "\"2,3\"", effectLists(top).join(" "));
  ok("prior centre and settings follow",
     centreOf(top, "segment_3_norm") == "1.05" && settingOf(top, "segment_3_norm").startsWith("1.092788 0.5 1.5 "));
  // Back to the order of step 3.
  w.getSegmentsTab()->moveDataSegment(0, 3);

  std::cout << "\n5. delete segment 1 (the only segment of SU2025)\n";
  w.statusBar()->clearMessage();
  w.getSegmentsTab()->deleteDataSegment(0);
  const QString notice = w.statusBar()->currentMessage();
  const QString deleted = saveAndCheckRoundTrip(w, path, "reopened and saved: the same file");
  ok("experiment SU2025 removed, and the status bar says so",
     !deleted.contains("experiment[SU2025]") && notice.contains("SU2025"), notice);
  ok("DIRECT renumbered", deleted.contains("experiment[DIRECT] segments=1,2   # kept with its comment"));
  ok("weight of the deleted segment dropped, the other one renumbered",
     deleted.contains("weight[1]=w3.dat") && !deleted.contains("w1.dat") && !deleted.contains("weight[2]"));
  ok("effect list 1,2 is 1", effectLists(deleted) == QStringList() << "\"1\"", effectLists(deleted).join(" "));
  ok("rows of the deleted segment dropped, the others renumbered",
     centreOf(deleted, "segment_1_norm") == "1.05" && centreOf(deleted, "segment_2_norm") == "1" &&
     centreOf(deleted, "segment_3_norm").isEmpty() &&
     settingOf(deleted, "segment_1_norm").startsWith("1.092788 0.5 1.5 ") &&
     settingOf(deleted, "segment_4_norm").isEmpty(),
     section(deleted, "parameterSettings").join(" | "));
  ok("the block's other keys kept", deleted.contains("# the THM run") && deleted.contains("kinematics=triple"));

  std::cout << "\n6. delete every segment the effect line names\n";
  w.statusBar()->clearMessage();
  w.getSegmentsTab()->deleteDataSegment(0);
  const QString noEffect = saveAndCheckRoundTrip(w, path, "reopened and saved: the same file");
  ok("the effect line is removed, and the status bar says so",
     effectLists(noEffect).isEmpty() && w.statusBar()->currentMessage().isEmpty() == false,
     effectLists(noEffect).join(" "));

  std::cout << "\n7. test segments: weightTest[k] follows\n";
  {
    ThmSettings s;
    s.weightTest[1] = "a.dat";
    s.weightTest[2] = "b.dat";
    s.weightTest[7] = "c.dat";  // beyond the table: left alone
    s.weight[1] = "d.dat";
    QVector<int> swap;
    swap << 2 << 1;
    ok("no experiment removed", s.renumberSegments(swap, true).isEmpty());
    ok("swapped", s.weightTest.value(1) == "b.dat" && s.weightTest.value(2) == "a.dat" &&
                  s.weightTest.value(7) == "c.dat" && s.weight.value(1) == "d.dat");
    QVector<int> drop;
    drop << 0 << 1;
    s.renumberSegments(drop, true);
    ok("deleted line dropped", s.weightTest.size() == 2 && s.weightTest.value(1) == "a.dat");
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
