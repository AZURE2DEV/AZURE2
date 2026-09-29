// Headless round trip of <parameterSettings> through the GUI.
//
// A level row is named "Level N Energy (MeV)" or "Level N Channel k Width
// (...)".  AZURE2 reads N as the N-th level it builds from <levels> (J-groups in
// order of first appearance, active levels only) and k as the k-th channel of
// that level (ParameterLimitsManager::SettingNameForMinuitName on energy_N and
// width_N_k).  The GUI used to name rows by the row of the levels model and the
// row of the channels model instead.  The model keeps the order of the file that
// was opened while the GUI writes its levels sorted by J, parity and energy, so
// a row saved from an unsorted project described one parameter in the file just
// written and another after reopening it: the second save moved limits, errors
// and nuisance priors between levels.  Channel rows beyond the first level
// never meant for AZURE2 what they meant in the GUI.
//
// Checked here: open -> save -> open -> save gives the same file twice; every
// level row names the parameter AZURE2 gives that name to (by the value it
// records); settings written in the engine's numbering of an unsorted file land
// on the intended parameters; and a file written by an earlier GUI, in its
// model-row numbering, keeps its settings on the parameters it meant.
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QDir>
#include <QFile>
#include <QRegExp>
#include <QString>
#include <QStringList>
#include <QTemporaryDir>
#include <QTextStream>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include "AZUREParams.h"
#include "AZURESetup.h"
#include "CNuc.h"
#include "Config.h"
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

// Read as LF text: a Windows checkout gives the test projects CRLF endings,
// and the insertions below anchor on "\n"-terminated lines.
static QString slurp(const QString& path) {
  QFile f(path);
  if(!f.open(QIODevice::ReadOnly)) return QString();
  return QString::fromUtf8(f.readAll()).remove('\r');
}
static void spit(const QString& path, const QString& text) {
  QFile f(path);
  if(f.open(QIODevice::WriteOnly | QIODevice::Truncate)) f.write(text.toUtf8());
}

// open() then save() through the main window; returns the saved text.
static QString openAndSave(AZURESetup& w, const QString& in, const QString& out) {
  QFile::remove(out);
  QFile::copy(in, out);
  w.open(out);
  w.saveProject();
  return slurp(out);
}

// The data rows of the <parameterSettings> section.
static QStringList settingRows(const QString& text) {
  QStringList rows;
  bool in = false;
  for(const QString& raw : text.split('\n')) {
    const QString line = raw.trimmed();
    if(line == "<parameterSettings>") { in = true; continue; }
    if(line.startsWith("<")) { in = false; continue; }
    if(in && !line.isEmpty() && !line.startsWith("#")) rows << line;
  }
  return rows;
}

// Replace (or add, after </targetInt>) the <parameterSettings> section.
static QString withSettings(QString text, const QStringList& rows) {
  const QString block = "<parameterSettings>\n" + rows.join("\n") + "\n</parameterSettings>\n";
  const int a = text.indexOf("<parameterSettings>");
  if(a >= 0) {
    const QString end = "</parameterSettings>\n";
    const int b = text.indexOf(end, a);
    return text.replace(a, b + end.size() - a, block);
  }
  const QString anchor = "</targetInt>\n";
  const int at = text.indexOf(anchor);
  if(at < 0) return QString();   // no anchor: fail the checks, never write mid-line
  return text.insert(at + anchor.size(), block);
}

// AZURE2's own names for the level parameters of a file, with the value it
// reads for each -- straight from CNuc, the way the fit builds them.
static std::map<std::string, double> engineLevelParameters(const QString& path) {
  std::map<std::string, double> result;
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = path.toStdString();
  CNuc compound;
  if(compound.Fill(cfg) == -1) return result;
  AZUREParams params;
  compound.FillMnParams(params.GetMinuitParams(), &cfg);
  const ROOT::Minuit2::MnUserParameters& p = params.GetMinuitParams();
  // energy_N and width_N_k, spelled as ParameterLimitsManager::
  // SettingNameForMinuitName spells them when it looks a setting up.
  QRegExp energy("^energy_(\\d+)$"), width("^width_(\\d+)_(\\d+)$");
  for(unsigned int i = 0; i < p.Params().size(); i++) {
    const QString minuitName = QString::fromStdString(p.GetName(i));
    QString name;
    if(energy.indexIn(minuitName) != -1)
      name = QString("Level %1 Energy (MeV)").arg(energy.cap(1));
    else if(width.indexIn(minuitName) != -1)
      name = QString("Level %1 Channel %2 Width (eV)").arg(width.cap(1)).arg(width.cap(2));
    else
      continue;
    result[name.toStdString()] = p.Value(i);
  }
  return result;
}

// Every level row of the saved file names the parameter AZURE2 gives that name
// to, as the value the row records shows.  Returns the offending rows.
static QStringList rowsNotNamingEngineParameter(const QString& savedPath) {
  const std::map<std::string, double> engine = engineLevelParameters(savedPath);
  QStringList bad;
  if(engine.empty()) return QStringList() << "(the engine could not read the file)";
  for(const QString& row : settingRows(slurp(savedPath))) {
    if(!row.startsWith("Level ")) continue;
    const int nameWords = row.section(' ', 2, 2) == "Channel" ? 6 : 4;
    QString name = row.section(' ', 0, nameWords - 1);
    const double value = row.section(' ', nameWords, nameWords).toDouble();
    // AZURE2 spells every width "(eV)"; the GUI names the unit the value is in.
    name.replace("Width (MeV^(1/2))", "Width (eV)");
    auto it = engine.find(name.toStdString());
    if(it == engine.end() ||
       std::fabs(it->second - value) > 1e-5 * std::max(std::fabs(it->second), 1e-30))
      bad << row;
  }
  return bad;
}

// The row for the parameter whose recorded value is `value` (first match).
static QString rowWithValue(const QStringList& rows, const QString& kind, double value) {
  for(const QString& row : rows) {
    if(!row.contains(kind)) continue;
    const int n = row.section(' ', 2, 2) == "Channel" ? 6 : 4;
    const double v = row.section(' ', n, n).toDouble();
    if(std::fabs(v - value) <= 1e-5 * std::fabs(value)) return row;
  }
  return QString();
}

// The settings columns of a row: everything after the value, less the minuit index.
static QString settingsOf(const QString& row) {
  const int n = row.section(' ', 2, 2) == "Channel" ? 6 : (row.startsWith("Level ") ? 4 : 1);
  return row.section(' ', n + 1, n + 5);
}

int main(int argc, char** argv) {
  // Keep the recent-file list the GUI writes out of the user's settings.
  QTemporaryDir settingsDir;
  qputenv("XDG_CONFIG_HOME", settingsDir.path().toUtf8());
  QApplication app(argc, argv);
  QCoreApplication::setOrganizationName("AZURE2-tests");
  QCoreApplication::setApplicationName("parameter_settings_test");

  QTemporaryDir work;
  AZURESetup w;
  const QString src7Li = QString(AZURE2_SOURCE_DIR) + "/tests/7Li_p_a/7Li_p_a.azr";
  const QString src13N = QString(AZURE2_SOURCE_DIR) + "/tests/13N/13N.azr";
  const QString plain7Li = slurp(src7Li);
  const QString plain13N = slurp(src13N);
  ok("found tests/7Li_p_a and tests/13N", !plain7Li.isEmpty() && !plain13N.isEmpty());

  // 1. 7Li_p_a: no <parameterSettings>, levels not in the order the GUI writes
  //    them (the J=0+ levels come first, the 0- one near the end).
  std::cout << "\n1. tests/7Li_p_a, no settings of its own\n";
  {
    spit(work.filePath("li.in"), plain7Li);
    const QString first = openAndSave(w, work.filePath("li.in"), work.filePath("li.azr"));
    spit(work.filePath("li2.in"), first);
    const QString second = openAndSave(w, work.filePath("li2.in"), work.filePath("li.azr"));
    ok("GUI wrote the rows", settingRows(first).size() > 100,
       QString::number(settingRows(first).size()));
    ok("second save is byte-identical to the first", second == first);
    spit(work.filePath("li3.in"), first);
    const QStringList bad = rowsNotNamingEngineParameter(work.filePath("li3.in"));
    ok("every level row names AZURE2's parameter of that name", bad.isEmpty(), bad.mid(0, 3).join(" | "));
  }

  // 2. The same project with settings on a level energy, a width on channel
  //    5 of the third level AZURE2 builds, a width on channel 2 of a later
  //    level and a normalization, named in AZURE2's numbering of this
  //    (unsorted) file: 20.2 MeV 0+ is its level 1, 27.494 MeV 0+ its level 2,
  //    22 MeV 1- its level 6.
  std::cout << "\n2. tests/7Li_p_a with settings in AZURE2's numbering of the unsorted file\n";
  {
    const QStringList rows = QStringList()
        << "Level 1 Energy (MeV) 20.2 19.9 20.5 0.1 0 1 level 0"
        << "Level 2 Channel 5 Width (eV) 3.9622e+07 3e+07 5e+07 2e+06 0 1 level 1"
        << "Level 6 Channel 2 Width (eV) 546170 100000 900000 50000 0 1 level 2"
        << "segment_1_norm 0.000597139 0.0004 0.0008 0 0 1 norm 3";
    // What AZURE2 says those names are, in the file as given.
    spit(work.filePath("probe.in"), withSettings(plain7Li, rows));
    const std::map<std::string, double> engine = engineLevelParameters(work.filePath("probe.in"));
    ok("AZURE2 reads them as intended in the input",
       engine.count("Level 1 Energy (MeV)") && engine.at("Level 1 Energy (MeV)") == 20.2 &&
       engine.count("Level 2 Channel 5 Width (eV)") && engine.at("Level 2 Channel 5 Width (eV)") == 39621970 &&
       engine.count("Level 6 Channel 2 Width (eV)") && engine.at("Level 6 Channel 2 Width (eV)") == 546170);

    const QString first = openAndSave(w, work.filePath("probe.in"), work.filePath("probe.azr"));
    spit(work.filePath("probe2.in"), first);
    const QString second = openAndSave(w, work.filePath("probe2.in"), work.filePath("probe.azr"));
    ok("second save is byte-identical to the first", second == first);

    const QStringList saved = settingRows(first);
    const QString energy = rowWithValue(saved, "Energy", 20.2);
    const QString width5 = rowWithValue(saved, "Width", 39621970);
    const QString width2 = rowWithValue(saved, "Width", 546170);
    ok("energy limits, error and prior kept on the 20.2 MeV level",
       settingsOf(energy) == "19.9 20.5 0.1 0 1", energy);
    ok("width settings kept on the 27.494 MeV level, channel 5",
       settingsOf(width5) == "3e+07 5e+07 2e+06 0 1", width5);
    ok("width settings kept on the 22 MeV level, channel 2",
       settingsOf(width2) == "100000 900000 50000 0 1", width2);
    ok("renamed to AZURE2's numbering of the saved (sorted) file",
       energy.startsWith("Level 2 Energy (MeV) ") && width5.startsWith("Level 3 Channel 5 Width (eV) ") &&
       width2.startsWith("Level 6 Channel 2 Width (eV) "),
       energy + " | " + width5 + " | " + width2);
    ok("normalization settings kept",
       saved.contains("segment_1_norm 0.000597139 0.0004 0.0008 0 0 1 norm 203"),
       saved.filter("segment_1_norm").join(" | "));
    int limited = 0;
    for(const QString& row : saved)
      if(row.startsWith("Level ") && settingsOf(row).section(' ', 0, 1) != "0 0") limited++;
    ok("no other level row picked up limits", limited == 3, QString::number(limited));
    spit(work.filePath("probe3.in"), first);
    const QStringList bad = rowsNotNamingEngineParameter(work.filePath("probe3.in"));
    ok("every level row names AZURE2's parameter of that name", bad.isEmpty(), bad.mid(0, 3).join(" | "));
  }

  // 3. tests/13N was saved by an earlier GUI, which numbered channels by their
  //    row in the channels model: its "Level 3 Channel 5" is the first channel
  //    of the third level, AZURE2's "Level 3 Channel 1".  Give that row and the
  //    level-4 one distinctive settings; they must stay on those parameters.
  std::cout << "\n3. tests/13N, settings written by an earlier GUI\n";
  {
    QString legacy = plain13N;
    legacy.replace("Level 3 Channel 5 Width (eV) 34000 0 0 3400 0 0 level 2",
                   "Level 3 Channel 5 Width (eV) 34000 1000 90000 777 0 1 level 2");
    legacy.replace("Level 4 Channel 8 Width (eV) 5400 0 0 540 0 0 level 4",
                   "Level 4 Channel 8 Width (eV) 5400 10 9000 55 0 1 level 4");
    ok("edited the legacy rows", legacy != plain13N);
    spit(work.filePath("n.in"), legacy);
    const QString first = openAndSave(w, work.filePath("n.in"), work.filePath("n.azr"));
    spit(work.filePath("n2.in"), first);
    const QString second = openAndSave(w, work.filePath("n2.in"), work.filePath("n.azr"));
    ok("second save is byte-identical to the first", second == first);
    const QStringList saved = settingRows(first);
    const QString w34000 = rowWithValue(saved, "Width", 34000);
    const QString w5400 = rowWithValue(saved, "Width", 5400);
    ok("settings stay on the 34000 eV width", settingsOf(w34000) == "1000 90000 777 0 1", w34000);
    ok("settings stay on the 5400 eV width", settingsOf(w5400) == "10 9000 55 0 1", w5400);
    ok("both renamed to AZURE2's numbering",
       w34000.startsWith("Level 3 Channel 1 Width (eV) ") && w5400.startsWith("Level 4 Channel 2 Width (eV) "),
       w34000 + " | " + w5400);
    ok("segment rows as before", saved.filter("segment_").size() == settingRows(plain13N).filter("segment_").size());
    spit(work.filePath("n3.in"), first);
    const QStringList bad = rowsNotNamingEngineParameter(work.filePath("n3.in"));
    ok("every level row names AZURE2's parameter of that name", bad.isEmpty(), bad.mid(0, 3).join(" | "));
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
