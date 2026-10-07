// Headless check that the GUI tells particle pairs apart the way AZURE2 does.
//
// AZURE2 knows a pair by its key, field 6 of a <levels> line (CNuc::IsPair
// compares keys only): the first line with a key defines that pair, two pairs
// may be alike in everything else.  examples/c12c12_tumino2018 relies on that.
// Its pairs 1 and 6 are both 12C+12C and differ only in the THM binding energy
// of field 32 (pair 1, the THM entrance, has B = 10.272312 MeV; pair 6, the
// entrance of the direct data, has none).  The GUI used to look a pair up by
// its physics without B, so pair 6 was taken for a duplicate of pair 1: a
// modal "Duplicate Pair" box during the read, no sixth pair, and channels
// pointing past the end of the pair list.  Its <segmentsTest> lines come in
// twins alike but for the active flag, which the segment reader likewise
// refused as duplicates (another modal box); a line read from the file is
// now kept, as the engine keeps it.
//
// Checked here, on examples/c12c12_tumino2018:
//  1. open: six pairs, 1 and 6 equal but for B; every channel's pair key is
//     the key of its line in the file, in file order;
//  2. save: every <levels> number of the file comes back exactly (lines
//     compared as a set, the level-id field left out), the sections the GUI
//     carries verbatim are verbatim, and the GUI's own output saved again is
//     byte-identical;
//  3. every <parameterSettings> level row the GUI (Fitting tab) wrote names
//     AZURE2's parameter of that name, as the value it records;
//  4. the engine's total chi-squared of the GUI-saved file is the original's
//     (112.0758: data 111.651 + the direct norms' priors centred on 1, kept by
//     the GUI as prior_centre rows of <parameterSettings>);
//  5. a trivial edit (one width of a pair-6 channel) round-trips: written,
//     read back on the pair-6 channel, nothing else moved, saved again
//     byte-identical;
//  6. the THM workspace opened and accepted untouched leaves the project byte
//     for byte, and its Channels page shows pair 1 (the THM entrance) only;
//  7. an interactive add refuses a true duplicate, not a pair that differs
//     in B.
//
// Files are read with CR removed (a Windows checkout gives CRLF endings).
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QDir>
#include <QFile>
#include <QLineEdit>
#include <QProcess>
#include <QRegExp>
#include <QString>
#include <QStringList>
#include <QTableWidget>
#include <QTemporaryDir>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include "AZUREParams.h"
#include "AZURESetup.h"
#include "CNuc.h"
#include "ChannelsModel.h"
#include "FittingTab.h"
#include "Config.h"
#include "LevelsTab.h"
#include "PairsModel.h"
#include "PairsTab.h"
#include "SegmentsTab.h"
#include "SegmentsTestModel.h"
#include "ThmChannelsPage.h"
#include "ThmSettings.h"
#include "ThmWorkspace.h"
#include "ThmExperimentsPage.h"
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
#ifndef AZURE2_BINARY
#error "AZURE2_BINARY must be defined"
#endif

static int fails = 0;
static void ok(const char* what, bool cond, const QString& detail = QString()) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what;
  if(!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if(!cond) fails++;
}

static QString slurp(const QString& path) {
  QFile f(path);
  if(!f.open(QIODevice::ReadOnly)) return QString();
  return QString::fromUtf8(f.readAll()).remove('\r');
}
static void spit(const QString& path, const QString& text) {
  QFile f(path);
  if(f.open(QIODevice::WriteOnly | QIODevice::Truncate)) f.write(text.toUtf8());
}

// The text between <tag> and </tag>, or a null string.
static QString section(const QString& text, const QString& tag) {
  const QString a = "<" + tag + ">\n", b = "</" + tag + ">";
  const int i = text.indexOf(a), j = text.indexOf(b);
  return i < 0 || j < i ? QString() : text.mid(i + a.size(), j - i - a.size());
}

// The non-blank lines of <levels>, split into fields.
static QList<QStringList> levelLines(const QString& text) {
  QList<QStringList> lines;
  for(const QString& line : section(text, "levels").split('\n')) {
    const QStringList f = line.simplified().split(' ', Qt::SkipEmptyParts);
    if(!f.isEmpty()) lines << f;
  }
  return lines;
}

// Every <levels> line as doubles, level-id field (9th) left out, fields 32-33
// (optional, 0 when absent) filled in; sorted, so the GUI's order is moot.
typedef std::vector<double> Row;
static std::vector<Row> levelRows(const QString& text, bool& parsed) {
  std::vector<Row> rows;
  parsed = true;
  for(const QStringList& f : levelLines(text)) {
    Row row;
    for(int i = 0; i < f.size(); i++) {
      bool good = false;
      const double v = f[i].toDouble(&good);
      if(!good) parsed = false;
      if(i != 8) row.push_back(v);
    }
    while(row.size() < 32) row.push_back(0.);
    rows.push_back(row);
  }
  std::sort(rows.begin(), rows.end());
  return rows;
}

// The data rows of <parameterSettings>.
static QStringList settingRows(const QString& text) {
  QStringList rows;
  for(const QString& raw : section(text, "parameterSettings").split('\n')) {
    const QString line = raw.trimmed();
    if(!line.isEmpty() && !line.startsWith("#")) rows << line;
  }
  return rows;
}

// AZURE2's names for the level parameters of a file, with their values,
// straight from CNuc as the fit builds them.
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

// Level rows of the saved file that do not name AZURE2's parameter of that
// name with the value the row records; count = level rows seen.
static QStringList rowsNotNamingEngineParameter(const QString& savedPath, int& count) {
  const std::map<std::string, double> engine = engineLevelParameters(savedPath);
  QStringList bad;
  count = 0;
  if(engine.empty()) return QStringList() << "(the engine could not read the file)";
  for(const QString& row : settingRows(slurp(savedPath))) {
    if(!row.startsWith("Level ")) continue;
    count++;
    const int nameWords = row.section(' ', 2, 2) == "Channel" ? 6 : 4;
    QString name = row.section(' ', 0, nameWords - 1);
    const double value = row.section(' ', nameWords, nameWords).toDouble();
    name.replace("Width (MeV^(1/2))", "Width (eV)");
    auto it = engine.find(name.toStdString());
    if(it == engine.end() || std::fabs(it->second - value) > 1e-5 * std::max(std::fabs(it->second), 1e-30))
      bad << row;
  }
  return bad;
}

// "Calculate Segments From Data" on dir/name; the stdout total, or "".
static QString engineChi2(const QString& dir, const QString& name) {
  QDir(dir).mkpath("output");
  QDir(dir).mkpath("checks");
  QProcess p;
  p.setWorkingDirectory(dir);
  p.setProcessChannelMode(QProcess::MergedChannels);
  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  if(!env.contains("OMP_NUM_THREADS")) env.insert("OMP_NUM_THREADS", "1");
  p.setProcessEnvironment(env);
  p.start(AZURE2_BINARY, QStringList() << "--no-gui" << "--no-readline" << name);
  if(!p.waitForStarted(30000)) return QString();
  p.write("1\n\n\n7\n");
  p.closeWriteChannel();
  if(!p.waitForFinished(600000)) { p.kill(); return QString(); }
  const QString out = QString::fromUtf8(p.readAll());
  QRegExp rx("Total Chi-Squared: ([-+0-9.eE]+)");
  return rx.indexIn(out) >= 0 ? rx.cap(1) : QString();
}

static void copyDir(const QString& from, const QString& to) {
  QDir().mkpath(to);
  for(const QString& f : QDir(from).entryList(QDir::Files)) QFile::copy(from + "/" + f, to + "/" + f);
}

// The pair key of every channel of the GUI's model, in model (= read) order.
static QStringList channelKeys(AZURESetup& w) {
  QStringList keys;
  for(const ChannelsData& c : w.getLevelsTab()->getChannelsModel()->getChannels()) keys << QString::number(c.pairIndex + 1);
  return keys;
}

static bool samePairBut(const PairsData& a, const PairsData& b) {
  return a.lightJ == b.lightJ && a.lightPi == b.lightPi && a.lightZ == b.lightZ && a.lightM == b.lightM &&
         a.lightG == b.lightG && a.heavyJ == b.heavyJ && a.heavyPi == b.heavyPi && a.heavyZ == b.heavyZ &&
         a.heavyM == b.heavyM && a.heavyG == b.heavyG && a.excitationEnergy == b.excitationEnergy &&
         a.seperationEnergy == b.seperationEnergy && a.channelRadius == b.channelRadius &&
         a.pairType == b.pairType && a.ecMultMask == b.ecMultMask;
}

int main(int argc, char** argv) {
  QTemporaryDir settingsDir;
  qputenv("XDG_CONFIG_HOME", settingsDir.path().toUtf8());
  QApplication app(argc, argv);
  ThmExperimentsPage::setDerivedDelay(0);  // derived values at once, as the checks read them
  QCoreApplication::setOrganizationName("AZURE2-tests");
  QCoreApplication::setApplicationName("pair_identity_test");

  const QString srcDir = QString(AZURE2_SOURCE_DIR) + "/examples/c12c12_tumino2018";
  const QString original = slurp(srcDir + "/c12c12_tumino2018.azr");
  ok("found examples/c12c12_tumino2018", original.contains("</levels>"));

  QTemporaryDir work;
  const QString origDir = work.filePath("orig"), guiDir = work.filePath("gui");
  copyDir(srcDir + "/data", origDir + "/data");
  copyDir(srcDir + "/data", guiDir + "/data");
  spit(origDir + "/run.azr", original);
  spit(guiDir + "/run.azr", original);
  const QString path = guiDir + "/run.azr";

  // The file's own pair keys, line by line, and its pairs 1 and 6.
  QStringList fileKeys;
  for(const QStringList& f : levelLines(original)) fileKeys << (f.size() > 5 ? f[5] : QString("?"));
  ok("the example has pair-6 lines", fileKeys.contains("6") && fileKeys.contains("1"));

  AZURESetup w;
  PairsModel* pairs = w.getPairsTab()->getPairsModel();

  // 1. Open.
  std::cout << "\n1. open\n";
  w.open(path);
  {
    const QList<PairsData> p = pairs->getPairs();
    ok("six pairs", p.size() == 6, QString::number(p.size()));
    if(p.size() == 6) {
      ok("pairs 1 and 6: both 12C+12C, alike but for B",
         p[0].lightZ == 6 && p[0].heavyZ == 6 && p[0].lightM == 11.99671 && samePairBut(p[0], p[5]));
      ok("pair 1 has the THM binding energy, pair 6 none",
         p[0].bindingEnergy == 10.272312 && p[5].bindingEnergy == 0.0,
         QString("%1 %2").arg(p[0].bindingEnergy).arg(p[5].bindingEnergy));
      ok("pairs 2-5: p+23Na, p+23Na*, a+20Ne, a+20Ne*",
         p[1].lightZ == 1 && p[2].lightZ == 1 && p[2].excitationEnergy == 0.4402 && p[3].lightZ == 2 &&
             p[4].lightZ == 2 && p[4].excitationEnergy == 1.6337);
    }
    // <segmentsTest> has each 12C+12C extrapolation twice, off and on, lines
    // alike but for the active flag; the engine reads both.
    int testLines = 0;
    for(const QString& line : section(original, "segmentsTest").split('\n')) testLines += !line.trimmed().isEmpty();
    ok("every <segmentsTest> line read, twins alike but for the active flag included",
       w.getSegmentsTab()->getSegmentsTestModel()->getLines().size() == testLines && testLines == 12,
       QString("%1 of %2").arg(w.getSegmentsTab()->getSegmentsTestModel()->getLines().size()).arg(testLines));
    const QStringList keys = channelKeys(w);
    ok("every channel's pair key is its line's", keys == fileKeys,
       keys.mid(0, 12).join(",") + " | " + fileKeys.mid(0, 12).join(","));
  }

  // 2. Save untouched.
  std::cout << "\n2. save\n";
  w.saveProject();
  const QString first = slurp(path);
  {
    bool p1 = false, p2 = false;
    const std::vector<Row> a = levelRows(original, p1), b = levelRows(first, p2);
    ok("both <levels> sections parse", p1 && p2);
    ok("every <levels> number comes back exactly", a == b && !a.empty(),
       QString("%1 vs %2 lines").arg(a.size()).arg(b.size()));
    // The GUI lays out the segment and <targetInt> columns its own way.
    auto tokens = [](const QString& text) {
      QStringList lines;
      for(const QString& line : text.split('\n'))
        if(!line.trimmed().isEmpty()) lines << line.simplified();
      return lines;
    };
    QString which;
    for(const QString& tag : QStringList({"segmentsData", "segmentsTest", "targetInt"}))
      if(tokens(section(first, tag)) != tokens(section(original, tag))) which += tag + " ";
    ok("segments and <targetInt>: the same lines, field for field", which.isEmpty(), which);
    ok("<thm> verbatim", section(first, "thm") == section(original, "thm"));
    ok("keys of the saved lines are the file's (as a multiset)", [&] {
      QStringList x, y = fileKeys;
      for(const QStringList& f : levelLines(first)) x << f[5];
      x.sort();
      y.sort();
      return x == y;
    }());
    spit(guiDir + "/again.in", first);
    QFile::remove(path);
    QFile::copy(guiDir + "/again.in", path);
    w.open(path);
    ok("reopened: six pairs, same channel keys (in the saved order)", pairs->numPairs() == 6 && [&] {
      QStringList saved;
      for(const QStringList& f : levelLines(first)) saved << f[5];
      return channelKeys(w) == saved;
    }());
    w.saveProject();
    ok("the GUI's own output saved again is byte-identical", slurp(path) == first);
  }
  if(qEnvironmentVariableIsSet("PAIR_IDENTITY_KEEP")) spit(qEnvironmentVariable("PAIR_IDENTITY_KEEP"), first);

  // 3. The Fitting tab's rows name the engine's parameters.
  std::cout << "\n3. parameter names\n";
  {
    int count = 0;
    const QStringList bad = rowsNotNamingEngineParameter(path, count);
    ok("every level row names AZURE2's parameter of that name", bad.isEmpty() && count > 200,
       QString("%1 rows; ").arg(count) + bad.mid(0, 3).join(" | "));
  }

  // 4. The engine sees the same project.
  std::cout << "\n4. chi-squared\n";
  {
    const QString c0 = engineChi2(origDir, "run.azr");
    const QString c1 = engineChi2(guiDir, "run.azr");
    std::cout << "  chi2 original " << c0.toStdString() << "   GUI-saved " << c1.toStdString() << std::endl;
    ok("engine chi-squared of the GUI-saved file equals the original's", !c0.isEmpty() && c0 == c1);
    ok("and is the example's 112.0758 (data 111.651 + priors 0.4248)",
       std::fabs(c0.toDouble() - 112.0758) < 0.01, c0);
  }

  // 5. A trivial edit: the width of the first pair-6 channel.
  std::cout << "\n5. edit\n";
  {
    w.open(path);
    ChannelsModel* cm = w.getLevelsTab()->getChannelsModel();
    const QList<ChannelsData> before = cm->getChannels();
    int target = -1;
    for(int i = 0; i < before.size() && target < 0; i++)
      if(before[i].pairIndex == 5) target = i;
    ok("found a pair-6 channel", target >= 0);
    if(target >= 0) {
      const double edited = 1.25e-17;
      cm->setData(cm->index(target, 6), edited, Qt::EditRole);
      // The Fitting tab, whose <parameterSettings> rows record the values,
      // takes them from the levels when it is shown; as if it had been.
      w.getFittingTab()->populateFromCurrentGUIState();
      w.saveProject();
      const QString saved = slurp(path);
      int lines = 0;
      for(const QStringList& f : levelLines(saved))
        if(f[5] == "6" && f[11].toDouble() == edited) lines++;
      ok("written on one pair-6 line", lines == 1, QString::number(lines));
      bool p1 = false, p2 = false;
      std::vector<Row> a = levelRows(first, p1), b = levelRows(saved, p2);
      int differing = 0;
      for(size_t i = 0; i < std::min(a.size(), b.size()); i++) differing += a[i] != b[i];
      ok("no other <levels> line moved", a.size() == b.size() && differing <= 2, QString::number(differing));
      w.open(path);
      const QList<ChannelsData> after = cm->getChannels();
      bool found = false;
      for(const ChannelsData& c : after) found = found || (c.pairIndex == 5 && c.reducedWidth == edited);
      ok("read back on a pair-6 channel, six pairs", found && pairs->numPairs() == 6);
      w.saveProject();
      ok("edited project saved again is byte-identical", slurp(path) == saved);
      if(qEnvironmentVariableIsSet("PAIR_IDENTITY_KEEP")) {
        spit(qEnvironmentVariable("PAIR_IDENTITY_KEEP") + ".edit1", saved);
        spit(qEnvironmentVariable("PAIR_IDENTITY_KEEP") + ".edit2", slurp(path));
      }
    }
  }

  // 6. The THM workspace, untouched.
  std::cout << "\n6. THM workspace\n";
  {
    spit(path, first);
    w.open(path);
    w.saveProject();
    const QString before = slurp(path);
    ThmSettings s;
    QString err;
    ok("the <thm> block opens", w.thmSettings(s, &err), err);
    {
      ThmWorkspace ws(&w, s);
      QTableWidget* t = ws.channelsPage->pairTable;
      QStringList shown;
      for(int r = 0; r < t->rowCount(); r++) shown << t->item(r, 0)->text();
      ok("Channels page: pair 1 only (the THM entrance, B set)", shown == QStringList({"1"}), shown.join(","));
      ok("Channels page: pair 1's binding energy", qobject_cast<QLineEdit*>(t->cellWidget(0, 3)) &&
                                                       qobject_cast<QLineEdit*>(t->cellWidget(0, 3))->text().toDouble() == 10.272312);
      ok("untouched: nothing refused", ws.validate().isEmpty(), ws.validate());
      ws.accept();
    }
    w.saveProject();
    ok("untouched workspace: byte-identical save", slurp(path) == before);
  }

  // 7. Interactive duplicates.
  std::cout << "\n7. duplicates\n";
  {
    PairsData p6 = pairs->getPairs().value(5);
    ok("pair 6 is found as itself", pairs->isPair(p6) == 5);
    p6.bindingEnergy = 1.0;
    ok("a pair differing in B is not a duplicate", pairs->isPair(p6) == -1);
  }

  std::cout << (fails ? "FAILED" : "PASSED") << std::endl;
  return fails ? 1 : 0;
}
