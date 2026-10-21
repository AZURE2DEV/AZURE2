// "Level N ..." entries of <parameterSettings> must name the engine's parameter energy_N / width_N_M of the
// file they are written into.  The GUI named them by its model row and wrote the levels sorted by J, parity and
// energy: after the first save of a file whose levels were in another order, or with an inactive level ahead of
// active ones, every limit and nuisance prior sat on another level (a prior on the 3.52 MeV level moved to the
// 5.0 MeV one and added 21,904 to chi2; 2026-10-21).  The project below has both: the levels are not in the
// written order and an inactive 1- level comes first.  Each entry carries a distinct error, and after a GUI
// read/write the error must sit on the same physical level and channel, judged by the engine itself (CNuc::Fill
// of the written file).  The older row-numbered block is read too, and a second read/write changes nothing.
#include <QApplication>
#include <QDir>
#include <QFile>
#include <QSettings>
#include <QStringList>
#include <QTemporaryDir>
#include <QTextStream>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include "AZURESetup.h"
#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "Config.h"
#include "JGroup.h"
#include "PPair.h"
struct SegPairs {int firstPair; int secondPair;};
Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

namespace {
int failures = 0;
void ok(const std::string &what, bool cond) {
  std::cout << (cond ? "  ok    " : "  FAIL  ") << what << std::endl;
  if (!cond) failures++;
}

// p + 7Be (pair 1) and p + 7Be*(0.429) (pair 2).  Level lines: J pi E fixE aa ir 2s 2l id active chfix width ...
const char *kPair1 = "0.5 1 1.5 -1 0 1.00783 7.01693 1 4 0.1375 0.1375 0 0 0.0 0 4.3 0 0 0";
const char *kPair2 = "0.5 1 0.5 -1 0.42908 1.00783 7.01693 1 4 0.1375 0.1375 0 0 0.0 0 4.3 0 0 0";
QString levelLine(const char *head, int ir, int twoS, int id, int active, double width) {
  return QString("%1 1 %2 %3 0 %4 %5 0 %6 %7").arg(head).arg(ir).arg(twoS).arg(id).arg(active).arg(width, 0, 'g', 10)
      .arg(ir == 1 ? kPair1 : kPair2);
}
QString project(const QString &settings) {
  QStringList l;
  l << "<config>" << "true #Perform A-Matrix Calculation" << "output/ #Full Path to Output Directory"
    << "checks/ #Full Path to Checks Directory";
  for (int i = 0; i < 8; i++) l << "none #check";
  l << "</config>" << "<levels>";
  l << levelLine("1 -1 2.0 0", 1, 2, 1, 0, 100000.0) << "";                 // inactive, first
  l << levelLine("2 -1 3.52 0", 1, 4, 2, 1, 7120152.86) << "";
  l << levelLine("1 -1 9.0 0", 1, 2, 3, 1, 9780292.75) << levelLine("1 -1 9.0 0", 2, 2, 3, 1, -15363065.0) << "";
  l << levelLine("1 -1 5.0 0", 1, 2, 4, 1, 842304.94) << levelLine("1 -1 5.0 0", 2, 2, 4, 1, 285549.81) << "";
  l << "</levels>" << "<segmentsData>" << "</segmentsData>" << "<segmentsTest>" << "</segmentsTest>"
    << "<targetInt>" << "</targetInt>" << "<parameterSettings>" << settings << "</parameterSettings>"
    << "<lastRun>" << "9997" << "\"\"" << "\"\"" << "0 0" << "0 \"\"" << "-1. -1. -1." << "</lastRun>";
  return l.join("\n") + "\n";
}
QString entry(const QString &name, double value, double error) {
  return QString("%1 %2 0 0 %3 0 1 level 0").arg(name).arg(value).arg(error);
}

// The physical owner of each error: level energy (MeV) and, for a width, the channel's pair.
const std::map<std::pair<double, int>, double> kExpected = {
    {{3.52, 0}, 0.011}, {{9.0, 0}, 0.022}, {{5.0, 0}, 0.033},
    {{3.52, 1}, 111}, {{9.0, 1}, 221}, {{9.0, 2}, 222}, {{5.0, 1}, 331}, {{5.0, 2}, 332}};

bool writeText(const QString &path, const QString &text) {
  QFile f(path);
  if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
  QTextStream(&f) << text;
  return true;
}
QString readText(const QString &path) {
  QFile f(path);
  if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) return QString();
  return QTextStream(&f).readAll();
}
QStringList settingsLines(const QString &text) {
  QStringList out;
  bool in = false;
  for (const QString &line : text.split('\n')) {
    if (line.trimmed() == "<parameterSettings>") in = true;
    else if (line.trimmed() == "</parameterSettings>") in = false;
    else if (in && line.startsWith("Level ")) out << line.trimmed();
  }
  return out;
}

// Every "Level N ..." entry of the written file against the engine's reading of that file.
void checkAgainstEngine(const std::string &label, const QString &path) {
  Config cfg(std::cout);
  cfg.configfile = path.toStdString();
  CNuc compound;
  ok(label + ": the engine reads the written file", compound.Fill(cfg) != -1);
  // energy_N -> (E, -), width_N_M -> (E, pair key)
  std::map<std::string, std::pair<double, int>> owner;
  int n = 1;
  for (int j = 1; j <= compound.NumJGroups(); j++) {
    JGroup *g = compound.GetJGroup(j);
    for (int la = 1; la <= g->NumLevels(); la++, n++) {
      const double e = g->GetLevel(la)->GetE();
      owner["Level " + std::to_string(n) + " Energy (MeV)"] = {e, 0};
      for (int ch = 1; ch <= g->NumChannels(); ch++)
        owner["Level " + std::to_string(n) + " Channel " + std::to_string(ch) + " Width (eV)"] =
            {e, compound.GetPair(g->GetChannel(ch)->GetPairNum())->GetPairKey()};
    }
  }
  const QStringList lines = settingsLines(readText(path));
  ok(label + ": 8 level entries written", lines.size() == 8);
  for (const QString &line : lines) {
    const QStringList t = line.split(' ', Qt::SkipEmptyParts);
    const bool width = line.contains(" Width ");
    const QString name = t.mid(0, width ? 6 : 4).join(' ');
    const double error = t.at(width ? 9 : 7).toDouble();  // name, value, lower, upper, error
    auto it = owner.find(name.toStdString());
    bool good = it != owner.end();
    if (good) {
      double e = std::round(it->second.first * 100.0) / 100.0;
      auto exp = kExpected.find({e, it->second.second});
      good = exp != kExpected.end() && std::fabs(exp->second - error) <= 1e-9 * std::max(1.0, exp->second);
    }
    ok(label + ": " + line.toStdString() + "  -> on its own level", good);
  }
}
}  // namespace

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  QTemporaryDir tmp;
  QSettings::setDefaultFormat(QSettings::IniFormat);
  QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, tmp.path() + "/settings");
  QSettings::setPath(QSettings::NativeFormat, QSettings::UserScope, tmp.path() + "/settings");
  const QString dir = tmp.path();
  QDir().mkpath(dir + "/output");
  QDir().mkpath(dir + "/checks");

  // In the engine numbering of this file as read: active levels group by group (2-: 3.52; 1-: 9.0, 5.0), channels
  // within the group (2-: p+7Be; 1-: p+7Be, p+7Be*).
  QStringList marked;
  marked << "# level numbering: engine"
         << entry("Level 1 Energy (MeV)", 3.52, 0.011) << entry("Level 2 Energy (MeV)", 9.0, 0.022)
         << entry("Level 3 Energy (MeV)", 5.0, 0.033) << entry("Level 1 Channel 1 Width (eV)", 7120152.86, 111)
         << entry("Level 2 Channel 1 Width (eV)", 9780292.75, 221) << entry("Level 2 Channel 2 Width (eV)", -15363065.0, 222)
         << entry("Level 3 Channel 1 Width (eV)", 842304.94, 331) << entry("Level 3 Channel 2 Width (eV)", 285549.81, 332);
  // The older GUI's numbering: model row (the inactive level is row 1) and the channel's row in the channel table.
  QStringList legacy;
  legacy << entry("Level 2 Energy (MeV)", 3.52, 0.011) << entry("Level 3 Energy (MeV)", 9.0, 0.022)
         << entry("Level 4 Energy (MeV)", 5.0, 0.033) << entry("Level 2 Channel 2 Width (eV)", 7120152.86, 111)
         << entry("Level 3 Channel 3 Width (eV)", 9780292.75, 221) << entry("Level 3 Channel 4 Width (eV)", -15363065.0, 222)
         << entry("Level 4 Channel 5 Width (eV)", 842304.94, 331) << entry("Level 4 Channel 6 Width (eV)", 285549.81, 332);

  for (const auto &kind : {std::make_pair(std::string("engine-numbered block"), marked.join("\n")),
                           std::make_pair(std::string("older row-numbered block"), legacy.join("\n"))}) {
    const QString in = dir + "/in.azr", a = dir + "/a.azr", b = dir + "/b.azr";
    ok(kind.first + ": fixture written", writeText(in, project(kind.second)));
    {
      AZURESetup setup;
      ok(kind.first + ": read", setup.readProjectFile(in));
      ok(kind.first + ": written", setup.writeProjectFile(a));
    }
    ok(kind.first + ": marked as engine numbering", readText(a).contains("# level numbering: engine"));
    checkAgainstEngine(kind.first, a);
    {
      AZURESetup setup;
      ok(kind.first + ": written file read back", setup.readProjectFile(a) && setup.writeProjectFile(b));
    }
    ok(kind.first + ": a second read/write leaves the entries as they were",
       settingsLines(readText(a)) == settingsLines(readText(b)));
  }

  // "Load from .sav": the engine's names of the file on disk (here the unsorted one, inactive level first),
  // matched to the GUI's rows by J, parity, energy and channel identity.
  {
    const QString in = dir + "/in.azr";
    writeText(in, project(marked.join("\n")));
    AZURESetup setup;
    setup.readProjectFile(in);
    const QMap<QPair<int, int>, QString> names = setup.EngineParameterNamesOnDisk();
    // rows: 0 inactive 1- 2.0 (channel row 0), 1 2- 3.52 (1), 2 1- 9.0 (2, 3), 3 1- 5.0 (4, 5)
    ok("on disk: no name for the inactive level", !names.contains(qMakePair(0, -1)) && !names.contains(qMakePair(0, 0)));
    ok("on disk: 3.52 MeV = energy_1, its width width_1_1",
       names.value(qMakePair(1, -1)) == "energy_1" && names.value(qMakePair(1, 1)) == "width_1_1");
    ok("on disk: 9.0 MeV = energy_2, widths width_2_1, width_2_2",
       names.value(qMakePair(2, -1)) == "energy_2" && names.value(qMakePair(2, 2)) == "width_2_1" &&
           names.value(qMakePair(2, 3)) == "width_2_2");
    ok("on disk: 5.0 MeV = energy_3, widths width_3_1, width_3_2",
       names.value(qMakePair(3, -1)) == "energy_3" && names.value(qMakePair(3, 4)) == "width_3_1" &&
           names.value(qMakePair(3, 5)) == "width_3_2");
  }

  if (failures) {
    std::cout << failures << " failure(s)" << std::endl;
    return 1;
  }
  std::cout << "all passed" << std::endl;
  return 0;
}
