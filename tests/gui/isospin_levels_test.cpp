// Isospin levels in the GUI (headless): a level ID whose lines span several
// compound nuclei (declared in the <isospin> block) is read as one GUI level per
// nucleus, and LevelsTab::writeNuclearFile writes those back consecutively under
// one ID; the block is kept and its level IDs follow the new numbering.
//
// Model: tests/isospin_levels/models/iso_mix.azr (4Li, 4H, 4He; a 2- T=1 level
// with lines in all three and a mix amplitude, a 2- T=0 level in 4He only).
//
// Runs without a display; the CMake target passes QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QFile>
#include <QMap>
#include <QRegularExpression>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <cstdlib>
#include <iostream>
#include "Config.h"
#include "LevelsTab.h"
#include "PairsModel.h"
#include "PairsTab.h"
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
  if (!cond && !detail.isEmpty()) std::cout << "  -- " << detail.toStdString();
  std::cout << std::endl;
  if (!cond) fails++;
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);

  QFile file(ISO_MODEL);
  if (!file.open(QIODevice::ReadOnly)) {
    std::cout << "cannot open " << ISO_MODEL << std::endl;
    return 1;
  }
  QTextStream in(&file);

  PairsTab pairsTab;
  PairsModel& pairs = *pairsTab.getPairsModel();
  LevelsTab levelsTab;
  levelsTab.setPairsModel(&pairs);
  QObject::connect(&levelsTab, SIGNAL(readNewPair(PairsData, int, bool)), &pairsTab, SLOT(addPair(PairsData, int, bool)));
  QObject::connect(&levelsTab, SIGNAL(readExistingPair(PairsData, int, bool)), &pairsTab, SLOT(addPair(PairsData, int, bool)));

  QString line;
  while (line.trimmed() != "<levels>" && !in.atEnd()) line = in.readLine();
  ok("levels block read", levelsTab.readNuclearFile(in));
  while (line.trimmed() != "<isospin>" && !in.atEnd()) line = in.readLine();
  ok("isospin block read", levelsTab.readIsospinBlock(in) && levelsTab.hasIsospinBlock());
  ok("three compound nuclei", pairs.nuclei().size() == 3, QString::number(pairs.nuclei().size()));

  QString levelsOut, blockOut;
  QTextStream lo(&levelsOut), bo(&blockOut);
  ok("levels written", levelsTab.writeNuclearFile(lo));
  ok("isospin block written", levelsTab.writeIsospinBlock(bo));
  lo.flush();
  bo.flush();

  // Group the written lines by level ID: an ID's lines must be consecutive.
  QMap<int, QList<int>> pairsOf;      // level ID -> pair keys of its lines
  QMap<int, double> energyOf;         // level ID -> energy of its first line
  QList<int> idOrder;
  bool consecutive = true;
  for (const QString& raw : levelsOut.split('\n', Qt::SkipEmptyParts)) {
    const QStringList t = raw.simplified().split(' ');
    if (t.size() < 12) continue;
    const int id = t.at(8).toInt();
    if (idOrder.isEmpty() || idOrder.last() != id) {
      if (idOrder.contains(id)) consecutive = false;
      idOrder.append(id);
      energyOf.insert(id, t.at(2).toDouble());
    }
    pairsOf[id].append(t.at(5).toInt());
  }
  ok("every level ID's lines are consecutive", consecutive, idOrder.size() ? QString::number(idOrder.size()) : QString());
  int idT1 = -1, idT0 = -1;
  for (auto it = pairsOf.begin(); it != pairsOf.end(); ++it) {
    QList<int> p = it.value();
    std::sort(p.begin(), p.end());
    if (p == QList<int>({1, 2, 3, 4, 5})) idT1 = it.key();
    if (p == QList<int>({3, 4, 5})) idT0 = it.key();
  }
  ok("the T=1 level is one ID with lines in 4Li, 4H and 4He", idT1 > 0, QString(levelsOut).left(200));
  ok("the T=0 level is one ID with lines in 4He only", idT0 > 0 && idT0 != idT1);
  ok("two level IDs in all", pairsOf.size() == 2, QString::number(pairsOf.size()));

  // The block's IDs follow the new numbering.
  QMap<QString, QStringList> kind;
  for (const QString& raw : blockOut.split('\n', Qt::SkipEmptyParts)) {
    const QStringList t = raw.section('#', 0, 0).simplified().split(' ', Qt::SkipEmptyParts);
    if (t.size() >= 2) kind[t.at(0)].append(t.mid(1).join(' '));
  }
  ok("reference and shift lines kept", kind.value("reference") == QStringList({"3 4 1"}) && kind.value("shift").size() == 2,
     kind.value("reference").join("|") + " / " + kind.value("shift").join("|"));
  ok("level lines name the new IDs",
     kind.value("level").contains(QString("%1 1").arg(idT1)) && kind.value("level").contains(QString("%1 0").arg(idT0)),
     kind.value("level").join("|"));
  ok("the mix line names the T=1 level's new ID", kind.value("mix").size() == 1 && kind.value("mix").first().startsWith(QString("%1 5 1 1").arg(idT1)),
     kind.value("mix").join("|"));

  // ISO_GUI_OUT=<path>: write the project with the GUI's levels and block spliced
  // in, so the engine can compare it with the original (tests/isospin_levels/check.sh).
  if (const char* path = std::getenv("ISO_GUI_OUT")) {
    file.seek(0);
    const QString text = QString::fromUtf8(file.readAll());
    const int a = text.indexOf("<levels>"), b = text.indexOf("</isospin>");
    QFile out(path);
    if (a >= 0 && b >= 0 && out.open(QIODevice::WriteOnly)) {
      QTextStream os(&out);
      os << text.left(a) << "<levels>\n" << levelsOut << "</levels>\n<isospin>\n" << blockOut
         << text.mid(b);
    }
  }

  std::cout << (fails ? "FAILED: " : "all ok: ") << fails << " failure(s)" << std::endl;
  return fails ? 1 : 0;
}
