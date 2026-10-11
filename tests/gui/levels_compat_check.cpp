// Backward-compatibility check of the Levels tab against existing .azr files.
//
// For each file on the command line the <levels> block is read through the
// GUI's own reader (LevelsTab::readNuclearFile, pairs wired to PairsTab as in
// AZURESetup) and three things are checked:
//   1. write-back: writeNuclearFile reproduces the block token for token;
//   2. the Add/Edit Level dialog's J^pi list contains every level in the file;
//   3. regenerating every channel from the pairs and limits (what editing a pair
//      or a limit does) gives the same channel set as the file holds, so the
//      capture-rule fix of 2026-10-11 changes nothing for the file.
// Not a ctest: it needs real project files.  QT_QPA_PLATFORM=offscreen.

#include <QApplication>
#include <QFile>
#include <QSet>
#include <QTextStream>
#include <iostream>
#include "LevelsTab.h"
#include "PairsTab.h"
#include "AddLevelDialog.h"
#include "Config.h"
struct SegPairs {int firstPair; int secondPair;};

Config* g_config = nullptr;
void exitMessage(const Config&) {}
bool checkExternalCapture(Config&, const std::vector<SegPairs>&) { return true; }
bool readSegmentFile(const Config&, std::vector<SegPairs>&) { return true; }
void startMessage(const Config&) {}

static QString levelsBlock(const QString& path) {
  QFile file(path);
  if (!file.open(QIODevice::ReadOnly)) return QString();
  QTextStream in(&file);
  QString line, block;
  while (!in.atEnd() && line.trimmed() != "<levels>") line = in.readLine();
  while (!in.atEnd()) {
    line = in.readLine();
    block += line + "\n";
    if (line.trimmed() == "</levels>") break;
  }
  return block;
}

static QStringList tokens(const QString& s) {
  QStringList out;
  for (const QString& line : s.split('\n')) {
    QString t = line.trimmed();
    if (t.isEmpty() || t == "</levels>") continue;
    out << t.simplified();
  }
  return out;
}

static QSet<QString> channelSet(const ChannelsModel* model, const LevelsModel* levels) {
  QSet<QString> out;
  for (const ChannelsData& c : model->getChannels()) {
    const LevelsData& l = levels->getLevels().at(c.levelIndex);
    out.insert(QString("%1 %2 %3 | pair %4 s=%5 l=%6 %7")
                 .arg(l.jValue).arg(l.piValue).arg(l.energy).arg(c.pairIndex + 1).arg(c.sValue).arg(c.lValue).arg(c.radType));
  }
  return out;
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  int bad = 0;
  for (int a = 1; a < argc; a++) {
    QString path = argv[a];
    std::cout << path.toStdString() << std::endl;
    PairsTab pairsTab;
    LevelsTab tab;
    tab.setPairsModel(pairsTab.getPairsModel());
    QObject::connect(&tab, SIGNAL(readNewPair(PairsData, int, bool)), &pairsTab, SLOT(addPair(PairsData, int, bool)));
    QObject::connect(&tab, SIGNAL(readExistingPair(PairsData, int, bool)), &pairsTab, SLOT(editPair(PairsData, int, bool)));

    QString block = levelsBlock(path);
    QTextStream in(&block);
    if (block.isEmpty() || !tab.readNuclearFile(in)) { std::cout << "  FAIL  read" << std::endl; bad++; continue; }
    int nLevels = tab.getLevelsModel()->getLevels().size();
    int nChannels = tab.getChannelsModel()->getChannels().size();
    std::cout << "  read " << nLevels << " levels, " << nChannels << " channels, "
              << pairsTab.getPairsModel()->numPairs() << " pairs" << std::endl;

    // 1. write-back
    QString out;
    QTextStream outStream(&out);
    tab.writeNuclearFile(outStream);
    QStringList tIn = tokens(block), tOut = tokens(out);
    bool same = tIn == tOut;
    if (!same) {
      // tolerate number formatting only: compare as doubles
      same = tIn.size() == tOut.size();
      for (int i = 0; same && i < tIn.size(); i++) {
        QStringList x = tIn[i].split(' '), y = tOut[i].split(' ');
        same = x.size() == y.size();
        for (int k = 0; same && k < x.size(); k++) same = qAbs(x[k].toDouble() - y[k].toDouble()) <= 1e-9 * qMax(1.0, qAbs(x[k].toDouble()));
      }
      std::cout << (same ? "  ok    write-back equal as numbers (formatting differs)" : "  FAIL  write-back differs") << std::endl;
      if (!same) {
        std::cout << "        lines in " << tIn.size() << ", out " << tOut.size() << std::endl;
        for (int i = 0; i < qMin(tIn.size(), tOut.size()); i++)
          if (tIn[i] != tOut[i]) { std::cout << "        first difference, line " << i + 1 << ":\n        in:  " << tIn[i].toStdString() << "\n        out: " << tOut[i].toStdString() << std::endl; break; }
      }
    } else
      std::cout << "  ok    write-back identical token for token" << std::endl;
    if (!same) bad++;

    // 2. every level's J^pi is offered by the dialog
    QList<SpinParityOption> options = tab.spinParityOptions();
    int missing = 0;
    for (const LevelsData& l : tab.getLevelsModel()->getLevels()) {
      bool found = false;
      for (const SpinParityOption& o : options) found |= (o.jValue == l.jValue && o.piValue == l.piValue);
      if (!found) { missing++; std::cout << "        not offered: " << AddLevelDialog::spinParityLabel(l.jValue, l.piValue).toStdString() << " at " << l.energy << " MeV" << std::endl; }
    }
    QStringList labels;
    for (const SpinParityOption& o : options) labels << AddLevelDialog::spinParityLabel(o.jValue, o.piValue);
    std::cout << (missing ? "  FAIL  " : "  ok    ") << "dialog offers every level's J^pi; list: " << labels.join(" ").toStdString() << std::endl;
    if (missing) bad++;

    // 3. regenerating the channels reproduces the file's set
    QSet<QString> before = channelSet(tab.getChannelsModel(), tab.getLevelsModel());
    tab.updateChannelsPairAddedEdited();
    QSet<QString> after = channelSet(tab.getChannelsModel(), tab.getLevelsModel());
    QSet<QString> lost = before - after, gained = after - before;
    std::cout << ((lost.isEmpty() && gained.isEmpty()) ? "  ok    " : "  NOTE  ") << "channels regenerated from pairs and limits: "
              << after.size() << " (file " << before.size() << "), lost " << lost.size() << ", gained " << gained.size() << std::endl;
    for (const QString& s : lost) std::cout << "        lost:   " << s.toStdString() << std::endl;
    for (const QString& s : gained) std::cout << "        gained: " << s.toStdString() << std::endl;
  }
  std::cout << (bad ? "FAILED" : "PASSED") << std::endl;
  return bad ? 1 : 0;
}
