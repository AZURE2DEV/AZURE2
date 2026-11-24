// The Levels tab's <levels> rows (docs/upgrade/qt_gui.md section 2): read from the engine's reader
// (azure2::front::Project::Levels, the columns of io::LevelColumnNames) into the pairs, levels and channels, and made
// back into rows for the engine's writer (Project::SetLevels: an unchanged row keeps its line, an edited one its
// layout).  A value the GUI did not change keeps the token it was read with, so that a number is never rounded by
// the GUI; a new value is written in the shortest form that reads back to the same double.  (Before 2026-11-18 the
// GUI rewrote the whole block sorted by nucleus, J, parity and energy, writtenLevelOrder.)
#include <QLocale>
#include <QMap>
#include <QSet>

#include "LevelsTab.h"

namespace {

// The columns of a row by their names (azure2::front::LevelColumnNames, the engine's io::LevelColumnNames): the GUI
// types no column position.
struct Columns {
  static int At(const char *name) { return azure2::front::LevelColumnIndex(name); }
  const int kLevelJ = At("levelJ"), kLevelPi = At("levelPi"), kLevelE = At("levelE"), kLevelFix = At("levelFix"),
            kAa = At("aa"), kIr = At("ir"), kS2 = At("s2"), kL2 = At("l2"), kLevelId = At("levelID"),
            kIsActive = At("isActive"), kChannelFix = At("channelFix"), kGamma = At("gamma"), kJ1 = At("j1"),
            kPi1 = At("pi1"), kJ2 = At("j2"), kPi2 = At("pi2"), kE2 = At("e2"), kM1 = At("m1"), kM2 = At("m2"),
            kZ1 = At("z1"), kZ2 = At("z2"), kEntranceSepE = At("entranceSepE"), kSepE = At("sepE"), kJ3 = At("j3"),
            kPi3 = At("pi3"), kE3 = At("e3"), kPType = At("pType"), kChRad = At("chRad"), kG1 = At("g1"),
            kG2 = At("g2"), kEcMultMask = At("ecMultMask"), kBindingE = At("bindingE"), kGammaIsRWA = At("gammaIsRWA"),
            kE1 = At("e1");
  const int kColumns = (int)azure2::front::LevelColumnNames().size();
};

QString Shortest(double v) { return QString::number(v, 'g', QLocale::FloatingPointShortest); }

double Number(const std::vector<std::string> &t, int c) {
  return c < (int)t.size() && !t[c].empty() ? QString::fromStdString(t[c]).toDouble() : 0.;
}

}  // namespace

bool LevelsTab::readLevels(const std::vector<azure2::front::LevelRow> &rows, QString &why) {
  const Columns col;
  const int kE1 = col.kE1, kIr = col.kIr, kPType = col.kPType, kL2 = col.kL2, kLevelPi = col.kLevelPi, kJ1 = col.kJ1,
            kPi1 = col.kPi1, kZ1 = col.kZ1, kM1 = col.kM1, kG1 = col.kG1, kJ2 = col.kJ2, kPi2 = col.kPi2,
            kZ2 = col.kZ2, kM2 = col.kM2, kG2 = col.kG2, kE2 = col.kE2, kSepE = col.kSepE, kChRad = col.kChRad,
            kEcMultMask = col.kEcMultMask, kIsActive = col.kIsActive, kLevelFix = col.kLevelFix, kLevelJ = col.kLevelJ,
            kLevelE = col.kLevelE, kLevelId = col.kLevelId, kChannelFix = col.kChannelFix, kS2 = col.kS2,
            kGamma = col.kGamma;
  maxLSpin->blockSignals(true);
  maxMultSpin->blockSignals(true);
  maxNumMultSpin->blockSignals(true);
  maxLSpin->setValue(0);
  maxMultSpin->setValue(0);
  maxNumMultSpin->setValue(0);
  int maxLValue = 0, maxMultValue = 0, maxNumMultValue = 0;
  int lastPair = -1, thisNumMult = 0;
  int rowNumber = 0;
  for (const azure2::front::LevelRow &row : rows) {
    rowNumber++;
    if (row.tokens.empty()) continue;
    if ((int)row.tokens.size() < kE1) {
      why = tr("<levels> row %1 has %2 columns").arg(rowNumber).arg(row.tokens.size());
      return false;
    }
    const std::vector<std::string> &t = row.tokens;
    const int ir = (int)Number(t, kIr), pairType = (int)Number(t, kPType), channelL = (int)Number(t, kL2);
    const double levelPi = Number(t, kLevelPi);
    thisNumMult = (ir == lastPair) ? thisNumMult + 1 : 1;
    lastPair = ir;
    if (pairType == 0) maxLValue = std::max(maxLValue, channelL / 2);
    else if (pairType == 10) {
      maxMultValue = std::max(maxMultValue, channelL / 2);
      maxNumMultValue = std::max(maxNumMultValue, thisNumMult);
    }
    PairsData pair = {Number(t, kJ1), (int)Number(t, kPi1), (int)Number(t, kZ1), Number(t, kM1), Number(t, kG1),
                      Number(t, kJ2), (int)Number(t, kPi2), (int)Number(t, kZ2), Number(t, kM2), Number(t, kG2),
                      Number(t, kE2), Number(t, kSepE), Number(t, kChRad), pairType, (int)Number(t, kEcMultMask),
                      Number(t, kE1)};
    const int pairIndex = ir - 1;
    if (pairIndex < 0) {
      why = tr("<levels> row %1: pair key %2").arg(rowNumber).arg(ir);
      return false;
    }
    // The pair goes to its own key, whatever the other keys hold: a project may name the same pair under two keys
    // (c12c12_tumino2018: key 1 for the THM data, key 6 for the direct data, each with its own widths).  Comparing
    // with every pair (isPair) refused such a key with a box per row and left its channels without a pair.
    if (pairsModel->numPairs() < ir) emit readNewPair(pair, pairIndex, true);
    else if (!pairsModel->pairAt(pairIndex, pair)) emit readExistingPair(pair, pairIndex, true);
    LevelsData level = {(int)Number(t, kIsActive), (int)Number(t, kLevelFix), Number(t, kLevelJ), (int)levelPi,
                        Number(t, kLevelE)};
    const Nucleus nucleus = PairsModel::nucleusOf(pair);  // the engine's rule: a level belongs to its channels' nucleus
    level.nucZ = nucleus.Z;
    level.nucA = nucleus.A;
    level.fileId = (int)Number(t, kLevelId);
    int levelIndex = levelsModel->isLevel(level);
    if (levelIndex == -1) {
      addLevel(level, true);
      levelIndex = levelsModel->isLevel(level);
    }
    QChar radType = 'P';
    if (pairType == 20)
      radType = channelL == 0 ? 'F' : 'G';
    else if (pairType == 10) {
      const int parityChange = (int)Number(t, kPi2) * (int)levelPi;
      const bool odd = (channelL / 2) % 2 != 0;
      radType = odd ? (parityChange == -1 ? 'E' : 'M') : (parityChange == -1 ? 'M' : 'E');
    }
    const int channelIndex = channelsModel->getChannels().size();
    channelsModel->insertRows(channelIndex, 1, QModelIndex());
    channelsModel->setData(channelsModel->index(channelIndex, 0), (int)Number(t, kChannelFix), Qt::EditRole);
    channelsModel->setData(channelsModel->index(channelIndex, 1), levelIndex, Qt::EditRole);
    channelsModel->setData(channelsModel->index(channelIndex, 2), pairIndex, Qt::EditRole);
    channelsModel->setData(channelsModel->index(channelIndex, 3), Number(t, kS2) / 2.0, Qt::EditRole);
    channelsModel->setData(channelsModel->index(channelIndex, 4), channelL / 2, Qt::EditRole);
    channelsModel->setData(channelsModel->index(channelIndex, 5), radType, Qt::EditRole);
    channelsModel->setData(channelsModel->index(channelIndex, 6), Number(t, kGamma), Qt::EditRole);
    QStringList tokens;
    for (const std::string &s : t) tokens << QString::fromStdString(s);
    channelsModel->setFileTokens(channelIndex, tokens);
  }
  maxLSpin->setValue(maxLValue);
  maxMultSpin->setValue(maxMultValue);
  maxNumMultSpin->setValue(maxNumMultValue);
  updateMissingSpinParities();
  maxLSpin->blockSignals(false);
  maxMultSpin->blockSignals(false);
  maxNumMultSpin->blockSignals(false);
  channelsView->resizeRowsToContents();
  levelsView->resizeRowsToContents();
  return true;
}

azure2::front::LevelRow LevelsTab::channelRow(const ChannelsData &ch, const LevelsData &level, const PairsData &pair,
                                               const PairsData &entrance, int id) const {
  const Columns col;
  const int kColumns = col.kColumns, kLevelJ = col.kLevelJ, kLevelPi = col.kLevelPi, kLevelE = col.kLevelE,
            kLevelFix = col.kLevelFix, kAa = col.kAa, kIr = col.kIr, kS2 = col.kS2, kL2 = col.kL2,
            kLevelId = col.kLevelId, kIsActive = col.kIsActive, kChannelFix = col.kChannelFix, kGamma = col.kGamma,
            kJ1 = col.kJ1, kPi1 = col.kPi1, kJ2 = col.kJ2, kPi2 = col.kPi2, kE2 = col.kE2, kM1 = col.kM1,
            kM2 = col.kM2, kZ1 = col.kZ1, kZ2 = col.kZ2, kEntranceSepE = col.kEntranceSepE, kSepE = col.kSepE,
            kJ3 = col.kJ3, kPi3 = col.kPi3, kE3 = col.kE3, kPType = col.kPType, kChRad = col.kChRad, kG1 = col.kG1,
            kG2 = col.kG2, kEcMultMask = col.kEcMultMask, kBindingE = col.kBindingE, kGammaIsRWA = col.kGammaIsRWA,
            kE1 = col.kE1;
  QVector<QString> tok(kColumns);
  // a number: the token read when it holds the same value, else the shortest form of the value
  auto number = [&ch, &tok](int c, double v) {
    const QString was = c < ch.fileTokens.size() ? ch.fileTokens.at(c) : QString();
    bool ok = false;
    const double old = was.toDouble(&ok);
    tok[c] = (ok && old == v) ? was : Shortest(v);
  };
  number(kLevelJ, level.jValue);
  number(kLevelPi, level.piValue);
  number(kLevelE, level.energy);
  number(kLevelFix, level.isFixed);
  number(kAa, 1);
  number(kIr, ch.pairIndex + 1);
  number(kS2, int(ch.sValue * 2));
  number(kL2, ch.lValue * 2);
  number(kLevelId, id);
  number(kIsActive, level.isActive);
  number(kChannelFix, ch.isFixed);
  number(kGamma, ch.reducedWidth);
  number(kJ1, pair.lightJ);
  number(kPi1, pair.lightPi);
  number(kJ2, pair.heavyJ);
  number(kPi2, pair.heavyPi);
  number(kE2, pair.excitationEnergy);
  number(kM1, pair.lightM);
  number(kM2, pair.heavyM);
  number(kZ1, pair.lightZ);
  number(kZ2, pair.heavyZ);
  number(kEntranceSepE, entrance.seperationEnergy);
  number(kSepE, pair.seperationEnergy);
  // j3 pi3 e3: not used by the engine; kept as read
  for (int c : {kJ3, kPi3}) tok[c] = c < ch.fileTokens.size() && !ch.fileTokens.at(c).isEmpty() ? ch.fileTokens.at(c) : "0";
  tok[kE3] = kE3 < ch.fileTokens.size() && !ch.fileTokens.at(kE3).isEmpty() ? ch.fileTokens.at(kE3) : "0.0";
  number(kPType, pair.pairType);
  number(kChRad, pair.channelRadius);
  number(kG1, pair.lightG);
  number(kG2, pair.heavyG);
  number(kEcMultMask, pair.ecMultMask);
  // the THM binding energy and the amplitude flag (columns 32, 33): as read, not edited here
  for (int c : {kBindingE, kGammaIsRWA}) tok[c] = c < ch.fileTokens.size() ? ch.fileTokens.at(c) : QString();
  // the light particle's excitation energy (column 34): written when set or when the row had it
  const bool hadE1 = kE1 < ch.fileTokens.size() && !ch.fileTokens.at(kE1).isEmpty();
  if (pair.lightExcitation != 0. || hadE1) number(kE1, pair.lightExcitation);
  azure2::front::LevelRow row;
  for (const QString &s : tok) row.tokens.push_back(s.toStdString());
  return row;
}

std::vector<azure2::front::LevelRow> LevelsTab::levelRows() {
  const QList<PairsData> pairs = pairsModel->getPairs();
  const QList<LevelsData> levels = levelsModel->getLevels();
  const QList<ChannelsData> channels = channelsModel->getChannels();
  // The levels in the order of the model: a file's levels in the file's order (they were added as they first
  // appeared), new ones after them, so that an edit moves no other row.  The members of an isospin level (one ID over
  // several nuclei) are written together.
  QList<int> order;
  {
    QSet<int> placed;
    for (int r = 0; r < levels.size(); r++) {
      if (placed.contains(r)) continue;
      order << r;
      placed.insert(r);
      if (levels.at(r).fileId <= 0) continue;
      for (int q = r + 1; q < levels.size(); q++)
        if (levels.at(q).fileId == levels.at(r).fileId && !placed.contains(q)) {
          order << q;
          placed.insert(q);
        }
    }
  }
  // The ID of a level: the one it was read with, else the next free one; the members of an isospin level (one ID
  // over several nuclei, writtenLevelOrder keeps them together) share it.
  QSet<int> used;
  for (const LevelsData &l : levels)
    if (l.fileId > 0) used.insert(l.fileId);
  int next = 1;
  std::vector<azure2::front::LevelRow> rows;
  for (int k = 0; k < order.size(); k++) {
    const int levelRow = order.at(k);
    const LevelsData &level = levels.at(levelRow);
    int id = level.fileId;
    if (id <= 0) {
      while (used.contains(next)) next++;
      id = next;
      used.insert(id);
    }
    for (const ChannelsData &ch : channels) {
      if (ch.levelIndex != levelRow || ch.pairIndex < 0 || ch.pairIndex >= pairs.size()) continue;
      const PairsData &pair = pairs.at(ch.pairIndex);
      rows.push_back(channelRow(ch, level, pair, pairs.at(0), id));
    }
    // the members of an isospin level are one level: no blank line between them
    if (k + 1 < order.size() && level.fileId > 0 && levels.at(order.at(k + 1)).fileId == level.fileId) continue;
    rows.push_back(azure2::front::LevelRow());
  }
  return rows;
}

void LevelsTab::addExtraPanel(const QString &title, QWidget *panel) {
  extras->addTab(panel, title);
  extras->setVisible(true);
}
