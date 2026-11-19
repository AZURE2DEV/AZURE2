// The Levels tab on the go (docs/upgrade/qt_gui.md section 11.2): the edits the tab checks (the add bar, a cell
// edited in place), removal of several levels with one confirmation, duplicates, bulk fix and free, background poles,
// the undo stack of the Pairs and Levels tabs, and the channels' quantities from the facade (section 11.4).
#include <QLocale>
#include <QMessageBox>
#include <QTimer>
#include <QToolTip>
#include <QUndoStack>
#include <algorithm>
#include <cmath>

#include "LevelsTab.h"
#include "ProjectHolder.h"

namespace {

/// A number as typed: the C locale first (the file's), then the user's.
bool ParseNumber(const QString &text, double &out) {
  bool ok = false;
  out = QLocale::c().toDouble(text.trimmed(), &ok);
  if (!ok) out = QLocale().toDouble(text.trimmed(), &ok);
  return ok && std::isfinite(out);
}

}  // namespace

// ---- the project, selection, status ----

ProjectHolder *LevelsTab::holder() const {
  QWidget *w = window();
  return w ? w->findChild<ProjectHolder *>() : nullptr;
}

QString LevelsTab::levelName(const LevelsData &level) const {
  return AddLevelDialog::spinParityLabel(level.jValue, level.piValue) + QString(" at %1 MeV").arg(level.energy);
}

void LevelsTab::setStatus(const QString &text, bool error) {
  status_->setText(error ? QString("<font color='#b00020'>%1</font>").arg(text.toHtmlEscaped()) : text.toHtmlEscaped());
}

void LevelsTab::showRefusal(const QModelIndex &index, const QString &why) {
  setStatus(tr("Refused: %1").arg(why), true);
  QTableView *view = index.model() == levelsModelProxy ? levelsView : channelsView;
  const QRect r = view->visualRect(index);
  QToolTip::showText(view->viewport()->mapToGlobal(r.bottomLeft()), why, view, r, 6000);
}

QList<int> LevelsTab::selectedLevelRows() const {
  QList<QModelIndex> shown = levelsView->selectionModel()->selectedRows();
  std::sort(shown.begin(), shown.end(), [](const QModelIndex &a, const QModelIndex &b) { return a.row() < b.row(); });
  QList<int> rows;
  for (const QModelIndex &i : shown) rows << levelsModelProxy->mapToSource(i).row();
  return rows;
}

int LevelsTab::currentLevelRow() const {
  const QModelIndex c = levelsView->currentIndex();
  if (c.isValid() && levelsView->selectionModel()->isRowSelected(c.row(), QModelIndex()))
    return levelsModelProxy->mapToSource(c).row();
  const QList<int> rows = selectedLevelRows();
  return rows.isEmpty() ? -1 : rows.first();
}

void LevelsTab::selectLevelRows(const QList<int> &rows) {
  QItemSelection selection;
  QModelIndex first;
  for (int r : rows) {
    if (r < 0 || r >= levelsModel->rowCount(QModelIndex())) continue;
    const QModelIndex p = levelsModelProxy->mapFromSource(levelsModel->index(r, 4));
    if (!p.isValid()) continue;  // filtered out
    selection.select(p.sibling(p.row(), 0), p.sibling(p.row(), levelsModelProxy->columnCount() - 1));
    if (!first.isValid()) first = p;
  }
  QItemSelectionModel *sm = levelsView->selectionModel();
  if (!first.isValid()) {
    sm->clearSelection();
    currentLevelChanged();
    return;
  }
  sm->setCurrentIndex(first, QItemSelectionModel::NoUpdate);
  sm->select(selection, QItemSelectionModel::ClearAndSelect);
  levelsView->scrollTo(first);
  currentLevelChanged();
}

// ---- the undo stack ----

LevelsState LevelsTab::currentState() const {
  LevelsState s;
  if (pairsModel) s.pairs = pairsModel->getPairs();
  s.levels = levelsModel->getLevels();
  s.channels = channelsModel->getChannels();
  s.maxL = maxLSpin->value();
  s.maxMult = maxMultSpin->value();
  s.maxNumMult = maxNumMultSpin->value();
  s.includeAllJPi = includeAllJPiCheck->isChecked();
  s.selectedLevels = selectedLevelRows();
  return s;
}

void LevelsTab::undoLater() {
  if (loading_ || restoring_ || undoPending_) return;
  undoPending_ = true;
  QTimer::singleShot(0, this, &LevelsTab::recordUndo);
}

void LevelsTab::recordUndo() {
  undoPending_ = false;
  if (loading_ || restoring_) return;
  const LevelsState now = currentState();
  if (now.SameContent(lastState_)) return;
  const QString text = undoText_.isEmpty() ? tr("edit") : undoText_;
  undoText_.clear();
  undo_->push(new LevelsUndoCommand(text, lastState_, now, [this](const LevelsState &s) { applyState(s); }));
  lastState_ = now;
}

void LevelsTab::applyState(const LevelsState &state) {
  restoring_ = true;
  // the pairs first (the nucleus column follows them), then the levels and their channels
  if (pairsModel) pairsModel->setPairs(state.pairs);
  levelsModel->setLevels(state.levels);
  channelsModel->setChannels(state.channels);
  for (QSpinBox *b : {maxLSpin, maxMultSpin, maxNumMultSpin}) b->blockSignals(true);
  maxLSpin->setValue(state.maxL);
  maxMultSpin->setValue(state.maxMult);
  maxNumMultSpin->setValue(state.maxNumMult);
  for (QSpinBox *b : {maxLSpin, maxMultSpin, maxNumMultSpin}) b->blockSignals(false);
  includeAllJPiCheck->blockSignals(true);
  includeAllJPiCheck->setChecked(state.includeAllJPi);
  includeAllJPiCheck->blockSignals(false);
  levelsView->resizeRowsToContents();
  updateMissingSpinParities();
  updateJpiFilter();
  selectLevelRows(state.selectedLevels);
  lastState_ = currentState();
  restoring_ = false;
}

void LevelsTab::finishLoading() {
  loading_ = false;
  undoPending_ = false;
  undo_->clear();
  updateJpiFilter();
  updateCount();
  // the first level selected: its channels are on screen at once (docs/upgrade/qt_gui.md section 11.2)
  if (selectedLevelRows().isEmpty() && levelsModelProxy->rowCount() > 0)
    selectLevelRows({levelsModelProxy->mapToSource(levelsModelProxy->index(0, 0)).row()});
  lastState_ = currentState();
}

// ---- the edits the tab checks ----

bool LevelsTab::addLevelFromBar(QString &why) {
  if (!pairsModel || pairsModel->getPairs().isEmpty()) {
    why = tr("there is no particle pair yet: add one in the tab Particle Pairs first");
    return false;
  }
  double energy = 0;
  if (!ParseNumber(addE_->text(), energy)) {
    why = addE_->text().trimmed().isEmpty() ? tr("type the level's energy E_x (MeV)")
                                            : tr("'%1' is not a number: the energy E_x in MeV").arg(addE_->text());
    return false;
  }
  LevelsData level{1, 0, addJ_->value(), addPi_->parity(), energy};
  const QList<Nucleus> nuclei = pairsModel->nuclei();
  if (nuclei.size() > 1 && addNucleus_->currentIndex() >= 0 && addNucleus_->currentIndex() < nuclei.size()) {
    level.nucZ = nuclei.at(addNucleus_->currentIndex()).Z;
    level.nucA = nuclei.at(addNucleus_->currentIndex()).A;
  } else if (nuclei.size() == 1) {
    level.nucZ = nuclei.first().Z;
    level.nucA = nuclei.first().A;
  }
  if (calculateChannels(level, -1).isEmpty()) {
    why = tr("%1 has no channel with these pairs and limits (maximum l %2, multipolarity %3): raise a limit or choose "
             "another J^pi").arg(AddLevelDialog::spinParityLabel(level.jValue, level.piValue)).arg(maxLSpin->value())
              .arg(maxMultSpin->value());
    return false;
  }
  const QList<LevelsData> levels = levelsModel->getLevels();
  for (const LevelsData &l : levels)
    if (l.jValue == level.jValue && l.piValue == level.piValue && std::fabs(l.energy - energy) < 1.e-3 &&
        (!l.hasNucleus() || !level.hasNucleus() || l.nucleus() == level.nucleus())) {
      why = tr("%1 exists (the engine takes two levels of one J^pi within 1 keV as one)").arg(levelName(l));
      return false;
    }
  undoText_ = tr("add level %1").arg(levelName(level));
  addLevel(level, false);
  selectLevelRows({(int)levelsModel->getLevels().size() - 1});
  setStatus(tr("Added %1 with %2 channels.").arg(levelName(level)).arg(calculateChannels(level, -1).size()));
  return true;
}

bool LevelsTab::setLevelEnergy(int row, const QString &text, QString &why) {
  const QList<LevelsData> levels = levelsModel->getLevels();
  if (row < 0 || row >= levels.size()) return false;
  double energy = 0;
  if (!ParseNumber(text, energy)) {
    why = tr("'%1' is not a number: the energy E_x in MeV").arg(text);
    return false;
  }
  const LevelsData &level = levels.at(row);
  if (energy == level.energy) return true;
  for (int i = 0; i < levels.size(); i++) {
    const LevelsData &l = levels.at(i);
    if (i != row && l.jValue == level.jValue && l.piValue == level.piValue && std::fabs(l.energy - energy) < 1.e-3 &&
        l.nucleus() == level.nucleus()) {
      why = tr("%1 exists (the engine takes two levels of one J^pi within 1 keV as one)").arg(levelName(l));
      return false;
    }
  }
  undoText_ = tr("energy of %1").arg(levelName(level));
  levelsModel->setData(levelsModel->index(row, 4), energy, Qt::EditRole);
  setStatus(tr("%1: E_x = %2 MeV.").arg(AddLevelDialog::spinParityLabel(level.jValue, level.piValue)).arg(energy));
  refreshChannelInfo();
  return true;
}

bool LevelsTab::setLevelJpi(int row, double j, int parity, QString &why) {
  const QList<LevelsData> levels = levelsModel->getLevels();
  if (row < 0 || row >= levels.size()) return false;
  const LevelsData &level = levels.at(row);
  if (j == level.jValue && parity == level.piValue) return true;
  LevelsData changed = level;
  changed.jValue = j;
  changed.piValue = parity;
  if (calculateChannels(changed, row).isEmpty()) {
    why = tr("%1 has no channel with these pairs and limits (maximum l %2): raise a limit or choose another J^pi")
              .arg(AddLevelDialog::spinParityLabel(j, parity)).arg(maxLSpin->value());
    return false;
  }
  for (int i = 0; i < levels.size(); i++) {
    const LevelsData &l = levels.at(i);
    if (i != row && l.jValue == j && l.piValue == parity && std::fabs(l.energy - level.energy) < 1.e-3 &&
        l.nucleus() == level.nucleus()) {
      why = tr("%1 exists").arg(levelName(l));
      return false;
    }
  }
  undoText_ = tr("J^pi of %1").arg(levelName(level));
  levelsModel->setData(levelsModel->index(row, 2), j, Qt::EditRole);
  levelsModel->setData(levelsModel->index(row, 3), parity, Qt::EditRole);
  // the channels of the new J^pi; those it keeps keep their widths and fix flags
  updateChannelsLevelEdited(row);
  updateMissingSpinParities();
  updateFilter(QItemSelection());
  setStatus(tr("%1 is now %2.").arg(levelName(level), AddLevelDialog::spinParityLabel(j, parity)));
  return true;
}

bool LevelsTab::setChannelValue(int row, const QString &text, QString &why) {
  const QList<ChannelsData> channels = channelsModel->getChannels();
  if (row < 0 || row >= channels.size()) return false;
  double value = 0;
  const azure2::front::ChannelQuantities q = channelQuantities(row);
  const QString what = QString::fromStdString(q.symbol + " [" + q.unit + "]");
  if (!ParseNumber(text, value)) {
    why = tr("'%1' is not a number: the value is %2, %3").arg(text, what, QString::fromStdString(q.meaning));
    return false;
  }
  if (value == channels.at(row).reducedWidth) return true;
  undoText_ = tr("value of a channel");
  channelsModel->setData(channelsModel->index(row, ChannelsModel::ValueColumn), value, Qt::EditRole);
  refreshChannelInfo();
  updateDetails(channelsView->selectionModel()->selection());
  return true;
}

bool LevelsTab::removeLevels(QList<int> rows, bool ask) {
  const QList<LevelsData> levels = levelsModel->getLevels();
  const QList<ChannelsData> channels = channelsModel->getChannels();
  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
  rows.erase(std::remove_if(rows.begin(), rows.end(), [&levels](int r) { return r < 0 || r >= levels.size(); }), rows.end());
  if (rows.isEmpty()) return false;
  if (ask) {
    QStringList names;
    for (int r : rows) {
      int n = 0;
      for (const ChannelsData &c : channels) n += c.levelIndex == r;
      if (names.size() < 15) names << tr("%1 (%2 channels)").arg(levelName(levels.at(r))).arg(n);
    }
    if (rows.size() > names.size()) names << tr("... and %1 more").arg(rows.size() - names.size());
    const auto answer = QMessageBox::question(this, tr("Remove Levels"),
                                              tr("Remove %n level(s) with their channels?", "", rows.size()) + "\n\n" +
                                                  names.join("\n") + "\n\n" + tr("Undo with Ctrl+Z."),
                                              QMessageBox::Yes | QMessageBox::Cancel, QMessageBox::Yes);
    if (answer != QMessageBox::Yes) return false;
  }
  // the neighbour of the first removed level (in the order shown) takes the selection
  const QModelIndex firstShown = levelsModelProxy->mapFromSource(levelsModel->index(rows.first(), 0));
  const int shownRow = firstShown.isValid() ? firstShown.row() : 0;
  undoText_ = rows.size() == 1 ? tr("remove %1").arg(levelName(levels.at(rows.first())))
                               : tr("remove %1 levels").arg(rows.size());
  levelsView->selectionModel()->clearSelection();
  for (int k = rows.size() - 1; k >= 0; k--) {
    levelsModel->removeRows(rows.at(k), 1, QModelIndex());
    updateChannelsLevelDeleted(rows.at(k));
  }
  updateMissingSpinParities();
  updateJpiFilter();
  const int count = levelsModelProxy->rowCount();
  if (count > 0)
    selectLevelRows({levelsModelProxy->mapToSource(levelsModelProxy->index(std::min(shownRow, count - 1), 0)).row()});
  setStatus(tr("Removed %n level(s); Ctrl+Z brings them back.", "", rows.size()));
  return true;
}

void LevelsTab::removeSelectedLevels() {
  const QList<int> rows = selectedLevelRows();
  if (!rows.isEmpty()) removeLevels(rows, true);
}

void LevelsTab::copyChannelValues(int from, int to) {
  const QList<ChannelsData> channels = channelsModel->getChannels();
  for (int i = 0; i < channels.size(); i++) {
    const ChannelsData &c = channels.at(i);
    if (c.levelIndex != to) continue;
    for (const ChannelsData &o : channels) {
      if (o.levelIndex != from || o.pairIndex != c.pairIndex || o.sValue != c.sValue || o.lValue != c.lValue ||
          o.radType != c.radType)
        continue;
      channelsModel->setData(channelsModel->index(i, 0), o.isFixed, Qt::EditRole);
      channelsModel->setData(channelsModel->index(i, ChannelsModel::ValueColumn), o.reducedWidth, Qt::EditRole);
      break;
    }
  }
}

void LevelsTab::duplicateSelectedLevels() {
  const QList<int> rows = selectedLevelRows();
  if (rows.isEmpty()) return;
  undoText_ = tr("duplicate %n level(s)", "", rows.size());
  QList<int> copies;
  for (int r : rows) {
    LevelsData copy = levelsModel->getLevels().at(r);
    copy.fileId = 0;
    // 10 keV above, or above that until no level of the J^pi is within 1 keV (the engine's rule for one level)
    bool clash = true;
    while (clash) {
      copy.energy = std::round((copy.energy + 0.01) * 1e9) / 1e9;
      clash = false;
      for (const LevelsData &l : levelsModel->getLevels())
        clash |= l.jValue == copy.jValue && l.piValue == copy.piValue && l.nucleus() == copy.nucleus() &&
                 std::fabs(l.energy - copy.energy) < 1.e-3;
    }
    addLevel(copy, false);
    const int to = (int)levelsModel->getLevels().size() - 1;
    copyChannelValues(r, to);
    copies << to;
  }
  selectLevelRows(copies);
  setStatus(tr("Duplicated %n level(s) 10 keV above; type the new energy.", "", rows.size()));
  const QModelIndex p = levelsModelProxy->mapFromSource(levelsModel->index(copies.first(), 4));
  if (p.isValid() && isVisible()) levelsView->edit(p);
}

void LevelsTab::setFixed(const QList<int> &rows, bool energies, bool fixed) {
  undoText_ = (fixed ? tr("fix ") : tr("free ")) + (energies ? tr("energies") : tr("widths"));
  if (energies) {
    for (int r : rows) levelsModel->setData(levelsModel->index(r, 1), fixed ? 1 : 0, Qt::EditRole);
    return;
  }
  const QList<ChannelsData> channels = channelsModel->getChannels();
  for (int i = 0; i < channels.size(); i++)
    if (rows.contains(channels.at(i).levelIndex)) channelsModel->setData(channelsModel->index(i, 0), fixed ? 1 : 0, Qt::EditRole);
}

void LevelsTab::addBackgroundPole() {
  const int row = currentLevelRow();
  if (row < 0) return;
  const QList<LevelsData> levels = levelsModel->getLevels();
  LevelsData pole = levels.at(row);
  double top = pole.energy;
  for (const LevelsData &l : levels)
    if (l.jValue == pole.jValue && l.piValue == pole.piValue && l.nucleus() == pole.nucleus()) top = std::max(top, l.energy);
  // a start value to be edited: 10 MeV above the J^pi's highest level, energy fixed, widths 0 and free
  pole.energy = std::round((top + 10.0) * 1e6) / 1e6;
  pole.isFixed = 1;
  pole.isActive = 1;
  pole.fileId = 0;
  undoText_ = tr("background pole %1").arg(levelName(pole));
  addLevel(pole, false);
  const int to = (int)levelsModel->getLevels().size() - 1;
  selectLevelRows({to});
  setStatus(tr("Added a background pole %1 (energy fixed, widths 0): set its energy and free the widths to fit.")
                .arg(levelName(pole)));
}

// ---- the channels' quantities (the facade) ----

azure2::front::ChannelQuantities LevelsTab::channelQuantities(int row) const {
  const QList<ChannelsData> channels = channelsModel->getChannels();
  const QList<LevelsData> levels = levelsModel->getLevels();
  const QList<PairsData> pairs = pairsModel ? pairsModel->getPairs() : QList<PairsData>();
  azure2::front::ChannelQuantities none;
  if (row < 0 || row >= channels.size()) return none;
  const ChannelsData &ch = channels.at(row);
  if (ch.levelIndex < 0 || ch.levelIndex >= levels.size() || ch.pairIndex < 0 || ch.pairIndex >= pairs.size()) return none;
  const LevelsData &level = levels.at(ch.levelIndex);
  const PairsData &pair = pairs.at(ch.pairIndex);
  const azure2::front::LevelRow r = channelRow(ch, level, pair, pairs.first(), std::max(1, level.fileId));
  double radius = 0;
  bool relativistic = false;
  if (ProjectHolder *h = holder()) {
    const QString rel = h->configValue("relativistic");
    relativistic = rel == "1" || rel == "true";
    if (pair.pairType == 0)
      radius = azure2::front::ChannelRadius(h->project(), ch.pairIndex + 1, ch.lValue, (int)std::lround(2 * ch.sValue),
                                            (int)std::lround(2 * level.jValue), pair.channelRadius)
                   .radius;
  }
  return azure2::front::ChannelAtLevel(r, radius, relativistic);
}

void LevelsTab::refreshChannelInfo() {
  const int level = currentLevelRow();
  const QList<ChannelsData> channels = channelsModel->getChannels();
  for (int i = 0; i < channels.size(); i++) {
    if (channels.at(i).levelIndex != level) continue;
    const azure2::front::ChannelQuantities q = channelQuantities(i);
    ChannelsModel::Info info;
    if (q.error.empty()) {
      info.meaning = QString::fromStdString(q.unit.empty() ? q.symbol : q.symbol + " [" + q.unit + "]");
      info.tip = QString::fromStdString(q.meaning);
      if (q.particle && std::isfinite(q.wignerRatio)) info.ratio = q.wignerRatio;
    } else {
      info.tip = QString::fromStdString(q.error);
    }
    channelsModel->setInfo(i, info);
  }
  channelsView->viewport()->update();
}
