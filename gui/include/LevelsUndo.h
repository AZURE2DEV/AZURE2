#ifndef LEVELSUNDO_H
#define LEVELSUNDO_H

#include <functional>
#include <QList>
#include <QString>
#include <QUndoCommand>

#include "ChannelsModel.h"
#include "LevelsModel.h"
#include "PairsModel.h"

/*!
 * What the Particle Pairs and Levels tabs hold (docs/upgrade/qt_gui.md section 11.2): the three models, the channel
 * limits and the hard-sphere option, and the levels selected.  An entry of the tabs' undo stack is the state before
 * and after one action; the project text is not touched (the tabs write <levels> on a save or a run).
 */
struct LevelsState {
  QList<PairsData> pairs;
  QList<LevelsData> levels;
  QList<ChannelsData> channels;
  int maxL = 0, maxMult = 0, maxNumMult = 0;
  bool includeAllJPi = false;
  QList<int> selectedLevels;  ///< model rows; not compared
  /// The same pairs, levels, channels, limits and option (the selection is not part of the comparison).
  bool SameContent(const LevelsState &o) const;
};

/// One undoable action: the state before and after it; the first redo (when the command is pushed) does nothing, the
/// action is already done.
class LevelsUndoCommand : public QUndoCommand {
 public:
  LevelsUndoCommand(const QString &text, const LevelsState &before, const LevelsState &after,
                    std::function<void(const LevelsState &)> apply);
  void undo() override;
  void redo() override;

 private:
  LevelsState before_, after_;
  std::function<void(const LevelsState &)> apply_;
  bool first_ = true;
};

#endif
