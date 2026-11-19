// The undo stack of the Particle Pairs and Levels tabs (docs/upgrade/qt_gui.md section 11.2).
#include "LevelsUndo.h"

namespace {

bool SamePair(const PairsData &a, const PairsData &b) {
  return a.lightJ == b.lightJ && a.lightPi == b.lightPi && a.lightZ == b.lightZ && a.lightM == b.lightM &&
         a.lightG == b.lightG && a.heavyJ == b.heavyJ && a.heavyPi == b.heavyPi && a.heavyZ == b.heavyZ &&
         a.heavyM == b.heavyM && a.heavyG == b.heavyG && a.excitationEnergy == b.excitationEnergy &&
         a.seperationEnergy == b.seperationEnergy && a.channelRadius == b.channelRadius && a.pairType == b.pairType &&
         a.ecMultMask == b.ecMultMask && a.lightExcitation == b.lightExcitation;
}

bool SameLevel(const LevelsData &a, const LevelsData &b) {
  return a.isActive == b.isActive && a.isFixed == b.isFixed && a.jValue == b.jValue && a.piValue == b.piValue &&
         a.energy == b.energy && a.nucZ == b.nucZ && a.nucA == b.nucA && a.fileId == b.fileId;
}

bool SameChannel(const ChannelsData &a, const ChannelsData &b) {
  return a.isFixed == b.isFixed && a.levelIndex == b.levelIndex && a.pairIndex == b.pairIndex && a.sValue == b.sValue &&
         a.lValue == b.lValue && a.radType == b.radType && a.reducedWidth == b.reducedWidth &&
         a.fileTokens == b.fileTokens;
}

template <class T, class F>
bool SameList(const QList<T> &a, const QList<T> &b, F same) {
  if (a.size() != b.size()) return false;
  for (int i = 0; i < a.size(); i++)
    if (!same(a.at(i), b.at(i))) return false;
  return true;
}

}  // namespace

bool LevelsState::SameContent(const LevelsState &o) const {
  return maxL == o.maxL && maxMult == o.maxMult && maxNumMult == o.maxNumMult && includeAllJPi == o.includeAllJPi &&
         SameList(pairs, o.pairs, SamePair) && SameList(levels, o.levels, SameLevel) &&
         SameList(channels, o.channels, SameChannel);
}

LevelsUndoCommand::LevelsUndoCommand(const QString &text, const LevelsState &before, const LevelsState &after,
                                     std::function<void(const LevelsState &)> apply) :
  QUndoCommand(text), before_(before), after_(after), apply_(std::move(apply)) {}

void LevelsUndoCommand::undo() { apply_(before_); }

void LevelsUndoCommand::redo() {
  if (first_) {
    first_ = false;
    return;
  }
  apply_(after_);
}
