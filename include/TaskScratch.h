#ifndef AZURE2_TASKSCRATCH_H
#define AZURE2_TASKSCRATCH_H

#include <memory>
#include <vector>

#include "ChannelSpace.h"

class AMatrixFunc;
class CNuc;
class Config;
class LevelState;
class Site;
class THMMatrixFunc;

namespace azure2 {

/*!
 * The scratch of the energy tasks one thread runs: the level-matrix solvers with their buffers, the channel-space
 * route's per-J-group values and the list of a task's units.  Owned by the caller's buffer (`Results::scratch`, one
 * per thread of a pass), so it lives as long as the evaluation's buffer and is never shared by two evaluations
 * (docs/upgrade/evaluation.md: no hidden state).  It holds no value an evaluation reads back: each solver is
 * re-pointed (Reset) at the task's compound nucleus and configuration before use, and the channel-space values are
 * keyed by the levels they were built from (ChannelSpace::Memo).  Before 2026-11-12 these were `thread_local`
 * objects of the passes.
 */
class TaskScratch {
 public:
  TaskScratch();
  ~TaskScratch();
  TaskScratch(const TaskScratch &) = delete;
  TaskScratch &operator=(const TaskScratch &) = delete;
  /// The point pass's solver (SolveAtOneEnergy), re-pointed at `compound`, the levels `state` and `configure`.
  AMatrixFunc &PointSolver(CNuc *compound, const LevelState &state, const Config &configure);
  /// The point pass's level-space solver of sigma_HOES, re-pointed likewise.
  THMMatrixFunc &PointThmSolver(CNuc *compound, const LevelState &state, const Config &configure);
  /// The gradient pass's solver (GradEnergyTasks.cpp), re-pointed likewise.
  AMatrixFunc &GradientSolver(CNuc *compound, const LevelState &state, const Config &configure);
  /// The gradient pass's THM solver (HoesTask, ThmPointModelGradient).  Not re-pointed: its caller does.
  THMMatrixFunc &GradientThmSolver(CNuc *compound, const LevelState &state, const Config &configure);
  /// The units of the task being solved (EnergyPass.cpp).
  std::vector<Site *> units;

 private:
  ChannelSpace::Memo channel_;
  std::unique_ptr<AMatrixFunc> point_, gradient_;
  std::unique_ptr<THMMatrixFunc> pointThm_, gradientThm_;
};

/*!
 * One TaskScratch per thread of a parallel pass, in a buffer that owns them.  A slot is made on its thread's first
 * request and written only by that thread (OpenMP thread number inside the pass's team; 0 outside a parallel region,
 * and inside an evaluation that runs within another parallel loop, whose own buffer it then is).  A copy of the
 * buffer starts empty: scratch carries nothing an evaluation needs.
 */
class TaskScratchSet {
 public:
  TaskScratchSet();
  ~TaskScratchSet();
  TaskScratchSet(const TaskScratchSet &);
  TaskScratchSet &operator=(const TaskScratchSet &);
  /// The calling thread's scratch.  Throws std::logic_error beyond kMaxThreads threads.
  TaskScratch &ForThisThread();
  static constexpr int kMaxThreads = 256;

 private:
  std::vector<std::unique_ptr<TaskScratch>> slots_;
};

}  // namespace azure2

#endif
