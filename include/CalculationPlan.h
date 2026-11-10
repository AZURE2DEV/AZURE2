#ifndef AZURE2_CALCULATIONPLAN_H
#define AZURE2_CALCULATIONPLAN_H

#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "Results.h"

class EData;
class EPoint;

namespace azure2 {

/*!
 * The calculation plan of a data set (docs/upgrade/data_model.md section 2.4, step S5): who calculates which point.
 * One task per point that calculates itself (a host); the points mapped onto it (EData::MapData: the same entrance,
 * energies, target effect and sub-point grid, and for the components of a SUM, RATIO or total-capture segment the
 * base point's level matrix) are the host's to calculate, with their sub-points, in the order listed -- the
 * component points grouped by exit, so each exit's T-matrix is built once (EPoint::Calculate).  A point holds the
 * index of the task that calculates it (EPoint::PlanTask), never a pointer to another point; the plan holds points by
 * position (PointRef) and is resolved in the data set that owns it, so a copy of the data set copies its plan as it
 * is.
 *
 * What a task writes: the fixed quantities and the angular-distribution coefficients of the record (Results.h) of
 * every point it calculates and of every sub-point of those points; nothing else.  The natural-term slots of a
 * record belong to the model pass (ModelPass.h).  What a task reads: the set-up (tables, grids, kernels) and the
 * slots it wrote itself (a fold reads its sub-points; an automatic tolerance its probes).  Built when the mapping is
 * built (at set-up, and after a change of an energy shift or of a sub-point grid, EData::FillEnergyShiftsFromParams),
 * and checked when built: every point of the data set is calculated by exactly one task, a mapped point is not a host,
 * and a host with sub-points maps only points with as many sub-points (Check).  Under AZURE_PLAN_CHECK=1 the point
 * pass also checks at run time that each task writes only its own slots and no slot is written twice
 * (EData::CalculatePointsInParallel).
 */
struct PointRef {
  enum List : int { kSegments = 0, kComponents = 1 };
  int list = kSegments;  ///< which list of the data set: its data/test segments, or its component segments
  int segment = 0;       ///< 0-based position of the segment in that list
  int point = 0;         ///< 0-based position of the point in its segment
};

struct CalculationTask {
  PointRef host;                ///< the point that calculates itself and the points below
  std::vector<PointRef> mapped;  ///< the points the host calculates after itself, in this order
};

class CalculationPlan {
 public:
  /// The plan of `data` from its mapping (EData::MapData): one task per point that is not mapped, in data order (the
  /// segments, then the component segments), and each (host, mapped) pair of `mapping` appended to its host's task in
  /// the order given; every point's task index set (EPoint::SetPlanTask).  Then Check: a plan that breaks it throws
  /// std::logic_error.
  void Build(EData &data, const std::vector<std::pair<PointRef, PointRef>> &mapping);
  int NumTasks() const { return (int)tasks_.size(); }
  const CalculationTask &Task(int task) const { return tasks_[task]; }

  /// The build-time assertion (see above), on the data set the plan was built for: an empty string when it holds,
  /// otherwise what is wrong (the first point that breaks it, by segment key and lab energy).
  std::string Check(EData &data) const;
  /// The slot ranges [begin, end) of Results::slot that task `task` writes, for the data set's current layout.
  void WriteRanges(EData &data, int task, std::vector<std::pair<std::size_t, std::size_t>> &ranges) const;

  /// AZURE_PLAN_CHECK=1: the point pass checks each task's writes at run time (slow; a debugging aid).
  static bool RunTimeCheckRequested();
  /// The point pass under AZURE_PLAN_CHECK=1: each host of `hosts` (in this order, serially) runs its task through
  /// `run`; every slot that changes must lie in that task's write ranges and change in one task only.  A breach
  /// throws std::logic_error naming the task, the point and the slot.
  static void RunChecked(EData &data, const std::vector<EPoint *> &hosts, Results &results,
                         const std::function<void(EPoint *)> &run);

 private:
  std::vector<CalculationTask> tasks_;
};

/// The points a host calculates besides itself, as its data set's plan lists them (empty for a point that is not
/// the host of a task: a mapped point, a sub-point, a probe copy, a point outside a data set).
class MappedPoints {
 public:
  MappedPoints() = default;
  MappedPoints(EData *data, const std::vector<PointRef> *refs) : data_(data), refs_(refs) {}
  int size() const { return refs_ ? (int)refs_->size() : 0; }
  /// The i-th mapped point (0-based).
  EPoint *operator[](int i) const;

 private:
  EData *data_ = nullptr;
  const std::vector<PointRef> *refs_ = nullptr;
};

}  // namespace azure2

#endif
