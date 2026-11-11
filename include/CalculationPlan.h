#ifndef AZURE2_CALCULATIONPLAN_H
#define AZURE2_CALCULATIONPLAN_H

#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "Results.h"

class WorkingCopy;
class Unit;
class CNuc;
class Config;

namespace azure2 {

/*!
 * The calculation plan of a data set (docs/upgrade/data_model.md section 2.4, step S5): who calculates which point.
 * One task per point that calculates itself (a host); the points mapped onto it (WorkingCopy::MapData: the same entrance,
 * energies, target effect and sub-point grid, and for the components of a SUM, RATIO or total-capture segment the
 * base point's level matrix) are the host's to calculate, with their sub-points, in the order listed -- the
 * component points grouped by exit, so each exit's T-matrix is built once (Unit::Calculate).  A point holds the
 * index of the task that calculates it (Unit::PlanTask), never a pointer to another point; the plan names points
 * by their index in the data set's units (PointRef, WorkingCopy::Unit; the points keep their indices after the set-up), so
 * a copy of the data set copies its plan as it is.
 *
 * What a task writes: the fixed quantities and the angular-distribution coefficients of the record (Results.h) of
 * every point it calculates and of every sub-point of those points; nothing else.  The natural-term slots of a
 * record belong to the model pass (ModelPass.h).  What a task reads: the set-up (tables, grids, kernels) and the
 * slots it wrote itself (a fold reads its sub-points; an automatic tolerance its probes).  Built when the mapping is
 * built (at set-up, and after a change of an energy shift or of a sub-point grid, WorkingCopy::FillEnergyShiftsFromParams),
 * and checked when built: every point of the data set is calculated by exactly one task, a mapped point is not a host,
 * and a host with sub-points maps only points with as many sub-points (Check).  Under AZURE_PLAN_CHECK=1 the point
 * pass also checks at run time that each task writes only its own slots and no slot is written twice
 * (WorkingCopy::CalculatePointsInParallel).
 */
using PointRef = int;  ///< a point by its index in the data set's units (WorkingCopy::Unit)

struct CalculationTask {
  PointRef host;                ///< the point that calculates itself and the points below
  std::vector<PointRef> mapped;  ///< the points the host calculates after itself, in this order
};

using UnitRef = int;  ///< a unit (a point or a sub-point) by its index in the data set's units (WorkingCopy::UnitAt)

/*!
 * One compound energy of the point pass (docs/upgrade/energy_grid.md section 8, first branch): every unit whose level
 * matrix is solved at that energy -- a point without a sub-point fold, or a sub-point of one -- with the same compound
 * nucleus, the same route (A-matrix, THM, R-matrix) and the same bits of E_x = AMatrixFunc::InEnergy.  The level
 * matrix is filled and factored once, with the first unit; the units are listed by T-matrix group -- (entrance, exit)
 * pair and the source of the external-capture amplitudes (ECAmplitudeSource) -- so each group's T-matrix is built
 * once, from its first unit; then each unit's observable is formed in turn
 * (SolveAtOneEnergy, PointPass.cpp).  AZURE_PLAN_CHECK=1 shows that this gives every unit the level matrix, T-matrix
 * and record its own solve gives, bit for bit (CheckEnergyTasks).
 */
struct EnergyTask {
  int nucleus = 0;             ///< compound nucleus of the units (PPair::GetNucleus of their entrance pairs, 1-based)
  char route = 'A';            ///< 'A' A-matrix, 'T' THM (THMMatrixFunc), 'R' R-matrix (RMatrixFunc)
  double energy = 0.0;         ///< E_x, MeV (AMatrixFunc::InEnergy), the same bits for every unit
  int numTMatrices = 0;        ///< T-matrix groups: the T-matrices the task builds
  int begin = 0;               ///< the task's units are [begin, begin + size) of the arrays of its EnergyTaskList
  int size = 0;
};

/*!
 * The energy tasks of a pass with their units in four parallel arrays, the units of one task side by side (a task is a
 * range of them, EnergyTask::begin and size): no allocation per task, so a rebuild (every evaluation that moves an
 * energy shift or a sub-point grid) costs one sort of the units, and a task reads its units' references from one
 * stretch of memory.
 */
struct EnergyTaskList {
  std::vector<EnergyTask> tasks;   ///< by decreasing cost (1 + T-matrices + units), ties in energy order
  std::vector<UnitRef> units;      ///< in the order calculated: by T-matrix group, then in plan order
  std::vector<int> tGroups;        ///< the T-matrix group of each unit (0, 1, ... inside its task; non-decreasing)
  std::vector<UnitRef> providers;  ///< the unit whose T-matrix each unit's observable is formed from (TaskUnits)
  std::vector<int> owners;         ///< the plan task (CalculationTask) of each unit: its point's
  int size() const { return (int)tasks.size(); }
  const EnergyTask &operator[](int k) const { return tasks[k]; }
};

/*!
 * The source of the external-capture amplitudes the T-matrix of `provider` is built from (GetECAmplitudeWithShift):
 * 0 without external capture; with the EC table 1 for amplitudes computed at the unit's energy (one function of the
 * energy), else the hash of the amplitudes the unit holds there (the set-up's: points at one energy need not hold the
 * same bits; Unit::ECAmplitudesHeldAt); without the table 2 + the segment key for the interpolation of
 * ECAmplitudeCache.  The units of a T-matrix group share it (a hash collision is refused when the groups are built,
 * Unit::SameECAmplitudesAs), so the T-matrix built from the first is the one each would build itself
 * (AZURE_PLAN_CHECK).
 */
long long ECAmplitudeSource(const Unit &provider, const Config &configure);

class CalculationPlan {
 public:
  /// The plan of `data` from its mapping (WorkingCopy::MapData): one task per point that is not mapped, in data order (the
  /// segments' own points, then the leaves' (WorkingCopy::PlanUnits)), and each (host, mapped) pair of `mapping` appended to its host's task in
  /// the order given; every point's task index set (Unit::SetPlanTask).  Then Check: a plan that breaks it throws
  /// std::logic_error.
  void Build(WorkingCopy &data, const std::vector<std::pair<PointRef, PointRef>> &mapping);
  int NumTasks() const { return (int)tasks_.size(); }
  const CalculationTask &Task(int task) const { return tasks_[task]; }
  /// The host point of task `task` in `data`.
  Unit *Host(WorkingCopy &data, int task) const;

  /// The build-time assertion (see above), on the data set the plan was built for: an empty string when it holds,
  /// otherwise what is wrong (the first point that breaks it, by segment key and lab energy).
  std::string Check(WorkingCopy &data) const;

  /*!
   * The energy tasks of the plan tasks `work` (indices of Task, in the order given): the units of each task's points --
   * a point itself, or its sub-points when it carries a sub-point fold -- grouped by (nucleus, route, E_x bits), and
   * inside a task into T-matrix groups by (entrance, exit, ECAmplitudeSource of the unit's provider), each in plan
   * order.  Kept between passes and rebuilt only when the work, a unit, a unit's E_x or an amplitude source has
   * changed (the energies move with shifts and re-placed sub-point grids, a copy drops the set-up's amplitudes; the
   * check walks the units, it does not sort).  Ordered by decreasing cost (1 + T-matrices + units), so the largest
   * tasks start first on a dynamic schedule.
   */
  const EnergyTaskList &EnergyTasks(WorkingCopy &data, CNuc *cnuc, const Config &configure,
                                    const std::vector<int> &work);
  /// The units of the plan task `task` in plan order (the host, then each mapped point; a point with a sub-point
  /// fold by its sub-points), appended to `units`; and to `providers`, if given, the unit each one's T-matrix comes
  /// from where a host calculates its mapped points (Unit::Calculate): the first unit at the same sub-point index
  /// (or point) of the same (entrance, exit) pair -- itself for the host, its host's for a mapped point of the host's
  /// exit.
  void TaskUnits(WorkingCopy &data, int task, std::vector<UnitRef> &units,
                 std::vector<UnitRef> *providers = nullptr) const;

  /// AZURE_PLAN_CHECK=1: the point pass checks each task's writes at run time (slow; a debugging aid).
  static bool RunTimeCheckRequested();
  /// The point pass under AZURE_PLAN_CHECK=1: steps 0 .. numSteps - 1 (the energy tasks, then the folds) run one
  /// after the other through `run`, which calculates step k and gives the slot ranges [begin, end) it may write; every
  /// slot that changes must lie in those ranges and change in one step only, and no two steps' ranges overlap.  The
  /// slots are compared after each step, or after each batch of steps in a pass of more than 2,000 (see the .cpp).
  /// A breach throws std::logic_error naming the step (`name`) and the slot.
  static void RunChecked(Results &results, int numSteps,
                         const std::function<void(int, std::vector<std::pair<std::size_t, std::size_t>> &)> &run,
                         const std::function<std::string(int)> &name);
  /// The slots of `u`'s own record apart from the natural terms (the model pass's), appended to `ranges`.
  static void AddRecordRange(const Unit &u, std::vector<std::pair<std::size_t, std::size_t>> &ranges);

 private:
  std::vector<CalculationTask> tasks_;
  // The energy tasks of the last pass and what they were built from: the work, every unit with its E_x.
  EnergyTaskList energy_tasks_;
  std::vector<int> energy_work_;
  std::vector<UnitRef> energy_units_;
  std::vector<double> energy_keys_;
  std::vector<long long> energy_sources_;
};

/// One energy task's calculation (EnergyPass.cpp): the level matrix filled and factored with units[0], then for each
/// unit in turn a T-matrix built from it when it starts a T-matrix group, and its observable into `results`.  The
/// units share the nucleus, route and E_x bits of units[0] (an EnergyTask, or a host and its mapped points).  The
/// groups are `tGroups` (one per unit, non-decreasing), or, if null, runs of one (entrance, exit) pair.  `progress`,
/// if given, holds the index of the unit being calculated (so a caller that catches an exception knows whose it was:
/// the fill's and factorisation's are units[0]'s).
void SolveAtOneEnergy(CNuc *cnuc, const Config &configure, Results &results, Unit *const *units, int n,
                      const int *tGroups = nullptr, int *progress = nullptr);
/// AZURE_PLAN_CHECK=1: every unit of `tasks` solved alone (fill and factorisation with it, the T-matrix of its
/// provider, its record) and compared bit for bit with the shared solve of its task (the bilinear forms g^T A g of
/// every J-group, the route, the T-matrix, the record), written into a copy of `results`.  Prints the counts and the first differences to stderr; returns
/// the number of units that differ.
long CheckEnergyTasks(WorkingCopy &data, CNuc *cnuc, const Config &configure, const Results &results,
                      const EnergyTaskList &tasks);

/// Whether `point` carries a sub-point fold (a target effect whose kernel is evaluated on sub-points): its sub-points
/// are the units of the energy pass, and the fold runs after it.
bool HasSubPointFold(const WorkingCopy &data, const Unit &point);

/// The points a host calculates besides itself, as its data set's plan lists them (empty for a point that is not
/// the host of a task: a mapped point, a sub-point, a probe copy, a point outside a data set).
class MappedPoints {
 public:
  MappedPoints() = default;
  MappedPoints(WorkingCopy *data, const std::vector<PointRef> *refs) : data_(data), refs_(refs) {}
  int size() const { return refs_ ? (int)refs_->size() : 0; }
  /// The i-th mapped point (0-based).
  Unit *operator[](int i) const;

 private:
  WorkingCopy *data_ = nullptr;
  const std::vector<PointRef> *refs_ = nullptr;
};

}  // namespace azure2

#endif
