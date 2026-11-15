// The point pass of an evaluation (WorkingCopy::CalculatePointsInParallel); the calculation at one compound energy (CalculationPlan.h, EnergyTask): one fill and factorisation of the level
// matrix, one T-matrix per group (entrance, exit, amplitude source), the observable of every unit; and the check (AZURE_PLAN_CHECK=1) that
// this shared solve gives every unit what its own solve gives, bit for bit.
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <set>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AMatrixFunc.h"
#include "CompoundSystem.h"
#include "CalculationPlan.h"
#include "Config.h"
#include "JGroup.h"
#include "NodeTable.h"
#include "Profile.h"
#include "Segment.h"
#include "RMatrixFunc.h"
#include "Results.h"
#include "THMMatrixFunc.h"
#include "Unit.h"
#include "WorkingCopy.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace azure2 {

namespace {

// Whether unit k starts another T-matrix group than unit `before`, the last one that read the T-matrix: by
// `tGroups`, or, if null, by a change of (entrance, exit) pair.
bool StartsTMatrix(Site *const *units, int before, int k, const int *tGroups) {
  if (tGroups) return tGroups[k] != tGroups[before];
  return units[k]->GetEntranceKey() != units[before]->GetEntranceKey() ||
         units[k]->GetExitKey() != units[before]->GetExitKey();
}

}  // namespace

void SolveAtOneEnergy(const CompoundSystem *compound, const LevelState &state, const std::vector<double> &thmCoherentValues, const Config &configure,
                      Results &results, Site *const *units, int n, const int *tGroups, int *progress) {
  if (progress) *progress = 0;
  if (n <= 0) return;
  // The readers of the T-matrix (every kind but kHoes) and of the level matrix's elements (kHoes, the HOES formula).
  int first = -1, firstHoes = -1;
  for (int k = 0; k < n; k++) {
    int &at = units[k]->Observable() == SiteObservable::kHoes ? firstHoes : first;
    if (at < 0) at = k;
  }
  // The solvers are the caller's scratch for this thread (Results::scratch), reused across the energies of its tasks:
  // their per-J-group buffers are sized once, and they are re-pointed at this compound and configuration first, so
  // nothing of another evaluation is read; the R-matrix route builds its own.
  TaskScratch &scratch = results.scratch.ForThisThread();
  AMatrixFunc *amatrix = nullptr;
  std::unique_ptr<RMatrixFunc> rmatrix;
  GenMatrixFunc *func = nullptr;
  if (first >= 0) {
    if (configure.paramMask & Config::USE_AMATRIX) {
      amatrix = &scratch.PointSolver(compound, state, configure);
      func = amatrix;
    } else {
      rmatrix.reset(new RMatrixFunc(compound, state, configure));
      func = rmatrix.get();
    }
    if (progress) *progress = first;
    func->ClearMatrices();
    func->FillMatrices(units[first]);
    func->InvertMatrices();
  }
  // sigma_HOES reads the level matrix's elements in level space: the fill above where it took no channel route (the
  // same matrices bit for bit), else a level-space fill of its own (THMMatrixFunc.h).
  AMatrixFunc *levels = nullptr;
  if (firstHoes >= 0) {
    if (amatrix && func == amatrix && amatrix->LevelSpaceThroughout()) {
      levels = amatrix;
    } else {
      THMMatrixFunc &thm = scratch.PointThmSolver(compound, state, configure);
      if (progress) *progress = firstHoes;
      thm.ClearMatrices();
      thm.FillMatrices(units[firstHoes]);
      thm.InvertMatrices();
      levels = &thm;
    }
  }
  // The first T-matrix after a fill is built directly, every later one through RecalculateTMatrix (which drops the
  // elements held).
  int last = -1;  // the last unit that read the T-matrix
  for (int k = 0; k < n; k++) {
    if (progress) *progress = k;
    if (units[k]->Observable() == SiteObservable::kHoes) {
      HoesCrossSection(*levels, units[k], thmCoherentValues, results);
      continue;
    }
    if (last < 0)
      func->CalculateTMatrix(units[k]);
    else if (StartsTMatrix(units, last, k, tGroups))
      func->RecalculateTMatrix(units[k]);
    last = k;
    func->CalculateCrossSection(units[k], results);
  }
}

namespace {

std::string Describe(const Site &u, const EnergyTask &task) {
  std::ostringstream s;
  s.precision(17);
  s << "segment key " << u.GetSegmentKey() << ", entrance " << u.GetEntranceKey() << ", exit " << u.GetExitKey()
    << ", E_cm " << u.GetCMEnergy() << " MeV (E_x " << task.energy << " MeV, nucleus " << task.nucleus << ", route "
    << task.route << ")";
  return s.str();
}

}  // namespace

long CheckEnergyTasks(WorkingCopy &data, const CompoundSystem *compound, const LevelState &state, const Config &configure, const Results &results,
                      const EnergyTaskList &tasks) {
  Results scratch = results;  // the check writes here, never into the evaluation's values
  std::unique_ptr<AMatrixFunc> ownA;
  if (configure.paramMask & Config::USE_AMATRIX) ownA.reset(new AMatrixFunc(compound, state, configure));
  std::unique_ptr<THMMatrixFunc> ownT(new THMMatrixFunc(compound, state, configure));
  long units = 0, shared = 0, differ = 0;
  long hoesTasks = 0, mixedTasks = 0, oneFill = 0;  // tasks with THM sites, with other readers too, solved by one fill
  std::vector<std::string> reports;
  std::vector<Site *> sites;
  std::vector<std::vector<double>> sharedRecords;
  for (const EnergyTask &task : tasks.tasks) {
    const int n = task.size;
    const UnitRef *taskUnits = tasks.units.data() + task.begin;
    units += n;
    if (n > 1) shared += n;
    // The task as the pass solves it (SolveAtOneEnergy), the records kept; then each unit's own solve.
    sites.resize((std::size_t)n);
    for (int k = 0; k < n; k++) sites[(std::size_t)k] = &data.SiteAt(taskUnits[k]);
    SolveAtOneEnergy(compound, state, data.ThmCoherentValues(), configure, scratch, sites.data(), n,
                     tasks.tGroups.data() + task.begin, nullptr);
    int hoes = 0, direct = -1;
    for (int k = 0; k < n; k++)
      if (sites[(std::size_t)k]->Observable() == SiteObservable::kHoes)
        hoes++;
      else if (direct < 0)
        direct = k;
    if (hoes > 0) hoesTasks++;
    if (hoes > 0 && direct >= 0) {
      mixedTasks++;
      if (ownA) {
        ownA->ClearMatrices();
        ownA->FillMatrices(sites[(std::size_t)direct]);
        if (ownA->LevelSpaceThroughout()) oneFill++;
      }
    }
    sharedRecords.assign((std::size_t)n, {});
    for (int k = 0; k < n; k++) {
      const Site *u = sites[(std::size_t)k];
      if (u->ResultSlot() == Unit::kNoRecord) continue;
      const std::size_t slot = u->ResultSlot(), size = Results::RecordSize(u->ResultAngularSlots(), 0);
      sharedRecords[(std::size_t)k].assign(scratch.slot.begin() + slot, scratch.slot.begin() + slot + size);
    }
    for (int k = 0; k < n; k++) {
      Site *u = sites[(std::size_t)k];
      std::string why;
      if (AMatrixFunc::InEnergy(compound, u) != task.energy) why = "E_x is not the task's";
      // The unit's own solve: the HOES formula in level space, every other kind through its own T-matrix.
      std::unique_ptr<RMatrixFunc> ownR;
      GenMatrixFunc *o = ownA.get();
      if (u->Observable() == SiteObservable::kHoes)
        o = ownT.get();
      else if (!o)
        ownR.reset(new RMatrixFunc(compound, state, configure)), o = ownR.get();
      o->ClearMatrices();
      o->FillMatrices(u);
      o->InvertMatrices();
      if (o == ownT.get()) {
        HoesCrossSection(*ownT, u, data.ThmCoherentValues(), scratch);
      } else {
        o->CalculateTMatrix(u);
        o->CalculateCrossSection(u, scratch);
      }
      const std::vector<double> &mine = sharedRecords[(std::size_t)k];
      if (why.empty() && u->ResultSlot() != Unit::kNoRecord &&
          std::memcmp(mine.data(), scratch.slot.data() + u->ResultSlot(), mine.size() * sizeof(double)) != 0)
        why = "record";
      if (why.empty()) continue;
      differ++;
      if (reports.size() < 10) reports.push_back(why + " differs at unit " + std::to_string(taskUnits[k]) + " (" +
                                                 Describe(*u, task) + ") from its own solve");
    }
  }
  std::cerr << "AZURE_PLAN_CHECK: " << units << " units in " << tasks.size() << " energy tasks (" << shared
            << " units at a shared energy) against their own solves: " << differ << " differ" << std::endl;
  if (hoesTasks > 0)
    std::cerr << "AZURE_PLAN_CHECK: THM sites in " << hoesTasks << " energy tasks, " << mixedTasks
              << " of them with other readers, " << oneFill << " of those solved by one fill" << std::endl;
  for (const std::string &r : reports) std::cerr << "AZURE_PLAN_CHECK:   " << r << std::endl;
  return differ;
}

}  // namespace azure2

std::vector<int> WorkingCopy::PassWork(const std::vector<const Segment *> &segments) const {
  // The plan tasks of the segments (one per point), in data order.
  std::vector<int> work;
  std::set<const Segment *> seen;
  for (const Segment *seg : segments) {
    if (!seg || !seen.insert(seg).second) continue;
    for (int t = 0; t < seg->NumTerms(); t++)
      for (const Unit &pt : seg->TermPoints(*this, t))
        work.push_back(pt.PlanTask());
  }
  return work;
}

const azure2::EnergyTaskList &WorkingCopy::EnergyTasksOf(const std::vector<const Segment *> &segments,
                                                         const Config &configure) {
  return plan_.EnergyTasks(*this, configure, PassWork(segments));
}

/*!
 * One pass over the points of an evaluation, shared by the CLI (AZURECalc::PrecalculatePoints) and the API's
 * calculations (AZUREAPI): parallel over the compound energies, or serial inside a caller's parallel region.
 */

int WorkingCopy::CalculatePointsInParallel(const std::vector<const Segment *> &segments, const CompoundSystem *compound,
                                           const LevelState &state, const Config &configure,
                                           azure2::Results &results) {
  // Called from inside a parallel region (Minuit2's numerical gradient, the sampler's walkers), the pass runs on the
  // calling thread alone: the energy tasks are shared as ever (one solve per compound energy), only no second team
  // is opened.  Every record it writes is the evaluation's own (`results`, this working copy), so evaluations
  // running at once share nothing here.  Without OpenMP the same, serially.
#ifdef _OPENMP
  const bool ownTeam = !omp_in_parallel();
#else
  const bool ownTeam = false;
#endif
  (void)ownTeam;
  if (configure.environment.serialPointPass) return -1;  // AZURE_SERIAL_CHI2
  const std::vector<int> work = PassWork(segments);
  if (azure2::NodeCheckRequested()) {
    const long differ = CheckNodeTable(compound, work);
    if (differ) throw std::logic_error("AZURE_NODE_CHECK: " + std::to_string(differ) + " differ (listed above)");
  }
  AZUREProfile::Timer profileTimer(configure.profile, AZUREProfile::kPoints);
  // A plan task whose calculation threw: the reason (empty: none).  Its point gets NaN (Unit::MarkFailed) after the
  // pass, and its folds do not run.
  std::vector<std::string> failure(plan_.NumTasks());
  for (int t : work)
    if (plan_.Point(*this, t)->InjectedFailure()) failure[t] = "injected failure (AZURE_TEST_FAIL_E).";
  // The points with a fold, largest first (sub-points times points), so a large fold does not end up last on a thread.
  std::vector<int> folds;
  for (int t : work)
    if (azure2::HasSubPointFold(*this, *plan_.Point(*this, t))) folds.push_back(t);
  std::stable_sort(folds.begin(), folds.end(), [this](int a, int b) {
    auto cost = [this](int t) { return 1.0 + plan_.Point(*this, t)->NumSubPoints(); };
    return cost(a) > cost(b);
  });
  // Before the energy pass: the sub-points' energy-dependent values (external capture), point by point.
  if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
#pragma omp parallel for schedule(dynamic, 1) if (ownTeam)
    for (int k = 0; k < (int)folds.size(); k++) plan_.Point(*this, folds[k])->PrepareSubPoints(this, compound, configure);
  }
  // The energy pass: one task per compound energy, each writing only its units' records.  A unit whose calculation
  // throws fails its plan task; the units after it are calculated with a new solve.
  AZUREProfile::Timer setupTimer(configure.profile, AZUREProfile::kEnergyTaskSetup);
  const azure2::EnergyTaskList &tasks = plan_.EnergyTasks(*this, configure, work);
  setupTimer.Stop();
  AZUREProfile::CountMany(configure.profile, AZUREProfile::kEnergyTasks, (long long)tasks.size());
  std::vector<std::vector<std::pair<int, std::string>>> taskFailures(tasks.size());
  auto runEnergyTask = [&](int k) {
    const azure2::EnergyTask &task = tasks[k];
    const int n = task.size;
    AZUREProfile::CountMany(configure.profile, AZUREProfile::kEnergyUnits, n);
    std::vector<Site *> &units = results.scratch.ForThisThread().units;
    units.resize(n);
    for (int i = 0; i < n; i++) units[i] = &SiteAt(tasks.units[task.begin + i]);
    int start = 0;
    while (start < n) {
      int progress = 0;
      try {
        azure2::SolveAtOneEnergy(compound, state, thmCoherentValues_, configure, results, units.data() + start, n - start,
                                 tasks.tGroups.data() + task.begin + start, &progress);
        return;
      } catch (const std::exception &e) {
        taskFailures[k].push_back({tasks.owners[task.begin + start + progress], e.what()});
      } catch (...) {
        taskFailures[k].push_back({tasks.owners[task.begin + start + progress], "unknown exception."});
      }
      start += progress + 1;
    }
  };
  // The folds, after every sub-point is calculated: each writes only its points' records.
  auto runFold = [&](int k) {
    const int t = folds[k];
    if (!failure[t].empty()) return;
    try {
      plan_.Point(*this, t)->FoldSubPoints(*this, configure, results);
    } catch (const std::exception &e) {
      failure[t] = e.what();
    } catch (...) {
      failure[t] = "unknown exception.";
    }
  };
  auto collectFailures = [&]() {
    std::vector<char> solveFailed(plan_.NumTasks(), 0);
    bool any = false;
    for (const auto &list : taskFailures)
      for (const auto &f : list) {
        if (failure[f.first].empty()) failure[f.first] = f.second;
        solveFailed[f.first] = 1, any = true;
      }
    // A slot whose solve failed fails every point that reads it, not only the one that writes it (SlotGroup.h).
    if (any)
      for (int t : folds) {
        const Unit &point = *plan_.Point(*this, t);
        for (int i = 1; i <= point.NumSubPoints() && failure[t].empty(); i++) {
          const int writer = UnitAt(SlotWriter(point.SubPointFirst() + i - 1)).PlanTask();
          if (solveFailed[writer]) failure[t] = failure[writer];
        }
      }
  };
  if (azure2::CalculationPlan::RunTimeCheckRequested()) {
    // One step after the other, each checked to write only its own records (the energy tasks: their units'; a fold:
    // its points'); then every unit's shared solve against its own.  The row sets' values after the last task.
    const int numTasks = (int)tasks.size();
    if (numTasks == 0) ContractRowSets(results);
    azure2::CalculationPlan::RunChecked(
        results, numTasks + (int)folds.size(),
        [&](int k, std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
          if (k < numTasks) {
            runEnergyTask(k);
            for (int i = 0; i < tasks[k].size; i++)
              azure2::CalculationPlan::AddRecordRange(SiteAt(tasks.units[tasks[k].begin + i]), ranges);
            if (k == numTasks - 1) collectFailures(), ContractRowSets(results);
          } else {
            // A fold writes its point's record and reads its slots'.
            runFold(k - numTasks);
            azure2::CalculationPlan::AddRecordRange(*plan_.Point(*this, folds[k - numTasks]), ranges);
          }
        },
        [&](int k) {
          return k < numTasks ? "energy task " + std::to_string(k) : "the fold of plan task " + std::to_string(folds[k - numTasks]);
        });
    if (numTasks == 0) collectFailures();
    const long differ = azure2::CheckEnergyTasks(*this, compound, state, configure, results, tasks);
    if (differ) throw std::logic_error("AZURE_PLAN_CHECK: " + std::to_string(differ) +
                                       " units differ from their own solve (listed above)");
  } else {
#pragma omp parallel for schedule(dynamic, 1) if (ownTeam)
    for (int k = 0; k < (int)tasks.size(); k++) runEnergyTask(k);
    collectFailures();
    ContractRowSets(results);
#pragma omp parallel for schedule(dynamic, 1) if (ownTeam)
    for (int k = 0; k < (int)folds.size(); k++) runFold(k);
  }
  // NaN on the point of a failed plan task, reported once (Unit::MarkFailed); the CLI's objective then recalculates
  // point by point (PrecalculatePoints).
  int failed = 0;
  for (int t : work)
    if (!failure[t].empty()) {
      plan_.Point(*this, t)->MarkFailed(results, "point", configure, failure[t].c_str());
      failed++;
    }
  return failed;
}
