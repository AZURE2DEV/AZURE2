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
#include "CNuc.h"
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

// Whether unit k starts a T-matrix group: by `tGroups`, or, if null, by a change of (entrance, exit) pair.
bool StartsTMatrix(Site *const *units, int k, const int *tGroups) {
  if (k == 0) return true;
  if (tGroups) return tGroups[k] != tGroups[k - 1];
  return units[k]->GetEntranceKey() != units[k - 1]->GetEntranceKey() ||
         units[k]->GetExitKey() != units[k - 1]->GetExitKey();
}

}  // namespace

void SolveAtOneEnergy(CNuc *cnuc, const Config &configure, Results &results, Site *const *units, int n,
                      const int *tGroups, int *progress) {
  if (progress) *progress = 0;
  if (n <= 0) return;
  if (units[0]->IsTHM()) {
    // THM (modified R-matrix) units: the shared interior (the level matrix) and the half-off-shell cross section of
    // each.  The object is reused across the energies of this thread: its per-J-group buffers are sized once (Reset
    // repoints it at this compound and configuration first, so nothing of an earlier session is read).
    static thread_local THMMatrixFunc thm(cnuc, configure);
    thm.Reset(cnuc, configure);
    thm.ClearMatrices();
    thm.FillMatrices(units[0]);
    thm.InvertMatrices();
    for (int k = 0; k < n; k++) {
      if (progress) *progress = k;
      thm.CalculateTHMCrossSection(units[k], results);
    }
    return;
  }
  // The A-matrix object is reused in the same way; the R-matrix route builds its own.
  static thread_local AMatrixFunc amatrix(cnuc, configure);
  std::unique_ptr<RMatrixFunc> rmatrix;
  GenMatrixFunc *func;
  if (configure.paramMask & Config::USE_AMATRIX) {
    amatrix.Reset(cnuc, configure);
    func = &amatrix;
  } else {
    rmatrix.reset(new RMatrixFunc(cnuc, configure));
    func = rmatrix.get();
  }
  func->ClearMatrices();
  func->FillMatrices(units[0]);
  func->InvertMatrices();
  // The first T-matrix after a fill is built directly, every later one through RecalculateTMatrix (which drops the
  // elements held).
  for (int k = 0; k < n; k++) {
    if (progress) *progress = k;
    if (k == 0)
      func->CalculateTMatrix(units[k]);
    else if (StartsTMatrix(units, k, tGroups))
      func->RecalculateTMatrix(units[k]);
    func->CalculateCrossSection(units[k], results);
  }
}

namespace {

// What differs between the shared solve `a` and a unit's own solve `b` at one energy: the route and the bilinear forms
// g_c^T A g_c' of every J-group of the nucleus (the only way the level matrix enters a T-matrix or a THM amplitude).
std::string CompareLevelMatrices(CNuc *cnuc, int nucleus, AMatrixFunc &a, AMatrixFunc &b) {
  for (int j = 1; j <= cnuc->NumJGroups(); j++) {
    JGroup *jg = cnuc->GetJGroup(j);
    if (!jg->IsInRMatrix() || jg->GetNucleus() != nucleus) continue;
    if (a.IsChannelRoute(j) != b.IsChannelRoute(j)) return "route of J-group " + std::to_string(j);
    for (int c = 1; c <= jg->NumChannels(); c++)
      for (int cp = 1; cp <= jg->NumChannels(); cp++) {
        const complex x = a.GetUBilinear(j, c, cp), y = b.GetUBilinear(j, c, cp);
        if (std::memcmp(&x, &y, sizeof x) != 0)
          return "level matrix of J-group " + std::to_string(j) + " (g^T A g at channels " + std::to_string(c) + ", " +
                 std::to_string(cp) + ")";
      }
  }
  return std::string();
}

std::string Describe(const Site &u, const EnergyTask &task) {
  std::ostringstream s;
  s.precision(17);
  s << "segment key " << u.GetSegmentKey() << ", entrance " << u.GetEntranceKey() << ", exit " << u.GetExitKey()
    << ", E_cm " << u.GetCMEnergy() << " MeV (E_x " << task.energy << " MeV, nucleus " << task.nucleus << ", route "
    << task.route << ")";
  return s.str();
}

}  // namespace

long CheckEnergyTasks(WorkingCopy &data, CNuc *cnuc, const Config &configure, const Results &results,
                      const EnergyTaskList &tasks) {
  Results scratch = results;  // the check writes here, never into the evaluation's values
  std::unique_ptr<AMatrixFunc> sharedA, ownA;
  if (configure.paramMask & Config::USE_AMATRIX) {
    sharedA.reset(new AMatrixFunc(cnuc, configure));
    ownA.reset(new AMatrixFunc(cnuc, configure));
  }
  std::unique_ptr<THMMatrixFunc> sharedT(new THMMatrixFunc(cnuc, configure)), ownT(new THMMatrixFunc(cnuc, configure));
  long units = 0, shared = 0, differ = 0;
  std::vector<std::string> reports;
  std::vector<double> own;
  for (const EnergyTask &task : tasks.tasks) {
    const int n = task.size;
    const UnitRef *taskUnits = tasks.units.data() + task.begin;
    const int *taskGroups = tasks.tGroups.data() + task.begin;
    units += n;
    if (n > 1) shared += n;
    const bool thm = task.route == 'T';
    std::unique_ptr<RMatrixFunc> sharedR, ownR;
    GenMatrixFunc *s, *o;
    if (thm)
      s = sharedT.get(), o = ownT.get();
    else if (task.route == 'A')
      s = sharedA.get(), o = ownA.get();
    else {
      sharedR.reset(new RMatrixFunc(cnuc, configure)), ownR.reset(new RMatrixFunc(cnuc, configure));
      s = sharedR.get(), o = ownR.get();
    }
    Site *first = &data.SiteAt(taskUnits[0]);
    s->ClearMatrices();
    s->FillMatrices(first);
    s->InvertMatrices();
    for (int k = 0; k < n; k++) {
      Site *u = &data.SiteAt(taskUnits[k]);
      std::string why;
      if (AMatrixFunc::InEnergy(cnuc, u) != task.energy) why = "E_x is not the task's";
      // The unit's own solve.
      o->ClearMatrices();
      o->FillMatrices(u);
      o->InvertMatrices();
      if (why.empty() && task.route != 'R')
        why = CompareLevelMatrices(cnuc, task.nucleus, *static_cast<AMatrixFunc *>(s), *static_cast<AMatrixFunc *>(o));
      if (!thm) {
        o->CalculateTMatrix(u);
        if (k == 0)
          s->CalculateTMatrix(u);
        else if (taskGroups[k] != taskGroups[k - 1])
          s->RecalculateTMatrix(u);
        if (why.empty() && !s->SameTMatrix(*o)) why = "T-matrix";
      }
      // The unit's record from its own solve, then from the shared one.
      const std::size_t slot = u->ResultSlot();
      const std::size_t size = Results::RecordSize(u->ResultAngularSlots(), 0);
      if (thm)
        static_cast<THMMatrixFunc *>(o)->CalculateTHMCrossSection(u, scratch);
      else
        o->CalculateCrossSection(u, scratch);
      if (slot != Unit::kNoRecord) own.assign(scratch.slot.begin() + slot, scratch.slot.begin() + slot + size);
      if (thm)
        static_cast<THMMatrixFunc *>(s)->CalculateTHMCrossSection(u, scratch);
      else
        s->CalculateCrossSection(u, scratch);
      if (why.empty() && slot != Unit::kNoRecord &&
          std::memcmp(own.data(), scratch.slot.data() + slot, size * sizeof(double)) != 0)
        why = "record";
      if (why.empty()) continue;
      differ++;
      if (reports.size() < 10) reports.push_back(why + " differs at unit " + std::to_string(taskUnits[k]) + " (" +
                                                 Describe(*u, task) + ") from unit " + std::to_string(taskUnits[0]));
    }
  }
  std::cerr << "AZURE_PLAN_CHECK: " << units << " units in " << tasks.size() << " energy tasks (" << shared
            << " units at a shared energy) against their own solves: " << differ << " differ" << std::endl;
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
 * One parallel pass over the points of an evaluation, shared by the CLI
 * (AZURECalc::PrecalculatePoints) and the API's calculations (AZUREAPI).
 */

int WorkingCopy::CalculatePointsInParallel(const std::vector<const Segment *> &segments, CNuc *theCNuc, const Config &configure,
                                     azure2::Results &results) {
#ifdef _OPENMP
  if (omp_in_parallel()) return -1;
#else
  return -1;
#endif
  if (std::getenv("AZURE_SERIAL_CHI2")) return -1;
  const std::vector<int> work = PassWork(segments);
  if (azure2::NodeCheckRequested()) {
    const long differ = CheckNodeTable(theCNuc, work);
    if (differ) throw std::logic_error("AZURE_NODE_CHECK: " + std::to_string(differ) + " differ (listed above)");
  }
  AZUREProfile::Timer profileTimer(AZUREProfile::kPoints);
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
#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < (int)folds.size(); k++) plan_.Point(*this, folds[k])->PrepareSubPoints(this, theCNuc, configure);
  }
  // The energy pass: one task per compound energy, each writing only its units' records.  A unit whose calculation
  // throws fails its plan task; the units after it are calculated with a new solve.
  AZUREProfile::Timer setupTimer(AZUREProfile::kEnergyTaskSetup);
  const azure2::EnergyTaskList &tasks = plan_.EnergyTasks(*this, configure, work);
  setupTimer.Stop();
  AZUREProfile::CountMany(AZUREProfile::kEnergyTasks, (long long)tasks.size());
  std::vector<std::vector<std::pair<int, std::string>>> taskFailures(tasks.size());
  auto runEnergyTask = [&](int k) {
    const azure2::EnergyTask &task = tasks[k];
    const int n = task.size;
    AZUREProfile::CountMany(AZUREProfile::kEnergyUnits, n);
    static thread_local std::vector<Site *> units;
    units.resize(n);
    for (int i = 0; i < n; i++) units[i] = &SiteAt(tasks.units[task.begin + i]);
    int start = 0;
    while (start < n) {
      int progress = 0;
      try {
        azure2::SolveAtOneEnergy(theCNuc, configure, results, units.data() + start, n - start,
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
      plan_.Point(*this, t)->FoldSubPoints(this, theCNuc, configure, results);
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
    // its points'); then every unit's shared solve against its own.
    const int numTasks = (int)tasks.size();
    azure2::CalculationPlan::RunChecked(
        results, numTasks + (int)folds.size(),
        [&](int k, std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
          if (k < numTasks) {
            runEnergyTask(k);
            for (int i = 0; i < tasks[k].size; i++)
              azure2::CalculationPlan::AddRecordRange(SiteAt(tasks.units[tasks[k].begin + i]), ranges);
            if (k == numTasks - 1) collectFailures();
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
    const long differ = azure2::CheckEnergyTasks(*this, theCNuc, configure, results, tasks);
    if (differ) throw std::logic_error("AZURE_PLAN_CHECK: " + std::to_string(differ) +
                                       " units differ from their own solve (listed above)");
  } else {
#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < (int)tasks.size(); k++) runEnergyTask(k);
    collectFailures();
#pragma omp parallel for schedule(dynamic, 1)
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
