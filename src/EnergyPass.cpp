// The point pass of an evaluation (WorkingCopy::CalculatePointsInParallel); the calculation at one compound energy (CalculationPlan.h, EnergyTask): one fill and factorisation of the level
// matrix, one T-matrix per (entrance, exit) pair, the observable of every unit; and the check (AZURE_PLAN_CHECK=1) that
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

// The T-matrix of `u`'s pair held by `func` (`held`: the pair of the T-matrix held now; {0, 0} none): built anew when
// the pair changes, else kept.  The first T-matrix after a fill is built directly, every later one through
// RecalculateTMatrix (which drops the elements held), as a host and its mapped points always did.
void TMatrixFor(GenMatrixFunc &func, Unit *u, std::pair<int, int> &held) {
  const std::pair<int, int> pair{u->GetEntranceKey(), u->GetExitKey()};
  if (pair == held) return;
  if (held == std::pair<int, int>{0, 0})
    func.CalculateTMatrix(u);
  else
    func.RecalculateTMatrix(u);
  held = pair;
}

}  // namespace

void SolveAtOneEnergy(CNuc *cnuc, const Config &configure, Results &results, Unit *const *units, int n,
                      int *progress) {
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
  std::pair<int, int> held{0, 0};
  for (int k = 0; k < n; k++) {
    if (progress) *progress = k;
    TMatrixFor(*func, units[k], held);
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

std::string Describe(const Unit &u, const EnergyTask &task) {
  std::ostringstream s;
  s.precision(17);
  s << "segment key " << u.GetSegmentKey() << ", entrance " << u.GetEntranceKey() << ", exit " << u.GetExitKey()
    << ", E_cm " << u.GetCMEnergy() << " MeV (E_x " << task.energy << " MeV, nucleus " << task.nucleus << ", route "
    << task.route << ")";
  return s.str();
}

}  // namespace

long CheckEnergyTasks(WorkingCopy &data, CNuc *cnuc, const Config &configure, const Results &results,
                      const std::vector<EnergyTask> &tasks) {
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
  for (const EnergyTask &task : tasks) {
    const int n = (int)task.units.size();
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
    Unit *first = &data.UnitAt(task.units[0]);
    s->ClearMatrices();
    s->FillMatrices(first);
    s->InvertMatrices();
    std::pair<int, int> held{0, 0};
    for (int k = 0; k < n; k++) {
      Unit *u = &data.UnitAt(task.units[k]);
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
        TMatrixFor(*s, u, held);
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
      if (reports.size() < 10) reports.push_back(why + " differs at unit " + std::to_string(task.units[k]) + " (" +
                                                 Describe(*u, task) + ") from unit " + std::to_string(task.units[0]));
    }
  }
  std::cerr << "AZURE_PLAN_CHECK: " << units << " units in " << tasks.size() << " energy tasks (" << shared
            << " units at a shared energy) against their own solves: " << differ << " differ" << std::endl;
  for (const std::string &r : reports) std::cerr << "AZURE_PLAN_CHECK:   " << r << std::endl;
  return differ;
}

}  // namespace azure2

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
  if (ScopedPointMappingOff::Active()) return -1;
  if (std::getenv("AZURE_SERIAL_CHI2")) return -1;
  std::vector<Unit *> work;
  std::set<const Segment *> seen;
  for (const Segment *seg : segments) {
    if (!seg || !seen.insert(seg).second) continue;
    for (int t = 0; t < seg->NumTerms(); t++)
      for (Unit &pt : seg->TermPoints(*this, t))
        if (!pt.IsMapped()) work.push_back(&pt);
  }
  // A host also calculates every point mapped onto it (WorkingCopy::MapData), so
  // its work grows with that number: one Fowler-type angular distribution at
  // a single energy is one item.  Largest items first, one at a time, so such
  // an item does not end up last on one thread.  Each item is calculated the
  // same way whatever the order.
  std::vector<double> cost(work.size());
  for (size_t k = 0; k < work.size(); k++)
    cost[k] = (1.0 + work[k]->NumSubPoints()) * (1.0 + work[k]->PlanMapped(this).size());
  std::vector<size_t> order(work.size());
  for (size_t k = 0; k < order.size(); k++) order[k] = k;
  std::stable_sort(order.begin(), order.end(), [&cost](size_t a, size_t b) { return cost[a] > cost[b]; });
  std::vector<Unit *> sorted(work.size());
  for (size_t k = 0; k < order.size(); k++) sorted[k] = work[order[k]];
  work.swap(sorted);
  int failed = 0;
  AZUREProfile::Timer profileTimer(AZUREProfile::kPoints);
  // One item per task of the plan: its host calculates itself and the points the plan maps onto it.
  auto runTask = [&](Unit *host) {
    try {
      host->Calculate(this, theCNuc, configure, results);
    } catch (const std::exception &e) {
      // NaN on the point and the points mapped onto it, reported once (Unit::MarkFailed); 0 on the point alone,
      // silently, before 2026-11-01 -- and the CLI's objective takes these values (PrecalculatePoints).
      host->MarkFailed(this, results, "point", configure, e.what());
      return 1;
    } catch (...) {
      host->MarkFailed(this, results, "point", configure, "unknown exception.");
      return 1;
    }
    return 0;
  };
  if (azure2::CalculationPlan::RunTimeCheckRequested() && !work.empty()) {
    azure2::CalculationPlan::RunChecked(
        results, (int)work.size(),
        [&](int k, std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
          failed += runTask(work[k]);
          plan_.WriteRanges(*this, work[k]->PlanTask(), ranges);
        },
        [&](int k) { return "task " + std::to_string(work[k]->PlanTask()); });
    std::vector<int> tasks;
    for (Unit *h : work) tasks.push_back(h->PlanTask());
    std::sort(tasks.begin(), tasks.end());
    const long differ =
        azure2::CheckEnergyTasks(*this, theCNuc, configure, results, plan_.EnergyTasks(*this, theCNuc, configure, tasks));
    if (differ) throw std::logic_error("AZURE_PLAN_CHECK: " + std::to_string(differ) +
                                       " units differ from their own solve (listed above)");
    return failed;
  }
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : failed)
  for (int k = 0; k < (int)work.size(); k++) failed += runTask(work[k]);
  return failed;
}
