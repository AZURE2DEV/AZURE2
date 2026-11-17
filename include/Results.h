#ifndef AZURE2_RESULTS_H
#define AZURE2_RESULTS_H

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

#include "TaskScratch.h"
#include "azure2/inference/LinearProfile.h"
#include "azure2/support/Environment.h"

namespace azure2 {

/*!
 * What one evaluation of the model produces (docs/upgrade/data_model.md section 2.1, decision M2 (c)): the values
 * every calculation unit (a data or test point, and each of its sub-points) computes for the current level
 * parameters and nuisances.  A unit owns one record, a fixed range of `slot` at an offset assigned when the units
 * are laid out (WorkingCopy::LayoutResults); nothing is looked up or allocated while the points are calculated.
 *
 * A record holds, in this order:
 *   kSigma    the unit's observable: sigma (b), dsigma/dOmega_c.m. (b/sr), a phase shift (deg), A_y (1), a
 *             Legendre coefficient a_k (1) or A_k (b/sr), the THM cross section; after the target fold for a
 *             point that carries one (never the segment's combination: that is the model value, below)
 *   kSigmaE1  the E1 part of an angle-integrated capture cross section (b)
 *   kSigmaE2  the E2 part (b)
 *   kAy       the analysing power A_y (1), kept beside the cross section, which weights its target average
 *   kA0       A_0 = sigma_int / 4 pi (b/sr) of an angular-distribution unit
 *   then, for a unit with an angular-distribution capacity K > 0: the number of coefficients written, then K values
 *   a_0 .. a_{K-1}; then, for the point of a sample (natural) segment with N > 0 nuclei: the number of terms written,
 *   then N abundance-weighted terms, one per compound nucleus (written with the model value, EvaluateModelPoint).
 *
 * Who holds what.  The layout -- where each unit's record starts and how many angular and natural slots it has -- is
 * a property of the data set and of every copy of it (WorkingCopy::LayoutResults; a unit keeps its offset).  The storage
 * belongs to whoever starts the evaluation: the CLI's calculation, each pooled working copy of a fit or of the MCMC,
 * the Python session, a probe or a reaction-rate point.  A unit holds no pointer to a buffer, so a copy of a unit or
 * of a data set cannot write into its original's values.
 *
 * The one way to reach a result: every accessor and every function that reads or writes one takes the buffer as an
 * argument, `const Results &` to read and `Results &` to write -- `point->GetFitCrossSection(results)`,
 * `point->Calculate(compound, state, config, results)`, `azure2::EvaluateModelPoint(segment, i, compound, state, config, results,
 * leaves)`.  There is no default buffer, no global and no thread-local one.  The owner hands its buffer to
 * WorkingCopy::FillEnergyShiftsFromParams, the last step of the hand-over (WorkingCopy::BeginEvaluation), which sizes it for the
 * data set's layout and, under AZURE_POISON_RESULTS=1, fills it with NaN.
 *
 * The struct is plain (no engine headers; inference/'s LinearProfile.h is plain too), so it can move to the experiment
 * library.
 */
struct Results {
  enum Quantity : std::size_t { kSigma = 0, kSigmaE1 = 1, kSigmaE2 = 2, kAy = 3, kA0 = 4, kNumFixed = 5 };

  std::vector<double> slot;

  /// The values the readers of angle-free slots fold: per row set (RowSet.h) and slot, its row contracted with the
  /// slot's record (WorkingCopy::ContractRowSets, Site::SlotValue); laid out by WorkingCopy::BuildRowSets.
  std::vector<double> contracted;

  /*!
   * What the data set predicts, apart from the leaf values above (docs/upgrade/data_model.md section 2.1, step S4):
   * one record per point of a data or test segment -- not of a leaf or a sub-point --, laid out with the
   * leaf records (WorkingCopy::LayoutResults; the point keeps its offset, Unit::SetModelRecord).  A record holds the model
   * value m_i, the segment's combination (SUM, RATIO, sample, transform, the E1 or E2 part) of its points' leaf
   * values, which the residual compares with the datum; then, for a point with an angular-distribution capacity
   * K > 0, the number of coefficients and K coefficients a_0 .. a_{K-1}.  Written by azure2::EvaluateModelPoint
   * (ModelPass.h) only.
   */
  std::vector<double> model;

  /*!
   * What the objective (Objective.h) writes, index = segment number (1-based; [0] unused): the normalization the
   * segment's data were compared with (the parameter's, or the profiled n* of a THM segment), its data chi^2 (a THM
   * experiment's segments: each one's share of the profile), the profile of each THM experiment (WorkingCopy::thmGroups_
   * order), and the totals of the terms that are not per segment, by slot.  The writers and the API read them here.
   */
  struct ObjectiveValues {
    std::vector<double> norm, segmentChi2;
    /// The covariance seed (M^-1 r)_i of every point in a data-covariance group (WorkingCopy::CovarianceChi2Correction):
    /// the point's share of d(chi^2)/d(r_i) / 2, which the gradient reads.  Indexed as Results::model (the point's
    /// model record, Unit::SetModelRecord); sized only when the data set has a covariance group.
    std::vector<double> covarianceSeed;
    std::vector<inference::ThmProfile> thmProfiles;
    /*!
     * The total of each term that keeps one, index = the term's slot (its place in Objective::Terms(),
     * ObjectiveTerm::Slot), and whether it was stored: written by the term's Add when its condition holds, read by
     * its field of chiSquared.out (ObjectiveTerm::WriteTotal).  Grown on the first store and never cleared: a total
     * keeps the value of the last evaluation that stored it.  0 and false for a slot never stored.
     */
    std::vector<double> totals;
    std::vector<char> stored;
    void StoreTotal(int slot, double value) {
      if (slot < 0) return;
      if ((std::size_t)slot >= totals.size()) {
        totals.resize((std::size_t)slot + 1, 0.0);
        stored.resize((std::size_t)slot + 1, 0);
      }
      totals[(std::size_t)slot] = value;
      stored[(std::size_t)slot] = 1;
    }
    double Total(int slot) const { return slot >= 0 && (std::size_t)slot < totals.size() ? totals[(std::size_t)slot] : 0.0; }
    bool Stored(int slot) const { return slot >= 0 && (std::size_t)slot < stored.size() && stored[(std::size_t)slot]; }
  } objective;

  /// The scratch of the energy tasks that write this buffer, one per thread (TaskScratch.h): solvers and their
  /// buffers, never a value of the evaluation.  A copy of the buffer starts with empty scratch.
  TaskScratchSet scratch;

  /// Length of a model record with angular capacity nAng.
  static std::size_t ModelRecordSize(std::size_t nAng) { return 1 + (nAng ? 1 + nAng : 0); }

  /// Length of a record with angular capacity nAng and nNat natural terms.
  static std::size_t RecordSize(std::size_t nAng, std::size_t nNat) {
    return kNumFixed + (nAng ? 1 + nAng : 0) + (nNat ? 1 + nNat : 0);
  }
  /// Offset, inside a record, of the angular-distribution count (the values follow it).
  static std::size_t AngularOffset() { return kNumFixed; }
  /// Offset, inside a record, of the natural-term count (the terms follow it).
  static std::size_t NaturalOffset(std::size_t nAng) { return kNumFixed + (nAng ? 1 + nAng : 0); }

  /// Every slot, model value and per-segment objective value quiet NaN: a value read before the evaluation writes it
  /// then shows as NaN.
  void Poison() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::vector<double> *v : {&slot, &contracted, &model, &objective.norm, &objective.segmentChi2, &objective.covarianceSeed})
      std::fill(v->begin(), v->end(), nan);
  }

  /// AZURE_POISON_RESULTS=1: reused result storage is poisoned on entry to every evaluation.
  static bool PoisonRequested() {
    static const bool on = azure2::support::ProcessEnvironment().On(azure2::support::Switch::kPoisonResults);
    return on;
  }
};

/// Copy a record laid out with (fromAng, fromNat) slots into one with (toAng, toNat): the fixed quantities, then
/// as much of the angular and natural parts as both hold (the counts included).
inline void CopyRecord(const double *from, std::size_t fromAng, std::size_t fromNat, double *to, std::size_t toAng,
                       std::size_t toNat) {
  std::copy(from, from + Results::kNumFixed, to);
  if (fromAng && toAng) {
    const std::size_t n = std::min(fromAng, toAng);
    std::copy(from + Results::AngularOffset(), from + Results::AngularOffset() + 1 + n, to + Results::AngularOffset());
    if (to[Results::AngularOffset()] > (double)toAng) to[Results::AngularOffset()] = (double)toAng;
  }
  if (fromNat && toNat) {
    const std::size_t n = std::min(fromNat, toNat);
    const double *f = from + Results::NaturalOffset(fromAng);
    double *t = to + Results::NaturalOffset(toAng);
    std::copy(f, f + 1 + n, t);
    if (t[0] > (double)toNat) t[0] = (double)toNat;
  }
}

}  // namespace azure2

#endif
