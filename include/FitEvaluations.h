#ifndef FITEVALUATIONS_H
#define FITEVALUATIONS_H

// The evaluations of the objective while a fitter runs (docs/upgrade/fit.md section 3): counted, the best point kept,
// each handed to the caller's FitControl as a progress record.

#include <mutex>
#include <vector>

#include "FitControl.h"

class AZURECalc;

/// Hands one record to `control`: its progress sink, then its checkpoint every AZURECalc::kOutputInterval evaluations
/// of the run (an evaluation record outside MINOS).  The one rule of the program's intermediate files for evaluations.
void HandOverEvaluation(const FitControl &control, const FitProgress &record);

/// While it lives, the observer of `objective`: every evaluation is counted, the lowest objective and its point are
/// kept, and the record goes to HandOverEvaluation with control's stop request in force.  The objective's previous
/// observer is put back at the end.
class FitEvaluations {
 public:
  /// `stage`: "fit", or "hesse" / "minos" for the error analysis after it (FitProgress::stage).
  FitEvaluations(AZURECalc &objective, const FitControl &control, const char *fitter, const char *stage = "fit");
  ~FitEvaluations();
  FitEvaluations(const FitEvaluations &) = delete;
  FitEvaluations &operator=(const FitEvaluations &) = delete;

  /// Evaluations seen so far.
  long Count();
  /// Whether a finite objective was seen, the lowest and its full parameter vector.
  bool HasBest();
  double BestObjective();
  std::vector<double> BestPoint();
  /// The caller's stop request is set.
  bool Stopped() const;

 private:
  void Evaluated(long evaluation, double objective, const std::vector<double> &p, bool errorAnalysis);

  AZURECalc &objective_;
  const FitControl &control_;
  const char *fitter_;
  const char *stage_;
  const EvaluationObserver *previous_;
  EvaluationObserver observer_;
  std::mutex mutex_;  ///< the counters and the best point: evaluations arrive from several threads
  long count_ = 0;
  bool hasBest_ = false;
  double best_ = 0.0;
  std::vector<double> bestPoint_;
};

#endif
