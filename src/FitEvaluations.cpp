// The evaluations of the objective while a fitter runs (FitEvaluations.h).
#include "FitEvaluations.h"

#include <cmath>
#include <string>

#include "AZURECalc.h"
#include "WorkingCopy.h"
#include "CompoundSystem.h"

void HandOverEvaluation(const FitControl &control, const FitProgress &record) {
  if (control.progress) control.progress(record);
  if (control.checkpoint && record.kind == FitProgress::Kind::kEvaluation && !record.errorAnalysis &&
      record.evaluation != 0 && record.evaluation % AZURECalc::kOutputInterval == 0 && record.x)
    control.checkpoint(*record.x);
}

FitEvaluations::FitEvaluations(AZURECalc &objective, const FitControl &control, const char *fitter,
                               const char *stage) :
  objective_(objective), control_(control), fitter_(fitter), stage_(stage), previous_(objective.Observer()) {
  observer_.evaluated = [this](long evaluation, double value, const std::vector<double> &p, bool errorAnalysis) {
    Evaluated(evaluation, value, p, errorAnalysis);
  };
  observer_.stop = control.stop;
  objective_.Observe(&observer_);
}

FitEvaluations::~FitEvaluations() { objective_.Observe(previous_); }

void FitEvaluations::Evaluated(long evaluation, double value, const std::vector<double> &p, bool errorAnalysis) {
  FitProgress record;
  record.kind = FitProgress::Kind::kEvaluation;
  record.fitter = fitter_;
  record.stage = stage_;
  record.evaluation = evaluation;
  record.objective = value;
  // The working copy's flag is the evaluation's own (a pooled copy made before MINOS began does not carry it): the
  // stage says it as well.
  record.errorAnalysis = errorAnalysis || std::string(stage_) == "minos";
  record.x = &p;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    record.calls = ++count_;
    if (std::isfinite(value) && (!hasBest_ || value < best_)) {
      hasBest_ = true;
      best_ = value;
      bestPoint_ = p;
    }
    if (hasBest_) record.best = best_;
  }
  HandOverEvaluation(control_, record);
}

long FitEvaluations::Count() {
  std::lock_guard<std::mutex> lock(mutex_);
  return count_;
}

bool FitEvaluations::HasBest() {
  std::lock_guard<std::mutex> lock(mutex_);
  return hasBest_;
}

double FitEvaluations::BestObjective() {
  std::lock_guard<std::mutex> lock(mutex_);
  return best_;
}

std::vector<double> FitEvaluations::BestPoint() {
  std::lock_guard<std::mutex> lock(mutex_);
  return bestPoint_;
}

bool FitEvaluations::Stopped() const { return control_.stop && control_.stop->load(std::memory_order_relaxed); }
