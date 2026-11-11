#ifndef AZURE2_EVALUATIONCOUNT_H
#define AZURE2_EVALUATIONCOUNT_H

#include <atomic>

/// Objective evaluations, incremented from every thread of a parallel gradient; a copy starts at its source's count.
struct EvaluationCount {
  std::atomic<int> n{0};
  EvaluationCount() = default;
  EvaluationCount(const EvaluationCount &o) : n(o.n.load(std::memory_order_relaxed)) {}
  EvaluationCount &operator=(const EvaluationCount &o) {
    n.store(o.n.load(std::memory_order_relaxed), std::memory_order_relaxed);
    return *this;
  }
};

#endif
