#ifndef STATEPOOL_H
#define STATEPOOL_H
// The idle values of the evaluations that run at the same time (a fit's threads, a sampler's walkers): a LevelState,
// a walker's state and inputs, or a working copy with its buffer (WorkingData.h).  Owned by whoever runs the fitter
// or the sampler, as its azure2::Results are, and lent to it, so that the evaluations (const methods) take and give
// back values that are not the evaluator's own.
#include <memory>
#include <mutex>
#include <stack>
#include <utility>

template <class Values>
class StatePool {
 public:
  /// An idle value, or a copy of `seed` when none is idle (the caller's values at this moment, as a copy of the
  /// compound copied them).
  std::unique_ptr<Values> Take(const Values &seed) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!idle_.empty()) {
      std::unique_ptr<Values> values = std::move(idle_.top());
      idle_.pop();
      return values;
    }
    return std::make_unique<Values>(seed);
  }
  /// An idle value, or what `make` returns when none is idle (a value that is made, not copied: a working copy).
  /// `make` runs under the pool's lock, so two threads never make one at the same time.
  template <class Make>
  std::unique_ptr<Values> TakeOrMake(Make make) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!idle_.empty()) {
      std::unique_ptr<Values> values = std::move(idle_.top());
      idle_.pop();
      return values;
    }
    return make();
  }
  /// Back to the idle ones; the next Take returns it with what the last evaluation left in it.
  void Give(std::unique_ptr<Values> values) {
    if (!values) return;
    std::lock_guard<std::mutex> lock(mutex_);
    idle_.push(std::move(values));
  }

 private:
  std::stack<std::unique_ptr<Values>> idle_;
  std::mutex mutex_;  ///< Take and Give run on every evaluating thread
};

#endif
