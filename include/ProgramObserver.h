#ifndef PROGRAMOBSERVER_H
#define PROGRAMOBSERVER_H

// What a front end that runs the program in its own process hands to the run (docs/upgrade/qt_gui.md section 3): the
// stop request it may set from another thread, and the sinks of the fit's progress and of the final chi^2.  Owned by
// the caller, outliving the run; the run holds a pointer (Config::observer, null for the command line, which prints
// these on standard output only).  Every member may be empty.  The sinks are called from the run's threads (the
// progress from several at once): each must be thread-safe and must not call back into the run.

#include <atomic>
#include <functional>

#include "FitControl.h"

struct ProgramObserver {
  /// Set by the caller to stop the run: a fit ends at its best point and the run writes its outputs there; the
  /// set-up's channel rows and external-capture integrals end early and the run is aborted.  A calculation without a
  /// fit runs to its end once its set-up is done.
  const std::atomic<bool> *stop = nullptr;
  /// Every progress record of the run's objective and fitter (the records the console's "Iteration:" lines print).
  std::function<void(const FitProgress &)> progress;
  /// The final evaluation's chi^2 of one data set: its segment number (the "Segment #" of the console and the first
  /// column of chiSquared.out), its chi^2 and its number of points.
  std::function<void(int segment, double chiSquared, int points)> dataSetChiSquared;
  /// The final evaluation's total chi^2 (data and priors).
  std::function<void(double chiSquared)> totalChiSquared;
  /// After each step of the sampler (command mcmc): the step (from 1), the steps asked, the walkers' mean
  /// log-posterior, log-likelihood and log-prior after it (the numbers of the console's progress line).
  std::function<void(int step, int total, double logPosterior, double logLikelihood, double logPrior)> samplerStep;
};

#endif
