#ifndef FITERRORS_H
#define FITERRORS_H

// The error analysis on a FitResult of a MIGRAD row (docs/upgrade/fit.md section 3): HESSE refines the covariance,
// MINOS gives the asymmetric intervals.  The program (src/ProgramFit.cpp) and a pyazr session
// (api/src/AZUREAPIFit.cpp) call the same functions.

#include <ostream>

#include "FitControl.h"

class AZURECalc;
class AZUREParams;

/// Whether `result` has a free parameter for HESSE: Minuit2 crashes in MnHesse on a state with none.
bool HasFreeParameterForHesse(const FitResult &result);

/// HESSE on result.minimum (updated in place).  True when it gives a valid covariance, which then replaces
/// result.covariance (covarianceFromHesse).  Needs result.minimum and a free parameter (HasFreeParameterForHesse).
bool RunHesse(AZURECalc &objective, FitResult &result);

/*!
 * MINOS on result.minimum with Up = `up` (the objective's error definition is set to it): for every parameter of
 * `params` in order, "\tParameter <n>..." on `log`, then its (lower, upper) interval, (0, 0) for a fixed one, into
 * result.minos, and the validity of each side (Minuit2's MinosError) into result.minosValidity.  The working copy is
 * put in error-analysis mode; the caller ends it.  A caller that reports progress observes the evaluations with a
 * FitEvaluations of stage "minos", which marks their records as MINOS's.
 */
void RunMinos(AZURECalc &objective, const AZUREParams &params, double up, FitResult &result, std::ostream &log);

#endif
