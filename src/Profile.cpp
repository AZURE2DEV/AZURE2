#include "Profile.h"
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>

namespace AZUREProfile {

namespace {
const char *const kNames[kNumStages] = {"evaluation",     "points",          "E-dependent values",
                                        "matrix fill",    "factor/solve",    "T-matrix",
                                        "channel solves", "cross sections",  "gradient",
                                        "grad forward",   "grad T-cotangent", "grad back-propagation",
                                        "grad level matrix", "grad rows",
                                        "shift cols (segment)", "shift cols (full)",
                                        "energy tasks: set-up"};
const char *const kCounterNames[kNumCounters] = {"channel-space solves", "level route: size",
                                                 "level route: no formal set-up", "level route: near pole",
                                                 "level route: pivot", "level route: forced",
                                                 "channel space, deflated", "final-state W: memo hits",
                                                 "final-state W: computed", "energy tasks (one solve each)",
                                                 "units of the energy tasks", "energy tasks rebuilt"};
}  // namespace

void Account::Report(std::ostream &out, const char *label) {
  std::ostringstream s;
  s << "[profile] " << label << " (ms, thread-summed; calls):";
  for (int k = 0; k < kNumStages; k++) {
    const long long ns = ns_[k].exchange(0);
    const long long calls = calls_[k].exchange(0);
    if (calls == 0) continue;
    s << "\n[profile]   " << std::left << std::setw(22) << kNames[k] << std::right << std::fixed
      << std::setprecision(1) << std::setw(10) << ns / 1.0e6 << "  " << calls;
  }
  for (int k = 0; k < kNumCounters; k++) {
    const long long n = count_[k].exchange(0);
    if (n == 0) continue;
    s << "\n[profile]   " << std::left << std::setw(32) << kCounterNames[k] << std::right << n;
  }
  out << s.str() << std::endl;
}

ScopedReport::ScopedReport(Account *account, const char *label) : account_(account), label_(label), outer_(false) {
  if (account_) outer_ = account_->Enter();
}

ScopedReport::~ScopedReport() {
  if (!account_) return;
  account_->Leave();
  if (outer_) account_->Report(std::cerr, label_);
}

}  // namespace AZUREProfile
