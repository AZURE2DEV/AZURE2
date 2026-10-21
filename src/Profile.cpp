#include "Profile.h"
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>

namespace AZUREProfile {

namespace {
std::atomic<long long> g_ns[kNumStages];
std::atomic<long long> g_calls[kNumStages];
const char *const kNames[kNumStages] = {"evaluation",     "points",          "E-dependent values",
                                        "matrix fill",    "factor/solve",    "T-matrix",
                                        "channel solves", "cross sections",  "gradient",
                                        "grad forward",   "grad T-cotangent", "grad back-propagation",
                                        "grad level matrix", "grad rows"};
}  // namespace

void Add(Stage stage, long long ns) {
  g_ns[stage].fetch_add(ns, std::memory_order_relaxed);
  g_calls[stage].fetch_add(1, std::memory_order_relaxed);
}

void Report(std::ostream &out, const char *label) {
  if (!Enabled()) return;
  std::ostringstream s;
  s << "[profile] " << label << " (ms, thread-summed; calls):";
  for (int k = 0; k < kNumStages; k++) {
    const long long ns = g_ns[k].exchange(0);
    const long long calls = g_calls[k].exchange(0);
    if (calls == 0) continue;
    s << "\n[profile]   " << std::left << std::setw(22) << kNames[k] << std::right << std::fixed
      << std::setprecision(1) << std::setw(10) << ns / 1.0e6 << "  " << calls;
  }
  out << s.str() << std::endl;
}

namespace {
std::atomic<int> g_depth{0};
}

ScopedReport::ScopedReport(const char *label) : label_(label), outer_(false) {
  if (Enabled()) outer_ = (g_depth.fetch_add(1) == 0);
}

ScopedReport::~ScopedReport() {
  if (!Enabled()) return;
  g_depth.fetch_sub(1);
  if (outer_) Report(std::cerr, label_);
}

}  // namespace AZUREProfile
