// The program's choice of fitter (docs/upgrade/fit.md section 5): Config::fitter names a row of Fitters(); with several
// switches the earliest row of the table among them runs, whatever their order and however often each is given (the
// rule of the old bits of Config::paramMask); no switch is the row without one.  The expected row is computed here
// from the table's order, not through ChooseFitter.
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Config.h"
#include "Fitters.h"

namespace {

int failures = 0;
long checks = 0;
void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    failures++;
    std::printf("  FAIL  %s\n", what.c_str());
  }
}

std::string Joined(const std::vector<std::string> &words) {
  std::string out;
  for (const std::string &w : words) out += (out.empty() ? "" : " ") + w;
  return out;
}

}  // namespace

int main() {
  std::ostringstream sink;
  Config configure(sink);
  const std::vector<FitterRow> &rows = Fitters();

  // The rows with a switch, by their place in the table; and the row without one.
  std::vector<std::size_t> switched;
  const FitterRow *byDefault = nullptr;
  for (std::size_t k = 0; k < rows.size(); k++) {
    if (*rows[k].option)
      switched.push_back(k);
    else if (!byDefault)
      byDefault = &rows[k];
  }
  check("the table has a row without a switch", byDefault != nullptr);
  check("the table has a row with a switch", !switched.empty());
  if (!byDefault) return 1;
  check("no switch: the row without one (" + std::string(byDefault->name) + ")",
        configure.fitter.empty() && FitterOfConfig(configure).name == byDefault->name);

  // Every non-empty set of switched rows, in every order, each switch given once and given twice in a row.
  const std::size_t n = switched.size();
  long cases = 0;
  for (unsigned mask = 1; mask < (1u << n); mask++) {
    std::vector<std::size_t> set;
    for (std::size_t i = 0; i < n; i++)
      if (mask & (1u << i)) set.push_back(switched[i]);
    const std::size_t earliest = *std::min_element(set.begin(), set.end());
    std::sort(set.begin(), set.end());
    do {
      for (int repeat = 1; repeat <= 2; repeat++) {
        configure.Reset();
        std::vector<std::string> given;
        for (std::size_t k : set)
          for (int r = 0; r < repeat; r++) {
            given.push_back(rows[k].option);
            check(std::string(rows[k].option) + ": a fitter switch", ApplyFitterSwitch(configure, rows[k].option));
          }
        const std::string ran = FitterOfConfig(configure).name;
        check(Joined(given) + ": runs " + rows[earliest].name + " (got " + ran + ")", ran == rows[earliest].name);
        cases++;
      }
    } while (std::next_permutation(set.begin(), set.end()));
  }
  std::printf("  %ld sets and orders of %zu fitter switches\n", cases, n);

  // A word that is no fitter switch is not taken and leaves the choice alone.
  configure.Reset();
  ApplyFitterSwitch(configure, rows[switched.back()].option);
  const std::string before = configure.fitter;
  check("--covariance-band: no fitter switch", !ApplyFitterSwitch(configure, "--covariance-band"));
  check("--covariance-band: the choice unchanged", configure.fitter == before);

  // The THM rule (src/AZURE2.cpp): migrad-gradient under the same precedence.
  if (FitterNamed("migrad-gradient")) {
    configure.Reset();
    ChooseFitter(configure, "migrad-gradient");
    check("THM alone: migrad-gradient", std::string(FitterOfConfig(configure).name) == "migrad-gradient");
    for (std::size_t k : switched) {
      configure.Reset();
      ApplyFitterSwitch(configure, rows[k].option);
      ChooseFitter(configure, "migrad-gradient");
      const std::string expected = k < (std::size_t)(FitterNamed("migrad-gradient") - &rows[0]) ? rows[k].name
                                                                                                 : "migrad-gradient";
      check(std::string("THM after ") + rows[k].option + ": " + expected,
            std::string(FitterOfConfig(configure).name) == expected);
    }
  } else {
    check("the row migrad-gradient, which the THM rule chooses, exists", false);
  }

  // Reset forgets the choice; a name that is no row is refused.
  configure.Reset();
  check("Reset: no fitter chosen", configure.fitter.empty());
  bool threw = false;
  try {
    ChooseFitter(configure, "no-such-fitter");
  } catch (const std::logic_error &) {
    threw = true;
  }
  check("a name that is no row: std::logic_error, the choice unchanged", threw && configure.fitter.empty());

  // A fitter this build leaves out is refused with how to build it.
  for (const AbsentFitterRow &absent : AbsentFitters()) {
    bool refused = false;
    try {
      ApplyFitterSwitch(configure, absent.option);
    } catch (const std::invalid_argument &) {
      refused = true;
    }
    check(std::string(absent.option) + " (not in this build): refused", refused && configure.fitter.empty());
  }

  // A removed fitter (gsl-lm, 2026-11-18): no row, and its switch is refused with when and what instead, never run
  // by another fitter in its place.
  check("gsl-lm is a removed fitter, not a row", RemovedFitterNamed("gsl-lm") && RemovedFitterNamed("--use-gsl-lm") &&
                                                     !FitterNamed("gsl-lm"));
  for (const RemovedFitterRow &removed : RemovedFitters()) {
    std::string message;
    try {
      ApplyFitterSwitch(configure, removed.option);
    } catch (const std::invalid_argument &e) {
      message = e.what();
    }
    check(std::string(removed.option) + " (removed): refused with 'removed', the date and --use-lm",
          message.find("removed") != std::string::npos && message.find(removed.removed) != std::string::npos &&
              message.find("--use-lm") != std::string::npos && configure.fitter.empty());
  }

  std::printf("%s: %ld checks, %d failed\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
