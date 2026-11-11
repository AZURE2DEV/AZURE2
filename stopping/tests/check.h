// Shared by the stopping tests: relative checks, a refusal check, and the reference rows of srim_reference.csv
// (stopping/tools/make_srim_reference.py, LunaSRIM's values).
#ifndef AZURE2_STOPPING_TESTS_CHECK_H
#define AZURE2_STOPPING_TESTS_CHECK_H

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace check {

inline int failures = 0;

inline void Rel(const std::string &what, double got, double want, double rel) {
  const double scale = std::fabs(want) > 0.0 ? std::fabs(want) : 1.0;
  const double dev = std::fabs(got - want) / scale;
  const bool ok = dev <= rel;
  std::printf("%-78s %s  got %.12g want %.12g (rel %.2e, tol %.1e)\n", what.c_str(), ok ? "PASS" : "FAIL", got, want, dev, rel);
  if (!ok) ++failures;
}

inline void True(const std::string &what, bool ok) {
  std::printf("%-78s %s\n", what.c_str(), ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}

/// Passes if f() throws an exception of type E.
template <class E, class F>
void Throws(const std::string &what, F f) {
  bool thrown = false;
  try {
    f();
  } catch (const E &e) {
    thrown = true;
    std::printf("    (refused: %s)\n", e.what());
  } catch (...) {
  }
  True(what, thrown);
}

/// The rows of the reference file whose first field is \p kind, as text fields (the kind dropped).
inline std::vector<std::vector<std::string>> Rows(const std::string &path, const std::string &kind) {
  std::vector<std::vector<std::string>> out;
  std::ifstream in(path);
  if (!in) {
    std::printf("cannot read %s\n", path.c_str());
    ++failures;
    return out;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream fields(line);
    std::string field;
    std::getline(fields, field, ',');
    if (field != kind) continue;
    std::vector<std::string> row;
    while (std::getline(fields, field, ',')) row.push_back(field);
    out.push_back(row);
  }
  if (out.empty()) {
    std::printf("no '%s' rows in %s\n", kind.c_str(), path.c_str());
    ++failures;
  }
  return out;
}

/// The largest relative deviation seen by Track, per label, printed by Summary (the report's table).
struct Largest {
  std::vector<std::pair<std::string, double>> worst;
  void Track(const std::string &label, double got, double want) {
    const double dev = std::fabs(got - want) / (std::fabs(want) > 0.0 ? std::fabs(want) : 1.0);
    for (auto &w : worst) {
      if (w.first == label) {
        if (dev > w.second) w.second = dev;
        return;
      }
    }
    worst.emplace_back(label, dev);
  }
  void Summary(const std::string &title) const {
    std::printf("largest relative deviation, %s:\n", title.c_str());
    for (const auto &w : worst) std::printf("    %-40s %.2e\n", w.first.c_str(), w.second);
  }
};

/// A quiet relative check for long tables: counts a failure and prints only failing rows.
inline void RelQuiet(const std::string &what, double got, double want, double rel) {
  const double dev = std::fabs(got - want) / (std::fabs(want) > 0.0 ? std::fabs(want) : 1.0);
  if (!(dev <= rel)) {
    std::printf("%-78s FAIL  got %.15g want %.15g (rel %.2e, tol %.1e)\n", what.c_str(), got, want, dev, rel);
    ++failures;
  }
}

inline int Done() {
  std::printf("%s: %d failure(s)\n", failures ? "FAILED" : "OK", failures);
  return failures ? 1 : 0;
}

}  // namespace check

#endif
