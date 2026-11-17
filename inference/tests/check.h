// The tests' one helper: a named check that prints its line and counts the failures.
#ifndef AZURE2_INFERENCE_TESTS_CHECK_H
#define AZURE2_INFERENCE_TESTS_CHECK_H

#include <cstdio>
#include <string>

namespace {
int failures = 0;
void check(const std::string &name, bool ok) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", name.c_str());
  if (!ok) failures++;
}
int result() {
  std::printf("%s\n", failures ? "FAILED" : "passed");
  return failures ? 1 : 0;
}
}  // namespace

#endif
