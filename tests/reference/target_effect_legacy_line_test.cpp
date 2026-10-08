/*!
 * A <targetInt> line written by an older AZURE2 stops after the attenuation (Q)
 * coefficients: it has no convolution-coefficient block.  The reader extracts
 * those fields anyway, the extraction fails at the end of the line, and a
 * failed extraction leaves its target untouched -- so the coefficient count
 * used to be read from uninitialized memory and drive the read loop.  Depending
 * on the stack contents the line then parsed, was rejected, or allocated until
 * std::bad_alloc (or the OOM killer) ended the run.  Seen on the 13C(p,g)
 * project, whose .azr was saved by an older GUI.
 *
 * Run:  tests/reference/target_effect_legacy_line_test   (ctest: target_effect_legacy_line)
 */
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "Config.h"
#include "TargetEffect.h"

Config *g_config = nullptr;

namespace {

int failures = 0;

void check(const std::string &name, bool ok) {
  std::cout << "  " << (ok ? "ok  " : "FAIL") << "  " << name << std::endl;
  if (!ok) failures++;
}

// Leave a recognisable non-zero pattern where the constructor's locals will
// live, so an uninitialized count is large rather than whatever happened to be
// there.  Not a guarantee, but it turns the old behaviour into a reliable
// failure on the usual x86-64 builds.
__attribute__((noinline)) void dirtyStack() {
  volatile unsigned char junk[16384];
  for (unsigned i = 0; i < sizeof(junk); i++) junk[i] = 0x3f;
}

void legacyLine(const Config &configure, const std::string &label, const std::string &line) {
  dirtyStack();
  std::istringstream in(line);
  TargetEffect effect(in, configure);
  bool parsed = !(in.rdstate() & (std::stringstream::failbit | std::stringstream::badbit));
  check(label + ": the line is accepted", parsed);
  check(label + ": it is active", effect.IsActive());
  check(label + ": no convolution", !effect.IsConvolution());
  check(label + ": no target integration", !effect.IsTargetIntegration());
  check(label + ": Q coefficients are on", effect.IsQCoefficients());
  check(label + ": Q coefficients are read",
        std::fabs(effect.GetQCoefficient(0) - 1.0) < 1e-12 &&
            std::fabs(effect.GetQCoefficient(1) - 1.0) < 1e-12 &&
            std::fabs(effect.GetQCoefficient(2) - 0.88) < 1e-12);
  check(label + ": no convolution coefficients", !effect.IsConvCoefficients());
  std::vector<int> segments = effect.GetSegmentsList();
  check(label + ": it applies to segment 19", segments.size() == 1 && segments[0] == 19);
}

}  // namespace

int main() {
  std::ostringstream sink;
  Config configure(sink);
  g_config = &configure;

  std::cout << "legacy <targetInt> line" << std::endl;
  const std::string legacy =
      "1              \"19\"           10             0              0              0              0"
      "               \"\" 0               1              3 1 1 0.88";
  legacyLine(configure, "trailing space", legacy + " ");
  legacyLine(configure, "no trailing space", legacy);

  // The current format, with an empty convolution block, still reads.
  dirtyStack();
  std::istringstream in(legacy + " 0 \"\" 0");
  TargetEffect current(in, configure);
  check("current format: the line is accepted",
        !(in.rdstate() & (std::stringstream::failbit | std::stringstream::badbit)));
  check("current format: Q coefficients are read", std::fabs(current.GetQCoefficient(2) - 0.88) < 1e-12);

  std::cout << (failures ? "FAILED" : "passed") << std::endl;
  return failures ? 1 : 0;
}
