// A <targetInt> line that ends inside its first nine fields (active, segments, points, isConvolution, sigma,
// isTargetIntegration, density, stopping, nParams): TargetEffect reads nothing more (TargetIntLine.complete is false)
// and, before step 3a, left its members uninitialised -- IsActive(), GetSigma() and the rest read indeterminate values.
// Each member now has an in-class initialiser; this test reads every accessor of such a line and checks the
// initialised values, so that a sanitizer or valgrind run of it sees no indeterminate read, and checks that a
// complete line still sets them (the initialisers are values the complete path never relies on).
//   target_effect_line_test
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

#include "Config.h"
#include "TargetEffect.h"
#include "azure2/convert/SpecsFromOldLines.h"

namespace {

int failures = 0;
void check(const std::string &what, bool ok) {
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

TargetEffect Read(const Config &c, const std::string &line) {
  std::istringstream in(line);
  bool readFailed = false;
  return TargetEffect(azure2::convert::TargetEffectSpecFromOldLine(in, readFailed), c);
}

}  // namespace

int main() {
  std::ostringstream log;
  Config configure(log);
  for (const std::string line : {"", "1", "1 \"1\" 20", "1 \"1\" 20 1 0.001", "1 \"1\" 20 1 0.001 1 5e17 \"a0\""}) {
    const TargetEffect e = Read(configure, line);
    const std::string w = "truncated line '" + line + "': ";
    check(w + "inactive", !e.IsActive());
    check(w + "no convolution", !e.IsConvolution() && !e.IsConvCoefficients() && e.GetSigma() == 0.0);
    check(w + "no target integration", !e.IsTargetIntegration() && e.GetDensity() == 0.0);
    check(w + "no Q coefficients", !e.IsQCoefficients());
    check(w + "no sub-points", e.NumSubPoints() == 0);
    check(w + "no straggling", !e.IsStraggling() && e.GetStragglingCoefficient() == 0.0);
    check(w + "grid defaults", e.GetResonanceWidthMultiplier() > 0.0 && e.GetPointsPerWidth() > 0.0);
  }
  const TargetEffect full = Read(configure, "1 \"1\" 20 1 0.001 1 5e17 \"a0\" 1 1e-20 0 0 0 \"\" 0 1 0.6 10 20");
  check("complete line: active, convolution, integration", full.IsActive() && full.IsConvolution() && full.IsTargetIntegration());
  check("complete line: sigma, density, sub-points", full.GetSigma() == 0.001 && full.GetDensity() == 5e17 && full.NumSubPoints() == 20);
  check("complete line: straggling and grid fields", full.IsStraggling() && full.GetStragglingCoefficient() == 0.6 &&
                                                         full.GetResonanceWidthMultiplier() == 10 && full.GetPointsPerWidth() == 20);
  std::printf("target_effect_line: %d failure(s)\n", failures);
  return failures ? 1 : 0;
}
