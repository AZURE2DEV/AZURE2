// AZUREParams::ReadUserParameters on the files a fit can write.
//
// A fit whose covariance is invalid writes "nan" in the error column of
// param.sav.  The reader used operator>>, which cannot parse "nan": the stream
// went fail() without eof() and the loop never ended (std::bad_alloc).  Checked
// here: nan and inf errors, a blank line, a line without a number, a last line
// without a newline, and the #parametrization line.  ctest gives it a 10 s
// timeout, so a hang fails.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "AZUREParams.h"
#include "Config.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

static int failures = 0;
static void check(bool ok, const char *what) {
  if (!ok) {
    failures++;
    std::printf("  FAIL %s\n", what);
  }
}

int main() {
  const char *file = "parameter_file_reference_test.sav";
  {
    std::ofstream out(file);
    out << "    #parametrization       1.0000000e+00       0.0000000e+00\n"
        << "            energy_1      -6.7442187e+00                 nan\n"
        << "           width_1_1      -4.1960559e+00                -nan\n"
        << "\n"
        << "            energy_2       2.5000000e+00                 inf\n"
        << "           not_a_row                 abc                 1.0\n"
        << "           width_2_1       1.2500000e-01       3.0000000e-03";  // no newline
  }
  AZUREParams params;
  auto &mp = params.GetMinuitParams();
  mp.Add("energy_1", 0.0, 0.5);
  mp.Add("width_1_1", 1.0, 0.25);
  mp.Add("energy_2", 1.0, 0.75);
  mp.Add("width_2_1", 1.0, 0.125);
  params.ReadUserParameters(std::string(file));
  std::remove(file);

  check(params.BasisTag() == AZUREParams::kBasisBrune, "#parametrization read as 1 (Brune)");
  check(mp.Value(0) == -6.7442187, "energy_1 value read past a nan error");
  check(mp.Error(0) == 0.5, "nan error leaves the step size");
  check(mp.Value(1) == -4.1960559, "width_1_1 value read past a -nan error");
  check(mp.Error(1) == 0.25, "-nan error leaves the step size");
  check(mp.Value(2) == 2.5, "energy_2 after a blank line");
  check(mp.Error(2) == 0.75, "inf error leaves the step size");
  check(mp.Value(3) == 0.125, "last line without a newline");
  check(mp.Error(3) == 3.0e-3, "a finite error is applied");

  if (failures) {
    std::printf("FAIL: %d parameter-file checks\n", failures);
    return 1;
  }
  std::printf("PASS: parameter file with nan errors, blank and malformed lines\n");
  return 0;
}
