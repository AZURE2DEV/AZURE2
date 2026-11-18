// AZUREParams::Add refuses a second parameter of a name it holds.  Minuit2's MnUserParameters::Add refuses one with a
// warning only, and the record was pushed anyway: the vector then held one parameter fewer than its records, and every
// later Record(i) failed with "N parameters but N+1 records" far from the cause (2026-11-18: two stopping scales named
// stopping_scale.C.1H, scale_names_test.py).  The message names the parameter and both owners.
#include <cstdio>
#include <stdexcept>
#include <string>

#include "AZUREParams.h"

namespace {
int failures = 0;
void Check(const std::string &what, bool ok, const std::string &detail = "") {
  std::printf("  %s  %s%s\n", ok ? "ok  " : "FAIL", what.c_str(), detail.empty() ? "" : ("  (" + detail + ")").c_str());
  if (!ok) failures++;
}
}  // namespace

int main() {
  AZUREParams p;
  p.Add("stopping_scale.C.1H", 1.0, 0.01, ParameterRecord::Apparatus(1, "stopping_scale", "1", {{"target", "tgt_A"}}));
  p.Add("stopping_scale.C_solid.1H", 1.0, 0.01, ParameterRecord::Apparatus(2, "stopping_scale", "1", {{"target", "tgt_B"}}));
  Check("two names, two parameters and two records",
        p.GetMinuitParams().Params().size() == 2 && p.Records().size() == 2);
  std::string message;
  try {
    p.Add("stopping_scale.C.1H", 1.0, 0.01, ParameterRecord::Apparatus(3, "stopping_scale", "1", {{"target", "tgt_C"}}));
  } catch (const std::logic_error &e) {
    message = e.what();
  }
  Check("a second parameter of one name is refused", !message.empty());
  Check("the refusal names the parameter and both owners",
        message.find("'stopping_scale.C.1H'") != std::string::npos && message.find("tgt_A") != std::string::npos &&
            message.find("tgt_C") != std::string::npos,
        message);
  Check("the vector is left as it was", p.GetMinuitParams().Params().size() == 2 && p.Records().size() == 2);
  std::printf(failures ? "duplicate_parameter: %d failure(s)\n" : "duplicate_parameter: all checks passed\n", failures);
  return failures ? 1 : 0;
}
