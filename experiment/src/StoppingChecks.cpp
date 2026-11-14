// The set-up refusals of a target's stopping from the stopping/ library (StoppingChecks.h).
#include "azure2/experiment/StoppingChecks.h"

#include <sstream>

namespace azure2 {
namespace experiment {

namespace {
Diagnostic Base(const std::string &target, const char *code) {
  Diagnostic d;
  d.code = code;
  d.section = "targets";
  d.recordKind = "target";
  d.record = target;
  d.field = "stopping";
  return d;
}
}  // namespace

Diagnostic StoppingDataMissing(const std::string &target, const std::string &why) {
  Diagnostic d = Base(target, codes::kStoppingData);
  d.message = "the stopping cannot be evaluated: " + why +
              " (an evaluation needs SRIM-2013's SCOEF03.DAT and SNUC03.dat, compiled in from stopping/data/srim/ or in "
              "the directory --srim-data-dir=DIR (the switch AZURE2_SRIM_DATA_DIR) names; a table needs its file)";
  return d;
}

Diagnostic StoppingOutOfRange(const std::string &target, const std::string &projectile, double atLab, const std::string &why) {
  Diagnostic d = Base(target, codes::kStoppingRange);
  std::ostringstream s, f;
  s.precision(10);
  f.precision(10);
  f << atLab;
  s << "the stopping of " << projectile << " is asked for at E_lab = " << atLab
    << " MeV, which a point's integration over the target reaches, outside its range: " << why
    << "; give a table that covers every energy the beam has in the target, or data inside the range";
  d.hasFound = true;
  d.found = f.str();
  d.message = s.str();
  return d;
}

}  // namespace experiment
}  // namespace azure2
