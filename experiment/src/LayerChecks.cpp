// The set-up refusals of a yield on a target layer (LayerChecks.h).
#include "azure2/experiment/LayerChecks.h"

#include <sstream>

namespace azure2 {
namespace experiment {

namespace {
Diagnostic Base(const LayerPoint &p, const char *code, const std::string &field) {
  Diagnostic d;
  d.code = code;
  d.section = "targets";
  d.recordKind = "target";
  d.record = p.target;
  d.field = field;
  return d;
}
std::string Point(const LayerPoint &p) {
  std::ostringstream s;
  s.precision(10);
  s << "the yield at E_lab = " << p.beamLab << " MeV" << (p.dataFile.empty() ? "" : " of " + p.dataFile);
  return s.str();
}
}  // namespace

Diagnostic LayerStoppingUndefined(const LayerPoint &p, double atLab, double value) {
  Diagnostic d = Base(p, codes::kLayerStopping, "stopping");
  std::ostringstream s;
  s.precision(10);
  s << Point(p) << ": the stopping is " << value << " at E_lab = " << atLab
    << " MeV, which the integration over the layer reaches (not positive and finite): give a stopping that is positive "
       "at every energy the beam has in the layer (from E_0 down to E_0 minus its energy loss, and a little below), or "
       "a smaller areal density";
  std::ostringstream f;
  f.precision(10);
  f << value;
  d.hasFound = true;
  d.found = f.str();
  d.message = s.str();
  return d;
}

Diagnostic LayerExitBelowModel(const LayerPoint &p, double density, double lowestLab) {
  Diagnostic d = Base(p, codes::kLayerExitBelowModel, "areal_density");
  std::ostringstream s;
  s.precision(10);
  s << Point(p) << ": the beam does not leave a layer of " << density << " atoms/cm2 above E_lab = " << lowestLab
    << " MeV, the lowest energy the model is evaluated at; it stops in the target: use a target of kind infinite "
       "(the thick-target yield), or the areal density the beam crosses";
  std::ostringstream f;
  f.precision(10);
  f << density;
  d.hasFound = true;
  d.found = f.str();
  d.message = s.str();
  return d;
}

}  // namespace experiment
}  // namespace azure2
