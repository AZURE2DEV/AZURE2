#include "azure2/experiment/ValueText.h"

#include <charconv>
#include <cmath>
#include <system_error>

#include "azure2/experiment/Units.h"

namespace azure2 {
namespace experiment {

std::string NuisanceText(const NuisanceWords &w) {
  if (w.value.empty()) return w.mode;
  std::string s = w.value + (w.unit.empty() ? "" : " " + w.unit) + " " + w.mode;
  if (!w.priorCentre.empty() || !w.priorWidth.empty()) s += " prior " + w.priorCentre + " " + w.priorWidth;
  if (!w.gridMax.empty()) s += " grid_max " + w.gridMax;
  return s;
}

std::string ShortestNumberText(double x) {
  char buffer[64];
  const std::to_chars_result r = std::to_chars(buffer, buffer + sizeof(buffer), x);
  return r.ec == std::errc() ? std::string(buffer, r.ptr) : std::string("nan");
}

bool ParameterValueEdit(const ValidationResult &validated, const ParameterSource &source, double engineValue,
                        ValueEdit &edit, std::string &why) {
  const ValidatedRecord *record = validated.Find(source.record);
  if (!record) {
    why = "no record '" + source.record + "'";
    return false;
  }
  const std::vector<const Value *> fields = record->GetAll(source.field);
  if (source.occurrence < 0 || source.occurrence >= (int)fields.size()) {
    why = source.record + " has no field " + source.field + " (occurrence " + std::to_string(source.occurrence) + ")";
    return false;
  }
  const Value *v = source.part.empty() ? fields[source.occurrence] : fields[source.occurrence]->Get(source.part);
  if (!v || v->type != ValueType::Nuisance || !v->nuisance.hasValue) {
    why = source.record + "." + source.field + " holds no nuisance with a value";
    return false;
  }
  const Nuisance &n = v->nuisance;
  // the inverse of the reader's ToBase (Units.h): the engine's value in the unit the nuisance is written in
  const double written = !n.unit                     ? engineValue
                         : source.engineUnit.empty() ? Rescale(engineValue, 1.0 / n.unit->factor)
                                                     : engineValue / ToUnit(1.0, n.unit, source.engineUnit);
  if (!std::isfinite(written)) {
    why = "the value is not finite";
    return false;
  }
  edit.line = v->line;
  edit.col = n.value.col;
  edit.old = n.value.text;
  edit.text = ShortestNumberText(written);
  return true;
}

}  // namespace experiment
}  // namespace azure2
