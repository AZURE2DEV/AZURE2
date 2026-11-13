// The base a fitted nuisance parameter's grids and limits are laid out from (LayoutBase.h).
#include "azure2/experiment/LayoutBase.h"

#include <sstream>

namespace azure2 {
namespace experiment {

double LayoutSpan::Highest(double base, double sigma, double gridMax) const {
  if (gridMax > 0.0) return gridMax;
  if (!(sigma > 0.0)) return widenWithout * base;
  const double prior = base + 3.0 * sigma, widened = widenWithPrior * base;
  return prior > widened ? prior : widened;
}

LayoutBase ChooseLayoutBase(const LayoutSpan &span, double start, const PriorCentre &centre, bool vary, double sigma,
                            double gridMax, const std::string &name, const std::string &unit) {
  LayoutBase out;
  out.value = start;
  if (!vary || !centre.written || centre.value == start) return out;
  const double c = centre.value, lowest = span.Lowest(c), highest = span.Highest(c, sigma, gridMax);
  if (c > 0.0 && start >= lowest && start <= highest) {
    out.value = c;
    return out;
  }
  std::ostringstream note;
  const std::string u = unit.empty() || unit == "1" ? "" : " " + unit;
  note << "Note: " << name << " starts at " << start << u << ", outside [" << lowest << ", " << highest << "]" << u
       << ", the range laid out from its prior's centre " << c << u
       << "; its grids and limits are laid out from the start.";
  out.note = note.str();
  return out;
}

}  // namespace experiment
}  // namespace azure2
