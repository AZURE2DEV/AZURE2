// A tabulated reference cross section (ReferenceTable.h).
#include "azure2/experiment/ReferenceTable.h"

#include <algorithm>
#include <cmath>
#include <istream>
#include <limits>
#include <sstream>

namespace azure2 {
namespace experiment {

std::string ReadReferenceRows(std::istream &in, std::vector<double> &energy, std::vector<double> &value) {
  energy.clear();
  value.clear();
  std::string line;
  int number = 0;
  while (std::getline(in, line)) {
    number++;
    const std::string text = line.substr(0, line.find('#'));
    std::istringstream s(text);
    double e = 0.0, v = 0.0;
    if (!(s >> e)) continue;  // blank or comment
    std::ostringstream where;
    where << "line " << number << ": ";
    if (!(s >> v)) return where.str() + "fewer than two numbers";
    if (!std::isfinite(e) || !std::isfinite(v) || !(v > 0.0)) {
      std::ostringstream m;
      m << where.str() << "the cross section " << v << " is not positive and finite";
      return m.str();
    }
    if (!energy.empty() && !(e > energy.back())) {
      std::ostringstream m;
      m << where.str() << "the energy " << e << " does not increase (previous " << energy.back() << ")";
      return m.str();
    }
    energy.push_back(e);
    value.push_back(v);
  }
  if (energy.size() < 2) return "fewer than two rows";
  return "";
}

double ReferenceAt(const ReferenceTable &table, double eLab_MeV) {
  const std::vector<double> &x = table.energyLab_MeV, &y = table.sigma_b;
  if (x.size() < 2 || !(eLab_MeV >= x.front() && eLab_MeV <= x.back())) return std::numeric_limits<double>::quiet_NaN();
  std::size_t k = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), eLab_MeV) - x.begin());
  k = std::min(std::max<std::size_t>(k, 1), x.size() - 1);  // the interval [x[k-1], x[k]], the last for E = x.back()
  if (eLab_MeV == x[k - 1]) return y[k - 1];
  if (eLab_MeV == x[k]) return y[k];
  const double t = (eLab_MeV - x[k - 1]) / (x[k] - x[k - 1]);
  return std::exp((1.0 - t) * std::log(y[k - 1]) + t * std::log(y[k]));
}

}  // namespace experiment
}  // namespace azure2
