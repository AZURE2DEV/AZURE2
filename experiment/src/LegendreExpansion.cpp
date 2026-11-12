// The normalisations of Legendre coefficients (LegendreExpansion.h).
#include "azure2/experiment/LegendreExpansion.h"

namespace azure2 {
namespace experiment {

LegendreExpansion LegendreExpansionOf(const std::string &word) {
  if (word == "absolute") return LegendreExpansion::Absolute;
  if (word == "rsl") return LegendreExpansion::Rsl;
  if (word == "2l2") return LegendreExpansion::TwoL2;
  if (word == "4pi") return LegendreExpansion::FourPi;
  return LegendreExpansion::AOverA0;
}

const char *LegendreExpansionName(LegendreExpansion e) {
  switch (e) {
    case LegendreExpansion::AOverA0: return "a_over_a0";
    case LegendreExpansion::Absolute: return "absolute";
    case LegendreExpansion::Rsl: return "rsl";
    case LegendreExpansion::TwoL2: return "2l2";
    case LegendreExpansion::FourPi: return "4pi";
  }
  return "a_over_a0";
}

const char *LegendreExpansionFormula(LegendreExpansion e) {
  switch (e) {
    case LegendreExpansion::AOverA0: return "A_L / A_0";
    case LegendreExpansion::Absolute: return "A_L (b/sr)";
    case LegendreExpansion::Rsl: return "A_L / ((2L+1) A_0)";
    case LegendreExpansion::TwoL2: return "2 A_L / (2L+1) (b/sr)";
    case LegendreExpansion::FourPi: return "4 pi A_L (b)";
  }
  return "A_L / A_0";
}

bool LegendreExpansionNormalised(LegendreExpansion e) {
  return e == LegendreExpansion::AOverA0 || e == LegendreExpansion::Rsl;
}

double LegendreExpansionFactor(LegendreExpansion e, int order) {
  const double twoLPlusOne = 2.0 * order + 1.0;
  switch (e) {
    case LegendreExpansion::AOverA0:
    case LegendreExpansion::Absolute: return 1.0;
    case LegendreExpansion::Rsl: return 1.0 / twoLPlusOne;
    case LegendreExpansion::TwoL2: return 2.0 / twoLPlusOne;
    case LegendreExpansion::FourPi: return 4.0 * 3.14159265358979323846;
  }
  return 1.0;
}

}  // namespace experiment
}  // namespace azure2
