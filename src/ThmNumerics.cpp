#include "ThmNumerics.h"

#include <sstream>

std::string ThmNumberText(double x) {
  std::ostringstream s;
  s.precision(6);
  s << x;
  return s.str();
}
