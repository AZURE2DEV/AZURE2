#ifndef AZURE2_CONVERT_PRIORCENTREROW_H
#define AZURE2_CONVERT_PRIORCENTREROW_H

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

namespace azure2 {
namespace convert {

/*!
 * A prior-centre row of the old <parameterSettings> section (docs/source/reference/data_formats.rst, nominal-norm):
 *
 *   segment_<N>_norm prior_centre <c>          the centre of segment N's normalization prior
 *   segment_<N>_energy_shift prior_centre <c>  the centre of its energy-shift prior
 *
 * Only rows of exactly three tokens whose second token is "prior_centre" are such rows; every other row is not one.
 * The grammar of WorkingCopy::ReadPriorCentres, moved here so that the engine (until format 3 writes the centre in the
 * field, decision F14) and the converter read the row with one parser.
 */
struct PriorCentreRow {
  enum Kind {
    NotARow,      ///< not a prior-centre row (left to the other readers of the section)
    Norm,         ///< segment_<N>_norm prior_centre <c>
    EnergyShift,  ///< segment_<N>_energy_shift prior_centre <c>
    Malformed,    ///< three tokens with prior_centre, but not one of the two forms with a finite number
  };
  Kind kind = NotARow;
  int segment = 0;  ///< N, the segment key (the line's place among the <segmentsData> lines)
  double centre = 0.0;
  std::string centreText;  ///< the number as written
};

/// `trimmed`: the row without leading and trailing blanks.
inline PriorCentreRow ReadPriorCentreRow(const std::string &trimmed) {
  PriorCentreRow out;
  std::istringstream stm(trimmed);
  std::vector<std::string> tok;
  std::string t;
  while (stm >> t) tok.push_back(t);
  if (tok.size() != 3 || tok[1] != "prior_centre") return out;

  const std::string &name = tok[0];
  bool isNorm = false;
  std::string number;
  const std::string normTail = "_norm", shiftTail = "_energy_shift";
  if (name.compare(0, 8, "segment_") == 0 && name.size() > 8 + normTail.size() &&
      name.compare(name.size() - normTail.size(), normTail.size(), normTail) == 0) {
    isNorm = true;
    number = name.substr(8, name.size() - 8 - normTail.size());
  } else if (name.compare(0, 8, "segment_") == 0 && name.size() > 8 + shiftTail.size() &&
             name.compare(name.size() - shiftTail.size(), shiftTail.size(), shiftTail) == 0) {
    number = name.substr(8, name.size() - 8 - shiftTail.size());
  }
  const bool digits = !number.empty() && number.find_first_not_of("0123456789") == std::string::npos;
  char *end = nullptr;
  const double centre = std::strtod(tok[2].c_str(), &end);
  if (!digits || end == tok[2].c_str() || *end != '\0' || !std::isfinite(centre)) {
    out.kind = PriorCentreRow::Malformed;
    return out;
  }
  out.kind = isNorm ? PriorCentreRow::Norm : PriorCentreRow::EnergyShift;
  out.segment = std::atoi(number.c_str());
  out.centre = centre;
  out.centreText = tok[2];
  return out;
}

}  // namespace convert
}  // namespace azure2

#endif
