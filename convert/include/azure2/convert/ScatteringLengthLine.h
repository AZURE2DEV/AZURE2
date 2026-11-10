#ifndef AZURE2_CONVERT_SCATTERINGLENGTHLINE_H
#define AZURE2_CONVERT_SCATTERINGLENGTHLINE_H

#include <sstream>
#include <string>

namespace azure2 {
namespace convert {

/*!
 * One line of the old <scatteringLengths> section (include/ScatteringLength.h):
 *
 *   <active> <pair key> <s> <Re a (fm)> <sigma (fm)> [<Im a (fm)> <sigma (fm)>]
 *
 * s < 0 is the coherent combination; '#' starts a comment.  The grammar of ScatteringLength's ReadConstraints,
 * moved here so that the engine (until the format-3 reader replaces the section, step 3b) and the converter to
 * format 3 read the line with one parser.  What it means (the pair, its channels) stays with the reader.
 */
struct ScatteringLengthLine {
  enum Kind { Blank, Constraint, Malformed };
  Kind kind = Blank;
  std::string error;  ///< for Malformed: what is wrong, in the engine's words
  int active = 0;
  int pairKey = 0;
  double s = 0.0;
  double value = 0.0, sigma = 0.0;
  bool hasImaginary = false;
  double imValue = 0.0, imSigma = 0.0;
};

inline ScatteringLengthLine ReadScatteringLengthLine(std::string line) {
  ScatteringLengthLine out;
  const std::size_t hash = line.find('#');
  if (hash != std::string::npos) line = line.substr(0, hash);
  std::istringstream stm(line);
  if (!(stm >> out.active)) return out;  // blank or comment
  out.kind = ScatteringLengthLine::Malformed;
  if (!(stm >> out.pairKey >> out.s >> out.value >> out.sigma)) {
    out.error = "expected 'active pairKey s value sigma [imValue imSigma]'";
    return out;
  }
  double im, imSigma;
  if (stm >> im) {
    if (!(stm >> imSigma)) {
      out.error = "an imaginary part needs its uncertainty too";
      return out;
    }
    out.hasImaginary = true;
    out.imValue = im;
    out.imSigma = imSigma;
  }
  std::string extra;
  if (stm >> extra) {
    out.error = "unexpected '" + extra + "'";
    return out;
  }
  out.kind = ScatteringLengthLine::Constraint;
  return out;
}

}  // namespace convert
}  // namespace azure2

#endif
