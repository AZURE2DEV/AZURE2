#ifndef AZURE2_CONVERT_ABUNDANCELINE_H
#define AZURE2_CONVERT_ABUNDANCELINE_H

#include <sstream>
#include <string>

namespace azure2 {
namespace convert {

/*!
 * One line of the old <abundances> section (docs/source/reference/data_formats.rst, abundances block):
 *
 *   [<sample>:] <Z> <A> <abundance> <vary> <prior sigma>     an abundance of the compound nucleus (Z, A)
 *   [<sample>:] sum <sigma>                                   the width of the sample's sum-to-one prior
 *
 * '#' starts a comment.  The grammar of CNuc::ReadAbundances, moved here so that the engine (until the format-3
 * reader replaces the section, step 3b) and the converter to format 3 read the line with one parser.  What a line
 * means (which compound nuclei exist, the defaults of a nucleus without a line) stays with the reader.
 */
struct AbundanceLine {
  enum Kind {
    Blank,      ///< empty, a comment, or a sample prefix alone: nothing to read
    Sum,        ///< [<sample>:] sum <sigma>
    Entry,      ///< [<sample>:] <Z> <A> <abundance> <vary> <sigma>
    BadSum,     ///< 'sum' without a number
    BadEntry,   ///< neither form
  };
  Kind kind = Blank;
  std::string sample;  ///< "" for the unnamed sample
  double sumSigma = 0.0;
  int z = 0, a = 0, vary = 0;
  double value = 0.0, sigma = 0.0;
};

inline AbundanceLine ReadAbundanceLine(const std::string &line) {
  AbundanceLine out;
  std::string entry = line;
  const std::size_t hash = entry.find('#');
  if (hash != std::string::npos) entry = entry.substr(0, hash);
  if (entry.find_first_not_of(" \t\r") == std::string::npos) return out;
  // optional "<name>:" prefix (one token ending in a colon)
  std::istringstream probe(entry);
  std::string first;
  probe >> first;
  if (!first.empty() && first.back() == ':') {
    out.sample = first.substr(0, first.size() - 1);
    std::string rest;
    std::getline(probe, rest);
    entry = rest;
  }
  std::istringstream stm(entry);
  std::string head;
  if (!(stm >> head)) return out;
  if (head == "sum") {
    double s;
    if (!(stm >> s)) {
      out.kind = AbundanceLine::BadSum;
      return out;
    }
    out.kind = AbundanceLine::Sum;
    out.sumSigma = s;
    return out;
  }
  std::istringstream stm2(entry);
  if (!(stm2 >> out.z >> out.a >> out.value >> out.vary >> out.sigma)) {
    out.kind = AbundanceLine::BadEntry;
    return out;
  }
  out.kind = AbundanceLine::Entry;
  return out;
}

}  // namespace convert
}  // namespace azure2

#endif
