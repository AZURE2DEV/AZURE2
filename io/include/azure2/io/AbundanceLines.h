#ifndef AZURE2_IO_ABUNDANCELINES_H
#define AZURE2_IO_ABUNDANCELINES_H

#include <istream>
#include <string>
#include <vector>

namespace azure2 {
namespace io {

/*!
 * One line of the optional <abundances> section.  The grammar, after an optional "<sample>:" prefix (one token ending
 * in a colon; none = the unnamed sample) and with "#" starting a comment:
 *
 *   <Z> <A> <abundance> <vary> <prior sigma>     the abundance of one compound nucleus (sigma 0 = no prior)
 *   sum <sigma>                                  the width of the sample's sum-to-one penalty
 *
 * Further tokens after the numbers are ignored.  A line holding only a prefix is skipped.
 */
struct AbundanceLine {
  std::string sample;     ///< "" for the unnamed sample
  bool isSum = false;     ///< a "sum <sigma>" line
  double sumSigma = 0.0;  ///< its sigma
  int Z = 0, A = 0;       ///< an abundance line: the compound nucleus
  double value = 0.0;     ///< its abundance
  int vary = 0;           ///< 1: a fit parameter
  double sigma = 0.0;     ///< its prior width; 0 = no prior
  std::string text;       ///< the line as written
};

/// The lines of an <abundances> section, in file order.
struct AbundanceSection {
  enum class Problem {
    None,
    BadSum,   ///< a "sum" line without a number (badLine)
    BadEntry  ///< a line that is neither a sum nor five numbers (badLine)
  };
  std::vector<AbundanceLine> lines;  ///< every line up to the first malformed one
  Problem problem = Problem::None;   ///< the lines after a malformed one are not read
  std::string badLine;
  bool complete = false;             ///< the section ended with its </abundances> line
};

/// Reads the <abundances> section as Config::SectionStream gives it: the first line is the <abundances> tag.
AbundanceSection ReadAbundanceSection(std::istream &in);

}  // namespace io
}  // namespace azure2

#endif
