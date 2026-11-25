#ifndef AZURE2_CONVERT_OLDERA_H
#define AZURE2_CONVERT_OLDERA_H

#include <string>
#include <utility>
#include <vector>

#include "azure2/io/Document.h"

namespace azure2 {
namespace convert {

/*!
 * The forms of an old project that the readers of later eras no longer read, each brought to the form that has the
 * same meaning (docs/upgrade/project_format3.md section 8.5; the history in tests/convert/fixtures/history/README.md).
 * Each appends (construct, note) pairs for the report; a non-empty return is the refusal.
 */

/// An AZURE2 before 2020: external capture in an <externalCapture> section ("isActive exitPairKey minJ maxJ multMask")
/// and <levels> rows of 30 columns.  Folded into the 31st column ecMultMask as the 2020 setup program imported such a
/// file; the section is removed.  An inactive line, or two masks for one pair, is refused.
std::string FoldOldExternalCapture(io::Document &doc, std::vector<std::pair<std::string, std::string>> &notes);

/// A <lastRun> section that is not closed at the end of the file (no engine read it): the closing tag added.
std::string CloseOldLastRun(std::string &text, std::vector<std::pair<std::string, std::string>> &notes);

/// A <targetInt> line of 2025-01-26 .. 01-30 whose convolution coefficients have no equation token: the line with the
/// polynomial equation that gives that era's sigma(E) inserted (`note` says so); the line itself otherwise.
std::string TargetLineOfFirstEra(const std::string &line, std::string &note);

}  // namespace convert
}  // namespace azure2

#endif
