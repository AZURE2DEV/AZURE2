#ifndef AZURE2_CONVERT_LINEREPAIR_H
#define AZURE2_CONVERT_LINEREPAIR_H

#include <string>
#include <vector>

namespace azure2 {
namespace convert {

/*!
 * A line of an old project's <segmentsData>, <segmentsTest> or <targetInt> that the old grammar cannot read
 * (docs/upgrade/qt_gui.md section 19.5): where it is, what it holds, what the grammar expects, and the readings the
 * converter can propose.  A reading is only a proposal: the converter never applies one by itself (the line may hold a
 * normalization, an energy shift or a data type, and a guessed value would change the fit without a word).  The user
 * takes one, or edits the line, in a copy of the project and converts the copy; AZURE2 --convert prints the same
 * diagnosis as text, and the Qt program offers the readings as buttons with the line editable in place.
 */
struct LineReading {
  std::string label;  ///< what the reading assumes, in words ("the normalization's error is missing: 0, no prior")
  std::string line;   ///< the line with that assumption written out (the old grammar reads it)
};

struct LineProblem {
  std::string section;  ///< "segmentsData", "segmentsTest" or "targetInt"
  int fileLine = 0;     ///< 1-based line number in the file as given; 0: the refusal is not about one such line
  std::string text;     ///< the line as written in the file
  std::string problem;  ///< the refusal (what the old grammar could not read)
  std::string expected; ///< the fields the old grammar reads on such a line, in order
  std::vector<LineReading> readings;
  bool Found() const { return fileLine > 0; }
  /// The diagnosis as text, a few lines: the line and its number, the fields expected, each reading.
  std::string Text() const;
};

/// Whether `token` is a number written with a decimal comma ("0,005", "-1,5e-3").  On a data line such a token is
/// read by the old grammar as the data file's name (the numbers before it end there), so the line is refused.
bool DecimalCommaNumber(const std::string &token);

/// Whether `token` is two numbers run together by a fixed-width writer: a number with more than 17 significant digits
/// (more than any double needs, so no writer of a single number produced it) or with a leading zero before a digit
/// ("013.0057").  The setup program of 2026-09-25 .. 2026-10-07 wrote a round-trip number as wide as its column with no
/// space after it (fixed by 2887d991 and 3e7bd976): "114.918600000000001" for J = 1 and E = 14.918600000000001.
bool GluedNumbers(const std::string &token);

/// The fields the old grammar reads on a line of `section`, in order, optional parts in brackets.
std::string ExpectedFields(const std::string &section);

/// Whether the old grammar reads `line` of `section` whole: no failed read, nothing left unread, a data line naming a
/// data file.
bool OldLineReads(const std::string &section, const std::string &line);

/// The readings the converter can propose for a line of `section` that the old grammar cannot read; each proposed line
/// reads (OldLineReads).  Empty when no reading is evident (the user then edits the line).
std::vector<LineReading> ProposeReadings(const std::string &section, const std::string &line);

}  // namespace convert
}  // namespace azure2

#endif
