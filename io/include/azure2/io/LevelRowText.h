#ifndef AZURE2_IO_LEVELROWTEXT_H
#define AZURE2_IO_LEVELROWTEXT_H

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace azure2 {
namespace io {

/*!
 * The text of a <levels> row (the grammar is LevelRows.h): its tokens by column name, and the row with some of them
 * replaced.  The writer of a row: an edit keeps every byte of the line but the tokens it replaces, each new token
 * right-aligned in the width of the old one (a wider one widens the line), so one edit does not reflow the block.
 *
 * Columns, file order: the 31 of every row (levelJ levelPi levelE levelFix aa ir s2 l2 levelID isActive channelFix
 * gamma j1 pi1 j2 pi2 e2 m1 m2 z1 z2 entranceSepE sepE j3 pi3 e3 pType chRad g1 g2 ecMultMask; s2 and l2 are the
 * doubled channel spin and l as written), then the optional ones: bindingE and gammaIsRWA (written together), e1
 * (written last).  The names are LevelRow's accessors.
 */
const std::vector<std::string> &LevelColumnNames();
/// The index of a column name in LevelColumnNames(), or -1.
int LevelColumnIndex(const std::string &name);

/// One row as text.
struct LevelRowText {
  std::string line;  ///< the line as written (without its ending)
  /// (first byte, one past the last) of the token of each column; (npos, npos) for an optional column not written
  std::vector<std::pair<std::size_t, std::size_t>> spans;
  bool hasThm = false;  ///< bindingE and gammaIsRWA are written
  bool hasE1 = false;   ///< e1 is written
  /// The token of a column as written; "" for an optional column the row does not write (its value is 0).
  std::string Token(int column) const;
};

/// Splits a row into its columns: false with `why` when it has fewer than 31 or more than 34 tokens (the reader,
/// LevelRow, decides whether the tokens are numbers).
bool SplitLevelRow(const std::string &line, LevelRowText &out, std::string &why);

/*!
 * The row with the tokens of `tokens` (column index -> new token) replaced.  An optional column the row does not write
 * is appended with its group: bindingE and gammaIsRWA after the 31st token (the other one "0" when not given), e1 at
 * the end (two blanks before the THM pair, one before e1 alone, as the converter and the GUI write them).  The tokens
 * are written as given (the caller formats numbers).
 */
std::string EditLevelRow(const LevelRowText &row, const std::map<int, std::string> &tokens);

}  // namespace io
}  // namespace azure2

#endif
