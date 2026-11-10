#ifndef AZURE2_EXPERIMENT_PARSERS_H
#define AZURE2_EXPERIMENT_PARSERS_H

#include <cstddef>
#include <string>
#include <vector>

#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/Schema.h"
#include "azure2/experiment/Values.h"
#include "azure2/io/Record.h"

namespace azure2 {
namespace experiment {

/*!
 * The parsers of the value types (Values.h): exactly one per type, used by Schema::Validate and by whoever needs a
 * value read the same way (step 3's Build reads the ValidatedRecord they produce; pyazr calls the validator).  Each
 * reads tokens from a cursor and reports problems into the context; the first problem of a field ends its parsing, so
 * one wrong token gives one diagnostic.
 *
 *   ValueParsers.cpp   numbers, integers, units, nuisances, enums, names, spins, nuclei, pairs, compositions,
 *                      formulas, number lists, windows, and the structured fields (Form)
 *   ReactionParser.cpp reaction expressions (section 3.5)
 *   ColumnsGrid.cpp    column roles (section 3.4.1) and extrapolation grids (section 3.8)
 */
struct TokenCursor {
  const std::vector<io::Token> *tokens = nullptr;
  std::size_t pos = 0;
  int line = 0;
  int endCol = 0;  ///< the column just after the field's last token (where a missing token is reported)

  bool AtEnd() const { return pos >= tokens->size(); }
  const io::Token &Peek() const { return (*tokens)[pos]; }
  const io::Token &Next() { return (*tokens)[pos++]; }
};

/// Where the diagnostics of one field go, with the record's context filled in.
struct ParseContext {
  std::vector<Diagnostic> *out = nullptr;
  Diagnostic base;                    ///< section, record kind, record, field, line
  const Quantity *quantity = nullptr; ///< for units "@quantity" (nullptr: the quantity is unknown, any unit passes)
  bool failed = false;

  void Error(const char *code, int col, const std::string &found, std::vector<std::string> allowed,
             const std::string &message, const Fix *fix = nullptr);
  /// "dataset 's01', field 'norm'" -- the start of every message.
  std::string Where() const;
};

/// Parse a value of `spec` from the cursor.  The caller checks that the cursor is at its end afterwards.
bool ParseValue(const ValueSpec &spec, TokenCursor &cur, ParseContext &ctx, Value &out);

// Pure helpers: text -> value, no diagnostics.
bool ParseNumberText(const std::string &text, Number &out);
bool IsName(const std::string &text);                  ///< [A-Za-z][A-Za-z0-9_-]* (section 3.3)
bool ParseSpinText(const std::string &text, bool withParity, Spin &out);
/// <A><Symbol> (13C) or, when `light` is set, a light-particle letter (n p d t h a g e- b- b+); `error` says why not.
bool ParseParticleText(const std::string &text, bool light, Particle &out, std::string &error, const char *&code);

// The structured and composite types, each in its own file.
bool ParseReaction(TokenCursor &cur, ParseContext &ctx, ReactionExpr &out);
bool ParseColumns(TokenCursor &cur, ParseContext &ctx, Columns &out);
bool ParseGrid(const ValueSpec &spec, TokenCursor &cur, ParseContext &ctx, Grid &out);

/// The column roles with a unit in brackets (energy or angle), and their dimension.
const std::vector<std::pair<std::string, std::string>> &AxisRoles();
/// The column roles without a unit.
const std::vector<std::string> &PlainRoles();

}  // namespace experiment
}  // namespace azure2

#endif
