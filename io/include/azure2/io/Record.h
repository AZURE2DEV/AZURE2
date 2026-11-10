#ifndef AZURE2_IO_RECORD_H
#define AZURE2_IO_RECORD_H

#include <string>
#include <vector>

namespace azure2 {
namespace io {

/*!
 * Records of a format-3 section (docs/upgrade/project_format3.md section 3.1, decisions F1, F2): a header line
 * `<kind> <name>` followed by field lines `<key> <value tokens>`, up to the next header or the end of the section.
 * Plain values only, header-only, so that experiment/ (the table that knows kinds and fields) can take them without
 * linking io/.  RecordReader.h reads and writes them.
 *
 * Lines are the file's line numbers (1-based); columns are 1-based byte offsets in the line.
 */
struct Token {
  std::string text;     ///< the token as written, without the quotes of a quoted token
  int col = 0;          ///< column of its first character (of the opening quote, for a quoted token)
  bool quoted = false;  ///< written as "..." (a path or formula with blanks, or any token its writer quoted)
};

struct Field {
  std::string key;
  std::vector<Token> tokens;          ///< the value tokens after the key, in order (may be empty)
  int line = 0;
  int col = 0;                        ///< column of the key
  std::vector<std::string> comments;  ///< whole-line comments directly above the field ("# ..."), in order
  std::string trailingComment;        ///< "# ..." after the value on the same line, or empty
};

struct Record {
  std::string kind;
  std::string name;                   ///< empty when the header has no name (the reader reports E_HEADER_NAME)
  int line = 0;                       ///< line of the header
  int col = 0;                        ///< column of the kind
  int nameCol = 0;                    ///< column of the name (0 without a name)
  std::vector<Field> fields;          ///< in file order; a key given twice appears twice (the validator reports it)
  std::vector<std::string> comments;  ///< whole-line comments directly above the header
  std::string trailingComment;        ///< "# ..." after the header on the same line
};

/// The records of one section, with the comments that belong to no record.
struct RecordSet {
  std::vector<Record> records;
  std::vector<std::string> trailingComments;  ///< whole-line comments after the last field of the last record
};

/// A syntax problem of a section's text, found by the reader.  `code` is one of the codes below; the schema
/// (experiment/, `Schema::ToJson`) lists them with the validator's.
struct RecordDiagnostic {
  std::string code;
  int line = 0;
  int col = 0;
  std::string found;    ///< the offending text
  std::string message;  ///< names the line's content and says what is allowed
};

/// Codes of the reader's diagnostics.
namespace record_codes {
constexpr const char *kUnterminatedQuote = "E_UNTERMINATED_QUOTE";  ///< a '"' without its closing '"' on the line
constexpr const char *kQuoteInToken = "E_QUOTE_IN_TOKEN";           ///< a '"' inside a token or right after one
constexpr const char *kUnknownKind = "E_UNKNOWN_KIND";              ///< a line before any header that is no header
constexpr const char *kHeaderName = "E_HEADER_NAME";                ///< a header without a name, or with more tokens
constexpr const char *kKeySyntax = "E_KEY_SYNTAX";                  ///< a field key that is not [A-Za-z][A-Za-z0-9_]*
}  // namespace record_codes

}  // namespace io
}  // namespace azure2

#endif
