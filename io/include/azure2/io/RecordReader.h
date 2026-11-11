#ifndef AZURE2_IO_RECORDREADER_H
#define AZURE2_IO_RECORDREADER_H

#include <string>
#include <vector>

#include "azure2/io/Document.h"
#include "azure2/io/Record.h"

namespace azure2 {
namespace io {

/*!
 * The generic reader of format-3 record sections (docs/upgrade/project_format3.md sections 3.1-3.3, 7.3).  It knows
 * no field and no kind of its own: the caller passes the record kinds of the section (experiment/'s table), and a
 * line is a header when its first token is one of them, a field line otherwise (kinds and keys are disjoint, F1).
 *
 * Lexical rules: blanks (spaces, tabs) separate tokens and indentation is not significant; '#' outside quotes starts
 * a comment that runs to the end of the line; a token that starts with '"' runs to the next '"' (no escapes) and may
 * hold blanks and '#'; blank lines are ignored.
 *
 * What is preserved: every record, field and token in order with its line and column; whole-line comments above a
 * header or a field, attached to it; a comment after a header or a value on the same line; whole-line comments after
 * the last field of a section.  Not preserved: blank lines, indentation, the blanks between tokens and the column of
 * a comment -- WriteRecords writes them in the canonical form, so the bytes of a hand-written section change but its
 * records do not: Read(Write(Read(x))) == Read(x) up to line and column numbers, and Write(Read(Write(r))) ==
 * Write(r) byte for byte.  A section whose text must stay byte-identical is edited through Document::ReplaceLines
 * instead (the lines of one field), as pyazr edits sections today.
 *
 * Problems are returned, all of them, never thrown or printed; a line with a problem is still read as far as it can
 * be (an unterminated quote ends at the end of the line), so the validator sees the rest of the record.
 */
struct RecordReadResult {
  RecordSet set;
  std::vector<RecordDiagnostic> diagnostics;
};

/// Read `lines` (without line endings) as records; `firstLine` is the file line number of lines[0] (1-based).
RecordReadResult ReadRecords(const std::vector<std::string> &lines, int firstLine, const std::vector<std::string> &kinds);

/// Read the body of a section of `doc` (the lines between its tags).
RecordReadResult ReadRecords(const Document &doc, const Document::Section &section, const std::vector<std::string> &kinds);

/// The canonical text of `set`, one line per element (no endings), ready for Document::ReplaceLines: records separated
/// by one blank line, fields indented by four spaces, the values of a record in one column (column 20, or one blank
/// after the record's longest key when that key has 15 characters or more), comments kept where they were attached, tokens quoted when they were or when they hold a blank, '#' or
/// nothing.  A token holding '"' cannot be written so that it reads back the same (the reader has no escapes).
std::vector<std::string> WriteRecords(const RecordSet &set);

/// The value tokens of one field as WriteRecords writes them: separated by one blank, quoted as said above.
std::string WriteTokens(const std::vector<Token> &tokens);

/// Split one line into tokens and a trailing comment (the lexical rules above), for callers that read one field.
/// Appends problems (E_UNTERMINATED_QUOTE, E_QUOTE_IN_TOKEN) with `line` as their line number.
std::vector<Token> Tokenize(const std::string &line, int lineNumber, std::string &comment,
                            std::vector<RecordDiagnostic> &problems);

}  // namespace io
}  // namespace azure2

#endif
