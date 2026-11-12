#ifndef AZURE2_IO_PROJECTFORMAT_H
#define AZURE2_IO_PROJECTFORMAT_H

#include <string>

#include "azure2/io/Document.h"

namespace azure2 {
namespace io {

/// The format of a project: 2 or 3 with a "format=2" or "format=3" line in <config>, else 1
/// (docs/upgrade/project.md section 4, docs/upgrade/project_format3.md section 8.4).
int FormatOf(const Document &doc);

/// `key=value` of a <config> line (a comment after '#' dropped, both sides trimmed); false if the line has no '='.
bool ConfigKeyValue(const std::string &line, std::string &key, std::string &value);

/*!
 * The refusal of docs/upgrade/project_format3.md section 8.4 for a text that is not a format-3 project ("" for a
 * format-3 text): the format found and the first old section with its line, the conversion command (CLI and pyazr),
 * and that the converter keeps the file; a format-3 header over an old section, and a text without <config> (empty,
 * or not a project), have their own lines.  The one text of every refusal of an old file: the CLI's runs, --validate,
 * --describe, --print-predictions, the API session (Config::ReadConfigFile) and pyazr (binding not_format3_message).
 * The command line exits with 3 after it.
 */
std::string NotFormat3Message(const std::string &text, const std::string &fileName);

}  // namespace io
}  // namespace azure2

#endif
