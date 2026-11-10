#ifndef AZURE2_IO_CONVERT_H
#define AZURE2_IO_CONVERT_H

#include <string>
#include <vector>

#include "azure2/io/Document.h"

namespace azure2 {
namespace io {

/// The format of a project: 2 or 3 with a "format=2" or "format=3" line in <config>, else 1
/// (docs/upgrade/project.md section 4, docs/upgrade/project_format3.md section 8.4).
int FormatOf(const Document &doc);

/*!
 * Converts a format-1 project to format 2 in place (docs/upgrade/project.md section 4) and returns what it changed,
 * one line per change (empty for a file already in format 2 or 3):
 *  - <lastRun> removed (the Qt GUI's runtime options, read by nothing);
 *  - <mcmc>'s run settings removed (the CLI asks for them), and <mcmc> with them when its <parameters> are empty;
 *  - useAdaptiveGrid moved from <potential> to <config>, and <potential> removed unless useHybridPotential=1 for
 *    the default or for any pair;
 *  - "format=2" added as the last line of <config>.
 * Every other byte is kept.  Converting twice is converting once.
 */
std::vector<std::string> ConvertToFormat2(Document &doc);

}  // namespace io
}  // namespace azure2

#endif
