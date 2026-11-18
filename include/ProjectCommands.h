#ifndef PROJECTCOMMANDS_H
#define PROJECTCOMMANDS_H

#include <string>
#include <vector>

#include "CommandLineValues.h"
#include "azure2/experiment/Diagnostics.h"

/*!
 * The project-file commands of AZURE2's command line, handled before the options of a calculation are parsed
 * (docs/source/reference/command_line.rst):
 *   --convert OLD NEW [--to 3 [--json-report R] [--answers A] [--exact-grids] [--strict] [--allow-questions]]
 *   --validate FILE [--json]
 *   --schema
 * Returns true when argv holds one of them; `exitCode` is then the code the program exits with.  `threads`: the
 * threads in force (--describe prints them).
 */
bool RunProjectCommand(int argc, char **argv, const azure2::ThreadsInForce &threads, int &exitCode);

/// What --validate finds in the bytes of a project (the same for the program and a front end, docs/upgrade/qt_gui.md):
/// the refusal of a file that is not format 3 (then nothing else is set), or every diagnostic of the reader, the
/// validator, the build checks and the targets' stopping, the stopping's report, the count of build checks skipped.
struct ProjectValidation {
  std::string notFormat3;
  std::vector<azure2::experiment::Diagnostic> diagnostics;
  std::string stoppingText;
  int checksSkipped = 0;
  bool ok = false;  ///< no error
};
/// `file` names the project (messages, and the directory its data files are relative to); `bytes` its text.
ProjectValidation ValidateProjectBytes(const std::string &bytes, const std::string &file);

#endif
