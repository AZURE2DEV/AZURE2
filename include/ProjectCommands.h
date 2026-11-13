#ifndef PROJECTCOMMANDS_H
#define PROJECTCOMMANDS_H

#include "CommandLineValues.h"

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

#endif
