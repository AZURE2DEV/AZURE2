#ifndef PROGRAMMAIN_H
#define PROGRAMMAIN_H

// The program's run of one command on a project (src/ProgramMain.cpp; docs/upgrade/qt_gui.md section 3).

#include <string>
#include <vector>

#include "CommandLineValues.h"

class Config;

/*!
 * What the program does after its process's start-up (the environment record, the threads, GSL's error handler and
 * the project-file commands, all in main): the command and options of `arguments` (argv without the program's name,
 * as `AZURE2 <command> [options] FILE`), the project read, the welcome lines, the run, the closing line, every line on
 * configure.outStream.  `configure` is fresh (Config(stream)); its `observer` may be set (ProgramObserver.h).
 * `threads` is printed ("Threads: N (source)").  Returns main's exit value: 0, 1 (aborted or refused project), 2 (a
 * refused command line), 3 (not a format-3 project), -1 (an unreadable project or input file, a failed set-up).
 */
int RunProgram(Config &configure, const std::vector<std::string> &arguments, const azure2::ThreadsInForce &threads);

#endif
