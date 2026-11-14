#ifndef RUNPREPARATION_H
#define RUNPREPARATION_H

// A run request turned into the run's Config (include/RunRequest.h): the one path of both readers.

#include <sstream>
#include <string>

#include "RunRequest.h"

class Config;
class RunMenu;

/// The section <name> of the project as a stream over its text (Config::SectionStream: io::Document, whole-line
/// tags), or a stream in the fail state when there is no such section.
std::istringstream SectionText(const Config &configure, const std::string &name);

/// The refusal of a file the request names that cannot be read (relative to the project's directory), in three
/// parts; empty when every named file can be read.  For a request from the command line, once the project is read.
std::string RequestFileRefusal(const Config &configure, const RunRequest &request);

/*!
 * Sets `configure` for `request`: the command's calculation, the band, MINOS's chi-squared rise, the parameter file,
 * the rate's pairs and temperatures, the external capture (whether the run has it, from the project; the amplitude
 * file).  The questions are walked in the table's order; with a `menu`, each question the command's menu asks is
 * asked as the walk reaches it (the external capture only once the project says the run has it), so the menu's
 * prompts and the project's messages come in the menu's order; without one, the request is complete (the command
 * line).  A project whose levels cannot be read stops the run (exit 1).
 */
void PrepareRun(Config &configure, RunRequest &request, RunMenu *menu);

#endif
