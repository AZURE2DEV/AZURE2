#ifndef RUNPREPARATION_H
#define RUNPREPARATION_H

// A run request turned into the run's Config (include/RunRequest.h).

#include <sstream>
#include <string>

#include "RunRequest.h"

class Config;

/// The section <name> of the project as a stream over its text (Config::SectionStream: io::Document, whole-line
/// tags), or a stream in the fail state when there is no such section.
std::istringstream SectionText(const Config &configure, const std::string &name);

/// The refusal of a file the request names that cannot be read (relative to the project's directory), in three
/// parts; empty when every named file can be read.  Once the project is read.
std::string RequestFileRefusal(const Config &configure, const RunRequest &request);

/*!
 * Sets `configure` for `request`: the command's calculation, the band, MINOS's chi-squared rise, the parameter file,
 * the rate's pairs and temperatures, the external capture (whether the run has it, from the project; the amplitude
 * file, with a notice when the run has none).  A project whose levels cannot be read stops the run (exit 1).
 */
void PrepareRun(Config &configure, RunRequest &request);

#endif
