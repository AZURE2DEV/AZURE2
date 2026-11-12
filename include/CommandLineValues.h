#ifndef AZURE2_COMMANDLINEVALUES_H
#define AZURE2_COMMANDLINEVALUES_H

#include <string>

/*!
 * The values of the command line's value-taking flags (parseOptions in src/AZURE2.cpp) and of the diagnostic
 * environment variable AZURE_REFINE_QUADRATURE.  A flag given without its value, or with a value that does not parse,
 * is an error (std::invalid_argument, its message naming the flag and the form it needs), never a warning and a
 * default: a flag ignored silently is a trap for a person and for an agent.  The CLI prints the message with the
 * usage line and exits 2 (OptionError); pyazr raises it.
 */
namespace azure2 {

/// Whether `option` is the flag `flag` (given as "--name="): "--name" alone or "--name=VALUE".
bool OptionIs(const std::string &option, const std::string &flag);
/// The whole number N of `option` = "<flag>N" (`flag` given as "--name="), low <= N <= high.
long OptionWholeNumber(const std::string &option, const std::string &flag, long low, long high);
/// The text of `option` = "<flag>TEXT" (`flag` given as "--name="), not empty.
std::string OptionText(const std::string &option, const std::string &flag);
/// The refinement factor k of the sub-point grids from AZURE_REFINE_QUADRATURE (Config::refineQuadrature): 1 when the
/// variable is unset, else a whole number from 1 to 64.
int RefineQuadratureFromEnvironment();
/// The CLI's refusal of a flag or variable: "ERROR: <message>" and the usage line on standard output; returns 2, the
/// exit status.
int OptionError(const std::string &message);

}  // namespace azure2

#endif
