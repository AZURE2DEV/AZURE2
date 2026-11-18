#ifndef AZURE2_COMMANDLINEVALUES_H
#define AZURE2_COMMANDLINEVALUES_H

#include <ostream>
#include <string>
#include <vector>

namespace azure2 {
namespace support {
class Environment;
}
}  // namespace azure2

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
/// The refinement factor k of the sub-point grids (Config::refineQuadrature) from the process's record:
/// --refine-quadrature=K or AZURE_REFINE_QUADRATURE, 1 when neither is given, else a whole number from 1 to 64.
int RefineQuadratureFromEnvironment();

/// The largest N of --threads=N.
constexpr long kMaxThreads = 4096;
/// The OpenMP threads of the process and where the count came from: "--threads", "--threads, over
/// OMP_NUM_THREADS=M", "OMP_NUM_THREADS", "the OpenMP runtime's default".
struct ThreadsInForce {
  int count = 1;
  std::string source;
};
/// --threads=N (the last one given; 1 <= N <= kMaxThreads, else std::invalid_argument): sets the OpenMP thread count
/// of the process.  Called by main before any parallel region; without the option, OMP_NUM_THREADS or the runtime's
/// default decide, as they did.  Returns the count in force.
ThreadsInForce ApplyThreadsOption(const std::vector<std::string> &arguments,
                                  const azure2::support::Environment &environment);
/// {"value": N, "source": "..."} for --describe --json.
std::string ThreadsJson(const ThreadsInForce &threads);
/// The CLI's refusal of a flag or variable: "ERROR: <message>" and the usage line on standard output; returns 2, the
/// exit status.
int OptionError(const std::string &message);
/// The same on `out` (the run's stream of a front end in its own process, ProgramMain.h).
int OptionError(std::ostream &out, const std::string &message);

}  // namespace azure2

#endif
