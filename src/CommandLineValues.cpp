// The values of the command line's value-taking flags, --threads and AZURE_REFINE_QUADRATURE (include/CommandLineValues.h).
#include "CommandLineValues.h"

#include "azure2/support/Environment.h"

#include <iostream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace azure2 {

namespace {

/// The flag's name without its "=" ("--max-calls" of "--max-calls=").
std::string NameOf(const std::string &flag) { return flag.substr(0, flag.size() - 1); }

/// The value of `option` = "<flag>VALUE"; an error when the flag is given without "=".
std::string ValueOf(const std::string &option, const std::string &flag, const std::string &form) {
  if (option.compare(0, flag.size(), flag) != 0)
    throw std::invalid_argument(NameOf(flag) + " needs a value: " + form + " (got '" + option + "')");
  return option.substr(flag.size());
}

/// The whole number of `text` with low <= N <= high, or an error naming `what` (the one checked parse, the switches').
long WholeNumber(const std::string &text, long low, long high, const std::string &what) {
  return azure2::support::CheckedWholeNumber(text, low, high, what);
}

}  // namespace

bool OptionIs(const std::string &option, const std::string &flag) {
  return option == NameOf(flag) || option.compare(0, flag.size(), flag) == 0;
}

long OptionWholeNumber(const std::string &option, const std::string &flag, long low, long high) {
  return WholeNumber(ValueOf(option, flag, flag + "N"), low, high, NameOf(flag));
}

std::string OptionText(const std::string &option, const std::string &flag) {
  const std::string text = ValueOf(option, flag, flag + "FILE");
  if (text.empty()) throw std::invalid_argument(NameOf(flag) + "= needs a value: " + flag + "FILE");
  return text;
}

int RefineQuadratureFromEnvironment() {
  // A diagnostic that silently did not refine would pass a gate it never ran: a bad value is an error, whether the
  // option (--refine-quadrature=K) or the variable gave it (the row's range, azure2/support/Environment.h).
  return static_cast<int>(azure2::support::SwitchWholeNumber(azure2::support::ProcessEnvironment(),
                                                             azure2::support::Switch::kRefineQuadrature, 1));
}

ThreadsInForce ApplyThreadsOption(const std::vector<std::string> &arguments,
                                  const azure2::support::Environment &environment) {
  ThreadsInForce threads;
  const std::string &variable = environment.OmpNumThreads();
  std::string given;
  for (const std::string &argument : arguments)
    if (OptionIs(argument, "--threads=")) given = argument;  // the last one given
  if (!given.empty()) {
    const long n = OptionWholeNumber(given, "--threads=", 1, kMaxThreads);
#ifdef _OPENMP
    omp_set_num_threads(static_cast<int>(n));
#endif
    threads.source = variable.empty() ? "--threads" : "--threads, over OMP_NUM_THREADS=" + variable;
  } else {
    threads.source = variable.empty() ? "the OpenMP runtime's default" : "OMP_NUM_THREADS";
  }
#ifdef _OPENMP
  threads.count = omp_get_max_threads();
#else
  threads.count = 1;
  if (!given.empty()) threads.source += " (a build without OpenMP: one thread)";
#endif
  return threads;
}

std::string ThreadsJson(const ThreadsInForce &threads) {
  return "{\"value\": " + std::to_string(threads.count) + ", \"source\": \"" + threads.source + "\"}";
}

int OptionError(const std::string &message) {
  std::cout << "ERROR: " << message << std::endl
            << "\tSyntax: AZURE2 <command> [options] FILE  (AZURE2 --help lists the commands and options)" << std::endl;
  return 2;
}

}  // namespace azure2
