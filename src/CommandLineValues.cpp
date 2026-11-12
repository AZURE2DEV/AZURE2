// The values of the command line's value-taking flags and of AZURE_REFINE_QUADRATURE (include/CommandLineValues.h).
#include "CommandLineValues.h"

#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

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

/// The whole number of `text` with low <= N <= high, or an error naming `what`.
long WholeNumber(const std::string &text, long low, long high, const std::string &what) {
  char *end = nullptr;
  errno = 0;
  const long n = text.empty() ? 0 : std::strtol(text.c_str(), &end, 10);
  if (text.empty() || *end != '\0' || errno == ERANGE || n < low || n > high)
    throw std::invalid_argument(what + " needs a whole number from " + std::to_string(low) + " to " +
                                std::to_string(high) + ", not '" + text + "'");
  return n;
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
  // A diagnostic that silently did not refine would pass a gate it never ran: a bad value is an error.
  const char *v = std::getenv("AZURE_REFINE_QUADRATURE");
  if (!v) return 1;
  return static_cast<int>(WholeNumber(v, 1, 64, "AZURE_REFINE_QUADRATURE"));
}

int OptionError(const std::string &message) {
  std::cout << "ERROR: " << message << std::endl
            << "\tSyntax: AZURE2 <options> configfile  (AZURE2 --help lists the options)" << std::endl;
  return 2;
}

}  // namespace azure2
